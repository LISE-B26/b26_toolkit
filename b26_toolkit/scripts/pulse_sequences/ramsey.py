"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    b26_toolkit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    b26_toolkit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with b26_toolkit.  If not, see <http://www.gnu.org/licenses/>.
"""

from b26_toolkit.scripts.pulse_sequences.rabi import Rabi
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_base_script import PulsedExperimentBaseScript
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, MicrowaveGenerator, Pulse
from b26_toolkit.plotting.plots_1d import plot_pulses, update_pulse_plot, plot_1d_simple_timetrace_ns, update_1d_simple, plot_psd
from pylabcontrol.core import Parameter
from b26_toolkit.data_processing.fit_functions import fit_rabi_decay, cose_with_decay
import numpy as np
import scipy.signal as spig
from copy import deepcopy
import itertools

MAX_AVERAGES_PER_SCAN = 1000000  # 1E6, the max number of loops per point allowed at one time (true max is ~4E6 since
                                 #pulseblaster stores this value in 22 bits in its register
                                # DS 20191216: changed from 1e5 to 1e6 since loop register is 20 bits. PB will throw error
                                # if too large

MAX_NUM_PULSES = 4096


class Ramsey(Rabi):
    """
    This script executes a Ramsey measurement. It uses a single_init scheme, i.e. there's only one laser pulse, where
    we collect the signal fluor. at the beginning and ref fluor. at the end.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('3pi_half_pulse_time', 75, float, 'time duration of a 3pi/2 pulse (in ns)')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'minimum time between pi pulses'),
            Parameter('max_time', 10000, float, 'maximum time between pi pulses'),
            Parameter('time_step', 5, float, 'time step increment of time between pi pulses (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 1000, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}


    def _create_pulse_sequences(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """
        pulse_sequences = []

        tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],self.settings['tau_times']['time_step'])
        tau_list = np.ndarray.tolist(tau_list) # 20180731 ER convert to list

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        # MM: updated to min_pulse_dur
        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        wait_time = 0

        for tau in tau_list:

            pulse_sequence = \
                [
                    Pulse(microwave_channel, laser_off_time + wait_time, pi_half_time),
                    Pulse(microwave_channel, laser_off_time + wait_time + pi_half_time + tau, pi_half_time),
                ]

            end_of_first_ramsey = laser_off_time + wait_time + tau + 2 * pi_half_time

            pulse_sequence += [
                Pulse('laser', end_of_first_ramsey + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', end_of_first_ramsey + delay_mw_readout + delay_readout, meas_time),
            ]

            start_of_second_ramsey = end_of_first_ramsey + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence += \
                [
                    Pulse(microwave_channel, start_of_second_ramsey + wait_time, pi_half_time),
                    Pulse(microwave_channel, start_of_second_ramsey + wait_time + pi_half_time + tau, three_pi_half_time)
                ]

            end_of_second_ramsey = start_of_second_ramsey + wait_time + pi_half_time + tau + three_pi_half_time

            pulse_sequence += [
                Pulse('laser', end_of_second_ramsey + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', end_of_second_ramsey + delay_mw_readout + delay_readout, meas_time)
            ]
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _plot(self, axislist, data = None):
        '''
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional) dataset to plot (dictionary that contains keys counts, tau, fits), if not provided use self.data
        '''

        if data is None:
            data = self.data

        if 1==2 and 'fits' in data.keys() and data['fits'] is not None:
            counts = data['counts'][:,1]/ data['counts'][:,0]
            tau = data['tau']
            fits = data['fits'] # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            #axislist[0].hold(True) ER 20181012

            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            #pi_time = 2*np.pi / fits[1] / 2
            pi_time = (np.pi - fits[2])/fits[1]
            pi_half_time = (np.pi/2 - fits[2])/fits[1]
            three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
            rabi_freq = 1000*fits[1]/(2*np.pi)
         #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            super(Rabi, self)._plot(axislist)
            axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9))
            axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)


class RamseyGarbageDoNotUse(PulsedExperimentBaseScript): # DS 20191122
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns)')
        ]),
        Parameter('tau', 500, float, 'free precession time (in ns)'),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 1000, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('wait_time', 10, float, 'wait time in ns'),
        Parameter('num_averages', 100000, int, 'number of averages'),
        Parameter('num_points', 100, int, 'number of points after average')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _create_pulse_sequences(self):
        '''

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''
        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        wait_time = self.settings['wait_time']
        tau = self.settings['tau']

        pulse_sequence = \
        [
            Pulse(microwave_channel, laser_off_time + wait_time, pi_half_time),
            Pulse(microwave_channel, laser_off_time + wait_time + pi_half_time + tau, pi_half_time),
        ]

        end_of_first_ramsey = laser_off_time + wait_time + tau + 2 * pi_half_time

        pulse_sequence += [
             Pulse('laser', end_of_first_ramsey + delay_mw_readout, nv_reset_time),
             Pulse('apd_readout', end_of_first_ramsey + delay_mw_readout + delay_readout, meas_time),
        ]

        start_of_second_ramsey = end_of_first_ramsey + delay_mw_readout + nv_reset_time + laser_off_time

        pulse_sequence += \
        [
            Pulse(microwave_channel, start_of_second_ramsey + wait_time, pi_half_time),
            Pulse(microwave_channel, start_of_second_ramsey + wait_time + pi_half_time + tau, three_pi_half_time)
        ]

        end_of_second_ramsey = start_of_second_ramsey + wait_time + pi_half_time + tau + three_pi_half_time

        pulse_sequence += [
            Pulse('laser', end_of_second_ramsey + delay_mw_readout, nv_reset_time),
            Pulse('apd_readout', end_of_second_ramsey + delay_mw_readout + delay_readout, meas_time)
        ]

        return [pulse_sequence], [tau], meas_time

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True}) # ER 20181018
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        ####################################SUPER _function###############################################
            # Set DAQ
        if self.settings['daq_type'] == 'PCI':
            self._daq = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self._daq = self.instruments['NI9402']['instance']

        # make sure the microwave_switch is turned off so that we don't burn any steel cables. ER 20181017
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

        # retrieve initial mw carrier frequency to protect against bad fits in NV ESR tracking
        last_mw = self.scripts['esr'].instruments['microwave_generator']['instance'].frequency
        pulse_ampl = self.scripts['esr'].instruments['microwave_generator']['instance'].amplitude

        # ER 20181214 retrieve modulation on or off for main experiment
        mod_flag = self.scripts['esr'].instruments['microwave_generator']['instance'].enable_modulation

        # Keeps track of index of current pulse sequence for plotting
        self.sequence_index = 0

        # self.is_valid and create pulses
        self.pulse_sequence, self.tau, self.measurement_gate_width = self.create_pulse_sequences()
        self.pulse_sequence = self.pulse_sequence[0]
        self.tau = self.tau[0]

        self.num_averages = self.settings['num_averages']
        self.num_points = self.settings['num_points']
        self.num_measurements = self.num_averages * self.num_points
        self.num_daq_reads_per_seq = self._count_daq_reads()

        self.sequence_duration = np.max([p.end_time for p in self.pulse_sequence])

        if self.num_measurements > MAX_AVERAGES_PER_SCAN:
            num_cmd_per_seq = len(self.pulse_sequence)
            max_num_seq_repeats = MAX_NUM_PULSES // num_cmd_per_seq

            if self.num_measurements > max_num_seq_repeats * MAX_AVERAGES_PER_SCAN:
                self.log('Error: Number of measurements exceeds hardware limits of {} measurements for this sequence.'.format(max_num_seq_repeats * MAX_AVERAGES_PER_SCAN))
                return

            # repeat sequences
            num_seq_repeats = self.num_measurements // MAX_AVERAGES_PER_SCAN
            self.num_measurements = num_seq_repeats * MAX_AVERAGES_PER_SCAN
            self.num_averages = self.num_measurements // self.num_points
            self._repeat_pulse_sequences(num_seq_repeats)

        print('Pulse sequence: ', self.pulse_sequence)
        print('Pulse sequence length: ', len(self.pulse_sequence))
        print('Num averages: ', self.num_averages)
        print('Num points: ', self.num_points)
        print('Num measurements: ', self.num_measurements)

        self.time_step = self.num_averages * self.sequence_duration
        # calculates the number of daq reads per loop requested in the pulse sequence by asking how many apd reads are
        # called for. if this is not calculated properly, daq will either end too early (number too low) or hang since it
        # never receives the rest of the counts (number too high)
        num_daq_reads = self._count_daq_reads()

        if self._abort:
            self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
            self.log('aborted pulseblaster script during loop')



        self.count_data = np.array(self._run_single_sequence(self.pulse_sequence,
                                                    np.minimum(MAX_AVERAGES_PER_SCAN, self.num_measurements),
                                                    num_daq_reads))
        self.data['counts'] = self._normalize_to_kCounts(self.count_data,
                                                         self.measurement_gate_width,
                                                         self.num_averages)
        self.data['time_step'] = self.time_step

        if (len(self.data['counts'][0]) == 1) and not self._abort:
            self.data['counts'] = np.array([item for sublist in self.data['counts'] for item in sublist])

        # save data on the fly so that we can start to analyze it while the experiment is running!
        if self.settings['save']:
            #     self.save_b26()
            self.save_data()

        ##################################################################################################

    def _count_daq_reads(self):
        num_daq_reads = 0
        for pulse in self.pulse_sequence:
            if pulse.channel_id == 'apd_readout':
                num_daq_reads += 1
        return num_daq_reads

    def _repeat_pulse_sequences(self, num_seq_repeats):
        print('repeated')
        repeated_sequence = deepcopy(self.pulse_sequence)
        for j in range(1, num_seq_repeats):
            sequence = deepcopy(self.pulse_sequence)
            for k in range(len(sequence)):
                sequence[k].start_time += self.sequence_duration * j
                sequence[k].end_time += self.sequence_duration * j
            repeated_sequence.extend(sequence)
        self.pulse_sequence = repeated_sequence

    def _sum_measurements(self, daq_results, num_daq_reads):
        summed_results = []
        for i in range(self.num_daq_reads_per_seq):
            daqslice = list(itertools.islice(daq_results, i, None, self.num_daq_reads_per_seq))
            groups = [sum(daqslice[j:(j + self.num_averages)]) for j in range(0, len(daqslice), self.num_averages)]
            summed_results.append(groups)
        return summed_results

    def _plot(self, axes_list, data=None):
            """
            Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
            received for each time
            Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
            performed

            Args:
                axes_list: list of axes to write plots to (uses first 2)
                data (optional) dataset to plot (dictionary that contains keys counts, tau), if not provided use self.data
            """

            #        if self.scripts['find_nv'].is_running: #self.scripts['find_nv'].data['maximum_point']:
            #            # self.scripts['find_nv'].plot(axes_list)
            #             self.scripts['find_nv']._plot(axes_list)
            #             self._plot_refresh = True
            #        else:
            if data is None:
                data = self.data

            if 'counts' in data.keys():
                # The following does not work for pulsedelays; you need to comment out the 'if' for it to work.
                # if counts != []:
                #     plot_1d_simple_timetrace_ns(axes_list[0], data['tau'], [data['counts'])
                #  print('plotting with tau values: ', data['tau'])
                counts = (data['counts'][0, :] - data['counts'][1, :]) / (data['counts'][0, :] + data['counts'][1, :])

                time = np.arange(0, self.num_points * self.time_step, self.time_step)
                plot_1d_simple_timetrace_ns(axes_list[0], time, [counts])
                freq, psd = spig.periodogram(counts, 1. / self.time_step)
                plot_psd(freq, psd, axes_list[1])
