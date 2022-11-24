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

import numpy as np
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_base_script import PulsedExperimentBaseScript
from b26_toolkit.instruments import NI6259, NI9402, NI6229, B22PulseBlaster, MicrowaveGenerator, Pulse
from b26_toolkit.plotting.plots_1d import plot_pulses, update_pulse_plot, plot_1d_simple_timetrace_ns, update_1d_simple
from pylabcontrol.core import Parameter
from b26_toolkit.data_processing.fit_functions import fit_rabi_decay, cose_with_decay


class Rabi(PulsedExperimentBaseScript):  # ER 5.25.2017
    """
This script applies a microwave pulse at fixed power for varying durations to measure Rabi oscillations.
Uses a double_init scheme
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 15, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 200, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 5., [2.5, 4., 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
    ]

    _INSTRUMENTS = {'NI6229': NI6229, 'NI6259': NI6259, 'NI9402': NI9402, 'PB': B22PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(Rabi, self)._function()

        if 'counts' in self.data.keys() and 'tau' in self.data.keys():
            counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
            print('shape')
            print(np.shape(self.count_data))
            tau = self.data['tau']

            try:
                fits = fit_rabi_decay(tau, counts, variable_phase=True)
                self.data['fits'] = fits
            except:
                self.data['fits'] = None
                self.log('rabi fit failed')

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
        # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
        #                  self.settings['tau_times']['time_step'])
        # JG 16-08-25 changed (15ns min spacing is taken care of later):
       # tau_list = list(range(int(self.settings['tau_times']['min_time']),
                            #  int(self.settings['tau_times']['max_time']),
                            #  self.settings['tau_times']['time_step']))

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)

        #MM: update 15 to min_pulse_duration
        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]


        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = [Pulse('laser', laser_off_time + tau + 2*40, nv_reset_time),
                              Pulse('apd_readout', laser_off_time + tau + 2*40 + delay_readout, meas_time)]

            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence.append(Pulse(microwave_channel,
                                            laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time,
                                            tau))

            pulse_sequence.append(Pulse('laser',
                                        laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time + tau + 2*40 + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time + tau + 2*40 + delay_mw_readout + delay_readout,
                                        meas_time))
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
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

        if 'fits' in data.keys() and data['fits'] is not None:
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

class RabiDoublePi(Rabi):
    """
    Runs Rabi, but instead of sweeping the duration of a single MW pulse, sweeps the duration of two back to back MW
    pulses with same durations. Used to calibrate pulses such that two pi pulses bring the Bloch vector back to the
    original state
    """

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
        # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
        #                  self.settings['tau_times']['time_step'])
        # JG 16-08-25 changed (15ns min spacing is taken care of later):
       # tau_list = list(range(int(self.settings['tau_times']['min_time']),
                            #  int(self.settings['tau_times']['max_time']),
                            #  self.settings['tau_times']['time_step']))

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)

        #MM: update 15 to min_pulse_duration
        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]


        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = [Pulse('laser', laser_off_time + tau + 2*40, nv_reset_time),
                              Pulse('apd_readout', laser_off_time + tau + 2*40 + delay_readout, meas_time)]

            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence.append(Pulse(microwave_channel,
                                            laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time,
                                            tau))
                pulse_sequence.append(Pulse(microwave_channel,
                                            laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2*40,
                                            tau))

            pulse_sequence.append(Pulse('laser',
                                        laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time + tau + 2*40 + tau + 2*40 + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time + tau + 2*40 + tau + 2*40 + delay_mw_readout + delay_readout,
                                        meas_time))
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

class RabiPowerSweepSingleTau(PulsedExperimentBaseScript):
    """
This script applies a microwave pulse at fixed power for varying durations to measure Rabi Oscillations
todo(emma): (write as a double_init scheme)
    """
    _DEFAULT_SETTINGS = [
        Parameter('min_mw_power', -45.0, float, 'minimum microwave power in dB'),
        Parameter('max_mw_power', -45.0, float, 'maximum microwave power in dB'),
        Parameter('mw_power_step', 1.0, float, 'power to step by in dB'),
        Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
        Parameter('mw_time', 200, float, 'total time of rabi oscillations (in ns)'),
        Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('reset_time', 10000, int, 'time with laser on at the beginning to reset state'),
    ]

    _INSTRUMENTS = {'daq': NI6229, 'PB': B22PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        # COMMENT_ME
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_frequency']})
        mw_power_values = np.arange(self.settings['min_mw_power'],
                                    self.settings['max_mw_power'] + self.settings['mw_power_step'],
                                    self.settings['mw_power_step'])

     #   print(mw_power_values)
        self.data = {'mw_power_values': mw_power_values, 'counts_for_mw': np.zeros(len(mw_power_values))}
        for index, power in enumerate(mw_power_values):
            self.instruments['mw_gen']['instance'].update({'amplitude': float(power)})
            super(RabiPowerSweepSingleTau, self)._function(self.data)
            self.data['counts_for_mw'][index] = self.data['counts'][0]

    def _create_pulse_sequences(self):
        '''

        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''
        pulse_sequences = []
        reset_time = self.settings['reset_time']
        mw_time = self.settings['mw_time']
        pulse_sequences.append([Pulse('laser', 0, reset_time),
                                Pulse('microwave_i', reset_time + 200, mw_time),
                                Pulse('laser', reset_time + mw_time + 300, self.settings['meas_time']),
                                Pulse('apd_readout', reset_time + mw_time + 300, self.settings['meas_time'])
                                ])

        end_time_max = 0
        for pulse_sequence in pulse_sequences:
            for pulse in pulse_sequence:
                end_time_max = max(end_time_max, pulse.start_time + pulse.duration)
        for pulse_sequence in pulse_sequences:
            pulse_sequence.append(
                Pulse('laser', end_time_max + 1850, 15))  # Jan Feb 1st 2017: what is 1850??? Need to comment!

        return pulse_sequences, [mw_time], self.settings['meas_time']

    def _plot(self, axes_list, data=None):
        '''
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional) dataset to plot (dictionary that contains keys counts_for_mw, mw_power_values), if not provided use self.data
        '''
        if data is None:
            data = self.data

        counts = data['counts_for_mw']
        x_data = data['mw_power_values']
        axis1 = axes_list[0]
        if not counts == []:
            plot_1d_simple_timetrace_ns(axis1, x_data, [counts], x_label='microwave power (dBm)')
        axis2 = axes_list[1]
        plot_pulses(axis2, self.pulse_sequences[self.sequence_index])

    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''
        counts = self.data['counts_for_mw']
        x_data = self.data['mw_power_values']
        axis1 = axes_list[0]
        if not counts == []:
            update_1d_simple(axis1, x_data, [counts])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])