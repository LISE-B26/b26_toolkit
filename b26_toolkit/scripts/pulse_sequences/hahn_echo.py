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
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_generic import PulsedExperimentGeneric
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, MicrowaveGenerator, Pulse, Commander
from pylabcontrol.core import Parameter, Script
from b26_toolkit.data_processing.fit_functions import fit_exp_decay, exp_offset
from b26_toolkit.scripts import FindNvPulsed, AutoFocusDAQPulsed, Esr


class HahnEcho(PulsedExperimentGeneric):
    """
    This script runs a Hahn echo on the NV to find the Hahn echo T2. To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', '+i', ['+i', '-i', '+q', '-q'], 'Channel to use for mw pulses'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            # We keep this but this is for legacy.
            Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns)')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'minimum time between pi pulses'),
            Parameter('max_time', 10000, float, 'maximum time between pi pulses'),
            Parameter('time_step', 5, float,
                  'time step increment of time between pi pulses (in ns)')
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

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'commander': Commander}
    from b26_toolkit.scripts.find_nv import FindNvStrobe
    _SCRIPTS = {'find_nv': FindNvStrobe}

    def __init__(self, instruments, scripts, name=None, settings=None, log_function=None, data_path=None):
        """
        Standard script initialization
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_generic import PulsedExperimentGeneric
        self._DEFAULT_SETTINGS += PulsedExperimentGeneric._DEFAULT_SETTINGS
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.ref_index = 0

    def _configure_instruments_start_of_script(self):
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

    def _function(self):
        self.data['fits'] = None
        super(HahnEcho, self)._function(self.data)

        counts = (- self.data['counts'][:, 1] + self.data['counts'][:,0]) / (self.data['counts'][:,1] + self.data['counts'][:, 0])
        tau = self.data['tau']

        try:
            fits = fit_exp_decay(tau, counts, offset = True, verbose = True)
            self.data['fits'] = fits
        except:
            self.data['fits'] = None
            self.log('t2 fit failed')

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

        if microwave_channel[-2] == '+':
            negative_channel = microwave_channel[:-2] + '-' + microwave_channel[-1]
        elif microwave_channel[-2] == '-':
            negative_channel = microwave_channel[:-2] + '+' + microwave_channel[-1]

        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = \
            [
                Pulse(microwave_channel, laser_off_time, pi_half_time),
                Pulse(microwave_channel, laser_off_time + pi_half_time + tau, pi_time),
                Pulse(microwave_channel, laser_off_time + pi_half_time + tau + pi_time + tau, pi_half_time)
            ]

            end_of_first_HE = laser_off_time + pi_half_time + tau + pi_time + tau + pi_half_time

            pulse_sequence += [
                 Pulse('laser', end_of_first_HE + delay_mw_readout, nv_reset_time),
                 Pulse('apd_readout', end_of_first_HE + delay_mw_readout + delay_readout, meas_time),
                 ]

            start_of_second_HE = end_of_first_HE + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence += \
            [
                Pulse(microwave_channel, start_of_second_HE, pi_half_time),
                Pulse(microwave_channel, start_of_second_HE + pi_half_time + tau, pi_time),
                Pulse(negative_channel, start_of_second_HE + pi_half_time + tau + pi_time + tau, pi_half_time)
            ]

            #end_of_second_HE = start_of_second_HE + pi_half_time / 2. + tau + tau - pi_half_time / 2. + pi_half_time
            end_of_second_HE = start_of_second_HE + pi_half_time + tau + pi_time + tau + pi_half_time

            pulse_sequence += [
                Pulse('laser', end_of_second_HE + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', end_of_second_HE + delay_mw_readout + delay_readout, meas_time)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _create_pulse_sequences_old(self):
        '''
        The old sequence where tau is counted from the middle of the pulse
        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''
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

        for tau in tau_list:
            pulse_sequence = [
                Pulse(microwave_channel, laser_off_time, pi_half_time),
                Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau - pi_time/2., pi_time),
                Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2., pi_half_time)
            ]

            end_of_first_HE = laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time

            pulse_sequence += [
                 Pulse('laser', end_of_first_HE + delay_mw_readout, nv_reset_time),
                 Pulse('apd_readout', end_of_first_HE + delay_mw_readout + delay_readout, meas_time),
                 ]

            start_of_second_HE = end_of_first_HE + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence += [
                Pulse(microwave_channel, start_of_second_HE, pi_half_time),
                Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau - pi_time/2., pi_time),
                Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau + tau - pi_half_time/2., three_pi_half_time)
            ]

            #end_of_second_HE = start_of_second_HE + pi_half_time / 2. + tau + tau - pi_half_time / 2. + pi_half_time
            end_of_second_HE = start_of_second_HE + pi_half_time / 2. + tau + tau - pi_half_time / 2. + three_pi_half_time

            pulse_sequence += [
                Pulse('laser', end_of_second_HE + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', end_of_second_HE + delay_mw_readout + delay_readout, meas_time)
            ]
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

        if data['fits'] is not None and 1==0:  # I don't want to see the fit, it's usually shit
            counts = (-data['counts'][:,1] + data['counts'][:,0])/ (data['counts'][:,0] + data['counts'][:,1])
            tau = data['tau']
            fits = data['fits']

            axislist[0].plot(tau, counts, 'b')
           # axislist[0].hold(True)

            axislist[0].plot(tau, exp_offset(tau, fits[0], fits[1], fits[2]))
            axislist[0].set_title('T2 decay time (simple exponential, p = 1): {:2.1f} ns'.format(fits[1]))
        else:
            super(HahnEcho, self)._plot(axislist)
            axislist[0].set_title('Hahn Echo mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9))
            axislist[0].legend(labels=('Ref Fluorescence', 'T2 Data'), fontsize=8)

class XY8k(HahnEcho): # ER 5.25.2017
    """
    Frankie's version of XY8k; the original one by Emma seems to have a mistake right before the last pi/2 or 3pi/2 pulse.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('pi_pulse_time_i', 50.0, float, 'time duration of a pi pulse (in ns) for the mw_channel'),
            Parameter('pi_pulse_time_q', 50.0, float, 'time duration of a pi pulse (in ns)'),
            Parameter('pi_half_pulse_time_q', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('three_pi_half_pulse_time_q',  25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('k', 1, int, 'number of pi pulse blocks of 8 in the XY8-k sequence')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'minimum time between pi pulses'),
            Parameter('max_time', 10000, float, 'maximum time between pi pulses'),
            Parameter('time_step', 200, float,
                  'time step increment of time between pi pulses (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 1000, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'commander': Commander}
    from b26_toolkit.scripts.find_nv import FindNvStrobe
    _SCRIPTS = {'find_nv': FindNvStrobe}

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
        tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],
                             self.settings['tau_times']['time_step'])
        tau_list = np.ndarray.tolist(tau_list)  # 20180731 ER convert to list

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= 15]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel_i = 'microwave_+i'
        microwave_channel_q = 'microwave_+q'
        microwave_channel_q_neg = 'microwave_-q'
        pi_time_i = self.settings['mw_pulses']['pi_pulse_time_i']
        pi_time_q = self.settings['mw_pulses']['pi_pulse_time_q']
        pi_half_time_q = self.settings['mw_pulses']['pi_half_pulse_time_q']
        # three_pi_half_time_q = self.settings['mw_pulses']['three_pi_half_pulse_time_q']
        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        k = int(self.settings['mw_pulses']['k'])

        for tau in tau_list:
            current_time = 0
            pulse_sequence = []
            for part in range(2):
                pulse_sequence += [Pulse(microwave_channel_q, current_time + laser_off_time, pi_half_time_q)]
                current_time += laser_off_time + pi_half_time_q

                for k_index in range(k):
                    for i in range(8):
                        if i == 0:
                            pulse_sequence += [Pulse(microwave_channel_i, current_time + tau, pi_time_i)]
                            current_time += tau + pi_time_i
                        elif i in [1, 3, 4, 6]:
                            pulse_sequence += [Pulse(microwave_channel_q, current_time + 2 * tau, pi_time_q)]
                            current_time += 2 * tau + pi_time_q
                        elif i in [2, 5, 7]:
                            pulse_sequence += [Pulse(microwave_channel_i, current_time + 2 * tau, pi_time_i)]
                            current_time += 2 * tau + pi_time_i
                    current_time += tau

                if part == 0:
                    pulse_sequence += [Pulse(microwave_channel_q, current_time, pi_half_time_q)]
                    pulse_sequence += [Pulse('laser', current_time + pi_half_time_q + delay_mw_readout, nv_reset_time)]
                    pulse_sequence += [Pulse('apd_readout', current_time + pi_half_time_q + delay_mw_readout + delay_readout, meas_time)]
                    current_time += pi_half_time_q + delay_mw_readout + nv_reset_time

                elif part == 1:
                    pulse_sequence += [Pulse(microwave_channel_q_neg, current_time, pi_half_time_q)]
                    pulse_sequence += [Pulse('laser', current_time + pi_half_time_q + delay_mw_readout, nv_reset_time)]
                    pulse_sequence += [Pulse('apd_readout', current_time + pi_half_time_q + delay_mw_readout + delay_readout, meas_time)]
                    current_time += pi_half_time_q + delay_mw_readout + nv_reset_time

            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, meas_time
