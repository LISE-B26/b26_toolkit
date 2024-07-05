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
from b26_toolkit.scripts import FindNvPulsed, Esr
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, MicrowaveGenerator, Pulse
from pylabcontrol.core import Parameter
from b26_toolkit.data_processing.fit_functions import cose_with_decay, fit_exp_decay
from b26_toolkit.plotting.plots_1d import plot_1d_simple_timetrace, plot_pulses


class ReadoutStartTime(PulsedExperimentGeneric):  # ER 10.21.2017
    """
    This script sweeps the start time of the APD readout window, using a fixed-duration readout pulse.
    The goal is to figure out when to turn on the readout window to maximize SNR.
    To maximize sensitivity and protect against laser power drift noise, we use the 'double-init' scheme.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulse', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', '+i', ['+i', '-i', '+q', '-q'], 'Channel to use for mw pulses'),
            Parameter('pi_time', 30.0, float, 'pi time in ns')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 15, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 200, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 5, [5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000, 500000],
                      'time step increment of readout pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('readout_window', 300, int, 'length of readout window')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        # COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulse']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulse']['mw_frequency']})
        super(ReadoutStartTime, self)._function(self.data)


        counts = self.data['counts'][:, 1] # / self.data['counts'][:, 0]
        tau = self.data['tau']

        try:
            fits = fit_exp_decay(tau, counts, varibale_phase=True)
            self.data['fits'] = fits
        except:
            self.data['fits'] = None
            self.log('fit failed')

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
        pulse_sequences = []
        # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
        #                  self.settings['tau_times']['time_step'])
        # JG 16-08-25 changed (15ns min spacing is taken care of later):
        tau_list = list(range(int(self.settings['tau_times']['min_time']), int(self.settings['tau_times']['max_time']),
                         self.settings['tau_times']['time_step']))

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= 15]
        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulse']['microwave_channel']

        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        pi_time = self.settings['mw_pulse']['pi_time']
        meas_time = self.settings['read_out']['readout_window']


        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', laser_off_time, nv_reset_time),
                 Pulse('apd_readout', laser_off_time+ delay_readout + tau, meas_time),
                 ]


            pulse_sequence += [Pulse(microwave_channel, laser_off_time + nv_reset_time + laser_off_time, pi_time)]

            pulse_sequence += [
                Pulse('laser', laser_off_time + nv_reset_time + laser_off_time + pi_time + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', laser_off_time + nv_reset_time + laser_off_time + pi_time + delay_mw_readout + delay_readout + tau, meas_time)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _plot(self, axislist, data=None):
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

        if data['fits'] is not None:
            counts = data['counts'][:, 1] / data['counts'][:, 0]
            tau = data['tau']
            fits = data['fits']  # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            axislist[0].hold(True)

            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            # pi_time = 2*np.pi / fits[1] / 2
            pi_time = (np.pi - fits[2]) / fits[1]
            pi_half_time = (np.pi / 2 - fits[2]) / fits[1]
            three_pi_half_time = (3 * np.pi / 2 - fits[2]) / fits[1]
            rabi_freq = 1000 * fits[1] / (2 * np.pi)
            #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
            axislist[0].set_title('Readout pulse width counts')
        else:
            super(ReadoutStartTime, self)._plot(axislist)
            axislist[0].set_title('Readout pulse width counts')
            axislist[0].legend(labels=('Ref Fluorescence', 'Pi pulse Data'), fontsize=8)

class ReadoutStartTimeResonant(ReadoutStartTime):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulse', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', '+i', ['+i', '-i', '+q', '-q'], 'Channel to use for mw pulses'),
            Parameter('pi_time', 30.0, float, 'pi time in ns')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 15, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 200, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 5, [5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000, 500000],
                      'time step increment of readout pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('readout_window', 300, int, 'length of readout window'),
            Parameter('red_on_time', 1000, int, 'red laser on time [ns]'),
            Parameter('red_off_time', 1000, int, 'red laser off time [ns]')
        ]),
        Parameter('with_pi_seq', True, bool, 'add mw pulse sequence'),
        Parameter('num_averages', 100000, int, 'number of averages'),
    ]


    def _function(self):
        # COMMENT_ME

        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulse']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulse']['mw_frequency']})
        super(ReadoutStartTime, self)._function(self.data)

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
        pulse_sequences = []
        # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
        #                  self.settings['tau_times']['time_step'])
        # JG 16-08-25 changed (15ns min spacing is taken care of later):
        tau_list = list(range(int(self.settings['tau_times']['min_time']), int(self.settings['tau_times']['max_time']),
                         self.settings['tau_times']['time_step']))

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= 15]
        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulse']['microwave_channel']

        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        pi_time = self.settings['mw_pulse']['pi_time']
        meas_time = self.settings['read_out']['readout_window']
        red_on_time = self.settings['read_out']['red_on_time']
        red_off_time = self.settings['read_out']['red_off_time']

        for tau in tau_list:
            pulse_sequence = [
                Pulse('laser', red_off_time, nv_reset_time),
                Pulse('red_laser', red_off_time + nv_reset_time + laser_off_time, red_on_time),
                Pulse('apd_readout', red_off_time + nv_reset_time + laser_off_time + delay_readout + tau, meas_time),
            ]

            end_of_first_meas = red_off_time + nv_reset_time + laser_off_time + red_on_time + red_off_time

            if self.settings['with_pi_seq']:
                pulse_sequence += [
                    Pulse('laser', end_of_first_meas, nv_reset_time),
                    Pulse(microwave_channel, end_of_first_meas + nv_reset_time + laser_off_time, pi_time),
                    Pulse('red_laser', end_of_first_meas + nv_reset_time + laser_off_time + pi_time + delay_mw_readout, red_on_time),
                    Pulse('apd_readout', end_of_first_meas + nv_reset_time + laser_off_time + pi_time + delay_mw_readout + delay_readout + tau, meas_time)
                ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _plot(self, axislist, data=None):
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

        super(ReadoutStartTime, self)._plot(axislist)
        axislist[0].set_title('Readout pulse width counts')
        axislist[0].legend(labels=('Ref Fluorescence', 'Pi pulse Data'), fontsize=8)

class ReadoutDuration(PulsedExperimentGeneric):
    """
  This script sweeps the readout pulse duration. To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time
      """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulse', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_time', 30.0, float, 'pi time in ns')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 15, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 200, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 5, [5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000, 500000],
                      'time step increment of readout pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        # COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulse']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulse']['mw_frequency']})
        super(ReadoutDuration, self)._function(self.data)

        counts = self.data['counts'][:, 1]  # / self.data['counts'][:, 0]
        tau = self.data['tau']

        try:
            fits = fit_exp_decay(tau, counts, varibale_phase=True)
            self.data['fits'] = fits
        except:
            self.data['fits'] = None
            self.log('fit failed')

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
        pulse_sequences = []
        # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
        #                  self.settings['tau_times']['time_step'])
        # JG 16-08-25 changed (15ns min spacing is taken care of later):
        tau_list = range(int(self.settings['tau_times']['min_time']), int(self.settings['tau_times']['max_time']),
                         self.settings['tau_times']['time_step'])

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= 15]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulse']['microwave_channel']

        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        pi_time = self.settings['mw_pulse']['pi_time']
        meas_time = 100 # Placeholder

        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', laser_off_time, nv_reset_time),
                 Pulse('apd_readout', laser_off_time + delay_readout, tau),
                 ]

            pulse_sequence += [Pulse(microwave_channel, laser_off_time + nv_reset_time + laser_off_time, pi_time)]

            pulse_sequence += [
                Pulse('laser', laser_off_time + nv_reset_time + laser_off_time + pi_time + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', laser_off_time + nv_reset_time + laser_off_time + pi_time + delay_mw_readout + delay_readout,
                      tau)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _plot(self, axislist, data=None):
        """
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed
        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional) dataset to plot (dictionary that contains keys counts, tau, fits), if not provided use self.data
        """

        if data is None:
            data = self.data

        if data['fits'] is not None:
            counts = data['counts'][:, 1] / data['counts'][:, 0]
            tau = data['tau']
            fits = data['fits']  # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            axislist[0].hold(True)

            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            # pi_time = 2*np.pi / fits[1] / 2
            pi_time = (np.pi - fits[2]) / fits[1]
            pi_half_time = (np.pi / 2 - fits[2]) / fits[1]
            three_pi_half_time = (3 * np.pi / 2 - fits[2]) / fits[1]
            rabi_freq = 1000 * fits[1] / (2 * np.pi)
            #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
            axislist[0].set_title('Readout pulse width counts')
        else:
            super(ReadoutDuration, self)._plot(axislist)
            axislist[0].set_title('Readout pulse width counts')
            axislist[0].legend(labels=('Ref Fluorescence', 'Pi pulse Data'), fontsize=8)

class ChoppedInit(PulsedExperimentGeneric):
    """
  This script sweeps the number of chopped laser pulses to find the optimum for initialization.
      """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulse', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', '+i', ['+i', '-i', '+q', '-q'], 'Channel to use for mw pulses'),
            Parameter('pi_time', 30.0, float, 'pi time in ns')
        ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 80, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 80, float, 'spacing between each chopped laser pulse'),
            Parameter('min_n', 1, int, 'min num of initialization laser pulses'),
            Parameter('max_n', 20, int, 'max num of initialization laser pulses'),
            Parameter('n_step', 1, int, 'step increment of number of initialization laser pulses')
        ]),
        Parameter('compare_to_square_pulse', False, bool, 'compare to initialization with square, continuous pulses'),
        Parameter('read_out', [
            Parameter('meas_time', 300, int, 'measurement time after initialization sequence'),
            Parameter('nv_reset_time', 1750, int, 'time with CW laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        # COMMENT_ME
        self.data['fits'] = None

        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulse']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulse']['mw_frequency']})
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})
        super(ChoppedInit, self)._function(self.data)

        counts = self.data['counts'][:, 1]  # / self.data['counts'][:, 0]
        tau = self.data['tau']

        try:
            fits = fit_exp_decay(tau, counts, varibale_phase=True)
            self.data['fits'] = fits
        except:
            self.data['fits'] = None
            self.log('fit failed')

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
        pulse_sequences = []

        tau_list = range(int(self.settings['initialization_laser']['min_n']),
                         int(self.settings['initialization_laser']['max_n']),
                         int(self.settings['initialization_laser']['n_step']))
        n_list = tau_list

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulse']['microwave_channel']

        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        pi_time = self.settings['mw_pulse']['pi_time']
        meas_time = self.settings['read_out']['meas_time']

        for n in list(n_list):
            pulse_sequence = []
            current_time = laser_off_time
            for i in range(int(n)):
                pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

            pulse_sequence += [
                Pulse('laser', current_time + laser_off_time, nv_reset_time),
                Pulse('apd_readout', current_time + laser_off_time + delay_readout, meas_time)]

            pulse_sequence += [Pulse(microwave_channel, current_time + laser_off_time + nv_reset_time + laser_off_time, pi_time)]

            current_time = current_time + laser_off_time + nv_reset_time + laser_off_time + pi_time + delay_mw_readout
            for i in range(int(n)):
                pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

            pulse_sequence += [
                Pulse('laser', current_time + laser_off_time, nv_reset_time),
                Pulse('apd_readout', current_time + laser_off_time + delay_readout, meas_time)]

            if self.settings['compare_to_square_pulse']:
                current_time = current_time + laser_off_time + nv_reset_time + laser_off_time
                if n > 0:
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed*n))
                current_time = current_time + nv_reset_time_pulsed*n

                pulse_sequence += [
                    Pulse('laser', current_time + laser_off_time, nv_reset_time),
                    Pulse('apd_readout', current_time + laser_off_time + delay_readout, meas_time)]

                pulse_sequence += [Pulse(microwave_channel, current_time + laser_off_time + nv_reset_time + laser_off_time, pi_time)]

                current_time = current_time + laser_off_time + nv_reset_time + laser_off_time + pi_time + delay_mw_readout
                if n > 0:
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed*n))
                current_time = current_time + nv_reset_time_pulsed*n

                pulse_sequence += [
                    Pulse('laser', current_time + laser_off_time, nv_reset_time),
                    Pulse('apd_readout', current_time + laser_off_time + delay_readout, meas_time)]




            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _plot(self, axislist, data=None):
        """
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed
        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional) dataset to plot (dictionary that contains keys counts, tau, fits), if not provided use self.data
        """

        if data is None:
            data = self.data

        if data['fits'] is not None:
            counts = data['counts'][:, 1] / data['counts'][:, 0]
            tau = data['tau']
            fits = data['fits']  # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            axislist[0].hold(True)

            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            # pi_time = 2*np.pi / fits[1] / 2
            pi_time = (np.pi - fits[2]) / fits[1]
            pi_half_time = (np.pi / 2 - fits[2]) / fits[1]
            three_pi_half_time = (3 * np.pi / 2 - fits[2]) / fits[1]
            rabi_freq = 1000 * fits[1] / (2 * np.pi)
            #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
            axislist[0].set_title('Readout pulse width counts')
        else:
            super(ChoppedInit, self)._plot(axislist)
            axislist[0].set_title('Readout pulse width counts')
            axislist[0].legend(labels=('Ref Fluorescence', 'Pi pulse Data'), fontsize=8)

class ReadoutStartTimeWithoutMW(PulsedExperimentGeneric):
    """

    This script sweeps the start time of the APD readout window, using a fixed-duration readout pulse.
    The goal is to figure out when to turn on the readout window to maximize SNR.
    Does not include a MW pulse (so no double-init scheme)

    """
    _DEFAULT_SETTINGS = [
        Parameter('count_source_pulse_width', 10000, int, 'How long to pulse the count source (in ns)'),
        Parameter('measurement_gate_pulse_width', 15, int, 'How long to have the DAQ acquire data (in ns)'),
        Parameter('min_delay', 0, int, 'minimum delay over which to scan'),
        Parameter('max_delay', 1000, int, 'maximum delay over which to scan'),
        Parameter('delay_interval_step_size', 15, int, 'Amount delay is increased for each new run'),
        Parameter('num_averages', 1000, int, 'number of times to average for each delay'),
        Parameter('reset_time', 10000, int, 'How long to wait for laser to turn off and reach steady state'),
    ]

    _SCRIPTS = {'find_nv': FindNvPulsed, 'esr': Esr}
    def _create_pulse_sequences(self):
        '''
        Creates a pulse sequence with no pulses for a reset time, then the laser on for a count_source_pulse_width time.
        The daq measurement window is then swept across the laser pulse window to measure counts from min_delay
        (can be negative) to max_delay, which are both defined with the laser on at 0.

        Returns: pulse_sequences, num_averages, tau_list, measurement_gate_width
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            measurement_gate_width: the width (in ns) of the daq measurement


        '''
        pulse_sequences = []
        gate_delays = list(range(self.settings['min_delay'], self.settings['max_delay'], self.settings['delay_interval_step_size']))
        reset_time = self.settings['reset_time']
        for delay in gate_delays:
            pulse_sequences.append([Pulse('laser', reset_time, self.settings['count_source_pulse_width']),
                                    Pulse('apd_readout', delay + reset_time,
                                           self.settings['measurement_gate_pulse_width'])
                                    ])
        return pulse_sequences, gate_delays, self.settings['measurement_gate_pulse_width']

    def _plot(self, axes_list, data=None):
        """
        Very similar to plot for PulseBlasterBaseScript but here deals with case where first plot has counts=[], using
        PulseBlasterBaseScript plot leads to errors
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional) dataset to plot (dictionary that contains keys counts, tau), if not provided use self.data
        """
        if data is None:
            data = self.data

        counts = data['counts']
        x_data = data['tau']
        axis1 = axes_list[0]
        plot_1d_simple_timetrace(axis1, x_data, [counts])
        axis2 = axes_list[1]
        plot_pulses(axis2, self.pulse_sequences[self.sequence_index])

# if __name__ == '__main__':

    # from pylabcontrol.core.scripts import Script # import script, AS and ER 20180426
    #
    # script = {}
    # instr = {}
    # script, failed, instr = Script.load_and_append({'PulseDelays': 'PulseDelays'}, script, instr)
    #
    # print(script)
    # print(('failed', failed))
    # print(instr)

