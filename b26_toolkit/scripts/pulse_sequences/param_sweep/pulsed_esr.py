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
from timeit import timeit as timeit
from b26_toolkit.scripts.pulse_sequences.param_sweep.param_sweep_generic import ParamSweepGeneric, ParamSweepFastGeneric
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, Pulse, MicrowaveGenerator, Commander
from b26_toolkit.plotting.plots_1d import  plot_pulses, update_pulse_plot, update_1d_simple, plot_1d_simple_freq
from pylabcontrol.core import Parameter, Script
from b26_toolkit.data_processing.fit_functions import fit_pulsed_odmr, pulsed_odmr_double, pulsed_odmr_single
import time as t


laser_pulse_end_delay = 100  # Time of end of PB pulse to AOM minus time of end of laser pulse


class PulsedEsr(ParamSweepGeneric):
    """
    Pulsed version of ESR. This script applies a microwave pulse at fixed power and duration while varying its frequency.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.87e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 1e8, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 50, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 340, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 800, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 500, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 40, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 260, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
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

    def define_sweep_parameters(self):
        """
        Define name of sweep parameters. Since this script is generic, it is coded with variables like 'param_start'.
        This function redirects the code to look for the corresponding settings in _DEFAULT_SETTINGS
        :return:
        """
        self.sweep_params = {'param_start': self.settings['freq_start'],
                             'param_stop': self.settings['freq_stop'],
                             'param_points': self.settings['freq_points'],
                             'param_switching_time': self.settings['mw_generator_switching_time']}

    def _configure_instruments_start_of_script(self):
        """
        Configure instruments at the beginning of a script, e.g. enable IQ modulation
        :return: None
        """
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_power']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

    def _configure_instruments_end_of_script(self):
        """
        Configure instruments right before the script finishes, e.g. turn off function generators
        :return: None
        """
        pass

    def _configure_instruments_start_of_sweep(self, param_current):
        """
        Configure instruments before running a sweep. For example, change MW frequency
        :return: None
        """
        self.instruments['mw_gen']['instance'].update({'frequency': float(param_current)})

    def _configure_param_array(self):
        """
        Contruct the frequency array and store it in a variable called 'self.params'. Despite the naming, it's just a list of parameters to be swept;
        it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        :return:
        """

        if 'freq_points' not in self.settings:
            self.params = [self.settings['freq_start']]
        elif self.settings['range_type'] == 'start_stop':
            if self.settings['freq_start'] > self.settings['freq_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True

            if self.settings['freq_start'] < 0 or self.settings['freq_stop'] > 4.05E9:
                self.log('start or stop frequency out of bounds')
                self._abort = True

            self.params = np.linspace(self.settings['freq_start'], self.settings['freq_stop'], self.settings['freq_points'])

        elif self.settings['range_type'] == 'center_range':
            if self.settings['freq_start'] < 2 * self.settings['freq_stop']:
                self.log('end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True
            self.params = np.linspace(self.settings['freq_start'] - self.settings['freq_stop'] / 2,
                                      self.settings['freq_start'] + self.settings['freq_stop'] / 2, self.settings['freq_points'])

    @staticmethod
    def dbm_to_vpp(self, dbm):
        dbm = float(dbm)
        return np.sqrt(10 ** (dbm / 10) / 1000 * 50) * 2 * np.sqrt(2)

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

        tau = self.settings['tau_mw']
        pulse_sequences = []
        tau_list = [tau]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['microwave_channel']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', laser_off_time + tau + 2 * 40, nv_reset_time),
                 Pulse('apd_readout',
                       laser_off_time + tau + 2 * 40 + delay_readout, meas_time)]
               #  Pulse('apd_readout', laser_off_time + tau + 2 * 40 + delay_readout, meas_time),
               #  ]
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence += [
                    Pulse(microwave_channel, laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time, tau)]

            pulse_sequence += [
                Pulse('laser',
                      laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2 * 40 + delay_mw_readout,
                      nv_reset_time),
                Pulse('apd_readout',
                      laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2 * 40 + delay_mw_readout + delay_readout,
                      meas_time)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, self.settings['read_out']['meas_time']


class PulsedEsrChained(PulsedEsr):
    """
    Faster version of PulsedESR. PulsedESR is the "proper" sequence for it has separate readout windows for reference and signal fluorescence, but the DAQ read
    speed limits how quickly one can repeat one single pulse sequence, e.g.  despite only requiring ~600 ns to reinitialize, PulsedESR requires 1.5 us initialization
    time to artificially slow down the pulse sequence repetition to avoid a DAQ crash.

    This fast version chains together multiple PulsedESR sequences while only using one long readout window. There is loss of contrast from the readout window
    being on during initialization and laser off times, but the effective pulse sequence repetition rate is now much higher, leading to an overall reduction in
    averaging time needed.
    """

    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('repetitions', 4, int, 'number of repetitions of Pulsed ESR sequence consisting of MW pi-pulse and reinitialization'),
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
    ]

    def _create_pulse_sequences(self):
        """
        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement
        """

        tau = self.settings['tau_mw']
        pulse_sequences = []
        tau_list = [tau]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = []
            #pulse_sequence += [Pulse('laser', 0, 500)]
            current_time = 0
            for part in range(int(self.settings['repetitions'])):  # Repeat pulse sequence 4 times because DAQ crashes if we run really short pulse sequences (~ several us long)
                # if tau is 0 there is actually no mw pulse
                if tau > 0:
                    pulse_sequence += [
                        Pulse(microwave_channel, current_time + laser_off_time, tau)]

                pulse_sequence += [
                    Pulse('laser',
                          current_time + laser_off_time + tau + delay_mw_readout,
                          nv_reset_time)
                ]
                current_time += laser_off_time + tau + delay_mw_readout + nv_reset_time

            meas_time_long = current_time - (laser_off_time + tau + delay_mw_readout + delay_readout) - nv_reset_time + delay_readout + meas_time
            pulse_sequence += [Pulse('apd_readout',
                                     laser_off_time + tau + delay_mw_readout + delay_readout, meas_time_long)]

            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, meas_time_long


class PulsedEsrFast(ParamSweepFastGeneric):
    """
    Based on PulsedESRFast. Here we start the DAQ and PB pulse sequence before each freq sweep, and the read the DAQ after the freq sweep
    This reduces the overhead from the typical way of starting and reading the DAQ for each freq value.
    """

    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.87e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 1e8, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'center_range', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 40, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 340, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 800, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 80, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 40, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 260, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('spacer', 0, int, 'delay (ns) added to the beginning of the pulse sequence'),
        Parameter('pulse_order', 'laser_before_MW', ['laser_before_MW', 'MW_before_laser'],
                  'whether spacer -> laser -> MW, or spacer -> MW -> laser; this is only relevant for non-zero spacer duration, because this'
                  'then determines whether the spacer is between the MW and laser pulses (in order), or between laser and MW pulses'),
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
    ]

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

    def define_sweep_parameters(self):
        """
        Define name of sweep parameters. Since this script is generic, it is coded with variables like 'param_start'.
        This function redirects the code to look for the corresponding settings in _DEFAULT_SETTINGS
        :return:
        """
        self.sweep_params = {'param_start': self.settings['freq_start'],
                             'param_stop': self.settings['freq_stop'],
                             'param_points': self.settings['freq_points'],
                             'param_switching_time': self.settings['mw_generator_switching_time']}

    def _configure_instruments_start_of_script(self):
        """
        Configure instruments at the beginning of a script, e.g. enable IQ modulation
        :return: None
        """
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_power']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

    def _configure_instruments_end_of_script(self):
        """
        Configure instruments right before the script finishes, e.g. turn off function generators
        :return: None
        """
        pass

    def _configure_instruments_start_of_sweep(self, param_current):
        """
        Configure instruments before running a sweep. For example, change MW frequency
        :return: None
        """
        self.instruments['mw_gen']['instance'].update({'frequency': float(param_current)})

    def _configure_param_array(self):
        """
        Contruct the frequency array and store it in a variable called 'self.params'. Despite the naming, it's just a list of parameters to be swept;
        it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        :return:
        """

        if 'freq_points' not in self.settings:
            self.params = [self.settings['freq_start']]
        elif self.settings['range_type'] == 'start_stop':
            if self.settings['freq_start'] > self.settings['freq_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True

            if self.settings['freq_start'] < 0 or self.settings['freq_stop'] > 4.05E9:
                self.log('start or stop frequency out of bounds')
                self._abort = True

            self.params = np.linspace(self.settings['freq_start'], self.settings['freq_stop'], self.settings['freq_points'])

        elif self.settings['range_type'] == 'center_range':
            if self.settings['freq_start'] < 2 * self.settings['freq_stop']:
                self.log('end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True
            self.params = np.linspace(self.settings['freq_start'] - self.settings['freq_stop'] / 2,
                                      self.settings['freq_start'] + self.settings['freq_stop'] / 2, self.settings['freq_points'])

    def _create_pulse_sequences(self):

        '''
        New version where the MW pulse is after the laser pulse.

        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''

        tau = self.settings['tau_mw']
        pulse_sequences = []
        tau_list = [tau]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        spacer = self.settings['spacer']

        for tau in tau_list:
            pulse_sequence = []
            if self.settings['pulse_order'] == 'laser_before_MW':

                pulse_sequence += [Pulse('laser',
                                         delay_mw_readout + spacer, nv_reset_time)]
                pulse_sequence += [Pulse('apd_readout',
                                         delay_mw_readout + spacer + delay_readout, meas_time)]
                pulse_sequence += [Pulse(microwave_channel,
                                         delay_mw_readout + spacer + nv_reset_time + laser_off_time, tau)]
            elif self.settings['pulse_order'] == 'MW_before_laser':
                pulse_sequence += [Pulse(microwave_channel, spacer + laser_off_time, tau)]
                pulse_sequence += [Pulse('laser', spacer + laser_off_time + tau + delay_mw_readout, nv_reset_time)]

                pulse_sequence += [Pulse('apd_readout',
                                         spacer + laser_off_time + tau + delay_mw_readout + delay_readout, meas_time)]

            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, meas_time


class PulsedEsrUpperLower(PulsedEsrFast):
    """
    Based on PulsedESRFast. Here we start the DAQ and PB pulse sequence before each freq sweep, and the read the DAQ after the freq sweep
    This reduces the overhead from the typical way of starting and reading the DAQ for each freq value.

    This script allows scanning over two frequency ranges, e.g. upper and lower ESR frequencies for +-1 states.
    """

    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start_lower', 2.87e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop_lower', 1e8, float, 'end frequency of scan in Hz'),
        Parameter('freq_start_upper', 2.87e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop_upper', 1e8, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'center_range', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 40, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 34, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 800, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 80, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 40, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 260, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('spacer', 0, int, 'delay (ns) added to the beginning of the pulse sequence'),
        Parameter('pulse_order', 'MW_before_laser', ['laser_before_MW', 'MW_before_laser'],
                  'whether spacer -> laser -> MW, or spacer -> MW -> laser; this is only relevant for non-zero spacer duration, because this'
                  'then determines whether the spacer is between the MW and laser pulses (in order), or between laser and MW pulses'),
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
    ]

    def define_sweep_parameters(self):
        """
        Define name of sweep parameters. Since this script is generic, it is coded with variables like 'param_start'.
        This function redirects the code to look for the corresponding settings in _DEFAULT_SETTINGS
        :return:
        """
        self.sweep_params = {'param_start': self.settings['freq_start_lower'],  # Placeholder, _configure_frequency_array will override this
                             'param_stop': self.settings['freq_stop_upper'],  # Placeholder, _configure_frequency_array will override this
                             'param_points': self.settings['freq_points']*2,
                             'param_switching_time': self.settings['mw_generator_switching_time']}

    def _configure_param_array(self):
        """
        Contruct the frequency array and store it in a variable called 'self.params'. Despite the naming, it's just a list of parameters to be swept;
        it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        :return:
        """

        freq_start_lower = self.settings['freq_start_lower']
        freq_stop_lower = self.settings['freq_stop_lower']
        freq_start_upper = self.settings['freq_start_upper']
        freq_stop_upper = self.settings['freq_stop_upper']
        freq_points = self.settings['freq_points']

        if self.settings['range_type'] == 'start_stop':
            if freq_start_lower > freq_stop_lower or freq_start_upper > freq_stop_upper:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True

            if freq_start_lower < 0 or freq_stop_lower > 4.05E9 or freq_start_upper < 0 or freq_stop_upper > 4.05E9:
                self.log('start or stop frequency out of bounds')
                self._abort = True

            self.mw_frequencies_upper = np.linspace(freq_start_upper, freq_stop_upper, freq_points)
            self.mw_frequencies_lower = np.linspace(freq_start_lower, freq_stop_lower, freq_points)

        elif self.settings['range_type'] == 'center_range':
            if freq_start_lower < 2 * freq_stop_lower or freq_start_upper < 2 * freq_stop_upper:
                self.log('end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True
            self.mw_frequencies_lower = np.linspace(freq_start_lower - freq_stop_lower / 2,
                                                    freq_start_lower + freq_stop_lower / 2, freq_points)
            self.mw_frequencies_upper = np.linspace(freq_start_upper - freq_stop_upper / 2,
                                                    freq_start_upper + freq_stop_upper / 2, freq_points)
        if np.max(self.mw_frequencies_lower) > np.min(self.mw_frequencies_upper):
            self.log('WARNING: Lower frequency scan overlaps with higher frequency scan')

        self.params = np.concatenate((self.mw_frequencies_lower, self.mw_frequencies_upper))
        self.settings['freq_points'] *= 2

    def get_axes_layout(self, figure_list):
        """
        returns the axes objects the script needs to plot its data
        the default creates a single axes object on each figure
        This can/should be overwritten in a child script if more axes objects are needed
        Args:
            figure_list: a list of figure objects
        Returns:
            axes_list: a list of axes objects

        """

        axes_list = []
        figure_list[0].clf()
        figure_list[1].clf()
        axes_list.append(figure_list[0].add_subplot(121))
        axes_list.append(figure_list[1].add_subplot(111))
        axes_list.append(figure_list[0].add_subplot(122))

        return axes_list


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

        freq_points_subplot = int(len(self.params) / 2)
        fig = axes_list[0].get_figure()

        if data is None:
            data = self.data

        if 'count_data' in data.keys():
            plot_1d_simple_freq(axes_list[0], self.mw_frequencies_lower, [data['count_data'][:freq_points_subplot]])
            plot_1d_simple_freq(axes_list[2], self.mw_frequencies_upper, [data['count_data'][freq_points_subplot:]])
            plot_pulses(axes_list[1], self.pulse_sequences[self.sequence_index])
        if 'fits' in data.keys() and data['fits'] is not None:
            if len(data['fits']) == 6:
                # rabi_freq, baseline, contrast_1, contrast_2, freq_1, freq_2
                title = 'Pulsed ESR f1 = {:0.5e} Hz, f2 = {:0.5e} Hz, wo = {:0.2e} Hz,\n contrast1 = {:.1%}, contrast2 = {:.1%}' \
                    .format(data['fits'][4], data['fits'][5], data['fits'][0],
                            data['fits'][2], data['fits'][3])


            elif len(data['fits']) == 4:
                # rabi_freq, baseline, contrast_1, freq_1
                title = 'Pulsed ESR f1 = {:0.5e} Hz, wo = {:0.2e} Hz,\n contrast1 = {:.1%}' \
                    .format(data['fits'][3], data['fits'][0],
                            data['fits'][2])

            fig.suptitle(title, bbox=dict(facecolor='white', alpha=0.7))
            freq_fine_lower = np.linspace(np.min(self.mw_frequencies_lower), np.max(self.mw_frequencies_lower), len(self.mw_frequencies_lower) * 8)
            freq_fine_upper = np.linspace(np.min(self.mw_frequencies_upper), np.max(self.mw_frequencies_upper), len(self.mw_frequencies_upper) * 8)

        if 'fits' in data.keys() and data['fits'] is not None:
            if len(data['fits']) == 6:
                #plot_1d_simple_freq(axes_list[0], x_data_fine, [pulsed_odmr_double(x_data_fine, *data['fits'])], alpha=0.7, title=title)
                plot_1d_simple_freq(axes_list[0], freq_fine_lower, [pulsed_odmr_double(freq_fine_lower, *data['fits'])], alpha=0.7)
                plot_1d_simple_freq(axes_list[2], freq_fine_upper, [pulsed_odmr_double(freq_fine_upper, *data['fits'])], alpha=0.7)
            elif len(data['fits']) == 4:
                #plot_1d_simple_freq(axes_list[0], x_data_fine, [pulsed_odmr_single(x_data_fine, *data['fits'])], alpha=0.7, title=title)
                plot_1d_simple_freq(axes_list[0], freq_fine_lower, [pulsed_odmr_single(freq_fine_lower, *data['fits'])], alpha=0.7)
                plot_1d_simple_freq(axes_list[2], freq_fine_upper, [pulsed_odmr_single(freq_fine_upper, *data['fits'])], alpha=0.7)

    def _update_plot_garbage(self, axes_list):
        """
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)
        """

        data = self.data
        counts = data['count_data']
        freq_points_subplot = int(len(self.params) / 2)
        freq_fine_lower = np.linspace(np.min(self.mw_frequencies_lower), np.max(self.mw_frequencies_lower), len(self.mw_frequencies_lower) * 8)
        freq_fine_upper = np.linspace(np.min(self.mw_frequencies_upper), np.max(self.mw_frequencies_upper), len(self.mw_frequencies_upper) * 8)

        # If fit is found and fit has not been plotted, plot both data and fit
        fit_in_plot = len(axes_list[0].lines) == (len(np.transpose(counts)) + 1)
        #update_1d_simple(axes_list[0], x_data, counts, fit_in_plot=fit_in_plot)

        print(np.shape(data['count_data']))
        print(freq_points_subplot)
        print(axes_list)
        print()
        update_1d_simple(axes_list[0], self.mw_frequencies_lower, data['count_data'][:freq_points_subplot], fit_in_plot=fit_in_plot)
        update_1d_simple(axes_list[2], self.mw_frequencies_upper, data['count_data'][freq_points_subplot:], fit_in_plot=fit_in_plot)

        if 'fits' in data.keys() and data['fits'] is not None:
            if len(data['fits']) == 6:
                fit_fn = pulsed_odmr_double
                #rabi_freq, baseline, contrast_1, contrast_2, freq_1, freq_2
                title = 'Pulsed ESR f1 = {:0.6e} Hz, f2 = {:0.6e} Hz, wo = {:0.2e},\n contrast1 = {:.1%}, contrast2 = {:.1%}' \
                    .format(data['fits'][4], data['fits'][5], data['fits'][0],
                            data['fits'][2], data['fits'][3])
            elif len(data['fits']) == 4:
                fit_fn = pulsed_odmr_single
                # rabi_freq, baseline, contrast_1, freq_1
                title = 'Pulsed ESR f1 = {:0.6e} Hz, wo = {:0.2e},\n contrast1 = {:.1%}' \
                    .format(data['fits'][3], data['fits'][0],
                            data['fits'][2])

            axes_list[0].set_title(title)
            if fit_in_plot:
                for index, counts in enumerate([fit_fn(freq_fine_lower, *data['fits'])]):
                    axes_list[0].lines[-1 - index].set_ydata(data['count_data'][:freq_points_subplot])
                for index, counts in enumerate([fit_fn(freq_fine_upper, *data['fits'])]):
                    axes_list[2].lines[-1 - index].set_ydata(data['count_data'][freq_points_subplot:])

            else:
                plot_1d_simple_freq(axes_list[0], freq_fine_lower, [fit_fn(freq_fine_lower, *data['fits'])], alpha=0.7)
                plot_1d_simple_freq(axes_list[0], freq_fine_upper, [fit_fn(freq_fine_upper, *data['fits'])], alpha=0.7)

        update_pulse_plot(axes_list[1], self.pulse_sequences[self.sequence_index])

    def _update_plot(self, axes_list):
        """
        I cannot get update plot to work properly with multiple axes, so I just force a replot of everything. Should be OK because the update
        freq is quite low (at least every few s), and the number of data points is usually at most several hundred.
        """
        for axis in axes_list:
            axis.clear()
        self._plot(axes_list)
