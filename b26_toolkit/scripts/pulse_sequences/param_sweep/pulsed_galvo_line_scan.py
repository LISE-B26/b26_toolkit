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

"""
RENAME THIS TO pusled_line_scan WHEN YOU GET THE CHANCE!!
"""
import numpy as np
import scipy as sp
from timeit import timeit as timeit
from b26_toolkit.scripts.pulse_sequences.param_sweep.param_sweep_generic import ParamSweepFastGeneric
from b26_toolkit.scripts.set_laser import SetLaserSingleAxis
from b26_toolkit.instruments.piezo_controller import MDT693A
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, Pulse, Commander
from pylabcontrol.core import Parameter, Script
from b26_toolkit.plotting.plots_1d import plot_pulses, plot_1d_simple_freq
import time as t


laser_pulse_end_delay = 100  # Time of end of PB pulse to AOM minus time of end of laser pulse


class PulsedGalvoLineScan(ParamSweepFastGeneric):
    """
    Varies the galvo voltage along one axis (using SetLaserSingleAxis) and stroboscopically measures the fluorescence.
    """

    _DEFAULT_SETTINGS = [
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('axis', 'x', ['x', 'y'], 'galvo scan axis'),
        Parameter('voltage_start', 0, float, 'start voltage (V) of scan (excuse the misleading variable name)'),
        Parameter('voltage_stop', 0.035, float, 'end voltage (V) of scan (excuse the misleading variable name)'),
        Parameter('range_type', 'center_range', ['start_stop', 'center_range'],
                  'start_stop: voltage range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('voltage_points', 20, int, 'number of voltages in scan'),
        Parameter('read_out', [
            Parameter('meas_time', 500, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 500, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 80, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_readout', 260, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('spacer', 0, int, 'delay (ns) added to the beginning of the pulse sequence'),
        Parameter('settle_time', .01, float, 'time wait after switching galvo voltages (s)')
    ]

    _SCRIPTS = {'set_laser_single_axis': SetLaserSingleAxis}

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
        self.sweep_params = {'param_start': self.settings['voltage_start'],
                             'param_stop': self.settings['voltage_stop'],
                             'param_points': self.settings['voltage_points'],
                             'param_switching_time': self.settings['settle_time']}



    def _create_pulse_sequences(self, get_duration=False):

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

        tau = 200  # Placeholder, since we need to generate a tau list for the generic script to "sweep" through
        pulse_sequences = []
        tau_list = [tau]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        spacer = self.settings['spacer']

        for tau in tau_list:
            pulse_sequence = []

            pulse_sequence += [Pulse('laser', spacer + laser_off_time, nv_reset_time)]
            pulse_sequence += [Pulse('apd_readout', spacer + laser_off_time + delay_readout, meas_time)]
            pulse_sequences.append(pulse_sequence)

        if get_duration:
            duration_1 = delay_readout + meas_time
            duration_2 = nv_reset_time
            return np.min([duration_1, duration_2])
        else:
            return pulse_sequences, tau_list, meas_time

    def _configure_instruments_start_of_script(self):
        """
        Configure instruments at the beginning of a script, e.g. enable IQ modulation
        :return: None
        """
        pass

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
        self.scripts['set_laser_single_axis'].settings['axis'] = self.settings['axis']
        self.scripts['set_laser_single_axis'].settings['point'] = float(param_current)
        self.scripts['set_laser_single_axis'].run(verbose=False)

    def _configure_param_array(self):
        """
        Contruct the frequency array and store it in a variable called 'self.params'. Despite the naming, it's just a list of parameters to be swept;
        it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        :return:
        """

        if 'voltage_points' not in self.settings:
            self.params = [self.settings['voltage_start']]
        elif self.settings['range_type'] == 'start_stop':
            if self.settings['voltage_start'] > self.settings['voltage_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True
            self.params = np.linspace(self.settings['voltage_start'], self.settings['voltage_stop'], self.settings['voltage_points'])

        elif self.settings['range_type'] == 'center_range':
            self.params = np.linspace(self.settings['voltage_start'] - self.settings['voltage_stop'] / 2,
                                      self.settings['voltage_start'] + self.settings['voltage_stop'] / 2, self.settings['voltage_points'])

        if np.min(self.params) < -2 or np.max(self.params) > 2:
            self.log('start or stop voltage out of bounds')
            self._abort = True

    def _plot(self, axes_list, data=None):

        super(PulsedGalvoLineScan, self)._plot(axes_list, data)
        axes_list[0].set_title('Galvo Single-axis Scan with Pulsed Readout')
        axes_list[0].set_xlabel('Galvo voltage (V)')


class AutofocusPulsedLineScan(PulsedGalvoLineScan):
    """
    Finds focus of NV along Z-axis by monitoring fluorescence while varying focus piezo voltage (with fixed galvo X&Y).
    This is NOT based on AutofocusGeneric: it does not scan the galvo X and Y voltages at each Z value. It is therefore much faster.
    However, this might not find the focus properly if there's significant coupling between XY and the Z axis.
    """

    _DEFAULT_SETTINGS = [
        Parameter('time_per_pt', .05, float),
        Parameter('fitting', 'Gaussian', ['max', 'Gaussian'], 'fitting algorithm to find the focus'),
        Parameter('voltage_start', 50, float, 'start voltage (V) of scan (excuse the misleading variable name)'),
        Parameter('voltage_stop', 20, float, 'end voltage (V) of scan (excuse the misleading variable name)'),
        Parameter('range_type', 'center_range', ['start_stop', 'center_range'],
                  'start_stop: voltage range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('voltage_points', 15, int, 'number of voltages in scan'),
        Parameter('read_out', [
            Parameter('meas_time', 500, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 500, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 80, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_readout', 260, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('spacer', 0, int, 'delay (ns) added to the beginning of the pulse sequence'),
        Parameter('settle_time', .01, float,
                  'time wait after switching galvo voltages (s)')
    ]

    _INSTRUMENTS = {'z_piezo': MDT693A, 'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'commander': Commander}

    def __init__(self, instruments, scripts, name=None, settings=None, log_function=None, data_path=None):
        """
        Standard script initialization
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_generic import PulsedExperimentGeneric
        self._DEFAULT_SETTINGS += PulsedExperimentGeneric._DEFAULT_SETTINGS
        self._DEFAULT_SETTINGS += [Parameter('num_averages', 1000000, int, 'number of averages')]  # This is needed, but hidden in the displayed settings
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.ref_index = 0

    def _configure_instruments_start_of_script(self):
        """
        Configure instruments at the beginning of a script, e.g. enable IQ modulation
        :return: None
        """

        sequence_duration = self._create_pulse_sequences(get_duration=True) / 1e9
        self.settings['num_averages'] = int(self.settings['time_per_pt'] / sequence_duration)
        self.settings['averaging_block_size'] = self.settings['num_averages']
        print('Running %i averages per pixel' % self.settings['num_averages'])

    def _configure_instruments_end_of_script(self):
        """
        Configure instruments right before the script finishes, e.g. turn off function generators
        :return: None
        """
        focus_voltage, fit_params = self.fit_focus()
        try:
            self.log('Setting piezo voltage to %.1f V' % float(focus_voltage))
            self.instruments['z_piezo']['instance'].set_voltage(float(focus_voltage))
        except ValueError:
            self.log('Error setting piezo voltage')
            self._abort = True

    def _configure_instruments_start_of_sweep(self, param_current):
        """
        Updates z-piezo voltage
        :return: None
        """

        try:
            self.instruments['z_piezo']['instance'].set_voltage(float(param_current))
        except ValueError:
            self.log('Error setting piezo voltage')
            self._abort = True

    def _configure_param_array(self):
        """
        Contruct the frequency array and store it in a variable called 'self.params'. Despite the naming, it's just a list of parameters to be swept;
        it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        :return:
        """

        if self.settings['range_type'] == 'start_stop':
            if self.settings['voltage_start'] > self.settings['voltage_stop']:
                self.log('end voltage must be larger than start voltage when range_type is start_stop. Abort script')
                self._abort = True
            self.params = np.linspace(self.settings['voltage_start'], self.settings['voltage_stop'], self.settings['voltage_points'])

        elif self.settings['range_type'] == 'center_range':
            self.params = np.linspace(self.settings['voltage_start'] - self.settings['voltage_stop'] / 2,
                                      self.settings['voltage_start'] + self.settings['voltage_stop'] / 2, self.settings['voltage_points'])

    def fit_focus(self):
        """
        fit the data and set piezo to focus spot

        if fails return None otherwise it returns the voltage for the piezo
        :return: voltage of focus, fit parameters
        """

        self.data['fit_parameters'] = [0, 0, 0, 0]
        self.data['fit_focus'] = []
        counts = np.array(self.data['count_data'])[:, 0]

        if self.settings['fitting'] == 'Gaussian':
            noise_guess = np.min(counts)
            amplitude_guess = np.max(counts) - noise_guess
            center_guess = np.mean(self.data['params'])
            width_guess = 0.8

            p2 = [noise_guess, amplitude_guess, center_guess, width_guess]

            return_voltage = None
            try:
                p2, success = sp.optimize.curve_fit(self.gaussian, self.data['params'],
                                                    counts, p0=p2,
                                                    bounds=([0, [np.inf, np.inf, 100., 100.]]), max_nfev=2000)

                return_voltage = p2[2]
                self.log('Gaussian fit parameters: ' + str(p2))
            except(ValueError, RuntimeError):
                self.log('Could not converge to fit parameters, keeping piezo at final position')
            finally:
                # even if there is an exception we want to script to continue
                sweep_voltages = self.data['params']
                if return_voltage is not None:
                    if return_voltage > np.max(sweep_voltages):
                        return_voltage = np.max(sweep_voltages)
                        self.log('Best fit found center to be above max sweep range, setting voltage to max, {:0.3f} V'.format(return_voltage))
                    elif return_voltage < np.min(sweep_voltages):
                        return_voltage = float(np.min(sweep_voltages))
                        self.log('Best fit found center to be below min sweep range, setting voltage to min, {:0.3f} V'.format(return_voltage))

        elif self.settings['fitting'] == 'max':
            return_voltage = self.data['params'][np.argmax(counts)]
            p2 = [0, 0, 0, 0]

        self.data['fit_focus'] = [return_voltage]
        self.data['fit_parameters'] = p2
        return return_voltage, p2

    def gaussian(self, x, noise, amp, center, width):
        return noise + amp * np.exp(-1.0 * (np.square((x - center)) / (2 * (width ** 2))))

    def _plot(self, axes_list, data=None):

        super(AutofocusPulsedLineScan, self)._plot(axes_list, data)

        if data is None:
            data = self.data
        if 'fit_parameters' in data.keys() and data['fit_parameters'] is not None and not (np.array_equal(data['fit_parameters'], [0, 0, 0, 0])):
            sweep_voltages = self.data['params']
            plot_1d_simple_freq(axes_list[0], sweep_voltages, [self.gaussian(sweep_voltages, *self.data['fit_parameters'])], alpha=0.7)

        axes_list[0].set_title('Stroboscopic Autofocus')
        axes_list[0].set_xlabel('Focus piezo voltage (Vz)')


class PulsedEsrLineScan(PulsedGalvoLineScan):
    """
    Based on PulsedESRFast. Here we start the DAQ and PB pulse sequence before each freq sweep, and the read the DAQ after the freq sweep
    This reduces the overhead from the typical way of starting and reading the DAQ for each freq value.
    """

    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('axis', 'x', ['x', 'y'], 'galvo scan axis'),
        Parameter('freq', 2.87e9, float, 'start frequency of scan in Hz'),

        Parameter('voltage_start', 0, float, 'start voltage (V) of scan (excuse the misleading variable name)'),
        Parameter('voltage_stop', 0.035, float, 'end voltage (V) of scan (excuse the misleading variable name)'),
        Parameter('range_type', 'center_range', ['start_stop', 'center_range'],
                  'start_stop: voltage range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('voltage_points', 20, int, 'number of voltages in scan'),

        Parameter('read_out', [
            Parameter('meas_time', 500, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 500, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 80, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 40, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 260, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('spacer', 0, int, 'delay (ns) added to the beginning of the pulse sequence'),
        Parameter('settle_time', .01, float,
                  'time wait after switching galvo voltages (s)')
    ]
    from b26_toolkit.scripts.set_laser import SetAttoSingleAxis

    _SCRIPTS = {'set_atto_single_axis': SetAttoSingleAxis}

    def _create_pulse_sequences(self, get_duration=False):
        """
        New version where the MW pulse is after the laser pulse.
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
        spacer = self.settings['spacer']

        for tau in tau_list:
            pulse_sequence = []
            pulse_sequence += [Pulse(microwave_channel, spacer + laser_off_time, tau)]
            pulse_sequence += [Pulse('laser', spacer + laser_off_time + tau + delay_mw_readout, nv_reset_time)]
            pulse_sequence += [Pulse('apd_readout', spacer + laser_off_time + tau + delay_mw_readout + delay_readout, meas_time)]

            pulse_sequences.append(pulse_sequence)

        if get_duration:
            duration_1 = delay_readout + meas_time
            duration_2 = nv_reset_time
            return np.min([duration_1, duration_2])
        else:
            return pulse_sequences, tau_list, meas_time

    def _configure_instruments_start_of_script(self):
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['freq']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

    def _configure_instruments_start_of_sweep(self, param_current):
        """
        Configure instruments before running a sweep. For example, change MW frequency
        :return: None
        """
        self.scripts['set_atto_single_axis'].settings['axis'] = self.settings['axis']
        self.scripts['set_atto_single_axis'].settings['point'] = float(param_current)
        self.scripts['set_atto_single_axis'].run(verbose=True)

    def _configure_param_array(self):
        """
        Contruct the frequency array and store it in a variable called 'self.params'. Despite the naming, it's just a list of parameters to be swept;
        it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        :return:
        """

        if self.settings['range_type'] == 'start_stop':
            if self.settings['voltage_start'] > self.settings['voltage_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True
            self.params = np.linspace(self.settings['voltage_start'], self.settings['voltage_stop'], self.settings['voltage_points'])

        elif self.settings['range_type'] == 'center_range':
            self.params = np.linspace(self.settings['voltage_start'] - self.settings['voltage_stop'] / 2,
                                      self.settings['voltage_start'] + self.settings['voltagestop'] / 2, self.settings['voltage_points'])

        if np.min(self.params) < 0 or np.max(self.params) > 150:  # Set the limit to be 100 V. Might be able to increase this a bit in cryo
            self.log('start or stop voltage out of bounds')
            self._abort = True

    def _plot(self, axes_list, data=None):
        super(PulsedEsrLineScan, self)._plot(axes_list, data)
        axes_list[0].set_title('Attocube Line Scan with PulsedEsr')
        axes_list[0].set_xlabel('Attocube voltage (V)')
