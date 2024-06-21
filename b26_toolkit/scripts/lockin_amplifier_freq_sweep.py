from pylabcontrol.core import Script, Parameter

# import standard libraries
import numpy as np
from b26_toolkit.instruments import RFGenerator, NI6259, NI9402, MokuLockInAmplifier
from b26_toolkit.plotting.plots_1d import plot_esr


from b26_toolkit.data_processing.esr_signal_processing import fit_esr

import time
import random


class LockInAmpliferFreqSweep(Script):

    _DEFAULT_SETTINGS = [
        Parameter('voltage', 0.1, float, 'voltage [Vpp]'),
        Parameter('freq_start', 10000., float, 'start frequency of scan'),
        Parameter('freq_stop', 20000., float, 'end frequency of scan'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'], 'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan'),
        Parameter('integration_time', 0.05, float, 'The TOTAL integration time'),
        Parameter('sample_rate', 100, int, 'sample rate in [Hz]'),
        Parameter('switching_time', .01, float, 'time wait after switching center frequencies on generator (s)'),
        Parameter('turn_off_after', False, bool, 'if true MW output is turned off after the measurement'),
        Parameter('randomize', True, bool, 'check to randomize esr frequencies'),
    ]

    _INSTRUMENTS = {
        'lia': MokuLockInAmplifier
    }

    _SCRIPTS = {}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        self._DEFAULT_SETTINGS += LockInAmpliferFreqSweep._DEFAULT_SETTINGS
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path = data_path)
        self._lia = self.instruments['lia']['instance']

    def _setup_instruments(self):
        raise NotImplementedError

    def get_freq_array(self):
        '''
        Construct a list of values through which we will sweep a parameter.
        Function is called get_freq_array, but the array does not have to be frequency. Can be e.g. voltage values for a piezo

        Returns:
            freq_values: array of the parameter values to be tested
            freq_range: the range of the parameter values to be tested, i.e. maximum frequency - minumum frequency
        '''

        if self.settings['range_type'] == 'start_stop':
            if self.settings['freq_start'] > self.settings['freq_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True

            if self.settings['freq_start'] < 0.001 or self.settings['freq_stop'] > 200e6: # freq range of the SRS
                self.log('start or stop frequency out of bounds')
                self._abort = True

            freq_values = np.linspace(self.settings['freq_start'], self.settings['freq_stop'], self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

        elif self.settings['range_type'] == 'center_range':
            if self.settings['freq_start'] < 2 * self.settings['freq_stop']:
                self.log('end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True
            freq_values = np.linspace(self.settings['freq_start']-self.settings['freq_stop']/2,
                                      self.settings['freq_start']+self.settings['freq_stop']/2, self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

            if self.settings['freq_stop'] > 1e9:
                self.log('freq_stop (range) is quite large --- did you mean to set \'range_type\' to \'start_stop\'? ')
        else:
            self.log('unknown range parameter. Abort script')
            self._abort = True

        return freq_values, freq_range

    def run_sweep(self, freq_values):
        raise NotImplementedError

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        # setup the microwave generator
        self._setup_instruments()


        # get the frequencices of the sweep
        freq_values, freq_range = self.get_freq_array()
        self.data = {'frequency': freq_values}

        # get the data for a single sweep. These are raw data.
        self.run_sweep(freq_values)

        if self.settings['turn_off_after']:
            self._lia.output_aux = None


class LockInAmplifierFreqSweepInternal(LockInAmpliferFreqSweep):
    _DEFAULT_SETTINGS = [
        Parameter('output_options', [
            Parameter('output_1', 'R', ['X', 'Y', 'R', 'Theta'], 'output for monitor 1'),
            Parameter('output_2', None, ['Y', 'Theta', None], 'output for monitor 2')
        ])
    ]

    def _setup_instruments(self):
        self._lia.output_aux_amplitude = self.settings['voltage']
        self._lia.output_aux = 'Demod'
        self._lia.output_main = None
        self._lia.set_monitor(1, 'MainOutput')

    def run_sweep(self, freq_values):
        output_1 = self.settings['output_options']['output_1']
        output_2 = self.settings['output_options']['output_2']

        outputs = [output_1]
        if output_2 is not None and output_1 != output_2:
            outputs.append(output_2)

        # initialize data arrays
        self.data['output'] = np.zeros([len(outputs), len(freq_values)])
        indices = list(range(len(freq_values)))

        if self.settings['randomize']:
            random.shuffle(indices)

        for freq_index in range(len(indices)):
            if self._abort:
                break

            freq = freq_values[indices[freq_index]]

            # change MW frequency
            self._lia.demod_frequency = freq
            time.sleep(self.settings['switching_time'])

            for i in range(len(outputs)):
                self._lia.output_main = outputs[i]
                self._lia.start_data_streaming(duration=self.settings['integration_time'],
                                               sample_rate=self.settings['sample_rate'])
                time.sleep(self.settings['integration_time'])
                self.data['output'][i, indices[freq_index]] = np.mean(self._lia.stream_data['ch1'])

            self.updateProgress((freq_index + 1) / len(indices) * 100.)

class LockInAmplifierFreqSweepExternal(LockInAmpliferFreqSweep):
    _DEFAULT_SETTINGS = [
        Parameter('output_options', 'R,Theta', ['R,Theta', 'X,Y'], 'output options for monitors')
    ]

    _INSTRUMENTS = {
        'lia': MokuLockInAmplifier,
        'rf_gen': RFGenerator
    }

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path = None):
        super().__init__(instruments, scripts=scripts, name=name, settings=settings, log_function=log_function, data_path=data_path)
        # defines which daqs contain the input and output based on user selection of daq interface
        self._rf_gen = self.instruments['rf_gen']['instance']

    def _setup_instruments(self):
        def _vpp_to_dBm(vpp):
            return 30. * np.log10(vpp ** 2 / (8 * 50))

        self._lia.set_monitor(1, 'MainOutput')
        self._lia.set_monitor(2, 'AuxOutput')
        self._lia.demod = 'ExternalPLL'
        self._lia.set_pll()
        

        self._lia.output_main, self._lia.output_aux = self.settings['output_options'].split(',')

        self._rf_gen.update({'enable_modulation': False,
                             'amplitude_rf': _vpp_to_dBm(self.settings['voltage'])})

    def run_sweep(self, freq_values):

        # initialize data arrays
        self.data['output'] = np.zeros([2, len(freq_values)])
        indices = list(range(len(freq_values)))

        if self.settings['randomize']:
            random.shuffle(indices)

        for freq_index in range(len(indices)):
            if self._abort:
                break

            freq = freq_values[indices[freq_index]]

            # change MW frequency
            self._rf_gen.update({'frequency': float(freq)})
            time.sleep(self.settings['switching_time'])

            self._lia.start_data_streaming(duration=self.settings['integration_time'],
                                           sample_rate=self.settings['sample_rate'])
            time.sleep(self.settings['integration_time'])
            self.data['output'][0, indices[freq_index]] = np.mean(self._lia.stream_data['ch1'])
            self.data['output'][1, indices[freq_index]] = np.mean(self._lia.stream_data['ch2'])

            self.updateProgress((freq_index + 1) / len(indices) * 100.)

class PhaseLockedLoop(Script):
    _DEFAULT_SETTINGS = [
        Parameter('voltage', 0.1, float, 'voltage [Vpp]'),
        Parameter('starting_freq', 10000., float, 'starting frequency of PLL [Hz]'),
        Parameter('phase', 0., float, 'set phase to lock to')
    ]

    _INSTRUMENTS = {
        'lia': MokuLockInAmplifier,
        'rf_gen': RFGenerator
    }

    _SCRIPTS = {}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function,
                        data_path=data_path)
        self._lia = self.instruments['lia']['instance']
        self._rf_gen = self.instruments['rf_gen']['instance']

    def _setup_instruments(self):
        def _vpp_to_dBm(vpp):
            return 30. * np.log10(vpp ** 2 / (8 * 50))
        
        self._rf_gen.update({'enable_modulation': True,
                             'modulation_type': 'FM',
                             'modulation_function': 'External',
                             'frequency': self.settings['starting_freq'],
                             'ampltitude_rf': _vpp_to_dBm(self.settings['voltage'])})
        
        self._lia.set_monitor(1, 'MainOutput')
        self._lia.set_monitor(2, 'AuxOutput')
        self._lia.output_main = 'R'
        self._lia.output_aux = 'Theta'
        self._lia.demod = 'ExternalPLL'
        self._lia.phase_shift = self.settings['phase']
        self._lia.set_pll()
        self._lia.pid = 'Aux'

    def _function(self):
        self._setup_instruments()
        
        
