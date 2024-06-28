from pylabcontrol.core import Script, Parameter

# import standard libraries
import numpy as np
from scipy.signal import periodogram
from b26_toolkit.instruments import RFGenerator, MokuLockInAmplifier
from b26_toolkit.tools.utils import get_param_array


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
        # Parameter('daq_type', 'Moku', ['Moku', 'PCI'], 'daq type for data acquisition'),
        # Parameter('daq_channels', [
        #     Parameter('ai_main', 'ai0', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4'], 'ai channel for main output'),
        #     Parameter('ai_aux', 'ai1', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4'], 'ai channel for main output'),
        # ])
    ]

    _INSTRUMENTS = {
        'lia': MokuLockInAmplifier,
    }

    _SCRIPTS = {}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        self._DEFAULT_SETTINGS += LockInAmpliferFreqSweep._DEFAULT_SETTINGS
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path = data_path)
        self._lia = self.instruments['lia']['instance']

    def _setup_instruments(self):
        raise NotImplementedError

    def _run_sweep(self, freq_values):
        raise NotImplementedError

    # def _setup_daq(self):
    #     self.daq = self.instruments['NI6259']['instance']
    #     self.tasks = []
    #     for chan in ['ai_main', 'ai_aux']:
    #         self.daq.settings['analog_input'][self.settings['daq_channels'][chan]]['sample_rate'] = self.settings['sample_rate']
    #         self.tasks.append(self.daq.setup_AI(self.settings['daq_channels'][chan],
    #                           self.settings['integration_time'] * self.settings['sample_rate'],
    #                           continuous=False))

    # def _read_data(self, channel='1'):
    #     if self.settings['daq_type'] == 'Moku':
    #         self._lia.start_data_streaming(duration=self.settings['integration_time'],
    #                                        rate=self.settings['sample_rate'])
    #         return self._lia.get_stream_data(self.settings['integration_time'] * self.settings['sample_rate'], channel='1')
    #     elif self.settings['daq_type'] == 'PCI':
    #         if channel == 'both':
    #             self.daq.run(self.tasks)
    #             self.data.read(self.tasks)
    #         else:
    #             self.daq.run(self.tasks[int(channel) - 1])



    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """
        # if self.settings['daq_type'] != 'Moku':
        #     self._setup_daq()

        # setup the microwave generator
        self._setup_instruments()


        # get the frequencices of the sweep

        freq_values = get_param_array(self.settings['freq_start'],
                                     self.settings['freq_stop'],
                                     self.settings['freq_points'],
                                     self.settings['range_type'])

        self.data = {'frequency': freq_values}

        # get the data for a single sweep. These are raw data.
        self._run_sweep(freq_values)

        if self.settings['turn_off_after']:
            self._lia.output_aux = 'None'
    
    def _plot(self, axes_list, data=None):
        if data is None:
            data = self.data

        if len(self.outputs) == 2 and self.outputs[1] == 'Theta':
            for i in range(data['output'].shape[0]):
                axes_list[i].set_xlabel('frequency [Hz]')
                axes_list[i].set_ylabel('signal [V]')
                axes_list[i].plot(data['frequency'], data['output'][i, :], label=self.outputs[i])
                axes_list[i].legend()
        else:
            axes_list[0].set_xlabel('frequency [Hz]')
            axes_list[0].set_ylabel('signal [V]')
            for i in range(data['output'].shape[0]):
                axes_list[0].plot(data['frequency'], data['output'][i, :], label=self.outputs[i])
            axes_list[0].legend()

    def _update_plot(self, axes_list, data=None):
        if data is None:
            data = self.data

        if len(self.outputs) == 2 and self.outputs[1] == 'Theta':
            for i in range(data['output'].shape[0]):
                axes_list[i].lines[0].set_ydata(data['output'][i, :])
                axes_list[i].relim()
                axes_list[i].autoscale_view()
        else:
            for i in range(data['output'].shape[0]):
                axes_list[0].lines[i].set_ydata(data['output'][i, :])
            axes_list[0].relim()
            axes_list[0].autoscale_view()


class LockInAmplifierFreqSweepInternal(LockInAmpliferFreqSweep):
    _DEFAULT_SETTINGS = [
        Parameter('output_options', [
            Parameter('output_1', 'R', ['X', 'Y', 'R', 'Theta'], 'output for monitor 1'),
            Parameter('output_2', 'None', ['Y', 'Theta', 'None'], 'output for monitor 2')
        ])
    ]

    def _setup_instruments(self):
        output_1 = self.settings['output_options']['output_1']
        output_2 = self.settings['output_options']['output_2']

        self.outputs = [output_1]
        if output_2 != 'None' and output_1 != output_2:
            self.outputs.append(output_2)


        self._lia.output_aux_amplitude = self.settings['voltage']
        self._lia.output_aux = 'Demod'
        self._lia.output_main = 'None'
        self._lia.set_monitor(1, 'MainOutput')

    def _run_sweep(self, freq_values):


        # initialize data arrays
        self.data['output'] = np.zeros([len(self.outputs), len(freq_values)])
        indices = list(range(len(freq_values)))

        if self.settings['randomize']:
            random.shuffle(indices)

        for freq_index in range(len(indices)):
            if self._abort:
                break

            freq = freq_values[indices[freq_index]]

            # change MW frequency
            self._lia.demod_frequency = float(freq)
            time.sleep(self.settings['switching_time'])

            for i in range(len(self.outputs)):
                self._lia.output_main = self.outputs[i]
                self._lia.start_data_streaming(duration=self.settings['integration_time'],
                                               rate=self.settings['sample_rate'])
                self.data['output'][i, indices[freq_index]] = self._lia.get_stream_data(self.settings['integration_time'] * self.settings['sample_rate'],
                                                                                        channel='1')

            self.progress = int((freq_index + 1) / len(indices) * 100.)
            self.updateProgress.emit(self.progress)

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
        self.outputs = []

    def _setup_instruments(self):
        def _vpp_to_dBm(vpp):
            return 30. * np.log10(vpp ** 2 / (8 * 50))

        self.outputs = self.settings['output_options'].split(',')

        self._lia.set_monitor(1, 'MainOutput')
        self._lia.set_monitor(2, 'AuxOutput')
        self._lia.demod = 'ExternalPLL'
        self._lia.set_pll()

        self._lia.output_main, self._lia.output_aux = self.outputs

        self._rf_gen.update({'enable_output': True,
                             'enable_modulation': False,
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
        Parameter('start_or_stop', True, bool, 'start or stop PLL'),
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

    def start_pll(self):
        def _vpp_to_dBm(vpp):
            return 30. * np.log10(vpp ** 2 / (8 * 50))
        
        self._rf_gen.update({'enable_modulation': True,
                             'enable_output': True,
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

    def stop_pll(self):
        self._lia.output_aux = 'None'

    def _function(self):
        if self.settings['start_or_stop']:
            self.start_pll()
        else:
            self.stop_pll()

class TempMeasurePll(Script):
    _DEFAULT_SETTINGS = [
        Parameter('integration_time', 0.05, float, 'The TOTAL integration time'),
        Parameter('sample_rate', 100, int, 'sample rate in [Hz]'),
    ]

    _SCRIPTS = {'pll': PhaseLockedLoop}

    def _function(self):
        self.scripts['pll'].start_pll()
        lia = self.scripts['pll'].instruments['lia']['instance']

        lia.start_data_streaming(duration=self.settings['integration_time'],
                                 rate=self.settings['sample_rate'])
        self.data = {'voltage': lia.get_stream_data(self.settings['integration_time'] * self.settings['sample_rate'],
                                                    channel='1', mean=False)}

        self.scripts['pll'].stop_pll()

    def _plot(self, axes_list, data=None):
        if data is None:
            data = self.data

        f, psd = periodogram(data['voltage'], fs=self.settings['sample_rate'])
        axes_list[0].semilogy(f, psd)
        axes_list[0].set_x_label('frequency [Hz]')
        axes_list[0].set_y_label('psd [V^2 / Hz]')

