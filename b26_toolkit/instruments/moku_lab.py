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

from pylabcontrol.core import Parameter, Instrument
from moku import instruments

class MokuLab(Instrument):
    _DEFAULT_SETTINGS = Parameter([
        Parameter('ipv6_address', '[fe80::72b3:d5ff:fe87:b4f9%11]', str, 'ip address of Moku'),
    ])

    _PROBES = {}

    def read_probes(self, key = None):
        pass

    # def __del__(self):
    #     # self.end_data_streaming()
    #     self._instrument.relinquish_ownership()

    def start_data_streaming(self, duration, mode='Normal', rate=1000):
        if self._is_connected:
            return self._instrument.start_streaming(duration=duration, mode=mode, rate=rate)
        else:
            raise ConnectionError('instrument not connected yet')

    def end_data_streaming(self):
        if self._is_connected:
            try:
                self._instrument.stop_streaming()
            except:
                self.log('no stream to stop')
        else:
            raise ConnectionError('instrument not connected yet')

    def start_data_logging(self, duration=1., mode='Normal', rate=1000):
        if self._is_connected:
            self._instrument.start_logging(duration, mode=mode, rate=rate)
        else:
            raise ConnectionError('instrument not connected yet')

    def stop_data_logging(self):
        if self._is_connected:
            self._instrument.stop_streaming()
        else:
            raise ConnectionError('instrument not connected yet')

    @property
    def stream_data(self):
        if self._is_connected:
            return self._instrument.get_stream_data()
        else:
            raise ConnectionError('instrument not connected yet')


class MokuLockInAmplifier(MokuLab):
    _LIA_INPUT_CHANNEL = 1
    _OUTPUT_MAIN_OPTIONS = ['X', 'Y', 'R', 'Theta', 'Offset', 'None']
    _OUTPUT_AUX_OPTIONS = ['Y', 'Theta', 'Demod', 'Aux', 'Offset', 'None']
    _MONITOR_CHANNELS = [1, 2]
    _MONITOR_SOURCES = ['None', 'Input1', 'Input2', 'ISignal', 'QSignal', 'MainOutput', 'AuxOutput', 'Demod']
    _DEMOD_MODES = ['Internal', 'External', 'ExternalPLL', 'None']
    _PID_CHANNELS = ['Main', 'Aux', 'None']

    _DEFAULT_SETTINGS = Parameter([
        Parameter('input', [
            Parameter('impedance', '1MOhm', ['50Ohm', '1MOhm'], 'input impedance'),
            Parameter('coupling', 'DC', ['DC', 'AC'], 'input coupling'),
            Parameter('attenuation', '0dB', ['0dB', '20dB'], 'input attenuation'),
        ]),
        Parameter('output', [
            Parameter('main', 'R', _OUTPUT_MAIN_OPTIONS, 'Source for the Main LIA output'),
            Parameter('aux', 'Demod', _OUTPUT_AUX_OPTIONS, 'Source for the Auxiliary LIA output'),
            Parameter('main_offset', 0., float, 'Main output DC offset [V]'),
            Parameter('aux_offset', 0., float, 'Aux output DC offset [V]'),
        ]),
        Parameter('monitor', [
            Parameter('source_1', 'MainOutput', _MONITOR_SOURCES, 'Monitor channel 1 source'),
            Parameter('source_2', 'AuxOutput', _MONITOR_SOURCES, 'Monitor channel 2 source'),
        ]),
        Parameter('output_aux_amplitude', 0.1, float, 'aux output amplitude'),
        Parameter('demodulation', [
            Parameter('mode', 'Internal', _DEMOD_MODES, 'demodulation mode'),
            Parameter('frequency', 100000., float, 'demodulation frequency [Hz]'),
            Parameter('phase', 0., float, 'phase of demodulation [degrees]'),
        ]),
        Parameter('low_pass_filter', [
            Parameter('corner_freq', 10., float, 'corner frequency of filter [Hz]'),
            Parameter('slope', 'Slope6dB', ['Slope6dB', 'Slope12dB', 'Slope18dB', 'Slope24dB'], 'slope per octave'),
        ]),
        Parameter('gain', [
            Parameter('main', 0., float, 'Main output gain [dB]'),
            Parameter('aux', 0., float, 'Auxiliary output gain [dB]'),
            Parameter('main_invert', False, bool, 'Invert main channel gain'),
            Parameter('aux_invert', False, bool, 'Invert auxiliary channel gain'),
            Parameter('main_gain_range', '0dB', ['0dB', '14dB'], 'Main output gain range'),
            Parameter('aux_gain_range', '0dB', ['0dB', '14dB'], 'aux output gain range'),
        ]),
        Parameter('pid', [
            Parameter('channel', 'Aux', _PID_CHANNELS, 'channel to pid on'),
            Parameter('prop_gain', 1., float, 'proprotional gain in PID, [dB]'),
            Parameter('int_crossover', 100., float, 'int crossover frequency [Hz]'),
        ])
    ])

    def __init__(self, name=None, settings=None):
        super().__init__(name, settings)
        self._instrument = instruments.LockInAmp(super()._DEFAULT_SETTINGS['ipv6_address'], force_connect=True)
        self._is_connected = True

    def update(self, settings):
        super().update(settings)

        for key, value in settings.items():
            if key == 'input':
                self.settings[key].update(value)
                self._instrument.set_frontend(self._LIA_INPUT_CHANNEL, **self.settings[key])
            elif key == 'output':
                self.settings[key].update(value)
                self._instrument.set_outputs(**self.settings[key])
            elif key == 'demodulation':
                self.settings[key].update(value)
                self._instrument.set_demodulation(**self.settings[key])
            elif key == 'low_pass_filter':
                self.settings[key].update(value)
                self._instrument.set_filter(**self.settings[key])
            elif key == 'gain':
                self.settings[key].update(value)
                self._instrument.set_gain(**self.settings[key])
            elif key == 'output_aux_amplitude':
                # frequency doesnt matter when aux output is not set to oscillator
                self._instrument.set_aux_output(frequency=self.settings['demodulation']['frequency'],
                                                amplitude=self.settings[key])
            elif key == 'monitor':
                self._instrument.set_monitor(1, self.settings['monitor']['source_1'])
                self._instrument.set_monitor(2, self.settings['monitor']['source_2'])
            elif key == 'pid':
                pid_settings = self.settings[key]
                self._instrument.use_pid(pid_settings['channel'])
                self._instrument.set_by_frequency(**pid_settings)

    @property
    def demod_frequency(self):
        return self._instrument.get_demodulation()['frequency']

    @demod_frequency.setter
    def demod_frequency(self, frequency):
        self.update({'demodulation': {'frequency': frequency}})

    @property
    def demod(self):
        return self._instrument.get_demodulation()['mode']

    @demod.setter
    def demod(self, mode):
        if mode not in self._DEMOD_MODES:
            raise KeyError('{} is not an allowed. Please choose one of {}'.format(mode, self._DEMOD_MODES))
        self.update({'demodulation': {'mode': mode}})

    @property
    def output_main(self):
        return self._instrument.get_outputs()['main']

    @output_main.setter
    def output_main(self, om):
        if om not in self._OUTPUT_MAIN_OPTIONS:
            raise KeyError('{} is not an allowed. Please choose one of {}'.format(om, self._OUTPUT_MAIN_OPTIONS))
        self.update({'output': {'main': om}})

    @property
    def output_aux(self):
        return self._instrument.get_outputs()['aux']

    @output_aux.setter
    def output_aux(self, om):
        if om not in self._OUTPUT_AUX_OPTIONS:
            raise KeyError('{} is not an allowed. Please choose one of {}'.format(om, self._OUTPUT_AUX_OPTIONS))
        self.update({'output': {'aux': om}})

    @property
    def output_aux_amplitude(self):
        return self._instrument.get_aux_output()['amplitude']

    @output_aux_amplitude.setter
    def output_aux_amplitude(self, oaa):
        self.update({'output_aux_amplitude': oaa})

    def set_monitor(self, channel, source):
        if channel in self._MONITOR_CHANNELS and source in self._MONITOR_SOURCES:
            self._instrument.set_monitor(channel, source)
        else:
            raise KeyError('not allowed channel or source')

    def set_pll(self):
        self._instrument.set_pll()
    @property
    def pid(self):
        return self._instrument.get_pid()
    
    @pid.setter
    def pid(self, pid_channel):
        self.update({'pid': {'channel': pid_channel}})

    @property
    def phase_shift(self):
        return self._instrument.get_demodulation()['phase']
    
    @phase_shift.setter
    def phase_shift(self, ps):
        self._instrument.set_demodulation(self.settings['demodulation']['mode'], phase=ps)
