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

class Moku(Instrument):
    _DEFAULT_SETTINGS = Parameter([
        Parameter('ip_address', 'Moku', str, 'ip address of Moku'),
    ])

    def __del__(self):
        self.end_data_acquisiton()
        self._instrument.relinquish_ownership()

    def start_data_streaming(self, duration, mode='Normal', rate=1000):
        if self._is_connected:
            return self._instrument.start_streaming(duration=duration, mode=mode, rate=rate)
        else:
            raise ConnectionError('instrument not connected yet')

    def end_data_streaming(self):
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

class MokuLockInAmplifier(Moku):
    _LIA_INPUT_CHANNEL = 1
    OUTPUT_MAIN_OPTIONS = ['X', 'Y', 'R', 'Theta', 'Offset', None]
    OUTPUT_AUX_OPTIONS = ['Y', 'Theta', 'Demod', 'Aux', 'Offset', None]

    _DEFAULT_SETTINGS = Parameter([
        Parameter('input', [
            Parameter('impedance', '1MOhm', ['50Ohm', '1MOhm'], 'input impedance'),
            Parameter('coupling', 'DC', ['DC', 'AC'], 'input coupling'),
            Parameter('attenuation', '0dB', ['0dB', '20dB'], 'input attenuation'),
            Parameter('strict', True, bool, 'Disable all implicit conversions and coercions')
        ]),
        Parameter('output', [
            Parameter('main', 'R', OUTPUT_MAIN_OPTIONS, 'Source for the Main LIA output'),
            Parameter('aux', 'Demod', OUTPUT_AUX_OPTIONS, 'Source for the Auxiliary LIA output'),
            Parameter('main_offset', 0., float, 'Main output DC offset [V]'),
            Parameter('aux_offset', 0., float, 'Aux output DC offset [V]'),
            Parameter('strict', True, bool, 'Disable all implicit conversions and coercions')
        ]),
        Parameter('output_aux_amplitude', 0.1, float, 'aux output amplitude'),
        Parameter('demodulation', [
            Parameter('mode', 'Internal', ['Internal', 'External', 'ExternalPLL', None]),
            Parameter('frequency', 100000., float, 'demodulation frequency [Hz]'),
            Parameter('phase', 0., float, 'phase of demodulation [degrees]'),
            Parameter('strict', True, bool, 'Disable all implicit conversions and coercions')
        ]),
        Parameter('low_pass_filter', [
            Parameter('corner_freq', 10., float, 'corner frequency of filter [Hz]'),
            Parameter('slope', 'Slope6dB', ['Slope6dB', 'Slope12dB', 'Slope18dB', 'Slope24dB'], 'slope per octave'),
            Parameter('strict', True, bool, 'Disable all implicit conversions and coercions')
        ]),
        Parameter('gain', [
            Parameter('main', 0., float, 'Main output gain [dB]'),
            Parameter('aux', 0., float, 'Auxiliary output gain [dB]'),
            Parameter('main_invert', False, bool, 'Invert main channel gain'),
            Parameter('aux_invert', False, bool, 'Invert auxiliary channel gain'),
            Parameter('main_gain_range', '0dB', ['0dB', '14dB'], 'Main output gain range'),
            Parameter('aux_gain_range', '0dB', ['0dB', '14dB'], 'aux output gain range'),
            Parameter('strict', True, bool, 'Disable all implicit conversions and coercions')
        ])
    ])

    def __init__(self, name=None, settings=None):
        super().__init__(name, settings)
        self._instrument = instruments.LockInAmp(super().settings['ip_address'])
        self._is_connected = True

    def update(self, settings):
        if not self._settings_initialized:
            return

        for key, value in settings.items():
            if key == 'input':
                self.settings[key].update(value)
                self._instrument.set_frontend(self._LIA_INPUT_CHANNEL, **self.settings[key])
                # self._lia.set_frontend(self._LIA_INPUT_CHANNEL,
                #                        impedance=value['impedance'],
                #                        coupling=value['coupling'],
                #                        attenuation=value['attenuation'],
                #                        strict=value['strict'])
            elif key == 'output':
                self.settings[key].update(value)
                self._instrument.set_outputs(**self.settings[key])
                # self._lia.set_outputs(main=value['main'],
                #                       aux=value['aux'],
                #                       main_offset=value['main_offset'],
                #                       aux_offset=value['aux_offset'],
                #                       strict=value['strict'])
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

    @property
    def demod_frequency(self):
        return self._instrument.get_demodulation()['frequency']

    @demod_frequency.setter
    def demod_frequency(self, frequency):
        self.update({'demodulation': {'frequency': frequency}})

    @property
    def output_main(self):
        return self._instrument.get_outputs()['main']

    @output_main.setter
    def output_main(self, om):
        if om not in self.OUTPUT_MAIN_OPTIONS:
            raise KeyError('{} is not an allowed. Please choose one of {}'.format(om, self.OUTPUT_MAIN_OPTIONS))
        self.update({'output': {'main': om}})

    @property
    def output_aux(self):
        return self._instrument.get_outputs()['aux']

    @output_aux.setter
    def output_aux(self, om):
        if om not in self.OUTPUT_AUX_OPTIONS:
            raise KeyError('{} is not an allowed. Please choose one of {}'.format(om, self.OUTPUT_AUX_OPTIONS))
        self.update({'output': {'aux': om}})

    @property
    def output_aux_amplitude(self):
        return self._instrument.get_aux_output()['amplitude']

    @output_aux_amplitude.setter
    def output_aux_amplitude(self, oaa):
        self.update({'output_aux_amplitude': oaa})

