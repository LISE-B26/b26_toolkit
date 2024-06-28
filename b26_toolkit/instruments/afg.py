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

import pyvisa
import pyvisa.errors
from pylabcontrol.core import Parameter, Instrument

# RANGE_MIN = 2025000000 #2.025 GHz

RANGE_MAX = 10. # Vpp, maximum power for the SRS IQ
FUNC_TYPE_TO_INTERNAL = {
    'Sine': 'SIN',
    'Square': 'SQU',
    'Pulse': 'PULS',
    'Sin(x)/x': 'SINC',
    'Noise': 'PRN',
    'DC': 'DC',
    'Gaussian': 'GAUS',
    'Lorentz': 'LOR',
    'Exponential Rise': 'ERIS',
    'Exponential Decay': 'EDEC',
    'Haversine': 'HAV'
}


class AFG3021C(Instrument):
    """
    This class implements the Tektronix AFG3021C. The class commuicates with the
    device over GPIB using pyvisa.
    """

    FREQ_MIN = 0
    FREQ_MAX = 25e6

    _DEFAULT_SETTINGS = Parameter([
        Parameter('USB_num', 0, int, 'USB device on which to connect'),
        Parameter('enable_output', False, bool, 'enable output'),
        # type of AWG function: should be Arb for
        Parameter('function', 'Square',
                  ['Sine', 'Square', 'Ramp', 'Pulse', 'Sin(x)/x', 'Noise', 'DC', 'Gaussian', 'Lorentz',
                   'Exponential Rise', 'Exponential Decay', 'Haversine'],
                  'Function'),

        # if a non-arb function is selected, program properties of the waveform
        # NOTE: needed for pulses too
        Parameter('frequency', 1.0e6, float, 'carrier frequency [Hz]'),
        Parameter('amplitude', 3.3, float, 'output amplitude peak to peak [Vpp]'),
        Parameter('phase', 0, float, 'output phase, in degrees'),
        Parameter('offset', 1.65, float, 'DC offset of waveform, [V]')
    ])

    MANUFACTURER_ID = '0x0699'
    MODEL_CODE = '0x0350'
    SERIAL_NUMBER = 'C020718'

    SIGNAL_MAX = 16382
    POINTS_MAX = 131027

    def __init__(self, name=None, settings=None):
        super().__init__(name, settings)
        # Issue where visa.ResourceManager() takes 4 minutes no longer happens after using pdb to debug (??? not sure why???)
        try:
            self._connect()
        except pyvisa.errors.VisaIOError:
            print('No AWG Detected!. Check that you are using the correct communication type')
            raise
        except Exception as e:
            raise e

    def __del__(self):
        self.afg.close()

    def _connect(self):
        rm = pyvisa.ResourceManager()
        self.afg = rm.open_resource(
            '::'.join(['USB' + str(self.settings['USB_num']), self.MANUFACTURER_ID, self.MODEL_CODE, self.SERIAL_NUMBER, 'INSTR']))
        self.afg.query('*IDN?')
        self.afg.write('SOUR1:FREQ:MODE CW')
        self.afg.write('SOUR1:VOLT:UNIT VPP')
        # at this point the analog output channels are still not enabled

    def update(self, settings):
        """
        Updates the internal settings of the MicrowaveGenerator, and then also updates physical parameters such as
        frequency, amplitude, modulation type, etc in the hardware
        Args:
            settings: a dictionary in the standard settings format
        """
        super().update(settings)
        # ===========================================
        for key, value in settings.items():
            if key != 'USB_num':
                if self.settings.valid_values[key] == bool: #converts booleans, which are more natural to store for on/off, to
                    arg = int(value)                #the integers used internally in the SRS
                elif key == 'function':
                    # if value == 'Arb':
                    #     self._waveform_to_memory(int(key[-1]))
                    arg = self._func_type_to_internal(value)
                else:
                    arg = value
                # elif (key == 'run_mode'):
                #     arg = self._run_mode_type_to_internal(key)
                cmd = self._param_to_internal(key)
                # only send update to instrument if connection to instrument has been established
                if self._settings_initialized:
                    self.settings[key] = value
                    self.afg.write(cmd + ' ' + str(arg)) # frequency change operation timed using timeit.timeit and
                                                           # completion confirmed by query('*OPC?'), found delay of <10ms

    @property
    def _PROBES(self):
        return {
            'enable_output': 'if output is enabled',
            'function': 'function of output signal',
            'frequency': 'carrier frequency of output[Hz]',
            'amplitude': 'amplitude of output [Vpp]',
            'phase': 'phase of output [degrees]',
            'offset': 'offset on output [V]'
        }

    def read_probes(self, key):
        assert(self._settings_initialized) #will cause read_probes to fail if settings (and thus also connection) not yet initialized
        assert key in list(self._PROBES.keys())

        #query always returns string, need to cast to proper return type
        key_internal = self._param_to_internal(key)
        value = self.afg.query(key_internal + '?')
        if key == 'enable_output':
            return bool(value)
        elif key == 'function':
            return self._internal_to_func_type(value)
        else:
            return float(value)

    @property
    def is_connected(self):
        try:
            self.afg.query('*IDN?') # arbitrary call to check connection, throws exception on failure to get response
            return True
        except pyvisa.errors.VisaIOError:
            return False

    def _func_type_to_internal(self, param):
        return FUNC_TYPE_TO_INTERNAL[param]

    def _internal_to_func_type(self, param):
        param = param.strip()
        for key in FUNC_TYPE_TO_INTERNAL:
            if FUNC_TYPE_TO_INTERNAL[key] == param:
                return key

    def _run_mode_type_to_internal(self, param):
        if param == 'run_mode_ch1':
            return 'SOUR1:BURS:STAT'

    def _param_to_internal(self, param):
        """
        Converts settings parameters to the corresponding key used for GPIB commands in the SRS.
        Args:
            param: settings parameter, ex. enable_output

        Returns: GPIB command, ex. ENBR
        """
        if param == 'enable_output':
            return 'OUTP1:STAT'
        elif param == 'function':
            return 'SOUR1:FUNC:SHAP'
        elif param == 'frequency':
            return 'SOUR1:FREQ'
        elif param == 'amplitude':
            return 'SOUR1:VOLT:LEV:IMM:AMPL'
        elif param == 'phase':
            return 'SOUR1:PHAS:ADJ'
        elif param == 'offset':
            return 'SOUR1:VOLT:LEV:IMM:OFFS'

        # Arbitrary waveform
        else:
            raise KeyError


if __name__ == '__main__':
    pass
