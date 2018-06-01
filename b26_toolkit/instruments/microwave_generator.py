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

import visa
import pyvisa.errors

from pylabcontrol.core import Parameter, Instrument

# RANGE_MIN = 2025000000 #2.025 GHz
RANGE_MIN = 1012500000
RANGE_MAX = 4050000000 #4.050 GHZ

class MicrowaveGenerator(Instrument):
    """
    This class implements the Stanford Research Systems SG384 microwave generator. The class commuicates with the
    device over GPIB using pyvisa.
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('connection_type', 'RS232', ['GPIB', 'RS232'], 'type of connection to open to controller'),
        Parameter('port', 5, list(range(0, 31)), 'GPIB or COM port on which to connect'),
        Parameter('GPIB_num', 0, int, 'GPIB device on which to connect'),
        Parameter('enable_output', False, bool, 'Type-N output enabled'),
        Parameter('frequency', 3e9, float, 'frequency in Hz, or with label in other units ex 300 MHz'),
        Parameter('amplitude', -60, float, 'Type-N amplitude in dBm'),
        Parameter('phase', 0, float, 'output phase'),
        Parameter('enable_modulation', True, bool, 'enable modulation'),
        Parameter('modulation_type', 'FM', ['AM', 'FM', 'PhaseM', 'Freq sweep', 'Pulse', 'Blank', 'IQ'],
                  'Modulation Type: 0= AM, 1=FM, 2= PhaseM, 3= Freq sweep, 4= Pulse, 5 = Blank, 6=IQ'),
        Parameter('modulation_function', 'External', ['Sine', 'Ramp', 'Triangle', 'Square', 'Noise', 'External'],
                  'Modulation Function: 0=Sine, 1=Ramp, 2=Triangle, 3=Square, 4=Noise, 5=External'),
        Parameter('pulse_modulation_function', 'External', ['Square', 'Noise(PRBS)', 'External'], 'Pulse Modulation Function: 3=Square, 4=Noise(PRBS), 5=External'),
        Parameter('dev_width', 32e6, float, 'Width of deviation from center frequency in FM')
    ])

    def __init__(self, name=None, settings=None):

        super(MicrowaveGenerator, self).__init__(name, settings)

        # XXXXX MW ISSUE = START
        #===========================================
        # Issue where visa.ResourceManager() takes 4 minutes no longer happens after using pdb to debug (??? not sure why???)
        try:
            self._connect()
        except pyvisa.errors.VisaIOError:
            print('No Microwave Controller Detected!. Check that you are using the correct communication type')
            raise
        except Exception as e:
            raise(e)
        #XXXXX MW ISSUE = END
        #===========================================

    def _connect(self):
        rm = visa.ResourceManager()
        if self.settings['connection_type'] == 'GPIB':
            self.srs = rm.open_resource(
                'GPIB' + str(self.settings['GPIB_num']) + '::' + str(self.settings['port']) + '::INSTR')
        elif self.settings['connection_type'] == 'RS232':
            self.srs = rm.open_resource('COM' + str(self.settings['port']))
            self.srs.baud_rate = 115200
        self.srs.query('*IDN?')

    #Doesn't appear to be necessary, can't manually make two sessions conflict, rms may share well
    # def __del__(self):
    #     self.srs.close()

    def update(self, settings):
        """
        Updates the internal settings of the MicrowaveGenerator, and then also updates physical parameters such as
        frequency, amplitude, modulation type, etc in the hardware
        Args:
            settings: a dictionary in the standard settings format
        """
        super(MicrowaveGenerator, self).update(settings)
        # XXXXX MW ISSUE = START
        # ===========================================
        for key, value in settings.items():
            if key == 'connection_type':
                self._connect()
            elif not (key == 'port' or key == 'GPIB_num'):
                if self.settings.valid_values[key] == bool: #converts booleans, which are more natural to store for on/off, to
                    value = int(value)                #the integers used internally in the SRS
                elif key == 'modulation_type':
                    value = self._mod_type_to_internal(value)
                elif key == 'modulation_function':
                    value = self._mod_func_to_internal(value)
                elif key == 'pulse_modulation_function':
                    value = self._pulse_mod_func_to_internal
                # elif key == 'frequency':
                #     if value > RANGE_MAX or value < RANGE_MIN:
                #         raise ValueError("Invalid frequency. All frequencies must be between 2.025 GHz and 4.050 GHz.")
                key = self._param_to_internal(key)

                # only send update to instrument if connection to instrument has been established
                if self._settings_initialized:
                    self.srs.write(key + ' ' + str(value)) # frequency change operation timed using timeit.timeit and
                                                           # completion confirmed by query('*OPC?'), found delay of <10ms
                    # print(self.srs.query('*OPC?'))

        # XXXXX MW ISSUE = END
        # ===========================================

    @property
    def _PROBES(self):
        return{
            'enable_output': 'if type-N output is enabled',
            'frequency': 'frequency of output in Hz',
            'amplitude': 'type-N amplitude in dBm',
            'phase': 'phase',
            'enable_modulation': 'is modulation enabled',
            'modulation_type': 'Modulation Type: 0= AM, 1=FM, 2= PhaseM, 3= Freq sweep, 4= Pulse, 5 = Blank, 6=IQ',
            'modulation_function': 'Modulation Function: 0=Sine, 1=Ramp, 2=Triangle, 3=Square, 4=Noise, 5=External',
            'pulse_modulation_function': 'Pulse Modulation Function: 3=Square, 4=Noise(PRBS), 5=External',
            'dev_width': 'Width of deviation from center frequency in FM'
        }

    def read_probes(self, key):
        # assert hasattr(self, 'srs') #will cause read_probes to fail if connection not yet established, such as when called in init
        assert(self._settings_initialized) #will cause read_probes to fail if settings (and thus also connection) not yet initialized
        assert key in list(self._PROBES.keys())

        #query always returns string, need to cast to proper return type
        if key in ['enable_output', 'enable_modulation']:
            key_internal = self._param_to_internal(key)
            value = int(self.srs.query(key_internal + '?'))
            if value == 1:
                value = True
            elif value == 0:
                value = False
        elif key in ['modulation_type','modulation_function','pulse_modulation_function']:
            key_internal = self._param_to_internal(key)
            value = int(self.srs.query(key_internal + '?'))
            if key == 'modulation_type':
                value = self._internal_to_mod_type(value)
            elif key == 'modulation_function':
                value = self._internal_to_mod_func(value)
            elif key == 'pulse_modulation_function':
                value = self._internal_to_pulse_mod_func(value)
        else:
            key_internal = self._param_to_internal(key)
            value = float(self.srs.query(key_internal + '?'))

        return value

    @property
    def is_connected(self):
        try:
            self.srs.query('*IDN?') # arbitrary call to check connection, throws exception on failure to get response
            return True
        except pyvisa.errors.VisaIOError:
            return False

    def _param_to_internal(self, param):
        """
        Converts settings parameters to the corresponding key used for GPIB commands in the SRS.
        Args:
            param: settings parameter, ex. enable_output

        Returns: GPIB command, ex. ENBR

        """
        if param == 'enable_output':
            return 'ENBR'
        elif param == 'frequency':
            return 'FREQ'
        elif param == 'amplitude':
            return 'AMPR'
        elif param == 'phase':
            return 'PHAS'
        elif param == 'enable_modulation':
            return 'MODL'
        elif param == 'modulation_type':
            return 'TYPE'
        elif param == 'modulation_function':
            return 'MFNC'
        elif param == 'pulse_modulation_function':
            return 'PFNC'
        elif param == 'dev_width':
            return 'FDEV'
        else:
            raise KeyError

    def _mod_type_to_internal(self, value):
        #COMMENT_ME
        if value == 'AM':
            return 0
        elif value == 'FM':
            return 1
        elif value == 'PhaseM':
            return 2
        elif value == 'Freq sweep':
            return 3
        elif value == 'Pulse':
            return 4
        elif value == 'Blank':
            return 5
        elif value == 'IQ':
            return 6
        else:
            raise KeyError

    def _internal_to_mod_type(self, value):
        #COMMENT_ME
        if value == 0:
            return 'AM'
        elif value == 1:
            return 'FM'
        elif value == 2:
            return 'PhaseM'
        elif value == 3:
            return 'Freq sweep'
        elif value == 4:
            return 'Pulse'
        elif value == 5:
            return 'Blank'
        elif value == 6:
            return 'IQ'
        else:
            raise KeyError

    def _mod_func_to_internal(self, value):
        #COMMENT_ME
        if value == 'Sine':
            return 0
        elif value == 'Ramp':
            return 1
        elif value == 'Triangle':
            return 2
        elif value == 'Square':
            return 3
        elif value == 'Noise':
            return 4
        elif value == 'External':
            return 5
        else:
            raise KeyError

    def _internal_to_mod_func(self, value):
        #COMMENT_ME
        if value == 0:
            return 'Sine'
        elif value == 1:
            return 'Ramp'
        elif value == 2:
            return 'Triangle'
        elif value == 3:
            return 'Square'
        elif value == 4:
            return 'Noise'
        elif value == 5:
            return 'External'
        else:
            raise KeyError

    def _pulse_mod_func_to_internal(self, value):
        #COMMENT_ME
        if value == 'Square':
            return 3
        elif value == 'Noise(PRBS)':
            return 4
        elif value == 'External':
            return 5
        else:
            raise KeyError

    def _internal_to_pulse_mod_func(self, value):
        #COMMENT_ME
        if value == 3:
            return 'Square'
        elif value == 4:
            return 'Noise(PRBS)'
        elif value == 5:
            return 'External'
        else:
            raise KeyError

if __name__ == '__main__':
    # from pylabcontrol.core import Instrument
    #
    # instruments =        {"MicrowaveGenerator": {
    #         "class": "MicrowaveGenerator",
    #         "settings": {
    #             "enable_output": False,
    #             "enable_modulation": True,
    #             "amplitude": -60,
    #             "GPIB_num": 0,
    #             "frequency": 3000000000.0,
    #             "dev_width": 32000000.0,
    #             "pulse_modulation_function": "External",
    #             "phase": 0,
    #             "modulation_function": "External",
    #             "port": 27,
    #             "modulation_type": "FM"
    #         }
    #     }}
    #
    # instr = {}
    #
    # instr, loaded_failed = Instrument.load_and_append(
    #     {name: instruments[name] for name in instruments.keys()}, instr)
    #
    # print(instr)
    # print(loaded_failed)
    # import pylabcontrol.instruments.microwave_generator MicrowaveGenerator

    mw = MicrowaveGenerator()
    # mw = MicrowaveGenerator(settings={'enable_modulation': True, 'frequency': 3000000000.0, 'dev_width': 32000000.0, 'pulse_modulation_function': 'External', 'phase': 0, 'port': 27, 'modulation_type': 'FM', 'enable_output': False, 'GPIB_num': 0, 'amplitude': -60, 'modulation_function': 'External'})

    print((mw.srs))

    # instrument_name= 'MicrowaveGenerator'
    # instrument_settings = {'enable_modulation': True, 'frequency': 3000000000.0, 'dev_width': 32000000.0, 'pulse_modulation_function': 'External', 'phase': 0, 'port': 27, 'modulation_type': 'FM', 'enable_output': False, 'GPIB_num': 0, 'amplitude': -60, 'modulation_function': 'External'}
    # class_of_instrument = MicrowaveGenerator
    # instrument_instance = class_of_instrument(name=instrument_name, settings=instrument_settings)
