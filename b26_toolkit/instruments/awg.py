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
import numpy as np
import matplotlib.pyplot as plt

from pylabcontrol.core import Parameter, Instrument
#from pulse_blaster import Pulse
from b26_toolkit.instruments import Pulse


# RANGE_MIN = 2025000000 #2.025 GHz
RANGE_MIN = -0.500 # V, minimum voltage for the SRS IQ
RANGE_MAX = 0.500 #V, maximum voltage for the SRS IQ

RANGE_MAX = 10. #Vpp, maximum power for the SRS IQ

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

class AWG(Instrument): # Emma Rosenfeld 20170822
    """
    This class implements the Tektronix AFG3022C. The class commuicates with the
    device over GPIB using pyvisa.
    """

    _DEFAULT_SETTINGS = Parameter([
      #  Parameter('connection_type', 'GPIB', ['GPIB', 'RS232'], 'type of connection to open to controller'),
        #Parameter('port', 11, list(range(0,31)), 'GPIB port on which to connect'),
        Parameter('USB_num', 0, int, 'GPIB device on which to connect'),
        Parameter('enable_output_ch1', False, bool, 'enable output CH1'),
        Parameter('enable_output_ch2', False, bool, 'enable output CH2'),

        # type of AWG function: should be Arb for
        Parameter('function_ch1', 'Arb',
                  ['Sine', 'Square', 'Ramp', 'Pulse', 'Arb', 'Sin(x)/x', 'Noise', 'DC', 'Gaussian', 'Lorentz',
                   'Exponential Rise', 'Exponential Decay', 'Haversine'],
                  'Function CH1'),
        Parameter('function_ch2', 'Arb',
                  ['Sine', 'Square', 'Ramp', 'Pulse', 'Arb', 'Sin(x)/x', 'Noise', 'DC', 'Gaussian', 'Lorentz',
                   'Exponential Rise', 'Exponential Decay', 'Haversine'],
                  'Function CH2'),

        # if a non-arb function is selected, program properties of the waveform
        # NOTE: needed for pulses too
        Parameter('default_waveform_ch1', [
            Parameter('frequency', 1e6, float, 'frequency in Hz'),
            Parameter('amplitude', 0.5, float, 'output amplitude peak to peak in volts'),
            Parameter('phase', 0, float, 'output phase, in degrees'),
            Parameter('offset', 0, float, 'DC offset of waveform, in volts')
        ]),
        Parameter('default_waveform_ch2', [
            Parameter('frequency', 1e6, float, 'frequency in Hz'),
            Parameter('amplitude', 0.5, float, 'output amplitude peak to peak in volts'),
            Parameter('phase', 0, float, 'output phase'),
            Parameter('offset', 0, float, 'DC offset of waveform, in volts')
        ]),

        # # Arbitrary waveforms
        # Parameter('arbitrary_waveform_ch1', [
        #     Parameter('time', np.zeros([1]), np.ndarray, '1D array of time values in seconds'),
        #     Parameter('amplitude', np.zeros([1]), np.ndarray, '1D array of amplitude values in volts')
        # ]),
        #
        # Parameter('arbitrary_waveform_ch2', [
        #     Parameter('time', np.zeros([1]), np.ndarray, '1D array of time values in seconds'),
        #     Parameter('amplitude', np.zeros([1]), np.ndarray, '1D array of amplitude values in volts')
        # ]),

        # Parameter('pulse_sequence_ch1', [
        #     Parameter('period', 1e6, float, 'time between  pulses in seconds'),
        #     Parameter('pi_pulse_width', 1e6, float, 'pulse width in seconds for pi pulse'),
        #     Parameter('half_pi_pulse_width', 1e6, float, 'pulse width in seconds for half pi pulse'),
        #     Parameter('x_amplitude', 0.5, float, 'pulse height in volts for x pulse'),
        #     Parameter('y_amplitude', 0.5, float, 'pulse height in volts for y pulse'),
        #     Parameter('sequence', [], )
        # ])

        # Parameter('run_mode_ch1', 'Burst', ['Continuous', 'Burst'], 'Run Mode: continuous or burst once'),
        # Parameter('run_mode_ch2', 'Burst', ['Continuous', 'Burst'], 'Run Mode: continuous or burst once'),
    ])

    MANUFACTURER_ID = '0x0699'
    MODEL_CODE = '0x034A'
    SERIAL_NUMBER = 'C020007'

    SIGNAL_MAX = 16382
    POINTS_MAX = 131027

    def __init__(self, name=None, settings=None, pulse_sequence_ch1=None, pulse_sequence_ch2=None):

        if pulse_sequence_ch1 is not None:
            self.pulse_sequence_ch1 = pulse_sequence_ch1
        if pulse_sequence_ch2 is not None:
            self.pulse_sequence_ch2 = pulse_sequence_ch2
        super(AWG, self).__init__(name, settings)

        #===========================================
        # Issue where visa.ResourceManager() takes 4 minutes no longer happens after using pdb to debug (??? not sure why???)
        try:
            self._connect()
        except pyvisa.errors.VisaIOError:
            print('No AWG Detected!. Check that you are using the correct communication type')
            raise
        except Exception as e:
            raise(e)
        #===========================================

        #self.update(self.settings)

    def _connect(self):
        rm = visa.ResourceManager()
        self.awg = rm.open_resource(
            '::'.join(['USB' + str(self.settings['USB_num']), self.MANUFACTURER_ID, self.MODEL_CODE, self.SERIAL_NUMBER, 'INSTR']))
        self.awg.query('*IDN?')
        self.awg.write('*RST')  # reset AWG
        self.awg.write('SOUR1:BURS:MODE TRIG')
        self.awg.write('SOUR2:BURS:MODE TRIG')
        self.awg.write('SOUR1:BURS:NCYC 1') # prepare settings for burst mode
        self.awg.write('SOUR2:BURS:NCYC 1')
        self.awg.write('OUTP1:TRIG:MODE EXT')
        self.awg.write('OUTP2:TRIG:MODE EXT')
        self.awg.write('SOUR1:VOLT:UNIT VPP')
        self.awg.write('SOUR2:VOLT:UNIT VPP')
        self.awg.write('SOUR1:VOLT:LEV:IMM:AMPL 1V')
        self.awg.write('SOUR2:VOLT:LEV:IMM:AMPL 1V')
        # at this point the analog output channels are still not enabled

    def update(self, settings):
        """
        Updates the internal settings of the MicrowaveGenerator, and then also updates physical parameters such as
        frequency, amplitude, modulation type, etc in the hardware
        Args:
            settings: a dictionary in the standard settings format
        """
        #print(self.settings)
        super(AWG, self).update(settings)
        # ===========================================
        for key, value in settings.items():
            if key != 'USB_num' and type(value) is not dict:
                if self.settings.valid_values[key] == bool: #converts booleans, which are more natural to store for on/off, to
                    arg = int(value)                #the integers used internally in the SRS
                elif (key == 'function_ch1' or key == 'function_ch2'):
                    if value == 'Arb':
                        self._waveform_to_memory(int(key[-1]))
                    arg = self._func_type_to_internal(value, key)
                elif (key == 'run_mode_ch1' or key == 'run_mode_ch2'):
                    arg = self._run_mode_type_to_internal(key)
                elif key == 'amplitude':
                    if value > RANGE_MAX or value < RANGE_MIN:
                        raise ValueError("Invalid amplitude. All amplitudes must be between -0.5 and +0.5V.")
                cmd = self._param_to_internal(key)
                # only send update to instrument if connection to instrument has been established
                if self._settings_initialized:
                    self.settings[key] = value
                    self.awg.write(cmd + ' ' + str(arg)) # frequency change operation timed using timeit.timeit and
                                                           # completion confirmed by query('*OPC?'), found delay of <10ms
            elif type(value) is dict: # if nested dictionary keep iterating
                for key2, value2 in value.items():
                    if key == 'default_waveform_ch1' or key == 'default_waveform_ch2':
                        if key2 == 'phase':
                            value2 = np.deg2rad(np.mod(value2, 360) - 180)
                        elif key2 == 'amplitude':
                            if value2 > RANGE_MAX or value2 < RANGE_MIN:
                                raise ValueError("Invalid amplitude. All amplitudes must be between -0.5 and +0.5V.")
                        elif key2 == 'offset':
                            if (value2 > 0.5 or value2 < -0.5):
                                raise ValueError("All voltages programmed on the AWG must be between -0.5 and +0.5V")

                    key2 = self._param_to_internal(key2, key)
                    # only send update to instrument if connection to instrument has been established
                    if self._settings_initialized:
                        self.awg.write(
                            key2 + ' ' + str(value2))  # frequency change operation timed using timeit.timeit and
                        # completion confirmed by query('*OPC?'), found delay of <10ms

                        # print(self.awg.query('*OPC?'))
                    # elif (key == 'arbitrary_waveform_ch1' and settings['function_ch1'] == 'Arb') \
                    #         or (key == 'arbitrary_waveform_ch2' and settings['function_ch2'] == 'Arb'):
                    #     if key2 == 'time' or key2 == 'amplitude' and (type(value2) is not np.ndarray or len(value2.shape) != 1):
                    #         raise ValueError('Time is not a 1D array.')

        # ===========================================

    @property
    def _PROBES(self):
        return{
        }

    def read_probes(self, key):

        # assert hasattr(self, 'awg') #will cause read_probes to fail if connection not yet established, such as when called in init
        assert(self._settings_initialized) #will cause read_probes to fail if settings (and thus also connection) not yet initialized
        assert key in list(self._PROBES.keys())

        #query always returns string, need to cast to proper return type
        if key in ['enable_output_ch1']:
            key_internal = self._param_to_internal(key, 0)
            value = int(self.awg.query(key_internal + '?'))
            if value == 1:
                value = True
            elif value == 0:
                value = False
        else:
            key_internal = self._param_to_internal(key, 0)
            value = float(self.awg.query(key_internal + '?'))

        return value

    @property
    def is_connected(self):
        try:
            self.awg.query('*IDN?') # arbitrary call to check connection, throws exception on failure to get response
            return True
        except pyvisa.errors.VisaIOError:
            return False

    def _func_type_to_internal(self, param, key0=0):
        if (key0 == 'function_ch1' and param == 'Arb'):
            return 'USER1'
        elif (key0 == 'function_ch2' and param == 'Arb'):
            return 'USER2'
        elif param == 'Sine':
            return 'SIN'
        elif param == 'Square':
            return 'SQU'
        elif param == 'Pulse':
            return 'PULS'
        elif param == 'Sin(x)/x':
            return 'SINC'
        elif param == 'Noise':
            return 'PRN'
        elif param == 'DC':
            return 'DC'
        elif param == 'Gaussian':
            return 'GAUS'
        elif param == 'Lorentz':
            return 'LOR'
        elif param == 'Exponential Rise':
            return 'ERIS'
        elif param == 'Exponential Decay':
            return 'EDEC'
        elif param == 'Haversine':
            return 'HAV'

    def _run_mode_type_to_internal(self, param):
        if param == 'run_mode_ch1':
            return 'SOUR1:BURS:STAT'
        elif param == 'run_mode_ch2':
            return 'SOUR2:BURS:STAT'

    def _param_to_internal(self, param, key0=0):
        """
        Converts settings parameters to the corresponding key used for GPIB commands in the SRS.
        Args:
            param: settings parameter, ex. enable_output

        Returns: GPIB command, ex. ENBR

        """
        if param == 'enable_output_ch1':
            return 'OUT1:STAT'
        elif param == 'enable_output_ch2':
            return 'OUT2:STAT'
        elif param == 'function_ch1':
            return 'SOUR1:FUNC:SHAP'
        elif param == 'function_ch2':
            return 'SOUR1:FUNC:SHAP'
        elif param == 'frequency':
            if key0 == 'default_waveform_ch1':
                return 'SOUR1:FREQ'
            elif key0 == 'default_waveform_ch2':
                return 'SOUR2:FREQ'
        elif param == 'amplitude':
            if key0 == 'default_waveform_ch1':
                return 'SOUR1:VOLT:LEV:IMM:AMPL'
            elif key0 == 'default_waveform_ch2':
                return 'SOUR2:VOLT:LEV:IMM:AMPL'
        elif param == 'phase':
            if key0 == 'default_waveform_ch1':
                return 'SOUR1:PHAS:ADJ'
            elif key0 == 'default_waveform_ch2':
                return 'SOUR2:PHAS:ADJ'
        elif param == 'offset':
            if key0 == 'default_waveform_ch1':
                return 'SOUR1:VOLT:LEV:IMM:OFFS'
            elif key0 == 'default_waveform_ch2':
                return 'SOUR2:VOLT:LEV:IMM:OFFS'

        # Arbitrary waveform
        else:
            raise KeyError

    def _waveform_to_memory(self, channel):
        if self._settings_initialized:
            if channel == 1:
                pulse_sequence = self.__pulse_sequence_ch1
            elif channel == 2:
                pulse_sequence = self.__pulse_sequence_ch2
            else:
                raise ValueError('There are only two channels for the Tektronix AFG3022C.')

            self.awg.write('DATA:DEF EMEM,' + str(self.POINTS_MAX))  # Reset edit memory
            for i in range(1, pulse_sequence.shape[1]):
                self.awg.write('DATA:DATA:LINE EMEM,' + ','.join([str(el) for el in [pulse_sequence[0, i - 1],
                                                                  pulse_sequence[1, i - 1],
                                                                  pulse_sequence[0, i],
                                                                  pulse_sequence[1, i]]]))

            self.awg.write('DATA:COPY USER' + str(channel) + ',EMEM')

    def _pulse_to_points(self, pulse_sequence):
        time = np.array([0.])
        amplitude = np.array([0.])

        for pulse in pulse_sequence:
            if type(pulse) is not Pulse:
                raise ValueError('List is not a sequence of Pulse objects.')
            if pulse.amplitude is None:
                raise ValueError('Amplitude not defined.')
            time = np.append(time, [pulse.start_time, pulse.start_time, pulse.end_time, pulse.end_time])
            amplitude = np.append(amplitude, [0., pulse.amplitude, pulse.amplitude, 0.])

        amplitude *= self.SIGNAL_MAX / np.max(amplitude)
        time *= self.POINTS_MAX / np.max(time)

        return np.array([np.floor(time).astype(int), np.floor(amplitude).astype(int)])

    @property
    def pulse_sequence_ch1(self):
        return self.__pulse_sequence_ch1

    @pulse_sequence_ch1.setter
    def pulse_sequence_ch1(self, pulse_sequence):
        self.__pulse_sequence_ch1 = self._pulse_to_points(pulse_sequence)

    @property
    def pulse_sequence_ch2(self):
        return self.__pulse_sequence_ch2

    @pulse_sequence_ch2.setter
    def pulse_sequence_ch2(self, pulse_sequence):
        self.__pulse_sequence_ch2 = self._pulse_to_points(pulse_sequence)


class AFG3022C_GARBAGE(Instrument):
    """
    This class implements the Tektronix AFG3021C. The class commuicates with the
    device over GPIB using pyvisa.
    """

    FREQ_MIN = 0
    FREQ_MAX = 25e6

    _DEFAULT_SETTINGS = Parameter([
      #  Parameter('connection_type', 'GPIB', ['GPIB', 'RS232'], 'type of connection to open to controller'),
        #Parameter('port', 11, list(range(0,31)), 'GPIB port on which to connect'),
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


    """ DaLi settings
    MANUFACTURER_ID = '0x0699'
    MODEL_CODE = '0x0350'
    SERIAL_NUMBER = 'C020718' 
    """

    MANUFACTURER_ID = '0x0699'
    MODEL_CODE = '0x034A'
    SERIAL_NUMBER = 'C020007'

    SIGNAL_MAX = 16382
    POINTS_MAX = 131027


    def __init__(self, name=None, settings=None):


        super().__init__(name, settings)

        #===========================================
        # Issue where visa.ResourceManager() takes 4 minutes no longer happens after using pdb to debug (??? not sure why???)
        try:
            self._connect()
        except pyvisa.errors.VisaIOError:
            print('No AWG Detected!. Check that you are using the correct communication type')
            raise
        except Exception as e:
            raise(e)
        #===========================================

        #self.update(self.settings)

    def __del__(self):
        self.afg.close()

    def _connect(self):
        rm = pyvisa.ResourceManager()
        self.afg = rm.open_resource(
            '::'.join(['USB' + str(self.settings['USB_num']), self.MANUFACTURER_ID, self.MODEL_CODE, self.SERIAL_NUMBER, 'INSTR']))
        self.afg.query('*IDN?')
        self.afg.write('SOUR2:FREQ:MODE CW')
        self.afg.write('SOUR1:FREQ:MODE CW')
        self.afg.write('SOUR2:VOLT:UNIT VPP')
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

        # assert hasattr(self, 'awg') #will cause read_probes to fail if connection not yet established, such as when called in init
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
        # if (key0 == 'function_ch1' and param == 'Arb'):
        #     return 'USER1'
        # elif (key0 == 'function_ch2' and param == 'Arb'):
        #     return 'USER2'
        return FUNC_TYPE_TO_INTERNAL[param]

    def _internal_to_func_type(self, param):
        # if (key0 == 'function_ch1' and param == 'Arb'):
        #     return 'USER1'
        # elif (key0 == 'function_ch2' and param == 'Arb'):
        #     return 'USER2'
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


class AFG3022C(Instrument): # Emma Rosenfeld 20170822
    """
    This class implements the Tektronix AFG3022C. The class commuicates with the
    device over GPIB using pyvisa.
    """

    _DEFAULT_SETTINGS = Parameter([
      #  Parameter('connection_type', 'GPIB', ['GPIB', 'RS232'], 'type of connection to open to controller'),
        #Parameter('port', 11, list(range(0,31)), 'GPIB port on which to connect'),
        Parameter('USB_num', 0, int, 'GPIB device on which to connect'),
        Parameter('enable_output_ch1', False, bool, 'enable output CH1'),
        Parameter('enable_output_ch2', False, bool, 'enable output CH2'),

        # type of AWG function: should be Arb for
        Parameter('function_ch1', 'Arb',
                  ['Sine', 'Square', 'Ramp', 'Pulse', 'Arb', 'Sin(x)/x', 'Noise', 'DC', 'Gaussian', 'Lorentz',
                   'Exponential Rise', 'Exponential Decay', 'Haversine'],
                  'Function CH1'),
        Parameter('function_ch2', 'Arb',
                  ['Sine', 'Square', 'Ramp', 'Pulse', 'Arb', 'Sin(x)/x', 'Noise', 'DC', 'Gaussian', 'Lorentz',
                   'Exponential Rise', 'Exponential Decay', 'Haversine'],
                  'Function CH2'),

        # if a non-arb function is selected, program properties of the waveform
        # NOTE: needed for pulses too
        Parameter('default_waveform_ch1', [
            Parameter('frequency', 1e6, float, 'frequency in Hz'),
            Parameter('amplitude', 0.5, float, 'output amplitude peak to peak in volts'),
            Parameter('phase', 0, float, 'output phase, in degrees'),
            Parameter('offset', 0, float, 'DC offset of waveform, in volts')
        ]),
        Parameter('default_waveform_ch2', [
            Parameter('frequency', 1e6, float, 'frequency in Hz'),
            Parameter('amplitude', 0.5, float, 'output amplitude peak to peak in volts'),
            Parameter('phase', 0, float, 'output phase'),
            Parameter('offset', 0, float, 'DC offset of waveform, in volts')
        ]),

        # # Arbitrary waveforms
        # Parameter('arbitrary_waveform_ch1', [
        #     Parameter('time', np.zeros([1]), np.ndarray, '1D array of time values in seconds'),
        #     Parameter('amplitude', np.zeros([1]), np.ndarray, '1D array of amplitude values in volts')
        # ]),
        #
        # Parameter('arbitrary_waveform_ch2', [
        #     Parameter('time', np.zeros([1]), np.ndarray, '1D array of time values in seconds'),
        #     Parameter('amplitude', np.zeros([1]), np.ndarray, '1D array of amplitude values in volts')
        # ]),

        # Parameter('pulse_sequence_ch1', [
        #     Parameter('period', 1e6, float, 'time between  pulses in seconds'),
        #     Parameter('pi_pulse_width', 1e6, float, 'pulse width in seconds for pi pulse'),
        #     Parameter('half_pi_pulse_width', 1e6, float, 'pulse width in seconds for half pi pulse'),
        #     Parameter('x_amplitude', 0.5, float, 'pulse height in volts for x pulse'),
        #     Parameter('y_amplitude', 0.5, float, 'pulse height in volts for y pulse'),
        #     Parameter('sequence', [], )
        # ])

        # Parameter('run_mode_ch1', 'Burst', ['Continuous', 'Burst'], 'Run Mode: continuous or burst once'),
        # Parameter('run_mode_ch2', 'Burst', ['Continuous', 'Burst'], 'Run Mode: continuous or burst once'),
    ])

    MANUFACTURER_ID = '0x0699'
    MODEL_CODE = '0x034A'
    SERIAL_NUMBER = 'C020007'

    SIGNAL_MAX = 16382
    POINTS_MAX = 131027


    def __init__(self, name=None, settings=None):

        super().__init__(name, settings)

        #===========================================
        # Issue where visa.ResourceManager() takes 4 minutes no longer happens after using pdb to debug (??? not sure why???)
        try:
            self._connect()
        except pyvisa.errors.VisaIOError:
            print('No AWG Detected!. Check that you are using the correct communication type')
            raise
        except Exception as e:
            raise(e)
        #===========================================

        #self.update(self.settings)

    def _connect(self):
        rm = visa.ResourceManager()
        self.afg = rm.open_resource(
            '::'.join(['USB' + str(self.settings['USB_num']), self.MANUFACTURER_ID, self.MODEL_CODE, self.SERIAL_NUMBER, 'INSTR']))
        self.afg.query('*IDN?')
        self.afg.write('SOUR2:FREQ:MODE CW')
        self.afg.write('SOUR1:FREQ:MODE CW')
        self.afg.write('SOUR2:VOLT:UNIT VPP')
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
            if key != 'USB_num' and type(value) is not dict:
                if self.settings.valid_values[key] == bool: #converts booleans, which are more natural to store for on/off, to
                    arg = int(value)                #the integers used internally in the SRS
                elif (key == 'function_ch1' or key == 'function_ch2'):
                    if value == 'Arb':
                        self._waveform_to_memory(int(key[-1]))
                    arg = self._func_type_to_internal(value, key)
                elif (key == 'run_mode_ch1' or key == 'run_mode_ch2'):
                    arg = self._run_mode_type_to_internal(key)
                elif key == 'amplitude':
                    if value > RANGE_MAX or value < RANGE_MIN:
                        raise ValueError("Invalid amplitude. All amplitudes must be between -0.5 and +0.5V.")
                cmd = self._param_to_internal(key)
                # only send update to instrument if connection to instrument has been established
                if self._settings_initialized:
                    self.settings[key] = value
                    self.afg.write(cmd + ' ' + str(arg)) # frequency change operation timed using timeit.timeit and
                                                           # completion confirmed by query('*OPC?'), found delay of <10ms
            elif type(value) is dict: # if nested dictionary keep iterating
                for key2, value2 in value.items():
                    if key == 'default_waveform_ch1' or key == 'default_waveform_ch2':
                        if key2 == 'phase':
                            value2 = np.deg2rad(np.mod(value2, 360) - 180)
                        elif key2 == 'amplitude':
                            if value2 > RANGE_MAX or value2 < RANGE_MIN:
                                raise ValueError("Invalid amplitude. All amplitudes must be between -0.5 and +0.5V.")
                        elif key2 == 'offset':
                            if (value2 > 0.5 or value2 < -0.5):
                                raise ValueError("All voltages programmed on the AWG must be between -0.5 and +0.5V")

                    key2 = self._param_to_internal(key2, key)
                    # only send update to instrument if connection to instrument has been established
                    if self._settings_initialized:
                        self.afg.write(
                            key2 + ' ' + str(value2))  # frequency change operation timed using timeit.timeit and
                        # completion confirmed by query('*OPC?'), found delay of <10ms

                        # print(self.awg.query('*OPC?'))
                    # elif (key == 'arbitrary_waveform_ch1' and settings['function_ch1'] == 'Arb') \
                    #         or (key == 'arbitrary_waveform_ch2' and settings['function_ch2'] == 'Arb'):
                    #     if key2 == 'time' or key2 == 'amplitude' and (type(value2) is not np.ndarray or len(value2.shape) != 1):
                    #         raise ValueError('Time is not a 1D array.')

        # ===========================================

    @property
    def _PROBES(self):
        return{
        }

    def read_probes(self, key):

        # assert hasattr(self, 'awg') #will cause read_probes to fail if connection not yet established, such as when called in init
        assert(self._settings_initialized) #will cause read_probes to fail if settings (and thus also connection) not yet initialized
        assert key in list(self._PROBES.keys())

        #query always returns string, need to cast to proper return type
        if key in ['enable_output_ch1']:
            key_internal = self._param_to_internal(key, 0)
            value = int(self.afg.query(key_internal + '?'))
            if value == 1:
                value = True
            elif value == 0:
                value = False
        else:
            key_internal = self._param_to_internal(key, 0)
            value = float(self.afg.query(key_internal + '?'))

        return value

    @property
    def is_connected(self):
        try:
            self.afg.query('*IDN?') # arbitrary call to check connection, throws exception on failure to get response
            return True
        except pyvisa.errors.VisaIOError:
            return False

    def _func_type_to_internal(self, param, key0=0):
        if (key0 == 'function_ch1' and param == 'Arb'):
            return 'USER1'
        elif (key0 == 'function_ch2' and param == 'Arb'):
            return 'USER2'
        elif param == 'Sine':
            return 'SIN'
        elif param == 'Square':
            return 'SQU'
        elif param == 'Pulse':
            return 'PULS'
        elif param == 'Sin(x)/x':
            return 'SINC'
        elif param == 'Noise':
            return 'PRN'
        elif param == 'DC':
            return 'DC'
        elif param == 'Gaussian':
            return 'GAUS'
        elif param == 'Lorentz':
            return 'LOR'
        elif param == 'Exponential Rise':
            return 'ERIS'
        elif param == 'Exponential Decay':
            return 'EDEC'
        elif param == 'Haversine':
            return 'HAV'

    def _run_mode_type_to_internal(self, param):
        if param == 'run_mode_ch1':
            return 'SOUR1:BURS:STAT'
        elif param == 'run_mode_ch2':
            return 'SOUR2:BURS:STAT'

    def _param_to_internal(self, param, key0=0):
        """
        Converts settings parameters to the corresponding key used for GPIB commands in the SRS.
        Args:
            param: settings parameter, ex. enable_output

        Returns: GPIB command, ex. ENBR

        """
        if param == 'enable_output_ch1':
            return 'OUT1:STAT'
        elif param == 'enable_output_ch2':
            return 'OUT2:STAT'
        elif param == 'function_ch1':
            return 'SOUR1:FUNC:SHAP'
        elif param == 'function_ch2':
            return 'SOUR1:FUNC:SHAP'
        elif param == 'frequency':
            if key0 == 'default_waveform_ch1':
                return 'SOUR1:FREQ'
            elif key0 == 'default_waveform_ch2':
                return 'SOUR2:FREQ'
        elif param == 'amplitude':
            if key0 == 'default_waveform_ch1':
                return 'SOUR1:VOLT:LEV:IMM:AMPL'
            elif key0 == 'default_waveform_ch2':
                return 'SOUR2:VOLT:LEV:IMM:AMPL'
        elif param == 'phase':
            if key0 == 'default_waveform_ch1':
                return 'SOUR1:PHAS:ADJ'
            elif key0 == 'default_waveform_ch2':
                return 'SOUR2:PHAS:ADJ'
        elif param == 'offset':
            if key0 == 'default_waveform_ch1':
                return 'SOUR1:VOLT:LEV:IMM:OFFS'
            elif key0 == 'default_waveform_ch2':
                return 'SOUR2:VOLT:LEV:IMM:OFFS'

        # Arbitrary waveform
        else:
            raise KeyError

    def _waveform_to_memory(self, channel):
        if self._settings_initialized:
            if channel == 1:
                pulse_sequence = self.__pulse_sequence_ch1
            elif channel == 2:
                pulse_sequence = self.__pulse_sequence_ch2
            else:
                raise ValueError('There are only two channels for the Tektronix AFG3022C.')

            self.afg.write('DATA:DEF EMEM,' + str(self.POINTS_MAX))  # Reset edit memory
            for i in range(1, pulse_sequence.shape[1]):
                self.afg.write('DATA:DATA:LINE EMEM,' + ','.join([str(el) for el in [pulse_sequence[0, i - 1],
                                                                  pulse_sequence[1, i - 1],
                                                                  pulse_sequence[0, i],
                                                                  pulse_sequence[1, i]]]))

            self.afg.write('DATA:COPY USER' + str(channel) + ',EMEM')

    def _pulse_to_points(self, pulse_sequence):
        time = np.array([0.])
        amplitude = np.array([0.])

        for pulse in pulse_sequence:
            if type(pulse) is not Pulse:
                raise ValueError('List is not a sequence of Pulse objects.')
            if pulse.amplitude is None:
                raise ValueError('Amplitude not defined.')
            time = np.append(time, [pulse.start_time, pulse.start_time, pulse.end_time, pulse.end_time])
            amplitude = np.append(amplitude, [0., pulse.amplitude, pulse.amplitude, 0.])

        amplitude *= self.SIGNAL_MAX / np.max(amplitude)
        time *= self.POINTS_MAX / np.max(time)

        return np.array([np.floor(time).astype(int), np.floor(amplitude).astype(int)])

    @property
    def pulse_sequence_ch1(self):
        return self.__pulse_sequence_ch1

    @pulse_sequence_ch1.setter
    def pulse_sequence_ch1(self, pulse_sequence):
        self.__pulse_sequence_ch1 = self._pulse_to_points(pulse_sequence)

    @property
    def pulse_sequence_ch2(self):
        return self.__pulse_sequence_ch2

    @pulse_sequence_ch2.setter
    def pulse_sequence_ch2(self, pulse_sequence):
        self.__pulse_sequence_ch2 = self._pulse_to_points(pulse_sequence)


if __name__ == '__main__':
    arb = AFG3022C()
    print(arb.afg)
    arb.update({'default_waveform_ch1': {'frequency': (2.9294e6-3e3)}})
