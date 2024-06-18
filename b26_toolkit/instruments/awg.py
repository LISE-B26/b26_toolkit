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

import pyvisa as visa
import pyvisa.errors
import numpy as np
from pylabcontrol.core import Parameter, Instrument
from b26_toolkit.instruments import Pulse


# RANGE_MIN = 2025000000 #2.025 GHz
RANGE_MIN = -3.4/2-0.1 # V, minimum voltage for the SRS IQ (should be -1, but we now have attenuators)
RANGE_MAX = 3.4/2+0.1 #V, maximum voltage for the SRS IQ (should be +1, but we now have attenuators)

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

class AFG3022C(Instrument):
    """
    This class implements the Tektronix AFG3022C. The class commuicates with the
    device over GPIB using pyvisa.

    It seems like pylabcontrol doesn't like nested instrument settings, so I flattened everything like in Grozny
    """

    _DEFAULT_SETTINGS = Parameter([
      #  Parameter('connection_type', 'GPIB', ['GPIB', 'RS232'], 'type of connection to open to controller'),
        #Parameter('port', 11, list(range(0,31)), 'GPIB port on which to connect'),
        Parameter('USB_num', 0, int, 'GPIB device on which to connect'),

        # if a non-arb function is selected, program properties of the waveform
        # NOTE: needed for pulses too

        Parameter('ch1_enable', False, bool, 'enable output CH1'),
        Parameter('ch1_function', 'Arb',
                  ['Sine', 'Square', 'Ramp', 'Pulse', 'Arb', 'Sin(x)/x', 'Noise', 'DC', 'Gaussian', 'Lorentz',
                   'Exponential Rise', 'Exponential Decay', 'Haversine'],
                  'Function CH1'),
        Parameter('ch1_run_mode', 'Continuous', ['Continuous', 'Burst', 'Unimplemented'],
                  'run mode of CH1. "Unimplemented" means the AFG is either in sweep or modulation mode; '
                  'choosing "unimplemented" will simply set it to "continuous"'),
        Parameter('ch1_frequency', 1e6, float, 'frequency in Hz'),
        Parameter('ch1_amplitude', 0.5, float, 'output amplitude peak to peak in volts'),
        Parameter('ch1_phase', 0.0, float, 'output phase in radians'),
        Parameter('ch1_offset', 0.0, float, 'DC offset of waveform, in volts'),
        Parameter('ch1_invert_polarity', False, bool, 'invert polarity of output'),
        Parameter('ch1_burst_mode', 'Trigger', ['Trigger', 'Gate'], 'Select burst mode, only applies if "Burst" is selected in ch1_mode. '
                                                                    'AFG is triggered on positive edge of external signal.'),
        Parameter('ch1_burst_ncycles', 1, int, 'number of cycles in each burst, anything above 1e6 will be interpreted as "infinite" by the AFG'),
        Parameter('ch1_pulse_duty_cycle', 50, float, 'duty cycle (%%) of pulse waveform setting range is 0.001%% to 99.999%% in increments of 0.001'),

        Parameter('ch2_enable', False, bool, 'enable output CH1'),
        Parameter('ch2_function', 'Arb',
                  ['Sine', 'Square', 'Ramp', 'Pulse', 'Arb', 'Sin(x)/x', 'Noise', 'DC', 'Gaussian', 'Lorentz',
                   'Exponential Rise', 'Exponential Decay', 'Haversine'],
                  'Function CH1'),
        Parameter('ch2_run_mode', 'Continuous', ['Continuous', 'Burst', 'Unimplemented'],
                  'run mode of CH2. "Unimplemented" means the AFG is either in sweep or modulation mode; '
                  'choosing "unimplemented" will simply set it to "continuous"'),
        Parameter('ch2_frequency', 1e6, float, 'frequency in Hz'),
        Parameter('ch2_amplitude', 0.5, float, 'output amplitude peak to peak in volts'),
        Parameter('ch2_phase', 0.0, float, 'output phase in radians'),
        Parameter('ch2_offset', 0.0, float, 'DC offset of waveform, in volts'),
        Parameter('ch2_invert_polarity', False, bool, 'invert polarity of output'),
        Parameter('ch2_burst_mode', 'Trigger', ['Trigger', 'Gate'], 'Select burst mode, only applies if "Burst" is selected in ch1_mode. '
                                                                    'AFG is triggered on positive edge of external signal.'),
        Parameter('ch2_burst_ncycles', 1, int, 'number of cycles in each burst, anything above 1e6 will be interpreted as "infinite" by the AFG'),
        Parameter('ch2_pulse_duty_cycle', 50, float, 'duty cycle (%%) of pulse waveform setting range is 0.001%% to 99.999%% in increments of 0.001'),
        ])

    MANUFACTURER_ID = '0x0699'
    MODEL_CODE = '0x034A'
    SERIAL_NUMBER = 'C020934'

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
        #self.afg.write('SOUR2:FREQ:MODE CW')
        #self.afg.write('SOUR1:FREQ:MODE CW')
        self.afg.write('SOUR2:VOLT:UNIT VPP')
        self.afg.write('SOUR1:VOLT:UNIT VPP')

        # Configure burst mode (not the same as activating it)
        self.afg.write('TRIG:SOUR EXT')
        self.afg.write('TRIG:SLOP POS')
        # at this point the analog output channels are still not enabled

    def update(self, settings):
        """
        Updates the internal settings of the MicrowaveGenerator, and then also updates physical parameters such as
        frequency, amplitude, modulation type, etc in the hardware
        Args:
            settings: a dictionary in the standard settings format
        """
        super(AFG3022C, self).update(settings)
        # ===========================================
        for key, value in settings.items():
            if key != 'USB_num':
                channel = key[2]
                if key[4:] == 'enable':
                    value = int(value)
                elif key[4:] == 'function':
                    value = self._func_type_to_internal(value, key)
                elif key[4:] == 'frequency':
                    if value > 25e6 or value < 0:
                        raise ValueError("Invalid frequency. All frequencies must be between 0 and 25 MHz.")
                    value = str(value)+'Hz'
                elif key[4:] == 'phase':
                    value = value % (2*np.pi)
                elif key[4:] == 'amplitude':
                    if value > RANGE_MAX or value < RANGE_MIN:
                        raise ValueError("Invalid amplitude. All amplitudes must be between -0.5 and +0.5V.")
                elif key[4:] == 'offset':
                    if (value > RANGE_MAX or value < RANGE_MIN):
                        raise ValueError("All voltages programmed on the AWG must be between -0.5 and +0.5V")
                elif key[4:] == 'invert_polarity':
                    if value:
                        value = 'INV'
                    else:
                        value = 'NORM'
                elif key[4:] == 'burst_mode':
                    if value == 'Trigger':
                        value = 'TRIG'
                    elif value == 'Gate':
                        value = 'GAT'
                elif key[4:] == 'burst_ncycles':
                    if value > 1e6:
                        value = 'INF'
                elif key[4:] == 'pulse_duty_cycle':
                    if value > 99.999:
                        value = 99.999
                    elif value < 0.001:
                        value = 0.001
                elif key[4:] == 'run_mode':
                    # We need a different code structure here bc there's no simple command in the API to change the run mode
                    # Have to write a different command depending on the requested mode
                    if value == 'Continuous':
                        key = key[:4]+'set_to_continuous'
                        value = 'CW'
                    elif value == 'Burst':
                        key = key[:4]+'set_to_burst'
                        value = 'ON'
                    elif value == 'Sweep':
                        key = key[:4]+'set_to_sweep'
                        value = 'SWE'
                    elif value == 'Unimplemented':
                        key = key[:4] + 'set_to_continuous'
                        value = 'CW'


                key = self._param_to_internal(key[4:], int(channel))

                # only send update to instrument if connection to instrument has been established
                if self._settings_initialized:
                    self.afg.write(
                        key + ' ' + str(value))  # frequency change operation timed using timeit.timeit and
                        # completion confirmed by query('*OPC?'), found delay of <10ms

    def align_phase(self):
        self.afg.write('SOUR1:PHAS:INIT')

    @property
    def _PROBES(self):
        return{
            'ch1_enable':'enable output CH1',
            'ch1_function': 'function CH1',
            'ch1_run_mode': 'run mode of CH1',
            'ch1_frequency': 'frequency in Hz',
            'ch1_amplitude': 'output amplitude peak to peak in volts',
            'ch1_phase': 'output phase in radians',
            'ch1_offset': 'DC offset of waveform, in volts',
            'ch1_invert_polarity': 'invert polarity of output',
            'ch1_burst_mode': '',
            'ch1_burst_ncycles': '',
            'ch1_pulse_duty_cycle': '',
            'ch2_enable': 'enable output CH2',
            'ch2_function': 'Function CH1',
            'ch2_frequency': 'frequency in Hz',
            'ch2_amplitude': 'output amplitude peak to peak in volts',
            'ch2_phase': 'output phase in radians',
            'ch2_offset': 'DC offset of waveform, in volts',
            'ch2_invert_polarity': 'invert polarity of output'}

    def read_probes(self, key):
        # assert hasattr(self, 'awg') #will cause read_probes to fail if connection not yet established, such as when called in init
        assert self._settings_initialized #will cause read_probes to fail if settings (and thus also connection) not yet initialized
        assert key in list(self._PROBES.keys())

        # Query always returns string, need to cast to proper return type
        channel = int(key[2])
        if key[4:] == 'run_mode':
            # Check each mode one by one
            def check_mode(param):
                # Query the AFG for a given parameter
                key_internal = self._param_to_internal(param, channel)
                value = self.afg.query(key_internal + '?').strip()
                return value
            if check_mode('set_to_continuous') ==  'CW':
                value = 'Continuous'
            elif bool(int(check_mode('set_to_burst'))):
                value = 'Burst'
            else:
                value = 'Unimplemented'
        else:
            key_internal = self._param_to_internal(key[4:], channel)
            if key[4:] == 'enable':
                #print(self.afg.query(key_internal + '?').strip())
                value = bool(int(self.afg.query(key_internal + '?').strip()))
            elif key[4:] in ['function', 'invert_polarity', 'burst_mode']:
                value = self._internal_to_param(self.afg.query(key_internal + '?').strip())
            elif key[4:] == 'burst_ncycles':
                value = float(self.afg.query(key_internal + '?').strip())
                if value > 1000000:
                    value = 1000001
                else:
                    value = int(value)
            else:
                value = float(self.afg.query(key_internal + '?').strip())
        return value

    @property
    def is_connected(self):
        try:
            self.afg.query('*IDN?') # arbitrary call to check connection, throws exception on failure to get response
            return True
        except pyvisa.errors.VisaIOError:
            return False

    def _func_type_to_internal(self, param, channel=1):
        if channel == 1 and param == 'Arb':
            return 'USER1'
        elif channel == 2 and param == 'Arb':
            return 'USER2'
        elif param == 'Sine':
            return 'SIN'
        elif param == 'Square':
            return 'SQU'
        elif param == 'Ramp':
            return 'RAMP'
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

    def _internal_to_param(self, value):
        conversion_dict = {
            'USER1': 'Arb',
            'USER2': 'Arb',
            'SIN': 'Sine',
            'SQU': 'Square',
            'RAMP': 'Ramp',
            'PULS': 'Pulse',
            'SINC': 'Sin(x)/x',
            'PRN': 'Noise',
            'DC': 'DC',
            'GAUS': 'Gaussian',
            'LOR': 'Lorentz',
            'ERIS': 'Exponential Rise',
            'EDEC': 'Exponential Decay',
            'HAV': 'Haversine',
            'NORM': False,
            'INV': True,
            'TRIG': 'Trigger',
            'GAT': 'Gate'
        }
        return conversion_dict[value]

    def _param_to_internal(self, param, channel=1):
        """
        Converts settings parameters to the corresponding key used for GPIB commands in the SRS.
        Args:
            param: settings parameter, ex. enable_output

        Returns: GPIB command, ex. ENBR

        """
        #print(param)
        if param == 'enable':
            if channel == 1:
                return 'OUTP1:STAT'
            elif channel == 2:
                return 'OUTP2:STAT'
        elif param == 'function':
            if channel == 1:
                return 'SOUR1:FUNC:SHAP'
            elif channel == 2:
                return 'SOUR2:FUNC:SHAP'
        elif param == 'frequency':
            if channel == 1:
                return 'SOUR1:FREQ:CW'
            elif channel == 2:
                return 'SOUR2:FREQ:CW'
        elif param == 'amplitude':
            if channel == 1:
                return 'SOUR1:VOLT:LEV:IMM:AMPL'
            elif channel == 2:
                return 'SOUR2:VOLT:LEV:IMM:AMPL'
        elif param == 'phase':
            if channel == 1:
                return 'SOUR1:PHAS:ADJ'
            elif channel == 2:
                return 'SOUR2:PHAS:ADJ'
        elif param == 'offset':
            if channel == 1:
                return 'SOUR1:VOLT:LEV:IMM:OFFS'
            elif channel == 2:
                return 'SOUR2:VOLT:LEV:IMM:OFFS'
        elif param == 'invert_polarity':
            if channel == 1:
                return 'OUTP1:POL'
            elif channel == 2:
                return 'OUTP2:POL'
        elif param == 'set_to_continuous' or param == 'set_to_sweep':
            return 'SOUR%i:FREQ:MODE' % channel
        elif param == 'set_to_burst':
            return 'SOUR%i:BURS:STAT' % channel
        elif param == 'burst_mode':
            return 'SOUR%i:BURS:MODE' % channel
        elif param == 'burst_ncycles':
            return 'SOUR%i:BURS:NCYC' % channel
        elif param == 'pulse_duty_cycle':
            return 'SOUR%i:PULS:DCYC' % channel

        # Arbitrary waveform
        else:
            raise KeyError

    def configure_iq(self, angle, pulse_length, angle_idle=0):
        """
        Configures the pulse amplitudes of the AFG to modify the IQ inputs to the SRS. We assume that channel 1&2 are connected to the SRS i&q respectively.
        Args:
            angle: rotation angle in radians of the MW, defined as arctan(q amplitude/i amplitude)
            pulse_length: length of pulse output in ns by AFG
            angle_idle: rotation angle in radians when AFG is not triggered (AFG can output an idle voltage when pulse waveform is selected)
        Returns:
        """
        magnitude = 3.4/2  # Match the PB high-level so that we don't have to recalibrate pi-pulses when switching between using just the PB to using AFG+PB
        voltages = magnitude*np.cos(angle), magnitude*np.sin(angle)
        voltages_idle = magnitude * np.cos(angle_idle), magnitude * np.sin(angle_idle)
        duty_cycle = 50

        for channel, voltage, voltage_idle in zip([1, 2], voltages, voltages_idle):
            self.update({'ch%i_function' % channel: 'Pulse'})
            self.update({'ch%i_run_mode' % channel: 'Burst'})
            self.update({'ch%i_frequency' % channel: float(1 / (pulse_length*1e-9 / (duty_cycle / 100)))})
            self.update({'ch%i_burst_ncycles' % channel: 1})
            self.update({'ch%i_pulse_duty_cycle' % channel: duty_cycle})

            if voltage > voltage_idle:
                self.update({'ch%i_invert_polarity' % channel: False})
            else:
                self.update({'ch%i_invert_polarity' % channel: True})
            self.update({'ch%i_offset' % channel: float(voltage + voltage_idle)/2})

            amp = float(np.abs(voltage - voltage_idle))
            if amp > 2e-4: # AFG cannot output less than 0.1 mV
                self.update({'ch%i_amplitude' % channel: float(np.abs(voltage - voltage_idle))})
            else:
                self.update({'ch%i_function' % channel: 'DC'})


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

class AFG3022C_02(AFG3022C):
    """
    Copy of AFG3022C, just so that we can use two AFGs at the same time (fun quirk of pylabcontrol)
    """



    MANUFACTURER_ID = '0x0699'
    MODEL_CODE = '0x034A'
    SERIAL_NUMBER = 'C020007'

    SIGNAL_MAX = 16382
    POINTS_MAX = 131027


class AFG3021C(AFG3022C):
    """
    Copy of AFG3022C, just so that we can use two AFGs at the same time (fun quirk of pylabcontrol)
    """

    MANUFACTURER_ID = '0x0699'
    MODEL_CODE = '0x0350'
    SERIAL_NUMBER = 'C020718'

    SIGNAL_MAX = 16382
    POINTS_MAX = 131027

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


if __name__ == '__main__':
    #arb_iq = AFG3022C_02()
    #print(arb_iq.afg)
    #print(arb_iq)
    #arb.update({'ch1_function': 'Pulse'})
    #arb.update({'ch1_run_mode': 'Continuous'})
    #arb.update({'ch1_burst_mode': 'Trigger'})
    #arb.update({'ch1_burst_ncycles': 1})
    #arb_iq.configure_iq(angle=np.pi, pulse_length=1e-6, angle_idle=0)
    #for angle in np.linspace(0,2*np.pi,100):
    #    arb_iq.configure_iq(angle=angle, pulse_length=286, angle_idle=0)
    #print(AFG3022C_02._DEFAULT_SETTINGS.valid_values['ch1_run_mode'])

    arb = AFG3021C()
    #print(arb)
    #print(arb.afg.query('*RST').strip())
