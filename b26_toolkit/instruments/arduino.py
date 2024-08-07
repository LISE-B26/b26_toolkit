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
from pylabcontrol.core import Instrument, Parameter
import time
import serial


class ArduinoUno(Instrument):
    _DEFAULT_SETTINGS = Parameter([
        Parameter('port', 'COM14', str, 'COM port that arduino board is on'),
        Parameter('baudrate', 57600, int, 'baudrate over serial port'),
        Parameter('timeout', 1., float, 'connection timeout in seconds')
    ])

    def __init__(self, name=None, settings=None):
        super().__init__(name, settings)
        self._is_connected = False

        try:
            self._connect(self.settings['port'],
                            self.settings['baudrate'],
                            self.settings['timeout'])
        except Exception as e:
            print(('Attocube not detected. Check connection.', UserWarning))
            print(e)
            raise e

    def __del__(self):
        if self._is_connected:
            self.ser.close()

    def _connect(self, port, baudrate, timeout):
        self.ser = serial.Serial(port=port,
                                 baudrate=baudrate,
                                 timeout=timeout)
        self._is_connected = True

    def flip(self, servoNumber, up):
        self.ser.write(bytearray([servoNumber, int(up)]))

    @property
    def _PROBES(self):
        return {}

    def read_probes(self, key):
        assert True

"""
A library to interface Arduino through serial connection
"""
import serial

class ArduinoZero(Instrument):
    _DEFAULT_SETTINGS = Parameter([
        Parameter('port', 'COM14', str, 'COM port that arduino board is on'),
        Parameter('baudrate', 9600, int, 'baudrate over serial port'),
        Parameter('timeout', 1., float, 'connection timeout in seconds'),
        Parameter('LB1005 integrator hold', [
            Parameter('channel', 13, int, 'channel to which LB1005 integrator hold toggle is connected; set true to disengage feedback loop'),
            Parameter('status', False, bool, 'True if voltage is high, false otherwise')
        ]),
        Parameter('pellicles', [
            Parameter('channel', 12, int, 'channel to which camera and LED pellicles are connected'),
            Parameter('status', False, bool, 'True if voltage is high, false otherwise')
        ]),
    ])

    def __init__(self, name=None, settings=None):
        super().__init__(name, settings)
        self._is_connected = False

        try:
            self._connect(self.settings['port'],
                            self.settings['baudrate'],
                            self.settings['timeout'])
        except Exception as e:
            print(('Arduino not detected. Check connection.', UserWarning))
            print(e)
            raise e

    def __del__(self):
        if self._is_connected:
            self.ser.close()

    def _connect(self, port, baudrate, timeout):
        self.ser = serial.Serial(port=port,
                                 baudrate=baudrate,
                                 timeout=timeout)
        self._is_connected = True

    def set_pin_mode(self, pin_number, mode):
        """
        Performs a pinMode() operation on pin_number
        Internally sends b'M{mode}{pin_number} where mode could be:
        - I for INPUT
        - O for OUTPUT
        - P for INPUT_PULLUP
        """
        command = (''.join(('M', mode, str(pin_number), '\n'))).encode()
        self.ser.write(command)
        time.sleep(0.5) # Wouldn't work without this

    def digital_write(self, pin_number, digital_value):
        """
        Writes the digital_value on pin_number
        """
        command = (''.join((str(1), ',', str(pin_number), ',', str(digital_value), '\n'))).encode()
        self.ser.write(command)


    @property
    def _PROBES(self):
        return {}

    def read_probes(self, key):
        assert True

    def update(self, settings):
        """
        Updates the internal settings of the MicrowaveGenerator, and then also updates physical parameters such as
        frequency, amplitude, modulation type, etc in the hardware
        Args:
            settings: a dictionary in the standard settings format
        """
        super(ArduinoZero, self).update(settings)

        for key, value in settings.items():
            if self._settings_initialized:
                #if key != 'port' and key != 'baudrate' and key != 'timeout':
                if type(value) is dict and 'status' in value.keys():
                    self.digital_write(int(self.settings[key]['channel']), int(value['status']))


if __name__ == '__main__':
        # instruments, failed = Instrument.load_and_append(instrument_dict={'GaugeController': PumpLinePressureGauge})


        # print((instruments['GaugeController']))
        # print(('PumpLinePressureGauge', instruments['GaugeController'].pressure))

        #instruments, failed = Instrument.load_and_append(instrument_dict={'arduino': ArduinoUno})
        a = ArduinoZero('hi')


        #a.digital_write(3, 0)
        for i in range(20):
            a.digital_write(13, 1)
            time.sleep(100e-3)
            a.digital_write(13, 0)
            time.sleep(100e-3)
