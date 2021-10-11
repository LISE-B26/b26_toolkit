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
import serial
from pylabcontrol.core import Instrument, Parameter

class ArduinoUno(Instrument):
    _DEFAULT_SETTINGS = Parameter([
        Parameter('port', 'COM5', str, 'COM port that arduino board is on'),
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
