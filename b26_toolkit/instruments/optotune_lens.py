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

MAX_CURRENT = 292.84 #mA, used as a scaling factor for setting current


class OptotuneLens(Instrument):
    """
    Instrument class to control an Optotune Lens Driver 4. Tested with an EL-10-30-TC, but any elens compatable with
    the Lens Driver 4 should work.
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('port', 'COM20', str, 'serial port on which to connect'),
        Parameter('baudrate', 115200, int, 'baudrate of connection'),
        Parameter('timeout', .1, float, 'connection timeout'),
        Parameter('current', 0.0, float, 'current applied to lens')
    ])

    def __init__(self, name = None, settings = None):
        """
        Initializes connection to optotune lens. If none found, raises exception.
        Args:
            name: instrument name
            settings: dictionary of settings to override defaults
        """
        super(OptotuneLens, self).__init__(name, settings)
        self._is_connected = False
        try:
            self.connect(port = self.settings['port'], baudrate = self.settings['baudrate'], timeout = self.settings['timeout'])
        except Exception:
            print('No ELens Detected')
            raise

    def connect(self, port, baudrate, timeout):
        """
        Connects to elens using the serial port
        Args:
            port: COM port on which to connect
            baudrate: baud rate of connection. Check value required by device
            timeout: time to wait before abandoning communication

        Poststate: self._is_connected is True

        """
        self.ser = serial.Serial(port = port, baudrate = baudrate, timeout = timeout)
        self._is_connected = True

    def update(self, settings):
        """
        Updates internal settings, and sets ELens current on hardware if that is changed
        Args:
            settings: dictionary of settings to update

        Poststate: changes current on ELens if it is updated

        """
        super(OptotuneLens, self).update(settings)
        for key, value in settings.items():
            if self._settings_initialized:
                if key == 'current':
                    self.set_current(value)

    @property
    def _PROBES(self):
        """
        Returns: a dictionary that contains the values that can be read from the instrument
        the key is the name of the value and the value of the dictionary is an info
        """
        return {
            'current': 'the current applied to the ELens'
        }

    def read_probes(self, key):
        """
        requests value from the instrument and returns it
        Args:
            key: name of requested value

        Returns: reads values from instrument

        """
        assert key in list(self._PROBES.keys())
        assert isinstance(key, str)

        if key in ['current']:
            self.ser.write(bytearray.fromhex('41 72 00 00 b4 27'))
            response = bytearray(self.ser.readline())
            current_scaled = int.from_bytes(response[1:3], 'big')
            return((current_scaled/4096) * 292.84)

    @property
    def is_connected(self):
        """

        Returns: True if currently connected to piezo controller, False otherwise

        """
        try:
            self.current
            return True
        except serial.serialutil.SerialTimeoutException:
            return False

    def __del__(self):
        """
        Ensures that connection to hardware is closed on deletion of PiezoController object, to prevent a stale
        connection from a closed object from blocking further connections from new objects
        """
        if self._is_connected:
            self.ser.close()

    def set_current(self, current_mA):
        """
        Sets the voltage on the piezo.
        Args:
            voltage: voltage (in V) to set

        """
        # self.ser.write((self.settings['axis'] + 'voltage=' + str(voltage) + '\r').encode('ascii'))
        if(current_mA > MAX_CURRENT):
            self.log('cannot set current above ' + str(MAX_CURRENT))
            return
        current_scaled = int(current_mA / MAX_CURRENT * 4096)
        writebytes = bytearray((b'Aw' + current_scaled.to_bytes(2, 'big')))
        checksum = CRC16IBM().compute_checksum_bytes(writebytes)
        self.ser.write(writebytes + checksum)
        error = self.ser.readline()
        if error:
            self.log('setting of current failed with error code ' + error)

class CRC16IBM():
    polynomial = 0xA001
    table = [None] * 256

    def __init__(self):
        for i in range(len(self.table)):
            value = 0
            temp = i
            for j in range(8):
                if(((value ^ temp) & 0x0001) != 0):
                    value = ((value >> 1) ^ self.polynomial)
                else:
                    value >>= 1
                temp >>= 1
            self.table[i] = value

    def compute_checksum(self, bytes):
        crc = 0
        for i in range(len(bytes)):
            index = (crc^bytes[i]) & 0xFF
            crc = (crc >> 8) ^ self.table[index]
        return crc

    def compute_checksum_bytes(self, bytes):
        crc = self.compute_checksum(bytes)
        return crc.to_bytes(2, 'little')

# crc = CRC16IBM()
# temp = crc.compute_checksum_bytes(bytearray.fromhex('41 72 00 00'))
# print(hex(temp[0]))
# print(hex(temp[1]))

# lens = OptotuneLens()
# print(lens.current)
# lens.set_current(0)
# print(lens.current)