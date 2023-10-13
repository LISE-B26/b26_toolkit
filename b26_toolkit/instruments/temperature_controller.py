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


class TemperatureController(Instrument):
    """
    This class implements the Lakeshore Model 335 Temperature controller. The class communicates with the device over RS232 using pyserial.
    """


    # ASCII Characters used for controller communication
    ETX = chr(3)
    CR = chr(13)
    LF = chr(10)
    ENQ = chr(5)
    ACK = chr(6)
    NAK = chr(21)

    _possible_com_ports = ['COM' + str(i) for i in range(0, 256)]

    _DEFAULT_SETTINGS = Parameter([
            Parameter('port', 'COM5', _possible_com_ports, 'com port to which the gauge controller is connected'),
            Parameter('timeout', 2.0, float, 'amount of time to wait for a response '
                                             'from the gauge controller for each query'),
            Parameter('baudrate', 57600, int, 'baudrate of serial communication with gauge'),
            Parameter('temperature_A', 300, float, 'temperature reading A; changing this does nothing'),
            Parameter('temperature_B', 300, float, 'temperature reading B; changing this does nothing'),
            Parameter('heater_1_range', 'off', ['off', 'low', 'medium', 'high'], 'heater power range (each higher level -> 10x power)'),
            Parameter('heater_1_setpoint', 0, float, 'temperature setpoint (K) for PID feedback'),
            Parameter('heater_1_ramp_rate', 0, float, 'ramp rate (0.1 - 100 K/min) of temperature setpoint, to avoid temperature shocks during cool down/warm up'),
            Parameter('heater_1_output', 0, float, 'percent of full scale power for heater 1; changing this does nothing')
        ])

    #serial_connection = serial.Serial(port=_DEFAULT_SETTINGS['port'], baudrate=_DEFAULT_SETTINGS['baudrate'],
    #                                           timeout=_DEFAULT_SETTINGS['timeout'])
    def __init__(self, name='TemperatureController', settings=None):
        """
        The serial connection should be setup with the following parameters:
        1 start bit, 8 data bits, No parity bit, 1 stop bit, no hardware
        handshake. These are all default for Serial and therefore not input
        below
        """
        super(TemperatureController, self).__init__(name, settings)
        self.serial_connection = serial.Serial(port=self.settings['port'], baudrate=self.settings['baudrate'],
                                               timeout=self.settings['timeout'], parity=serial.PARITY_ODD, bytesize=serial.SEVENBITS)

    @property
    def _PROBES(self):
        """

        Returns: A dictionary of key-value string-string pairs. keys: probe names, values: probe descriptions

        """
        return {
            'temperature_A': 'numerical temperature read from temperature controller',
            'temperature_B': 'numerical temperature read from temperature controller',
            'heater_1_range': 'heater power range (each higher level -> 10x power)',
            'heater_1_setpoint': 'temperature setpoint for PID feedback',
            'heater_1_ramp_rate': 'temperature setpoint for PID feedback',
            'heater_1_output': 'percent of full scale power for heater 1'
        }

    def update(self, settings):
        super(TemperatureController, self).update(settings)

        for key, value in settings.items():
            if key == 'heater_1_range':
                value_dict = {'off': 0, 'low': 1, 'medium': 2, 'high': 3}
                self.serial_connection.write(('RANGE 1,%i \r\n' % value_dict[value]).encode())
            elif key == 'heater_1_setpoint':
                self.serial_connection.write(('SETP 1,%.3f \r\n' % float(value)).encode())
            elif key == 'heater_1_ramp_rate':
                self.serial_connection.write(('RAMP 1,1,%.3f \r\n' % float(value)).encode())


    def read_probes(self, probe_name):
        """
        Args:
            probe_name: Name of the probe to get the value of from the Pressure Gauge (e.g., 'pressure')

        Returns:
            value of the probe from the Pressure Gauge
        """

        assert probe_name in list(self._PROBES.keys())

        if probe_name == 'temperature_A':
            return self._get_temperature('A')
        if probe_name == 'temperature_B':
            return self._get_temperature('B')
        if probe_name == 'heater_1_range':
            self.serial_connection.write('RANGE? 1 \r\n'.encode())
            response = int(self.serial_connection.readline().strip())
            response_dict = {0: 'off', 1: 'low', 2: 'medium', 3: 'high'}
            return response_dict[response]
        if probe_name == 'heater_1_setpoint':
            self.serial_connection.write('SETP? 1 \r\n'.encode())
            return float(self.serial_connection.readline())
        if probe_name == 'heater_1_output':
            self.serial_connection.write('HTR? 1 \r\n'.encode())
            return float(self.serial_connection.readline())
        if probe_name == 'heater_1_ramp_rate':
            self.serial_connection.write('RAMP? 1 \r\n'.encode())
            response = (self.serial_connection.readline().strip()).decode("utf-8")
            status, rate = response.split(',')
            if status == 0:
                return 0
            else:
                return float(rate)

        else:
            message = '\'{0}\' not found as a probe in the class. ' \
                      'Expected either \'pressure\', \'units\', or \'model\''.format(probe_name)
            raise AttributeError(message)

    def _get_temperature(self, channel):
        """
        Returns the pressure currently read by the guage controller.

        :return: pressure
        """
        assert self.serial_connection.isOpen()
        self.serial_connection.write(('KRDG? %s\r\n' % str(channel)).encode())

        # response returns string with temperature values interspaced with \x characters
        response = self.serial_connection.readline()
        temperature = float(response[1:7])

        # ======================================
        # output when cold (>100K)
        # '\xab\xb6\xb0\xae\xb57\xb5\r\x8a'
        # '\xab\xb6\xb0\xae\xb57\xb3\r\x8a'
        # '\xab\xb6\xb0\xae\xb571\r\x8a'
        # '\xab\xb6\xb0\xae\xb5\xb6\xb9\r\x8a'
        # '\xab\xb6\xb0\xae\xb5\xb67\r\x8a'
        # ======================================

        # ======================================
        # output when cold (>100K)
        # '\xab\xb6\xb0\xae\xb57\xb5\r\x8a'
        # '\xab\xb6\xb0\xae\xb57\xb3\r\x8a'
        # '\xab\xb6\xb0\xae\xb571\r\x8a'
        # '\xab\xb6\xb0\xae\xb5\xb6\xb9\r\x8a'
        # '\xab\xb6\xb0\xae\xb5\xb67\r\x8a'
        # ======================================

        # repr(response) forces python not to interpret \x as a hex value
        # temperature = float(''.join([s for s in repr(response) if s.isdigit()])[0:-1])/1000.0

        return temperature

    def _set_heater_range(self, heater_range):
        #self.serial_connection.write('RANGE 1,%i \r\n'%heater_range.encode())
        pass

    def _set_heater_ramp_rate(self, ramp_rate):
        #self.serial_connection.write('RAMP 1,1,%.2f \r\n' % float(ramp_rate).encode())
        pass

    def _set_heater_setpoint(self, setpoint):
        #self.serial_connection.write('SETP 1,%.2f \r\n' % float(setpoint).encode())
        pass

    def is_connected(self):
        """
        checks if serial connection is still open with instrument.

        :return: boolean connection status
        """
        return self.serial_connection.isOpen()

    def __del__(self):
        """
        Destructor, to close the serial connection when the instance is this class is garbage collected
        """
        self.serial_connection.close()

if __name__ == '__main__':
        instruments, failed = Instrument.load_and_append(instrument_dict={'TemperatureController': TemperatureController})
        print(instruments['TemperatureController'].read_probes('temperature_A'))
        print(instruments['TemperatureController'].read_probes('heater_1_range'))
        print(instruments['TemperatureController'].read_probes('heater_1_ramp_rate'))

