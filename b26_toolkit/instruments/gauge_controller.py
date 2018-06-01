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


class PressureGauge(Instrument):
    """
    This class implements the AGC100 pressure gauge. The class communicates with the device over RS232 using pyserial.
    """

    # Translations of the controller's status messages
    MEASUREMENT_STATUS = {
        '0': 'Measurement data okay',
        '1': 'Underrange',
        '2': 'Overrange',
        '3': 'Sensor error',
        '4': 'Sensor off',
        '5': 'No sensor',
        '6': 'Identification error',
        '7': 'Error FRG-720, FRG-730'
    }

    # Translation of the controller's units check  messages
    MEASUREMENT_UNITS = {
        '0': 'mbar/bar',
        '1': 'Torr',
        '2': 'Pascal',
        '3': 'Micron'
    }

    # ASCII Characters used for controller communication
    ETX = chr(3)
    CR = chr(13)
    LF = chr(10)
    ENQ = chr(5)
    ACK = chr(6)
    NAK = chr(21)

    _possible_com_ports = ['COM' + str(i) for i in range(0, 256)]

    _DEFAULT_SETTINGS = Parameter([
            Parameter('port', 'COM7', _possible_com_ports, 'com port to which the gauge controller is connected'),
            Parameter('timeout', 1.0, float, 'amount of time to wait for a response '
                                             'from the gauge controller for each query'),
            Parameter('baudrate', 9600, int, 'baudrate of serial communication with gauge')
        ])

    #serial_connection = serial.Serial(port=_DEFAULT_SETTINGS['port'], baudrate=_DEFAULT_SETTINGS['baudrate'],
    #                                           timeout=_DEFAULT_SETTINGS['timeout'])
    def __init__(self, name='PressureGauge', settings=None):
        """
        The serial connection should be setup with the following parameters:
        1 start bit, 8 data bits, No parity bit, 1 stop bit, no hardware
        handshake. These are all default for Serial and therefore not input
        below
        """

        super(PressureGauge, self).__init__(name, settings)
        self.serial_connection = serial.Serial(port=self.settings['port'], baudrate=self.settings['baudrate'],
                                               timeout=self.settings['timeout'])

    @property
    def _PROBES(self):
        """

        Returns: A dictionary of key-value string-string pairs. keys: probe names, values: probe descriptions

        """
        return {
            'pressure': 'numerical pressure read from Pressure Gauge',
            'units': 'Units used by pressure gauge',
            'model': 'Model of the pressure gauge'
        }

    def update(self, settings):
        super(PressureGauge, self).update(settings)

    def read_probes(self, probe_name):
        """
        Args:
            probe_name: Name of the probe to get the value of from the Pressure Gauge (e.g., 'pressure')

        Returns:
            value of the probe from the Pressure Gauge
        """

        probe_name = probe_name.lower()  # making sure the probe is lowercase

        if probe_name == 'pressure':
            return self._get_pressure()
        elif probe_name == 'units':
            return self._get_units()
        elif probe_name == 'model':
            return self._get_model()
        else:
            message = '\'{0}\' not found as a probe in the class. ' \
                      'Expected either \'pressure\', \'units\', or \'model\''.format(probe_name)
            raise AttributeError(message)

    def _check_acknowledgement(self, response):
        """
        _check_acknowledgement raises an error if the response passed in indicates an negatice response from the guage.

        :param response: the string response from the Guage Controller
        """

        if response == self.NAK + self.CR + self.LF:
            message = 'Serial communication returned negative acknowledge (NAK). ' \
                      'Check AGC100 documentation for more details.'
            raise IOError(message)

        elif response != self.ACK + self.CR + self.LF:
            message = 'Serial communication returned unknown response:\n{}' \
                ''.format(repr(response))
            raise AssertionError(message)

    def _get_pressure(self):
        """
        Returns the pressure currently read by the guage controller.

        :return: pressure
        """
        assert self.serial_connection.isOpen()

        self.serial_connection.write('PR1' + self.CR + self.LF)
        acknowledgement = self.serial_connection.readline()
        self._check_acknowledgement(acknowledgement)

        self.serial_connection.write(self.ENQ)
        err_msg_and_pressure = self.serial_connection.readline().rstrip(self.LF).rstrip(self.CR)

        err_msg = err_msg_and_pressure[0]
        pressure = float(err_msg_and_pressure[3:])

        if err_msg != '0':
            print(('xx', err_msg, pressure))
            message = 'Pressure query resulted in an error: ' + self.MEASUREMENT_STATUS[err_msg]
            # raise IOError(message) # JG: don't raise the error because this crashes the programm, rather we want to return an invalid value

        self.serial_connection.write(self.CR + self.LF)
        return pressure

    def _get_model(self):
        """
        Returns the model of the connected gauge controller.
        :return: model name
        """
        assert self.serial_connection.isOpen()

        self.serial_connection.write('TID' + self.CR + self.LF)
        acknowledgement = self.serial_connection.readline(25)
        self._check_acknowledgement(acknowledgement)

        self.serial_connection.write(self.ENQ)
        model = self.serial_connection.readline().rstrip(self.LF).rstrip(self.CR)

        self.serial_connection.write(self.CR + self.LF)

        return model

    def _get_units(self):
        """
        Returns the units that are in use by the guage controller.

        :return: gauge units (either bar, Torr, Pascal, or Micron)
        """
        #assert self.ser.isOpen()

        self.serial_connection.write('UNI' + self.CR + self.LF)
        acknowledgement = self.serial_connection.readline()
        self._check_acknowledgement(acknowledgement)

        self.serial_connection.write(self.ENQ)
        unit = self.MEASUREMENT_UNITS[self.serial_connection.readline().rstrip(self.LF).rstrip(self.CR)]

        self.serial_connection.write(self.CR + self.LF)

        return unit

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

class PumpLinePressureGauge(PressureGauge):
    """
    This class implements the AGC100 pressure gauge. The class communicates with the device over RS232 using pyserial.
    """

    _possible_com_ports = ['COM' + str(i) for i in range(0, 256)]

    _DEFAULT_SETTINGS = Parameter([
            Parameter('port', 'COM6', _possible_com_ports, 'com port to which the gauge controller is connected'),
            Parameter('timeout', 1.0, float, 'amount of time to wait for a response '
                                             'from the gauge controller for each query'),
            Parameter('baudrate', 9600, int, 'baudrate of serial communication with gauge')
        ])

class ChamberPressureGauge(PressureGauge):
    """
    This class implements the AGC100 pressure gauge. The class communicates with the device over RS232 using pyserial.
    """

    _possible_com_ports = ['COM' + str(i) for i in range(0, 256)]

    _DEFAULT_SETTINGS = Parameter([
            Parameter('port', 'COM7', _possible_com_ports, 'com port to which the gauge controller is connected'),
            Parameter('timeout', 1.0, float, 'amount of time to wait for a response '
                                             'from the gauge controller for each query'),
            Parameter('baudrate', 9600, int, 'baudrate of serial communication with gauge')
        ])



if __name__ == '__main__':
        instruments, failed = Instrument.load_and_append(instrument_dict={'GaugeController': PumpLinePressureGauge})


        print((instruments['GaugeController']))
        print(('PumpLinePressureGauge', instruments['GaugeController'].pressure))

        instruments, failed = Instrument.load_and_append(instrument_dict={'GaugeController': ChamberPressureGauge})
        print(('ChamberPressureGauge', instruments['GaugeController'].pressure))