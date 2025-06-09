"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    pylabcontrol is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    pylabcontrol is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with pylabcontrol.  If not, see <http://www.gnu.org/licenses/>.
"""

from pylabcontrol.core import Instrument, Parameter
from struct import unpack
import pyvisa
import numpy as np
import time

import matplotlib.pyplot as plt


# class shutterlloscope(Instrument):

class IRShutter(Instrument):
    """
    Code for a Thorlabs shutter controller. This is connected to the computer via USB, and the Instrument
    interacts with the controller using PySerial and sending commands as defined in the controller documentation.
    """


    # String returned by spectrum analyzer upon querying it with '*IDN?'
    _INSTRUMENT_IDENTIFIER = 'THORLABS SC10 VERSION 1.01'

    _DEFAULT_SETTINGS = Parameter([
        Parameter('visa_resource', 'ASRL29::INSTR', str,
                  'pyVisa instrument identifier, to make a connection using the pyVisa package.'),
        Parameter('baudrate', 9600, int, 'baudrate of connection'),
        Parameter('timeout', 5, float, 'connection timeout'),
        Parameter('shutter_status', False, bool, "status of shutter")
    ])

    _COMMANDS = {
        'shutter_status': 'ens'
    }

    def __init__(self, name=None, settings=None):
        """

        Args:
            name (str): optional name of instance of class
            settings (list): list of other values to initialize class with

        """

        super(IRShutter, self).__init__(name, settings)
        print('the shutter is initialized')
        # keep track of when the instrument was updated last to prevent sending requests to frequently
        # self._last_update_time = time.time()

        rm = pyvisa.ResourceManager()

        # todo: JG 20170623 implement proper error handling when insturment is not connected.
        self.shutter = rm.open_resource(self.settings['visa_resource'])
        # self.shutter.query_delay = 1000

        self.shutter.write_termination = '\r'
        self.shutter.read_termination = '\r'

        self.shutter.timeout = self.settings['timeout']
        self.shutter.baud_rate = self.settings['baudrate']  # 9600

        self.update(self.settings)

    def _query(self, command):
        self.shutter.query(command, delay = 0.5)
        return self.shutter.read()

    def update(self, settings):
        """
        updates the instrument parameters

        Args:
            settings: dictionary that contains the parameter indentifiers (keys) and the new parameters values (value)

        """
        # Places the shutter in the factory default setup state.

        super(IRShutter, self).update(settings)
        if not self._settings_initialized:
            return

        for key, value in settings.items():
            # print(key)
            if key in ['shutter_status']:
                if value:
                    self.turn_on()
                if not value:
                    self.turn_off()
                # command = self._COMMANDS[key]
                # print(command)
                # self.shutter.write((command + ' ' + str(value)).encode())
                # self.shutter.query(command, delay = 0.5)
                # self.shutter.read()#, termination = '\r', encoding = 'ascii')


    @property
    def _PROBES(self):
        return {'shutter_status': 'status of shutter'}

    def read_probes(self, probe_name):

        # self._wait_for_shutter()
        # assert(False)
        # print(probe_name)
        assert self._settings_initialized  # will cause read_probes to fail if settings (and thus also connection) not yet initialized
        if probe_name == 'visa_resource':
            return self.settings['visa_resource']
        if probe_name == 'timeout':
            return self.settings['timeout']
        if probe_name == 'baudrate':
            return self.settings['baudrate']
        # print(probe_name)
        assert probe_name in list(self._PROBES.keys())

        output = self.get_status()

        return output

    def get_status(self):
        status = self._query('ens?')
        if status == '0':
            return False
        if status == '1':
            return True
        else:
            raise AssertionError

    def turn_on(self):
        """
        Checks if the shutter is open and opens it if not

        """
        status = self.get_status()
        if status:
            return
        else:
            print('turning on')
            self.shutter.query('ens', delay = 0.5)
        print(self.get_status())
        return

    def turn_off(self):
        """
        Checks if the shutter is open and opens it if not

        """
        status = self.get_status()
        if not status:
            return
        else:
            print('turning off')
            self.shutter.query('ens', delay = 0.5)
        print(self.get_status())
        return


    @property
    def is_connected(self):
        """
        Checks if the instrument is connected.
        Returns: True if connected, False otherwise.

        """
        identification = self.shutter.query("id?", delay = 0.5)
        identification = self.shutter.read()
        # identification = self.shutter.read_raw().decode()[:-1]
        # print(identification)
        return identification == self._INSTRUMENT_IDENTIFIER

if __name__ == '__main__':

        print('create shutter instance:')
        shut = IRShutter()
        if shut.is_connected:
            print('is connected')
        shut.turn_on()

        # print('=============')
        # # oscil.update({'waveform':{'CHAN1':{'vert_scale':0.5}}})
        #
        # # print((oscil.settings))
        #
        # # oscil.osci.write(':SINGLE')  # start a single acquisition
        # # oscil.osci.write(':TFORce')
        #
        # data, preambleBlock = oscil.get_timetrace()
        #
        # dt = preambleBlock['dt']
        # time = dt*np.arange(len(data))
        # print(('data', data))
        # plt.plot(time, data, '-x')
        #
        # plt.show()

