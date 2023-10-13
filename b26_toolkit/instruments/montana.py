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
import time, datetime
import pandas as pd
import os, tailer

TEMP_INDICES = {'platform_temp': 3, 'stage_1_temp': 5, 'stage_2_temp': 6}

class CryoStation(Instrument):
    """
    instrument class to talk to get infos from Montana Cryostation
    Now this doesn't actually communicate with the Cryostation but only reads data from a log-file
    """
    _DEFAULT_SETTINGS = Parameter([
        Parameter('path', 'C:/Cryostation/Temperature Data/', str, 'path to log file of cryostation'),
    ])


    def __init__(self, name='CryoStation', settings=None):

        super().__init__(name, settings)
        self._today = time.strftime('%m_%d_%Y')

        # create available probes dynamically from headers of logfile
        self._filepath = "{:s}/MI_DiagnosticsDataLog {:s}.csv".format(self.settings['path'], self._today)
        self._is_connected = os.path.isfile(self._filepath)

    @property
    def _PROBES(self):
        '''

        Returns: a dictionary that contains the values that can be read from the instrument
        the key is the name of the value and the value of the dictionary is an info

        '''
        return {
            'platform_temp': 'temperature of platform',
            'stage_1_temp': 'temperature of stage 1',
            'stage_2_temp': 'temperature of stage 2'
        }

    def read_probes(self, key):
        '''

        requests value from the instrument and returns it
        Args:
            key: name of requested value

        Returns: reads values from instrument

        '''

        key = key.lower()
        assert key in list(self._PROBES.keys()), "key assertion failed {:s}".format(str(key))

        todayTemp = time.strftime('%m_%d_%Y')
        newFilepath = "{:s}/MI_DiagnosticsDataLog {:s}.csv".format(self.settings['path'], todayTemp)
        if todayTemp != self._today and os.path.isfile(newFilepath):
                self._today = todayTemp
                self._filepath = newFilepath

        try:
            with open(self._filepath) as file:
                row = tailer.tail(file, 1)[1]
                return float(row.split(', ')[TEMP_INDICES[key]])
        except:
            raise;

if __name__ == '__main__':
    instruments, failed = Instrument.load_and_append(instrument_dict={'CryoStation': CryoStation})
    print((instruments['CryoStation'].platform_temp))