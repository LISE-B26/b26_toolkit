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
import os
class CryoStation(Instrument):
    """
    instrument class to talk to get infos from Montana Cryostation
    Now this doesn't actually communicate with the Cryostation but only reads data from a log-file
    """
    _DEFAULT_SETTINGS = Parameter(
        Parameter('path', 'C:/Cryostation/Temperature Data/', str, 'path to log file of cryostation'),
    )


    def __init__(self, name = None, settings = None):

        super(CryoStation, self).__init__(name, settings)
        # apply all settings to instrument
        self.update(self.settings)

        try:
            # create available probes dynamically from headers of logfile
            filepath = "{:s}/MI_DiagnosticsDataLog {:s}.csv".format(self.settings['path'], time.strftime('%m_%d_%Y'))
            data = pd.read_csv(filepath)
            self._dynamic_probes = {
                elem.lstrip().lower().replace(' ', '_').replace('.', '').replace( ')', '').replace( '(', '').replace('/', '-'): elem
                for elem in data.columns}
            self._is_connected = True

        except IOError:
            self._is_connected = False
        except:
            raise ImportError

    def update(self, settings):
        '''
        updates the internal dictionary, just call function of super class

        '''
        super(CryoStation, self).update(settings)


    @property
    def is_connected(self):
        '''
        check if instrument is active and connected and return True in that case
        :return: bool
        '''
        return self._is_connected

    @property
    def _PROBES(self):
        '''

        Returns: a dictionary that contains the values that can be read from the instrument
        the key is the name of the value and the value of the dictionary is an info

        '''
        user_specific_probes = {
            'platform_temp': 'temperature of platform',
            'stage_1_temp': 'temperature of stage 1',
            'stage_2_temp': 'temperature of stage 2'
        }

        user_specific_probes.update(self._dynamic_probes)

        return user_specific_probes

    def read_probes(self, key):
        '''

        requests value from the instrument and returns it
        Args:
            key: name of requested value

        Returns: reads values from instrument

        '''

        key = key.lower()
        assert key in list(self._PROBES.keys()), "key assertion failed {:s}".format(str(key))



        # catch change of date
        time_tag = datetime.datetime.now()
        filepath = "{:s}/MI_DiagnosticsDataLog {:s}.csv".format(self.settings['path'],time_tag.strftime('%m_%d_%Y'))
        while os.path.exists(filepath) == False:
            time_tag -= datetime.timedelta(hours=1)
            filepath = "{:s}/MI_DiagnosticsDataLog {:s}.csv".format(self.settings['path'],time_tag.strftime('%m_%d_%Y'))

        data = pd.read_csv(filepath)

        # create dictionary with last row as values
        data = dict(data.iloc[-1])


        print(('xxxx', self._dynamic_probes))

        # since we striped some characters when defining the probes we have to find the right key,
        # which is give by the valeu in self._dynamic_probes
        key = self._dynamic_probes[key]



        return data[key]
