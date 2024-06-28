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

from pylabcontrol.core import Parameter, Instrument


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


class Commander(Instrument):
    """
    This class does not implement a physical instrument. However, we can tell scripts to check the settings of this ''instrument'' and conditionally do things
    like FindNV, run autofocus, etc.

    For example, in every iteration of the loop inside ESR, it can check if the Commander's setting for ''run_findnv'' is checked, and it will then run it
    before the next iteration.

    It's basically a hack, might be better to implement a class separate from instruments for this purpose.
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('find_nv', False, bool, 'run FindNV'),
        Parameter('autofocus', False, bool, 'run AutoFocus')
        ])

    def __init__(self, name=None, settings=None):
        '''
        Connects to attocube controller ANC300 through USB using PySerial. Throws exception if unable to connect.
        baudrate was not specified in manual but 9600 seems to work.
        :param name: name of instruments
        :param settings: settings to update instrument with
        '''
        super().__init__(name, settings)
        self._is_connected = True


    def update(self, settings):
        '''
        Updates the internal settings, updating voltage or frequency
        Args:
            settings: a dictionary in the same form as settings with the new values
        '''
        super().update(settings)


    @property
    def _PROBES(self):
        return {}

    def read_probes(self, key):
        assert key in list(self._PROBES.keys())
        assert isinstance(key, str)


if __name__ == '__main__':
    commander = Commander()
    print(commander.settings)
    commander.settings.update({'find_nv': True})
    print(commander.settings)


