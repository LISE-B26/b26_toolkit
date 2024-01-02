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


from b26_toolkit.instruments import ANC300
from pylabcontrol.core import Parameter, Script


class AttoOffsetZ(Script):
    """
    Sets the offset of Z-axis Attocube.
    Deliberately left out the other axes to avoid user error
    """
    _DEFAULT_SETTINGS = [
        Parameter('z_offset', 0, float, 'Offset voltage for Z-axis Attocube')
    ]

    _INSTRUMENTS = {'ANC300': ANC300}
    _SCRIPTS = {}

    def _function(self):
        """
        Performs a multiple attocube step with the voltage and frequency specified in instrument,
        and the direction and number specified in settings
        """
        attocube = self.instruments['ANC300']['instance']
        attocube._set_offset(3, self.settings['z_offset'])
        self.log('Z-axis Attocube offset changed to %.1f V' % attocube.read_probes('z_offset'))
