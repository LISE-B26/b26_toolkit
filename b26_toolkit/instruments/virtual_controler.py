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

	#########
#
#p=PID(3.0,0.4,1.2)
#p.setPoint(5.0)
#while True:
#     pid = p.update(measurement_value)
#
#
from pylabcontrol.core import Instrument,Parameter



class Plant(Instrument):
    _DEFAULT_SETTINGS = Parameter([
        Parameter('set_point', 0.0, float, 'setpoint to which to stabilize'),
        Parameter('gains', [
            Parameter('proportional', 0.0, float, 'proportional gain'),
            Parameter('integral', 0.0, float, 'integral gain')
        ]),
        Parameter('time_step', 1.0, float, 'time_step of loop'),
        Parameter('output_range', [
            Parameter('min', -10000, float, 'min allowed value for PI-loop output'),
            Parameter('max', 10000, float, 'max allowed value for PI-loop output')
        ]),
    ])
    _PROBES = {}

    def __init__(self, name=None, settings=None):
        super(PIControler, self).__init__(name, settings)
        self.reset()

    def update(self, settings):
        super(PIControler, self).update(settings)

if __name__ == '__main__':
    pi = PIControler()
