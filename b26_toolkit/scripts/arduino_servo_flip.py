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

from pylabcontrol.core import Script, Parameter
from b26_toolkit.instruments import ArduinoZero, B26PulseBlaster


class ArduinoServoFlip(Script):
    _DEFAULT_SETTINGS = [Parameter('servo_num', 0, int, 'servo id'),
                         Parameter('position', True, bool, 'checked for up, unchecked for down')]
    
    _INSTRUMENTS = {'arduino': ArduinoZero}
    _SCRIPTS = {}

    def _function(self):
        self.instruments['arduino']['instance'].flip(self.settings['servo_num'], self.settings['position'])

class ToggleCameraView(Script):
    _INSTRUMENTS = {'arduino': ArduinoZero, 'PB': B26PulseBlaster}
    _DEFAULT_SETTINGS = [Parameter('camera_view', False, bool, 'Toggle camera view')]

    def _function(self):
        state = self.settings['camera_view']
        #self.instruments['arduino']['instance'].update({'LED_pellicle': {'status': state}})
        #self.instruments['arduino']['instance'].update({'camera_pellicle': {'status': state}})

        self.instruments['arduino']['instance']['LED_pellicle']['status'] = state
        self.instruments['arduino']['instance']['camera_pellicle']['status'] = state

        self.instruments['PB']['instance'].update({'laser': {'status': not state}})

