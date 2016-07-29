"""
    This file is part of b26_toolkit, a PyLabControl add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    Foobar is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""

from b26_toolkit.src.instruments import MaestroLightControl
from PyLabControl.src.core import Script, Parameter


class CameraOn(Script):
    # COMMENT_ME

    _DEFAULT_SETTINGS = [
        Parameter('On', True, bool, '')
    ]

    _INSTRUMENTS = {
        'light_control' : MaestroLightControl
    }
    _SCRIPTS = {}

    def __init__(self, instruments, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that makes use of an instrument
        Args:
            instruments: instruments the script will make use of
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """

        # call init of superclass
        Script.__init__(self, name, settings, instruments, log_function= log_function, data_path = data_path)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in _DEFAULT_SETTINGS
        for this dummy example we just implement a counter
        """

        if self.settings['On'] == True:
            # fluorescence filter
            self.instruments['white_light'].update({'open': False})
            self.instruments['filter_wheel'].update({'current_position': 'position_3'})
            self.instruments['block_ir'].update({'open': True})
            self.instruments['block_green'].update({'open': True})

            self.log('camera on')
        else:
            # high ND
            self.instruments['filter_wheel'].update({'current_position': 'position_1'})
            self.instruments['block_ir'].update({'open': False})
            self.instruments['block_green'].update({'open': False})
            self.instruments['white_light'].update({'open': True})



            self.log('camera off')