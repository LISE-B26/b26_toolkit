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

import numpy as np
from matplotlib import patches

from b26_toolkit.src.instruments import MagnetCoils
from PyLabControl.src.core import Script, Parameter


class SetMagneticCoils(Script):
    """
This script sets the magnetic field coils to the given magnetic field values
    """

    _DEFAULT_SETTINGS = [
        Parameter('magnetic_fields',
                  [
                      Parameter('x_field', 0, float, 'x field in Gauss'),
                      Parameter('y_field', 0, float, 'y field in Gauss'),
                      Parameter('z_field', 0, float, 'z field in Gauss')
                  ])
    ]

    _INSTRUMENTS = {'MagnetCoils':  MagnetCoils}

    _SCRIPTS = {}


    def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_function= log_function, data_path = data_path)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """
        new_x_field = self.settings['magnetic_fields']['x_field']
        new_y_field = self.settings['magnetic_fields']['y_field']
        new_z_field = self.settings['magnetic_fields']['z_field']

        self.instruments['MagnetCoils']['instance'].update({'magnetic_fields': {'x_field': new_x_field, 'y_field': new_y_field, 'z_field': new_z_field}})

        self.log('Magnetic Field set to Bx={:.4}G, By={:.4}G, Bz={:.4}G'.format(self.settings['point']['x'], self.settings['point']['y']))


if __name__ == '__main__':
    from PyLabControl.src.core import Instrument

    instruments, instruments_failed = Instrument.load_and_append({'daq':  'DAQ'})

    script, failed, instruments = Script.load_and_append(script_dict={'SetLaser':'SetLaser'}, instruments = instruments)

    print(script)
    print(failed)
    # print(instruments)