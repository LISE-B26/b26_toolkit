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

import numpy as np
from matplotlib import patches

from b26_toolkit.instruments import MagnetCoils
from pylabcontrol.core import Script, Parameter
from b26_toolkit.data_processing.coordinate_conversions import spherical_to_cartesian


class SetMagneticCoils(Script):
    """
This script sets the magnetic field coils to the given magnetic field values
    """

    _DEFAULT_SETTINGS = [
        Parameter('magnetic_fields',
                  [
                      Parameter('x/r_field', 0.0, float, 'x or r component of field vector in Gauss'),
                      Parameter('y/theta_field', 0.0, float, 'y or theta component of field vector in Gauss'),
                      Parameter('z/phi_field', 0.0, float, 'z or phi component of field vector in Gauss'),
                      Parameter('coordinate_system', 'Cartesian', ['Cartesian', 'Spherical'], 'Chooses what coordinate system to use to interpret field parameters')
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

        self.data = {}

        input_vector = [self.settings['magnetic_fields']['x/r_field'], self.settings['magnetic_fields']['y/theta_field'], self.settings['magnetic_fields']['z/phi_field']]

        if self.settings['magnetic_fields']['coordinate_system'] == 'Spherical':
            [new_x_field, new_y_field, new_z_field] = spherical_to_cartesian(input_vector)
        elif self.settings['magnetic_fields']['coordinate_system'] == 'Cartesian':
            [new_x_field, new_y_field, new_z_field] = input_vector

        # new_x_field = 0
        # new_y_field = 0
        # new_z_field = 10

        # try:
        #     self.instruments['MagnetCoils']['instance'].calc_voltages_for_fields([new_x_field, new_y_field, new_z_field])
        # except ValueError:
        #     raise
        #     self.log('Could not set magnetic field. Reverting to previous value.')

        #need to cast from float64 in the 'Spherical' case

        # take settings defined in the script and update with settings for the fields
        dictator = self.instruments['MagnetCoils']['settings']
        dictator.update({'magnetic_fields': {'x_field': float(new_x_field), 'y_field': float(new_y_field), 'z_field': float(new_z_field)}})

        self.instruments['MagnetCoils']['instance'].update(dictator)
        #
        self.log('Magnetic Field set to Bx={:.4}G, By={:.4}G, Bz={:.4}G'.format(self.instruments['MagnetCoils']['instance'].settings['magnetic_fields']['x_field'],
                                                                                self.instruments['MagnetCoils']['instance'].settings['magnetic_fields']['y_field'],
                                                                                self.instruments['MagnetCoils']['instance'].settings['magnetic_fields']['z_field'])
                 )

        print(('requested fields', self.instruments['MagnetCoils']['instance'].requested_fields))
        print(('applied fields', self.instruments['MagnetCoils']['instance'].applied_fields))

        self.data['new_voltages'] = self.instruments['MagnetCoils']['instance'].new_voltages
        self.data['requested_fields'] = self.instruments['MagnetCoils']['instance'].requested_fields
        self.data['applied_fields'] = self.instruments['MagnetCoils']['instance'].applied_fields

if __name__ == '__main__':
    from pylabcontrol.core import Instrument

    instruments, instruments_failed = Instrument.load_and_append({'MagnetCoils': MagnetCoils })

    script, failed, instruments = Script.load_and_append(script_dict={'SetMagneticCoils': SetMagneticCoils}, instruments = instruments)

    script['SetMagneticCoils']._function()

    print(script)
    print(failed)
    # print(instruments)