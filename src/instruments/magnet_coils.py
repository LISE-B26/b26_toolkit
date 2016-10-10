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

from PyLabControl.src.core import Instrument, Parameter
from b26_toolkit.src.instruments import DAQ


class MagnetCoils(DAQ):
    """

    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('device', 'Dev1', (str), 'Name of DAQ device'),
        Parameter('analog_output',
                  [
                      Parameter('ao0',
                        [
                            Parameter('channel', 0, [0, 1, 2, 3], 'output channel'),
                            Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                            Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                            Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                        ]
                                ),
                      Parameter('ao1',
                        [
                            Parameter('channel', 1, [0, 1, 2, 3], 'output channel'),
                            Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                            Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                            Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                        ]
                                ),
                      Parameter('ao2',
                        [
                            Parameter('channel', 2, [0, 1, 2, 3], 'output channel'),
                            Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                            Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                            Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                        ]
                                ),
                      Parameter('ao3',
                        [
                            Parameter('channel', 3, [0, 1, 2, 3], 'output channel'),
                            Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                            Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                            Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                        ]
                                )
                  ]
                  ),
        Parameter('magnet_channels',
                  [
                      Parameter('x_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'output channel for x field'),
                      Parameter('y_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'output channel for x field'),
                      Parameter('z_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'], 'output channel for x field')
                  ]),
        Parameter('magnetic_fields',
                  [
                      Parameter('x_field', 0, float, 'x field in Gauss'),
                      Parameter('y_field', 0, float, 'y field in Gauss'),
                      Parameter('z_field', 0, float, 'z field in Gauss')
                  ]),
        Parameter('field_calibration',
                  [
                      Parameter('x',
                                [
                                    Parameter('max_field', 150.0, float, 'maximum allowed magnetic field on each axis'),
                                    Parameter('voltage_at_max_field', 1.0, float, 'input voltage corresponding to max field')
                                ]
                                ),
                      Parameter('y',
                                [
                                    Parameter('max_field', 150.0, float, 'maximum allowed magnetic field on each axis'),
                                    Parameter('voltage_at_max_field', 1.0, float,
                                              'input voltage corresponding to max field')
                                ]
                                ),
                      Parameter('z',
                                [
                                    Parameter('max_field', 150.0, float, 'maximum allowed magnetic field on each axis'),
                                    Parameter('voltage_at_max_field', 1.0, float,
                                              'input voltage corresponding to max field')
                                ]
                                )
                  ])

    ])

    def update(self, settings):
        """
        Updates the settings in software and, if applicable, takes an action to modify the hardware, such as opening
        a beamblock or spinning a filterwheel
        Args:
            settings: a dictionary in the standard settings format
        """
        def calc_voltage_for_field(field, axis):
            """
            Calculates the voltage to apply to the current generating circuit that will result in the inputted field
            Args:
                field: magnetic field for which to find corresponding voltage
                axis: axis on which this voltage will be applied, must be one of ['x', 'y', 'z']

            Returns: voltage

            """
            if axis not in ['x', 'y', 'z']:
                raise ValueError('invalid axis given')
            max_field = self.settings['field_calibration'][axis]['max_field']
            if field > max_field:
                raise ValueError('given field exceeds maximum possible field')
            voltage_at_max_field = self.settings['field_calibration'][axis]['voltage_at_max_field']
            return((field / max_field * voltage_at_max_field))

        # call the update_parameter_list to update the parameter list
        super(MagnetCoils, self).update(settings)
        # now we actually apply these newsettings to the hardware
        for key, value in settings.iteritems():
            if key == 'magnetic_fields':
                if 'x_field' in value.keys():
                    new_field = calc_voltage_for_field(value['x_field'], 'x')
                    new_field =[new_field, new_field]
                    self.AO_init([self.settings['magnet_channels']['x_channel']], new_field)
                elif 'y_field' in value.keys():
                    new_field = calc_voltage_for_field(value['y_field'], 'y')
                    new_field =[new_field, new_field]
                    self.AO_init([self.settings['magnet_channels']['y_channel']], new_field)
                elif 'z_field' in value.keys():
                    new_field = calc_voltage_for_field(value['z_field'], 'z')
                    new_field =[new_field, new_field]
                    self.AO_init([self.settings['magnet_channels']['z_channel']], new_field)
                else:
                    raise('Invalid field axis, must be either x_field, y_field, or z_field')

                self.AO_run()
                self.AO_waitToFinish()
                self.AO_stop()


if __name__ == '__main__':
    instruments, failed = Instrument.load_and_append(instrument_dict={'MagnetCoils': MagnetCoils})

    print(instruments['MagnetCoils'].is_connected())
    # instruments['MagnetCoils'].update({'magnetic_fields':{'x_field': 1}})