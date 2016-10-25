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
from b26_toolkit.src.instruments import NI9263


class MagnetCoils(NI9263):
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
                      Parameter('x_coil',
                                [
                                    Parameter('voltage_at_max_field', .368, float, 'input voltage corresponding to max field'),
                                    Parameter('x_field_max', 31.1, float, 'x field at max input voltage'),
                                    Parameter('y_field_max', 1.8, float, 'y field at max input voltage'),
                                    Parameter('z_field_max', -13.1, float, 'z field at max input voltage')
                                ]
                                ),
                      Parameter('y_coil',
                                [
                                    Parameter('voltage_at_max_field', .368, float, 'input voltage corresponding to max field'),
                                    Parameter('x_field_max', -6.9, float, 'x field at max input voltage'),
                                    Parameter('y_field_max', -32.6, float, 'y field at max input voltage'),
                                    Parameter('z_field_max', -13.4, float, 'z field at max input voltage')
                                ]
                                ),
                      Parameter('z_coil',
                                [
                                    Parameter('voltage_at_max_field', .260, float, 'input voltage corresponding to max field'),
                                    Parameter('x_field_max', .2, float, 'x field at max input voltage'),
                                    Parameter('y_field_max', -7.0, float, 'y field at max input voltage'),
                                    Parameter('z_field_max', 46.7, float, 'z field at max input voltage')
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
        def calc_voltages_for_fields(fields):
            """
            Calculates the voltage to apply to the current generating circuit that will result in the inputted field
            Args:
                field: magnetic field for which to find corresponding voltage
                axis: axis on which this voltage will be applied, must be one of ['x', 'y', 'z']

            Returns: voltage

            """
            max_voltages = np.array([self.settings['field_calibration']['x_coil']['voltage_at_max_field'], self.settings['field_calibration']['y_coil']['voltage_at_max_field'], self.settings['field_calibration']['z_coil']['voltage_at_max_field']])
            Cinv = self.calc_conversion_matrix()
            relative_voltages = Cinv * np.transpose(fields)
            print('relative', relative_voltages)
            new_voltages = np.dot(np.transpose(relative_voltages), max_voltages)
            print('new', new_voltages )

            if (new_voltages > max_voltages).any():
                raise ValueError('given field exceeds maximum possible field')

            return(new_voltages)

        # call the update_parameter_list to update the parameter list
        super(MagnetCoils, self).update(settings)
        # now we actually apply these newsettings to the hardware
        # if any of the settings updated are the fields...
        for key, value in settings.iteritems():
            if key == 'magnetic_fields':
                if any(x in value.keys() for x in ['x_field', 'y_field', 'z_field']):
                    new_field_x = self.settings['magnetic_fields']['x_field']
                    new_field_y = self.settings['magnetic_fields']['y_field']
                    new_field_z = self.settings['magnetic_fields']['z_field']
                    new_voltages = calc_voltages_for_fields(np.matrix([new_field_x, new_field_y, new_field_z]))
                    new_voltages = [new_voltages] * 2 #convert to form required for daq output
                    print('output voltages', new_voltages)
                    # self.AO_init([self.settings['magnet_channels']['x_channel'], self.settings['magnet_channels']['y_channel'],
                    #               self.settings['magnet_channels']['z_channel']], new_voltages)
                    # self.AO_run()
                    # self.AO_waitToFinish()
                    # self.AO_stop()

                    # even if multiple fields updated in the same pass, this will update all of them, so run this
                    # at most once
                    break


            # if key == 'magnetic_fields':
            #     if 'x_field' in value.keys():
            #         new_field = calc_voltage_for_field(value['x_field'], 'x')
            #         new_field =[new_field, new_field]
            #         self.AO_init([self.settings['magnet_channels']['x_channel']], new_field)
            #     elif 'y_field' in value.keys():
            #         new_field = calc_voltage_for_field(value['y_field'], 'y')
            #         new_field =[new_field, new_field]
            #         self.AO_init([self.settings['magnet_channels']['y_channel']], new_field)
            #     elif 'z_field' in value.keys():
            #         new_field = calc_voltage_for_field(value['z_field'], 'z')
            #         new_field =[new_field, new_field]
            #         self.AO_init([self.settings['magnet_channels']['z_channel']], new_field)
            #     else:
            #         raise('Invalid field axis, must be either x_field, y_field, or z_field')
            #
            #     self.AO_run()
            #     self.AO_waitToFinish()
            #     self.AO_stop()

    def calc_conversion_matrix(self):
        """
        One can map input voltages V to output B fields through the matrix equation B = C V, where C is a 3x3 conversion
        matrix and, for example, the first row of the matrix equation gives B_x = C_xx V_x + C_yx V_y + C_zx V_z. We
        want to go the other way, converting from B fields to voltages, so this equation inverts C.
        Returns:

        """
        cxx = self.settings['field_calibration']['x_coil']['x_field_max']
        cxy = self.settings['field_calibration']['x_coil']['y_field_max']
        cxz = self.settings['field_calibration']['x_coil']['z_field_max']
        cyx = self.settings['field_calibration']['y_coil']['x_field_max']
        cyy = self.settings['field_calibration']['y_coil']['y_field_max']
        cyz = self.settings['field_calibration']['y_coil']['z_field_max']
        czx = self.settings['field_calibration']['z_coil']['x_field_max']
        czy = self.settings['field_calibration']['z_coil']['y_field_max']
        czz = self.settings['field_calibration']['z_coil']['z_field_max']

        C = np.matrix([[cxx, cyx, czx], [cxy, cyy, czy], [cxz, cyz, czz]])
        print('C', C)

        print('Cinv', np.linalg.inv(C))
        return(np.linalg.inv(C))



if __name__ == '__main__':
    instruments, failed = Instrument.load_and_append(instrument_dict={'MagnetCoils': MagnetCoils})

    instruments['MagnetCoils'].update({'magnetic_fields': {'z_field': 1}})

    # print(instruments['MagnetCoils'].is_connected())
    # instruments['MagnetCoils'].update({'magnetic_fields':{'x_field': 1}})