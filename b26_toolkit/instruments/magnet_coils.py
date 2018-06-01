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

from pylabcontrol.core import Instrument, Parameter
from b26_toolkit.instruments import NI9263


class MagnetCoils(NI9263):
    """

    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('device', 'cDAQ9184-1BA7633Mod4', (str), 'Name of DAQ device'),
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
                      Parameter('x_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'output channel for x field'),
                      Parameter('y_channel', 'ao2', ['ao0', 'ao1', 'ao2', 'ao3'], 'output channel for x field'),
                      Parameter('z_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'output channel for x field')
                  ]),
        Parameter('magnetic_fields',
                  [
                      Parameter('x_field', 0, float, 'x field in Gauss'),
                      Parameter('y_field', 0, float, 'y field in Gauss'),
                      Parameter('z_field', 0, float, 'z field in Gauss')
                  ]),
        # Parameter('field_calibration',
        #           [
        #               Parameter('x_coil',
        #                         [
        #                             Parameter('voltage_at_max_field', .368, float, 'input voltage corresponding to max field'),
        #                             Parameter('x_field_max', 32.8, float, 'x field at max input voltage'),
        #                             Parameter('y_field_max', 0.4, float, 'y field at max input voltage'),
        #                             Parameter('z_field_max', -11.4, float, 'z field at max input voltage')
        #                         ]
        #                         ),
        #               Parameter('y_coil',
        #                         [
        #                             Parameter('voltage_at_max_field', .278, float, 'input voltage corresponding to max field'),
        #                             Parameter('x_field_max', -5.7, float, 'x field at max input voltage'),
        #                             Parameter('y_field_max', -24.0, float, 'y field at max input voltage'),
        #                             Parameter('z_field_max', -7.9, float, 'z field at max input voltage')
        #                         ]
        #                         ),
        #               Parameter('z_coil',
        #                         [
        #                             Parameter('voltage_at_max_field', .245, float, 'input voltage corresponding to max field'),
        #                             Parameter('x_field_max', 1.9, float, 'x field at max input voltage'),
        #                             Parameter('y_field_max', -5.6, float, 'y field at max input voltage'),
        #                             Parameter('z_field_max', 47.7, float, 'z field at max input voltage')
        #                         ]
        #                         )
        #           ]),
        # Parameter('field_calibration',
        #           [
        #               Parameter('x_coil',
        #                         [
        #                             Parameter('voltage_at_max_field', .736, float,
        #                                       'input voltage corresponding to max field'),
        #                             Parameter('x_field_max', 65.6, float, 'x field in Gauss at max input voltage'),
        #                             Parameter('y_field_max', .8, float, 'y field at max input voltage'),
        #                             Parameter('z_field_max', -22.8, float, 'z field at max input voltage'),
        #                             Parameter('flip_polarity', False, bool, 'if true, signs of the field as defined above above are flipped')
        #                         ]
        #                         ),
        #               Parameter('y_coil',
        #                         [
        #                             Parameter('voltage_at_max_field', .556, float,
        #                                       'input voltage corresponding to max field'),
        #                             Parameter('x_field_max', -11.4, float, 'x field at max input voltage'),
        #                             Parameter('y_field_max', -48.0, float, 'y field at max input voltage'),
        #                             Parameter('z_field_max', -15.8, float, 'z field at max input voltage'),
        #                             Parameter('flip_polarity', False, bool, 'if true, signs of the field as defined above above are flipped')
        #                         ]
        #                         ),
        #               Parameter('z_coil',
        #                         [
        #                             Parameter('voltage_at_max_field', .49, float,
        #                                       'input voltage corresponding to max field'),
        #                             Parameter('x_field_max', 3.8, float, 'x field at max input voltage'),
        #                             Parameter('y_field_max', -11.2, float, 'y field at max input voltage'),
        #                             Parameter('z_field_max', 95.4, float, 'z field at max input voltage'),
        #                             Parameter('flip_polarity', False, bool, 'if true, signs of the field as defined above above are flipped')
        #                         ]
        #                         )
        #           ]),
        Parameter('field_calibration',
                  [
                      Parameter('x_coil',
                                [
                                    Parameter('voltage_at_max_field', .736, float,
                                              'input voltage corresponding to max field'),
                                    Parameter('x_field_max', 0, float, 'x field at max input voltage'),
                                    Parameter('y_field_max', 0, float, 'y field at max input voltage'),
                                    Parameter('z_field_max', 0, float, 'z field at max input voltage'),
                                    Parameter('flip_polarity', False, bool,
                                              'if true, signs of the field as defined above above are flipped')
                                ]
                                ),
                      Parameter('y_coil',
                                [
                                    Parameter('voltage_at_max_field', .556, float,
                                              'input voltage corresponding to max field'),
                                    Parameter('x_field_max', 0, float, 'x field at max input voltage'),
                                    Parameter('y_field_max', 0, float, 'y field at max input voltage'),
                                    Parameter('z_field_max', 0, float, 'z field at max input voltage'),
                                    Parameter('flip_polarity', False, bool,
                                              'if true, signs of the field as defined above above are flipped')
                                ]
                                ),
                      Parameter('z_coil',
                                [
                                    Parameter('voltage_at_max_field', .408, float,
                                              'input voltage corresponding to max field'),
                                    Parameter('x_field_max', 0, float, 'x field at max input voltage'),
                                    Parameter('y_field_max', 0, float, 'y field at max input voltage'),
                                    Parameter('z_field_max', 152.88, float, 'z field at max input voltage'),
                                    Parameter('flip_polarity', False, bool,
                                              'if true, signs of the field as defined above above are flipped')
                                ]
                                )
                  ]),
        Parameter('use_approximate_fields', True, bool, 'if out of bounds fields given, uses closest allowed fields'),
        Parameter('magnetic fields input type', 'Gauss', ['Gauss', 'Volts'])
    ])

    def update(self, settings):
        """
        Updates the settings in software and, if applicable, takes an action to modify the hardware, such as opening
        a beamblock or spinning a filterwheel
        Args:
            settings: a dictionary in the standard settings format
        """
        # call the update_parameter_list to update the parameter list
        super(MagnetCoils, self).update(settings)
        # now we actually apply these newsettings to the hardware
        # if any of the settings updated are the fields...
        for key, value in settings.items():
            if key == 'magnetic_fields':
                if any(x in list(value.keys()) for x in ['x_field', 'y_field', 'z_field']):
                    new_field_x = self.settings['magnetic_fields']['x_field']
                    new_field_y = self.settings['magnetic_fields']['y_field']
                    new_field_z = self.settings['magnetic_fields']['z_field']
                    new_voltages = self.calc_voltages_for_fields(np.array([new_field_x, new_field_y, new_field_z]))
                    # convert to form required for daq output
                    # new_voltages = np.transpose(np.column_stack((new_voltages[0], new_voltages[1], new_voltages[2])))
                    # new_voltages = (np.repeat(new_voltages, 2, axis=1))
                    # print('nv', new_voltages)
                    self.set_analog_voltages({self.settings['magnet_channels']['x_channel']: new_voltages[0], self.settings['magnet_channels']['y_channel']: new_voltages[1],
                                   self.settings['magnet_channels']['z_channel']: new_voltages[2]})

                    # even if multiple fields updated in the same pass, this will update all of them, so run this
                    # at most once
                    break

    def calc_conversion_matrix(self):
        """
        One can map input voltages V to output B fields through the matrix equation B = C V, where C is a 3x3 conversion
        matrix and, for example, the first row of the matrix equation gives B_x = C_xx V_x + C_yx V_y + C_zx V_z. We
        want to go the other way, converting from B fields to voltages, so this equation inverts C.
        Returns: conversion matrix to map from V to B, and its inverse to map from B to V

        """
        polarity = {'x': 1, 'y': 1, 'z': 1}
        if self.settings['field_calibration']['x_coil']['flip_polarity']:
            polarity['x'] = -1
        if self.settings['field_calibration']['y_coil']['flip_polarity']:
            polarity['y'] = -1
        if self.settings['field_calibration']['z_coil']['flip_polarity']:
            polarity['z'] = -1
        cxx = self.settings['field_calibration']['x_coil']['x_field_max'] * polarity['x']
        cxy = self.settings['field_calibration']['x_coil']['y_field_max'] * polarity['x']
        cxz = self.settings['field_calibration']['x_coil']['z_field_max'] * polarity['x']
        cyx = self.settings['field_calibration']['y_coil']['x_field_max'] * polarity['y']
        cyy = self.settings['field_calibration']['y_coil']['y_field_max'] * polarity['y']
        cyz = self.settings['field_calibration']['y_coil']['z_field_max'] * polarity['y']
        czx = self.settings['field_calibration']['z_coil']['x_field_max'] * polarity['z']
        czy = self.settings['field_calibration']['z_coil']['y_field_max'] * polarity['z']
        czz = self.settings['field_calibration']['z_coil']['z_field_max'] * polarity['z']

        Cmat = np.matrix([[cxx, cyx, czx], [cxy, cyy, czy], [cxz, cyz, czz]])
        return(Cmat, np.linalg.inv(Cmat))

    def calc_voltages_for_fields(self, fields):
        """
        Calculates the voltage to apply to the current generating circuit that will result in the inputted field
        Args:
            field: magnetic field for which to find corresponding voltage
            axis: axis on which this voltage will be applied, must be one of ['x', 'y', 'z']

        Returns: voltage

        """
        used_approx_fields = False
        max_voltages = np.array([self.settings['field_calibration']['x_coil']['voltage_at_max_field'],
                                 self.settings['field_calibration']['y_coil']['voltage_at_max_field'],
                                 self.settings['field_calibration']['z_coil']['voltage_at_max_field']])
        #first solve exactly for voltages to apply for the requested fields
        Cmat, Cinv = self.calc_conversion_matrix()
        relative_voltages = np.matmul(Cinv, fields)
        new_voltages = np.array([relative_voltages[0, 0] * max_voltages[0], relative_voltages[0, 1] * max_voltages[1],
                                 relative_voltages[0, 2] * max_voltages[2]])

        # self.log('Calculated voltages to reach desired field: ' + str(new_voltages))

        if self.settings['use_approximate_fields']:
            if (np.abs(new_voltages) > max_voltages).any():
                relative_voltages = self.optimize_fields(fields, Cmat, relative_voltages)
                new_voltages = np.array(
                    [relative_voltages[0] * max_voltages[0], relative_voltages[1] * max_voltages[1],
                     relative_voltages[2] * max_voltages[2]])
                self.log('Given field exceeds maximum possible field, outputting best approximation with these voltages: ' + str(new_voltages))
                self.log('Corresponding best possible field (in G): ' + str(np.dot(Cmat, relative_voltages)))
                used_approx_fields = True
            elif (new_voltages < 0).any():
                relative_voltages = self.optimize_fields(fields, Cmat, relative_voltages)
                new_voltages = np.array(
                    [relative_voltages[0] * max_voltages[0], relative_voltages[1] * max_voltages[1],
                     relative_voltages[2] * max_voltages[2]])
                print('Could not match desired fields, outputting best approximation')
                self.log('Given field exceeds maximum possible field, outputting best approximation with these voltages: ' + str(new_voltages))
                self.log('Corresponding best possible field (in G): ' + str(np.dot(Cmat, relative_voltages)))
                used_approx_fields = True
        else:
            if (np.abs(new_voltages) > max_voltages).any():
                self.log('The requested field exceeds the maximum possible field')
                raise ValueError('Given field exceeds maximum possible field')
            elif (new_voltages < 0).any():
                mask = (new_voltages < 0)
                negative_axes = [x[1] for x in zip(*(mask, ['x', 'y', 'z'])) if x[0]]
                self.log('Field not reachable, polarity switch required on ' + ','.join('{}'.format(k) for k in negative_axes))
                raise ValueError('Field not reachable, polarity switch required on ' + ','.join('{}'.format(k) for k in negative_axes))

        #to allow scripts to access data
        self.new_voltages = new_voltages
        self.requested_fields = fields
        if used_approx_fields:
            self.applied_fields = np.transpose(np.dot(Cmat, relative_voltages))
        else:
            self.applied_fields = fields
        return (new_voltages)

    def optimize_fields(self, fields, Cmat, proposed_voltage):
        """
        Finds closest approximation to fields based on the reachable space allowed by Cmat. In particular,
        finds the output voltage V that minimizes |Cmat * V - fields|^2, or equivalently
        |output_fields - requested_fields|^2
        Args:
            fields: field to match
            Cmat: mapping from voltages to output fields
            proposed_voltage: initial guess for fields. Cinv * fields is a good place to start.

        Returns:

        """
        import scipy.optimize as optimize

        def f(x):
            y = np.squeeze(np.asarray(np.dot(Cmat, x) - fields)) #convert from matrix to array
            return np.dot(y, y)

        bounds = ((0, 1), (0, 1), (0, 1))
        opt = optimize.minimize(f, proposed_voltage, bounds = bounds)

        return opt.x


if __name__ == '__main__':
    instruments, failed = Instrument.load_and_append(instrument_dict={'MagnetCoils': MagnetCoils})

    # instruments['MagnetCoils'].update({'magnetic_fields': {'x_field': 10, 'y_field': -10, 'z_field': 40}})
    instruments['MagnetCoils'].update({'magnetic_fields': {'x_field': 10, 'y_field': -1, 'z_field': 0}})


    # print(instruments['MagnetCoils'].is_connected())
    # instruments['MagnetCoils'].update({'magnetic_fields':{'x_field': 1}})