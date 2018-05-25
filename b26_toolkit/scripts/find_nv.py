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

from copy import deepcopy

import numpy as np
import trackpy as tp
from matplotlib import patches

from b26_toolkit.plotting.plots_2d import plot_fluorescence_new
from pylabcontrol.core import Script, Parameter
from b26_toolkit.scripts import GalvoScan, SetLaser


class FindNV(Script):
    """
GalvoScan uses the apd, daq, and galvo to sweep across voltages while counting photons at each voltage,
resulting in an image in the current field of view of the objective.

Known issues:
    1.) if fits are poor, check  sweep_range. It should extend significantly beyond end of NV on both sides.
    """

    _DEFAULT_SETTINGS = [
        Parameter('initial_point',
                  [Parameter('x', 0, float, 'x-coordinate'),
                   Parameter('y', 0, float, 'y-coordinate')
                   ]),
        Parameter('sweep_range', .03, float, 'voltage range to sweep over to find a max'),
        Parameter('num_points', 61, int, 'number of points to sweep in the sweep range'),
        Parameter('nv_size', 11, int, 'TEMP: size of nv in pixels - need to be refined!!'),
        Parameter('min_mass', 180, int, 'TEMP: brightness of nv - need to be refined!!'),
        Parameter('number_of_attempts', 1, int, 'Number of times to decrease min_mass if an NV is not found'),
        Parameter('center_on_current_location', False, bool, 'Check to use current galvo location rather than ')
    ]

    _INSTRUMENTS = {}

    _SCRIPTS = {'take_image': GalvoScan, 'set_laser': SetLaser}

    def __init__(self, scripts, name = None, settings = None, log_function = None, timeout = 1000000000, data_path = None):

        Script.__init__(self, name, scripts = scripts, settings=settings, log_function=log_function, data_path = data_path)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        attempt_num = 1
        print((self.scripts['take_image'].instruments))
        if self.settings['center_on_current_location']:
            # fixed for cold setup and new DAQ ER 6/4/17
            #daq_pt = self.scripts['take_image'].instruments['daq']['instance'].get_analog_voltages([self.scripts['take_image'].settings['DAQ_channels']['x_ao_channel'], self.scripts['take_image'].settings['DAQ_channels']['y_ao_channel']])
            daq_pt = self.scripts['take_image'].instruments['NI6259']['instance'].get_analog_voltages([self.scripts['take_image'].settings['DAQ_channels']['x_ao_channel'], self.scripts['take_image'].settings['DAQ_channels']['y_ao_channel']])
            self.settings['initial_point'].update({'x': daq_pt[0], 'y': daq_pt[1]})
        initial_point = self.settings['initial_point']
        nv_size = self.settings['nv_size']
        min_mass = self.settings['min_mass']

        self.data = {'maximum_point': None,
                     'initial_point': initial_point,
                     'image_data': [],
                     'extent': [],
                     'fluorescence': None
                     }

        def pixel_to_voltage(pt, extent, image_dimensions):
            """"
            pt: point in pixels
            extent: [xVmin, Vmax, Vmax, yVmin] in volts
            image_dimensions: dimensions of image in pixels

            Returns: point in volts
            """

            image_x_len, image_y_len = image_dimensions
            image_x_min, image_x_max, image_y_max, image_y_min = extent

            assert image_x_max > image_x_min
            assert image_y_max > image_y_min

            volt_per_px_x = (image_x_max - image_x_min) / (image_x_len-1)
            volt_per_px_y = (image_y_max - image_y_min) / (image_y_len-1)

            V_x = volt_per_px_x * pt[0] + image_x_min
            V_y = volt_per_px_y * pt[1] + image_y_min

            return [V_x, V_y]

        def min_mass_adjustment(min_mass):
            #COMMENT_ME
            return (min_mass - 40)

        self.scripts['take_image'].settings['point_a'].update({'x': self.settings['initial_point']['x'], 'y': self.settings['initial_point']['y']})
        self.scripts['take_image'].settings['point_b'].update({'x': self.settings['sweep_range'], 'y': self.settings['sweep_range']})
        self.scripts['take_image'].update({'RoI_mode': 'center'})
        self.scripts['take_image'].settings['num_points'].update({'x': self.settings['num_points'], 'y': self.settings['num_points']})

        self.scripts['take_image'].run()

        self.data['image_data'] = deepcopy(self.scripts['take_image'].data['image_data'])
        self.data['extent'] = deepcopy(self.scripts['take_image'].data['extent'])
        while True: # modified ER 5/27/2017 to implement tracking
            locate_info = tp.locate(self.data['image_data'], nv_size, minmass=min_mass)
            po = [self.data['initial_point']['x'], self.data['initial_point']['y']]
            if len(locate_info) == 0:
                self.data['maximum_point'] = {'x': float(po[0]), 'y': float(po[1])}
            else:

                # all the points that have been identified as valid NV centers
                pts = [pixel_to_voltage(p, self.data['extent'], np.shape(self.data['image_data'])) for p in
                       locate_info[['x', 'y']].as_matrix()]
                #print(pts)
                if len(pts) > 1:
                    self.log('FindNV found more than one NV in the scan image. Selecting the one closest to initial point.')
                # pick the one that is closest to the original one
                pm = pts[np.argmin(np.array([np.linalg.norm(p - np.array(po)) for p in pts]))]
                self.data['maximum_point'] = {'x': float(pm[0]), 'y': float(pm[1])}
                counter = 0
                for p in pts: # record maximum counts = fluorescence
                    if p[1] == self.data['maximum_point']['y']:
                        self.data['fluorescence'] = 2*locate_info[['signal']].as_matrix()[counter]
                        print('fluorescence of the NV, kCps:')
                        print((self.data['fluorescence']))
                        counter += 1
                break

            if attempt_num <= self.settings['number_of_attempts']:
                min_mass = min_mass_adjustment(min_mass)
                attempt_num += 1
            else:
                self.log('FindNV did not find an NV --- setting laser to initial point instead, setting fluorescence to zero')
                self.data['fluorescence'] = 0.0
                break

        self.scripts['set_laser'].settings['point'].update(self.data['maximum_point'])
        self.scripts['set_laser'].run()


    @staticmethod
    def plot_data(axes_list, data):
        plot_fluorescence_new(data['image_data'], data['extent'], axes_list[0])

        initial_point = data['initial_point']
        patch = patches.Circle((initial_point['x'], initial_point['y']), .001, ec='g', fc='none', ls='dashed')
        axes_list[0].add_patch(patch)
        axes_list[0].text(initial_point['x'], initial_point['y'] - .002, 'initial point', color='g', fontsize=8)

        # plot marker
        if data['maximum_point']:
            maximum_point = data['maximum_point']
            patch = patches.Circle((maximum_point['x'], maximum_point['y']), .001, ec='r', fc='none', ls='dashed')
            axes_list[0].add_patch(patch)
            axes_list[0].text(maximum_point['x'], maximum_point['y'] - .002, 'found NV', color='r', fontsize=8)

    def _plot(self, axes_list, data = None):
        """
        plotting function for find_nv
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
            data: data (dictionary that contains keys image_data, extent, initial_point, maximum_point) if not provided use self.data
        """
        if data is None:
            data = self.data


        if self._current_subscript_stage['current_subscript'] == self.scripts['take_image']:
            self.scripts['take_image']._plot(axes_list)
        else:
            self.plot_data(axes_list, data)


    def _update_plot(self, axes_list):
        """
        update plotting function for find_nv
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
        """

        if self._current_subscript_stage['current_subscript'] == self.scripts['take_image']:
            self.scripts['take_image']._update_plot(axes_list)

        if self.data['maximum_point']:
            maximum_point = self.data['maximum_point']
            patch = patches.Circle((maximum_point['x'], maximum_point['y']), .001, ec='r', fc='none', ls='dashed')
            axes_list[0].add_patch(patch)
            axes_list[0].text(maximum_point['x'], maximum_point['y'] - .002, 'found NV', color='r', fontsize=8)



    def get_axes_layout(self, figure_list):
        """
        returns the axes objects the script needs to plot its data
        the default creates a single axes object on each figure
        This can/should be overwritten in a child script if more axes objects are needed
        Args:
            figure_list: a list of figure objects
        Returns:
            axes_list: a list of axes objects

        """

        # create a new figure list that contains only figure 1, this assures that the super.get_axes_layout doesn't
        # empty the plot contained on figure 2
        return super(FindNV, self).get_axes_layout([figure_list[0]])

# class FindNV_cDAQ(FindNV):
#     """
# GalvoScan uses the apd, daq, and galvo to sweep across voltages while counting photons at each voltage,
# resulting in an image in the current field of view of the objective.
#
# Known issues:
#     1.) if fits are poor, check  sweep_range. It should extend significantly beyond end of NV on both sides.
#     """
#
#     _SCRIPTS = {'take_image': GalvoScan_cDAQ, 'set_laser': SetLaser_cDAQ}
#
#
#
#     if __name__ == '__main__':
#         script, failed, instruments = Script.load_and_append(script_dict={'FindMaxCounts': 'FindMaxCounts'})
#
#         print(script)
#         print(failed)
#         print(instruments)
