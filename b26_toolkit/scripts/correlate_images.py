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

from b26_toolkit.data_processing import correlation, shift_NVs
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from pylabcontrol.core import Script, Parameter
from b26_toolkit.scripts import GalvoScan
from copy import deepcopy
import numpy as np
from matplotlib import patches


class Take_And_Correlate_Images(Script):
    '''
    Takes a galvo scan, compares it to a previous galvo scan to find the relative shift, and then updates a list of
    nvs based on this shift so that they will give the current coordinates of those nvs
    '''

    _DEFAULT_SETTINGS = [
        Parameter('use_trackpy', False, bool, 'Use trackpy to create artificial nv-only images to filter out background')
    ]

    _INSTRUMENTS = {}
    _SCRIPTS = {'GalvoScan': GalvoScan}
    # _SCRIPTS = {}

    def __init__(self, instruments = None, name = None, settings = None, scripts = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_function= log_function, data_path = data_path)

        self.data = {'baseline_image': [], 'new_image': [], 'image_extent': [], 'correlation_image': []}

    def _function(self):
        """
        # Takes a new image, and correlates this with the image provided to baseline_image in self.data. It uses the
        determined pixel shift to calculate a shift for each of the nvs in the old_nv_list, which is given to it by
        a superscript, then store it as new_NV_list in data
        """

        if self.data['baseline_image'] == []:
            self.log('No baseline image avaliable. Taking baseline.')
        elif self.data['image_extent'] == []:
            self.log('No image extent avaliable. Script may have been run in error.')

        if not self.data['baseline_image'] == []:
            # use same settings as initial image
            scan = self.scripts['GalvoScan']
            scan.settings['point_a']['x'] = self.data['image_extent'][0]
            scan.settings['point_b']['x'] = self.data['image_extent'][1]
            scan.settings['point_a']['y'] = self.data['image_extent'][3]
            scan.settings['point_b']['y'] = self.data['image_extent'][2]
            scan.settings['RoI_mode'] = 'corner'

            self.scripts['GalvoScan'].run()

            self.data['new_image'] = self.scripts['GalvoScan'].data['image_data']

            dx_voltage, dy_voltage, self.data['correlation_image'] = correlation(self.data['baseline_image'],
                                                   self.data['image_extent'], self.data['new_image'],
                                                   self.data['image_extent'], use_trackpy=self.settings['use_trackpy'])

            self.data['shift'] = [dx_voltage, dy_voltage]
            print((self.data['shift']))

        else:
            self.scripts['GalvoScan'].run()
            self.data['baseline_image'] = self.scripts['GalvoScan'].data['image_data']
            self.data['image_extent'] = self.scripts['GalvoScan'].data['extent']

    def _plot(self, axes_list):
        '''
        Plots the newly taken galvo scan to axis 2, and the correlation image to axis 1
        Args:
            axes_list: list of axes to plot to. Uses two axes.

        '''
        if self.scripts['GalvoScan'].is_running:
            self.scripts['GalvoScan']._plot(axes_list)
        else:
            if not self.data['new_image'] == [] and not self.data['image_extent'] == []:
                data = self.data['new_image']
                extent = self.data['image_extent']
                plot_fluorescence_new(data, extent, axes_list[1])
                if not self.data['correlation_image'] == []:
                    axes_list[0].imshow(self.data['correlation_image'])
            else:
                self.scripts['GalvoScan']._plot(axes_list)

    def _update_plot(self, axes_list):
        '''
        Plots the newly taken galvo scan to axis 2, and the correlation image to axis 1
        Args:
            axes_list: list of axes to plot to. Uses two axes.

        '''
        if self.scripts['GalvoScan'].is_running:
            self.scripts['GalvoScan']._update_plot(axes_list)
        else:
            if not self.data['new_image'] == [] and not self.data['image_extent'] == []:
                data = self.data['new_image']
                update_fluorescence(data, axes_list[1])
                if not self.data['correlation_image'] == []:
                    axes_list[0].imshow(self.data['correlation_image'])
            else:
                self.scripts['GalvoScan']._update_plot(axes_list)


class Track_Correlate_Images(Script):
    '''
Track_Correlate_Images:
1.) Reads the current position of the galvo mirror: pt_0.
2.) Takes a galvo scan, compares it to a previous galvo scan to find the relative shift: dp.
3.) Sets the position of the galvo mirror to its initial position plus the shift: pt_0 + dp
    '''

    _DEFAULT_SETTINGS = [
        Parameter('mode', 'plain', ['plain', 'edge_detection', 'trackpy'], 'mode for correlation algorithm: plain images, identify points using trackpy to filter out background from nv-images or edge detection'),
        Parameter('baseline_update_frequency', 0, int, 'Use the last acquired image as the baseline for the next run after x executions of the script. x = 0 never update. Tip: use x=1 to update baseline'),
        Parameter('display_processed_images', False, bool, 'Show processed images used for correlation insead of raw images')
    ]
    _INSTRUMENTS = {}
    _SCRIPTS = {'take_baseline_image': GalvoScan, 'take_new_image': GalvoScan}

    def __init__(self, instruments = None, name = None, settings = None, scripts = None, log_function = None, data_path = None):
        Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_function= log_function, data_path = data_path)

        self.data = {'baseline_image': [],
                     'baseline_extent':[],
                     'new_image': [],
                     'new_image_extent': [],
                     'initial_galvo_position':[],
                     'shift': [],
                     'correlation_image': []
                     }

        self.baseline_processed_image = self.data['baseline_image']
        self.new_processed_image = self.data['new_image']
        self.count_executions = 0 # counts how often the script has been updated

    def _function(self):
        """
        see class doc string
        """
        def update_baseline():
            """
            update baseline image from the subscript data
            """
            self.scripts['take_baseline_image'].run()
            self.data['baseline_image'] = deepcopy(self.scripts['take_baseline_image'].data['image_data'])
            self.data['baseline_extent'] = deepcopy(self.scripts['take_baseline_image'].data['extent'])

        baseline_update_frequency = self.settings['baseline_update_frequency']
        self.count_executions += 1

        if self.settings['mode'] == 'plain':
            use_trackpy = False
            use_edge_detection = False
        elif self.settings['mode'] == 'edge_detection':
            use_trackpy = False
            use_edge_detection = True
        elif self.settings['mode'] == 'trackpy':
            use_trackpy = True
            use_edge_detection = False


        # get galvo position
        self.data['initial_galvo_position'] = np.array(self.scripts['take_new_image'].get_galvo_location())
        self.log('galvo at to x={:0.3f} y={:0.3f}'.format(self.data['initial_galvo_position'][0], self.data['initial_galvo_position'][1]))

        if self.data['baseline_image'] == []:
            update_baseline()

        else:
            self.scripts['take_new_image'].run()
            self.data['new_image'] = deepcopy(self.scripts['take_new_image'].data['image_data'])
            self.data['new_image_extent'] = deepcopy(self.scripts['take_new_image'].data['extent'])


            dx_voltage, dy_voltage, self.data['correlation_image'], self.baseline_processed_image, self.new_processed_image = correlation(self.data['baseline_image'],
                                                                                 self.data['baseline_extent'], self.data['new_image'],
                                                                                 self.data['new_image_extent'], use_trackpy=use_trackpy,
                                                                                 use_edge_detection=use_edge_detection
                                                                                 )
            self.data['shift'] = np.array((dx_voltage, dy_voltage))

            if baseline_update_frequency > 0 and self.count_executions % baseline_update_frequency == 0:
                self.log('updating baseline image')
                update_baseline()

            #need to set new galvo position after taking new baseline
            final_galvo_position = self.data['initial_galvo_position'] + self.data['shift']

            self.log('setting galvo to x={:0.3f} y={:0.3f}'.format(final_galvo_position[0], final_galvo_position[1]))
            # self.scripts['take_new_image'].set_galvo_location(final_galvo_position)

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
        axes_list = []
        if self._plot_refresh is True:
            figure_list[0].clf()
            axes_list.append(figure_list[0].add_subplot(111))
            figure_list[1].clf()
            axes_list.append(figure_list[1].add_subplot(121))
            axes_list.append(figure_list[1].add_subplot(122))
        else:
            for fig in figure_list:
                axes_list.append(fig.axes[0])

        return axes_list

    def _plot(self, axes_list):
        '''
        Plots the newly taken galvo scan to axis 2, and the correlation image to axis 1
        Args:
            axes_list: list of axes to plot to. Uses two axes.

        '''
        if self.scripts['take_baseline_image'].is_running:
            self.scripts['take_baseline_image']._plot(axes_list)
        elif self.scripts['take_new_image'].is_running:
            self.scripts['take_new_image']._plot(axes_list)
        else:
            #when clicking on script to select it for the first time, don't attempt to plot as no data avaliable
            if self.data['baseline_image'] == []:
                return

            if not self.settings['display_processed_images'] or self.baseline_processed_image == []:
                plot_fluorescence_new(self.data['baseline_image'], self.data['baseline_extent'], axes_list[0])
            else:
                plot_fluorescence_new(self.baseline_processed_image, self.data['baseline_extent'], axes_list[0])

            #for the first run when only the baseline is taken, plot only that image
            if self.data['new_image'] == []:
                return

            # add patch that marks the region of interest
            x, y = self.data['new_image_extent'][0], self.data['new_image_extent'][3]
            w, h = (self.data['new_image_extent'][1] - self.data['new_image_extent'][0]), (self.data['new_image_extent'][2] - self.data['new_image_extent'][3])
            patch = patches.Rectangle((x,y), w,h, ec='c', fc='none', ls='dashed')
            axes_list[0].add_patch(patch)

            # add patch that marks the region of interest
            x, y = self.data['new_image_extent'][0] - self.data['shift'][0], self.data['new_image_extent'][3] - self.data['shift'][1]
            w, h = (self.data['new_image_extent'][1] - self.data['new_image_extent'][0]), (self.data['new_image_extent'][2] - self.data['new_image_extent'][3])
            patch = patches.Rectangle((x,y), w,h, ec='r', fc='none', ls='dashed')
            axes_list[0].add_patch(patch)

            # plot correlation image
            if not self.data['correlation_image'] == []:
                # axes_list[1].imshow(self.data['correlation_image'], interpolation="nearest")
                plot_fluorescence_new(self.data['correlation_image'], [0,self.data['correlation_image'].shape[0], self.data['correlation_image'].shape[1], 0], axes_list[1])
                axes_list[1].set_title('correlation image')

            # plot new image
            if not self.data['new_image'] == []:
                if not self.settings['display_processed_images']:
                    plot_fluorescence_new(self.data['new_image'], self.data['new_image_extent'], axes_list[2])
                else:
                    plot_fluorescence_new(self.new_processed_image, self.data['new_image_extent'], axes_list[2])
                axes_list[2].set_title('new image shifted by dx={:0.3f} dy={:0.3f}'.format(self.data['shift'][0], self.data['shift'][1]))


    def _update_plot(self, axes_list):
        '''
        Plots the newly taken galvo scan to axis 2, and the correlation image to axis 1
        Args:
            axes_list: list of axes to plot to. Uses two axes.

        '''
        if self.scripts['take_baseline_image'].is_running:
            self.scripts['take_baseline_image']._update_plot(axes_list)
        elif self.scripts['take_new_image'].is_running:
            self.scripts['take_new_image']._update_plot(axes_list)
        else:
            self._plot(axes_list)

if __name__ == '__main__':
    script, failed, instr = Script.load_and_append({'Correlate_Images': Take_And_Correlate_Images})

    print(script)
    print(failed)
    print(instr)