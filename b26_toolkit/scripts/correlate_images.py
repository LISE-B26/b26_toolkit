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

from b26_toolkit.data_processing import correlation, shift_NVs, compare_galvos
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from pylabcontrol.core import Script, Parameter
from b26_toolkit.scripts import GalvoScan
from copy import deepcopy
import numpy as np
from matplotlib import patches
from b26_toolkit.instruments import PiezoController



class TakeAndCorrelateImages(Script):
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


class TrackCorrelateImages(Script):
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


class ImageCorrelation(Script):
    '''
    Takes a galvo scan, compares it to a previous galvo scan to find the relative shift
    '''

    _DEFAULT_SETTINGS = [
        Parameter('scan', 'baseline', ['baseline', 'comparison'],
                  'Baseline: take baseline galvo scan; Comparison: take comparison galvo scan and calculate correlation')
    ]

    _INSTRUMENTS = {}
    _SCRIPTS = {'GalvoScan': GalvoScan}

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

        if self.settings['scan'] == 'baseline':
            self.data['new_image'] = []
            self.data['correlation_image'] = []
            self.log('Taking baseline galvo scan')
            self.scripts['GalvoScan'].run()
            self.data['baseline_image'] = self.scripts['GalvoScan'].data['image_data']
            self.data['image_extent'] = self.scripts['GalvoScan'].data['extent']

        elif self.settings['scan'] == 'comparison':
            if self.data['baseline_image'] == []:
                self.log('No baseline image avaliable!')
            else:
                self.log('Taking comparison galvo scan')
                scan = self.scripts['GalvoScan']
                scan.settings['point_a']['x'] = self.data['image_extent'][0]
                scan.settings['point_b']['x'] = self.data['image_extent'][1]
                scan.settings['point_a']['y'] = self.data['image_extent'][3]
                scan.settings['point_b']['y'] = self.data['image_extent'][2]
                scan.settings['RoI_mode'] = 'corner'

                self.scripts['GalvoScan'].run()

                self.data['new_image'] = self.scripts['GalvoScan'].data['image_data']
                dx_voltage, dy_voltage, self.data['correlation_image'] = compare_galvos(self.data['baseline_image'],
                                                                                        self.data['image_extent'],
                                                                                        self.data['new_image'])

                self.data['shift'] = [dx_voltage, dy_voltage]
                self.log('Shift of Vx={:.3f},Vy={:.3f}V or Vx={:.3f},Vy={:.3f}um'
                         .format(*self.data['shift'],*np.array(self.data['shift'])*90))
                print((self.data['shift']))


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
            elif not self.data['baseline_image'] == []:
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
                # update_fluorescence(data, axes_list[1])
                update_fluorescence(data, axes_list[0])
                if not self.data['correlation_image'] == []:
                    axes_list[0].imshow(self.data['correlation_image'])
            elif not self.data['baseline_image'] == []:
                self.scripts['GalvoScan']._update_plot(axes_list)

class atto_compen(Script):
    '''
        Takes a galvo scan, compares it to a previous galvo scan to find the relative shift
        Optionally, also moves the attocube to compensate for any discovered shifts if they are above a threshold
        '''

    _DEFAULT_SETTINGS = [
        Parameter('toggle Attocube compensation', [
            Parameter('on/off', False, bool, 'If checked, apply DC voltage to Attocube to compensate for drifts'),
            Parameter('error_threshold', 0.0005, float, 'Minimum shift (V on galvo) for compensation to occur'),
            Parameter('adaptive_calibration', [
                Parameter('on/off', False, bool, 'If checked, adjust calibration by comparing images '
                                                           'before/after moving Attocubes'),
                Parameter('threshold', 0.3, float, 'Max allowed relative adjustment in calibration')]),
            Parameter('attocube_calibration', [
                Parameter('x_calib', .00101, float, 'V_y_galvo/V_x_attocube'),
                Parameter('y_calib', -.00090, float, 'V_x_galvo/V_y_attocube')]),
            Parameter('max_step_size', [
                Parameter('x', 15, float, 'Max step size for X attocube in V'),
                Parameter('y', 15, float, 'Max step size for Y attocube in V')])
            ])]

    _INSTRUMENTS = {'attocube_DC_controller': PiezoController}
    _SCRIPTS = {'ImageCorrelation': ImageCorrelation}

    def __init__(self, instruments=None, name=None, settings=None, scripts=None, log_function=None, data_path=None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings=settings, instruments=instruments, scripts=scripts,
                        log_function=log_function, data_path=data_path)

        self.data = {'baseline_image': [], 'new_image': [], 'image_extent': [], 'correlation_image': [],
                     'attocube_displacement': [], 'adjusted_calibration_x': None, 'adjusted_calibration_y': None}


    def _function(self):
        # Perform image correlation and obtain shifts in galvo voltage

        if self.scripts['ImageCorrelation'].data['baseline_image'] == []:
            self.log('No baseline image found. Take a baseline measurement first by running Image Correlation')

        elif self.scripts['ImageCorrelation'].settings['scan'] == 'baseline':
            self.log('Change ImageCorrelation''s scan type to ''comparison''!')
        else:
            self.scripts['ImageCorrelation'].run()
            self.data['shift'] = self.scripts['ImageCorrelation'].data['shift']
            self.data['new_image'] = self.scripts['ImageCorrelation'].data['new_image']
            self.data['image_extent'] = self.scripts['ImageCorrelation'].data['image_extent']
            self.data['baseline_image'] = self.scripts['ImageCorrelation'].data['baseline_image']
            self.data['correlation_image'] = self.scripts['ImageCorrelation'].data['correlation_image']

            if self.settings['toggle Attocube compensation']['on/off'] and not self.data['baseline_image'] == [] \
                    and np.max(np.abs(self.data['shift'])) > self.settings['toggle Attocube compensation']['error_threshold']:
                current_x, current_y = self.read_piezo()

                # Note that the galvo XY axis corresponds to the Attocube YX axis
                x_calib = self.settings['toggle Attocube compensation']['attocube_calibration']['x_calib']
                y_calib = self.settings['toggle Attocube compensation']['attocube_calibration']['y_calib']

                if not self.data['adjusted_calibration_x'] is None and not self.data['adjusted_calibration_y'] is None:
                    x_calib_new = self.data['adjusted_calibration_x']
                    y_calib_new = self.data['adjusted_calibration_y']

                    if self.settings['toggle Attocube compensation']['adaptive_calibration']['on/off']:
                        if max(np.abs(x_calib_new / x_calib - 1.), np.abs(y_calib_new / y_calib - 1.)) < \
                                self.settings['toggle Attocube compensation']['adaptive_calibration']['threshold']:
                            x_calib = x_calib_new
                            y_calib = y_calib_new
                        else:
                            print('Adjusted calibration differs too much from old one. Will revert to old calibration.')

                print('Current calibration: {:.5f},{:.5f}'.format(x_calib,y_calib))
                step_x = self.data['shift'][1]/x_calib
                step_y = self.data['shift'][0]/y_calib
                destination_x = float(current_x - step_x)
                destination_y = float(current_y - step_y)

                if 0 < destination_x < 60 and 0 < destination_y < 60 \
                        and step_x < self.settings['toggle Attocube compensation']['max_step_size']['x'] \
                        and step_y < self.settings['toggle Attocube compensation']['max_step_size']['y']:
                    # Ensure that the new voltage setting is within the Attocube's range, and that the step size is
                    # within reason (0.5 um), in case the correlation algorithm calculated an incorrectly big shift.

                    self.log('Target position (Attocube voltage): {:.2f},{:.2f}'.format(destination_x,destination_y))
                    self.move_piezo(destination_x, destination_y)
                else:
                    self.log('Target position out of range! Attocube will not be moved.')

                # Run correlation again to see if shift was eliminated
                self.scripts['ImageCorrelation'].run()
                self.data['shift_after_compensation'] = self.scripts['ImageCorrelation'].data['shift']

                self.data['adjusted_calibration_x'] = -(self.data['shift_after_compensation'][1]-self.data['shift'][1])/step_x
                self.data['adjusted_calibration_y'] = -(self.data['shift_after_compensation'][0]-self.data['shift'][0])/step_y



    def move_piezo(self, x, y):
        print(x, y)
        self.piezo = self.instruments['attocube_DC_controller']['instance']
        self.piezo.axis = 'x'
        self.piezo.voltage = x
        self.piezo.axis = 'y'
        self.piezo.voltage = y

    def read_piezo(self):
        self.piezo = self.instruments['attocube_DC_controller']['instance']
        self.piezo.axis = 'x'
        x_posn = self.piezo.read_probes('voltage')
        self.piezo.axis = 'y'
        y_posn = self.piezo.read_probes('voltage')

        return x_posn, y_posn

    def _plot(self, axes_list):
        '''
        Plots the newly taken galvo scan to axis 2, and the correlation image to axis 1
        Args:
            axes_list: list of axes to plot to. Uses two axes.

        '''
        if self.scripts['ImageCorrelation'].is_running:
            self.scripts['ImageCorrelation']._plot(axes_list)
        else:
            if not self.data['new_image'] == [] and not self.data['image_extent'] == []:
                data = self.data['new_image']
                extent = self.data['image_extent']
                plot_fluorescence_new(data, extent, axes_list[1])
                if not self.data['correlation_image'] == []:
                    axes_list[0].imshow(self.data['correlation_image'])
            elif not self.data['baseline_image'] == []:
                self.scripts['ImageCorrelation']._plot(axes_list)

    def _update_plot(self, axes_list):
        '''
        Plots the newly taken galvo scan to axis 2, and the correlation image to axis 1
        Args:
            axes_list: list of axes to plot to. Uses two axes.

        '''
        if self.scripts['ImageCorrelation'].is_running:
            self.scripts['ImageCorrelation']._update_plot(axes_list)
        else:
            if not self.data['new_image'] == [] and not self.data['image_extent'] == []:
                data = self.data['new_image']
                # update_fluorescence(data, axes_list[1])
                update_fluorescence(data, axes_list[0])
                if not self.data['correlation_image'] == []:
                    axes_list[0].imshow(self.data['correlation_image'])
            elif not self.data['baseline_image'] == []:
                self.scripts['ImageCorrelation']._update_plot(axes_list)

if __name__ == '__main__':
    script, failed, instr = Script.load_and_append({'Correlate_Images': TakeAndCorrelateImages})

    print(script)
    print(failed)
    print(instr)