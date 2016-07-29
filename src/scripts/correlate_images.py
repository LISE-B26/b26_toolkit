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

from b26_toolkit.src.data_processing import correlation, shift_NVs
from b26_toolkit.src.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from PyLabControl.src.core import Script, Parameter
from b26_toolkit.src.scripts import GalvoScan


class Take_And_Correlate_Images_2(Script):
    '''
    Takes a galvo scan, compares it to a previous galvo scan to find the relative shift, and then updates a list of
    nvs based on this shift so that they will give the current coordinates of those nvs
    '''

    _DEFAULT_SETTINGS = [
        Parameter('use_trackpy', False, bool, 'Use trackpy to create artificial nv-only images to filter out background')
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

        self.data = {'baseline_image': [], 'new_image': [], 'image_extent': [], 'old_nv_list':[], 'new_NV_list': [], 'correlation_image': []}

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
        elif self.data['old_nv_list'] == []:
            self.log('No nv list avaliable. Scipt may have been run in error.')

        if not self.data['baseline_image'] == []:
            #use same settings as initial image
            scan = self.scripts['GalvoScan'].scripts['acquire_image']
            scan.settings['point_a']['x'] = self.data['image_extent'][0]
            scan.settings['point_b']['x'] = self.data['image_extent'][1]
            scan.settings['point_a']['y'] = self.data['image_extent'][3]
            scan.settings['point_b']['y'] = self.data['image_extent'][2]

            self.scripts['GalvoScan'].run()

            self.data['new_image'] = self.scripts['GalvoScan'].scripts['acquire_image'].data['image_data']

            dx_voltage, dy_voltage, self.data['correlation_image'] = correlation(self.data['baseline_image'],
                                                   self.data['image_extent'], self.data['new_image'],
                                                   self.data['image_extent'], use_trackpy=self.settings['use_trackpy'])

            self.data['new_NV_list'] = shift_NVs(dx_voltage, dy_voltage, self.data['old_nv_list'])

        else:
            self.scripts['GalvoScan'].run()
            self.data['baseline_image'] = self.scripts['GalvoScan'].data['image_data']

    def _plot(self, axes_list):
        '''
        Plots the newly taken galvo scan to axis 2, and the correlation image to axis 1
        Args:
            axes_list: list of axes to plot to. Uses two axes.

        '''
        data = self.scripts['GalvoScan'].data['image_data']
        extent = self.scripts['GalvoScan'].data['extent']
        plot_fluorescence_new(data, extent, axes_list[1])
        if not self.data['correlation_image'] == []:
            axes_list[0].imshow(self.data['correlation_image'])

    def _update_plot(self, axes_list):
        '''
        Plots the newly taken galvo scan to axis 2, and the correlation image to axis 1
        Args:
            axes_list: list of axes to plot to. Uses two axes.

        '''
        data = self.scripts['GalvoScan'].data['image_data']
        extent = self.scripts['GalvoScan'].data['extent']
        update_fluorescence(data, axes_list[1])
        if not self.data['correlation_image'] == []:
            axes_list[0].imshow(self.data['correlation_image'])


if __name__ == '__main__':
    script, failed, instr = Script.load_and_append({'Correlate_Images': 'Correlate_Images'})

    print(script)
    print(failed)
    print(instr)
