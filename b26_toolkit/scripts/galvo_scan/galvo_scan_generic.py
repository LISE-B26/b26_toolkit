"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    pylabcontrol is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    pylabcontrol is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with pylabcontrol.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import time
from b26_toolkit.instruments import NI6259
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from pylabcontrol.core import Script, Parameter

class GalvoScanGeneric(Script):
    """
    GalvoScan uses the apd, daq, and galvo to sweep across voltages while counting photons at each voltage,
    resulting in an image in the current field of view of the objective.
    """

    _DEFAULT_SETTINGS = [
        Parameter('point_a',
                  [Parameter('x', 0, float, 'x-coordinate'),
                   Parameter('y', 0, float, 'y-coordinate')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 1.0, float, 'x-coordinate'),
                   Parameter('y', 1.0, float, 'y-coordinate')
                   ]),
        Parameter('RoI_mode', 'center', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                                           corner: pta and ptb are diagonal corners of rectangle.\n \
                                                           center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 25, int, 'number of x points to scan'),
                   Parameter('y', 25, int, 'number of y points to scan')
                   ]),
        Parameter('time_per_pt', .002, [.0001,.001, .002, .005, .01, .015, .02], 'time in s to measure at each point'),
        Parameter('settle_time', .0002, [.0002], 'wait time between points to allow galvo to settle in seconds'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        # Parameter('DAQ_channels',
        #            [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
        #             Parameter('y_ao_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
        #             Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1'], 'Daq channel used for counter')
        #           ]),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'], 'return to the corn')
    ]

    _INSTRUMENTS = {}
    _SCRIPTS = {}

    _ACQ_TYPE = 'line' #this defines if the galvo acquisition is line by line or point by point, the default is line

    def __init__(self, instruments, name=None, settings=None, log_function=None, data_path=None):
        '''
        Initializes GalvoScan script for use in gui

        Args:
            instruments: list of instrument objects
            name: name to give to instantiated script object
            settings: dictionary of new settings to pass in to override defaults
            log_function: log function passed from the gui to direct log calls to the gui log
            data_path: path to save data

        '''
        Script.__init__(self, name, settings=settings, instruments=instruments, log_function=log_function,
                        data_path = data_path)


    def setup_scan(self):
        """
        prepares the scan
        Returns:

        """
        pass


    def _function(self):
        """
        Executes threaded galvo scan
        """




        self.data = {'image_data': np.zeros((self.settings['num_points']['y'], self.settings['num_points']['x']))}
        self.data['extent'] = self.pts_to_extent(self.settings['point_a'], self.settings['point_b'], self.settings['RoI_mode'])

        [xVmin, xVmax, yVmax, yVmin] = self.data['extent']
        self.x_array = np.linspace(xVmin, xVmax, self.settings['num_points']['x'], endpoint=True)
        self.y_array = np.linspace(yVmin, yVmax, self.settings['num_points']['y'], endpoint=True)

        if self._ACQ_TYPE == 'point':
            self.data['point_data'] = [] # stores the complete data acquired at each point, image_data holds only a scalar at each point

        self.setup_scan()

        # initial_position = self.instruments['daq']['instance'].get_analog_out_voltages([self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']])
        initial_position = self.get_galvo_location()
        if initial_position == []:
            print('WARNING!! GALVO POSITION COULD NOT BE DETERMINED. SET ENDING ending_behavior TO leave_at_corner')
            self.settings['ending_behavior'] = 'leave_at_corner'

        Nx, Ny = self.settings['num_points']['x'], self.settings['num_points']['y']

        for yNum in range(0, Ny):

            if self._ACQ_TYPE == 'line':
                if self._abort:
                    break
                line_data = self.read_line(self.y_array[yNum])
                self.data['image_data'][yNum] = line_data
                self.progress = float(yNum + 1) / Ny * 100
                self.updateProgress.emit(int(self.progress))

            elif self._ACQ_TYPE == 'point':
                for xNum in range(0, Nx):
                    if self._abort:
                        break

                    point_data = self.read_point(self.x_array[xNum], self.y_array[yNum])
                    self.data['image_data'][yNum, xNum] = np.mean(point_data)

                    self.data['point_data'].append(point_data)
                    self.progress = float(yNum * Nx + 1 + xNum) / (Nx * Ny) * 100

                    # JG: tmp print info about progress
                    print(('current acquisition {:02d}/{:02d} ({:0.2f}%)'.format(yNum * Nx + xNum, Nx * Ny, self.progress)))

                    self.updateProgress.emit(int(self.progress))

                # fill the rest of the array with the mean of the data up to now (otherwise it's zero and the data is not visible in the plot)
                if yNum<Ny:
                    self.data['image_data'][yNum + 1:, :] = np.mean(self.data['image_data'][0:yNum, :].flatten())

        #set point after scan based on ending_behavior setting
        if self.settings['ending_behavior'] == 'leave_at_corner':
            return
        elif self.settings['ending_behavior'] == 'return_to_start':
            self.set_galvo_location(initial_position)
        elif self.settings['ending_behavior'] == 'return_to_origin':
            self.set_galvo_location([0,0])

    def get_galvo_location(self):
        """
        returns the current position of the galvo
        Returns: list with two floats, which give the x and y position of the galvo mirror
        """
        raise NotImplementedError
        return galvo_position

    def set_galvo_location(self, galvo_position):
        """
        sets the current position of the galvo
        galvo_position: list with two floats, which give the x and y position of the galvo mirror
        """
        raise NotImplementedError

    def read_line(self, y_pos):
        """
        reads a line of data from the DAQ, this function is used if _ACQ_TYPE = 'line'
        Args:
            y_pos: y position of the scan

        Returns:

        """
        raise NotImplementedError

    def read_point(self, x_pos, y_pos):
        """
        reads a line of data from the DAQ, this function is used if _ACQ_TYPE = 'point'
        Args:
            x_pos: x position of the scan
            y_pos: y position of the scan
        Returns:

        """
        raise NotImplementedError

    @staticmethod
    def pts_to_extent(pta, ptb, roi_mode):
        """

        Args:
            pta: point a
            ptb: point b
            roi_mode:   mode how to calculate region of interest
                        corner: pta and ptb are diagonal corners of rectangle.
                        center: pta is center and ptb is extend or rectangle

        Returns: extend of region of interest [xVmin, xVmax, yVmax, yVmin]

        """
        if roi_mode == 'corner':
            xVmin = min(pta['x'], ptb['x'])
            xVmax = max(pta['x'], ptb['x'])
            yVmin = min(pta['y'], ptb['y'])
            yVmax = max(pta['y'], ptb['y'])
        elif roi_mode == 'center':
            xVmin = pta['x'] - float(ptb['x']) / 2.
            xVmax = pta['x'] + float(ptb['x']) / 2.
            yVmin = pta['y'] - float(ptb['y']) / 2.
            yVmax = pta['y'] + float(ptb['y']) / 2.
        return [xVmin, xVmax, yVmax, yVmin]

    def _plot(self, axes_list, data = None):
        """
        Plots the galvo scan image
        Args:
            axes_list: list of axes objects on which to plot the galvo scan on the first axes object
            data: data (dictionary that contains keys image_data, extent) if not provided use self.data
        """

        if data is None:
            data = self.data

        plot_fluorescence_new(data['image_data'], data['extent'], axes_list[0], max_counts=self.settings['max_counts_plot'])

    def _update_plot(self, axes_list):
        """
        updates the galvo scan image
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
        """
        update_fluorescence(self.data['image_data'], axes_list[0], self.settings['max_counts_plot'])

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

        # only pick the first figure from the figure list, this avoids that get_axes_layout clears all the figures
        return super(GalvoScanGeneric, self).get_axes_layout([figure_list[0]])


# uncommented by Jan March 12th 2018. I think this is not used anymore
class GalvoScanDAQ(GalvoScanGeneric):
    """
    GalvoScan uses the apd, daq, and galvo to sweep across voltages while counting photons at each voltage,
    resulting in an image in the current field of view of the objective.
    """

    _DEFAULT_SETTINGS = [
        Parameter('DAQ_channels',
                   [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
                    Parameter('y_ao_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
                    Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1'], 'Daq channel used for counter')
                  ])
    ]

    _INSTRUMENTS = {'daq':  NI6259}

    def __init__(self, instruments, name=None, settings=None, log_function=None, data_path=None):
        '''
        Initializes GalvoScan script for use in gui

        Args:
            instruments: list of instrument objects
            name: name to give to instantiated script object
            settings: dictionary of new settings to pass in to override defaults
            log_function: log function passed from the gui to direct log calls to the gui log
            data_path: path to save data

        '''
        self._DEFAULT_SETTINGS = self._DEFAULT_SETTINGS + GalvoScanGeneric._DEFAULT_SETTINGS
        Script.__init__(self, name, settings=settings, instruments=instruments, log_function=log_function, data_path = data_path)

    def setup_scan(self):

        self.clockAdjust = int((self.settings['time_per_pt'] + self.settings['settle_time']) / self.settings['settle_time'])

        [xVmin, xVmax, yVmax, yVmin] = self.pts_to_extent(self.settings['point_a'],self.settings['point_b'],self.settings['RoI_mode'])


        self.x_array = np.repeat(np.linspace(xVmin, xVmax, self.settings['num_points']['x'], endpoint=True),self.clockAdjust)
        self.y_array = np.linspace(yVmin, yVmax, self.settings['num_points']['y'], endpoint=True)
        sample_rate = float(1) / self.settings['settle_time']
        self.instruments['daq']['instance'].settings['analog_output'][
            self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
        self.instruments['daq']['instance'].settings['analog_output'][
            self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate
        self.instruments['daq']['instance'].settings['digital_input'][
            self.settings['DAQ_channels']['counter_channel']]['sample_rate'] = sample_rate

    def get_galvo_location(self):
        """
        returns the current position of the galvo
        Returns: list with two floats, which give the x and y position of the galvo mirror
        """
        galvo_position = self.instruments['daq']['instance'].get_analog_voltages([
            self.settings['DAQ_channels']['x_ao_channel'],
            self.settings['DAQ_channels']['y_ao_channel']]
        )
        return galvo_position

    def set_galvo_location(self, galvo_position):
        """
        sets the current position of the galvo
        galvo_position: list with two floats, which give the x and y position of the galvo mirror
        """
        if galvo_position[0] > 1 or galvo_position[0] < -1 or galvo_position[1] > 1 or galvo_position[1] < -1:
            raise ValueError('The script attempted to set the galvo position to an illegal position outside of +- 1 V')

        pt = galvo_position
        daq = self.instruments['daq']['instance']
        # daq API only accepts either one point and one channel or multiple points and multiple channels
        pt = np.transpose(np.column_stack((pt[0],pt[1])))
        pt = (np.repeat(pt, 2, axis=1))

        daq.setup_AO([self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], pt)
        daq.AO_run()
        daq.AO_waitToFinish()
        daq.AO_stop()

    def read_line(self, y_pos):
        """
        reads a line of data from the DAQ
        Args:
            y_pos: y position of the scan

        Returns:

        """
        # initialize APD thread
        clk_source = self.instruments['daq']['instance'].DI_init(self.settings['DAQ_channels']['counter_channel'],
                                                                 len(self.x_array) + 1)
        self.initPt = np.transpose(np.column_stack((self.x_array[0],y_pos)))
        self.initPt = (np.repeat(self.initPt, 2, axis=1))

        # move galvo to first point in line
        self.instruments['daq']['instance'].setup_AO(
            [self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']],
            self.initPt, "")
        self.instruments['daq']['instance'].AO_run()
        self.instruments['daq']['instance'].AO_waitToFinish()
        self.instruments['daq']['instance'].AO_stop()
        self.instruments['daq']['instance'].setup_AO([self.settings['DAQ_channels']['x_ao_channel']], self.x_array,
                                                     clk_source)
        # start counter and scanning sequence
        self.instruments['daq']['instance'].AO_run()
        self.instruments['daq']['instance'].DI_run()
        self.instruments['daq']['instance'].AO_waitToFinish()
        self.instruments['daq']['instance'].AO_stop()
        xLineData, _ = self.instruments['daq']['instance'].DI_read()
        self.instruments['daq']['instance'].DI_stop()
        diffData = np.diff(xLineData)

        summedData = np.zeros(len(self.x_array) / self.clockAdjust)
        for i in range(0, int((len(self.x_array) / self.clockAdjust))):
            summedData[i] = np.sum(diffData[(i * self.clockAdjust + 1):(i * self.clockAdjust + self.clockAdjust - 1)])
        # also normalizing to kcounts/sec
        return summedData * (.001 / self.settings['time_per_pt'])



if __name__ == '__main__':
    from pylabcontrol.core import Instrument
    # from b26_toolkit.pylabcontrol.instruments import NI7845RMain
    #
    # fpga = NI7845RMain()
    #
    #
    # g = GalvoScanFPGA(instruments={'NI7845RMain':fpga}, name='test_fpga_scan', settings=None, log_function=None, data_path=None)
    # print(fpga)


    # instruments, failed =  Instrument.load_and_append(instrument_dict ={'NI7845RMain': 'NI7845RMain'}, raise_errors=True )

    script, failed, instruments = Script.load_and_append(script_dict={'GalvoScanFPGA': 'GalvoScanFPGA'}, raise_errors=True)
    #
    print(script)
    print(failed)
# # print(instruments)

