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

from b26_toolkit.src.instruments import NI6259, NI9263, NI9402
from b26_toolkit.src.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from PyLabControl.src.core import Script, Parameter


class GalvoScan(Script):

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
                  [Parameter('x', 126, int, 'number of x points to scan'),
                   Parameter('y', 126, int, 'number of y points to scan')
                   ]),
        Parameter('time_per_pt', .002, [.0005, .001, .002, .005, .01, .015, .02], 'time in s to measure at each point'),
        Parameter('settle_time', .0002, [.0002], 'wait time between points to allow galvo to settle'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('DAQ_channels',
                   [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
                    Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
                    Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'], 'Daq channel used for counter')
                  ]),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'], 'return to the corn'),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    _INSTRUMENTS = {'NI6259':  NI6259, 'NI9263': NI9263, 'NI9402': NI9402}

    _SCRIPTS = {}

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
                        data_path=data_path)
        # defines which daqs contain the input and output based on user selection of daq interface
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['NI6259']['instance']
            self.daq_out = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_in = self.instruments['NI9402']['instance']
            self.daq_out = self.instruments['NI9263']['instance']

    def _function(self):
        """
        Executes threaded galvo scan
        """

        # update_time = datetime.datetime.now()

        # self._plot_refresh = True

        # self._plotting = True

        def init_scan():
            self._recording = False

            self.clockAdjust = int(
                (self.settings['time_per_pt'] + self.settings['settle_time']) / self.settings['settle_time'])

            [self.xVmin, self.xVmax, self.yVmax, self.yVmin] = self.pts_to_extent(self.settings['point_a'],
                                                                                  self.settings['point_b'],
                                                                                  self.settings['RoI_mode'])

            self.x_array = np.repeat(
                np.linspace(self.xVmin, self.xVmax, self.settings['num_points']['x'], endpoint=True),
                self.clockAdjust)
            self.y_array = np.linspace(self.yVmin, self.yVmax, self.settings['num_points']['y'], endpoint=True)
            sample_rate = float(1) / self.settings['settle_time']
            self.daq_out.settings['analog_output'][
                self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
            self.daq_out.settings['analog_output'][
                self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate
            self.daq_in.settings['digital_input'][
                self.settings['DAQ_channels']['counter_channel']]['sample_rate'] = sample_rate
            self.data = {'image_data': np.zeros((self.settings['num_points']['y'], self.settings['num_points']['x'])),
                         'bounds': [self.xVmin, self.xVmax, self.yVmin, self.yVmax]}

        if self.settings['daq_type'] == 'PCI':
            initial_position = self.daq_out.get_analog_voltages(
                [self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']])

        init_scan()
        self.data['extent'] = [self.xVmin, self.xVmax, self.yVmax, self.yVmin]

        for yNum in xrange(0, len(self.y_array)):

            if self._abort:
                break
            # set galvo to initial point of next line
            self.initPt = [self.x_array[0], self.y_array[yNum]]
            self.daq_out.set_analog_voltages(
                {self.settings['DAQ_channels']['x_ao_channel']: self.initPt[0],
                 self.settings['DAQ_channels']['y_ao_channel']: self.initPt[1]})

            # initialize APD thread
            ctrtask = self.daq_in.setup_counter(
                self.settings['DAQ_channels']['counter_channel'],
                len(self.x_array) + 1)
            aotask = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel']],
                                                                  self.x_array, ctrtask)

            # start counter and scanning sequence
            self.daq_out.run(aotask)
            self.daq_in.run(ctrtask)
            self.daq_out.waitToFinish(aotask)
            self.daq_out.stop(aotask)
            xLineData, _ = self.daq_in.read(ctrtask)
            self.daq_in.stop(ctrtask)
            diffData = np.diff(xLineData)

            summedData = np.zeros(len(self.x_array) / self.clockAdjust)
            for i in range(0, int((len(self.x_array) / self.clockAdjust))):
                summedData[i] = np.sum(
                    diffData[(i * self.clockAdjust + 1):(i * self.clockAdjust + self.clockAdjust - 1)])
            # also normalizing to kcounts/sec
            self.data['image_data'][yNum] = summedData * (.001 / self.settings['time_per_pt'])

            self.progress = float(yNum + 1) / len(self.y_array) * 100
            self.updateProgress.emit(int(self.progress))

        # set point after scan based on ending_behavior setting
        if self.settings['ending_behavior'] == 'leave_at_corner':
            return

        elif self.settings['ending_behavior'] == 'return_to_start':
            if self.settings['daq_type'] == 'PCI':
                self.set_galvo_location(initial_position)
            else:
                self.log('Could not determine initial position with this daq. Instead using leave_at_corner behavior')
                return

        elif self.settings['ending_behavior'] == 'return_to_origin':
            self.set_galvo_location([0, 0])

    def get_galvo_location(self):
        """
        Returns the current position of the galvo. Requires a daq with analog inputs internally routed to the analog
        outputs (ex. NI6259. Note that the cDAQ does not have this capability).
        Returns: list with two floats, which give the x and y position of the galvo mirror
        """
        galvo_position = self.daq_out.get_analog_voltages([
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
        # daq API only accepts either one point and one channel or multiple points and multiple channels
        pt = np.transpose(np.column_stack((pt[0], pt[1])))
        pt = (np.repeat(pt, 2, axis=1))

        task = self.daq_out.setup_AO(
            [self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], pt)
        self.daq_out.run(task)
        self.daq_out.waitToFinish(task)
        self.daq_out.stop(task)


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


    def _plot(self, axes_list, data=None):
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
        return super(GalvoScan, self).get_axes_layout([figure_list[0]])

class GalvoScan_cDAQ(GalvoScan):
    _INSTRUMENTS = {'daq_out': NI9263, 'daq_in': NI9402}

    def _function(self):
        """
        Executes threaded galvo scan
        """

        # update_time = datetime.datetime.now()

        # self._plot_refresh = True

        # self._plotting = True

        def init_scan():
            self._recording = False
            # self._abort = False

            self.clockAdjust = int(
                (self.settings['time_per_pt'] + self.settings['settle_time']) / self.settings['settle_time'])

            [self.xVmin, self.xVmax, self.yVmax, self.yVmin] = self.pts_to_extent(self.settings['point_a'],
                                                                                  self.settings['point_b'],
                                                                                  self.settings['RoI_mode'])

            self.x_array = np.repeat(
                np.linspace(self.xVmin, self.xVmax, self.settings['num_points']['x'], endpoint=True),
                self.clockAdjust)
            self.y_array = np.linspace(self.yVmin, self.yVmax, self.settings['num_points']['y'], endpoint=True)
            sample_rate = float(1) / self.settings['settle_time']
            self.instruments['daq_out']['instance'].settings['analog_output'][
                self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
            self.instruments['daq_out']['instance'].settings['analog_output'][
                self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate
            self.instruments['daq_in']['instance'].settings['digital_input'][
                self.settings['DAQ_channels']['counter_channel']]['sample_rate'] = sample_rate
            self.data = {'image_data': np.zeros((self.settings['num_points']['y'], self.settings['num_points']['x'])),
                         'bounds': [self.xVmin, self.xVmax, self.yVmin, self.yVmax]}

        # not possible for cDAQ
        # initial_position = self.instruments['daq']['instance'].get_analog_voltages(
        #     [self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']])

        init_scan()
        self.data['extent'] = [self.xVmin, self.xVmax, self.yVmax, self.yVmin]

        for yNum in xrange(0, len(self.y_array)):

            if self._abort:
                break
            # set galvo to initial point of next line
            self.initPt = [self.x_array[0], self.y_array[yNum]]
            self.instruments['daq_out']['instance'].set_analog_voltages(
                {self.settings['DAQ_channels']['x_ao_channel']: self.initPt[0],
                 self.settings['DAQ_channels']['y_ao_channel']: self.initPt[1]})

            # initialize APD thread
            ctrtask = self.instruments['daq_in']['instance'].setup_counter(
                self.settings['DAQ_channels']['counter_channel'],
                len(self.x_array) + 1)
            aotask = self.instruments['daq_out']['instance'].setup_AO([self.settings['DAQ_channels']['x_ao_channel']],
                                                                  self.x_array, ctrtask)
            # aotask = self.instruments['daq_out']['instance'].setup_AO([self.settings['DAQ_channels']['x_ao_channel']],
            #                                                           self.x_array)

            # start counter and scanning sequence
            self.instruments['daq_out']['instance'].run(aotask)
            self.instruments['daq_in']['instance'].run(ctrtask)
            self.instruments['daq_out']['instance'].waitToFinish(aotask)
            self.instruments['daq_out']['instance'].stop(aotask)
            xLineData, _ = self.instruments['daq_in']['instance'].read(ctrtask)
            self.instruments['daq_in']['instance'].stop(ctrtask)
            diffData = np.diff(xLineData)

            summedData = np.zeros(len(self.x_array) / self.clockAdjust)
            for i in range(0, int((len(self.x_array) / self.clockAdjust))):
                summedData[i] = np.sum(
                    diffData[(i * self.clockAdjust + 1):(i * self.clockAdjust + self.clockAdjust - 1)])
            # also normalizing to kcounts/sec
            self.data['image_data'][yNum] = summedData * (.001 / self.settings['time_per_pt'])

            self.progress = float(yNum + 1) / len(self.y_array) * 100
            self.updateProgress.emit(int(self.progress))

        # set point after scan based on ending_behavior setting
        if self.settings['ending_behavior'] == 'leave_at_corner':
            return
        elif self.settings['ending_behavior'] == 'return_to_start':
            # pt = (self.settings['point_a']['x'], self.settings['point_a']['y'])
            # pt = np.transpose(np.column_stack((pt[0], pt[1])))
            # pt = (np.repeat(initial_position, 2, axis=1))
            # self.set_galvo_location(initial_position)
            return

        elif self.settings['ending_behavior'] == 'return_to_origin':
            self.set_galvo_location([0, 0])
            # pt = (0, 0)
            # pt = np.transpose(np.column_stack((pt[0], pt[1])))
            # pt = (np.repeat(pt, 2, axis=1))

            # self.instruments['daq']['instance'].setup_AO(
            #     [self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], pt)
            # self.instruments['daq']['instance'].AO_run()
            # self.instruments['daq']['instance'].AO_waitToFinish()
            # self.instruments['daq']['instance'].AO_stop()

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
        daq = self.instruments['daq_out']['instance']
        # daq API only accepts either one point and one channel or multiple points and multiple channels
        pt = np.transpose(np.column_stack((pt[0], pt[1])))
        pt = (np.repeat(pt, 2, axis=1))

        task = daq.setup_AO(
            [self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], pt)
        daq.run(task)
        daq.waitToFinish(task)
        daq.stop(task)





if __name__ == '__main__':
    script, failed, instruments = Script.load_and_append(script_dict={'GalvoScan': 'GalvoScan'})

    print(script)
    print(failed)
    # print(instruments)

