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
import time
from b26_toolkit.scripts.galvo_scan.galvo_scan import GalvoScan
from b26_toolkit.instruments.piezo_controller import PiezoController
from b26_toolkit.instruments import B26KDC001x, B26KDC001z, B26KDC001y, NI9219
from b26_toolkit.scripts.set_laser import SetLaser, SetLaserSingleAxis
from b26_toolkit.scripts.daq_read_ai import Daq_Read_Analog
from b26_toolkit.scripts.pulse_sequences.laser_pulses import PulsedReadout
from b26_toolkit.scripts.pulse_sequences.param_sweep.pulsed_esr import PulsedEsrFast
from pylabcontrol.core import Parameter, Script
from b26_toolkit.scripts import FindNv, EsrSimple
from b26_toolkit.scripts.find_nv import FindNvStrobe
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from b26_toolkit.plotting.plots_1d import plot_1d_simple_freq
from b26_toolkit.scripts.set_laser import SetAtto
from b26_toolkit.data_analysis.nv_optical_response import B_field_from_esr

nv_gyro = 2.8025e6

class AttoGridScanGeneric(Script):

    """
    Base script for moving a scanner (e.g. Attocube) in a grid and taking a measurement at each point.
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
        Parameter('RoI_mode', 'corner', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                               corner: pta and ptb are diagonal corners of rectangle.\n \
                                               center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 5, int, 'number of x points to scan, if 1 then perform a line scan along the other axis'),
                   Parameter('y', 5, int, 'number of y points to scan, if 1 then perform a line scan along the other axis')
                   ]),
        Parameter('settle_time', 0.5, float, 'time (s) to wait after moving to a new point, before taking data'),
        Parameter('Tracking', [
            Parameter('on/off', True, bool, 'used to turn on tracking'),
            Parameter('every_N', 1, int, 'track every n points')]),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn')
    ]

    _SCRIPTS = {'find_nv': FindNv, 'ESR_simple': EsrSimple}
    _INSTRUMENTS = {}

    def _function(self):

        try:
            self.setup_scan()
        except AttributeError:
            return

        Nx, Ny = self.settings['num_points']['x'], self.settings['num_points']['y']
        self.data['actual_Vx'] = np.zeros(Nx * Ny)
        self.data['actual_Vy'] = np.zeros(Nx * Ny)

        if self.settings['num_points']['y'] == 1:
            self.linescan_axis = 'x'
        elif self.settings['num_points']['x'] == 1:
            self.linescan_axis = 'y'
        else:
            self.linescan_axis = '2d'

        self._initialize_data()  # Initialize self.data

        [xVmin, xVmax, yVmax, yVmin] = self.data['extent']
        if self.linescan_axis == 'x' and self.settings['RoI_mode'] == 'corner' and self.settings['point_a']['x'] > self.settings['point_b']['x']:
            xVmin, xVmax = self.settings['point_a']['x'],  self.settings['point_b']['x']
        elif self.linescan_axis == 'y' and self.settings['RoI_mode'] == 'corner' and self.settings['point_a']['y'] > self.settings['point_b']['y']:
            yVmin, yVmax = self.settings['point_a']['y'], self.settings['point_b']['y']

        self.x_array = np.linspace(xVmin, xVmax, self.settings['num_points']['x'], endpoint=True)
        self.y_array = np.linspace(yVmin, yVmax, self.settings['num_points']['y'], endpoint=True)

        self.data['x_array'] = self.x_array
        self.data['y_array'] = self.y_array

        try:
            self.check_bounds()
        except AttributeError:
            self._abort = True

        self.x_array = self.x_array.tolist()
        self.y_array = self.y_array.tolist()

        self.point_index = 0

        xNum_array = list(range(0, Nx))
        pixels_completed = 0
        for yNum in range(0, Ny):
            self.move_piezo(0, self.y_array[yNum], direction='y')
            for xNum in xNum_array:
                if self._abort:
                    break

                # print('Moving piezo to %.1f, %.1f' % ((self.x_array[xNum], self.y_array[yNum])))
                self.move_piezo(self.x_array[xNum], self.y_array[yNum], direction='x')
                time.sleep(self.settings['settle_time'])

                self.progress = float(yNum * Nx + xNum) / (Nx * Ny) * 100
                self.updateProgress.emit(int(self.progress))
                self.point_index += 1
                if 'Tracking' in self.settings:
                    if self.settings['Tracking']['on/off'] and int(self.point_index) % self.settings['Tracking']['every_N'] == 0:
                        self.scripts['find_nv'].run(verbose=False)

                        # Move piezo to the set location again. For some reason, FindNv seems to affect the DAQ outputs to the Attocube
                        print('Moving piezo (again) to %.1f, %.1f' % (self.x_array[xNum], self.y_array[yNum]))
                        self.move_piezo(self.x_array[xNum], self.y_array[yNum])
                        time.sleep(.5)

                point_value = self.read_point(yNum, xNum)
                self.data['point_value'][yNum, xNum] = point_value

                pixels_completed += 1
                self.progress = float(pixels_completed) / (Nx * Ny) * 100
                self.updateProgress.emit(int(self.progress))

            if 'scan_mode' in self.settings and self.settings['scan_mode'] == 'meander':
                xNum_array = xNum_array[::-1]

        # set end position after scan based on ending_behavior setting
        if self.settings['ending_behavior'] == 'leave_at_corner':
            return
        elif self.settings['ending_behavior'] == 'return_to_start':
            self.move_piezo(self.settings['point_a']['x'], self.settings['point_a']['y'])
        elif self.settings['ending_behavior'] == 'return_to_origin':
            self.move_piezo(0, 0)

    def _initialize_data(self):
        self.data = {'point_value': np.zeros((self.settings['num_points']['y'], self.settings['num_points']['x'])),
                     'extent': self.pts_to_extent(self.settings['point_a'], self.settings['point_b'], self.settings['RoI_mode'])}

    def check_bounds(self):
        # Check if there are negative voltages in the scan range
        negative_voltage_err = 'Piezo voltage cannot be < 0 V'
        if np.any(self.x_array < 0):
            raise ValueError(negative_voltage_err)
        elif np.any(self.y_array < 0):
            raise ValueError(negative_voltage_err)

    def setup_scan(self):
        pass

    def _read_point(self, yNum, xNum):
        """
        Wrapper function for read_point(). Normally just runs read_point(), but for things like gradient measurements, we can use
        the wrapper function to run read_point with different parameters and extract a single value based on the multiple measurements
        :param yNum: index number along y axis in the scan grid
        :param xNum: index number along x axis in the scan grid
        :return: None
        """
        return self.read_point(yNum, xNum)

    def read_point(self, yNum, xNum):
        """
        Replace this with some data-taking script, e.g. read the counts (to perform a galvoscan with the attocubes
        instead) or take ESR (for a 2D scan of B field)
        Returns: measured value at given pt, with other data that won't be plotted (e.g. for ESRs, we want to plot
        only the on-axis field, and save but not plot the off-axis field)
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
        label = ['B field scan w/ Attocube', r'V$_x$ [V]', r'V$_y$ [V]', 'On-axis field (G)']

        if self.linescan_axis == '2d':
            plot_fluorescence_new(data['point_value'], data['extent'], axes_list[0], max_counts=-1, labels=label)
        elif self.linescan_axis == 'x' or self.linescan_axis == 'y':
            raise NotImplementedError


    def _update_plot(self, axes_list):
        """
        updates the galvo scan image
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
        """
        update_fluorescence(self.data['point_value'], axes_list[0], -1)

    def move_piezo(self, x, y, direction='both'):
        """
        Action to move the scanner
        :param x: location along x that scanner should move to
        :param y: location along y that scanner should move to
        :return: None
        """
        raise NotImplementedError


class AttoGridScanEsr(AttoGridScanGeneric):
    """
    Moves Attocubes in a grid and takes EsrSimple (CW) at each point, resulting in a magnetic field scan
    """
    _DEFAULT_SETTINGS = [
        Parameter('zero_field_splitting', 2.8707e9, float, 'Dip in ESR with no external B field applied'),
        Parameter('point_a',
                  [Parameter('x', 0, float, 'x-coordinate'),
                   Parameter('y', 0, float, 'y-coordinate')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 1.0, float, 'x-coordinate'),
                   Parameter('y', 1.0, float, 'y-coordinate')
                   ]),
        Parameter('RoI_mode', 'corner', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                                   corner: pta and ptb are diagonal corners of rectangle.\n \
                                                   center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 5, int, 'number of x points to scan'),
                   Parameter('y', 5, int, 'number of y points to scan')
                   ]),
        Parameter('settle_time', 0.5, float, 'time (s) to wait after moving to a new point, before taking data'),
        Parameter('Tracking', [
            Parameter('on/off', True, bool, 'used to turn on tracking'),
            Parameter('every_N', 1, int, 'track every n points')]),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn')
    ]

    _SCRIPTS = {'find_nv': FindNv, 'ESR_simple': EsrSimple, 'SetAttoANC300': SetAtto}

    def _initialize_data(self):
        self.data = {'point_value': np.zeros((self.settings['num_points']['y'], self.settings['num_points']['x'])),
                     'extent': self.pts_to_extent(self.settings['point_a'], self.settings['point_b'], self.settings['RoI_mode']),
                     'point_freq': np.zeros((self.settings['num_points']['y'], self.settings['num_points']['x']))}
        if 'measure_gradient' in self.settings and self.settings['measure_gradient']['enable']:
            self.data['point_freq_2'] = np.zeros((self.settings['num_points']['y'], self.settings['num_points']['x']))

    def move_piezo(self, x, y, direction='both'):
        self.scripts['SetAttoANC300'].settings['point'].update({'x': float(x), 'y': float(y)})
        self.scripts['SetAttoANC300'].run()

    def read_point(self, yNum, xNum):
        """
        Take ESR (for a 2D scan of B field)
        Returns: measured value at given pt, with other data that won't be plotted
        """
        self.scripts['ESR_simple'].settings['tag'] = 'esr_simple_%i,%i' % (yNum, xNum)

        if self.settings['track_esr']:
            if yNum == 0 and xNum == 0:
                # For origin (where scan begins), simply use ESR subscript settings to set center freq
                pass
            elif xNum > 0:
                # Center current ESR scan range around extracted frequency of the left adjacent (previous) point
                predicted_freq = self.data['point_freq'][yNum][xNum - 1]
                if predicted_freq != 0:
                    print('Predicted freq: %.3e using pt: %i, %i' % (predicted_freq, yNum, xNum - 1))
                    self.scripts['ESR_simple'].settings['range_type'] = 'center_range'
                    self.scripts['ESR_simple'].settings['freq_start'] = float(predicted_freq)
            elif xNum == 0:
                # Center current ESR scan range around extracted frequency of the upper adjacent (previous) point
                predicted_freq = self.data['point_freq'][yNum - 1][xNum]
                if predicted_freq != 0:
                    print('Predicted freq: %.3e using pt: %i, %i' % (predicted_freq, yNum - 1, xNum))
                    self.scripts['ESR_simple'].settings['range_type'] = 'center_range'
                    self.scripts['ESR_simple'].settings['freq_start'] = float(predicted_freq)

        self.scripts['ESR_simple'].run()

        try:
            point_data = self.scripts['ESR_simple'].data['fits']
        except:
            point_data = None

        if point_data is not None:

            if len(point_data) == 6:
                fp = point_data[5]
                fn = point_data[4]
                self.data['point_freq'][yNum][xNum] = point_data[4]
                point_value = B_field_from_esr(fp, fn, D=self.settings['zero_field_splitting'],
                                               gamma=27.969e9, angular_freq=False, verbose=False)[0] * 1e4
            elif len(point_data) == 4:
                self.data['point_freq'][yNum][xNum] = point_data[3]
                point_value = (point_data[3] - self.settings['zero_field_splitting']) / nv_gyro
        else:
            self.data['point_freq'][yNum][xNum] = 0
            point_value = 0

        return point_value

        return point_value


class AttoGridScanPulsedEsr(AttoGridScanEsr):
    """
    Moves Attocubes in a grid and takes PulsedEsr at each point, resulting in a magnetic field scan
    """
    _DEFAULT_SETTINGS = [
        Parameter('zero_field_splitting', 2.8707e9, float, 'Dip in ESR with no external B field applied'),
        Parameter('point_a',
                  [Parameter('x', 0, float, 'x-coordinate'),
                   Parameter('y', 0, float, 'y-coordinate')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 1.0, float, 'x-coordinate'),
                   Parameter('y', 1.0, float, 'y-coordinate')
                   ]),
        Parameter('RoI_mode', 'corner', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                                   corner: pta and ptb are diagonal corners of rectangle.\n \
                                                   center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 5, int, 'number of x points to scan'),
                   Parameter('y', 5, int, 'number of y points to scan')
                   ]),
        Parameter('settle_time', 0.5, float, 'time (s) to wait after moving to a new point, before taking data'),
        Parameter('track_esr', False, bool, 'center ESR scan around the frequency measured previously at a neighboring grid point'),
        Parameter('measure_gradient', [
            Parameter('enable', False, bool, 'measure field gradient by taking ESR at two Z values at each lateral grid point'),
            Parameter('z1', 0, float, 'take first ESR at this Attocube Z voltage (V)'),
            Parameter('z2', 0, float, 'take second ESR at this Attocube Z voltage (V)'),
        ]),
        Parameter('Tracking', [
            Parameter('on/off', True, bool, 'used to turn on tracking'),
            Parameter('every_N', 1, int, 'track every n points')]),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn')
    ]

    _SCRIPTS = {'find_nv': FindNvStrobe, 'PulsedEsrFast': PulsedEsrFast, 'SetAttoANC300': SetAtto}

    def move_piezo(self, x, y, direction='both'):
        self.scripts['SetAttoANC300'].settings['point'].update({'x': float(x), 'y': float(y)})
        self.scripts['SetAttoANC300'].run(verbose=False)

    def read_point(self, yNum, xNum):
        """
        Take ESR (for a 2D scan of B field)
        Returns: measured value at given pt, with other data that won't be plotted
        """

        self.scripts['PulsedEsrFast'].settings['tag'] = 'pulsedesrfast_%i,%i' % (yNum, xNum)

        if self.settings['track_esr']:
            if yNum == 0 and xNum == 0:
                # For origin (where scan begins), simply use ESR subscript settings to set center freq
                pass
            elif xNum > 0:
                # Center current ESR scan range around extracted frequency of the left adjacent (previous) point
                predicted_freq = self.data['point_freq'][yNum][xNum-1]
                if predicted_freq != 0:
                    print('Predicted freq: %.3e using pt: %i, %i' % (predicted_freq, yNum, xNum-1))
                    self.scripts['PulsedEsrFast'].settings['range_type'] = 'center_range'
                    self.scripts['PulsedEsrFast'].settings['freq_start'] = float(predicted_freq)
            elif xNum == 0:
                # Center current ESR scan range around extracted frequency of the upper adjacent (previous) point
                predicted_freq = self.data['point_freq'][yNum - 1][xNum]
                if predicted_freq != 0:
                    print('Predicted freq: %.3e using pt: %i, %i' % (predicted_freq, yNum-1, xNum))
                    self.scripts['PulsedEsrFast'].settings['range_type'] = 'center_range'
                    self.scripts['PulsedEsrFast'].settings['freq_start'] = float(predicted_freq)

        self.scripts['PulsedEsrFast'].run(verbose=False)

        try:
            point_data = self.scripts['PulsedEsrFast'].data['fits']
        except:
            point_data = None

        if point_data is not None:

            if len(point_data) == 6:
                fp = point_data[5]
                fn = point_data[4]
                self.data['point_freq'][yNum][xNum] = point_data[4]
                point_value = B_field_from_esr(fp, fn, D=self.settings['zero_field_splitting'],
                                               gamma=27.969e9, angular_freq=False, verbose=False)[0] * 1e4
            elif len(point_data) == 4:
                self.data['point_freq'][yNum][xNum] = point_data[3]
                point_value = (point_data[3] - self.settings['zero_field_splitting']) / nv_gyro
        else:
            self.data['point_freq'][yNum][xNum] = 0
            point_value = 0

        return point_value

    def is_valid(self):
        return self.scripts['PulsedEsrFast'].is_valid()

    def _plot(self, axes_list, data=None):

        if data is None:
            data = self.data
        label = ['B field scan w/ Attocube', r'V$_x$ [V]', r'V$_y$ [V]', 'On-axis field (G)']

        if self.linescan_axis == '2d':
            plot_fluorescence_new(np.abs(data['point_value']), data['extent'], axes_list[0], max_counts=-1, labels=label)
        elif self.linescan_axis == 'x' or self.linescan_axis == 'y':
            raise NotImplementedError

        axes_list[1].clear()

        current_esr_data = None
        if self._current_subscript_stage['current_subscript'] is self.scripts['PulsedEsrFast'] and 'count_data' in self.scripts['PulsedEsrFast'].data:
            current_esr_data = self.scripts['PulsedEsrFast'].data['count_data']
            esr_freqs = self.scripts['PulsedEsrFast'].data['params']

        if current_esr_data is not None:
            plot_1d_simple_freq(axes_list[1], esr_freqs, [current_esr_data])

    def _update_plot(self, axes_list):
        update_fluorescence(np.abs(self.data['point_value']), axes_list[0], -1)

        axes_list[1].clear()

        current_esr_data = None
        if self._current_subscript_stage['current_subscript'] is self.scripts['PulsedEsrFast'] and 'count_data' in self.scripts['PulsedEsrFast'].data:
            current_esr_data = self.scripts['PulsedEsrFast'].data['count_data']
            esr_freqs = self.scripts['PulsedEsrFast'].data['params']

        if current_esr_data is not None:
            plot_1d_simple_freq(axes_list[1], esr_freqs, [current_esr_data])


class AttoGridScanGalvoScan(AttoGridScanEsr):
    """
    Moves Attocubes in a grid and takes a galvo scan at each point. Can be used to calibrate the displacement of the Attocube.
    """
    _DEFAULT_SETTINGS = [
        Parameter('zero_field_splitting', 2.8707e9, float, 'Dip in ESR with no external B field applied'),
        Parameter('point_a',
                  [Parameter('x', 0, float, 'x-coordinate'),
                   Parameter('y', 0, float, 'y-coordinate')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 1.0, float, 'x-coordinate'),
                   Parameter('y', 1.0, float, 'y-coordinate')
                   ]),
        Parameter('RoI_mode', 'corner', ['corner'], 'mode to calculate region of interest.\n \
                                                   corner: pta and ptb are diagonal corners of rectangle.\n \
                                                   center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 5, int, 'number of x points to scan'),
                   Parameter('y', 5, int, 'number of y points to scan')
                   ]),
        Parameter('settle_time', 0.5, float, 'time (s) to wait after moving to a new point, before taking data'),
        Parameter('Tracking', [
            Parameter('on/off', True, bool, 'used to turn on tracking'),
            Parameter('every_N', 1, int, 'track every n points')]),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn')
    ]

    _SCRIPTS = {'find_nv': FindNv, 'take_image': GalvoScan, 'SetAttoANC300': SetAtto}

    def move_piezo(self, x, y, direction='both'):
        self.scripts['SetAttoANC300'].settings['point'].update({'x': float(x), 'y': float(y)})
        self.scripts['SetAttoANC300'].run()
        self.scripts['take_image'].settings['tag'] = 'take_image_%.1f,%.1f' % (x, y)

    def read_point(self, yNum, xNum):
        """
        Take ESR (for a 2D scan of B field)
        Returns: measured value at given pt, with other data that won't be plotted
        """

        self.scripts['take_image'].run()
        point_value = 0
        return point_value

    def _plot(self, axes_list, data = None):
        """
        Plots the galvo scan image
        Args:
            axes_list: list of axes objects on which to plot the galvo scan on the first axes object
            data: data (dictionary that contains keys image_data, extent) if not provided use self.data
        """

        pass


    def _update_plot(self, axes_list):
        """
        updates the galvo scan image
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
        """
        pass


class AttoGridScanEsrPiezoController(AttoGridScanGeneric):
    """
    Moves Attocubes in a grid and takes EsrSimple (CW) at each point, resulting in a magnetic field scan. This script moves the Attocubes by communicating
    digitally with a PiezoController connected to the Attocubes, instead of using DAQ outputs to control the DC-inputs of an Attocube/piezo controller
    """
    _SCRIPTS = {'find_nv': FindNv, 'ESR_simple': EsrSimple}
    _INSTRUMENTS = {'piezo_controller': PiezoController}

    def move_piezo(self, x, y, direction='both'):
        piezo = self.instruments['piezo_controller']['instance']

        piezo.axis = 'x'
        piezo.voltage = x

        piezo.axis = 'y'
        piezo.voltage = y

    def read_point(self, yNum, xNum):
        """
        Take ESR (for a 2D scan of B field)
        Returns: measured value at given pt, with other data that won't be plotted
        """

        self.scripts['ESR_simple'].run()
        point_data = self.scripts['ESR_simple'].data['fit_params']
        print(point_data)
        if len(point_data) == 5:
            fp = point_data[5]
            fn = point_data[4]
            point_value = B_field_from_esr(fp, fn, D=self.settings['zero_field_splitting'],
                                           gamma=27.969e9, angular_freq=False, verbose=False)[0] * 1e4
        else:
            point_value = 0

        return point_value


class GalvoGridScanPulsed(AttoGridScanGeneric):
    """
    Test script. Very slow because it has to program the PB at each point. Use GalvoScanStrobe instead.
    Runs set laser at each point in a defined grid, and get the counts by running StroboscopicReadout.
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
                  [Parameter('x', 5, int, 'number of x points to scan, if 1 then perform a line scan along the other axis'),
                   Parameter('y', 5, int, 'number of y points to scan, if 1 then perform a line scan along the other axis')
                   ]),
        Parameter('settle_time', 0.01, float, 'time (s) to wait after moving to a new point, before taking data'),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn')
    ]

    _SCRIPTS = {'set_laser': SetLaser, 'set_laser_single_axis': SetLaserSingleAxis, 'stroboscopic_readout': PulsedReadout}

    def read_point(self, yNum, xNum):
        """
        Replace this with some data-taking script, e.g. read the counts (to perform a galvoscan with the attocubes
        instead) or take ESR (for a 2D scan of B field)
        Returns: measured value at given pt, with other data that won't be plotted (e.g. for ESRs, we want to plot
        only the on-axis field, and save but not plot the off-axis field)
        """
        self.scripts['stroboscopic_readout'].run(verbose=False)
        print(self.scripts['stroboscopic_readout'].data['counts'])
        return self.scripts['stroboscopic_readout'].data['counts'][0]

    def move_piezo(self, x, y, direction='both'):
        if direction == 'x':
            self.scripts['set_laser_single_axis'].settings['axis'] = 'x'
            self.scripts['set_laser_single_axis'].settings['point'] = x
            self.scripts['set_laser_single_axis'].run(verbose=False)
        elif direction == 'y':
            self.scripts['set_laser_single_axis'].settings['axis'] = 'y'
            self.scripts['set_laser_single_axis'].settings['point'] = y
            self.scripts['set_laser_single_axis'].run(verbose=False)
        elif direction == 'both':
            self.scripts['set_laser'].settings['point']['x'] = x
            self.scripts['set_laser'].settings['point']['y'] = y
            self.scripts['set_laser'].run(verbose=False)

    def check_bounds(self):
        pass


class ServoGridScanFringes(AttoGridScanGeneric):
    _DEFAULT_SETTINGS = [
        Parameter('point_a',
                  [Parameter('x', 0, float, 'x-coordinate'),
                   Parameter('y', 0, float, 'y-coordinate')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 1.0, float, 'x-coordinate'),
                   Parameter('y', 1.0, float, 'y-coordinate')
                   ]),
        Parameter('RoI_mode', 'corner', ['center', 'corner'], 'mode to calculate region of interest.\n \
                                                   corner: pta and ptb are diagonal corners of rectangle.\n \
                                                   center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 10, int, 'number of x points to scan, if 1 then perform a line scan along the other axis'),
                   Parameter('y', 10, int, 'number of y points to scan, if 1 then perform a line scan along the other axis')
                   ]),
        Parameter('scan_mode', 'meander', ['meander', 'book'], 'Meander: scan from left to right, and then back to left for the next line, etc; '
                                                               'Book: always scan from left to right'),
        Parameter('statistics', 'min-max', ['min-max', 'mean', 'std']),
        Parameter('settle_time', 0.05, float, 'time (s) to wait after moving to a new point, before taking data'),
        Parameter('Tracking', [
            Parameter('on/off', False, bool, 'used to turn on tracking'),
            Parameter('every_N', 1, int, 'track every n points')]),
        Parameter('ending_behavior', 'leave_at_corner', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn')
    ]

    _INSTRUMENTS = {'XServo': B26KDC001x, 'YServo': B26KDC001y}
    _SCRIPTS = {'daq_read_ai': Daq_Read_Analog}

    def check_bounds(self):
        pass

    def setup_scan(self):
        servo_x = self.instruments['XServo']['instance']
        servo_y = self.instruments['YServo']['instance']
        servo_x.settings['velocity'] = 4
        print(servo_x.get_velocity())
        servo_x.set_velocity()
        servo_y.settings['velocity'] = 2.2
        servo_y.set_velocity()

    def move_piezo(self, x, y, direction='both'):
        servo_x = self.instruments['XServo']['instance']
        servo_y = self.instruments['YServo']['instance']
        if direction == 'x':
            servo_x.settings['position'] = x
            t_start = time.time()
            servo_x.set_position()
            print('Time to move piezo in x: %.3f'%(time.time()-t_start))
        elif direction == 'y':
            servo_y.settings['position'] = y
            servo_y.set_position()
        elif direction == 'both':
            servo_x.settings['position'] = x  # update the position setting of the instrument
            servo_y.settings['position'] = y
            servo_x.set_position()  # actually move the instrument to that location. If this is not within the safety
            servo_y.set_position()  # limits of the instruments, it will not actually move and say so in the log
        else:
            raise ValueError

    def read_point(self, yNum, xNum):
        """
        Take ESR (for a 2D scan of B field)
        Returns: measured value at given pt, with other data that won't be plotted
        """

        self.scripts['daq_read_ai'].run(verbose=False)
        time.sleep(self.settings['settle_time'])

        data = self.scripts['daq_read_ai'].data['voltage']
        self.daq_data = data
        if data:  # Need this conditional because data might be empty list if abort in the middle of Daq Read Analog
            if self.settings['statistics'] == 'mean':
                value_from_data = np.mean(data)
            elif self.settings['statistics'] == 'std':
                value_from_data = np.std(data)
            elif self.settings['statistics'] == 'min-max':
                value_from_data = np.max(data) - np.min(data)

        return value_from_data

    def read_point_disabled(self):
        """
        Take ESR (for a 2D scan of B field)
        Returns: measured value at given pt, with other data that won't be plotted
        """

        return 0


if __name__ == '__main__':
    script, failed, instr = Script.load_and_append({'AttoStep': 'AttoStep'})

    print(script)
    print(failed)
    print(instr)
    # fp = Find_Points(settings={'path': 'Z:/Lab/Cantilever/Measurements/__tmp__', 'tag':'nvs'})
    # fp.run()

    # plt.pcolor(fp.data['image'])
    # print(fp.data['image_gaussian'].shape)
    # plt.pcolor(fp.data['image'])
    # plt.imshow(fp.data['image'], cmap = 'pink', interpolation = 'nearest')
    #
    #
    # for x in fp.data['NV_positions']:
    #     plt.plot(x[0],x[1],'ro')
    #
    # plt.show()

    # plt.figure()
    # plt.imshow(fp.data['image_gaussian'])
    # Axes3D.plot(fp.data['image_gaussian'])
    # plt.show()
    # print(max(fp.data['image']))
    # print(max(fp.data['image_gaussian'].flatten()))
    # print('NV_positions', fp.data['NV_positions'])
