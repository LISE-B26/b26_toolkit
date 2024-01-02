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
from b26_toolkit.instruments import PiezoController
from b26_toolkit.scripts.galvo_scan.galvo_scan import GalvoScan
from b26_toolkit.scripts.set_laser import SetLaser, SetLaserSingleAxis
from b26_toolkit.scripts.pulse_sequences.laser_pulses import PulsedReadout
from b26_toolkit.scripts.pulse_sequences.param_sweep.pulsed_esr import PulsedEsrFast
from pylabcontrol.core import Parameter, Script
from b26_toolkit.scripts import FindNv, EsrSimple
from b26_toolkit.scripts.find_nv import FindNvStrobe
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from b26_toolkit.plotting.plots_1d import plot_1d_simple_freq
from b26_toolkit.scripts.set_laser import SetAtto
from b26_toolkit.data_analysis.nv_optical_response import B_field_from_esr


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
        Parameter('RoI_mode', 'corner', ['corner'], 'mode to calculate region of interest.\n \
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

        self.data = {'point_value': np.zeros((self.settings['num_points']['y'], self.settings['num_points']['x'])),
                     'extent': self.pts_to_extent(self.settings['point_a'], self.settings['point_b'],
                                                 self.settings['RoI_mode'])}

        [xVmin, xVmax, yVmax, yVmin] = self.data['extent']
        if self.linescan_axis == 'x' and self.settings['RoI_mode'] == 'corner' and self.settings['point_a']['x'] > self.settings['point_b']['x']:
            xVmin, xVmax = self.settings['point_a']['x'],  self.settings['point_b']['x']
        elif self.linescan_axis == 'y' and self.settings['RoI_mode'] == 'corner' and self.settings['point_a']['y'] > self.settings['point_b']['y']:
            yVmin, yVmax = self.settings['point_a']['y'], self.settings['point_b']['y']

        print([xVmin, xVmax, yVmax, yVmin])
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
        for yNum in range(0, Ny):
            self.move_piezo(0, self.y_array[yNum], direction='y')
            for xNum in range(0, Nx):
                if self._abort:
                    break
                print('Moving piezo to %.1f, %.1f' % ((self.x_array[xNum], self.y_array[yNum])))
                self.move_piezo(self.x_array[xNum], self.y_array[yNum], direction='x')
                time.sleep(self.settings['settle_time'])
                #self.data['actual_Vx'][yNum*Nx+xNum], self.data['actual_Vy'][yNum*Nx+xNum] = self.read_piezo()

                self.progress = float(yNum * Nx + xNum) / (Nx * Ny) * 100
                self.updateProgress.emit(int(self.progress))
                self.point_index += 1
                if 'tracking' in self.settings:
                    if self.settings['Tracking']['on/off'] and int(self.point_index) % self.settings['Tracking']['every_N'] == 0:
                        self.scripts['find_nv'].run(verbose=False)

                        # Move piezo to the set location again. For some reason, FindNv seems to affect the DAQ outputs to the Attocube
                        print('Moving piezo (again) to %.1f, %.1f' % ((self.x_array[xNum], self.y_array[yNum])))
                        self.move_piezo(self.x_array[xNum], self.y_array[yNum])
                        time.sleep(.5)

                point_value = self.read_point()
                self.data['point_value'][yNum, xNum] = point_value

                self.progress = float(yNum * Nx + 1 + xNum) / (Nx * Ny) * 100

                #print('Current acquisition {:02d}/{:02d}'.format(yNum * Nx + xNum, Nx * Ny))

                self.updateProgress.emit(int(self.progress))

        # set end position after scan based on ending_behavior setting
        if self.settings['ending_behavior'] == 'leave_at_corner':
            return
        elif self.settings['ending_behavior'] == 'return_to_start':
            self.move_piezo(self.settings['point_a']['x'], self.settings['point_a']['y'])
        elif self.settings['ending_behavior'] == 'return_to_origin':
            self.move_piezo(0, 0)

    def check_bounds(self):
        # Check if there are negative voltages in the scan range
        negative_voltage_err = 'Piezo voltage cannot be < 0 V'
        if np.any(self.x_array < 0):
            raise ValueError(negative_voltage_err)
        elif np.any(self.y_array < 0):
            raise ValueError(negative_voltage_err)

    def setup_scan(self):
        pass

    def read_point(self):
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

    _SCRIPTS = {'find_nv': FindNv, 'ESR_simple': EsrSimple, 'SetAttoANC300': SetAtto}

    def move_piezo(self, x, y, direction='both'):
        self.scripts['SetAttoANC300'].settings['point'].update({'x': float(x), 'y': float(y)})
        self.scripts['SetAttoANC300'].run()
        time.sleep(.1)  # Wait for Find NV to settle

        self.scripts['ESR_simple'].settings['tag'] = 'esr_simple_%.1f,%.1f' % (x, y)

    def read_point(self):
        """
        Take ESR (for a 2D scan of B field)
        Returns: measured value at given pt, with other data that won't be plotted
        """

        self.scripts['ESR_simple'].run()
        point_data = self.scripts['ESR_simple'].data['fit_params']
        print(point_data)
        if point_data is not None and len(point_data) == 6:
            fp = point_data[5]
            fn = point_data[4]
            point_value = B_field_from_esr(fp, fn, D=self.settings['zero_field_splitting'],
                                           gamma=27.969e9, angular_freq=False, verbose=False)[0] * 1e4
            print(point_value)
            print()
        else:
            point_value = 0

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

    _SCRIPTS = {'find_nv': FindNvStrobe, 'PulsedEsrFaster': PulsedEsrFast, 'SetAttoANC300': SetAtto}

    def move_piezo(self, x, y, direction='both'):
        self.scripts['SetAttoANC300'].settings['point'].update({'x': float(x), 'y': float(y)})
        self.scripts['SetAttoANC300'].run()
        time.sleep(self.settings['settle_time'])  # Wait for Find NV to settle

        self.scripts['PulsedEsrFaster'].settings['tag'] = 'pulsedesrfaster_%.1f,%.1f' % (x, y)

    def read_point(self):
        """
        Take ESR (for a 2D scan of B field)
        Returns: measured value at given pt, with other data that won't be plotted
        """

        self.scripts['PulsedEsrFaster'].run()

        try:
            point_data = self.scripts['PulsedEsrFaster'].data['fits']
        except:
            point_data = None

        if point_data is not None and len(point_data) == 6:
            fp = point_data[5]
            fn = point_data[4]
            point_value = B_field_from_esr(fp, fn, D=self.settings['zero_field_splitting'],
                                           gamma=27.969e9, angular_freq=False, verbose=False)[0] * 1e4
        else:
            point_value = 0

        return point_value

    def is_valid(self):
        return self.scripts['PulsedEsrFaster'].is_valid()

    def _plot(self, axes_list, data=None):

        if data is None:
            data = self.data
        label = ['B field scan w/ Attocube', r'V$_x$ [V]', r'V$_y$ [V]', 'On-axis field (G)']

        if self.linescan_axis == '2d':
            plot_fluorescence_new(data['point_value'], data['extent'], axes_list[0], max_counts=-1, labels=label)
        elif self.linescan_axis == 'x' or self.linescan_axis == 'y':
            raise NotImplementedError

        axes_list[1].clear()

        current_esr_data = None
        if self._current_subscript_stage['current_subscript'] is self.scripts['PulsedEsrFaster'] and 'count_data' in self.scripts['PulsedEsrFaster'].data:
            current_esr_data = self.scripts['PulsedEsrFaster'].data['count_data']
            esr_freqs = self.scripts['PulsedEsrFaster'].data['params']

        if current_esr_data is not None:
            plot_1d_simple_freq(axes_list[1], esr_freqs, [current_esr_data])

    def _update_plot(self, axes_list):
        update_fluorescence(self.data['point_value'], axes_list[0], -1)

        axes_list[1].clear()

        current_esr_data = None
        if self._current_subscript_stage['current_subscript'] is self.scripts['PulsedEsrFaster'] and 'count_data' in self.scripts['PulsedEsrFaster'].data:
            current_esr_data = self.scripts['PulsedEsrFaster'].data['count_data']
            esr_freqs = self.scripts['PulsedEsrFaster'].data['params']

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
        time.sleep(.1)  # Wait for Find NV to settle

        self.scripts['take_image'].settings['tag'] = 'take_image_%.1f,%.1f' % (x, y)

    def read_point(self):
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

    def read_point(self):
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

    def read_point(self):
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
        time.sleep(self.settings['settle_time'])

    def check_bounds(self):
        pass


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
