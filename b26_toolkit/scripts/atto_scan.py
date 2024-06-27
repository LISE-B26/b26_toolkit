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
from b26_toolkit.instruments import ANC300, ANC350, PiezoController, NI6259, NI9263, NI9402
from b26_toolkit.scripts.daq_read_counter import DaqReadCounter, Daq_Read_Counter_NI6259
from b26_toolkit.scripts.galvo_scan.galvo_scan import GalvoScan
from b26_toolkit.plotting.plots_1d import plot_counts_vs_pos, update_counts_vs_pos
from b26_toolkit.plotting.plots_2d import plot_fluorescence_pos, update_fluorescence
from pylabcontrol.core import Parameter, Script
from b26_toolkit.scripts import FindNV, ESR_simple, SetAtto
from b26_toolkit.scripts.autofocus import AutoFocusDAQ
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from b26_toolkit.scripts.set_laser import SetAttoANC300
from b26_toolkit.data_analysis.nv_optical_response import B_field_from_esr



class AttoStep(Script):
    """Steps Attocube in number of steps specificed for each axis. Note that running this script will disable any offset afterwards!"""
    _DEFAULT_SETTINGS = [
        Parameter('controller_type', 'ANC350', ['ANC300', 'ANC350'], 'attocube controller model'),
        Parameter('num_steps', [
            Parameter('x', 0, int, 'num steps along x-axis'),
            Parameter('y', 0, int, 'num steps along y-axis'),
            Parameter('z', 0, int, 'num steps along z-axis')
        ]
    )]

    _INSTRUMENTS = {'ANC300': ANC300, 'ANC350': ANC350}
    _SCRIPTS = {}

    def _function(self):
        """
        Performs a multiple attocube step with the voltage and frequency specified in instrument,
        and the direction and number specified in settings
        """
        attocube = self.instruments[self.settings['controller_type']]['instance']
        for axis in self.settings['num_steps']:
            attocube.multistep(axis, self.settings['num_steps'][axis])


class AttoSetAndStep(AttoStep):
    """Sets DC input to zero first (with 16 Hz filter on to prevent any slippage), then steps Attocube in number of steps specificed for each axis.
    Note that running this script will disable any offset afterwards!"""
    _DEFAULT_SETTINGS = [
        Parameter('controller_type', 'ANC300', ['ANC300', 'ANC350'], 'attocube controller model'),
        Parameter('num_steps', [
            Parameter('x', 0, int, 'num steps along x-axis'),
            Parameter('y', 0, int, 'num steps along y-axis'),
            Parameter('z', 0, int, 'num steps along z-axis')
        ]
    )]

    _INSTRUMENTS = {'ANC300': ANC300, 'ANC350': ANC350}
    _SCRIPTS = {'SetAttoANC300': SetAttoANC300}

    def _function(self):
        """
        Performs a multiple attocube step with the voltage and frequency specified in instrument,
        and the direction and number specified in settings
        """
        attocube = self.instruments[self.settings['controller_type']]['instance']

        self.scripts['SetAttoANC300'].run()

        for axis in self.settings['num_steps']:
            attocube.multistep(axis, self.settings['num_steps'][axis])


class AttoStepXY(AttoStep):
    # Constrains Attostep so that only X & Y moves
    _DEFAULT_SETTINGS = [
        Parameter('controller_type', 'ANC350', ['ANC300', 'ANC350'], 'attocube controller model'),
        Parameter('num_steps', [
            Parameter('x', 0, int, 'num steps along x-axis'),
            Parameter('y', 0, int, 'num steps along y-axis'),
        ]
    )]


class AttoScanOpenLoop(Script):
    _DEFAULT_SETTINGS = [
        Parameter('controller_type', 'ANC350', ['ANC300', 'ANC350'], 'attocube controller model'),
        Parameter('scan_axis', 'y', ['x', 'y'], 'axis to scan on'),
        Parameter('num_steps', 100, int, 'number of points in the scan'),
        Parameter('num_points', 10, int, 'number of DAQ counter measurements'),
    ]

    _SCRIPTS = {'daq_read_counter': Daq_Read_Counter}
    _INSTRUMENTS = {'ANC300': ANC300, 'ANC350': ANC350}

    def _function(self):
        self.attocube = self.instruments[self.settings['controller_type']]['instance']

        scan_axis = self.settings['scan_axis']
        steps_per_meas = self.settings['num_steps'] // self.settings['num_points']

        self.data = {
            'counts': [],
            'positions': np.arange(0, self.settings['num_steps'], steps_per_meas)
        }

        for i in range(len(self.data['positions'])):
            if self._abort:
                break

            # run daq_read_counter or the relevant script to get fluorescence
            self.scripts['daq_read_counter'].run()

            # add to output structures which will be plotted
            data = self.scripts['daq_read_counter'].data['counts']
            self.data['counts'].append(np.mean(data))

            self.attocube.multistep(scan_axis, steps_per_meas)

            self.progress = i * 100. / self.settings['num_points']
            self.updateProgress.emit(int(self.progress))

        # clean up data, as in daq_read_counter
        self.data['counts'] = list(self.data['counts'])
        self.data['positions'] = list(self.data['positions'])

    def _plot(self, axes_list, data=None):
        if data is None:
            data = self.data

        if data:
            axes_list[0].set_xlabel('steps')
            axes_list[0].set_ylabel('kCounts/sec')

            axes_list[0].plot(data['positions'][:len(data['counts'])], data['counts'], linewidth=2.0)

    def _update_plot(self, axes_list):
        update_counts_vs_pos(axes_list[0], self.data['counts'], self.data['positions'][:len(self.data['counts'])])


class AttoScanOpenLoopTrackNI6259(Script):

    ''''

    ER 20200911

    Same as AttoScanOpenLoop,
    but adds findNV and AutofocusDAQ between the attoscan and DAQ read counter scripts,
    in order to make sure we are on the NV

    Takes the daq read counter script corresponding to NI6259 DAQ (Alice)

    '''
    _DEFAULT_SETTINGS = [
        Parameter('controller_type', 'ANC350', ['ANC300', 'ANC350'], 'attocube controller model'),
        Parameter('scan_axis', 'y', ['x', 'y'], 'axis to scan on'),
        Parameter('num_steps', 100, int, 'number of points in the scan'),
        Parameter('num_points', 10, int, 'number of DAQ counter measurements'),
    ]

    _SCRIPTS = {'daq_read_counter': Daq_Read_Counter_NI6259, 'find_nv': FindNV, 'autofocus': AutoFocusDAQ}
    _INSTRUMENTS = {'attocube': ANC300}

    def _function(self):
        self.attocube = self.instruments[self.settings['controller_type']]['instance']

        scan_axis = self.settings['scan_axis']
        steps_per_meas = self.settings['num_steps'] // self.settings['num_points']

        self.data = {
            'counts': [],
            'positions': np.arange(0, self.settings['num_steps'], steps_per_meas)
        }

        for i in range(len(self.data['positions'])):
            if self._abort:
                break

            # track to the NV
            # ER 20200911
            self.scripts['autofocus'].run()
            self.scripts['find_nv'].run()

            # run daq_read_counter or the relevant script to get fluorescence
            self.scripts['daq_read_counter'].run()

            # add to output structures which will be plotted
            data = self.scripts['daq_read_counter'].data['counts']
            self.data['counts'].append(np.mean(data))

            self.attocube.multistep(scan_axis, steps_per_meas)

            self.progress = i * 100. / self.settings['num_points']
            self.updateProgress.emit(int(self.progress))

        # clean up data, as in daq_read_counter
        self.data['counts'] = list(self.data['counts'])
        self.data['positions'] = list(self.data['positions'])

    def _plot(self, axes_list, data=None):
        if data is None:
            data = self.data

        if data:
            axes_list[0].set_xlabel('steps')
            axes_list[0].set_ylabel('kCounts/sec')

            axes_list[0].plot(data['positions'][:len(data['counts'])], data['counts'], linewidth=2.0)

    def _update_plot(self, axes_list):
        update_counts_vs_pos(axes_list[0], self.data['counts'], self.data['positions'][:len(self.data['counts'])])


class AttoScanOpenLoop2D(Script):
    _DEFAULT_SETTINGS = [
        Parameter('controller_type', 'ANC350', ['ANC300', 'ANC350'], 'attocube controller model'),
        Parameter('outer_loop', [
            Parameter('scan_axis', 'x', ['x', 'y'], 'axis to scan on'),
            Parameter('num_steps', 100, int, 'number of points in the scan'),
            Parameter('num_points', 10, int, 'number of DAQ counter measurements')
        ]
    )]

    _SCRIPTS = {'inner_scan': AttoScanOpenLoop}
    _INSTRUMENTS = {}

    def _function(self):
        self.inner_scan = self.scripts['inner_scan']
        self.attocube = self.instruments[self.settings['controller_type']]['instance']

        outer_settings = self.settings['outer_loop']
        inner_settings = self.inner_scan.settings

        if outer_settings['scan_axis'] == inner_settings['scan_axis']:
            raise ValueError('scan axes cannot be the same')

        steps_per_outermeas = outer_settings['num_steps'] // outer_settings['num_points']
        outer_pos = np.arange(0, outer_settings['num_steps'], steps_per_outermeas)

        self.data = {
            'counts': np.zeros((outer_settings['num_points'], inner_settings['num_points'])),
            'positions': {
                'inner': np.arange(0, inner_settings['num_steps'],
                                   inner_settings['num_steps'] // inner_settings['num_points']),
                'outer': outer_pos
            }
        }

        for i_out in range(len(outer_pos)):

            if self._abort:
                break

            self.inner_scan.run()
            if i_out % 2 == 0:
                self.data['counts'][i_out, :] = self.inner_scan.data['counts']
            else:
                self.data['counts'][i_out, :] = np.flip(self.inner_scan.data['counts'], 0)

            self.inner_scan.settings['num_steps'] = -self.inner_scan.settings['num_steps']

            self.attocube.multistep(outer_settings['scan_axis'], steps_per_outermeas)

            self.progress = i_out * 100. / self.settings['outer_loop']['num_points']
            self.updateProgress.emit(int(self.progress))

        self.data['counts'] = list(self.data['counts'])

    def _plot(self, axes_list, data=None):
        # COMMENT_ME

        if data is None:
            data = self.data

        extent = [0, self.data['positions']['inner'][-1],
                  0, self.data['positions']['outer'][-1]]

        if data:
            plot_fluorescence_pos(data['counts'], extent, axes_list[0])
            axes_list[0].set_xlabel(r'pos$_{outer}$ [steps]')
            axes_list[0].set_ylabel(r'pos$_{inner}$ [steps]')

    def _update_plot(self, axes_list):
        update_fluorescence(self.data['counts'], axes_list[0])


class AttoScanClosedLoop(Script):

    # 20211012 FF: Not finished
    _DEFAULT_SETTINGS = [
        Parameter('scan_axis', 'x', ['x', 'y'], 'axis to scan on'),
        Parameter('direction', 'positive', ['positive', 'negative'], 'direction to scan'),
        Parameter('num_points', 100, int, 'number of points in the scan'),
        Parameter('min_pos', 0., float, 'minimum position of scan (um)'),
        Parameter('max_pos', 5., float, 'maximum position of scan (um)'),
        Parameter('time_per_pt', 0.5, float, 'time to wait at each point (s)'),
    ]

    _SCRIPTS = {'daq_read_counter': Daq_Read_Counter, "attocube": ANC300}

    def __init__(self, name=None, settings=None, instruments=None, scripts=None, log_function=None, data_path=None):
        """
        Default script initialization
        """
        self.attocube = self.instruments['attocube']['instance']
        super().__init__(name=None, settings=None, instruments=None, scripts=None, log_function=None, data_path=None)

    def _get_scan_positions(self, params):
        return np.linspace(params['min_pos'], params['max_pos'], params['num_points'])

        Script.__init__(self, name, settings = settings, instruments = instruments, log_function= log_function, data_path = data_path)
        
    def _function(self):
        scan_axis = self.settings['scan_axis']

        self.positions = self._get_scan_positions(self.settings)

        self.data = {
            'counts': np.zeros(self.settings['num_points']),
            'position': self.positions
        }

        indices = range(self.settings['num_points'])
        if self.settings['direction'] == 'negative':
            indices = reversed(indices)

        for i in indices:
            if self._abort:
                break

            self.attocube.move_absolute(scan_axis, self.positions[i])

            # run daq_read_counter or the relevant script to get fluorescence
            self.scripts['daq_read_counter'].run()
            time.sleep(self.settings['time_per_pt'])
            self.scripts['daq_read_counter'].stop()

            # add to output structures which will be plotted
            data = self.scripts['daq_read_counter'].data['counts']
            self.data['counts'][i] = np.mean(data)

            self.progress = i * 100. / self.settings['num_points']
            self.updateProgress.emit(int(self.progress))

        # clean up data, as in daq_read_counter
        self.data['counts'] = list(self.data['counts'])

    def _plot(self, axes_list, data=None):
        if data is None:
            data = self.data

        if data:
            plot_counts_vs_pos(axes_list[0], data['counts'], self.positions)

    def _update_plot(self, axes_list):
        update_counts_vs_pos(axes_list[0], self.data['counts'], self.positions)


class AttoScanClosedLoop2D(Script):
    _DEFAULT_SETTINGS = [
        Parameter('outer_loop',
                  [
                      Parameter('scan_axis', 'x', ['x', 'y'], 'axis to scan on'),
                      Parameter('direction', 'positive', ['positive', 'negative'], 'direction to scan'),
                      Parameter('num_points', 100, int, 'number of points in the scan'),
                      Parameter('min_pos', 0., float, 'minimum position of scan (um)'),
                      Parameter('max_pos', 5., float, 'maximum position of scan (um)'),
                      Parameter('time_per_pt', 0.5, float, 'time to wait at each point (s)'),
                  ]),
        Parameter('inner_loop',
                  [
                      Parameter('scan_axis', 'x', ['x', 'y'], 'axis to scan on'),
                      Parameter('direction', 'positive', ['positive', 'negative'], 'direction to scan'),
                      Parameter('num_points', 100, int, 'number of points in the scan'),
                      Parameter('min_pos', 0., float, 'minimum position of scan (um)'),
                      Parameter('max_pos', 5., float, 'maximum position of scan (um)'),
                      Parameter('time_per_pt', 0.5, float, 'time to wait at each point (s)'),
                  ]),
        Parameter('time_per_pt', 0.5, float, 'time to wait at each point (s)'),
    ]

    _SCRIPTS = {'attoscan': AttoScanClosedLoop}

    def __init__(self, name=None, settings=None, instruments=None, scripts=None, log_function=None, data_path=None):
        """
        Default script initialization
        """
        self.attocube = self.scripts['attoscan']['instruments']['instance']

        self.inner_scan = self.scripts['attoscan']
        self.inner_scan.settings = self.settings['inner_loop']

        super().__init__(name=None, settings=None, instruments=None, scripts=None, log_function=None, data_path=None)

    def _switch_direction(self, direction):
        if direction == 'positive':
            return 'negative'
        return 'positive'

    def _function(self):
        outer_settings = self.settings['outer_loop']
        inner_settings = self.settings['inner_loop']

        if outer_settings['scan_axis'] == inner_settings['scan_axis']:
            raise ValueError('scan axes cannot be the same')

        self.outer_pos = self._get_scan_positions(outer_settings),

        self.data['counts'] = np.zeros((outer_settings['num_points'], inner_settings['num_points']))

        out_range = range(outer_settings['num_points'])
        if outer_settings['direction'] == 'negative':
            out_range = reversed(out_range)

        for i_out in out_range:
            self.attocube.move_absolute(outer_settings['scan_axis'], self.positions['outer_loop'][i_out])

            self.inner_scan.run()

            self.data['counts'][i_out, :] = self.inner_scan.data['counts']

            self.inner_scan['direction'] = self._switch_direction(self.inner_scan['direction'])

            self.progress = i_out * 100. / self.settings['outer_loop']['num_points']
            self.updateProgress.emit(int(self.progress))

        self.data['counts'] = list(self.data['counts'])

    def _plot(self, axes_list, data=None):
        # COMMENT_ME

        extent = [self.settings['inner_loop']['min_pos'], self.settings['inner_loop']['max_pos'],
                  self.settings['outer_loop']['min_pos'], self.settings['outer_loop']['max_pos']]

        if data is None:
            data = self.data

        if data:
            plot_fluorescence_pos(data['counts'], extent, axes_list[0])

    def _update_plot(self, axes_list):
        update_fluorescence(self.data['counts'], axes_list[0])


class AttoScanGridGeneric(Script):

    ''''

    FF 20211012

    Base script for taking measurements while scanning Attocube. This uses SetAttoANC300 to move the Attocube at each point in the scan grid.

    '''

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
        Parameter('Tracking', [
            Parameter('on/off', True, bool, 'used to turn on tracking'),
            Parameter('every_N', 1, int, 'track every n points')]),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn')
    ]

    _SCRIPTS = {'find_nv': FindNV, 'ESR_simple': ESR_simple}
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

        print(self.x_array)
        print(self.y_array)

        self.data['x_array'] = self.x_array
        self.data['y_array'] = self.y_array

        # Check if there are negative voltages in the scan range
        negative_voltage_err = 'Piezo voltage cannot be < 0 V'
        if np.any(self.x_array < 0):
            raise ValueError(negative_voltage_err)
        elif np.any(self.y_array < 0):
            raise ValueError(negative_voltage_err)

        self.x_array = self.x_array.tolist()
        self.y_array = self.y_array.tolist()

        self.point_index = 0
        for yNum in range(0, Ny):
            for xNum in range(0, Nx):
                if self._abort:
                    break
                print('Moving piezo')
                self.move_piezo(self.x_array[xNum], self.y_array[yNum])
                time.sleep(.5)
                print('Finished moving piezo')
                #self.data['actual_Vx'][yNum*Nx+xNum], self.data['actual_Vy'][yNum*Nx+xNum] = self.read_piezo()

                self.point_index += 1
                if self.settings['Tracking']['on/off'] and int(self.point_index) % self.settings['Tracking']['every_N'] == 0:
                    self.scripts['find_nv'].run()

                point_value = self.read_point()
                self.data['point_value'][yNum, xNum] = point_value

                self.progress = float(yNum * Nx + 1 + xNum) / (Nx * Ny) * 100

                #print(('Current acquisition {:02d}/{:02d} ({:0.2f}%)'.format(yNum * Nx + xNum, Nx * Ny,
                #                                                             self.progress)))
                print('Current acquisition {:02d}/{:02d}'.format(yNum * Nx + xNum, Nx * Ny))

                #self.updateProgress.emit(int(self.progress))

        # set end position after scan based on ending_behavior setting
        if self.settings['ending_behavior'] == 'leave_at_corner':
            return
        elif self.settings['ending_behavior'] == 'return_to_start':
            self.move_piezo(self.settings['point_a']['x'], self.settings['point_a']['y'])
        elif self.settings['ending_behavior'] == 'return_to_origin':
            self.move_piezo(0, 0)

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

    def move_piezo(self, x, y):
        raise NotImplementedError

class AttoScanGridEsr(AttoScanGridGeneric):
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
        Parameter('Tracking', [
            Parameter('on/off', True, bool, 'used to turn on tracking'),
            Parameter('every_N', 1, int, 'track every n points')]),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn')
    ]

    _SCRIPTS = {'find_nv': FindNV, 'ESR_simple': ESR_simple, 'SetAttoANC300': SetAttoANC300}

    def move_piezo(self, x, y):
        self.scripts['SetAttoANC300'].settings['point'].update({'x': float(x), 'y': float(y)})
        self.scripts['SetAttoANC300'].run()
        time.sleep(.1) # Wait for Find NV to settle

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



class AttoScanGridGalvoScan(AttoScanGridEsr):
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
        Parameter('Tracking', [
            Parameter('on/off', True, bool, 'used to turn on tracking'),
            Parameter('every_N', 1, int, 'track every n points')]),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn')
    ]

    _SCRIPTS = {'find_nv': FindNV, 'take_image': GalvoScan, 'SetAttoANC300': SetAttoANC300}

    def move_piezo(self, x, y):
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

class AttoScanGridEsrPiezoController(AttoScanGridGeneric):
    _SCRIPTS = {'find_nv': FindNV, 'ESR_simple': ESR_simple}
    _INSTRUMENTS = {'piezo_controller': PiezoController}

    def move_piezo(self, x, y):
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


class AttoOffsetZ(Script):
    """
    Sets the offset of Z-axis Attocube.
    Deliberately left out the other axes to avoid user error
    """
    _DEFAULT_SETTINGS = [
        Parameter('z_offset', 0, float, 'Offset voltage for Z-axis Attocube')
    ]

    _INSTRUMENTS = {'ANC300': ANC300}
    _SCRIPTS = {}

    def _function(self):
        """
        Performs a multiple attocube step with the voltage and frequency specified in instrument,
        and the direction and number specified in settings
        """
        attocube = self.instruments['ANC300']['instance']
        attocube._set_offset(3, self.settings['z_offset'])
        self.log('Z-axis Attocube offset changed to %.1f'%attocube.read_probes('z_offset'))

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
