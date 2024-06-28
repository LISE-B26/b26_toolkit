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
from b26_toolkit.instruments import ANC300, ANC350
from b26_toolkit.scripts.daq_read_counter import DaqReadCounterOld, DaqReadCounterNi6259

from b26_toolkit.plotting.plots_1d import plot_counts_vs_pos, update_counts_vs_pos
from b26_toolkit.plotting.plots_2d import plot_fluorescence_pos
from pylabcontrol.core import Parameter, Script
from b26_toolkit.scripts import FindNv
from b26_toolkit.scripts.autofocus import AutoFocusDAQ
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from b26_toolkit.scripts.set_laser import SetAtto


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

    _INSTRUMENTS = {'ANC300': ANC300}
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

    _INSTRUMENTS = {'ANC300': ANC300}
    #_INSTRUMENTS = {'ANC300': ANC300, 'ANC350': ANC350}
    _SCRIPTS = {'SetAttoANC300': SetAtto}

    def _function(self):
        """
        Performs a multiple attocube step with the voltage and frequency specified in instrument,
        and the direction and number specified in settings
        """
        attocube = self.instruments[self.settings['controller_type']]['instance']

        self.scripts['SetAttoANC300'].run()

        for axis in self.settings['num_steps']:
            attocube.multistep(axis, self.settings['num_steps'][axis])


class AttoStepXy(AttoStep):
    # Constrains Attostep so that only X & Y moves
    _DEFAULT_SETTINGS = [
        Parameter('controller_type', 'ANC350', ['ANC300', 'ANC350'], 'attocube controller model'),
        Parameter('num_steps', [
            Parameter('x', 0, int, 'num steps along x-axis'),
            Parameter('y', 0, int, 'num steps along y-axis'),
        ]
    )]


class AttoStepScanOpenLoop(Script):
    _DEFAULT_SETTINGS = [
        Parameter('controller_type', 'ANC350', ['ANC300', 'ANC350'], 'attocube controller model'),
        Parameter('scan_axis', 'y', ['x', 'y'], 'axis to scan on'),
        Parameter('num_steps', 100, int, 'number of points in the scan'),
        Parameter('num_points', 10, int, 'number of DAQ counter measurements'),
    ]

    _SCRIPTS = {'daq_read_counter': DaqReadCounterOld}
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
            axes_list[0].set_xlabel('[steps]')
            axes_list[0].set_ylabel('[kCt/s]')

            axes_list[0].plot(data['positions'][:len(data['counts'])], data['counts'], linewidth=2.0)

    def _update_plot(self, axes_list):
        update_counts_vs_pos(axes_list[0], self.data['counts'], self.data['positions'][:len(self.data['counts'])])


class AttoStepScanOpenLoopTrackNI6259(Script):

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

    _SCRIPTS = {'daq_read_counter': DaqReadCounterNi6259, 'find_nv': FindNv, 'autofocus': AutoFocusDAQ}
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


class AttoStepScanOpenLoop2D(Script):
    _DEFAULT_SETTINGS = [
        Parameter('controller_type', 'ANC350', ['ANC300', 'ANC350'], 'attocube controller model'),
        Parameter('outer_loop', [
            Parameter('scan_axis', 'x', ['x', 'y'], 'axis to scan on'),
            Parameter('num_steps', 100, int, 'number of points in the scan'),
            Parameter('num_points', 10, int, 'number of DAQ counter measurements')
        ]
    )]

    _SCRIPTS = {'inner_scan': AttoStepScanOpenLoop}
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


class AttoStepScanClosedLoop(Script):

    # 20211012 FF: Not finished
    _DEFAULT_SETTINGS = [
        Parameter('scan_axis', 'x', ['x', 'y'], 'axis to scan on'),
        Parameter('direction', 'positive', ['positive', 'negative'], 'direction to scan'),
        Parameter('num_points', 100, int, 'number of points in the scan'),
        Parameter('min_pos', 0., float, 'minimum position of scan (um)'),
        Parameter('max_pos', 5., float, 'maximum position of scan (um)'),
        Parameter('time_per_pt', 0.5, float, 'time to wait at each point (s)'),
    ]

    _SCRIPTS = {'daq_read_counter': DaqReadCounterOld, "attocube": ANC300}

    def __init__(self, name=None, settings=None, instruments=None, scripts=None, log_function=None, data_path=None):
        """
        Default script initialization
        """
        self.attocube = self.instruments['attocube']['instance']
        super().__init__(name=None, settings=None, instruments=None, scripts=None, log_function=None, data_path=None)

    def _get_scan_positions(self, params):
        return np.linspace(params['min_pos'], params['max_pos'], params['num_points'])

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


class AttoStepScanClosedLoop2d(Script):
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

    _SCRIPTS = {'attoscan': AttoStepScanClosedLoop}

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