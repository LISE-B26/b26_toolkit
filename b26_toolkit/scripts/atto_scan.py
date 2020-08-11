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
from b26_toolkit.instruments import Attocube, AttocubeXY
from b26_toolkit.scripts.daq_read_counter import Daq_Read_Counter
from b26_toolkit.plotting.plots_1d import plot_counts_vs_pos, update_counts_vs_pos
from b26_toolkit.plotting.plots_2d import plot_fluorescence_pos, update_fluorescence
from pylabcontrol.core import Parameter, Script


class AttoStep(Script):
    # COMMENT_ME
    _DEFAULT_SETTINGS = [
        Parameter('axis', 'z', ['x', 'y', 'z'], 'Axis to step on'),
        Parameter('direction', 'Up', ['Up', 'Down'], 'step direction, up or down in voltage (or on physical switch)')
    ]

    _INSTRUMENTS = {'attocube': Attocube}
    _SCRIPTS = {}

    def __init__(self, instruments = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Default script initialization
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, log_function= log_function, data_path = data_path)

    def _function(self):
        """
        Performs a single attocube step with the voltage and frequency, and in the direction, specified in settings
        """
        attocube = self.instruments['attocube']['instance']
        attocube_voltage = self.instruments['attocube']['settings'][self.settings['axis']]['voltage']
        attocube.update({self.settings['axis']: {'voltage': attocube_voltage}})
        attocube_freq = self.instruments['attocube']['settings'][self.settings['axis']]['freq']
        attocube.update({self.settings['axis']: {'freq': attocube_freq}})
        if self.settings['direction'] == 'Up':
            dir = 0
        elif self.settings['direction'] == 'Down':
            dir = 1
        self.instruments['attocube']['instance'].step(self.settings['axis'], dir)

class AttoStepXY(AttoStep):
    _DEFAULT_SETTINGS = [
        Parameter('axis', 'x', ['x', 'y'], 'Axis to step on'),
        Parameter('direction', 'Up', ['Up', 'Down'], 'step direction, up or down in voltage (or on physical switch)')
    ]

    _INSTRUMENTS = {'attocube': AttocubeXY}

class AttoScanOpenLoop(Script):

    _DEFAULT_SETTINGS = [
        Parameter('scan_axis', 'x', ['x', 'y'], 'axis to scan on'),
        Parameter('direction', 'positive', ['positive', 'negative'], 'direction to scan'),
        Parameter('num_steps', 100, int, 'number of points in the scan'),
        Parameter('num_points', 10, int, 'number of DAQ counter measurements'),
        Parameter('time_per_pt', 0.5, float, 'time to wait at each point (s)'),
    ]

    _SCRIPTS = {'daq_read_counter': Daq_Read_Counter, 'attocube': AttocubeXY}

    def __init__(self, name=None, settings=None, instruments=None, scripts=None, log_function=None, data_path=None):
        """
        Default script initialization
        """
        self.attocube = self.instruments['attocube']['instance']
        super().__init__(name=None, settings=None, instruments=None, scripts=None, log_function=None, data_path=None)

    def _function(self):
        scan_axis = self.settings['scan_axis']
        
        self.data = {
            'counts': [],
            'position': np.linspace(0, self.settings['num_steps'], self.settings['num_points'])
        }

        if self.settings['direction'] == 'negative':
            self.data['position'] = -self.data['position']

        for i in range(self.settings['num_steps']):
            if self._abort:
                break

            # run daq_read_counter or the relevant script to get fluorescence
            if i % (self.settings['num_steps'] // self.settings['num_points']) == 0:
                self.scripts['daq_read_counter'].run()
                time.sleep(self.settings['time_per_pt'])
                self.scripts['daq_read_counter'].stop()

                # add to output structures which will be plotted
                data = self.scripts['daq_read_counter'].data['counts']
                self.data['counts'].append(np.mean(data))

            self.attocube.step(scan_axis, self.settings['direction'] == 'negative')

            self.progress = i * 100. / self.settings['num_points']
            self.updateProgress.emit(int(self.progress))
        
        # clean up data, as in daq_read_counter
        self.data['counts'] = list(self.data['counts'])

    def _plot(self, axes_list, data=None):
        if data is None:
            data = self.data

        if data:
            axes_list[0].set_xlabel('steps')
            axes_list[0].set_ylabel('kCounts/sec')

            axes_list[0].plot(self.positions, data['counts'], linewidth=2.0)

    def _update_plot(self, axes_list):
        update_counts_vs_pos(axes_list[0], self.data['counts'], self.positions)

class AttoScanClosedLoop(Script):

    _DEFAULT_SETTINGS = [
        Parameter('scan_axis', 'x', ['x', 'y'], 'axis to scan on'),
        Parameter('direction', 'positive', ['positive', 'negative'], 'direction to scan'),
        Parameter('num_points', 100, int, 'number of points in the scan'),
        Parameter('min_pos', 0., float, 'minimum position of scan (um)'),
        Parameter('max_pos', 5., float, 'maximum position of scan (um)'),
        Parameter('time_per_pt', 0.5, float, 'time to wait at each point (s)'),
    ]

    _SCRIPTS = {'daq_read_counter': Daq_Read_Counter, "attocube": AttocubeXY}

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

class AttoScanOpenLoop_2D(Script):
    
    _DEFAULT_SETTINGS = [
        Parameter('outer_loop',
                [
                    Parameter('scan_axis', 'x', ['x', 'y'], 'axis to scan on'),
                    Parameter('direction', 'positive', ['positive', 'negative'], 'direction to scan'),
                    Parameter('num_points', 10, int, 'number of points in the scan'),
                    Parameter('num_steps', 100, int, 'minimum position of scan (um)'),
                    Parameter('time_per_pt', 0.5, float, 'time to wait at each point (s)'),
                ]),
        Parameter('inner_loop',
                [
                    Parameter('scan_axis', 'x', ['x', 'y'], 'axis to scan on'),
                    Parameter('direction', 'positive', ['positive', 'negative'], 'direction to scan'),
                    Parameter('num_points', 10, int, 'number of points in the scan'),
                    Parameter('num_steps', 100, int, 'minimum position of scan (um)'),
                    Parameter('time_per_pt', 0.5, float, 'time to wait at each point (s)'),
                ]),
        Parameter('time_per_pt', 0.5, float, 'time to wait at each point (s)'),
    ]

    _SCRIPTS = {'attoscan': AttoScanOpenLoop}

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

        self.outer_pos = np.linspace(0, outer_settings['num_steps'], outer_settings['num_points'])
        

        self.data = {
            'counts': np.zeros((outer_settings['num_points'], inner_settings['num_points'])),
            'positions': {
                'outer': np.linspace(0, outer_settings['num_steps'], outer_settings['num_points']),
                'inner': np.linspace(0, inner_settings['num_steps'], inner_settings['num_points'])
            }
        }

        if outer_settings['direction'] == 'negative':
            self.data['positions']['outer'] = -self.data['positions']['outer']

        if inner_settings['direction'] == 'negative':
            self.data['positions']['inner'] = -self.data['positions']['inner']


        for i_out in range(outer_settings['num_points']):

            if self._abort:
                break

            if i_out % (outer_settings['num_steps'] // outer_settings['num_points']) == 0:
                self.inner_scan.run()
                self.data['counts'][i_out, :] = self.inner_scan.data['counts']                
                self.inner_scan['direction'] = self._switch_direction(self.inner_scan['direction'])

            self.attocube.step(outer_settings['scan_axis'], outer_settings['direction'])

            self.progress = i_out * 100. / self.settings['outer_loop']['num_points']
            self.updateProgress.emit(int(self.progress))

        self.data['counts'] = list(self.data['counts'])

    def _plot(self, axes_list, data=None):
        # COMMENT_ME

        extent = [0, self.data['positions']['inner'][-1], 0, self.data['positions']['outer'][-1]]

        if data is None:
            data = self.data

        if data:
            plot_fluorescence_pos(data['counts'], extent, axes_list[0])
            axes_list[0].set_xlabel(r'pos$_{outer}$ [steps]')
            axes_list[0].set_ylabel(r'pos$_{inner}$ [steps]')

    def _update_plot(self, axes_list):
        extent = [0, self.data['positions']['inner'][-1], 0, self.data['positions']['outer'][-1]]

        update_fluorescence(self.data['counts'], axes_list[0])

class AttoScanClosedLoop_2D(Script):

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

        extent = [self.settings['inner_loop']['min_pos'], self.settings['inner_loop']['max_pos'], self.settings['outer_loop']['min_pos'], self.settings['outer_loop']['max_pos']]

        if data is None:
            data = self.data

        if data:
            plot_fluorescence_pos(data['counts'], extent, axes_list[0])

    def _update_plot(self, axes_list):
        update_fluorescence(self.data['counts'], axes_list[0])

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
