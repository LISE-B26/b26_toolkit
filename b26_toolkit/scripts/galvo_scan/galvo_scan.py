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

from b26_toolkit.instruments import NI6259, NI9263, NI9402, PiezoController, MicrowaveGenerator, ANC300, B26PulseBlaster
from b26_toolkit.scripts.galvo_scan.galvo_scan_generic import GalvoScanGeneric
from b26_toolkit.scripts.daq_read_counter_timetrace import DaqTimeTraceNi6259
from b26_toolkit.scripts.set_laser import SetLaser
from pylabcontrol.core import Parameter, Script
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence


class GalvoScan(GalvoScanGeneric):
    """
    Uses the DAQ analog outputs to sweep the voltages sent to the galvo mirrors, and collect the photon count at each voltage.
    This results in an image in the current field of view of the objective.
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
                  [Parameter('x', 64, int, 'number of x points to scan'),
                   Parameter('y', 64, int, 'number of y points to scan')
                   ]),
        Parameter('time_per_pt', .002, float, 'time in s to measure at each point'),
        Parameter('settle_time', .0002, float, 'wait time between points to allow galvo to settle'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('DAQ_channels',
                   [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
                    Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
                    Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'], 'Daq channel used for counter')
                    ]),
        Parameter('ending_behavior', 'leave_at_corner', ['return_to_start', 'return_to_origin', 'leave_at_corner'], 'return to the corn'),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    _INSTRUMENTS = {'NI6259':  NI6259, 'NI9263': NI9263, 'NI9402': NI9402}

    _SCRIPTS = {}

    def __init__(self, instruments, scripts = None, name=None, settings=None, log_function=None, data_path=None):
        """
        Initializes GalvoScan script for use in gui

        Args:
            instruments: list of instrument objects
            name: name to give to instantiated script object
            settings: dictionary of new settings to pass in to override defaults
            log_function: log function passed from the gui to direct log calls to the gui log
            data_path: path to save data
        """

        Script.__init__(self, name, settings=settings, instruments=instruments, scripts=scripts, log_function=log_function,
                        data_path=data_path)

        device_list = NI6259.get_connected_devices()
        if not (self.instruments['NI6259']['instance'].settings['device'] in device_list):
            self.settings['daq_type'] = 'cDAQ'
        # # defines which daqs contain the input and output based on user selection of daq interface
        # if self.settings['daq_type'] == 'PCI':
        #     self.daq_in = self.instruments['NI6259']['instance']
        #     self.daq_out = self.instruments['NI6259']['instance']
        # elif self.settings['daq_type'] == 'cDAQ':
        #     self.daq_in = self.instruments['NI9402']['instance']
        #     self.daq_out = self.instruments['NI9263']['instance']

    def setup_scan(self):
        """
        setup the scan, i.e. identify the instruments and set up sample rate and such


        :return:
        """
        # defines which daqs contain the input and output based on user selection of daq interface
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['NI6259']['instance']
            self.daq_out = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_in = self.instruments['NI9402']['instance']
            self.daq_out = self.instruments['NI9263']['instance']

        # checks that requested daqs are actually physically present in the system
        device_list = NI6259.get_connected_devices()
        if not(self.daq_in.settings['device'] in device_list):
            self.log('The requested input daq ' + self.daq_in.settings['device'] + ' is not connected to this computer. Possible daqs are '
                     + str(device_list) + '. Please choose one of these and try again.')
            raise AttributeError

        if not (self.daq_out.settings['device'] in device_list):
            self.log('The requested output daq ' + self.daq_out.settings[
                'device'] + ' is not connected to this computer. Possible daqs are '
                     + str(device_list) + '. Please choose one of these and try again.')
            raise AttributeError


        self.clockAdjust = int(
            (self.settings['time_per_pt'] + self.settings['settle_time']) / self.settings['settle_time'])

        if self._ACQ_TYPE == 'line':
            # overwrites existing self.x_array to include clockAdjust
            self.x_array = np.repeat(self.x_array, self.clockAdjust)
        elif self._ACQ_TYPE == 'point':
            pass

        sample_rate = float(1) / self.settings['settle_time']
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate
        self.daq_in.settings['digital_input'][
            self.settings['DAQ_channels']['counter_channel']]['sample_rate'] = sample_rate

    def check_bounds(self):
        if np.any(np.abs(self.x_array) > 2) or np.any(np.abs(self.y_array) > 2):
            self.log('GalvoScan range exceeded!')
            raise AttributeError

    def before_line_scan(self):
        """
        Action before scanning each line
        :return: None
        """
        pass

    def read_line(self, y_pos):
        self.before_line_scan()

        self.initPt = [self.x_array[0], y_pos]
        timer_start = time.time()

        self.daq_out.settings['analog_output'][self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = 10000
        self.daq_out.settings['analog_output'][self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = 10000

        self.daq_out.set_analog_voltages(
            {self.settings['DAQ_channels']['x_ao_channel']: self.initPt[0],
             self.settings['DAQ_channels']['y_ao_channel']: self.initPt[1]})

        #print('Time after set outputs: %.2e' % (time.time()-timer_start))

        # initialize APD thread
        ctrtask = self.daq_in.setup_counter(
            self.settings['DAQ_channels']['counter_channel'],
            len(self.x_array) + 1)

        #print('Time after setup counter: %.2e' % (time.time() - timer_start))
        aotask = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel']],
                                       self.x_array, ctrtask)

        #print('Time after setup AO: %.2e' % (time.time() - timer_start))

        # start counter and scanning sequence
        self.daq_out.run(aotask)
        #print('Time after run AO task: %.2e' % (time.time() - timer_start))
        self.daq_in.run(ctrtask)
        #print('Time after run counter task: %.2e' % (time.time() - timer_start))
        self.daq_out.waitToFinish(aotask)
        self.daq_out.stop(aotask)
        #print('Time after stop AO task: %.2e' % (time.time() - timer_start))
        xLineData, _ = self.daq_in.read(ctrtask)
        #print('Time after read counter task: %.2e' % (time.time() - timer_start))
        self.daq_in.stop(ctrtask)
        #print('Time after stop counter task: %.2e' % (time.time() - timer_start))
        diffData = np.diff(xLineData)

        # Create a list of 0s with the same length as the number of points in x dir requested
        summedData = np.zeros(int(len(self.x_array) / self.clockAdjust))
        for i in range(len(summedData)):
            summedData[i] = np.sum(
                diffData[(i * self.clockAdjust + 1):((i+1) * self.clockAdjust - 1)])

        # also normalizing to kcounts/sec
        normalizedData = (summedData * (.001 / self.time_per_pt_actual()))
        if 'safety_threshold' in self.settings:
            safety_threshold = self.settings['safety_threshold']
            if safety_threshold != -1 and any(pixel > safety_threshold for pixel in normalizedData):
                self.safety_threshold_exceeded = True
                self.safety_exceeded_behavior()
        return normalizedData

    def time_per_pt_actual(self):
        # Actual time spent measuring counts at each point. Normally the same as time_per_pt setting, but might be much shorter if using pulse sequences
        return self.settings['time_per_pt']


    def get_galvo_location(self):
        """
        Returns the current position of the galvo. Requires a daq with analog inputs internally routed to the analog
        outputs (ex. NI6259. Note that the cDAQ does not have this capability).
        Returns: list with two floats, which give the x and y position of the galvo mirror
        """
        if self.settings['daq_type'] == 'PCI':
            initial_position = self.daq_out.get_analog_voltages([self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']])
        else:
            initial_position = []
        return initial_position


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


class GalvoScanTimetrace(GalvoScanGeneric):
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
                  [Parameter('x', 64, int, 'number of x points to scan'),
                   Parameter('y', 64, int, 'number of y points to scan')
                   ]),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    _INSTRUMENTS = {}

    _SCRIPTS = {'Daq_timetrace': DaqTimeTraceNi6259, 'SetLaser': SetLaser}

    def __init__(self, name=None, settings=None, instruments=None, scripts=None, log_function=None, data_path=None):
        '''
        Initializes GalvoScan script for use in gui

        Args:
            instruments: list of instrument objects
            name: name to give to instantiated script object
            settings: dictionary of new settings to pass in to override defaults
            log_function: log function passed from the gui to direct log calls to the gui log
            data_path: path to save data

        '''
        Script.__init__(self, name=name, settings=settings, instruments=instruments, scripts=scripts, log_function=log_function, data_path=data_path)

    def setup_scan(self):
        """
        setup the scan, i.e. identify the instruments and set up sample rate and such
        """
        self.scripts['SetLaser'].settings['daq_type'] = 'PCI'
        self.data['counts'] = []

    def read_line(self, y_pos):
        set_laser_script = self.scripts['SetLaser']
        set_laser_script.settings['point']['y'] = y_pos

        line_data = np.zeros(len(self.x_array))

        self.data['counts'].append([])

        for i in range(len(self.x_array)):
            if self._abort:
                break

            set_laser_script.settings['point']['x'] = self.x_array[i]
            set_laser_script.run()

            self.scripts['Daq_timetrace'].run()
            counts = self.scripts['Daq_timetrace'].data['counts']

            line_data[i] = np.mean(counts)
            self.data['counts'][-1].append(counts)

        return line_data

    def get_galvo_location(self):
        """
        Returns the current position of the galvo. Requires a daq with analog inputs internally routed to the analog
        outputs (ex. NI6259. Note that the cDAQ does not have this capability).
        Returns: list with two floats, which give the x and y position of the galvo mirror
        """
        point = self.scripts['SetLaser'].settings['point']
        return [point['x'], point['y']]

    def set_galvo_location(self, galvo_position):
        """
        sets the current position of the galvo
        galvo_position: list with two floats, which give the x and y position of the galvo mirror
        """
        if galvo_position[0] > 1 or galvo_position[0] < -1 or galvo_position[1] > 1 or galvo_position[1] < -1:
            raise ValueError('The script attempted to set the galvo position to an illegal position outside of +- 1 V')

        self.scripts['SetLaser'].settings['point']['x'] = galvo_position[0]
        self.scripts['SetLaser'].settings['point']['y'] = galvo_position[1]


class GalvoScanSafe(GalvoScan):
    """
    Uses the DAQ analog outputs to sweep the voltages sent to the galvo mirrors, and collect the photon count at each voltage.
    This results in an image in the current field of view of the objective.
    This script will turn off the laser using the AOM if the photon count exceeds a safety threshold (e.g. accidentally hitting a micromagnet in vacuum)
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
                  [Parameter('x', 64, int, 'number of x points to scan'),
                   Parameter('y', 64, int, 'number of y points to scan')
                   ]),
        Parameter('time_per_pt', .002, float, 'time in s to measure at each point'),
        Parameter('settle_time', .0002, float, 'wait time between points to allow galvo to settle'),
        Parameter('safety_threshold', 10, float, 'PulseBlaster will turn off the laser when a line contains a pixel exceeding this threshold (kCt/s); '
                                                 'set to -1 to disable this feature'),
        Parameter('turn_off_laser_after', True, bool, 'turn off laser after scan'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('DAQ_channels',
                   [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
                    Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
                    Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'], 'Daq channel used for counter')
                  ]),
        Parameter('ending_behavior', 'leave_at_corner', ['return_to_start', 'return_to_origin', 'leave_at_corner'], 'return to the corn'),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9263': NI9263, 'NI9402': NI9402, 'PB': B26PulseBlaster}

    def safety_exceeded_behavior(self):
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.log('Line contains a pixel exceeding the safety threshold! Turned off laser with PulseBlaster.')
        self._abort = True

    def before_scan(self):
        """
        Runs something before starting the scan.
        This can be something like moving the laser to a certain position, running findNV, or turning on MW
        Returns:

        """
        self.instruments['PB']['instance'].update({'laser': {'status': True}})


    def after_scan(self):
        """
        Runs something after finishing the scan.
        This can be something like moving the laser to a certain position, running findNV, or turning off MW
        Returns:

        """
        if self.settings['turn_off_laser_after']:
            self.instruments['PB']['instance'].update({'laser': {'status': False}})


if __name__ == '__main__':
    script, failed, instruments = Script.load_and_append(script_dict={'GalvoScan': 'GalvoScan'})

    print(script)
    print(failed)
    # print(instruments)

