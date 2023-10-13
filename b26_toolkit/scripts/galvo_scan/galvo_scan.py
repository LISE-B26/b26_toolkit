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
from b26_toolkit.scripts.set_laser import SetLaser
from pylabcontrol.core import Script, Parameter
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence

class GalvoScan(GalvoScanGeneric):

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
        Parameter('time_per_pt', .002, [.0005, .001, .002, .005, .01, .015, .02, .05, .08, .1], 'time in s to measure at each point'),
        Parameter('settle_time', .0002, [.0002, .0005], 'wait time between points to allow galvo to settle'),
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
        '''
        Initializes GalvoScan script for use in gui

        Args:
            instruments: list of instrument objects
            name: name to give to instantiated script object
            settings: dictionary of new settings to pass in to override defaults
            log_function: log function passed from the gui to direct log calls to the gui log
            data_path: path to save data

        '''
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
            self.x_array = np.repeat(self.x_array,self.clockAdjust)
        elif self._ACQ_TYPE == 'point':
            pass

        sample_rate = float(1) / self.settings['settle_time']
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate
        self.daq_in.settings['digital_input'][
            self.settings['DAQ_channels']['counter_channel']]['sample_rate'] = sample_rate

    def read_line(self, y_pos):
        self.initPt = [self.x_array[0], y_pos]
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

        # Create a list of 0s with the same length as the number of points in x dir requested
        summedData = np.zeros(int(len(self.x_array) / self.clockAdjust))
        for i in range(len(summedData)):
            summedData[i] = np.sum(
                diffData[(i * self.clockAdjust + 1):((i+1) * self.clockAdjust - 1)])

        # also normalizing to kcounts/sec
        normalizedData = (summedData * (.001 / self.settings['time_per_pt']))
        if 'safety_threshold' in self.settings and any(pixel > self.settings['safety_threshold'] for pixel in normalizedData):
            self.safety_threshold_exceeded = True
            self.safety_exceeded_behavior()
        return normalizedData

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



class GalvoScanSafe(GalvoScan):

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
        Parameter('time_per_pt', .002, [.0005, .001, .002, .005, .01, .015, .02, .05, .08, .1], 'time in s to measure at each point'),
        Parameter('settle_time', .0002, [.0002, .0005], 'wait time between points to allow galvo to settle'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('DAQ_channels',
                   [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
                    Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
                    Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'], 'Daq channel used for counter')
                  ]),
        Parameter('ending_behavior', 'leave_at_corner', ['return_to_start', 'return_to_origin', 'leave_at_corner'], 'return to the corn'),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan'),
        Parameter('safety_threshold', 10, float, 'The PulseBlaster will turn off the laser when a line contains a pixel exceeding this threshold (kCt/s)')
    ]

    _INSTRUMENTS = {'NI6259':  NI6259, 'NI9263': NI9263, 'NI9402': NI9402, 'PB': B26PulseBlaster}

    _SCRIPTS = {}

    def safety_exceeded_behavior(self):
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.log('Line contains a pixel exceeding the safety threshold! Turned off laser with PulseBlaster.')
        self._abort = True


class AttoScanAO(GalvoScan):
    """
    Basically the same as a galvo scan, but sends the analog voltages from the DAQ to a piezo controller which then goes to the Attocube
    Requires the piezo controller instrument to read the voltage gain setting, and adjust the analog voltage accordingly
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
                  [Parameter('x', 64, int, 'number of x points to scan'),
                   Parameter('y', 64, int, 'number of y points to scan')
                   ]),
        Parameter('time_per_pt', .01, [.01, .015, .02, .025, .03, .035, .04, .05, .1],
                  'time in s to measure at each point'),
        Parameter('settle_time', .01, [.001, .005, .01], 'wait time between points to allow galvo to settle, \n'
                                                         'time_per_pt must be divisible by settle_time!'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('DAQ_channels',
                  [Parameter('x_ao_channel', 'ao2', ['ao0', 'ao1', 'ao2', 'ao3'],
                             'Daq channel used for x voltage analog output'),
                   Parameter('y_ao_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'],
                             'Daq channel used for y voltage analog output'),
                   Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'],
                             'Daq channel used for counter')
                   ]),
        Parameter('ending_behavior', 'leave_at_corner', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn'),
        #Parameter('ending_behavior', 'return_to_start', ['return_to_start'],
        #          'sets Attocube location after scan is complete'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan'),

        Parameter('microwaves',
                  [Parameter('enable', False, bool, 'enable microwaves'),
                   Parameter('power_out', -45.0, float, 'output power (dBm)'),
                   Parameter('freq', 2.82e9, float, 'frequency of MW'),
                   Parameter('normalize_scan', False, bool, 'Normalize scan by taking each line twice, with MW off and on'),
                   Parameter('turn_off_after', False, bool, 'If true MW output is turned off after the measurement'),
                   Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')])
    ]
    _INSTRUMENTS = {'piezo_controller': PiezoController, 'microwave_generator': MicrowaveGenerator,
                    'NI6259':  NI6259, 'NI9263': NI9263, 'NI9402': NI9402}

    #_SCRIPTS = {'set_atto': SetAtto}

    def check_bounds(self):
        if np.any(self.x_array < 0) or np.any(self.y_array < 0):
            self.log('Attocube cannot accept negative voltages!')
            raise AttributeError

        # Spec says not to send in >1V/ms to DC-in, which is equiv to asking ANC300 to output 15V/ms
        # In our code that's equiv to speed = 15000
        self.speed_max = 1*self.scale()/0.001
        self.speed = (np.max(self.x_array) - np.min(self.x_array)) * self.scale() / float(self.settings['num_points']['x']) / self.settings['time_per_pt']
        self.speed_ratio = self.speed/self.speed_max

        if self.speed > self.speed_max * .5:
            self.log('Scan speed of %.1f V/s exceeds safety limit!'%self.speed)
            raise AttributeError
        else:
            self.log('Scan speed: %.1f V/s'%self.speed)

    def scale(self):
        voltage_limit = int(self.instruments['piezo_controller']['instance'].read_probes('voltage_limit'))
        if voltage_limit == 75:
            scale = 7.5
        elif voltage_limit == 100:
            scale = 10
        elif voltage_limit == 150:
            scale = 15
        return scale

    def before_scan(self):
        if self.settings['microwaves']['enable']:
            self.instruments['microwave_generator']['instance'].update({'frequency': float(self.settings['microwaves']['freq'])})
            self.instruments['microwave_generator']['instance'].update({'amplitude': self.settings['microwaves']['power_out']})
            self.instruments['microwave_generator']['instance'].update({'enable_modulation': False})
            self.instruments['microwave_generator']['instance'].update({'enable_output': True})
        else:
            self.instruments['microwave_generator']['instance'].update({'enable_output': False})

    def after_scan(self):
        if self.settings['microwaves']['turn_off_after']:
            self.instruments['microwave_generator']['instance'].update({'enable_output': False})

    def read_line_wrapper(self, y_pos):
        """
        If normalization is on, call read_line twice, with MW off and on.
        Args:
            y_pos: y position of the scan

        Returns:

        """
        if self.settings['microwaves']['normalize_scan'] and self.settings['microwaves']['enable']:
            self.instruments['microwave_generator']['instance'].update(
                {'frequency': float(self.settings['microwaves']['freq'])})
            #self.instruments['microwave_generator']['instance'].update({'enable_output': True})
            time.sleep(self.settings['settle_time'] * self.settings['num_points']['x'])
            data_MW_on = self.read_line(y_pos)
            self.instruments['microwave_generator']['instance'].update(
                {'frequency': float(2.87e9)})
            #self.instruments['microwave_generator']['instance'].update({'enable_output': False})
            time.sleep(self.settings['settle_time']*self.settings['num_points']['x'])
            data_MW_off = self.read_line(y_pos)

            line_data = (data_MW_on - data_MW_off)/(data_MW_off)
        else:
            line_data = self.read_line(y_pos)
        return line_data

    #def set_galvo_location(self, galvo_position):
    #    self.scripts['set_atto'].settings['initial_point'] = galvo_position
    #    self.scripts['set_atto'].run()

    def before_line_scan(self):
        pass

    def read_line(self, y_pos):
        self.initPt = [self.x_array[0], y_pos]
        self.daq_out.set_analog_voltages(
            {self.settings['DAQ_channels']['x_ao_channel']: self.initPt[0],
             self.settings['DAQ_channels']['y_ao_channel']: self.initPt[1]})

        self.before_line_scan()

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

        # Create a list of 0s with the same length as the number of points in x dir requested
        summedData = np.zeros(int(len(self.x_array) / self.clockAdjust))
        for i in range(len(summedData)):
            summedData[i] = np.sum(
                diffData[(i * self.clockAdjust + 1):((i+1) * self.clockAdjust - 1)])

        #################
        # If forward scan speed is x of the max speed, do the reverse scan at 0.1 max speed
        if self.speed_ratio < 0.04:
            reverse_speedup = int(1/self.speed_ratio*0.04)
        else:
            reverse_speedup = 1

        sample_rate = float(1) / self.settings['settle_time'] * reverse_speedup

        if sample_rate > 90000:
            # Lower speed up rate if the required sampling rate is too high for the DAQ to handle
            sample_rate = 90000
            reverse_speedup = sample_rate * self.settings['settle_time']

        #print("Moving back to start of x_line with %i X forward speed" % reverse_speedup)

        # Change sample rate for speed up
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate
        self.daq_in.settings['digital_input'][
            self.settings['DAQ_channels']['counter_channel']]['sample_rate'] = sample_rate


        # Return to beginning of the line
        # initialize APD thread
        self.x_array_reverse = self.x_array[::-1]
        clktask = self.daq_in.setup_clock(
            self.settings['DAQ_channels']['counter_channel'],
            len(self.x_array_reverse) + 1)

        # move in the opposite direction

        _aotask = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel']],
                                        self.x_array_reverse, clktask)
        # start counter and scanning sequence
        self.daq_out.run(_aotask)
        self.daq_in.run(clktask)
        self.daq_out.waitToFinish(_aotask)
        self.daq_out.stop(_aotask)
        #_xLineData_reverse, _rev = self.daq_in.read(_aitask)  # read the ai data for this scan
        #self.daq_in.stop(_aitask)
        self.daq_in.stop(clktask)

        # Revert change to sample rate
        sample_rate = float(1) / self.settings['settle_time']
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate
        self.daq_in.settings['digital_input'][
            self.settings['DAQ_channels']['counter_channel']]['sample_rate'] = sample_rate
        #################

        # also normalizing to kcounts/sec
        return summedData * (.001 / self.time_per_pt_actual())

    def time_per_pt_actual(self):
        # Actual time spent measuring counts at each point. Normally the same as time_per_pt setting, but might be much shorter if using pulse sequences
        return self.settings['time_per_pt']

    def _plot(self, axes_list, data = None):
        """
        Plots the galvo scan image
        Args:
            axes_list: list of axes objects on which to plot the galvo scan on the first axes object
            data: data (dictionary that contains keys image_data, extent) if not provided use self.data
        """

        if data is None:
            data = self.data

        labels = ['B-field Contour Scan (f = %.3e Hz)' % self.settings['microwaves']['freq'], r'V$_x$ [V]', r'V$_y$ [V]', 'kcounts/sec']
        plot_fluorescence_new(data['image_data'], data['extent'], axes_list[0], max_counts=self.settings['max_counts_plot'], labels=labels)


class AttoScanANC300(AttoScanAO):
    _INSTRUMENTS = {'microwave_generator': MicrowaveGenerator, 'NI6259': NI6259, 'NI9263': NI9263, 'NI9402': NI9402, 'ANC300': ANC300}

    def scale(self):
        # DC-in functionality has a built-in gain of 15
        return 15

    def before_scan(self):
        self.instruments['ANC300']['instance']._set_mode(1, 'input')
        self.instruments['ANC300']['instance']._set_mode(2, 'input')
        super(AttoScanANC300, self).before_scan()



if __name__ == '__main__':
    script, failed, instruments = Script.load_and_append(script_dict={'GalvoScan': 'GalvoScan'})

    print(script)
    print(failed)
    # print(instruments)

