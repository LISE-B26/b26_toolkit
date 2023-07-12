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

from b26_toolkit.instruments import NI6259, NI9263, NI9402, NI9215, NI9263_02
from pylabcontrol.core import Script, Parameter
from b26_toolkit.scripts.galvo_scan.galvo_scan_generic import GalvoScanGeneric


class GalvoScanPhotodiode(GalvoScanGeneric):

    '''

    This GalvoScan records the analog voltage coming out of a photodiode, and which is assumed to be connected to an
    analog input voltage on a daq.

    The original purpose of this scan was to set up the 'green interferometer' in a confocal, for the strings project.
    We would like to see the string motion in situ, and want to use the interface between the sample and the diamond
    as an interference pathway.

    However, we need to be sure to place the galvo position to the iron on the resonator to get the maximum signal.
    To do so we need to setup this galvo scan to make sure the photodiode is 'seeing' the same galvo xy voltages for the
    center of the iron piece as the APD does.

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
        Parameter('RoI_mode', 'center', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                                           corner: pta and ptb are diagonal corners of rectangle.\n \
                                                           center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 126, int, 'number of x points to scan'),
                   Parameter('y', 126, int, 'number of y points to scan')
                   ]),
        Parameter('settle_time', .001, float, 'wait time between points to allow galvo to settle'),
        Parameter('samples_per_point', 24, int, 'number of samples per pt, each sample takes as long as the settle time'),
        #Parameter('time_per_pt', .002, [.0005, .001, .002, .005, .01, .015, .02], 'time in s to measure at each point'),
        Parameter('statistics', 'mean', ['mean', 'std', 'min-max'], 'value to plot based on the samples'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        #Parameter('scanner', 'kinematic_mount', ['kinematic_mount', 'galvo'], 'choose whether to use galvo or kinematic mount to scan the sample'),
        Parameter('DAQ_channels',
                   [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
                    Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
                    Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'], 'Daq channel counter to use as a clock for ao and ai'),
                    Parameter('ai_channel', 'ai0', ['ai0', 'ai1', 'ai2', 'ai3'], 'Daq channel used for photodiode voltage')
                  ]),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'], 'return to the corn'),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    _INSTRUMENTS = {'NI9215': NI9215, 'NI9402': NI9402, 'NI9263': NI9263}

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

    def setup_scan(self):
        """
        setup the scan, i.e. identify the instruments and set up sample rate and such


        :return:
        """


        # ER commented out 20181221 - still need to merge daq code with other scripts correctly

        # defines which daqs contain the input and output based on user selection of daq interface
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['daq']['instance']
            self.daq_out = self.instruments['daq']['instance']

        elif self.settings['daq_type'] == 'cDAQ': # ER 20181221 these lines will fail on other computers since key for daq instrument is now just 'daq'
            self.daq_in = self.instruments['NI9215']['instance']

            if 'NI9263' in self.instruments:
                self.daq_out = self.instruments['NI9263']['instance']
            elif 'NI9263_02' in self.instruments:
                self.daq_out = self.instruments['NI9263_02']['instance']

            self.daq_clk = self.instruments['NI9402']['instance']

        #checks that requested daqs are actually physically present in the system
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

        #self.clockAdjust = int(
        #    (self.settings['time_per_pt'] + self.settings['settle_time']) / self.settings['settle_time'])
        self.clockAdjust = int(self.settings['samples_per_point'] + 1)

        # overwrites existing self.x_array to include clockAdjust
        self.x_array = np.repeat(self.x_array, self.clockAdjust)
        sample_rate = float(1) / self.settings['settle_time']
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate
        self.daq_in.settings['analog_input'][
            self.settings['DAQ_channels']['ai_channel']]['sample_rate'] = sample_rate

    def read_line(self, y_pos):
        self.initPt = [self.x_array[0], y_pos]
        self.daq_out.set_analog_voltages(
            {self.settings['DAQ_channels']['x_ao_channel']: self.initPt[0],
             self.settings['DAQ_channels']['y_ao_channel']: self.initPt[1]})

        # initialize APD thread
        clktask = self.daq_clk.setup_clock(
            self.settings['DAQ_channels']['counter_channel'],
            len(self.x_array) + 1)

        # initialize the photodiode channel
        aitask = self.daq_in.setup_AI(
            self.settings['DAQ_channels']['ai_channel'], len(self.x_array) + 1, continuous =False, clk_source=clktask)
        aotask = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel']],
                                       self.x_array, clktask)

        # start counter and scanning sequence
        self.daq_in.run(aitask)
        self.daq_out.run(aotask)
        self.daq_clk.run(clktask)
        self.daq_out.waitToFinish(aotask)
        self.daq_out.stop(aotask)
        xLineData, _ = self.daq_in.read(aitask) # read the ai data for this scan
        self.daq_in.stop(aitask)
        self.daq_clk.stop(clktask)
      #  diffData = np.diff(xLineData)

        processed_data = np.zeros(int(len(self.x_array) / self.clockAdjust))
        for i in range(0, int((len(self.x_array) / self.clockAdjust))):
            xLineData_block = xLineData[(i * self.clockAdjust + 1):(i * self.clockAdjust + self.clockAdjust - 1)]
            if self.settings['statistics'] == 'mean':
                processed_data[i] = np.mean(xLineData_block)
            elif self.settings['statistics'] == 'std':
                processed_data[i] = np.std(xLineData_block)
            elif self.settings['statistics'] == 'min-max':
                processed_data[i] = np.max(xLineData_block)-np.min(xLineData_block)
        return processed_data

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

    def scale(self):
        return 1

    def check_bounds(self):
        if np.any(np.abs(self.x_array) > 2) or np.any(np.abs(self.y_array) > 2):
            self.log('Exceeded galvo range!')
            raise AttributeError

class KinematicMountScanPhotodiode(GalvoScanPhotodiode):

    '''

    Same as GalvoScanPhotodiode, but uses a piezo-actuated kinematic mount instead of galvo mirrors.
    Sets different voltage scales (15x for piezo controller) and voltage bounds (non-negative)
    '''

    _DEFAULT_SETTINGS = [
        Parameter('point_a',
                  [Parameter('x', 75, float, 'x-coordinate'),
                   Parameter('y', 75, float, 'y-coordinate')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 150, float, 'x-coordinate'),
                   Parameter('y', 150, float, 'y-coordinate')
                   ]),
        Parameter('RoI_mode', 'center', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                                           corner: pta and ptb are diagonal corners of rectangle.\n \
                                                           center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 64, int, 'number of x points to scan'),
                   Parameter('y', 64, int, 'number of y points to scan')
                   ]),
        Parameter('settle_time', .001, float, 'wait time between points to allow galvo to settle'),
        Parameter('samples_per_point', 24, int, 'number of samples per pt, each sample takes as long as the settle time'),
        Parameter('statistics', 'mean', ['mean', 'std', 'min-max'], 'value to plot based on the samples'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('DAQ_channels',
                   [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
                    Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
                    Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'], 'Daq channel counter to use as a clock for ao and ai'),
                    Parameter('ai_channel', 'ai0', ['ai0', 'ai1', 'ai2', 'ai3'], 'Daq channel used for photodiode voltage')
                  ]),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'], 'return to the corn'),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    _INSTRUMENTS = {'NI9215': NI9215, 'NI9402': NI9402, 'NI9263_02': NI9263_02}

    def check_bounds(self):
        if np.any(self.x_array < 0) or np.any(self.y_array < 0):
            self.log('Piezo cannot accept negative voltages!')
            raise AttributeError
        if np.any(np.abs(self.x_array) > 150) or np.any(np.abs(self.y_array) > 150):
            self.log('Exceeded galvo range!')
            raise AttributeError

    def scale(self):
        return 15

if __name__ == '__main__':
    script, failed, instruments = Script.load_and_append(script_dict={'GalvoScanPhotodiode': 'GalvoScanPhotodiode'})

    print(script)
    print(failed)
    # print(instruments)

