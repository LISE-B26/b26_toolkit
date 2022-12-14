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
import time

import numpy as np

from b26_toolkit.instruments import NI6229, NI6259, NI9263, NI9402
from pylabcontrol.core import Script, Parameter
from b26_toolkit.scripts.galvo_scan.galvo_scan_generic import GalvoScanGeneric

'''
FG: file copied from galvo_scan_photodiode.py and modified to work as AFM scanning script for B22 scanning probe setup
'''


class AFMScanClassic(GalvoScanGeneric): # ER 20181221

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
        Parameter('RoI_mode', 'center', ['corner'], 'mode to calculate region of interest.\n \
                                                           corner: pta and ptb are diagonal corners of rectangle.\n \
                                                           center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 126, int, 'number of x points to scan'),
                   Parameter('y', 126, int, 'number of y points to scan')
                   ]),
        Parameter('time_per_pt', .002,float, 'time in s to measure at each point'),
        Parameter('settle_time', .002,float, 'wait time between points to allow galvo to settle'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('DAQ_channels',
                   [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
                    Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
                    Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'], 'Daq channel counter to use as a clock for ao and ai'),
                    Parameter('ai_channel', 'ai2', ['ai0', 'ai1', 'ai2', 'ai3'], 'Daq channel used for photodiode voltage')
                  ]),
        #Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'], 'return to the corn'),
        Parameter('ending_behavior', 'leave_at_corner', ['leave_at_corner'], 'return to the corn'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

   # _INSTRUMENTS = {'NI6259':  NI6259, 'NI9263': NI9263, 'NI9402': NI9402} when pulse_blaster_base_script uses 'daq' as the key for daq, we can't use these settings
    _INSTRUMENTS = {'daq': NI6229} # change daq names for when using this script on computers that aren't Alice

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

        # ER commented out 20181221 - still need to merge daq code with other scripts correctly
        #device_list = NI6259.get_connected_devices()
        #if not (self.instruments['NI6259']['instance'].settings['device'] in device_list):
        #    self.settings['daq_type'] = 'cDAQ'
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


        # ER commented out 20181221 - still need to merge daq code with other scripts correctly

        # defines which daqs contain the input and output based on user selection of daq interface
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['daq']['instance']
            self.daq_out = self.instruments['daq']['instance']
        elif self.settings['daq_type'] == 'cDAQ': # ER 20181221 these lines will fail on other computers since key for daq instrument is now just 'daq'
            self.daq_in = self.instruments['NI9402']['instance']
            self.daq_out = self.instruments['NI9263']['instance']

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

        self.clockAdjust = int(
            (self.settings['time_per_pt'] + self.settings['settle_time']) / self.settings['settle_time'])

        # overwrites existing self.x_array to include clockAdjust
        self.x_array = np.repeat(self.x_array,self.clockAdjust)
        sample_rate = float(1) / self.settings['settle_time']
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate
        self.daq_in.settings['analog_input'][
            self.settings['DAQ_channels']['ai_channel']]['sample_rate'] = sample_rate

    def after_scan(self):
        print('point a x')
        print(self.settings['point_a']['x'])
        print(self.x_array)
        x = np.linspace(self.x_array[0], self.settings['point_a']['x'] / self.scale(), 10)
        x = np.repeat(x, self.clockAdjust)
        print('x')
        print(x)
        clktask = self.daq_in.setup_clock(
            self.settings['DAQ_channels']['counter_channel'],len(x) + 1)
        aotask = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel']],x, clktask)
        self.daq_out.run(aotask)
        self.log('moving back to start along x')
        self.daq_in.run(clktask)
        self.daq_out.waitToFinish(aotask)
        self.daq_out.stop(aotask)
        self.daq_in.stop(clktask)

        print(self.y_current)
        y = np.linspace(self.y_current,self.settings['point_a']['y']/self.scale(), 10)
        y = np.repeat(y, self.clockAdjust)
        print('y')
        print(y)
        clktask = self.daq_in.setup_clock(
            self.settings['DAQ_channels']['counter_channel'],len(y) + 1)
        aotask = self.daq_out.setup_AO([self.settings['DAQ_channels']['y_ao_channel']],y, clktask)
        self.daq_out.run(aotask)
        self.log('moving back to start along y')
        self.daq_in.run(clktask)
        self.daq_out.waitToFinish(aotask)
        self.daq_out.stop(aotask)
        self.daq_in.stop(clktask)
        self.log("returned to start")

    def read_line(self, y_pos):
        # reverse direction of scan for uneven line numbers
        #linenumber = np.where(self.y_array==y_pos)[0]
        #print("scanning linenumber: %i"%linenumber)
        #if linenumber%2 !=0:
        #    isinstance(self.x_array,np.ndarray)
        self.x_array = self.x_array
        #else:
        #self.x_array = np.concatenate((self.x_array,self.x_array[::-1]))

        print("y_pos: %f"%y_pos)
        print("x:")
        print(self.x_array)

        self.initPt = [self.x_array[0], y_pos]
        self.daq_out.set_analog_voltages(
            {self.settings['DAQ_channels']['x_ao_channel']: self.initPt[0],
             self.settings['DAQ_channels']['y_ao_channel']: self.initPt[1]})

        # initialize APD thread
        clktask = self.daq_in.setup_clock(
            self.settings['DAQ_channels']['counter_channel'],
            len(self.x_array) + 1)

        # initialize the photodiode channel
        aitask = self.instruments['daq']['instance'].setup_AI(
            self.settings['DAQ_channels']['ai_channel'], len(self.x_array) + 1, continuous =False, clk_source=clktask)
        aotask = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel']],
                                       self.x_array, clktask)

        # start counter and scanning sequence
        print("moving right")
        self.daq_in.run(aitask)
        self.daq_out.run(aotask)
        self.daq_in.run(clktask)
        self.daq_out.waitToFinish(aotask)
        self.daq_out.stop(aotask)
     #   xLineData, _ = self.daq_in.read(ctrtask)
        xLineData, _ = self.daq_in.read(aitask) # read the ai data for this scan
        #print("xlineData")
        #print(np.asarray(xLineData))
        self.daq_in.stop(aitask)
        self.daq_in.stop(clktask)
      #  diffData = np.diff(xLineData)

        meanData = np.zeros(int(len(self.x_array) / self.clockAdjust))
        for i in range(0, int((len(self.x_array) / self.clockAdjust))):
            meanData[i] = np.mean(
               xLineData[(i * self.clockAdjust + 1):(i * self.clockAdjust + self.clockAdjust - 1)])
               # diffData[(i * self.clockAdjust + 1):(i * self.clockAdjust + self.clockAdjust - 1)])
        # also normalizing to kcounts/sec

        # also flip output values
        #if linenumber%2 !=0:
        #    meanData = meanData[::-1]

        # move tip back to starting position
        self.finPt = [self.x_array[-1], y_pos]
        self.daq_out.set_analog_voltages(
            {self.settings['DAQ_channels']['x_ao_channel']: self.finPt[0],
             self.settings['DAQ_channels']['y_ao_channel']: self.finPt[1]})
        # initialize APD thread
        clktask = self.daq_in.setup_clock(
            self.settings['DAQ_channels']['counter_channel'],
            len(self.x_array) + 1)

        # initialize the photodiode channel
        _aitask = self.instruments['daq']['instance'].setup_AI(
            self.settings['DAQ_channels']['ai_channel'], len(self.x_array) + 1, continuous =False, clk_source=clktask)
        # move in the opposite direction
        _aotask = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel']],
                                       self.x_array[::-1], clktask)

        # start counter and scanning sequence
        print("moving left")
        self.daq_in.run(_aitask)
        self.daq_out.run(_aotask)
        self.daq_in.run(clktask)
        self.daq_out.waitToFinish(_aotask)
        self.daq_out.stop(_aotask)
        _xLineData_reverse, _rev = self.daq_in.read(_aitask) # read the ai data for this scan
        #print("xlineData")
        #print(np.asarray(xLineData))
        self.daq_in.stop(_aitask)
        self.daq_in.stop(clktask)
#
        #meanData = np.split(meanData,2)[0]
        return meanData

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
        if galvo_position[0] > 10 or galvo_position[0] < 0 or galvo_position[1] > 10 or galvo_position[1] < 0:
            raise ValueError('The script attempted to set the nanomax to an illegal position outside of 0-10 V')

        pt = galvo_position
        # daq API only accepts either one point and one channel or multiple points and multiple channels
        pt = np.transpose(np.column_stack((pt[0], pt[1])))
        pt = (np.repeat(pt, 2, axis=1))

        task = self.daq_out.setup_AO(
            [self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], pt)
        self.daq_out.run(task)
        self.daq_out.waitToFinish(task)
        self.daq_out.stop(task)

    def check_bounds(self):
        if np.any(self.x_array > 10) or np.any(self.y_array > 10):
            self.log('Piezo control voltage for NanoMax are limited to +10 volts!!!')
            raise AttributeError
        if np.any(self.x_array < 0) or np.any(self.y_array < 0):
            self.log('Piezo control voltage for NanoMax need to be positive!!!')
            raise AttributeError

    def scale(self):
        # todo: ask DaLi what the units of the scale are!
        #  backpiezo input: amplification: 7.5
        #  so when driving on the back, the motion is 2um/V = 2nm/mV
        # Nanomax scale:
        var = self.instruments
        print(var)
        #voltage_limit_x = int(self.instruments['daq']['analog_output'][])
        #print(voltage_limit_x)
        # for NanoMax the Voltage limit is 75V
        voltage_limit = 75
        if voltage_limit == 75:
            scale = 7.5
        elif voltage_limit == 100:
            scale = 10
        elif voltage_limit == 150:
            scale = 15
        return scale

if __name__ == '__main__':
    script, failed, instruments = Script.load_and_append(script_dict={'AFMScanClassic': 'AFMScanClassic'})
    print(script)
    print(failed)
    # print(instruments)

