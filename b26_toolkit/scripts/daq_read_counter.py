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
from collections import deque
import numpy as np

from b26_toolkit.instruments import NI6259
from b26_toolkit.plotting.plots_1d import plot_counts, update_1d_simple, update_counts_vs_pos
from pylabcontrol.core import Parameter, Script


class Daq_Read_Counter(Script):
    """
This script reads the Counter input from the DAQ and plots it. Only implemented for the PCI DAQ!!!!
    """
    _DEFAULT_SETTINGS = [
        Parameter('integration_time', .25, float, 'Time per data point (s)'),
        Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1'], 'Daq channel used for counter'),
        Parameter('total_int_time', 3.0, float, 'Total time to integrate (s) (if -1 then it will go indefinitely)'), # added by ER 20180606
        Parameter('track_laser_power_photodiode1',
                  [
                      Parameter('on/off', False, bool,
                                'If true, measure and normalize out laser power drifts during daq_read_counter'),
                      Parameter('ai_channel', 'ai2', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4'],
                                'channel to use for analog input, to which the photodiode is connected')
                  ]),
        Parameter('track_laser_power_photodiode2',
                  [
                      Parameter('on/off', False, bool, 'If true, measure and save laser power drifts during daq_read_counter on this photodiode. Cant use both simultaneously'),
                      Parameter('ai_channel', 'ai4', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4'], 'channel to use for photodiode 2, cant be the same as the track_laser_power photodiode')
                  ])
    ]

    _INSTRUMENTS = {'daq': NI6259}

    _SCRIPTS = {
    }

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.data = {'counts': deque(), 'laser_power': deque(), 'normalized_counts': deque(), 'laser_power2': deque()}


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        if self.settings['track_laser_power_photodiode1']['on/off'] and self.settings['track_laser_power_photodiode2']['on/off']:
            print('cant use both photodiodes at the same time - only use one AI channel at a time, unfortunately :-(')
            return

        sample_rate = float(2) / self.settings['integration_time']
        normalization = self.settings['integration_time']/.001
        self.instruments['daq']['instance'].settings['digital_input'][self.settings['counter_channel']]['sample_rate'] = sample_rate
        self.data = {'counts': deque(), 'laser_power': deque(), 'normalized_counts': deque(), 'laser_power2': deque()}
        self.last_value = 0
        sample_num = 2

        task = self.instruments['daq']['instance'].setup_counter("ctr0", sample_num, continuous_acquisition=True)

        if self.settings['track_laser_power_photodiode1']['on/off'] == True:
            aitask = self.instruments['daq']['instance'].setup_AI(self.settings['track_laser_power_photodiode1']['ai_channel'], sample_num,
                                          continuous=True, # continuous sampling still reads every clock tick, here set to the clock of the counter
                                          clk_source=task)

        if self.settings['track_laser_power_photodiode2']['on/off'] == True:
            aitask2 = self.instruments['daq']['instance'].setup_AI(self.settings['track_laser_power_photodiode2']['ai_channel'], sample_num,
                                          continuous=True, # continuous sampling still reads every clock tick, here set to the clock of the counter
                                          clk_source=task)
            print('aitask2: ', aitask2)

        # maximum number of samples if total_int_time > 0
        if self.settings['total_int_time'] > 0:
            max_samples = np.floor(self.settings['total_int_time']/self.settings['integration_time'])

        # start counter and scanning sequence
        if (self.settings['track_laser_power_photodiode1']['on/off'] and not self.settings['track_laser_power_photodiode2']['on/off']):
            self.instruments['daq']['instance'].run(aitask)
        elif (self.settings['track_laser_power_photodiode2']['on/off'] and not self.settings['track_laser_power_photodiode1']['on/off']):
            self.instruments['daq']['instance'].run(aitask2)

        self.instruments['daq']['instance'].run(task)

        # ER 20180827 wait for at least one clock tick to go by to start with a full clock tick of acquisition time for the first bin
        time.sleep(self.settings['integration_time'])

        sample_index = 0 # keep track of samples made to know when to stop if finite integration time

        while True:
            if self._abort:
                break

            # TODO: this is currently a nonblocking read so we add a time.sleep at the end so it doesn't read faster
            # than it acquires, this should be replaced with a blocking read in the future
            if self.settings['track_laser_power_photodiode1']['on/off'] == True:
                raw_data_laser, num_read_laser = self.instruments['daq']['instance'].read(aitask)
            if self.settings['track_laser_power_photodiode2']['on/off'] == True:
                raw_data_laser2, num_read_laser2 = self.instruments['daq']['instance'].read(aitask2)

            raw_data, num_read = self.instruments['daq']['instance'].read(task)
            #skip first read, which gives an anomolous value
            if num_read.value == 1:
                self.last_value = raw_data[0] #update running value to last measured value to prevent count spikes
                time.sleep(2.0 / sample_rate)
                continue

            tmp_count = 0
            for value in raw_data:
                new_val = ((float(value) - self.last_value) / normalization)
                self.data['counts'].append(new_val)
                self.last_value = value
                if self.settings['track_laser_power_photodiode1']['on/off'] == True:
                    self.data['laser_power'].append(raw_data_laser[tmp_count])
                if self.settings['track_laser_power_photodiode2']['on/off'] == True:
                    self.data['laser_power2'].append(raw_data_laser2[tmp_count])

                tmp_count = tmp_count + 1

            if self.settings['total_int_time'] > 0:
                self.progress = sample_index/max_samples
            else:
                self.progress = 50.
            self.updateProgress.emit(int(self.progress))

            time.sleep(2.0 / sample_rate)
            sample_index = sample_index + 1
            if self.settings['total_int_time'] > 0. and sample_index >= max_samples: # if the maximum integration time is hit
                self._abort = True # tell the script to abort

        # clean up APD tasks
        self.instruments['daq']['instance'].stop(task)
        if self.settings['track_laser_power_photodiode1']['on/off'] == True:
            self.instruments['daq']['instance'].stop(aitask)
        if self.settings['track_laser_power_photodiode2']['on/off'] == True:
            self.instruments['daq']['instance'].stop(aitask2)

        self.data['counts'] = list(self.data['counts'])

        if self.settings['track_laser_power_photodiode1']['on/off'] == True:
            self.data['laser_power'] = list(self.data['laser_power'])
            self.data['normalized_counts'] = list(np.divide(np.multiply(self.data['counts'], np.mean(self.data['laser_power'])), self.data['laser_power']))
        if self.settings['track_laser_power_photodiode2']['on/off'] == True:
            self.data['laser_power2'] = list(self.data['laser_power2'])

    def plot(self, figure_list):
        super(Daq_Read_Counter, self).plot([figure_list[1]])

    def _plot(self, axes_list, data = None):
        # COMMENT_ME

        if data is None:
            data = self.data

        if data:
            if self.settings['track_laser_power_photodiode1']['on/off'] == True:
                array_to_plot = np.delete(np.divide(np.multiply(self.data['counts'], np.mean(self.data['laser_power'])), self.data['laser_power']),0)
            else:
                array_to_plot = np.delete(data['counts'], 0)

            plot_counts(axes_list[0], array_to_plot)

    def _update_plot(self, axes_list, data = None):
        if data is None:
            data = self.data

        if data:
            if self.settings['track_laser_power_photodiode1']['on/off'] == True:
                array_to_plot = np.delete(np.divide(np.multiply(self.data['counts'], np.mean(self.data['laser_power'])), self.data['laser_power']), 0)
            else:
                array_to_plot = np.delete(data['counts'], 0)

            update_counts_vs_pos(axes_list[0], array_to_plot, np.linspace(0, len(array_to_plot), len(array_to_plot)))

if __name__ == '__main__':
    script = {}
    instr = {}
    script, failed, instr = Script.load_and_append({'Daq_Read_Cntr': 'Daq_Read_Cntr'}, script, instr)

    print(script)
    print(failed)
    print(instr)