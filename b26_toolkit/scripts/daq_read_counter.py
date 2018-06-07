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
from b26_toolkit.plotting.plots_1d import plot_counts
from pylabcontrol.core import Parameter, Script


class Daq_Read_Counter(Script):
    """
This script reads the Counter input from the DAQ and plots it.
    """
    _DEFAULT_SETTINGS = [
        Parameter('integration_time', .25, float, 'Time per data point (s)'),
        Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1'], 'Daq channel used for counter'),
        Parameter('total_int_time', 3.0, float, 'Total time to integrate (s) (if -1 then it will go indefinitely)') # added by ER 20180606
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

        self.data = {'counts': deque()}


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

      #  print(('settings in instrument', self.instruments['daq']['settings']))

        sample_rate = float(1) / self.settings['integration_time']
        normalization = self.settings['integration_time']/.001
        self.instruments['daq']['instance'].settings['digital_input'][self.settings['counter_channel']]['sample_rate'] = sample_rate
      #  print('setting sample rate')

        self.data = {'counts': deque()}

        self.last_value = 0

        sample_num = 2

      #  print(('settings in instrument 2', self.instruments['daq']['settings']))

     #   print(('here', self.instruments['daq']))

        task = self.instruments['daq']['instance'].setup_counter("ctr0", sample_num, continuous_acquisition=True)

        # maximum number of samples if total_int_time > 0
        if self.settings['total_int_time'] > 0:
            max_samples = np.floor(self.settings['total_int_time']/self.settings['integration_time'])

        # start counter and scanning sequence
        self.instruments['daq']['instance'].run(task)

        sample_index = 0 # keep track of samples made to know when to stop if finite integration time

        while True:
            if self._abort:
                break

            # TODO: this is currently a nonblocking read so we add a time.sleep at the end so it doesn't read faster
            # than it acquires, this should be replaced with a blocking read in the future
            raw_data, num_read = self.instruments['daq']['instance'].read(task)
            #skip first read, which gives an anomolous value
            if num_read.value == 1:
                self.last_value = raw_data[0] #update running value to last measured value to prevent count spikes
                time.sleep(2.0 / sample_rate)
                continue
           # print(('raw data length: ', len(raw_data)))
            for value in raw_data:
                self.data['counts'].append(((float(value) - self.last_value) / normalization))
                self.last_value = value
            self.progress = 50.
            self.updateProgress.emit(int(self.progress))

            time.sleep(2.0 / sample_rate)
            sample_index = sample_index + 1
            if self.settings['total_int_time'] > 0. and sample_index >= max_samples: # if the maximum integration time is hit
                self._abort = True # tell the script to abort

        # clean up APD tasks
        self.instruments['daq']['instance'].stop(task)
        self.data['counts'] = list(self.data['counts'])

    def plot(self, figure_list):
        super(Daq_Read_Counter, self).plot([figure_list[1]])

    def _plot(self, axes_list, data = None):
        # COMMENT_ME

        if data is None:
            data = self.data

        if data:
            plot_counts(axes_list[0], data['counts'])

if __name__ == '__main__':
    script = {}
    instr = {}
    script, failed, instr = Script.load_and_append({'Daq_Read_Cntr': 'Daq_Read_Cntr'}, script, instr)

    print(script)
    print(failed)
    print(instr)