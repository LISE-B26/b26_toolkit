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

from b26_toolkit.instruments import NI6259, NI9215, B26PulseBlaster
from b26_toolkit.plotting.plots_1d import plot_voltage
from pylabcontrol.core import Parameter, Script


class Daq_Read_Analog(Script):
    """
This script reads the analog input from the DAQ and plots it. ER 20180626

ER made changes to the daq 20190325 without testing it -- new code (commented with ER and the date) needs to be tested!!

    """
    _DEFAULT_SETTINGS = [
        Parameter('sampling_rate', 100, float, 'rate (Hz) at which samples are read by DAQ'),
        Parameter('ai_channel', 'ai0', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4'], 'Daq channel used for counter'),
        Parameter('total_int_time', 3.0, float, 'Total time to integrate (s) (if -1 then it will go indefinitely)'),
        Parameter('max_len_to_plot', 1000, int, 'plots the last n samples'),
        Parameter('daq_read_rate', 2, float, 'rate (Hz) at which samples are requested from the DAQ and plotted, default = 2'),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'daq to be used for pulse sequence'), # ER 20190325
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9215': NI9215}

    _SCRIPTS = {}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.data_to_plot = {'voltage': deque(maxlen = self.settings['max_len_to_plot'])}
        self.data['voltage'] = list()


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

      #  print(('settings in instrument', self.instruments['daq']['settings']))

        # new ER 20190325
        if self.settings['daq_type'] == 'PCI':
            daq = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            daq = self.instruments['NI9215']['instance']

        self.data_to_plot = {'voltage': deque(maxlen=self.settings['max_len_to_plot'])}
        self.data['voltage'] = list()

        # maximum number of samples if total_int_time > 0
        if self.settings['total_int_time'] > 0:
            max_samples = int(np.floor(self.settings['total_int_time'] * self.settings['sampling_rate']))
            if max_samples == 0:
                self.log('Sampling rate too low to even get one sample during total measurement time')
                return
        else:
            max_samples = np.inf

        sample_num = int(self.settings['sampling_rate']/self.settings['daq_read_rate'])
        if sample_num == 0:
            sample_num = 1
        elif sample_num > max_samples:
            sample_num = max_samples

        daq.settings['analog_input'][self.settings['ai_channel']]['sample_rate'] = self.settings['sampling_rate']
        task = daq.setup_AI(self.settings['ai_channel'], sample_num, continuous=True)
        daq.run(task)

        sample_index = 0  # keep track of samples made to know when to stop if finite integration time

        while True:
            if self._abort:
                break

            # TODO: this is currently a nonblocking read so we add a time.sleep at the end so it doesn't read faster
            # than it acquires, this should be replaced with a blocking read in the future
            raw_data, num_read = daq.read(task)
            for value in raw_data:
                self.data_to_plot['voltage'].append(float(value))
                self.data['voltage'].append(float(value))
            if self.settings['total_int_time'] > 0:
                self.progress = sample_index/max_samples
            else:
                self.progress = 50.
            self.updateProgress.emit(int(self.progress))

            time.sleep(sample_num / self.settings['sampling_rate'])
            sample_index += sample_num
            if self.settings['total_int_time'] > 0. and sample_index >= max_samples: # if the maximum integration time is hit
                self._abort = True

        # clean up APD tasks
        daq.stop(task)

        self.data_to_plot['voltage'] = list(self.data_to_plot['voltage'])

    def plot(self, figure_list):
        super(Daq_Read_Analog, self).plot([figure_list[1]])

    def _plot(self, axes_list, data = None):
        # COMMENT_ME
#        axes_list[0].hold(False)
        if data is None:
            data = self.data_to_plot

        if data:
            axes_list[0].clear()
            plot_voltage(axes_list[0], data['voltage'])

if __name__ == '__main__':
    script = {}
    instr = {}

    print(script)
    print(instr)