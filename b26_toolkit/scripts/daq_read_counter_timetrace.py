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

from b26_toolkit.instruments import NI6259, NI9402
from b26_toolkit.plotting.plots_1d import plot_counts, update_counts
from pylabcontrol.core import Parameter, Script


class Daq_Read_Counter_TimeTrace(Script):
    """
This script reads the Counter input from the DAQ for a give duration and plots the time trace.
Future: plot also the PSD
    """
    _DEFAULT_SETTINGS = [
        Parameter('integration_time', .25, float, 'Time per data point (s)'),
        Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1'], 'Daq channel used for counter'),
        Parameter('total_int_time', 3.0, float, 'Total time to integrate (s)'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for counting')
    ]

    _INSTRUMENTS = {'NI6259':  NI6259, 'NI9402': NI9402}

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

        self.data = {'counts':[]}


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        if self.settings['daq_type'] == 'PCI':
            self.daq = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq = self.instruments['NI9402']['instance']


        sample_rate = float(1) / self.settings['integration_time']
        # normalization = self.settings['integration_time']/.001
        self.daq.settings['digital_input'][self.settings['counter_channel']]['sample_rate'] = sample_rate

        # maximum number of samples if total_int_time > 0
        if self.settings['total_int_time'] > 0:
            number_of_samples = int(np.floor(self.settings['total_int_time']/self.settings['integration_time']))
        else:
            self.log('total measurement time must be positive. Abort script')
            return


        # initialize APD thread
        ctrtask = self.daq.setup_counter(
            self.settings['counter_channel'], number_of_samples + 1)

        # start counter and scanning sequence
        self.daq.run(ctrtask)

        print('JG asdadad I am here')
        data, _ = self.daq.read(ctrtask)

        print('JG data', data)
        self.daq.stop(ctrtask)
        counts = np.diff(data)  # counter gives the accumulated counts, thus the diff gives the counts per interval

        self.data['counts'] = counts

    def plot(self, figure_list):
        super(Daq_Read_Counter_TimeTrace, self).plot([figure_list[1]])

    def _plot(self, axes_list, data = None):
        # COMMENT_ME

        if data is None:
            data = self.data

        if data:
            plot_counts(axes_list[0], data['counts'])

    def _update_plot(self, axes_list):
        if self.data:
            update_counts(axes_list[0], self.data['counts'])


if __name__ == '__main__':
    script = {}
    instr = {}
    script, failed, instr = Script.load_and_append({'Daq_Read_Cntr': 'Daq_Read_Cntr'}, script, instr)

    print(script)
    print(failed)
    print(instr)