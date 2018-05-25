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
from matplotlib import patches

from b26_toolkit.instruments import TemperatureController
from pylabcontrol.core import Script, Parameter
from collections import deque
import time

from b26_toolkit.plotting.plots_1d import plot_temperature
class ReadTemperatureLakeshore(Script):
    """
This script reads the temperature from the Lakeshore controller.
    """

    _DEFAULT_SETTINGS = [
        Parameter('sample_rate', 1, float, 'Rate at which temperature is recorded in seconds'),
        Parameter('max_samples', 10, int, 'Number of samples'),
    ]

    _INSTRUMENTS = {'temp_controller': TemperatureController}

    _SCRIPTS = {

    }

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        """
        init script function see superclass for details on the parameters
        """
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.data = None


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        sample_rate = self.settings['sample_rate']
        max_samples = self.settings['max_samples']

        self.data = {'temperature': deque()}
        c = 0
        while True:
            if self._abort:
                self.data['temperature'] = np.array(self.data['temperature'])
                break

            temperature = self.instruments['temp_controller']['instance'].temperature

            if len(temperature) ==2:
                temperature = temperature[0]
            else:
                print('warning! temperature reading failed. expected a list of length 2 received the following:')
                print(temperature)
                temperature = 0
            self.data['temperature'].append(temperature)
            self.progress = 50.
            self.updateProgress.emit(int(self.progress))


            c +=1
            if c>=max_samples and max_samples>0:
                self._abort = True
                self.data['temperature'] = np.array(self.data['temperature'])

            time.sleep(1.0 / sample_rate)

    def _plot(self, axes_list, data = None):

        if data is None:
            data = self.data
        if data:
            plot_temperature(axes_list[0], data['temperature'], self.settings['sample_rate'])

if __name__ == '__main__':
    pass