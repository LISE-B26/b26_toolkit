"""
    This file is part of b26_toolkit, a PyLabControl add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    Foobar is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""
import time
import numpy as np

from b26_toolkit.src.plotting.plots_1d import plot_1d_simple_timetrace_ns, update_1d_simple

from PyLabControl.src.core import Script, Parameter
from b26_toolkit.src.instruments import PressureGauge


class RecordPressures(Script):
    # COMMENT_ME

    _DEFAULT_SETTINGS = [
        Parameter('time_interval', 60.0, float, 'Time between points (s)'),
    ]

    _INSTRUMENTS = {
        'pressure_gauge' : PressureGauge
    }

    _SCRIPTS = {}

    def __init__(self, instruments = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, log_function= log_function, data_path=data_path)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        gauge = self.instruments['pressure_gauge']['instance']
        self.data['time'] = []
        self.data['pressures'] = []
        time_index = 0

        while not self._abort:
            self.data['pressures'].append(gauge._get_pressure())
            self.data['time'].append(time_index * self.settings['time_interval'])
            time_index += 1
            time.sleep(self.settings['time_interval'])
            self.force_update()
            self.progress = 50
            self.updateProgress.emit(int(self.progress))


    def _plot(self, axes_list):
        '''
        Args:
            axes_list: list of axes objects on which to plot the keyseight spectrun on the first axes object
            data: data (dictionary that contains keys amplitudes, frequencies) if not provided use self.data
        '''

        axes_list[0].plot(self.data['time'], self.data['pressures'])
        axes_list[0].set_xlabel('time (s)')
        axes_list[0].set_ylabel('pressure (Torr)')
