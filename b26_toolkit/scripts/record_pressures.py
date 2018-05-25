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

from b26_toolkit.plotting.plots_1d import plot_1d_simple_timetrace_ns, update_1d_simple

from pylabcontrol.core import Script, Parameter
from b26_toolkit.instruments import ChamberPressureGauge, PumpLinePressureGauge, TemperatureController
from b26_toolkit.instruments import CryoStation

class RecordPressures(Script):
    # COMMENT_ME

    _DEFAULT_SETTINGS = [
        Parameter('time_interval', 60.0, float, 'Time between points (s)'),
    ]

    _INSTRUMENTS = {
        'chamber_pressure_gauge' : ChamberPressureGauge,
        'pump_line_pressure_gauge': PumpLinePressureGauge,
        'temp_controller': TemperatureController,
        'cryo_station': CryoStation
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

        chamber_gauge = self.instruments['chamber_pressure_gauge']['instance']
        pump_line_gauge = self.instruments['pump_line_pressure_gauge']['instance']
        temp_controller = self.instruments['temp_controller']['instance']
        cryo_station = self.instruments['cryo_station']['instance']

        self.data['time'] = []
        self.data['chamber_pressures'] = []
        self.data['pump_line_pressures'] = []
        self.data['temperatures'] = []
        self.data['temperatures_raw'] = []
        self.data['Platform_Temp'] = []
        self.data['Stage_1_Temp'] = []
        self.data['Stage_2_Temp'] = []

        time_index = 0

        while not self._abort:
            self.data['chamber_pressures'].append(chamber_gauge.pressure)
            self.data['pump_line_pressures'].append(pump_line_gauge.pressure)
            temp, raw = temp_controller.temperature
            self.data['temperatures'].append(temp)
            self.data['temperatures_raw'].append(raw)

            self.data['Platform_Temp'].append(cryo_station.Platform_Temp)
            self.data['Stage_1_Temp'].append(cryo_station.stage_1_temp)
            self.data['Stage_2_Temp'].append(cryo_station.stage_2_temp)


            self.data['time'].append(time_index * self.settings['time_interval'])
            time_index += 1
            self.force_update()
            self.progress = 50
            self.updateProgress.emit(int(self.progress))
            time.sleep(self.settings['time_interval'])


    def _plot(self, axes_list):
        '''
        Args:
            axes_list: list of axes objects on which to plot the keyseight spectrun on the first axes object
            data: data (dictionary that contains keys amplitudes, frequencies) if not provided use self.data
        '''

        time = self.data['time']

        if max(time)>3600:
            time = np.array(self.data['time'])/3600
            time_label = 'time (h)'
        elif max(time)>60:
            time = np.array(self.data['time'])/60
            time_label = 'time (min)'
        else:
            time = np.array(self.data['time'])
            time_label = 'time (s)'


        axes_list[1].plot(time, self.data['Platform_Temp'],
                          time, self.data['Stage_1_Temp'],
                          time, self.data['Stage_2_Temp']
                          )
        axes_list[1].set_xlabel(time_label)
        axes_list[1].set_ylabel('temparatures (K)')

        axes_list[1].legend(labels=('Platform', 'Stage 1', 'Stage 2'), fontsize=8)


        #10/13/16 AK: pump line connection was broken, temporarily comment out
        axes_list[0].plot(time, self.data['chamber_pressures'],
                          time, self.data['pump_line_pressures']
                          )
        # axes_list[0].plot(time, self.data['chamber_pressures'])
        axes_list[0].set_xlabel(time_label)
        axes_list[0].set_ylabel('pressure (Torr)')

        ax2 = axes_list[0].twinx()
        ax2.plot(time, self.data['temperatures'], 'r')
        ax2.set_ylabel('Temperature (K)', color='r')
        axes_list[0].legend(labels=('chamber', 'pump line'), fontsize=8)


