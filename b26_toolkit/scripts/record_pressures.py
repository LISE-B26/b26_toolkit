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
from datetime import datetime
import time
from pylabcontrol.core import Script, Parameter
from b26_toolkit.instruments import ChamberPressureGauge, TemperatureController

class RecordPressures(Script):
    # COMMENT_ME
    _DEFAULT_SETTINGS = [
        Parameter('time_interval', 60.0, float, 'Time between points (s)'),
    ]
    _INSTRUMENTS = {'chamber_pressure_gauge': ChamberPressureGauge}
    _SCRIPTS = {}

    def __init__(self, instruments=None, name=None, settings=None, log_function=None, data_path=None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings=settings, instruments=instruments, log_function=log_function, data_path=data_path)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """
        chamber_gauge = self.instruments['chamber_pressure_gauge']['instance']
        # pump_line_gauge = self.instruments['pump_line_pressure_gauge']['instance']
        # temp_controller = self.instruments['temp_controller']['instance']
        # cryo_station = self.instruments['cryo_station']['instance']
        self.data['time'] = []
        self.data['chamber_pressures'] = []
        # self.data['pump_line_pressures'] = []
        self.data['temperatures'] = []
        self.data['temperatures_raw'] = []
        self.data['Platform_Temp'] = []
        self.data['Stage_1_Temp'] = []
        self.data['Stage_2_Temp'] = []
        time_index = 0
        while not self._abort:
            self.data['chamber_pressures'].append(chamber_gauge.pressure)
            # self.data['pump_line_pressures'].append(pump_line_gauge.pressure)
            # temp, raw = temp_controller.temperature
            # self.data['temperatures'].append(temp)
            # self.data['temperatures_raw'].append(raw)
            # self.data['Platform_Temp'].append(cryo_station.Platform_Temp)
            # self.data['Stage_1_Temp'].append(cryo_station.stage_1_temp)
            # self.data['Stage_2_Temp'].append(cryo_station.stage_2_temp)
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
        if max(time) > 3600:
            time = np.array(self.data['time'])/3600
            time_label = 'time (h)'
        elif max(time) > 60:
            time = np.array(self.data['time'])/60
            time_label = 'time (min)'
        else:
            time = np.array(self.data['time'])
            time_label = 'time (s)'

        #10/13/16 AK: pump line connection was broken, temporarily comment out
        axes_list[0].semilogy(time, self.data['chamber_pressures'])

        axes_list[0].set_xlabel(time_label)
        axes_list[0].set_ylabel('pressure (Torr)')


class RecordPressuresTemperature(Script):
    # COMMENT_ME

    _DEFAULT_SETTINGS = [
        Parameter('time_interval', 5.0, float, 'Time between points (s)'),
    ]

    _INSTRUMENTS = {
        'chamber_pressure_gauge' : ChamberPressureGauge,
         'temp_controller': TemperatureController
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
        temp_controller = self.instruments['temp_controller']['instance']

        self.data['time'] = []
        self.data['time_index'] = []
        self.data['temperaturesA'] = []
        self.data['temperaturesB'] = []
        self.data['chamber_pressures'] = []

        time_index = 0

        while not self._abort:


            self.data['chamber_pressures'].append(chamber_gauge.pressure)
            tempA, tempB = temp_controller.read_probes('temperature_A'), temp_controller.read_probes('temperature_B')
            self.data['temperaturesA'].append(tempA)
            self.data['temperaturesB'].append(tempB)

            self.data['time_index'].append(time_index * self.settings['time_interval'])
            time_index += 1
            if len(self.data['time']) == 0:
                time_start = datetime.now()
            time_elapsed = (datetime.now()-time_start).total_seconds()
            self.data['time'].append(time_elapsed)

            self.force_update()
            self.progress = 50
            self.updateProgress.emit(int(self.progress))

            if self.settings['save'] and time_index % 10 == 0:
                self.save_data()
            time.sleep(self.settings['time_interval'])

    def _plot_old(self, axes_list):
        '''
        Args:
            axes_list: list of axes objects on which to plot the keyseight spectrun on the first axes object
            data: data (dictionary that contains keys amplitudes, frequencies) if not provided use self.data
        '''

        time = np.array(self.data['time'])
        if max(time)>3600:
            time /= 3600
            time_label = 'time (h)'
        elif max(time) > 60:
            time /= 60
            time_label = 'time (min)'
        else:
            time_label = 'time (s)'

        lns1 = axes_list[0].plot(time, self.data['chamber_pressures'], color='black',label='Chamber pressure')

        axes_list[0].set_yscale('log')
        axes_list[0].set_xlabel(time_label)
        axes_list[0].set_ylabel('Pressure (Torr)')

        ax2 = axes_list[0].twinx()
        lns2 = ax2.plot(time, self.data['temperaturesA'], color ='red',label='Temp A')
        lns3 = ax2.plot(time, self.data['temperaturesB'], color = 'blue',label='Temp B')
        ax2.set_ylabel('Temperature (K)')

        lns = lns1 + lns2 + lns3
        labs = [l.get_label() for l in lns]
        ax2.legend(lns, labs, loc=0)

        ax2.text(0.5, 0.05, 'Temp A: %.3f Temp B: %.3f Pressure: %.3e'%(self.data['temperaturesA'][-1], self.data['temperaturesB'][-1], self.data['chamber_pressures'][-1]),
             horizontalalignment='center', verticalalignment='bottom', transform=ax2.transAxes)
        ax2.text(0.5, 0.0, 'Temp A std: %.3f Temp B std: %.3f Pressure std: %.3e' % (
        np.std(self.data['temperaturesA'][-10:-1]), np.std(self.data['temperaturesB'][-10:-1]), np.std(self.data['chamber_pressures'][-10:-1])),
                 horizontalalignment='center', verticalalignment='bottom', transform=ax2.transAxes)

    def get_axes_layout(self, figure_list):
        """
        returns the axes objects the script needs to plot its data
        the default creates a single axes object on each figure
        This can/should be overwritten in a child script if more axes objects are needed
        Args:
            figure_list: a list of figure objects
        Returns:
            axes_list: a list of axes objects

        """

        axes_list = []
        if self._plot_refresh is True:
            fig = figure_list[0]
            fig.clf()
            axes_list.append(fig.add_subplot(311))
            axes_list.append(fig.add_subplot(312))
            axes_list.append(fig.add_subplot(313))

        else:
            for fig in figure_list:
                axes_list.append(fig.axes[0])

        return axes_list


    def _plot(self, axes_list):
        '''
        Args:
            axes_list: list of axes objects on which to plot the keyseight spectrun on the first axes object
            data: data (dictionary that contains keys amplitudes, frequencies) if not provided use self.data
        '''

        time = np.array(self.data['time'])
        if max(time) > 3600:
            time /= 3600
            time_label = 'time (h)'
        elif max(time) > 60:
            time /= 60
            time_label = 'time (min)'
        else:
            time_label = 'time (s)'

        lns1 = axes_list[0].plot(time, self.data['chamber_pressures'], color='black',label='Chamber pressure')
        lns2 = axes_list[1].plot(time, self.data['temperaturesA'], color ='red',label='Temp A')
        lns3 = axes_list[2].plot(time, self.data['temperaturesB'], color = 'blue',label='Temp B')

        axes_list[0].set_yscale('log')

        axes_list[0].set_ylabel('Pressure (Torr)')
        axes_list[1].set_ylabel('Temperature A (K)')
        axes_list[2].set_ylabel('Temperature B (K)')
        axes_list[-1].set_xlabel(time_label)

        lns = lns1 + lns2 + lns3
        labs = [l.get_label() for l in lns]

        axes_list[0].set_title('Pressure: %.3e Torr; Temp A: %.3f ${{\pm}}$%.3f K; Temp B: %.3f ${{\pm}}$%.3f K'%
                               (self.data['chamber_pressures'][-1],
                                self.data['temperaturesA'][-1], np.std(self.data['temperaturesA'][-10:-1]),
                                self.data['temperaturesB'][-1], np.std(self.data['temperaturesB'][-10:-1])))
