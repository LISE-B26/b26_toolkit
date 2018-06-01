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


from b26_toolkit.instruments import NI6259,PiezoController
from pylabcontrol.core import Parameter, Script
from b26_toolkit.plotting.plots_1d import plot_1d_simple_timetrace_ns, update_1d_simple

# JG: NOT FINISHED JUST COPIED CODE AND DIDN"T ADAPT IT YET
class SimplePiezoSweep(Script):
    """
SimplePiezoSweep: Reads analog input (e.g. from photodiode) at different piezo voltages
    """

    _DEFAULT_SETTINGS = [
        Parameter('voltages',
                  [Parameter('min', 0., float, 'min voltage'),
                   Parameter('max', 100, float, 'max voltage'),
                   Parameter('N', 25, int, 'voltage steps')
                  ]),
        Parameter('dt', 0.1, float, 'time between voltage steps (seconds)'),
        Parameter('ai_channel', 'ai0', ['ai0', 'ai1', 'ai2', 'ai3'], 'Daq channel used for voltage analog input')
    ]

    _SCRIPTS = {
    }

    _INSTRUMENTS = {
        'daq':NI6259,
        'z_piezo': PiezoController
    }
    def __init__(self, instruments, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings, instruments, scripts, log_function= log_function, data_path = data_path)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """
        self.data['signal'] = []
        self.data['voltages'] = []

        dt = self.settings['dt']
        sweep_voltages = np.linspace(self.settings['voltages']['min'], self.settings['voltages']['max'], self.settings['voltages']['N'])

        for index, voltage in enumerate(sweep_voltages):
            self._step_piezo(voltage, dt)
            signal = self._get_voltage()
            self.data['signal'].append(signal)

            self.data['voltages'].append(voltage)
            self.progress = 100. * index / len(sweep_voltages)
            self.updateProgress.emit(self.progress if self.progress < 100 else 99)



    def _step_piezo(self, voltage, wait_time):
        """
        steps the piezo.  Has to be overwritten specifically for each different hardware realization
        voltage: target piezo voltage
        wait_time: settle time after voltage step
        """
        z_piezo = self.instruments['z_piezo']['instance']
        # set the voltage on the piezo
        z_piezo.voltage = float(voltage)
        time.sleep(wait_time)

    def _get_voltage(self):
        """
        returns the current position of the galvo
        Returns: list with two floats, which give the x and y position of the galvo mirror

        """
        voltage = self.instruments['daq']['instance'].get_analog_voltages([
            self.settings['ai_channel']
        ])

        return voltage


    def _plot(self, axes_list, data=None):
        # fit the data and set piezo to focus spot
        if data is None:
            data = self.data

        axis_focus, axis_image = axes_list

        axis_focus.plot(data['voltages'],data['signal'], 'o-')

        axis_focus.hold(False)

        axis_focus.set_xlabel('Piezo Voltage [V]')
        axis_focus.set_ylabel('Detector Voltage [V]')

class DAQ_Timetrace(Script):
    """
SimplePiezoSweep: Reads analog input (e.g. from photodiode) at different piezo voltages
    """

    _DEFAULT_SETTINGS = [
        Parameter('sample_rate', 1000, float, 'sample acquisition rate (Hz)'),
        Parameter('ai_channel', 'ai2', ['ai0', 'ai1', 'ai2', 'ai3'], 'Daq channel used for voltage analog input'),
        Parameter('acquisition_time', 10, float, 'time to acquire (s)')
    ]

    _SCRIPTS = {
    }

    _INSTRUMENTS = {
        'daq':NI6259
    }
    def __init__(self, instruments, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings, instruments, scripts, log_function= log_function, data_path = data_path)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """
        self.data['times'] = np.linspace(0, self.settings['acquisition_time']*1e9, self.settings['acquisition_time']*self.settings['sample_rate'])
        self.data['voltages'] = []

        sample_rate = self.settings['sample_rate']
        channel = self.settings['ai_channel']
        self.instruments['daq']['instance'].settings['analog_input'][channel]['sample_rate'] = sample_rate

        task = self.instruments['daq']['instance'].setup_AI(channel, 1000, continuous = True)
        self.instruments['daq']['instance'].run(task)

        start_time = time.time()
        while ((time.time() - start_time) < self.settings['acquisition_time']):
            if self._abort:
                break

            raw_data, num_read = self.instruments['daq']['instance'].read(task)
            self.data['voltages'] = self.data['voltages'] + list(raw_data)

            self.progress = 50.
            self.updateProgress.emit(int(self.progress))

        # clean up APD tasks
        self.instruments['daq']['instance'].stop(task)


    def _plot(self, axes_list, data = None):
        # COMMENT_ME

        if data is None:
            data = self.data

        axes_list[0].hold(False)

        if data:
            # plot_1d_simple_timetrace_ns(axes_list[0], self.data['times'][0:len(data['voltages'])], np.array(data['voltages']))
            axes_list[0].plot(data['voltages'])