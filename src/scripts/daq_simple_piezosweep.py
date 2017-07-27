"""
    This file is part of b26_toolkit, a PyLabControl add-on for experiments in Harvard LISE B26.
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


from b26_toolkit.src.instruments import NI6259,PiezoController
from PyLabControl.src.core import Parameter, Script


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
        Parameter('dt', 0.1, float, 'time between voltage steps (seconds)')
    ]

    _SCRIPTS = {
    }

    _INSTRUMENTS = {
        'daq':NI6259,
        'z_piezo': PiezoController
    }
    def __init__(self, scripts, instruments = None, name = None, settings = None, log_function = None, data_path = None):
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

        dt = self.settings['dt']
        sweep_voltages = np.linspace(self.settings['voltages']['min'], self.settings['voltages']['max'], self.settings['voltages']['N'])

        for voltage in sweep_voltages:
            self._step_piezo(voltage, dt)
            signal = self._get_voltage()
            self.data['signal'].append(signal)

        self.data['voltages'] = sweep_voltages



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
            self.settings['ao_channel']
        ])

        return voltage

