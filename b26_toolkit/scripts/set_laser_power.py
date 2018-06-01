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

from b26_toolkit.instruments import NI6259, NI9263
from pylabcontrol.core import Script, Parameter


class SetLaserPower(Script):
    """
This script sets the laser power by changing the analog voltage on the fiber coupled attenuator. ER 20180308

Spec sheet at https://www.thorlabs.com/drawings/3df25d052232c1a0-4F4D2098-BE71-279F-FAD6F90E43F9341F/V450A-SpecSheet.pdf
V450A
    """

    _DEFAULT_SETTINGS = [
        Parameter('attenuation', 0.0, float, 'attenuation factor, bewteen 1-1000'),
        Parameter('DAQ_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for voltage analog output'),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    _INSTRUMENTS = {'NI6259_laser_power':  NI6259, 'NI9263_laser_power': NI9263}

    _SCRIPTS = {}


    def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_function= log_function, data_path = data_path)
        if self.settings['daq_type'] == 'PCI':
            self.daq_out = self.instruments['NI6259_laser_power']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_out = self.instruments['NI9263_laser_power']['instance']

    def _function(self):

        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        # parameters / specs of the fiber couple attenuator
        min_volts = 0.0
        max_volts = 5.0

        pt = self.get_voltage(self.settings['attenuation'])

        if pt < min_volts or pt > max_volts:
            raise ValueError("Invalid voltage. Must be between 0 and 5 volts")

        task = self.daq_out.setup_AO(self.settings['DAQ_channels'], pt)
        self.daq_out.run(task)
        self.daq_out.waitToFinish(task)
        self.daq_out.stop(task)
        self.log('laser set to attenuation'.format(self.settings['attenuation']))

    def get_voltage(self, attenuation):
        # this function returns the voltage needed for a given attenuation
        # fit to a quartic polynomial

        voltage = a4*attenuation^4 + a3*attenuation^3 + a2*attenuation^2 + a1*attenuation + a0
        return voltage

if __name__ == '__main__':
    from pylabcontrol.core import Instrument

    # instruments, instruments_failed = Instrument.load_and_append({'daq':  'NI6259'})

    script, failed, instruments = Script.load_and_append(script_dict={'SetLaserPower': 'SetLaserPower'})

    print(script)
    print(failed)
    # print(instruments)