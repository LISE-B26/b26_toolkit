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

from b26_toolkit.instruments import ZIHF2
from pylabcontrol.core import Script, Parameter
import time

class SetZIOutput(Script):
    """

This script sets the ZI frequency, channel, amplitude, and offset, as well as turns the output on.

ER 20181120

    """

    _DEFAULT_SETTINGS = [
        Parameter('channel', 0, [0, 1], 'signal output channel'),
        Parameter('freq', 1e6, float, 'frequency of output channel'),
        Parameter('amp', 0.1, float, 'amplitude of output channel (V)'),
        Parameter('aux',
                  [
                      Parameter('channel', 0, [0, 1], 'auxilary channel'),
                      Parameter('offset', 1.0, float, 'offset in volts')
                  ]
                  )
    ]

    _INSTRUMENTS = {'ZI': ZIHF2}

    _SCRIPTS = {}


    def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_function= log_function, data_path = data_path)


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        # setup the ZI
        instr_stngs = self.instruments['ZI']['instance'].settings
        instr_stngs['sigouts']['on'] = False
        instr_stngs['sigouts']['range'] = 10.
        instr_stngs['sigouts']['add'] = True
        instr_stngs['aux']['channel'] = self.settings['aux']['channel']
        instr_stngs['sigouts']['channel'] = self.settings['channel']

        self.instruments['ZI']['instance'].update(instr_stngs)


        # check that the minimum voltage is positive
        vmin = self.settings['aux']['offset'] - self.settings['amp']
        if vmin < 0:
            self.log('need to make sure the minimum voltage is > 0 V. Fix parameters please!')
            return

        vmax = self.settings['aux']['offset'] + self.settings['amp']

        range = 0.01

        if vmax > 0.01 and vmax < 0.1:
            range = 0.1
        elif vmax > 0.1 and vmax < 1.0:
            range = 1.0
        elif vmax > 1.0 and vmax < 10.0:
            range = 10.0
        else:
            self.log('need to make sure the maximum voltage is < 10 V. Fix parameters please!')

        # if OK, let's go ahead and update the parameters on the ZI

        instr_stngs['freq'] = self.settings['freq']
        instr_stngs['amp'] = self.settings['amp']
        instr_stngs['sigouts']['range'] = range
        instr_stngs['aux']['offset'] = self.settings['aux']['offset']
        instr_stngs['sigouts']['on'] = True

        self.instruments['ZI']['instance'].update(instr_stngs)

        self.log('updated ZI instrument settings.')
        self.log('frequency set to {:0.3f} MHz'.format(self.settings['freq']*1e-6))
        self.log('amplitude set to {:0.3f} V'.format(self.settings['amp']))
        self.log('offset set to {:0.3f} V'.format(self.settings['aux']['offset']))

        time.sleep(0.5)


