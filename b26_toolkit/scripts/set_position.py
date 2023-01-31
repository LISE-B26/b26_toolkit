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

from b26_toolkit.instruments import NI6259, NI9263, NI6229
from pylabcontrol.core import Script, Parameter
from b26_toolkit.scripts import AFMScan


class SetPosition(AFMScan):
    """
This script positions the sample to a point
    """
    _DEFAULT_SETTINGS = [
        Parameter('point_a',
                  [Parameter('x', 0, float, 'x-coordinate'),
                   Parameter('y', 0, float, 'y-coordinate')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 75.0, [75.0], 'x-coordinate'),
                   Parameter('y', 75.0, [75.0], 'y-coordinate')
                   ]),
        Parameter('speed_limit [V/s]',1,float,'speed at which the probe is moved to start'),
        Parameter('settle_time', .0005,[.0005], 'wait time between points to allow galvo to settle'),
        Parameter('time_per_pt', .001,[.001], 'time in s to measure at each point'),
        Parameter('num_points',
                  [Parameter('x', 1, [1], 'number of x points to scan'),
                   Parameter('y', 1, [1], 'number of y points to scan')
                   ]),
        Parameter('RoI_mode', 'corner', ['corner'], 'mode to calculate region of interest.\n \
                                                           corner: pta and ptb are diagonal corners of rectangle.\n \
                                                          center: pta is center and pta is extend or rectangle'),
        Parameter('ending_behavior', 'return_to_a', ['return_to_a'], 'return to the corn'),
        Parameter('DAQ_channels',
                  [Parameter('x_ao_channel', 'ao2', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
                   Parameter('y_ao_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
                   Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'], 'Daq channel counter to use as a clock for ao and ai'),
                   Parameter('ai_channel', 'ai0', ['ai0', 'ai1', 'ai2', 'ai3'], 'Daq channel used for photodiode voltage')
                   ]),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]


    # _INSTRUMENTS = {'NI6259':  NI6259, 'NI9263': NI9263, 'NI9402': NI9402} when pulse_blaster_base_script uses 'daq' as the key for daq, we can't use these settings
    _INSTRUMENTS = {'daq': NI6229} # change daq names for when using this script on computers that aren't Alice

    _SCRIPTS = {}

    def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        #super(SetPosition,self).__init__(instruments,name=None,settings=None,log_function=None,data_path=None)
        Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_function= log_function, data_path = data_path)

        # set some dummy variables to enable GalvoScan
        #self.settings['point_b']={'x':75,'y':75}
        #self.settings['settle_time']=0.0005


    def after_scan(self):
        print("after scan")
        pass

    def read_line(self, y_pos):
        print("read line, y_pos: %f" %y_pos)
        pass

    def _plot(self, axes_list, data = None):
        print("do not plot anything")
        pass

    def _update_plot(self, axes_list):
        pass


if __name__ == '__main__':
    from pylabcontrol.core import Instrument

    # instruments, instruments_failed = Instrument.load_and_append({'daq':  'NI6259'})

    script, failed, instruments = Script.load_and_append(script_dict={'SetPosition': 'SetPosition'})

    print(script)
    print(failed)
    # print(instruments)