"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    pylabcontrol is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    pylabcontrol is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with pylabcontrol.  If not, see <http://www.gnu.org/licenses/>.
"""

from b26_toolkit.instruments import Oscilloscope
from pylabcontrol.core import Script
import numpy as np

from pylabcontrol.data_processing.signal_processing import power_spectral_density

class KeysightOsciGetTimeTrace(Script):
    # COMMENT_ME

    _DEFAULT_SETTINGS = [
        # Parameter('start_frequency', 2.7e9, float, 'start frequency of spectrum'),
        # Parameter('stop_frequency', 3e9, float, 'end frequency of spectrum'),
        # Parameter('output_power',0.0, float, 'output power (dBm)'),
        # Parameter('output_on',True, bool, 'enable output'),
    ]

    _INSTRUMENTS = {
        'osci' : Oscilloscope
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



        instrument = self.instruments['osci']['instance']
        settings = self.instruments['osci']['settings']

        instrument.reset()

        instrument.update(settings)
        trace, preamble = instrument.get_timetrace()

        self.data = {'voltage': trace, 'meta_data': preamble}
        print('acquired spectrum')



    def _plot(self, axes_list, data = None):
        '''
        Plots the galvo scan image
        Args:
            axes_list: list of axes objects on which to plot the keyseight spectrun on the first axes object
            data: data (dictionary that contains keys amplitudes, frequencies) if not provided use self.data
        '''
        if data is None:
            data = self.data

        dt = self.data['meta_data']['xincrement']
        data = self.data['voltage']
        time = dt*np.arange(len(data))
        axes_list[0].plot(time, data, '-')
        axes_list[0].set_xlabel('time (s)')
        axes_list[0].set_ylabel('signal (arb.)')


        F, P = power_spectral_density(data, dt)

        print(('JG adasd', data, dt))

        axes_list[1].plot(F, P, '-')
        axes_list[1].set_xlabel('freq (Hz)')
        axes_list[1].set_ylabel('signal (arb.)')

        axes_list[1].set_xscale("log")
        # JG: try to display on a log scale, this doesn't work if the psd is negative or zero (which might happen if the oscilloscope is out of range)
        if np.mean(data) > 0:
            axes_list[1].set_yscale("log")




if __name__ == '__main__':
    from pylabcontrol.core import Instrument
    # from b26_toolkit.pylabcontrol.instruments import NI7845RMain
    #
    # fpga = NI7845RMain()
    #
    #
    # g = GalvoScanFPGA(instruments={'NI7845RMain':fpga}, name='test_fpga_scan', settings=None, log_function=None, data_path=None)
    # print(fpga)


    # instruments, failed =  Instrument.load_and_append(instrument_dict ={'NI7845RMain': 'NI7845RMain'}, raise_errors=True )

    script, failed, instruments = Script.load_and_append(script_dict={'GalvoScanFPGA': 'GalvoScanFPGA'}, raise_errors=True)
    #
    print(script)
    print(failed)
    # # print(instruments)

