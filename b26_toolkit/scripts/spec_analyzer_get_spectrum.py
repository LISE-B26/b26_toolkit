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

from b26_toolkit.instruments import SpectrumAnalyzer
    # , MicrowaveGenerator, CryoStation
from pylabcontrol.core import Script


class SpecAnalyzerGetSpectrum(Script):
    # COMMENT_ME

    _DEFAULT_SETTINGS = [
        # Parameter('start_frequency', 2.7e9, float, 'start frequency of spectrum'),
        # Parameter('stop_frequency', 3e9, float, 'end frequency of spectrum'),
        # Parameter('output_power',0.0, float, 'output power (dBm)'),
        # Parameter('output_on',True, bool, 'enable output'),
    ]

    _INSTRUMENTS = {
        'spectrum_analyzer' : SpectrumAnalyzer
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



        instrument = self.instruments['spectrum_analyzer']['instance']
        settings = self.instruments['spectrum_analyzer']['settings']

        instrument.update(settings)
        trace = instrument.trace

        self.data = trace
        print('acquired spectrume')



    def _plot(self, axes_list, data = None):
        '''
        Plots the galvo scan image
        Args:
            axes_list: list of axes objects on which to plot the keyseight spectrun on the first axes object
            data: data (dictionary that contains keys amplitudes, frequencies) if not provided use self.data
        '''
        if data is None:
            data = self.data

        axes_list[0].plot(data['frequencies'], data['amplitudes'])
        axes_list[0].set_xlabel('frequencies')
        axes_list[0].set_xlabel('spectrum (??)')