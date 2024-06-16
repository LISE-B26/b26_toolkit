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
from pylabcontrol.core import Parameter, Script
from b26_toolkit.plotting.plots_1d import plot_counts, update_counts,  plot_psd
from b26_toolkit.instruments.oscilloscope import RigolOscilloscope
from scipy.signal import periodogram

class OscilloscopeTimetrace(Script):
    """
    This script reads the Counter input from the DAQ for a give duration and plots the time trace.
        """
    _DEFAULT_SETTINGS = [
        Parameter('integration_time', 1.0, float, 'Total time per trace (s)'),
        Parameter('num_points', '1M', ['AUTO', '1k', '10k', '100k', '1M', '5M', '10M', '25M', '50M', '100M', '200M'],
                  'Number of points in each trace'),
        Parameter('offset', 0.0, float, 'voltage offset [V]'),
        Parameter('vert_scale', 0.001, float, 'voltage scale [V]'),
        Parameter('acquisition_channel', 'CHAN1', ['CHAN1', 'CHAN2', 'CHAN3', 'CHAN4'], 'Daq channel used for counter'),
        Parameter('num_traces', 1, int, 'Number of traces that should be added together in Fourier space'),
    ]

    _INSTRUMENTS = {'oscope': RigolOscilloscope}

    _SCRIPTS = {}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.oscope = self.instruments['oscope']['instance']

    def setup_scope(self):
        """
        Here we set up the scope,

        Returns:

        """

        integration_time = self.settings['integration_time']
        num_points = self.settings['num_points']
        acquisition_channel = self.settings['acquisition_channel']

        self.oscope.settings['acq_memory_depth'] = num_points
        self.oscope.update(self.oscope.settings)

        self.oscope.settings['timebase_format'] = 'MAIN'
        self.oscope.settings['timebase'] = integration_time / 10.0
        self.oscope.settings['channel'] = acquisition_channel
        self.oscope.update(self.oscope.settings)
        self.oscope.settings['vert_scale'] = self.settings['vert_scale']
        self.oscope.update(self.oscope.settings)
        self.oscope.settings['offset'] = self.settings['offset']
        self.oscope.update(self.oscope.settings)

    def read_scope_data(self):

        num_traces = self.settings['num_traces']

        trace_data = []
        for i in range(num_traces):
            trace_data.append(self.oscope.get_timetrace())

        return trace_data

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        self.setup_scope()

        self.progress = 50
        self.updateProgress.emit(self.progress)

        self.data['voltage'] = self.read_scope_data()

    def _plot(self, axes_list, data=None):
        # COMMENT_ME

        if data is None:
            data = self.data

        if len(data['voltage']) > 0:
            psd_sum = np.zeros(int(len(data['voltage'][0]) / 2) + 1)
            for scope_data in data['voltage']:
                if len(scope_data) > 0:
                    freq, psd = periodogram(scope_data, fs=(len(scope_data) / (self.settings['integration_time'])))
                    psd_sum += psd

            psd_sum /= len(data['voltage'])

            if np.any(data['voltage']):
                plot_counts(axes_list[0], data['voltage'][0], int_time=self.settings['integration_time'] / len(data[0]))
                plot_psd(freq[1:], psd_sum[1:], axes_list[1], y_scaling='log', x_scaling='lin')  # remove dc component ER 20190129

    # def _update_plot(self, axes_list):
    #     if self.data:
    #         axis_timetrace, axis_fft = axes_list
    #
    #         update_counts(axes_list[0], self.data['counts'])