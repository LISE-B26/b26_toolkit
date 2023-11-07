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
from instruments.oscilloscope import rigol_Oscilloscope
from pylabcontrol.data_processing.signal_processing import power_spectral_density


class OScope_Timetrace(Script):
    """
    This script reads the Counter input from the DAQ for a give duration and plots the time trace.
        """
    _DEFAULT_SETTINGS = [
        Parameter('integration_time', 1.0, float, 'Total time per trace (s)'),
        Parameter('num_points', '1M', ['AUTO', '1k', '10k', '100k', '1M', '5M', '10M', '25M', '50M', '100M', '200M'],
                  'Number of points in each trace'),
        Parameter('acquisition_channel', 'CHAN1', ['CHAN1', 'CHAN2', 'CHAN3', 'CHAN4'], 'Daq channel used for counter'),
        Parameter('num_traces', 1, int, 'Number of traces that should be added together in Fourier space'),
        # TODO: delete these as they were daq specific
        # Parameter('ai_channel', 'ai0', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4'], 'Daq channel used for analog in'),
        # Parameter('counter_only', True, bool, 'use counter input only')
    ]

    _INSTRUMENTS = {'oscope': rigol_Oscilloscope}

    _SCRIPTS = {

    }

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.data = {'psd': [], 'freq': [], 'voltage': []}

    # TODO: might not need this task setup formally, unlike the daq
    def setup_scope_tasks(self, number_of_samples):
        """
        Set up the tasks and return a list of tasks

        be careful to return the right order of tasks!

        Args:
            number_of_samples:

        Returns: list of daq, daq_task pairs

        """

        oscope = self.instruments['scope']['instance']
        # initialize APD thread
        # scope_task = oscope.setup_scope(
        #     self.settings['counter_channel'], number_of_samples + 1,
        #     continuous_acquisition=False
        # )

        # daq_ai = self.instruments['daq_ai']['instance']
        # aitask = daq_ai.setup_AI(self.settings['ai_channel'], number_of_samples,
        #                               continuous=False, # continuous sampling still reads every clock tick, here set to the clock of the counter
        #                               clk_source=ctrtask)
        scope_task = oscope.setup_scope()
        return [oscope, scope_task]

    def setup_scope(self):
        """
        Here we set up the scope,

        Returns:

        """

        scope = self.instruments['oscope']['instance']

        integration_time = self.settings['integration_time']
        num_points = self.settings['num_points']
        acquisition_channel = self.settings['acquisition_channel']

        scope.settings['waveform']['timebase'] = integration_time / 10.0
        scope.settings['acquisition']['memory_depth'] = num_points
        scope.settings['waveform']['channel'] = acquisition_channel
        scope.update(scope.settings)

    def read_scope_data(self):

        data = {}
        # for daq, task in reversed(scope_tasks):
        #     if daq is not None:
        # scope, task = scope_tasks
        # task_data, samples_per_channel_read = scope.read(task)

        num_traces = self.settings['num_traces']

        trace_data = []
        for i in range(num_traces):
            # print(self.instruments['oscope']['instance'].query(':BUS1:RS232:BAUD?'))
            scope_data, preambleBlock = self.instruments['oscope']['instance'].get_timetrace()
            trace_data.append((scope_data, preambleBlock))

        # integration_time = self.settings['integration_time']

        # TODO: can we assume the freqs are all the same?
        psd_outputs = []
        for scope_data, preambleBlock in trace_data:
            freq, psd = power_spectral_density(scope_data, preambleBlock['dt'])
            psd_outputs.append(psd)

        psds_sum = np.sum(psd_outputs, axis=0)

        data['voltage'] = scope_data
        data['psd'] = psds_sum
        data['freq'] = freq

        return data
    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        # sample_rate = float(1) / self.settings['integration_time']
        self.setup_scope()

        # # maximum number of samples if total_int_time > 0
        # if self.settings['acquisition_time'] > 0:
        #     number_of_samples = int(np.floor(self.settings['acquisition_time'] / self.settings['integration_time']))
        # else:
        #     self.log('total measurement time must be positive. Abort script')
        #     return

        self.progress = 50
        self.updateProgress.emit(self.progress)

        # scope_tasks = self.setup_scope_tasks(number_of_samples)
        #
        # for scope, task in scope_tasks:
        #     # start task
        #     if scope is not None:
        #         scope.run(task)

        self.data = self.read_scope_data()

        # # clean up
        # for scope, task in scope_tasks:
        #     # start task
        #     if scope is not None:
        #         scope.stop(task)

    def plot(self, figure_list):
        super(OScope_Timetrace, self).plot(figure_list)

    def _plot(self, axes_list, data=None):
        # COMMENT_ME

        if data is None:
            data = self.data

        # This is how you could plot the freq and psd from a signal it seems, but instead we already have these
        # plotting_data = data['data']
        #
        # for signal in [plotting_data]: #ER 20190130
        #     if len(signal) > 0:
        #         # 20191105 ER get rid of normalization
        #         #plot_counts(axes_list[0], signal/np.mean(signal))
        #         plot_counts(axes_list[0], signal)
        #         #freq, psd = power_spectral_density(signal/np.mean(signal), self.settings['integration_time'])
        #         freq, psd = power_spectral_density(signal, self.settings['integration_time'])
        #
        #        # plot_psd(freq, psd, axes_list[1], y_scaling='log', x_scaling='log')
        #         plot_psd(freq[1:], psd[1:], axes_list[1], y_scaling='log', x_scaling='lin') # remove dc component ER 20190129

        plotting_psd = data['psd']
        plotting_freq = data['freq']
        if np.any(data['psd']):
            plot_psd(plotting_freq[1:], plotting_psd[1:], axes_list[1], y_scaling='log', x_scaling='lin')  # remove dc component ER 20190129

    # def _update_plot(self, axes_list):
    #     if self.data:
    #         axis_timetrace, axis_fft = axes_list
    #
    #         update_counts(axes_list[0], self.data['counts'])