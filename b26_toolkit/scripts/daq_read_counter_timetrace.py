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
from collections import deque
import numpy as np
import time
from copy import deepcopy

from b26_toolkit.instruments import NI6259, NI9402, NI9219
from b26_toolkit.plotting.plots_1d import plot_counts, update_counts,  plot_psd
from pylabcontrol.core import Parameter, Script
from pylabcontrol.data_processing.signal_processing import power_spectral_density


class DaqReadCounterTimeTrace(Script):
    """
    This script reads the Counter input from the DAQ for a given duration and plots the time trace. This is a more accurately timed version of DaqReadCounter.
    DaqReadCounter continuously runs a DAQ task over and over again, so there's some timing inaccuracy from the DAQ overhead and computer/DAQ clock differences.
    This script runs a single DAQ task and waits until it finishes, so the timing is entirely handled by the DAQ clock and should be reliable.
    Downside is there's no realtime info on counts.
    """
    _DEFAULT_SETTINGS = [
        Parameter('integration_time', .25, float, 'Time per data point (s)'),
        Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1'], 'Daq channel used for counter'),
        Parameter('total_int_time', 3.0, float, 'Total time to integrate (s)'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for counting')
    ]

    _INSTRUMENTS = {'NI6259':  NI6259, 'NI9402': NI9402}

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

        self.data = {'counts': []}


    def get_daq_task(self, number_of_samples):

        # initialize APD thread
        ctrtask = self.daq.setup_counter(
            self.settings['counter_channel'], number_of_samples + 1)
        return ctrtask

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        if self.settings['daq_type'] == 'PCI':
            self.daq = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq = self.instruments['NI9402']['instance']


        sample_rate = float(1) / self.settings['integration_time']
        # normalization = self.settings['integration_time']/.001
        self.daq.settings['digital_input'][self.settings['counter_channel']]['sample_rate'] = sample_rate

        # maximum number of samples if total_int_time > 0
        if self.settings['total_int_time'] > 0:
            number_of_samples = int(np.floor(self.settings['total_int_time']/self.settings['integration_time']))
        else:
            self.log('total measurement time must be positive. Abort script')
            return

        self.progress = 0
        self.updateProgress.emit(self.progress)

        daq_task = self.get_daq_task(number_of_samples)
        # start counter and scanning sequence
        self.daq.run(daq_task)

        if self.settings['total_int_time']/100 < 1:
            time_divisions = int(self.settings['total_int_time'] / 1)
            if time_divisions == 0:
                time_divisions = 1
        else:
            time_divisions = 100

        print(time_divisions)
        for t in range(int(time_divisions)):
            print(self.settings['total_int_time']/time_divisions)
            time.sleep(self.settings['total_int_time']/time_divisions)
            self.progress = float(t)/time_divisions*100
            self.updateProgress.emit(self.progress)

        data, _ = self.daq.read(daq_task)

        self.daq.stop(daq_task)
        counts = np.diff(data)  # counter gives the accumulated counts, thus the diff gives the counts per interval

        self.data['counts'] = counts * sample_rate/1000  # multiply by the sample rate to get kcounts /second

    def plot(self, figure_list):
        super(DaqReadCounterTimeTrace, self).plot(figure_list)

    def _plot(self, axes_list, data = None):
        # COMMENT_ME

        if data is None:
            data = self.data


        if len(data['counts'])>0:
            plot_counts(axes_list[0], data['counts'])
            freq, psd = power_spectral_density(data['counts'], self.settings['integration_time'])
            plot_psd(freq, psd, axes_list[1], y_scaling='log', x_scaling='log')

    def _update_plot(self, axes_list):
        pass
        if self.data and 1 != 1:
            axis_timetrace, axis_fft = axes_list

            update_counts(axes_list[0], self.data['counts'])


# work in progress (JG) write generic timetrace class
class DaqTimeTraceGeneric(Script):
    """
This script reads the Counter input from the DAQ for a give duration and plots the time trace.
Future: plot also the PSD
    """
    _DEFAULT_SETTINGS = [
        Parameter('integration_time', .25, float, 'Time per data point (s)'),
        Parameter('acquisition_time', 3.0, float, 'Total acquisition time (s)')
    ]

    _INSTRUMENTS = {}

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

        self.data = {'counts':[]}


    def get_daq_task(self, number_of_samples):
        """
         setup the daq taks using the daq instances created with setup_daq()

         overwrite this function

        Args:
            number_of_samples:

        Returns:

        """

        # overwite this, create the daq count
        raise NotImplementedError
        return ctrtask

    def setup_daq_ai(self, sample_rate, daq_ai, ai_channel):
        """

        generic function to setup ai task

        Args:
            sample_rate:
            daq_ai:
            ai_channel:

        Returns:

        """

    def setup_daq_counter(self, sample_rate, daq_counter, counter_channel):
        """

        generic function to setup a counter task

        Args:
            sample_rate:
            daq_counter:
            counter_channel:

        Returns:

        """
        # overwite this, create the daq instance references

        # # normalization = self.settings['integration_time']/.001
        daq_counter.settings['digital_input'][counter_channel]['sample_rate'] = sample_rate

    def setup_daq(self):
        """
        Here we set up the daqs that we need,
        use the functions setup_daq_counter() and setup_daq_ai() from the generic class
        to implement digital and analog channels


        Returns:

        """
        raise NotImplementedError
        setup_daq_counter(self, sample_rate, daq_counter, counter_channel)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        sample_rate = float(1) / self.settings['integration_time']


        self.setup_daq()  # setup the daq


        # maximum number of samples if total_int_time > 0
        if self.settings['acquisition_time'] > 0:
            number_of_samples = int(np.floor(self.settings['acquisition_time']/self.settings['integration_time']))
        else:
            self.log('total measurement time must be positive. Abort script')
            return

        self.progress = 50
        self.updateProgress.emit(self.progress)

        daq_task = self.get_daq_task(number_of_samples)


        # start counter and scanning sequence
        self.daq.run(daq_task)

        data, _ = self.daq.read(daq_task)

        self.daq.stop(daq_task)
        counts = np.diff(data)  # counter gives the accumulated counts, thus the diff gives the counts per interval

        self.data['counts'] = counts * sample_rate/1000  # multiply by the sample rate to get kcounts /second

    def plot(self, figure_list):
        super(DaqReadCounterTimeTrace, self).plot(figure_list)

    def _plot(self, axes_list, data = None):
        # COMMENT_ME

        if data is None:
            data = self.data


        if len(data['counts'])>0:
            plot_counts(axes_list[0], data['counts'])
            freq, psd = power_spectral_density(data['counts'], self.settings['integration_time'])
            plot_psd(freq, psd, axes_list[1], y_scaling='log', x_scaling='log')

    def _update_plot(self, axes_list):
        if self.data:
            axis_timetrace, axis_fft = axes_list

            update_counts(axes_list[0], self.data['counts'])


# work in progress (JG) write generic timetrace class
class DaqTimeTraceNi9402Ni9219(Script):
    """
    This script reads the Counter input from the DAQ for a give duration and plots the time trace.
        """
    _DEFAULT_SETTINGS = [
        Parameter('integration_time', .25, float, 'Time per data point (s)'),
        Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1'], 'Daq channel used for counter'),
        Parameter('acquisition_time', 3.0, float, 'Total acquisition time (s)'),
        Parameter('ai_channel', 'ai0', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4'], 'Daq channel used for analog in'),
        Parameter('ctr_or_ai', 'ctr', ['ctr', 'ai'], 'select ctr or ai'),
        Parameter('num_averages', 1, int, 'number of averages for psd')
    ]

    _INSTRUMENTS = {'daq_ai': NI9219, 'daq_counter': NI9402}

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


    def setup_daq_tasks(self, number_of_samples):
        """
        setup the tasks and return a list of tasks

        be carefull to return the right order of tasks!

        Args:
            number_of_samples:

        Returns: list of daq, daq_task pairs

        """

        daq_counter = self.instruments['daq_counter']['instance']
        # initialize APD thread
        ctrtask = daq_counter.setup_counter(
            self.settings['counter_channel'], number_of_samples + 1,
            continuous_acquisition=False
        )

        if self.settings['ctr_or_ai'] == 'ctr':
            return [[daq_counter, ctrtask]]

        daq_ai = self.instruments['daq_ai']['instance']
        aitask = daq_ai.setup_AI(self.settings['ai_channel'], number_of_samples,
                                      continuous=False, # continuous sampling still reads every clock tick, here set to the clock of the counter
                                      clk_source=ctrtask)

        return [[daq_ai, aitask], [daq_counter, ctrtask]]

    def setup_daq(self, sample_rate):
        """
        Here we set up the daqs that we need,

        Returns:

        """
        daq_counter = self.instruments['daq_counter']['instance']

        counter_channel = self.settings['counter_channel']


        daq_counter.settings['digital_input'][counter_channel]['sample_rate'] = sample_rate

    def read_daq_data(self, daq_tasks, sample_rate):

        data = {}
        for daq, task in reversed(daq_tasks):
            if daq is not None:
                task_data, samples_per_channel_read = daq.read(task)

                if daq == self.instruments['daq_counter']['instance']:
                    counts = np.diff(task_data)  # counter gives the accumulated counts, thus the diff gives the counts per interval
                    data['counts'] = counts * sample_rate / 1000  # multiply by the sample rate to get kcounts /second
                elif daq == self.instruments['daq_ai']['instance']:
                    data['ai'] = np.array(task_data)

                else:
                    raise KeyError('unknown daq type in read_daq_data()')

        return data
    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        # if self.settings['daq_type'] == 'PCI':
        #     self.daq = self.instruments['NI6259']['instance']
        # elif self.settings['daq_type'] == 'cDAQ':
        #     self.daq = self.instruments['NI9402']['instance']
        self.data = {'counts': [], 'ai':[]}

        sample_rate = float(1) / self.settings['integration_time']
        self.setup_daq(sample_rate)


        # maximum number of samples if total_int_time > 0
        if self.settings['acquisition_time'] > 0:
            number_of_samples = int(np.floor(self.settings['acquisition_time'] / self.settings['integration_time']))
        else:
            self.log('total measurement time must be positive. Abort script')
            return

        # self.progress = 50
        # self.updateProgress.emit(self.progress)
        self.data['start_time'] = time.time() # in seconds

        for i in range(self.settings['num_averages']):
            daq_tasks = self.setup_daq_tasks(number_of_samples)

            # DS20191016 Record start time
            # print('Start time: {}'.format(time.time()))

            for daq, task in daq_tasks:
                # start task
                if daq is not None:
                    daq.run(task)

            data = self.read_daq_data(daq_tasks, sample_rate)

            if self.settings['ctr_or_ai'] == 'ai':
                self.data['ai'].append(data['ai'])
            else:
                self.data['counts'].append(data['counts'])


            for daq, task in daq_tasks:
                # start task
                if daq is not None:
                    daq.stop(task)

        # clean up
        if self.settings['ctr_or_ai'] == 'ai':
            self.data['ai'] = np.array(self.data['ai']).transpose()
        else:
            self.data['counts'] = np.array(self.data['counts']).transpose()
    def plot(self, figure_list):
        super(DaqTimeTraceNi9402Ni9219, self).plot(figure_list)

    def _plot(self, axes_list, data=None):
        # COMMENT_ME

        if data is None:
            data = self.data

        if self.settings['ctr_or_ai'] == 'ctr':
            plotting_data = np.array(data['counts'])
        else:
            plotting_data = np.array(data['ai'])
        if len(plotting_data) == 0:
            return
        print(plotting_data.shape)
        # plot_counts(axes_list[0], plotting_data[0], int_time=self.settings['integration_time'])
        psdAvg = np.zeros(int(plotting_data.shape[0] / 2 + 1))

        for i in range(plotting_data.shape[1]): #ER 20190130
            if len(plotting_data[:, i]) > 0:
                # 20191105 ER get rid of normalization
                #plot_counts(axes_list[0], signal/np.mean(signal))
                #freq, psd = power_spectral_density(signal/np.mean(signal), self.settings['integration_time'])
                freq, psd = power_spectral_density(plotting_data[:, i], self.settings['integration_time'])
                psdAvg += psd

               # plot_psd(freq, psd, axes_list[1], y_scaling='log', x_scaling='log')
        plot_psd(freq[1:], psdAvg[1:] / plotting_data.shape[1], axes_list[1], y_scaling='log', x_scaling='lin') # remove dc component ER 20190129


    # def _update_plot(self, axes_list):
    #     if self.data:
    #         axis_timetrace, axis_fft = axes_list
    #
    #         update_counts(axes_list[0], self.data['counts'])


class DaqTimeTraceNi6259(DaqTimeTraceNi9402Ni9219):
    _INSTRUMENTS = {'daq_ai': NI6259, 'daq_counter': NI6259}


class DaqTimeTraceNi6259Pd(DaqTimeTraceNi9402Ni9219): # ER 20200731

    '''


    Daq timetrace for NI6259 - plots both the counter APD data and the photodiode AI data.

    '''

    _INSTRUMENTS = {'daq_ai': NI6259, 'daq_counter': NI6259}

    def _plot(self, axes_list, data=None):
        # COMMENT_ME

        if data is None:
            data = self.data

        # for signal in [data['counts']]: #, data['ai']]: ER 20190130
        #     if len(signal) > 0:
        #         # 20191105 ER get rid of normalization
        #         #plot_counts(axes_list[0], signal/np.mean(signal))
        #         plot_counts(axes_list[0], signal)
        #         #freq, psd = power_spectral_density(signal/np.mean(signal), self.settings['integration_time'])
        #         freq, psd = power_spectral_density(signal, self.settings['integration_time'])
        #        # plot_psd(freq, psd, axes_list[1], y_scaling='log', x_scaling='log')
        #         plot_psd(freq[1:], psd[1:], axes_list[1], y_scaling='log', x_scaling='lin') # remove dc component ER 20190129

        for signal in data['ai']:
            if len(signal) > 0:
                # 20191105 ER get rid of normalization
                #plot_counts(axes_list[0], signal/np.mean(signal))
                plot_counts(axes_list[0], signal)
                #freq, psd = power_spectral_density(signal/np.mean(signal), self.settings['integration_time'])
                freq, psd = power_spectral_density(signal, self.settings['integration_time'])
               # plot_psd(freq, psd, axes_list[1], y_scaling='log', x_scaling='log')
                plot_psd(freq[1:], psd[1:], axes_list[1], y_scaling='log', x_scaling='lin', color='r') # remove dc component ER 20190129

# if __name__ == '__main__':
#     script = {}
#     instr = {}
#     script, failed, instr = Script.load_and_append({'Daq_Read_Cntr': 'Daq_Read_Cntr'}, script, instr)
#
#     print(script)
#     print(failed)
#     print(instr)