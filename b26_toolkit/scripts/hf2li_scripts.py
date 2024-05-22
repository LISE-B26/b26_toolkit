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

from collections import deque
import numpy as np
import time

import zhinst.utils

from b26_toolkit.plotting.plots_1d import plot_psd
from pylabcontrol.core import Script, Parameter
from b26_toolkit.instruments import Hf2Li
from b26_toolkit.plotting.plots_1d import plot_1d_simple_timetrace
from scipy.signal import periodogram as periodogram

class LockInDaqRead(Script):
    """
This script performs a frequency sweep with the Zurich Instrument HF2 Series Lock-in amplifier
    """
    _DEFAULT_SETTINGS = [
        Parameter('demods',
                  [Parameter('0',
                             [Parameter('filter_order', 4, [1, 2, 3, 4, 5, 6, 7, 8], 'low pass filter order'),
                              Parameter('bandwidth', 10., float, 'low pass filter bandwidth (Hz)'),
                              Parameter('data',
                                        [Parameter('x', False, bool, 'collect data on x quadrature for chosen demodulators'),
                                         Parameter('y', False, bool, 'collect data on y quadrature for chosen demodulators'),
                                         Parameter('freq', False, bool, 'collect data on frequency for chosen demodulators')])
                              ]),
                   Parameter('1',
                             [Parameter('filter_order', 4, [1, 2, 3, 4, 5, 6, 7, 8], 'low pass filter order'),
                              Parameter('bandwidth', 10., float, 'low pass filter bandwidth (Hz)'),
                              Parameter('data',
                                        [Parameter('x', False, bool, 'collect data on x quadrature for chosen demodulators'),
                                         Parameter('y', False, bool, 'collect data on y quadrature for chosen demodulators'),
                                         Parameter('freq', False, bool, 'collect data on frequency for chosen demodulators')])
                              ]),
                   Parameter('2',
                             [Parameter('filter_order', 4, [1, 2, 3, 4, 5, 6, 7, 8], 'low pass filter order'),
                              Parameter('bandwidth', 10., float, 'low pass filter bandwidth (Hz)'),
                              Parameter('data',
                                        [Parameter('x', False, bool, 'collect data on x quadrature for chosen demodulators'),
                                         Parameter('y', False, bool, 'collect data on y quadrature for chosen demodulators'),
                                         Parameter('freq', False, bool, 'collect data on frequency for chosen demodulators')])
                              ]),
                   Parameter('3',
                             [Parameter('filter_order', 4, [1, 2, 3, 4, 5, 6, 7, 8], 'low pass filter order'),
                              Parameter('bandwidth', 10., float, 'low pass filter bandwidth (Hz)'),
                              Parameter('data',
                                        [Parameter('x', False, bool, 'collect data on x quadrature for chosen demodulators'),
                                         Parameter('y', False, bool, 'collect data on y quadrature for chosen demodulators'),
                                         Parameter('freq', False, bool, 'collect data on frequency for chosen demodulators')])
                              ]),
                   Parameter('4',
                             [Parameter('filter_order', 4, [1, 2, 3, 4, 5, 6, 7, 8], 'low pass filter order'),
                              Parameter('bandwidth', 10., float, 'low pass filter bandwidth (Hz)'),
                              Parameter('data',
                                        [Parameter('x', False, bool, 'collect data on x quadrature for chosen demodulators'),
                                         Parameter('y', False, bool, 'collect data on y quadrature for chosen demodulators'),
                                         Parameter('freq', False, bool, 'collect data on frequency for chosen demodulators')])
                              ]),
                   Parameter('5',
                             [Parameter('filter_order', 4, [1, 2, 3, 4, 5, 6, 7, 8], 'low pass filter order'),
                              Parameter('bandwidth', 10., float, 'low pass filter bandwidth (Hz)'),
                              Parameter('data',
                                        [Parameter('x', False, bool, 'collect data on x quadrature for chosen demodulators'),
                                         Parameter('y', False, bool, 'collect data on y quadrature for chosen demodulators'),
                                         Parameter('freq', False, bool, 'collect data on frequency for chosen demodulators')])
                              ])
                   ]),
        Parameter('sampling_rate', 1000, float, 'number of samples per second (Hz)'),
        Parameter('segment_duration', 1, float, 'duration (s) of one sampling segment'),
        Parameter('segment_num', 1, int, 'number of segments'),
        Parameter('fft', False, bool, 'plot FFT of time series data (using Hanning window); this might take a while for large datasets')
    ]

    _INSTRUMENTS = {'hf2li': Hf2Li}

    _SCRIPTS = {}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path=data_path)

    def configure_lock_in(self):
        self.device = self.instruments['hf2li']['instance'].device
        self.session = self.instruments['hf2li']['instance'].session

        self.daq_module = self.session.modules.daq
        self.daq_module.device(self.device)
        self.daq_module.type(0)  # continuous acquisition
        self.daq_module.grid.mode(4)  # data points will be exactly on-grid
        self.daq_module.count(self.settings['segment_num'])  # number of "runs" of data to collect
        self.daq_module.duration(self.settings['segment_duration'])  # duration of each run
        num_cols = int(np.ceil(self.settings['sampling_rate'] * self.settings['segment_duration']))
        self.daq_module.grid.cols(num_cols)

        # Make a list of samples to collect (e.g. demod 4 x, demod 5 y, etc)
        # Also set sampling rate for all demods, enable demods
        self.sample_nodes = []
        self.data_labels = []

        for demod in self.settings['demods'].keys():
            demod_num = int(demod)
            if self.settings['demods'][demod]['data']['x']:
                self.sample_nodes.append(self.device.demods[demod_num].sample.x)
                self.data_labels.append('demod%i_x' % demod_num)
            if self.settings['demods'][demod]['data']['y']:
                self.sample_nodes.append(self.device.demods[demod_num].sample.y)
                self.data_labels.append('demod%i_y' % demod_num)
            if self.settings['demods'][demod]['data']['freq']:
                self.sample_nodes.append(self.device.demods[demod_num].sample.frequency)
                self.data_labels.append('demod%i_freq' % demod_num)
            if any(self.settings['demods'][demod]['data'].values()):  # Only configure demod if data from it will be collected
                self.device.demods[demod_num].enable(True)  # Enable demods
                self.device.demods[demod_num].trigger(0)  # Continuous trigger
                self.device.demods[demod_num].rate(self.settings['sampling_rate'])  # Set sampling rate for all enabled demods
                order = self.settings['demods'][demod]['filter_order']  # Set demod filter order
                self.device.demods[demod_num].order(order)
                tc = zhinst.utils.bw2tc(self.settings['demods'][demod]['bandwidth'], order)
                self.device.demods[demod_num].timeconstant(tc)  # Set demod time constant based on user-defined bandwidth

        if self.settings['fft']:
            self.data_labels_fft = ['%s_fft' % data_label for data_label in self.data_labels]

        if not self.sample_nodes:
            self.log('Error: No data channels requested! Aborting script.')
            self._abort = True

        for node in self.sample_nodes:
            self.daq_module.subscribe(node)


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        self.configure_lock_in()  # Configure lock-in amplifier settings (e.g. demod filter BWs, sampling freq)
        # Initialize fields in self.data
        self.data = {'time': [], 'sampling_rate_actual': [], 'bandwidth_actual': []}
        for label in self.data_labels:
            self.data[label] = []
        if self.settings['fft']:
            for label in self.data_labels_fft:
                self.data[label] = []

        # Get actual device settings after configuration; might differ from input settings b/c of rounding
        for demod in self.settings['demods'].keys():
            demod_num = int(demod)
            self.data['sampling_rate_actual'].append(self.device.demods[demod_num].rate())  # Get sampling rate
            order = self.device.demods[demod_num].order()  # Get demod filter order
            self.data['bandwidth_actual'].append(zhinst.utils.bw2tc(self.device.demods[demod_num].timeconstant(), order))

        # Run daq
        self.daq_module.execute()

        while True:
            if self.daq_module.raw_module.finished():
                self.log('DAQ has finished collecting data, now processing...')
                break
            elif self._abort:
                self.log('Script aborted.')
                break

        # Read from daq
        daq_data = self.daq_module.read(raw=False, clk_rate=self.instruments['hf2li']['instance'].clockbase)
        for node, data_label in zip(self.sample_nodes, self.data_labels):
            if node in daq_data.keys():
                for i, sig_burst in enumerate(daq_data[node]):
                    self.data['time'] = sig_burst.time
                    self.data[data_label].append(sig_burst.value[0])
                    if self.settings['fft']:
                        freqs, sxx = periodogram(sig_burst.value[0], fs=self.data['sampling_rate_actual'][0], window='hann',
                                                 return_onesided=True, scaling='density')
                        self.data['fft_freq'] = freqs
                        self.data['%s_fft' % data_label].append(sxx)

    def _plot(self, axes_list, data=None):
        """

        :param axes_list:
        :param data:
        :return:
        """
        if data is None:
            data = self.data

        if len(data['time']) > 0:  # Make sure some data was written before plotting
            data_plt = []
            plt_labels = []
            for data_label in self.data_labels:
                for i, segment in enumerate(self.data[data_label]):
                    data_plt.append(segment)
                    plt_labels.append('%s, chunk %i' % (data_label, i))

            plot_1d_simple_timetrace(axes_list[0], data['time'], data_plt, y_label='V', time_unit='s')
            axes_list[0].legend(plt_labels)

            if self.settings['fft']:
                fft_plt = []
                fft_labels = []
                for data_label_fft in self.data_labels_fft:
                    for i, segment in enumerate(self.data[data_label_fft]):
                        fft_plt.append(segment)
                        fft_labels.append('%s, chunk %i' % (data_label_fft, i))

                        plot_psd(data['fft_freq'], segment, axes_list[1])
                axes_list[1].set_ylabel(r'V^2/Hz')
                axes_list[1].legend(fft_labels)

        else:
            print("Data not yet collected! No data plotted.")

