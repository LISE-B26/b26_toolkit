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
from copy import deepcopy

import numpy as np

from b26_toolkit.plotting.plots_1d import plot_psd
from pylabcontrol.core import Script, Parameter
from b26_toolkit.scripts import ZISweeper


class ZISweeperHighResolution(Script):
    """
This script takes a high resolution frequency sweep with the Zurich Instrument HF2 Lock-in amplifier.
First it acquires a sweep over a larger frequecy range. Then it finds the maximum signal and performs a second sweep with high resolution around that maximum.
    """
    _DEFAULT_SETTINGS = [
        Parameter('high_res_df', 1000, float, 'frequency step of high res. scan (Hz)'),
        Parameter('high_res_N', 21, int, 'number of data points of high res. scan'),
    ]

    _INSTRUMENTS = {}

    _SCRIPTS = {'zi sweep' : ZISweeper}

    def __init__(self, scripts, name = None, settings = None, log_function = None, timeout = 1000000000, data_path = None):
        self._recording = False
        self._timeout = timeout

        Script.__init__(self, name, settings, scripts = scripts, log_function= log_function, data_path = data_path)

        self.data = deque()

        self._sweep_values =  list({'frequency' : [], 'x' : [], 'y' : [], 'phase': [], 'r':[]}.keys())


    # def _receive_signal(self, progess_sub_script):
        # # calculate progress of this script based on progress in subscript
        #
        # if self.current_subscript == 'quick scan':
        #     progress = int(self.weights['quick scan'] * progess_sub_script)
        # elif self.current_subscript == 'high res scan':
        #     progress = int(self.weights['quick scan']*100 + self.weights['high res scan'] * progess_sub_script)
        # else:
        #     progress = None
        # # if calculated progress is 100 force it to 99, because we still have to save before script is finished
        # if progress>= 100:
        #     progress = 99
        #
        # if progress is not None:
        #     self.updateProgress.emit(progress)
        #
        # if not progess_sub_script.is_running:
        #     self.current_subscript = None

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """



        def calculate_weights():
            """
            calculate a weight inversely proportional to the expected to duration of the two steps in the
            script

            Returns: weights as a dictionary for the two steps

            """
            weights = {}


            # estimate run time of step 1 (fast sweep)
            f_range = sweeper_script.settings['stop'] - sweeper_script.settings['start']
            N_samples = sweeper_script.settings['samplecount']
            df = f_range / N_samples

            t = N_samples / df

            weights['quick scan'] = t

            # estimate run time of step 2 (high res sweep)
            df = self.settings['high_res_df']
            N_samples = self.settings['high_res_N']

            t = N_samples / df

            weights['high res scan'] = t


            total_time = sum([v for k, v in weights.items()])

            weights = {k: v/total_time for k, v in weights.items()}

            print(('weights',weights))

            return weights

        # initializes data
        self.data = {'low_res_r':None, 'low_res_freq': None, 'high_res_r':None, 'high_res_freq':None}

        sweeper_script = self.scripts['zi sweep']
        #save initial settings, so that we can reset at the end of the script
        initial_settings = deepcopy(sweeper_script.settings)
        self.weights = calculate_weights()

        # run low resolution scan
        print('run low resolution scan')
        sweeper_script.run()
        # get data from sweeper script
        self.data['low_res_r'] = deepcopy(sweeper_script.data[-1]['r'])
        self.data['low_res_freq'] = deepcopy(sweeper_script.data[-1]['frequency'])

        # find max
        fo = self.data['low_res_freq'][np.argmax(self.data['low_res_r'])]

        # calc new range
        # make sure that we convert back to native python types (numpy file types don't pass the Parameter validation)
        df = self.settings['high_res_df']
        N = int(self.settings['high_res_N'])
        f_start, f_end = float(fo - N / 2 * df), float(fo + N / 2 * df)
        print(('f_start, f_end', f_start, f_end))


        self.log('found peak at {:1.2e}'.format(fo))

        # update sweeper
        sweeper_script.update({
            'start' : f_start,
            'stop' : f_end,
            'samplecount' : N
        })

        # run high resolution scan
        print('run high resolution scan')
        sweeper_script.run()
        # get data from sweeper script
        self.data['high_res_r'] = deepcopy(sweeper_script.data[-1]['r'])
        self.data['high_res_freq'] = deepcopy(sweeper_script.data[-1]['frequency'])
        # self.data = sweeper_script.data[-1]

        # set the sweeper script back to initial settings
        sweeper_script.update(initial_settings)


    def _plot(self, axes_list, data = None):
        """
        plots the zi instrument frequency sweep

        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional): dataset to plot (dictionary that contains keys r, frequency), if not provided use self.data
        """

        if self.scripts['zi sweep'].is_running:
            print('sweeper is running')

            data = self.scripts['zi sweep'].data[-1]


            if self.data['high_res_r'] is None:
                # this means that we run the low res scan

                low_res_r = data['r']
                low_res_freq = data['frequency']
                low_res_freq = low_res_freq[np.isfinite(low_res_r)]
                low_res_r = low_res_r[np.isfinite(low_res_r)]

                high_res_freq = None
                high_res_r = None
            else:
                # this means that we run the high res scan
                high_res_r = data['r']
                high_res_freq = data['frequency']
                high_res_freq = high_res_freq[np.isfinite(high_res_r)]
                high_res_r = high_res_r[np.isfinite(high_res_r)]

                low_res_freq = None
                low_res_r = None
        else:
            if data is None:
                data = self.data
            high_res_freq = data['high_res_freq']
            high_res_r = data['high_res_r']
            low_res_freq = data['low_res_freq']
            low_res_r = data['low_res_r']

        if low_res_r is not None:
            axes = axes_list[1]
            plot_psd(low_res_freq, low_res_r, axes, y_scaling = 'lin', x_scaling = 'lin')
            axes.set_title('low resolution scan')


        if high_res_r is not None:
            axes = axes_list[0]
            plot_psd(high_res_freq, high_res_r, axes, y_scaling = 'lin', x_scaling = 'lin')
            axes.set_title('high resolution scan')
