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


from b26_toolkit.plotting.plots_1d import plot_psd, update_1d_simple
from pylabcontrol.core import Script, Parameter

from b26_toolkit.instruments import ZIHF2


class ZISweeper(Script):
    """
This script performs a frequency sweep with the Zurich Instrument HF2 Series Lock-in amplifier
    """
    _DEFAULT_SETTINGS = [
        Parameter('start', 1.8e6, float, 'start value of sweep'),
        Parameter('stop', 1.9e6, float, 'end value of sweep'),
        Parameter('samplecount', 101, int, 'number of data points'),
        # Parameter('gridnode', 'oscs/0/freq', ['oscs/0/freq', 'oscs/1/freq'], 'output channel =not 100% sure, double check='),
        # Parameter('output', 0, [0, 1],'output channel =not 100% sure, double check='),
        Parameter('xmapping', 'linear', ['linear', 'logarithmic'], 'mapping 0 = linear, 1 = logarithmic'),
        Parameter('ymapping', 'linear', ['linear', 'logarithmic'], 'display of y-axis'),
        Parameter('bandwidthcontrol', 'automatic', ['automatic'], '2 = automatic bandwidth control'),
        Parameter('scan', 'sequential', ['sequential', 'binary', 'bidirecctional'], 'scan direction 0 = sequential, 1 = binary (non-sequential, each point once), 2 = bidirecctional (forward then reverse)'),
        Parameter('loopcount', 1, int, 'number of times it sweeps'),
        Parameter('averaging/sample', 1, int, 'number of samples to average over')
    ]

    _INSTRUMENTS = {'zihf2' : ZIHF2}

    _SCRIPTS = {}

    def __init__(self, instruments, name = None, settings = None, log_function = None, timeout = 1000000000, data_path = None):

        self._recording = False
        self._timeout = timeout

        Script.__init__(self, name, settings, instruments, log_function= log_function, data_path = data_path)

        self.sweeper = self.instruments['zihf2']['instance'].daq.sweep(self._timeout)
        self.sweeper.set('sweep/device', self.instruments['zihf2']['instance'].device)

        self.data = deque()

        # todo: clean this up! and plot data in gui!
        self._sweep_values =  list({'frequency' : [], 'x' : [], 'y' : [], 'phase': [], 'r':[]}.keys())


    def settings_to_commands(self, settings):
        '''
        converts dictionary to list of  setting, which can then be passed to the zi controler
        :param dictionary = dictionary that contains the commands
        :return: commands = list of commands, which can then be passed to the zi controler
        '''
        # create list that is passed to the ZI controler

        commands = []
        for key, val in settings.items():
            if isinstance(val, dict) and 'value' in val:
                commands.append(['sweep/%s' % (key), val['value']])
            elif key in ('start', 'stop', 'samplecount', 'gridnode', 'loopcount', 'averaging/sample'):
                commands.append(['sweep/%s' % (key), val])
            elif key in ('xmapping'):
                if val.lower() == 'linear':
                    val = 0
                elif val.lower() == 'logarithmic':
                    val = 1
                else:
                    raise ValueError
                commands.append(['sweep/%s' % (key), val])
            elif key in ('bandwidthcontrol'):
                if val.lower() == 'automatic':
                    val = 2
                else:
                    raise ValueError
                commands.append(['sweep/%s' % (key), val])
            elif key in ('scan'):
                if val.lower() == 'sequential':
                    val = 0
                elif val.lower() == 'binary':
                    val = 1
                elif val.lower() == 'bidirecctional':
                    val = 2
                else:
                    raise ValueError
                commands.append(['sweep/%s' % (key), val])

        # settings of instrument
        out_channel = self.instruments['zihf2']['settings']['sigouts']['channel']
        commands.append(['sweep/gridnode', 'oscs/%s/freq' % (out_channel)])

        return commands


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        self.instruments['zihf2']['instance'].update(self.instruments['zihf2']['settings'])

        self.data.clear() # clear data queue
        commands = self.settings_to_commands(self.settings)
        self.sweeper.set(commands)

        path = '/%s/demods/%d/sample' % (self.instruments['zihf2']['instance'].device, self.instruments['zihf2']['instance'].settings['demods']['channel'])
        self.sweeper.subscribe(path)
        self.sweeper.execute()

        N_loops = self.settings['loopcount']
        last_progress = 0
        loopcount = 0
        while not self.sweeper.finished():
            time.sleep(1)

            new_progress = self.sweeper.progress()
            # poor mans way of keeping track of the repetitions, there should be a command to get this directly from the ZI API
            if new_progress < last_progress:
                loopcount +=1
            last_progress = new_progress

            self.progress = float(100.*(self.sweeper.progress()+loopcount) / N_loops)
            # self.sweeper.finished
            data = self.sweeper.read(True)# True: flattened dictionary

            #  ensures that first point has completed before attempting to read data
            if (path not in data) or not (data[path][0]):
                continue

            data = data[path][0][0] # the data is nested, we remove the outer brackets with [0][0]
            # now we only want a subset of the data porvided by ZI
            data = {k : data[k] for k in self._sweep_values}

            start = time.time()
            self.data.append(data)

            if (time.time() - start) > self._timeout:
                # If for some reason the sweep is blocking, force the end of the
                # measurement
                print("\nSweep still not finished, forcing finish...")
                self.sweeper.finish()
                self._recording = False

            print(("Individual sweep %.2f%% complete. \n" % (self.progress)))

            self.updateProgress.emit(int(self.progress))

        if self.sweeper.finished():
            self._recording = False



    def _plot(self, axes_list, data = None, trace_only = False):
        """
        plots the zi instrument frequency sweep

        Args:
            axes_list: list of axes to write plots to (uses first)
            data (optional): dataset to plot (dictionary that contains keys r, frequency, phase), if not provided use self.data
        """

        if data is None:
            data = self.data

        if isinstance(data, deque):
            data = data[-1]

        r = data['r']
        freq = data['frequency']
        freq = freq[np.isfinite(r)]
        phase = data['phase']
        phase = phase[np.isfinite(r)]
        r = r[np.isfinite(r)]

        x_scaling = self.settings['xmapping'][0:3]
        y_scaling = self.settings['ymapping'][0:3]

        # plot amplitude
        axes = axes_list[0]
        axes.hold(False)
        plot_psd(freq, r, axes, x_scaling = x_scaling, y_scaling = y_scaling)

        # axes.set_xlim([min(freq), max(freq)])
        axes.set_ylim([min(r), max(r)])
        axes.set_ylabel('amplitude (Vrms)')

        # plot phase
        if not trace_only:
            axes = axes_list[1]
            axes.hold(False)
            plot_psd(freq, phase, axes, x_scaling=x_scaling, y_scaling='lin')
            # axes.set_xlim([min(freq), max(freq)])
            axes.set_ylim([min(phase), max(phase)])
            axes.set_ylabel('phase (rad)')


if __name__ == '__main__':
    script, failed, instruments = Script.load_and_append(script_dict={'ZISweeper': 'ZISweeper'})


