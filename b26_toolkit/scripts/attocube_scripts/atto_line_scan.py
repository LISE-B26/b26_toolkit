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
import time
from b26_toolkit.scripts.pulse_sequences.param_sweep.pulsed_esr import PulsedEsrFast, PulsedEsrUpperLower
from pylabcontrol.core import Parameter, Script
from b26_toolkit.scripts import FindNv, EsrSimple
from b26_toolkit.scripts.find_nv import FindNvStrobe
from b26_toolkit.scripts.attocube_scripts.atto_offset_z import AttoOffsetZ
from b26_toolkit.data_analysis.nv_optical_response import B_field_from_esr
from b26_toolkit.plotting.plots_1d import plot_1d_simple_freq, update_1d_simple


class AttoLineScanGeneric(Script):

    """
    Base script for moving a scanner (e.g. Attocube) in a grid and taking a measurement at each point.
    """

    _DEFAULT_SETTINGS = [
        Parameter('voltage_start', 50, float, 'start voltage (V) of scan (excuse the misleading variable name)'),
        Parameter('voltage_stop', 20, float, 'end voltage (V) of scan (excuse the misleading variable name)'),
        Parameter('range_type', 'center_range', ['start_stop', 'center_range'],
                  'start_stop: voltage range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('voltage_points', 15, int, 'number of voltages in scan'),
        Parameter('settle_time', .5, float,
                  'time wait after changing piezo voltages (s)')
    ]

    _SCRIPTS = {'find_nv': FindNv, 'ESR_simple': EsrSimple}
    _INSTRUMENTS = {}

    def _configure_param_array(self):
        """
        Contruct the frequency array and store it in a variable called 'self.params'. Despite the naming, it's just a list of parameters to be swept;
        it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        :return:
        """

        if self.settings['range_type'] == 'start_stop':
            if self.settings['voltage_start'] > self.settings['voltage_stop']:
                self.log('end voltage must be larger than start voltage when range_type is start_stop. Abort script')
                self._abort = True
            self.params = np.linspace(self.settings['voltage_start'], self.settings['voltage_stop'], self.settings['voltage_points'])

        elif self.settings['range_type'] == 'center_range':
            self.params = np.linspace(self.settings['voltage_start'] - self.settings['voltage_stop'] / 2,
                                      self.settings['voltage_start'] + self.settings['voltage_stop'] / 2, self.settings['voltage_points'])

    def setup_scan(self):
        pass

    def check_bounds(self):
        pass

    def _function(self):

        try:
            self.setup_scan()
            self._configure_param_array()
        except:
            return

        self.data = {'point_value': np.zeros(self.settings['voltage_points']), 'params': self.params}

        self.data['params'] = self.params

        try:
            self.check_bounds()
        except AttributeError:
            self._abort = True

        self.params = self.params.tolist()

        self.point_index = 0
        pixels_completed = 0
        for i in range(0, len(self.params)):
            if self._abort:
                break
            self.move_piezo(self.params[i])
            time.sleep(self.settings['settle_time'])
            point_value = self.read_point()
            self.data['point_value'][i] = point_value
            pixels_completed += 1
            self.progress = float(pixels_completed) / len(self.params) * 100
            self.updateProgress.emit(int(self.progress))



    def read_point(self):
        """
        Replace this with some data-taking script, e.g. read the counts (to perform a galvoscan with the attocubes
        instead) or take ESR (for a 2D scan of B field)
        Returns: measured value at given pt, with other data that won't be plotted (e.g. for ESRs, we want to plot
        only the on-axis field, and save but not plot the off-axis field)
        """
        raise NotImplementedError

    def move_piezo(self, voltage):
        """
        Action to move the scanner
        :param x: location along x that scanner should move to
        :param y: location along y that scanner should move to
        :return: None
        """
        raise NotImplementedError

    def _plot(self, axes_list, data=None):
        """
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional) dataset to plot (dictionary that contains keys counts, tau), if not provided use self.data
        """

        if data is None:
            data = self.data

        if 'point_value' in data.keys():
            plot_1d_simple_freq(axes_list[0], data['params'], [data['point_value']], x_label='Attocube voltage (V)', y_label='Field (G)')

    def _update_plot(self, axes_list):
        """
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)
        """

        data = self.data
        counts = data['point_value']
        x_data = data['params']

        # If fit is found and fit has not been plotted, plot both data and fit
        fit_in_plot = len(axes_list[0].lines) == len(np.transpose(counts)) + 1
        update_1d_simple(axes_list[0], x_data, counts, fit_in_plot=fit_in_plot)


class AttoLineScanPulsedEsr(AttoLineScanGeneric):
    """
    Base script for moving a scanner (e.g. Attocube) in a grid and taking a measurement at each point.
    """

    _DEFAULT_SETTINGS = [
        Parameter('voltage_start', 50, float, 'start voltage (V) of scan (excuse the misleading variable name)'),
        Parameter('voltage_stop', 20, float, 'end voltage (V) of scan (excuse the misleading variable name)'),
        Parameter('range_type', 'center_range', ['start_stop', 'center_range'],
                  'start_stop: voltage range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('voltage_points', 15, int, 'number of voltages in scan'),
        Parameter('settle_time', .01, float,
                  'time wait after changing piezo voltages (s)'),
        Parameter('zero_field_splitting', 2.8707e9, float, 'Dip in ESR with no external B field applied'),
        Parameter('smart_scan', False, bool, 'If enabled, automatically runs PulsedEsrUpperLower when splitting becomes larger than the single scan range')
    ]

    _SCRIPTS = {'find_nv': FindNvStrobe, 'PulsedEsrFaster': PulsedEsrFast, 'PulsedEsrUpperLower': PulsedEsrUpperLower, 'atto_offset_z': AttoOffsetZ}


    def setup_scan(self):
        if self.settings['smart_scan']:
            assert self.scripts['PulsedEsrFaster'].settings['range_type'] == 'center_range'
        self.scan_script = self.scripts['PulsedEsrFaster']
        self.scripts['PulsedEsrUpperLower'].settings['freq_stop_lower'] = self.scripts['PulsedEsrFaster'].settings['freq_stop']
        self.scripts['PulsedEsrUpperLower'].settings['freq_stop_upper'] = self.scripts['PulsedEsrFaster'].settings['freq_stop']
        self.scripts['PulsedEsrUpperLower'].settings['averaging_block_size'] = self.scripts['PulsedEsrFaster'].settings['averaging_block_size']
        self.scripts['PulsedEsrUpperLower'].settings['num_averages'] = self.scripts['PulsedEsrFaster'].settings['num_averages']

        # Copy over settings from PulsedEsrFaster to UpperLower script
        for key in self.scripts['PulsedEsrFaster'].settings.keys():
            if key in self.scripts['PulsedEsrUpperLower'].settings.keys() and key != 'tag':
                print(key)
                self.scripts['PulsedEsrUpperLower'].settings[key] = self.scripts['PulsedEsrFaster'].settings[key]

    def move_piezo(self, voltage):
        self.scripts['atto_offset_z'].settings['z_offset'] = voltage
        self.scripts['atto_offset_z'].run()

    def read_point(self):
        """
        Take ESR (for a 2D scan of B field)
        Returns: measured value at given pt, with other data that won't be plotted
        """
        self.scan_script.run()

        try:
            # self.scan_script.settings['tag'] = ""%self.scan_script.settings['tag']
            point_data = self.scan_script.data['fits']
        except:
            point_data = None

        if point_data is not None and len(point_data) == 6:
            fp = float(point_data[5])
            fn = float(point_data[4])
            if self.settings['smart_scan']:
                if self.scan_script == self.scripts['PulsedEsrFaster'] and np.abs(fp-fn) > self.scripts['PulsedEsrFaster'].settings['freq_stop']*0.6:
                    self.scan_script = self.scripts['PulsedEsrUpperLower']
                    self.log('Switching to PulsedEsrUpperlower with fn=%.4f, fp=%.4f GHz' % (fn*1e-9, fp*1e-9))
                    print('Switching to PulsedEsrUpperlower with fn=%.4f, fp=%.4f GHz' % (fn * 1e-9, fp * 1e-9))
                elif self.scan_script == self.scripts['PulsedEsrUpperLower']:
                    self.log('Adjusting PulsedEsrUpperlower with fn=%.4f, fp=%.4f GHz' % (fn * 1e-9, fp * 1e-9))

                fn_new = fn + (fn - self.scripts['PulsedEsrUpperLower'].settings['freq_start_lower']) - self.scripts['PulsedEsrUpperLower'].settings['freq_stop_lower']*0.15
                fp_new = fp + (fp - self.scripts['PulsedEsrUpperLower'].settings['freq_start_upper']) + self.scripts['PulsedEsrUpperLower'].settings['freq_stop_upper']*0.15
                self.scripts['PulsedEsrUpperLower'].settings['freq_start_lower'] = fn_new
                self.scripts['PulsedEsrUpperLower'].settings['freq_start_upper'] = fp_new

            point_value = B_field_from_esr(fp, fn, D=self.settings['zero_field_splitting'],
                                           gamma=27.969e9, angular_freq=False, verbose=False)[0] * 1e4
        else:
            point_value = 0

        return point_value


