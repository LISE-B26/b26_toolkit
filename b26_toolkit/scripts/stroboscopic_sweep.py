# """
#     This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
#     Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell
#
#     b26_toolkit is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     b26_toolkit is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with b26_toolkit.  If not, see <http://www.gnu.org/licenses/>.
# """

import numpy as np
import time, random
from pylabcontrol.core import Script, Parameter
from b26_toolkit.instruments import NI6259, NI9402, AFG3021C

TTL = 3.3  # V


class StroboscopicSweep(Script):
    """
    DaLi will never read this because he never writes doc strings for his code
    """

    _DEFAULT_SETTINGS = [
        Parameter('sweep_avg', 1, int, 'number of sweep averages'),
        Parameter('freq_start', 1e3, float, 'start frequency of scan'),
        Parameter('freq_stop', 5e4, float, 'end frequency of scan'),
        Parameter('freq_ref', 25e3, float, 'reference frequency'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan'),
        Parameter('integration_time', 0.05, float,
                  'The TOTAL integration time. The integration time per daq bins is integration_time / num_samps_per_pt'),
        Parameter('num_samps_per_pt', 100, int,
                  'Number of samples within one DAQ buffer, for each frequency. This should be large because the first measurement is thrown away (may start in the middle of a clock tick)'),
        Parameter('afg_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)'),
        Parameter('turn_off_after', False, bool, 'if true MW output is turned off after the measurement'),
        Parameter('save_full_spectrum', True, bool, 'If true save all the esr traces individually'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan'),
        Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1'], 'Daq channel used for counter'),
        Parameter('sweep_type', 'up', ['up', 'down', 'randomize'], 'type of sweep'),
        Parameter('save_timetrace', True, bool,
                  'check to save the measured fluorescence over time. This is identical to the full esr when the freq. are not randomized')
    ]

    _INSTRUMENTS = {
        'afg': AFG3021C,
        'NI6259': NI6259,  # PCI
        'NI9402': NI9402,  # cDAQ
    }

    _SCRIPTS = {}

    def __init__(self, instruments, scripts = None, name=None, settings=None, log_function=None, data_path = None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path = data_path)
        # defines which daqs contain the input and output based on user selection of daq interface
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_in = self.instruments['NI9402']['instance']

        self.afg = self.instruments['afg']['instance']

    def setup_afg(self):
        '''

        set relevant parameters on the MW generator.

        '''
        self.afg.update({'frequency': self.settings['freq_start'],
                         'amplitude': TTL,
                         'offset': TTL / 2,
                         'enable_output': False,
                         'phase': 0,
                         'function': 'Square'})

    def setup_daq(self):
        '''

        Initialize the relevant DAQ and set the sample rate.

        '''

        sample_rate = 1.0 / (self.settings['integration_time'] / self.settings['num_samps_per_pt']) # DAQ minimum buffer size is 2
        self.daq_in.settings['digital_input']['ctr0']['sample_rate'] = sample_rate

    def get_freq_array(self):
        '''

        Construct the frequency array based on the setting parameters.

        Returns:
            freq_values: array of the frequencies to be tested
            freq_range: the range of the frequency to be tested, i.e. maximum frequency - minumum frequency
        '''

        freq_start = self.settings['freq_start']
        freq_stop = self.settings['freq_stop']
        # contruct the frequency array
        if self.settings['range_type'] == 'start_stop':
            if freq_start > freq_stop:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True

            if freq_start < self.afg.FREQ_MIN or freq_stop > self.afg.FREQ_MAX: # freq range of the SRS
                self.log('start or stop frequency out of bounds')
                self._abort = True

            freq_values = np.linspace(freq_start, freq_stop, self.settings['freq_points'])

        elif self.settings['range_type'] == 'center_range':
            if freq_start < 2 * freq_stop:
                self.log('end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True

            freq_values = np.linspace(freq_start - freq_stop / 2,
                                      freq_start + freq_stop / 2, self.settings['freq_points'])

            if self.settings['freq_stop'] > 1e9:
                self.log('freq_stop (range) is quite large --- did you mean to set \'range_type\' to \'start_stop\'? ')
        else:
            self.log('unknown range parameter. Abort script')
            self._abort = True

        freq_range = np.max(freq_values) - np.min(freq_values)

        return freq_values, freq_range

    def run_sweep(self, freq_values):
        '''

        Actually runs the ESR sweep, for a single average.

        Returns:
            esr_data
            laser_data
            normalized_data

        '''

        start_run_sweep = time.time()

        num_samps = self.settings['num_samps_per_pt'] + 1 # acquire this many samples per point. The counter starts in the middle of a
                                                          # clock tick, so throw out the first sample and take num_samps_per_pt + 1 samples

        # initialize data arrays
        single_sweep_data = np.zeros(len(freq_values))

        indices = list(range(len(freq_values)))
        if self.settings['sweep_type'] == 'down':
            indices.reverse()
        elif self.settings['sweep_type'] == 'randomize':
            random.shuffle(indices)

        self.afg.update({'enable_output': True})

        t = np.linspace(0, self.settings['integration_time'], self.settings['num_samps_per_pt'])
        for freq_index in indices:

            single_freq_start_t = time.time()

            if self._abort:
                break

            # change MW frequency
            self.afg.update({'frequency': float(freq_values[freq_index])})# + self.settings['freq_ref'])})
            time.sleep(self.settings['afg_switching_time'])

            # fourierComp = np.exp(-2 * 1j * np.pi * self.settings['freq_ref'] * t) * self.measure_signal(num_samps)
            fourierComp = np.exp(-2 * 1j * np.pi * freq_values[freq_index]  * t) * self.measure_signal(num_samps)


            single_sweep_data[freq_index] = np.abs(np.sum(fourierComp)) ** 2


        # normalize  single sweep data to kcounts/sec
        single_sweep_data = single_sweep_data * (.001 / self.settings['integration_time'])**2

        end_run_sweep = time.time()


        return single_sweep_data, indices

    def measure_signal(self, num_samps):
        """

        measure the signal with the APD

        Args:
            num_samps:

        Returns:

        """

        start_meas = time.time()

        # setup the tasks
        ctrtask = self.daq_in.setup_counter(self.settings['counter_channel'], num_samps)

        t2 = time.time()-start_meas

        self.daq_in.run(ctrtask)  # the counter clock turns on and starts the AI task


        time.sleep(self.settings['integration_time'])

        t4 = time.time()-start_meas

        # read the data
        raw_data, _ = self.daq_in.read_counter(ctrtask)

        signal = np.diff(raw_data)  # take the total counts, neglecting the first element

        self.daq_in.stop(ctrtask)  # stop the clock task last

        t6 = time.time() - start_meas

        return signal

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        self.lines = []

        start_time = time.time()

        # setup the daq
        self.setup_daq()

        # setup the microwave generator
        self.setup_afg()

        # intialize some of the fields in self.data
        self.data = {'frequency': [], 'data': [], 'avrg_counts': []}

        # get the frequencices of the sweep
        freq_values, freq_range = self.get_freq_array()
        self.data.update({'frequency': freq_values})
        # initialize data arrays
        sweep_data = np.zeros((self.settings['sweep_avg'], len(freq_values))) # for the raw esr data

        # to get the timetrace we keep track of the indecies
        if self.settings['save_timetrace']:
            index_data = np.zeros((self.settings['sweep_avg'], len(freq_values))).astype(int)  # for the raw esr data


        # run sweeps
        for scan_num in range(self.settings['sweep_avg']):

            self.log('starting average ' + str(scan_num) + ', time elapsed: ' + str(np.round(time.time()-start_time, 1)) + 's')

            # get the data for a single sweep. These are raw data.
            single_sweep_data, indices = self.run_sweep(freq_values)

            if self._abort:
                break

            # save the single sweep data
            sweep_data[scan_num] = single_sweep_data

            if self.settings['save_timetrace']:
                index_data[scan_num] = indices

            # current non-normalized averaged data to plot and fit to if laser power tracking is off
            sweep_avg = np.mean(sweep_data[:(scan_num + 1)], axis=0)

            # # fit to the data
            # fit_params = fit_esr(freq_values, esr_avg, min_counts = self.settings['fit_constants']['minimum_counts'],
            #                     contrast_factor=self.settings['fit_constants']['contrast_factor'])

            self.data.update({'data': sweep_avg})

            progress = self._calc_progress(scan_num)
            self.updateProgress.emit(progress)

        if self.settings['save_full_spectrum']:
            self.data.update({'sweep_data': sweep_data})

        if self.settings['save_timetrace']:
            # using the indeces, we retrieve the right time ordering for the esr data and save it as the timetrace data
            self.data.update({'index_data':index_data})
            # self.data.update({'time_trace':[esr[idx] for esr, idx in zip(esr_data, index_data)]})


        if self.settings['turn_off_after']:
            self.afg.update({'enable_output': False})

    def _calc_progress(self, scan_num):
        #COMMENT_ME

        progress = float(scan_num - 1) / self.settings['sweep_avg'] * 100.
        self.progress = progress
        return int(progress)

    def _plot(self, axes_list, data=None):
        """
        plotting function for esr
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
            data: data (dictionary that contains keys frequency, data and fit_params) if not provided use self.data
        Returns:

        """

        if data is None:
            data = self.data

        axis = axes_list[0]
        plot_marker_data = 'b'
        linestyle = '-'
        marker = '.'

        axis.clear()  # ER 20181012 - matplotlib axes.hold() removed in update to 3.0.0
        if np.shape(data['frequency']) == np.shape(data['data']):  # ER 20190129
            axis.semilogy(data['frequency'], data['data'], plot_marker_data, linestyle=linestyle, marker=marker)
        # axes.hold(True) #ER 20181012

        title = 'Strobe'

        axis.set_title(title)
        axis.set_xlabel('Frequency (Hz)')
        axis.set_ylabel('(Kcounts/s)^2')

    def get_axes_layout(self, figure_list):
        """
        returns the axes objects the script needs to plot its data
        the default creates a single axes object on each figure
        This can/should be overwritten in a child script if more axes objects are needed
        Args:
            figure_list: a list of figure objects
        Returns:
            axes_list: a list of axes objects

        """
        new_figure_list = [figure_list[1]]
        return super().get_axes_layout(new_figure_list)