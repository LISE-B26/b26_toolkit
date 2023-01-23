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
from copy import deepcopy, copy

from pylabcontrol.core import Script, Parameter

# import standard libraries
import numpy as np
from b26_toolkit.scripts import FindNV
from b26_toolkit.instruments import MicrowaveGenerator, NI6259, NI9263, NI9402, B26PulseBlaster, RFGenerator
from b26_toolkit.scripts.spec_analyzer_get_spectrum import SpecAnalyzerGetSpectrum
from b26_toolkit.plotting.plots_1d import plot_esr


# from b26_toolkit.plotting.plots_1d import plot_diff_freq_vs_freq
from b26_toolkit.data_processing.esr_signal_processing import fit_esr
import time
import random


class ESR_DAQ_FM(Script):
    """
    This class runs ESR on an NV center, outputing microwaves using a MicrowaveGenerator and reading in NV counts using
    a DAQ. It uses FM using AO2 on the DAQ, which is off by a few MHz but may be faster than the other ESR script.

    Does not turn the PB card microwave channel off after the experiment!!

    """

    _DEFAULT_SETTINGS = [
        Parameter('power_out', -45.0, float, 'output power (dBm)'),
        Parameter('esr_avg', 50, int, 'number of esr averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'], 'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan'),
        Parameter('integration_time', 0.01, float, 'measurement time of fluorescent counts (must be a multiple of settle time)'),
        Parameter('settle_time', .0002, float, 'time wait after changing frequencies using daq (s)'),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)'),
        Parameter('turn_off_after', False, bool, 'if true MW output is turned off after the measurement'),
        Parameter('take_ref', True, bool, 'If true normalize each frequency sweep by the average counts. This should be renamed at some point because now we dont take additional data for the reference.'),
        Parameter('save_full_esr', True, bool, 'If true save all the esr traces individually'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan'),
        Parameter('fit_constants',
                  [
                      Parameter('minimum_counts', .5, float, 'minumum counts for an ESR to not be considered noise'),
                      Parameter('contrast_factor', 1.5, float, 'minimum contrast for an ESR to not be considered noise')
                  ]),
        Parameter('track_laser_power',
                  [
                      Parameter('on/off', False, bool, 'If true, measure and normalize out laser power drifts during esr'),
                      Parameter('ai_channel', 'ai4', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4'], 'channel to use for analog input, to which the photodiode is connected')
                  ]),
    ]

    _INSTRUMENTS = {
        'microwave_generator': MicrowaveGenerator,
        'NI6259': NI6259,
        'NI9263': NI9263,
        'NI9402': NI9402
    }

    _SCRIPTS = {}

    def __init__(self, instruments, scripts = None, name=None, settings=None, log_function=None, data_path = None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path = data_path)
        # defines which daqs contain the input and output based on user selection of daq interface
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['NI6259']['instance']
            self.daq_out = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_in = self.instruments['NI9402']['instance']
            self.daq_out = self.instruments['NI9263']['instance']

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """


        def get_frequency_voltages(freq_values, sec_num, dev_width, freq_array):
            """

            Args:
                freq_values: frequency values of the whole scan, including repeats for each frequency for the number of clock ticks in integration time
                sec_num: number of frequency section
                dev_width: width of frequency section


            Returns:

            """

            # calculate the minimum ad and max frequency of current section
            sec_min = min(freq_values) +  dev_width* 2 * sec_num
            sec_max = sec_min + dev_width * 2

            # make freq. array for current section
            freq_section_array = freq_array[np.where(np.logical_and(freq_array >= sec_min,freq_array < sec_max))]

            # if section is empty skip
            if len(freq_section_array) == 0:
                center_frequency = None
                freq_voltage_array = None

            else:
                center_frequency = (sec_max + sec_min) / 2.0
                freq_voltage_array = ((freq_section_array - sec_min) / (dev_width * 2)) * 2 - 1  # normalize voltages to +-1 range

            return freq_voltage_array, center_frequency

        def read_freq_section(freq_voltage_array, center_freq, clock_adjust):
            """
            reads a frequency section from the DAQ

            Args:
                freq_voltage_array: voltages corresponding to the frequency section to be measured (see get_frequency_voltages())
                center_freq:  center frequency corresponding to the frequency section to be measured (see get_frequency_voltages())
                clock_adjust: factor that specifies how many samples+1 go into the duration of the integration time in
                    order to allow for settling time. For example, if the settle time is .0002 and the integration time
                    is .01, the clock adjust is (.01+.0002)/.01 = 51, so 50 samples fit into the originally requested
                    .01 seconds, and each .01 seconds has a 1 sample (.0002 second) rest time.

            Returns: data from daq

            """
            self.instruments['microwave_generator']['instance'].update({'frequency': float(center_freq)})

            time.sleep(self.settings['mw_generator_switching_time'])

            ctrtask = self.daq_in.setup_counter("ctr0", len(freq_voltage_array) + 1)
            aotask = self.daq_out.setup_AO(["ao2"], freq_voltage_array, ctrtask)

            if self.settings['track_laser_power']['on/off'] == True and self.settings['daq_type'] == 'PCI':
                aitask = self.daq_in.setup_AI(self.settings['track_laser_power']['ai_channel'], len(freq_voltage_array), continuous=False, clk_source=ctrtask) # for optional arguments spell out every one if there are multiple
            elif self.settings['track_laser_power']['on/off'] == True and self.settings['daq_type'] != 'PCI':
                raise(NotImplementedError('cant use laser power tracking without the PCI daq!!'))

            # start counter and scanning sequence
            if self.settings['track_laser_power']['on/off'] == True and self.settings['daq_type'] == 'PCI':
                self.daq_in.run(aitask) # AI is actually tied to the clock, when this runs the clock actually starts
            self.daq_out.run(aotask) # AO is actually tied to the clock, when this runs the clock actually starts
            self.daq_in.run(ctrtask) # the counter clock turns on and starts the AO task

            self.daq_out.waitToFinish(aotask) # should only have to wait for one task to finish since ai and ao task are locked to the same clock
            self.daq_out.stop(aotask) # don't stop the ai task yet!

            raw_data, _ = self.daq_in.read_counter(ctrtask)
            if self.settings['track_laser_power']['on/off'] == True and self.settings['daq_type'] == 'PCI':
                raw_data_laser, num_read = self.daq_in.read(aitask)

            # raw_data = sweep_mw_and_count_APD(freq_voltage_array, dt)
            # counter counts continiously so we take the difference to get the counts per time interval
            diff_data = np.diff(raw_data)
            laser_data = np.zeros(int(len(freq_voltage_array) / clock_adjust))
            summed_data = np.zeros(int(len(freq_voltage_array) / clock_adjust))
            normalized_data = np.zeros(int(len(freq_voltage_array) / clock_adjust))

            for i in range(0, int((len(freq_voltage_array) / clock_adjust))):
                summed_data[i] = np.sum(diff_data[(i * clock_adjust + 1):(i * clock_adjust + clock_adjust - 1)])
                if self.settings['track_laser_power']['on/off'] == True and self.settings['daq_type'] == 'PCI':
                    laser_data[i] = np.sum(raw_data_laser[(i * clock_adjust + 1):(i * clock_adjust + clock_adjust - 1)])
            if self.settings['track_laser_power']['on/off'] == True and self.settings['daq_type'] == 'PCI':
                normalized_data = np.divide(np.multiply(summed_data, np.mean(laser_data)), laser_data)

            # clean up APD tasks
            if self.settings['track_laser_power']['on/off'] == True and self.settings['daq_type'] == 'PCI':
                self.daq_in.stop(aitask) # only stop teh ai task when you've extracted the data you need!!

            self.daq_in.stop(ctrtask) # stop the clock task last

            return summed_data, normalized_data, laser_data

        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['NI6259']['instance']
            self.daq_out = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_in = self.instruments['NI9402']['instance']
            self.daq_out = self.instruments['NI9263']['instance']

        self.lines = []

        take_ref = self.settings['take_ref']

        # contruct the frequency array
        if self.settings['range_type'] == 'start_stop':
            if self.settings['freq_start']>self.settings['freq_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True

            if self.settings['freq_start'] < 0 or self.settings['freq_stop'] > 4.05E9:
                self.log('start or stop frequency out of bounds')
                self._abort = True

            freq_values = np.linspace(self.settings['freq_start'], self.settings['freq_stop'], self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)
        elif self.settings['range_type'] == 'center_range':
            if self.settings['freq_start'] < 2 * self.settings['freq_stop']:
                self.log('end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True
            freq_values = np.linspace(self.settings['freq_start']-self.settings['freq_stop']/2,
                                      self.settings['freq_start']+self.settings['freq_stop']/2, self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

            if self.settings['freq_stop'] > 1e9:
                self.log('freq_stop (range) is quite large --- did you mean to set \'range_type\' to \'start_stop\'? ')
        else:
            self.log('unknown range parameter. Abort script')
            self._abort = True


        num_freq_sections = int(freq_range) / int(self.instruments['microwave_generator']['instance'].settings['dev_width']*2) + 1
        clock_adjust = int((self.settings['integration_time'] + self.settings['settle_time']) / self.settings['settle_time'])
        freq_array = np.repeat(freq_values, clock_adjust)
        self.instruments['microwave_generator']['instance'].update({'amplitude': self.settings['power_out']})
        self.instruments['microwave_generator']['instance'].update({'modulation_type': 'FM'})
      #  self.instruments['microwave_generator']['instance'].update({'dev_width': 3.2E7})

        # ER 20171128
        self.instruments['microwave_generator']['instance'].update({'dev_width': self.instruments['microwave_generator']['instance'].settings['dev_width']})
        self.instruments['microwave_generator']['instance'].update({'enable_modulation': True})

        sample_rate = float(1) / self.settings['settle_time']
        self.daq_out.settings['analog_output']['ao2']['sample_rate'] = sample_rate
        self.daq_in.settings['digital_input']['ctr0']['sample_rate'] = sample_rate

        self.instruments['microwave_generator']['instance'].update({'enable_output': True})

        esr_data = np.zeros((self.settings['esr_avg'], len(freq_values)))
        full_laser_data = np.zeros((self.settings['esr_avg'], len(freq_values)))
        full_normalized_data = np.zeros((self.settings['esr_avg'], len(freq_values)))
        avrg_counts = np.zeros(self.settings['esr_avg']) # here we save the avrg of the esr scan which we will use to normalize
        norm_data = np.zeros(len(freq_values))
        self.data = {'frequency': [], 'data': [], 'fit_params': [], 'avrg_counts' : avrg_counts}

        # run sweeps
        for scan_num in range(0, self.settings['esr_avg']):
            if self._abort:
                break
            esr_data_pos = 0

            for sec_num in range(0, int(num_freq_sections)):

                freq_voltage_array, center_freq = get_frequency_voltages(freq_values,
                                                                         sec_num,
                                                                         self.instruments['microwave_generator']['instance'].settings['dev_width'],
                                                                         freq_array)
                # if section is empty skip
                if len(freq_voltage_array) is None:
                    continue

                summed_data, normalized_data, laser_data = read_freq_section(freq_voltage_array, center_freq, clock_adjust)

                # also normalizing to kcounts/sec
                esr_data[scan_num, esr_data_pos:(esr_data_pos + len(summed_data))] = summed_data * (.001 / self.settings['integration_time'])
                full_laser_data[scan_num, esr_data_pos:(esr_data_pos + len(laser_data))] = laser_data
                full_normalized_data[scan_num, esr_data_pos:(esr_data_pos + len(laser_data))] = normalized_data * (.001 / self.settings['integration_time'])
                esr_data_pos += len(summed_data)


            avrg_counts[scan_num] = np.mean(esr_data[scan_num])
            norm_data = np.mean(full_normalized_data[0:(scan_num + 1)] , axis=0)

            if take_ref is True:
                esr_data[scan_num] /=avrg_counts[scan_num]

            esr_avg = np.mean(esr_data[0:(scan_num + 1)] , axis=0)
            if not self.settings['track_laser_power']['on/off']:
                fit_params = fit_esr(freq_values, esr_avg, min_counts = self.settings['fit_constants']['minimum_counts'],
                                    contrast_factor=self.settings['fit_constants']['contrast_factor'])
            elif self.settings['track_laser_power']['on/off'] == True and self.settings['daq_type'] == 'PCI':
                fit_params = fit_esr(freq_values, norm_data, min_counts = self.settings['fit_constants']['minimum_counts'],
                                    contrast_factor=self.settings['fit_constants']['contrast_factor'])
            self.data.update({'frequency': freq_values, 'data': esr_avg, 'fit_params': fit_params})

            if self.settings['track_laser_power']['on/off'] == True and self.settings['daq_type'] == 'PCI':
                self.data.update({'full_normalized_data': full_normalized_data, 'full_laser_data': full_laser_data})
                self.data.update({'norm_data': norm_data})

            if self.settings['save_full_esr']:
                self.data.update({'esr_data':esr_data})


            progress = self._calc_progress(scan_num)
            self.updateProgress.emit(progress)

        if self.settings['turn_off_after']:
            self.instruments['microwave_generator']['instance'].update({'enable_output': False})


    def _calc_progress(self, scan_num):
        #COMMENT_ME

        progress = float(scan_num) / self.settings['esr_avg'] * 100.
        self.progress = progress
        return int(progress)

    def _plot(self, axes_list, data = None):
        """
        plotting function for esr
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
            data: data (dictionary that contains keys frequency, data and fit_params) if not provided use self.data
        Returns:

        """
        if data is None:
            data = self.data
        if not self.settings['track_laser_power']['on/off']:
            plot_esr(axes_list[0], data['frequency'], data['data'], data['fit_params'])
        elif self.settings['track_laser_power']['on/off'] and self.settings['daq_type'] == 'PCI':
            plot_esr(axes_list[0], data['frequency'], data['norm_data'], data['fit_params'])


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
        return super(ESR_DAQ_FM, self).get_axes_layout(new_figure_list)

class ESR_tracking(Script):
    """
    This class runs ESR on an NV center, outputing microwaves using a MicrowaveGenerator and reading in NV counts using
    a DAQ. Each frequency is set explicitly on the SRS, instead of using FM.

    This script includes an option to track using findnv

    MM 072619
    """

    _DEFAULT_SETTINGS = [
        Parameter('power_out', -25.0, float, 'output power (dBm)'),
        Parameter('esr_avg', 20, int, 'number of esr averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan'),
        Parameter('integration_time', 0.05, float,
                  'The TOTAL integration time. The integration time per daq bins is integration_time / num_samps_per_pt'),
        Parameter('num_samps_per_pt', 100, int,
                  'Number of samples within one DAQ buffer, for each frequency. This should be large because the first measurement is thrown away (may start in the middle of a clock tick)'),
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)'),
        Parameter('turn_off_after', False, bool, 'if true MW output is turned off after the measurement'),
        Parameter('save_full_esr', True, bool, 'If true save all the esr traces individually'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan'),
        Parameter('fit_constants',
                  [
                      Parameter('minimum_counts', .5, float, 'minumum counts for an ESR to not be considered noise'),
                      Parameter('contrast_factor', 1.5, float, 'minimum contrast for an ESR to not be considered noise')
                  ]),
        Parameter('randomize', True, bool, 'check to randomize esr frequencies'),
        Parameter('save_timetrace', True, bool,
                  'check to save the measured fluorescence over time. This is identical to the full esr when the freq. are not randomized'),
        Parameter('Tracking', [
            Parameter('on/off', True, bool, 'used to turn on tracking'),
            Parameter('n_loops', 0.85, float, 'track once per n loops')
        ])
    ]

    _INSTRUMENTS = {
        'microwave_generator': MicrowaveGenerator,
        'NI6259': NI6259,  # PCI
        'NI9402': NI9402,  # cDAQ
    }

    _SCRIPTS = {'find_nv': FindNV}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)
        # defines which daqs contain the input and output based on user selection of daq interface
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_in = self.instruments['NI9402']['instance']

    def setup_microwave_gen(self):
        '''

        set relevant parameters on the MW generator.

        '''
        self.instruments['microwave_generator']['instance'].update({'amplitude': self.settings['power_out']})
        self.instruments['microwave_generator']['instance'].update({'enable_modulation': False})
        self.instruments['microwave_generator']['instance'].update({'enable_output': True})

    def setup_daq(self):
        '''

        Initialize the relevant DAQ and set the sample rate.

        '''
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_in = self.instruments['NI9402']['instance']

        sample_rate = float(1) / (self.settings['integration_time'] / self.settings[
            'num_samps_per_pt'])  # DAQ minimum buffer size is 2
        self.daq_in.settings['digital_input']['ctr0']['sample_rate'] = sample_rate

    def get_freq_array(self):
        '''

        Construct the frequency array based on the setting parameters.

        Returns:
            freq_values: array of the frequencies to be tested
            freq_range: the range of the frequency to be tested, i.e. maximum frequency - minumum frequency
        '''

        # contruct the frequency array
        if self.settings['range_type'] == 'start_stop':
            if self.settings['freq_start'] > self.settings['freq_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True

            if self.settings['freq_start'] < 950E3 or self.settings['freq_stop'] > 4.05E9:  # freq range of the SRS
                self.log('start or stop frequency out of bounds')
                self._abort = True

            freq_values = np.linspace(self.settings['freq_start'], self.settings['freq_stop'],
                                      self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

        elif self.settings['range_type'] == 'center_range':
            if self.settings['freq_start'] < 2 * self.settings['freq_stop']:
                self.log(
                    'end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True
            freq_values = np.linspace(self.settings['freq_start'] - self.settings['freq_stop'] / 2,
                                      self.settings['freq_start'] + self.settings['freq_stop'] / 2,
                                      self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

            if self.settings['freq_stop'] > 1e9:
                self.log('freq_stop (range) is quite large --- did you mean to set \'range_type\' to \'start_stop\'? ')
        else:
            self.log('unknown range parameter. Abort script')
            self._abort = True

        return freq_values, freq_range

    def run_sweep(self, freq_values, time_tracking=False):
        '''

        Actually runs the ESR sweep, for a single average.

        Returns:
            esr_data
            laser_data
            normalized_data

        '''

        if time_tracking:
            start_run_sweep = time.time()
        num_samps = self.settings[
                        'num_samps_per_pt'] + 1  # acquire this many samples per point. The counter starts in the middle of a
        # clock tick, so throw out the first sample and take num_samps_per_pt + 1 samples

        # initialize data arrays
        single_sweep_data = np.zeros(len(freq_values))

        indeces = list(range(len(freq_values)))
        if self.settings['randomize']:
            random.shuffle(indeces)

        for freq_index in indeces:

            single_freq_start_t = time.time()

            if self._abort:
                break

            freq = freq_values[freq_index]

            # change MW frequency
            self.instruments['microwave_generator']['instance'].update({'frequency': float(freq)})
            time.sleep(self.settings['mw_generator_switching_time'])

            if time_tracking:
                print('time for a single SRS frequency change: ' + str(time.time() - single_freq_start_t) + ' s')

            single_sweep_data[freq_index] = self.measure_signal(num_samps)

            if time_tracking:
                print('time for a single frequency data point: ' + str(time.time() - single_freq_start_t) + ' s')

        # normalize  single sweep data to kcounts/sec
        single_sweep_data = single_sweep_data * (.001 / self.settings['integration_time'])

        end_run_sweep = time.time()

        if time_tracking:
            print('run sweep took: ' + str(end_run_sweep - start_run_sweep) + ' s')

        return single_sweep_data, indeces

    def measure_signal(self, num_samps, time_tracking=False):
        """

        measure the signal with the APD

        Args:
            num_samps:

        Returns:

        """
        if time_tracking:
            start_meas = time.time()

        # setup the tasks
        ctrtask = self.daq_in.setup_counter("ctr0", num_samps)

        if time_tracking:
            t2 = time.time() - start_meas
            print('time to setup the task: ' + str(t2) + ' s')

        self.daq_in.run(ctrtask)  # the counter clock turns on and starts the AI task

        time.sleep(self.settings['integration_time'])

        if time_tracking:
            t4 = time.time() - start_meas
            print('time elapsed up to sleep: ' + str(t4) + ' s')

        # read the data
        raw_data, _ = self.daq_in.read_counter(ctrtask)

        signal = np.sum(np.diff(raw_data))  # take the total counts, neglecting the first element

        self.daq_in.stop(ctrtask)  # stop the clock task last

        if time_tracking:
            t6 = time.time() - start_meas
            print('time elapsed up to stopping the task: ' + str(t6) + ' s')
            print('diff-ing all times: ')
            print('setup: ', t2)
            print('running: ', t4 - t2)
            print('stopping: ', t6 - t4)
            print('total time: ', time.time() - start_meas)

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
        self.setup_microwave_gen()

        # intialize some of the fields in self.data
        self.data = {'frequency': [], 'data': [], 'fit_params': [], 'avrg_counts': []}

        # get the frequencices of the sweep
        freq_values, freq_range = self.get_freq_array()
        self.data.update({'frequency': freq_values})
        # initialize data arrays
        esr_data = np.zeros((self.settings['esr_avg'], len(freq_values)))  # for the raw esr data

        # to get the timetrace we keep track of the indecies
        if self.settings['save_timetrace']:
            index_data = np.zeros((self.settings['esr_avg'], len(freq_values))).astype(int)  # for the raw esr data

        # run sweeps
        for scan_num in range(0, self.settings['esr_avg']):

            if self._abort:
                break

            if self.settings['Tracking']['on/off'] and scan_num % self.settings['Tracking']['n_loops'] == 0:
                self.scripts['find_nv'].run()
                self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']
                self.setup_daq()


            self.log('starting average ' + str(scan_num) + ', time elapsed: ' + str(
                np.round(time.time() - start_time, 1)) + 's')

            # get the data for a single sweep. These are raw data.
            single_sweep_data, indeces = self.run_sweep(freq_values)

            # save the single sweep data
            esr_data[scan_num] = single_sweep_data

            if self.settings['save_timetrace']:
                index_data[scan_num] = indeces

            # current non-normalized averaged data to plot and fit to if laser power tracking is off
            esr_avg = np.mean(esr_data[0:(scan_num + 1)], axis=0)

            # fit to the data
            fit_params = fit_esr(freq_values, esr_avg, min_counts=self.settings['fit_constants']['minimum_counts'],
                                 contrast_factor=self.settings['fit_constants']['contrast_factor'])

            self.data.update({'data': esr_avg, 'fit_params': fit_params})

            progress = self._calc_progress(scan_num)
            self.updateProgress.emit(progress)

        if self.settings['save_full_esr']:
            self.data.update({'esr_data': esr_data})

        if self.settings['save_timetrace']:
            # using the indeces, we retrieve the right time ordering for the esr data and save it as the timetrace data
            self.data.update({'index_data': index_data})
            # self.data.update({'time_trace':[esr[idx] for esr, idx in zip(esr_data, index_data)]})


        if self.settings['turn_off_after']:
            self.instruments['microwave_generator']['instance'].update({'enable_output': False})

    def _calc_progress(self, scan_num):
        # COMMENT_ME

        progress = float(scan_num) / self.settings['esr_avg'] * 100.
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

        plot_esr(axes_list[0], data['frequency'], data['data'], data['fit_params'])

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
        return super(ESR_tracking, self).get_axes_layout(new_figure_list)

class ESR_simple(Script):

    """
    This class runs ESR on an NV center, outputing microwaves using a MicrowaveGenerator and reading in NV counts using
    a DAQ. Each frequency is set explicitly on the SRS, instead of using FM.

    This is the simple, i.e. minimal version of it
    """

    _DEFAULT_SETTINGS = [
        Parameter('power_out', -45.0, float, 'output power (dBm)'),
        Parameter('esr_avg', 50, int, 'number of esr averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'], 'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan'),
        Parameter('integration_time', 0.05, float, 'The TOTAL integration time. The integration time per daq bins is integration_time / num_samps_per_pt'),
        Parameter('num_samps_per_pt', 100, int, 'Number of samples within one DAQ buffer, for each frequency. This should be large because the first measurement is thrown away (may start in the middle of a clock tick)'),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)'),
        Parameter('turn_off_after', False, bool, 'if true MW output is turned off after the measurement'),
        Parameter('save_full_esr', True, bool, 'If true save all the esr traces individually'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan'),
        Parameter('fit_constants',
                  [
                      Parameter('minimum_counts', .5, float, 'minumum counts for an ESR to not be considered noise'),
                      Parameter('contrast_factor', 1.5, float, 'minimum contrast for an ESR to not be considered noise')
                  ]),
        Parameter('randomize', True, bool, 'check to randomize esr frequencies'),
        Parameter('save_timetrace', True, bool, 'check to save the measured fluorescence over time. This is identical to the full esr when the freq. are not randomized')
    ]

    _INSTRUMENTS = {
        'microwave_generator': MicrowaveGenerator,
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

    def setup_microwave_gen(self):
        '''

        set relevant parameters on the MW generator.

        '''
        self.instruments['microwave_generator']['instance'].update({'amplitude': self.settings['power_out']})
        self.instruments['microwave_generator']['instance'].update({'enable_modulation': False})
        self.instruments['microwave_generator']['instance'].update({'enable_output': True})

    def setup_daq(self):
        '''

        Initialize the relevant DAQ and set the sample rate.

        '''
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_in = self.instruments['NI9402']['instance']

        sample_rate = float(1) / (self.settings['integration_time']/self.settings['num_samps_per_pt']) # DAQ minimum buffer size is 2
        self.daq_in.settings['digital_input']['ctr0']['sample_rate'] = sample_rate

    def get_freq_array(self):
        '''

        Construct the frequency array based on the setting parameters.

        Returns:
            freq_values: array of the frequencies to be tested
            freq_range: the range of the frequency to be tested, i.e. maximum frequency - minumum frequency
        '''

        # contruct the frequency array
        if self.settings['range_type'] == 'start_stop':
            if self.settings['freq_start']>self.settings['freq_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True

            if self.settings['freq_start'] < 950E3 or self.settings['freq_stop'] > 4.05E9: # freq range of the SRS
                self.log('start or stop frequency out of bounds')
                self._abort = True

            freq_values = np.linspace(self.settings['freq_start'], self.settings['freq_stop'], self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

        elif self.settings['range_type'] == 'center_range':
            if self.settings['freq_start'] < 2 * self.settings['freq_stop']:
                self.log('end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True
            freq_values = np.linspace(self.settings['freq_start']-self.settings['freq_stop']/2,
                                      self.settings['freq_start']+self.settings['freq_stop']/2, self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

            if self.settings['freq_stop'] > 1e9:
                self.log('freq_stop (range) is quite large --- did you mean to set \'range_type\' to \'start_stop\'? ')
        else:
            self.log('unknown range parameter. Abort script')
            self._abort = True

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

        indeces = list(range(len(freq_values)))
        if self.settings['randomize']:
            random.shuffle(indeces)

        for freq_index in indeces:

            single_freq_start_t = time.time()

            if self._abort:
                break

            freq = freq_values[freq_index]

            # change MW frequency
            self.instruments['microwave_generator']['instance'].update({'frequency': float(freq)})
            time.sleep(self.settings['mw_generator_switching_time'])


            single_sweep_data[freq_index] = self.measure_signal(num_samps)


        # normalize  single sweep data to kcounts/sec
        single_sweep_data = single_sweep_data * (.001 / self.settings['integration_time'])

        end_run_sweep = time.time()


        return single_sweep_data, indeces

    def measure_signal(self, num_samps):
        """

        measure the signal with the APD

        Args:
            num_samps:

        Returns:

        """

        start_meas = time.time()

        # setup the tasks
        ctrtask = self.daq_in.setup_counter("ctr0", num_samps)

        t2 = time.time()-start_meas

        self.daq_in.run(ctrtask)  # the counter clock turns on and starts the AI task


        time.sleep(self.settings['integration_time'])

        t4 = time.time()-start_meas

        # read the data
        raw_data, _ = self.daq_in.read_counter(ctrtask)

        signal = np.sum(np.diff(raw_data))  # take the total counts, neglecting the first element

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
        self.setup_microwave_gen()

        # intialize some of the fields in self.data
        self.data = {'frequency': [], 'data': [], 'fit_params': [], 'avrg_counts': []}

        # get the frequencices of the sweep
        freq_values, freq_range = self.get_freq_array()
        self.data.update({'frequency': freq_values})
        # initialize data arrays
        esr_data = np.zeros((self.settings['esr_avg'], len(freq_values))) # for the raw esr data

        # to get the timetrace we keep track of the indecies
        if self.settings['save_timetrace']:
            index_data = np.zeros((self.settings['esr_avg'], len(freq_values))).astype(int)  # for the raw esr data


        # run sweeps
        for scan_num in range(0, self.settings['esr_avg']):

            self.log('starting average ' + str(scan_num) + ', time elapsed: ' + str(np.round(time.time()-start_time, 1)) + 's')

            # get the data for a single sweep. These are raw data.
            single_sweep_data, indeces = self.run_sweep(freq_values)

            if self._abort:
                break

            # save the single sweep data
            esr_data[scan_num] = single_sweep_data

            if self.settings['save_timetrace']:
                index_data[scan_num] = indeces

            # current non-normalized averaged data to plot and fit to if laser power tracking is off
            esr_avg = np.mean(esr_data[0:(scan_num + 1)], axis=0)

            # fit to the data
            fit_params = fit_esr(freq_values, esr_avg, min_counts = self.settings['fit_constants']['minimum_counts'],
                                contrast_factor=self.settings['fit_constants']['contrast_factor'])

            self.data.update({'data': esr_avg, 'fit_params': fit_params})

            progress = self._calc_progress(scan_num)
            self.updateProgress.emit(progress)

        if self.settings['save_full_esr']:
            self.data.update({'esr_data': esr_data})

        if self.settings['save_timetrace']:
            # using the indeces, we retrieve the right time ordering for the esr data and save it as the timetrace data
            self.data.update({'index_data':index_data})
            # self.data.update({'time_trace':[esr[idx] for esr, idx in zip(esr_data, index_data)]})


        if self.settings['turn_off_after']:
            self.instruments['microwave_generator']['instance'].update({'enable_output': False})

    def _calc_progress(self, scan_num):
        #COMMENT_ME

        progress = float(scan_num) / self.settings['esr_avg'] * 100.
        self.progress = progress
        return int(progress)

    def _plot(self, axes_list, data = None):
        """
        plotting function for esr
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
            data: data (dictionary that contains keys frequency, data and fit_params) if not provided use self.data
        Returns:

        """

        if data is None:
            data = self.data

        plot_esr(axes_list[0], data['frequency'], data['data'], data['fit_params'])

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
        return super(ESR_simple, self).get_axes_layout(new_figure_list)

class ESR_simple_lowerupper(Script):

    _DEFAULT_SETTINGS = [
        Parameter('frequency_1', 2.8e9, float, 'first center frequency (Hz)'),
        Parameter('frequency_2', 3.1e9, float, 'second center frequency (Hz)'),
        Parameter('frequency_width', 100.0e6, float, 'frequency width of ESR scan')
    ]

    _INSTRUMENTS = {}

    _SCRIPTS = {'esr_simple': ESR_simple}


    def _function(self):

        freqs = [self.settings['frequency_1'], self.settings['frequency_2']]
        width = self.settings['frequency_width']

        for freq in freqs:
            self.data = copy(self.scripts['esr_simple'].data)

            self.scripts['esr_simple'].settings['freq_start'] = freq
            self.scripts['esr_simple'].settings['freq_stop'] = width
            self.scripts['esr_simple'].settings['range_type'] = 'center_range'
            self.scripts['esr_simple'].run()
            self.data = copy(self.scripts['esr_simple'].data)

    def _plot(self, axes_list, data = None):
        """
        plotting function for find_nv
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
            data: data (dictionary that contains keys image_data, extent, initial_point, maximum_point) if not provided use self.data
        """

        print('inside _plot!')
        if data is None:
            data = self.data
      #  print(self._current_subscript_stage['current_subscript'])

        if self._current_subscript_stage['current_subscript'] == self.scripts['esr_simple']:
            print('got to here!! ')
            self.scripts['esr_simple']._plot(axes_list)


    def _update_plot(self, axes_list):
        """
        update plotting function for find_nv
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
        """
        print('inside _update_plot!')
    #    print(self._current_subscript_stage['current_subscript'])
        if self._current_subscript_stage['current_subscript'] == self.scripts['esr_simple']:
            print('got to update plot!! ')
            self.scripts['esr_simple']._update_plot(axes_list)

# re-written by ER 20180831 to enable esr without frequency modulation
class ESR(Script):
    """
    This class runs ESR on an NV center, outputing microwaves using a MicrowaveGenerator and reading in NV counts using
    a DAQ. Each frequency is set explicitly on the SRS, instead of using FM.
    """

    _DEFAULT_SETTINGS = [
        Parameter('power_out', -45.0, float, 'output power (dBm)'),
        Parameter('esr_avg', 1, int, 'number of esr averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'], 'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan'),
        Parameter('integration_time', 0.05, float, 'The TOTAL integration time. The integration time per daq bins is integration_time / num_samps_per_pt'),
        Parameter('num_samps_per_pt', 100, int, 'Number of samples within one DAQ buffer, for each frequency. This should be large because the first measurement is thrown away (may start in the middle of a clock tick)'),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)'),
        Parameter('turn_off_after', False, bool, 'if true MW output is turned off after the measurement'),
        Parameter('norm_to_ref', True, bool, 'If true normalize each frequency sweep by the average counts.'),
        Parameter('save_full_esr', True, bool, 'If true save all the esr traces individually'),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan'),
        Parameter('fit_constants',
                  [
                      Parameter('minimum_counts', .5, float, 'minumum counts for an ESR to not be considered noise'),
                      Parameter('contrast_factor', 1.5, float, 'minimum contrast for an ESR to not be considered noise')
                  ]),
        Parameter('track_laser_power',
                  [
                      Parameter('on/off', False, bool, 'If true, measure and normalize out laser power drifts during esr'),
                      Parameter('ai_channel', 'ai4', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4'], 'channel to use for analog input, to which the photodiode is connected')
                  ]),
        Parameter('randomize', True, bool, 'check to randomize esr frequencies'),
    ]

    _INSTRUMENTS = {
        'microwave_generator': MicrowaveGenerator,
        'NI6259': NI6259,
        'NI9263': NI9263,
        'NI9402': NI9402,
        'PB': B26PulseBlaster
    }

    _SCRIPTS = {}

    def __init__(self, instruments, scripts = None, name=None, settings=None, log_function=None, data_path = None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path = data_path)
        # defines which daqs contain the input and output based on user selection of daq interface
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_in = self.instruments['NI9402']['instance']

    def setup_microwave_gen(self):
        '''

        set relevant parameters on the MW generator.

        '''
        self.instruments['microwave_generator']['instance'].update({'amplitude': self.settings['power_out']})
        self.instruments['microwave_generator']['instance'].update({'enable_modulation': False})
        self.instruments['microwave_generator']['instance'].update({'enable_output': True})

    def setup_daq(self):
        '''

        Initialize the relevant DAQ and set the sample rate.

        '''
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_in = self.instruments['NI9402']['instance']

        sample_rate = float(1) / (self.settings['integration_time']/self.settings['num_samps_per_pt']) # DAQ minimum buffer size is 2, so we break the integration time in half
        self.daq_in.settings['digital_input']['ctr0']['sample_rate'] = sample_rate

    def setup_pb(self): # ER 20181017
        '''

        Setup the channels on the PB card.

        '''
        if self.instruments['microwave_generator']['instance'].amplitude < -10.0:
            self.instruments['PB']['instance'].update({'microwave_switch': {'status': True}})
        else:
            self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

    def get_freq_array(self):
        '''

        Construct the frequency array based on the setting parameters.

        Returns:
            freq_values: array of the frequencies to be tested
            freq_range: the range of the frequency to be tested, i.e. maximum frequency - minumum frequency
        '''

        # contruct the frequency array
        if self.settings['range_type'] == 'start_stop':
            if self.settings['freq_start']>self.settings['freq_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True

            if self.settings['freq_start'] < 950E3 or self.settings['freq_stop'] > 4.05E9: # freq range of the SRS
                self.log('start or stop frequency out of bounds')
                self._abort = True

            freq_values = np.linspace(self.settings['freq_start'], self.settings['freq_stop'], self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

        elif self.settings['range_type'] == 'center_range':
            if self.settings['freq_start'] < 2 * self.settings['freq_stop']:
                self.log('end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True
            freq_values = np.linspace(self.settings['freq_start']-self.settings['freq_stop']/2,
                                      self.settings['freq_start']+self.settings['freq_stop']/2, self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

            if self.settings['freq_stop'] > 1e9:
                self.log('freq_stop (range) is quite large --- did you mean to set \'range_type\' to \'start_stop\'? ')
        else:
            self.log('unknown range parameter. Abort script')
            self._abort = True

        return freq_values, freq_range

    def run_sweep(self, freq_values):
        '''

        Actually runs the ESR sweep, for a single average.

        Returns:
            esr_data
            laser_data
            normalized_data

        '''

        num_samps = self.settings['num_samps_per_pt'] + 1 # acquire this many samples per point. The counter starts in the middle of a
                                                          # clock tick, so throw out the first sample and take num_samps_per_pt + 1 samples

        # initialize data arrays
        single_sweep_data = np.zeros(len(freq_values))
        single_sweep_laser_data = np.zeros(len(freq_values))

        indeces = list(range(len(freq_values)))
        if self.settings['randomize']:
            random.shuffle(indeces)

        print("executing run_sweep!")
        for freq_index in indeces:

            if self._abort:
                break

            freq = freq_values[freq_index]

            # change MW frequency
            self.instruments['microwave_generator']['instance'].update({'frequency': float(freq)})
            time.sleep(self.settings['mw_generator_switching_time'])

            # setup the tasks
            ctrtask = self.daq_in.setup_counter("ctr0", num_samps)
            if self.settings['track_laser_power']['on/off']:
                aitask = self.daq_in.setup_AI(self.settings['track_laser_power']['ai_channel'], num_samps,
                                              continuous=False, clk_source=ctrtask)

            if self.settings['track_laser_power']['on/off']:
                self.daq_in.run(aitask) # AI is actually tied to the clock, when this runs the clock actually starts

            self.daq_in.run(ctrtask) # the counter clock turns on and starts the AI task
            time.sleep(self.settings['integration_time'])
           # if self.settings['track_laser_power']['on/off']:
           #     self.daq_out.waitToFinish(aitask) # wait for the tasks to be done

            # read the data
            raw_data, _ = self.daq_in.read_counter(ctrtask)
            single_sweep_data[freq_index] = np.sum(np.diff(raw_data)) # take the total counts, neglecting the first element

            if self.settings['track_laser_power']['on/off']:
                raw_data_laser, _ = self.daq_in.read(aitask)
                single_sweep_laser_data[freq_index] = np.mean(raw_data_laser)

            # clean up APD tasks
            if self.settings['track_laser_power']['on/off']:
                self.daq_in.stop(aitask)  # only stop teh ai task when you've extracted the data you need!!
            self.daq_in.stop(ctrtask)  # stop the clock task last

        return single_sweep_data, single_sweep_laser_data

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        start_time = time.time()

        # if tracking laser power drifts, check for PCI daq
        if self.settings['track_laser_power']['on/off'] and not self.settings['daq_type'] == 'PCI':
            print("tracking laser power drifts only enabled for PCI daq")
            self._abort = True

        self.lines = []
        take_ref = self.settings['norm_to_ref']

        # setup the daq
        self.setup_daq()

        # setup the microwave generator
        self.setup_microwave_gen()

        # setup the pulseblaster card (i.e., turn mw_switch on)
        self.setup_pb()

        # intialize some of the fields in self.data
        self.data = {'frequency': [], 'data': [], 'fit_params': [], 'avrg_counts': []}

        # get the frequencices of the sweep
        freq_values, freq_range = self.get_freq_array()
        self.data.update({'frequency': freq_values})
        # initialize data arrays
        esr_data = np.zeros((self.settings['esr_avg'], len(freq_values))) # for the raw esr data

        # to get the timetrace we keep track of the indecies WRONG? 20221115
        if self.settings['save_timetrace']:
            index_data = np.zeros((self.settings['esr_avg'], len(freq_values))).astype(int)  # for the raw esr data

        #WRONG?20221115
        laser_data = np.zeros((self.settings['esr_avg'], len(freq_values))) # for the raw photodiode data
        avrg_counts = np.zeros(self.settings['esr_avg']) # average counts for EACH average to normalize the plot if take_ref is true

        # run sweeps
        for scan_num in range(0, self.settings['esr_avg']):


            esr_data_pos = 0

            if self._abort:
                break

            # get the data for a single sweep. These are raw data.
            self.log('starting average number: ' + str(scan_num) + ' time elapsed: ' + str(time.time()-start_time))
            single_sweep_data, single_sweep_laser_data = self.run_sweep(freq_values)


            # save the single sweep data and normalize to kcounts/sec
            esr_data[scan_num, esr_data_pos:(esr_data_pos + len(single_sweep_data))] = single_sweep_data * (.001 / self.settings['integration_time'])
            laser_data[scan_num, esr_data_pos:(esr_data_pos + len(single_sweep_laser_data))] = single_sweep_laser_data
            esr_data_pos += len(single_sweep_data)

            # average counts of the expt
            avrg_counts[scan_num] = np.mean(esr_data[scan_num])

            if take_ref is True:
                esr_data[scan_num] /=avrg_counts[scan_num]

            # normalize the data if track_laser
            if self.settings['track_laser_power']['on/off']:
                laser_norm_data = np.divide(esr_data, laser_data)

                # average of the normalized data for the number of averages completed so far, to plot and fit to if laser power tracking is on
                data_laser_norm = (np.mean(laser_norm_data[0:(scan_num+1)], axis=0)) #*tmp_laser

            # current non-normalized averaged data to plot and fit to if laser power tracking is off
            esr_avg = np.mean(esr_data[0:(scan_num + 1)] , axis=0)

            if not self.settings['track_laser_power']['on/off']:
                # fit to the data
                fit_params = fit_esr(freq_values, esr_avg, min_counts = self.settings['fit_constants']['minimum_counts'],
                                    contrast_factor=self.settings['fit_constants']['contrast_factor'])
            elif self.settings['track_laser_power']['on/off']:
                # fit to the data
                fit_params = fit_esr(freq_values, data_laser_norm, min_counts = self.settings['fit_constants']['minimum_counts'],
                                    contrast_factor=self.settings['fit_constants']['contrast_factor'])

                # save the data
                self.data.update({'laser_data': laser_data})
                self.data.update({'data_laser_norm': data_laser_norm})

            self.data.update({'data': esr_avg, 'fit_params': fit_params})
       #     self.data.update({'data': avrg_counts})  # ER 20181022

            if self.settings['save_full_esr']:
                self.data.update({'esr_data':esr_data})

            progress = self._calc_progress(scan_num)
            self.updateProgress.emit(progress)

        # turn off the PB card channel: this is so that it's not left on, you run a high amplitude experiment next (e.g. Rabi), and burn the cables or CPW
        # ER 20181017
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

        if self.settings['turn_off_after']:
            self.instruments['microwave_generator']['instance'].update({'enable_output': False})

    def _calc_progress(self, scan_num):
        #COMMENT_ME

        progress = float(scan_num) / self.settings['esr_avg'] * 100.
        self.progress = progress
        return int(progress)

    def _plot(self, axes_list, data = None):
        """
        plotting function for esr
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
            data: data (dictionary that contains keys frequency, data and fit_params) if not provided use self.data
        Returns:

        """

        if data is None:
            data = self.data
        if not self.settings['track_laser_power']['on/off']:
            plot_esr(axes_list[0], data['frequency'], data['data'], data['fit_params'])
        elif self.settings['track_laser_power']['on/off'] and self.settings['daq_type'] == 'PCI':
            plot_esr(axes_list[0], data['frequency'], data['data_laser_norm'], data['fit_params'])

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
        return super(ESR, self).get_axes_layout(new_figure_list)


class ESR_simple_lowerupper(Script):

    _DEFAULT_SETTINGS = [
        Parameter('frequency_1', 2.8e9, float, 'first center frequency (Hz)'),
        Parameter('frequency_2', 3.1e9, float, 'second center frequency (Hz)'),
        Parameter('frequency_width', 100.0e6, float, 'frequency width of ESR scan')
    ]

    _INSTRUMENTS = {}

    _SCRIPTS = {'esr_simple': ESR_simple}


    def _function(self):

        freqs = [self.settings['frequency_1'], self.settings['frequency_2']]
        width = self.settings['frequency_width']

        for freq in freqs:
            self.data = copy(self.scripts['esr_simple'].data)

            self.scripts['esr_simple'].settings['freq_start'] = freq
            self.scripts['esr_simple'].settings['freq_stop'] = width
            self.scripts['esr_simple'].settings['range_type'] = 'center_range'
            self.scripts['esr_simple'].run()
            self.data = copy(self.scripts['esr_simple'].data)

    def _plot(self, axes_list, data = None):
        """
        plotting function for find_nv
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
            data: data (dictionary that contains keys image_data, extent, initial_point, maximum_point) if not provided use self.data
        """

        print('inside _plot!')
        if data is None:
            data = self.data
      #  print(self._current_subscript_stage['current_subscript'])

        if self._current_subscript_stage['current_subscript'] == self.scripts['esr_simple']:
            print('got to here!! ')
            self.scripts['esr_simple']._plot(axes_list)


    def _update_plot(self, axes_list):
        """
        update plotting function for find_nv
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
        """
        print('inside _update_plot!')
    #    print(self._current_subscript_stage['current_subscript'])
        if self._current_subscript_stage['current_subscript'] == self.scripts['esr_simple']:
            print('got to update plot!! ')
            self.scripts['esr_simple']._update_plot(axes_list)


# new version of ESR (JG)
# a
class ESR_2(ESR):
    """
    This class runs ESR on an NV center, outputing microwaves using a MicrowaveGenerator and
    reading in NV counts using a DAQ.
    Each frequency is set explicitly on the SRS, instead of using FM.

    JG (work in progrss)
    This version builds on ESR and allows to track the fluorescence in time

    """

    def run_sweep(self, freq_values):
        '''

        Actually runs the ESR sweep, for a single average.

        Returns:
            esr_data
            laser_data
            normalized_data

        '''

        num_samps = self.settings['num_samps_per_pt'] + 1 # acquire this many samples per point. The counter starts in the middle of a
                                                          # clock tick, so throw out the first sample and take num_samps_per_pt + 1 samples

        # initialize data arrays
        single_sweep_data = np.zeros(len(freq_values))
        single_sweep_laser_data = np.zeros(len(freq_values))

        indeces = list(range(len(freq_values)))
        if self.settings['randomize']:
            random.shuffle(indeces)

        for freq_index in indeces:

            if self._abort:
                break

            freq = freq_values[freq_index]

            # change MW frequency
            self.instruments['microwave_generator']['instance'].update({'frequency': float(freq)})
            time.sleep(self.settings['mw_generator_switching_time'])

            # setup the tasks
            ctrtask = self.daq_in.setup_counter("ctr0", num_samps)
            if self.settings['track_laser_power']['on/off']:
                aitask = self.daq_in.setup_AI(self.settings['track_laser_power']['ai_channel'], num_samps,
                                              continuous=False, clk_source=ctrtask)

            if self.settings['track_laser_power']['on/off']:
                self.daq_in.run(aitask) # AI is actually tied to the clock, when this runs the clock actually starts

            self.daq_in.run(ctrtask) # the counter clock turns on and starts the AI task
            time.sleep(self.settings['integration_time'])

            # read the data
            raw_data, _ = self.daq_in.read_counter(ctrtask)
            single_sweep_data[freq_index] = np.sum(np.diff(raw_data)) # take the total counts, neglecting the first element

            if self.settings['track_laser_power']['on/off']:
                raw_data_laser, _ = self.daq_in.read(aitask)
                single_sweep_laser_data[freq_index] = np.mean(raw_data_laser)

            # clean up APD tasks
            if self.settings['track_laser_power']['on/off']:
                self.daq_in.stop(aitask)  # only stop teh ai task when you've extracted the data you need!!
            self.daq_in.stop(ctrtask)  # stop the clock task last

        return single_sweep_data, single_sweep_laser_data

# re-written by ER 20180904 to check the ESR code.
class ESR_Spec_Ana(Script):
    """
    This class runs ESR, while taking a trace on the spectrum analyzer and finding the peak to check the ESR code.

    """

    _DEFAULT_SETTINGS = [
        Parameter('power_out', -45.0, float, 'output power (dBm)'),
        Parameter('esr_avg', 50, int, 'number of esr averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'], 'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan'),
        Parameter('integration_time', 0.05, float, 'measurement time of fluorescent counts (must be a multiple of settle time)'),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)'),
        Parameter('turn_off_after', False, bool, 'if true MW output is turned off after the measurement'),
        Parameter('norm_to_ref', True, bool, 'If true normalize each frequency sweep by the average counts.'),
        Parameter('save_full_esr', True, bool, 'If true save all the esr traces individually'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan'),
        Parameter('fit_constants',
                  [
                      Parameter('minimum_counts', .5, float, 'minumum counts for an ESR to not be considered noise'),
                      Parameter('contrast_factor', 1.5, float, 'minimum contrast for an ESR to not be considered noise')
                  ]),
        Parameter('track_laser_power',
                  [
                      Parameter('on/off', False, bool, 'If true, measure and normalize out laser power drifts during esr'),
                      Parameter('ai_channel', 'ai4', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4'], 'channel to use for analog input, to which the photodiode is connected')
                  ]),
    ]

    _INSTRUMENTS = {
        'microwave_generator': MicrowaveGenerator,
        'NI6259': NI6259,
        'NI9402': NI9402,
    }

    _SCRIPTS = {'spec_anal_get_trace': SpecAnalyzerGetSpectrum}

    def __init__(self, instruments, scripts = None, name=None, settings=None, log_function=None, data_path = None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path = data_path)
        # defines which daqs contain the input and output based on user selection of daq interface
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_in = self.instruments['NI9402']['instance']

    def setup_microwave_gen(self):
        '''

        set relevant parameters on the MW generator.

        '''
        self.instruments['microwave_generator']['instance'].update({'amplitude': self.settings['power_out']})
        self.instruments['microwave_generator']['instance'].update({'enable_modulation': False})
        self.instruments['microwave_generator']['instance'].update({'enable_output': True})

    def setup_daq(self):
        '''

        Initialize the relevant DAQ and set the sample rate.

        '''
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_in = self.instruments['NI9402']['instance']

        sample_rate = float(2) / self.settings['integration_time'] # DAQ minimum buffer size is 2, so we break the integration time in half
        self.daq_in.settings['digital_input']['ctr0']['sample_rate'] = sample_rate

    def get_freq_array(self):
        '''

        Construct the frequency array based on the setting parameters.

        Returns:
            freq_values: array of the frequencies to be tested
            freq_range: the range of the frequency to be tested, i.e. maximum frequency - minumum frequency
        '''

        # contruct the frequency array
        if self.settings['range_type'] == 'start_stop':
            if self.settings['freq_start']>self.settings['freq_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True

            if self.settings['freq_start'] < 950E3 or self.settings['freq_stop'] > 4.05E9:
                self.log('start or stop frequency out of bounds')
                self._abort = True

            freq_values = np.linspace(self.settings['freq_start'], self.settings['freq_stop'], self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

        elif self.settings['range_type'] == 'center_range':
            if self.settings['freq_start'] < 2 * self.settings['freq_stop']:
                self.log('end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True
            freq_values = np.linspace(self.settings['freq_start']-self.settings['freq_stop']/2,
                                      self.settings['freq_start']+self.settings['freq_stop']/2, self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

            if self.settings['freq_stop'] > 1e9:
                self.log('freq_stop (range) is quite large --- did you mean to set \'range_type\' to \'start_stop\'? ')
        else:
            self.log('unknown range parameter. Abort script')
            self._abort = True

        return freq_values, freq_range

    def run_sweep(self, freq_values):
        '''

        Actually runs the ESR sweep, for a single average.

        Returns:
            esr_data
            laser_data
            normalized_data

        '''

        num_samps = 2 # acquire 1 sample per point

        # initialize data arrays
        single_sweep_data = np.zeros(len(freq_values))
        single_sweep_laser_data = np.zeros(len(freq_values))
        single_sweep_sa_data = np.zeros(len(freq_values))

        freq_index = 0

        for freq in freq_values:

            if self._abort:
                break

            # change MW frequency
            self.instruments['microwave_generator']['instance'].update({'frequency': float(freq)})
            time.sleep(self.settings['mw_generator_switching_time'])

            # setup the tasks
            ctrtask = self.daq_in.setup_counter("ctr0", num_samps)
            if self.settings['track_laser_power']['on/off']:
                aitask = self.daq_in.setup_AI(self.settings['track_laser_power']['ai_channel'], num_samps,
                                              continuous=False, clk_source=ctrtask)

            if self.settings['track_laser_power']['on/off']:
                self.daq_in.run(aitask) # AI is actually tied to the clock, when this runs the clock actually starts

            self.daq_in.run(ctrtask) # the counter clock turns on and starts the AI task

            #todo: talk to aaron about waitToFinish timeout duration
            #if self.settings['track_laser_power']['on/off']:
                #self.daq_in.waitToFinish(aitask) # wait for the tasks to be done
            # read the data
            raw_data, _ = self.daq_in.read_counter(ctrtask)
            single_sweep_data[freq_index] = raw_data[1]/2 # take last element since counter is cumulative, then divide by 2 to get mean

            if self.settings['track_laser_power']['on/off']:
                raw_data_laser, _ = self.daq_in.read(aitask)
                single_sweep_laser_data[freq_index] = np.mean(raw_data_laser)

            # clean up APD tasks
            if self.settings['track_laser_power']['on/off']:
                self.daq_in.stop(aitask)  # only stop teh ai task when you've extracted the data you need!!
            self.daq_in.stop(ctrtask)  # stop the clock task last

            self.scripts['spec_anal_get_trace'].run()
            single_sweep_sa_data[freq_index] = self.scripts['spec_anal_get_trace'].data['peak_freq']

            freq_index += 1

        return single_sweep_data, single_sweep_laser_data, single_sweep_sa_data

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        # if tracking laser power drifts, check for PCI daq
        if self.settings['track_laser_power']['on/off'] and not self.settings['daq_type'] == 'PCI':
            print("tracking laser power drifts only enabled for PCI daq")
            self._abort = True

        self.lines = []
        take_ref = self.settings['norm_to_ref']

        # setup the daq
        self.setup_daq()

        # setup the microwave generator
        self.setup_microwave_gen()

        # intialize some of the fields in self.data
        self.data = {'frequency': [], 'data': [], 'fit_params': [], 'avrg_counts' : [], 'peaks': [], 'temp_peaks': []}

        # get the frequencices of the sweep
        freq_values, freq_range = self.get_freq_array()
        self.data.update({'frequency': freq_values})
        # initialize data arrays
        esr_data = np.zeros((self.settings['esr_avg'], len(freq_values))) # for the raw esr data
        laser_data = np.zeros((self.settings['esr_avg'], len(freq_values))) # for the raw photodiode data
        spec_anal_peaks = np.zeros((self.settings['esr_avg'], len(freq_values)))
        avrg_counts = np.zeros(self.settings['esr_avg']) # average counts for EACH average to normalize the plot if take_ref is true

        # run sweeps
        for scan_num in range(0, self.settings['esr_avg']):

            esr_data_pos = 0

            if self._abort:
                break

            # get the data for a single sweep
            single_sweep_data, single_sweep_laser_data, single_sweep_sa_data = self.run_sweep(freq_values)

            # save the single sweep data and normalize to kcounts/sec
            esr_data[scan_num, esr_data_pos:(esr_data_pos + len(single_sweep_data))] = single_sweep_data * (.001 / self.settings['integration_time'])
            laser_data[scan_num, esr_data_pos:(esr_data_pos + len(single_sweep_laser_data))] = single_sweep_laser_data
            spec_anal_peaks[scan_num, esr_data_pos:(esr_data_pos + len(single_sweep_laser_data))] = single_sweep_sa_data
            esr_data_pos += len(single_sweep_data)

            # average counts of the expt
            avrg_counts[scan_num] = np.mean(esr_data[scan_num])

            if take_ref is True:
                esr_data[scan_num] /=avrg_counts[scan_num]

            # normalize the data if track_laser
            if self.settings['track_laser_power']['on/off']:
                laser_norm_data = np.divide(esr_data, laser_data)

                # average of the normalized data for the number of averages completed so far, to plot and fit to if laser power tracking is on
                tmp_laser = np.mean(np.mean(laser_data[scan_num])) # instantaneous average of the laser power
                data_laser_norm = (np.mean(laser_norm_data[0:(scan_num+1)], axis=0))*tmp_laser

            # current non-normalized averaged data to plot and fit to if laser power tracking is off
            esr_avg = np.mean(esr_data[0:(scan_num + 1)] , axis=0)

            if not self.settings['track_laser_power']['on/off']:
                # fit to the data
                fit_params = fit_esr(freq_values, esr_avg, min_counts = self.settings['fit_constants']['minimum_counts'],
                                    contrast_factor=self.settings['fit_constants']['contrast_factor'])
            elif self.settings['track_laser_power']['on/off']:
                # fit to the data
                fit_params = fit_esr(freq_values, data_laser_norm, min_counts = self.settings['fit_constants']['minimum_counts'],
                                    contrast_factor=self.settings['fit_constants']['contrast_factor'])

                # save the data
                self.data.update({'laser_data': laser_data})
                self.data.update({'data_laser_norm': data_laser_norm})

            self.data.update({'data': esr_avg, 'fit_params': fit_params})
            self.data.update({'peaks': spec_anal_peaks})
            self.data.update({'temp_peaks': single_sweep_sa_data})
            self.data.update({'avrg_counts': avrg_counts})

            if self.settings['save_full_esr']:
                self.data.update({'esr_data': esr_data})

            progress = self._calc_progress(scan_num)
            self.updateProgress.emit(progress)

        if self.settings['turn_off_after']:
            self.instruments['microwave_generator']['instance'].update({'enable_output': False})

    def _calc_progress(self, scan_num):
        #COMMENT_ME

        progress = float(scan_num) / self.settings['esr_avg'] * 100.
        self.progress = progress

        return int(progress)

    def _plot(self, axes_list, data = None):
        """
        plotting function for esr
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
            data: data (dictionary that contains keys frequency, data and fit_params) if not provided use self.data
        Returns:

        """

        if data is None:
            data = self.data
 #       if not self.settings['track_laser_power']['on/off']:
 #           plot_esr(axes_list[0], data['frequency'], data['data'], data['fit_params'])
 #       elif self.settings['track_laser_power']['on/off'] and self.settings['daq_type'] == 'PCI':
 #           plot_esr(axes_list[0], data['frequency'], data['norm_data'], data['fit_params'])
        plot_diff_freq_vs_freq(axes_list[0], data['frequency'] - data['temp_peaks'], data['frequency'])

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
        return super(ESR_Spec_Ana, self).get_axes_layout(new_figure_list)

class NMR_CW(ESR_simple):
    """
    This class runs ESR on an NV center, outputing microwaves using a MicrowaveGenerator and reading in NV counts using
    a DAQ. Each frequency is set explicitly on the SRS, instead of using FM.

    This is the simple, i.e. minimal version of it
    """

    _DEFAULT_SETTINGS = [
        Parameter('power_out', -45.0, float, 'output power (dBm)'),
        Parameter('esr_avg', 1, int, 'number of esr averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'], 'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan'),
        Parameter('integration_time', 0.05, float, 'The TOTAL integration time. The integration time per daq bins is integration_time / num_samps_per_pt'),
        Parameter('num_samps_per_pt', 100, int, 'Number of samples within one DAQ buffer, for each frequency. This should be large because the first measurement is thrown away (may start in the middle of a clock tick)'),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)'),
        Parameter('turn_off_after', False, bool, 'if true MW output is turned off after the measurement'),
        Parameter('save_full_esr', True, bool, 'If true save all the esr traces individually'),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan'),
        Parameter('fit_constants',
                  [
                      Parameter('minimum_counts', .5, float, 'minumum counts for an ESR to not be considered noise'),
                      Parameter('contrast_factor', 1.5, float, 'minimum contrast for an ESR to not be considered noise')
                  ]),
        Parameter('randomize', True, bool, 'check to randomize esr frequencies'),
        Parameter('save_timetrace', True, bool, 'check to save the measured fluorescence over time. This is identical to the full esr when the freq. are not randomized')
    ]

    _INSTRUMENTS = {
        'rf_generator': RFGenerator,
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

    def setup_microwave_gen(self):
        '''

        set relevant parameters on the MW generator.

        '''
        self.instruments['rf_generator']['instance'].update({'amplitude_rf': self.settings['power_out']})
        self.instruments['rf_generator']['instance'].update({'enable_modulation': False})
        self.instruments['rf_generator']['instance'].update({'enable_rf_output': True})

    def setup_daq(self):
        '''

        Initialize the relevant DAQ and set the sample rate.

        '''
        if self.settings['daq_type'] == 'PCI':
            self.daq_in = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_in = self.instruments['NI9402']['instance']

        sample_rate = float(1) / (self.settings['integration_time']/self.settings['num_samps_per_pt']) # DAQ minimum buffer size is 2
        self.daq_in.settings['digital_input']['ctr0']['sample_rate'] = sample_rate

    def get_freq_array(self):
        '''

        Construct the frequency array based on the setting parameters.

        Returns:
            freq_values: array of the frequencies to be tested
            freq_range: the range of the frequency to be tested, i.e. maximum frequency - minumum frequency
        '''

        # contruct the frequency array
        if self.settings['range_type'] == 'start_stop':
            if self.settings['freq_start']>self.settings['freq_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True

            #if self.settings['freq_start'] < 950E3 or self.settings['freq_stop'] > 4.05E9: # freq range of the SRS
            #    self.log('start or stop frequency out of bounds')
            #    self._abort = True

            freq_values = np.linspace(self.settings['freq_start'], self.settings['freq_stop'], self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

        elif self.settings['range_type'] == 'center_range':
            if self.settings['freq_start'] < 2 * self.settings['freq_stop']:
                self.log('end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True
            freq_values = np.linspace(self.settings['freq_start']-self.settings['freq_stop']/2,
                                      self.settings['freq_start']+self.settings['freq_stop']/2, self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

            if self.settings['freq_stop'] > 1e9:
                self.log('freq_stop (range) is quite large --- did you mean to set \'range_type\' to \'start_stop\'? ')
        else:
            self.log('unknown range parameter. Abort script')
            self._abort = True

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

        indeces = list(range(len(freq_values)))
        if self.settings['randomize']:
            random.shuffle(indeces)

        for freq_index in indeces:

            single_freq_start_t = time.time()

            if self._abort:
                break

            freq = freq_values[freq_index]

            # change MW frequency
            self.instruments['rf_generator']['instance'].update({'frequency': float(freq)})
            time.sleep(self.settings['mw_generator_switching_time'])


            single_sweep_data[freq_index] = self.measure_signal(num_samps)


        # normalize  single sweep data to kcounts/sec
        single_sweep_data = single_sweep_data * (.001 / self.settings['integration_time'])

        end_run_sweep = time.time()


        return single_sweep_data, indeces

    def measure_signal(self, num_samps):
        """

        measure the signal with the APD

        Args:
            num_samps:

        Returns:

        """

        start_meas = time.time()

        # setup the tasks
        ctrtask = self.daq_in.setup_counter("ctr0", num_samps)

        t2 = time.time()-start_meas

        self.daq_in.run(ctrtask)  # the counter clock turns on and starts the AI task


        time.sleep(self.settings['integration_time'])

        t4 = time.time()-start_meas

        # read the data
        raw_data, _ = self.daq_in.read_counter(ctrtask)

        signal = np.sum(np.diff(raw_data))  # take the total counts, neglecting the first element

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
        self.setup_microwave_gen()

        # intialize some of the fields in self.data
        self.data = {'frequency': [], 'data': [], 'fit_params': [], 'avrg_counts': []}

        # get the frequencices of the sweep
        freq_values, freq_range = self.get_freq_array()
        self.data.update({'frequency': freq_values})
        # initialize data arrays
        esr_data = np.zeros((self.settings['esr_avg'], len(freq_values))) # for the raw esr data

        # to get the timetrace we keep track of the indecies
        if self.settings['save_timetrace']:
            index_data = np.zeros((self.settings['esr_avg'], len(freq_values))).astype(int)  # for the raw esr data

        # run sweeps
        for scan_num in range(0, self.settings['esr_avg']):

            self.log('starting average ' + str(scan_num) + ', time elapsed: ' + str(np.round(time.time()-start_time, 1)) + 's')

            # get the data for a single sweep. These are raw data.
            single_sweep_data, indeces = self.run_sweep(freq_values)

            if self._abort:
                break

            # save the single sweep data
            esr_data[scan_num] = single_sweep_data

            if self.settings['save_timetrace']:
                index_data[scan_num] = indeces

            # current non-normalized averaged data to plot and fit to if laser power tracking is off
            esr_avg = np.mean(esr_data[0:(scan_num + 1)], axis=0)

            # fit to the data
            fit_params = fit_esr(freq_values, esr_avg, min_counts = self.settings['fit_constants']['minimum_counts'],
                                contrast_factor=self.settings['fit_constants']['contrast_factor'])

            self.data.update({'data': esr_avg, 'fit_params': fit_params})

            progress = self._calc_progress(scan_num)
            self.updateProgress.emit(progress)

        if self.settings['save_full_esr']:
            self.data.update({'esr_data': esr_data})

        if self.settings['save_timetrace']:
            # using the indeces, we retrieve the right time ordering for the esr data and save it as the timetrace data
            self.data.update({'index_data':index_data})
            # self.data.update({'time_trace':[esr[idx] for esr, idx in zip(esr_data, index_data)]})


        if self.settings['turn_off_after']:
            self.instruments['rf_generator']['instance'].update({'enable_rf_output': False})


if __name__ == '__main__':
    script = {}
    instr = {}
    script, failed, instr = Script.load_and_append({'ESR': 'ESR'}, script, instr)

    print(script)
    print(failed)
    print(instr)