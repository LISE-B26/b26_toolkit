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

from pylabcontrol.core import Script, Parameter

# import standard libraries
import numpy as np
from b26_toolkit.instruments import MicrowaveGenerator, NI6259, NI9263, NI9402
from b26_toolkit.plotting.plots_1d import plot_esr

from b26_toolkit.data_processing.esr_signal_processing import fit_esr
import time

class ESR_FM_Dither(Script):
    """
    This class runs ESR on an NV center, outputing microwaves using a MicrowaveGenerator and reading in NV counts using
    a DAQ. It uses FM using AO2 on the DAQ, which is off by a few MHz but may be faster than the other ESR script.
    """

    _DEFAULT_SETTINGS = [
        Parameter('power_out', -45.0, float, 'output power (dBm)'),
        Parameter('esr_avg', 50, int, 'number of esr averages'),
        Parameter('freq_one', 2.82e9, float, 'start frequency of scan'),
        Parameter('freq_two', 2.92e9, float, 'end frequency of scan'),
        Parameter('reps_per_average', 100, int, 'number of points per frequency per average'),
        Parameter('integration_time', 90, float, 'measurement time of fluorescent counts (must be a multiple of settle time) in us'),
        Parameter('settle_time', 10, float, 'time wait after changing frequencies using daq (us)'),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)'),
        Parameter('turn_off_after', False, bool, 'if true MW output is turned off after the measurement'),
        Parameter('take_ref', True, bool, 'If true normalize each frequency sweep by the average counts. This should be renamed at some point because now we dont take additional data for the reference.'),
        Parameter('save_full_esr', True, bool, 'If true save all the esr traces individually'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan'),
        Parameter('FM_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3']),
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
            aotask = self.daq_out.setup_AO([self.settings['FM_channel']], freq_voltage_array, ctrtask)

            if self.settings['track_laser_power']['on/off'] == True and self.settings['daq_type'] == 'PCI':
                aitask = self.daq_in.setup_AI(self.settings['track_laser_power']['ai_channel'], len(freq_voltage_array), continuous=False, clk_source=ctrtask) # for optional arguments spell out every one if there are multiple
            elif self.settings['track_laser_power']['on/off'] == True and self.settings['daq_type'] != 'PCI':
                raise(NotImplementedError('cant use laser power tracking without the PCI daq!!'))

            # start counter and scanning sequence
            if self.settings['track_laser_power']['on/off'] == True and self.settings['daq_type'] == 'PCI':
                self.daq_in.run(aitask) # AI is actually tied to the clock, this runs when the clock actually starts
            self.daq_out.run(aotask) # AO is actually tied to the clock, this runs when the clock actually starts
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

        freq_values = [self.settings['freq_one'], self.settings['freq_two']]

        if(np.abs(self.settings['freq_one'] - self.settings['freq_two']) > 64e6):
            print('Two frequencies used must be between 64 MHz')
            return
        clock_adjust = int((self.settings['integration_time'] + self.settings['settle_time']) / self.settings['settle_time'])
        freq_array = np.tile(np.repeat(freq_values, clock_adjust), self.settings['reps_per_average'])
        self.instruments['microwave_generator']['instance'].update({'amplitude': self.settings['power_out']})
        self.instruments['microwave_generator']['instance'].update({'modulation_type': 'FM'})

        # ER 20171128
        self.instruments['microwave_generator']['instance'].update({'dev_width': self.instruments['microwave_generator']['instance'].settings['dev_width']})
        self.instruments['microwave_generator']['instance'].update({'enable_modulation': True})

        sample_rate = float(1) / (self.settings['settle_time']*(1e-6))
        self.daq_out.settings['analog_output'][self.settings['FM_channel']]['sample_rate'] = sample_rate
        self.daq_in.settings['digital_input']['ctr0']['sample_rate'] = sample_rate

        self.instruments['microwave_generator']['instance'].update({'enable_output': True})

        esr_data = np.zeros((self.settings['esr_avg'], len(freq_values)*self.settings['reps_per_average']))
        full_laser_data = np.zeros((self.settings['esr_avg'], len(freq_values)*self.settings['reps_per_average']))
        full_normalized_data = np.zeros((self.settings['esr_avg'], len(freq_values)*self.settings['reps_per_average']))
        avrg_counts = np.zeros(self.settings['esr_avg']) # here we save the avrg of the esr scan which we will use to normalize
        norm_data = np.zeros(len(freq_values))
        self.data = {'frequency': [], 'data': [], 'fit_params': [], 'avrg_counts' : avrg_counts}

        # run sweeps
        for scan_num in range(0, self.settings['esr_avg']):
            if self._abort:
                break
            esr_data_pos = 0

            freq_voltage_array, center_freq = get_frequency_voltages(freq_values,
                                                                     0,
                                                                     self.instruments['microwave_generator'][
                                                                         'instance'].settings['dev_width'],
                                                                     freq_array)

            summed_data, normalized_data, laser_data = read_freq_section(freq_voltage_array, center_freq, clock_adjust)

            # also normalizing to kcounts/sec
            esr_data[scan_num, esr_data_pos:(esr_data_pos + len(summed_data))] = summed_data * (.001 / (self.settings['integration_time'] * (1e-6)))
            full_laser_data[scan_num, esr_data_pos:(esr_data_pos + len(laser_data))] = laser_data
            full_normalized_data[scan_num, esr_data_pos:(esr_data_pos + len(laser_data))] = normalized_data * (.001 / (self.settings['integration_time'] * (1e-6)))
            esr_data_pos += len(summed_data)


            avrg_counts[scan_num] = np.mean(esr_data[scan_num])
            norm_data = np.mean(full_normalized_data[0:(scan_num + 1)] , axis=0)

            esr_avg = np.mean(esr_data[0:(scan_num + 1)] , axis=0)
            # if not self.settings['track_laser_power']['on/off']:
            #     fit_params = fit_esr(freq_values, esr_avg, min_counts = self.settings['fit_constants']['minimum_counts'],
            #                         contrast_factor=self.settings['fit_constants']['contrast_factor'])
            # elif self.settings['track_laser_power']['on/off'] == True and self.settings['daq_type'] == 'PCI':
            #     fit_params = fit_esr(freq_values, norm_data, min_counts = self.settings['fit_constants']['minimum_counts'],
            #                         contrast_factor=self.settings['fit_constants']['contrast_factor'])
            # self.data.update({'frequency': freq_values, 'data': esr_avg, 'fit_params': fit_params})

            self.data.update({'frequency': np.tile(freq_values, self.settings['reps_per_average']), 'data': esr_avg})


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
            plot_esr(axes_list[0], data['frequency'], data['data'], linestyle = 'None', marker = 'o')
        elif self.settings['track_laser_power']['on/off'] and self.settings['daq_type'] == 'PCI':
            plot_esr(axes_list[0], data['frequency'], data['norm_data'], linestyle = 'None', marker = 'o')


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
        return super(ESR_FM_Dither, self).get_axes_layout(new_figure_list)

if __name__ == '__main__':
    script = {}
    instr = {}
    script, failed, instr = Script.load_and_append({'ESR': 'ESR_dithering'}, script, instr)

    print(script)
    print(failed)
    print(instr)