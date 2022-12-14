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
from scipy.ndimage.filters import uniform_filter1d
import matplotlib.pyplot as plt
from b26_toolkit.instruments import NI6259, NI9402,NI6229, NI9219, MicrowaveGenerator, PiezoController
from b26_toolkit.plotting.plots_1d import plot_counts, update_1d_simple, update_counts_vs_pos, update_counts
from pylabcontrol.core import Parameter, Script
from b26_toolkit.scripts import FindNV


class Daq_Read_Counter(Script):
    """
This script reads the Counter input from the DAQ and plots it.

WARNING: Only implemented either for the PCI DAQ (NI6259) or cDAQ (NI9402) !!!!

If you want to use it make sure that the right instrument is defined in _INSTRUMENTS = {'daq': NI9402} in the python code.

    """
    _DEFAULT_SETTINGS = [
        Parameter('integration_time', .25, float, 'Time per data point (s)'),
        Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr2'], 'Daq channel used for counter'),
        Parameter('total_int_time', 3.0, float, 'Total time to integrate (s) (if -1 then it will go indefinitely)'), # added by ER 20180606
        Parameter('track_laser_power_photodiode1',
                  [
                      Parameter('on/off', False, bool,
                                'If true, measure and normalize out laser power drifts during daq_read_counter'),
                      Parameter('ai_channel', 'ai2', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4'],
                                'channel to use for analog input, to which the photodiode is connected')
                  ]),
        Parameter('track_laser_power_photodiode2',
                  [
                      Parameter('on/off', False, bool, 'If true, measure and save laser power drifts during daq_read_counter on this photodiode. Cant use both simultaneously'),
                      Parameter('ai_channel', 'ai4', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4'], 'channel to use for photodiode 2, cant be the same as the track_laser_power photodiode')
                  ])
    ]
    _INSTRUMENTS = {'daq': NI6229}
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

        self.data = {'counts': deque(), 'laser_power': deque(), 'normalized_counts': deque(), 'laser_power2': deque()}

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        if self.settings['track_laser_power_photodiode1']['on/off'] and self.settings['track_laser_power_photodiode2']['on/off']:
            print('cant use both photodiodes at the same time - only use one AI channel at a time, unfortunately :-(')
            return

        sample_rate = float(2) / self.settings['integration_time']
        normalization = self.settings['integration_time']/.001
        self.instruments['daq']['instance'].settings['digital_input'][self.settings['counter_channel']]['sample_rate'] = sample_rate
        self.data = {'counts': deque(), 'laser_power': deque(), 'normalized_counts': deque(), 'laser_power2': deque()}
        self.last_value = 0
        sample_num = 2

        task = self.instruments['daq']['instance'].setup_counter(self.settings['counter_channel'], sample_num, continuous_acquisition=True)

        if self.settings['track_laser_power_photodiode1']['on/off'] == True:
            aitask = self.instruments['daq']['instance'].setup_AI(self.settings['track_laser_power_photodiode1']['ai_channel'], sample_num,
                                          continuous=True, # continuous sampling still reads every clock tick, here set to the clock of the counter
                                          clk_source=task)

        if self.settings['track_laser_power_photodiode2']['on/off'] == True:
            aitask2 = self.instruments['daq']['instance'].setup_AI(self.settings['track_laser_power_photodiode2']['ai_channel'], sample_num,
                                          continuous=True, # continuous sampling still reads every clock tick, here set to the clock of the counter
                                          clk_source=task)
            print('aitask2: ', aitask2)

        # maximum number of samples if total_int_time > 0
        if self.settings['total_int_time'] > 0:
            max_samples = np.floor(self.settings['total_int_time']/self.settings['integration_time'])

        # start counter and scanning sequence
        if (self.settings['track_laser_power_photodiode1']['on/off'] and not self.settings['track_laser_power_photodiode2']['on/off']):
            self.instruments['daq']['instance'].run(aitask)
        elif (self.settings['track_laser_power_photodiode2']['on/off'] and not self.settings['track_laser_power_photodiode1']['on/off']):
            self.instruments['daq']['instance'].run(aitask2)

        self.instruments['daq']['instance'].run(task)

        # ER 20180827 wait for at least one clock tick to go by to start with a full clock tick of acquisition time for the first bin
        time.sleep(self.settings['integration_time'])

        sample_index = 0 # keep track of samples made to know when to stop if finite integration time

        while True:
            if self._abort:
                break

            # TODO: this is currently a nonblocking read so we add a time.sleep at the end so it doesn't read faster
            # than it acquires, this should be replaced with a blocking read in the future
            if self.settings['track_laser_power_photodiode1']['on/off'] == True:
                raw_data_laser, num_read_laser = self.instruments['daq']['instance'].read(aitask)
            if self.settings['track_laser_power_photodiode2']['on/off'] == True:
                raw_data_laser2, num_read_laser2 = self.instruments['daq']['instance'].read(aitask2)

            raw_data, num_read = self.instruments['daq']['instance'].read(task)

            #skip first read, which gives an anomolous value
            if num_read.value == 1:
                self.last_value = raw_data[0] #update running value to last measured value to prevent count spikes
                time.sleep(2.0 / sample_rate)
                continue

            tmp_count = 0
            for value in raw_data:
                new_val = ((float(value) - self.last_value) / normalization)
                self.data['counts'].append(new_val)
                self.last_value = value
                if self.settings['track_laser_power_photodiode1']['on/off'] == True:
                    self.data['laser_power'].append(raw_data_laser[tmp_count])
                if self.settings['track_laser_power_photodiode2']['on/off'] == True:
                    self.data['laser_power2'].append(raw_data_laser2[tmp_count])

                tmp_count = tmp_count + 1

            if self.settings['total_int_time'] > 0:
                self.progress = sample_index/max_samples
            else:
                self.progress = 50.
            self.updateProgress.emit(int(self.progress))

            time.sleep(2.0 / sample_rate)
            sample_index = sample_index + 1
            if self.settings['total_int_time'] > 0. and sample_index >= max_samples: # if the maximum integration time is hit
                self._abort = True # tell the script to abort

        # clean up APD tasks
        self.instruments['daq']['instance'].stop(task)
        if self.settings['track_laser_power_photodiode1']['on/off'] == True:
            self.instruments['daq']['instance'].stop(aitask)
        if self.settings['track_laser_power_photodiode2']['on/off'] == True:
            self.instruments['daq']['instance'].stop(aitask2)

        self.data['counts'] = list(self.data['counts'])

        if self.settings['track_laser_power_photodiode1']['on/off'] == True:
            self.data['laser_power'] = list(self.data['laser_power'])
            self.data['normalized_counts'] = list(np.divide(np.multiply(self.data['counts'], np.mean(self.data['laser_power'])), self.data['laser_power']))
        if self.settings['track_laser_power_photodiode2']['on/off'] == True:
            self.data['laser_power2'] = list(self.data['laser_power2'])

    def plot(self, figure_list):
        super(Daq_Read_Counter, self).plot([figure_list[1]])

    def _plot(self, axes_list, data = None):
        # COMMENT_ME

        if data is None:
            data = self.data

        if len(data['counts']) > 0:
            if self.settings['track_laser_power_photodiode1']['on/off'] == True:
                array_to_plot = np.delete(np.divide(np.multiply(self.data['counts'], np.mean(self.data['laser_power'])), self.data['laser_power']),0)
            else:
                array_to_plot = np.delete(data['counts'], 0)

            plot_counts(axes_list[0], array_to_plot)

    def _update_plot(self, axes_list, data = None):
        if data is None:
            data = self.data

        if data:
            if self.settings['track_laser_power_photodiode1']['on/off'] == True:
                array_to_plot = np.delete(np.divide(np.multiply(self.data['counts'], np.mean(self.data['laser_power'])), self.data['laser_power']), 0)
            else:
                array_to_plot = np.delete(data['counts'], 0)

            update_counts_vs_pos(axes_list[0], array_to_plot, np.linspace(0, len(array_to_plot), len(array_to_plot)))

class Daq_Read_Counter_2_Channel(Script):
    """
This script reads the Counter input from the DAQ and plots it.

WARNING: Only implemented either for the PCI DAQ (NI6259) or cDAQ (NI9402) !!!!

If you want to use it make sure that the right instrument is defined in _INSTRUMENTS = {'daq': NI9402} in the python code.

    """
    _DEFAULT_SETTINGS = [
        Parameter('integration_time', .25, float, 'Time per data point (s)'),
        Parameter('counter_channel_a', 'ctr0', ['ctr0', 'ctr2'], 'Daq channel used for counter'),
        Parameter('counter_channel_b', 'ctr2', ['ctr0', 'ctr2'], 'Daq channel used for counter'),
        Parameter('total_int_time', -1, float, 'Total time to integrate (s) (if -1 then it will go indefinitely)'),
        Parameter('trim_plot', -1, int, 'Keep only the last n data points to keep plot decluttered (saved data unaffected), (if -1 then no trim)'),
        Parameter('microwaves',
                  [Parameter('power_out', -45.0, float, 'output power (dBm) (NOT implemented)'),
                   Parameter('freq_mod_rate', 1e3, float, 'frequency modulation rate (NOT implemented)'),
                   Parameter('freq_mod_dev', 1e6, float, 'frequency modulation width (NOT implemented)'),
                   Parameter('tracking',
                             [Parameter('tracking_on', False, bool, 'feedback on normalized count difference to locate ESR dip'),
                              Parameter('freq_start', 2.87e9, float, 'initial guess for ESR dip'),
                              Parameter('P', 6e6, float, 'initial guess for ESR dip'),
                              Parameter('I', 0, float, 'initial guess for ESR dip'),
                              Parameter('low-pass cutoff', -1, float, 'frequency of low-pass filter on error signal (if -1 then disable filter)')
                              ]),
                   Parameter('attocube_feedback',
                             [Parameter('feedback_on', False, bool,
                                        'move Attocube to eliminate count difference and maintain specified ESR frequency'),
                              Parameter('freq_target', 2.87e9, float, 'target ESR frequency'),
                              Parameter('gradient', 6e6, float, 'gradient at current point (T/m)'),
                              Parameter('attocube_calibration', 35, float, 'Attocube calibration (nm/V)')
                              ])
                   ]),
        Parameter('nv_tracking', [
            Parameter('on/off', False, bool, 'used to turn on tracking'),
            Parameter('threshold', 0.85, float, 'threshold for tracking'),
            Parameter('init_fluor', 20., float, 'initial fluorescence of the NV to compare to, in kcps')
        ])
    ]

    _INSTRUMENTS = {'daq': NI9402,  'microwave_generator': MicrowaveGenerator, 'piezo_controller': PiezoController}

    _SCRIPTS = {'find_nv': FindNV}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.data = {'counts_a': deque(), 'counts_b': deque(), 'counts_diff': deque(), 'mw_freq': deque(), 'attocube_voltage': deque()}

        plt.ioff()

    def low_pass(self, x_new, y_old, cutoff=20):
        # Simple low-pass filter
        alpha = self.settings['integration_time'] / (self.settings['integration_time'] + 1 / (2 * np.pi * cutoff))
        y_new = x_new * alpha + (1 - alpha) * y_old
        return y_new

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        sample_num = 2
        sample_rate = float(sample_num) / self.settings['integration_time']
        normalization = self.settings['integration_time']/.001/float(sample_num)

        self.instruments['daq']['instance'].settings['digital_input'][self.settings['counter_channel_a']]['sample_rate'] = sample_rate
        self.instruments['daq']['instance'].settings['digital_input'][self.settings['counter_channel_b']]['sample_rate'] = sample_rate

        mw_gen = self.instruments['microwave_generator']['instance']
        self.instruments['piezo_controller']['instance'].axis = 'z'
        self.attocube_voltage = self.instruments['piezo_controller']['instance'].read_probes('voltage')
        self.attocube_voltage_initial = self.attocube_voltage

        self.data = {'counts_a': deque(), 'counts_b': deque(), 'counts_diff': deque(), 'mw_freq': deque(), 'attocube_voltage': deque()}
        self.last_value_a = 0
        self.last_value_b = 0

        # maximum number of samples if total_int_time > 0
        if self.settings['total_int_time'] > 0:
            max_samples = np.floor(self.settings['total_int_time']/self.settings['integration_time'])

        self.task_a = self.instruments['daq']['instance'].setup_counter(self.settings['counter_channel_a'], sample_num,
                                                                   continuous_acquisition=True)
        # Task b only requires a counter and uses the clock from task a. .run() is used to start the clock; the
        # counter starts with just setup_counter(), so there's no need to call .run() here
        existing_clock_channel = self.instruments['daq']['instance'].settings['digital_input'][self.settings['counter_channel_a']]['clock_counter_channel']
        self.task_b = self.instruments['daq']['instance'].setup_counter(self.settings['counter_channel_b'], sample_num,
                                                                   continuous_acquisition=True,existing_clock_channel=existing_clock_channel)
        self.instruments['daq']['instance'].run(self.task_a)


        # ER 20180827 wait for at least one clock tick to go by to start with a full clock tick of acquisition time for the first bin
        time.sleep(self.settings['integration_time'])

        err_sum = 0
        mw_freq = self.settings['microwaves']['tracking']['freq_start']
        if self.settings['microwaves']['tracking']['tracking_on']:
            mw_gen.update({'enable_modulation': True})
            mw_gen.update({'frequency': mw_freq})

        self.threshold_counts = self.settings['nv_tracking']['init_fluor'] * self.settings['nv_tracking']['threshold']
        sample_index = 0  # keep track of samples made to know when to stop if finite integration time
        self.loop_iteration = 0

        while True:
            if self._abort:
                break
            self.loop_iteration += 1
            # TODO: this is currently a nonblocking read so we add a time.sleep at the end so it doesn't read faster
            # than it acquires, this should be replaced with a blocking read in the future

            raw_data_a, num_read_a = self.instruments['daq']['instance'].read(self.task_a)
            raw_data_b, num_read_b = self.instruments['daq']['instance'].read(self.task_b)
            #print(num_read_a.value, list(raw_data_a))
            #print(num_read_b.value, list(raw_data_b))
            #print()
            #skip first read, which gives an anomolous value

            if num_read_a.value < sample_num:
                #self.last_value_a = raw_data_a[0] #update running value to last measured value to prevent count spikes
                time.sleep(float(sample_num) / sample_rate)
                continue

            if num_read_b.value < sample_num:
                #self.last_value_b = raw_data_b[0] #update running value to last measured value to prevent count spikes
                time.sleep(float(sample_num) / sample_rate)
                continue

            if self.last_value_a == 0 or self.last_value_b == 0:
                self.last_value_a = raw_data_a[0]
                raw_data_a = raw_data_a[1:]
                self.last_value_b = raw_data_b[0]
                raw_data_b = raw_data_b[1:]


            for value in raw_data_a:
                self.new_val_a = ((float(value) - self.last_value_a) / normalization)
                #print('New_value:'+str(float(value) - self.last_value_a))
                #print('Last_value:' + str(value))
                self.data['counts_a'].append(self.new_val_a)
                self.last_value_a = value

            for value in raw_data_b:
                self.new_val_b = ((float(value) - self.last_value_b) / normalization)
                self.data['counts_b'].append(self.new_val_b)
                self.last_value_b = value

            new_val_diff = 2 * float(self.new_val_a - self.new_val_b) / float(self.new_val_a + self.new_val_b)
            if len(self.data['counts_diff']) > 1 and self.settings['microwaves']['tracking']['low-pass cutoff'] != -1:
                new_val_diff = self.low_pass(new_val_diff, self.data['counts_diff'][-1],
                                             cutoff=self.settings['microwaves']['tracking']['low-pass cutoff'])
            self.data['counts_diff'].append(new_val_diff)

            # Adjust applied MW frequency to track ESR frequency
            if self.settings['microwaves']['tracking']['tracking_on']:
                err_sum += new_val_diff
                mw_freq_adjust = self.settings['microwaves']['tracking']['P'] * new_val_diff + \
                                 self.settings['microwaves']['tracking']['I'] * err_sum
                if np.abs(mw_freq_adjust) < 6e6:
                    mw_freq = float(mw_freq - mw_freq_adjust)
                    mw_gen.update({'frequency': mw_freq})
            self.data['mw_freq'].append(mw_freq)

            # Adjust Attocube to eliminate deviation in ESR frequency from target
            if self.settings['microwaves']['attocube_feedback']['feedback_on'] and self.loop_iteration % 4 == 0:
                err_freq = mw_freq - self.settings['microwaves']['attocube_feedback']['freq_target'] # Error in ESR frequency (Hz)
                err_field = err_freq / 2.8e6*1e-4 # Error in B field (T)
                err_displacement = err_field / self.settings['microwaves']['attocube_feedback']['gradient'] * 1e9 # Error in displacement (nm)
                err_voltage = err_displacement / self.settings['microwaves']['attocube_feedback']['attocube_calibration'] # Error in Attocube voltage

                # Change sign of voltage adjustment based on whether it's the lower or upper dip
                # e.g. for the lower dip, a frequency that is below target means that the field is too large, requiring
                # us to move the Attocube down, vice versa for the upper dip
                if mw_freq < 2.87e9:
                    err_voltage = err_voltage * -1
                    print('Err voltage: '+ str(err_voltage))

                    print(np.all(np.abs(np.array(self.data['counts_diff'])[-8:]) < 0.03))
                    #print(np.abs(self.attocube_voltage_initial - self.attocube_voltage) < 8)
                    print(np.average(np.array(self.data['counts_diff'])[-10:]) < .006)


                if np.abs(err_voltage) < .5 and np.all(np.abs(np.array(self.data['counts_diff'])[-6:]) < 0.03) and np.average(np.array(self.data['counts_diff'])[-6:]) < .006:
                    self.instruments['piezo_controller']['instance'].axis = 'z'
                    self.attocube_voltage = self.instruments['piezo_controller']['instance'].read_probes('voltage')
                    self.attocube_voltage = float(self.attocube_voltage - err_voltage)
                    if np.abs(self.attocube_voltage_initial - self.attocube_voltage) < 8 and np.abs(err_freq) > 2e6:  # Restrict total range of Attocube movement
                        self.instruments['piezo_controller']['instance'].voltage = self.attocube_voltage
                        self.data['attocube_voltage'].append(self.attocube_voltage)
                        print('New voltage: ' + str(self.attocube_voltage))
                        print('Initial voltage: ' + str(self.attocube_voltage_initial))

            if self.settings['total_int_time'] > 0:
                self.progress = sample_index/max_samples
            else:
                self.progress = 50.
            self.updateProgress.emit(int(self.progress))

            time.sleep(float(sample_num) / sample_rate)
            sample_index = sample_index + 1
            if self.settings['total_int_time'] > 0. and sample_index >= max_samples: # if the maximum integration time is hit
                self._abort = True # tell the script to abort

            def restart_loop():
                self.task_a = self.instruments['daq']['instance'].setup_counter(self.settings['counter_channel_a'],
                                                                           sample_num, continuous_acquisition=True)
                self.task_b = self.instruments['daq']['instance'].setup_counter(self.settings['counter_channel_b'],
                                                                           sample_num, continuous_acquisition=True,
                                                                           existing_clock_channel=existing_clock_channel)
                self.instruments['daq']['instance'].run(self.task_a)
                time.sleep(self.settings['integration_time'])
                self.last_value_a = 0
                self.last_value_b = 0

            # Script crashes after about 500 loops if the sampling rate is set to 1 kS/s
            if self.loop_iteration > 200000 / sample_rate:
                print(self.loop_iteration)
                self.loop_iteration = 0
                self.instruments['daq']['instance'].stop(self.task_a)
                self.instruments['daq']['instance'].stop(self.task_b)
                restart_loop()

            if self.settings['nv_tracking']['on/off'] and \
                    (self.new_val_a < self.threshold_counts/2 or
                     self.new_val_b < self.threshold_counts/2):
                print(self.loop_iteration)
                self.loop_iteration = 0
                self.instruments['daq']['instance'].stop(self.task_a)
                self.instruments['daq']['instance'].stop(self.task_b)
                self.scripts['find_nv'].run()
                self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']
                restart_loop()

        # clean up APD tasks
        self.instruments['daq']['instance'].stop(self.task_a)
        self.instruments['daq']['instance'].stop(self.task_b)

        self.data['counts_a'] = list(self.data['counts_a'])
        self.data['counts_b'] = list(self.data['counts_b'])
        self.data['counts_diff'] = list(self.data['counts_diff'])
        self.data['mw_freq'] = list(self.data['mw_freq'])
        self.data['attocube_voltage'] = list(self.data['attocube_voltage'])

        self.instruments['microwave_generator']['instance'].update({'enable_modulation': False})
        self.settings['microwaves']['tracking']['freq_start'] = self.data['mw_freq'][-1]
        print('Next MW start freq: ' + str(self.settings['microwaves']['tracking']['freq_start']))


    def plot(self, figure_list):
        super(Daq_Read_Counter_2_Channel, self).plot([figure_list[1]])

    def _plot(self, axes_list, data = None):
        # COMMENT_ME

        if data is None:
            data = self.data

        if len(data['counts_diff']) > 0:
            array_to_plot = np.array([np.delete(data['counts_diff'], 0), np.delete(data['mw_freq'], 0)])
            #if self.settings['trim_plot'] != -1 and len(data['counts_diff']) > self.settings['trim_plot']:
            #        array_to_plot = array_to_plot[:, -self.settings['trim_plot']:]
            self.lns1 = axes_list[0].plot(array_to_plot[0], color='black', label='Normalized count diff')
            axes_list[0].axhline(y=0.0, color='grey', linestyle='-', alpha = .2)
            self.ax2 = axes_list[0].twinx()
            self.lns2 = self.ax2.plot(array_to_plot[1], color='red', label='MW freq')
            self.ax2.set_ylabel('[Hz]')

            lns = self.lns1 + self.lns2
            labs = [l.get_label() for l in lns]
            self.ax2.legend(lns, labs, loc=0)

    def _update_plot(self, axes_list, data=None):
        if data is None:
            data = self.data

        if data:
            array_to_plot = np.array([np.delete(data['counts_diff'], 0), np.delete(data['mw_freq'], 0)])

            if self.settings['trim_plot'] != -1 and len(data['counts_diff']) > self.settings['trim_plot']:
                    array_to_plot = array_to_plot[:, -self.settings['trim_plot']:]
            self.lns1[0].set_ydata(array_to_plot[0])
            self.lns1[0].set_xdata(range(0, len(array_to_plot[0])))
            self.lns2[0].set_ydata(array_to_plot[1])
            self.lns2[0].set_xdata(range(0, len(array_to_plot[1])))

            axes_list[0].relim()
            axes_list[0].autoscale_view()
            self.ax2.relim()
            self.ax2.autoscale_view()



"""def _update_plot(self, axes_list, data = None):
        if data is None:
            data = self.data

        if data:
            array_to_plot = np.array([np.delete(data['counts_diff'], 0), np.delete(data['mw_freq'], 0)])
            #array_to_plot = np.delete(data['counts_diff'], 0)

            # Decimate data to speed up plotting once data set becomes too large


            if self.settings['decimate_plot']:
                if len(data['counts_diff']) > self.max_pts:
                    array_to_plot = array_to_plot[:, ::int(len(data['counts_a']) / self.max_pts)]
                    #array_to_plot = array_to_plot[::int(len(data['counts_a']) / self.max_pts)]

            #update_counts_vs_pos(axes_list[0], array_to_plot, np.linspace(0, len(array_to_plot), len(array_to_plot)))
            update_counts(axes_list[0], array_to_plot)"""


class Daq_Read_Counter_NI6259(Daq_Read_Counter):
    _INSTRUMENTS = {'daq': NI6259}

class Daq_Read_Counter_NI6229(Daq_Read_Counter):
    _INSTRUMENTS = {'daq': NI6229}

if __name__ == '__main__':
    script = {}
    instr = {}
    script, failed, instr = Script.load_and_append({'Daq_Read_Counter': 'Daq_Read_Counter'}, script, instr)

    print(script)
    print(failed)
    print(instr)