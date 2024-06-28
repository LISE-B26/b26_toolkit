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
import datetime
import time
from pylabcontrol.core import Script, Parameter
from b26_toolkit.scripts.hf2li_scripts import LockInDaqRead
from b26_toolkit.scripts.pulse_sequences.pi_pulse_train import PiPulseTrainBackaction
from b26_toolkit.instruments import MicrowaveGenerator, Hf2Li, B26PulseBlaster
from b26_toolkit.plotting.plots_1d import plot_1d_simple_freq, update_1d_simple

class BackactionGeneric(Script):
    """
    Runs a pulse sequence while varying some parameter, and measures the mechanical response using the lock-in DAQ
    """
    _DEFAULT_SETTINGS = [
        Parameter('param_start', 2.87e9, float, 'start value of parameter sweep'),
        Parameter('param_stop', 1e8, float, 'end value of parameter sweep'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('param_points', 50, int, 'number of parameter points in scan'),
        Parameter('param_switching_time', .01, float, 'time wait after changing parameter value (s)'),
        Parameter('randomize', True, bool, 'check to randomize runs of the pulse sequence'),
    ]

    _INSTRUMENTS = {'mw_gen': MicrowaveGenerator, 'hf2li': Hf2Li, 'PB': B26PulseBlaster}
    _SCRIPTS = {'backaction_pulse_seq': PiPulseTrainBackaction, 'lock_in_daq_read': LockInDaqRead}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path=data_path)
        self.params = []

        self.sweep_params = {}
        self.sweep_params['param_start'] = self.settings['param_start']
        self.sweep_params['param_stop'] = self.settings['param_stop']
        self.sweep_params['param_points'] = self.settings['param_points']
        self.sweep_params['param_switching_time'] = self.settings['param_switching_time']

    def _configure_instruments_start_of_script(self):
        """
        Configure instruments at the beginning of a script, e.g. enable IQ modulation
        :return: None
        """
        pass

    def _configure_instruments_end_of_script(self):
        """
        Configure instruments right before the script finishes, e.g. turn off function generators
        :return: None
        """
        pass

    def _configure_instruments_for_param(self, param_current):
        """
        Configure instruments before running a sweep. For example, change MW frequency
        :return: None
        """
        raise NotImplementedError

    def _configure_subscripts_for_param(self, param_current):
        """

        :param param_current:
        :return:
        """

    def _configure_param_array(self):
        """
        Contruct the frequency array and store it in a variable called 'self.params'. Despite the naming, it's just a list of parameters to be swept;
        it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        :return:
        """

        raise NotImplementedError


    def _initialize_data(self):
        """
        Initialize arrays that will contain data
        :return: None
        """
        self.data['params'] = self.params
        self.data_labels = []
        for demod in self.scripts['lock_in_daq_read'].settings['demods'].keys():
            demod_num = int(demod)
            if self.scripts['lock_in_daq_read'].settings['demods'][demod]['data']['x']:
                self.data_labels.append('demod%i_x' % demod_num)
            if self.scripts['lock_in_daq_read'].settings['demods'][demod]['data']['y']:
                self.data_labels.append('demod%i_y' % demod_num)
            if self.scripts['lock_in_daq_read'].settings['demods'][demod]['data']['freq']:
                self.data_labels.append('demod%i_freq' % demod_num)

        print('Found %i data labels: ' % len(self.data_labels))
        self.data['value'] = np.zeros((len(self.data_labels), len(self.params)))

    def _function(self, in_data=None):
        """
        The core function can be configured in two ways. In the first, it iterates over each MW freq/param and calls the _run_sweep() function
        from pulsed_experiment_base_script. In this _run_sweep() function, we tell the PB to run the same pulse sequence for one averaging blcok
        and then read out the counts (it's a "sweep" over a list of one pulse seq, as opposed to the usual sweep over many pulse seqs with different
        pulse durations (tau).

        In the second way, we program and start the PB with averaging block size * num of freq points sequences, and read the DAQ counts AFTER iterating
        over all the freq values. This eliminates a lot of the overhead from having to program/read the DAQ repeatedly for each freq point, at the cost
        of less frequent updates. To synchronize the pulse seqs with the MW freq sweep, we tell the PB to "WAIT" between each averaging block for each freq
        and tell it to resume after adjusting the MW freq

        The second way is significantly faster for a slow DAQ (e.g. USB DAQ). It shouldn't do much for a PCI DAQ which has much less read/write overhead
        :param in_data: I don't know what this is
        :return: Nothing
        """

        self._configure_param_array()
        self._configure_instruments_start_of_script()
        self._initialize_data()

        # self.num_averages = self.settings['num_averages']
        self.num_averages = 1

        # Keeps track of index of current pulse sequence for plotting
        self.sequence_index = 0

        breaker = False
        for loop_index in range(self.num_averages):
            time_start = time.time()
            if breaker:
                break

            param_indices = list(range(len(self.params)))
            if self.settings['randomize']:
                np.random.shuffle(param_indices)

            for i in param_indices:
                if self._abort:
                    print('Aborting!!')
                    if self.instruments['PB']['instance'].settings['PB_type'] == 'USB':
                        self.instruments['PB']['instance'].stop_pulse_seq()
                    elif self.instruments['PB']['instance'].settings['PB_type'] == 'PCI':
                        try:
                            if self.instruments['PB']['instance'].pb.pb_read_status() & 0b100 == 0b100:  # Check if PB is running
                                self.instruments['PB']['instance'].pb.pb_close()
                            print('Sent close command to PB')
                        except:
                            print('Unable to read PB status. Did not send close command.')

                    self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
                    self.log('Aborted pulseblaster script during loop')
                    breaker = True
                    break

                self.param_current = self.params[i]
                self._configure_instruments_for_param(self.param_current)
                self._configure_subscripts_for_param(self.param_current)
                time.sleep(self.sweep_params['param_switching_time'])

                self.scripts['backaction_pulse_seq'].run()
                self.scripts['lock_in_daq_read'].run()
                self.data['value'][:,i] = self.process_daq()

            #  Save data after each averaging block so that we can analyze data externally on the go
            if self.settings['save']:
                self.save_data()

            time_elapsed = time.time() - time_start
            self.log("Completed loop %i of %i in %s" %
                     (loop_index + 1, self.num_averages, str(datetime.timedelta(seconds=time_elapsed))[:-7]))

            # This is needed to update the plotting. We put it here after performing the fitting
            self.updateProgress.emit(float(loop_index + 1) / int(self.num_averages) * 100)

        self._configure_instruments_end_of_script()

    def process_daq(self):
        """
        Process the data collected by the daq after running self.scripts['lock_in_daq_read']
        :return: 2D array (M*N): M = no. of demod channels (e.g. demod 4, quadrature X), N = length of data for each demod channel
        """
        raise NotImplementedError


class BackactionMwFreqSweep(BackactionGeneric):
    """
    Sweeps the MW freq of the pi-pulse train and measures the mechanical response using the lock-in DAQ
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_sweep', [
            Parameter('freq_start', 2.87e9, float, 'start frequency of scan in Hz'),
            Parameter('freq_stop', 1e8, float, 'end frequency of scan in Hz'),
            Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                      'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
            Parameter('freq_points', 50, int, 'number of frequencies in scan in Hz'),
            Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)')]),
        Parameter('fft_integration', [
            Parameter('freq_start', 10, float, 'center frequency (Hz) of integration window over FFT obtained from lock-in DAQ'),
            Parameter('freq_stop', 1, float, 'frequency span (Hz) of integration window over FFT obtained from lock-in DAQ'),
            Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                      'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop')
        ]),
        Parameter('randomize', True, bool, 'check to randomize runs of the pulse sequence')
    ]

    _INSTRUMENTS = {'mw_gen': MicrowaveGenerator, 'hf2li': Hf2Li, 'PB': B26PulseBlaster}
    _SCRIPTS = {'backaction_pulse_seq': PiPulseTrainBackaction, 'lock_in_daq_read': LockInDaqRead}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path=data_path)
        self.params = []

        self.sweep_params = {'param_start': self.settings['mw_sweep']['freq_start'],
                             'param_stop': self.settings['mw_sweep']['freq_stop'],
                             'param_points': self.settings['mw_sweep']['freq_points'],
                             'param_switching_time': self.settings['mw_sweep']['mw_generator_switching_time']}

    def _configure_instruments_start_of_script(self):
        """
        Configure instruments at the beginning of a script, e.g. enable IQ modulation
        :return: None
        """
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})


    def _configure_instruments_end_of_script(self):
        """
        Configure instruments right before the script finishes, e.g. turn off function generators
        :return: None
        """

        # Normally this is just done for safety. Here it also reprograms the pulseblaster so that it stops responding to the HW trigger
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

    def _configure_instruments_for_param(self, param_current):
        """
        Configure instruments before running a sweep. For example, change MW frequency
        :return: None
        """
        # self.instruments['mw_gen']['instance'].update({'frequency': float(param_current)})
        self.scripts['backaction_pulse_seq'].settings['mw_pulses']['mw_frequency'] = float(param_current)

    def _configure_subscripts_for_param(self, param_current):
        """

        :param param_current:
        :return:
        """
        # Enable FFT setting so that we have frequency domain data to process later
        self.scripts['lock_in_daq_read'].settings['fft'] = True

    def _configure_param_array(self):
        """
        Contruct the frequency array and store it in a variable called 'self.params'. Despite the naming, it's just a list of parameters to be swept;
        it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        :return:
        """

        mw_settings = self.settings['mw_sweep']
        if mw_settings['range_type'] == 'start_stop':
            if mw_settings['freq_start'] > mw_settings['freq_stop']:
                self.log('Error: end freq must be larger than start freq when range_type is start_stop..', flag='error')
                self._abort = True

            if mw_settings['freq_start'] < 0 or mw_settings['freq_stop'] > 4.05E9:
                self.log('Error: start or stop frequency out of bounds.', flag='error')
                self._abort = True

            self.params = np.linspace(mw_settings['freq_start'], mw_settings['freq_stop'], mw_settings['freq_points'])
        elif mw_settings['range_type'] == 'center_range':
            if mw_settings['freq_start'] < 2 * mw_settings['freq_stop']:
                self.log('End freq(range) must be smaller than 2x start freq (center) when range_type is center_range.', flag='error')
                self._abort = True
            self.params = np.linspace(mw_settings['freq_start'] - mw_settings['freq_stop'] / 2,
                                      mw_settings['freq_start'] + mw_settings['freq_stop'] / 2, mw_settings['freq_points'])

    def process_daq(self):
        integration_settings = self.settings['fft_integration']
        if integration_settings['range_type'] == 'start_stop':
            if integration_settings['freq_start'] > integration_settings['freq_stop']:
                self.log('Error: end freq must be larger than start freq when range_type is start_stop.', flag='error')
                self._abort = True
            int_freq_min, int_freq_max = integration_settings['freq_start'], integration_settings['freq_stop']
        elif integration_settings['range_type'] == 'center_range':
            if integration_settings['freq_start'] < 2 * integration_settings['freq_stop']:
                self.log('End freq(range) must be smaller than 2x start freq (center) when range_type is center_range.', flag='error')
                self._abort = True
            int_freq_min, int_freq_max = (integration_settings['freq_start'] - integration_settings['freq_stop'] / 2,
                                          integration_settings['freq_start'] + integration_settings['freq_stop'] / 2)
        else:
            raise KeyError('Unrecognized key for FFT integration settings')

        fft_freq = self.scripts['lock_in_daq_read'].data['fft_freq']
        int_freq_min_i, int_freq_max_i = np.argmin(np.abs(fft_freq-int_freq_min)), np.argmin(np.abs(fft_freq-int_freq_max))

        fft_integrals = []

        if self.scripts['lock_in_daq_read'].settings['segment_num'] > 1:
            self.log('Multiple segments detected in script "lock_in_daq_read". Data will be averaged over segments.')
        for data_label in self.data_labels:
            data = np.array(self.scripts['lock_in_daq_read'].data[data_label])
            fft_integrals.append(np.sum(np.average(data,axis=0)[int_freq_min_i:int_freq_max_i+1]))
        print(fft_integrals)
        return fft_integrals

    def _plot(self, axes_list, data=None):
        """

        :param axes_list:
        :param data:
        :return:
        """
        if data is None:
            data = self.data
        print('data[value]')
        print(data['value'])
        if len(data['value']) > 0:  # Make sure some data was written before plotting
            data_plt = []
            plt_labels = []
            plot_1d_simple_freq(axes_list[0], data['params'], data['value'])
            # axes_list[0].legend(plt_labels)

        else:
            print("Data not yet collected! No data plotted.")

    def _update_plot(self, axes_list):
        """
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)
        """

        update_1d_simple(axes_list[0], self.data['params'], self.data['value'])

