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
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_base_script import PulsedExperimentBaseScript
from copy import deepcopy
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, Pulse, MicrowaveGenerator, MicrowaveGenerator2, RFGenerator
from b26_toolkit.plotting.plots_1d import plot_1d_simple_timetrace_ns, plot_pulses, update_pulse_plot, update_1d_simple, plot_1d_simple_freq, plot_pulsedesr
from pylabcontrol.core import Script, Parameter
from b26_toolkit.data_processing.esr_signal_processing import fit_esr, fit_double_lorentzian, double_lorentzian
import time as t
import datetime


#MAX_AVERAGES_PER_SCAN = 5e4

class PulsedESR(PulsedExperimentBaseScript):
    """
This script applies a microwave pulse at fixed power and durations for varying frequencies.
Uses double_init scheme.

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _configure_instruments_start_of_script(self):
        """
        Configure instruments at the beginning of a script. For example, enable IQ modulation
        :return: None
        """
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_power']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

    def _configure_instruments_start_of_sweep(self):
        """
        Configure instruments before running a sweep. For example, change MW frequency
        :return: None
        """
        self.instruments['mw_gen']['instance'].update({'frequency': float(self.mw_frequency_current)})

    def _function(self, in_data=None):
        self.data['fits'] = None

        self._configure_instruments_start_of_script()

        # contruct the frequency array
        if self.settings['range_type'] == 'start_stop':
            if self.settings['freq_start'] > self.settings['freq_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True

            if self.settings['freq_start'] < 0 or self.settings['freq_stop'] > 4.05E9:
                self.log('start or stop frequency out of bounds')
                self._abort = True

            self.mw_frequencies = np.linspace(self.settings['freq_start'], self.settings['freq_stop'], self.settings['freq_points'])

        elif self.settings['range_type'] == 'center_range':
            if self.settings['freq_start'] < 2 * self.settings['freq_stop']:
                self.log('end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True
            self.mw_frequencies = np.linspace(self.settings['freq_start'] - self.settings['freq_stop'] / 2,
                                              self.settings['freq_start'] + self.settings['freq_stop'] / 2, self.settings['freq_points'])


        # divides the total number of averages requested into a number of slices of MAX_AVERAGES_PER_SCAN and a remainder.
        # This is required because the pulseblaster won't accept more than ~4E6 loops (22 bits available to store loop
        # number) so need to break it up into smaller chunks (use 1E6 so initial results display faster)
        self.pulse_sequences, self.tau_list, self.measurement_gate_width = self.create_pulse_sequences()
        self.num_averages = self.settings['num_averages']
        (self.num_1E5_avg_pb_programs, remainder) = divmod(self.num_averages, self.settings['averaging_block_size'])

        # retrieve initial mw carrier frequency to protect against bad fits in NV ESR tracking
        last_mw = self.scripts['esr'].instruments['microwave_generator']['instance'].frequency
        pulse_ampl = self.scripts['esr'].instruments['microwave_generator']['instance'].amplitude

        self.laser_status_before_script = self.instruments['PB']['instance'].settings['laser']['status']

        # ER 20181214 retrieve modulation on or off for main experiment
        mod_flag = self.scripts['esr'].instruments['microwave_generator']['instance'].enable_modulation

        # Keeps track of index of current pulse sequence for plotting
        self.sequence_index = 0

        if in_data is None:
            in_data = {}

        # calculates the number of daq reads per loop requested in the pulse sequence by asking how many apd reads are
        # called for. if this is not calculated properly, daq will either end too early (number too low) or hang since it
        # never receives the rest of the counts (number too high)
        num_daq_reads = 0

        for pulse in self.pulse_sequences[0]:
            if pulse.channel_id == 'apd_readout':
                num_daq_reads += 1
        self._initialize_data(num_daq_reads, in_data)

        # run find_nv if tracking is on ER 5/30/2017
        if self.settings['Tracking']['on/off']:
            self.scripts['find_nv'].run()
            if self.scripts['find_nv'].data[
                'fluorescence'] == 0.0:  # if it doesn't find an NV, abort the experiment
                self.log('Could not find an NV in FindNV.')
                self._abort = True
                return  # exit function in case no NV is found

        #self.log("Averaging over %i blocks of %.1e" % (self.num_1E5_avg_pb_programs, self.settings['averaging_block_size']))

        breaker = 0
        for average_loop in range(int(self.num_1E5_avg_pb_programs)):
            time_start = t.time()
            if breaker:
                break
            #self.log("Running average block {0} of {1}".format(average_loop + 1, int(self.num_1E5_avg_pb_programs)))

            mw_freq_indices = list(range(len(self.mw_frequencies)))
            if self.settings['randomize']:
                np.random.shuffle(mw_freq_indices)

            for i in mw_freq_indices:
                if self._abort:
                    print('aborting!!')
                    # ER 20200828 stop the pulseblaster
                    if self.instruments['PB']['instance'].settings['PB_type'] == 'USB':
                        print('stopping pulse seq: abort!! ')
                        self.instruments['PB']['instance'].stop_pulse_seq()

                    self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

                    self.log('Aborted pulseblaster script during loop')
                    breaker = 1
                    break

                # MM 20190621
                # if self.settings['Autofocus_Tracking']['on/off'] and average_loop % self.settings['Autofocus_Tracking']['track_every_N'] == 0:
                #    self.scripts['autofocus'].run()
                # Retrack NV afterwards.
                #    self.scripts['find_nv'].run()
                #    self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']

                time.sleep(self.settings['mw_generator_switching_time'])
                self.mw_frequency_current = self.mw_frequencies[i]
                self._configure_instruments_start_of_sweep()

                # ER 20181028
                if self.settings['ESR_Tracking']['on/off'] and average_loop % self.settings['ESR_Tracking'][
                    'track_every_N'] == 0:
                    self.scripts['esr'].run()

                    # retrieve the new mw frequency: if there are two frequencies in the fit, pick the one closest to the old frequency
                    fit_params = self.scripts['esr'].data['fit_params']

                    # default update flag to false
                    update_mw = False

                    if fit_params is not None and len(fit_params) and fit_params[0] != -1:  # check if fit valid
                        if len(fit_params) == 4:
                            # single peak
                            if (fit_params[2] - last_mw) ** 2 < (self.settings['ESR_Tracking'][
                                                                     'allowed_delta_freq'] * 1e6) ** 2:  # check if new value is within range allowed
                                update_mw = True
                            new_mw = fit_params[2]
                        elif len(fit_params) == 6:
                            # double peak, don't update the frequency - the fit may be bad
                            update_mw = False

                    if update_mw:
                        # self.instruments['mw_gen'].update({'frequency': new_mw})
                        self.scripts['esr'].instruments['microwave_generator']['instance'].update(
                            {'frequency': float(new_mw)})
                        self.log('Updated mw carrier frequency to: {}'.format(new_mw))
                        self.scripts['esr'].instruments['microwave_generator']['instance'].update(
                            {'amplitude': float(pulse_ampl)})
                        self.scripts['esr'].instruments['microwave_generator']['instance'].update(
                            {'enable_modulation': bool(mod_flag)})

                        last_mw = new_mw
                    else:
                        # self.instruments['mw_gen'].update({'frequency': last_mw})
                        self.scripts['esr'].instruments['microwave_generator']['instance'].update(
                            {'frequency': float(last_mw)})
                        self.log(
                            'Not updating the mw carrier frequency. SRS carrier frequency kept at {0} Hz'.format(
                                last_mw))
                        self.scripts['esr'].instruments['microwave_generator']['instance'].update(
                            {'amplitude': float(pulse_ampl)})
                        self.scripts['esr'].instruments['microwave_generator']['instance'].update(
                            {'enable_modulation': bool(mod_flag)})

                self.current_averages = (average_loop + 1) * self.settings['averaging_block_size']
                #    print('tau sequences running: ', self.tau_list)

                self._run_sweep(self.pulse_sequences, self.settings['averaging_block_size'], num_daq_reads)
                self.esr_counts[i][average_loop] = self.result_current
                self.data['esr_counts'][i] = np.average(self.esr_counts[i][0:average_loop+1], axis=0)

                # save data on the fly so that we can start to analyze it while the experiment is running!
                if self.settings['save']:
                    self.save_data()

            time_elapsed = t.time() - time_start
            self.log("Completed average block %i of %i in %s" %
                     (average_loop + 1, int(self.num_1E5_avg_pb_programs), str(datetime.timedelta(seconds=time_elapsed))[:-7]))

            if 'esr_counts' in self.data.keys() and 'mw_frequencies' in self.data.keys():
                try:
                    freq = self.data['mw_frequencies']
                    counts = self.data['esr_counts'][:, 0]
                    freq_spacing = freq[1]-freq[0]
                    est_contrast = (min(counts)-max(counts))
                    fits = fit_double_lorentzian(freq, counts,
                                                 starting_params = [max(counts),1e6,est_contrast,est_contrast,
                                                                    min(freq)+(max(freq)-min(freq))*0.2,max(freq)-(max(freq)-min(freq))*0.2])
                    self.data['fits'] = fits
                    print(fits)
                except:
                    self.data['fits'] = None
                    self.log('ESR fit failed')

        if remainder != 0 and not self._abort:
            self.current_averages = self.num_averages
            self._run_sweep(self.pulse_sequences, remainder, num_daq_reads)

        if (len(self.data['counts'][0]) == 1) and not self._abort:
            print("(len(self.data['counts'][0]) == 1) and not self._abort")
            self.data['counts'] = np.array([item for sublist in self.data['counts'] for item in sublist])

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

        #        if self.scripts['find_nv'].is_running: #self.scripts['find_nv'].data['maximum_point']:
        #            # self.scripts['find_nv'].plot(axes_list)
        #             self.scripts['find_nv']._plot(axes_list)
        #             self._plot_refresh = True
        #        else:
        if data is None:
            data = self.data

        if 'esr_counts' in data.keys():
            plot_1d_simple_freq(axes_list[0], data['mw_frequencies'], [data['esr_counts']])
            plot_pulses(axes_list[1], self.pulse_sequences[self.sequence_index])
        if 'fits' in data.keys() and data['fits'] is not None:
            title = 'Pulsed ESR f1 = {:0.6e} Hz, f2 = {:0.6e} Hz, wo = {:0.2e},\n contrast1 = {:.1%}, contrast2 = {:.1%}' \
                .format(data['fits'][4], data['fits'][5], data['fits'][1],
                        np.abs(data['fits'][2]) / data['fits'][0], np.abs(data['fits'][3]) / data['fits'][0])
            axes_list[0].set_title(title)

            x_data = data['mw_frequencies']
            x_data_fine = np.linspace(np.min(x_data), np.max(x_data), len(x_data) * 2)
            plot_1d_simple_freq(axes_list[0], x_data_fine, [double_lorentzian(x_data_fine, *data['fits'])], alpha=0.5, title=title)

    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''
        #        if self.scripts['find_nv'].is_running:
        #            self.scripts['find_nv']._update_plot(axes_list)
        #        else:

        data = self.data
        counts = data['esr_counts']
        x_data = data['mw_frequencies']
        x_data_fine = np.linspace(np.min(x_data), np.max(x_data), len(x_data)*2)

        # If fit is found and fit has not been plotted, plot both data and fit
        fit_in_plot = len(axes_list[0].lines) == len(np.transpose(counts)) + 1
        update_1d_simple(axes_list[0], x_data, counts, fit_in_plot=fit_in_plot)
        if 'fits' in data.keys() and data['fits'] is not None:
            title = 'Pulsed ESR f1 = {:0.6e} Hz, f2 = {:0.6e} Hz, wo = {:0.2e},\n contrast1 = {:.1%}, contrast2 = {:.1%}'\
                .format(data['fits'][4], data['fits'][5], data['fits'][1],
                        np.abs(data['fits'][2]) / data['fits'][0], np.abs(data['fits'][3]) / data['fits'][0])
            axes_list[0].set_title(title)
            if fit_in_plot:
                for index, counts in enumerate([double_lorentzian(x_data_fine, *data['fits'])]):
                    axes_list[0].lines[-1 - index].set_ydata(counts)

            else:
                plot_1d_simple_freq(axes_list[0], x_data_fine, [double_lorentzian(x_data_fine, *data['fits'])], alpha=0.5)

        update_pulse_plot(axes_list[1], self.pulse_sequences[self.sequence_index])


    def _create_pulse_sequences(self):

        '''

        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''

        tau = self.settings['tau_mw']
        pulse_sequences = []
        tau_list = [tau]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', laser_off_time + tau + 2 * 40, nv_reset_time),
                 Pulse('apd_readout',
                       laser_off_time + tau + 2 * 40 + delay_readout, meas_time)]
               #  Pulse('apd_readout', laser_off_time + tau + 2 * 40 + delay_readout, meas_time),
               #  ]
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence += [
                    Pulse(microwave_channel, laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time, tau)]

            pulse_sequence += [
                Pulse('laser',
                      laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2 * 40 + delay_mw_readout,
                      nv_reset_time),
                Pulse('apd_readout',
                      laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2 * 40 + delay_mw_readout + delay_readout,
                      meas_time)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, self.settings['read_out']['meas_time']

    def _initialize_data(self, num_daq_reads, in_data):
        signal = [0.0]
        norms = np.repeat([0.0], (num_daq_reads - 1))
        self.count_data = np.repeat([np.append(signal, norms)], len(self.pulse_sequences), axis=0)
        self.data = in_data
        self.data['tau'] = np.array(self.tau_list)
        self.data['counts'] = deepcopy(self.count_data)
        self.data['mw_frequencies'] = self.mw_frequencies
        self.esr_counts = np.zeros((len(self.mw_frequencies), int(self.num_1E5_avg_pb_programs), num_daq_reads))
        self.data['esr_counts'] = np.zeros((len(self.mw_frequencies),
                                            num_daq_reads))
        if self.settings['save_full']: # ER 20210331
            print('num avgs in initialize ', self.num_averages)
            self.data['full_contrast'] = np.zeros((int(self.num_averages/self.settings['averaging_block_size']), len(self.pulse_sequences)))


class PulsedESRFast(PulsedESR):
    def _create_pulse_sequences(self):

        '''

        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''

        tau = self.settings['tau_mw']
        pulse_sequences = []
        tau_list = [tau]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = []
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence += [
                    Pulse(microwave_channel, laser_off_time, tau)]

            pulse_sequence += [
                Pulse('laser',
                      laser_off_time + tau + 2 * 40 + delay_mw_readout,
                      nv_reset_time),
                Pulse('apd_readout',
                      laser_off_time + tau + 2 * 40 + delay_mw_readout + delay_readout,
                      meas_time),
                Pulse('apd_readout',
                      laser_off_time + tau + 2 * 40 + delay_mw_readout + nv_reset_time - meas_time,
                      meas_time)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, self.settings['read_out']['meas_time']


class PulsedESRFast_MWGen2(PulsedESR):
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator2}


class PulsedNMR(PulsedESR):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.82e9, float,
                      'frequency of hyperfine transition corresponding to desired nuclear spin polarization'),
            Parameter('tau_mw', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('tau_rf', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('cooldown_time', 5000, float,
                  'cooldown time between running each sequence, to minimize heating effects'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator,
                    'rf_gen': RFGenerator}

    def _create_pulse_sequences(self):
        '''

        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_tau = self.settings['mw_pulses']['tau_mw']
        microwave_channel = 'microwave_i'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        tau = self.settings['mw_pulses']['tau_mw']
        rf_tau = self.settings['rf_pulses']['tau_rf']
        cooldown_time = self.settings['cooldown_time']
        pulse_sequences = []
        tau_list = [tau]

        for tau in tau_list:
            pulse_sequence = [Pulse('laser', cooldown_time + laser_off_time, nv_reset_time),
                              #Pulse('apd_readout', tau + laser_off_time + delay_readout, meas_time),
                              Pulse(microwave_channel, cooldown_time + laser_off_time + nv_reset_time + laser_off_time, mw_tau)
                              ]
            # if tau is 0 there is actually no mw pulse
            if rf_tau > 0:
                pulse_sequence.append(
                    Pulse('rf_switch', cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle, rf_tau))
                pulse_sequence.append(
                    Pulse(rf_channel, cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle - 60, rf_tau + 120))

            pulse_sequence.append(Pulse(microwave_channel,
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + rf_tau + rf_settle,
                                        mw_tau))
            pulse_sequence.append(Pulse('laser',
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + rf_tau + rf_settle + mw_tau + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + rf_tau + rf_settle + mw_tau + delay_mw_readout + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + rf_tau + rf_settle + mw_tau + delay_mw_readout + nv_reset_time - meas_time,
                                        meas_time))
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:

            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, self.settings['read_out']['meas_time']

    def _function(self, in_data=None):
        self.instruments['PB']['instance'].update({'rf_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['rf_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['rf_gen']['instance'].update({'enable_modulation': True})
        self.instruments['rf_gen']['instance'].update({'amplitude_rf': self.settings['rf_pulses']['rf_power']})
        #self.instruments['rf_gen']['instance'].update({'frequency': self.settings['rf_frequency']})
        #self.instruments['rf_gen']['instance'].update({'enable_output': True})

        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

        assert self.settings['freq_start'] < self.settings['freq_stop']
        self.rf_frequencies = np.linspace(self.settings['freq_start'], self.settings['freq_stop'],
                                          self.settings['freq_points'])

        # divides the total number of averages requested into a number of slices of MAX_AVERAGES_PER_SCAN and a remainder.
        # This is required because the pulseblaster won't accept more than ~4E6 loops (22 bits available to store loop
        # number) so need to break it up into smaller chunks (use 1E6 so initial results display faster)
        self.pulse_sequences, self.tau_list, self.measurement_gate_width = self.create_pulse_sequences()
        self.num_averages = self.settings['num_averages']
        (self.num_1E5_avg_pb_programs, remainder) = divmod(self.num_averages, self.settings['averaging_block_size'])

        # retrieve initial mw carrier frequency to protect against bad fits in NV ESR tracking
        last_mw = self.scripts['esr'].instruments['microwave_generator']['instance'].frequency
        pulse_ampl = self.scripts['esr'].instruments['microwave_generator']['instance'].amplitude

        # ER 20181214 retrieve modulation on or off for main experiment
        mod_flag = self.scripts['esr'].instruments['microwave_generator']['instance'].enable_modulation

        # Keeps track of index of current pulse sequence for plotting
        self.sequence_index = 0

        if in_data is None:
            in_data = {}

        # calculates the number of daq reads per loop requested in the pulse sequence by asking how many apd reads are
        # called for. if this is not calculated properly, daq will either end too early (number too low) or hang since it
        # never receives the rest of the counts (number too high)
        num_daq_reads = 0

        for pulse in self.pulse_sequences[0]:
            if pulse.channel_id == 'apd_readout':
                num_daq_reads += 1
        self._initialize_data(num_daq_reads, in_data)

        # run find_nv if tracking is on ER 5/30/2017
        if self.settings['Tracking']['on/off']:
            self.scripts['find_nv'].run()
            if self.scripts['find_nv'].data[
                'fluorescence'] == 0.0:  # if it doesn't find an NV, abort the experiment
                self.log('Could not find an NV in FindNV.')
                self._abort = True
                return  # exit function in case no NV is found

        self.log("Averaging over %i blocks of %.1e" % (self.num_1E5_avg_pb_programs, self.settings['averaging_block_size']))

        breaker = 0
        for average_loop in range(int(self.num_1E5_avg_pb_programs)):
            self.log("Running average block {0} of {1}".format(average_loop + 1, int(self.num_1E5_avg_pb_programs)))

            rf_freq_indices = list(range(len(self.rf_frequencies)))
            if self.settings['randomize']:
                np.random.shuffle(rf_freq_indices)
            print(rf_freq_indices)
            for i in rf_freq_indices:

                if breaker:
                    break
                time.sleep(self.settings['mw_generator_switching_time'])
                self.instruments['rf_gen']['instance'].update({'frequency': float(self.rf_frequencies[i])})
                print('Switching to %.2e Hz' % float(self.rf_frequencies[i]))
                print(i)

                if self._abort:
                    print('aborting!!')
                    # ER 20200828 stop the pulseblaster
                    if self.instruments['PB']['instance'].settings['PB_type'] == 'USB':
                        print('stopping pulse seq: abort!! ')
                        self.instruments['PB']['instance'].stop_pulse_seq()

                    self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

                    self.log('Aborted pulseblaster script during loop')
                    breaker = 1
                    break


                # MM 20190621
                # if self.settings['Autofocus_Tracking']['on/off'] and average_loop % self.settings['Autofocus_Tracking']['track_every_N'] == 0:
                #    self.scripts['autofocus'].run()
                # Retrack NV afterwards.
                #    self.scripts['find_nv'].run()
                #    self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']

                # ER 20181028
                if self.settings['ESR_Tracking']['on/off'] and average_loop % self.settings['ESR_Tracking'][
                    'track_every_N'] == 0:
                    self.scripts['esr'].run()

                    # retrieve the new mw frequency: if there are two frequencies in the fit, pick the one closest to the old frequency
                    fit_params = self.scripts['esr'].data['fit_params']

                    # default update flag to false
                    update_mw = False

                    if fit_params is not None and len(fit_params) and fit_params[0] != -1:  # check if fit valid
                        if len(fit_params) == 4:
                            # single peak
                            if (fit_params[2] - last_mw) ** 2 < (self.settings['ESR_Tracking'][
                                                                     'allowed_delta_freq'] * 1e6) ** 2:  # check if new value is within range allowed
                                update_mw = True
                            new_mw = fit_params[2]
                        elif len(fit_params) == 6:
                            # double peak, don't update the frequency - the fit may be bad
                            update_mw = False

                    if update_mw:
                        # self.instruments['mw_gen'].update({'frequency': new_mw})
                        self.scripts['esr'].instruments['microwave_generator']['instance'].update(
                            {'frequency': float(new_mw)})
                        self.log('Updated mw carrier frequency to: {}'.format(new_mw))
                        self.scripts['esr'].instruments['microwave_generator']['instance'].update(
                            {'amplitude': float(pulse_ampl)})
                        self.scripts['esr'].instruments['microwave_generator']['instance'].update(
                            {'enable_modulation': bool(mod_flag)})

                        last_mw = new_mw
                    else:
                        # self.instruments['mw_gen'].update({'frequency': last_mw})
                        self.scripts['esr'].instruments['microwave_generator']['instance'].update(
                            {'frequency': float(last_mw)})
                        self.log(
                            'Not updating the mw carrier frequency. SRS carrier frequency kept at {0} Hz'.format(
                                last_mw))
                        self.scripts['esr'].instruments['microwave_generator']['instance'].update(
                            {'amplitude': float(pulse_ampl)})
                        self.scripts['esr'].instruments['microwave_generator']['instance'].update(
                            {'enable_modulation': bool(mod_flag)})

                self.current_averages = (average_loop + 1) * self.settings['averaging_block_size']
                #    print('tau sequences running: ', self.tau_list)

                self._run_sweep(self.pulse_sequences, self.settings['averaging_block_size'], num_daq_reads)
                self.esr_counts[i][average_loop] = self.result_current
                #print(self.esr_counts[i][0:average_loop + 1])
                self.data['esr_counts'][i] = np.average(self.esr_counts[i][0:average_loop + 1], axis=0)

            # save data on the fly so that we can start to analyze it while the experiment is running!
            if self.settings['save']:
                #     self.save_b26()
                self.save_data()
            #     self.save_log()
            #     self.save_image_to_disk()

        if remainder != 0 and not self._abort:
            self.current_averages = self.num_averages
            self._run_sweep(self.pulse_sequences, remainder, num_daq_reads)

        if (len(self.data['counts'][0]) == 1) and not self._abort:
            print("(len(self.data['counts'][0]) == 1) and not self._abort")
            self.data['counts'] = np.array([item for sublist in self.data['counts'] for item in sublist])




    def _initialize_data(self, num_daq_reads, in_data):
        signal = [0.0]
        norms = np.repeat([0.0], (num_daq_reads - 1))
        self.count_data = np.repeat([np.append(signal, norms)], len(self.pulse_sequences), axis=0)
        self.data = in_data
        self.data['tau'] = np.array(self.tau_list)
        self.data['counts'] = deepcopy(self.count_data)
        self.data['rf_frequencies'] = self.rf_frequencies
        self.esr_counts = np.zeros((len(self.rf_frequencies), int(self.num_1E5_avg_pb_programs), num_daq_reads))
        self.data['esr_counts'] = np.zeros((len(self.rf_frequencies),
                                            num_daq_reads))
        if self.settings['save_full']: # ER 20210331
            print('num avgs in initialize ', self.num_averages)
            self.data['full_contrast'] = np.zeros((int(self.num_averages/self.settings['averaging_block_size']), len(self.pulse_sequences)))

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

        #        if self.scripts['find_nv'].is_running: #self.scripts['find_nv'].data['maximum_point']:
        #            # self.scripts['find_nv'].plot(axes_list)
        #             self.scripts['find_nv']._plot(axes_list)
        #             self._plot_refresh = True
        #        else:
        if data is None:
            data = self.data

        if 'esr_counts' in data.keys():
            # The following does not work for pulsedelays; you need to comment out the 'if' for it to work.
            # if counts != []:
            #     plot_1d_simple_timetrace_ns(axes_list[0], data['tau'], [data['cousants'])
          #  print('plotting with tau values: ', data['tau'])
            plot_1d_simple_timetrace_ns(axes_list[0], data['rf_frequencies'], [data['esr_counts']])
            plot_pulses(axes_list[1], self.pulse_sequences[self.sequence_index])

    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''
        #        if self.scripts['find_nv'].is_running:
        #            self.scripts['find_nv']._update_plot(axes_list)
        #        else:
        counts = self.data['esr_counts']
        x_data = self.data['rf_frequencies']
        axis1 = axes_list[0]
        if not counts == []:
            update_1d_simple(axis1, x_data, [counts])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])


class PulsedESRPolarized(PulsedESR):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power_1', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency_1', 2.82e9, float,
                      'frequency of hyperfine transition corresponding to desired nuclear spin polarization'),
            Parameter('tau_mw_1', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_power_2', -45.0, float, 'microwave power in dBm'),
            Parameter('tau_mw_2', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('tau_rf', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('cooldown_time', 5000, float, 'cooldown time between running each sequence, to minimize heating effects'),
        Parameter('polarization_iterations', 2, int, 'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int, 'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            #Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int, 'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'mw_gen_2': MicrowaveGenerator2,
                    'rf_gen': RFGenerator}

    def _create_pulse_sequences(self):
        '''

        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''

        #nv_reset_time = self.settings['read_out']['nv_reset_time']
        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_tau_1 = self.settings['mw_pulses']['tau_mw_1']
        mw_tau_2 = self.settings['mw_pulses']['tau_mw_2']
        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        tau = self.settings['mw_pulses']['tau_mw_1']
        rf_tau = self.settings['rf_pulses']['tau_rf']
        cooldown_time = self.settings['cooldown_time']
        pulse_sequences = []
        tau_list = [tau]

        for tau in tau_list:

            pulse_sequence = [Pulse('laser', laser_off_time + cooldown_time, nv_reset_time_long)]
            current_time = laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_tau_1))
                if rf_tau > 0:
                    pulse_sequence.append(
                        Pulse('rf_switch', current_time + mw_tau_1 + rf_settle, rf_tau))
                    pulse_sequence.append(
                        Pulse(rf_channel, current_time + mw_tau_1 + rf_settle - 60, rf_tau + 120))
                    current_time = current_time + mw_tau_1 + rf_settle + rf_tau + rf_settle
                #if self.settings['initialization_laser']['chopped']:
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait
                #else:
                #    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time))
                #    current_time = current_time + nv_reset_time



            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + laser_off_time,
                                        mw_tau_2))
            pulse_sequence.append(Pulse('laser',
                                        current_time + laser_off_time + mw_tau_2 + delay_mw_readout,
                                        nv_reset_time_long))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + laser_off_time + mw_tau_2 + delay_mw_readout + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + laser_off_time + mw_tau_2 + delay_mw_readout + nv_reset_time_long - meas_time,
                                        meas_time))
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:

            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, self.settings['read_out']['meas_time']

    def _configure_instruments_start_of_script(self):
        self.instruments['PB']['instance'].update({'rf_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        #self.instruments['rf_gen']['instance'].update({'modulation_type': 'IQ'})
        #self.instruments['rf_gen']['instance'].update({'enable_modulation': True})
        #self.instruments['rf_gen']['instance'].update({'amplitude_rf': self.settings['rf_pulses']['rf_power']})

        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power_1']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency_1']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power_2']})
        self.instruments['mw_gen_2']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency_1']})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})

    def _configure_instruments_start_of_sweep(self):
        self.instruments['mw_gen_2']['instance'].update({'frequency': float(self.mw_frequency_current)})

    def _function_GARBAGE(self, in_data=None):

        self.instruments['PB']['instance'].update({'rf_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['rf_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['rf_gen']['instance'].update({'enable_modulation': True})
        self.instruments['rf_gen']['instance'].update({'amplitude_rf': self.settings['rf_pulses']['rf_power']})

        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power_1']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency_1']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power_2']})
        self.instruments['mw_gen_2']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency_1']})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})

        assert self.settings['freq_start'] < self.settings['freq_stop']
        self.mw_frequencies = np.linspace(self.settings['freq_start'], self.settings['freq_stop'],
                                          self.settings['freq_points'])

        # divides the total number of averages requested into a number of slices of MAX_AVERAGES_PER_SCAN and a remainder.
        # This is required because the pulseblaster won't accept more than ~4E6 loops (22 bits available to store loop
        # number) so need to break it up into smaller chunks (use 1E6 so initial results display faster)
        self.pulse_sequences, self.tau_list, self.measurement_gate_width = self.create_pulse_sequences()
        self.num_averages = self.settings['num_averages']
        (self.num_1E5_avg_pb_programs, remainder) = divmod(self.num_averages, self.settings['averaging_block_size'])

        # retrieve initial mw carrier frequency to protect against bad fits in NV ESR tracking
        last_mw = self.scripts['esr'].instruments['microwave_generator']['instance'].frequency
        pulse_ampl = self.scripts['esr'].instruments['microwave_generator']['instance'].amplitude

        # ER 20181214 retrieve modulation on or off for main experiment
        mod_flag = self.scripts['esr'].instruments['microwave_generator']['instance'].enable_modulation

        # Keeps track of index of current pulse sequence for plotting
        self.sequence_index = 0

        if in_data is None:
            in_data = {}

        # calculates the number of daq reads per loop requested in the pulse sequence by asking how many apd reads are
        # called for. if this is not calculated properly, daq will either end too early (number too low) or hang since it
        # never receives the rest of the counts (number too high)
        num_daq_reads = 0

        for pulse in self.pulse_sequences[0]:
            if pulse.channel_id == 'apd_readout':
                num_daq_reads += 1
        self._initialize_data(num_daq_reads, in_data)

        # run find_nv if tracking is on ER 5/30/2017
        if self.settings['Tracking']['on/off']:
            self.scripts['find_nv'].run()
            if self.scripts['find_nv'].data[
                'fluorescence'] == 0.0:  # if it doesn't find an NV, abort the experiment
                self.log('Could not find an NV in FindNV.')
                self._abort = True
                return  # exit function in case no NV is found

        self.log("Averaging over %i blocks of %.1e" % (self.num_1E5_avg_pb_programs, self.settings['averaging_block_size']))
        for average_loop in range(int(self.num_1E5_avg_pb_programs)):
            self.log("Running average block {0} of {1}".format(average_loop + 1, int(self.num_1E5_avg_pb_programs)))

            mw_freq_indices = list(range(len(self.mw_frequencies)))
            if self.settings['randomize']:
                np.random.shuffle(mw_freq_indices)
            for i in mw_freq_indices:
                time.sleep(self.settings['mw_generator_switching_time'])
                self.instruments['mw_gen_2']['instance'].update({'frequency': float(self.mw_frequencies[i])})
                print('Switching to %.2e Hz' % float(self.mw_frequencies[i]))
                print(i)

                if self._abort:
                    print('aborting!!')
                    # ER 20200828 stop the pulseblaster
                    if self.instruments['PB']['instance'].settings['PB_type'] == 'USB':
                        print('stopping pulse seq: abort!! ')
                        self.instruments['PB']['instance'].stop_pulse_seq()

                    self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

                    self.log('Aborted pulseblaster script during loop')
                    break


                self.current_averages = (average_loop + 1) * self.settings['averaging_block_size']
                #    print('tau sequences running: ', self.tau_list)

                self._run_sweep(self.pulse_sequences, self.settings['averaging_block_size'], num_daq_reads)
                self.esr_counts[i][average_loop] = self.result_current
                #print(self.esr_counts[i][0:average_loop + 1])
                self.data['esr_counts'][i] = np.average(self.esr_counts[i][0:average_loop + 1], axis=0)

        if remainder != 0 and not self._abort:
            self.current_averages = self.num_averages
            self._run_sweep(self.pulse_sequences, remainder, num_daq_reads)

        if (len(self.data['counts'][0]) == 1) and not self._abort:
            print("(len(self.data['counts'][0]) == 1) and not self._abort")
            self.data['counts'] = np.array([item for sublist in self.data['counts'] for item in sublist])

        # save data on the fly so that we can start to analyze it while the experiment is running!
        if self.settings['save']:
            #     self.save_b26()
            self.save_data()
        #     self.save_log()
        #     self.save_image_to_disk()


class PulsedEsrRfHeating(PulsedESRPolarized):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power_1', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency_1', 2.82e9, float,
                      'frequency of hyperfine transition corresponding to desired nuclear spin polarization'),
            Parameter('tau_mw_1', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_power_2', -45.0, float, 'microwave power in dBm'),
            Parameter('tau_mw_2', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('tau_rf', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('cooldown_time', 5000, float, 'cooldown time between running each sequence, to minimize heating effects'),
        Parameter('polarization_iterations', 2, int, 'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int, 'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            #Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int, 'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'mw_gen_2': MicrowaveGenerator2,
                    'rf_gen': RFGenerator}

    def _create_pulse_sequences(self):
        '''

        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''

        #nv_reset_time = self.settings['read_out']['nv_reset_time']
        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_tau_1 = self.settings['mw_pulses']['tau_mw_1']
        mw_tau_2 = self.settings['mw_pulses']['tau_mw_2']
        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        tau = self.settings['mw_pulses']['tau_mw_1']
        rf_tau = self.settings['rf_pulses']['tau_rf']
        cooldown_time = self.settings['cooldown_time']
        pulse_sequences = []
        tau_list = [tau]

        for tau in tau_list:



            pulse_sequence = []
            current_time = laser_off_time + cooldown_time
            if rf_tau > 0:
                pulse_sequence.append(
                    Pulse(rf_channel, current_time, rf_tau))
                pulse_sequence.append(
                    Pulse('rf_switch', current_time, rf_tau))

            current_time += rf_tau + rf_settle
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time,
                                        mw_tau_2))
            pulse_sequence.append(Pulse('laser',
                                        current_time + mw_tau_2 + delay_mw_readout,
                                        nv_reset_time_long))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + mw_tau_2 + delay_mw_readout + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + mw_tau_2 + delay_mw_readout + nv_reset_time_long - meas_time,
                                        meas_time))


            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:

            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, self.settings['read_out']['meas_time']


class NuclearPolarizationLaserDuration(PulsedExperimentBaseScript):
    # Please update this with structure from PulsedESRPolarized, FF
    _DEFAULT_SETTINGS = [
        Parameter('mw_power_1', -45.0, float, 'microwave power in dBm'),
        Parameter('mw_frequency_1', 2.82e9, float, 'frequency of hyperfine transition that we want to depopulate'),
        Parameter('mw_frequency_2', 2.82e9, float, 'frequency of hyperfine transition for which we want to measure the population'),
        Parameter('tau_mw_1', 80, float, 'the time duration of the microwaves in ns'),
        Parameter('mw_power_2', -45.0, float, 'microwave power in dBm'),
        Parameter('tau_mw_2', 80, float, 'the time duration of the microwaves in ns'),
        Parameter('rf_power', -10, float, 'RF power in dBm'),
        Parameter('tau_rf', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
        Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns'),
        Parameter('cooldown_time', 5000, float, 'cooldown time between running each sequence'),
        Parameter('polarization_iterations', 2, int, 'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('initialization_laser', [
            Parameter('chopped', False, bool, 'whether to use chopped laser pulses for initialization, overrides nv_reset_time if True'),
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('min_n', 1, int, 'min num of initialization laser pulses'),
            Parameter('max_n', 200, int, 'max num of initialization laser pulses'),
            Parameter('n_step', 1, int, 'step increment of number of initialization laser pulses')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int, 'time with laser on to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ])
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'mw_gen_2': MicrowaveGenerator2,
                    'rf_gen': RFGenerator}

    def _function(self, in_data=None):
        self.instruments['PB']['instance'].update({'rf_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['rf_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['rf_gen']['instance'].update({'enable_modulation': True})
        self.instruments['rf_gen']['instance'].update({'amplitude_rf': self.settings['rf_power']})
        #self.instruments['rf_gen']['instance'].update({'frequency': self.settings['rf_frequency']})
        #self.instruments['rf_gen']['instance'].update({'enable_output': True})

        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_power_1']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_frequency_1']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_power_2']})
        self.instruments['mw_gen_2']['instance'].update({'frequency': self.settings['mw_frequency_2']})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})

        super(NuclearPolarizationLaserDuration, self)._function()

        if 'counts' in self.data.keys() and 'tau' in self.data.keys():
            counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
            tau = self.data['tau']


            self.data['fits'] = None
            self.log('rabi fit failed')

    def _create_pulse_sequences(self):
        '''
        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement
        '''
        pulse_sequences = []

        if self.settings['initialization_laser']['chopped']:
            tau_list = range(int(self.settings['initialization_laser']['min_n']),
                             int(self.settings['initialization_laser']['max_n']),
                             int(self.settings['initialization_laser']['n_step']))
        else:
            print('Not implemented for CW')

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_tau_1 = self.settings['tau_mw_1']
        mw_tau_2 = self.settings['tau_mw_2']
        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_settle']

        #tau = self.settings['tau_mw_1']
        rf_tau = self.settings['tau_rf']
        pulse_sequences = []

        tau_list = range(int(self.settings['initialization_laser']['min_n']),
                         int(self.settings['initialization_laser']['max_n']),
                         int(self.settings['initialization_laser']['n_step']))

        for tau in tau_list:
            # cooldown_time = int(rf_tau/2. / 2.) * 2 # time between each pulse sequence, as a multiple of tau, such that we have a constant duty cycle of rf pulses to avoid overheating
            cooldown_time = self.settings['cooldown_time']
            pulse_sequence = [Pulse('laser', laser_off_time + cooldown_time, nv_reset_time_long)]
            current_time = laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_tau_1))
                if rf_tau > 0:
                    pulse_sequence.append(
                        Pulse('rf_switch', current_time + mw_tau_1 + rf_settle, rf_tau))
                    pulse_sequence.append(
                        Pulse(rf_channel, current_time + mw_tau_1 + rf_settle - 60, rf_tau + 120))
                    current_time = current_time + mw_tau_1 + rf_settle + rf_tau + rf_settle
                if self.settings['initialization_laser']['chopped']:
                    for i in range(int(tau)):
                        pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                        current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait
                else:
                    raise NotImplementedError
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time))
                    current_time = current_time + nv_reset_time

            pulse_sequence.append(Pulse(microwave_channel,
                                        current_time + laser_off_time,
                                        mw_tau_1))
            pulse_sequence.append(Pulse('laser',
                                        current_time + laser_off_time + mw_tau_1 + delay_mw_readout,
                                        nv_reset_time_long))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + laser_off_time + mw_tau_1 + delay_mw_readout + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + laser_off_time + mw_tau_1 + delay_mw_readout + nv_reset_time_long - meas_time,
                                        meas_time))
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time