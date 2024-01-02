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
import time, datetime, itertools
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_generic import PulsedExperimentGeneric
from copy import deepcopy
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, MicrowaveGenerator, Commander
from b26_toolkit.plotting.plots_1d import plot_pulses, update_pulse_plot, update_1d_simple, plot_1d_simple_freq
from pylabcontrol.core import Parameter
from b26_toolkit.data_processing.fit_functions import fit_pulsed_odmr, pulsed_odmr_double, pulsed_odmr_single
import time as t


class ParamSweepGeneric(PulsedExperimentGeneric):
    """
    Pulsed version of ESR. This script applies a microwave pulse at fixed power and duration while varying its frequency.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('param_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('param_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('param_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('param_switching_time', .01, float, 'time wait after changing parameter (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'commander': Commander}

    def define_sweep_parameters(self):
        """
        Define name of sweep parameters. Since this script is generic, it is coded with variables like 'param_start'.
        This function redirects the code to look for the corresponding settings in _DEFAULT_SETTINGS
        :return:
        """
        self.sweep_params = {'param_start': self.settings['param_start'],
                             'param_stop': self.settings['param_stop'],
                             'param_points': self.settings['param_points'],
                             'param_switching_time': self.settings['param_switching_time']}

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

    def _configure_instruments_start_of_sweep(self, param_current):
        """
        Configure instruments before running a sweep. For example, change MW frequency
        :return: None
        """
        raise NotImplementedError

    def _configure_param_array(self):
        """
        Contruct the parameter array and store it in a variable called 'mw_frequencies'. Despite the naming, it's just a list of parameters to be swept;
        it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        :return:
        """

        raise NotImplementedError

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

        self.define_sweep_parameters()
        self.data['fits'] = None

        self._configure_instruments_start_of_script()

        self._configure_param_array()

        # Divides the total number of averages requested into a number of slices of MAX_AVERAGES_PER_SCAN and a remainder.
        # This is required because the pulseblaster won't accept more than ~4E6 loops (22 bits available to store loop
        # number) so need to break it up into smaller chunks (use 1E6 so initial results display faster)

        self.pulse_sequences, self.tau_list, self.measurement_gate_width = self.create_pulse_sequences()

        self.num_averages = self.settings['num_averages']
        (self.num_1E5_avg_pb_programs, remainder) = divmod(self.num_averages, self.settings['averaging_block_size'])

        self.laser_status_before_script = self.instruments['PB']['instance'].settings['laser']['status']

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

        if num_daq_reads > 0:
            self._initialize_data(num_daq_reads, in_data)
        else:
            self._initialize_data(1, in_data)

        breaker = 0
        for loop_index, average_loop in enumerate(list(range(int(self.num_1E5_avg_pb_programs)))):
            time_start = t.time()
            if breaker:
                break

            param_indices = list(range(len(self.params)))
            if self.settings['randomize']:
                np.random.shuffle(param_indices)

            self._track_nv(scenario='before_block')

            self._daq_action_before_param_sweep(num_daq_reads)
            self.first_time_program_PB = True
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
                self._configure_instruments_start_of_sweep(self.param_current)
                time.sleep(self.sweep_params['param_switching_time'])

                self.current_averages = (average_loop + 1) * self.settings['averaging_block_size']
                self.result_current = None  # run_sweep will give a new value/array for self.result_current

                # Function is called "run_sweep" but we're not sweeping anything. We just want to use the generic function from
                # pulsed_experiment_base_script. In a typical pulsed experiment script we ask the PB to sweep through different taus, running
                # one averaging block for each tau. Here we just ask the PB to run one averaging block for one "tau" (i.e. same pulse seq)
                self._run_sweep(self.pulse_sequences, self.settings['averaging_block_size'], num_daq_reads)

                if self.result_current is not None and len(np.shape(self.result_current)) == 1 and len(self.result_current) == num_daq_reads:
                    # If there's a DAQ read inside _run_sweep(), result_current will be a 1D array, holding the values for each
                    # DAQ read for the current MW freq/param value. e.g. if there are two readout windows, result_current will have len 2
                    self.esr_counts[i][average_loop] = self.result_current
                    self.data['count_data'][i] = np.average(self.esr_counts[i][0:average_loop+1], axis=0)

            if breaker:
                self._daq_action_abort()
                break
            else:
                self._daq_action_after_param_sweep(num_daq_reads)

            # If the DAQ is read outside of _run_sweep, it must be run inside _daq_action_after_freq_sweep
            # The DAQ data will then have values for ALL frequency values. It will be a 2D array, with 1st dim being num of freq points,
            # and 2nd dim being num of daq reads
            if len(np.shape(self.result_current)) == 2 and np.shape(self.result_current)[1] == num_daq_reads:
                for i, index_random in enumerate(param_indices):
                    self.esr_counts[index_random][average_loop] = self.result_current[i]
                    self.data['count_data'][index_random] = np.average(self.esr_counts[index_random][0:average_loop + 1], axis=0)

            if 'count_data' in self.data.keys() and 'params' in self.data.keys():
                try:
                    freq = self.data['params']
                    counts = self.data['count_data'][:, np.shape(self.data['count_data'])[1]-1]
                    fits = fit_pulsed_odmr(freq, counts,
                                           rabi_freq=1/(self.settings['tau_mw'] *1e-9 * 2), contrast_1=0.05, contrast_2=0.05)
                    self.data['fits'] = fits
                except:
                    self.data['fits'] = None

            # Save data after each averaging block so that we can analyze data externally on the go
            if self.settings['save']:
                self.save_data()

            time_elapsed = t.time() - time_start
            self.log("Completed average block %i of %i in %s" %
                     (average_loop + 1, int(self.num_1E5_avg_pb_programs), str(datetime.timedelta(seconds=time_elapsed))[:-7]))

            # This is needed to update the plotting. We put it here after performing the fitting
            self.updateProgress.emit(float(loop_index + 1) / int(self.num_1E5_avg_pb_programs) * 100)

        if remainder != 0 and not self._abort:
            # Not gonna bother implementing this. Don't see why this will be useful.
            self.log('Total num of averages not a multiple of averaging block size. The remaining %i averages will be skipped over.' %
                     int(remainder))

        self._configure_instruments_end_of_script()
        self.updateProgress.emit(float(loop_index + 1) / int(self.num_1E5_avg_pb_programs) * 100)

        if (len(self.data['counts'][0]) == 1) and not self._abort:
            self.data['counts'] = np.array([item for sublist in self.data['counts'] for item in sublist])

    def _daq_action_before_param_sweep(self, num_daq_reads):
        """
        DAQ action before sweeping through all the MW freq/param values. Usually we want to setup and start the DAQ here and
        tell the PB to start running the sequences
        :param num_daq_reads: number of DAQ reads, or the number of readout windows in the pulse sequence
        :return:
        """
        pass

    def _daq_action_after_param_sweep(self, num_daq_reads):
        """
        DAQ action after sweeping through all the MW freq/param values. Usually we want to read the DAQ here and
        stop the PB if it's still running
        :param num_daq_reads: number of DAQ reads, or the number of readout windows in the pulse sequence
        :return:
        """
        pass

    def _daq_action_abort(self):
        """
        DAQ action when aborting the loop
        :return:
        """
        pass

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

        if 'count_data' in data.keys():
            plot_1d_simple_freq(axes_list[0], data['params'], [data['count_data']])
            plot_pulses(axes_list[1], self.pulse_sequences[self.sequence_index])
        if 'fits' in data.keys() and data['fits'] is not None:
            if len(data['fits']) == 6:
                # rabi_freq, baseline, contrast_1, contrast_2, freq_1, freq_2
                title = 'Pulsed ESR f1 = {:0.6e} Hz, f2 = {:0.6e} Hz, wo = {:0.2e},\n contrast1 = {:.1%}, contrast2 = {:.1%}' \
                    .format(data['fits'][4], data['fits'][5], data['fits'][0],
                            data['fits'][2], data['fits'][3])
            elif len(data['fits']) == 4:
                # rabi_freq, baseline, contrast_1, freq_1
                title = 'Pulsed ESR f1 = {:0.6e} Hz, wo = {:0.2e},\n contrast1 = {:.1%}' \
                    .format(data['fits'][3], data['fits'][0],
                            data['fits'][2])
            axes_list[0].set_title(title)

            x_data = data['params']
            x_data_fine = np.linspace(np.min(x_data), np.max(x_data), len(x_data) * 8)

        if 'fits' in data.keys() and data['fits'] is not None:
            if len(data['fits']) == 6:
                plot_1d_simple_freq(axes_list[0], x_data_fine, [pulsed_odmr_double(x_data_fine, *data['fits'])], alpha=0.7, title=title)
            elif len(data['fits']) == 4:
                plot_1d_simple_freq(axes_list[0], x_data_fine, [pulsed_odmr_single(x_data_fine, *data['fits'])], alpha=0.7, title=title)

    def _update_plot(self, axes_list):
        """
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)
        """

        data = self.data
        counts = data['count_data']
        x_data = data['params']
        x_data_fine = np.linspace(np.min(x_data), np.max(x_data), len(x_data)*8)

        # If fit is found and fit has not been plotted, plot both data and fit
        fit_in_plot = len(axes_list[0].lines) == len(np.transpose(counts)) + 1
        update_1d_simple(axes_list[0], x_data, counts, fit_in_plot=fit_in_plot)
        if 'fits' in data.keys() and data['fits'] is not None:
            if len(data['fits']) == 6:
                #rabi_freq, baseline, contrast_1, contrast_2, freq_1, freq_2
                title = 'Pulsed ESR f1 = {:0.6e} Hz, f2 = {:0.6e} Hz, wo = {:0.2e},\n contrast1 = {:.1%}, contrast2 = {:.1%}' \
                    .format(data['fits'][4], data['fits'][5], data['fits'][0],
                            data['fits'][2], data['fits'][3])
            elif len(data['fits']) == 4:
                # rabi_freq, baseline, contrast_1, freq_1
                title = 'Pulsed ESR f1 = {:0.6e} Hz, wo = {:0.2e},\n contrast1 = {:.1%}' \
                    .format(data['fits'][3], data['fits'][0],
                            data['fits'][2])
            axes_list[0].set_title(title)
            if fit_in_plot:
                if len(data['fits']) == 6:
                    for index, counts in enumerate([pulsed_odmr_double(x_data_fine, *data['fits'])]):
                        axes_list[0].lines[-1 - index].set_ydata(counts)
                elif len(data['fits']) == 4:
                    for index, counts in enumerate([pulsed_odmr_single(x_data_fine, *data['fits'])]):
                        axes_list[0].lines[-1 - index].set_ydata(counts)

            else:
                if len(data['fits']) == 6:
                    plot_1d_simple_freq(axes_list[0], x_data_fine, [pulsed_odmr_double(x_data_fine, *data['fits'])], alpha=0.7)
                elif len(data['fits']) == 4:
                    plot_1d_simple_freq(axes_list[0], x_data_fine, [pulsed_odmr_single(x_data_fine, *data['fits'])], alpha=0.7)

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

        raise NotImplementedError

    def _initialize_data(self, num_daq_reads, in_data):
        signal = [0.0]
        norms = np.repeat([0.0], (num_daq_reads - 1))
        self.count_data = np.repeat([np.append(signal, norms)], len(self.pulse_sequences), axis=0)
        self.data = in_data
        self.data['tau'] = np.array(self.tau_list)
        self.data['counts'] = deepcopy(self.count_data)
        self.data['params'] = self.params
        self.esr_counts = np.zeros((len(self.params), int(self.num_1E5_avg_pb_programs), num_daq_reads))
        self.data['count_data'] = np.zeros((len(self.params), num_daq_reads))
        if self.settings['save_full']: # ER 20210331
            print('num avgs in initialize ', self.num_averages)
            self.data['full_contrast'] = np.zeros((int(self.num_averages/self.settings['averaging_block_size']), len(self.pulse_sequences)))


class ParamSweepFastGeneric(ParamSweepGeneric):
    """
    Based on PulsedESRFast. Here we start the DAQ and PB pulse sequence before each freq sweep, and the read the DAQ after the freq sweep
    This reduces the overhead from the typical way of starting and reading the DAQ for each freq value.
    """

    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('param_start', 2.87e9, float, 'start frequency of scan in Hz'),
        Parameter('param_stop', 1e8, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'center_range', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('param_points', 40, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 340, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 800, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 80, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 40, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 260, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('spacer', 0, int, 'delay (ns) added to the beginning of the pulse sequence'),
        Parameter('param_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
    ]

    def _run_sweep(self, pulse_sequences, num_loops_sweep, num_daq_reads, verbose=False):
        # We just take the first pulse sequence, because they should be all the same for PulsedESR-type scripts
        if self.first_time_program_PB:
            self.instruments['PB']['instance'].program_pb(self.pulse_sequences[0],
                                                          num_loops=self.settings['averaging_block_size'], num_loops_w_pause=self.sweep_params['param_points'])
            self.first_time_program_PB = False

        # We use a lower-level start function because start_pulse_seq() closes the connection the PB after starting the pulse seq
        # Here we program the PB such that it waits for a start command after measuring an averaging block for each frequency, so we need to keep the comm open
        self.instruments['PB']['instance'].pb.pb_start()

        # Allocate some time for the PB to run before moving on to the next MW freq/parameter
        # This time is the estimated time of one averaging block with 4.5% safety margin
        # With inadequate sleep time, MW freq will switch in the middle of an averaging block! (Saw shift in ESR freqs with bad wait time)
        time.sleep(self.instruments['PB']['instance'].estimated_runtime/1000*1.045)

    def _daq_action_before_param_sweep(self, num_daq_reads):
        if self.settings['daq_type'] == 'PCI':
            daq = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            daq = self.instruments['NI9402']['instance']

        # Tell the DAQ to expect enough DAQ reads for an entire frequency sweep
        if num_daq_reads != 0:
            self.task = daq.setup_gated_counter('ctr0', int(self.settings['averaging_block_size'] * self.sweep_params['param_points'] * num_daq_reads))
            daq.run(self.task)

    def _daq_action_after_param_sweep(self, num_daq_reads):
        # The way create_commands in the pulse_blaster instrument works, each iteration (if WAIT feature is used) will end with a WAIT command
        # At the end of the freq sweep we need to explicitly close the PB connection otherwise the PB will wait for further instructions forever
        if self.instruments['PB']['instance'].settings['PB_type'] == 'USB':
            print('stopping pulse seq: ')
            self.instruments['PB']['instance'].stop_pulse_seq()
        elif self.instruments['PB']['instance'].settings['PB_type'] == 'PCI':
            self.instruments['PB']['instance'].pb.pb_close()

        num_loops_sweep = self.settings['averaging_block_size']*num_daq_reads
        if self.settings['daq_type'] == 'PCI':
            daq = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            daq = self.instruments['NI9402']['instance']

        result = []

        if num_daq_reads == 0:
            result = np.zeros((self.sweep_params['param_points'], 0))

        if num_daq_reads != 0:
            result = np.zeros((self.sweep_params['param_points'], num_daq_reads))
            result_array, temp = daq.read(self.task)  # thread waits on DAQ getting the right number of gates
            daq.stop(self.task)
            result_array = np.array_split(result_array, self.sweep_params['param_points'])

            for param_index, param_block in enumerate(result_array):
                for i in range(num_daq_reads):
                    result[param_index][i] = sum(itertools.islice(param_block, i, None, num_daq_reads))

        self.result_current = np.array([self._normalize_to_kCounts(np.array(param_block), self.measurement_gate_width,
                                                                   num_loops_sweep) for param_block in result])

    def _daq_action_abort(self):
        if self.settings['daq_type'] == 'PCI':
            daq = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            daq = self.instruments['NI9402']['instance']
        daq.stop(self.task)

        if self.instruments['PB']['instance'].settings['PB_type'] == 'PCI':
            try:
                self.instruments['PB']['instance'].pb.pb_close()
                print("Sent close command to PB")
            except:
                print("Could not send close command to PB")