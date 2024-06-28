"""
This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

pylabcontrol is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pylabcontrol is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pylabcontrol.  If not, see <http://www.gnu.org/licenses/>.
"""

import random
import datetime
import time as t
import itertools
from copy import deepcopy
import numpy as np
import time
from pylabcontrol.core import Script, Parameter
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, Pulse, MicrowaveGenerator
from b26_toolkit.plotting.plots_1d import plot_1d_simple_timetrace, plot_pulses, update_pulse_plot, update_1d_simple


class PulsedExperimentGeneric(Script):
    """
    This class is a base class that should be inherited by all classes that utilize the pulseblaster for experiments. The
    _function part of this class takes care of high-level interaction with the pulseblaster for experiment control and optionally
    the daq for reading counter input (usually from the APD). It also provides all of the functionality needed to run a
    standard Script such as plotting.
    To use this class, the inheriting class need only overwrite _create_pulse_sequences to create the proper pulse sequence
    for a given experiment
    """
    _DEFAULT_SETTINGS = [
        Parameter('averaging_block_size', 50000, int, 'number of averages in each averaging block, '
                                                      'too small of a block size probably introduces a larger relative overhead from reading DAQ/updating plots'),
        Parameter('randomize', True, bool, 'check to randomize runs of the pulse sequence'),
        Parameter('mw_switch', [
            Parameter('add', True, bool,  'check to add mw switch to every i/q pulse and to use switch to carve out pulses. Note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 60, int, 'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'mw_switch', ['mw_switch', 'mw_iq'],'determines if mw pulses are carved out by mw-switch or by i and q channels of mw source '),
            Parameter('no_iq_overlap', True, bool,'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ]),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'daq to be used for pulse sequence'),
        Parameter('save_full', False, bool, 'save every average'),
        Parameter('track_nv', [
            Parameter('below_threshold', False, bool, 'run FindNv whenever counts fall below a given threshold'),
            Parameter('threshold', 0.85, float, 'run FindNv whenever counts fall below this threshold'),
            Parameter('init_fluor', 20., float, 'initial fluorescence of the NV to compare to, in kcps'),
            Parameter('before_block', False, bool, 'run FindNv before each averaging block')]),
        Parameter('track_focus', [
            Parameter('threshold', 0.85, float, 'run FindNv whenever counts fall below this threshold'),
            Parameter('init_fluor', 20., float, 'initial fluorescence of the NV to compare to, in kcps'),
            Parameter('before_block', False, bool, 'run FindNv before each averaging block')])
    ]
    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    # Leave _SCRIPTS = {}. To enable tracking, use pulsed_experiment_tracking instead!
    _SCRIPTS = {}

    def __init__(self, instruments, name=None, settings=None, log_function=None, data_path=None):
        '''
        Initializes GalvoScan script for use in gui

        Args:
            instruments: list of instrument objects
            name: name to give to instantiated script object
            settings: dictionary of new settings to pass in to override defaults
            log_function: log function passed from the gui to direct log calls to the gui log
            data_path: path to save data

        '''

        self._DEFAULT_SETTINGS += PulsedExperimentGeneric._DEFAULT_SETTINGS
        Script.__init__(self, name, settings=settings, instruments=instruments, log_function=log_function,
                        data_path = data_path)

        # Index of the daq read bin that will be used to check if FindNv is needed
        # e.g. Typically, if this is 0, that means we will check the reference fluor (from the first daq read)
        self.ref_index = 0

    def _calc_progress(self, index):
        """
        Calculates progress of inner loop (in _run_sweep)
        :param index:
        :return: Script progress in percentage
        """
        progress_inner = index / len(self.pulse_sequences)

        # progress of outer loop (in _function)
        if self.current_averages >= self.settings['averaging_block_size']:
            progress = (self.current_averages + (1.0 - progress_inner) * self.settings['averaging_block_size']) / self.num_averages
        else:
            # if self.current_averages < self.settings['averaging_block_size'], there is only a single run and
            # therefore the progress of the inner loop equals the outer loop
            progress = progress_inner

        self.progress = 100.0 * progress
        return int(round(self.progress))

    def _configure_instruments_start_of_script(self):
        """
        Configure instruments at the beginning of a script, e.g. enable IQ modulation
        :return: None
        """
        pass

    def _configure_instruments_for_tau(self, tau):
        """
        Configure instruments before taking measurements for a new tau value.
        :return: None
        """
        pass

    def _configure_instruments_end_of_script(self):
        """
        Configure instruments right before the script finishes, e.g. turn off function generators
        :return: None
        """
        pass

    def _function(self, in_data=None):
        """
        This is the core loop in which the desired experiment specified by the inheriting script's pulse sequence
        is performed.

        in_data: input data dictionary, caution 'tau' and 'counts' will be overwritten here!

        Poststate: self.data contains two key/value pairs, 'tau' and 'counts'
            'tau': a list of the times tau that are scanned over for the relative experiment (ex wait times between pulses)
             'counts': the counts received from the experiment. This is a len('tau') list, with each element being a list
             of length 1, 2, or 3 corresponding to the sum over all trials for a single tau time. In this sublist, the
             first value is the signal, the second (optional) value is counts in the |0> state (the maximum counts for
             normalization), and the third (optional) value is the counts in the |1> state (the minimum counts for
             normalization)

        """

        self._configure_instruments_start_of_script()

        # make sure the microwave_switch is turned off so that we don't burn any steel cables. ER 20181017
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # Remember if laser was on prior to this script
        self.laser_status_before_script = self.instruments['PB']['instance'].settings['laser']['status']

        # Keeps track of index of current pulse sequence for plotting
        self.sequence_index = 0

        # self.is_valid and create pulses
        self.pulse_sequences, self.tau_list, self.measurement_gate_width = self.create_pulse_sequences()
        self.num_averages = self.settings['num_averages']

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

        # divides the total number of averages requested into a number of slices of MAX_AVERAGES_PER_SCAN and a remainder.
        # This is required because the pulseblaster won't accept more than ~4E6 loops (22 bits available to store loop
        # number) so need to break it up into smaller chunks (use 1E6 so initial results display faster)
        (num_1E5_avg_pb_programs, remainder) = divmod(self.num_averages, self.settings['averaging_block_size']) # Name is not accurate anymore, block size is no longer fixed to 1e5, FF

        self.log("Averaging over %i blocks of %.1e" % (num_1E5_avg_pb_programs, self.settings['averaging_block_size']))

        for average_loop in range(int(num_1E5_avg_pb_programs)):
            time_start = t.time()
            if self._abort:
                print('Aborting!!')
                # ER 20200828 stop the pulseblaster
                if self.instruments['PB']['instance'].settings['PB_type'] == 'USB':
                    print('Stopping pulse seq: abort!! ')
                    self.instruments['PB']['instance'].stop_pulse_seq()

                self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

                self.log('Aborted pulseblaster script during loop')
                break

            self._track_nv(scenario='before_block')
            self.current_averages = (average_loop + 1) * self.settings['averaging_block_size']
            self._run_sweep(self.pulse_sequences, self.settings['averaging_block_size'], num_daq_reads)

            # save data on the fly so that we can start to analyze it while the experiment is running!
            if self.settings['save']:
                self.save_data()

            time_elapsed = t.time() - time_start
            self.log("Completed average block %i of %i in %s" %
                     (average_loop + 1, int(num_1E5_avg_pb_programs), str(datetime.timedelta(seconds=time_elapsed))[:-7]))
            if 'loop_delay' in self.settings:
                time.sleep(self.settings['loop_delay'])

        if remainder != 0 and not self._abort:
            self.current_averages = self.num_averages
            self._run_sweep(self.pulse_sequences, remainder, num_daq_reads)

        if (len(self.data['counts'][0]) == 1) and not self._abort:
            self.data['counts'] = np.array([item for sublist in self.data['counts'] for item in sublist])

        if self.instruments['PB']['instance'].settings['constant_microwave_heat_load']['enable']:
            self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
            self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
            self.instruments['mw_gen']['instance'].update({'amplitude': self.instruments['PB']['instance'].settings['constant_microwave_heat_load']['mw_power']})
            self.instruments['mw_gen']['instance'].update({'frequency': self.instruments['PB']['instance'].settings['constant_microwave_heat_load']['mw_frequency']})
            self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
            self.instruments['mw_gen']['instance'].update({'enable_output': True})
            self.instruments['PB']['instance'].mw_duty_cycle_loop()

    def _initialize_data(self, num_daq_reads, in_data):
        signal = [0.0]
        norms = np.repeat([0.0], (num_daq_reads - 1))
        self.count_data = np.repeat([np.append(signal, norms)], len(self.pulse_sequences), axis=0)
        self.data = in_data
        self.data['tau'] = np.array(self.tau_list)
        self.data['counts'] = deepcopy(self.count_data)
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

        if data is None:
            data = self.data

        if 'counts' in data.keys():
            # The following does not work for pulsedelays; you need to comment out the 'if' for it to work.
            # if counts != []:
            #     plot_1d_simple_timetrace_ns(axes_list[0], data['tau'], [data['cousants'])
            plot_1d_simple_timetrace(axes_list[0], data['tau'], [data['counts']])
            plot_pulses(axes_list[1], self.pulse_sequences[self.sequence_index])

    def _update_plot(self, axes_list):
        """
        Updates plots specified in _plot above

        Args:
            :param axes_list: list of axes to write plots to (uses first 2)
            :return: None
        """

        counts = self.data['counts']
        x_data = self.data['tau']
        axis1 = axes_list[0]
        if not counts is []:
            update_1d_simple(axis1, x_data, [counts])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])

    def _run_sweep(self, pulse_sequences, num_loops_sweep, num_daq_reads, verbose=False):
        """
        Each pulse sequence specified in pulse_sequences is run num_loops_sweep consecutive times.

        Args:
            :param pulse_sequences: a list of pulse sequences to run, each corresponding to a different value of tau. Each
                             sequence is a list of Pulse objects specifying a given pulse sequence
            :param num_loops_sweep: number of times to repeat each sequence before moving on to the next one
            :param num_daq_reads: number of times the daq must read for each sequence (generally 1, 2, or 3)

        Poststate: self.data['counts'] is updated with the acquired data

        """

        rand_indexes = list(range(len(pulse_sequences)))

        # if 'track_nv' in self.settings and 'find_nv' in self.scripts:
        #     threshold = self.settings['track_nv']['threshold']
        #     init_fluor = self.settings['track_nv']['init_fluor']

        if self.settings['randomize']:
            random.shuffle(rand_indexes)
        if verbose:
            print(('_run_sweep number of pulse sequences', len(pulse_sequences)))

        breaker = False
        for index, sequence in enumerate(pulse_sequences):
            if breaker or self._abort:
                break
            if verbose:
                print(('_run_sweep index', index, len(pulse_sequences)))
            rand_index = rand_indexes[index]

            if self._abort:
                print('aborting in run sweep!!')

                if self.instruments['PB']['instance'].settings['PB_type'] == 'USB':
                    print('stopping pulse seq in run sweep: abort!! ')
                    self.instruments['PB']['instance'].stop_pulse_seq()

                self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
                breaker = True
                break

            self._configure_instruments_for_tau(self.tau_list[rand_index])

            timer_start = time.time()
            result = self._run_single_sequence(pulse_sequences[rand_index], num_loops_sweep, num_daq_reads)  # keep entire array

            #print('Time for running single seq: %.2e' % float(time.time() - timer_start))

            # self.result_current is added here for pulsed_esr. In a usual pulse sequence (e.g. Rabi), there are 2
            # loops: tau and averaging blocks. In pulsed_esr, there are 3 loops: MW frequency, tau (with 1 value), and
            # averaging blocks. Because of this different structure I added self.result_current so that I can do the
            # averaging in _function()
            self.result_current = self._normalize_to_kCounts(np.array(result), self.measurement_gate_width, num_loops_sweep)

            # for tracking ER 20210331
            counts_to_check = self._normalize_to_kCounts(np.array(result), self.measurement_gate_width, num_loops_sweep)
            counts_temp = counts_to_check[self.ref_index]  # ref_index is set as 0 by default when script is initialized

            if 'commander' in self.instruments:
                if self.instruments['commander']['instance'].settings['find_nv']:
                    self._track_nv(scenario='force')
                    self.instruments['commander']['instance'].settings['find_nv'] = False

                if self.instruments['commander']['instance'].settings['autofocus'] and 'autofocus' in self.scripts:
                    print('Autofocus not implemented')

            # counts_unsatisfactory = (1 + (1 - threshold)) * init_fluor < counts_temp or threshold * init_fluor > counts_temp
            abort_script = self._track_nv(scenario='check_threshold', counts_temp=counts_temp)  # If FindNv fails a few times to locate NV, abort script
            if abort_script:
                self._abort = True
                breaker = True

            # Added another breaker here so that data does not get saved if the sweep is aborted.
            # Not the most elegant solution but this works
            if breaker:
                break

            # ER 20210331
            if self.settings['save_full']:
                self.data['full_contrast'][int(self.current_averages/self.settings['averaging_block_size'])-1][rand_index] = \
                    (counts_to_check[1]-counts_to_check[0])/np.mean(counts_to_check)

            self.count_data[rand_index] = self.count_data[rand_index] + result
            self.data['counts'][rand_index] = self._normalize_to_kCounts(self.count_data[rand_index], self.measurement_gate_width,
                                                                         self.current_averages)

            self.sequence_index = rand_index
            if index+1 == len(rand_indexes):
                self.sequence_index_next = rand_indexes[index] # We can't show the upcoming sequence because we haven't shuffled the next sweep yet
            else:
                self.sequence_index_next = rand_indexes[index+1]

            self.updateProgress.emit(self._calc_progress(index))

    def _run_single_sequence(self, pulse_sequence, num_loops, num_daq_reads):
        '''
        Runs a single pulse sequence, num_loops consecutive times
        Args:
            pulse_sequence: a list of Pulse objects specifying a pulse sequence
            num_loops: number of times to repeat the pulse sequence
            num_daq_reads: number of times sequence requires that the

        Returns: a list containing, 1, 2, or 3 values depending on the pulse sequence
        counts, the second is the number of

        '''

        if self.settings['daq_type'] == 'PCI':
            daq = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            daq = self.instruments['NI9402']['instance']

        timer_start = time.time()
        try:
            self.instruments['PB']['instance'].program_pb(pulse_sequence, num_loops=num_loops)
        except AssertionError:
            self.log('Error programming PB, aborting script')
            self._abort = True
            return np.zeros(num_daq_reads)

        #print('Time after program PB: %.2e' % (time.time() - timer_start))

        # TODO(AK): figure out if timeout is actually needed
        timeout = 2 * self.instruments['PB']['instance'].estimated_runtime

        if num_daq_reads != 0:
            task = daq.setup_gated_counter('ctr0', int(num_loops * num_daq_reads))
            #print('Time after setup counter: %.2e' % (time.time() - timer_start))
            daq.run(task)
            #print('Time after run daq: %.2e' % (time.time() - timer_start))

        self.instruments['PB']['instance'].start_pulse_seq()
        #print('Time after start seq: %.2e' % (time.time() - timer_start))

        result = []

        if num_daq_reads == 0:
            result = [0]

        if num_daq_reads != 0:
            result_array, temp = daq.read(task)   # thread waits on DAQ getting the right number of gates
            #print('Time after read daq: %.2e' % (time.time() - timer_start))
            #t4 = t.perf_counter()
            for i in range(num_daq_reads):
                result.append(sum(itertools.islice(result_array, i, None, num_daq_reads)))
            #print('Time after slice: %.2e' % (time.time() - timer_start))

        # clean up APD tasks
        if num_daq_reads != 0:
            daq.stop(task)
            #print('Time after stop task: %.2e' % (time.time() - timer_start))

        # if num_daq_reads == 0:
        #     # If we're reading DAQ samples, we know that an expected number of sequences has run when we've collected enough samples
        #     # If there are no DAQ samples to wait for, we need to manually wait for the estimated pulse seq duration
        #     time.sleep(self.instruments['PB']['instance'].estimated_runtime*1e-3)

        if self.instruments['PB']['instance'].settings['PB_type'] == 'USB':
            print('stopping pulse seq: ')
            self.instruments['PB']['instance'].stop_pulse_seq()

        return result

    def _sum_measurements(self, daq_results, num_daq_reads):

        summed_results = []
        for i in range(num_daq_reads):
            summed_results.append(sum(itertools.islice(daq_results, i, None, num_daq_reads)))
        return summed_results

    def _create_pulse_sequences(self):
        '''
        A function to create the pulse sequence. This must be overwritten in scripts inheriting from this script
        The pulse sequences are pulse-blaster friendly opposed to the settings which are human-readable!!

        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences


        '''
        raise NotImplementedError

    def _normalize(self, signal, baseline_max=0, baseline_min=0):
        """
        Normalizes the signal values given a maximum value (counts in |0>) and optionally minimum value (counts in |1>,
        set to 0 if none provided). If no maximum given, returns the bare signal.
        Args:
            signal:
            baseline_max:
            baseline_min:

        Returns:

        """
        if baseline_max == 0:
            return signal
        else:
            return (signal - baseline_min) / (baseline_max - baseline_min)

    def _normalize_to_kCounts(self, signal, gate_width=1, num_averages=1):
        """
        Converts the signal from counts/gate_width as returned by the PB to kcounts/s
        Args:
            signal: data to normalize
            gate_width: length of each bin
            num_averages: number of times each bin is read

        Returns:

        """
        return signal * 1E6 / (gate_width * num_averages)  # 1E6 is to convert from ns to ms

    def is_valid(self, pulse_sequences=None, verbose=True):
        """
        Checks if a pulse blaster script is valid

        Args:
            create_pulse: is true creates the pulses first before validation
            verbose: print more stuff if true

        Returns:

        """
        self.verbose = verbose
        # make sure that we have the pulse sequences that correspond to the current settings
        if pulse_sequences is None:
            pulse_sequences, __, __ = self.create_pulse_sequences()

        for pulse_sequence in pulse_sequences:
            if self._is_bad_pulse_sequence(pulse_sequence):
                print('bad pulse sequence is: ', pulse_sequence)
                return False

        # Calculate total length of pulse sequence
        # We only do this if there's a single sequence, because that's usually when we care about the duty cycle
        # e.g. For PulsedESR, we just run the same sequence repeatedly for different freqs
        # e.g. For HahnEcho, we run through a list of sequences with variable tau, and we don't usually care about the duty cycle

        self.sequence_durations = []
        self.laser_duties = []
        self.mw_duties = []

        for pulse_sequence in pulse_sequences:

            laser_duration = 0
            mw_duration = 0
            pulse_end_times = []

            for pulse in pulse_sequence:
                if pulse.duration is not None:
                    pulse_end_times.append(pulse.start_time + pulse.duration)
                    if pulse.channel_id == 'laser':
                        laser_duration += pulse.duration
                    if pulse.channel_id == 'microwave_switch':
                        mw_duration += pulse.duration

            sequence_duration = np.max(pulse_end_times)
            self.sequence_durations.append(sequence_duration)
            self.laser_duties.append(laser_duration/sequence_duration)
            self.mw_duties.append(mw_duration/sequence_duration)
            if 'mw_power' in self.settings:
                power_dbm = self.settings['mw_power']
            elif 'mw_pulses' in self.settings and 'mw_power' in self.settings['mw_pulses']:
                power_dbm = self.settings['mw_pulses']['mw_power']
            elif 'mw_pulse' in self.settings and 'mw_power' in self.settings['mw_pulses']:
                power_dbm = self.settings['mw_pulse']['mw_power']
            else:
                power_dbm = -100
                print('No microwave power setting found in script!')

            # power_dbm = self.instruments['mw_gen']['instance'].settings['amplitude']
            self.mw_avg_power = np.array(self.mw_duties) * 0.001 * 10**(power_dbm/10)

        if verbose:
            if len(pulse_sequences) == 1:
                if self.laser_duties[0] == 0:
                    laser_duty_str = 'None'
                else:
                    laser_duty_str = str('%.2e' % self.laser_duties[0])
                if self.mw_duties[0] == 0:
                    mw_duty_str = 'None'
                    mw_avg_power_str = 'None'
                else:
                    mw_duty_str = str('%.2e' % self.mw_duties[0])
                    mw_avg_power_str = str('%.2e' % self.mw_avg_power[0])

                log_str = (self.sequence_durations[0]*1e-9, laser_duty_str, mw_duty_str, mw_avg_power_str)
                self.log('Seq duration: %.2e, laser duty: %s, MW duty: %s, MW avg power: %s' % log_str)
            else:
                log_str = (np.min(self.sequence_durations) * 1e-9, np.max(self.sequence_durations) * 1e-9,
                          np.min(self.laser_duties), np.max(self.laser_duties),
                          np.min(self.mw_duties), np.max(self.mw_duties),
                           np.min(self.mw_avg_power), np.max(self.mw_avg_power))
                self.log('Seq duration: %.2e to %.2e s, laser duty: %.2e to %.2e, MW duty: %.2e to %.2e, MW avg power: %.2e to %.2e' %
                         log_str)

        return True

    def _get_overlapping_pulses(self, pulse_sequence, verbose=False):
        """
        Finds and returns a list of ordered pairs of pulses from pulse_sequence that currently overlap, thus likely
        not being the pulse sequence that was intended. (It would be strange to hope for this)
        Args:
            pulse_sequence:
            verbose: if True, calls extra print() statements during execution, helpful for debugging.

        Returns:
            (list) a list of pairs of overlapping pulses in pulse_sequence

        """
        if 'mw_switch' in self.settings and self.settings['mw_switch']['no_iq_overlap']:
            overlapping_pulses = B26PulseBlaster.find_overlapping_pulses(pulse_sequence,
                                                                         combine_channels=['microwave_i',
                                                                                           'microwave_q'])
        else:
            overlapping_pulses = B26PulseBlaster.find_overlapping_pulses(pulse_sequence)

        if verbose and overlapping_pulses:
                print('Found overlapping pulses:', overlapping_pulses)

        return overlapping_pulses

    def _get_commands_that_are_too_short(self, pulse_sequence, verbose=False):
        """
        Finds if the current pulse sequence would compile to commands that have intervals of less than 15 ns between
        them, which it not supported by the pulse blaster as it is currently implemented.

        Args:
            pulse_sequence(list): list of pulses to check for commands that are too close
            verbose: if True, calls extra print() statements during execution, helpful for debugging.

        Returns:
            (list) A list of pulseblaster commands resulting from pulse_sequence that are too short

        """
        verbose=True
        pulse_blaster = self.instruments['PB']['instance']

        delayed_pulse_collection = pulse_blaster.create_physical_pulse_seq(pulse_sequence)

        pb_state_changes = pulse_blaster.generate_pb_sequence(delayed_pulse_collection)
        pb_commands = pulse_blaster.create_commands(pb_state_changes, self.settings['num_averages'])
        short_pulses = [command for command in pb_commands if command.duration < pulse_blaster.settings['min_pulse_dur']]

        if verbose and short_pulses:
            print('Found short pulses: ', short_pulses)

        return short_pulses

    def _has_pulses_with_bad_start_time(self, pulse_sequence, verbose=False):
        """
        Checks for pulses that have a start_time between 0 and 1 (exclusive), indicating that the user may have
        used the wrong units for the pulses in the pulse sequence.

        Args:
            pulse_sequence(list): list of pulses to check for bad start_times
            verbose(bool): if True, calls extra print() statements during execution, helpful for debugging.

        Returns:
            (bool) True if the pulse sequence contains a pulse with start time between 0 and 1, else False.
        """
        for pulse in pulse_sequence:
            if 0 < pulse.start_time < 1:
                if verbose:
                    print('Found pulse with start time between 0 and 1 -- remember that the units are nanoseconds.')
                return True

            if np.mod(pulse.start_time, 1.e9/(1.0e6*self.instruments['PB']['instance'].settings['clock_speed'])) != 0:
                # print("hello!")
                # print(pulse.start_time)
                # print(self.instruments['PB']['instance'].settings['clock_speed'])
                # print(1.e9/(1.0e6*self.instruments['PB']['instance'].settings['clock_speed']))
                return True
            if np.mod(pulse.start_time + pulse.duration, 1.e9/(1.0e6*self.instruments['PB']['instance'].settings['clock_speed'])) != 0:
                # print("???")
                # print(pulse.start_time + pulse.duration)
                # print(self.instruments['PB']['instance'].settings['clock_speed'])
                # print(1.e9 / (1.0e6 * self.instruments['PB']['instance'].settings['clock_speed']))
                return True

        return False

    def _compiles_to_too_many_commands_for_pb(self, pulse_sequence, verbose=False):
        """
        Checks to see if the current pulse sequence would result in more than 4096 commands to the pulseblaster,
        which is too much for it to handle.

        Args:
            pulse_sequence: A list of pulse objects to check for compilation issue
            verbose(bool): if True, calls extra print() statements during execution, helpful for debugging.

        Returns:
            (bool) True if the pulse sequence results into too many (4096) commands, else False

        """
        pulse_blaster = self.instruments['PB']['instance']
        delayed_pulse_collection = pulse_blaster.create_physical_pulse_seq(pulse_sequence) # ER 20181027 don't have time to debug why these delays are giving errors today
        pb_state_changes = pulse_blaster.generate_pb_sequence(delayed_pulse_collection)
        pb_commands = pulse_blaster.create_commands(pb_state_changes, self.settings['num_averages'])

        if len(pb_commands) < 4096:
            return False
        else:
            if verbose:
                print('Unfortunately compiled a pulse_sequence into too many pulse blaser commands, cannot program it')
            return True

    def _is_bad_pulse_sequence(self, pulse_sequence, verbose=True):
        """

        validates the pulse sequences, i.e. checks if the pulse sequences are compatible with all the constrains of the pulseblaster,
        e.g.
            - that pulses have to be multiples of 2.5ns
            - pulses have to be spaced at least 15ns apart
            - pulse can not be overlapping
            - that pulse sequences are compiled to fewer than 4096 pulse blaster commands (too many for the pulse blaster)

        Args:
            pulse_sequence(list): a list of Pulse objects to check for syndromes
            verbose(bool): if True, calls extra print() statements during execution, helpful for debugging.

        Returns:
            (bool) True if the pulse sequence contains one of the above issues, else False
        """

        if self._get_overlapping_pulses(pulse_sequence, verbose=verbose):
            return True
        elif self._get_commands_that_are_too_short(pulse_sequence, verbose=verbose):
            return True
        elif self._has_pulses_with_bad_start_time(pulse_sequence, verbose=verbose):
            return True
        elif self._compiles_to_too_many_commands_for_pb(pulse_sequence, verbose=verbose):
            return True
        else:
            return False

    def _combine_pulses(self, pulse_sequence, channel_id, overlap_window):
        """

        combines overlapping pulses in pulse_sequence with channeld id "channeld_id"

        Args:
            pulse_sequence: pulse sequence, list of Pulse object
            channel_id: id of channel where we combine the pulses if they overlap
            overlap_window: if pulses are closer than the "overlap_window" time window, they are considered overlapping and will be combined

        Returns: new pulse sequence, where the overlapping pulses have been combined

        """

        # split the sequence in the channels that we want to combine and the ones that we keep
        sequence_id = [pulse for pulse in pulse_sequence if
                       pulse.channel_id == channel_id]  # the sequences belonging to channel_id
        sequence_remainder = [pulse for pulse in pulse_sequence if
                              pulse.channel_id != channel_id]  # the sequences not belonging to channel_id

        # sort
        sequence_id = sorted(sequence_id, key=lambda pulse: pulse.start_time)

        # since the length of sequence_id shrinks over time we check on every iteration if we reached the end of the list
        # instead of doing a for index in range(len(sequence_id)) loop
        index = 0
        while True:
            if index >= len(sequence_id) - 1:
                break
            first_pulse = sequence_id[index]
            second_pulse = sequence_id[index + 1]
            if ((second_pulse.start_time) - (first_pulse.start_time + first_pulse.duration)) < overlap_window:

                # the combined pulse duration is the max of either the first pulse duration or
                # the differnence between the end of the second pulse and the start of the first pulse
                # the first case seems a bit unusual but can happen, for mw switch pulses where both
                # the I and Q channel are mapped onto the same channel (mw_switch)
                combined_pulse_duration = max(
                    first_pulse.duration,
                    (second_pulse.start_time - first_pulse.start_time) + second_pulse.duration
                )

                combined_pulse = Pulse(channel_id, first_pulse.start_time, combined_pulse_duration)
                sequence_id.remove(first_pulse)
                sequence_id.remove(second_pulse)
                sequence_id.insert(index, combined_pulse)
            else:
                index += 1

        # "adding" list in python concatenates them!
        return sequence_id + sequence_remainder

    def _add_mw_switch_to_sequences(self, pulse_sequences):
        """
        Adds the microwave switch to a sequence by toggling it on/off for every microwave_i or microwave_q pulse,
        with a buffer given by mw_switch_extra_time
        Args:
            pulse_sequences: Pulse sequence without mw switch
        Returns: Pulse sequence with mw switch added in appropriate places
        """
        gating = self.settings['mw_switch']['gating']
        mw_switch_time = self.settings['mw_switch']['extra_time']

        pulse_sequences_with_mw_switch = []
        for pulse_sequence in pulse_sequences:
            mw_switch_pulses = []
            # add a switch pulse for each microwave pulse, pulses are carved with i and q channels and wide mw switch pulses are added to surpress leakage
            if gating == 'mw_iq':
                for pulse in pulse_sequence:
                    if pulse.channel_id in ['microwave_i', 'microwave_q', 'microwave_i_2', 'microwave_q_2']:
                        mw_switch_pulses.append(Pulse('microwave_switch', pulse.start_time - mw_switch_time, pulse.duration + 2 * mw_switch_time))

                # add the mw switch pulses to the pulse sequences
                pulse_sequence.extend(mw_switch_pulses)

                # combine overlapping pulses and those that are within 2*mw_switch_extra_time
                pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_switch', overlap_window= 2 * mw_switch_time)


            elif gating == 'mw_switch':
                # in the case gating == 'mw_switch', the pulse is carved with the mw switch
                # thus, we extend the duration of the i and q pulses by mw_switch_time before and after

                for index, pulse in enumerate(pulse_sequence):
                    if pulse.channel_id in ['microwave_i', 'microwave_q', 'microwave_i_2', 'microwave_q_2']:
                        mw_switch_pulses.append(Pulse('microwave_switch', pulse.start_time, pulse.duration))

                        # replace the i and q pulses with wider pulses
                        new_pulse = Pulse(pulse.channel_id, pulse.start_time - mw_switch_time, pulse.duration + 2 * mw_switch_time)
                        pulse_sequence.remove(pulse)
                        pulse_sequence.insert(index, new_pulse)

                # add the mw switch pulses to the pulse sequences
                pulse_sequence.extend(mw_switch_pulses)

                # combine overlapping pulses and those that are within 2*mw_switch_extra_time
                pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_i', overlap_window= 2 * mw_switch_time)
                pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_i_2', overlap_window=2 * mw_switch_time)
                pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_q', overlap_window=2 * mw_switch_time)
                pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_q_2',overlap_window=2 * mw_switch_time)
                #pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_switch', overlap_window=2 * mw_switch_time)
            pulse_sequences_with_mw_switch.append(pulse_sequence)

        return pulse_sequences_with_mw_switch

    def _add_rf_switch_to_sequences(self, pulse_sequences):
        """
        Adds the RF switch to a sequence by toggling it on/off for every rf_i or rf_q pulse,
        with a buffer given by rf_switch_extra_time
        Args:
            pulse_sequences: Pulse sequence without RF switch
        Returns: Pulse sequence with RF switch added in appropriate places
        """
        gating = self.settings['rf_switch']['gating']
        rf_switch_time = self.settings['rf_switch']['extra_time']

        pulse_sequences_with_rf_switch = []
        for pulse_sequence in pulse_sequences:
            rf_switch_pulses = []
            # add a switch pulse for each microwave pulse, pulses are carved with i and q channels and wide mw switch pulses are added to surpress leakage
            if gating == 'rf_iq':
                for pulse in pulse_sequence:
                    if pulse.channel_id in ['rf_i', 'rf_q']:
                        rf_switch_pulses.append(Pulse('rf_switch', pulse.start_time - rf_switch_time, pulse.duration + 2 * rf_switch_time))

                # add the mw switch pulses to the pulse sequences
                pulse_sequence.extend(rf_switch_pulses)

                # combine overlapping pulses and those that are within 2*mw_switch_extra_time
                pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='rf_switch', overlap_window= 2 * rf_switch_time)


            elif gating == 'rf_switch':
                # in the case gating == 'mw_switch', the pulse is carved with the mw switch
                # thus, we extend the duration of the i and q pulses by mw_switch_time before and after

                for index, pulse in enumerate(pulse_sequence):
                    if pulse.channel_id in ['rf_i', 'rf_q']:
                        rf_switch_pulses.append(Pulse('rf_switch', pulse.start_time, pulse.duration))

                        # replace the i and q pulses with wider pulses
                        new_pulse = Pulse(pulse.channel_id, pulse.start_time - rf_switch_time, pulse.duration + 2 * rf_switch_time)
                        pulse_sequence.remove(pulse)
                        pulse_sequence.insert(index, new_pulse)

                # add the mw switch pulses to the pulse sequences
                pulse_sequence.extend(rf_switch_pulses)

                # combine overlapping pulses and those that are within 2*mw_switch_extra_time
                pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='rf_i', overlap_window= 2 * rf_switch_time)
                pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='rf_q', overlap_window=2 * rf_switch_time)
                pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='rf_switch', overlap_window=2 * rf_switch_time)
            pulse_sequences_with_rf_switch.append(pulse_sequence)

        return pulse_sequences_with_rf_switch

    def _plot_validate(self, axes_list):
        """
        Preview pulse sequence by plotting first and last sequence to plots 1 and 0, respectively
        Args:
            axes_list: List containing axes to plot to

        Returns:

        """
        axis0 = axes_list[0]
        axis1 = axes_list[1]

        pulse_sequences, tau_list, _ = self.create_pulse_sequences(logging=False)

      #  if pulse_sequences[0]:
        if pulse_sequences:
            plot_pulses(axis0, pulse_sequences[0])
            axis0.set_title('Pulse Visualization for Minimum tau (tau = {:.0f} ns)'.format(tau_list[0]))

            plot_pulses(axis1, pulse_sequences[-1])
            axis1.set_title('Pulse Visualization for Maximum tau (tau = {:.0f} ns)'.format(tau_list[-1]))
        else:
            print('No pulse sequences passed validation!!!')

    def create_pulse_sequences(self, logging=True):
        """
        A function to create the pulse sequence.

        NOTE THAT this is different from _create_pulse_sequences(), which must be overwritten in scripts inheriting from this script
        in contrast this function here takes the output from _create_pulse_sequences() and cleans it up!

        Returns:
            pulse_sequences: a list of pulse sequences to be run
            tau_list: the list of tau times we vary
            measurement_gate_width: measurement window time
        """

        pulse_sequences, tau_list, measurement_gate_width = self._create_pulse_sequences()
        logging = self.verbose
        pulse_sequences = self.process_virtual_channels(pulse_sequences)

        """
        # Adding microwave switch
        if 'mw_switch' in self.settings and self.settings['mw_switch']['add']:
            if logging:
                self.log('Adding microwave switch to pulse sequences')
            pulse_sequences = self._add_mw_switch_to_sequences(pulse_sequences)
        if 'rf_switch' in self.settings and self.settings['rf_switch']['add']:
            if logging:
                self.log('Adding rf switch to pulse sequences')
            pulse_sequences = self._add_rf_switch_to_sequences(pulse_sequences)

        # look for bad pulses, i.e. that don't comply with the requirements, e.g. given the pulse-blaster specs or
        # requiring that pulses don't overlap
        """

        valid_pulse_sequences = []
        valid_tau_list = []
        invalid_tau_list = []
        for pulse_sequence, tau in zip(pulse_sequences, tau_list):
            if not self._is_bad_pulse_sequence(pulse_sequence):
                valid_pulse_sequences.append(pulse_sequence)
                valid_tau_list.append(tau)
            else:
                invalid_tau_list.append(tau)
                print('Invalid sequence is: ', pulse_sequence)

        if logging:
            if invalid_tau_list:
                self.log("The pulse sequences corresponding to the following tau's were *invalid*, thus will not be "
                         "included: " + str(invalid_tau_list), flag='reminder')
            else:
                self.log("All generated pulse sequences are valid. No tau times will be skipped.")

            self.log("{:d} different tau times have passed validation".format(len(valid_tau_list)))

        return valid_pulse_sequences, valid_tau_list, measurement_gate_width

    def stop(self):
        """
        Stop currently executed pulse blaster sequence
        NOT CURRENTLY WORKING, WILL CRASH PULSEBLASTER
        """


        super(PulsedExperimentGeneric, self).stop()

        # self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('Stopping PB generic')
        #
        # print('Pulsed exp generic end behavior')
        # self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        # self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        # # self.instruments['mw_gen']['instance'].update({'amplitude': self.instruments['PB']['instance'].settings['constant_microwave_heat_load']['mw_power']})
        # self.instruments['mw_gen']['instance'].update(
        #     {'frequency': self.instruments['PB']['instance'].settings['constant_microwave_heat_load']['mw_frequency']})
        # self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # self.instruments['mw_gen']['instance'].update({'enable_output': True})
        # self.instruments['PB']['instance'].mw_duty_cycle_loop()

    def _track_nv(self, scenario, counts_temp=None):
        """
        Performs action to track to NV
        :param scenario:
                force -> run FindNv without checking any condition
                before_block -> run FindNv if 'before_block' in settings is True
                check_threshold -> run FindNv if counts are below threshold
        :param counts_temp: counts from the pulse sequence, to be checked against the threshold fluorescence, only needed if scenario is 'check_threshold'
        :return: whether to abort script, if FindNv fails consecutively multiple times
        """
        abort_script = False

        if scenario == 'before_block' and self.settings['track_nv']['before_block']:
            self.log('Running FindNv before averaging block')
            self.scripts['find_nv'].run(verbose=True)
            self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']
        if scenario == 'force':
            self.log('Running FindNv on manual request')
            self.scripts['find_nv'].run(verbose=True)
            self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']
        elif scenario == 'check_threshold' and self.settings['track_nv']['below_threshold']:
            threshold = self.settings['track_nv']['threshold']
            init_fluor = self.settings['track_nv']['init_fluor']
            findnv_attempts = 0
            counts_unsatisfactory = (1 + (1 - threshold)) * init_fluor < counts_temp or threshold * init_fluor > counts_temp
            while counts_unsatisfactory:
                self.log('Counts are below threshold; running FindNv')
                self.scripts['find_nv'].run(verbose=True)
                self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']
                findnv_attempts += 1
                threshold *= 0.9
                counts_unsatisfactory = (1 + (1 - threshold)) * init_fluor < counts_temp or threshold * init_fluor > counts_temp
                if findnv_attempts >= 3:
                    self.log('FindNv was unsuccessful after 5 attempts. Aborting script', flag='error')
                    abort_script = True
                    break
                elif counts_unsatisfactory:
                    self.log('Counts still below threshold after FindNv, running FindNv again with lower threshold of %.1f kCt/s' % (threshold*init_fluor),
                             flag='reminder')

        return abort_script

    def process_virtual_channels(self, pulse_sequences):
        '''
        Helper function that converts virtual channels into physical channels while adding microwave switch gating
        properly. This function includes the processing step done in _add_mw_switch_to_sequences.
        TODO: incorporate _add_rf_switch_to_sequences in this function.
        The output is defaulted to microwave_+i.
        Outputting to -i means that we're sending in a pulse to microwave_polarity channel.
        Outputting to +q means that we're sending in a pulse to microwave_iq channel.
        Outputting to -q means that we're sending in a pulse to both microwave_polarity and microwave_iq.
        '''
        gating = self.settings['mw_switch']['gating']
        mw_switch_time = self.settings['mw_switch']['extra_time']

        processed_pulse_sequences = []
        for pulse_sequence in pulse_sequences:
            if 'mw_switch' in self.settings and self.settings['mw_switch']['add']:
                mw_switch_pulses = []
                for pulse in pulse_sequence:
                    if pulse.channel_id in ['microwave_+i', 'microwave_-i', 'microwave_+q', 'microwave_-q']:
                        if gating == 'mw_iq':
                            mw_switch_pulses.append(Pulse('microwave_switch', pulse.start_time - mw_switch_time,
                                                          pulse.duration + 2 * mw_switch_time))
                        if gating == 'mw_switch':
                            mw_switch_pulses.append(Pulse('microwave_switch', pulse.start_time, pulse.duration))
                            # replace the i and q pulses with wider pulses. Due to some flipping issue with the switch
                            # we're giving it more time in the initial flip to stabilize the output before mw_switch is turned on.
                            # new_pulse = Pulse(pulse.channel_id, pulse.start_time - 2 * mw_switch_time,
                            #                   pulse.duration + 3 * mw_switch_time)
                            # pulse.start_time = pulse.start_time - 2 * mw_switch_time
                            # pulse.duration = pulse.duration + 3 * mw_switch_time

                            new_pulse = Pulse(pulse.channel_id, pulse.start_time - 1 * mw_switch_time,
                                              pulse.duration + 2 * mw_switch_time)
                            pulse.start_time = pulse.start_time - 1 * mw_switch_time
                            pulse.duration = pulse.duration + 2 * mw_switch_time

                        pulse_sequence.extend(mw_switch_pulses)

                        # combine overlapping pulses and those that are within 2*mw_switch_extra_time
                        pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_+i',
                                                              overlap_window=2 * mw_switch_time)
                        pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_-i',
                                                              overlap_window=2 * mw_switch_time)
                        pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_+q',
                                                              overlap_window=2 * mw_switch_time)
                        pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_-q',
                                                              overlap_window=2 * mw_switch_time)
                        pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_switch',
                                                              overlap_window=2 * mw_switch_time)

            new_pulse_sequence = []
            for pulse in pulse_sequence:
                if pulse.channel_id == 'microwave_+i':
                    pass
                elif pulse.channel_id == 'microwave_-i':
                    new_pulse_sequence.append(Pulse('microwave_polarity', pulse.start_time, pulse.duration))
                elif pulse.channel_id == 'microwave_+q':
                    new_pulse_sequence.append(Pulse('microwave_iq', pulse.start_time, pulse.duration))
                elif pulse.channel_id == 'microwave_-q':
                    new_pulse_sequence.append(Pulse('microwave_iq', pulse.start_time, pulse.duration))
                    new_pulse_sequence.append(Pulse('microwave_polarity', pulse.start_time, pulse.duration))
                else:
                    new_pulse_sequence.append(pulse)
            processed_pulse_sequences.append(new_pulse_sequence)
        return processed_pulse_sequences


if __name__ == '__main__':
    """
    script = {}
    instr = {}
    script, failed, instr = Script.load_and_append({'ExecutePulseBlasterSequence': 'ExecutePulseBlasterSequence'},
                                                   script, instr)

    script, failed, instr = Script.load_and_append({'ExecutePulseBlasterSequence': {'class':'ExecutePulseBlasterSequence',
                                                                                    'package':'b26toolkit'}},
                                                   script, instr)


    print(script)
    print(failed)
    print(instr)
    """
    pass
