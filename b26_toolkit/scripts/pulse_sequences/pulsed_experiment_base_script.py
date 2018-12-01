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

import itertools
from copy import deepcopy

import numpy as np

from b26_toolkit.scripts import FindNV, ESR
from b26_toolkit.instruments import NI6259, B26PulseBlaster, Pulse, MicrowaveGenerator
from b26_toolkit.plotting.plots_1d import plot_1d_simple_timetrace_ns, plot_pulses, update_pulse_plot, update_1d_simple
from pylabcontrol.core import Script, Parameter
import random

MAX_AVERAGES_PER_SCAN = 100000  # 1E5, the max number of loops per point allowed at one time (true max is ~4E6 since
                                 #pulseblaster stores this value in 22 bits in its register


class PulsedExperimentBaseScript(Script):
    """
This class is a base class that should be inherited by all classes that utilize the pulseblaster for experiments. The
_function part of this class takes care of high-level interaction with the pulseblaster for experiment control and optionally
the daq for reading counter input (usually from the APD). It also provides all of the functionality needed to run a
standard Script such as plotting.
To use this class, the inheriting class need only overwrite _create_pulse_sequences to create the proper pulse sequence
for a given experiment
    """
    _DEFAULT_SETTINGS = [
        Parameter('Tracking', [
            Parameter('on/off', True, bool, 'used to turn on tracking'),
            Parameter('threshold', 0.85, float, 'threshold for tracking'),
            Parameter('init_fluor', 20., float, 'initial fluorescence of the NV to compare to, in kcps')
        ]),
        Parameter('ESR_Tracking', [
            Parameter('on/off', False, bool, 'turn on to track NV ESR'),
            Parameter('track_every_N', 10, int, 'track every N averages of the relevant sequence'),
            Parameter('allowed_delta_freq', 2., float, 'do not change the mw carrier frequency if the new ESR is different by more than this amount (MHz) (protects against bad fits)')
        ]),
        Parameter('randomize', True, bool, 'check to randomize runs of the pulse sequence'),
        Parameter('mw_switch', [
            Parameter('add', True, bool,  'check to add mw switch to every i and q pulse and to use switch to carve out pulses. note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 50, int, 'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'mw_switch', ['mw_switch', 'mw_iq'],'determines if mw pulses are carved out by mw-switch or by i and q channels of mw source '),
            Parameter('no_iq_overlap', True, bool,'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ])
    ]
    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster}

    _SCRIPTS = {'find_nv': FindNV, 'esr': ESR}

    def __init__(self, instruments, scripts, name=None, settings=None, log_function=None, data_path=None):
        """
        Standard script initialization
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        self._DEFAULT_SETTINGS += PulsedExperimentBaseScript._DEFAULT_SETTINGS

        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

    def _calc_progress(self, index):
        # progress of inner loop (in _run_sweep)
        progress_inner = index / len(self.pulse_sequences)

        # progress of outer loop (in _function)
        if self.current_averages >= MAX_AVERAGES_PER_SCAN:
            progress = (self.current_averages + (1.0 - progress_inner) * MAX_AVERAGES_PER_SCAN) / self.num_averages
        else:
            # if self.current_averages < MAX_AVERAGES_PER_SCAN, there is only a single run and
            # therefore the progress of the inner loop equals the outer loop
            progress = progress_inner

        self.progress = 100.0 * progress
        return int(round(self.progress))

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

        # make sure the microwave_switch is turned off so that we don't burn any steel cables. ER 20181017
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

        # retrieve initial mw carrier frequency to protect against bad fits in NV ESR tracking
        last_mw = self.scripts['esr'].instruments['microwave_generator']['instance'].frequency
        pulse_ampl = self.scripts['esr'].instruments['microwave_generator']['instance'].amplitude

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
        signal = [0.0]
        norms = np.repeat([0.0], (num_daq_reads - 1))
        self.count_data = np.repeat([np.append(signal, norms)], len(self.pulse_sequences), axis=0)
        self.data = in_data
        self.data['tau'] = np.array(self.tau_list)
        self.data['counts'] = deepcopy(self.count_data)

        # divides the total number of averages requested into a number of slices of MAX_AVERAGES_PER_SCAN and a remainder.
        # This is required because the pulseblaster won't accept more than ~4E6 loops (22 bits available to store loop
        # number) so need to break it up into smaller chunks (use 1E6 so initial results display faster)
        (num_1E5_avg_pb_programs, remainder) = divmod(self.num_averages, MAX_AVERAGES_PER_SCAN)
        # run find_nv if tracking is on ER 5/30/2017
        if self.settings['Tracking']['on/off']:
            self.scripts['find_nv'].run()
            if self.scripts['find_nv'].data['fluorescence'] == 0.0: # if it doesn't find an NV, abort the experiment
                self.log('Could not find an NV in FindNV.')
                self._abort = True
                return  # exit function in case no NV is found

        self.log("Averaging over {0} blocks of 1e5".format(num_1E5_avg_pb_programs))
        for average_loop in range(int(num_1E5_avg_pb_programs)):
            self.log("Running average block {0} of {1}".format(average_loop+1, int(num_1E5_avg_pb_programs)))
            if self._abort:
                self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
                self.log('aborted pulseblaster script during loop')
                break

            # ER 20181028
            if self.settings['ESR_Tracking']['on/off'] and average_loop % self.settings['ESR_Tracking']['track_every_N']==0:
                self.scripts['esr'].run()

                # retrieve the new mw frequency: if there are two frequencies in the fit, pick the one closest to the old frequency
                fit_params = self.scripts['esr'].data['fit_params']

                # default update flag to false
                update_mw = False

                if fit_params is not None and len(fit_params) and fit_params[0] != -1:  # check if fit valid
                    if len(fit_params) == 4:
                        # single peak
                        if (fit_params[2] - last_mw)**2 < (self.settings['ESR_Tracking']['allowed_delta_freq']*1e6)**2: # check if new value is within range allowed
                            update_mw = True
                        new_mw = fit_params[2]
                    elif len(fit_params) == 6:
                        # double peak, don't update the frequency - the fit may be bad
                        update_mw = False

                if update_mw:
                    #self.instruments['mw_gen'].update({'frequency': new_mw})
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'frequency': float(new_mw)})
                    self.log('updated mw carrier frequency to: {}'.format(new_mw))
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'amplitude': float(pulse_ampl)})
                    last_mw = new_mw
                else:
                    #self.instruments['mw_gen'].update({'frequency': last_mw})
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'frequency': float(last_mw)})
                    self.log('not updating the mw carrier frequency. SRS carrier frequency kept at {0} Hz'.format(last_mw))
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'amplitude': float(pulse_ampl)})

            self.current_averages = (average_loop + 1) * MAX_AVERAGES_PER_SCAN
            print('tau sequences running: ', self.tau_list)
            self._run_sweep(self.pulse_sequences, MAX_AVERAGES_PER_SCAN, num_daq_reads)
        if remainder != 0 and not self._abort:
            self.current_averages = self.num_averages
            self._run_sweep(self.pulse_sequences, remainder, num_daq_reads)

        if (len(self.data['counts'][0]) == 1) and not self._abort:
            self.data['counts'] = np.array([item for sublist in self.data['counts'] for item in sublist])

        # save data on the fly so that we can start to analyze it while the experiment is running!
        if self.settings['save']:
        #     self.save_b26()
            self.save_data()
        #     self.save_log()
        #     self.save_image_to_disk()

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

        if 'counts' in data.keys():
            # The following does not work for pulsedelays; you need to comment out the 'if' for it to work.
            # if counts != []:
            #     plot_1d_simple_timetrace_ns(axes_list[0], data['tau'], [data['counts'])
            print('plotting with tau values: ', data['tau'])
            plot_1d_simple_timetrace_ns(axes_list[0], data['tau'], [data['counts']])
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
        counts = self.data['counts']
        x_data = self.data['tau']
        axis1 = axes_list[0]
        if not counts == []:
            update_1d_simple(axis1, x_data, [counts])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])

    def _run_sweep(self, pulse_sequences, num_loops_sweep, num_daq_reads, verbose=False):
        """
        Each pulse sequence specified in pulse_sequences is run num_loops_sweep consecutive times.

        Args:
            pulse_sequences: a list of pulse sequences to run, each corresponding to a different value of tau. Each
                             sequence is a list of Pulse objects specifying a given pulse sequence
            num_loops_sweep: number of times to repeat each sequence before moving on to the next one
            num_daq_reads: number of times the daq must read for each sequence (generally 1, 2, or 3)

        Poststate: self.data['counts'] is updated with the acquired data

        """
        # ER 20180731 set init_fluor to zero
     #   self.data['init_fluor'] = 0.

        rand_indexes = list(range(len(pulse_sequences)))
        random.shuffle(rand_indexes)
        if verbose:
            print(('_run_sweep number of pulse sequences', len(pulse_sequences)))

        for index, sequence in enumerate(pulse_sequences):
            if verbose:
                print(('_run_sweep index', index, len(pulse_sequences)))

            rand_index = rand_indexes[index]
            if self._abort:
                self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
                break
            result = self._run_single_sequence(pulse_sequences[rand_index], num_loops_sweep, num_daq_reads)  # keep entire array
            self.count_data[rand_index] = self.count_data[rand_index] + result

            counts_to_check = self._normalize_to_kCounts(np.array(result), self.measurement_gate_width, num_loops_sweep)
            self.data['counts'][rand_index] = self._normalize_to_kCounts(self.count_data[rand_index], self.measurement_gate_width,
                                                                    self.current_averages)

            self.sequence_index = index
            counts_temp = counts_to_check[0]
            # track to the NV if necessary ER 5/31/17
            if self.settings['Tracking']['on/off']:
                if (1+(1-self.settings['Tracking']['threshold']))*self.settings['Tracking']['init_fluor'] < counts_temp or \
                        self.settings['Tracking']['threshold']*self.settings['Tracking']['init_fluor'] > counts_temp:
                    if verbose:
                        print('TRACKING TO NV...')
                    self.scripts['find_nv'].run()
                    self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']
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
        self.instruments['PB']['instance'].program_pb(pulse_sequence, num_loops=num_loops)
        # TODO(AK): figure out if timeout is actually needed
        timeout = 2 * self.instruments['PB']['instance'].estimated_runtime

        if num_daq_reads != 0:
            task = self.instruments['daq']['instance'].setup_gated_counter('ctr0', int(num_loops * num_daq_reads))
            self.instruments['daq']['instance'].run(task)

        self.instruments['PB']['instance'].start_pulse_seq()
        result = []
        if num_daq_reads != 0:
            result_array, _ = self.instruments['daq']['instance'].read(task)  # thread waits on DAQ getting the right number of gates
            for i in range(num_daq_reads):
                result.append(sum(itertools.islice(result_array, i, None, num_daq_reads)))
        # clean up APD tasks
        if num_daq_reads != 0:
            self.instruments['daq']['instance'].stop(task)

        return result

    # MUST BE IMPLEMENTED IN INHERITING SCRIPT
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

    def is_valid(self, pulse_sequences=None, verbose=False):
        """
        Checks if a pulse blaster script is valid

        Args:
            create_pulse: is true creates the pulses first before validation
            verbose: print more stuff if true

        Returns:

        """

        # make sure that we have the pulse sequences that correspond to the current settings

        if pulse_sequences is None:
            pulse_sequences, __, __ = self.create_pulse_sequences()

        for pulse_sequence in pulse_sequences:
            if self._is_bad_pulse_sequence(pulse_sequence):
                return False

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
        if self.settings['mw_switch']['no_iq_overlap']:
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
        pulse_blaster = self.instruments['PB']['instance']

        delayed_pulse_collection = pulse_blaster.create_physical_pulse_seq(pulse_sequence)

        pb_state_changes = pulse_blaster.generate_pb_sequence(delayed_pulse_collection)
        pb_commands = pulse_blaster.create_commands(pb_state_changes, self.settings['num_averages'])
        short_pulses = [command for command in pb_commands if command.duration < 15]

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
                return True
            if np.mod(pulse.start_time + pulse.duration, 1.e9/(1.0e6*self.instruments['PB']['instance'].settings['clock_speed'])) != 0:
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

    def _is_bad_pulse_sequence(self, pulse_sequence, verbose=False):
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
                    if pulse.channel_id in ['microwave_i', 'microwave_q']:
                        mw_switch_pulses.append(Pulse('microwave_switch', pulse.start_time - mw_switch_time, pulse.duration + 2 * mw_switch_time))

                # add the mw switch pulses to the pulse sequences
                pulse_sequence.extend(mw_switch_pulses)

                # combine overlapping pulses and those that are within 2*mw_switch_extra_time
                pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_switch', overlap_window= 2 * mw_switch_time)


            elif gating == 'mw_switch':
                # in the case gating == 'mw_switch', the pulse is carved with the mw switch
                # thus, we extend the duration of the i and q pulses by mw_switch_time before and after

                for index, pulse in enumerate(pulse_sequence):
                    if pulse.channel_id in ['microwave_i', 'microwave_q']:
                        mw_switch_pulses.append(Pulse('microwave_switch', pulse.start_time, pulse.duration))

                        # replace the i and q pulses with wider pulses
                        new_pulse = Pulse(pulse.channel_id, pulse.start_time - mw_switch_time, pulse.duration + 2 * mw_switch_time)
                        pulse_sequence.remove(pulse)
                        pulse_sequence.insert(index, new_pulse)

                # add the mw switch pulses to the pulse sequences
                pulse_sequence.extend(mw_switch_pulses)

                # combine overlapping pulses and those that are within 2*mw_switch_extra_time
                pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_i', overlap_window= 2 * mw_switch_time)
                pulse_sequence = self._combine_pulses(pulse_sequence, channel_id='microwave_q', overlap_window=2 * mw_switch_time)
            pulse_sequences_with_mw_switch.append(pulse_sequence)

        return pulse_sequences_with_mw_switch

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
            axis0.set_title('Pulse Visualization for Minimum tau (tau = {:f} ns)'.format(tau_list[0]))

            plot_pulses(axis1, pulse_sequences[-1])
            axis1.set_title('Pulse Visualization for Maximum tau (tau = {:f} ns)'.format(tau_list[-1]))
        else:
            print('no pulse sequences passed validation!!!')

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

        # Adding microwave switch
        if self.settings['mw_switch']['add']:
            if logging:
                self.log('Adding microwave switch to pulse sequences')
            pulse_sequences = self._add_mw_switch_to_sequences(pulse_sequences)

        # look for bad pulses, i.e. that don't comply with the requirements, e.g. given the pulse-blaster specs or
        # requiring that pulses don't overlap

        valid_pulse_sequences = []
        valid_tau_list = []
        invalid_tau_list = []
        for pulse_sequence, tau in zip(pulse_sequences, tau_list):
            if not self._is_bad_pulse_sequence(pulse_sequence):
                valid_pulse_sequences.append(pulse_sequence)
                valid_tau_list.append(tau)
            else:
                invalid_tau_list.append(tau)
                print('invalid sequence is: ', pulse_sequence)

        if logging:
            if invalid_tau_list:
                self.log("The pulse sequences corresponding to the following tau's were *invalid*, thus will not be "
                         "included when running this experiment: " + str(invalid_tau_list))
            else:
                self.log("All generated pulse sequences are valid. No tau times will be skipped in this experiment.")

            self.log("{:d} different tau times have passed validation".format(len(valid_tau_list)))

        return valid_pulse_sequences, valid_tau_list, measurement_gate_width

    def stop(self):
        """
        Stop currently executed pulse blaster sequence
        NOT CURRENTLY WORKING, WILL CRASH PULSEBLASTER
        """
        # self.instruments['PB']['instance'].stop()
        super(PulsedExperimentBaseScript, self).stop()


if __name__ == '__main__':
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
