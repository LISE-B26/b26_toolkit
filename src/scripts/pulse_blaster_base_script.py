"""
    This file is part of b26_toolkit, a PyLabControl add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    Foobar is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""

import itertools
from copy import deepcopy

import numpy as np

from b26_toolkit.src.scripts import FindNV
from b26_toolkit.src.instruments import NI6259, B26PulseBlaster, Pulse
from b26_toolkit.src.plotting.plots_1d import plot_1d_simple_timetrace_ns, plot_pulses, update_pulse_plot, update_1d_simple
from PyLabControl.src.core.scripts import Script, Parameter
import random

MAX_AVERAGES_PER_SCAN = 100000  # 1E5, the max number of loops per point allowed at one time (true max is ~4E6 since
                                 #pulseblaster stores this value in 22 bits in its register


class PulseBlasterBaseScript(Script):
    '''
This class is a base class that should be inherited by all classes that utilize the pulseblaster for experiments. The
_function part of this class takes care of high-level interaction with the pulseblaster for experiment control and optionally
the daq for reading counter input (usually from the APD). It also provides all of the functionality needed to run a
standard Script such as plotting.
To use this class, the inheriting class need only overwrite _create_pulse_sequences to create the proper pulse sequence
for a given experiment
    '''
    _DEFAULT_SETTINGS = [
        Parameter('Tracking', [
            Parameter('on/off', True, bool, 'used to turn on tracking'),
            Parameter('threshold', 0.85, float, 'threshold for tracking'),
        ]),
        Parameter('randomize', True, bool, 'check to randomize runs of the pulse sequence'),


    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster}

    _SCRIPTS = {'find_nv': FindNV}

    def __init__(self, instruments, scripts, name=None, settings=None, log_function=None, data_path=None):
        """
        Standard script initialization
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        self._DEFAULT_SETTINGS += PulseBlasterBaseScript._DEFAULT_SETTINGS

        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)


    def _calc_progress(self, index):
        # progress of inner loop (in _run_sweep)
        progress_inner = float(index) / len(self.pulse_sequences)

        # progress of outer loop (in _function)
        if self.current_averages >= MAX_AVERAGES_PER_SCAN:
            progress = float(self.current_averages + (progress_inner - 1.0) * MAX_AVERAGES_PER_SCAN) / self.num_averages
        else:
            progress = progress_inner
        self.progress = 100.0 * progress
        return int(progress)

    def _function(self, in_data={}):
        '''
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

        '''
        self.sequence_index = 0
        # self.pulse_sequences, self.num_averages, tau_list, self.measurement_gate_width = self._create_pulse_sequences()
        #wrapper

        self.pulse_sequences, self.num_averages, tau_list, self.measurement_gate_width = self._process_sequences()

        # print('number of sequences after validation ', len(self.pulse_sequences))

        #calculates the number of daq reads per loop requested in the pulse sequence by asking how many apd reads are
        #called for. if this is not calculated properly, daq will either end too early (number too low) or hang since it
        #never receives the rest of the counts (number too high)
        num_daq_reads = 0

        for pulse in self.pulse_sequences[0]:
            if pulse.channel_id == 'apd_readout':
                num_daq_reads += 1
        signal = [0.0]
        norms = np.repeat([0.0], (num_daq_reads - 1))
        self.count_data = np.repeat([np.append(signal, norms)], len(self.pulse_sequences), axis=0)
        self.data = in_data
        self.data.update({'tau': np.array(tau_list), 'counts': deepcopy(self.count_data)})
        #divides the total number of averages requested into a number of slices of MAX_AVERAGES_PER_SCAN and a remainer.
        #This is required because the pulseblaster won't accept more than ~4E6 loops (22 bits avaliable to store loop
        #number) so need to break it up into smaller chunks (use 1E6 so initial results display faster)
        (num_1E5_avg_pb_programs, remainder) = divmod(self.num_averages, MAX_AVERAGES_PER_SCAN)

        # run find_nv if tracking is on ER 5/30/2017
        if self.settings['Tracking']['on/off']:
            self.scripts['find_nv'].run()
            if self.scripts['find_nv'].data['fluorescence'] == 0.0: # if it doesn't find an NV, abort the experiment
                self._abort = 1
            #self._plot_refresh = True

        self.log("Averaging over {0} blocks of 1e5".format(num_1E5_avg_pb_programs))

        for average_loop in range(int(num_1E5_avg_pb_programs)):
            self.log("Running average block {0} of {1}".format(average_loop+1, int(num_1E5_avg_pb_programs)))
            if self._abort:
                break
            # print('loop ' + str(average_loop))
            self.current_averages = (average_loop + 1) * MAX_AVERAGES_PER_SCAN
            self._run_sweep(self.pulse_sequences, MAX_AVERAGES_PER_SCAN, num_daq_reads)

        if remainder != 0 and not self._abort:
            self.current_averages = self.num_averages
            self._run_sweep(self.pulse_sequences, remainder, num_daq_reads)

        if (len(self.data['counts'][0]) == 1) and not self._abort:
            self.data['counts'] = np.array([item for sublist in self.data['counts'] for item in sublist])

        # if self.settings['save']:
        #     self.save_b26()
        #     self.save_data()
        #     self.save_log()
        #     self.save_image_to_disk()

    def _plot(self, axes_list, data = None):
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
        counts = data['counts']
        x_data = data['tau']
        axis1 = axes_list[0]
    # The following does not work for pulsedelays; you need to comment out the 'if' for it to work.
    # if counts != []:
    #     plot_1d_simple_timetrace_ns(axis1, x_data, [counts])
        plot_1d_simple_timetrace_ns(axis1, x_data, [counts])
        axis2 = axes_list[1]
        plot_pulses(axis2, self.pulse_sequences[self.sequence_index])

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

    def _run_sweep(self, pulse_sequences, num_loops_sweep, num_daq_reads):
        '''
        Each pulse sequence specified in pulse_sequences is run num_loops_sweep consecutive times.

        Args:
            pulse_sequences: a list of pulse sequences to run, each corresponding to a different value of tau. Each
                             sequence is a list of Pulse objects specifying a given pulse sequence
            num_loops_sweep: number of times to repeat each sequence before moving on to the next one
            num_daq_reads: number of times the daq must read for each sequence (generally 1, 2, or 3)

        Poststate: self.data['counts'] is updated with the acquired data

        '''

        # randomize the indexes of the pulse sequences to run, to reduce heating. ER 5/25/2017
        rand_indexes = []
        for i in range(0, len(pulse_sequences)):
            rand_indexes.append(i)
        if self.settings['randomize']:
            random.shuffle(rand_indexes)
        for index, sequence in enumerate(pulse_sequences):
            rand_index = rand_indexes[index]
            if self._abort:
                break
            result = self._single_sequence(pulse_sequences[rand_index], num_loops_sweep, num_daq_reads)  # keep entire array
            self.count_data[rand_index] = self.count_data[rand_index] + result

            # emma 10/22/17: just for readout loop
            #self.measurement_gate_width = sequence[2][3]
            #print self.measurement_gate_width

            counts_to_check = self._normalize_to_kCounts(np.array(result), self.measurement_gate_width, num_loops_sweep)
            self.data['counts'][rand_index] = self._normalize_to_kCounts(self.count_data[rand_index], self.measurement_gate_width,
                                                                    self.current_averages)
            self.sequence_index = rand_index
            counts_temp = counts_to_check[0]
            # throw error if tracking is on and you haven't ran find nv ER 6/2/17
            if self.settings['Tracking']['on/off']:
                if self.scripts['find_nv'].data['fluorescence']:
                    self.data['init_fluor'] = deepcopy(self.scripts['find_nv'].data['fluorescence'])
                else:
                    raise AttributeError('need to run find NV first for tracking!')

            # track to the NV if necessary ER 5/31/17
            if (self.settings['Tracking']['on/off']):
                if (self.settings['Tracking']['threshold']*self.data['init_fluor'] > counts_temp or
                            (2-self.settings['Tracking']['threshold'])*self.data['init_fluor'] < counts_temp):
             #      self._plot_refresh = True
                    print('TRACKING TO NV...')
                    self.scripts['find_nv'].run()
              #     self._plot_refresh = True
                    self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']
            self.updateProgress.emit(self._calc_progress(index))

    def _single_sequence(self, pulse_sequence, num_loops, num_daq_reads):
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

        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences


        '''
        raise NotImplementedError



    def _normalize(self, signal, baseline_max=0, baseline_min=0):
        '''
        Normalizes the signal values given a maximum value (counts in |0>) and optionally minimum value (counts in |1>,
        set to 0 if none provided). If no maximum given, returns the bare signal.
        Args:
            signal:
            baseline_max:
            baseline_min:

        Returns:

        '''
        if baseline_max == 0:
            return signal
        else:
            return ((signal - baseline_min) / (baseline_max - baseline_min))

    def _normalize_to_kCounts(self, signal, gate_width=1, num_averages=1):
        """
        Converts the signal from counts/gate_width as returned by the PB to kcounts/s
        Args:
            signal: data to normalize
            gate_width: length of each bin
            num_averages: number of times each bin is read

        Returns:

        """
        return (1. * signal * (1E6 / (gate_width * num_averages)))  # 1E6 is to convert from ns to ms

    def validate(self):
        """
        Checks if a pulse sequence is a valid input to the pulseblaster, that is the sequence has no overlapping
        pulses, or time between pulses of <15ns
        """
        def add_mw_switch_to_sequences(pulse_sequences):
            """
            Adds the microwave switch to a sequence by toggling it on/off for every microwave_i or microwave_q pulse,
            with a buffer given by mw_switch_extra_time
            Args:
                pulse_sequences: Pulse sequence without mw switch

            Returns: Pulse sequence with mw switch added in appropriate places

            """
            if not 'mw_switch_extra_time' in self.settings.keys():
                #default to a 40 ns buffer
                mw_switch_time = 40
            for sequence in pulse_sequences:
                mw_switch_pulses = []
                # add a switch pulse for each microwave pulse
                for pulse in sequence:
                    if pulse.channel_id in ['microwave_i', 'microwave_q']:
                        mw_switch_pulses.append(Pulse('microwave_switch', pulse.start_time - mw_switch_time,
                                                      pulse.duration + 2 * mw_switch_time))
                # combine overlapping pulses and those that are within 2*mw_switch_extra_time
                mw_switch_pulses = sorted(mw_switch_pulses, key=lambda pulse: pulse.start_time)
                index = 0
                while True:
                    if index >= len(mw_switch_pulses) - 1:
                        break
                    first_pulse = mw_switch_pulses[index]
                    second_pulse = mw_switch_pulses[index + 1]
                    if ((second_pulse.start_time) - (
                        first_pulse.start_time + first_pulse.duration)) < 2 * mw_switch_time:
                        new_pulse = Pulse('microwave_switch', first_pulse.start_time, max(first_pulse.duration, (
                            second_pulse.start_time - first_pulse.start_time) + second_pulse.duration))
                        mw_switch_pulses.remove(first_pulse)
                        mw_switch_pulses.remove(second_pulse)
                        mw_switch_pulses.insert(index, new_pulse)
                    else:
                        index += 1
                sequence.extend(mw_switch_pulses)

            return pulse_sequences

        pulse_blaster = self.instruments['PB']['instance']
        self.pulse_sequences, num_averages, tau_list, measurement_gate_width = self._create_pulse_sequences()
        self.pulse_sequences = add_mw_switch_to_sequences(self.pulse_sequences) #UNCOMMENT TO ADD SWITCH
        failure_list = []
        for pulse_sequence in self.pulse_sequences:
            overlapping_pulses = B26PulseBlaster.find_overlapping_pulses(pulse_sequence)
            if not overlapping_pulses == []:
                failure_list.append(B26PulseBlaster.find_overlapping_pulses(pulse_sequence))
                break
            for pulse in pulse_sequence:
                assert pulse.start_time == 0 or pulse.start_time > 1, \
                    'found a start time that was between 0 and 1. Remember pulse times are in nanoseconds!'
                assert pulse.duration > 1, \
                    'found a pulse duration less than 1. Remember durations are in nanoseconds, and you can\'t have a 0 duration pulse'

            # process the pulse collection into a format that is designed to deal with the low-level spincore API
            delayed_pulse_collection = pulse_blaster.create_physical_pulse_seq(pulse_sequence)
            pb_state_changes = pulse_blaster.generate_pb_sequence(delayed_pulse_collection)
            pb_commands = pulse_blaster.create_commands(pb_state_changes, self.settings['num_averages'])

            assert len(pb_commands) < 4096, "Generated a number of commands too long for the pulseblaster!"

            short_pulses = [command for command in pb_commands if command.duration < 15]
            if short_pulses:
                failure_list.append(short_pulses[0])
            else:
                failure_list.append([])  # good sequence

        if any([isinstance(a, pulse_blaster.PBCommand) for a in failure_list]):
            self.log('Validation failed. At least one pulse in the sequence is invalid.')

        return self.pulse_sequences, num_averages, tau_list, measurement_gate_width, failure_list

    def _plot_validate(self, axes_list):
        """
        Preview pulse sequence by plotting first and last sequence to plots 1 and 2
        Args:
            axes_list: List containing axes to plot to

        Returns:

        """
        axis1 = axes_list[0]
        axis2 = axes_list[1]
        plot_pulses(axis2, self.pulse_sequences[0])
        plot_pulses(axis1, self.pulse_sequences[-1])

    def _process_sequences(self):
        """
        gets pulse sequences from self.validate
        Returns: pulse_sequences, num_averages, tau_list, measurement_gate_width
        """
        pulse_sequences, num_averages, tau_list, measurement_gate_width, failure_list = self.validate()
        delete_list = []
        for index, result in enumerate(failure_list):
            if not result == []:
                delete_list.append(index)
        if not delete_list == []:
            delete_list.reverse()
            for index in delete_list:
                pulse_sequences.pop(index)

        tau_list = np.delete(np.array(tau_list), delete_list).tolist()

        return pulse_sequences, num_averages, tau_list, measurement_gate_width

    def stop(self):
        """
        Stop currently executed pulse blaster sequence
        NOT CURRENTLY WORKING, WILL CRASH PULSEBLASTER
        """
        # self.instruments['PB']['instance'].stop()
        super(PulseBlasterBaseScript, self).stop()


if __name__ == '__main__':
    script = {}
    instr = {}
    script, failed, instr = Script.load_and_append({'ExecutePulseBlasterSequence': 'ExecutePulseBlasterSequence'},
                                                   script, instr)

    print(script)
    print(failed)
    print(instr)
