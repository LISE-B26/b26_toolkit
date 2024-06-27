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
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_generic import PulsedExperimentGeneric
from b26_toolkit.instruments import NI6259, NI9402, PulseBlasterHwTrig, MicrowaveGenerator, Pulse, Commander
from b26_toolkit.plotting.plots_1d import plot_pulses, update_pulse_plot, plot_1d_simple_timetrace, update_1d_simple
from pylabcontrol.core import Parameter, Script
import itertools, ctypes, time, warnings
from copy import deepcopy
import random
import datetime
import time as t


class PiPulseTrainBackaction(PulsedExperimentGeneric):
    """
    Runs a MW pi-pulse train interrutped by laser reinitialization pulses. No readout.
    """

    _DEFAULT_SETTINGS = [
        Parameter('train_freq', 1.5e6, float, 'freq of pulse train in Hz. To drive the mechanics with pi pulses, this should be twice the mech freq'),
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'frequency of hyperfine transition'),
            Parameter('microwave_channel', 'i', ['i', 'q', 'alternate i and q'], 'Channel to use for mw pulses, UNIMPLEMENTED'),
            Parameter('tau_mw', 80, float, 'the time duration of the microwaves in ns. Make sure that 1/2f < tau_laser < 1/f.'),
            Parameter('n', 201, int, 'num of pi pulses, must be odd number')
        ]),
        Parameter('laser_pulses', [
            Parameter('tau_laser', 500, float, 'laser pulse duration. Make sure that 1/2f < tau_laser < 1/f, where f = train_freq, to ensure that '
                                               'the PB does not get triggered twice during the LOW part of the square wave going into the HW trigger port'),
            Parameter('n', 3, int, 'num of laser pulses, must be odd number')
        ]),
        Parameter('loop_delay', 1, float, 'delay in s between loops, one loop consisting of the full train of laser and pi pulses'),
        Parameter('num_averages', 1, int, 'number of averages')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB_hw_trig': PulseBlasterHwTrig, 'mw_gen': MicrowaveGenerator, 'commander': Commander}
    from b26_toolkit.scripts.find_nv import FindNvStrobe
    _SCRIPTS = {'find_nv': FindNvStrobe}

    def __init__(self, instruments, scripts, name=None, settings=None, log_function=None, data_path=None):
        """
        Standard script initialization
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_generic import PulsedExperimentGeneric
        self._DEFAULT_SETTINGS += PulsedExperimentGeneric._DEFAULT_SETTINGS
        # self._DEFAULT_SETTINGS['averaging_block_size'] = 1
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.ref_index = 0
        self.instruments['PB'] = self.instruments['PB_hw_trig']

    def _configure_instruments_start_of_script(self):
        """
        Configure instruments at the beginning of a script, e.g. enable IQ modulation
        :return: None
        """

        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

    def _create_pulse_sequences(self):
        """
        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement
        """

        tau_mech = 1/self.settings['train_freq']*1e9
        mw_tau = self.settings['mw_pulses']['tau_mw']
        n_pi = self.settings['mw_pulses']['n']
        laser_tau = self.settings['laser_pulses']['tau_laser']
        n_reset = self.settings['laser_pulses']['n']

        if self.settings['mw_pulses']['microwave_channel'] == 'alternate i and q':
            alternate = True
            microwave_channel_1 = 'microwave_i'
            microwave_channel_2 = 'microwave_q'
        elif self.settings['mw_pulses']['microwave_channel'] == 'alternate i_2 and q_2':
            alternate = True
            microwave_channel_1 = 'microwave_i_2'
            microwave_channel_2 = 'microwave_q_2'
        else:
            alternate = False
            microwave_channel_1 = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        tau = self.settings['mw_pulses']['tau_mw']
        pulse_sequences = []

        if (n_pi % 2) != 1:
            self.log('Error: number of pi pulses must be an odd number.', flag='error')
            self._abort = True

        if (n_reset % 2) != 1:
            self.log('Error: number of laser pulses must be an odd number.', flag='error')
            self._abort = True

        pulse_sequence = [Pulse('atto_trig', 10, mw_tau), Pulse('rf_i', 10, int(tau_mech/2*1.2/10)*10)]
        pulse_sequences.append(pulse_sequence)

        # Return pulse sequence and dummy variables
        return pulse_sequences, [tau], tau

    def _create_pulse_sequences_segment(self):
        """
        :return: [[pulse_sequence_laser], [pulse_sequence_mw]]
        pulse_sequence_laser, pulse_sequence_mw: a list of pulse sequences, each corresponding to a different time 'tau' that is to be scanned over.
        Each pulse sequence is a list of pulse objects containing the desired pulses.
        """
        
        tau_train = 1/self.settings['train_freq']*1e9
        tau_mw = self.settings['mw_pulses']['tau_mw']
        tau_laser = self.settings['laser_pulses']['tau_laser']
        n_pi = self.settings['mw_pulses']['n']
        n_reset = self.settings['laser_pulses']['n']

        if self.settings['mw_pulses']['microwave_channel'] == 'alternate i and q':
            alternate = True
            microwave_channel_1 = 'microwave_i'
            microwave_channel_2 = 'microwave_q'
        elif self.settings['mw_pulses']['microwave_channel'] == 'alternate i_2 and q_2':
            alternate = True
            microwave_channel_1 = 'microwave_i_2'
            microwave_channel_2 = 'microwave_q_2'
        else:
            alternate = False
            microwave_channel_1 = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        if (n_pi % 2) != 1:
            self.log('Error: number of pi pulses must be an odd number.', flag='error')
            self._abort = True
        if (n_reset % 2) != 1:
            self.log('Error: number of laser pulses must be an odd number.', flag='error')
            self._abort = True
        if tau_mw > int(tau_train * .8):
            self.log('Error: MW pulse duration longer than 80% of the period of the pulse train.', flag='error')
            self._abort = True
        if tau_laser > int(tau_train * .8):
            self.log('Error: Laser pulse duration longer than 80% of the period of the pulse train.', flag='error')
            self._abort = True

        spacer_pulse = Pulse('rf_i', 10, int(tau_train/2*1.2/10)*10)
        laser_delay = int(self.instruments['PB_hw_trig']['instance'].settings['laser']['delay_time'])
        mw_delay = int(self.instruments['PB_hw_trig']['instance'].settings['microwave_switch']['delay_time'])
        pulse_sequence_mw = [Pulse(microwave_channel_1, 10 + laser_delay, tau_mw), spacer_pulse]
        pulse_sequence_laser = [Pulse('laser', 10 + mw_delay, tau_laser), spacer_pulse]

        return [[pulse_sequence_laser], [pulse_sequence_mw]]

    def create_pulse_sequences_segment(self, logging=True):
        """
        A function to create the pulse sequence for either the laser or MW pulses
        Unlike create_pulse_sequences, this function does not automatically use the pulse sequences from _create_pulse_sequences() and requires explicit input

        NOTE THAT this is different from _create_pulse_sequences(), which must be overwritten in scripts inheriting from this script
        in contrast this function here takes the output from _create_pulse_sequences() and cleans it up!

        :param pulse_sequences_segment: a list of lists of pulse sequences to be processed, first list for laser pulses, second list for microwave pulses
        :param logging: if logging is enabled
        :return:
            pulse_sequences: a list of pulse sequences to be run
            tau_list: the list of tau times we vary
            measurement_gate_width: measurement window time
        """

        logging = self.verbose
        # pulse_sequences_segment_raw = self._create_pulse_sequences_segment()

        pulse_sequences_segment = []
        for pulse_sequences in self._create_pulse_sequences_segment():
            if 'mw_switch' in self.settings and self.settings['mw_switch']['add']:
                if logging:
                    self.log('Adding microwave switch to pulse sequences')
                pulse_sequences = self._add_mw_switch_to_sequences(pulse_sequences)
            if 'rf_switch' in self.settings and self.settings['rf_switch']['add']:
                if logging:
                    self.log('Adding rf switch to pulse sequences')
                pulse_sequences = self._add_rf_switch_to_sequences(pulse_sequences)
            pulse_sequences_segment.append(pulse_sequences)
        print('Length of pulse sequences segment: %i'%len(pulse_sequences_segment))

        # look for bad pulses, i.e. that don't comply with the requirements, e.g. given the pulse-blaster specs or requiring that pulses don't overlap

        # valid_pulse_sequences_segment = []
        # for pulse_sequences in pulse_sequences_segment:
        #     valid_pulse_sequences = []
        #     invalid_pulse_sequences = []
        #     for pulse_sequence in pulse_sequences:
        #         if not self._is_bad_pulse_sequence(pulse_sequence):
        #             valid_pulse_sequences.append(pulse_sequence)
        #         else:
        #             invalid_pulse_sequences.append(pulse_sequence)
        #             print('Invalid sequence is: ', pulse_sequence)
        #
        #     if logging:
        #         if len(invalid_pulse_sequences) == 0:
        #             self.log("All generated pulse sequences are valid. No tau times will be skipped.")
        #         else:
        #             self.log("Not all pulse sequences are valid. %i sequences will be skipped"%len(invalid_pulse_sequences))
        #
        #         self.log("{:d} different tau times have passed validation".format(len(valid_pulse_sequences)))
        #     valid_pulse_sequences_segment.append(valid_pulse_sequences)
        #     print('Length of valid pulse sequences segment: %i' % len(valid_pulse_sequences_segment))
        print('Valid pulse sequences:')
        print(pulse_sequences_segment)
        return pulse_sequences_segment

    def _run_single_sequence(self, pulse_sequence, num_loops, num_daq_reads):
        """
        Runs a single pulse sequence, num_loops consecutive times
        Args:
            pulse_sequence: a list of Pulse objects specifying a pulse sequence
            num_loops: number of times to repeat the pulse sequence
            num_daq_reads: number of times sequence requires that the

        Returns: a list containing, 1, 2, or 3 values depending on the pulse sequence counts, the second is the number of
        """

        if self.settings['daq_type'] == 'PCI':
            daq = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            daq = self.instruments['NI9402']['instance']

        try:
            pulse_sequence_laser, pulse_sequence_mw = self.create_pulse_sequences_segment()
            self.instruments['PB']['instance'].program_pb(pulse_sequence_laser[0], pulse_sequence_mw[0],
                                                          num_loops_laser=self.settings['laser_pulses']['n'],
                                                          num_loops_mw=self.settings['mw_pulses']['n'])
        except AssertionError:
            self.log('Error programming PB, aborting script')
            self._abort = True
            return np.zeros(num_daq_reads)

        if num_daq_reads != 0:
            task = daq.setup_gated_counter('ctr0', int(num_loops * num_daq_reads))
            daq.run(task)

        self.instruments['PB']['instance'].start_pulse_seq()

        result = []

        if num_daq_reads == 0:
            result = [0]

        if num_daq_reads != 0:
            result_array, temp = daq.read(task)   # thread waits on DAQ getting the right number of gates

            for i in range(num_daq_reads):
                result.append(sum(itertools.islice(result_array, i, None, num_daq_reads)))

        # clean up APD tasks
        if num_daq_reads != 0:
            daq.stop(task)

        if self.instruments['PB']['instance'].settings['PB_type'] == 'USB':
            print('stopping pulse seq: ')
            self.instruments['PB']['instance'].stop_pulse_seq()

        return result

    def is_valid(self, pulse_sequences=None, verbose=True):
        """
        Checks if a pulse blaster script is valid

        Args:
            create_pulse: is true creates the pulses first before validation
            verbose: print more stuff if true

        Returns:

        """
        self.verbose = verbose
        self.log('Warning: Validation unavailable for pulse sequences with hardware trigger, please check pulses on scope')

        if verbose:
            total_duration = 1.0e9/self.settings['train_freq']*(self.settings['laser_pulses']['n']+self.settings['mw_pulses']['n'])
            laser_duty = self.settings['laser_pulses']['tau_laser']*self.settings['laser_pulses']['n']/total_duration
            mw_duty = self.settings['mw_pulses']['tau_mw'] * self.settings['mw_pulses']['n'] / total_duration

            log_str = (total_duration*1e-9, laser_duty, mw_duty)
            self.log('Sequence duration: %.2e, laser duty cycle: %.2e, MW duty cycle: %.2e' %
                     log_str)
        return True

    def _plot_validate(self, axes_list):
        """
        Preview pulse sequence by plotting first and last sequence to plots 1 and 0, respectively
        Args:
            axes_list: List containing axes to plot to

        Returns:

        """
        pass

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

        # valid_pulse_sequences = [pulse_sequences]
        # valid_tau_list = []
        # invalid_tau_list = []
        # for pulse_sequence, tau in zip(pulse_sequences, tau_list):
        #     if not self._is_bad_pulse_sequence(pulse_sequence):
        #         valid_pulse_sequences.append(pulse_sequence)
        #         valid_tau_list.append(tau)
        #     else:
        #         invalid_tau_list.append(tau)
        #         print('Invalid sequence is: ', pulse_sequence)
        #
        # if logging:
        #     if invalid_tau_list:
        #         self.log("The pulse sequences corresponding to the following tau's were *invalid*, thus will not be "
        #                  "included: " + str(invalid_tau_list), flag='reminder')
        #     else:
        #         self.log("All generated pulse sequences are valid. No tau times will be skipped.")
        #
        #     self.log("{:d} different tau times have passed validation".format(len(valid_tau_list)))

        return pulse_sequences, tau_list, measurement_gate_width

    def stop(self):
        """
        Stop currently executed pulse blaster sequence
        NOT CURRENTLY WORKING, WILL CRASH PULSEBLASTER
        """
        # self.instruments['PB']['instance'].stop()
        super(PulsedExperimentGeneric, self).stop()
        print('STOP NOT IMPLEMENTED!!')




