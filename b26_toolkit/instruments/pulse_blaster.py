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

from pylabcontrol.core import Instrument, Parameter
from collections import namedtuple
import numpy as np
import itertools, ctypes, datetime, time, warnings
from pylabcontrol.core.read_write_functions import get_config_value
import os

# Pulse = namedtuple('Pulse', ('channel_id', 'start_time', 'duration'))


class Pulse(object):
    """
    pulse object that defines a pulse that is sent to the microwave source

    """
    def __init__(self, channel_id, start_time, duration=None, end_time=None, amplitude=None):

        self.channel_id = channel_id
        self.start_time = start_time
        self.amplitude = amplitude

        if duration is None:
            assert end_time is not None
            self.end_time = end_time
        if end_time is None:
            assert duration is not None
            self.duration = duration

    def __str__(self):
        """


        Returns: a representation of the Pulse object as a string, e.g.  when calling print()

        """

        if('amplitude' in vars(self)):
            return 'Pulse(id = {:s}, start = {:0.1f}ns, end = {:0.1f}ns, duration = {:0.1f}ns, amplitude = {:0.1f}V)'.format(
                self.channel_id,
                self.start_time,
                self.end_time,
                self.duration,
                self.amplitude)
        else:
            return 'Pulse(id = {:s}, start = {:0.1f}ns, end = {:0.1f}ns, duration = {:0.1f}ns)'.format(
                self.channel_id,
                self.start_time,
                self.end_time,
                self.duration)

    def __repr__(self):
        """


        Returns: a representation of the Pulse object as a string, e.g.  when the default output

        """


        return self.__str__()

    @property
    def duration(self):
        return self._duration

    @duration.setter
    def duration(self, duration):
        assert duration > 0, 'pulse duration has to be of finite duration but is {:d}'.format(duration)

        self._duration = duration
        self._end_time = self.start_time + duration

    @property
    def end_time(self):
        return self._end_time

    @end_time.setter
    def end_time(self, end_time):
        assert end_time > self.start_time, 'pulse end time has to be later than pulse start time'

        self._end_time = end_time
        self._duration = end_time - self.start_time

    @property
    def amplitude(self):
        return self._amplitude

    @amplitude.setter
    def amplitude(self, amplitude):
        self._amplitude = amplitude

    @staticmethod
    def is_overlapping(pulse1, pulse2, dead_time = 0):
        """
        Returns True if pulse overlap, barely touching pulses are **not** overlapping

        Args:
            pulse1: pulse 1 (Pulse object)
            pulse1: pulse 2 (Pulse object)
            dead_time: additional time inbetween pulses, if pulses are closer than dead time, they are considered overlapping

        Returns:
            boolean

        """

        # this is the overlap condition for two pulses
#         is_overlapping = True if pulse1.start_time < pulse2.end_time and pulse2.start_time < pulse1.end_time else False
        if pulse1.start_time < pulse2.end_time+dead_time and pulse2.start_time < pulse1.end_time+dead_time:
            is_overlapping = True
        else:
            is_overlapping = False
        return is_overlapping


class PulseBlaster(Instrument):
    """
    This Instrument controls a SpinCore Pulseblaster
    """
    PBStateChange = namedtuple('PBStateChange', ('channel_bits', 'time'))
    PBCommand = namedtuple('PBCommand', ('channel_bits', 'duration', 'command', 'command_arg'))
    MAX_COMMANDS = 4096

    _PROBES = {}

    _DEFAULT_SETTINGS = Parameter([
        Parameter('output_0', [
            Parameter('channel', 0, int, 'channel to which laser is connected'),
            Parameter('status', False, bool, 'True if voltage is high to the laser, false otherwise'),
            Parameter('delay_time', 0.2, float, 'delay time between pulse sending time and laser switch on [ns]')
        ]),
        Parameter('output_1', [
            Parameter('channel', 1, int, 'channel to which the daq is connected to'),
            Parameter('status', True, bool, 'True if voltage is high to the daq, false otherwise'),
            Parameter('delay_time', 0.2, float, 'delay time between pulse sending time and daq acknowledgement [ns]')
        ]),
        Parameter('output_2', [
            Parameter('channel', 2, int, 'channel to which the the microwave p trigger is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the microwave p trigger, false otherwise'),
            Parameter('delay_time', 0.2, float, 'delay time between pulse sending time and microwave p trigger [ns]')
        ]),
        Parameter('output_3', [
            Parameter('channel', 100, int, 'channel to which the the microwave q trigger is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the microwave q trigger, false otherwise'),
            Parameter('delay_time', 0.2, float, 'delay time between pulse sending time and microwave q trigger [ns]')
        ]),
        Parameter('output_4', [
            Parameter('channel', 3, int, 'channel to which the microwave switch is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the microwave switch, false otherwise'),
            Parameter('delay_time', 0.2, float, 'delay time between pulse sending time and microwave switch [ns]')
        ]),
        Parameter('clock_speed', 500.0, [100.0, 250.0, 400.0, 500.], 'Clock speed of the pulse blaster [MHz]'),
        Parameter('min_pulse_dur', 2, [15, 20, 50, 2, 12], 'Minimum allowed pulse duration (ns)'),
        # Min pulse duration setting can be removed. Min pulse duration is 5 clock cycles. Can be even shorter with short pulse feature, but doesn't seem implemented in code
        Parameter('PB_type', 'PCI', ['PCI', 'USB'], 'Type of pulseblaster used')
    ])

    PULSE_PROGRAM = ctypes.c_int(0)
    LONG_DELAY_THRESHOLD = 500
    # LONG_DELAY_THRESHOLD = 640
    # Had some issues with setting it to 640 (why was it even 640??). We only need the long delay feature when 32 bits aren't enough to specify the length of
    # the pulse; for a 500 MHz clock, this is 8.59s (just as website specified: https://www.spincore.com/products/PulseBlasterESR-PRO/)


    PB_INSTRUCTIONS = {
        'CONTINUE': ctypes.c_int(0),
        'STOP': ctypes.c_int(1),
        'LOOP': ctypes.c_int(2),
        'END_LOOP': ctypes.c_int(3),
        'BRANCH': ctypes.c_int(6),
        'LONG_DELAY': ctypes.c_int(7),
        'WAIT': ctypes.c_int(8)
    }

    def __init__(self, name=None, settings=None):
        try:
            self.dll_path = get_config_value('PULSEBLASTER_DLL_PATH', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.txt'))
        except IOError:
            warnings.warn(" ** Pulseblaster DLL not found. If it should be present, check the path.")
            dll_path = None
            print(('Expected dll_path: ', self.dll_path))
            self.is_conneted = False
        try:
          #  self.pb = ctypes.windll.LoadLibrary(dll_path) commented AS and ER 20180503
            self.pb = ctypes.CDLL('spinapi64')
            self.prepare_function_calls()
        except WindowsError:
            self.is_conneted = False
            warnings.warn("Pulseblaster DLL not found. If it should be present, check the path:")
            print(('Expected dll_path: ', self.dll_path)) # ER 20201009

        self.is_conneted = False

        super(PulseBlaster, self).__init__(name, settings)
        self.estimated_runtime = None
        self.sequence_start_time = None


    def min_pulse_dur(self):
        """
        :return: Minimum pulse duration in ns (5 clock cycles)
        """
        return int(1e3/self.settings['clock_speed']*5)

    def prepare_function_calls(self):

        '''

        todo(arthur): add docstring

        ER attempt:

        prepares fields in self.pb for the input arguments of the ctypes functions used in this instrument file.

        '''
        self.pb.pb_get_version.restype = ctypes.c_char_p
        self.pb.pb_get_error.restype = ctypes.c_char_p
        self.pb.pb_read_status.restype = ctypes.c_int

        self.pb.pb_count_boards.restype = ctypes.c_int

        self.pb.pb_init.restype = ctypes.c_int

        self.pb.pb_select_board.argtype = ctypes.c_int
        self.pb.pb_select_board.restype = ctypes.c_int

        self.pb.pb_set_debug.argtype = ctypes.c_int
        self.pb.pb_set_debug.restype = ctypes.c_int

        self.pb.pb_set_defaults.restype = ctypes.c_int

        self.pb.pb_core_clock.argtype = ctypes.c_double
        self.pb.pb_core_clock.restype = ctypes.c_int

        self.pb.pb_write_register.argtype = (ctypes.c_int, ctypes.c_int)
        self.pb.pb_write_register.restype = ctypes.c_int

        self.pb.pb_start_programming.argtype = ctypes.c_int
        self.pb.pb_start_programming.restype = ctypes.c_int

        self.pb.pb_stop_programming.restype = ctypes.c_int

        self.pb.pb_start.restype = ctypes.c_int
        self.pb.pb_stop.restype = ctypes.c_int
        self.pb.pb_reset.restype = ctypes.c_int
        self.pb.pb_close.restype = ctypes.c_int

        self.pb.pb_inst_pbonly.argtype = (
            ctypes.c_uint,  # flags
            ctypes.c_int,  # inst
            ctypes.c_int,  # inst_data
            ctypes.c_double,  # length
        )

    def update(self, settings):
        # call the update_parameter_list to update the parameter list
        super(PulseBlaster, self).update(settings)
        assert self.pb.pb_init() == 0, 'Could not initialize the PulseBlaster on pb_init() command.'
        self.pb.pb_core_clock(ctypes.c_double(self.settings['clock_speed']))
        self.pb.pb_start_programming(self.PULSE_PROGRAM)
        start = self.pb.pb_inst_pbonly(ctypes.c_int((self.settings2bits() | 0xE00000)),
                                       self.PB_INSTRUCTIONS['CONTINUE'], ctypes.c_int(0), ctypes.c_double(int(5e8)))
        address = self.pb.pb_inst_pbonly(ctypes.c_int((self.settings2bits() | 0xE00000)),
                                         self.PB_INSTRUCTIONS['BRANCH'], ctypes.c_int(1), ctypes.c_double(int(5e8)))
        # print("address: {}".format(address))+

        # need three 1's (0xE00000) at the start of the bitstring sent to the PB card to get the card to work

        # See here for explanation: http://www.spincore.com/support/spinapi/using_spin_api_pb.shtml
        # assert address >= 0, "An invalid parameter was sent to the pulseblaster"
        self.pb.pb_stop_programming()
        self.pb.pb_start()
        assert (self.pb.pb_read_status() & 0b100) == 0b100, 'PulseBlaster did not begin running after start() called.'
        self.pb.pb_close()

    def settings2bits(self):
        # Convert settings of PB outputs to bit string (e.g. all outputs off correspond to a string of 0s)
        bits = 0
        for output, output_params in self.settings.items():
            if isinstance(output_params, dict) and 'channel' in output_params and 'status' in output_params:
                if output_params['status']:
                    bits |= 1 << output_params['channel']

        return bits

    def get_delay(self, channel_id):
        """
        Gets the delay for an inputted channel id, a name or number.
        Args:
            channel_id: channel for which you want the delay. Must be an integer or string

        Returns:
            the delay, in ns, for a given channel id

        """
        if isinstance(channel_id, str):
            channel_name = channel_id
            return self.settings[channel_name]['delay_time']

        elif isinstance(channel_id, int):
            channel_num = channel_id
            for key, value in self.settings.items():
                if isinstance(value, dict) and 'channel' in list(value.keys()) and value['channel'] == channel_num:
                    return self.settings[key]['delay_time']

        raise AttributeError('Could not find delay of channel name or number: {0}'.format(str(channel_id)))

    @staticmethod
    def find_overlapping_pulses(pulses, combine_channels=None):
        """
        Finds all overlapping pulses in a collection of pulses, and returns the clashing pulses.

        Args:
            pulses: An iterable collection of Pulse objects
            combine_channels: list of *two* channels that should not be overlapping, if empty list only pulses
                with the same channel_id can be overlapping.

        Returns:
            A list of length-2 tuples of overlapping pulses. Each pair of pulses has the earlier pulse in the first
            position

        """

        if combine_channels is None:
            combine_channels = set()

        def get_overlapping_pulses(pulse1, pulse2):
            """
            Returns overlapping pulses as a tuple if time pulse 1 and 2 overlap in time

            Args:
                pulse1: pulse 1 (Pulse object)
                pulse2: pulse 2 (Pulse object)

            Returns:
                overlapping pulses as a tuple if time intervals 1 and 2 overlap in time

            """

            if Pulse.is_overlapping(pulse1, pulse2):
                overlapping_pulses = [tuple(sorted([pulse1, pulse2], key=lambda pulse: pulse.start_time))]
            else:
                overlapping_pulses = []

            return overlapping_pulses

        overlapping_pulses = []
        channel_ids = set([p.channel_id for p in pulses])  # get all channel ids as a list

        # check for overlapping pulses in the combined channels
        pulse_list = [p for p in pulses if p.channel_id in combine_channels]  # get all the pulses with any of the channel id
        for pulse1, pulse2 in itertools.combinations(pulse_list, 2):
            if Pulse.is_overlapping(pulse1, pulse2):
                overlapping_pulses.append(tuple(sorted([pulse1, pulse2], key=lambda pulse: pulse.start_time)))

        # check for overlapping pulses in the each of the channels that are not combined
        for channel_id in (channel_ids - set(combine_channels)):
            pulse_list = [pulse for pulse in pulses if pulse.channel_id == channel_id]  # get all the pulses with channel id
            for pulse1, pulse2 in itertools.combinations(pulse_list, 2):
                if Pulse.is_overlapping(pulse1, pulse2):
                    overlapping_pulses.append(tuple(sorted([pulse1, pulse2], key=lambda pulse: pulse.start_time)))

        return overlapping_pulses

    def create_physical_pulse_seq(self, pulse_collection):
        """
        Creates the physical pulse sequence from a pulse_collection, adding delays to each pulse, and ensuring the first
        pulse starts at time = 0.

        Args:
            pulse_collection: An iterable collection of pulses, named tuples with (name, start_time, pulse_duration)

        Returns:
            A list of pulses.

        """

        # Start time of earliest pulse
        min_pulse_time = np.min([pulse.start_time for pulse in pulse_collection])

        assert min_pulse_time >= 0, 'pulse with negative start time detected, that is not a valid pulse'

        # Round pulse durations to be integer multiples of clock period
        clock_period = 1.0e3 / self.settings['clock_speed']  # clock period in ns

        # Add delays to each pulse
        delayed_pulse_collection = [Pulse(pulse.channel_id,
                                          np.round((pulse.start_time - self.get_delay(pulse.channel_id))/clock_period)*clock_period,
                                          np.round(pulse.duration/clock_period)*clock_period)
                                    for pulse in pulse_collection]

        # Make sure the pulses start at same time as min_pulse_time
        delayed_min_pulse_time = np.min([pulse.start_time for pulse in delayed_pulse_collection])
        if delayed_min_pulse_time < min_pulse_time:
            delayed_pulse_collection = [Pulse(pulse.channel_id,
                                                   np.round((pulse.start_time - delayed_min_pulse_time + min_pulse_time)/clock_period)*clock_period,
                                                   np.round(pulse.duration/clock_period)*clock_period)
                                        for pulse in delayed_pulse_collection]

        # Return the sorted list of pulses, sorted by when they start
        return delayed_pulse_collection

    def generate_pb_sequence(self, pulses):
        """
        Creates a (ordered) list of PBStateChange objects for use with the pulseblaster API from a collection
        of Pulse's. Specifically, generate_pb_sequence generates a corresponding sequence of (bitstring, duration)
        objects that indicate the channels to keep on and for how long before the next instruction.

        Args:
            pulses: An iterable collection of Pulse's

        Returns:
            A (ordered) list of PBStateChange objects that indicate what bitstrings to turn on at what times.

        """

        # Create a dictionary with key=time and val=list of channels to toggle at that time
        # Note: we do not specifically keep track of whether we toggle on or off at a given time (is not necessary)
        pb_command_dict = {}
        for pulse in pulses:
            pulse_channel = self._get_channel(pulse.channel_id)
            pb_command_dict.setdefault(pulse.start_time, []).append(
                1 << pulse_channel)  # bitshifts by channel number to create array of 1 (0) for each channel on (off)
            pulse_end_time = pulse.start_time + pulse.duration
            pb_command_dict.setdefault(pulse_end_time, []).append(1 << pulse_channel)

        # Make sure we have a command at time=0, the command to have nothing on.
        if 0 not in list(pb_command_dict.keys()):
            pb_command_dict[0] = 0

        # For each time, combine all of the channels we need to toggle into a single bit string, and add it to a
        # command list of PBStateChange objects
        pb_command_list = []
        for instruction_time, bit_strings in pb_command_dict.items():
            channel_bits = np.bitwise_xor.reduce(bit_strings)
            if channel_bits != 0 or instruction_time == 0:
                pb_command_list.append(self.PBStateChange(channel_bits, instruction_time))

        # sort the list by the time a command needs to be placed
        pb_command_list.sort(key=lambda x: x.time)

        def change_to_propagating_signal(state_change_collection):
            # COMMENT_ME

            propagating_state_changes = []
            for idx in range(0, len(state_change_collection) - 1):

                # for the first command, just take the bitstring of the first element
                if idx == 0:
                    new_channel_bits = state_change_collection[0].channel_bits

                # otherwise, xor with the previous bitstring we computed (in propagating_state_changes)
                else:
                    new_channel_bits = propagating_state_changes[idx - 1].channel_bits ^ state_change_collection[
                        idx].channel_bits

                time_between_change = state_change_collection[idx + 1].time - state_change_collection[idx].time
                propagating_state_changes.append(self.PBStateChange(new_channel_bits, time_between_change))

            return propagating_state_changes

        # change this list so that instead of absolute times, they are durations, and the bitstrings properly propagate
        # i.e., if we want to keep channel 0 on at t = 1000 but now want to turn on channel 4, our bitstring would be
        # 10001 = 17 for the command at this time.
        pb_command_list = change_to_propagating_signal(pb_command_list)

        return pb_command_list

    def _get_long_delay_breakdown(self, pb_state_change, command='CONTINUE', command_arg=0):
        """
        Returns a list of PBCommand objects that correspond to a given state change and command, properly splitting up
        the state change into LONG_DELAY commands if necessary. The outputted list will always have the given command as
        the *last* command in the list.

        Args:
            pb_state_change: a PBStateChange object to find the long_delay breakdown of
            command: the command you want to initiate for this time interval
            command_arg: the argument to that command.

        Returns:
            list of PBCommand objects corresponding to a given state change.

        """

        # if the state duration is less than the LONG_DELAY threshold, return the normally formatted command
        if pb_state_change.time < self.LONG_DELAY_THRESHOLD:
            instruction_list = [
                self.PBCommand(pb_state_change.channel_bits, pb_state_change.time, command, command_arg)]
            return instruction_list

        instruction_list = []
        (num_long_delays, remainder) = divmod(pb_state_change.time, self.LONG_DELAY_THRESHOLD)

        # If the remaining time after filling with LONG_DELAYs is smaller than the minimum instruction duration,
        # split one LONG_DELAY command into one 'CONTINUE' command and one command of the type passed in
        if remainder < self.min_pulse_dur() and num_long_delays > 2:
            num_long_delays -= 1
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits, self.LONG_DELAY_THRESHOLD,
                               command='LONG_DELAY', command_arg=num_long_delays))
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits, self.LONG_DELAY_THRESHOLD / 2,
                               command='CONTINUE', command_arg=0))
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits,
                               self.LONG_DELAY_THRESHOLD / 2 + remainder, command, command_arg))
            return instruction_list

        # if there aren't any more LONG_DELAYs after doing this, just send in two commands
        elif remainder < self.min_pulse_dur() and num_long_delays == 2:
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits, self.LONG_DELAY_THRESHOLD / 2,
                               command='CONTINUE', command_arg=0))
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits, self.LONG_DELAY_THRESHOLD / 2,
                               command='CONTINUE', command_arg=0))
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits, self.LONG_DELAY_THRESHOLD / 2,
                               command='CONTINUE', command_arg=0))
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits,
                               self.LONG_DELAY_THRESHOLD / 2 + remainder, command, command_arg))
            return instruction_list

        # if there aren't any more LONG_DELAYs after doing this, just send in two commands
        elif remainder < self.min_pulse_dur() and num_long_delays == 1:
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits, self.LONG_DELAY_THRESHOLD / 2,
                               command='CONTINUE', command_arg=0))
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits,
                               self.LONG_DELAY_THRESHOLD / 2 + remainder, command, command_arg))
            return instruction_list

        # you cannot call LONG_DELAY once, so you have to call two CONTINUE instead of you were going to
        # have it loop only once.
        elif num_long_delays == 1:
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits, self.LONG_DELAY_THRESHOLD / 2,
                               command='CONTINUE', command_arg=0))
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits, self.LONG_DELAY_THRESHOLD / 2,
                               command='CONTINUE', command_arg=0))
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits, remainder, command, command_arg))
            return instruction_list

        # Otherwise, just use the LONG_DELAYs command followed by the given command for a duration given by remainder
        else:
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits, self.LONG_DELAY_THRESHOLD,
                               command='LONG_DELAY', command_arg=num_long_delays))
            instruction_list.append(
                self.PBCommand(pb_state_change.channel_bits, remainder, command, command_arg))
            return instruction_list

    def create_commands(self, pb_state_changes, num_loops=1):
        """
        Creates a list of commands to program the pulseblaster with, assuming that the user wants to loop over the
        state changes indicated in pb_state_changes for num_loops number of times. This function properly figures out
        when to use the LONG_DELAY pulseblaster command vs. CONTINUE, and also leaves the pulseblaster back in its
        steady-state conditions (following self.instruments.settings) when finished.


        Args:
            pb_state_changes: An ordered collection of state changes for the pulseblaster
            num_loops: the number of times the user wants to loop through the state changes (inner loop)
        Returns:
            An ordered list of PBComnmand objects to program the pulseblaster with.

        """
        if len(pb_state_changes) == 0:
            return []

        pb_commands = []

        pb_commands += list(reversed(self._get_long_delay_breakdown(pb_state_changes[0], command='LOOP', command_arg=num_loops)))

        for index in range(1, len(pb_state_changes)-1):
            pb_commands += list(reversed(self._get_long_delay_breakdown(pb_state_changes[index], command='CONTINUE')))

        pb_commands += self._get_long_delay_breakdown(pb_state_changes[-1], command='END_LOOP', command_arg=0)

        # Revert PB outputs to original settings after finishing commands
        # e.g. if laser channel was ON before running the sequence, turn it back ON after sequence is done
        pb_commands.append(self.PBCommand(self.settings2bits(), 1000, command='BRANCH', command_arg=len(pb_commands)))

        return pb_commands


    @staticmethod
    def estimate_runtime(pulses, num_loops=1):
        """
        Estimates the number of milliseconds required to complete the given pulse_collection for a given number of loops

        Args:
            pulses: a collection of Pulse objects
            num_loops: the number of iterations of pulse_collection

        Returns: estimated milliseconds to complete the pulse sequence.

        """
        return float(num_loops * max([pulse.start_time + pulse.duration for pulse in pulses])) / 1E6

    def program_pb(self, pulse_collection, num_loops=1):
        """
        Programs the pulseblaster to perform the pulses in the given pulse_collection on the next time start_pulse_seq()
        is called. The pulse collection must contain at least 2 pulses. Currently, we do not support time resolution below
        15 ns (FF: is this still true?).

        Args:
            pulse_collection: A collection of Pulse objects
            num_loops: The number of times to perform the given pulse collection

        Returns:

        """

        # Check for errors in the given pulse_collection
        # assert len(pulse_collection) > 1, 'pulse program must have at least 2 pulses'
        assert num_loops < (1 << 20), 'cannot have more than 2^20 (approx 1 million) loop iterations'
        if self.find_overlapping_pulses(pulse_collection):
            print('Pulses sent to the pulseblaster: ', pulse_collection)
            raise AttributeError('Found overlapping pulses in given pulse collection')

        for pulse in pulse_collection:
            assert pulse.start_time == 0 or pulse.start_time > 1, \
                'Found a start time between 0 and 1. Remember pulse times are in nanoseconds!'
            assert pulse.duration > 1, \
                'Found a pulse duration less than 1. Remember durations are in nanoseconds, and you can\'t have a 0 duration pulse'

        # Process the pulse collection into a format compatible with the low-level spincore API
        delayed_pulse_collection = self.create_physical_pulse_seq(pulse_collection)
        self.estimated_runtime = self.estimate_runtime(delayed_pulse_collection, num_loops)
        pb_state_changes = self.generate_pb_sequence(delayed_pulse_collection)
        pb_commands = self.create_commands(pb_state_changes, num_loops)

        assert len(pb_commands) < PulseBlaster.MAX_COMMANDS, "Generated a number of commands too long for the pulseblaster!"

        for command in pb_commands:
            if command.duration < self.min_pulse_dur():
                print('less than min pulse duration!!')
                print(command)
                print(command.duration)
                print(pb_commands)
                raise RuntimeError("Detected command with duration < min_pulse_dur.")

        assert self.pb.pb_init() == 0, 'Could not initialize the pulseblaster on pb_init() command.'
        self.pb.pb_core_clock(ctypes.c_double(self.settings['clock_speed']))
        self.pb.pb_start_programming(self.PULSE_PROGRAM)

        for pb_instruction in pb_commands:
            # note tha
            # t change types to the appropriate c type, as well as set certain bits in the channel bits to 1 in
            # order to properly output the signal
            return_value = self.pb.pb_inst_pbonly(ctypes.c_uint(pb_instruction.channel_bits | 0xE00000),
                                                  self.PB_INSTRUCTIONS[pb_instruction.command],
                                                  ctypes.c_int(int(pb_instruction.command_arg)),
                                                  ctypes.c_double(pb_instruction.duration))

            # print(bin(pb_instruction.channel_bits | 0xE00000))
            # print(pb_instruction.duration)
            # print(pb_instruction.command_arg)
            # print(pb_instruction.duration)
            assert return_value >=0, 'There was an error while programming the pulseblaster'
        self.pb.pb_stop_programming()

    def start_pulse_seq(self):
        """
        Starts the pulse sequence programmed in program_pb, and confirms that the pulseblaster is running.

        Returns:

        """

        self.pb.pb_start()
        assert self.pb.pb_read_status() & 0b100 == 0b100, 'pulseblaster did not begin running after start() called.'

        if self.settings['PB_type'] == 'PCI':  # leave USB PB connection open to stop it later
            self.pb.pb_close()
        self.sequence_start_time = datetime.datetime.now()

    def stop_pulse_seq(self):
        """
        Stops pulse sequence. Only needs to be called for USB pulseblasters. Requires that the PB-computer connection
        is still open.

        Returns:

        """
        self.pb.pb_stop()
        self.pb.pb_close()

    def wait(self):
        # COMMENT_ME
        if self.estimated_runtime is None or self.sequence_start_time is None:
            return
        if (datetime.datetime.now() - self.sequence_start_time).microseconds/1000.0 < self.estimated_runtime:
            time.sleep(self.estimated_runtime/1000.0 - (datetime.datetime.now() - self.sequence_start_time).total_seconds())

        self.estimated_runtime = None
        self.sequence_start_time = None

    def _get_channel(self, channel_id):
        # COMMENT_ME
        if isinstance(channel_id, (int, float)):
            return channel_id
        elif isinstance(channel_id, str):
            if channel_id in list(self.settings.keys()) and isinstance(self.settings[channel_id], dict) and 'channel' in \
                    list(self.settings[channel_id].keys()):
                return self.settings[channel_id]['channel']
            else:
                raise AttributeError('Could not find channel with the following id: {0}'.format(channel_id))

        raise AttributeError(
            'channel id must be either an integer or a string. Instead, this was passed in: {0}'.format(channel_id))


class B26PulseBlaster(PulseBlaster):
    # COMMENT_ME
    _DEFAULT_SETTINGS = Parameter([
        Parameter('laser', [
            Parameter('channel', 0, int, 'channel to which laser is connected'),
            Parameter('status', False, bool, 'True if voltage is high to the laser, false otherwise'),
            Parameter('delay_time', 350.2, float, 'delay time between pulse sending time and laser switch on [ns]')
        ]),
        Parameter('apd_readout', [
            Parameter('channel', 1, int, 'channel to which the daq is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the daq, false otherwise'),
            Parameter('delay_time', 0.2, float, 'delay time between pulse sending time and daq acknowledgement [ns]')
        ]),
        Parameter('microwave_i', [
            Parameter('channel', 2, int, 'channel to which the microwave p trigger is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the microwave p trigger, false otherwise'),
            Parameter('delay_time', 0.2, float, 'delay time between pulse sending time and microwave p trigger [ns]')
        ]),
        Parameter('microwave_q', [
            Parameter('channel', 3, int, 'channel to which the microwave q trigger is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the microwave q trigger, false otherwise'),
            Parameter('delay_time', 0.2, float, 'delay time between pulse sending time and microwave q trigger [ns]')
        ]),
        Parameter('microwave_switch', [
            Parameter('channel', 4, int, 'channel to which the microwave switch is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the microwave switch, false otherwise'),
            Parameter('delay_time', 0.2, float, 'delay time between pulse sending time and microwave switch [ns]')
        ]),
        Parameter('atto_trig', [
            Parameter('channel', 11, int, 'channel used to trigger function generator which in turn is connected to an Attocube'),
            Parameter('status', False, bool, 'True if voltage is high to the off-channel, false otherwise'),
            Parameter('delay_time', 0, float, 'delay time between pulse sending time and atto_trig on [ns]')
        ]),
        Parameter('rf_i', [
            Parameter('channel', 6, int, 'channel to which the rf p trigger is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the off-channel, false otherwise'),
            Parameter('delay_time', 0, float, 'delay time between pulse sending time and off channel on [ns]')
        ]),
        Parameter('rf_switch', [
            Parameter('channel', 7, int, 'channel to which the rf switch is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the off-channel, false otherwise'),
            Parameter('delay_time', 0, float, 'delay time between pulse sending time and off channel on [ns]')
        ]),
        Parameter('microwave_i_2', [
            Parameter('channel', 8, int, 'channel to which the microwave p trigger is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the off-channel, false otherwise'),
            Parameter('delay_time', 0, float, 'delay time between pulse sending time and off channel on [ns]')
        ]),
        Parameter('microwave_q_2', [
            Parameter('channel', 9, int, 'channel to which the microwave q trigger is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the off-channel, false otherwise'),
            Parameter('delay_time', 0, float, 'delay time between pulse sending time and off channel on [ns]')
        ]),
        Parameter('spacer', [
            Parameter('channel', 10, int, 'placeholder channel to control the length of pulse sequences, do not actually use the physical output!'),
            Parameter('status', False, bool, 'True if voltage is high to the off-channel, false otherwise'),
            Parameter('delay_time', 0, float, 'delay time between pulse sending time and off channel on [ns]')
        ]),
        Parameter('red_laser', [
            Parameter('channel', 5, int,
                      'placeholder channel to control the length of pulse sequences, do not actually use the physical output!'),
            Parameter('status', False, bool, 'True if voltage is high to the off-channel, false otherwise'),
            Parameter('delay_time', 0, float, 'delay time between pulse sending time and off channel on [ns]')
        ]),
        Parameter('clock_speed', 500, [100, 250, 400, 500], 'Clock speed of the pulse blaster [MHz]'),
        Parameter('min_pulse_dur', 2, [15, 20, 50, 2, 12], 'Minimum allowed pulse duration (ns)'),
        Parameter('PB_type', 'PCI', ['PCI', 'USB'], 'Type of pulseblaster used')
    ])

    _PROBES = {}

    # def __init__(self, name=None, settings=None):
    #     #COMMENT_ME
    #     super(B26PulseBlaster, self).__init__(name, settings)


    # def update(self, settings):
    #     # call the update_parameter_list to update the parameter list
    #     # oh god this is confusing
    #     super(B26PulseBlaster, self).update(settings)

    def read_probes(self, key):
        pass

    def get_name(self, channel):
        #COMMENT_ME
        for key, value in self.settings:
            if 'channel' in list(value.keys()) and value['channel'] == channel:
                return key

        raise AttributeError('Could not find instrument name attached to channel {s}'.format(channel))


class PulseBlasterHwTrig(B26PulseBlaster):
    _DEFAULT_SETTINGS = Parameter([
        Parameter('laser', [
            Parameter('channel', 0, int, 'channel to which laser is connected'),
            Parameter('status', False, bool, 'True if voltage is high to the laser, false otherwise'),
            Parameter('delay_time', 350.2, float, 'delay time between pulse sending time and laser switch on [ns]')
        ]),
        Parameter('apd_readout', [
            Parameter('channel', 1, int, 'channel to which the daq is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the daq, false otherwise'),
            Parameter('delay_time', 0.2, float, 'delay time between pulse sending time and daq acknowledgement [ns]')
        ]),
        Parameter('microwave_i', [
            Parameter('channel', 2, int, 'channel to which the microwave p trigger is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the microwave p trigger, false otherwise'),
            Parameter('delay_time', 0.2, float, 'delay time betweepulse sending time and microwave p trigger [ns]')
        ]),
        Parameter('microwave_q', [
            Parameter('channel', 8, int, 'channel to which the microwave q trigger is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the microwave q trigger, false otherwise'),
            Parameter('delay_time', 0.2, float, 'delay time between pulse sending time and microwave q trigger [ns]')
        ]),
        Parameter('microwave_switch', [
            Parameter('channel', 3, int, 'channel to which the microwave switch is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the microwave switch, false otherwise'),
            Parameter('delay_time', 0.2, float, 'delay time between pulse sending time and microwave switch [ns]')
        ]),
        Parameter('atto_trig', [
            Parameter('channel', 4, int, 'channel used to trigger function generator which in turn is connected to an Attocube'),
            Parameter('status', False, bool, 'True if voltage is high to the off-channel, false otherwise'),
            Parameter('delay_time', 0, float, 'delay time between pulse sending time and atto_trig on [ns]')
        ]),
        Parameter('rf_i', [
            Parameter('channel', 5, int, 'channel to which the rf p trigger is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the off-channel, false otherwise'),
            Parameter('delay_time', 0, float, 'delay time between pulse sending time and off channel on [ns]')
        ]),
        Parameter('rf_switch', [
            Parameter('channel', 6, int, 'channel to which the rf switch is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the off-channel, false otherwise'),
            Parameter('delay_time', 0, float, 'delay time between pulse sending time and off channel on [ns]')
        ]),
        Parameter('microwave_i_2', [
            Parameter('channel', 7, int, 'channel to which the microwave p trigger is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the off-channel, false otherwise'),
            Parameter('delay_time', 0, float, 'delay time between pulse sending time and off channel on [ns]')
        ]),
        Parameter('microwave_q_2', [
            Parameter('channel', 9, int, 'channel to which the microwave q trigger is connected to'),
            Parameter('status', False, bool, 'True if voltage is high to the off-channel, false otherwise'),
            Parameter('delay_time', 0, float, 'delay time between pulse sending time and off channel on [ns]')
        ]),
        Parameter('spacer', [
            Parameter('channel', 10, int, 'placeholder channel to control the length of pulse sequences, do not actually use the physical output!'),
            Parameter('status', False, bool, 'True if voltage is high to the off-channel, false otherwise'),
            Parameter('delay_time', 0, float, 'delay time between pulse sending time and off channel on [ns]')
        ]),
        Parameter('clock_speed', 500, [100, 250, 400, 500], 'Clock speed of the pulse blaster [MHz]'),
        Parameter('min_pulse_dur', 2, [15, 20, 50, 2, 12], 'Minimum allowed pulse duration (ns)'),
        Parameter('PB_type', 'PCI', ['PCI', 'USB'], 'Type of pulseblaster used')
    ])

    _PROBES = {}

    def program_pb(self, pulse_collection_laser, pulse_collection_mw, num_loops_laser, num_loops_mw):
        """
        Programs the pulseblaster to perform the pulses in the given pulse_collection on the next time start_pulse_seq()
        is called. The pulse collection must contain at least 2 pulses. Currently, we do not support time resolution below
        15 ns (FF: is this still true?).

        Args:
            pulse_collection: A collection of Pulse objects
            num_loops: The number of times to perform the given pulse collection

        Returns:
        """

        # Check for errors in the given pulse_collection
        # assert len(pulse_collection) > 1, 'pulse program must have at least 2 pulses'
        assert num_loops_laser < (1 << 20), 'cannot have more than 2^20 (approx 1 million) loop iterations'
        assert num_loops_mw < (1 << 20), 'cannot have more than 2^20 (approx 1 million) loop iterations'

        for pulse_collection in [pulse_collection_laser, pulse_collection_mw]:
            if self.find_overlapping_pulses(pulse_collection):
                print('Pulses sent to the pulseblaster: ', pulse_collection)
                raise AttributeError('Found overlapping pulses in given pulse collection')
            for pulse in pulse_collection:
                assert pulse.start_time == 0 or pulse.start_time > 1, \
                    'Found a start time between 0 and 1. Remember pulse times are in nanoseconds!'
                assert pulse.duration > 1, \
                    'Found a pulse duration less than 1. Remember durations are in nanoseconds, and you can\'t have a 0 duration pulse'

        # Process the pulse collection into a format compatible with the low-level spincore API
        delayed_pulse_collection_laser = self.create_physical_pulse_seq(pulse_collection_laser)
        delayed_pulse_collection_mw = self.create_physical_pulse_seq(pulse_collection_mw)

        self.estimated_runtime = self.estimate_runtime(pulse_collection_laser, num_loops_laser) + self.estimate_runtime(pulse_collection_mw, num_loops_mw)
        pb_state_changes_laser = self.generate_pb_sequence(pulse_collection_laser)
        pb_state_changes_mw = self.generate_pb_sequence(pulse_collection_mw)
        pb_commands = self.create_commands(pb_state_changes_laser, pb_state_changes_mw, num_loops_laser, num_loops_mw)

        assert len(pb_commands) < 4096, "Generated a number of commands too long (>= 4096) for the pulseblaster!"

        for command in pb_commands:
            if command.duration < self.min_pulse_dur():
                print('less than min pulse duration!!')
                print(command)
                print(command.duration)
                print(pb_commands)
                raise RuntimeError("Detected command with duration < min_pulse_dur.")

        assert self.pb.pb_init() == 0, 'Could not initialize the pulseblaster on pb_init() command.'
        self.pb.pb_core_clock(ctypes.c_double(self.settings['clock_speed']))
        self.pb.pb_start_programming(self.PULSE_PROGRAM)

        for pb_instruction in pb_commands:
            # return_value = self.pb.pb_inst_pbonly(ctypes.c_int(pb_instruction.channel_bits | 0xE00000),
            #                                       self.PB_INSTRUCTIONS[pb_instruction.command],
            #                                       ctypes.c_int(int(pb_instruction.command_arg)),
            #                                       ctypes.c_double(int(pb_instruction.duration)))
            return_value = self.pb.pb_inst_pbonly(ctypes.c_int(pb_instruction.channel_bits),
                                                  self.PB_INSTRUCTIONS[pb_instruction.command],
                                                  ctypes.c_int(int(pb_instruction.command_arg)),
                                                  ctypes.c_double(int(pb_instruction.duration)))
            if return_value < 0:
                print('Pulseblaster instruction error!')
                print(pb_commands)
                print(pb_instruction.duration)
                print(pb_instruction)
            assert return_value >= 0, 'There was an error while programming the pulseblaster'
        self.pb.pb_stop_programming()

    def create_commands(self, pb_state_changes_laser, pb_state_changes_mw, num_loops_laser, num_loops_mw, verbose=False):
        """
        Creates a list of commands to program the pulseblaster with. In PulseBlasterHwTrig the code has a different structure than create_commands in
        PulseBlaster. Here we first loop through pb_state_changes_laser for num_loops_laser and then loop through pb_states_changes_mw for num_loops_mw.
        In both loops the PB waits for a hardware trigger before sending out pulses. Therefore, we end up with a synchronized pulse train of laser initialization
        pulses and then pi pulses.

        Not implemented yet: The conventional way of dealing with pulse delays (e.g. AOM delay) is to simply shift the pulse forward by the delay time. However,
        the PB cannot wait for a hardware trigger while sending out a pulse simultaneously. The current workaround is to simply shift forward the laser pulses
        by a num of pulse train periods longer than the delay time.

        Args:
            pb_state_changes_laser: An ordered collection of state changes for the PB for controlling the laser; this needs to fit inside one period of the pulse train
            pb_state_changes_mw: An ordered collection of state changes for the PB for controlling the laser; this needs to fit inside one period of the pulse train
            num_loops_laser: number of iterations of the laser pulse collection
            num_loops_mw: number of iterations of the MW pulse collection
            verbose: if True, prints out all the PB commands

        Returns:
            An ordered list of PBComnmand objects to program the pulseblaster with.
        """

        if len(pb_state_changes_laser) + len(pb_state_changes_mw) == 0:
            print('0 state changes requested; create_commands() did not produce any PB commands.')
            return []

        pb_commands = []

        # Briefly turn off all PB outputs before running sequence
        pb_commands.append(self.PBCommand(0, 1000, command='CONTINUE', command_arg=0))
        for j in range(num_loops_laser):
            pb_commands.append(self.PBCommand(0, 10, command='WAIT', command_arg=0))
            for index in range(0, len(pb_state_changes_laser) - 1):
                pb_commands += list(reversed(self._get_long_delay_breakdown(pb_state_changes_laser[index], command='CONTINUE')))
            pb_commands += self._get_long_delay_breakdown(pb_state_changes_laser[-1], command='CONTINUE', command_arg=0)

        loop_begin_index = len(pb_commands)
        pb_commands.append(self.PBCommand(0, 10, command='LOOP', command_arg=num_loops_mw))
        # Wait command needs a minimum of 5 clock cycles. Command arg is unused.
        pb_commands.append(self.PBCommand(0, 10, command='WAIT', command_arg=len(pb_commands)))
        for index in range(0, len(pb_state_changes_mw)-1):
            pb_commands += list(reversed(self._get_long_delay_breakdown(pb_state_changes_mw[index], command='CONTINUE')))
        pb_commands += self._get_long_delay_breakdown(pb_state_changes_mw[-1], command='END_LOOP', command_arg=loop_begin_index)

        pb_commands.append(self.PBCommand(0, 10, command='BRANCH', command_arg=1))

        if verbose:
            print('Output from create_commands:')
            for command in pb_commands:
                print(command)
        return pb_commands

    def start_pulse_seq(self):
        """
        Starts the pulse sequence programmed in program_pb, and confirms that the PulseBlaster is running.
        Unlike in the original PulseBlaster class, this runs pb_reset() before running pb_start() to ensure that the PB starts from the beginning of the seq

        Returns:

        """
        print('You are using PB HW Trig')
        self.pb.pb_reset()
        super(PulseBlasterHwTrig, self).start_pulse_seq()
        #
        #
        # self.pb.pb_start()
        # assert self.pb.pb_read_status() & 0b100 == 0b100, 'pulseblaster did not begin running after start() called.'
        #
        # if self.settings['PB_type'] == 'PCI':  # leave USB PB connection open to stop it later
        #     self.pb.pb_close()
        # self.sequence_start_time = datetime.datetime.now()


if __name__ == '__main__':
    inst = PulseBlasterHwTrig()
    inst.pb.pb_init()
    inst.pb.pb_stop()
    print('{:04b}'.format(inst.pb.pb_read_status()))
    inst.pb.pb_close()

if __name__ == '__main2__':
    for i in range(1):
        import time
        # time.sleep(.1)
        from b26_toolkit.instruments.awg import AFG3022C
        arb = AFG3022C()
        arb.update({'ch2_enable': True})
        arb.update({'ch2_frequency': 4e6})
        arb.update({'ch2_amplitude': 1.0})
        arb.update({'ch2_run_mode': 'Burst'})
        arb.update({'ch2_burst_ncycles': int(10000000)})
        arb.afg.write('TRIG:SOUR EXT')


        inst = B26PulseBlaster()
        inst.update({'atto_trig': {'status': False}})
        # inst.update({'spacer': {'delay_time': 0.0}})
        # inst.update({'atto_trig': {'status': False}})

        # del inst
        # inst = B26PulseBlaster()
        # inst.update(inst.settings)



        # inst.pb.pb_init()
        # print('{:04b}'.format(inst.pb.pb_read_status()))
        # inst.pb.pb_core_clock(ctypes.c_double(inst.settings['clock_speed']))
        # inst.pb.pb_start_programming(inst.PULSE_PROGRAM)
        #
        # inst.pb.pb_inst_pbonly(ctypes.c_int(16), inst.PB_INSTRUCTIONS['CONTINUE'], ctypes.c_int(0), ctypes.c_double(int(1000)))
        # inst.pb.pb_inst_pbonly(ctypes.c_int(16), inst.PB_INSTRUCTIONS['BRANCH'], ctypes.c_int(1), ctypes.c_double(int(1000)))
        #
        # inst.pb.pb_stop_programming()
        # inst.pb.pb_start()
        # inst.pb.pb_close()



        pulse_collection = [Pulse('atto_trig', 10, int(80)), Pulse('rf_i', 10, int(160))]

        inst.program_pb(pulse_collection, num_loops=1)

        print('{:04b}'.format(inst.pb.pb_read_status()))
        inst.pb.pb_reset()
        inst.pb.pb_start()
        inst.pb.pb_close()

        time.sleep(.1)

        arb.afg.write('TRIG:SEQ:IMM')



        # inst.pb.pb_init()
        # print('{:04b}'.format(inst.pb.pb_read_status()))
        # inst.pb.pb_core_clock(ctypes.c_double(inst.settings['clock_speed']))
        # inst.pb.pb_start_programming(inst.PULSE_PROGRAM)
        #
        # inst.pb.pb_inst_pbonly(ctypes.c_int(0), inst.PB_INSTRUCTIONS['CONTINUE'], ctypes.c_int(0), ctypes.c_double(int(1000)))
        # inst.pb.pb_inst_pbonly(ctypes.c_int(0), inst.PB_INSTRUCTIONS['WAIT'], ctypes.c_int(0), ctypes.c_double(int(1000)))
        # inst.pb.pb_inst_pbonly(ctypes.c_int(16), inst.PB_INSTRUCTIONS['LOOP'], ctypes.c_int(7), ctypes.c_double(int(100)))
        # inst.pb.pb_inst_pbonly(ctypes.c_int(0), inst.PB_INSTRUCTIONS['CONTINUE'], ctypes.c_int(0), ctypes.c_double(int(1000)))
        # # inst.pb.pb_inst_pbonly(ctypes.c_int(16), inst.PB_INSTRUCTIONS['CONTINUE'], ctypes.c_int(0), ctypes.c_double(int(1e8)))
        # inst.pb.pb_inst_pbonly(ctypes.c_int(0), inst.PB_INSTRUCTIONS['END_LOOP'], ctypes.c_int(2), ctypes.c_double(int(100)))
        # inst.pb.pb_inst_pbonly(ctypes.c_int(0), inst.PB_INSTRUCTIONS['BRANCH'], ctypes.c_int(5), ctypes.c_double(int(1000)))
        #
        # inst.pb.pb_stop_programming()
        # inst.pb.pb_reset()
        # inst.pb.pb_start()
        # inst.pb.pb_close()


        # inst.pb.pb_init()
        # print('{:04b}'.format(inst.pb.pb_read_status()))
        # inst.pb.pb_core_clock(ctypes.c_double(inst.settings['clock_speed']))
        # inst.pb.pb_start_programming(inst.PULSE_PROGRAM)
        #
        # inst.pb.pb_inst_pbonly(ctypes.c_int(0), inst.PB_INSTRUCTIONS['CONTINUE'], ctypes.c_int(0), ctypes.c_double(int(1000)))
        # inst.pb.pb_inst_pbonly(ctypes.c_int(0), inst.PB_INSTRUCTIONS['LOOP'], ctypes.c_int(7), ctypes.c_double(int(1000)))
        # inst.pb.pb_inst_pbonly(ctypes.c_int(0), inst.PB_INSTRUCTIONS['WAIT'], ctypes.c_int(0), ctypes.c_double(int(100)))
        # inst.pb.pb_inst_pbonly(ctypes.c_int(0), inst.PB_INSTRUCTIONS['CONTINUE'], ctypes.c_int(0), ctypes.c_double(int(1000)))
        # # inst.pb.pb_inst_pbonly(ctypes.c_int(16), inst.PB_INSTRUCTIONS['CONTINUE'], ctypes.c_int(0), ctypes.c_double(int(1e8)))
        # inst.pb.pb_inst_pbonly(ctypes.c_int(16), inst.PB_INSTRUCTIONS['END_LOOP'], ctypes.c_int(1), ctypes.c_double(int(100)))
        # inst.pb.pb_inst_pbonly(ctypes.c_int(0), inst.PB_INSTRUCTIONS['BRANCH'], ctypes.c_int(5), ctypes.c_double(int(1000)))
        #
        # inst.pb.pb_stop_programming()
        # inst.pb.pb_reset()
        # inst.pb.pb_start()
        # inst.pb.pb_close()