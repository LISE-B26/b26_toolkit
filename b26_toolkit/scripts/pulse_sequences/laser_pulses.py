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

from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_generic import PulsedExperimentGeneric
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, Pulse, MicrowaveGenerator
from pylabcontrol.core import Parameter
import numpy as np


class LaserPulses(PulsedExperimentGeneric):
    """
    Sends a simple train of laser pulses without any readout window.
    """
    _DEFAULT_SETTINGS = [
        Parameter('laser_pulse_duration', 200, float, 'Length of each laser pulse (ns)'),
        Parameter('laser_off_time', 200, float, 'Time between each laser pulse (ns)'),
        Parameter('repetitions', 50, int, 'Number of repetitions'),
        Parameter('num_averages', 100000, int, 'number of averages')
    ]

    def _create_pulse_sequences(self, get_duration=False):
        """
        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement
        """

        tau = 200  # Placeholder, since we need to generate a tau list for the generic script to "sweep" through
        pulse_sequences = []
        tau_list = [tau]

        laser_pulse_duration = self.settings['laser_pulse_duration']
        laser_off_time = self.settings['laser_off_time']

        for tau in tau_list:
            pulse_sequence = []
            current_time = 0
            for part in range(int(self.settings['repetitions'])):  # Repeat pulse sequence 4 times because DAQ crashes if we run really short pulse sequences (~ several us long)
                # if tau is 0 there is actually no mw pulse
                if tau > 0:
                    pulse_sequence += [
                        Pulse('laser', current_time + laser_off_time, laser_pulse_duration)]
                current_time += laser_off_time + laser_pulse_duration

            pulse_sequences.append(pulse_sequence)

        if get_duration:
            # Return total sequence duration and total laser duration
            return current_time
        else:
            return pulse_sequences, tau_list, 100


class PulsedReadout(PulsedExperimentGeneric):
    """
    Single laser pulse with readout window
    """
    _DEFAULT_SETTINGS = [
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster}

    def _create_pulse_sequences(self, get_duration=False):

        '''

        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''

        tau = 200  # Placeholder, since we need to generate a tau list for the generic script to "sweep" through
        pulse_sequences = []
        tau_list = [tau]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']

        for tau in tau_list:
            pulse_sequence = []

            pulse_sequence += [Pulse('laser', laser_off_time, nv_reset_time)]
            pulse_sequence += [Pulse('apd_readout', laser_off_time + delay_readout, meas_time)]

            pulse_sequences.append(pulse_sequence)

        if get_duration:
            # Return total sequence duration and total laser duration
            duration = laser_off_time + np.min([nv_reset_time, delay_readout + meas_time])
            return duration
        else:
            return pulse_sequences, tau_list, self.settings['read_out']['meas_time']


class LaserPiPulses(PulsedExperimentGeneric):
    """
    Sends a MW pi-pulse followed by a laser reset pulse, with no readout window.
    Typically, this script runs in the background of another script (e.g. magnetic contour scan), the latter of which handles the readout
    """

    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('read_out', [
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
        ]),
        Parameter('repetitions', 100, int, 'number of repetitions of Pulsed ESR sequence consisting of MW pi-pulse and reinitialization'),
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _configure_instruments_start_of_script(self):
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_power']})
        # self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

    def _create_pulse_sequences(self, get_duration=False):
        """
        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement
        """

        tau = self.settings['tau_mw']
        pulse_sequences = []
        tau_list = [tau]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        microwave_channel = 'microwave_+i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        current_time = 0
        for tau in tau_list:
            pulse_sequence = []
            current_time = 0
            for part in range(int(self.settings['repetitions'])):  # Repeat pulse sequence 4 times because DAQ crashes if we run really short pulse sequences (~ several us long)
                # if tau is 0 there is actually no mw pulse
                if tau > 0:
                    if tau > 0:
                        pulse_sequence += [Pulse(microwave_channel, current_time + laser_off_time, tau)]

                if nv_reset_time > 0:
                    pulse_sequence += [Pulse('laser', current_time + laser_off_time + tau + delay_mw_readout, nv_reset_time)]
                current_time += laser_off_time + tau + delay_mw_readout + nv_reset_time

            pulse_sequences.append(pulse_sequence)

        if get_duration:
            # Return total sequence duration and total laser duration
            return current_time
        else:
            return pulse_sequences, tau_list, 100
