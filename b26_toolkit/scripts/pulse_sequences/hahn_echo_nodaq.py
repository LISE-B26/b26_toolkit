from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_generic import PulsedExperimentGenericNoDAQ
from b26_toolkit.instruments import B26PulseBlaster, Pulse #, MicrowaveGenerator
from pylabcontrol.core import Parameter
import numpy as np


class HahnEchoNoMeasurement(PulsedExperimentGenericNoDAQ):

    _DEFAULT_SETTINGS = [
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'minimum time between pi pulses'),
            Parameter('max_time', 10000, float, 'maximum time between pi pulses'),
            Parameter('time_step', 5, float,
                      'time step increment of time between pi pulses (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
        ]),
        Parameter('num_averages', 100000, int, 'number of averages')
    ]

    _INSTRUMENTS = {'PB': B26PulseBlaster}
    _SCRIPTS = {}

    def _create_pulse_sequences(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement
        """

        pulse_sequences = []
        tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],self.settings['tau_times']['time_step'])
        tau_list = np.ndarray.tolist(tau_list) # 20180731 ER convert to list

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        # MM: updated to min_pulse_dur
        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time = self.settings['read_out']['nv_reset_time']

        laser_off_time = self.settings['read_out']['laser_off_time']

        for tau in tau_list:

            end_of_first_HE = laser_off_time + tau + tau

            pulse_sequence = [
                 Pulse('laser', end_of_first_HE, nv_reset_time)
                 ]

            start_of_second_HE = end_of_first_HE + nv_reset_time + laser_off_time

            end_of_second_HE = start_of_second_HE + tau + tau

            pulse_sequence += [
                Pulse('laser', end_of_second_HE, nv_reset_time)
            ]

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, 100000