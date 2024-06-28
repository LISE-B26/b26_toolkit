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
from b26_toolkit.instruments import Pulse
from pylabcontrol.core import Parameter
from b26_toolkit.scripts import Esr, FindNvPulsed
from b26_toolkit.tools.utils import get_param_array


class HeraldedCoolingNSweep(PulsedExperimentGeneric):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('tau_mech', 1000, float, 'the time duration of the microwaves (in ns)'),
        ]),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('n_params',[
                      Parameter('n_start', 1, int, 'start frequency of scan in Hz'),
                      Parameter('n_stop', 100, int, 'end frequency of scan in Hz'),
                      Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                                'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
                      Parameter('n_points', 100, int, 'number of frequencies in scan in Hz'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mech_readout', [
            Parameter('pulse_width', 100, int, 'mech trigger pulse width [ns]'),
            Parameter('meas_time', 100000, int, 'wait for mechanics readout [ns]')
        ]),
        Parameter('repetitions', 1, int,
                  'number of repetitions of Pulsed ESR sequence consisting of MW pi-pulse and reinitialization'),
    ]

    _SCRIPTS = {'find_nv': FindNvPulsed, 'esr': Esr}

    def define_sweep_parameters(self):
        """
        Define name of sweep parameters. Since this script is generic, it is coded with variables like 'param_start'.
        This function redirects the code to look for the corresponding settings in _DEFAULT_SETTINGS
        :return:
        """
        self.sweep_params = {'param_start': self.settings['n_start'],
                             'param_stop': self.settings['n_stop'],
                             'param_points': self.settings['n_points'],
                             'param_switching_time': 0}

    def _configure_instruments_start_of_script(self):
        """
        Configure instruments at the beginning of a script, e.g. enable IQ modulation
        :return: None
        """
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

    def _configure_instruments_start_of_sweep(self, param_current):
        pass

    def _configure_param_array(self):
        self.params = get_param_array(self.sweep_params['param_start'],
                                      self.sweep_params['param_stop'],
                                      self.sweep_params['param_points'],
                                      self.settings['range_type'])

    def _create_pulse_sequences(self):
        readout_settings = self.settings['read_out']
        meas_time = readout_settings['meas_time']
        laser_on_time = readout_settings['nv_reset_time']
        laser_off_time = readout_settings['laser_off_time']
        delay_readout = readout_settings['delay_readout']
        delay_mw_readout = readout_settings['delay_mw_readout']

        mw_settings = self.settings['mw_pulses']
        pi_time = mw_settings['pi_pulse_time']
        pi_half_time = mw_settings['pi_half_pulse_time']
        tau_mech = mw_settings['tau_mech']

        mech_settings = self.settings['mech_readout']
        mech_pulse_width = mech_settings['pulse_width']
        mech_wait_time = mech_settings['meas_time']

        pulse_sequences = []
        for n in self.params:
            pulse_sequence = [Pulse('laser', mech_wait_time, laser_on_time),
                              Pulse('apd_readout', mech_wait_time + delay_readout, meas_time),
                              Pulse('microwave_i', mech_wait_time + laser_on_time + laser_off_time, pi_half_time)]

            train_start = mech_wait_time + laser_on_time + laser_off_time + pi_half_time / 2
            for i in range(n):
                pulse_sequence.append(Pulse('microwave_i', train_start + (i + 1) * tau_mech - pi_time / 2, pi_time))


            pulse_sequence.extend([Pulse('microwave_i', train_start + n * tau_mech - pi_half_time / 2, pi_half_time),
                                   Pulse('laser', train_start + n * tau_mech + delay_mw_readout, laser_on_time),
                                   Pulse('apd_readout', train_start + n * tau_mech + delay_mw_readout + delay_readout, meas_time),
                                   Pulse('mech', train_start + n * tau_mech + delay_mw_readout + laser_on_time + laser_off_time, mech_pulse_width)])

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, self.params, meas_time
