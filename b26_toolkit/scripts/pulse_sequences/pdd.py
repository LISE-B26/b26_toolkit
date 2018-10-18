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

from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_base_script import PulsedExperimentBaseScript
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, MicrowaveGenerator, Pulse
from pylabcontrol.core import Parameter, Script

class PDD(PulsedExperimentBaseScript):
    """
This script runs a PDD ( Periodic Dynamical Decoupling) sequence for different number of pi pulses.
For a single pi-pulse this is a Hahn-echo sequence.
For zero pulses this is a Ramsey sequence.

The sequence is pi/2 - tau/4 - (tau/4 - pi  - tau/4)^n - tau/4 - pi/2

Tau/2 is the time between the center of the pulses!

todo(emma): (make double_init sheme)


    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -2, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('pi_pulse_time', 50, float, 'time duration of pi-pulse (in ns)'),
            Parameter('number_of_pi_pulses', 1, list(range(0, 17)), 'number of pi pulses')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 15, float, 'min value for tau, the free evolution time in between pulses (in ns)'),
            Parameter('max_time', 30, float, 'max value for tau, the free evolution time in between pulses (in ns)'),
            Parameter('time_step', 5, [5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000], 'step size for tau, the free evolution time in between pulses (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after CPMG sequence (in ns)'),
            Parameter('nv_reset_time', 3000, int, 'time duration of the green laser to reset the spin state'),
            Parameter('ref_meas_off_time', 1000, int, 'laser off time before taking reference measurement at the end of init (ns)'),
            Parameter('delay_mw_init', 1000, int, 'delay between initialization and mw (in ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)')
        ]),
        Parameter('num_averages', 1000, int, 'number of averages (should be less than a million)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    _SCRIPTS = {}


    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

    def _function(self):
        #COMMENT_ME
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(PDD, self)._function()


    def _create_pulse_sequences(self):
        '''
        creates the pulse sequence for the Hahn echo /
        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''
        pulse_sequences = []


        tau_list = list(range(int(self.settings['tau_times']['min_time']),
                         int(self.settings['tau_times']['max_time'] + self.settings['tau_times']['time_step']),
                         self.settings['tau_times']['time_step']))



        reset_time = self.settings['read_out']['nv_reset_time']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = pi_time/2.0

        ref_meas_off_time = self.settings['read_out']['ref_meas_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_init = self.settings['read_out']['delay_mw_init']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        number_of_pi_pulses = self.settings['mw_pulses']['number_of_pi_pulses']


        for tau in tau_list:

            # pulse_sequence = [Pulse('laser', 0, reset_time - ref_meas_off_time - 15 - meas_time),
            #                   Pulse('apd_readout', reset_time - 15 - meas_time, meas_time),
            #                   Pulse('laser', reset_time - 15 - meas_time, meas_time),
            #                   Pulse('microwave_i', reset_time + delay_mw_init, pi_half_time)
            #                   ]
            # # 16-08-25 JG: changed :
            pulse_sequence = [Pulse('laser', 0, reset_time - ref_meas_off_time - 15 - meas_time),
                              Pulse('apd_readout', reset_time - 15 - meas_time, meas_time),
                              Pulse('laser', reset_time - 15 - meas_time, meas_time),
                              Pulse('microwave_i', reset_time + delay_mw_init-pi_half_time/2, pi_half_time)
                              ]


            # next_pi_pulse_time = reset_time + delay_mw_init + pi_half_time + tau
            # # 16-08-19 JG: changed :
            next_pi_pulse_time = reset_time + delay_mw_init
            # # 16-08-25 JG: changed :
            # next_pi_pulse_time = reset_time + delay_mw_init - pi_half_time / 2 + tau / 2

            for n in range(1, number_of_pi_pulses + 1):
                next_pi_pulse_time += tau/2
                pulse_sequence.extend([Pulse('microwave_q', next_pi_pulse_time - pi_time/2, pi_time)])
                # next_pi_pulse_time += tau*2 + pi_time
                # 16-08-19 JG: changed:
                # next_pi_pulse_time += tau
                # 16 - 08 -24 JG: changed
                next_pi_pulse_time += tau/2

            if number_of_pi_pulses == 0:
                next_pi_pulse_time += tau

            # pulse_sequence.extend([Pulse('microwave_i', next_pi_pulse_time-tau, pi_half_time),
            #                         Pulse('laser', next_pi_pulse_time-tau + delay_mw_readout + pi_half_time, meas_time),
            #                         Pulse('apd_readout',next_pi_pulse_time-tau + delay_mw_readout + pi_half_time, meas_time)
            #                         ])
            # 16-08-19 JG: changed:
            # pulse_sequence.extend([Pulse('microwave_i', next_pi_pulse_time-tau/2 + pi_half_time, pi_half_time),
            #                         Pulse('laser',      next_pi_pulse_time-tau/2 + pi_time + delay_mw_readout, meas_time),
            #                         Pulse('apd_readout',next_pi_pulse_time-tau/2 + pi_time + delay_mw_readout, meas_time)
            #                         ])
            # pulse_sequences.append(pulse_sequence)
            # 16 - 08 -24 JG: changed
            # pulse_sequence.extend([Pulse('microwave_i', next_pi_pulse_time + pi_half_time, pi_half_time),
            #                        Pulse('laser', next_pi_pulse_time + pi_time + delay_mw_readout, meas_time),
            #                        Pulse('apd_readout', next_pi_pulse_time + pi_time + delay_mw_readout,
            #                              meas_time)
            #                        ])

            # # 16-08-25 JG: changed :
            pulse_sequence.extend([Pulse('microwave_i', next_pi_pulse_time - pi_half_time/2, pi_half_time),
                                   Pulse('laser', next_pi_pulse_time + pi_half_time + delay_mw_readout, meas_time),
                                   Pulse('apd_readout', next_pi_pulse_time + pi_half_time + delay_mw_readout, meas_time)
                                   ])

            pulse_sequences.append(pulse_sequence)

        # TEMPORATTY: THIS IS TO SEE IF THE OVERALL TIME OF A SEQUENCE SHOULD ALWAYS BE THE SAME
        # IF WE WANT TO KEEP THIS ADD ADDITIONAL PARAMETER TO THE SCRIPT SETTINGS
        # end_time_max = 0
        # for pulse_sequence in pulse_sequences:
        #     for pulse in pulse_sequence:
        #         end_time_max = max(end_time_max, pulse.start_time + pulse.duration)
        # for pulse_sequence in pulse_sequences:
        #     pulse_sequence.append(Pulse('laser', end_time_max + 1850, 15))

        return pulse_sequences, tau_list, meas_time
