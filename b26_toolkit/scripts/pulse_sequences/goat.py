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

# import numpy as np
# from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_base_script import PulsedExperimentBaseScript
# from pylabcontrol.core import Parameter
# from b26_toolkit.scripts import ESR, FindNVPulsed
#
#
# #
# # class Goat(PulsedExperimentBaseScript):
# #     _DEFAULT_SETTINGS = [
# #         Parameter('mw_power', -45.0, float, 'microwave power in dB'),
# #         Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
# #         Parameter('num_averages', 1000000, int, 'number of averages'),
# #         Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
# #         Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
# #         Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
# #                   'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
# #         Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
# #         Parameter('read_out', [
# #             Parameter('meas_time', 250, float, 'measurement time (in ns)'),
# #             Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
# #             Parameter('laser_off_time', 1000, int, 'minimum laser off time before taking measurements (ns)'),
# #             Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
# #             Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
# #         ]),
# #         Parameter('repetitions', 4, int,
# #                   'number of repetitions of Pulsed ESR sequence consisting of MW pi-pulse and reinitialization'),
# #         Parameter('mw_generator_switching_time', .01, float,
# #                   'time wait after switching center frequencies on generator (s)')
# #     ]
# #
# #     _SCRIPTS = {'find_nv': FindNVPulsed, 'esr': ESR}