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
from pylabcontrol.core import Parameter
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_base_script import PulsedExperimentBaseScript
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, Pulse
from b26_toolkit.plotting.plots_1d import plot_pulses, update_pulse_plot, update_1d_simple

class StroboscopicReadout(PulsedExperimentBaseScript):

    _DEFAULT_SETTINGS = [
        Parameter('cycle_period', 10000, int, 'laser and readout off time'),
        Parameter('laser_time', 1000, int, 'time laser is on [ns]'),
        Parameter('meas_time', 1000, int, 'measurement time [ns]'),
        Parameter('delay_readout', 0, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
        Parameter('pulses_per_seq', 10, int, 'number of laser pulses per sequence'),
        Parameter('num_averages', 100000, float, 'acquisition time for each frequency'),
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster}

    def _create_pulse_sequences(self):
        if self.settings['pulses_per_seq'] > self.instruments['PB']['instance'].MAX_COMMANDS / 2:
            raise Exception('pulses_per_seq exceeds pulseblaster maximum of {}'.format(self.instruments['PB']['instance'].MAX_COMMANDS))
        pulseSequences = []

        period = self.settings['cycle_period']
        meas_time = self.settings['meas_time']
        laser_time = self.settings['laser_time']
        delay_readout = self.settings['delay_readout']

        for i in range(self.settings['pulses_per_seq']):
            pulseSequences.extend([Pulse('laser', period * (i + 1) - laser_time, laser_time),
                                Pulse('apd_readout', period * (i + 1) - meas_time + delay_readout, meas_time)])
        # print(pulseSequences)
        return [pulseSequences], [period], meas_time

    # def _plot(self, axislist, data=None):
    #
    #     if data is None:
    #         data = self.data
    #
    #     plot_pulses(axislist[1], self.pulse_sequences[self.sequence_index])
    #
    #     axislist[0].plot(1.0e9 / self.settings['cycle_period'], data['counts'])
    #     axislist[0].set_xlabel('strobing frequency [Hz]')
    #     axislist[0].set_ylabel('counts [kcounts/s]')
    #     print(data['counts'])
    #
    # def _update_plot(self, axes_list):
    #
    #     counts = self.data['counts']
    #     taus = self.data['tau']
    #     update_1d_simple(axes_list[0], taus, [counts])
    #     update_pulse_plot(axes_list[1], self.pulse_sequences[self.sequence_index])
    #     print(counts)