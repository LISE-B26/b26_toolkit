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


###############PROBABLY DOESNT WORK#############################
import numpy as np
from pylabcontrol.core import Parameter
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_base_script import PulsedExperimentBaseScript
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, Pulse
from b26_toolkit.plotting.plots_1d import plot_pulses, update_pulse_plot, update_1d_simple

class StroboscopicMechanics(PulsedExperimentBaseScript):

    _DEFAULT_SETTINGS = [
        Parameter('taus_times', [
            Parameter('min_time', 1e4, float, 'minimum strobing frequency [Hz]'),
            Parameter('max_time', 1e6, float, 'maximum strobing frequency [Hz]'),
            Parameter('time_step', 50000, float, 'frequency step [Hz]')
        ]),
        Parameter('read_out', [
            Parameter('laser_time', 10, int, 'time laser is on [ns]'),
            Parameter('meas_time', 10, int, 'measurement time [ns]'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, float, 'acquisition time for each frequency'),
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster}

    def _create_pulse_sequences(self):

        pulseSequences = []

        strobeFreqs = np.arange(self.settings['strobe_freqs']['min_freq'],
                                  self.settings['strobe_freqs']['max_freq'],
                                  self.settings['strobe_freqs']['freq_step'])

        taus = 1.0e9 / strobeFreqs
        measTime = self.settings['read_out']['meas_time']
        laserTime = self.settings['read_out']['laser_time']

        for tau in taus:
            pulseSequences.append([Pulse('laser', tau - laserTime, laserTime), Pulse('apd_readout', tau - measTime, measTime)])

        return pulseSequences, taus, measTime

    def _plot(self, axislist, data=None):

        if data is None:
            data = self.data

        plot_pulses(axislist[1], self.pulse_sequences[self.sequence_index])

        axislist[0].plot(1.0e9 / data['tau'], data['counts'])
        axislist[0].set_xlabel('strobing frequency [Hz]')
        axislist[0].set_ylabel('counts [kcounts/s]')

    def _update_plot(self, axes_list):

        counts = self.data['counts']
        taus = self.data['tau']
        update_1d_simple(axes_list[0], taus, [counts])
        update_pulse_plot(axes_list[1], self.pulse_sequences[self.sequence_index])