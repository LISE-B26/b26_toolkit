"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
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
import numpy as np
from b26_toolkit.src.scripts.pulse_sequences.pulsed_experiment_base_script import PulsedExperimentBaseScript
from b26_toolkit.src.instruments import NI6259, B26PulseBlaster, MicrowaveGenerator, Pulse
from b26_toolkit.src.plotting.plots_1d import plot_esr, plot_pulses
from pylabcontrol.src.core import Parameter

class PulsedESR(PulsedExperimentBaseScript): # ER 20170616 - wrote for symmetry between 0 and -1 state
    """
This script applies a microwave pulse at fixed power and durations for varying frequencies.
Uses double_init scheme.

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_power']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

        assert self.settings['freq_start'] < self.settings['freq_stop']

        self.data = {'mw_frequencies': np.linspace(self.settings['freq_start'], self.settings['freq_stop'],
                                                   self.settings['freq_points']), 'esr_counts': []}

        for i, mw_frequency in enumerate(self.data['mw_frequencies']):
            self._loop_count = i
            self.instruments['mw_gen']['instance'].update({'frequency': float(mw_frequency)})
            super(PulsedESR, self)._function(self.data)
            self.data['esr_counts'].append(self.data['counts'])

    # def _calc_progress(self):
    #     #COMMENT_ME
    #     # todo: change to _calc_progress(self, index):
    #     progress = int(100. * (self._loop_count) / self.settings['freq_points'])
    #     return progress

    def _plot(self, axes_list, data = None):
        '''
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional) dataset to plot, if not provided use self.data
        '''
        if data is None:
            data = self.data

        mw_frequencies = data['mw_frequencies']
        esr_counts = data['esr_counts']
        axis1 = axes_list[0]
        if not esr_counts == []:
            counts = esr_counts
            plot_esr(axis1, mw_frequencies[0:len(counts)], counts)
            axis1.hold(False)
        axis2 = axes_list[1]
        plot_pulses(axis2, self.pulse_sequences[0])

    def _update_plot(self, axes_list):
        mw_frequencies = self.data['mw_frequencies']
        esr_counts = self.data['esr_counts']
        axis1 = axes_list[0]
        if not esr_counts == []:
            counts = esr_counts
            plot_esr(axis1, mw_frequencies[0:len(counts)], counts)
            axis1.hold(False)
            # axis2 = axes_list[1]
            # update_pulse_plot(axis2, self.pulse_sequences[0])

    def _create_pulse_sequences(self):

        '''

        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''

        tau = self.settings['tau_mw']
        pulse_sequences = []
        tau_list = [tau]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', laser_off_time + tau + 2 * 40, nv_reset_time)]
               #  Pulse('apd_readout', laser_off_time + tau + 2 * 40 + delay_readout, meas_time),
               #  ]
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence += [
                    Pulse(microwave_channel, laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time, tau)]

            pulse_sequence += [
                Pulse('laser',
                      laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2 * 40 + delay_mw_readout,
                      nv_reset_time),
                Pulse('apd_readout',
                      laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2 * 40 + delay_mw_readout + delay_readout,
                      meas_time)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, self.settings['read_out']['meas_time']