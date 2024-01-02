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

from b26_toolkit.scripts.pulse_sequences.rabi import Rabi
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, MicrowaveGenerator, Pulse
from pylabcontrol.core import Parameter
from b26_toolkit.data_processing.fit_functions import fit_rabi_decay, cose_with_decay
import numpy as np

MAX_AVERAGES_PER_SCAN = 1000000  # 1E6, the max number of loops per point allowed at one time (true max is ~4E6 since
                                 #pulseblaster stores this value in 22 bits in its register
                                # DS 20191216: changed from 1e5 to 1e6 since loop register is 20 bits. PB will throw error
                                # if too large

MAX_NUM_PULSES = 4096


class Ramsey(Rabi):
    """
    This script executes a Ramsey measurement. It uses a single_init scheme, i.e. there's only one laser pulse, where
    we collect the signal fluor. at the beginning and ref fluor. at the end.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('3pi_half_pulse_time', 75, float, 'time duration of a 3pi/2 pulse (in ns)')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'minimum time between pi pulses'),
            Parameter('max_time', 10000, float, 'maximum time between pi pulses'),
            Parameter('time_step', 5, float, 'time step increment of time between pi pulses (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 1000, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}


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
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        wait_time = 0

        for tau in tau_list:

            pulse_sequence = \
                [
                    Pulse(microwave_channel, laser_off_time + wait_time, pi_half_time),
                    Pulse(microwave_channel, laser_off_time + wait_time + pi_half_time + tau, pi_half_time),
                ]

            end_of_first_ramsey = laser_off_time + wait_time + tau + 2 * pi_half_time

            pulse_sequence += [
                Pulse('laser', end_of_first_ramsey + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', end_of_first_ramsey + delay_mw_readout + delay_readout, meas_time),
            ]

            start_of_second_ramsey = end_of_first_ramsey + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence += \
                [
                    Pulse(microwave_channel, start_of_second_ramsey + wait_time, pi_half_time),
                    Pulse(microwave_channel, start_of_second_ramsey + wait_time + pi_half_time + tau, three_pi_half_time)
                ]

            end_of_second_ramsey = start_of_second_ramsey + wait_time + pi_half_time + tau + three_pi_half_time

            pulse_sequence += [
                Pulse('laser', end_of_second_ramsey + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', end_of_second_ramsey + delay_mw_readout + delay_readout, meas_time)
            ]
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _plot(self, axislist, data = None):
        '''
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional) dataset to plot (dictionary that contains keys counts, tau, fits), if not provided use self.data
        '''

        if data is None:
            data = self.data

        if 1==2 and 'fits' in data.keys() and data['fits'] is not None:
            counts = data['counts'][:,1]/ data['counts'][:,0]
            tau = data['tau']
            fits = data['fits'] # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            #axislist[0].hold(True) ER 20181012

            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            #pi_time = 2*np.pi / fits[1] / 2
            pi_time = (np.pi - fits[2])/fits[1]
            pi_half_time = (np.pi/2 - fits[2])/fits[1]
            three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
            rabi_freq = 1000*fits[1]/(2*np.pi)
         #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            super(Rabi, self)._plot(axislist)
            axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9))
            axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)


