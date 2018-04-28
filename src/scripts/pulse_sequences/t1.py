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
from pylabcontrol.src.core import Parameter, Script
from pylabcontrol.src.scripts import SelectPoints
from b26_toolkit.src.data_processing.fit_functions import fit_exp_decay
from b26_toolkit.src.scripts import ESR
from .rabi import Rabi

class T1(PulsedExperimentBaseScript):  # ER 10.21.2017
    """
This script sweeps the readout pulse duration. Uses a double_init scheme
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulse', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_time', 30.0, float, 'pi time in ns')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 15, float, 'minimum time for T1 (in ns)'),
            Parameter('max_time', 200, float, 'total time for T1 (in ns)'),
            Parameter('time_step', 5, [5, 10, 20, 50, 100, 200, 500, 1000, 10000, 50000, 100000, 500000],
                      'time step increment of readout pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('nv_reset_time', 7000, int, 'time with laser on to reset state'),
            Parameter('meas_time', 1000, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        # COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulse']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulse']['mw_frequency']})
        super(T1, self)._function(self.data)


        counts = (self.data['counts'][:, 0]  - self.data['counts'][:,1]) / (self.data['counts'][:, 0] + self.data['counts'][:,1])
        tau = self.data['tau']

        try:
            fits = fit_exp_decay(tau, counts, varibale_phase=True)
            self.data['fits'] = fits
        except:
            self.data['fits'] = None
            self.log('fit failed')

    def _create_pulse_sequences(self):
        '''

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''
        pulse_sequences = []
        # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
        #                  self.settings['tau_times']['time_step'])
        # JG 16-08-25 changed (15ns min spacing is taken care of later):
        tau_list = list(range(int(self.settings['tau_times']['min_time']), int(self.settings['tau_times']['max_time']),
                         self.settings['tau_times']['time_step']))

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= 15]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulse']['microwave_channel']

        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        pi_time = self.settings['mw_pulse']['pi_time']
        meas_time = self.settings['read_out']['meas_time']


        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', laser_off_time + tau, nv_reset_time),
                 Pulse('apd_readout', laser_off_time+ delay_readout + tau, meas_time),
                 ]
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence += [Pulse(microwave_channel, laser_off_time + nv_reset_time + laser_off_time + tau, pi_time)]

            pulse_sequence += [
                Pulse('laser', laser_off_time + nv_reset_time + laser_off_time + delay_mw_readout + 2*tau, nv_reset_time),
                Pulse('apd_readout', laser_off_time + nv_reset_time + laser_off_time + delay_mw_readout + delay_readout + 2*tau, meas_time)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _plot(self, axislist, data=None):
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

        if data['fits'] is not None:
            counts = (data['counts'][:, 0] - data['counts'][:, 1]) /(data['counts'][:, 1] + data['counts'][:, 0])
            tau = data['tau']
            fits = data['fits']  # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            axislist[0].hold(True)

            axislist[0].plot(tau, fit_exp_decay(tau, *fits), 'k', lw=3)
            #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
            axislist[0].set_title('T1 counts')
        else:
            super(T1, self)._plot(axislist)
            axislist[0].set_title('Readout pulse width counts')
            axislist[0].legend(labels=('Ref Fluorescence', 'Pi pulse Data'), fontsize=8)

class T1SingleInit(PulsedExperimentBaseScript): # ER 5.25.2017
    """
This script takes a T1 by measuring the decay of the ms = 0 population,into +/-1. To avoid needing a pi pulse we only do this once on ms = 0
    """
    _DEFAULT_SETTINGS = [
        Parameter('tau_times', [
            Parameter('min_time', 15, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 200, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 5, [5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000, 500000],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None

        super(T1SingleInit, self)._function(self.data)
        counts = self.data['counts']
        tau = self.data['tau']


        try:
            fits = fit_exp_decay(tau, counts, varibale_phase=True)
            self.data['fits'] = fits
        except:
            self.data['fits'] = None
            self.log('exp fit failed')

    def _create_pulse_sequences(self):
        '''

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''
        pulse_sequences = []
        # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
        #                  self.settings['tau_times']['time_step'])
        # JG 16-08-25 changed (15ns min spacing is taken care of later):
        tau_list = list(range(int(self.settings['tau_times']['min_time']), int(self.settings['tau_times']['max_time']),self.settings['tau_times']['time_step']))

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= 15]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', laser_off_time + tau + 2*40, nv_reset_time),
                 Pulse('apd_readout', laser_off_time + tau + 2*40 + delay_readout, meas_time),
                 ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
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

        if data['fits'] is not None:
            counts = data['counts'][:,1]/ data['counts'][:,0]
            tau = data['tau']
            fits = data['fits'] # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            axislist[0].hold(True)

            axislist[0].plot(tau, fit_exp_decay(tau, *fits), 'k', lw=3)

         #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
            axislist[0].set_title('T1 decay of ms = 0 population')
        else:
            super(T1SingleInit, self)._plot(axislist)
            axislist[0].set_title('T1 decay of ms = 0 population')

class T1_double_init_many_NVs(Script):
    _DEFAULT_SETTINGS = [
        Parameter('esr_peak', 'upper', ['upper', 'lower', 'both'], 'if ESR fits two peaks, defines which one to use')
    ]
    _INSTRUMENTS = {}
    _SCRIPTS = {'select_NVs': SelectPoints, 'ESR': ESR, 'Rabi': Rabi, 'T1': T1}

    def _function(self):
        for num, nv_loc in enumerate(self.scripts['select_NVs'].data['nv_locations']):
            if self._abort:
                break
            find_NV_rabi = self.scripts['Rabi'].scripts['find_nv']
            find_NV_rabi.settings['initial_point']['x'] = nv_loc[0]
            find_NV_rabi.settings['initial_point']['y'] = nv_loc[1]
            find_NV_rabi.run()
            self.scripts['ESR'].settings['tag'] = 'esr_NV' + str(num)
            self.scripts['ESR'].run()
            fit_params = self.scripts['ESR'].data['fit_params']
            if fit_params is None:
                continue
            if len(fit_params) == 4:
                freqs = [fit_params[2]]
            elif len(fit_params == 6):
                if self.settings['esr_peak'] == 'lower':
                    freqs = [fit_params[4]]
                elif self.settings['esr_peak'] == 'upper':
                    freqs = [fit_params[5]]
                elif self.settings['esr_peak'] == 'both':
                    freqs = [fit_params[4], fit_params[5]]
            for freq in freqs:
                if self._abort:
                    break
                rabi = self.scripts['Rabi']
                rabi.settings['tag'] = 'rabi_NV' + str(num)
                rabi.settings['mw_pulses']['mw_frequency'] = float(freq)
                rabi.run()
                rabi_fit = rabi.data['fits']
                if rabi_fit is None:
                    continue
                pi_time = abs((np.pi - rabi_fit[2])/rabi_fit[1])
                pi_time = min(max(np.round(pi_time / 2.5) * 2.5, 15.), 300.) #round to nearest 2.5 ns
                find_NV_T1 = self.scripts['T1'].scripts['find_nv']
                find_NV_T1.settings['initial_point']['x'] = find_NV_rabi.data['maximum_point']['x']
                find_NV_T1.settings['initial_point']['y'] = find_NV_rabi.data['maximum_point']['y']
                T1 = self.scripts['T1']
                T1.settings['mw_pulse']['mw_frequency'] = float(freq)
                T1.settings['mw_pulse']['pi_time'] = float(pi_time)
                T1.settings['tag'] = 't1_' + '_NV' + str(num)
                T1.run()

    def plot(self, figure_list):
        if self._current_subscript_stage is not None:
            if self._current_subscript_stage['current_subscript'] is not None:
                self._current_subscript_stage['current_subscript'].plot(figure_list)


    def skip_next(self):
        for script in self.scripts.values():
            script.stop()
