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
import time
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_base_script import PulsedExperimentBaseScript
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, MicrowaveGenerator, Pulse, RFGenerator, MicrowaveGenerator2, AFG3022C, Commander
from b26_toolkit.plotting.plots_1d import plot_pulses, update_pulse_plot, plot_1d_simple_timetrace_ns, update_1d_simple
from pylabcontrol.core import Parameter
from b26_toolkit.scripts import FindNV, ESR
from b26_toolkit.data_processing.fit_functions import fit_rabi_decay, cose_with_decay
import random
from copy import deepcopy

class NuclearRabiSRS(PulsedExperimentBaseScript): # FF_20220725
    """
This script measures Rabi oscillations of a nuclear spin coupled to an NV center. This uses the SRS to drive the RF coil, as opposed to the AFG3022C.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('mw_tau', 300, float, 'microwave pi-pulse duration in ns'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -45.0, float, 'microwave power in dB'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_channel', 'i', ['i', 'q'], 'Channel to use for rf pulses'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns'),
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            #Parameter('time_step', 5., [2.5, 4., 5., 10., 20., 50., 100., 200., 500., 1000., 2000., 5000., 10000., 20000., 100000., 500000.],
            #      'time step increment of rabi pulse duration (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                      Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'], 'use a constant cooldown time or set it based on length of entire sequence'),
                      Parameter('cooldown_time', 5000, float, 'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                                                              'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 4000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 250, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
        Parameter('rf_switch', [
            Parameter('add', True, bool,
                      'check to add mw switch to every i and q pulse and to use switch to carve out pulses. note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 50, int,
                      'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'rf_switch', ['rf_switch', 'rf_iq'],
                      'determines if rf pulses are carved out by mw-switch or by i and q channels of rf source '),
            Parameter('no_iq_overlap', True, bool,
                      'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ]),
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator,
                    'rf_gen': RFGenerator}

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})

        self.instruments['rf_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['rf_gen']['instance'].update({'enable_modulation': True})
        self.instruments['rf_gen']['instance'].update({'amplitude_rf': self.settings['rf_pulses']['rf_power']})
        self.instruments['rf_gen']['instance'].update({'frequency': self.settings['rf_pulses']['rf_frequency']})

        super(NuclearRabiSRS, self)._function()

        if 'counts' in self.data.keys() and 'tau' in self.data.keys():
            counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]

            tau = self.data['tau']

            #try:
            #    fits = fit_rabi_decay(tau, counts, variable_phase=True)
            #    self.data['fits'] = fits
            #except:
            #    self.data['fits'] = None
            #    self.log('rabi fit failed')
            self.data['fits'] = None
            #self.log('rabi fit failed')

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

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_tau = self.settings['mw_pulses']['mw_tau']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
        rf_channel = 'rf_' + self.settings['rf_pulses']['rf_channel']


        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        for tau in tau_list:
            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = tau/duty_cycle - tau
                cooldown_time = int(int(cooldown_time/100)*100)

            pulse_sequence = [Pulse('laser', cooldown_time + laser_off_time, nv_reset_time),
                              #Pulse('apd_readout', cooldown_time + laser_off_time + delay_readout, meas_time),
                              Pulse(microwave_channel, cooldown_time + laser_off_time + nv_reset_time + laser_off_time, mw_tau)
                              ]
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence.append(
                    Pulse(rf_channel, cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle, tau))


            pulse_sequence.append(Pulse(microwave_channel,
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + tau + rf_settle,
                                        mw_tau))
            pulse_sequence.append(Pulse('laser',
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + tau + rf_settle + mw_tau + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + tau + rf_settle + mw_tau + delay_mw_readout + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + tau + rf_settle + mw_tau + delay_mw_readout + nv_reset_time - meas_time,
                                        meas_time))
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

        if 'fits' in data.keys() and data['fits'] is not None and 1 == 2:
            counts = data['counts'][:, 1] / data['counts'][:, 0]
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
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            super(NuclearRabi, self)._plot(axislist)
            axislist[0].set_title('Nuclear Rabi rf_power:{:0.1f}dBm, rf_freq:{:0.3f} MHz'.format(self.settings['rf_pulses']['rf_power'], self.settings['rf_pulses']['rf_frequency']*1e-6))
            axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)


class NuclearRabi(PulsedExperimentBaseScript): # FF_20220725
    """
This script measures Rabi oscillations of a nuclear spin coupled to an NV center.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('mw_tau', 300, float, 'microwave pi-pulse duration in ns'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -45.0, float, 'microwave power in dB'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_channel', 'i', ['i', 'q'], 'Channel to use for rf pulses'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns, for pulses like nuclear polarization')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            #Parameter('time_step', 5., [2.5, 4., 5., 10., 20., 50., 100., 200., 500., 1000., 2000., 5000., 10000., 20000., 100000., 500000.],
            #      'time step increment of rabi pulse duration (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                      Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'], 'use a constant cooldown time or set it based on length of entire sequence'),
                      Parameter('cooldown_time', 5000, float, 'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                                                              'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 4000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 250, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
        Parameter('rf_switch', [
            Parameter('add', True, bool,
                      'check to add mw switch to every i and q pulse and to use switch to carve out pulses. note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 50, int,
                      'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'rf_switch', ['rf_switch', 'rf_iq'],
                      'determines if rf pulses are carved out by mw-switch or by i and q channels of rf source '),
            Parameter('no_iq_overlap', True, bool,
                      'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ]),
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator,
                    'afg': AFG3022C, 'commander': Commander}
    _SCRIPTS = {'find_nv': FindNV, 'esr': ESR}

    def dbm_to_vpp(self, dbm):
        dbm = float(dbm)
        return np.sqrt(10 ** (dbm / 10) / 1000 * 50) * 2 * np.sqrt(2)

    def _configure_instruments_start_of_script(self):
        self.instruments['PB']['instance'].update({'rf_switch': {'status': False}})
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

        self.instruments['afg']['instance'].update({'ch1_invert_polarity': False})
        self.instruments['afg']['instance'].update({'ch2_invert_polarity': False})
        self.instruments['afg']['instance'].update({'ch1_phase': 0})
        self.instruments['afg']['instance'].update({'ch2_phase': 3.1415926})

        for channel in [1, 2]:
            self.instruments['afg']['instance'].update({'ch%i_amplitude' % channel: float(self.dbm_to_vpp(self.settings['rf_pulses']['rf_power']))})
            self.instruments['afg']['instance'].update({'ch%i_frequency' % channel: float(self.settings['rf_pulses']['rf_frequency'])})
            self.instruments['afg']['instance'].update({'ch%i_function' % channel: 'Sine'})
            #self.instruments['afg']['instance'].update({'ch%i_phase' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_offset' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_enable' % channel: True})

        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None
        self._configure_instruments_start_of_script()


        super(NuclearRabi, self)._function()

        if 'counts' in self.data.keys() and 'tau' in self.data.keys():
            counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]

            tau = self.data['tau']

            #try:
            #    fits = fit_rabi_decay(tau, counts, variable_phase=True)
            #    self.data['fits'] = fits
            #except:
            #    self.data['fits'] = None
            #    self.log('rabi fit failed')
            self.data['fits'] = None
            #self.log('rabi fit failed')

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

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_tau = self.settings['mw_pulses']['mw_tau']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
        rf_channel = 'rf_' + self.settings['rf_pulses']['rf_channel']


        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        for tau in tau_list:
            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = tau/duty_cycle - tau
                cooldown_time = int(int(cooldown_time/100)*100)

            pulse_sequence = [Pulse('laser', cooldown_time + laser_off_time, nv_reset_time),
                              #Pulse('apd_readout', cooldown_time + laser_off_time + delay_readout, meas_time),
                              Pulse(microwave_channel, cooldown_time + laser_off_time + nv_reset_time + laser_off_time, mw_tau)
                              ]
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence.append(
                    Pulse(rf_channel, cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle, tau))


            pulse_sequence.append(Pulse(microwave_channel,
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + tau + rf_settle,
                                        mw_tau))
            pulse_sequence.append(Pulse('laser',
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + tau + rf_settle + mw_tau + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + tau + rf_settle + mw_tau + delay_mw_readout + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + tau + rf_settle + mw_tau + delay_mw_readout + nv_reset_time - meas_time,
                                        meas_time))
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

        super(NuclearRabi, self)._plot(axislist)
        axislist[0].set_title('Nuclear Rabi rf_power:{:0.1f}dBm, rf_freq:{:0.3f} MHz'.format(self.settings['rf_pulses']['rf_power'], self.settings['rf_pulses']['rf_frequency']*1e-6))
        axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)


class RabiPolarized(NuclearRabi):

    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('microwave_channel_1', 'i', ['i', 'q'], 'Channel on mw_gen_1 to use for mw pulses'),
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float, 'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('microwave_channel_2', 'i', ['i', 'q'], 'Channel on mw_gen_2 to use for mw pulses'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float, 'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            # Parameter('time_step', 5., [2.5, 4., 5., 10., 20., 50., 100., 200., 500., 1000., 2000., 5000., 10000., 20000., 100000., 500000.],
            #      'time step increment of rabi pulse duration (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'],
                          'use a constant cooldown time or set it based on length of entire sequence'),
                Parameter('cooldown_time', 1, float,
                          'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                          'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int,
                      'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            # Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int,
                      'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('rf_switch', [
            Parameter('add', True, bool,
                      'check to add mw switch to every i and q pulse and to use switch to carve out pulses. note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 50, int,
                      'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'rf_switch', ['rf_switch', 'rf_iq'],
                      'determines if rf pulses are carved out by mw-switch or by i and q channels of rf source '),
            Parameter('no_iq_overlap', True, bool,
                      'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ])
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator,
                    'afg': AFG3022C, 'mw_gen_2': MicrowaveGenerator2}
    _SCRIPTS = {'find_nv': FindNV, 'esr': ESR}

    def _configure_instruments_start_of_script(self):
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_1_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_1_frequency']})

        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_2_power']})
        self.instruments['mw_gen_2']['instance'].update({'frequency': self.settings['mw_pulses']['mw_2_frequency']})

        self.instruments['afg']['instance'].update({'ch1_invert_polarity': False})
        self.instruments['afg']['instance'].update({'ch2_invert_polarity': False})

        for channel in [1, 2]:
            self.instruments['afg']['instance'].update({'ch%i_amplitude' % channel: float(self.dbm_to_vpp(self.settings['rf_pulses']['rf_power']))})
            self.instruments['afg']['instance'].update({'ch%i_frequency' % channel: float(self.settings['rf_pulses']['rf_frequency'])})
            self.instruments['afg']['instance'].update({'ch%i_function' % channel: 'Sine'})
            #self.instruments['afg']['instance'].update({'ch%i_phase' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_offset' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_enable' % channel: True})

        self.instruments['afg']['instance'].align_phase()
        self.instruments['afg']['instance'].update({'ch1_phase': 0})
        self.instruments['afg']['instance'].update({'ch2_phase': 3.1415926})

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

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])


        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']

        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel_1']
        microwave_channel_2 = 'microwave_' + self.settings['mw_pulses']['microwave_channel_2']+'_2'
        #microwave_channel = 'microwave_i'
        #microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        pulse_sequences = []

        for tau in tau_list:
            t_rotation = tau

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = rf_pi_time/duty_cycle - rf_pi_time
                cooldown_time = int(int(cooldown_time/100)*100)

            pulse_sequence = []
            start_time = 0
            for part in range(2):
                # Initialize nuclear and electronic spin at the beginning of each sequence
                pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
                current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
                for iteration in range(self.settings['polarization_iterations']):
                    pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                    if rf_pi_time > 0:
                        pulse_sequence.append(Pulse(rf_channel, current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                        current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                    for i in range(self.settings['initialization_laser']['n']):
                        pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                        current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

                if part == 0:
                    pulse_sequence.append(Pulse(microwave_channel_2, current_time + laser_off_time, t_rotation))

                elif part == 1:
                    pulse_sequence.append(Pulse(microwave_channel, current_time + laser_off_time, t_rotation))

                current_time_2 = current_time + laser_off_time + t_rotation + delay_mw_readout

                pulse_sequence.append(Pulse('laser', current_time_2,
                                            nv_reset_time_long))
                pulse_sequence.append(Pulse('apd_readout',
                                            current_time_2 + delay_readout,
                                            meas_time))

                if part == 1:
                    pulse_sequence.append(Pulse('apd_readout',
                                                current_time_2 + nv_reset_time_long - meas_time,
                                                meas_time))

                # Update start time for the part 2 of the sequence (where we look at the other nuclear spin population)
                start_time = current_time_2 + nv_reset_time_long

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

        super(NuclearRabi, self)._plot(axislist)
        axislist[0].set_title('Nuclear Rabi rf_power:{:0.1f}dBm, rf_freq:{:0.3f} MHz'.format(self.settings['rf_pulses']['rf_power'], self.settings['rf_pulses']['rf_frequency']*1e-6))
        axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)



class TransportDebugSingleSWAP(RabiPolarized):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'],
                          'use a constant cooldown time or set it based on length of entire sequence'),
                Parameter('cooldown_time', 1, float,
                          'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                          'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('t_movement', 20000, float, 'momvement/wait time in ns before reading state of the nuclear spin'),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int,
                      'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
            # Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int,
                      'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 250, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 10, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 500, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 100, float, 'short wait time between each readout laser pulse')
            ])
        ]),
        Parameter('rf_switch', [
            Parameter('add', True, bool,
                      'check to add mw switch to every i and q pulse and to use switch to carve out pulses. note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 50, int,
                      'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'rf_switch', ['rf_switch', 'rf_iq'],
                      'determines if rf pulses are carved out by mw-switch or by i and q channels of rf source '),
            Parameter('no_iq_overlap', True, bool,
                      'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ])
    ]

    def _create_pulse_sequences_garbage(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        self.num_ref_reads = 3
        pulse_sequences = []
        tau_transport = 0 # Remove this when we do actual transport stuff!

        for tau in tau_list:

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = rf_pi_time / duty_cycle - rf_pi_time
                cooldown_time = int(int(cooldown_time / 100) * 100)

            pulse_sequence = []
            start_time = 0

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
            for i in range(self.num_ref_reads):
                pulse_sequence.append(
                    Pulse('apd_readout', start_time + laser_off_time + cooldown_time + nv_reset_time_long - (i + 1) * (
                                meas_time + delay_readout),
                          meas_time))
            current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse(rf_channel, current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

            # Finished nuclear polarization, now onto the main sequence
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + laser_off_time,
                                        mw_2_pi_time))
            pulse_sequence.append(Pulse(rf_channel,
                                        current_time + laser_off_time + mw_2_pi_time + rf_settle,
                                        rf_pi_half_time))
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + laser_off_time + mw_2_pi_time + rf_settle + rf_pi_half_time +
                                        rf_settle,
                                        mw_2_pi_time))
            # Wait for tau to accumulate phase (Ramsey)
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + laser_off_time + mw_2_pi_time + rf_settle + rf_pi_half_time +
                                        rf_settle + mw_2_pi_time + tau,
                                        mw_2_pi_time))
            # Phase is now stored in nuclear spin, nuclear and e-spin now disentangled
            # Wait for time t_transport, but for now we skip the wait time
            pulse_sequence.append(Pulse(rf_channel,
                                        current_time + laser_off_time + mw_2_pi_time + rf_settle + rf_pi_half_time +
                                        rf_settle + mw_2_pi_time + tau + mw_2_pi_time + tau_transport + rf_settle,
                                        rf_pi_half_time))



            current_time = current_time + laser_off_time + mw_2_pi_time + rf_settle + rf_pi_half_time + rf_settle + \
                           mw_2_pi_time + tau + mw_2_pi_time + tau_transport + rf_settle + rf_pi_half_time + rf_settle
            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                if m == 0:
                    pulse_sequence.append(Pulse(microwave_channel_2,
                                                current_time,
                                                mw_2_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                else:
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time,
                                                mw_1_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

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

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        self.num_ref_reads = 3
        pulse_sequences = []
        t_movement = self.settings['t_movement'] # Remove this when we do actual transport stuff!

        for tau in tau_list:

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = rf_pi_time / duty_cycle - rf_pi_time
                cooldown_time = int(int(cooldown_time / 100) * 100)

            pulse_sequence = []
            start_time = 0

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
            for i in range(self.num_ref_reads):
                pulse_sequence.append(
                    Pulse('apd_readout', start_time + laser_off_time + cooldown_time + nv_reset_time_long - (i + 1) * (
                                meas_time + delay_readout),
                          meas_time))
            current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse(rf_channel, current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

            # Finished nuclear polarization, now onto the main sequence
            if tau > 0:
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time + laser_off_time,
                                            tau))
            pulse_sequence.append(Pulse(rf_channel,
                                        current_time + laser_off_time + tau + delay_mw_readout,
                                        rf_pi_time))
            pulse_sequence.append(Pulse(microwave_channel,
                                        current_time + laser_off_time + tau + delay_mw_readout + rf_pi_time + rf_settle,
                                        mw_1_pi_time))

            # Phase is now stored in nuclear spin, nuclear and e-spin now disentangled
            # Wait for time t_movement
            current_time = current_time + laser_off_time + tau + delay_mw_readout + rf_pi_time + rf_settle + mw_1_pi_time + delay_mw_readout + t_movement

            # Read out state of the nuclear spin
            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                pulse_sequence.append(Pulse(microwave_channel,
                                            current_time,
                                            mw_1_pi_time))
                pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                            nv_reset_time_short))
                pulse_sequence.append(
                    Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                          meas_time))
                current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

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

        if 'fits' in data.keys() and data['fits'] is not None and 1 == 2:
            counts = data['counts'][:, 1] / data['counts'][:, 0]
            tau = data['tau']
            fits = data['fits']  # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            # axislist[0].hold(True) ER 20181012

            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            # pi_time = 2*np.pi / fits[1] / 2
            pi_time = (np.pi - fits[2]) / fits[1]
            pi_half_time = (np.pi / 2 - fits[2]) / fits[1]
            three_pi_half_time = (3 * np.pi / 2 - fits[2]) / fits[1]
            rabi_freq = 1000 * fits[1] / (2 * np.pi)
            axislist[0].set_title(
                'Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(
                    rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            if 'counts' in data.keys():
                avg_counts = np.transpose(np.array([np.average(self.data['counts'][:, self.num_ref_reads:], axis=1)]))
                ref_counts = np.transpose(np.array([np.average(self.data['counts'][:, 0:self.num_ref_reads], axis=1)]))
                first_counts = np.transpose(
                    np.array([np.average(self.data['counts'][:, self.num_ref_reads:self.num_ref_reads + 1], axis=1)]))
                # ref_counts = self.data['counts'][:, 0:1]
                # first_counts = self.data['counts'][:, 1:2]

                # plot_1d_simple_timetrace_ns(axislist[0], data['tau'], [data['counts']])
                plot_1d_simple_timetrace_ns(axislist[0], data['tau'], [ref_counts, first_counts, avg_counts])
                plot_pulses(axislist[1], self.pulse_sequences[self.sequence_index])
            axislist[0].set_title(
                'Nuclear Rabi rf_power:{:0.1f}dBm, rf_freq:{:0.3f} MHz'.format(self.settings['rf_pulses']['rf_power'],
                                                                               self.settings['rf_pulses'][
                                                                                   'rf_frequency'] * 1e-6))
            axislist[0].legend(labels=('Ref Fluorescence', 'First Readout', 'Avg Readout'), fontsize=8)

    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''
        #        if self.scripts['find_nv'].is_running:
        #            self.scripts['find_nv']._update_plot(axes_list)
        #        else:

        x_data = self.data['tau']
        # avg_counts = np.transpose(np.array([np.average(self.data['counts'][:, 1:], axis=1)]))
        # ref_counts = self.data['counts'][:, 0:1]
        # first_counts = self.data['counts'][:, 1:2]

        avg_counts = np.transpose(np.array([np.average(self.data['counts'][:, self.num_ref_reads:], axis=1)]))
        ref_counts = np.transpose(np.array([np.average(self.data['counts'][:, 0:self.num_ref_reads], axis=1)]))
        first_counts = np.transpose(
            np.array([np.average(self.data['counts'][:, self.num_ref_reads:self.num_ref_reads + 1], axis=1)]))

        axis1 = axes_list[0]
        if not self.data['counts'] == []:
            update_1d_simple(axis1, x_data, [ref_counts, first_counts, avg_counts])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])


class TransportDebugDoubleSWAP(RabiPolarized):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_1_pi_half_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_pi_half_time', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'],
                          'use a constant cooldown time or set it based on length of entire sequence'),
                Parameter('cooldown_time', 1, float,
                          'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                          'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('t_movement', 20000, float, 'momvement/wait time in ns before reading state of the nuclear spin'),
        Parameter('storage_axis', 'z', ['x', 'y', 'z'],
                  'rotate spin onto this axis before storage onto nuclear spin; for z, there is no rotation'),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int,
                      'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
            # Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int,
                      'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 250, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 10, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 500, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 100, float, 'short wait time between each readout laser pulse')
            ])
        ]),
        Parameter('rf_switch', [
            Parameter('add', True, bool,
                      'check to add mw switch to every i and q pulse and to use switch to carve out pulses. note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 50, int,
                      'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'rf_switch', ['rf_switch', 'rf_iq'],
                      'determines if rf pulses are carved out by mw-switch or by i and q channels of rf source '),
            Parameter('no_iq_overlap', True, bool,
                      'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ])
    ]

    def _create_pulse_sequences_garbage(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        self.num_ref_reads = 3
        pulse_sequences = []
        tau_transport = 0 # Remove this when we do actual transport stuff!

        for tau in tau_list:

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = rf_pi_time / duty_cycle - rf_pi_time
                cooldown_time = int(int(cooldown_time / 100) * 100)

            pulse_sequence = []
            start_time = 0

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
            for i in range(self.num_ref_reads):
                pulse_sequence.append(
                    Pulse('apd_readout', start_time + laser_off_time + cooldown_time + nv_reset_time_long - (i + 1) * (
                                meas_time + delay_readout),
                          meas_time))
            current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse(rf_channel, current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

            # Finished nuclear polarization, now onto the main sequence
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + laser_off_time,
                                        mw_2_pi_time))
            pulse_sequence.append(Pulse(rf_channel,
                                        current_time + laser_off_time + mw_2_pi_time + rf_settle,
                                        rf_pi_half_time))
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + laser_off_time + mw_2_pi_time + rf_settle + rf_pi_half_time +
                                        rf_settle,
                                        mw_2_pi_time))
            # Wait for tau to accumulate phase (Ramsey)
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + laser_off_time + mw_2_pi_time + rf_settle + rf_pi_half_time +
                                        rf_settle + mw_2_pi_time + tau,
                                        mw_2_pi_time))
            # Phase is now stored in nuclear spin, nuclear and e-spin now disentangled
            # Wait for time t_transport, but for now we skip the wait time
            pulse_sequence.append(Pulse(rf_channel,
                                        current_time + laser_off_time + mw_2_pi_time + rf_settle + rf_pi_half_time +
                                        rf_settle + mw_2_pi_time + tau + mw_2_pi_time + tau_transport + rf_settle,
                                        rf_pi_half_time))



            current_time = current_time + laser_off_time + mw_2_pi_time + rf_settle + rf_pi_half_time + rf_settle + \
                           mw_2_pi_time + tau + mw_2_pi_time + tau_transport + rf_settle + rf_pi_half_time + rf_settle
            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                if m == 0:
                    pulse_sequence.append(Pulse(microwave_channel_2,
                                                current_time,
                                                mw_2_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                else:
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time,
                                                mw_1_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _create_pulse_sequences_simple_rf(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_1_pi_half_time = self.settings['mw_pulses']['mw_1_pi_half_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']
        mw_2_pi_half_time = self.settings['mw_pulses']['mw_2_pi_half_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        self.num_ref_reads = 3
        pulse_sequences = []
        t_movement = self.settings['t_movement'] # Remove this when we do actual transport stuff!

        for tau in tau_list:

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = rf_pi_time / duty_cycle - rf_pi_time
                cooldown_time = int(int(cooldown_time / 100) * 100)

            pulse_sequence = []
            start_time = 0

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
            """
            for i in range(self.num_ref_reads):
                pulse_sequence.append(
                    Pulse('apd_readout', start_time + laser_off_time + cooldown_time + nv_reset_time_long - (i + 1) * (
                                meas_time + delay_readout),
                          meas_time))
          """
            current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse(rf_channel, current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time += nv_reset_time_pulsed + nv_reset_time_wait
            current_time += laser_off_time-nv_reset_time_wait
            # Finished nuclear polarization, now onto the main sequence

            if self.settings['storage_axis'] == 'z':
                pass
            elif self.settings['storage_axis'] == 'x':
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time,
                                            mw_2_pi_half_time))
                #pulse_sequence.append(Pulse(microwave_channel_2,
                #                            current_time + mw_2_pi_half_time + int(rf_pi_time/4) - int(mw_2_pi_time/2),
                #                            mw_2_pi_time))
                current_time += mw_2_pi_half_time + 0*delay_mw_readout

            elif self.settings['storage_axis'] == 'y':
                raise NotImplementedError

            pulse_sequence.append(Pulse(rf_channel,
                                        current_time,
                                        rf_pi_time))


            pulse_sequence.append(Pulse(microwave_channel,
                                        current_time + rf_pi_time + rf_settle,
                                        mw_1_pi_time))

            # Phase is now stored in nuclear spin, nuclear and e-spin now disentangled
            # Wait for time t_movement
            current_time += rf_pi_time + rf_settle + mw_1_pi_time + delay_mw_readout + t_movement

            # Retrieve from memory
            pulse_sequence.append(Pulse(microwave_channel,
                                        current_time,
                                        mw_1_pi_time))
            pulse_sequence.append(Pulse(rf_channel,
                                        current_time + mw_2_pi_time + delay_mw_readout,
                                        rf_pi_time))
            if tau > 0:
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time + mw_2_pi_time + delay_mw_readout + rf_pi_time + rf_settle,
                                            tau))
            current_time += mw_1_pi_time + delay_mw_readout + rf_pi_time + rf_settle + tau + delay_mw_readout

            # Readout
            pulse_sequence.append(Pulse('laser',
                                        current_time,
                                        nv_reset_time_long))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + nv_reset_time_long - meas_time,
                                        meas_time))
            # Read out state of the nuclear spin
            """
            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                pulse_sequence.append(Pulse(microwave_channel,
                                            current_time,
                                            mw_1_pi_time))
                pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                            nv_reset_time_short))
                pulse_sequence.append(
                    Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                          meas_time))
                current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
            """

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

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

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_1_pi_half_time = self.settings['mw_pulses']['mw_1_pi_half_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']
        mw_2_pi_half_time = self.settings['mw_pulses']['mw_2_pi_half_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        self.num_ref_reads = 3
        pulse_sequences = []
        t_movement = self.settings['t_movement']  # Remove this when we do actual transport stuff!

        for tau in tau_list:

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = rf_pi_time / duty_cycle - rf_pi_time
                cooldown_time = int(int(cooldown_time / 100) * 100)

            pulse_sequence = []
            start_time = 0

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
            """
            for i in range(self.num_ref_reads):
                pulse_sequence.append(
                    Pulse('apd_readout', start_time + laser_off_time + cooldown_time + nv_reset_time_long - (i + 1) * (
                                meas_time + delay_readout),
                          meas_time))
          """
            current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse(rf_channel, current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time += nv_reset_time_pulsed + nv_reset_time_wait
            current_time += laser_off_time - nv_reset_time_wait
            # Finished nuclear polarization, now onto the main sequence

            if self.settings['storage_axis'] == 'z':
                pass
            elif self.settings['storage_axis'] == 'x':
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time,
                                            mw_2_pi_half_time))
                # pulse_sequence.append(Pulse(microwave_channel_2,
                #                            current_time + mw_2_pi_half_time + int(rf_pi_time/4) - int(mw_2_pi_time/2),
                #                            mw_2_pi_time))
                current_time += mw_2_pi_half_time + 0 * delay_mw_readout

            elif self.settings['storage_axis'] == 'y':
                raise NotImplementedError

            pulse_sequence.append(Pulse(rf_channel,
                                        current_time,
                                        rf_pi_half_time))
            pulse_sequence.append(Pulse(microwave_channel,
                                        current_time + rf_pi_half_time + rf_settle,
                                        mw_1_pi_time))
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + rf_pi_half_time + rf_settle,
                                        mw_2_pi_time))
            pulse_sequence.append(Pulse(microwave_channel,
                                        current_time + rf_pi_half_time + rf_settle + mw_1_pi_time + rf_pi_half_time*2,
                                        mw_1_pi_time))
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + rf_pi_half_time + rf_settle + mw_2_pi_time + rf_pi_half_time*2,
                                        mw_2_pi_time))
            pulse_sequence.append(Pulse(rf_channel,
                                        current_time + rf_pi_half_time + rf_settle + mw_2_pi_time + rf_pi_half_time * 2 + mw_2_pi_time,
                                        rf_pi_half_time))


            pulse_sequence.append(Pulse(microwave_channel,
                                        current_time + rf_pi_half_time + rf_settle + mw_2_pi_time + rf_pi_half_time * 2 + mw_2_pi_time + rf_pi_half_time + rf_settle,
                                        mw_1_pi_time))

            # Phase is now stored in nuclear spin, nuclear and e-spin now disentangled
            # Wait for time t_movement
            current_time += rf_pi_half_time + rf_settle + mw_2_pi_time + rf_pi_half_time * 2 + mw_2_pi_time + rf_pi_half_time + rf_settle + mw_1_pi_time + delay_mw_readout + t_movement

            # Retrieve from memory
            pulse_sequence.append(Pulse(microwave_channel,
                                        current_time,
                                        mw_1_pi_time))

            pulse_sequence.append(Pulse(rf_channel,
                                        current_time + mw_1_pi_time,
                                        rf_pi_half_time))
            pulse_sequence.append(Pulse(microwave_channel,
                                        current_time + mw_1_pi_time + rf_pi_half_time + rf_settle,
                                        mw_1_pi_time))
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + mw_1_pi_time + rf_pi_half_time + rf_settle,
                                        mw_2_pi_time))
            pulse_sequence.append(Pulse(microwave_channel,
                                        current_time + mw_1_pi_time + rf_pi_half_time + rf_settle + mw_1_pi_time + rf_pi_half_time * 2,
                                        mw_1_pi_time))
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + mw_1_pi_time + rf_pi_half_time + rf_settle + mw_2_pi_time + rf_pi_half_time * 2,
                                        mw_2_pi_time))
            pulse_sequence.append(Pulse(rf_channel,
                                        current_time + mw_1_pi_time + rf_pi_half_time + rf_settle + mw_2_pi_time + rf_pi_half_time * 2 + mw_2_pi_time,
                                        rf_pi_half_time))

            if tau > 0:
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time + mw_1_pi_time + rf_pi_half_time + rf_settle + mw_2_pi_time + rf_pi_half_time * 2 + mw_2_pi_time + rf_pi_half_time + rf_settle,
                                            tau))
            current_time += mw_1_pi_time + rf_pi_half_time + rf_settle + mw_2_pi_time + rf_pi_half_time * 2 + mw_2_pi_time + rf_pi_half_time + rf_settle + tau + delay_mw_readout

            # Readout
            pulse_sequence.append(Pulse('laser',
                                        current_time,
                                        nv_reset_time_long))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + nv_reset_time_long - meas_time,
                                        meas_time))
            # Read out state of the nuclear spin
            """
            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                pulse_sequence.append(Pulse(microwave_channel,
                                            current_time,
                                            mw_1_pi_time))
                pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                            nv_reset_time_short))
                pulse_sequence.append(
                    Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                          meas_time))
                current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
            """

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time


class TransportDebugDoubleSWAPRepetitiveReadout(TransportDebugSingleSWAP):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_1_pi_half_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_pi_half_time', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'],
                          'use a constant cooldown time or set it based on length of entire sequence'),
                Parameter('cooldown_time', 1, float,
                          'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                          'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('t_movement', 20000, float, 'momvement/wait time in ns before reading state of the nuclear spin'),
        Parameter('storage_axis', 'z', ['x', 'y', 'z'],
                  'rotate spin onto this axis before storage onto nuclear spin; for z, there is no rotation'),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int,
                      'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
            # Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int,
                      'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 250, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 10, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 500, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 100, float, 'short wait time between each readout laser pulse')
            ])
        ]),
        Parameter('rf_switch', [
            Parameter('add', True, bool,
                      'check to add mw switch to every i and q pulse and to use switch to carve out pulses. note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 50, int,
                      'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'rf_switch', ['rf_switch', 'rf_iq'],
                      'determines if rf pulses are carved out by mw-switch or by i and q channels of rf source '),
            Parameter('no_iq_overlap', True, bool,
                      'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ])
    ]

    def _create_pulse_sequences_garbage(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        self.num_ref_reads = 3
        pulse_sequences = []
        tau_transport = 0 # Remove this when we do actual transport stuff!

        for tau in tau_list:

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = rf_pi_time / duty_cycle - rf_pi_time
                cooldown_time = int(int(cooldown_time / 100) * 100)

            pulse_sequence = []
            start_time = 0

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
            for i in range(self.num_ref_reads):
                pulse_sequence.append(
                    Pulse('apd_readout', start_time + laser_off_time + cooldown_time + nv_reset_time_long - (i + 1) * (
                                meas_time + delay_readout),
                          meas_time))
            current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse(rf_channel, current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

            # Finished nuclear polarization, now onto the main sequence
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + laser_off_time,
                                        mw_2_pi_time))
            pulse_sequence.append(Pulse(rf_channel,
                                        current_time + laser_off_time + mw_2_pi_time + rf_settle,
                                        rf_pi_half_time))
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + laser_off_time + mw_2_pi_time + rf_settle + rf_pi_half_time +
                                        rf_settle,
                                        mw_2_pi_time))
            # Wait for tau to accumulate phase (Ramsey)
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + laser_off_time + mw_2_pi_time + rf_settle + rf_pi_half_time +
                                        rf_settle + mw_2_pi_time + tau,
                                        mw_2_pi_time))
            # Phase is now stored in nuclear spin, nuclear and e-spin now disentangled
            # Wait for time t_transport, but for now we skip the wait time
            pulse_sequence.append(Pulse(rf_channel,
                                        current_time + laser_off_time + mw_2_pi_time + rf_settle + rf_pi_half_time +
                                        rf_settle + mw_2_pi_time + tau + mw_2_pi_time + tau_transport + rf_settle,
                                        rf_pi_half_time))



            current_time = current_time + laser_off_time + mw_2_pi_time + rf_settle + rf_pi_half_time + rf_settle + \
                           mw_2_pi_time + tau + mw_2_pi_time + tau_transport + rf_settle + rf_pi_half_time + rf_settle
            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                if m == 0:
                    pulse_sequence.append(Pulse(microwave_channel_2,
                                                current_time,
                                                mw_2_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                else:
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time,
                                                mw_1_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

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

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_1_pi_half_time = self.settings['mw_pulses']['mw_1_pi_half_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']
        mw_2_pi_half_time = self.settings['mw_pulses']['mw_2_pi_half_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        self.num_ref_reads = 3
        pulse_sequences = []
        t_movement = self.settings['t_movement'] # Remove this when we do actual transport stuff!

        for tau in tau_list:

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = rf_pi_time / duty_cycle - rf_pi_time
                cooldown_time = int(int(cooldown_time / 100) * 100)

            pulse_sequence = []
            start_time = 0

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))

            for i in range(self.num_ref_reads):
                pulse_sequence.append(
                    Pulse('apd_readout', start_time + laser_off_time + cooldown_time + nv_reset_time_long - (i + 1) * (
                                meas_time + delay_readout),
                          meas_time))

            current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse(rf_channel, current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time += nv_reset_time_pulsed + nv_reset_time_wait
            current_time += laser_off_time-nv_reset_time_wait
            # Finished nuclear polarization, now onto the main sequence

            if self.settings['storage_axis'] == 'z':
                pass
            elif self.settings['storage_axis'] == 'x':
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time,
                                            mw_2_pi_half_time))
                current_time += mw_2_pi_half_time + delay_mw_readout
            elif self.settings['storage_axis'] == 'y':
                raise NotImplementedError

            pulse_sequence.append(Pulse(rf_channel,
                                        current_time,
                                        rf_pi_time))
            pulse_sequence.append(Pulse(microwave_channel,
                                        current_time + rf_pi_time + rf_settle,
                                        mw_1_pi_time))

            # Phase is now stored in nuclear spin, nuclear and e-spin now disentangled
            # Wait for time t_movement
            current_time += rf_pi_time + rf_settle + mw_1_pi_time + delay_mw_readout + t_movement

            # Retrieve from memory
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time,
                                        mw_2_pi_time))
            pulse_sequence.append(Pulse(rf_channel,
                                        current_time + mw_2_pi_time + delay_mw_readout,
                                        rf_pi_time))
            if tau > 0:
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time + mw_2_pi_time + delay_mw_readout + rf_pi_time + rf_settle,
                                            tau))
            pulse_sequence.append(Pulse(rf_channel,
                                        current_time + mw_2_pi_time + delay_mw_readout + rf_pi_time + rf_settle + tau + delay_mw_readout,
                                        rf_pi_time))
            current_time += mw_2_pi_time + delay_mw_readout + rf_pi_time + rf_settle + tau + delay_mw_readout + rf_pi_time + rf_settle

            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                if m == 0:
                    pulse_sequence.append(Pulse('laser', current_time,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + delay_readout,
                              meas_time))
                    current_time = current_time + nv_reset_time_short + laser_off_time_short
                else:
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time,
                                                mw_1_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time


class TransportDebug02(RabiPolarized):
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

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])


        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        pulse_sequences = []

        for tau in tau_list:
            t_rotation = tau

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = rf_pi_time/duty_cycle - rf_pi_time
                cooldown_time = int(int(cooldown_time/100)*100)

            pulse_sequence = []
            start_time = 0
            for part in range(2):
                # Initialize nuclear and electronic spin at the beginning of each sequence
                pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
                current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
                for iteration in range(self.settings['polarization_iterations']):
                    pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                    if rf_pi_time > 0:
                        pulse_sequence.append(Pulse(rf_channel, current_time + mw_1_pi_time, rf_pi_time))
                        current_time = current_time + mw_1_pi_time + rf_pi_time + rf_settle
                    for i in range(self.settings['initialization_laser']['n']):
                        pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                        current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

                pulse_sequence.append(Pulse(microwave_channel_2, current_time + laser_off_time, mw_2_pi_time))
                if rf_pi_half_time > 0:
                    pulse_sequence.append(Pulse(rf_channel, current_time + laser_off_time + mw_2_pi_time, rf_pi_half_time))
                if part == 0:
                    pulse_sequence.append(Pulse(microwave_channel_2,
                                                current_time + laser_off_time + mw_2_pi_time + rf_pi_half_time + rf_settle,
                                                t_rotation))

                elif part == 1:
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time + laser_off_time + mw_2_pi_time + rf_pi_half_time + rf_settle,
                                                t_rotation))

                current_time_2 = current_time + laser_off_time + mw_2_pi_time + rf_pi_half_time + rf_settle + t_rotation + delay_mw_readout

                pulse_sequence.append(Pulse('laser', current_time_2,
                                            nv_reset_time_long))
                pulse_sequence.append(Pulse('apd_readout',
                                            current_time_2 + delay_readout,
                                            meas_time))

                if part == 1:
                    pulse_sequence.append(Pulse('apd_readout',
                                                current_time_2 + nv_reset_time_long - meas_time,
                                                meas_time))

                # Update start time for the part 2 of the sequence (where we look at the other nuclear spin population)
                start_time = current_time_2 + nv_reset_time_long

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time


class TransportDebug03(RabiPolarized):
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

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])


        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        pulse_sequences = []

        for tau in tau_list:

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = rf_pi_time/duty_cycle - rf_pi_time
                cooldown_time = int(int(cooldown_time/100)*100)

            pulse_sequence = []
            start_time = 0
            for part in range(2):
                # Initialize nuclear and electronic spin at the beginning of each sequence
                pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
                current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
                for iteration in range(self.settings['polarization_iterations']):
                    pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                    if rf_pi_time > 0:
                        pulse_sequence.append(Pulse(rf_channel, current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                        current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                    for i in range(self.settings['initialization_laser']['n']):
                        pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                        current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

                pulse_sequence.append(Pulse(microwave_channel_2, current_time + laser_off_time, mw_2_pi_time))
                if tau > 0:
                    pulse_sequence.append(Pulse(rf_channel, current_time + laser_off_time + mw_2_pi_time + rf_settle, tau))
                if part == 0:
                    pulse_sequence.append(Pulse(microwave_channel_2,
                                                current_time + laser_off_time + mw_2_pi_time + rf_settle + tau + rf_settle,
                                                mw_2_pi_time))
                    current_time_2 = current_time + laser_off_time + mw_2_pi_time + rf_settle + tau + rf_settle + mw_2_pi_time

                elif part == 1:
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time + laser_off_time + mw_2_pi_time + rf_settle + tau + rf_settle,
                                                mw_1_pi_time))

                    current_time_2 = current_time + laser_off_time + mw_2_pi_time + rf_settle + tau + rf_settle + mw_1_pi_time + delay_mw_readout

                pulse_sequence.append(Pulse('laser', current_time_2,
                                            nv_reset_time_long))
                pulse_sequence.append(Pulse('apd_readout',
                                            current_time_2 + delay_readout,
                                            meas_time))

                if part == 1:
                    pulse_sequence.append(Pulse('apd_readout',
                                                current_time_2 + nv_reset_time_long - meas_time,
                                                meas_time))

                # Update start time for the part 2 of the sequence (where we look at the other nuclear spin population)
                start_time = current_time_2 + nv_reset_time_long

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time


class NuclearRabiRepetitiveReadout(RabiPolarized):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            # Parameter('time_step', 5., [2.5, 4., 5., 10., 20., 50., 100., 200., 500., 1000., 2000., 5000., 10000., 20000., 100000., 500000.],
            #      'time step increment of rabi pulse duration (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'],
                          'use a constant cooldown time or set it based on length of entire sequence'),
                Parameter('cooldown_time', 1, float,
                          'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                          'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int,
                      'num of initialization laser pulses, set n=1 and wait_duration=0 for just one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            # Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int,
                      'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 10, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 500, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 100, float, 'short wait time between each readout laser pulse')
            ])
        ]),
        Parameter('rf_switch', [
            Parameter('add', True, bool,
                      'check to add mw switch to every i and q pulse and to use switch to carve out pulses. note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 50, int,
                      'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'rf_switch', ['rf_switch', 'rf_iq'],
                      'determines if rf pulses are carved out by mw-switch or by i and q channels of rf source '),
            Parameter('no_iq_overlap', True, bool,
                      'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ])
    ]

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

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])


        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        self.num_ref_reads = 3
        pulse_sequences = []

        for tau in tau_list:

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = rf_pi_time/duty_cycle - rf_pi_time
                cooldown_time = int(int(cooldown_time/100)*100)

            pulse_sequence = []
            start_time = 0

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
            for i in range(self.num_ref_reads):
                pulse_sequence.append(
                    Pulse('apd_readout', start_time + laser_off_time + cooldown_time + nv_reset_time_long - (i+1) * (meas_time+delay_readout),
                          meas_time))
            current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse(rf_channel, current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

            pulse_sequence.append(Pulse(microwave_channel_2, current_time + laser_off_time, mw_2_pi_time))
            if tau > 0:
                pulse_sequence.append(Pulse(rf_channel, current_time + laser_off_time + mw_2_pi_time + rf_settle, tau))

            current_time = current_time + laser_off_time + mw_2_pi_time + rf_settle + tau + rf_settle
            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                if m == 0:
                    pulse_sequence.append(Pulse(microwave_channel_2,
                                                current_time,
                                                mw_2_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                                            meas_time))
                    current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                else:
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time,
                                                mw_1_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                                                meas_time))
                    current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

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

        if 'fits' in data.keys() and data['fits'] is not None and 1 == 2:
            counts = data['counts'][:, 1] / data['counts'][:, 0]
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
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            if 'counts' in data.keys():
                avg_counts = np.transpose(np.array([np.average(self.data['counts'][:, self.num_ref_reads:], axis=1)]))
                ref_counts = np.transpose(np.array([np.average(self.data['counts'][:, 0:self.num_ref_reads], axis=1)]))
                first_counts = np.transpose(np.array([np.average(self.data['counts'][:, self.num_ref_reads:self.num_ref_reads + 1], axis=1)]))
                #ref_counts = self.data['counts'][:, 0:1]
                #first_counts = self.data['counts'][:, 1:2]

                #plot_1d_simple_timetrace_ns(axislist[0], data['tau'], [data['counts']])
                plot_1d_simple_timetrace_ns(axislist[0], data['tau'], [ref_counts, first_counts, avg_counts])
                plot_pulses(axislist[1], self.pulse_sequences[self.sequence_index])
            axislist[0].set_title('Nuclear Rabi rf_power:{:0.1f}dBm, rf_freq:{:0.3f} MHz'.format(self.settings['rf_pulses']['rf_power'], self.settings['rf_pulses']['rf_frequency']*1e-6))
            axislist[0].legend(labels=('Ref Fluorescence', 'First Readout', 'Avg Readout'), fontsize=8)

    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''
        #        if self.scripts['find_nv'].is_running:
        #            self.scripts['find_nv']._update_plot(axes_list)
        #        else:

        x_data = self.data['tau']
        #avg_counts = np.transpose(np.array([np.average(self.data['counts'][:, 1:], axis=1)]))
        #ref_counts = self.data['counts'][:, 0:1]
        #first_counts = self.data['counts'][:, 1:2]

        avg_counts = np.transpose(np.array([np.average(self.data['counts'][:, self.num_ref_reads:], axis=1)]))
        ref_counts = np.transpose(np.array([np.average(self.data['counts'][:, 0:self.num_ref_reads], axis=1)]))
        first_counts = np.transpose(np.array([np.average(self.data['counts'][:, self.num_ref_reads:self.num_ref_reads + 1], axis=1)]))

        axis1 = axes_list[0]
        if not self.data['counts'] == []:
            update_1d_simple(axis1, x_data, [ref_counts, first_counts, avg_counts])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])


class NuclearRabiConsecutive(NuclearRabiRepetitiveReadout):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns'),
            Parameter('n', 4, int, 'number of consecutive RF pulses')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'],
                          'use a constant cooldown time or set it based on length of entire sequence'),
                Parameter('cooldown_time', 1, float,
                          'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                          'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int,
                      'num of initialization laser pulses, set n=1 and wait_duration=0 for just one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            # Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int,
                      'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 10, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 500, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 100, float, 'short wait time between each readout laser pulse')
            ])
        ]),
        Parameter('rf_switch', [
            Parameter('add', True, bool,
                      'check to add mw switch to every i and q pulse and to use switch to carve out pulses. note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 50, int,
                      'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'rf_switch', ['rf_switch', 'rf_iq'],
                      'determines if rf pulses are carved out by mw-switch or by i and q channels of rf source '),
            Parameter('no_iq_overlap', True, bool,
                      'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ])
    ]

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

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])


        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        self.num_ref_reads = 3
        num_rf_pulses = self.settings['rf_pulses']['n']
        pulse_sequences = []

        for tau in tau_list:

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                rf_on_total_time = self.settings['polarization_iterations']*rf_pi_time + tau
                cooldown_time = rf_on_total_time*(1/duty_cycle - 1)
                cooldown_time = int(int(cooldown_time/100)*100)

            pulse_sequence = []
            start_time = 0

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
            for i in range(self.num_ref_reads):
                pulse_sequence.append(
                    Pulse('apd_readout', start_time + laser_off_time + cooldown_time + nv_reset_time_long - (i+1) * (meas_time+delay_readout),
                          meas_time))
            current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

            pulse_sequence.append(Pulse(microwave_channel_2, current_time + laser_off_time, mw_2_pi_time))
            if tau > 0:
                tau = int(tau/num_rf_pulses/4)*4
                for i in range(num_rf_pulses):
                    pulse_sequence.append(Pulse('rf_switch', current_time + laser_off_time + mw_2_pi_time + rf_settle, tau))
                    current_time += laser_off_time + mw_2_pi_time + rf_settle + tau + rf_settle

            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                if m == 0:
                    pulse_sequence.append(Pulse(microwave_channel_2,
                                                current_time,
                                                mw_2_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                                            meas_time))
                    current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                else:
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time,
                                                mw_1_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                                                meas_time))
                    current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

            rf_i_duration = int((current_time - start_time)/2)*2
            pulse_sequence.append(Pulse(rf_channel, start_time, rf_i_duration))

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _add_rf_switch_to_sequences(self, pulse_sequences):
        return pulse_sequences

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

        if 'fits' in data.keys() and data['fits'] is not None and 1 == 2:
            counts = data['counts'][:, 1] / data['counts'][:, 0]
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
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            if 'counts' in data.keys():
                avg_counts = np.transpose(np.array([np.average(self.data['counts'][:, self.num_ref_reads+1:], axis=1)]))
                ref_counts = np.transpose(np.array([np.average(self.data['counts'][:, 0:self.num_ref_reads], axis=1)]))
                first_counts = np.transpose(np.array([np.average(self.data['counts'][:, self.num_ref_reads:self.num_ref_reads + 1], axis=1)]))
                #ref_counts = self.data['counts'][:, 0:1]
                #first_counts = self.data['counts'][:, 1:2]

                #plot_1d_simple_timetrace_ns(axislist[0], data['tau'], [data['counts']])
                plot_1d_simple_timetrace_ns(axislist[0], data['tau'], [ref_counts, first_counts, avg_counts])
                plot_pulses(axislist[1], self.pulse_sequences[self.sequence_index])
            axislist[0].set_title('Nuclear Rabi rf_power:{:0.1f}dBm, rf_freq:{:0.3f} MHz'.format(self.settings['rf_pulses']['rf_power'], self.settings['rf_pulses']['rf_frequency']*1e-6))
            axislist[0].legend(labels=('Ref Fluorescence', 'First Readout', 'Avg Readout'), fontsize=8)

    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''
        #        if self.scripts['find_nv'].is_running:
        #            self.scripts['find_nv']._update_plot(axes_list)
        #        else:

        x_data = self.data['tau']
        #avg_counts = np.transpose(np.array([np.average(self.data['counts'][:, 1:], axis=1)]))
        #ref_counts = self.data['counts'][:, 0:1]
        #first_counts = self.data['counts'][:, 1:2]

        avg_counts = np.transpose(np.array([np.average(self.data['counts'][:, self.num_ref_reads+1:], axis=1)]))
        ref_counts = np.transpose(np.array([np.average(self.data['counts'][:, 0:self.num_ref_reads], axis=1)]))
        first_counts = np.transpose(np.array([np.average(self.data['counts'][:, self.num_ref_reads:self.num_ref_reads + 1], axis=1)]))

        axis1 = axes_list[0]
        if not self.data['counts'] == []:
            update_1d_simple(axis1, x_data, [ref_counts, first_counts, avg_counts])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])


class TransportDebug04(NuclearRabiRepetitiveReadout):
    """
    Create coherent superposition on the electron, transfer to nuclear, and then drive nuclear rotations
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_1_pi_half_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_pi_half_time', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'],
                          'use a constant cooldown time or set it based on length of entire sequence'),
                Parameter('cooldown_time', 1, float,
                          'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                          'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('t_movement', 20000, float, 'momvement/wait time in ns before reading state of the nuclear spin'),
        Parameter('storage_axis', 'z', ['x', 'y', 'z'],
                  'rotate spin onto this axis before storage onto nuclear spin; for z, there is no rotation'),
        Parameter('dynamic_decoupling', [
            Parameter('enable', True, bool, 'enable dynamically decoupled RF pulses'),
            Parameter('n', 1, int, 'number of echo units'),
            Parameter('rf_switch_offset', 0, float, 'extra fraction of a resonator period that leaks out of a gated pulse')
        ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int,
                      'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
            # Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int,
                      'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 250, int,
                      'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 10, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 500, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 100, float, 'short wait time between each readout laser pulse')
            ])
        ]),
        Parameter('rf_switch', [
            Parameter('add', True, bool,
                      'check to add mw switch to every i and q pulse and to use switch to carve out pulses. note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 50, int,
                      'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'rf_switch', ['rf_switch', 'rf_iq'],
                      'determines if rf pulses are carved out by mw-switch or by i and q channels of rf source '),
            Parameter('no_iq_overlap', True, bool,
                      'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ])
    ]

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

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_1_pi_half_time = self.settings['mw_pulses']['mw_1_pi_half_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']
        mw_2_pi_half_time = self.settings['mw_pulses']['mw_2_pi_half_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        microwave_channel_2_q = 'microwave_q_2'
        microwave_channel_1_q = 'microwave_q'

        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']
        rf_switch_offset = self.settings['dynamic_decoupling']['rf_switch_offset']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        self.num_ref_reads = 3
        pulse_sequences = []
        t_movement = self.settings['t_movement'] # Remove this when we do actual transport stuff!
        rf_period = 1/self.settings['rf_pulses']['rf_frequency']*1e9
        print(rf_period)

        for tau in tau_list:

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                rf_on_total_time = self.settings['polarization_iterations'] * rf_pi_time + tau + rf_pi_half_time
                cooldown_time = rf_on_total_time * (1 / duty_cycle - 1)
                cooldown_time = int(int(cooldown_time / 100) * 100)

            pulse_sequence = []
            start_time = 0

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))

            for i in range(self.num_ref_reads):
                pulse_sequence.append(
                    Pulse('apd_readout', start_time + laser_off_time + cooldown_time + nv_reset_time_long - (i + 1) * (
                                meas_time + delay_readout),
                          meas_time))

            current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time += nv_reset_time_pulsed + nv_reset_time_wait
            current_time += laser_off_time-nv_reset_time_wait
            # Finished nuclear polarization, now onto the main sequence

            if self.settings['storage_axis'] == 'z':
                if tau > 0:
                    pulse_sequence.append(Pulse(microwave_channel_2_q,
                                                current_time,
                                                tau))
                current_time += tau + delay_mw_readout
            elif self.settings['storage_axis'] == 'y' or self.settings['storage_axis'] == 'x':
                if tau > 0:
                    pulse_sequence.append(Pulse(microwave_channel_2_q,
                                                current_time,
                                                tau))
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time + tau + delay_mw_readout,
                                            mw_2_pi_half_time))
                current_time += tau + delay_mw_readout + mw_2_pi_half_time

            tau_half = int((int(rf_pi_time/2/rf_period)+rf_switch_offset)*rf_period/10)*10 # Make sure time is an integer multiple of RF freq and also 4 ns (to avoid PB errors)
            print(tau_half)
            if self.settings['dynamic_decoupling']['enable']:
                pulse_sequence.append(Pulse('rf_switch',
                                            current_time,
                                            tau_half))
                pulse_sequence.append(Pulse(microwave_channel,
                                            current_time + tau_half + rf_settle,
                                            mw_1_pi_time))
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time + tau_half + rf_settle + mw_1_pi_time + delay_mw_readout,
                                            mw_2_pi_time))
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time + tau_half + rf_settle + mw_1_pi_time + delay_mw_readout + mw_2_pi_time + tau_half * 2,
                                            mw_2_pi_time))
                pulse_sequence.append(Pulse(microwave_channel,
                                            current_time + tau_half + rf_settle + mw_1_pi_time + delay_mw_readout + mw_2_pi_time + tau_half * 2 + mw_2_pi_time + delay_mw_readout,
                                            mw_1_pi_time))
                pulse_sequence.append(Pulse('rf_switch',
                                            current_time + tau_half + rf_settle + mw_1_pi_time + delay_mw_readout + mw_2_pi_time + tau_half * 2 + mw_2_pi_time + delay_mw_readout + mw_1_pi_time + rf_settle,
                                            tau_half))

                current_time += tau_half + rf_settle + mw_1_pi_time + delay_mw_readout + mw_2_pi_time + tau_half * 2 + mw_2_pi_time + \
                                delay_mw_readout + mw_1_pi_time + rf_settle + tau_half + rf_settle

            else:
                pulse_sequence.append(Pulse('rf_switch',
                                            current_time,
                                            rf_pi_time))
                current_time += rf_pi_time + rf_settle

            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time,
                                        mw_2_pi_time))
            current_time += mw_2_pi_time

            if self.settings['storage_axis'] == 'x' or self.settings['storage_axis'] == 'y':
                pulse_sequence.append(Pulse('rf_switch',
                                            current_time,
                                            rf_pi_half_time))
                current_time += rf_pi_half_time + rf_settle + t_movement
            else:
                current_time += t_movement

            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                if m == 0:
                    pulse_sequence.append(Pulse(microwave_channel_2,
                                                current_time,
                                                mw_2_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                                            meas_time))
                    current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                else:
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time,
                                                mw_1_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                                                meas_time))
                    current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
            pulse_sequence.append(Pulse(rf_channel, 0, current_time))
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time


    def _add_rf_switch_to_sequences(self, pulse_sequences):
        return pulse_sequences


class NuclearRabiDecoherenceProtected(NuclearRabiConsecutive):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_1_pi_half_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_1_three_pi_half_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns'),
            Parameter('n', 4, int, 'number of consecutive RF pulses')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            # Parameter('time_step', 5., [2.5, 4., 5., 10., 20., 50., 100., 200., 500., 1000., 2000., 5000., 10000., 20000., 100000., 500000.],
            #      'time step increment of rabi pulse duration (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'],
                          'use a constant cooldown time or set it based on length of entire sequence'),
                Parameter('cooldown_time', 1, float,
                          'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                          'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int,
                      'num of initialization laser pulses, set n=1 and wait_duration=0 for just one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            # Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int,
                      'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 100, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 16, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 400, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 100, float, 'short wait time between each readout laser pulse')
            ])
        ]),
        Parameter('rf_switch', [
            Parameter('add', True, bool,
                      'check to add mw switch to every i and q pulse and to use switch to carve out pulses. note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 60, int,
                      'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'rf_switch', ['rf_switch', 'rf_iq'],
                      'determines if rf pulses are carved out by mw-switch or by i and q channels of rf source '),
            Parameter('no_iq_overlap', True, bool,
                      'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ])
    ]

    def _create_pulse_sequences_van_der_sar(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        #short_pulses = [x for x in tau_list if x < min_pulse_dur]
        #print('Found short pulses: ', short_pulses)
        #tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']

        pulse_sequences = []
        self.num_ref_reads = 3

        for n in tau_list:
            tau = int((6+0.5)/self.settings['rf_pulses']['rf_frequency']*1e9)
            tau = int(tau/4)*4

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = tau*4*n/duty_cycle - tau*4*n
                cooldown_time = int(int(cooldown_time/100)*100)

            pulse_sequence = []
            start_time = 0

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time,
                                        nv_reset_time_long))

            current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence

            pulse_sequence.append(Pulse(microwave_channel, current_time + laser_off_time, mw_1_pi_time))
            current_time += laser_off_time + mw_1_pi_time + delay_mw_readout

            if n > 0:
                pulse_sequence.append(Pulse(rf_channel,
                                            current_time,
                                            tau * 4 * int(n)))
                for iteration in range(int(n)):
                    pulse_sequence.append(Pulse(microwave_channel_2,
                                                current_time + tau,
                                                mw_2_pi_time))
                    pulse_sequence.append(Pulse(microwave_channel_2,
                                                current_time + tau + tau * 2,
                                                mw_2_pi_time))
                    current_time += tau*4

            pulse_sequence.append(Pulse(microwave_channel,
                                        current_time,
                                        mw_1_pi_time))
            pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                        nv_reset_time_long))
            pulse_sequence.append(Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_long - meas_time,
                                        meas_time))
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

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

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        #mw_1_pi_half_time = self.settings['mw_pulses']['mw_1_pi_half_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']
        #mw_2_pi_half_time = self.settings['mw_pulses']['mw_2_pi_half_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']
        self.num_ref_reads = 3
        pulse_sequences = []
        #t_movement = self.settings['t_movement'] # Remove this when we do actual transport stuff!
        rf_period = 1/self.settings['rf_pulses']['rf_frequency']*1e9

        for tau in tau_list:

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                rf_on_total_time = self.settings['polarization_iterations'] * rf_pi_time + tau + rf_pi_half_time
                cooldown_time = rf_on_total_time * (1 / duty_cycle - 1)
                cooldown_time = int(int(cooldown_time / 100) * 100)

            pulse_sequence = []
            start_time = 0

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))

            for i in range(self.num_ref_reads):
                pulse_sequence.append(
                    Pulse('apd_readout', start_time + laser_off_time + cooldown_time + nv_reset_time_long - (i + 1) * (
                                meas_time + delay_readout),
                          meas_time))

            current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time += nv_reset_time_pulsed + nv_reset_time_wait
            current_time += laser_off_time-nv_reset_time_wait
            # Finished nuclear polarization, now onto the main sequence

            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time,
                                        mw_2_pi_time))
            current_time += mw_2_pi_time + delay_mw_readout

            tau_half = int(int(tau/2/rf_period)*rf_period/10)*10 # Make sure time is an integer multiple of RF freq and also 4 ns (to avoid PB errors)

            if tau > 0:
                pulse_sequence.append(Pulse('rf_switch',
                                            current_time,
                                            tau_half))
                pulse_sequence.append(Pulse(microwave_channel,
                                            current_time + tau_half + rf_settle,
                                            mw_1_pi_time))
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time + tau_half + rf_settle + mw_1_pi_time + delay_mw_readout,
                                            mw_2_pi_time))
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time + tau_half + rf_settle + mw_1_pi_time + delay_mw_readout + mw_2_pi_time + tau_half * 2,
                                            mw_2_pi_time))
                pulse_sequence.append(Pulse(microwave_channel,
                                            current_time + tau_half + rf_settle + mw_1_pi_time + delay_mw_readout + mw_2_pi_time + tau_half*2 + mw_2_pi_time + delay_mw_readout,
                                            mw_1_pi_time))

                pulse_sequence.append(Pulse('rf_switch',
                                            current_time + tau_half + rf_settle + mw_1_pi_time + delay_mw_readout + mw_2_pi_time + tau_half * 2 + mw_2_pi_time + delay_mw_readout + mw_1_pi_time + rf_settle,
                                            tau_half))
                current_time += tau_half + rf_settle + mw_1_pi_time + delay_mw_readout + mw_2_pi_time + tau_half * 2 + mw_2_pi_time + delay_mw_readout + mw_1_pi_time + rf_settle + tau_half + rf_settle

            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                if m == 0:
                    pulse_sequence.append(Pulse(microwave_channel_2,
                                                current_time,
                                                mw_2_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                                            meas_time))
                    current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                else:
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time,
                                                mw_1_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                                                meas_time))
                    current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
            pulse_sequence.append(Pulse(rf_channel, start_time, int((current_time-start_time)/2)*2))
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time


    def _add_rf_switch_to_sequences(self, pulse_sequences):
        return pulse_sequences


class NuclearSwapDecoherenceProtected(NuclearRabiConsecutive):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_1_pi_half_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_1_three_pi_half_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_pi_half_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_three_pi_half_time', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_time_protected', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns'),
            Parameter('n', 4, int, 'number of consecutive RF pulses')
        ]),
        Parameter('t_entangled', 100, float, 'time to accumulate phase'),
        Parameter('storage_axis', 'z', ['x', 'y', 'z'],
                  'rotate spin onto this axis before storage onto nuclear spin; for z, there is no rotation'),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            # Parameter('time_step', 5., [2.5, 4., 5., 10., 20., 50., 100., 200., 500., 1000., 2000., 5000., 10000., 20000., 100000., 500000.],
            #      'time step increment of rabi pulse duration (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'],
                          'use a constant cooldown time or set it based on length of entire sequence'),
                Parameter('cooldown_time', 1, float,
                          'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                          'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int,
                      'num of initialization laser pulses, set n=1 and wait_duration=0 for just one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            # Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int,
                      'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 100, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 16, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 400, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 100, float, 'short wait time between each readout laser pulse')
            ])
        ]),
        Parameter('rf_switch', [
            Parameter('add', True, bool,
                      'check to add mw switch to every i and q pulse and to use switch to carve out pulses. note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 60, int,
                      'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'rf_switch', ['rf_switch', 'rf_iq'],
                      'determines if rf pulses are carved out by mw-switch or by i and q channels of rf source '),
            Parameter('no_iq_overlap', True, bool,
                      'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ])
    ]

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

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_1_pi_half_time = self.settings['mw_pulses']['mw_1_pi_half_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']
        mw_2_pi_half_time = self.settings['mw_pulses']['mw_2_pi_half_time']
        mw_2_three_pi_half_time = self.settings['mw_pulses']['mw_2_three_pi_half_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        microwave_channel_2_q = 'microwave_q_2'
        microwave_channel_1_q = 'microwave_q'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_time_protected = self.settings['rf_pulses']['rf_pi_time_protected']
        self.num_ref_reads = 3
        pulse_sequences = []
        t_entangled = self.settings['t_entangled'] # Remove this when we do actual transport stuff!
        rf_period = 1/self.settings['rf_pulses']['rf_frequency']*1e9

        for tau in tau_list:
            start_time = 0
            pulse_sequence = []
            for i in range(2):

                if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                    cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
                elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                    duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                    rf_on_total_time = self.settings['polarization_iterations'] * rf_pi_time + tau
                    cooldown_time = rf_on_total_time * (1 / duty_cycle - 1)
                    cooldown_time = int(int(cooldown_time / 100) * 100)

                # Initialize nuclear and electronic spin at the beginning of each sequence
                pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))

                current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
                for iteration in range(self.settings['polarization_iterations']):
                    pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                    if rf_pi_time > 0:
                        pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time + rf_settle, rf_pi_time))
                        current_time = current_time + mw_1_pi_time + rf_settle + rf_pi_time + rf_settle
                    for i in range(self.settings['initialization_laser']['n']):
                        pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                        current_time += nv_reset_time_pulsed + nv_reset_time_wait
                current_time += laser_off_time-nv_reset_time_wait
                # Finished nuclear polarization, now onto the main sequence


                if self.settings['storage_axis'] == 'z':
                    if i == 1:
                        pulse_sequence.append(Pulse(microwave_channel_2,
                                                    current_time,
                                                    mw_2_pi_time))
                        current_time += mw_2_pi_time + delay_mw_readout
                elif self.settings['storage_axis'] == 'x':
                    if i == 0:
                        pulse_sequence.append(Pulse(microwave_channel_2,
                                                    current_time,
                                                    mw_2_pi_half_time))
                        current_time += mw_2_pi_half_time + delay_mw_readout
                    elif i == 1:
                        pulse_sequence.append(Pulse(microwave_channel_2,
                                                    current_time,
                                                    mw_2_three_pi_half_time))
                        current_time += mw_2_three_pi_half_time + delay_mw_readout
                elif self.settings['storage_axis'] == 'y':
                    raise NotImplementedError

                tau_half = int(int(rf_pi_time_protected/2/rf_period)*rf_period/10)*10 # Make sure time is an integer multiple of RF freq and also 4 ns (to avoid PB errors)


                pulse_sequence.append(Pulse('rf_switch',
                                            current_time,
                                            tau_half))
                pulse_sequence.append(Pulse(microwave_channel,
                                            current_time + tau_half + rf_settle,
                                            mw_1_pi_time))
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time + tau_half + rf_settle + mw_1_pi_time + delay_mw_readout,
                                            mw_2_pi_time))
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time + tau_half + rf_settle + mw_1_pi_time + delay_mw_readout + mw_2_pi_time + tau_half * 2,
                                            mw_2_pi_time))
                pulse_sequence.append(Pulse(microwave_channel,
                                            current_time + tau_half + rf_settle + mw_1_pi_time + delay_mw_readout + mw_2_pi_time + tau_half*2 + mw_2_pi_time + delay_mw_readout,
                                            mw_1_pi_time))
                pulse_sequence.append(Pulse('rf_switch',
                                            current_time + tau_half + rf_settle + mw_1_pi_time + delay_mw_readout + mw_2_pi_time + tau_half * 2 + mw_2_pi_time + delay_mw_readout + mw_1_pi_time + rf_settle,
                                            tau_half))
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time + tau_half + rf_settle + mw_1_pi_time + delay_mw_readout + mw_2_pi_time + tau_half * 2 + mw_2_pi_time +
                                            delay_mw_readout + mw_1_pi_time + rf_settle + tau_half + rf_settle + t_entangled,
                                            mw_2_pi_time))
                current_time += tau_half + rf_settle + mw_1_pi_time + delay_mw_readout + mw_2_pi_time + tau_half * 2 + mw_2_pi_time + \
                                delay_mw_readout + mw_1_pi_time + rf_settle + tau_half + rf_settle + t_entangled + mw_2_pi_time + delay_mw_readout


                if tau > 0:
                    pulse_sequence.append(Pulse('rf_switch',
                                                current_time,
                                                tau))
                    current_time += tau + rf_settle

                # readout
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time,
                                            mw_2_pi_time))
                pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                            nv_reset_time_short))
                pulse_sequence.append(
                    Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                          meas_time))
                current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                start_time = current_time
            pulse_sequence.append(Pulse(rf_channel, 0, int((current_time-0)/2)*2))
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time


    def _add_rf_switch_to_sequences(self, pulse_sequences):
        return pulse_sequences

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

        if 1 == 2 and 'fits' in data.keys() and data['fits'] is not None:
            counts = data['counts'][:, 1] / data['counts'][:, 0]
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
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            super(NuclearRabi, self)._plot(axislist)
            axislist[0].set_title('Nuclear Rabi rf_power:{:0.1f}dBm, rf_freq:{:0.3f} MHz'.format(self.settings['rf_pulses']['rf_power'], self.settings['rf_pulses']['rf_frequency']*1e-6))
            axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)

    def _update_plot(self, axislist):
        super(NuclearRabi, self)._update_plot(axislist)


class RamseyNuclearPolarized(NuclearSwapDecoherenceProtected):


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

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_1_pi_half_time = self.settings['mw_pulses']['mw_1_pi_half_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']
        mw_2_pi_half_time = self.settings['mw_pulses']['mw_2_pi_half_time']
        mw_2_three_pi_half_time = self.settings['mw_pulses']['mw_2_three_pi_half_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_time_protected = self.settings['rf_pulses']['rf_pi_time_protected']
        self.num_ref_reads = 3
        pulse_sequences = []
        t_entangled = self.settings['t_entangled'] # Remove this when we do actual transport stuff!
        rf_period = 1/self.settings['rf_pulses']['rf_frequency']*1e9

        for tau in tau_list:
            start_time = 0
            pulse_sequence = []
            for i in range(2):

                if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                    cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
                elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                    duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                    rf_on_total_time = self.settings['polarization_iterations'] * rf_pi_time + tau
                    cooldown_time = rf_on_total_time * (1 / duty_cycle - 1)
                    cooldown_time = int(int(cooldown_time / 100) * 100)

                # Initialize nuclear and electronic spin at the beginning of each sequence
                pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))

                current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
                for iteration in range(self.settings['polarization_iterations']):
                    pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                    if rf_pi_time > 0:
                        pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time + delay_mw_readout, rf_pi_time))
                        current_time = current_time + mw_1_pi_time + delay_mw_readout + rf_pi_time + rf_settle
                    for i in range(self.settings['initialization_laser']['n']):
                        pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                        current_time += nv_reset_time_pulsed + nv_reset_time_wait
                current_time += laser_off_time-nv_reset_time_wait
                # Finished nuclear polarization, now onto the main sequence

                if i == 0:
                    pulse_sequence.append(Pulse(microwave_channel_2, current_time, mw_2_pi_half_time))
                    pulse_sequence.append(Pulse(microwave_channel_2, current_time + mw_2_pi_half_time + tau, mw_2_pi_half_time))
                    current_time = current_time + mw_2_pi_half_time + tau + mw_2_pi_half_time + delay_mw_readout

                elif i == 1:
                    pulse_sequence.append(Pulse(microwave_channel_2, current_time, mw_2_pi_half_time))
                    pulse_sequence.append(Pulse(microwave_channel_2, current_time + mw_2_pi_half_time + tau, mw_2_three_pi_half_time))
                    current_time = current_time + mw_2_pi_half_time + tau + mw_2_three_pi_half_time + delay_mw_readout

                pulse_sequence += [
                    Pulse('laser', current_time, nv_reset_time_long),
                    Pulse('apd_readout', current_time + delay_readout, meas_time),
                ]

                current_time += nv_reset_time_long + laser_off_time
                start_time = current_time

            pulse_sequence.append(Pulse(rf_channel, 0, int((current_time-0)/2)*2))
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time


class HahnEchoNuclearPolarized(NuclearSwapDecoherenceProtected):

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

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_1_pi_half_time = self.settings['mw_pulses']['mw_1_pi_half_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']
        mw_2_pi_half_time = self.settings['mw_pulses']['mw_2_pi_half_time']
        mw_2_three_pi_half_time = self.settings['mw_pulses']['mw_2_three_pi_half_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_time_protected = self.settings['rf_pulses']['rf_pi_time_protected']
        self.num_ref_reads = 3
        pulse_sequences = []
        t_entangled = self.settings['t_entangled']  # Remove this when we do actual transport stuff!
        rf_period = 1 / self.settings['rf_pulses']['rf_frequency'] * 1e9

        for tau in tau_list:
            start_time = 0
            pulse_sequence = []
            for i in range(2):

                if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                    cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
                elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                    duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                    rf_on_total_time = self.settings['polarization_iterations'] * rf_pi_time + tau
                    cooldown_time = rf_on_total_time * (1 / duty_cycle - 1)
                    cooldown_time = int(int(cooldown_time / 100) * 100)

                # Initialize nuclear and electronic spin at the beginning of each sequence
                pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))

                current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
                for iteration in range(self.settings['polarization_iterations']):
                    pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                    if rf_pi_time > 0:
                        pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time + delay_mw_readout, rf_pi_time))
                        current_time = current_time + mw_1_pi_time + delay_mw_readout + rf_pi_time + rf_settle
                    for i in range(self.settings['initialization_laser']['n']):
                        pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                        current_time += nv_reset_time_pulsed + nv_reset_time_wait
                current_time += laser_off_time - nv_reset_time_wait
                # Finished nuclear polarization, now onto the main sequence

                if i == 0:
                    pulse_sequence.append(Pulse(microwave_channel_2, current_time, mw_2_pi_half_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, current_time + mw_2_pi_half_time + tau, mw_2_pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, current_time + mw_2_pi_half_time + tau + mw_2_pi_time + tau, mw_2_pi_half_time))
                    current_time = current_time + mw_2_pi_half_time + tau +  mw_2_pi_time + tau + mw_2_pi_half_time + delay_mw_readout

                elif i == 1:
                    pulse_sequence.append(Pulse(microwave_channel_2, current_time, mw_2_pi_half_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, current_time + mw_2_pi_half_time + tau, mw_2_pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, current_time + mw_2_pi_half_time + tau + mw_2_pi_time + tau,
                              mw_2_three_pi_half_time))
                    current_time = current_time + mw_2_pi_half_time + tau + mw_2_pi_time + tau + mw_2_three_pi_half_time + delay_mw_readout

                pulse_sequence += [
                    Pulse('laser', current_time, nv_reset_time_long),
                    Pulse('apd_readout', current_time + delay_readout, meas_time),
                ]

                current_time += nv_reset_time_long + laser_off_time
                start_time = current_time

            pulse_sequence.append(Pulse(rf_channel, 0, int((current_time - 0) / 2) * 2))
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time


class TransportDqma(RabiPolarized):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns'),
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_time_i', 5000, float, 'time for RF pulse to drive nuclear spin in ns (when rf_i is on)'),
            Parameter('rf_pi_half_time_i', 5000, float, 'time for RF pulse to drive nuclear spin in ns (when rf_i is on)'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'],
                          'use a constant cooldown time or set it based on length of entire sequence'),
                Parameter('cooldown_time', 1, float,
                          'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                          'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('t_movement', 0, float, 'momvement/wait time in ns before reading state of the nuclear spin'),
        Parameter('dynamic_decoupling', [
            Parameter('sequence', 'none', ['none', 'echo', 'cpmg', 'xy'], 'dynamic decoupling sequence'),
            Parameter('k', 1, int, 'number of pulses in the decoupling sequence; this only applies for CPMG-k and XY-k')
            ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int, 'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time_long', 5000, int, 'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 100, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 8, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 400, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 250, float, 'short wait time between each readout laser pulse')
            ])
        ])
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator,
                    'afg': AFG3022C, 'mw_gen_2': MicrowaveGenerator2, 'commander': Commander}

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

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        microwave_channel, microwave_channel_2, microwave_channel_q, microwave_channel_2_q, rf_channel = \
            'microwave_i', 'microwave_i_2', 'microwave_q', 'microwave_q_2', 'rf_i'

        mw_1_pi_time, mw_2_pi_time = [self.settings['mw_pulses'][key] for key in ['mw_1_pi_time', 'mw_2_pi_time']]

        meas_time, nv_reset_time_long, laser_off_time, delay_mw_readout, delay_readout = \
            [self.settings['read_out'][key] for key in
             ['meas_time', 'nv_reset_time_long', 'laser_off_time', 'delay_mw_readout', 'delay_readout']]

        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']

        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']

        rf_pi_time, rf_pi_half_time, rf_pi_time_i, rf_pi_half_time_i, rf_settle = \
            [self.settings['rf_pulses'][key] for key in ['rf_pi_time', 'rf_pi_half_time', 'rf_pi_time_i', 'rf_pi_half_time_i', 'rf_settle']]

        t_movement = self.settings['t_movement']

        pulse_sequences = []
        for tau in tau_list:
            t_sense = tau

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = (rf_pi_time*self.settings['polarization_iterations']+rf_pi_half_time*2)*(1/duty_cycle - 1)
                cooldown_time = int(int(cooldown_time/100)*100)

            pulse_sequence = []
            start_time = 0

            for part in range(2):
                # Initialize nuclear and electronic spin at the beginning of each sequence
                pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
                current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
                for iteration in range(self.settings['polarization_iterations']):
                    pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                    if rf_pi_time > 0:
                        pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time, rf_pi_time))
                        current_time = current_time + mw_1_pi_time + rf_pi_time + rf_settle
                    for i in range(self.settings['initialization_laser']['n']):
                        pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                        current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

                # Flip e-spin to m_s=1 state so that we can drive RF
                pulse_sequence.append(Pulse(microwave_channel_2, current_time + laser_off_time, mw_2_pi_time))

                # Prepare nuclear spin in superposition
                pulse_sequence.append(Pulse('rf_switch', current_time + laser_off_time + mw_2_pi_time + delay_mw_readout, rf_pi_half_time))
                current_time += laser_off_time + mw_2_pi_time + delay_mw_readout + rf_pi_half_time + rf_settle
                # Access nuclear memory; choose to accumulate conditional on which nuclear state

                channel_odd, channel_even = microwave_channel_2_q, microwave_channel_q
                pi_time_odd, pi_time_even = mw_2_pi_time, mw_1_pi_time

                if self.settings['dynamic_decoupling']['sequence'] == 'none':

                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time,
                                                pi_time_odd))

                    # Wait for t_sense and disentangle from memory
                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time + pi_time_odd + t_sense,
                                                pi_time_odd))

                    current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout

                elif self.settings['dynamic_decoupling']['sequence'] == 'echo':
                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time,
                                                pi_time_odd))

                    # Wait for t_sense and disentangle from memory
                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time + pi_time_odd + t_sense,
                                                pi_time_odd))

                    current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout

                    pulse_sequence.append(Pulse(channel_even,
                                                current_time,
                                                pi_time_even))

                    # Wait for t_sense and disentangle from memory
                    pulse_sequence.append(Pulse(channel_even,
                                                current_time + pi_time_even + t_sense,
                                                pi_time_even))

                    current_time += pi_time_even + t_sense + pi_time_even + delay_mw_readout

                # Wait for t_movement and rotate nuclear spin to Z axis for readout
                if part == 1:
                    pulse_sequence.append(
                        Pulse(rf_channel,
                              current_time + t_movement - 60,
                              rf_pi_half_time_i + 2*60))
                pulse_sequence.append(
                    Pulse('rf_switch',
                          current_time + t_movement,
                          rf_pi_half_time_i))

                current_time += t_movement + rf_pi_half_time_i + rf_settle

                if part == 10:

                    for m in range(self.settings['read_out']['repetitive_readout']['m']):
                        if m == 0:
                            pulse_sequence.append(Pulse(microwave_channel_2,
                                                        current_time,
                                                        mw_2_pi_time))
                            pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                        nv_reset_time_short))
                            pulse_sequence.append(
                                Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                                      meas_time))
                            current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                        else:
                            pulse_sequence.append(Pulse(microwave_channel,
                                                        current_time,
                                                        mw_1_pi_time))
                            pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                        nv_reset_time_short))
                            pulse_sequence.append(
                                Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                                      meas_time))
                            current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

                    # Update start time for the part 2 of the sequence (where we look at the other nuclear spin population)
                    start_time = current_time

                else:
                    for m in range(self.settings['read_out']['repetitive_readout']['m']):
                        if m == 0:
                            pulse_sequence.append(Pulse(microwave_channel,
                                                        current_time,
                                                        mw_1_pi_time))
                            pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                        nv_reset_time_short))
                            pulse_sequence.append(
                                Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                                      meas_time))
                            current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                        else:
                            pulse_sequence.append(Pulse(microwave_channel_2,
                                                        current_time,
                                                        mw_2_pi_time))
                            pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                        nv_reset_time_short))
                            pulse_sequence.append(
                                Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                                      meas_time))
                            current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

                    start_time = current_time

            # Get the length of the whole sequence. This basically determines the frequency of the flexure stage momvement. We'd then put in this total duration
            # into the FieldProfile script to measure what the profile is under the same timing.
            self.log('Total length of pulse sequence: %i' % current_time)
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _add_rf_switch_to_sequences(self, pulse_sequences):
        return pulse_sequences

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

        if 'fits' in data.keys() and data['fits'] is not None and 1 == 2:
            counts = data['counts'][:, 1] / data['counts'][:, 0]
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
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            if 'counts' in data.keys():
                num_daq_reads = self.settings['read_out']['repetitive_readout']['m']
                avg_counts_1 = np.transpose(np.array([np.average(self.data['counts'][:, 0:num_daq_reads], axis=1)]))
                first_counts_1 = np.transpose(np.array([np.average(self.data['counts'][:, 0:1], axis=1)]))
                avg_counts_2 = np.transpose(np.array([np.average(self.data['counts'][:, num_daq_reads:], axis=1)]))
                first_counts_2 = np.transpose(np.array([np.average(self.data['counts'][:, num_daq_reads:num_daq_reads+1], axis=1)]))

                #plot_1d_simple_timetrace_ns(axislist[0], data['tau'], [first_counts_1, first_counts_2, avg_counts_1, avg_counts_2])
                plot_1d_simple_timetrace_ns(axislist[0], data['tau'],[avg_counts_1, avg_counts_2])
                plot_pulses(axislist[1], self.pulse_sequences[self.sequence_index])
            axislist[0].set_title('Coherent Transport w/ Direct Quantum Memory Access')
            #[0].legend(labels=('Nuclear State 0 (first readout)', 'Nuclear State 1 (first readout)', 'Nuclear State 0 (avg readout)', 'Nuclear State 1 (avg readout)'), fontsize=8)
            axislist[0].legend(labels=('Nuclear State 0 (avg readout)', 'Nuclear State 1 (avg readout)'), fontsize=8)

    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''
        #        if self.scripts['find_nv'].is_running:
        #            self.scripts['find_nv']._update_plot(axes_list)
        #        else:

        x_data = self.data['tau']


        num_daq_reads = self.settings['read_out']['repetitive_readout']['m']
        avg_counts_1 = np.transpose(np.array([np.average(self.data['counts'][:, 0:num_daq_reads], axis=1)]))
        first_counts_1 = np.transpose(np.array([np.average(self.data['counts'][:, 0:1], axis=1)]))
        avg_counts_2 = np.transpose(np.array([np.average(self.data['counts'][:, num_daq_reads:], axis=1)]))
        first_counts_2 = np.transpose(np.array([np.average(self.data['counts'][:, num_daq_reads:num_daq_reads + 1], axis=1)]))

        axis1 = axes_list[0]
        if not self.data['counts'] == []:
            update_1d_simple(axis1, x_data, [avg_counts_1, avg_counts_2])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])


class TransportDqmaCphase(TransportDqma):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_1_pi_half_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_1_three_pi_half_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_pi_half_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_three_pi_half_time', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'],
                          'use a constant cooldown time or set it based on length of entire sequence'),
                Parameter('cooldown_time', 1, float,
                          'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                          'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('t_sense', 20000, float, 'time to accumulate phase in the entangled state'),
        Parameter('t_movement', 20000, float, 'momvement/wait time in ns before reading state of the nuclear spin'),
        Parameter('echo', 'none', ['none', 'add', 'subtract']),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 3000, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 0, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int,
                      'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
            # Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int,
                      'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 250, int,
                      'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 10, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 500, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 100, float, 'short wait time between each readout laser pulse')
            ])
        ]),
        Parameter('rf_switch', [
            Parameter('add', True, bool,
                      'check to add mw switch to every i and q pulse and to use switch to carve out pulses. note that iq pulses become longer by 2*extra-time'),
            Parameter('extra_time', 50, int,
                      'extra time that is added before and after the time of the i/q pulses in ns'),
            Parameter('gating', 'rf_switch', ['rf_switch', 'rf_iq'],
                      'determines if rf pulses are carved out by mw-switch or by i and q channels of rf source '),
            Parameter('no_iq_overlap', True, bool,
                      'Toggle to check for overlapping i q output. In general i and q channels should not be on simultaneously.')
        ])
    ]

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

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])


        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_1_pi_time = self.settings['mw_pulses']['mw_1_pi_time']
        mw_2_pi_time = self.settings['mw_pulses']['mw_2_pi_time']
        mw_1_pi_half_time = self.settings['mw_pulses']['mw_1_pi_half_time']
        mw_2_pi_half_time = self.settings['mw_pulses']['mw_2_pi_half_time']
        mw_1_three_pi_half_time = self.settings['mw_pulses']['mw_1_three_pi_half_time']
        mw_2_three_pi_half_time = self.settings['mw_pulses']['mw_2_three_pi_half_time']

        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        microwave_channel_2_q = 'microwave_q_2'
        microwave_channel_1_q = 'microwave_q'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        rf_pi_time = self.settings['rf_pulses']['rf_pi_time']
        rf_pi_half_time = self.settings['rf_pulses']['rf_pi_half_time']

        t_sense = self.settings['t_sense']
        t_movement = self.settings['t_movement']
        pulse_sequences = []

        for tau in tau_list:
            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = (rf_pi_time*self.settings['polarization_iterations']+rf_pi_half_time*2)*(1/duty_cycle - 1)
                cooldown_time = int(int(cooldown_time/100)*100)

            pulse_sequence = []
            start_time = 0
            for part in range(2):
                # Initialize nuclear and electronic spin at the beginning of each sequence
                pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
                current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
                for iteration in range(self.settings['polarization_iterations']):
                    pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                    if rf_pi_time > 0:
                        pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time, rf_pi_time))
                        current_time = current_time + mw_1_pi_time + rf_pi_time + rf_settle
                    for i in range(self.settings['initialization_laser']['n']):
                        pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                        current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

                # Flip e-spin to m_s=1 state so that we can drive RF
                pulse_sequence.append(Pulse(microwave_channel_2, current_time + laser_off_time, mw_2_pi_time))

                # Prepare nuclear spin in superposition
                pulse_sequence.append(Pulse('rf_switch', current_time + laser_off_time + mw_2_pi_time + delay_mw_readout, rf_pi_half_time))
                current_time += laser_off_time + mw_2_pi_time + delay_mw_readout + rf_pi_half_time + rf_settle
                # Access nuclear memory; choose to accumulate conditional on which nuclear state

                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time,
                                            mw_2_pi_half_time))
                if tau > 0:
                    pulse_sequence.append(Pulse(microwave_channel_2_q,
                                                current_time + mw_2_pi_half_time + delay_mw_readout,
                                                tau))
                pulse_sequence.append(Pulse(microwave_channel_2,
                                            current_time + mw_2_pi_half_time + delay_mw_readout + tau + delay_mw_readout,
                                            mw_2_three_pi_half_time))

                # Wait for t_movement
                current_time += mw_2_pi_half_time + delay_mw_readout + tau + delay_mw_readout + mw_2_three_pi_half_time + delay_mw_readout + t_movement

                # If 'echo' is checked, perform another z rotation. The phases should add up (or cancel)
                if self.settings['echo'] == 'subtract':
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time,
                                                mw_1_pi_half_time))
                    if tau > 0:
                        pulse_sequence.append(Pulse(microwave_channel_1_q,
                                                    current_time + mw_1_pi_half_time + delay_mw_readout,
                                                    tau))
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time + mw_1_pi_half_time + delay_mw_readout + tau + delay_mw_readout,
                                                mw_1_three_pi_half_time))

                    current_time += mw_1_pi_half_time + delay_mw_readout + tau + delay_mw_readout + mw_1_three_pi_half_time + delay_mw_readout

                elif self.settings['echo'] == 'add':
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time,
                                                mw_1_three_pi_half_time))
                    if tau > 0:
                        pulse_sequence.append(Pulse(microwave_channel_1_q,
                                                    current_time + mw_1_three_pi_half_time + delay_mw_readout,
                                                    tau))
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time + mw_1_three_pi_half_time + delay_mw_readout + tau + delay_mw_readout,
                                                mw_1_pi_half_time))

                    current_time += mw_1_three_pi_half_time + delay_mw_readout + tau + delay_mw_readout + mw_1_pi_half_time + delay_mw_readout


                # Wait for t_movement and rotate nuclear spin to Z axis for readout
                pulse_sequence.append(
                    Pulse('rf_switch',
                          current_time + rf_settle,
                          rf_pi_half_time))

                current_time += rf_settle + rf_pi_half_time + rf_settle

                if part == 1:

                    for m in range(self.settings['read_out']['repetitive_readout']['m']):
                        if m == 0:
                            pulse_sequence.append(Pulse(microwave_channel_2,
                                                        current_time,
                                                        mw_2_pi_time))
                            pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                        nv_reset_time_short))
                            pulse_sequence.append(
                                Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                                      meas_time))
                            current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                        else:
                            pulse_sequence.append(Pulse(microwave_channel,
                                                        current_time,
                                                        mw_1_pi_time))
                            pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                        nv_reset_time_short))
                            pulse_sequence.append(
                                Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                                      meas_time))
                            current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

                    # Update start time for the part 2 of the sequence (where we look at the other nuclear spin population)
                    start_time = current_time

                elif part == 0:
                    for m in range(self.settings['read_out']['repetitive_readout']['m']):
                        if m == 0:
                            pulse_sequence.append(Pulse(microwave_channel,
                                                        current_time,
                                                        mw_1_pi_time))
                            pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                        nv_reset_time_short))
                            pulse_sequence.append(
                                Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                                      meas_time))
                            current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                        else:
                            pulse_sequence.append(Pulse(microwave_channel_2,
                                                        current_time,
                                                        mw_2_pi_time))
                            pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                        nv_reset_time_short))
                            pulse_sequence.append(
                                Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                                      meas_time))
                            current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

                    start_time = current_time

            pulse_sequence.append(Pulse(rf_channel,
                                        0,
                                        start_time))
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time


class TransportDqmaDetuningFringes(TransportDqma):
    """
    With a finite detuning of the RF from the nuclear spin frequency, we expect to see oscillations for long movmement times
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns'),
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_three_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_time_i', 5000, float, 'NOT USED: time for RF pulse to drive nuclear spin in ns (when rf_i is on)'),
            Parameter('rf_pi_half_time_i', 5000, float, 'NOT USED: time for RF pulse to drive nuclear spin in ns (when rf_i is on)'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns'),
            Parameter('rf_echo', False, bool, 'add RF pi pulse to echo out any detunings, turning this into a nuclear T2 measurement')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)'),
            Parameter('cooldown', [
                Parameter('timing', 'duty_cycle', ['constant', 'duty_cycle'],
                          'use a constant cooldown time or set it based on length of entire sequence'),
                Parameter('cooldown_time', 1, float,
                          'if timing=constant, the same cooldown time in ns will be added to each sequence; '
                          'if timing=duty_cycle, cooldown_time/tau = 1-duty_cycle')])
        ]),
        Parameter('t_sense', 0, float, 'momvement/wait time in ns before reading state of the nuclear spin'),
        Parameter('dynamic_decoupling', [
            Parameter('sequence', 'none', ['none', 'echo', 'cpmg', 'xy'], 'dynamic decoupling sequence'),
            Parameter('k', 1, int, 'number of pulses in the decoupling sequence; this only applies for CPMG-k and XY-k')
            ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int, 'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time_long', 5000, int, 'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 100, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 8, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 400, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 250, float, 'short wait time between each readout laser pulse')
            ])
        ])
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator,
                    'afg': AFG3022C, 'mw_gen_2': MicrowaveGenerator2, 'commander': Commander}

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

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        microwave_channel, microwave_channel_2, microwave_channel_q, microwave_channel_2_q, rf_channel = \
            'microwave_i', 'microwave_i_2', 'microwave_q', 'microwave_q_2', 'rf_i'

        mw_1_pi_time, mw_2_pi_time = [self.settings['mw_pulses'][key] for key in ['mw_1_pi_time', 'mw_2_pi_time']]

        meas_time, nv_reset_time_long, laser_off_time, delay_mw_readout, delay_readout = \
            [self.settings['read_out'][key] for key in
             ['meas_time', 'nv_reset_time_long', 'laser_off_time', 'delay_mw_readout', 'delay_readout']]

        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']

        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']

        rf_pi_time, rf_pi_half_time, rf_three_pi_half_time, \
        rf_pi_time_i, rf_pi_half_time_i, rf_settle = \
            [self.settings['rf_pulses'][key] for key in ['rf_pi_time', 'rf_pi_half_time', 'rf_three_pi_half_time',
                                                         'rf_pi_time_i', 'rf_pi_half_time_i', 'rf_settle']]

        t_sense = self.settings['t_sense']

        pulse_sequences = []
        for tau in tau_list:
            t_movement = tau

            if self.settings['tau_times']['cooldown']['timing'] == 'constant':
                cooldown_time = self.settings['tau_times']['cooldown']['cooldown_time']
            elif self.settings['tau_times']['cooldown']['timing'] == 'duty_cycle':
                duty_cycle = self.settings['tau_times']['cooldown']['cooldown_time']
                cooldown_time = (rf_pi_time*self.settings['polarization_iterations']+rf_pi_half_time*2)*(1/duty_cycle - 1)
                cooldown_time = int(int(cooldown_time/100)*100)

            pulse_sequence = []
            start_time = 0

            for part in range(2):
                # Initialize nuclear and electronic spin at the beginning of each sequence
                pulse_sequence.append(Pulse('laser', start_time + laser_off_time + cooldown_time, nv_reset_time_long))
                current_time = start_time + laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
                for iteration in range(self.settings['polarization_iterations']):
                    pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                    if rf_pi_time > 0:
                        pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time, rf_pi_time))
                        current_time = current_time + mw_1_pi_time + rf_pi_time + rf_settle
                    for i in range(self.settings['initialization_laser']['n']):
                        pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                        current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

                # Flip e-spin to m_s=1 state so that we can drive RF
                pulse_sequence.append(Pulse(microwave_channel_2, current_time + laser_off_time, mw_2_pi_time))

                # Prepare nuclear spin in superposition
                pulse_sequence.append(Pulse('rf_switch', current_time + laser_off_time + mw_2_pi_time + delay_mw_readout, rf_pi_half_time))
                current_time += laser_off_time + mw_2_pi_time + delay_mw_readout + rf_pi_half_time + rf_settle
                # Access nuclear memory; choose to accumulate conditional on which nuclear state

                channel_odd, channel_even = microwave_channel_2_q, microwave_channel_q
                pi_time_odd, pi_time_even = mw_2_pi_time, mw_1_pi_time

                if self.settings['dynamic_decoupling']['sequence'] == 'none':

                    if not self.settings['rf_pulses']['rf_echo']:
                        print('Echo disabled')
                        pulse_sequence.append(Pulse(channel_odd,
                                                    current_time,
                                                    pi_time_odd))

                        # Wait for t_sense and disentangle from memory
                        pulse_sequence.append(Pulse(channel_odd,
                                                    current_time + pi_time_odd + t_sense,
                                                    pi_time_odd))
                    else:
                        print('Echo mode, disabling entangling/disentangling gates')

                    current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout

                elif self.settings['dynamic_decoupling']['sequence'] == 'echo':
                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time,
                                                pi_time_odd))

                    # Wait for t_sense and disentangle from memory
                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time + pi_time_odd + t_sense,
                                                pi_time_odd))

                    current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout

                    pulse_sequence.append(Pulse(channel_even,
                                                current_time,
                                                pi_time_even))

                    # Wait for t_sense and disentangle from memory
                    pulse_sequence.append(Pulse(channel_even,
                                                current_time + pi_time_even + t_sense,
                                                pi_time_even))

                    current_time += pi_time_even + t_sense + pi_time_even + delay_mw_readout
                if self.settings['rf_pulses']['rf_echo']:
                    pulse_sequence.append(Pulse('rf_switch', int(int((current_time + t_movement/2 - int(rf_pi_time / 2))/2)*2), rf_pi_time))

                # Wait for t_movement and rotate nuclear spin to Z axis for readout
                if part == 0:
                    pulse_sequence.append(Pulse('rf_switch', current_time + t_movement, rf_pi_half_time))
                    current_time += t_movement + rf_pi_half_time + rf_settle
                elif part == 1:
                    #pulse_sequence.append(Pulse('rf_switch', current_time + t_movement, rf_three_pi_half_time))
                    pulse_sequence.append(Pulse('rf_switch', current_time + t_movement, rf_pi_half_time))
                    pulse_sequence.append(Pulse(rf_channel, current_time + t_movement - rf_settle, rf_pi_half_time + rf_settle*2))
                    #current_time += t_movement + rf_three_pi_half_time + rf_settle
                    current_time += t_movement + rf_pi_half_time + rf_settle

                for m in range(self.settings['read_out']['repetitive_readout']['m']):
                    if m == 0:
                        pulse_sequence.append(Pulse(microwave_channel,
                                                    current_time,
                                                    mw_1_pi_time))
                        pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                    nv_reset_time_short))
                        pulse_sequence.append(
                            Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                                  meas_time))
                        current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                    else:
                        pulse_sequence.append(Pulse(microwave_channel_2,
                                                    current_time,
                                                    mw_2_pi_time))
                        pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                    nv_reset_time_short))
                        pulse_sequence.append(
                            Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                                  meas_time))
                        current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

                start_time = current_time

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _function(self):
        super(TransportDqmaDetuningFringes, self)._function()


class TransportDqmaMovement(RabiPolarized):
    """
    More careful timing of the whole pulse sequence to match the timing of FieldProfile
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns'),
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_three_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns'),
            Parameter('rf_settle_entangle', 2000, float, 'settling time after RF pulse in ns, for entangling CNOT pulse'),
            Parameter('final_pi_half_phase', [
                Parameter('modulation', 'fixed', ['fixed', 'sweep'], 'phase modulation of final pi-half pulse'),
                Parameter('value', 0, float, 'modulation set to fixed: this is the phase in radians of the final pi/2 pulse rel. to the first; modulation set to '
                                             'sweep: rel. phase of final pi/2 pulse will be set as phi = tau*freq, where freq in Hz is the value entered here.')
            ])
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('movement_timing', [
            Parameter('total_time', 600000, float, 'length of each part of the pulse sequence (there are two parts)'),
            Parameter('trig_duration', 200, float, 'length of trigger pulse in ns to send to function generator controlling attocube'),
            Parameter('settle_time', 2000, float, 'time in ns to wait between each pulse sequence, use the same setting as in FieldProfilePulsed'),
            Parameter('movement_end', 500000, float, 'time in ns of the end of the movement, as determined from FieldProfile script; '
                                                     'the pulse sequence will be timed such that the readout pulses happen after this time'),
            Parameter('rf_echo', [
                Parameter('enable', False, bool, 'apply RF pi pulse to echo out phase accumulation on nucleus from field changes during movement'),
                Parameter('echo_time', 250000, float, 'pi pulse will be centered at this time')
            ])
        ]),
        Parameter('dynamic_decoupling', [
            Parameter('sequence', 'none', ['none', 'echo', 'cpmg', 'xy'], 'dynamic decoupling sequence'),
            Parameter('k', 1, int, 'number of pulses in the decoupling sequence; this only applies for CPMG-k and XY-k')
            ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int, 'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time_long', 5000, int, 'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 100, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 8, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 400, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 250, float, 'short wait time between each readout laser pulse')
            ])
        ])
    ]

    _SCRIPTS = {'find_nv': FindNV, 'esr': ESR}

    def _configure_instruments_for_tau(self, tau):
        num_tau = len(self.data['tau']) - 4  # Actual number of tau values
        norm1 = np.average(np.average(self.data['counts'][num_tau:num_tau + 2], axis=1))
        norm2 = np.average(np.average(self.data['counts'][num_tau + 2:num_tau + 4], axis=1))

        if self.settings['rf_pulses']['final_pi_half_phase']['modulation'] == 'sweep':
            phase = float(np.pi*2*self.settings['rf_pulses']['final_pi_half_phase']['value'] * tau/1e9)
            self.instruments['afg']['instance'].align_phase()
            self.instruments['afg']['instance'].update({'ch1_phase': 0})
            self.instruments['afg']['instance'].update({'ch2_phase': phase})

        if tau == self.settings['tau_times']['max_time'] + 1000 or tau == self.settings['tau_times']['max_time'] + 2000:
            print('Norm meas for min fluor')
            self.log('Updated fluor bounds: %.2f, %.2f' % (norm1, norm2))
            self.instruments['afg']['instance'].update({'ch2_phase': np.pi})
        elif tau == self.settings['tau_times']['max_time'] + 3000 or tau == self.settings['tau_times']['max_time'] + 4000:
            print('Norm meas for max fluor')
            self.log('Updated fluor bounds: %.2f, %.2f' % (norm1, norm2))
            self.instruments['afg']['instance'].update({'ch2_phase': 0})
        else:
            self.instruments['afg']['instance'].update({'ch2_phase': 0})

        time.sleep(0.01)


    def _configure_instruments_start_of_script(self):
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_1_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_1_frequency']})

        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_2_power']})
        self.instruments['mw_gen_2']['instance'].update({'frequency': self.settings['mw_pulses']['mw_2_frequency']})

        self.instruments['afg']['instance'].update({'ch1_invert_polarity': False})
        self.instruments['afg']['instance'].update({'ch2_invert_polarity': False})

        for channel in [1, 2]:
            self.instruments['afg']['instance'].update({'ch%i_amplitude' % channel: float(self.dbm_to_vpp(self.settings['rf_pulses']['rf_power']))})
            self.instruments['afg']['instance'].update({'ch%i_frequency' % channel: float(self.settings['rf_pulses']['rf_frequency'])})
            self.instruments['afg']['instance'].update({'ch%i_function' % channel: 'Sine'})
            #self.instruments['afg']['instance'].update({'ch%i_phase' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_offset' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_enable' % channel: True})

        self.instruments['afg']['instance'].align_phase()
        self.instruments['afg']['instance'].update({'ch1_phase': 0})

        if self.settings['rf_pulses']['final_pi_half_phase']['modulation'] == 'fixed':
            ch2_phase = float(self.settings['rf_pulses']['final_pi_half_phase']['value'])
            self.instruments['afg']['instance'].align_phase()
            self.instruments['afg']['instance'].update({'ch2_phase': ch2_phase})
        else:
            self.instruments['afg']['instance'].update({'ch2_phase': 0})


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

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])
        max_tau = self.settings['tau_times']['max_time']
        tau_list = np.concatenate((tau_list, np.array([1000, 2000, 3000, 4000]) + max_tau))
        print('Combined tau list')
        print(tau_list)


        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        microwave_channel, microwave_channel_2, microwave_channel_q, microwave_channel_2_q, rf_channel = \
            'microwave_i', 'microwave_i_2', 'microwave_q', 'microwave_q_2', 'rf_i'

        mw_1_pi_time, mw_2_pi_time = [self.settings['mw_pulses'][key] for key in ['mw_1_pi_time', 'mw_2_pi_time']]

        meas_time, nv_reset_time_long, laser_off_time, delay_mw_readout, delay_readout = \
            [self.settings['read_out'][key] for key in
             ['meas_time', 'nv_reset_time_long', 'laser_off_time', 'delay_mw_readout', 'delay_readout']]

        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']

        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']

        rf_pi_time, rf_pi_half_time, rf_three_pi_half_time, rf_settle, rf_settle_entangle = \
            [self.settings['rf_pulses'][key] for key in ['rf_pi_time', 'rf_pi_half_time', 'rf_three_pi_half_time', 'rf_settle', 'rf_settle_entangle']]

        total_time = self.settings['movement_timing']['total_time']
        trig_duration = self.settings['movement_timing']['trig_duration']

        pulse_sequences = []
        for tau in tau_list:
            t_sense = tau

            # Cooldown time set to 0 because the timing of the pulses need to be fixed, based on the field profile measurement
            # Heating shouldn't be an issue since the pulse sequence is expected to be >500 us anyway
            cooldown_time = 0

            pulse_sequence = []
            start_time = 0 + int(total_time)*0
            current_time = start_time

            # Send a pulse that doesn't do anything for the total movement time, to set the durations of part 0 and 1 of the pulse sequence
            pulse_sequence.append(Pulse('spacer', total_time-200, 200))
            current_time += self.settings['movement_timing']['settle_time']
            pulse_sequence.append(Pulse('atto_trig', current_time, trig_duration))

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', current_time + laser_off_time + cooldown_time, nv_reset_time_long))

            final_readout_time = current_time + laser_off_time + cooldown_time + nv_reset_time_long - meas_time - 500
            #for m in range(self.settings['read_out']['repetitive_readout']['m']):
                #pulse_sequence.append(Pulse('apd_readout', final_readout_time, meas_time))
                #final_readout_time -= (meas_time + delay_readout)

            # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            current_time += laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time
            self.log('First nuc pol. MW pulse at %i ns' % current_time)

            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

            # Flip e-spin to m_s=1 state so that we can drive RF
            pulse_sequence.append(Pulse(microwave_channel_2, current_time + laser_off_time, mw_2_pi_time))

            # Prepare nuclear spin in superposition
            pulse_sequence.append(Pulse('rf_switch', current_time + laser_off_time + mw_2_pi_time + delay_mw_readout, rf_pi_half_time))
            print('Center of first pi/2 pulse: %i'% (current_time + laser_off_time + mw_2_pi_time + delay_mw_readout+rf_pi_half_time/2))
            current_time += laser_off_time + mw_2_pi_time + delay_mw_readout + rf_pi_half_time + rf_settle_entangle

            channel_odd, channel_even = microwave_channel_2_q, microwave_channel_q

            pi_time_odd, pi_time_even = mw_2_pi_time, mw_1_pi_time

            # Changed channels on 20230123, reverted 20230124 11 min after midnight
            #channel_odd, channel_even = microwave_channel, microwave_channel_2
            #pi_time_odd, pi_time_even = mw_1_pi_time, mw_2_pi_time

            # Access nuclear memory; choose to accumulate conditional on which nuclear state

            if self.settings['dynamic_decoupling']['sequence'] == 'none':
                if tau <= max_tau:
                    print('curr time:%i'%current_time)
                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time,
                                                pi_time_odd))

                    # Wait for t_sense and disentangle from memory
                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time + pi_time_odd + t_sense,
                                                pi_time_odd))
                    self.log('Disentangling MW pulse at %i ns' % (current_time + pi_time_odd + t_sense))
                else:
                    print('tau > max tau')

                current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout

            elif self.settings['dynamic_decoupling']['sequence'] == 'echo':
                pulse_sequence.append(Pulse(channel_odd,
                                            current_time,
                                            pi_time_odd))

                # Wait for t_sense and disentangle from memory
                pulse_sequence.append(Pulse(channel_odd,
                                            current_time + pi_time_odd + t_sense,
                                            pi_time_odd))

                current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout

                pulse_sequence.append(Pulse(channel_even,
                                            current_time,
                                            pi_time_even))

                # Wait for t_sense and disentangle from memory
                pulse_sequence.append(Pulse(channel_even,
                                            current_time + pi_time_even + t_sense,
                                            pi_time_even))

                current_time += pi_time_even + t_sense + pi_time_even + delay_mw_readout

            if self.settings['movement_timing']['rf_echo']['enable']:
                # Add RF echo pulse, time it with reference to the center of the pulse, also round it to an even number to avoid PB issues
                echo_time = int(total_time)*0 + self.settings['movement_timing']['rf_echo']['echo_time']
                pulse_sequence.append(Pulse('rf_switch', echo_time - int(rf_pi_time/2/2)*2, rf_pi_time))

            # Move for a certain time and rotate nuclear spin to Z axis for readout
            t_movement_end = int(total_time)*0 + self.settings['movement_timing']['movement_end']

            pulse_sequence.append(Pulse('rf_switch', t_movement_end - (rf_pi_half_time + rf_settle), rf_pi_half_time))
            pulse_sequence.append(Pulse(rf_channel, t_movement_end - (rf_pi_half_time + rf_settle) - rf_settle/2, rf_pi_half_time + rf_settle))
            print('Center of second pi/2 pulse: %i' % (t_movement_end - (rf_pi_half_time + rf_settle) + rf_pi_half_time / 2))

            current_time = t_movement_end
            self.log('First readout MW pulse at %i ns' % current_time)

            channel_readout_odd, channel_readout_even = microwave_channel, microwave_channel_2
            pi_time_odd, pi_time_even = mw_1_pi_time, mw_2_pi_time

            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                if m == 0:
                    pulse_sequence.append(Pulse(channel_readout_odd,
                                                current_time,
                                                pi_time_odd))
                    pulse_sequence.append(Pulse('laser', current_time + pi_time_odd + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + pi_time_odd + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + pi_time_odd + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                else:
                    pulse_sequence.append(Pulse(channel_readout_even,
                                                current_time,
                                                pi_time_even))
                    pulse_sequence.append(Pulse('laser', current_time + pi_time_even + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + pi_time_even + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + pi_time_even + delay_mw_readout + nv_reset_time_short + laser_off_time_short

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _add_rf_switch_to_sequences(self, pulse_sequences):
        return pulse_sequences

    def _configure_instruments_for_tau_OBSOLETE(self, tau):
        print(tau)
        if tau > self.settings['tau_times']['max_time']:
            self.instruments['afg']['instance'].align_phase()
            self.instruments['afg']['instance'].update({'ch1_phase': 0})
            self.instruments['afg']['instance'].update({'ch2_phase': np.pi})
            print('set phase to pi')
        else:
            self.instruments['afg']['instance'].align_phase()
            self.instruments['afg']['instance'].update({'ch1_phase': 0})
            self.instruments['afg']['instance'].update({'ch2_phase': 0})
            print('set phase to 0')
        time.sleep(0.01)


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

        if 'fits' in data.keys() and data['fits'] is not None and 1 == 2:
            counts = data['counts'][:, 1] / data['counts'][:, 0]
            tau = data['tau']
            fits = data['fits'] # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            pi_time = (np.pi - fits[2])/fits[1]
            pi_half_time = (np.pi/2 - fits[2])/fits[1]
            three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
            rabi_freq = 1000*fits[1]/(2*np.pi)
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            if 'counts' in data.keys():
                num_daq_reads = self.settings['read_out']['repetitive_readout']['m']
                num_tau = len(self.data['tau'])-4 # Actual number of tau values, recall that we made some 'fake' taus to tell the program to run different pulse sequences
                x_data = self.data['tau'][:num_tau]
                avg_counts_0 = np.transpose(np.array([np.average(self.data['counts'][:num_tau, :num_daq_reads], axis=1)]))
                #avg_counts_1 = np.transpose(np.array([np.average(self.data['counts'][num_tau:, :num_daq_reads], axis=1)]))
                plot_1d_simple_timetrace_ns(axislist[0], x_data, [avg_counts_0])
                #plot_1d_simple_timetrace_ns(axislist[0], data['tau'],[avg_counts_1])
                plot_pulses(axislist[1], self.pulse_sequences[self.sequence_index])
            axislist[0].set_title('Coherent Transport w/ Direct Quantum Memory Access')
            axislist[0].legend(labels=('Nuclear State 0 (avg readout)', 'Nuclear State 1 (avg readout)'), fontsize=8)

    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''
        #        if self.scripts['find_nv'].is_running:
        #            self.scripts['find_nv']._update_plot(axes_list)
        #        else:

        num_daq_reads = self.settings['read_out']['repetitive_readout']['m']
        num_tau = len(self.data['tau'])-4  # Actual number of tau values
        x_data = self.data['tau'][:num_tau]
        avg_counts_0 = np.transpose(np.array([np.average(self.data['counts'][:num_tau, :num_daq_reads], axis=1)]))
        #avg_counts_1 = np.transpose(np.array([np.average(self.data['counts'][num_tau:, :num_daq_reads], axis=1)]))

        axis1 = axes_list[0]
        if not self.data['counts'] == []:
            update_1d_simple(axis1, x_data, [avg_counts_0])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])



class TransportDqmaMovementSingle(RabiPolarized):
    """
    More careful timing of the whole pulse sequence to match the timing of FieldProfile
    No reference measurement
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns'),
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_three_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns'),
            #Parameter('invert_final_pi_half_pulse', False, bool, 'if true, final RF pi/2 pulse will be 180 deg out of phase w/ first one'),
            Parameter('final_pi_half_pulse_phase_vary_rate', 1e6, float, 'phase of final pi/2 pulse will be set as phi = tau*freq, where freq in Hz is the value entered here.')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('movement_timing', [
            Parameter('total_time', 600000, float, 'length of each part of the pulse sequence (there are two parts)'),
            Parameter('trig_duration', 200, float, 'length of trigger pulse in ns to send to function generator controlling attocube'),
            Parameter('settle_time', 2000, float, 'time in ns to wait between each pulse sequence, use the same setting as in FieldProfilePulsed'),
            Parameter('movement_end', 500000, float, 'time in ns of the end of the movement, as determined from FieldProfile script; '
                                                     'the pulse sequence will be timed such that the readout pulses happen after this time'),
            Parameter('rf_echo', [
                Parameter('enable', False, bool, 'apply RF pi pulse to echo out phase accumulation on nucleus from field changes during movement'),
                Parameter('echo_time', 250000, float, 'pi pulse will be centered at this time')
            ])
        ]),
        Parameter('dynamic_decoupling', [
            Parameter('sequence', 'none', ['none', 'echo', 'cpmg', 'xy'], 'dynamic decoupling sequence'),
            Parameter('k', 1, int, 'number of pulses in the decoupling sequence; this only applies for CPMG-k and XY-k')
            ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int, 'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time_long', 5000, int, 'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 100, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 8, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 400, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 250, float, 'short wait time between each readout laser pulse')
            ])
        ])
    ]

    _SCRIPTS = {'find_nv': FindNV, 'esr': ESR}

    def _configure_instruments_for_tau(self, tau):
        phase = float(np.pi*2*self.settings['rf_pulses']['final_pi_half_pulse_phase_vary_rate'] * tau/1e9)
        self.instruments['afg']['instance'].align_phase()
        self.instruments['afg']['instance'].update({'ch1_phase': 0})
        self.instruments['afg']['instance'].update({'ch2_phase': phase})


    def _create_pulse_sequences_two_part(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        microwave_channel, microwave_channel_2, microwave_channel_q, microwave_channel_2_q, rf_channel = \
            'microwave_i', 'microwave_i_2', 'microwave_q', 'microwave_q_2', 'rf_i'

        mw_1_pi_time, mw_2_pi_time = [self.settings['mw_pulses'][key] for key in ['mw_1_pi_time', 'mw_2_pi_time']]

        meas_time, nv_reset_time_long, laser_off_time, delay_mw_readout, delay_readout = \
            [self.settings['read_out'][key] for key in
             ['meas_time', 'nv_reset_time_long', 'laser_off_time', 'delay_mw_readout', 'delay_readout']]

        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']

        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']

        rf_pi_time, rf_pi_half_time, rf_three_pi_half_time, rf_settle = \
            [self.settings['rf_pulses'][key] for key in ['rf_pi_time', 'rf_pi_half_time', 'rf_three_pi_half_time', 'rf_settle']]


        total_time = self.settings['movement_timing']['total_time']
        trig_duration = self.settings['movement_timing']['trig_duration']

        pulse_sequences = []
        for tau in tau_list:
            t_sense = tau

            # Cooldown time set to 0 because the timing of the pulses need to be fixed, based on the field profile measurement
            # Heating shouldn't be an issue since the pulse sequence is expected to be >500 us anyway
            cooldown_time = 0

            pulse_sequence = []
            for part in range(2):
                start_time = 0 + int(total_time * part)
                current_time = start_time

                # Send a pulse that doesn't do anything for the total movement time, to set the durations of part 0 and 1 of the pulse sequence
                pulse_sequence.append(Pulse('spacer', start_time, total_time))
                current_time += self.settings['movement_timing']['settle_time']
                pulse_sequence.append(Pulse('atto_trig', current_time, trig_duration))

                # Initialize nuclear and electronic spin at the beginning of each sequence
                pulse_sequence.append(Pulse('laser', current_time + laser_off_time + cooldown_time, nv_reset_time_long))
                # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
                current_time += laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time

                for iteration in range(self.settings['polarization_iterations']):
                    pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                    if rf_pi_time > 0:
                        pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time, rf_pi_time))
                        current_time = current_time + mw_1_pi_time + rf_pi_time + rf_settle
                    for i in range(self.settings['initialization_laser']['n']):
                        pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                        current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

                # Flip e-spin to m_s=1 state so that we can drive RF
                pulse_sequence.append(Pulse(microwave_channel_2, current_time + laser_off_time, mw_2_pi_time))

                # Prepare nuclear spin in superposition
                pulse_sequence.append(Pulse('rf_switch', current_time + laser_off_time + mw_2_pi_time + delay_mw_readout, rf_pi_half_time))
                current_time += laser_off_time + mw_2_pi_time + delay_mw_readout + rf_pi_half_time + rf_settle

                channel_odd, channel_even = microwave_channel_2_q, microwave_channel_q
                pi_time_odd, pi_time_even = mw_2_pi_time, mw_1_pi_time
                # Access nuclear memory; choose to accumulate conditional on which nuclear state

                if self.settings['dynamic_decoupling']['sequence'] == 'none':

                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time,
                                                pi_time_odd))

                    # Wait for t_sense and disentangle from memory
                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time + pi_time_odd + t_sense,
                                                pi_time_odd))

                    current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout

                elif self.settings['dynamic_decoupling']['sequence'] == 'echo':
                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time,
                                                pi_time_odd))

                    # Wait for t_sense and disentangle from memory
                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time + pi_time_odd + t_sense,
                                                pi_time_odd))

                    current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout

                    pulse_sequence.append(Pulse(channel_even,
                                                current_time,
                                                pi_time_even))

                    # Wait for t_sense and disentangle from memory
                    pulse_sequence.append(Pulse(channel_even,
                                                current_time + pi_time_even + t_sense,
                                                pi_time_even))

                    current_time += pi_time_even + t_sense + pi_time_even + delay_mw_readout

                if self.settings['movement_timing']['rf_echo']['enable']:
                    # Add RF echo pulse, time it with reference to the center of the pulse, also round it to an even number to avoid PB issues
                    echo_time = int(total_time * part) + self.settings['movement_timing']['rf_echo']['echo_time']
                    pulse_sequence.append(Pulse('rf_switch', echo_time - int(rf_pi_time/2/2)*2, rf_pi_time))

                # Move for a certain time and rotate nuclear spin to Z axis for readout
                # Part 1's 3pi/2 pulse ends precisely at t_movement_end; time part 0's pi/2 pulse such that its center coincides with that of part 1's 3pi/2
                t_movement_end = int(total_time * part) + self.settings['movement_timing']['movement_end']
                print(t_movement_end)
                if part == 0:
                    #pulse_sequence.append(Pulse('rf_switch', t_movement_end - int((rf_three_pi_half_time/2 + rf_pi_half_time/2 + rf_settle)/2)*2, rf_pi_half_time))
                    pulse_sequence.append(Pulse('rf_switch', t_movement_end - (rf_pi_half_time + rf_settle), rf_pi_half_time))
                elif part == 1:
                    #pulse_sequence.append(Pulse('rf_switch', t_movement_end - (rf_three_pi_half_time + rf_settle), rf_three_pi_half_time))
                    pulse_sequence.append(Pulse('rf_switch', t_movement_end - (rf_pi_half_time + rf_settle), rf_pi_half_time))
                    pulse_sequence.append(Pulse(rf_channel, t_movement_end - (rf_pi_half_time + rf_settle) - rf_settle/2, rf_pi_half_time + rf_settle))

                current_time = t_movement_end

                for m in range(self.settings['read_out']['repetitive_readout']['m']):
                    if m == 0:
                        pulse_sequence.append(Pulse(microwave_channel,
                                                    current_time,
                                                    mw_1_pi_time))
                        pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                    nv_reset_time_short))
                        pulse_sequence.append(
                            Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                                  meas_time))
                        current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                    else:
                        pulse_sequence.append(Pulse(microwave_channel_2,
                                                    current_time,
                                                    mw_2_pi_time))
                        pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                    nv_reset_time_short))
                        pulse_sequence.append(
                            Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                                  meas_time))
                        current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

                    #self.log('Total length of pulse sequence: %i' % (current_time-laser_off_time_short))
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

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

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        microwave_channel, microwave_channel_2, microwave_channel_q, microwave_channel_2_q, rf_channel = \
            'microwave_i', 'microwave_i_2', 'microwave_q', 'microwave_q_2', 'rf_i'

        mw_1_pi_time, mw_2_pi_time = [self.settings['mw_pulses'][key] for key in ['mw_1_pi_time', 'mw_2_pi_time']]

        meas_time, nv_reset_time_long, laser_off_time, delay_mw_readout, delay_readout = \
            [self.settings['read_out'][key] for key in
             ['meas_time', 'nv_reset_time_long', 'laser_off_time', 'delay_mw_readout', 'delay_readout']]

        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']

        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']

        rf_pi_time, rf_pi_half_time, rf_three_pi_half_time, rf_settle = \
            [self.settings['rf_pulses'][key] for key in ['rf_pi_time', 'rf_pi_half_time', 'rf_three_pi_half_time', 'rf_settle']]

        total_time = self.settings['movement_timing']['total_time']
        trig_duration = self.settings['movement_timing']['trig_duration']
        #part = int(self.settings['rf_pulses']['invert_final_pi_half_pulse'])
        part = 1

        pulse_sequences = []
        for tau in tau_list:
            t_sense = tau

            # Cooldown time set to 0 because the timing of the pulses need to be fixed, based on the field profile measurement
            # Heating shouldn't be an issue since the pulse sequence is expected to be >500 us anyway
            cooldown_time = 0

            pulse_sequence = []
            start_time = 0 + int(total_time)*0
            current_time = start_time

            # Send a pulse that doesn't do anything for the total movement time, to set the durations of part 0 and 1 of the pulse sequence
            pulse_sequence.append(Pulse('spacer', start_time, total_time))
            current_time += self.settings['movement_timing']['settle_time']
            pulse_sequence.append(Pulse('atto_trig', current_time, trig_duration))

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', current_time + laser_off_time + cooldown_time, nv_reset_time_long))
            # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            current_time += laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time
            self.log('First nuc pol. MW pulse at %i ns' % current_time)

            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

            # Flip e-spin to m_s=1 state so that we can drive RF
            pulse_sequence.append(Pulse(microwave_channel_2, current_time + laser_off_time, mw_2_pi_time))

            # Prepare nuclear spin in superposition
            pulse_sequence.append(Pulse('rf_switch', current_time + laser_off_time + mw_2_pi_time + delay_mw_readout, rf_pi_half_time))
            current_time += laser_off_time + mw_2_pi_time + delay_mw_readout + rf_pi_half_time + rf_settle

            channel_odd, channel_even = microwave_channel_2_q, microwave_channel_q
            pi_time_odd, pi_time_even = mw_2_pi_time, mw_1_pi_time
            # Access nuclear memory; choose to accumulate conditional on which nuclear state

            if self.settings['dynamic_decoupling']['sequence'] == 'none':

                pulse_sequence.append(Pulse(channel_odd,
                                            current_time,
                                            pi_time_odd))

                # Wait for t_sense and disentangle from memory
                pulse_sequence.append(Pulse(channel_odd,
                                            current_time + pi_time_odd + t_sense,
                                            pi_time_odd))
                self.log('Disentangling MW pulse at %i ns' % (current_time + pi_time_odd + t_sense))

                current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout

            elif self.settings['dynamic_decoupling']['sequence'] == 'echo':
                pulse_sequence.append(Pulse(channel_odd,
                                            current_time,
                                            pi_time_odd))

                # Wait for t_sense and disentangle from memory
                pulse_sequence.append(Pulse(channel_odd,
                                            current_time + pi_time_odd + t_sense,
                                            pi_time_odd))

                current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout

                pulse_sequence.append(Pulse(channel_even,
                                            current_time,
                                            pi_time_even))

                # Wait for t_sense and disentangle from memory
                pulse_sequence.append(Pulse(channel_even,
                                            current_time + pi_time_even + t_sense,
                                            pi_time_even))

                current_time += pi_time_even + t_sense + pi_time_even + delay_mw_readout

            if self.settings['movement_timing']['rf_echo']['enable']:
                # Add RF echo pulse, time it with reference to the center of the pulse, also round it to an even number to avoid PB issues
                echo_time = int(total_time)*0 + self.settings['movement_timing']['rf_echo']['echo_time']
                pulse_sequence.append(Pulse('rf_switch', echo_time - int(rf_pi_time/2/2)*2, rf_pi_time))

            # Move for a certain time and rotate nuclear spin to Z axis for readout
            # Part 1's 3pi/2 pulse ends precisely at t_movement_end; time part 0's pi/2 pulse such that its center coincides with that of part 1's 3pi/2
            t_movement_end = int(total_time)*0 + self.settings['movement_timing']['movement_end']
            print(t_movement_end)
            if part == 0:
                #pulse_sequence.append(Pulse('rf_switch', t_movement_end - int((rf_three_pi_half_time/2 + rf_pi_half_time/2 + rf_settle)/2)*2, rf_pi_half_time))
                pulse_sequence.append(Pulse('rf_switch', t_movement_end - (rf_pi_half_time + rf_settle), rf_pi_half_time))

            elif part == 1:
                #pulse_sequence.append(Pulse('rf_switch', t_movement_end - (rf_three_pi_half_time + rf_settle), rf_three_pi_half_time))
                pulse_sequence.append(Pulse('rf_switch', t_movement_end - (rf_pi_half_time + rf_settle), rf_pi_half_time))
                pulse_sequence.append(Pulse(rf_channel, t_movement_end - (rf_pi_half_time + rf_settle) - rf_settle/2, rf_pi_half_time + rf_settle))

            current_time = t_movement_end
            self.log('First readout MW pulse at %i ns' % current_time)
            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                if m == 0:
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time,
                                                mw_1_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_1_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + mw_1_pi_time + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + mw_1_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                else:
                    pulse_sequence.append(Pulse(microwave_channel_2,
                                                current_time,
                                                mw_2_pi_time))
                    pulse_sequence.append(Pulse('laser', current_time + mw_2_pi_time + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + mw_2_pi_time + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + mw_2_pi_time + delay_mw_readout + nv_reset_time_short + laser_off_time_short

                #self.log('Total length of pulse sequence: %i' % (current_time-laser_off_time_short))
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _add_rf_switch_to_sequences(self, pulse_sequences):
        return pulse_sequences

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

        if 'fits' in data.keys() and data['fits'] is not None and 1 == 2:
            counts = data['counts'][:, 1] / data['counts'][:, 0]
            tau = data['tau']
            fits = data['fits'] # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            pi_time = (np.pi - fits[2])/fits[1]
            pi_half_time = (np.pi/2 - fits[2])/fits[1]
            three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
            rabi_freq = 1000*fits[1]/(2*np.pi)
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            if 'counts' in data.keys():
                num_daq_reads = self.settings['read_out']['repetitive_readout']['m']
                avg_counts_1 = np.transpose(np.array([np.average(self.data['counts'][:, 0:num_daq_reads], axis=1)]))
                plot_1d_simple_timetrace_ns(axislist[0], data['tau'],[avg_counts_1])
                plot_pulses(axislist[1], self.pulse_sequences[self.sequence_index])
            axislist[0].set_title('Coherent Transport w/ Direct Quantum Memory Access')
            axislist[0].legend(labels=('Nuclear State 0 (avg readout)', 'Nuclear State 1 (avg readout)'), fontsize=8)

    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''
        #        if self.scripts['find_nv'].is_running:
        #            self.scripts['find_nv']._update_plot(axes_list)
        #        else:

        x_data = self.data['tau']


        num_daq_reads = self.settings['read_out']['repetitive_readout']['m']
        avg_counts_1 = np.transpose(np.array([np.average(self.data['counts'][:, 0:num_daq_reads], axis=1)]))

        axis1 = axes_list[0]
        if not self.data['counts'] == []:
            update_1d_simple(axis1, x_data, [avg_counts_1])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])


class TransportDqmaMovementRfPhase(RabiPolarized):
    """
    More careful timing of the whole pulse sequence to match the timing of FieldProfile
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns'),
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_three_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns'),
            Parameter('rf_settle_entangle', 2000, float, 'settling time after RF pulse in ns, for entangling CNOT pulse'),
            Parameter('final_pi_half_phase', [
                Parameter('modulation', 'fixed', ['fixed'], 'phase modulation of final pi-half pulse'),
                Parameter('value', 0, float, 'modulation set to fixed: this is the phase in radians of the final pi/2 pulse rel. to the first; modulation set to '
                                             'sweep: rel. phase of final pi/2 pulse will be set as phi = tau*freq, where freq in Hz is the value entered here.')
            ])
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('rf_phase', [
            Parameter('min_phase', 0, float, 'minimum phase in radians of final RF pi/2 pulse, rel. to first RF pi/2 pulse'),
            Parameter('max_phase', 6.28318531, float, 'maximum phase in radians of final RF pi/2 pulse, rel. to first RF pi/2 pulse'),
            Parameter('phase_points', 40, int, 'number of phase steps'),
        ]),
        Parameter('t_sense', 1000, float, 'sensing duration between entangling/disentangling CNOT gates in ns'),
        Parameter('movement_timing', [
            Parameter('total_time', 600000, float, 'length of each part of the pulse sequence (there are two parts)'),
            Parameter('trig_duration', 200, float, 'length of trigger pulse in ns to send to function generator controlling attocube'),
            Parameter('settle_time', 2000, float, 'time in ns to wait between each pulse sequence, use the same setting as in FieldProfilePulsed'),
            Parameter('movement_end', 500000, float, 'time in ns of the end of the movement, as determined from FieldProfile script; '
                                                     'the pulse sequence will be timed such that the readout pulses happen after this time'),
            Parameter('rf_echo', [
                Parameter('enable', False, bool, 'apply RF pi pulse to echo out phase accumulation on nucleus from field changes during movement'),
                Parameter('echo_time', 250000, float, 'pi pulse will be centered at this time')
            ])
        ]),
        Parameter('dynamic_decoupling', [
            Parameter('sequence', 'none', ['none', 'echo', 'cpmg', 'xy'], 'dynamic decoupling sequence'),
            Parameter('k', 1, int, 'number of pulses in the decoupling sequence; this only applies for CPMG-k and XY-k')
            ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int, 'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time_long', 5000, int, 'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 100, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 8, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 400, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 250, float, 'short wait time between each readout laser pulse')
            ])
        ])
    ]

    _SCRIPTS = {'find_nv': FindNV, 'esr': ESR}

    def _configure_instruments_for_tau(self, phase):
        num_tau = len(self.data['tau'])-4  # Actual number of tau values
        norm1 = np.average(np.average(self.data['counts'][num_tau:num_tau + 2], axis=1))
        norm2 = np.average(np.average(self.data['counts'][num_tau + 2:num_tau + 4], axis=1))

        if phase == self.settings['rf_phase']['max_phase'] + 3000 or phase == self.settings['rf_phase']['max_phase'] + 4000:
            print('Norm meas for min fluor')
            #print(norm1, norm2)
            self.log('Updated fluor bounds: %.2f, %.2f' % (norm1, norm2))
            self.instruments['afg']['instance'].update({'ch2_phase': np.pi})
            flag = 1
        elif phase == self.settings['rf_phase']['max_phase'] + 1000 or phase == self.settings['rf_phase']['max_phase'] + 2000:
            print('Norm meas for max fluor')
            self.log('Updated fluor bounds: %.2f, %.2f'%(norm1, norm2))
            self.instruments['afg']['instance'].update({'ch2_phase': 0})
        else:
            self.instruments['afg']['instance'].update({'ch2_phase': float(phase)})
        time.sleep(0.01)



    def _configure_instruments_start_of_script(self):
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_1_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_1_frequency']})

        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_2_power']})
        self.instruments['mw_gen_2']['instance'].update({'frequency': self.settings['mw_pulses']['mw_2_frequency']})

        self.instruments['afg']['instance'].update({'ch1_invert_polarity': False})
        self.instruments['afg']['instance'].update({'ch2_invert_polarity': False})

        for channel in [1, 2]:
            self.instruments['afg']['instance'].update({'ch%i_amplitude' % channel: float(self.dbm_to_vpp(self.settings['rf_pulses']['rf_power']))})
            self.instruments['afg']['instance'].update({'ch%i_frequency' % channel: float(self.settings['rf_pulses']['rf_frequency'])})
            self.instruments['afg']['instance'].update({'ch%i_function' % channel: 'Sine'})
            #self.instruments['afg']['instance'].update({'ch%i_phase' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_offset' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_enable' % channel: True})

        self.instruments['afg']['instance'].align_phase()
        self.instruments['afg']['instance'].update({'ch1_phase': 0})

        if self.settings['rf_pulses']['final_pi_half_phase']['modulation'] == 'fixed':
            ch2_phase = float(self.settings['rf_pulses']['final_pi_half_phase']['value'])
            self.instruments['afg']['instance'].align_phase()
            self.instruments['afg']['instance'].update({'ch2_phase': ch2_phase})
        else:
            self.instruments['afg']['instance'].update({'ch2_phase': 0})


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

        tau_list = np.linspace(self.settings['rf_phase']['min_phase'], self.settings['rf_phase']['max_phase'], self.settings['rf_phase']['phase_points'])
        max_tau = self.settings['rf_phase']['max_phase']
        tau_list = np.concatenate((tau_list, np.array([1000, 2000, 3000, 4000]) + max_tau))

        microwave_channel, microwave_channel_2, microwave_channel_q, microwave_channel_2_q, rf_channel = \
            'microwave_i', 'microwave_i_2', 'microwave_q', 'microwave_q_2', 'rf_i'

        mw_1_pi_time, mw_2_pi_time = [self.settings['mw_pulses'][key] for key in ['mw_1_pi_time', 'mw_2_pi_time']]

        meas_time, nv_reset_time_long, laser_off_time, delay_mw_readout, delay_readout = \
            [self.settings['read_out'][key] for key in
             ['meas_time', 'nv_reset_time_long', 'laser_off_time', 'delay_mw_readout', 'delay_readout']]

        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']

        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']

        rf_pi_time, rf_pi_half_time, rf_three_pi_half_time, rf_settle, rf_settle_entangle = \
            [self.settings['rf_pulses'][key] for key in ['rf_pi_time', 'rf_pi_half_time', 'rf_three_pi_half_time', 'rf_settle', 'rf_settle_entangle']]

        total_time = self.settings['movement_timing']['total_time']
        trig_duration = self.settings['movement_timing']['trig_duration']

        pulse_sequences = []
        for tau in tau_list:
            t_sense = self.settings['t_sense']

            # Cooldown time set to 0 because the timing of the pulses need to be fixed, based on the field profile measurement
            # Heating shouldn't be an issue since the pulse sequence is expected to be >500 us anyway
            cooldown_time = 0

            pulse_sequence = []
            start_time = 0 + int(total_time)*0
            current_time = start_time

            # Send a pulse that doesn't do anything for the total movement time, to set the durations of part 0 and 1 of the pulse sequence
            pulse_sequence.append(Pulse('spacer', total_time-200, 200))
            current_time += self.settings['movement_timing']['settle_time']
            pulse_sequence.append(Pulse('atto_trig', current_time, trig_duration))

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', current_time + laser_off_time + cooldown_time, nv_reset_time_long))

            final_readout_time = current_time + laser_off_time + cooldown_time + nv_reset_time_long - meas_time - 500
            #for m in range(self.settings['read_out']['repetitive_readout']['m']):
                #pulse_sequence.append(Pulse('apd_readout', final_readout_time, meas_time))
                #final_readout_time -= (meas_time + delay_readout)

            # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            current_time += laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time
            #self.log('First nuc pol. MW pulse at %i ns' % current_time)

            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

            # Flip e-spin to m_s=1 state so that we can drive RF
            pulse_sequence.append(Pulse(microwave_channel_2, current_time + laser_off_time, mw_2_pi_time))

            # Prepare nuclear spin in superposition
            pulse_sequence.append(Pulse('rf_switch', current_time + laser_off_time + mw_2_pi_time + delay_mw_readout, rf_pi_half_time))
            current_time += laser_off_time + mw_2_pi_time + delay_mw_readout + rf_pi_half_time + rf_settle_entangle

            channel_odd, channel_even = microwave_channel_2_q, microwave_channel_q
            pi_time_odd, pi_time_even = mw_2_pi_time, mw_1_pi_time

            # Changed channels on 20230123, reverted 20230124 11 min after midnight
            #channel_odd, channel_even = microwave_channel, microwave_channel_2
            #pi_time_odd, pi_time_even = mw_1_pi_time, mw_2_pi_time

            # Access nuclear memory; choose to accumulate conditional on which nuclear state

            if self.settings['dynamic_decoupling']['sequence'] == 'none':
                if tau <= max_tau:
                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time,
                                                pi_time_odd))
                    self.log('Entangling MW pulse at %i ns' % (current_time))
                    # Wait for t_sense and disentangle from memory
                    pulse_sequence.append(Pulse(channel_odd,
                                                current_time + pi_time_odd + t_sense,
                                                pi_time_odd))
                    #self.log('Disentangling MW pulse at %i ns' % (current_time + pi_time_odd + t_sense))

                current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout

            elif self.settings['dynamic_decoupling']['sequence'] == 'echo':
                pulse_sequence.append(Pulse(channel_odd,
                                            current_time,
                                            pi_time_odd))

                # Wait for t_sense and disentangle from memory
                pulse_sequence.append(Pulse(channel_odd,
                                            current_time + pi_time_odd + t_sense,
                                            pi_time_odd))

                current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout

                pulse_sequence.append(Pulse(channel_even,
                                            current_time,
                                            pi_time_even))

                # Wait for t_sense and disentangle from memory
                pulse_sequence.append(Pulse(channel_even,
                                            current_time + pi_time_even + t_sense,
                                            pi_time_even))

                current_time += pi_time_even + t_sense + pi_time_even + delay_mw_readout

            if self.settings['movement_timing']['rf_echo']['enable']:
                # Add RF echo pulse, time it with reference to the center of the pulse, also round it to an even number to avoid PB issues
                echo_time = int(total_time)*0 + self.settings['movement_timing']['rf_echo']['echo_time']
                pulse_sequence.append(Pulse('rf_switch', echo_time - int(rf_pi_time/2/2)*2, rf_pi_time))

            # Move for a certain time and rotate nuclear spin to Z axis for readout
            t_movement_end = int(total_time)*0 + self.settings['movement_timing']['movement_end']

            pulse_sequence.append(Pulse('rf_switch', t_movement_end - (rf_pi_half_time + rf_settle), rf_pi_half_time))
            pulse_sequence.append(Pulse(rf_channel, t_movement_end - (rf_pi_half_time + rf_settle) - rf_settle/2, rf_pi_half_time + rf_settle))

            current_time = t_movement_end
            #self.log('First readout MW pulse at %i ns' % current_time)

            channel_readout_odd, channel_readout_even = microwave_channel, microwave_channel_2
            pi_time_odd, pi_time_even = mw_1_pi_time, mw_2_pi_time

            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                if m == 0:
                    pulse_sequence.append(Pulse(channel_readout_odd,
                                                current_time,
                                                pi_time_odd))
                    pulse_sequence.append(Pulse('laser', current_time + pi_time_odd + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + pi_time_odd + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + pi_time_odd + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                else:
                    pulse_sequence.append(Pulse(channel_readout_even,
                                                current_time,
                                                pi_time_even))
                    pulse_sequence.append(Pulse('laser', current_time + pi_time_even + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + pi_time_even + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + pi_time_even + delay_mw_readout + nv_reset_time_short + laser_off_time_short

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time


    def _add_rf_switch_to_sequences(self, pulse_sequences):
        return pulse_sequences



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

        if 'fits' in data.keys() and data['fits'] is not None and 1 == 2:
            counts = data['counts'][:, 1] / data['counts'][:, 0]
            tau = data['tau']
            fits = data['fits'] # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            pi_time = (np.pi - fits[2])/fits[1]
            pi_half_time = (np.pi/2 - fits[2])/fits[1]
            three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
            rabi_freq = 1000*fits[1]/(2*np.pi)
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            if 'counts' in data.keys():
                num_daq_reads = self.settings['read_out']['repetitive_readout']['m']
                num_tau = len(self.data['tau'])-4 # Actual number of tau values, recall that we made some 'fake' taus to tell the program to run different pulse sequences
                x_data = self.data['tau'][:num_tau]
                avg_counts_0 = np.transpose(np.array([np.average(self.data['counts'][:num_tau, :num_daq_reads], axis=1)]))
                #avg_counts_1 = np.transpose(np.array([np.average(self.data['counts'][num_tau:, :num_daq_reads], axis=1)]))
                plot_1d_simple_timetrace_ns(axislist[0], x_data, [avg_counts_0])
                #plot_1d_simple_timetrace_ns(axislist[0], data['tau'],[avg_counts_1])
                plot_pulses(axislist[1], self.pulse_sequences[self.sequence_index])
            axislist[0].set_title('Coherent Transport w/ Direct Quantum Memory Access')
            axislist[0].legend(labels=('Nuclear State 0 (avg readout)', 'Nuclear State 1 (avg readout)'), fontsize=8)

    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''
        #        if self.scripts['find_nv'].is_running:
        #            self.scripts['find_nv']._update_plot(axes_list)
        #        else:

        num_daq_reads = self.settings['read_out']['repetitive_readout']['m']
        num_tau = len(self.data['tau'])-4  # Actual number of tau values
        x_data = self.data['tau'][:num_tau]
        avg_counts_0 = np.transpose(np.array([np.average(self.data['counts'][:num_tau, :num_daq_reads], axis=1)]))
        #avg_counts_1 = np.transpose(np.array([np.average(self.data['counts'][num_tau:, :num_daq_reads], axis=1)]))

        axis1 = axes_list[0]
        if not self.data['counts'] == []:
            update_1d_simple(axis1, x_data, [avg_counts_0])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])


class TransportDqmaMovementPiPulseTime(RabiPolarized):
    """
    Sweeps the timing of the decoupling RF pi pulse to echo out the phase accumulation from movement
    More careful timing of the whole pulse sequence to match the timing of FieldProfile
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_1_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_1_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 1, also the transition to be depopulated for initialization'),
            Parameter('mw_1_pi_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns'),
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_three_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns'),
            Parameter('rf_settle_entangle', 2000, float, 'settling time after RF pulse in ns, for entangling CNOT pulse'),
            Parameter('final_pi_half_phase', [
                Parameter('modulation', 'fixed', ['fixed', 'sweep'], 'phase modulation of final pi-half pulse'),
                Parameter('value', 0, float, 'modulation set to fixed: this is the phase in radians of the final pi/2 pulse rel. to the first; modulation set to '
                                             'sweep: rel. phase of final pi/2 pulse will be set as phi = tau*freq, where freq in Hz is the value entered here.')
            ])
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('tau_times', [
            Parameter('min_time', 100, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 2000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 500., float, 'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('t_sense', 1000, float, 'sensing duration between entangling/disentangling CNOT gates in ns, if t = -1, disable those gates entirely'),
        Parameter('movement_timing', [
            Parameter('total_time', 600000, float, 'length of each part of the pulse sequence (there are two parts)'),
            Parameter('trig_duration', 200, float, 'length of trigger pulse in ns to send to function generator controlling attocube'),
            Parameter('settle_time', 2000, float, 'time in ns to wait between each pulse sequence, use the same setting as in FieldProfilePulsed'),
            Parameter('movement_end', 500000, float, 'time in ns of the end of the movement, as determined from FieldProfile script; '
                                                     'the pulse sequence will be timed such that the readout pulses happen after this time')
        ]),
        Parameter('dynamic_decoupling', [
            Parameter('sequence', 'none', ['none', 'echo', 'cpmg', 'xy'], 'dynamic decoupling sequence'),
            Parameter('k', 1, int, 'number of pulses in the decoupling sequence; this only applies for CPMG-k and XY-k')
            ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int, 'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time_long', 5000, int, 'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 100, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('repetitive_readout', [
                Parameter('m', 8, int, 'number of repetitions'),
                Parameter('nv_reset_time_short', 400, float, 'short laser pulse to reset NV'),
                Parameter('laser_off_time_short', 250, float, 'short wait time between each readout laser pulse')
            ])
        ])
    ]

    _SCRIPTS = {'find_nv': FindNV, 'esr': ESR}

    def _configure_instruments_for_tau(self, tau):

        if self.settings['rf_pulses']['final_pi_half_phase']['modulation'] == 'sweep':
            phase = float(np.pi*2*self.settings['rf_pulses']['final_pi_half_phase']['value'] * tau/1e9)
            self.instruments['afg']['instance'].align_phase()
            self.instruments['afg']['instance'].update({'ch1_phase': 0})
            self.instruments['afg']['instance'].update({'ch2_phase': phase})
            time.sleep(0.01)

        elif self.settings['rf_pulses']['final_pi_half_phase']['modulation'] == 'fixed':
            if tau > self.settings['tau_times']['max_time']:
                phase = self.settings['rf_pulses']['final_pi_half_phase']['value']
                self.instruments['afg']['instance'].update({'ch2_phase': phase})
            else:
                self.instruments['afg']['instance'].update({'ch2_phase': 0})


    def _configure_instruments_start_of_script(self):
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_1_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_1_frequency']})

        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_2_power']})
        self.instruments['mw_gen_2']['instance'].update({'frequency': self.settings['mw_pulses']['mw_2_frequency']})

        self.instruments['afg']['instance'].update({'ch1_invert_polarity': False})
        self.instruments['afg']['instance'].update({'ch2_invert_polarity': False})

        for channel in [1, 2]:
            self.instruments['afg']['instance'].update({'ch%i_amplitude' % channel: float(self.dbm_to_vpp(self.settings['rf_pulses']['rf_power']))})
            self.instruments['afg']['instance'].update({'ch%i_frequency' % channel: float(self.settings['rf_pulses']['rf_frequency'])})
            self.instruments['afg']['instance'].update({'ch%i_function' % channel: 'Sine'})
            #self.instruments['afg']['instance'].update({'ch%i_phase' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_offset' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_enable' % channel: True})

        self.instruments['afg']['instance'].align_phase()
        self.instruments['afg']['instance'].update({'ch1_phase': 0})

        if self.settings['rf_pulses']['final_pi_half_phase']['modulation'] == 'fixed':
            ch2_phase = float(self.settings['rf_pulses']['final_pi_half_phase']['value'])
            self.instruments['afg']['instance'].align_phase()
            self.instruments['afg']['instance'].update({'ch2_phase': ch2_phase})
        else:
            self.instruments['afg']['instance'].update({'ch2_phase': 0})


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

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])
        max_tau = self.settings['tau_times']['max_time']
        tau_list = np.concatenate((tau_list, tau_list + max_tau))
        print('Combined tau list')
        print(tau_list)


        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        microwave_channel, microwave_channel_2, microwave_channel_q, microwave_channel_2_q, rf_channel = \
            'microwave_i', 'microwave_i_2', 'microwave_q', 'microwave_q_2', 'rf_i'

        mw_1_pi_time, mw_2_pi_time = [self.settings['mw_pulses'][key] for key in ['mw_1_pi_time', 'mw_2_pi_time']]

        meas_time, nv_reset_time_long, laser_off_time, delay_mw_readout, delay_readout = \
            [self.settings['read_out'][key] for key in
             ['meas_time', 'nv_reset_time_long', 'laser_off_time', 'delay_mw_readout', 'delay_readout']]

        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']

        laser_off_time_short = self.settings['read_out']['repetitive_readout']['laser_off_time_short']
        nv_reset_time_short = self.settings['read_out']['repetitive_readout']['nv_reset_time_short']

        rf_pi_time, rf_pi_half_time, rf_three_pi_half_time, rf_settle, rf_settle_entangle = \
            [self.settings['rf_pulses'][key] for key in ['rf_pi_time', 'rf_pi_half_time', 'rf_three_pi_half_time', 'rf_settle', 'rf_settle_entangle']]

        total_time = self.settings['movement_timing']['total_time']
        trig_duration = self.settings['movement_timing']['trig_duration']
        t_sense = self.settings['t_sense']

        pulse_sequences = []
        for tau in tau_list:
            if tau > max_tau:
                tau -= max_tau
            # Cooldown time set to 0 because the timing of the pulses need to be fixed, based on the field profile measurement
            # Heating shouldn't be an issue since the pulse sequence is expected to be >500 us anyway
            cooldown_time = 0

            pulse_sequence = []
            start_time = 0 + int(total_time)*0
            current_time = start_time

            # Send a pulse that doesn't do anything for the total movement time, to set the durations of part 0 and 1 of the pulse sequence
            pulse_sequence.append(Pulse('spacer', total_time-200, 200))
            current_time += self.settings['movement_timing']['settle_time']
            pulse_sequence.append(Pulse('atto_trig', current_time, trig_duration))

            # Initialize nuclear and electronic spin at the beginning of each sequence
            pulse_sequence.append(Pulse('laser', current_time + laser_off_time + cooldown_time, nv_reset_time_long))

            final_readout_time = current_time + laser_off_time + cooldown_time + nv_reset_time_long - meas_time - 500
            #for m in range(self.settings['read_out']['repetitive_readout']['m']):
                #pulse_sequence.append(Pulse('apd_readout', final_readout_time, meas_time))
                #final_readout_time -= (meas_time + delay_readout)

            # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            current_time += laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time
            #self.log('First nuc pol. MW pulse at %i ns' % current_time)

            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_1_pi_time))
                if rf_pi_time > 0:
                    pulse_sequence.append(Pulse('rf_switch', current_time + mw_1_pi_time, rf_pi_time))
                    current_time = current_time + mw_1_pi_time + rf_pi_time + rf_settle
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait

            # Flip e-spin to m_s=1 state so that we can drive RF
            pulse_sequence.append(Pulse(microwave_channel_2, current_time + laser_off_time, mw_2_pi_time))

            # Prepare nuclear spin in superposition
            pulse_sequence.append(Pulse('rf_switch', current_time + laser_off_time + mw_2_pi_time + delay_mw_readout, rf_pi_half_time))
            #print('Center of first pi/2 pulse: %i'% (current_time + laser_off_time + mw_2_pi_time + delay_mw_readout+rf_pi_half_time/2))
            current_time += laser_off_time + mw_2_pi_time + delay_mw_readout + rf_pi_half_time + rf_settle_entangle

            channel_odd, channel_even = microwave_channel_2_q, microwave_channel_q
            pi_time_odd, pi_time_even = mw_2_pi_time, mw_1_pi_time

            #channel_odd, channel_even = microwave_channel_2, microwave_channel
            #pi_time_odd, pi_time_even = mw_2_pi_time, mw_1_pi_time

            if tau <= max_tau and t_sense != -1:
                pulse_sequence.append(Pulse(channel_odd,
                                            current_time,
                                            pi_time_odd))
                self.log('Entangling MW pulse at %i ns' % (current_time))

                # Wait for t_sense and disentangle from memory
                pulse_sequence.append(Pulse(channel_odd,
                                            current_time + pi_time_odd + t_sense,
                                            pi_time_odd))
                #self.log('Disentangling MW pulse at %i ns' % (current_time + pi_time_odd + t_sense))


            # Add RF echo pulse, time it with reference to the center of the pulse, also round it to an even number to avoid PB issues
            echo_time = int(total_time)*0 + tau
            pulse_sequence.append(Pulse('rf_switch', echo_time - int(rf_pi_time/2/2)*2, rf_pi_time))

            # Move for a certain time and rotate nuclear spin to Z axis for readout
            t_movement_end = int(total_time)*0 + self.settings['movement_timing']['movement_end']

            pulse_sequence.append(Pulse('rf_switch', t_movement_end - (rf_pi_half_time + rf_settle), rf_pi_half_time))
            pulse_sequence.append(Pulse(rf_channel, t_movement_end - (rf_pi_half_time + rf_settle) - rf_settle/2, rf_pi_half_time + rf_settle))
            #print('Center of second pi/2 pulse: %i' % (t_movement_end - (rf_pi_half_time + rf_settle) + rf_pi_half_time / 2))

            current_time = t_movement_end
            #self.log('First readout MW pulse at %i ns' % current_time)

            channel_readout_odd, channel_readout_even = microwave_channel, microwave_channel_2
            pi_time_odd, pi_time_even = mw_1_pi_time, mw_2_pi_time

            for m in range(self.settings['read_out']['repetitive_readout']['m']):
                if m == 0:
                    pulse_sequence.append(Pulse(channel_readout_odd,
                                                current_time,
                                                pi_time_odd))
                    pulse_sequence.append(Pulse('laser', current_time + pi_time_odd + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + pi_time_odd + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + pi_time_odd + delay_mw_readout + nv_reset_time_short + laser_off_time_short
                else:
                    pulse_sequence.append(Pulse(channel_readout_even,
                                                current_time,
                                                pi_time_even))
                    pulse_sequence.append(Pulse('laser', current_time + pi_time_even + delay_mw_readout,
                                                nv_reset_time_short))
                    pulse_sequence.append(
                        Pulse('apd_readout', current_time + pi_time_even + delay_mw_readout + delay_readout,
                              meas_time))
                    current_time = current_time + pi_time_even + delay_mw_readout + nv_reset_time_short + laser_off_time_short

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _add_rf_switch_to_sequences(self, pulse_sequences):
        return pulse_sequences

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

        if 'fits' in data.keys() and data['fits'] is not None and 1 == 2:
            counts = data['counts'][:, 1] / data['counts'][:, 0]
            tau = data['tau']
            fits = data['fits'] # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            pi_time = (np.pi - fits[2])/fits[1]
            pi_half_time = (np.pi/2 - fits[2])/fits[1]
            three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
            rabi_freq = 1000*fits[1]/(2*np.pi)
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            if 'counts' in data.keys():
                num_daq_reads = self.settings['read_out']['repetitive_readout']['m']
                num_tau = int(len(self.data['tau']) / 2) # Actual number of tau values, recall that we made some 'fake' taus to tell the program to run different pulse sequences
                print(num_tau)
                x_data = self.data['tau'][:num_tau]
                avg_counts_0 = np.transpose(np.array([np.average(self.data['counts'][:num_tau, :num_daq_reads], axis=1)]))
                avg_counts_1 = np.transpose(np.array([np.average(self.data['counts'][num_tau:, :num_daq_reads], axis=1)]))
                plot_1d_simple_timetrace_ns(axislist[0], x_data, [avg_counts_0, avg_counts_1])
                #plot_1d_simple_timetrace_ns(axislist[0], data['tau'],[avg_counts_1])
                plot_pulses(axislist[1], self.pulse_sequences[self.sequence_index])
            axislist[0].set_title('Coherent Transport w/ Direct Quantum Memory Access')
            axislist[0].legend(labels=('Nuclear State 0 (avg readout)', 'Nuclear State 1 (avg readout)'), fontsize=8)

    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''
        #        if self.scripts['find_nv'].is_running:
        #            self.scripts['find_nv']._update_plot(axes_list)
        #        else:

        num_daq_reads = self.settings['read_out']['repetitive_readout']['m']
        num_tau = int(len(self.data['tau']) / 2)  # Actual number of tau values
        x_data = self.data['tau'][:num_tau]
        avg_counts_0 = np.transpose(np.array([np.average(self.data['counts'][:num_tau, :num_daq_reads], axis=1)]))
        avg_counts_1 = np.transpose(np.array([np.average(self.data['counts'][num_tau:, :num_daq_reads], axis=1)]))

        axis1 = axes_list[0]
        if not self.data['counts'] == []:
            update_1d_simple(axis1, x_data, [avg_counts_0,avg_counts_1])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])