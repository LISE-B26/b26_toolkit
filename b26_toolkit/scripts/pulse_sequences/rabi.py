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
from pylabcontrol.core import Parameter, Script
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_generic import PulsedExperimentGeneric
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, MicrowaveGenerator, Pulse, Commander
from b26_toolkit.data_processing.fit_functions import fit_rabi_decay, cose_with_decay


class Rabi(PulsedExperimentGeneric):
    """
    This script applies a microwave pulse at fixed power for varying durations to measure Rabi oscillations.
    Uses a double_init scheme
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', '+i', ['+i', '-i', '+q', '-q'], 'Channel to use for mw pulses')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 20, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 400, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 20., [2.5, 4., 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 340, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 1000000, int, 'number of averages'),
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'commander': Commander}
    #_SCRIPTS = {'find_nv': FindNV, 'esr': ESR}

    def _configure_instruments_start_of_script(self):
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

    def _function(self):
        # COMMENT_ME
        self.data['fits'] = None
        super(Rabi, self)._function()

        if 'counts' in self.data.keys() and 'tau' in self.data.keys():
            counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
            tau = self.data['tau']

            try:
                fits = fit_rabi_decay(tau, counts, variable_phase=True)
                self.data['fits'] = fits
            except:
                self.data['fits'] = None
                self.log('rabi fit failed')

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

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        #MM: update 15 to min_pulse_duration
        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]


        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = [Pulse('laser', laser_off_time + tau, nv_reset_time),
                              Pulse('apd_readout', laser_off_time + tau + delay_readout, meas_time)]

            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence.append(Pulse(microwave_channel,
                                            laser_off_time + tau + nv_reset_time + laser_off_time,
                                            tau))

            pulse_sequence.append(Pulse('laser',
                                        laser_off_time + tau + nv_reset_time + laser_off_time + tau + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        laser_off_time + tau + nv_reset_time + laser_off_time + tau + delay_mw_readout + delay_readout,
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

        if 'fits' in data.keys() and data['fits'] is not None:
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

class RabiResonant(Rabi):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', '+i', ['+i', '-i', '+q', '-q'], 'Channel to use for mw pulses')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 20, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 400, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 20., [2.5, 4., 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('red_on_time', 1000, int, 'time that red laser is on'),
            Parameter('red_off_time', 1000, int, 'time that red laser is off'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('measure_ref', True, bool, 'add sequence to measure ms=0 state as reference')
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
        pulse_sequences = []

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)

        #MM: update 15 to min_pulse_duration
        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]


        nv_reset_time = self.settings['read_out']['nv_reset_time']
        red_on_time = self.settings['read_out']['red_on_time']
        red_off_time = self.settings['read_out']['red_off_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = [Pulse('laser', red_off_time, nv_reset_time)]

            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence.append(Pulse(microwave_channel,
                                            red_off_time + nv_reset_time + laser_off_time,
                                            tau))

            pulse_sequence.append(Pulse('red_laser',
                                        red_off_time + nv_reset_time + laser_off_time + tau + delay_mw_readout,
                                        red_on_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        red_off_time + nv_reset_time + laser_off_time + tau + delay_mw_readout + delay_readout,
                                        meas_time))

            if self.settings['measure_ref']:
                end_first_seq = red_off_time + nv_reset_time + laser_off_time + tau + delay_mw_readout + red_on_time + red_off_time
                pulse_sequence += [
                    Pulse('laser', end_first_seq, nv_reset_time),
                    Pulse('red_laser', end_first_seq + nv_reset_time + laser_off_time, red_on_time),
                    Pulse('apd_readout', end_first_seq + nv_reset_time + laser_off_time + delay_readout, meas_time)
                ]

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

class RabiChoppedInit(PulsedExperimentGeneric):
    """
    This script applies a microwave pulse at fixed power for varying durations to measure Rabi oscillations.
    Uses a double_init scheme
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', '+i', ['+i', '-i', '+q', '-q'], 'Channel to use for mw pulses')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 20, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 400, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 20., [2.5, 4., 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 340, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('init_laser', [
            Parameter('pulse_duration', 80, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 80, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int, 'num of initialization laser pulses')
        ]),
        Parameter('num_averages', 1000000, int, 'number of averages'),
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'commander': Commander}
    #_SCRIPTS = {'find_nv': FindNV, 'esr': ESR}

    def _configure_instruments_start_of_script(self):
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

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

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        #MM: update 15 to min_pulse_duration
        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]


        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        init_pulse_time =  self.settings['init_laser']['pulse_duration']
        init_wait_time = self.settings['init_laser']['wait_duration']
        init_n = self.settings['init_laser']['n']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = []
            current_time = laser_off_time
            for i in range(init_n):
                pulse_sequence.append(Pulse('laser', current_time, init_pulse_time))
                current_time += init_pulse_time + init_wait_time

            current_time -= init_wait_time  # Remove the last, unnecessary wait time. If you want spacing between reset pulse(s) and MW, use laser_off_time

            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence.append(Pulse(microwave_channel, current_time + laser_off_time, tau))

            pulse_sequence.append(Pulse('laser', current_time + laser_off_time + tau + delay_mw_readout, nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout', current_time + laser_off_time + tau + delay_mw_readout + delay_readout,
                                        meas_time))
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

        super(RabiChoppedInit, self)._plot(axislist)
        axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'],
                                                                                     self.settings['mw_pulses']['mw_frequency'] * 1e-9))


class RabiDoublePi(Rabi):
    """
    Runs Rabi, but instead of sweeping the duration of a single MW pulse, sweeps the duration of two back to back MW
    pulses with same durations. Used to calibrate pulses such that two pi pulses bring the Bloch vector back to the
    original state
    """

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

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)

        #MM: update 15 to min_pulse_duration
        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]


        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = [Pulse('laser', laser_off_time + tau + 2*40, nv_reset_time),
                              Pulse('apd_readout', laser_off_time + tau + 2*40 + delay_readout, meas_time)]

            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence.append(Pulse(microwave_channel,
                                            laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time,
                                            tau))
                pulse_sequence.append(Pulse(microwave_channel,
                                            laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2*40,
                                            tau))

            pulse_sequence.append(Pulse('laser',
                                        laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time + tau + 2*40 + tau + 2*40 + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time + tau + 2*40 + tau + 2*40 + delay_mw_readout + delay_readout,
                                        meas_time))
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time


class PiPulseTrain(Rabi):
    """
    Runs a MW pi-pulse train. Measures the contrast as we vary the number of pi-pulses.
    """

    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.82e9, float, 'frequency of hyperfine transition'),
            Parameter('microwave_channel', '+i', ['+i', '-i', '+q', '-q',
                                                  '+i, +q', '+i, -i', '+q, -q',
                                                  '+i, -i, +q, -q',
                                                  '+i, +q, -i, -q'
                                                  ], 'Channel to use for mw pulses'),
            Parameter('tau_mw', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('spacing', 100, float, 'spacing in ns between consecutive pi-pulses')
        ]),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('n_start', 0, int, 'start num of pulses'),
        Parameter('n_stop', 5, int, 'end num of pulses'),
        Parameter('n_step', 1, int, 'step size for varying n'),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 80, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 80, float, 'wait time after MW pi pulse before sending in laser pulse'),
            Parameter('n_step', 1, int, 'send laser pulse every n pi pulses'),
            Parameter('n_start', 1, int, 'index for first laser pulse'),
        ]),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'commander': Commander}
    from b26_toolkit.scripts.find_nv import FindNvStrobe
    _SCRIPTS = {'find_nv': FindNvStrobe}

    def __init__(self, instruments, scripts, name=None, settings=None, log_function=None, data_path=None):
        """
        Standard script initialization
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_generic import PulsedExperimentGeneric
        self._DEFAULT_SETTINGS += PulsedExperimentGeneric._DEFAULT_SETTINGS
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.ref_index = 0

    def _create_pulse_sequences(self):

        """
        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement
        """

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_tau = self.settings['mw_pulses']['tau_mw']
        spacing = self.settings['mw_pulses']['spacing']

        if self.settings['mw_pulses']['microwave_channel'] == '+i, +q':
            microwave_channel_0, microwave_channel_1 = 'microwave_+i', 'microwave_+q'
        elif self.settings['mw_pulses']['microwave_channel'] == '+i, -i':
            microwave_channel_0, microwave_channel_1 = 'microwave_+i', 'microwave_-i'
        elif self.settings['mw_pulses']['microwave_channel'] == '+q, -q':
            microwave_channel_0, microwave_channel_1 = 'microwave_+q', 'microwave_-q'
        elif self.settings['mw_pulses']['microwave_channel'] == '+i, -i, +q, -q':
            microwave_channel_0 = 'microwave_+i'
            microwave_channel_1 = 'microwave_-i'
            microwave_channel_2 = 'microwave_+q'
            microwave_channel_3 = 'microwave_-q'
        elif self.settings['mw_pulses']['microwave_channel'] == '+i, +q, -i, -q':
            microwave_channel_0 = 'microwave_+i'
            microwave_channel_1 = 'microwave_+q'
            microwave_channel_2 = 'microwave_-i'
            microwave_channel_3 = 'microwave_-q'
        else:
            microwave_channel_0 = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        tau = self.settings['mw_pulses']['tau_mw']
        pulse_sequences = []

        if self.settings['n_start'] > self.settings['n_stop']:
            self.log('Error: end n must be larger than start n when range_type is start_stop.', flag='error')
            self._abort = True

        if self.settings['n_start'] < 0:
            self.log('Error: number of pulses cannot be negative', flag='error')
            self._abort = True

        n_list = np.arange(self.settings['n_start'], self.settings['n_stop'], self.settings['n_step'])

        n_list = np.array(n_list, dtype=int)

        reinit_pulse_duration = self.settings['initialization_laser']['pulse_duration']
        reinit_wait_duration = self.settings['initialization_laser']['wait_duration']
        reinit_n_step = self.settings['initialization_laser']['n_step']
        reinit_n_start = self.settings['initialization_laser']['n_start']

        for n in n_list:
            pulse_sequence = [Pulse('laser', laser_off_time + mw_tau, nv_reset_time),
                              Pulse('apd_readout', laser_off_time + mw_tau + delay_readout, meas_time)]
            current_time = laser_off_time + mw_tau + nv_reset_time + laser_off_time


            for i in range(n):
                if self.settings['mw_pulses']['microwave_channel'] == '+i, +q':
                    if i % 2 == 1:  # Alternate to MW chan 2 for odd-numbered pulses
                        pulse_sequence.append(Pulse(microwave_channel_1, current_time, mw_tau))
                    else:
                        pulse_sequence.append(Pulse(microwave_channel_0, current_time, mw_tau))
                elif self.settings['mw_pulses']['microwave_channel'] in ['+i, -i, +q, -q', '+i, +q, -i, -q']:
                    if i % 4 == 3:  # Alternate to MW chan 2 for odd-numbered pulses
                        pulse_sequence.append(Pulse(microwave_channel_3, current_time, mw_tau))
                    elif i % 4 == 2:
                        pulse_sequence.append(Pulse(microwave_channel_2, current_time, mw_tau))
                    elif i % 4 == 1:
                        pulse_sequence.append(Pulse(microwave_channel_1, current_time, mw_tau))
                    else:
                        pulse_sequence.append(Pulse(microwave_channel_0, current_time, mw_tau))
                else:
                    pulse_sequence.append(Pulse(microwave_channel_0, current_time, mw_tau))
                if (i + 1) % (reinit_n_step) == 0:
                    pulse_sequence.append(Pulse('laser', current_time + mw_tau + reinit_wait_duration, reinit_pulse_duration))
                current_time += mw_tau + spacing

            current_time += delay_mw_readout
            pulse_sequence.append(Pulse('laser', current_time, nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout', current_time + delay_readout, meas_time))
            pulse_sequence.append(Pulse('apd_readout', current_time + nv_reset_time - meas_time, meas_time))
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, n_list, self.settings['read_out']['meas_time']

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

        super(PiPulseTrain, self)._plot(axislist)
        axislist[0].set_title('Pi-pulse train')
        axislist[0].set_xlabel('Number of pulses')
        axislist[0].legend(labels=('Ref Fluorescence', 'Signal Fluorescence'), fontsize=8)
