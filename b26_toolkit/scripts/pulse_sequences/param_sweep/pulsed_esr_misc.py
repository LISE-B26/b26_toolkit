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
from b26_toolkit.scripts.pulse_sequences.param_sweep.pulsed_esr import PulsedEsr
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_generic import PulsedExperimentGeneric
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, Pulse, MicrowaveGenerator, MicrowaveGenerator2, RFGenerator, AFG3022C, Commander, AFG3022C_02
from b26_toolkit.plotting.plots_1d import plot_1d_simple_timetrace_ns, plot_pulses, update_pulse_plot, update_1d_simple, plot_1d_simple_freq
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from pylabcontrol.core import Parameter
import time as t


laser_pulse_end_delay = 100  # Time of end of PB pulse to AOM minus time of end of laser pulse

class PulsedEsrMwGen2(PulsedEsr):
    """
    Same as PulsedEsr but for a second microwave generator.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('microwave_channel', 'i_2', ['i_2', 'q_2'], 'Channel to use for mw pulses'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen_2': MicrowaveGenerator2, 'commander': Commander}

    def _configure_instruments_start_of_script(self):
        """
        Configure instruments at the beginning of a script, e.g. enable IQ modulation
        :return: None
        """
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_power']})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})

    def _configure_instruments_start_of_sweep(self, param_current):
        """
        Configure instruments before running a sweep. For example, change MW frequency
        :return: None
        """
        self.instruments['mw_gen_2']['instance'].update({'frequency': float(param_current)})


class PulsedEsrFastMwGen2(PulsedEsr):
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
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator2, 'commander': Commander}


class PulsedNmr(PulsedEsr):
    # Uses the AFG3022C to send MHz signals to the RF coil to drive the nuclear spin. Code is also greatly simplified compared to PulsedNMR_SRS
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.82e9, float,
                      'frequency of hyperfine transition corresponding to desired nuclear spin polarization'),
            Parameter('tau_mw', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('tau_rf', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('cooldown_time', 5000, float,
                  'cooldown time between running each sequence, to minimize heating effects'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator,
                    'afg': AFG3022C, 'commander': Commander}

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

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_tau = self.settings['mw_pulses']['tau_mw']
        microwave_channel = 'microwave_i'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        tau = self.settings['mw_pulses']['tau_mw']
        rf_tau = self.settings['rf_pulses']['tau_rf']
        cooldown_time = self.settings['cooldown_time']
        pulse_sequences = []
        tau_list = [tau]

        for tau in tau_list:
            pulse_sequence = [Pulse('laser', cooldown_time + laser_off_time, nv_reset_time),
                              #Pulse('apd_readout', cooldown_time + laser_off_time + delay_readout, meas_time),
                              Pulse(microwave_channel, cooldown_time + laser_off_time + nv_reset_time + laser_off_time, mw_tau)
                              ]
            # if tau is 0 there is actually no mw pulse
            if rf_tau > 0:
                pulse_sequence.append(
                    Pulse('rf_switch', cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle, rf_tau))
                pulse_sequence.append(
                    Pulse(rf_channel, cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle - 60, rf_tau + 120))

            pulse_sequence.append(Pulse(microwave_channel,
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + rf_tau + rf_settle,
                                        mw_tau))
            pulse_sequence.append(Pulse('laser',
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + rf_tau + rf_settle + mw_tau + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + rf_tau + rf_settle + mw_tau + delay_mw_readout + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        cooldown_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau + rf_settle + rf_tau + rf_settle + mw_tau + delay_mw_readout + nv_reset_time - meas_time - 200,
                                        meas_time))

            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:

            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, self.settings['read_out']['meas_time']

    def _configure_instruments_start_of_script(self):
        self.instruments['PB']['instance'].update({'rf_switch': {'status': False}})
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

        self.instruments['afg']['instance'].update({'ch1_invert_polarity': False})
        self.instruments['afg']['instance'].update({'ch2_invert_polarity': False})
        self.instruments['afg']['instance'].update({'ch1_phase': 0})
        self.instruments['afg']['instance'].update({'ch2_phase': 3.1415926})
        for channel in [1, 2]:
            self.instruments['afg']['instance'].update({'ch%i_amplitude' % channel: float(self.dbm_to_vpp(self.settings['rf_pulses']['rf_power']))})
            self.instruments['afg']['instance'].update({'ch%i_function' % channel: 'Sine'})
            #self.instruments['afg']['instance'].update({'ch%i_phase' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_offset' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_enable' % channel: True})

        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})


    def _configure_instruments_start_of_sweep(self, param_current):
        for channel in [1, 2]:
            self.instruments['afg']['instance'].update({'ch%i_frequency' % channel: float(param_current)})


class PulsedNmrSrs(PulsedNmr):
    """
    Uses the SRS to drive the RF coil to drive nuclear transitions. Not recommended; you should use PulsedNMR instead which uses the AFG3022C, which will be used
    to drive nuclear transitions in all other scripts anyway.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.82e9, float,
                      'frequency of hyperfine transition corresponding to desired nuclear spin polarization'),
            Parameter('tau_mw', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('tau_rf', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('cooldown_time', 5000, float,
                  'cooldown time between running each sequence, to minimize heating effects'),
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
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator,
                    'rf_gen': RFGenerator, 'commander': Commander}

    def _configure_instruments_start_of_script(self):
        self.instruments['PB']['instance'].update({'rf_switch': {'status': False}})
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

        # IQ modulation doesn't even work for MHz signals from the SRS!
        self.instruments['rf_gen']['instance'].update({'enable_modulation': False})
        self.instruments['rf_gen']['instance'].update({'enable_rf_output': True})
        self.instruments['rf_gen']['instance'].update({'amplitude_rf': self.settings['rf_pulses']['rf_power']})

        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})


    def _configure_instruments_start_of_sweep(self, param_current):
        self.instruments['rf_gen']['instance'].update({'frequency': param_current})


class PulsedEsrPolarized(PulsedNmr):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power_1', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency_1', 2.82e9, float,
                      'frequency of hyperfine transition corresponding to desired nuclear spin polarization'),
            Parameter('tau_mw_1', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_power_2', -45.0, float, 'microwave power in dBm'),
            Parameter('tau_mw_2', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('tau_rf', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('cooldown_time', 5000, float, 'cooldown time between running each sequence, to minimize heating effects'),
        Parameter('polarization_iterations', 2, int, 'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int, 'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            #Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int, 'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'mw_gen_2': MicrowaveGenerator2,
                    'afg': AFG3022C, 'commander': Commander}

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

        #nv_reset_time = self.settings['read_out']['nv_reset_time']
        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_tau_1 = self.settings['mw_pulses']['tau_mw_1']
        mw_tau_2 = self.settings['mw_pulses']['tau_mw_2']
        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        tau = self.settings['mw_pulses']['tau_mw_1']
        rf_tau = self.settings['rf_pulses']['tau_rf']
        cooldown_time = self.settings['cooldown_time']
        pulse_sequences = []
        tau_list = [tau]

        for tau in tau_list:

            pulse_sequence = [Pulse('laser', laser_off_time + cooldown_time, nv_reset_time_long)]
            current_time = laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_tau_1))
                if rf_tau > 0:
                    pulse_sequence.append(
                        Pulse('rf_switch', current_time + mw_tau_1 + delay_mw_readout, rf_tau))
                    pulse_sequence.append(
                        Pulse(rf_channel, current_time + mw_tau_1 + delay_mw_readout - 60, rf_tau + 120))
                    current_time = current_time + mw_tau_1 + delay_mw_readout + rf_tau + rf_settle
                #if self.settings['initialization_laser']['chopped']:
                for i in range(self.settings['initialization_laser']['n']):
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                    current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait
                #else:
                #    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time))
                #    current_time = current_time + nv_reset_time



            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time + laser_off_time,
                                        mw_tau_2))
            pulse_sequence.append(Pulse('laser',
                                        current_time + laser_off_time + mw_tau_2 + delay_mw_readout,
                                        nv_reset_time_long))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + laser_off_time + mw_tau_2 + delay_mw_readout + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + laser_off_time + mw_tau_2 + delay_mw_readout + nv_reset_time_long - meas_time - laser_pulse_end_delay,
                                        meas_time))
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:

            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, self.settings['read_out']['meas_time']

    def _configure_instruments_start_of_script(self):
        self.instruments['PB']['instance'].update({'rf_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')

        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power_1']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency_1']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power_2']})
        self.instruments['mw_gen_2']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency_1']})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})

        self.instruments['afg']['instance'].update({'ch1_invert_polarity': False})
        self.instruments['afg']['instance'].update({'ch2_invert_polarity': False})

        for channel in [1, 2]:
            self.instruments['afg']['instance'].update({'ch%i_amplitude' % channel: float(self.dbm_to_vpp(self.settings['rf_pulses']['rf_power']))})
            self.instruments['afg']['instance'].update({'ch%i_frequency' % channel: float(self.settings['rf_pulses']['rf_frequency'])})
            self.instruments['afg']['instance'].update({'ch%i_function' % channel: 'Sine'})
            #self.instruments['afg']['instance'].update({'ch%i_phase' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_offset' % channel: 0})
            self.instruments['afg']['instance'].update({'ch%i_enable' % channel: True})
            self.instruments['afg']['instance'].update({'ch%i_enable' % channel: True})

        self.instruments['afg']['instance'].align_phase()
        self.instruments['afg']['instance'].update({'ch1_phase': 0})
        self.instruments['afg']['instance'].update({'ch2_phase': 3.1415926})

    def _configure_instruments_start_of_sweep(self, param_current):
        self.instruments['mw_gen_2']['instance'].update({'frequency': float(param_current)})


class PulsedEsrRfHeating(PulsedEsrPolarized):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power_1', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency_1', 2.82e9, float,
                      'frequency of hyperfine transition corresponding to desired nuclear spin polarization'),
            Parameter('tau_mw_1', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_power_2', -45.0, float, 'microwave power in dBm'),
            Parameter('tau_mw_2', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('tau_rf', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('cooldown_time', 5000, float, 'cooldown time between running each sequence, to minimize heating effects'),
        Parameter('polarization_iterations', 2, int, 'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('initialization_laser', [
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('n', 1, int, 'num of initialization laser pulses, set n=1 and wait_duration=0 if you want to use one long rectangular laser pulse'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            #Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int, 'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'mw_gen_2': MicrowaveGenerator2,
                    'rf_gen': RFGenerator}

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

        #nv_reset_time = self.settings['read_out']['nv_reset_time']
        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_tau_1 = self.settings['mw_pulses']['tau_mw_1']
        mw_tau_2 = self.settings['mw_pulses']['tau_mw_2']
        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_pulses']['rf_settle']

        tau = self.settings['mw_pulses']['tau_mw_1']
        rf_tau = self.settings['rf_pulses']['tau_rf']
        cooldown_time = self.settings['cooldown_time']
        pulse_sequences = []
        tau_list = [tau]

        for tau in tau_list:



            pulse_sequence = []
            current_time = laser_off_time + cooldown_time
            if rf_tau > 0:
                pulse_sequence.append(
                    Pulse(rf_channel, current_time, rf_tau))
                pulse_sequence.append(
                    Pulse('rf_switch', current_time, rf_tau))

            current_time += rf_tau + rf_settle
            pulse_sequence.append(Pulse(microwave_channel_2,
                                        current_time,
                                        mw_tau_2))
            pulse_sequence.append(Pulse('laser',
                                        current_time + mw_tau_2 + delay_mw_readout,
                                        nv_reset_time_long))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + mw_tau_2 + delay_mw_readout + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + mw_tau_2 + delay_mw_readout + nv_reset_time_long - meas_time,
                                        meas_time))


            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:

            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, self.settings['read_out']['meas_time']


class NuclearPolarizationLaserDuration(PulsedExperimentGeneric):
    # Please update this with structure from PulsedESRPolarized, FF
    _DEFAULT_SETTINGS = [
        Parameter('mw_power_1', -45.0, float, 'microwave power in dBm'),
        Parameter('mw_frequency_1', 2.82e9, float, 'frequency of hyperfine transition that we want to depopulate'),
        Parameter('mw_frequency_2', 2.82e9, float, 'frequency of hyperfine transition for which we want to measure the population'),
        Parameter('tau_mw_1', 80, float, 'the time duration of the microwaves in ns'),
        Parameter('mw_power_2', -45.0, float, 'microwave power in dBm'),
        Parameter('tau_mw_2', 80, float, 'the time duration of the microwaves in ns'),
        Parameter('rf_power', -10, float, 'RF power in dBm'),
        Parameter('tau_rf', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
        Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns'),
        Parameter('cooldown_time', 5000, float, 'cooldown time between running each sequence'),
        Parameter('polarization_iterations', 2, int, 'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('initialization_laser', [
            Parameter('chopped', False, bool, 'whether to use chopped laser pulses for initialization, overrides nv_reset_time if True'),
            Parameter('pulse_duration', 6, float, 'duration of each chopped laser pulse'),
            Parameter('wait_duration', 30, float, 'spacing between each chopped laser pulse'),
            Parameter('min_n', 1, int, 'min num of initialization laser pulses'),
            Parameter('max_n', 200, int, 'max num of initialization laser pulses'),
            Parameter('n_step', 1, int, 'step increment of number of initialization laser pulses')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('nv_reset_time_long', 5000, int, 'time with laser on to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ])
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'mw_gen_2': MicrowaveGenerator2,
                    'rf_gen': RFGenerator}

    def _function(self, in_data=None):
        self.instruments['PB']['instance'].update({'rf_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['rf_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['rf_gen']['instance'].update({'enable_modulation': True})
        self.instruments['rf_gen']['instance'].update({'amplitude_rf': self.settings['rf_power']})
        #self.instruments['rf_gen']['instance'].update({'frequency': self.settings['rf_frequency']})
        #self.instruments['rf_gen']['instance'].update({'enable_output': True})

        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_power_1']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_frequency_1']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_power_2']})
        self.instruments['mw_gen_2']['instance'].update({'frequency': self.settings['mw_frequency_2']})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})

        super(NuclearPolarizationLaserDuration, self)._function()

        if 'counts' in self.data.keys() and 'tau' in self.data.keys():
            counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
            tau = self.data['tau']


            self.data['fits'] = None
            self.log('rabi fit failed')

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

        if self.settings['initialization_laser']['chopped']:
            tau_list = range(int(self.settings['initialization_laser']['min_n']),
                             int(self.settings['initialization_laser']['max_n']),
                             int(self.settings['initialization_laser']['n_step']))
        else:
            print('Not implemented for CW')

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        nv_reset_time_long = self.settings['read_out']['nv_reset_time_long']
        nv_reset_time_pulsed = self.settings['initialization_laser']['pulse_duration']
        nv_reset_time_wait = self.settings['initialization_laser']['wait_duration']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_tau_1 = self.settings['tau_mw_1']
        mw_tau_2 = self.settings['tau_mw_2']
        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'
        rf_channel = 'rf_i'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        rf_settle = self.settings['rf_settle']

        #tau = self.settings['tau_mw_1']
        rf_tau = self.settings['tau_rf']
        pulse_sequences = []

        tau_list = range(int(self.settings['initialization_laser']['min_n']),
                         int(self.settings['initialization_laser']['max_n']),
                         int(self.settings['initialization_laser']['n_step']))

        for tau in tau_list:
            # cooldown_time = int(rf_tau/2. / 2.) * 2 # time between each pulse sequence, as a multiple of tau, such that we have a constant duty cycle of rf pulses to avoid overheating
            cooldown_time = self.settings['cooldown_time']
            pulse_sequence = [Pulse('laser', laser_off_time + cooldown_time, nv_reset_time_long)]
            current_time = laser_off_time + cooldown_time + nv_reset_time_long + laser_off_time  # Begin building the pulse sequence at a given offset, which acts as a cooldown time between each sequence
            for iteration in range(self.settings['polarization_iterations']):
                pulse_sequence.append(Pulse(microwave_channel, current_time, mw_tau_1))
                if rf_tau > 0:
                    pulse_sequence.append(
                        Pulse('rf_switch', current_time + mw_tau_1 + rf_settle, rf_tau))
                    pulse_sequence.append(
                        Pulse(rf_channel, current_time + mw_tau_1 + rf_settle - 60, rf_tau + 120))
                    current_time = current_time + mw_tau_1 + rf_settle + rf_tau + rf_settle
                if self.settings['initialization_laser']['chopped']:
                    for i in range(int(tau)):
                        pulse_sequence.append(Pulse('laser', current_time, nv_reset_time_pulsed))
                        current_time = current_time + nv_reset_time_pulsed + nv_reset_time_wait
                else:
                    raise NotImplementedError
                    pulse_sequence.append(Pulse('laser', current_time, nv_reset_time))
                    current_time = current_time + nv_reset_time

            pulse_sequence.append(Pulse(microwave_channel,
                                        current_time + laser_off_time,
                                        mw_tau_1))
            pulse_sequence.append(Pulse('laser',
                                        current_time + laser_off_time + mw_tau_1 + delay_mw_readout,
                                        nv_reset_time_long))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + laser_off_time + mw_tau_1 + delay_mw_readout + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + laser_off_time + mw_tau_1 + delay_mw_readout + nv_reset_time_long - meas_time,
                                        meas_time))
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time


class CnotFidelity(PulsedEsrPolarized):
    """
    Executes a CNOT gate before a PulsedESR. Depending on the fidelity of this first CNOT gate, the ESR lineshape changes:
    A perfect CNOT gate one transition A will lead to transition A disappearing in the ESR while unaffecting other transitions.

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power_1', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency_1', 2.82e9, float,
                      'frequency of hyperfine transition corresponding to desired nuclear spin polarization'),
            Parameter('tau_mw_1', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_power_2', -45.0, float, 'microwave power in dBm'),
            Parameter('tau_mw_2', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'mw_gen_2': MicrowaveGenerator2,
                    'commander': Commander}

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

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_tau_1 = self.settings['mw_pulses']['tau_mw_1']
        mw_tau_2 = self.settings['mw_pulses']['tau_mw_2']
        microwave_channel = 'microwave_i'
        microwave_channel_2 = 'microwave_i_2'

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        tau = self.settings['mw_pulses']['tau_mw_1']
        pulse_sequences = []
        tau_list = [tau]

        current_time = 0
        for tau in tau_list:
            pulse_sequence = []
            for part in range(2):

                pulse_sequence += \
                    [Pulse('laser', current_time + laser_off_time, nv_reset_time),
                     Pulse('apd_readout', current_time + laser_off_time + delay_readout, meas_time)]
                if part == 0:
                    pulse_sequence += \
                        [Pulse(microwave_channel, current_time + laser_off_time + nv_reset_time + laser_off_time, mw_tau_1)]
                    current_time = current_time + laser_off_time + nv_reset_time + laser_off_time + mw_tau_1 + delay_mw_readout
                else:
                    current_time = current_time + laser_off_time + nv_reset_time + laser_off_time

                if tau > 0:
                    pulse_sequence += [
                        Pulse(microwave_channel_2, current_time, mw_tau_2)]

                pulse_sequence += [
                    Pulse('laser',
                          current_time + mw_tau_2 + delay_mw_readout,
                          nv_reset_time),
                    Pulse('apd_readout',
                          current_time + mw_tau_2 + delay_mw_readout + delay_readout,
                          meas_time)
                ]
                current_time = current_time + mw_tau_2 + delay_mw_readout + nv_reset_time
            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, self.settings['read_out']['meas_time']

    def _configure_instruments_start_of_script(self):
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')

        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power_1']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency_1']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power_2']})
        self.instruments['mw_gen_2']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency_1']})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})

    def _configure_instruments_start_of_sweep(self, param_current):
        self.instruments['mw_gen_2']['instance'].update({'frequency': float(param_current)})


class CnotPowerSweep(PulsedEsr):
    """
    Sweeps the power to optimize a pi-pulse given a fixed pi time.

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_frequency', 2.82e9, float, 'frequency of hyperfine transition'),
            Parameter('microwave_channel', 'i', ['i', 'q', 'i_2', 'q_2'], 'Channel to use for mw pulses'),
            Parameter('tau_mw', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('n', 0, int, 'number of pi-pulses, the larger n is the more accurate we can determine the correct power for a pi-pulse'),
            Parameter('spacing', 100, float, 'spacing in ns between consecutive pi-pulses')
        ]),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('power_start', -20, float, 'start frequency of scan in Hz'),
        Parameter('power_stop', -10, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('power_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset electronic spin'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator, 'mw_gen_2': MicrowaveGenerator2, 'commander': Commander}

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

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        mw_tau = self.settings['mw_pulses']['tau_mw']
        spacing = self.settings['mw_pulses']['spacing']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        n = int(self.settings['mw_pulses']['n'])

        tau = self.settings['mw_pulses']['tau_mw']
        pulse_sequences = []
        tau_list = [tau]

        for tau in tau_list:
            pulse_sequence = [Pulse('laser', laser_off_time + mw_tau, nv_reset_time),
                              Pulse('apd_readout', laser_off_time + mw_tau + delay_readout, meas_time)]

            current_time = laser_off_time + mw_tau + nv_reset_time + laser_off_time
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                for i in range(n):
                    pulse_sequence.append(Pulse(microwave_channel,
                                                current_time,
                                                mw_tau))
                    current_time += mw_tau + spacing

            pulse_sequence.append(Pulse('laser',
                                        current_time,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        current_time + delay_readout,
                                        meas_time))
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, self.settings['read_out']['meas_time']

    def _configure_instruments_start_of_script(self):
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')

        mw_channel = self.settings['mw_pulses']['microwave_channel']
        if mw_channel == 'i' or mw_channel == 'q':
            mw_gen_instrument = 'mw_gen'
        elif mw_channel == 'i_2' or mw_channel == 'q_2':
            mw_gen_instrument = 'mw_gen_2'

        self.instruments[mw_gen_instrument]['instance'].update({'modulation_type': 'IQ'})
        self.instruments[mw_gen_instrument]['instance'].update({'enable_modulation': True})
        self.instruments[mw_gen_instrument]['instance'].update({'amplitude': self.settings['power_start']})
        self.instruments[mw_gen_instrument]['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments[mw_gen_instrument]['instance'].update({'enable_output': True})



    def _configure_instruments_start_of_sweep(self, power_current):
        mw_channel = self.settings['mw_pulses']['microwave_channel']
        if mw_channel == 'i' or mw_channel == 'q':
            mw_gen_instrument = 'mw_gen'
        elif mw_channel == 'i_2' or mw_channel == 'q_2':
            mw_gen_instrument = 'mw_gen_2'

        self.instruments[mw_gen_instrument]['instance'].update({'amplitude': float(power_current)})

    def _configure_param_array(self):
        # Contruct the frequency array and store it in a variable called 'mw_frequencies'. Despite the naming, it's just a list of parameters to be swept;
        # it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode

        if self.settings['range_type'] == 'start_stop':
            if self.settings['power_start'] > self.settings['power_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True


            self.params = np.linspace(self.settings['power_start'], self.settings['power_stop'], self.settings['power_points'])

        elif self.settings['range_type'] == 'center_range':

            self.params = np.linspace(self.settings['power_start'] - self.settings['power_stop'] / 2,
                                      self.settings['power_start'] + self.settings['power_stop'] / 2, self.settings['power_points'])


class FieldProfilePulsed(PulsedEsr):
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 10, int, 'number of frequencies in scan'),
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'beginning of first readout window'),
            Parameter('max_time', 10000, float, 'beginning of last readout window'),
            Parameter('time_step', 1000, float, 'time step between beginning of readout windows')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 1000, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 1000, int, 'delay between mw and readout (in ns)'),
            #Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('num_cycles', 5, int, 'number of initialization + pi-pulse cycles for each readout; the time step must be longer than '
                                            'num_cycles*(nv_reset_time+laser_off_time+delay_mw_readout)')
        ]),
        Parameter('atto_trig', [
            Parameter('trig_duration', 200, float,
                      'length of trigger pulse in ns to send to function generator controlling attocube'),
            Parameter('settle_time', 200000, float,
                      'time in ns to wait between each pulse sequence'),
            Parameter('total_time', 600000, float, 'length of total pulse sequence'),
        ]),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)'),
        Parameter('num_averages', 100000, int, 'number of averages')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _create_pulse_sequences(self):
        '''

        Returns: pulse_sequence, num_averages, tau_list, meas_time
            pulse_sequences: a single pulse sequences, with daq reads at specified taus
            num_averages: the number of times to repeat each pulse sequence

        '''

        if self.settings['read_out']['meas_time'] > self.settings['tau_times']['time_step']:
            self.log('Readout window must be shorter between spacing between beginning of readout windows!')
            raise AttributeError

        subblock_tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],
                             self.settings['tau_times']['time_step'])
        subblock_tau_list = np.ndarray.tolist(subblock_tau_list)  # 20180731 ER convert to list
        self.subblock_tau_list = subblock_tau_list

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        subblock_tau_list = [x for x in subblock_tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        #delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_i'# + self.settings['mw_pulses']['microwave_channel']
        pi_time = self.settings['tau_mw']
        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        num_cycles = self.settings['read_out']['num_cycles']

        atto_trig_duration = self.settings['atto_trig']['trig_duration']
        atto_settle_time = self.settings['atto_trig']['settle_time']
        total_time = self.settings['atto_trig']['total_time']

        tau = self.settings['tau_mw']
        pulse_sequences = []
        tau_list = [tau]

        for tau in tau_list:
            pulse_sequence = [Pulse('atto_trig', atto_settle_time, atto_trig_duration)]
            pulse_sequence.append(Pulse('spacer', total_time - 200, 200))
            """
            for subblock_tau in subblock_tau_list:
                # if tau is 0 there is actually no mw pulse
                if subblock_tau > 0:
                    block_beginning = atto_settle_time + subblock_tau
                    pulse_sequence += [Pulse(microwave_channel, block_beginning - 30, pi_time + 30 * 2)]
                    pulse_sequence += [Pulse('microwave_switch', block_beginning, pi_time)]
                    pulse_sequence += [Pulse('laser', block_beginning + pi_time + delay_mw_readout, meas_time + nv_reset_time)]
                    pulse_sequence += [Pulse('apd_readout', block_beginning + pi_time + delay_mw_readout + delay_readout, meas_time)]
            """

            for subblock_tau in subblock_tau_list:
                # if tau is 0 there is actually no mw pulse
                if subblock_tau > 0:
                    block_beginning = atto_settle_time + subblock_tau
                    for i in range(num_cycles):
                        pulse_sequence += [Pulse(microwave_channel, block_beginning, pi_time)]
                        #pulse_sequence += [Pulse('microwave_switch', block_beginning, pi_time)]
                        pulse_sequence += [Pulse('laser', block_beginning + pi_time + delay_mw_readout, nv_reset_time)]
                        block_beginning += pi_time + delay_mw_readout + nv_reset_time + laser_off_time

                    pulse_sequence += [Pulse('apd_readout', atto_settle_time + subblock_tau, block_beginning - (atto_settle_time + subblock_tau))]
            pulse_sequences.append(pulse_sequence)

        #self.log('Total length of pulse sequence (including settling time): %i' % (block_beginning + pi_time + delay_mw_readout + delay_readout + meas_time))
        self.log('Total length of pulse sequence (including settling time): %i' %
                 (atto_settle_time + subblock_tau + block_beginning - (atto_settle_time + subblock_tau)))

        meas_time = block_beginning - (atto_settle_time + subblock_tau)
        return pulse_sequences, tau_list, meas_time

    def _plot(self, axes_list, data=None):
        """
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional) dataset to plot (dictionary that contains keys counts, tau), if not provided use self.data
        """

        if data is None:
            data = self.data

        if 'counts' in data.keys():
            plot_pulses(axes_list[1], self.pulse_sequences[self.sequence_index])
            counts = np.array(self.data['count_data'])
            label = ['Pulsed ESR while moving', 'time', 'freq', 'kcounts/s']
            extent = self.subblock_tau_list[0], self.subblock_tau_list[-1], data['params'][-1], data['params'][0]
            aspect = np.abs((extent[1]-extent[0])/(extent[3]-extent[2]))
            plot_fluorescence_new(counts, extent, axes_list[0], max_counts=-1, labels=label, aspect=aspect)


    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''

        counts = np.array(self.data['count_data'])

        update_fluorescence(counts, axes_list[0], -1)
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])


class FieldProfileCw(FieldProfilePulsed):
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 10, int, 'number of frequencies in scan'),
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'beginning of first readout window'),
            Parameter('max_time', 10000, float, 'beginning of last readout window'),
            Parameter('time_step', 1000, float, 'time step between beginning of readout windows'),
            Parameter('readout_window', 800, float, 'length of readout window, must be smaller than time_step')
        ]),

        Parameter('atto_trig', [
            Parameter('trig_duration', 200, float,
                      'length of trigger pulse in ns to send to function generator controlling attocube'),
            Parameter('settle_time', 200000, float,
                      'time in ns to wait between each pulse sequence'),
            Parameter('total_time', 600000, float, 'length of total pulse sequence'),
        ]),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)'),
        Parameter('num_averages', 100000, int, 'number of averages')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _create_pulse_sequences(self):
        '''

        Returns: pulse_sequence, num_averages, tau_list, meas_time
            pulse_sequences: a single pulse sequences, with daq reads at specified taus
            num_averages: the number of times to repeat each pulse sequence

        '''

        subblock_tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],
                             self.settings['tau_times']['time_step'])
        subblock_tau_list = np.ndarray.tolist(subblock_tau_list)  # 20180731 ER convert to list
        self.subblock_tau_list = subblock_tau_list

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        subblock_tau_list = [x for x in subblock_tau_list if x == 0 or x >= min_pulse_dur]

        microwave_channel = 'microwave_i'

        atto_trig_duration = self.settings['atto_trig']['trig_duration']
        atto_settle_time = self.settings['atto_trig']['settle_time']
        total_time = self.settings['atto_trig']['total_time']

        readout_window = self.settings['tau_times']['readout_window']

        tau = 100
        pulse_sequences = []
        tau_list = [tau]

        for tau in tau_list:
            pulse_sequence = [Pulse('atto_trig', atto_settle_time, atto_trig_duration)]
            pulse_sequence.append(Pulse('spacer', total_time - 200, 200))

            for subblock_tau in subblock_tau_list:
                # if tau is 0 there is actually no mw pulse
                if subblock_tau > 0:
                    block_beginning = atto_settle_time + subblock_tau
                    pulse_sequence += [Pulse(microwave_channel, block_beginning, readout_window)]
                    pulse_sequence += [Pulse('laser', block_beginning, readout_window)]
                    pulse_sequence += [Pulse('apd_readout', block_beginning, readout_window)]
            pulse_sequences.append(pulse_sequence)

        meas_time = readout_window
        return pulse_sequences, tau_list, meas_time


class TransportDqmaRfPhase(PulsedEsrPolarized):
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
            Parameter('mw_2_pi_time', 80, float, 'the time duration of the microwaves in ns')
        ]),
        Parameter('rf_pulses', [
            Parameter('rf_power', -10, float, 'RF power in dBm'),
            Parameter('rf_frequency', 3e6, float, 'frequency of hyperfine interaction in Hz'),
            Parameter('rf_pi_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_three_pi_half_time', 5000, float, 'time for RF pulse to drive nuclear spin in ns'),
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('rf_phase', [
            Parameter('min_phase', 0, float, 'minimum phase in radians of final RF pi/2 pulse, rel. to first RF pi/2 pulse'),
            Parameter('max_phase', 6.28318531, float, 'maximum phase in radians of final RF pi/2 pulse, rel. to first RF pi/2 pulse'),
            Parameter('phase_points', 40, int, 'number of phase steps'),
        ]),
        Parameter('t_sense', 0, float, 'sensing time in ns'),
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
        ]),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching the phase on RF generator (s)')
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

        tau_list = [self.settings['t_sense']]

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
            for part in range(1):
                start_time = 0 + int(total_time * part)
                current_time = start_time

                # Send a pulse that doesn't do anything for the total movement time, to set the durations of part 0 and 1 of the pulse sequence
                pulse_sequence.append(Pulse('spacer', total_time - 200, 200))
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
                    if part == 0:
                        pulse_sequence.append(Pulse(channel_odd,
                        #pulse_sequence.append(Pulse(microwave_channel,
                                                    current_time,
                                                    pi_time_odd))

                        # Wait for t_sense and disentangle from memory
                        pulse_sequence.append(Pulse(channel_odd,
                        #pulse_sequence.append(Pulse(microwave_channel,
                                                    current_time + pi_time_odd + t_sense,
                                                    pi_time_odd))

                        current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout
                    elif part == 1:
                        pulse_sequence.append(Pulse(channel_even,
                        #pulse_sequence.append(Pulse(microwave_channel_2,
                                                    current_time,
                                                    pi_time_even))

                        # Wait for t_sense and disentangle from memory
                        pulse_sequence.append(Pulse(channel_even,
                        #pulse_sequence.append(Pulse(microwave_channel_2,
                                                    current_time + pi_time_even + t_sense,
                                                    pi_time_even))

                        current_time += pi_time_even + t_sense + pi_time_even + delay_mw_readout

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


    def _configure_instruments_start_of_sweep(self, ch2_phase_current):
        if ch2_phase_current == self.settings['rf_phase']['max_phase'] + 2000:
            print('Norm meas for min fluor')
            self.instruments['afg']['instance'].update({'ch2_phase': np.pi})
        elif ch2_phase_current == self.settings['rf_phase']['max_phase'] + 1000:
            print('Norm meas for max fluor')
            self.instruments['afg']['instance'].update({'ch2_phase': 0})
        else:
            self.instruments['afg']['instance'].update({'ch2_phase': ch2_phase_current})

    def _configure_param_array(self):
        # Contruct the frequency array and store it in a variable called 'mw_frequencies'. Despite the naming, it's just a list of parameters to be swept;
        # it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        self.params = np.linspace(self.settings['rf_phase']['min_phase'], self.settings['rf_phase']['max_phase'], self.settings['rf_phase']['phase_points'])


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
            counts = data['count_data'][:, 1] / data['count_data'][:, 0]
            tau = data['tau']
            fits = data['fits'] # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            #axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            pi_time = (np.pi - fits[2])/fits[1]
            pi_half_time = (np.pi/2 - fits[2])/fits[1]
            three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
            rabi_freq = 1000*fits[1]/(2*np.pi)
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            if 'count_data' in data.keys():
                num_daq_reads = self.settings['read_out']['repetitive_readout']['m']
                avg_counts_1 = np.transpose(np.array([np.average(self.data['count_data'][:, 0:num_daq_reads], axis=1)]))
                first_counts_1 = np.transpose(np.array([np.average(self.data['count_data'][:, 0:1], axis=1)]))
                avg_counts_2 = np.transpose(np.array([np.average(self.data['count_data'][:, num_daq_reads:], axis=1)]))
                first_counts_2 = np.transpose(np.array([np.average(self.data['count_data'][:, num_daq_reads:num_daq_reads+1], axis=1)]))

                #plot_1d_simple_timetrace_ns(axislist[0], data['tau'], [first_counts_1, first_counts_2, avg_counts_1, avg_counts_2])
                plot_1d_simple_timetrace_ns(axislist[0], data['params'],[avg_counts_1, avg_counts_2])
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

        x_data = self.data['params']


        num_daq_reads = self.settings['read_out']['repetitive_readout']['m']
        avg_counts_1 = np.transpose(np.array([np.average(self.data['count_data'][:, 0:num_daq_reads], axis=1)]))
        first_counts_1 = np.transpose(np.array([np.average(self.data['count_data'][:, 0:1], axis=1)]))
        avg_counts_2 = np.transpose(np.array([np.average(self.data['count_data'][:, num_daq_reads:], axis=1)]))
        first_counts_2 = np.transpose(np.array([np.average(self.data['count_data'][:, num_daq_reads:num_daq_reads + 1], axis=1)]))

        axis1 = axes_list[0]
        if not self.data['count_data'] == []:
            update_1d_simple(axis1, x_data, [avg_counts_1, avg_counts_2])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])


class TransportDqmaMwPhase(PulsedEsrPolarized):
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
            Parameter('rf_settle', 2000, float, 'settling time before and after RF pulse in ns')
        ]),
        Parameter('polarization_iterations', 2, int,
                  'number of times to repeat the nuclear spin initialization sequence'),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('mw_phase', [
            Parameter('min_phase', 0, float, 'minimum phase in radians of final RF pi/2 pulse, rel. to first RF pi/2 pulse'),
            Parameter('max_phase', 6.28318531, float, 'maximum phase in radians of final RF pi/2 pulse, rel. to first RF pi/2 pulse'),
            Parameter('phase_points', 40, int, 'number of phase steps'),
        ]),
        Parameter('t_sense', 0, float, 'sensing time in ns'),
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
        ]),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching the phase on RF generator (s)')
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

        tau_list = [self.settings['t_sense']]

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
                    if part == 0 or part == 1:
                        pulse_sequence.append(Pulse(channel_odd,
                        #pulse_sequence.append(Pulse(microwave_channel,
                                                    current_time,
                                                    pi_time_odd))

                        # Wait for t_sense and disentangle from memory
                        pulse_sequence.append(Pulse(channel_odd,
                        #pulse_sequence.append(Pulse(microwave_channel,
                                                    current_time + pi_time_odd + t_sense,
                                                    pi_time_odd))

                        current_time += pi_time_odd + t_sense + pi_time_odd + delay_mw_readout
                    elif part == 99999:
                        pulse_sequence.append(Pulse(channel_even,
                        #pulse_sequence.append(Pulse(microwave_channel_2,
                                                    current_time,
                                                    pi_time_even))

                        # Wait for t_sense and disentangle from memory
                        pulse_sequence.append(Pulse(channel_even,
                        #pulse_sequence.append(Pulse(microwave_channel_2,
                                                    current_time + pi_time_even + t_sense,
                                                    pi_time_even))

                        current_time += pi_time_even + t_sense + pi_time_even + delay_mw_readout

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

                pulse_sequence.append(Pulse('rf_switch', t_movement_end - (rf_pi_half_time + rf_settle), rf_pi_half_time))
                if part == 1:
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


    def _configure_instruments_start_of_script(self):
        self.instruments['PB']['instance'].update({'rf_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')

        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_1_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_1_frequency']})
        self.instruments['mw_gen']['instance'].update({'phase': float(0)})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_2_power']})
        self.instruments['mw_gen_2']['instance'].update({'frequency': self.settings['mw_pulses']['mw_2_frequency']})
        self.instruments['mw_gen_2']['instance'].update({'phase': float(0)})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})

        for channel in [1, 2]:
            self.instruments['afg']['instance'].update({'ch%i_amplitude' % channel: float(self.dbm_to_vpp(self.settings['rf_pulses']['rf_power']))})
            self.instruments['afg']['instance'].update({'ch%i_frequency' % channel: float(self.settings['rf_pulses']['rf_frequency'])})
            self.instruments['afg']['instance'].update({'ch%i_run_mode' % channel: 'Continuous'})
            self.instruments['afg']['instance'].update({'ch%i_function' % channel: 'Sine'})
            self.instruments['afg']['instance'].update({'ch%i_invert_polarity' % channel: False})
            self.instruments['afg']['instance'].update({'ch%i_offset' % channel: 0})

        self.instruments['afg']['instance'].align_phase()
        self.instruments['afg']['instance'].update({'ch1_phase': 0})
        self.instruments['afg']['instance'].update({'ch2_phase': 3.1415926})

    def _configure_instruments_start_of_sweep(self, phase_current):
        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_2_power']})
        self.instruments['mw_gen_2']['instance'].update({'frequency': self.settings['mw_pulses']['mw_2_frequency']})
        #print(phase_current)
        self.instruments['mw_gen_2']['instance'].update({'phase': float(phase_current)})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})

    def _configure_param_array(self):
        # Contruct the frequency array and store it in a variable called 'mw_frequencies'. Despite the naming, it's just a list of parameters to be swept;
        # it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        self.params = np.linspace(self.settings['mw_phase']['min_phase'], self.settings['mw_phase']['max_phase'], self.settings['mw_phase']['phase_points'])


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
            counts = data['count_data'][:, 1] / data['count_data'][:, 0]
            tau = data['tau']
            fits = data['fits'] # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            #axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            pi_time = (np.pi - fits[2])/fits[1]
            pi_half_time = (np.pi/2 - fits[2])/fits[1]
            three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
            rabi_freq = 1000*fits[1]/(2*np.pi)
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            if 'count_data' in data.keys():
                num_daq_reads = self.settings['read_out']['repetitive_readout']['m']
                avg_counts_1 = np.transpose(np.array([np.average(self.data['count_data'][:, 0:num_daq_reads], axis=1)]))
                first_counts_1 = np.transpose(np.array([np.average(self.data['count_data'][:, 0:1], axis=1)]))
                avg_counts_2 = np.transpose(np.array([np.average(self.data['count_data'][:, num_daq_reads:], axis=1)]))
                first_counts_2 = np.transpose(np.array([np.average(self.data['count_data'][:, num_daq_reads:num_daq_reads+1], axis=1)]))

                #plot_1d_simple_timetrace_ns(axislist[0], data['tau'], [first_counts_1, first_counts_2, avg_counts_1, avg_counts_2])
                plot_1d_simple_timetrace_ns(axislist[0], data['params'],[avg_counts_1, avg_counts_2])
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

        x_data = self.data['params']


        num_daq_reads = self.settings['read_out']['repetitive_readout']['m']

        avg_counts_1 = np.transpose(np.array([np.average(self.data['count_data'][:, 0:num_daq_reads], axis=1)]))
        first_counts_1 = np.transpose(np.array([np.average(self.data['count_data'][:, 0:1], axis=1)]))
        avg_counts_2 = np.transpose(np.array([np.average(self.data['count_data'][:, num_daq_reads:], axis=1)]))
        first_counts_2 = np.transpose(np.array([np.average(self.data['count_data'][:, num_daq_reads:num_daq_reads + 1], axis=1)]))

        axis1 = axes_list[0]
        if not self.data['count_data'] == []:
            update_1d_simple(axis1, x_data, [avg_counts_1, avg_counts_2])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])


class MicrowaveRotationAxis(PulsedEsrPolarized):
    """
    With a finite detuning of the RF from the nuclear spin frequency, we expect to see oscillations for long movmement times
    """

    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_2_power', -45.0, float, 'microwave power in dBm'),
            Parameter('mw_2_frequency', 2.82e9, float,
                      'frequency of hyperfine transition 2, also the transition to be populated for initialization'),
            Parameter('mw_2_pi_half_time', 80, float, 'the time duration of the microwaves in ns'),
            Parameter('mw_pulse_spacing', 200, float, 'spacing between consecutive microwave pulses')
        ]),
        Parameter('num_averages', 5000000, int, 'number of averages'),
        Parameter('mw_phase', [
            Parameter('min_phase', 0, float, 'minimum phase in radians when phase control pulse is applied, I(Q) corresponds to 0(pi/2) rad'),
            Parameter('max_phase', 6.28318531, float, 'maximum phase in radians when phase control pulse is applied'),
            Parameter('phase_points', 40, int, 'number of phase steps'),
            Parameter('idle_phase', 0, float, 'idle phase in radians when phase control pulse is NOT applied')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 5000, int, 'duration of long laser pulse to reset both electronic and nuclear spin'),
            Parameter('laser_off_time', 500, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 250, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 100, int, 'delay between laser on and readout (given by spontaneous decay rate)')]),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching the phase on SRS (s)')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator,
                    'afg_iq': AFG3022C_02, 'mw_gen_2': MicrowaveGenerator2, 'commander': Commander}

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

        tau_list = [self.settings['mw_pulses']['mw_2_pi_half_time']]
        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        short_pulses = [x for x in tau_list if x < min_pulse_dur]
        print('Found short pulses: ', short_pulses)
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        microwave_channel, microwave_channel_2, microwave_channel_q, microwave_channel_2_q, rf_channel = \
            'microwave_i', 'microwave_i_2', 'microwave_q', 'microwave_q_2', 'rf_i'

        mw_2_pi_half_time = self.settings['mw_pulses']['mw_2_pi_half_time']
        mw_pulse_spacing = self.settings['mw_pulses']['mw_pulse_spacing']

        meas_time, nv_reset_time, laser_off_time, delay_mw_readout, delay_readout = \
            [self.settings['read_out'][key] for key in
             ['meas_time', 'nv_reset_time', 'laser_off_time', 'delay_mw_readout', 'delay_readout']]

        pulse_sequences = []
        for tau in tau_list:
            pulse_sequence = []
            if tau > 0:
                pulse_sequence.append(Pulse('laser', laser_off_time, nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout', laser_off_time + delay_readout, meas_time))
                pulse_sequence.append(Pulse('microwave_switch', laser_off_time + nv_reset_time + laser_off_time, mw_2_pi_half_time))
                pulse_sequence.append(Pulse('microwave_i_2', laser_off_time + nv_reset_time + laser_off_time - 30,
                                            mw_2_pi_half_time*2+400))
                pulse_sequence.append(Pulse('microwave_switch', laser_off_time + nv_reset_time + laser_off_time + mw_2_pi_half_time + mw_pulse_spacing,
                                            mw_2_pi_half_time))
                #pulse_sequence.append(Pulse('microwave_i_2', laser_off_time + nv_reset_time + laser_off_time + mw_2_pi_half_time + mw_pulse_spacing - 30,
                #                            250))

                pulse_sequence.append(Pulse('laser',
                                            laser_off_time + nv_reset_time + laser_off_time + mw_2_pi_half_time + mw_pulse_spacing + mw_2_pi_half_time + delay_mw_readout,
                                            nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                            laser_off_time + nv_reset_time + laser_off_time + mw_2_pi_half_time + mw_pulse_spacing + mw_2_pi_half_time + delay_mw_readout + delay_readout,
                                            meas_time))

            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _configure_instruments_start_of_script(self):
        self.instruments['PB']['instance'].update({'rf_switch': {'status': False}})
        # print('successfully turned off microwave switch')
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        # print('successfully turned off microwave switch')

        self.instruments['mw_gen_2']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_2']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_2']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_2_power']})
        self.instruments['mw_gen_2']['instance'].update({'frequency': self.settings['mw_pulses']['mw_2_frequency']})
        self.instruments['mw_gen_2']['instance'].update({'enable_output': True})

        self.instruments['afg_iq']['instance'].configure_iq(angle=0, pulse_length=self.settings['mw_pulses']['mw_2_pi_half_time'] + 30 * 2, angle_idle=0)
        for channel in [1, 2]:
            self.instruments['afg_iq']['instance'].update({'ch%i_enable' % channel: True})

    def _configure_instruments_start_of_sweep(self, phase_current):
        self.instruments['afg_iq']['instance'].configure_iq(angle=phase_current,
                                                            pulse_length=self.settings['mw_pulses']['mw_2_pi_half_time'] + 30 * 2 + 400,
                                                            angle_idle=self.settings['mw_phase']['idle_phase'])
        self.instruments['afg_iq']['instance'].align_phase()

    def _configure_param_array(self):
        # Contruct the frequency array and store it in a variable called 'mw_frequencies'. Despite the naming, it's just a list of parameters to be swept;
        # it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        self.params = np.linspace(self.settings['mw_phase']['min_phase'], self.settings['mw_phase']['max_phase'],
                                  self.settings['mw_phase']['phase_points'])

    def _add_mw_switch_to_sequences(self, pulse_sequences):
        return pulse_sequences


class RamseyFreqSweep(PulsedEsr):
    """
    Takes a Ramsey measurement with fixed duration of pi/2 pulse and phase accumulation time, while sweeping the frequency. This is to precisely find a
    detuning. For example, consider a 15N NV: with fast pi/2 pulses we will see Ramsey fringes due to the detunings from the two 15N hyperfine
    frequencies. If we fix the phase accumulation time to be 1/1.5MHz/2 = 330 ns, we can sweep the frequency of the pi/2 pulse until the readout signal is at
    an extrema, allowing us to precisely determine the frequency halfway between the two hyperfine frequencies (1.5 MHz is half the hyperfine splitting of 3 MHz).
    """

    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('pi_half_pulse_time', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('tau', 80, float, 'the time duration of the phase accumulation time in the Ramsey sequence (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('mw_generator_switching_time', .01, float,
                  'time wait after switching center frequencies on generator (s)')
    ]

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

        tau = self.settings['tau']
        pulse_sequences = []
        tau_list = [tau]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_i'
        pi_half_time = self.settings['pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', laser_off_time, nv_reset_time),
                 Pulse('apd_readout', laser_off_time + delay_readout, meas_time)]

            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence += [
                    Pulse(microwave_channel, laser_off_time + nv_reset_time + laser_off_time, pi_half_time)]
                pulse_sequence += [
                    Pulse(microwave_channel, laser_off_time + nv_reset_time + laser_off_time + pi_half_time + tau, pi_half_time)]

            pulse_sequence += [
                Pulse('laser',
                      laser_off_time + nv_reset_time + laser_off_time + pi_half_time + tau + pi_half_time + delay_mw_readout,
                      nv_reset_time),
                Pulse('apd_readout',
                      laser_off_time + nv_reset_time + laser_off_time + pi_half_time + tau + pi_half_time + delay_mw_readout + delay_readout,
                      meas_time)
            ]

            pulse_sequences.append(pulse_sequence)
        return pulse_sequences, tau_list, self.settings['read_out']['meas_time']