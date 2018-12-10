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
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_base_script import PulsedExperimentBaseScript
from b26_toolkit.instruments import NI6259, B26PulseBlaster, MicrowaveGenerator, Pulse
from pylabcontrol.core import Parameter, Script
from pylabcontrol.scripts import SelectPoints
from b26_toolkit.data_processing.fit_functions import fit_exp_decay, exp_offset
from b26_toolkit.scripts import ESR
from .rabi import Rabi

class HahnEcho(PulsedExperimentBaseScript): # ER 5.25.2017
    """
This script runs a Hahn echo on the NV to find the Hahn echo T2. To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns)')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'minimum time between pi pulses'),
            Parameter('max_time', 10000, float, 'maximum time between pi pulses'),
            Parameter('time_step', 5, [2.5, 5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000, 500000],
                  'time step increment of time between pi pulses (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 1000, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True}) # ER 20181018
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(HahnEcho, self)._function(self.data)

        counts = (- self.data['counts'][:, 1] + self.data['counts'][:,0]) / (self.data['counts'][:,1] + self.data['counts'][:, 0])
        tau = self.data['tau']

        try:
            fits = fit_exp_decay(tau, counts, offset = True, verbose = True)
            self.data['fits'] = fits
        except:
            self.data['fits'] = None
            self.log('t2 fit failed')

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
        tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],self.settings['tau_times']['time_step'])
        tau_list = np.ndarray.tolist(tau_list) # 20180731 ER convert to list

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= 15]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = \
            [
                Pulse(microwave_channel, laser_off_time, pi_half_time),
                Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau - pi_time/2., pi_time),
                Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2., pi_half_time)
            ]

            end_of_first_HE = laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time

            pulse_sequence += [
                 Pulse('laser', end_of_first_HE + delay_mw_readout, nv_reset_time),
                 Pulse('apd_readout', end_of_first_HE + delay_mw_readout + delay_readout, meas_time),
                 ]

            start_of_second_HE = end_of_first_HE + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence += \
            [
                Pulse(microwave_channel, start_of_second_HE, pi_half_time),
                Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau - pi_time/2., pi_time),
                Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau + tau - pi_half_time/2., three_pi_half_time)
            ]

            end_of_second_HE = start_of_second_HE + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time

            pulse_sequence += [
                Pulse('laser', end_of_second_HE + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', end_of_second_HE + delay_mw_readout + delay_readout, meas_time)
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
            counts = (-data['counts'][:,1] + data['counts'][:,0])/ (data['counts'][:,0] + data['counts'][:,1])
            tau = data['tau']
            fits = data['fits']

            axislist[0].plot(tau, counts, 'b')
           # axislist[0].hold(True)

            axislist[0].plot(tau, exp_offset(tau, fits[0], fits[1], fits[2]))
            axislist[0].set_title('T2 decay time (simple exponential, p = 1): {:2.1f} ns'.format(fits[1]))
        else:
            super(HahnEcho, self)._plot(axislist)
            axislist[0].set_title('Hahn Echo mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9))
            axislist[0].legend(labels=('Ref Fluorescence', 'T2 Data'), fontsize=8)

class HahnEcho_bothIQ(PulsedExperimentBaseScript): # ER 20181013
    """
This script runs a Hahn echo on the NV to find the Hahn echo T2. To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time.
The difference between this script and HahnEcho is that the phases of the pulses should be xyx, instead of xxx; i.e., we use BOTH IQ channels now.
The channel 'microwave_channel' in the settings corresponds to the pi/2 and 3pi/2 pulses. The other one the user chooses will be for the pi pulse.
We do this to check if we can extend T2 with something like this, which may help mitigate pulse errors.

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pi/2 pulses (other one will be for the pi pulses)'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns)')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'minimum time between pi pulses'),
            Parameter('max_time', 10000, float, 'maximum time between pi pulses'),
            Parameter('time_step', 5, [2.5, 5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000, 500000],
                  'time step increment of time between pi pulses (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 1000, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(HahnEcho_bothIQ, self)._function(self.data)

        counts = (- self.data['counts'][:, 1] + self.data['counts'][:,0]) / (self.data['counts'][:,1] + self.data['counts'][:, 0])
        tau = self.data['tau']

        try:
            fits = fit_exp_decay(tau, counts, offset = True, verbose = True)
            self.data['fits'] = fits
        except:
            self.data['fits'] = None
            self.log('t2 fit failed')

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
        tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],self.settings['tau_times']['time_step'])
        tau_list = np.ndarray.tolist(tau_list) # 20180731 ER convert to list

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= 15]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
        if self.settings['mw_pulses']['microwave_channel'] == 'i':
            mw_chan_pi = 'q'
        else:
            mw_chan_pi = 'i'
        microwave_channel_pi = 'microwave_' + mw_chan_pi
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = \
            [
                Pulse(microwave_channel, laser_off_time, pi_half_time),
                Pulse(microwave_channel_pi, laser_off_time + pi_half_time + tau - pi_time/2., pi_time),
                Pulse(microwave_channel, laser_off_time + pi_half_time + tau + tau, pi_half_time)
            ]

            end_of_first_HE = laser_off_time + pi_half_time + tau + tau + pi_half_time

            pulse_sequence += [
                 Pulse('laser', end_of_first_HE + delay_mw_readout, nv_reset_time),
                 Pulse('apd_readout', end_of_first_HE + delay_mw_readout + delay_readout, meas_time),
                 ]

            start_of_second_HE = end_of_first_HE + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence += \
            [
                Pulse(microwave_channel, start_of_second_HE, pi_half_time),
                Pulse(microwave_channel_pi, start_of_second_HE + pi_half_time + tau - pi_time/2., pi_time),
                Pulse(microwave_channel, start_of_second_HE + pi_half_time + tau + tau, three_pi_half_time)
            ]

            end_of_second_HE = start_of_second_HE + pi_half_time + tau + tau + three_pi_half_time

            pulse_sequence += [
                Pulse('laser', end_of_second_HE + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', end_of_second_HE + delay_mw_readout + delay_readout, meas_time)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

class HahnEchoManyNVs(Script):
    _DEFAULT_SETTINGS = [
        Parameter('esr_peak', 'upper', ['upper', 'lower', 'both'], 'if ESR fits two peaks, defines which one to use')
    ]
    _INSTRUMENTS = {}
    _SCRIPTS = {'select_NVs': SelectPoints, 'ESR': ESR, 'Rabi': Rabi, 'HE': HahnEcho}

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
                print('running rabi')
                rabi = self.scripts['Rabi']
                rabi.settings['tag'] = 'rabi_NV' + str(num)
                rabi.settings['mw_pulses']['mw_frequency'] = float(freq)
                print('about to run rabi')
                rabi.run()
                rabi_fit = rabi.data['fits']
                if rabi_fit is None:
                    continue
                pi_time = abs((np.pi - rabi_fit[2])/rabi_fit[1])
                pi_time = min(max(np.round(pi_time / 2.5) * 2.5, 15.), 300.) #round to nearest 2.5 ns
                pi_half_time = min(max((np.pi / 2 - rabi_fit[2]) / rabi_fit[1], 15.), 300)
                three_pi_half_time = min(max((3 * np.pi / 2 - rabi_fit[2]) / rabi_fit[1], 15.), 300)
                find_NV_HE = self.scripts['HE'].scripts['find_nv']
                find_NV_HE.settings['initial_point']['x'] = find_NV_rabi.data['maximum_point']['x']
                find_NV_HE.settings['initial_point']['y'] = find_NV_rabi.data['maximum_point']['y']
                HE = self.scripts['HE']
                HE.settings['mw_pulses']['mw_frequency'] = float(freq)
                HE.settings['mw_pulses']['pi_time'] = float(pi_time)
                HE.settings['mw_pulses']['pi_half_time'] = float(pi_half_time)
                HE.settings['mw_pulses']['3pi_half_time'] = float(three_pi_half_time)
                HE.settings['tag'] = 'HE' + '_NV' + str(num)
                HE.run()

    def plot(self, figure_list):
        if self._current_subscript_stage is not None:
            if self._current_subscript_stage['current_subscript'] is not None:
                self._current_subscript_stage['current_subscript'].plot(figure_list)


    def skip_next(self):
        for script in self.scripts.values():
            script.stop()
