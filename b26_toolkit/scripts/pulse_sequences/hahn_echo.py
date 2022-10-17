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
from b26_toolkit.instruments import NI6259, NI9402, NI9263_02, B26PulseBlaster, MicrowaveGenerator, Pulse
from pylabcontrol.core import Parameter, Script
from pylabcontrol.scripts import SelectPoints
from b26_toolkit.data_processing.fit_functions import fit_exp_decay, exp_offset
from b26_toolkit.plotting.plots_1d import plot_1d_simple_timetrace_ns, plot_pulses, update_pulse_plot, update_1d_simple
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from b26_toolkit.scripts import ESR
from b26_toolkit.scripts.esr import ESR_tracking
from .rabi import Rabi
from b26_toolkit.scripts.autofocus import AutoFocusDAQ
from scipy import signal
import random


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
            Parameter('time_step', 5, float,
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
        Parameter('num_averages', 100000, int, 'number of averages')
    ]

    #041619 MM added cDAQ
    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

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

        for tau in tau_list:
            pulse_sequence = \
            [
                Pulse(microwave_channel, laser_off_time, pi_half_time),
                Pulse(microwave_channel, laser_off_time + pi_half_time + tau, pi_time),
                Pulse(microwave_channel, laser_off_time + pi_half_time + tau + pi_time + tau, pi_half_time)
            ]

            end_of_first_HE = laser_off_time + pi_half_time + tau + pi_time + tau + pi_half_time

            pulse_sequence += [
                 Pulse('laser', end_of_first_HE + delay_mw_readout, nv_reset_time),
                 Pulse('apd_readout', end_of_first_HE + delay_mw_readout + delay_readout, meas_time),
                 ]

            start_of_second_HE = end_of_first_HE + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence += \
            [
                Pulse(microwave_channel, start_of_second_HE, pi_half_time),
                Pulse(microwave_channel, start_of_second_HE + pi_half_time + tau, pi_time),
                Pulse(microwave_channel, start_of_second_HE + pi_half_time + tau + pi_time + tau, three_pi_half_time)
            ]

            #end_of_second_HE = start_of_second_HE + pi_half_time / 2. + tau + tau - pi_half_time / 2. + pi_half_time
            end_of_second_HE = start_of_second_HE + pi_half_time + tau + pi_time + tau + three_pi_half_time

            pulse_sequence += [
                Pulse('laser', end_of_second_HE + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', end_of_second_HE + delay_mw_readout + delay_readout, meas_time)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _create_pulse_sequences_old(self):
        '''
        The old sequence where tau is counted from the middle of the pulse
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

            #end_of_second_HE = start_of_second_HE + pi_half_time / 2. + tau + tau - pi_half_time / 2. + pi_half_time
            end_of_second_HE = start_of_second_HE + pi_half_time / 2. + tau + tau - pi_half_time / 2. + three_pi_half_time

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

class HahnEcho_AttoTrig_GARBAGE(HahnEcho):
    """
This script runs a Hahn echo on the NV to find the Hahn echo T2. To symmetrize the sequence between the 0 and +/-1 state
we reinitialize every time. Additionally, immediately after the first pi/2 pulse, it sends a trigger pulse to a function
generator which moves the Attocube.

This script was written on 20220707 by FF for measuring the NV T2 while a magnet is moved away from it, to see if the T2 is
preserved during movement in a high magnetic field gradient.

We never took (or tried to) good data with this, so this hasn't been tested well. Please delete after the first strings paper
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power (pulsed) in dBm'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns)')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'minimum time between pi pulses'),
            Parameter('max_time', 10000, float, 'maximum time between pi pulses'),
            Parameter('time_step', 5, float,
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
        Parameter('atto_trig', [
            Parameter('trig_duration', 200, float,
                      'length of trigger pulse to send to function generator controlling attocube')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

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
        tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],
                             self.settings['tau_times']['time_step'])
        tau_list = np.ndarray.tolist(tau_list)  # 20180731 ER convert to list

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

        atto_trig_duration = self.settings['atto_trig']['trig_duration']

        for tau in tau_list:
            pulse_sequence = \
                [
                    Pulse(microwave_channel, laser_off_time, pi_half_time),
                    Pulse(microwave_channel, laser_off_time + pi_half_time / 2. + tau - pi_time / 2., pi_time),
                    Pulse(microwave_channel, laser_off_time + pi_half_time / 2. + tau + tau - pi_half_time / 2.,
                          pi_half_time),
                    Pulse('atto_trig', laser_off_time + pi_half_time, atto_trig_duration)
                ]

            end_of_first_HE = laser_off_time + pi_half_time / 2. + tau + tau - pi_half_time / 2. + pi_half_time

            pulse_sequence += [
                Pulse('laser', end_of_first_HE + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', end_of_first_HE + delay_mw_readout + delay_readout, meas_time),
            ]

            start_of_second_HE = end_of_first_HE + delay_mw_readout + nv_reset_time + laser_off_time

            if start_of_second_HE < 100000:
                start_of_second_HE = 100000

            pulse_sequence += \
                [
                    Pulse(microwave_channel, start_of_second_HE, pi_half_time),
                    Pulse(microwave_channel, start_of_second_HE + pi_half_time / 2. + tau - pi_time / 2., pi_time),
                    Pulse(microwave_channel, start_of_second_HE + pi_half_time / 2. + tau + tau - pi_half_time / 2.,
                          three_pi_half_time),
                    Pulse('atto_trig', start_of_second_HE + pi_half_time, atto_trig_duration)
                ]

            end_of_second_HE = start_of_second_HE + pi_half_time / 2. + tau + tau - pi_half_time / 2. + pi_half_time

            pulse_sequence += [
                Pulse('laser', end_of_second_HE + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', end_of_second_HE + delay_mw_readout + delay_readout, meas_time)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

class FieldProfile_Pulsed_GARBAGE(HahnEcho_AttoTrig_GARBAGE):
    """
    This script initializes the NV, then moves the Attocube for a certain time by sending a trigger pulse to the Attocube
    controller. At tau during the Attocube's path, a MW pi pulse is sent. If the NV transition at this point in time
    coincides with the MW frequency, the NV will flip its state, and less so if it's off resonant.

    Essentially this is doing pulsed ODMR at various times during a triggered movement sequence.
    """

    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'minimum time between pi pulses'),
            Parameter('max_time', 10000, float, 'maximum time between pi pulses'),
            Parameter('time_step', 5, [2.5, 5, 10, 20, 40, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 30000,
                                       50000, 100000, 500000],
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
        Parameter('atto_trig', [
            Parameter('trig_duration', 200, float,
                      'length of trigger pulse in ns to send to function generator controlling attocube'),
            Parameter('settle_time', 200000, float,
                      'time in ns to wait before moving Attocube again')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'NI9263_02': NI9263_02, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        # Same as default Hahn echo, but we also setup the trigger for a function generator connected to a high-power piezo controller

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True}) # ER 20181018
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})

        """
        sample_rate = 100000  # Use max sampling rate
        daq_out = self.instruments['NI9263_02']['instance']
        daq_out.settings['analog_output']['ao0']['sample_rate'] = sample_rate  # hard-coded for now for testing
        trig_source = ('/cDAQ1Mod2/PFI2').encode('ascii')  # hard-coded for now for testing

        amplitude_daq = self.settings['waveform']['amplitude'] / self.settings['waveform']['amplification_factor']
        ao_task = daq_out.setup_AO_triggered(channels=['ao0'],
                                             waveform=self.triangular_waveform(sample_rate, amplitude_daq,
                                                                               self.settings['waveform'][
                                                                                   'total_duration']),
                                             clk_source="", trig_source=trig_source)
        daq_out.run(ao_task)
        """

        super(HahnEcho, self)._function(self.data)

        counts = (- self.data['counts'][:, 1] + self.data['counts'][:,0]) / (self.data['counts'][:,1] + self.data['counts'][:, 0])
        tau = self.data['tau']

        try:
            fits = fit_exp_decay(tau, counts, offset = True, verbose = True)
            self.data['fits'] = fits
        except:
            self.data['fits'] = None
            self.log('t2 fit failed')

    def triangular_waveform(self, sample_rate, amplitude, period):
        """
        Returns a single period of a triangular wave
        Args:
            sample_rate: sample rate of analog output
            amplitude: max voltage of triangular wave (min is 0)
            period: time in ns for voltage to ramp from 0 to max and then back to 0

        Returns: a single period of a triangular wave

        """
        period = period * 1e-9
        t_end = period
        t_array = np.linspace(0, t_end, int(t_end * sample_rate))
        waveform = amplitude * (signal.sawtooth(2 * np.pi * t_array / period, 0.5) + 1)
        print(waveform)

        return waveform

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
        tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],
                             self.settings['tau_times']['time_step'])
        tau_list = np.ndarray.tolist(tau_list)  # 20180731 ER convert to list

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        # MM: updated to min_pulse_dur
        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        atto_trig_duration = self.settings['atto_trig']['trig_duration']
        atto_movement_time = self.settings['atto_trig']['settle_time']
        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', atto_movement_time + laser_off_time, nv_reset_time),
                 Pulse('apd_readout', atto_movement_time + laser_off_time + delay_readout, meas_time),
                 Pulse('atto_trig', atto_movement_time + laser_off_time + nv_reset_time + laser_off_time, atto_trig_duration)
                 ]
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence += [Pulse(microwave_channel, atto_movement_time + laser_off_time + nv_reset_time + laser_off_time + tau, pi_time)]

            pulse_sequence += [
                Pulse('laser', atto_movement_time + laser_off_time + nv_reset_time + laser_off_time + delay_mw_readout + tau, nv_reset_time),
                Pulse('apd_readout', atto_movement_time + laser_off_time + nv_reset_time + laser_off_time + delay_mw_readout + delay_readout + tau, meas_time)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

class FieldProfile_CW(PulsedExperimentBaseScript):
    """
    This script initializes the NV, then moves the Attocube for a certain time by sending a trigger pulse to the Attocube
    controller. At tau during the Attocube's path, a MW pi pulse is sent. If the NV transition at this point in time
    coincides with the MW frequency, the NV will flip its state, and less so if it's off resonant.

    Essentially this is doing pulsed ODMR at various times during a triggered movement sequence.
    """

    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB')]),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 10, int, 'number of frequencies in scan'),
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'beginning of first readout window'),
            Parameter('max_time', 10000, float, 'beginning of last readout window'),
            Parameter('time_step', 1000, [20, 50, 80, 100, 200, 1000, 2000, 5000, 10000, 20000, 40000, 50000, 100000, 200000],
                      'time step between beginning of readout windows')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 1000, float, 'measurement time after rabi sequence (in ns)'),
        ]),
        Parameter('atto_trig', [
            Parameter('trig_duration', 200, float,
                      'length of trigger pulse in ns to send to function generator controlling attocube'),
            Parameter('settle_time', 200000, float,
                      'time in ns to wait before moving Attocube again')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
        Parameter('averaging_block', 8000, int, 'size of averaging block (plot will be updated after each block')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'NI9263_02': NI9263_02, 'PB': B26PulseBlaster,
                    'mw_gen': MicrowaveGenerator}


    def _function(self, in_data=None):
        """
        This is the core loop in which the desired experiment specified by the inheriting script's pulse sequence
        is performed.

        in_data: input data dictionary, caution 'tau' and 'counts' will be overwritten here!

        Poststate: self.data contains two key/value pairs, 'tau' and 'counts'
            'tau': a list of the times tau that are scanned over for the relative experiment (ex wait times between pulses)
             'counts': the counts received from the experiment. This is a len('tau') list, with each element being a list
             of length 1, 2, or 3 corresponding to the sum over all trials for a single tau time. In this sublist, the
             first value is the signal, the second (optional) value is counts in the |0> state (the maximum counts for
             normalization), and the third (optional) value is the counts in the |1> state (the minimum counts for
             normalization)

        """
        self.MAX_AVERAGES_PER_SCAN = self.settings['averaging_block']
        # retrieve initial mw carrier frequency to protect against bad fits in NV ESR tracking
        last_mw = self.scripts['esr'].instruments['microwave_generator']['instance'].frequency
        pulse_ampl = self.scripts['esr'].instruments['microwave_generator']['instance'].amplitude
        mod_flag = self.scripts['esr'].instruments['microwave_generator']['instance'].enable_modulation

        # Keeps track of index of current pulse sequence for plotting
        self.sequence_index = 0

        # self.is_valid and create pulses
        self.pulse_sequences, self.tau_list, self.measurement_gate_width = self.create_pulse_sequences()
        self.num_averages = self.settings['num_averages']

        if in_data is None:
            in_data = {}

        # calculates the number of daq reads per loop requested in the pulse sequence by asking how many apd reads are
        # called for. if this is not calculated properly, daq will either end too early (number too low) or hang since it
        # never receives the rest of the counts (number too high)
        num_daq_reads = 0

        for pulse in self.pulse_sequences[0]:
            if pulse.channel_id == 'apd_readout':
                num_daq_reads += 1
        self._initialize_data(num_daq_reads, in_data)

        self.freq_values, self.freq_range = self.get_freq_array()

        self.instruments['mw_gen']['instance'].update({'enable_modulation': False})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

        # divides the total number of averages requested into a number of slices of MAX_AVERAGES_PER_SCAN and a remainder.
        # This is required because the pulseblaster won't accept more than ~4E6 loops (22 bits available to store loop
        # number) so need to break it up into smaller chunks (use 1E6 so initial results display faster)
        (num_1E5_avg_pb_programs, remainder) = divmod(self.num_averages, self.MAX_AVERAGES_PER_SCAN)
        # run find_nv if tracking is on ER 5/30/2017
        if self.settings['Tracking']['on/off']:
            self.scripts['find_nv'].run()
            if self.scripts['find_nv'].data['fluorescence'] == 0.0: # if it doesn't find an NV, abort the experiment
                self.log('Could not find an NV in FindNV.')
                self._abort = True
                return  # exit function in case no NV is found

        self.log("Averaging over {0} blocks of 1e5".format(num_1E5_avg_pb_programs))
        for average_loop in range(int(num_1E5_avg_pb_programs)):
            self.log("Running average block {0} of {1}".format(average_loop+1, int(num_1E5_avg_pb_programs)))
            if self._abort:
                print('aborting!!')
                # ER 20200828 stop the pulseblaster
                if self.instruments['PB']['instance'].settings['PB_type'] == 'USB':
                    print('stopping pulse seq: abort!! ')
                    self.instruments['PB']['instance'].stop_pulse_seq()

                self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

                self.log('Aborted pulseblaster script during loop')
                break

            if self.settings['ESR_Tracking']['on/off'] and average_loop % self.settings['ESR_Tracking']['track_every_N']==0:
                self.scripts['esr'].run()

                # retrieve the new mw frequency: if there are two frequencies in the fit, pick the one closest to the old frequency
                fit_params = self.scripts['esr'].data['fit_params']

                # default update flag to false
                update_mw = False

                if fit_params is not None and len(fit_params) and fit_params[0] != -1:  # check if fit valid
                    if len(fit_params) == 4:
                        # single peak
                        if (fit_params[2] - last_mw)**2 < (self.settings['ESR_Tracking']['allowed_delta_freq']*1e6)**2: # check if new value is within range allowed
                            update_mw = True
                        new_mw = fit_params[2]
                    elif len(fit_params) == 6:
                        # double peak, don't update the frequency - the fit may be bad
                        update_mw = False

                if update_mw:
                    #self.instruments['mw_gen'].update({'frequency': new_mw})
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'frequency': float(new_mw)})
                    self.log('Updated mw carrier frequency to: {}'.format(new_mw))
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'amplitude': float(pulse_ampl)})
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'enable_modulation': bool(mod_flag)})

                    last_mw = new_mw
                else:
                    #self.instruments['mw_gen'].update({'frequency': last_mw})
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'frequency': float(last_mw)})
                    self.log('Not updating the mw carrier frequency. SRS carrier frequency kept at {0} Hz'.format(last_mw))
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'amplitude': float(pulse_ampl)})
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'enable_modulation': bool(mod_flag)})

            self.current_averages = (average_loop + 1) * self.MAX_AVERAGES_PER_SCAN

            self._run_sweep(self.pulse_sequences, self.MAX_AVERAGES_PER_SCAN, num_daq_reads)

        if remainder != 0 and not self._abort:
            self.current_averages = self.num_averages
            self._run_sweep(self.pulse_sequences, remainder, num_daq_reads)

        if (len(self.data['counts'][0]) == 1) and not self._abort:
            self.data['counts'] = np.array([item for sublist in self.data['counts'] for item in sublist])

        # save data on the fly so that we can start to analyze it while the experiment is running!
        if self.settings['save']:
        #     self.save_b26()
            self.save_data()
        #     self.save_log()
        #     self.save_image_to_disk()

    def get_freq_array(self):
        return ESR.get_freq_array(self)

    def _run_sweep(self, pulse_sequences, num_loops_sweep, num_daq_reads, verbose=False):
        """
        Each pulse sequence specified in pulse_sequences is run num_loops_sweep consecutive times.

        Args:
            pulse_sequences: a list of pulse sequences to run, each corresponding to a different value of tau. Each
                             sequence is a list of Pulse objects specifying a given pulse sequence
            num_loops_sweep: number of times to repeat each sequence before moving on to the next one
            num_daq_reads: number of times the daq must read for each sequence (generally 1, 2, or 3)

        Poststate: self.data['counts'] is updated with the acquired data

        """

        rand_indexes = list(range(len(pulse_sequences)))

        if self.settings['randomize']:
            # Randomizing order w/ which to run signal and normalizing sequences'
            random.shuffle(rand_indexes)

        if verbose:
            print(('_run_sweep number of pulse sequences', len(pulse_sequences)))

        for index, sequence in enumerate(pulse_sequences):
            print(('_run_sweep index', index, len(pulse_sequences)))
            if verbose:
                print(('_run_sweep index', index, len(pulse_sequences)))

            rand_index = rand_indexes[index]
            if self._abort:
                print('aborting in run sweep!!')

                if self.instruments['PB']['instance'].settings['PB_type'] == 'USB':
                    print('stopping pulse seq in run sweep: abort!! ')
                    self.instruments['PB']['instance'].stop_pulse_seq()

                break


            self.instruments['mw_gen']['instance'].update({'frequency': float(self.freq_values[rand_index])})
            result = self._run_single_sequence(pulse_sequences[rand_index], num_loops_sweep,
                                               num_daq_reads)  # keep entire array
            # self.result_current is added here for pulsed_esr. In a usual pulse sequence (e.g. Rabi), there are 2
            # loops: tau and averaging blocks. In pulsed_esr, there are 3 loops: MW frequency, tau (with 1 value), and
            # averaging blocks. Because of this different structure I added self.result_current so that I can do the
            # averaging in _function()
            self.result_current = self._normalize_to_kCounts(np.array(result), self.measurement_gate_width,
                                                             num_loops_sweep)

            # self.result_current = np.array(result)

            # for tracking ER 20210331
            counts_to_check = self._normalize_to_kCounts(np.array(result), self.measurement_gate_width,
                                                         num_loops_sweep)
            # ER 20210331
            if self.settings['save_full']:
                self.data['full_contrast'][int(self.current_averages / self.MAX_AVERAGES_PER_SCAN) - 1][rand_index] = \
                    (counts_to_check[1] - counts_to_check[0]) / np.mean(counts_to_check)

            if self.settings['normalize_block']:

                result = [1, result[1] / result[0]]

                # do the averaging with existing data
                self.count_data[rand_index] = (self.count_data[rand_index] * self.current_averages + np.array(
                    result) * num_loops_sweep) \
                                              / (self.current_averages + num_loops_sweep)
                self.data['counts'][rand_index] = self.count_data[rand_index]

            else:
                self.count_data[rand_index] = self.count_data[rand_index] + result
                self.data['counts'][rand_index] = self._normalize_to_kCounts(self.count_data[rand_index],
                                                                             self.measurement_gate_width,
                                                                             self.current_averages)

            print('Measurement gate width: %.2e'%self.measurement_gate_width)

            self.sequence_index = rand_index
            counts_temp = counts_to_check[0]
            # track to the NV if necessary ER 5/31/17
            if self.settings['Tracking']['on/off']:
                if (1 + (1 - self.settings['Tracking']['threshold'])) * self.settings['Tracking'][
                    'init_fluor'] < counts_temp or \
                        self.settings['Tracking']['threshold'] * self.settings['Tracking']['init_fluor'] > counts_temp:
                    if verbose:
                        print('TRACKING TO NV...')
                    self.scripts['find_nv'].run()
                    self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data[
                        'maximum_point']
                    # MM 20190621
            #                    self.scripts['autofocus'].scripts['take_image'].settings['point_a'] = self.scripts['find_nv'].data['maximum_point']
            self.updateProgress.emit(self._calc_progress(index))

    def _create_pulse_sequences(self):
        '''

        Returns: pulse_sequence, num_averages, tau_list, meas_time
            pulse_sequences: a single pulse sequences, with daq reads at specified taus
            num_averages: the number of times to repeat each pulse sequence

        '''

        if self.settings['read_out']['meas_time'] > self.settings['tau_times']['time_step']:
            self.log('Readout window must be shorter between spacing between beginning of readout windows!')
            raise AttributeError

        tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],
                             self.settings['tau_times']['time_step'])
        tau_list = np.ndarray.tolist(tau_list)  # 20180731 ER convert to list

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        meas_time = self.settings['read_out']['meas_time']

        atto_trig_duration = self.settings['atto_trig']['trig_duration']
        atto_movement_time = self.settings['atto_trig']['settle_time']

        pulse_sequence = [Pulse('atto_trig', atto_movement_time, atto_trig_duration)]
        pulse_sequence += [Pulse('laser', 0, atto_movement_time + tau_list[-1] + meas_time)]
        pulse_sequence += [Pulse('microwave_switch', 0, atto_movement_time + tau_list[-1] + meas_time)]


        for tau in tau_list:
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence += [Pulse('apd_readout', atto_movement_time + tau, meas_time)]

        self.tau_list = tau_list

        pulse_sequences = []

        self.freq_values, self.freq_range = self.get_freq_array()
        for freq in self.freq_values:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def create_pulse_sequences(self, logging=True):
        """
        A function to create the pulse sequence.

        NOTE THAT this is different from _create_pulse_sequences(), which must be overwritten in scripts inheriting from this script
        in contrast this function here takes the output from _create_pulse_sequences() and cleans it up!

        Returns:
            pulse_sequences: a list of pulse sequences to be run
            tau_list: the list of tau times we vary
            measurement_gate_width: measurement window time
        """

        pulse_sequences, tau_list, measurement_gate_width = self._create_pulse_sequences()

        # Adding microwave switch
        if self.settings['mw_switch']['add']:
            if logging:
                self.log('Adding microwave switch to pulse sequences')
            pulse_sequences = self._add_mw_switch_to_sequences(pulse_sequences)
        if 'rf_switch' in self.settings and self.settings['rf_switch']['add']:
            if logging:
                self.log('Adding rf switch to pulse sequences')
            pulse_sequences = self._add_rf_switch_to_sequences(pulse_sequences)

        # look for bad pulses, i.e. that don't comply with the requirements, e.g. given the pulse-blaster specs or
        # requiring that pulses don't overlap

        valid_pulse_sequences = []
        for pulse_sequence in pulse_sequences:
            if not self._is_bad_pulse_sequence(pulse_sequence):
                valid_pulse_sequences.append(pulse_sequence)
            else:
                print('invalid sequence is: ', pulse_sequence)

        if logging:
            if valid_pulse_sequences:
                self.log("All generated pulse sequences are valid. No tau times will be skipped.")
            else:
                self.log("The pulse sequences is invalid")

        return valid_pulse_sequences, tau_list, measurement_gate_width

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

        #        if self.scripts['find_nv'].is_running: #self.scripts['find_nv'].data['maximum_point']:
        #            # self.scripts['find_nv'].plot(axes_list)
        #             self.scripts['find_nv']._plot(axes_list)
        #             self._plot_refresh = True
        #        else:
        if data is None:
            data = self.data

        if 'counts' in data.keys():
            plot_pulses(axes_list[1], self.pulse_sequences[self.sequence_index])
            counts = np.array(self.data['counts'])
            if len(data['counts']) == 2:
                #plot_1d_simple_timetrace_ns(axes_list[0], data['tau'], [counts[0]/counts[1]], y_label= 'Normalized kcounts/s')
                plot_1d_simple_timetrace_ns(axes_list[0], data['tau'], counts)
            elif len(data['counts']) == 1:
                plot_1d_simple_timetrace_ns(axes_list[0], data['tau'], counts)
                plot_pulses(axes_list[1], self.pulse_sequences[self.sequence_index])
            else:
                label = ['CW ESR while moving', 'time', 'freq', 'kcounts/s']
                extent = data['tau'][0], data['tau'][-1], self.freq_values[-1], self.freq_values[0]
                aspect = np.abs((extent[1]-extent[0])/(extent[3]-extent[2]))
                plot_fluorescence_new(counts, extent, axes_list[0], max_counts=-1, labels=label, aspect=aspect)


    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''
        #        if self.scripts['find_nv'].is_running:
        #            self.scripts['find_nv']._update_plot(axes_list)
        #        else:
        counts = np.array(self.data['counts'])
        x_data = self.data['tau']
        axis1 = axes_list[0]
        if not counts == []:
            if len(counts) == 2:
                #update_1d_simple(axis1, x_data, [counts[0]/counts[1]])
                update_1d_simple(axis1, x_data, counts)
            elif len(counts) == 1:
                update_1d_simple(axis1, x_data, counts)
            else:
                update_fluorescence(counts, axes_list[0], -1)
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])

class CoilProfile_CW(FieldProfile_CW):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB')]),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 10, int, 'number of frequencies in scan'),
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'beginning of first readout window'),
            Parameter('max_time', 10000, float, 'beginning of last readout window'),
            Parameter('time_step', 1000,
                      [20, 50, 80, 100, 200, 1000, 2000, 5000, 10000, 20000, 40000, 50000, 100000, 200000],
                      'time step between beginning of readout windows')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 1000, float, 'measurement time after rabi sequence (in ns)'),
        ]),
        Parameter('atto_trig', [
            Parameter('trig_duration', 200, float,
                      'length of trigger pulse in ns to send to function generator controlling attocube'),
            Parameter('settle_time', 200000, float,
                      'time in ns to wait before moving Attocube again')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
        Parameter('averaging_block', 8000, int, 'size of averaging block (plot will be updated after each block')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'NI9263_02': NI9263_02, 'PB': B26PulseBlaster,
                    'mw_gen': MicrowaveGenerator}

    def _create_pulse_sequences(self):
        '''

        Returns: pulse_sequence, num_averages, tau_list, meas_time
            pulse_sequences: a single pulse sequences, with daq reads at specified taus
            num_averages: the number of times to repeat each pulse sequence

        '''

        if self.settings['read_out']['meas_time'] > self.settings['tau_times']['time_step']:
            self.log('Readout window must be shorter between spacing between beginning of readout windows!')
            raise AttributeError

        tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],
                             self.settings['tau_times']['time_step'])
        tau_list = np.ndarray.tolist(tau_list)  # 20180731 ER convert to list

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        meas_time = self.settings['read_out']['meas_time']

        atto_trig_duration = self.settings['atto_trig']['trig_duration']
        atto_movement_time = self.settings['atto_trig']['settle_time']

        pulse_sequence = [Pulse('rf_i', atto_movement_time - 60, atto_trig_duration + 60*2)]
        pulse_sequence += [Pulse('rf_switch', atto_movement_time, atto_trig_duration)]
        pulse_sequence += [Pulse('laser', 0, atto_movement_time + tau_list[-1] + meas_time)]
        pulse_sequence += [Pulse('microwave_switch', 0, atto_movement_time + tau_list[-1] + meas_time)]


        for tau in tau_list:
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence += [Pulse('apd_readout', atto_movement_time + tau, meas_time)]

        self.tau_list = tau_list

        pulse_sequences = []

        self.freq_values, self.freq_range = self.get_freq_array()
        for freq in self.freq_values:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

class FieldProfile_Pulsed(FieldProfile_CW):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)')]),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 10, int, 'number of frequencies in scan'),
        Parameter('tau_times', [
            Parameter('min_time', 500, float, 'beginning of first readout window'),
            Parameter('max_time', 10000, float, 'beginning of last readout window'),
            Parameter('time_step', 1000, [1000, 2000, 5000, 10000, 20000, 40000, 50000, 100000, 200000],
                      'time step between beginning of readout windows')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 1000, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 1000, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('atto_trig', [
            Parameter('trig_duration', 200, float,
                      'length of trigger pulse in ns to send to function generator controlling attocube'),
            Parameter('settle_time', 200000, float,
                      'time in ns to wait before moving Attocube again')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
        Parameter('averaging_block', 8000, int, 'size of averaging block (plot will be updated after each block')
    ]


    def _function(self, in_data=None):
        """
        This is the core loop in which the desired experiment specified by the inheriting script's pulse sequence
        is performed.

        in_data: input data dictionary, caution 'tau' and 'counts' will be overwritten here!

        Poststate: self.data contains two key/value pairs, 'tau' and 'counts'
            'tau': a list of the times tau that are scanned over for the relative experiment (ex wait times between pulses)
             'counts': the counts received from the experiment. This is a len('tau') list, with each element being a list
             of length 1, 2, or 3 corresponding to the sum over all trials for a single tau time. In this sublist, the
             first value is the signal, the second (optional) value is counts in the |0> state (the maximum counts for
             normalization), and the third (optional) value is the counts in the |1> state (the minimum counts for
             normalization)

        """
        self.MAX_AVERAGES_PER_SCAN = self.settings['averaging_block']
        # retrieve initial mw carrier frequency to protect against bad fits in NV ESR tracking
        last_mw = self.scripts['esr'].instruments['microwave_generator']['instance'].frequency
        pulse_ampl = self.scripts['esr'].instruments['microwave_generator']['instance'].amplitude
        mod_flag = self.scripts['esr'].instruments['microwave_generator']['instance'].enable_modulation

        # Keeps track of index of current pulse sequence for plotting
        self.sequence_index = 0

        # self.is_valid and create pulses
        self.pulse_sequences, self.tau_list, self.measurement_gate_width = self.create_pulse_sequences()
        self.num_averages = self.settings['num_averages']

        if in_data is None:
            in_data = {}

        # calculates the number of daq reads per loop requested in the pulse sequence by asking how many apd reads are
        # called for. if this is not calculated properly, daq will either end too early (number too low) or hang since it
        # never receives the rest of the counts (number too high)
        num_daq_reads = 0

        for pulse in self.pulse_sequences[0]:
            if pulse.channel_id == 'apd_readout':
                num_daq_reads += 1
        self._initialize_data(num_daq_reads, in_data)

        self.freq_values, self.freq_range = self.get_freq_array()

        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

        # divides the total number of averages requested into a number of slices of MAX_AVERAGES_PER_SCAN and a remainder.
        # This is required because the pulseblaster won't accept more than ~4E6 loops (22 bits available to store loop
        # number) so need to break it up into smaller chunks (use 1E6 so initial results display faster)
        (num_1E5_avg_pb_programs, remainder) = divmod(self.num_averages, self.MAX_AVERAGES_PER_SCAN)
        # run find_nv if tracking is on ER 5/30/2017
        if self.settings['Tracking']['on/off']:
            self.scripts['find_nv'].run()
            if self.scripts['find_nv'].data['fluorescence'] == 0.0: # if it doesn't find an NV, abort the experiment
                self.log('Could not find an NV in FindNV.')
                self._abort = True
                return  # exit function in case no NV is found

        self.log("Averaging over {0} blocks of 1e5".format(num_1E5_avg_pb_programs))
        for average_loop in range(int(num_1E5_avg_pb_programs)):
            self.log("Running average block {0} of {1}".format(average_loop+1, int(num_1E5_avg_pb_programs)))
            if self._abort:
                print('aborting!!')
                # ER 20200828 stop the pulseblaster
                if self.instruments['PB']['instance'].settings['PB_type'] == 'USB':
                    print('stopping pulse seq: abort!! ')
                    self.instruments['PB']['instance'].stop_pulse_seq()

                self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

                self.log('Aborted pulseblaster script during loop')
                break

            if self.settings['ESR_Tracking']['on/off'] and average_loop % self.settings['ESR_Tracking']['track_every_N']==0:
                self.scripts['esr'].run()

                # retrieve the new mw frequency: if there are two frequencies in the fit, pick the one closest to the old frequency
                fit_params = self.scripts['esr'].data['fit_params']

                # default update flag to false
                update_mw = False

                if fit_params is not None and len(fit_params) and fit_params[0] != -1:  # check if fit valid
                    if len(fit_params) == 4:
                        # single peak
                        if (fit_params[2] - last_mw)**2 < (self.settings['ESR_Tracking']['allowed_delta_freq']*1e6)**2: # check if new value is within range allowed
                            update_mw = True
                        new_mw = fit_params[2]
                    elif len(fit_params) == 6:
                        # double peak, don't update the frequency - the fit may be bad
                        update_mw = False

                if update_mw:
                    #self.instruments['mw_gen'].update({'frequency': new_mw})
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'frequency': float(new_mw)})
                    self.log('Updated mw carrier frequency to: {}'.format(new_mw))
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'amplitude': float(pulse_ampl)})
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'enable_modulation': bool(mod_flag)})

                    last_mw = new_mw
                else:
                    #self.instruments['mw_gen'].update({'frequency': last_mw})
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'frequency': float(last_mw)})
                    self.log('Not updating the mw carrier frequency. SRS carrier frequency kept at {0} Hz'.format(last_mw))
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'amplitude': float(pulse_ampl)})
                    self.scripts['esr'].instruments['microwave_generator']['instance'].update({'enable_modulation': bool(mod_flag)})

            self.current_averages = (average_loop + 1) * self.MAX_AVERAGES_PER_SCAN

            self._run_sweep(self.pulse_sequences, self.MAX_AVERAGES_PER_SCAN, num_daq_reads)

        if remainder != 0 and not self._abort:
            self.current_averages = self.num_averages
            self._run_sweep(self.pulse_sequences, remainder, num_daq_reads)

        if (len(self.data['counts'][0]) == 1) and not self._abort:
            self.data['counts'] = np.array([item for sublist in self.data['counts'] for item in sublist])

        # save data on the fly so that we can start to analyze it while the experiment is running!
        if self.settings['save']:
        #     self.save_b26()
            self.save_data()
        #     self.save_log()
        #     self.save_image_to_disk()

    def _create_pulse_sequences(self):
        '''

        Returns: pulse_sequence, num_averages, tau_list, meas_time
            pulse_sequences: a single pulse sequences, with daq reads at specified taus
            num_averages: the number of times to repeat each pulse sequence

        '''

        if self.settings['read_out']['meas_time'] > self.settings['tau_times']['time_step']:
            self.log('Readout window must be shorter between spacing between beginning of readout windows!')
            raise AttributeError

        tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],
                             self.settings['tau_times']['time_step'])
        tau_list = np.ndarray.tolist(tau_list)  # 20180731 ER convert to list

        min_pulse_dur = self.instruments['PB']['instance'].settings['min_pulse_dur']
        tau_list = [x for x in tau_list if x == 0 or x >= min_pulse_dur]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        atto_trig_duration = self.settings['atto_trig']['trig_duration']
        atto_movement_time = self.settings['atto_trig']['settle_time']

        pulse_sequence = [Pulse('atto_trig', atto_movement_time, atto_trig_duration)]

        for tau in tau_list:
            if tau > 0:
                block_beginning = atto_movement_time + tau
                pulse_sequence += [Pulse(microwave_channel, block_beginning-30, pi_time + 30*2)]
                pulse_sequence += [Pulse('microwave_switch', block_beginning, pi_time)]
                pulse_sequence += [Pulse('laser', block_beginning + pi_time + delay_mw_readout, meas_time + nv_reset_time)]
                pulse_sequence += [Pulse('apd_readout', block_beginning + pi_time + delay_mw_readout + delay_readout, meas_time)]


        self.tau_list = tau_list

        pulse_sequences = []

        self.freq_values, self.freq_range = self.get_freq_array()
        for freq in self.freq_values:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

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

    _INSTRUMENTS = {'NI6259': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

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
    '''
    Assumes tracking in all subscripts.
    '''
    #Updated 072919 MM: for various clockspeeds, using ESR w/ find_nv tracking.
    _DEFAULT_SETTINGS = [
        Parameter('esr_peak', 'upper', ['upper', 'lower', 'both'], 'if ESR fits two peaks, defines which one to use'),
        Parameter('default_mw_frequency', 2.87e9, float, 'microwave frequency in Hz to use if no peak found on ESR'),
        Parameter('use_autofocus', False, bool, 'run autofocus in between every nv'),
        Parameter('default_continuum', 15, float, 'continuum counts in kcounts/sec to use if esr fit fails'),
        Parameter('rabi_attempts', 3, int, 'number of times to attempt Rabi fit before making rough guess.')
    ]
    _INSTRUMENTS = {}
    _SCRIPTS = {'select_NVs': SelectPoints, 'ESR': ESR_tracking, 'Rabi': Rabi, 'HE': HahnEcho, 'autofocus': AutoFocusDAQ}

    def _function(self):
        #To keep track of overall shift in galvoscan position between nvs.
        shift_x = 0
        shift_y = 0
        for num, nv_loc in enumerate(self.scripts['select_NVs'].data['nv_locations']):
            if self._abort:
                break


            #Update position to NV location + galvoscan shift; autofocus.
            self.scripts['ESR'].settings['tag'] = 'esr_NV' + str(num)
            find_NV_ESR = self.scripts['ESR'].scripts['find_nv']
            find_NV_ESR.settings['initial_point']['x'] = nv_loc[0] + shift_x
            find_NV_ESR.settings['initial_point']['y'] = nv_loc[1] + shift_y
            find_NV_ESR.run()
            if self.settings['use_autofocus']:
                autofocus = self.scripts['autofocus']
                autofocus.scripts['take_image'].settings['point_a'] = find_NV_ESR.data['maximum_point']
                autofocus.run()

            #ESR section: run, and extract frequency peaks from fit params.
            #Note: ESR starts by running find_nv, so no need to re-run before ESR called.
            #Making sure mw switch is turned on:
            self.scripts['Rabi'].instruments['PB']['instance'].update({'microwave_switch': {'status': True}})
            self.scripts['ESR'].run()
            fit_params = self.scripts['ESR'].data['fit_params']
            if fit_params is None:
                freqs = [self.settings['default_mw_frequency']]
                continuum = self.settings['default_continuum']
            elif len(fit_params) == 4:
                freqs = [fit_params[2]]
                continuum = self.settings['default_continuum']
            elif len(fit_params == 6):
                continuum = fit_params[0]
                if self.settings['esr_peak'] == 'lower':
                    freqs = [fit_params[4]]
                elif self.settings['esr_peak'] == 'upper':
                    freqs = [fit_params[5]]
                elif self.settings['esr_peak'] == 'both':
                    freqs = [fit_params[4], fit_params[5]]
            for freq in freqs:
                if self._abort:
                    break

                #Rabi section:

                #Update Rabi location and run findnv
                good_rabi = False
                tries_remaining = self.settings['rabi_attempts']
                find_NV_rabi = self.scripts['Rabi'].scripts['find_nv']
                find_NV_rabi.settings['initial_point'] = self.scripts['ESR'].scripts['find_nv'].data['maximum_point']

                #Try rabi fit 3 times.
                while (not good_rabi) and (tries_remaining > 0):
                    find_NV_rabi.run()
                    #Run Rabi and extract pulse durations:
                    print('running rabi')
                    rabi = self.scripts['Rabi']
                    rabi.settings['tag'] = 'rabi_NV' + str(num)
                    rabi.settings['mw_pulses']['mw_frequency'] = float(freq)
                    rabi.settings['Tracking']['init_fluor'] = continuum
                    print('about to run rabi')
                    rabi.run()
                    fits = rabi.data['fits']
                    tries_remaining -= 1
                    if fits is not None:
                        clockspeed = rabi.instruments['PB']['instance'].settings['clock_speed']
                        pb_round = 1 / clockspeed * 1e3
                        # Start w/ raw pi_time, pi/2 time, and round s.t. diff is divisible by 2 clock periods
                        pi_time = round((np.pi - fits[2]) / fits[1] / pb_round) * pb_round
                        pi_half_time = (np.pi / 2 - fits[2]) / fits[1]
                        pi_half_diff = round((pi_time - pi_half_time) / (2 * pb_round)) * 2 * pb_round
                        pi_half_time = pi_time - pi_half_diff
                        # Rounding 3/2 pi time:
                        three_pi_half_time = round((3 * np.pi / 2 - fits[2]) / fits[1] / pb_round) * pb_round
                        if pi_half_time > 12:
                            good_rabi = True
                        elif tries_remaining == 0:
                            norm_counts = rabi.data['counts'][:, 1] / rabi.data['counts'][:, 0]
                            tau = rabi.data['tau']

                            from scipy.signal import argrelmin
                            from scipy.optimize import curve_fit

                            #Find initial guess
                            min_ix = argrelmin(norm_counts)[0][0]
                            A_guess = (1 - norm_counts[min_ix]) / 2
                            p0 = [A_guess, np.pi / tau[min_ix], 0, 1 - A_guess]

                            #Fit
                            cos_fxn = lambda t, A, w, phi, c: A * np.cos(w * t + phi) + c
                            lb = [0, 0, 0, 0]
                            ub = [1, np.inf, np.pi, 1]
                            popt, pcov = curve_fit(cos_fxn, tau, norm_counts, p0=p0, bounds=(lb, ub))

                            #Extract pi time info
                            clockspeed = rabi.instruments['PB']['instance'].settings['clock_speed']
                            pb_round = 1 / clockspeed * 1e3

                            pi_time = round((np.pi - popt[2]) / popt[1] / pb_round)*pb_round
                            pi_half_time = (np.pi - 2 * popt[2]) / (2 * popt[1])
                            pi_half_diff = round((pi_time - pi_half_time) / (2 * pb_round)) * 2 * pb_round
                            pi_half_time = pi_time - pi_half_diff
                            three_pi_half_time = round((3 * np.pi - 2 * popt[2]) / (2 * popt[1])/pb_round)*pb_round

                            if pi_half_time > 12:
                                good_rabi = True

                if good_rabi:
                    #If successful Rabi fit, proceed to HE
                    #HE section:
                    find_NV_HE = self.scripts['HE'].scripts['find_nv']
                    find_NV_HE.settings['initial_point'] = find_NV_rabi.data['maximum_point']
                    find_NV_HE.run()
                    HE = self.scripts['HE']
                    HE.settings['mw_pulses']['mw_frequency'] = float(freq)
                    HE.settings['mw_pulses']['pi_pulse_time'] = float(pi_time)
                    HE.settings['mw_pulses']['pi_half_pulse_time'] = float(pi_half_time)
                    HE.settings['mw_pulses']['3pi_half_pulse_time'] = float(three_pi_half_time)
                    HE.settings['tag'] = 'HE' + '_NV' + str(num)
                    HE.settings['Tracking']['init_fluor'] = continuum
                    HE.run()

                    #update shifts to track long-term galvo drift
                    shift_x = find_NV_HE.data['maximum_point']['x'] - nv_loc[0]
                    shift_y = find_NV_HE.data['maximum_point']['y'] - nv_loc[1]
    def plot(self, figure_list):
        if self._current_subscript_stage is not None:
            if self._current_subscript_stage['current_subscript'] is not None:
                self._current_subscript_stage['current_subscript'].plot(figure_list)


    def skip_next(self):
        for script in self.scripts.values():
            script.stop()
