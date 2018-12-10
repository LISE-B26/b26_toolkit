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
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, MicrowaveGenerator, Pulse
from pylabcontrol.core import Parameter, Script
from b26_toolkit.data_processing.fit_functions import fit_exp_decay, exp_offset

class XY8_k(PulsedExperimentBaseScript): # ER 5.25.2017
    """
This script runs a Hahn echo on the NV to find the Hahn echo T2. To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pi pulses'),
            Parameter('microwave_channel_pi2', 'q', ['i', 'q'], 'Channel to use for the mw pi/2 pulses'),
            Parameter('pi_pulse_time_mwchan', 50.0, float, 'time duration of a pi pulse (in ns) for the mw_channel'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns)'),
            Parameter('pi_pulse_blocks_k', 1, int, 'number of pi pulse blocks of 8 in the XY8-k sequence')
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

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(XY8_k, self)._function(self.data)

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
        #tau_list = list(range(int(self.settings['tau_times']['min_time']), int(self.settings['tau_times']['max_time']),self.settings['tau_times']['time_step']))

        tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],
                             self.settings['tau_times']['time_step'])
        tau_list = np.ndarray.tolist(tau_list)  # 20180731 ER convert to list

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= 15]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
        microwave_channel_pi2 = 'microwave_' + self.settings['mw_pulses']['microwave_channel_pi2']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_time_mwchan = self.settings['mw_pulses']['pi_pulse_time_mwchan']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']


        for tau in tau_list:
            pulse_sequence = \
            [
                Pulse(microwave_channel_pi2, laser_off_time, pi_half_time), # pi/2 pulse
            ]

          #  next_pi_t = laser_off_time + pi_half_time/2. + tau/2 - pi_time/2.
            next_pi_t = laser_off_time + pi_half_time + tau/2 - pi_time/2.

            N = self.settings['mw_pulses']['pi_pulse_blocks_k']*8
            counter = 0
            for ind in range(0, N):
                if counter == 0 or counter == 2 or counter == 5 or counter == 7:
                    pulse_sequence += [
                        Pulse(microwave_channel_pi2, next_pi_t, pi_time) # pulses along x
                    ]
                else:
                    pulse_sequence += [
                        Pulse(microwave_channel, next_pi_t, pi_time_mwchan) # pulses along y
                    ]
                if counter == 7:
                    counter = -1
                next_pi_t = next_pi_t + tau
                counter += 1

            pulse_sequence += [
                Pulse(microwave_channel_pi2, next_pi_t - tau/2., pi_half_time)
                ]
            end_of_first_CPMG = next_pi_t - tau/2. + pi_half_time

            pulse_sequence += \
                [
                 Pulse('laser', end_of_first_CPMG + delay_mw_readout, nv_reset_time),
                 Pulse('apd_readout', end_of_first_CPMG + delay_mw_readout + delay_readout, meas_time)
                ]

            start_of_second_CPMG = end_of_first_CPMG + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence += \
            [
                Pulse(microwave_channel_pi2, start_of_second_CPMG, pi_half_time),
            ]

            next_pi_t = start_of_second_CPMG + pi_half_time + tau/2. - pi_time/2.
            counter = 0
            for ind in range(0, N):
                if counter == 0 or counter == 2 or counter == 5 or counter == 7:
                    pulse_sequence += [
                        Pulse(microwave_channel_pi2, next_pi_t, pi_time) # pulses along x
                    ]
                else:
                    pulse_sequence += [
                        Pulse(microwave_channel, next_pi_t, pi_time_mwchan) # pulses along y
                    ]
                if counter == 7:
                    counter = -1
                next_pi_t = next_pi_t + tau
                counter += 1

            pulse_sequence += [
                Pulse(microwave_channel_pi2, next_pi_t - tau/2., three_pi_half_time)
                ]

            end_of_second_CPMG = next_pi_t - tau/2. + three_pi_half_time

            pulse_sequence += [
                Pulse('laser', end_of_second_CPMG + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', end_of_second_CPMG + delay_mw_readout + delay_readout, meas_time)
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

        if 'fits'in data and data['fits'] is not None:
            counts = (-data['counts'][:,1] + data['counts'][:,0])/ (data['counts'][:,0] + data['counts'][:,1])
            tau = data['tau']
            fits = data['fits']

            axislist[0].plot(tau, counts, 'b')
         #   axislist[0].hold(True) ER 20181012

            axislist[0].plot(tau, exp_offset(tau, fits[0], fits[1], fits[2]))
            axislist[0].set_title('T2 decay time (simple exponential, p = 1): {:2.1f} ns'.format(fits[1]))
        else:
            super(XY8_k, self)._plot(axislist)
            axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9))
            axislist[0].legend(labels=('Ref Fluorescence', 'T2 Data'), fontsize=8)

class XY4(PulsedExperimentBaseScript):
    """
This script runs an XY4 pulse sequence.
todo(emma): (make double_init scheme)
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses',[
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            # Parameter('mw_switch_extra_time', 15, int, 'Time to add before and after microwave switch is turned on'),
            Parameter('pi_pulse_time', 50, float, 'time duration of pi-pulse (in ns)'),
            Parameter('number_of_pulse_blocks', 1, list(range(1, 17)), 'number of alternating x-y-x-y-y-x-y-x pulses'),
        ]),
        Parameter('tau_times',[
            Parameter('time_step', 5, [5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000],
                      'time step increment of time between pulses (in ns)'),
            Parameter('min_time', 100, float, 'minimum time between pulses (in ns)'),
            Parameter('max_time', 1000, float, 'maximum time between pulses (in ns)'),
        ]),
        Parameter('read_out',[
            Parameter('delay_mw_init', 1000, int, 'delay between initialization and mw (in ns)'),
            Parameter('delay_mw_readout', 200, int, 'delay between mw and readout (in ns)'),
            Parameter('meas_time', 250, float, 'measurement time after CPMG sequence (in ns)'),
            Parameter('nv_reset_time', 3000, int, 'time with laser on at the beginning to reset state'),
            Parameter('ref_meas_off_time', 1000, int,'laser off time before taking reference measurement at the end of init (ns)')
        ]),
        Parameter('num_averages', 1000, int, 'number of averages (should be less than a million)'),
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}
    _SCRIPTS = {}

    def _function(self):
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(XY4, self)._function()

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
        pulse_sequences = []
        # tau_list = range(int(max(15, self.settings['min_delay_time'])), int(self.settings['max_delay_time'] + 15),
        #                  self.settings['delay_time_step'])

        # JG: changed the previous because the 15ns is taken care of later
        tau_list = list(range(int(self.settings['tau_times']['min_time']),
                         int(self.settings['tau_times']['max_time']),
                         self.settings['tau_times']['time_step']))

        reset_time = self.settings['read_out']['nv_reset_time']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = pi_time/2.0

        ref_meas_off_time = self.settings['read_out']['ref_meas_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_init = self.settings['read_out']['delay_mw_init']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        number_of_pulse_blocks = self.settings['mw_pulses']['number_of_pulse_blocks']


        for tau in tau_list:

            pulse_sequence = []

            #initialize and pi/2 pulse
            pulse_sequence.extend([Pulse('laser', 0, reset_time - ref_meas_off_time - 15 - meas_time),
                                   Pulse('apd_readout', reset_time - 15 - meas_time, meas_time),
                                   Pulse('laser', reset_time - 15 - meas_time, meas_time),
                                   Pulse('microwave_i', reset_time + delay_mw_init, pi_half_time)
                                   ])

            #CPMG xyxyyxyx loops added number_of_pulse_blocks times
            section_begin_time = reset_time + delay_mw_init + pi_half_time - tau/2 #for the first pulse, only wait tau/2
            # JG 16-08-19 - begin changed to pi time instead of pi/2
            # section_begin_time = reset_time + delay_mw_init + pi_time
            # JG 16-08-19 - end

            for i in range(0, number_of_pulse_blocks):
                pulse_sequence.extend([Pulse('microwave_i', section_begin_time + 1*tau - pi_half_time, pi_time),
                                       Pulse('microwave_q', section_begin_time + 2*tau - pi_half_time, pi_time),
                                       Pulse('microwave_i', section_begin_time + 3*tau - pi_half_time, pi_time),
                                       Pulse('microwave_q', section_begin_time + 4*tau - pi_half_time, pi_time),
                                      ])
                section_begin_time += 4*tau

            #pi/2 and readout
            pulse_sequence.extend([Pulse('microwave_i', section_begin_time + tau/2, pi_half_time),
                                   Pulse('laser',       section_begin_time + tau/2 + pi_half_time + delay_mw_readout, meas_time),
                                   Pulse('apd_readout', section_begin_time + tau/2 + pi_half_time + delay_mw_readout, meas_time)])

            # JG 16-08-19 - begin changed to pi time instead of pi/2
            # pulse_sequence.extend([Pulse('microwave_i', section_begin_time + tau, pi_half_time),
            #                        Pulse('laser',       section_begin_time + tau + pi_half_time + delay_mw_readout, meas_time),
            #                        Pulse('apd_readout', section_begin_time + tau + pi_half_time + delay_mw_readout, meas_time)])
            # JG 16-08-19 - end


            pulse_sequences.append(pulse_sequence)


        # end_time_max = 0
        # for pulse_sequence in pulse_sequences:
        #     for pulse in pulse_sequence:
        #         end_time_max = max(end_time_max, pulse.start_time + pulse.duration)
        # for pulse_sequence in pulse_sequences:
        #     pulse_sequence.append(Pulse('laser', end_time_max + 1850, 15))

        return pulse_sequences, tau_list, meas_time


    def _plot(self, axislist, data = None):
        """
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional) dataset to plot (dictionary that contains keys counts, tau), if not provided use self.data
        """

        super(XY4, self)._plot(axislist, data)
        axislist[0].set_title('XY4')
        axislist[0].legend(labels=('Ref Fluorescence', 'XY4 data'), fontsize=8)

class XYXY(PulsedExperimentBaseScript): # ER 5.25.2017
    """
Pulse sequence is X Y X Y X Y X Y .... to accumulate pulse errors and calibrate phase - typically a very short tau should be used
Uses double_init scheme
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pi pulses'),
            Parameter('microwave_channel_pi2', 'q', ['i', 'q'], 'Channel to use for the mw pi/2 pulses'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns)'),
            Parameter('pi_pulse_blocks_k', 1, int, 'number of pi pulse blocks of 8 in the XY8-k sequence')
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

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(XYXY, self)._function(self.data)

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
        tau_list = list(range(int(self.settings['tau_times']['min_time']), int(self.settings['tau_times']['max_time']),self.settings['tau_times']['time_step']))

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= 15]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
        microwave_channel_pi2 = 'microwave_' + self.settings['mw_pulses']['microwave_channel_pi2']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']


        for tau in tau_list:
            pulse_sequence = \
            [
                Pulse(microwave_channel_pi2, laser_off_time, pi_half_time), # pi/2 pulse
            ]

            next_pi_t =  laser_off_time + pi_half_time/2. + tau/2 - pi_time/2.
            N = self.settings['mw_pulses']['pi_pulse_blocks_k']*8
            counter = 0
            for ind in range(0, N):
                print("counter")
                print(counter)
                if counter %2 == 0:
                    pulse_sequence += [
                        Pulse(microwave_channel, next_pi_t, pi_time) # pulses along x
                    ]
                else:
                    pulse_sequence += [
                        Pulse(microwave_channel_pi2, next_pi_t, pi_time) # pulses along y
                    ]
                next_pi_t = next_pi_t + tau
                counter += 1

            pulse_sequence += [
                Pulse(microwave_channel_pi2, next_pi_t - tau + pi_time/2. + tau/2 - pi_half_time/2., pi_half_time)
                ]
            end_of_first_CPMG = next_pi_t - tau + pi_time/2. + tau/2 - pi_half_time/2. + pi_half_time

            pulse_sequence += \
                [
                 Pulse('laser', end_of_first_CPMG + delay_mw_readout, nv_reset_time),
                 Pulse('apd_readout', end_of_first_CPMG + delay_mw_readout + delay_readout, meas_time)
                ]

            start_of_second_CPMG = end_of_first_CPMG + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence += \
            [
                Pulse(microwave_channel_pi2, start_of_second_CPMG, pi_half_time),
            ]

            next_pi_t =  start_of_second_CPMG + pi_half_time/2. + tau/2 - pi_time/2.
            counter = 0
            for ind in range(0, N):
                print("counter take 2")
                print(counter)
                if counter %2 == 0:
                    pulse_sequence += [
                        Pulse(microwave_channel, next_pi_t, pi_time) # pulses along x
                    ]
                else:
                    pulse_sequence += [
                        Pulse(microwave_channel_pi2, next_pi_t, pi_time) # pulses along y
                    ]
                next_pi_t = next_pi_t + tau
                counter += 1

            pulse_sequence += [
                Pulse(microwave_channel_pi2, next_pi_t - tau + pi_time/2. + tau/2 - three_pi_half_time/2., three_pi_half_time)
                ]

            end_of_second_CPMG = next_pi_t - tau + pi_time/2. + tau/2 - three_pi_half_time/2. + three_pi_half_time

            pulse_sequence += [
                Pulse('laser', end_of_second_CPMG + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', end_of_second_CPMG + delay_mw_readout + delay_readout, meas_time)
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
         #   axislist[0].hold(True) ER 20181012

            axislist[0].plot(tau, exp_offset(tau, fits[0], fits[1], fits[2]))
            axislist[0].set_title('T2 decay time (simple exponential, p = 1): {:2.1f} ns'.format(fits[1]))
        else:
            super(XYXY, self)._plot(axislist)
            axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9))
            axislist[0].legend(labels=('Ref Fluorescence', 'T2 Data'), fontsize=8)