"""
    This file is part of b26_toolkit, a PyLabControl add-on for experiments in Harvard LISE B26.
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
from b26_toolkit.src.scripts.pulse_blaster_base_script import PulseBlasterBaseScript
from b26_toolkit.src.instruments import NI6259, B26PulseBlaster, MicrowaveGenerator, Pulse
from b26_toolkit.src.plotting.plots_1d import plot_esr, plot_pulses
from PyLabControl.src.core import Parameter

class PulsedESR_double_init(PulseBlasterBaseScript): # ER 20170616 - wrote for symmetry between 0 and -1 state
    """
This script applies a microwave pulse at fixed power and durations for varying frequencies
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
            super(PulsedESR_double_init, self)._function(self.data)
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

        return pulse_sequences, self.settings['num_averages'], tau_list, self.settings['read_out']['meas_time']




class PulsedESR(PulseBlasterBaseScript):
    """
This script applies a microwave pulse at fixed power and durations for varying frequencies
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('tau_mw', 200, float, 'the time duration of the microwaves (in ns)'),
        Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('reset_time', 1000000, int, 'time with laser on at the beginning to reset state'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
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

        reset_time = self.settings['reset_time']
        tau = self.settings['tau_mw']
        pulse_sequences = [[Pulse('laser', 0, reset_time),
                            Pulse('microwave_i', reset_time, tau),
                            Pulse('laser', reset_time + tau, self.settings['meas_time']),
                            Pulse('apd_readout', reset_time + tau, self.settings['meas_time'])
                            ]]

        tau_list = [tau]
        end_time_max = 0
        for pulse_sequence in pulse_sequences:
            for pulse in pulse_sequence:
                end_time_max = max(end_time_max, pulse.start_time + pulse.duration)
        for pulse_sequence in pulse_sequences:
            pulse_sequence.append(Pulse('laser', end_time_max + 1850, 15))

        return pulse_sequences, self.settings['num_averages'], tau_list, self.settings['meas_time']


class PulsedESRSlow(PulseBlasterBaseScript):
    """
This script applies a microwave pulse at fixed power and durations for varying frequencies.
This is the CW version, where we apply the MW only for short times but still much longer than a pi/2 pulse, ie. a few micro seconds to avoid heating of the sample.
This is different from the actual pulsed ESR, where we apply pi/2 pulses to get the max contrast.
    """


    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
            Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
            Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
        ]),
        Parameter('read_out', [
            Parameter('mw_off_time', 4000, int, '1.) Time the MWs are off (us) and the measurement time. Laser is on and photons are counted during this time.'),
            Parameter('mw_on_time', 250, float, '2.) Time the MWs are on (us). If measure_ref is True Laser is on and photons are counted during time of duration mw_off_time (i.e. the measurement time)!!.'),
            Parameter('laser_off_time', 250, float, '3.) Laser is off during this time after mW_on and mw_off pulse (us). Set to zero to skip this.'),
            Parameter('measure_ref', True, bool, 'If true take reference measurement. In that case mw_on_time has to be larger than mw_off_time. If not mw_on_time is set equal to mw_off_time'),
            Parameter('delay_mw_readout', 100, int, 'delay between laser on and readout (in ns)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

        assert self.settings['mw_pulses']['freq_start'] < self.settings['mw_pulses']['freq_stop']

        self.data = {'mw_frequencies': np.linspace(self.settings['mw_pulses']['freq_start'], self.settings['mw_pulses']['freq_stop'],
                                                   self.settings['mw_pulses']['freq_points']), 'esr_counts': []}

        for i, mw_frequency in enumerate(self.data['mw_frequencies']):
            self._loop_count = i
            self.instruments['mw_gen']['instance'].update({'frequency': float(mw_frequency)})
            super(PulsedESRSlow, self)._function(self.data)

            self.data['esr_counts'].append(self.data['counts'][0])

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
        esr_counts = np.array(data['esr_counts'])

        if len(np.shape(esr_counts))== 3:
            esr_counts  = esr_counts[:,0,0]/esr_counts[:,0,1]

        # print('sXXadsdasda', np.shape(esr_counts)[-1])
        #
        # # if there is two measurement per run, the second serves as a normalization measurement
        # if len(esr_counts.T) == 2:
        #     print('sadsdasda', np.shape(esr_counts))
        #     esr_counts = esr_counts[:,0] / esr_counts[:,1]
        #     print('ggggg', len(esr_counts))

        axis1 = axes_list[0]
        if not esr_counts == []:
            counts = esr_counts
            plot_esr(axis1, mw_frequencies[0:len(counts)], counts)
            axis1.hold(False)
            # axis1.set_title('avrg count')
        axis2 = axes_list[1]
        plot_pulses(axis2, self.pulse_sequences[0])


    def _update_plot(self, axes_list):
        mw_frequencies = self.data['mw_frequencies']
        esr_counts = np.array(self.data['esr_counts'])


        # if there is two measurement per run, the second serves as a normalization measurement
        if len(np.shape(esr_counts))== 3:
            esr_counts  = esr_counts[:,0,0]/esr_counts[:,0,1]


        print((len(esr_counts)))
        axis1 = axes_list[0]
        if not esr_counts == []:
            counts = esr_counts
            plot_esr(axis1, mw_frequencies[0:len(counts)], counts)
            axis1.hold(False)
            # axis2 = axes_list[1]
            # update_pulse_plot(axis2, self.pulse_sequences[0])

    def _create_pulse_sequences(self):

        '''

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement(s)
        '''


        # on contrast to other script the following times are given in us and have to be converted to ns
        mw_on_time = self.settings['read_out']['mw_on_time']*1e3
        mw_off_time = self.settings['read_out']['mw_off_time']*1e3
        laser_off_time = self.settings['read_out']['laser_off_time']*1e3

        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        measure_ref = self.settings['read_out']['measure_ref']

        # minimum pulse length is 15ns, if less set to zero, i.e. don't turn laser off
        if laser_off_time <= 15:
            laser_off_time = 0
        # enforce minimum pulse length (15ns)
        if delay_mw_readout <=15:
            delay_mw_readout = 15

        # if measuring the reference fluorescence duing the MW off time enforce that the MW are off at least as long as it takes to take the measurement
        if (measure_ref is True) and (mw_off_time < mw_on_time):
            mw_off_time = mw_on_time

        # todo: JG - this is a quick fix and should be handled by pulse_blaster_script.validate, which adds a pulse for the microwave switch, this has a hard coded delay of 40ns /
        # since here the mw pulse is the first pulse this results in negative times
        mw_offset_time = 40


        if laser_off_time == 0:
            if (measure_ref is True):
                pulse_sequences = [[Pulse('laser', mw_offset_time+ 0, mw_on_time+mw_off_time+2*delay_mw_readout),
                                    Pulse('microwave_i', mw_offset_time+ 0, mw_on_time+delay_mw_readout),
                                    Pulse('apd_readout', mw_offset_time+ delay_mw_readout, mw_on_time),
                                    Pulse('apd_readout', mw_offset_time+ 2*delay_mw_readout+mw_on_time, mw_on_time), # the readout is actually on mw_on_time long but the mw are mw_off_time off
                                    Pulse('off_channel', mw_offset_time+ 2 * delay_mw_readout + mw_on_time + mw_off_time, 15)
                                    # Pulse('laser', mw_offset_time+ 2 * delay_mw_readout + mw_on_time + mw_off_time,15) # at the end we want mw and laser to be off, so we add this short laser pulse because 'off' doesn't exist
                                    ]]
            else:
                pulse_sequences = [[Pulse('laser', mw_offset_time+ 0, mw_on_time+2*delay_mw_readout),
                                    Pulse('microwave_i', mw_offset_time+ 0, mw_on_time+delay_mw_readout),
                                    Pulse('apd_readout', mw_offset_time+ delay_mw_readout, mw_on_time),
                                    Pulse('off_channel', mw_offset_time+ 2 * delay_mw_readout + mw_on_time + mw_off_time, 15)
                                    # Pulse('laser', mw_offset_time+ 2 * delay_mw_readout + mw_on_time + mw_off_time,15) # at the end we want mw and laser to be off, so we add this short laser pulse because 'off' doesn't exist
                                    ]]
        else:
            if (measure_ref is True):
                pulse_sequences = [[Pulse('laser', mw_offset_time+ 0, mw_on_time+delay_mw_readout),
                                    Pulse('microwave_i',mw_offset_time+  0, mw_on_time+delay_mw_readout),
                                    Pulse('apd_readout',mw_offset_time+  delay_mw_readout, mw_on_time),
                                    Pulse('laser', mw_offset_time+ mw_on_time + delay_mw_readout + laser_off_time, mw_off_time+delay_mw_readout),
                                    Pulse('apd_readout', mw_offset_time+ 2*delay_mw_readout+mw_on_time+ laser_off_time, mw_on_time), # the readout is actually on mw_on_time long but the mw are mw_off_time off
                                    Pulse('off_channel', mw_offset_time+ 2 * delay_mw_readout + mw_on_time + 2 * laser_off_time + mw_off_time,15)
                                    # at the end we want mw and laser to be off, so we add this short laser pulse because 'off' doesn't exist
                                    # Pulse('laser', 2*delay_mw_readout+mw_on_time+ 2*laser_off_time +mw_off_time, 15) # at the end we want mw and laser to be off, so we add this short laser pulse because 'off' doesn't exist
                                    ]]
            else:
                pulse_sequences = [[Pulse('laser', mw_offset_time+ 0, mw_on_time+delay_mw_readout),
                                    Pulse('microwave_i',mw_offset_time+  0, mw_on_time+delay_mw_readout),
                                    Pulse('apd_readout',mw_offset_time+  delay_mw_readout, mw_on_time),
                                    Pulse('off_channel', mw_offset_time+ 2 * delay_mw_readout + mw_on_time + 2 * laser_off_time + mw_off_time,15)
                                    # at the end we want mw and laser to be off, so we add this short laser pulse because 'off' doesn't exist
                                    # Pulse('laser', 2*delay_mw_readout+mw_on_time+ 2*laser_off_time +mw_off_time, 15) # at the end we want mw and laser to be off, so we add this short laser pulse because 'off' doesn't exist
                                    ]]


        tau_list = [mw_on_time]
        # end_time_max = 0
        # for pulse_sequence in pulse_sequences:
        #     for pulse in pulse_sequence:
        #         end_time_max = max(end_time_max, pulse.start_time + pulse.duration)
        # for pulse_sequence in pulse_sequences:
        #     pulse_sequence.append(Pulse('laser', end_time_max + 1850, 15))

        return pulse_sequences, self.settings['num_averages'], tau_list, mw_on_time



class Pulsed_ESR_Pulsed_Laser(PulseBlasterBaseScript):
    """
This script applies a microwave pulse at fixed power for varying durations to measure Rabi Oscillations
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
        Parameter('delay_until_mw', 100, float, 'total time of rabi oscillations (in ns)'),
        Parameter('mw_duration', 200, float, 'total time of rabi oscillations (in ns)'),
        Parameter('time_step', 15, float,
                  'time step increment of rabi pulse duration (in ns)'),
        Parameter('time', 400, float, 'total time of rabi oscillations (in ns)'),
        Parameter('meas_time', 15, float, 'measurement time after rabi sequence (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('reset_time', 1000000, int, 'time with laser on at the beginning to reset state')
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_frequency']})
        super(Pulsed_ESR_Pulsed_Laser, self)._function()

    def _create_pulse_sequences(self):
        """
        Returns:
            pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """

        pulse_sequences = []
        tau_list = list(range(int(max(15, self.settings['time_step'])), int(self.settings['time'] + 15),
                         int(self.settings['time_step'])))
        reset_time = self.settings['reset_time']
        for tau in tau_list:
            if tau < self.settings['delay_until_mw']:
                pulse_sequences.append([Pulse('laser', 0, reset_time),
                                        Pulse('microwave_i', reset_time + self.settings['delay_until_mw'],
                                              self.settings['mw_duration']),
                                        Pulse('apd_readout', reset_time + tau, self.settings['meas_time'])
                                        ])
            else:
                pulse_sequences.append([Pulse('laser', 0, reset_time + max(tau + self.settings['meas_time'],
                                                                           self.settings['delay_until_mw'] +
                                                                           self.settings[
                                                                               'mw_duration'])),
                                        Pulse('microwave_i', reset_time + self.settings['delay_until_mw'],
                                              self.settings['mw_duration']),
                                        Pulse('apd_readout', reset_time + tau, self.settings['meas_time'])
                                        ])
        end_time_max = 0
        for pulse_sequence in pulse_sequences:
            for pulse in pulse_sequence:
                end_time_max = max(end_time_max, pulse.start_time + pulse.duration)
        for pulse_sequence in pulse_sequences:
            pulse_sequence[0] = Pulse('laser', 0, end_time_max)

        return pulse_sequences, self.settings['num_averages'], tau_list, self.settings['meas_time']

class ESRSingleFreqCont(PulseBlasterBaseScript):
    """
This script applies a microwave pulse at fixed power and durations for varying frequencies.
This is the CW version, where we apply the MW only for short times but still much longer than a pi/2 pulse, ie. a few micro seconds to avoid heating of the sample.
This is different from the actual pulsed ESR, where we apply pi/2 pulses to get the max contrast.
    """


    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('frequency', 2.82e9, float, 'frequency (Hz)')
        ]),
        Parameter('read_out', [
            Parameter('integration_time', 4000, int, '1.) Time the MWs are off (us) and the measurement time. Laser is on and photons are counted during this time.'),
            Parameter('delay_mw_readout', 100, int, 'delay between laser on and readout (in ns)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
        Parameter('max_points', 100, int, 'number of points to display if 0 show all')
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})

        self.data = {'mw_frequencies': np.linspace(self.settings['mw_pulses']['freq_start'], self.settings['mw_pulses']['freq_stop'],
                                                   self.settings['mw_pulses']['freq_points']), 'esr_counts': []}

        self.instruments['mw_gen']['instance'].update({'frequency': float(self.settings['mw_pulses']['frequency'])})

        while self._abort is False:

            super(ESRSingleFreqCont, self)._function(self.data)

            self.data['esr_counts'].append(self.data['counts'][0])


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
        esr_counts = np.array(data['esr_counts'])

        # if there is two measurement per run, the second serves as a normalization measurement
        if len(np.shape(esr_counts))== 2:
            esr_counts  = esr_counts[:,0]/esr_counts[:,1]


        axis1 = axes_list[0]
        if not esr_counts == []:
            counts = esr_counts
            plot_esr(axis1, mw_frequencies[0:len(counts)], counts)
            axis1.hold(False)
            # axis1.set_title('avrg count')
        axis2 = axes_list[1]
        plot_pulses(axis2, self.pulse_sequences[0])


    def _update_plot(self, axes_list):
        mw_frequencies = self.data['mw_frequencies']
        esr_counts = np.array(self.data['esr_counts'])


        # if there is two measurement per run, the second serves as a normalization measurement
        if len(np.shape(esr_counts)) == 2:
            esr_counts  = esr_counts[:,0]/esr_counts[:,1]

        axis1 = axes_list[0]
        if not esr_counts == []:
            counts = esr_counts
            plot_esr(axis1, mw_frequencies[0:len(counts)], counts)
            axis1.hold(False)
            # axis2 = axes_list[1]
            # update_pulse_plot(axis2, self.pulse_sequences[0])

    def _create_pulse_sequences(self):

        '''

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement(s)
        '''


        # on contrast to other script the following times are given in us and have to be converted to ns
        integration_time = self.settings['read_out']['integration_time']*1e3
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']



        pulse_sequences = [[Pulse('laser', 0, 2*integration_time+2*delay_mw_readout),
                            Pulse('apd_readout', delay_mw_readout, integration_time), #read fluourescence
                            Pulse('microwave_i',2*delay_mw_readout + integration_time, integration_time),
                            Pulse('apd_readout',2*delay_mw_readout + integration_time, integration_time)
                            ]]

        tau_list = [integration_time]
        # end_time_max = 0
        # for pulse_sequence in pulse_sequences:
        #     for pulse in pulse_sequence:
        #         end_time_max = max(end_time_max, pulse.start_time + pulse.duration)
        # for pulse_sequence in pulse_sequences:
        #     pulse_sequence.append(Pulse('laser', end_time_max + 1850, 15))

        return pulse_sequences, self.settings['num_averages'], tau_list, mw_on_time
