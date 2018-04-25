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
from b26_toolkit.src.plotting.plots_1d import plot_esr, plot_pulses, update_pulse_plot, plot_1d_simple_timetrace_ns, update_1d_simple
from PyLabControl.src.core import Parameter, Script
from PyLabControl.src.scripts import SelectPoints
from b26_toolkit.src.data_processing.fit_functions import fit_rabi_decay, cose_with_decay, fit_exp_decay, exp_offset
from b26_toolkit.src.scripts import ESR

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

# NOT FINSISHED - TESTING PHASE!!!
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

class Rabi(PulseBlasterBaseScript):
    """
This script applies a microwave pulse at fixed power for varying durations to measure Rabi Oscillations
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses')
        ]),
        Parameter('tau_times', [
            Parameter('max_time', 200, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 5, [5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000, 500000],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 4000, int, 'time with laser on at the beginning to reset state'),
            Parameter('ref_meas_off_time', 1000, int,
                      'laser off time before taking reference measurement at the end of init (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(Rabi, self)._function(self.data)

        counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
        tau = self.data['tau']


        try:
            fits = fit_rabi_decay(tau, counts, variable_phase=True)
            self.data['fits'] = fits
        except:
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
        # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
        #                  self.settings['tau_times']['time_step'])
        # JG 16-08-25 changed (15ns min spacing is taken care of later):
        tau_list = list(range(0, int(self.settings['tau_times']['max_time']),self.settings['tau_times']['time_step']))

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= 15]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        ref_meas_off_time = self.settings['read_out']['ref_meas_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', 0, nv_reset_time - ref_meas_off_time - meas_time - ref_meas_off_time),
                 Pulse('apd_readout', nv_reset_time - ref_meas_off_time - meas_time, meas_time),
                 Pulse('laser', nv_reset_time - ref_meas_off_time - meas_time, meas_time)
                 ]
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence += [Pulse(microwave_channel, nv_reset_time, tau)]

            pulse_sequence += [
                Pulse('laser', nv_reset_time + tau + delay_mw_readout, meas_time),
                Pulse('apd_readout', nv_reset_time + tau + delay_mw_readout, meas_time)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)



            # pulse_sequences.append([Pulse('laser', 0,       nv_reset_time - ref_meas_off_time - meas_time - ref_meas_off_time),
            #                         Pulse('apd_readout',    nv_reset_time - ref_meas_off_time - meas_time, meas_time),
            #                         Pulse('laser',          nv_reset_time - ref_meas_off_time - meas_time, meas_time),
            #                         Pulse(microwave_channel,nv_reset_time,  tau),
            #                         Pulse('laser',          nv_reset_time + tau + delay_mw_readout, meas_time),
            #                         Pulse('apd_readout',    nv_reset_time + tau + delay_mw_readout, meas_time)
            #                         ])


        # end_time_max = 0
        # for pulse_sequence in pulse_sequences:
        #     for pulse in pulse_sequence:
        #         end_time_max = max(end_time_max, pulse.start_time + pulse.duration)
        # for pulse_sequence in pulse_sequences:
        #     pulse_sequence.append(Pulse('laser', end_time_max + 1850, 15))

        return pulse_sequences, self.settings['num_averages'], tau_list, meas_time



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
            fits = data['fits']

            axislist[0].plot(tau, counts, 'b')
            axislist[0].hold(True)

            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            pi_time = 2*np.pi / fits[1] / 2
            axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time))
        else:
            super(Rabi, self)._plot(axislist)
            axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9))
            axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)


class Rabi_double_init(PulseBlasterBaseScript): # ER 5.25.2017
    """
This script applies a microwave pulse at fixed power for varying durations to measure Rabi oscillations.
To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses')
        ]),
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
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(Rabi_double_init, self)._function()

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
        # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
        #                  self.settings['tau_times']['time_step'])
        # JG 16-08-25 changed (15ns min spacing is taken care of later):
        tau_list = list(range(int(self.settings['tau_times']['min_time']),
                              int(self.settings['tau_times']['max_time']),
                              self.settings['tau_times']['time_step']))

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= 15]

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

            pulse_sequence.append(Pulse('laser',
                                        laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time + tau + 2*40 + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time + tau + 2*40 + delay_mw_readout + delay_readout,
                                        meas_time))
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)


        return pulse_sequences, self.settings['num_averages'], tau_list, meas_time



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
            axislist[0].hold(True)

            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            #pi_time = 2*np.pi / fits[1] / 2
            pi_time = (np.pi - fits[2])/fits[1]
            pi_half_time = (np.pi/2 - fits[2])/fits[1]
            three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
            rabi_freq = 1000*fits[1]/(2*np.pi)
         #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
            axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
        else:
            super(Rabi_double_init, self)._plot(axislist)
            axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9))
            axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)

class readout_double_init(PulseBlasterBaseScript):  # ER 10.21.2017
    """
This script sweeps the readout pulse duration. To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulse', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_time', 30.0, float, 'pi time in ns')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 15, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 200, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 5, [5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000, 500000],
                      'time step increment of readout pulse duration (in ns)')
        ]),
        Parameter('read_out', [
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
        # COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulse']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulse']['mw_frequency']})
        super(readout_double_init, self)._function(self.data)


        counts = self.data['counts'][:, 1] # / self.data['counts'][:, 0]
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
        meas_time = 100


        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', laser_off_time, nv_reset_time),
                 Pulse('apd_readout', laser_off_time+ delay_readout, tau),
                 ]
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence += [Pulse(microwave_channel, laser_off_time + nv_reset_time + laser_off_time, pi_time)]

            pulse_sequence += [
                Pulse('laser', laser_off_time + nv_reset_time + laser_off_time + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', laser_off_time + nv_reset_time + laser_off_time + delay_mw_readout + delay_readout, tau)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, self.settings['num_averages'], tau_list, meas_time

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
            counts = data['counts'][:, 1] / data['counts'][:, 0]
            tau = data['tau']
            fits = data['fits']  # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            axislist[0].hold(True)

            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            # pi_time = 2*np.pi / fits[1] / 2
            pi_time = (np.pi - fits[2]) / fits[1]
            pi_half_time = (np.pi / 2 - fits[2]) / fits[1]
            three_pi_half_time = (3 * np.pi / 2 - fits[2]) / fits[1]
            rabi_freq = 1000 * fits[1] / (2 * np.pi)
            #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
            axislist[0].set_title('Readout pulse width counts')
        else:
            super(readout_double_init, self)._plot(axislist)
            axislist[0].set_title('Readout pulse width counts')
            axislist[0].legend(labels=('Ref Fluorescence', 'Pi pulse Data'), fontsize=8)

class T1_double_init(PulseBlasterBaseScript):  # ER 10.21.2017
    """
This script sweeps the readout pulse duration. To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time
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
        super(T1_double_init, self)._function(self.data)


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

        return pulse_sequences, self.settings['num_averages'], tau_list, meas_time

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
            super(T1_double_init, self)._plot(axislist)
            axislist[0].set_title('Readout pulse width counts')
            axislist[0].legend(labels=('Ref Fluorescence', 'Pi pulse Data'), fontsize=8)

class readout_T_double_init(PulseBlasterBaseScript):  # ER 10.21.2017
    """
This script sweeps the readout pulse rise time. To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulse', [
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_time', 30.0, float, 'pi time in ns')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 15, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 200, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 5, [5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000, 500000],
                      'time step increment of readout pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
            Parameter('readout_window', 300, int, 'length of readout window')
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
        super(readout_T_double_init, self)._function(self.data)


        counts = self.data['counts'][:, 1] # / self.data['counts'][:, 0]
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
        meas_time = self.settings['read_out']['readout_window']


        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', laser_off_time, nv_reset_time),
                 Pulse('apd_readout', laser_off_time+ delay_readout + tau, meas_time),
                 ]
            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                pulse_sequence += [Pulse(microwave_channel, laser_off_time + nv_reset_time + laser_off_time, pi_time)]

            pulse_sequence += [
                Pulse('laser', laser_off_time + nv_reset_time + laser_off_time + delay_mw_readout, nv_reset_time),
                Pulse('apd_readout', laser_off_time + nv_reset_time + laser_off_time + delay_mw_readout + delay_readout + tau, meas_time)
            ]
            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, self.settings['num_averages'], tau_list, meas_time

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
            counts = data['counts'][:, 1] / data['counts'][:, 0]
            tau = data['tau']
            fits = data['fits']  # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            axislist[0].hold(True)

            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            # pi_time = 2*np.pi / fits[1] / 2
            pi_time = (np.pi - fits[2]) / fits[1]
            pi_half_time = (np.pi / 2 - fits[2]) / fits[1]
            three_pi_half_time = (3 * np.pi / 2 - fits[2]) / fits[1]
            rabi_freq = 1000 * fits[1] / (2 * np.pi)
            #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
            axislist[0].set_title('Readout pulse width counts')
        else:
            super(readout_T_double_init, self)._plot(axislist)
            axislist[0].set_title('Readout pulse width counts')
            axislist[0].legend(labels=('Ref Fluorescence', 'Pi pulse Data'), fontsize=8)

class HahnEcho_double_init(PulseBlasterBaseScript): # ER 5.25.2017
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
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(HahnEcho_double_init, self)._function(self.data)

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

        return pulse_sequences, self.settings['num_averages'], tau_list, meas_time



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
            axislist[0].hold(True)

            axislist[0].plot(tau, exp_offset(tau, fits[0], fits[1], fits[2]))
            axislist[0].set_title('T2 decay time (simple exponential, p = 1): {:2.1f} ns'.format(fits[1]))
        else:
            super(HahnEcho_double_init, self)._plot(axislist)
            axislist[0].set_title('Hahn Echo mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9))
            axislist[0].legend(labels=('Ref Fluorescence', 'T2 Data'), fontsize=8)

class CPMG_double_init(PulseBlasterBaseScript): # ER 5.25.2017
    """
This script runs a Hahn echo on the NV to find the Hahn echo T2. To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time
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
            Parameter('Number of pi pulses N', 4, int, 'number of pi pulses in the CPMG-N sequence')
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
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(CPMG_double_init, self)._function(self.data)

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
                Pulse(microwave_channel_pi2, laser_off_time, pi_half_time),
                Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau/2 - pi_time/2., pi_time)
            ]

            next_pi_t =  laser_off_time + pi_half_time/2. + tau/2 - pi_time/2. + tau
            N = self.settings['mw_pulses']['Number of pi pulses N']

            for ind in range(0, N-1):
                pulse_sequence += [
                    Pulse(microwave_channel, next_pi_t, pi_time)
                ]
                next_pi_t = next_pi_t + tau

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
                Pulse(microwave_channel, start_of_second_CPMG + pi_half_time/2. + tau/2 - pi_time/2., pi_time)
            ]
            next_pi_t =  start_of_second_CPMG + pi_half_time/2. + tau/2 - pi_time/2. + tau
            for ind in range(0, N-1):
                pulse_sequence += [
                    Pulse(microwave_channel, next_pi_t, pi_time)
                    ]
                next_pi_t = next_pi_t + tau

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

        return pulse_sequences, self.settings['num_averages'], tau_list, meas_time



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
            axislist[0].hold(True)

            axislist[0].plot(tau, exp_offset(tau, fits[0], fits[1], fits[2]))
            axislist[0].set_title('T2 decay time (simple exponential, p = 1): {:2.1f} ns'.format(fits[1]))
        else:
            super(CPMG_double_init, self)._plot(axislist)
            axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9))
            axislist[0].legend(labels=('Ref Fluorescence', 'T2 Data'), fontsize=8)

class XY8_double_init(PulseBlasterBaseScript): # ER 5.25.2017
    """
This script runs a Hahn echo on the NV to find the Hahn echo T2. To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time
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

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(XY8_double_init, self)._function(self.data)

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
                if counter == 0 or counter == 2 or counter == 5 or counter == 7:
                    pulse_sequence += [
                        Pulse(microwave_channel_pi2, next_pi_t, pi_time) # pulses along x
                    ]
                else:
                    pulse_sequence += [
                        Pulse(microwave_channel, next_pi_t, pi_time) # pulses along y
                    ]
                if counter == 7:
                    counter = -1
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
                if counter == 0 or counter == 2 or counter == 5 or counter == 7:
                    pulse_sequence += [
                        Pulse(microwave_channel_pi2, next_pi_t, pi_time) # pulses along x
                    ]
                else:
                    pulse_sequence += [
                        Pulse(microwave_channel, next_pi_t, pi_time) # pulses along y
                    ]
                if counter == 7:
                    counter = -1
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


        return pulse_sequences, self.settings['num_averages'], tau_list, meas_time



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
            axislist[0].hold(True)

            axislist[0].plot(tau, exp_offset(tau, fits[0], fits[1], fits[2]))
            axislist[0].set_title('T2 decay time (simple exponential, p = 1): {:2.1f} ns'.format(fits[1]))
        else:
            super(XY8_double_init, self)._plot(axislist)
            axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9))
            axislist[0].legend(labels=('Ref Fluorescence', 'T2 Data'), fontsize=8)

# pulse sequence is X Y X Y X Y X Y .... to accumulate pulse errors and calibrate phase
class XYXY_double_init(PulseBlasterBaseScript): # ER 5.25.2017
    """
This script runs a Hahn echo on the NV to find the Hahn echo T2. To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time
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

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(XYXY_double_init, self)._function(self.data)

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

        return pulse_sequences, self.settings['num_averages'], tau_list, meas_time



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
            axislist[0].hold(True)

            axislist[0].plot(tau, exp_offset(tau, fits[0], fits[1], fits[2]))
            axislist[0].set_title('T2 decay time (simple exponential, p = 1): {:2.1f} ns'.format(fits[1]))
        else:
            super(XYXY_double_init, self)._plot(axislist)
            axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9))
            axislist[0].legend(labels=('Ref Fluorescence', 'T2 Data'), fontsize=8)

class Rabi_Power_Sweep_Single_Tau(PulseBlasterBaseScript):
    """
This script applies a microwave pulse at fixed power for varying durations to measure Rabi Oscillations
    """
    _DEFAULT_SETTINGS = [
        Parameter('min_mw_power', -45.0, float, 'minimum microwave power in dB'),
        Parameter('max_mw_power', -45.0, float, 'maximum microwave power in dB'),
        Parameter('mw_power_step', 1.0, float, 'power to step by in dB'),
        Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
        Parameter('mw_time', 200, float, 'total time of rabi oscillations (in ns)'),
        Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('reset_time', 10000, int, 'time with laser on at the beginning to reset state'),
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_frequency']})
        mw_power_values = np.arange(self.settings['min_mw_power'],
                                    self.settings['max_mw_power'] + self.settings['mw_power_step'],
                                    self.settings['mw_power_step'])

        print(mw_power_values)
        self.data = {'mw_power_values': mw_power_values, 'counts_for_mw': np.zeros(len(mw_power_values))}
        for index, power in enumerate(mw_power_values):
            self.instruments['mw_gen']['instance'].update({'amplitude': float(power)})
            super(Rabi_Power_Sweep_Single_Tau, self)._function(self.data)
            self.data['counts_for_mw'][index] = self.data['counts'][0]

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
        reset_time = self.settings['reset_time']
        mw_time = self.settings['mw_time']
        pulse_sequences.append([Pulse('laser', 0, reset_time),
                                Pulse('microwave_i', reset_time + 200, mw_time),
                                Pulse('laser', reset_time + mw_time + 300, self.settings['meas_time']),
                                Pulse('apd_readout', reset_time + mw_time + 300, self.settings['meas_time'])
                                ])

        end_time_max = 0
        for pulse_sequence in pulse_sequences:
            for pulse in pulse_sequence:
                end_time_max = max(end_time_max, pulse.start_time + pulse.duration)
        for pulse_sequence in pulse_sequences:
            pulse_sequence.append(Pulse('laser', end_time_max + 1850, 15)) # Jan Feb 1st 2017: what is 1850??? Need to comment!

        return pulse_sequences, self.settings['num_averages'], [mw_time], self.settings['meas_time']

    def _plot(self, axes_list, data = None):
        '''
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional) dataset to plot (dictionary that contains keys counts_for_mw, mw_power_values), if not provided use self.data
        '''
        if data is None:
            data = self.data

        counts = data['counts_for_mw']
        x_data = data['mw_power_values']
        axis1 = axes_list[0]
        if not counts == []:
            plot_1d_simple_timetrace_ns(axis1, x_data, [counts], x_label='microwave power (dBm)')
        axis2 = axes_list[1]
        plot_pulses(axis2, self.pulse_sequences[self.sequence_index])

    def _update_plot(self, axes_list):
        '''
        Updates plots specified in _plot above
        Args:
            axes_list: list of axes to write plots to (uses first 2)

        '''
        counts = self.data['counts_for_mw']
        x_data = self.data['mw_power_values']
        axis1 = axes_list[0]
        if not counts == []:
            update_1d_simple(axis1, x_data, [counts])
        axis2 = axes_list[1]
        update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])

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

class CalibrateMeasurementWindow(PulseBlasterBaseScript):
    """
This script find the optimal duration of the measurment window.
It applies a sliding measurement window with respect to a readout from the NV 0 state and the NV 1 state.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
        Parameter('pi_pulse_time', 50, float, 'time duration of pi-pulse (in ns)'),
        Parameter('readout_window_incremement', 10, [5, 10, 20, 50, 100], 'time step increment of measurement duration (in ns)'),
        Parameter('initial_readout_displacement', -80, int, 'min time of measurement duration (in ns)'),
        Parameter('final_readout_displacement', 450, int, 'max time of measurement duration (in ns)'),
        Parameter('reset_time', 3000, int, 'time with laser on at the beginning to reset state'),
        Parameter('delay_init_mw', 200, int, 'time delay before pi pulse after NV reset'),
        Parameter('delay_mw_readout', 200, int, 'time delay before readout after pi pulse'),
        Parameter('measurement_window_width', 20, int, 'the width of the sliding readout window'),
        Parameter('laser_on_time', 500, list(range(100, 1201, 100)), 'time laser is on for readout'),
        Parameter('ref_meas_off_time', 1000, int, 'time reset laser is turned off before reference measurement is made'),
        Parameter('num_averages', 1000000, int, 'number of averages')
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_frequency']})
        super(CalibrateMeasurementWindow, self)._function()

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
        pulse_sequences = []
        tau_list = list(range(self.settings['initial_readout_displacement'],
                         self.settings['final_readout_displacement'],
                         self.settings['readout_window_incremement']))
        reset_time = self.settings['reset_time']

        for tau in tau_list:
            pulse_sequences.append([Pulse('laser', 0, reset_time - self.settings['ref_meas_off_time'] - self.settings['laser_on_time']),
                                    Pulse('apd_readout', reset_time - self.settings['laser_on_time'] + tau, self.settings['measurement_window_width']),
                                    Pulse('laser', reset_time - self.settings['laser_on_time'], self.settings['laser_on_time']),
                                    Pulse('microwave_i', reset_time + self.settings['delay_init_mw'], self.settings['pi_pulse_time']),
                                    Pulse('laser', reset_time + self.settings['delay_init_mw'] + self.settings['pi_pulse_time'] + self.settings[
                                        'delay_mw_readout'], self.settings['laser_on_time']),
                                    Pulse('apd_readout', reset_time + self.settings['delay_init_mw'] + self.settings['pi_pulse_time'] +
                                          self.settings['delay_mw_readout'] + tau, self.settings['measurement_window_width'])
                                    ])

        return pulse_sequences, self.settings['num_averages'], tau_list, self.settings['measurement_window_width']

    def _plot(self, axes_list, data = None):
        """
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first)
            data (optional): dataset to plot (dictionary that contains keys counts, tau), if not provided use self.data
        """

        super(CalibrateMeasurementWindow, self)._plot(axes_list, data)
        axes_list[0].set_title('Measurement Calibration')
        axes_list[0].legend(labels=('|0> State Fluorescence', '|1> State Fluoresence'), fontsize=8)

class XY8(PulseBlasterBaseScript):
    """
This script runs an XY pulse sequence.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses',[
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            # Parameter('mw_switch_extra_time', 15, int, 'Time to add before and after microwave switch is turned on'),
            Parameter('pi_pulse_time', 50, float, 'time duration of pi-pulse (in ns)'),
            Parameter('number_of_pulse_blocks', 1, list(range(1, 17)), 'number of alternating x-y-x-y-y-x-y-x pulses'),
            Parameter('end_in_0', False, bool, 'end with 3pi/2 pulse so end state is |0> rather than |1>')
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

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}
    _SCRIPTS = {}

    def _function(self):
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(XY8, self)._function()

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
                                   Pulse('microwave_i', reset_time + delay_mw_init-pi_half_time/2, pi_half_time)
                                   ])

            #CPMG xyxyyxyx loops added number_of_pulse_blocks times
            section_begin_time = reset_time + delay_mw_init - tau/2 #for the first pulse, only wait tau/2
            # JG 16-08-19 - begin changed to pi time instead of pi/2
            # section_begin_time = reset_time + delay_mw_init + pi_time
            # JG 16-08-19 - end

            # for i in range(0, number_of_pulse_blocks):
            #     pulse_sequence.extend([Pulse('microwave_i', section_begin_time + 1*tau - pi_half_time, pi_time),
            #                            Pulse('microwave_q', section_begin_time + 2*tau - pi_half_time, pi_time),
            #                            Pulse('microwave_i', section_begin_time + 3*tau - pi_half_time, pi_time),
            #                            Pulse('microwave_q', section_begin_time + 4*tau - pi_half_time, pi_time),
            #                            Pulse('microwave_q', section_begin_time + 5*tau - pi_half_time, pi_time),
            #                            Pulse('microwave_i', section_begin_time + 6*tau - pi_half_time, pi_time),
            #                            Pulse('microwave_q', section_begin_time + 7*tau - pi_half_time, pi_time),
            #                            Pulse('microwave_i', section_begin_time + 8*tau - pi_half_time, pi_time)
            #                           ])
            #     section_begin_time += 8*tau

            # AK 17-02-28 - switched to yx rather than xy since we saw echo was better with rephasing pulses
            #               perpendicular to pi/2 pulses
            for i in range(0, number_of_pulse_blocks):
                pulse_sequence.extend([Pulse('microwave_q', section_begin_time + 1*tau - pi_half_time, pi_time),
                                       Pulse('microwave_i', section_begin_time + 2*tau - pi_half_time, pi_time),
                                       Pulse('microwave_q', section_begin_time + 3*tau - pi_half_time, pi_time),
                                       Pulse('microwave_i', section_begin_time + 4*tau - pi_half_time, pi_time),
                                       Pulse('microwave_i', section_begin_time + 5*tau - pi_half_time, pi_time),
                                       Pulse('microwave_q', section_begin_time + 6*tau - pi_half_time, pi_time),
                                       Pulse('microwave_i', section_begin_time + 7*tau - pi_half_time, pi_time),
                                       Pulse('microwave_q', section_begin_time + 8*tau - pi_half_time, pi_time)
                                      ])
                section_begin_time += 8*tau


            if self.settings['mw_pulses']['end_in_0']:
                # 3pi/2 and readout
                pulse_sequence.extend([Pulse('microwave_i', section_begin_time + tau / 2 - 3*pi_half_time/4, 3*pi_half_time),
                                       Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                             meas_time),
                                       Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                             meas_time)])
            else:
                #pi/2 and readout
                pulse_sequence.extend([Pulse('microwave_i', section_begin_time + tau/2 - pi_half_time/2, pi_half_time),
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

        return pulse_sequences, self.settings['num_averages'], tau_list, meas_time


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

        super(XY8, self)._plot(axislist, data)
        axislist[0].set_title('XY8')
        axislist[0].legend(labels=('Ref Fluorescence', 'XY8 data'), fontsize=8)

class XY4(PulseBlasterBaseScript):
    """
This script runs a CPMG pulse sequence.
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

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}
    _SCRIPTS = {}

    def _function(self):
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
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

        return pulse_sequences, self.settings['num_averages'], tau_list, meas_time


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

class PDD(PulseBlasterBaseScript):
    """
This script runs a PDD ( Periodic Dynamical Decoupling) sequence for different number of pi pulses.
For a single pi-pulse this is a Hahn-echo sequence.
For zero pulses this is a Ramsey sequence.

The sequence is pi/2 - tau/4 - (tau/4 - pi  - tau/4)^n - tau/4 - pi/2

Tau/2 is the time between the center of the pulses!


    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -2, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('pi_pulse_time', 50, float, 'time duration of pi-pulse (in ns)'),
            Parameter('number_of_pi_pulses', 1, list(range(0, 17)), 'number of pi pulses')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 15, float, 'min value for tau, the free evolution time in between pulses (in ns)'),
            Parameter('max_time', 30, float, 'max value for tau, the free evolution time in between pulses (in ns)'),
            Parameter('time_step', 5, [5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000], 'step size for tau, the free evolution time in between pulses (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 250, float, 'measurement time after CPMG sequence (in ns)'),
            Parameter('nv_reset_time', 3000, int, 'time duration of the green laser to reset the spin state'),
            Parameter('ref_meas_off_time', 1000, int, 'laser off time before taking reference measurement at the end of init (ns)'),
            Parameter('delay_mw_init', 1000, int, 'delay between initialization and mw (in ns)'),
            Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)')
        ]),
        Parameter('num_averages', 1000, int, 'number of averages (should be less than a million)')
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    _SCRIPTS = {}


    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

    def _function(self):
        #COMMENT_ME
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        super(PDD, self)._function()


    def _create_pulse_sequences(self):
        '''
        creates the pulse sequence for the Hahn echo /
        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''
        pulse_sequences = []


        tau_list = list(range(int(self.settings['tau_times']['min_time']),
                         int(self.settings['tau_times']['max_time'] + self.settings['tau_times']['time_step']),
                         self.settings['tau_times']['time_step']))



        reset_time = self.settings['read_out']['nv_reset_time']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = pi_time/2.0

        ref_meas_off_time = self.settings['read_out']['ref_meas_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_init = self.settings['read_out']['delay_mw_init']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        number_of_pi_pulses = self.settings['mw_pulses']['number_of_pi_pulses']


        for tau in tau_list:

            # pulse_sequence = [Pulse('laser', 0, reset_time - ref_meas_off_time - 15 - meas_time),
            #                   Pulse('apd_readout', reset_time - 15 - meas_time, meas_time),
            #                   Pulse('laser', reset_time - 15 - meas_time, meas_time),
            #                   Pulse('microwave_i', reset_time + delay_mw_init, pi_half_time)
            #                   ]
            # # 16-08-25 JG: changed :
            pulse_sequence = [Pulse('laser', 0, reset_time - ref_meas_off_time - 15 - meas_time),
                              Pulse('apd_readout', reset_time - 15 - meas_time, meas_time),
                              Pulse('laser', reset_time - 15 - meas_time, meas_time),
                              Pulse('microwave_i', reset_time + delay_mw_init-pi_half_time/2, pi_half_time)
                              ]


            # next_pi_pulse_time = reset_time + delay_mw_init + pi_half_time + tau
            # # 16-08-19 JG: changed :
            next_pi_pulse_time = reset_time + delay_mw_init
            # # 16-08-25 JG: changed :
            # next_pi_pulse_time = reset_time + delay_mw_init - pi_half_time / 2 + tau / 2

            for n in range(1, number_of_pi_pulses + 1):
                next_pi_pulse_time += tau/2
                pulse_sequence.extend([Pulse('microwave_q', next_pi_pulse_time - pi_time/2, pi_time)])
                # next_pi_pulse_time += tau*2 + pi_time
                # 16-08-19 JG: changed:
                # next_pi_pulse_time += tau
                # 16 - 08 -24 JG: changed
                next_pi_pulse_time += tau/2

            if number_of_pi_pulses == 0:
                next_pi_pulse_time += tau

            # pulse_sequence.extend([Pulse('microwave_i', next_pi_pulse_time-tau, pi_half_time),
            #                         Pulse('laser', next_pi_pulse_time-tau + delay_mw_readout + pi_half_time, meas_time),
            #                         Pulse('apd_readout',next_pi_pulse_time-tau + delay_mw_readout + pi_half_time, meas_time)
            #                         ])
            # 16-08-19 JG: changed:
            # pulse_sequence.extend([Pulse('microwave_i', next_pi_pulse_time-tau/2 + pi_half_time, pi_half_time),
            #                         Pulse('laser',      next_pi_pulse_time-tau/2 + pi_time + delay_mw_readout, meas_time),
            #                         Pulse('apd_readout',next_pi_pulse_time-tau/2 + pi_time + delay_mw_readout, meas_time)
            #                         ])
            # pulse_sequences.append(pulse_sequence)
            # 16 - 08 -24 JG: changed
            # pulse_sequence.extend([Pulse('microwave_i', next_pi_pulse_time + pi_half_time, pi_half_time),
            #                        Pulse('laser', next_pi_pulse_time + pi_time + delay_mw_readout, meas_time),
            #                        Pulse('apd_readout', next_pi_pulse_time + pi_time + delay_mw_readout,
            #                              meas_time)
            #                        ])

            # # 16-08-25 JG: changed :
            pulse_sequence.extend([Pulse('microwave_i', next_pi_pulse_time - pi_half_time/2, pi_half_time),
                                   Pulse('laser', next_pi_pulse_time + pi_half_time + delay_mw_readout, meas_time),
                                   Pulse('apd_readout', next_pi_pulse_time + pi_half_time + delay_mw_readout, meas_time)
                                   ])

            pulse_sequences.append(pulse_sequence)

        # TEMPORATTY: THIS IS TO SEE IF THE OVERALL TIME OF A SEQUENCE SHOULD ALWAYS BE THE SAME
        # IF WE WANT TO KEEP THIS ADD ADDITIONAL PARAMETER TO THE SCRIPT SETTINGS
        # end_time_max = 0
        # for pulse_sequence in pulse_sequences:
        #     for pulse in pulse_sequence:
        #         end_time_max = max(end_time_max, pulse.start_time + pulse.duration)
        # for pulse_sequence in pulse_sequences:
        #     pulse_sequence.append(Pulse('laser', end_time_max + 1850, 15))

        return pulse_sequences, self.settings['num_averages'], tau_list, meas_time

class XY(PulseBlasterBaseScript):
    """
This script runs a XY sequence for different number of pi pulses. Without pi-pulse this is a Ramsey sequence.
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
        Parameter('pi_half_pulse_time', 50, float, 'time duration of pi-pulse (in ns)'),
        Parameter('number_of__pi_pulses', 0, list(range(0,17)), 'number of pi pulses'),
        Parameter('tau', [
            Parameter('min', 15, float, 'min value for tau, the free evolution time in between pulses (in ns)'),
            Parameter('max', 30, float, 'max value for tau, the free evolution time in between pulses (in ns)'),
            Parameter('step', 5, float, 'step size for tau, the free evolution time in between pulses (in ns)'),
        ]),
        Parameter('meas_time', 300, float, 'measurement time after CPMG sequence (in ns)'),
        Parameter('num_averages', 1000, int, 'number of averages (should be less than a million)'),
        Parameter('reset_time', 1000, int, 'time duration of the green laser to reset the spin state'),
        Parameter('delay_init_mw', 100, int, 'delay between initialization and mw (in ns)'),
        Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
        Parameter('ref_meas_off_time', 1000, int,'laser off time before taking reference measurement at the end of init (ns)')
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    _SCRIPTS = {}


    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

    def _function(self):
        #COMMENT_ME
        self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_frequency']})
        super(XY, self)._function()


    def _create_pulse_sequences(self):
        '''
        creates the pulse sequence for the Hahn echo /
        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''
        pulse_sequences = []

        tau_list = list(range(int(max(15,self.settings['tau']['min'])), int(self.settings['tau']['max']),int(self.settings['tau']['step'])))
        reset_time = self.settings['reset_time']
        mw_delay_time = self.settings['delay_init_mw']
        delay_after_mw = self.settings['delay_mw_readout']
        pi_half_pulse_time = self.settings['pi_half_pulse_time']
        meas_time  = self.settings['meas_time']
        number_of__pi_pulses =  self.settings['number_of__pi_pulses']

        for tau in tau_list:
            # if number_of__pi_pulses == 0:
            #     pulse_sequences.append([Pulse('laser', 0, reset_time),
            #                             Pulse('microwave_i', reset_time+ mw_delay_time, pi_half_pulse_time),
            #                             Pulse('microwave_i', reset_time + mw_delay_time+ pi_half_pulse_time + tau, pi_half_pulse_time),
            #                             Pulse('laser', reset_time + mw_delay_time+ pi_half_pulse_time + tau + pi_half_pulse_time, meas_time),
            #                             Pulse('apd_readout', reset_time + mw_delay_time+ pi_half_pulse_time + tau + pi_half_pulse_time, meas_time)
            #                             ])
            # else:

            pulse_sequence = []

            pulse_sequence.extend([Pulse('laser', 0, reset_time - self.settings['ref_meas_off_time'] - 15 - self.settings['meas_time']),
                                    Pulse('apd_readout', reset_time - 15 - self.settings['meas_time'], self.settings['meas_time']),
                                    Pulse('laser', reset_time - 15 - self.settings['meas_time'], self.settings['meas_time']),
                                    Pulse('microwave_i', reset_time + mw_delay_time, pi_half_pulse_time)
                                    ])

            next_pi_pulse_time = reset_time + mw_delay_time + pi_half_pulse_time + tau

            for n in range(1, number_of__pi_pulses + 1):
                pulse_sequence.extend([Pulse('microwave_i', next_pi_pulse_time, 2*pi_half_pulse_time)])
                pulse_sequence.extend([Pulse('microwave_q', next_pi_pulse_time + 2*pi_half_pulse_time + 2*tau, 2*pi_half_pulse_time)])
                next_pi_pulse_time += tau*4 + 4*pi_half_pulse_time

            pulse_sequence.extend([Pulse('microwave_i', next_pi_pulse_time-tau,pi_half_pulse_time),
                                    Pulse('laser', next_pi_pulse_time-tau + delay_after_mw + pi_half_pulse_time, meas_time),
                                    Pulse('apd_readout',next_pi_pulse_time-tau + delay_after_mw + pi_half_pulse_time, meas_time)
                                    ])

            pulse_sequences.append(pulse_sequence)


        # TEMPORATTY: THIS IS TO SEE IF THE OVERALL TIME OF A SEQUENCE SHOULD ALWAYS BE THE SAME
        # IF WE WANT TO KEEP THIS ADD ADDITIONAL PARAMETER TO THE SCRIPT SETTINGS
        # end_time_max = 0
        # for pulse_sequence in pulse_sequences:
        #     for pulse in pulse_sequence:
        #         end_time_max = max(end_time_max, pulse.start_time + pulse.duration)
        # for pulse_sequence in pulse_sequences:
        #     pulse_sequence.append(Pulse('laser', end_time_max + 1850, 15))

        return pulse_sequences, self.settings['num_averages'], tau_list, self.settings['meas_time']

class T1(PulseBlasterBaseScript):
    """
This script measures the relaxation time of an NV center
    """
    _DEFAULT_SETTINGS = [
        Parameter('time_step', 1000, int, 'time step increment of T1 measurement (ns)'),
        Parameter('max_time', 200000, float, 'total time of T1 measurement (ns)'),
        Parameter('meas_time', 300, float, 'measurement time of fluorescence counts (ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter('nv_reset_time', 3000, int, 'time with laser on at the beginning to reset state (ns)'),
        Parameter('ref_meas_off_time', 1000, int,'laser off time before taking reference measurement at the end of init (ns)'),
        Parameter('tau_scale', 'linear', ['linear', 'logarithmic'])
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster}

    def _create_pulse_sequences(self):
        """
        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a  pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """

        pulse_sequences = []
        if self.settings['time_step'] % 5 != 0:
            raise AttributeError('given time_step is not a multiple of 5')

        tau_list = list(range(0, int(self.settings['max_time'] + self.settings['time_step']), self.settings['time_step']))
        reset_time = self.settings['nv_reset_time']

        # reduce the initialization time by 15 ns to avoid touching DAQ pulses
        # (they are problematic because the DAQ expects two pulse but get only one because they get merged by the pulse blaster)
        for tau in tau_list:
            pulse_sequences.append(
                [Pulse('laser', 0, reset_time - self.settings['ref_meas_off_time'] - 15 - self.settings['meas_time']),
                 Pulse('apd_readout', reset_time - 15 - self.settings['meas_time'], self.settings['meas_time']),
                 Pulse('laser', reset_time - 15 - self.settings['meas_time'], self.settings['meas_time']),
                 Pulse('apd_readout', reset_time + tau, self.settings['meas_time']),
                 Pulse('laser', reset_time + tau, self.settings['meas_time']),
                 ])

        return pulse_sequences, self.settings['num_averages'], tau_list, self.settings['meas_time']

    def _plot(self, axislist, data = None):
        """
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first)
            data (optional): dataset to plot (dictionary that contains keys counts, tau), if not provided use self.data
        """
        super(T1, self)._plot(axislist, data)
        axislist[0].set_title('T1')
        axislist[0].legend(labels=( 'Ref Fluorescence', 'T1 data'), fontsize=8)

class T1SpinFlip(PulseBlasterBaseScript):
    """
This script measures the relaxation time of an NV center.
Optionally a microwave pulse is applied as part of the initialization to prepare the system in a different state
    """
    _DEFAULT_SETTINGS = [
        Parameter(
            'tau_times',
            [
                Parameter('time_step', 1000, int, 'time step increment of T1 measurement (ns)'),
                Parameter('max_time', 200000, float, 'max time of T1 measurement (ns)'),
                Parameter('min_time', 0, float, 'min time of T1 measurement (ns)')
            ]
        ),
        Parameter('num_averages', 1000000, int, 'number of averages'),
        Parameter(
            'read_out',
            [
                Parameter('meas_time', 700, float, 'measurement time of fluorescence counts (ns)'),
                Parameter('nv_reset_time', 3000, int, 'time with laser on at the beginning to reset state (ns)'),
                Parameter('ref_meas_off_time', 1000, int, 'laser off time before taking reference measurement at the end of init (ns)')
            ]
        ),
        Parameter('apply mw-pulse', True, bool, 'if true a pi pulse is at the beginning of the measurement'),
        Parameter('mw-pulse',
                  [
                      Parameter('mw_frequency', 2.87e9, float, 'microwave frequency of pi pulse (Hz)'),
                      Parameter('mw_duration', 300, int, 'pi pulse duration (ns)'),
                      Parameter('mw_channel', 'i', ['i', 'q'], 'select i or q channel for i pulse'),
                      Parameter('mw_power', -2, float, 'microwave power (dB)')
                  ]
                  )
    ]

    _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}

    def _function(self):
        #COMMENT_ME
        self.instruments['mw_gen']['instance'].update({
            'modulation_type': 'IQ',
            'amplitude': self.settings['mw-pulse']['mw_power'],
            'frequency': self.settings['mw-pulse']['mw_frequency']
        })
        super(T1SpinFlip, self)._function()

    def _create_pulse_sequences(self):
        """
        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a  pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """

        pulse_sequences = []
        if self.settings['tau_times']['time_step'] % 5 != 0:
            raise AttributeError('given time_step is not a multiple of 5')

        tau_list = list(range(int(self.settings['tau_times']['min_time']),
                         int(self.settings['tau_times']['max_time'] + self.settings['tau_times']['time_step']),
                         self.settings['tau_times']['time_step']))


        # if self.settings['apply mw-pulse']:
        #     ref_meas_off_time = self.settings['read_out']['ref_meas_off_time'] + self.settings['mw-pulse']['mw_duration']
        # else:
        ref_meas_off_time = self.settings['read_out']['ref_meas_off_time']
        reset_time = self.settings['read_out']['nv_reset_time']
        meas_time = self.settings['read_out']['meas_time']

        microwave_channel = 'microwave_' + self.settings['mw-pulse']['mw_channel']
        microwave_duration = self.settings['mw-pulse']['mw_duration']

        # reduce the initialization time by 15 ns to avoid touching DAQ pulses
        # (they are problematic because the DAQ expects two pulse but get only one because they get merged by the pulse blaster)
        def build_sequence(tau):
            """
            builds the sequence for a given tau
            Args:
                tau: the time after the initialization at which to measure the population

            Returns: the sequence for tau

            """

            if self.settings['apply mw-pulse']:
                sequence = [
                    Pulse('laser', 0,       reset_time - ref_meas_off_time - meas_time - 15 - ref_meas_off_time- microwave_duration),
                    Pulse('apd_readout',    reset_time - ref_meas_off_time - meas_time - 15 - microwave_duration, meas_time),
                    Pulse('laser',          reset_time - ref_meas_off_time - meas_time - 15 - microwave_duration, meas_time),
                    Pulse(microwave_channel,reset_time - 15 - ref_meas_off_time/2.- microwave_duration, microwave_duration)
                ]
            else:
                sequence = [
                    Pulse('laser', 0, reset_time - ref_meas_off_time - meas_time - 15),
                    Pulse('apd_readout', reset_time - 15 - meas_time, meas_time),
                    Pulse('laser', reset_time - 15 - meas_time, meas_time)
                ]

            sequence += [
                Pulse('apd_readout', reset_time + tau, meas_time),
                Pulse('laser',       reset_time + tau, meas_time)
            ]

            return sequence

        pulse_sequences = [build_sequence(tau) for tau in tau_list]

        return pulse_sequences, self.settings['num_averages'], tau_list, self.settings['read_out']['meas_time']

    def _plot(self, axislist, data = None):
        """
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first)
            data (optional): dataset to plot (dictionary that contains keys counts, tau), if not provided use self.data
        """
        super(T1SpinFlip, self)._plot(axislist, data)
        axislist[0].set_title('T1')
        axislist[0].legend(labels=( 'Ref Fluorescence', 'T1 data'), fontsize=8)

class T1_single_init(PulseBlasterBaseScript): # ER 5.25.2017
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

        super(T1_single_init, self)._function(self.data)
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

        return pulse_sequences, self.settings['num_averages'], tau_list, meas_time



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
            super(T1_single_init, self)._plot(axislist)
            axislist[0].set_title('T1 decay of ms = 0 population')

class T1_double_init_many_NVs(Script):
    _DEFAULT_SETTINGS = [
        Parameter('esr_peak', 'upper', ['upper', 'lower', 'both'], 'if ESR fits two peaks, defines which one to use')
    ]
    _INSTRUMENTS = {}
    _SCRIPTS = {'select_NVs': SelectPoints, 'ESR': ESR, 'Rabi': Rabi_double_init, 'T1': T1_double_init}

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

class HE_double_init_many_NVs(Script):
    _DEFAULT_SETTINGS = [
        Parameter('esr_peak', 'upper', ['upper', 'lower', 'both'], 'if ESR fits two peaks, defines which one to use')
    ]
    _INSTRUMENTS = {}
    _SCRIPTS = {'select_NVs': SelectPoints, 'ESR': ESR, 'Rabi': Rabi_double_init, 'HE': HahnEcho_double_init}

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



if __name__ == '__main__':
    # ===================================== 1  ================================================================

    updated_scripts, load_failed, updated_instruments = Script.load_and_append({'XY8':'XY8_double_init'}, package='b26_toolkit')
    xy8 = updated_scripts['XY8']


    xy8.update({'Tracking':{'on/off':False}}) # turn off tracking because this will cause an error if we don't run findnv
    print(xy8)

    xy8.is_valid()















    # # ===================================== 1  ================================================================
    # import os
    # from PyLabControl.src.core.read_write_functions import load_b26_file
    # filename = os.path.normpath('Z:/Lab/Cantilever/Measurements/20180319_Sample_37_diamond_Y_NA0.9/180321-19_41_00_xy8_6_double_init_7Vpp_2MHz/180321-19_41_00_xy8_6_double_init_7Vpp_2MHz.b26')
    # print(filename)
    # script = load_b26_file(filename)['scripts']
    # print(script.keys())
    #
    #
    #
    #
    # updated_scripts, load_failed, updated_instruments = Script.load_and_append(script)
    # xy8 = updated_scripts['XY8_double_init']
    # print(xy8)
    #
    # pulse_sequences, num_averages, tau_list, measurement_gate_width, failure_list = xy8.validate()
    #
    # print('ddddd', len(failure_list))
    #
    # print('ddddd', failure_list)
    # script = {}
    # instr = {}
    # script, failed, instr = Script.load_and_append({'T1_double_init_many_NVs': 'T1_double_init_many_NVs'}, script, instr)
    # print(script)
    # print('failed', failed)
    # print(instr)