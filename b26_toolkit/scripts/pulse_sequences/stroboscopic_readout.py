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
from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, Pulse, MicrowaveGenerator
from b26_toolkit.plotting.plots_1d import plot_counts, update_counts_vs_pos
from b26_toolkit.scripts import Esr


class StroboscopicReadout(PulsedExperimentGeneric):

    _DEFAULT_SETTINGS = [
        Parameter('cycle_period', 10000, int, 'laser and readout off time'),
        Parameter('laser_time', 1000, int, 'time laser is on [ns]'),
        Parameter('meas_time', 1000, int, 'measurement time [ns]'),
        Parameter('delay_readout', 0, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
        Parameter('laser_off_time', 1000, int,
                  'minimum laser off time before taking measurements (ns)'),
        Parameter('red_on_time', 1000, int, 'time that red laser is on'),
        Parameter('pulses_per_seq', 1, int, 'number of laser pulses per sequence'),
        Parameter('num_averages', 100000, float, 'acquisition time for each frequency'),
        Parameter('resonant_readout', [
            Parameter('add_red_laser', False, bool, 'add resonant readout'),
            Parameter('add_pi_pulse', False, bool, 'add pi pulse before readout'),
            Parameter('delay_mw_readout', 210, int, 'delay after mw pulse before readout [ns]'),
            Parameter('mw_power', -45.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', '+i', ['+i', '-i', '+q', '-q'], 'Channel to use for mw pulses'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
        ])

    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}
    _SCRIPTS = {'esr': Esr}

    def _function(self):
        #COMMENT_ME no
        if self.settings['resonant_readout']['add_red_laser'] and self.settings['resonant_readout']['add_pi_pulse']:
            self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
            self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
            self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['resonant_readout']['mw_power']})
            self.instruments['mw_gen']['instance'].update({'frequency': self.settings['resonant_readout']['mw_frequency']})
        super()._function()

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        self._DEFAULT_SETTINGS += PulsedExperimentGeneric._DEFAULT_SETTINGS

        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.ref_index = 0

        # if self.settings['averaging_block_size'] > PulsedExperimentGeneric.MAX_BLOCK_SIZE:
        #     raise Exception('block size too large')

    def _create_pulse_sequences(self):
        # if self.settings['pulses_per_seq'] > self.instruments['PB']['instance'].MAX_COMMANDS / 2:
        #     raise Exception('pulses_per_seq exceeds pulseblaster maximum of {}'.format(self.instruments['PB']['instance'].MAX_COMMANDS))
        pulseSequences = []

        period = self.settings['cycle_period']
        meas_time = self.settings['meas_time']
        laser_time = self.settings['laser_time']
        green_off_time = self.settings['laser_off_time']
        red_on_time = self.settings['red_on_time']
        delay_readout = self.settings['delay_readout']
        microwave_channel = 'microwave_' + self.settings['resonant_readout']['microwave_channel']
        pi_time = self.settings['resonant_readout']['pi_pulse_time']
        delay_mw_readout = self.settings['resonant_readout']['delay_mw_readout']

        try:
            for i in range(self.settings['pulses_per_seq']):
                if not self.settings['resonant_readout']['add_red_laser']:
                    pulseSequence = [
                        Pulse('laser', period * (i + 1) - laser_time, laser_time),
                        Pulse('apd_readout', period * (i + 1) - laser_time + delay_readout, meas_time)
                    ]
                else:
                    if self.settings['resonant_readout']['add_pi_pulse']:
                        pulseSequence = [
                            Pulse('laser', period * (i + 1) - (laser_time + green_off_time + red_on_time + delay_mw_readout + pi_time), laser_time),
                            Pulse(microwave_channel, period * (i + 1) - (red_on_time + delay_mw_readout + pi_time), pi_time),
                            Pulse('red_laser', period * (i + 1) - red_on_time, red_on_time),
                            Pulse('apd_readout', period * (i + 1) - red_on_time + delay_readout, meas_time)
                        ]
                    else:
                        pulseSequence = [
                            Pulse('laser', period * (i + 1) - (laser_time + green_off_time + red_on_time), laser_time),
                            Pulse('red_laser', period * (i + 1) - red_on_time, red_on_time),
                            Pulse('apd_readout', period * (i + 1) - red_on_time + delay_readout, meas_time)
                        ]
                pulseSequences.extend(pulseSequence)
        except AssertionError:
            self.log('Pulse start times might be negative. Try again')

        return [pulseSequences], [period], meas_time

class StroboscopicReadoutRealtime(Script):
    _DEFAULT_SETTINGS = [
        Parameter('num_measurements', -1, int, 'Total number number of runs of strobe readout'),

    ]
    _INSTRUMENTS = {}

    _SCRIPTS = {'strobe_readout': StroboscopicReadout}

    def _function(self):

        num_meas = self.settings['num_measurements']
        readout = self.scripts['strobe_readout']
        self.data = {'counts': []}

        meas_count = 0
        while meas_count != num_meas and not self._abort:
            readout.run()
            counts = np.array(readout.data['counts'])
            if len(counts.shape) == 1:
                self.data['counts'].append(counts[0])

            self.updateProgress.emit(int(self.progress))
            meas_count += 1

    def _plot(self, axes_list, data=None):
        if data is None:
            data = self.data


        settings = self.scripts['strobe_readout'].settings

        counts = data['counts']
        print(counts)
        self._plot_line = axes_list[0].plot(settings['num_averages'] * settings['cycle_period'] * 1e-9 * np.arange(len(counts)), counts, linewidth=1.25)
        # axis.hold(False)

        axes_list[0].set_xlabel('time [s]')
        axes_list[0].set_ylabel('[kCounts/s]')

    def _update_plot(self, axes_list, data=None):
        if data is None:
            data = self.data

        if data and len(data) > 0:
            settings = self.scripts['strobe_readout'].settings
            counts = data['counts']

            self._plot_line[0].set_ydata(counts)
            self._plot_line[0].set_xdata(settings['num_averages'] * settings['cycle_period'] * np.arange(len(counts)) * 1e-9)
            axes_list[0].relim()
            axes_list[0].autoscale_view()

