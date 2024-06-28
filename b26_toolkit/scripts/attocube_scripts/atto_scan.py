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
import copy

import numpy as np
import time

from b26_toolkit.instruments import NI6259, NI9263, NI9402, PiezoController, MicrowaveGenerator, ANC300, B26PulseBlaster
from b26_toolkit.scripts.galvo_scan.galvo_scan import GalvoScan
from pylabcontrol.core import Parameter, Script
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from b26_toolkit.scripts.pulse_sequences.laser_pulses import LaserPiPulses
from b26_toolkit.scripts.set_laser import SetAttoSingleAxis
from b26_toolkit.scripts.attocube_scripts.atto_offset_z import AttoOffsetZ


class AttoScanPiezoController(GalvoScan):
    """
    Scans the Attocube and record the fluorescence while applying a fixed MW tone, to map out magnetic field contours
    Based on GalvoScan, but sends the analog voltages from the DAQ to a piezo controller which then goes to the Attocube
    Requires the piezo controller instrument to read the voltage gain setting, and adjust the analog voltage accordingly
    """

    _DEFAULT_SETTINGS = [
        Parameter('point_a',
                  [Parameter('x', 0, float, 'x-coordinate'),
                   Parameter('y', 0, float, 'y-coordinate')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 1.0, float, 'x-coordinate'),
                   Parameter('y', 1.0, float, 'y-coordinate')
                   ]),
        Parameter('RoI_mode', 'corner', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                                               corner: pta and ptb are diagonal corners of rectangle.\n \
                                                               center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 64, int, 'number of x points to scan'),
                   Parameter('y', 64, int, 'number of y points to scan')
                   ]),
        Parameter('time_per_pt', .002, float, 'time in s to measure at each point'),
        Parameter('settle_time', .0002, float, 'wait time between points to allow galvo to settle'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('DAQ_channels',
                  [Parameter('x_ao_channel', 'ao2', ['ao0', 'ao1', 'ao2', 'ao3'],
                             'Daq channel used for x voltage analog output'),
                   Parameter('y_ao_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'],
                             'Daq channel used for y voltage analog output'),
                   Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'],
                             'Daq channel used for counter')
                   ]),
        Parameter('ending_behavior', 'leave_at_corner', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan'),
        Parameter('microwaves',
                  [Parameter('enable', False, bool, 'enable microwaves'),
                   Parameter('power_out', -45.0, float, 'output power (dBm)'),
                   Parameter('freq', 2.82e9, float, 'frequency of MW'),
                   Parameter('normalize_scan', False, bool, 'Normalize scan by taking each line twice, with MW off and on'),
                   Parameter('turn_off_after', False, bool, 'If true MW output is turned off after the measurement'),
                   Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')])
    ]
    _INSTRUMENTS = {'piezo_controller': PiezoController, 'microwave_generator': MicrowaveGenerator,
                    'NI6259':  NI6259, 'NI9263': NI9263, 'NI9402': NI9402}

    def check_bounds(self):
        if np.any(self.x_array < 0) or np.any(self.y_array < 0):
            self.log('Attocube cannot accept negative voltages!')
            raise AttributeError

        # Spec says not to send in >1V/ms to DC-in, which is equiv to asking ANC300 to output 15V/ms
        # In our code that's equiv to speed = 15000
        self.speed_max = 1*self.scale()[0]/0.001
        self.speed = (np.max(self.x_array) - np.min(self.x_array)) * self.scale()[0] / float(self.settings['num_points']['x']) / self.settings['time_per_pt']
        self.speed_ratio = self.speed/self.speed_max

        if self.speed > self.speed_max * .5:
            self.log('Scan speed of %.1f V/s exceeds safety limit!'%self.speed)
            raise AttributeError
        else:
            self.log('Scan speed: %.1f V/s'%self.speed)

    def scale(self):
        voltage_limit = int(self.instruments['piezo_controller']['instance'].read_probes('voltage_limit'))
        if voltage_limit == 75:
            scale = [7.5, 7.5]
        elif voltage_limit == 100:
            scale = [10, 10]
        elif voltage_limit == 150:
            scale = [15, 15]
        return scale

    def before_scan(self):
        if self.settings['microwaves']['enable']:
            self.instruments['microwave_generator']['instance'].update({'frequency': float(self.settings['microwaves']['freq'])})
            self.instruments['microwave_generator']['instance'].update({'amplitude': self.settings['microwaves']['power_out']})
            self.instruments['microwave_generator']['instance'].update({'enable_modulation': True})
            self.instruments['microwave_generator']['instance'].update({'enable_output': True})
        else:
            self.instruments['microwave_generator']['instance'].update({'enable_output': False})

    def after_scan(self):
        if self.settings['microwaves']['turn_off_after']:
            self.instruments['microwave_generator']['instance'].update({'enable_output': False})

    def read_line_wrapper(self, y_pos):
        """
        If normalization is on, call read_line twice, with MW off and on.
        Args:
            y_pos: y position of the scan

        Returns:

        """
        if self.settings['microwaves']['normalize_scan'] and self.settings['microwaves']['enable']:
            self.instruments['microwave_generator']['instance'].update(
                {'frequency': float(self.settings['microwaves']['freq'])})
            #self.instruments['microwave_generator']['instance'].update({'enable_output': True})
            #time.sleep(self.settings['settle_time'] * self.settings['num_points']['x'])
            data_MW_on = self.read_line(y_pos)
            self.instruments['microwave_generator']['instance'].update(
                {'frequency': float(2.87e9)})
            #self.instruments['microwave_generator']['instance'].update({'enable_output': False})
            #time.sleep(self.settings['settle_time']*self.settings['num_points']['x'])
            data_MW_off = self.read_line(y_pos)

            line_data = (data_MW_on - data_MW_off)/data_MW_off
        else:
            line_data = self.read_line(y_pos)
        return line_data

    def before_line_scan(self):
        pass

    def read_line(self, y_pos):
        self.initPt = [self.x_array[0], y_pos]
        self.daq_out.set_analog_voltages(
            {self.settings['DAQ_channels']['x_ao_channel']: self.initPt[0],
             self.settings['DAQ_channels']['y_ao_channel']: self.initPt[1]})

        self.before_line_scan()

        # initialize APD thread
        ctrtask = self.daq_in.setup_counter(
            self.settings['DAQ_channels']['counter_channel'],
            len(self.x_array) + 1)
        aotask = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel']],
                                       self.x_array, ctrtask)

        # start counter and scanning sequence
        self.daq_out.run(aotask)
        self.daq_in.run(ctrtask)
        self.daq_out.waitToFinish(aotask)
        self.daq_out.stop(aotask)
        xLineData, _ = self.daq_in.read(ctrtask)
        self.daq_in.stop(ctrtask)
        diffData = np.diff(xLineData)

        # Create a list of 0s with the same length as the number of points in x dir requested
        summedData = np.zeros(int(len(self.x_array) / self.clockAdjust))
        for i in range(len(summedData)):
            summedData[i] = np.sum(
                diffData[(i * self.clockAdjust + 1):((i+1) * self.clockAdjust - 1)])

        #################
        # If forward scan speed is x of the max speed, do the reverse scan at 0.1 max speed
        if self.speed_ratio < 0.04:
            reverse_speedup = int(1/self.speed_ratio*0.04)
        else:
            reverse_speedup = 1

        sample_rate = float(1) / self.settings['settle_time'] * reverse_speedup

        if sample_rate > 90000:
            # Lower speed up rate if the required sampling rate is too high for the DAQ to handle
            sample_rate = 90000
            reverse_speedup = sample_rate * self.settings['settle_time']

        #print("Moving back to start of x_line with %i X forward speed" % reverse_speedup)

        # Change sample rate for speed up
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate
        self.daq_in.settings['digital_input'][
            self.settings['DAQ_channels']['counter_channel']]['sample_rate'] = sample_rate


        # Return to beginning of the line
        # initialize APD thread
        self.x_array_reverse = self.x_array[::-1]
        clktask = self.daq_in.setup_clock(
            self.settings['DAQ_channels']['counter_channel'],
            len(self.x_array_reverse) + 1)

        # move in the opposite direction
        _aotask = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel']],
                                        self.x_array_reverse, clktask)
        # start counter and scanning sequence
        self.daq_out.run(_aotask)
        self.daq_in.run(clktask)
        self.daq_out.waitToFinish(_aotask)
        self.daq_out.stop(_aotask)
        #_xLineData_reverse, _rev = self.daq_in.read(_aitask)  # read the ai data for this scan
        #self.daq_in.stop(_aitask)
        self.daq_in.stop(clktask)

        # Revert change to sample rate
        sample_rate = float(1) / self.settings['settle_time']
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate
        self.daq_in.settings['digital_input'][
            self.settings['DAQ_channels']['counter_channel']]['sample_rate'] = sample_rate
        #################

        # also normalizing to kcounts/sec
        return summedData * (.001 / self.time_per_pt_actual())

    def time_per_pt_actual(self):
        # Actual time spent measuring counts at each point. Normally the same as time_per_pt setting, but might be much shorter if using pulse sequences
        return self.settings['time_per_pt']

    def _plot(self, axes_list, data = None):
        """
        Plots the galvo scan image
        Args:
            axes_list: list of axes objects on which to plot the galvo scan on the first axes object
            data: data (dictionary that contains keys image_data, extent) if not provided use self.data
        """

        if data is None:
            data = self.data

        labels = ['B-field Contour Scan (f = %.3e Hz)' % self.instruments['microwave_generator']['instance'].settings['frequency'], r'V$_x$ [V]', r'V$_y$ [V]', '[kCt/s]']
        plot_fluorescence_new(data['image_data'], data['extent'], axes_list[0], max_counts=self.settings['max_counts_plot'], labels=labels)


class AttoScan(AttoScanPiezoController):
    """
    Scans the Attocube and record the fluorescence while applying a fixed MW tone, to map out magnetic field contours
    Based on GalvoScan, but instead of scanning the galvo using DAQ outputs, scans the Attocube using DAQ outouts
    """
    _INSTRUMENTS = {'microwave_generator': MicrowaveGenerator, 'NI6259': NI6259, 'NI9263': NI9263, 'NI9402': NI9402, 'ANC300': ANC300}

    def scale(self):
        # DC-in functionality has a built-in gain of 15
        return [15, 15]

    def before_scan(self):
        self.instruments['ANC300']['instance']._set_mode(1, 'input')
        self.instruments['ANC300']['instance']._set_mode(2, 'input')
        super(AttoScan, self).before_scan()


class AttoScanPulsed(AttoScan):
    """
    """
    _DEFAULT_SETTINGS = [
        Parameter('point_a',
                  [Parameter('x', 0, float, 'x-coordinate'),
                   Parameter('y', 0, float, 'y-coordinate')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 1.0, float, 'x-coordinate'),
                   Parameter('y', 1.0, float, 'y-coordinate')
                   ]),
        Parameter('RoI_mode', 'corner', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                                                   corner: pta and ptb are diagonal corners of rectangle.\n \
                                                                   center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 64, int, 'number of x points to scan'),
                   Parameter('y', 64, int, 'number of y points to scan')
                   ]),
        Parameter('time_per_pt', .002, float, 'time in s to measure at each point'),
        Parameter('settle_time', .0002, float, 'wait time between points to allow galvo to settle'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('DAQ_channels',
                  [Parameter('x_ao_channel', 'ao2', ['ao0', 'ao1', 'ao2', 'ao3'],
                             'Daq channel used for x voltage analog output'),
                   Parameter('y_ao_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'],
                             'Daq channel used for y voltage analog output'),
                   Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'],
                             'Daq channel used for counter')
                   ]),
        Parameter('ending_behavior', 'leave_at_corner', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn'),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan'),

        Parameter('microwaves',
                  [Parameter('enable', False, bool, 'enable microwaves'),
                   Parameter('power_out', 'configure in subscript', ['configure in subscript'], 'output power (dBm)'),
                   Parameter('freq', 'configure in subscript', ['configure in subscript'], 'frequency of MW'),
                   Parameter('normalize_scan', False, bool, 'Normalize scan by taking each line twice, with MW off and on'),
                   Parameter('turn_off_after', False, bool, 'NOT IMPLEMENTED; if true MW output is turned off after the measurement'),
                   Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')])
    ]

    _INSTRUMENTS = {'microwave_generator': MicrowaveGenerator, 'NI6259': NI6259, 'NI9263': NI9263, 'NI9402': NI9402, 'ANC300': ANC300, 'PB': B26PulseBlaster,
                    'mw_gen': MicrowaveGenerator}
    _SCRIPTS = {'PulsedESRSingleBlind': LaserPiPulses}
    _ACQ_TYPE = 'line'

    def before_line_scan(self):
        self.scripts['PulsedESRSingleBlind'].run(verbose=False)

    def setup_scan(self):
        sequence_duration = self.scripts['PulsedESRSingleBlind']._create_pulse_sequences(get_duration=True) / 1e9
        line_duration = (self.settings['time_per_pt'] + self.settings['settle_time']) * self.settings['num_points']['x']

        # Repeat sequence to cover entire duration of scanning a single line, with the addition of a hard-coded time (0.3) to account for any delays
        num_sequences = int((line_duration + 0.4) / sequence_duration)
        self.scripts['PulsedESRSingleBlind'].settings['averaging_block_size'] = num_sequences
        self.log('Running %i sequences over line duration of %.2f s' % (num_sequences, line_duration))

        self.scripts['PulsedESRSingleBlind'].settings['num_averages'] = self.scripts['PulsedESRSingleBlind'].settings['averaging_block_size']
        #self.scripts['PulsedESRSingleBlind'].settings['freq_start'] = self.settings['microwaves']['freq']
        #self.scripts['PulsedESRSingleBlind'].settings['mw_power'] = self.settings['microwaves']['power_out']

        mw_duration = self.scripts['PulsedESRSingleBlind'].settings['tau_mw'] * self.scripts['PulsedESRSingleBlind'].settings['repetitions'] / 1e9
        mw_duty_cycle = mw_duration / sequence_duration

        equiv_cw_power = self.scripts['PulsedESRSingleBlind'].settings['mw_power'] + np.log10(mw_duty_cycle)*10
        self.log('MW duty cycle: %.2e; avg power: %.1f dBm' % (mw_duty_cycle, equiv_cw_power))
        self.is_valid()
        super(AttoScanPulsed, self).setup_scan()


    def time_per_pt_actual(self):
        return self.settings['time_per_pt'] * self.scripts['PulsedESRSingleBlind'].laser_duties[0]

    def before_scan(self):
        self.instruments['ANC300']['instance']._set_mode(1, 'input')
        self.instruments['ANC300']['instance']._set_mode(2, 'input')

        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        self.instruments['PB']['instance'].update({'laser': {'status': False}})

        if self.settings['microwaves']['enable']:
            # Running the script sets the MW frequency again, so there's actually no need to configure it before the scan
            # However, the plotting needs to know the MW frequency to show the correct plot title
            self.instruments['microwave_generator']['instance'].update({'frequency': self.scripts['PulsedESRSingleBlind'].settings['freq_start']})
            self.instruments['microwave_generator']['instance'].update({'enable_output': True})
        else:
            self.instruments['microwave_generator']['instance'].update({'enable_output': False})

        # Temporarily disable constant heat load function, which sends in
        self.constant_microwave_heat_load = copy.deepcopy(self.instruments['PB']['instance'].settings['constant_microwave_heat_load']['enable'])
        self.instruments['PB']['instance'].settings['constant_microwave_heat_load']['enable'] = False

    def after_scan(self):
        self.instruments['PB']['instance'].settings['constant_microwave_heat_load']['enable'] = self.constant_microwave_heat_load
        if self.instruments['PB']['instance'].settings['constant_microwave_heat_load']['enable']:
            self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
            self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
            self.instruments['mw_gen']['instance'].update({'amplitude': self.instruments['PB']['instance'].settings['constant_microwave_heat_load']['mw_power']})
            self.instruments['mw_gen']['instance'].update({'frequency': self.instruments['PB']['instance'].settings['constant_microwave_heat_load']['mw_frequency']})
            self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
            self.instruments['mw_gen']['instance'].update({'enable_output': True})
            self.instruments['PB']['instance'].mw_duty_cycle_loop()

    def is_valid(self, pulse_sequences=None, verbose=True):
        """
        Validates the pulse sequence, also displays the laser/MW duty cycle
        Returns: None
        """
        return self.scripts['PulsedESRSingleBlind'].is_valid()

class AttoScanZSlicePulsed(AttoScanPulsed):
    """
    """
    _DEFAULT_SETTINGS = [
        Parameter('point_a',
                  [Parameter('x', 0, float, 'x-coordinate'),
                   Parameter('z', 0, float, 'z-coordinate')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 1.0, float, 'x-coordinate'),
                   Parameter('z', 1.0, float, 'z-coordinate')
                   ]),
        Parameter('RoI_mode', 'corner', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                                                   corner: pta and ptb are diagonal corners of rectangle.\n \
                                                                   center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 64, int, 'number of x points to scan'),
                   Parameter('y', 64, int, 'number of y points to scan')
                   ]),
        Parameter('time_per_pt', .002, float, 'time in s to measure at each point'),
        Parameter('settle_time', .0002, float, 'wait time between points to allow galvo to settle'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('DAQ_channels',
                  [Parameter('x_ao_channel', 'ao2', ['ao0', 'ao1', 'ao2', 'ao3'],
                             'Daq channel used for x voltage analog output'),
                   Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'],
                             'Daq channel used for counter')
                   ]),
        Parameter('ending_behavior', 'leave_at_corner', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn'),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan'),

        Parameter('microwaves',
                  [Parameter('enable', False, bool, 'enable microwaves'),
                   Parameter('power_out', 'configure in subscript', ['configure in subscript'], 'output power (dBm)'),
                   Parameter('freq', 'configure in subscript', ['configure in subscript'], 'frequency of MW'),
                   Parameter('normalize_scan', False, bool, 'Normalize scan by taking each line twice, with MW off and on'),
                   Parameter('turn_off_after', False, bool, 'NOT IMPLEMENTED; if true MW output is turned off after the measurement'),
                   Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')])
    ]

    _INSTRUMENTS = {'microwave_generator': MicrowaveGenerator, 'NI6259': NI6259, 'NI9263': NI9263, 'NI9402': NI9402, 'ANC300': ANC300, 'PB': B26PulseBlaster}
    _SCRIPTS = {'PulsedESRSingleBlind': LaserPiPulses, 'AttoOffsetZ': AttoOffsetZ, 'SetAttoSingleAxis': SetAttoSingleAxis}
    _ACQ_TYPE = 'line'

    def scale(self):
        """
        Custom scaling for voltages
        Since X is sent as a DC voltage to the input port of the Attocube controller, it needs to be divided by the amplification factor (15)
        Z does not need to be scaled, since we're not using the input port to control Z. We're using AttoOffsetZ which communicates digitally with the controller
        Returns:
        1 as default
        """
        return [15, 1]

    def before_scan(self):
        self.instruments['ANC300']['instance']._set_mode(1, 'input')

        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        self.instruments['PB']['instance'].update({'laser': {'status': False}})

        if self.settings['microwaves']['enable']:
            # Running the script sets the MW frequency again, so there's actually no need to configure it before the scan
            # However, the plotting needs to know the MW frequency to show the correct plot title
            self.instruments['microwave_generator']['instance'].update({'frequency': self.scripts['PulsedESRSingleBlind'].settings['freq_start']})
            self.instruments['microwave_generator']['instance'].update({'enable_output': True})
        else:
            self.instruments['microwave_generator']['instance'].update({'enable_output': False})

    def read_line(self, z_pos):
        self.daq_out.set_analog_voltages({self.settings['DAQ_channels']['x_ao_channel']: float(self.x_array[0])})
        self.scripts['AttoOffsetZ'].settings['z_offset'] = float(z_pos)
        # print('Running AttoOffsetZ with z = %.1f' % z_pos)
        self.scripts['AttoOffsetZ'].run(verbose=False)

        self.before_line_scan()

        # initialize APD thread
        ctrtask = self.daq_in.setup_counter(
            self.settings['DAQ_channels']['counter_channel'],
            len(self.x_array) + 1)
        aotask = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel']],
                                       self.x_array, ctrtask)

        # start counter and scanning sequence
        self.daq_out.run(aotask)
        self.daq_in.run(ctrtask)
        self.daq_out.waitToFinish(aotask)
        self.daq_out.stop(aotask)
        xLineData, _ = self.daq_in.read(ctrtask)
        self.daq_in.stop(ctrtask)
        diffData = np.diff(xLineData)

        # Create a list of 0s with the same length as the number of points in x dir requested
        summedData = np.zeros(int(len(self.x_array) / self.clockAdjust))
        for i in range(len(summedData)):
            summedData[i] = np.sum(
                diffData[(i * self.clockAdjust + 1):((i+1) * self.clockAdjust - 1)])

        #################
        # If forward scan speed is x of the max speed, do the reverse scan at 0.04 max speed
        if self.speed_ratio < 0.04:
            reverse_speedup = int(1/self.speed_ratio*0.04)
        else:
            reverse_speedup = 1

        sample_rate = float(1) / self.settings['settle_time'] * reverse_speedup

        if sample_rate > 90000:
            # Lower speed up rate if the required sampling rate is too high for the DAQ to handle
            sample_rate = 90000
            reverse_speedup = sample_rate * self.settings['settle_time']

        #print("Moving back to start of x_line with %i X forward speed" % reverse_speedup)

        # Change sample rate for speed up
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
        self.daq_in.settings['digital_input'][
            self.settings['DAQ_channels']['counter_channel']]['sample_rate'] = sample_rate

        # Return to beginning of the line
        # initialize APD thread
        self.x_array_reverse = self.x_array[::-1]
        clktask = self.daq_in.setup_clock(
            self.settings['DAQ_channels']['counter_channel'],
            len(self.x_array_reverse) + 1)

        # move in the opposite direction
        _aotask = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel']],
                                        self.x_array_reverse, clktask)
        # start counter and scanning sequence
        self.daq_out.run(_aotask)
        self.daq_in.run(clktask)
        self.daq_out.waitToFinish(_aotask)
        self.daq_out.stop(_aotask)
        #_xLineData_reverse, _rev = self.daq_in.read(_aitask)  # read the ai data for this scan
        #self.daq_in.stop(_aitask)
        self.daq_in.stop(clktask)

        # Revert change to sample rate
        sample_rate = float(1) / self.settings['settle_time']
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
        self.daq_in.settings['digital_input'][
            self.settings['DAQ_channels']['counter_channel']]['sample_rate'] = sample_rate

        # also normalizing to kcounts/sec
        return summedData * (.001 / self.time_per_pt_actual())

    @staticmethod
    def pts_to_extent(pta, ptb, roi_mode):
        """

        Args:
            pta: point a
            ptb: point b
            roi_mode:   mode how to calculate region of interest
                        corner: pta and ptb are diagonal corners of rectangle.
                        center: pta is center and ptb is extend or rectangle

        Returns: extend of region of interest [xVmin, xVmax, yVmax, yVmin]

        """
        if roi_mode == 'corner':
            xVmin = min(pta['x'], ptb['x'])
            xVmax = max(pta['x'], ptb['x'])
            yVmin = min(pta['z'], ptb['z'])
            yVmax = max(pta['z'], ptb['z'])
        elif roi_mode == 'center':
            xVmin = pta['x'] - float(ptb['x']) / 2.
            xVmax = pta['x'] + float(ptb['x']) / 2.
            yVmin = pta['z'] - float(ptb['z']) / 2.
            yVmax = pta['z'] + float(ptb['z']) / 2.
        return [xVmin, xVmax, yVmax, yVmin]

    def _plot(self, axes_list, data = None):
        """
        Plots the galvo scan image
        Args:
            axes_list: list of axes objects on which to plot the galvo scan on the first axes object
            data: data (dictionary that contains keys image_data, extent) if not provided use self.data
        """

        if data is None:
            data = self.data

        labels = ['B-field Contour Scan Z Slice (f = %.3e Hz)' % self.instruments['microwave_generator']['instance'].settings['frequency'], r'V$_x$ [V]', r'V$_z$ [V]',
                  '[kCt/s]']
        plot_fluorescence_new(data['image_data'], data['extent'], axes_list[0], max_counts=self.settings['max_counts_plot'], labels=labels, aspect='auto')