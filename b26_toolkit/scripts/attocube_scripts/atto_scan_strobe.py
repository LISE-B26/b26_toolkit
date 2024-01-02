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

from b26_toolkit.instruments import NI6259, NI9263, NI9402, MicrowaveGenerator, ANC300, B26PulseBlaster
from b26_toolkit.scripts.galvo_scan.galvo_scan_strobe import GalvoScanStrobe
from pylabcontrol.core import Parameter, Script
from b26_toolkit.scripts.pulse_sequences.param_sweep.pulsed_galvo_line_scan import PulsedEsrLineScan


class AttoScanAnc300Pulsed(GalvoScanStrobe):
    """
    Use AttoScanStrobe instead! 'Pulsed' version sends a train of laser and pi-pulses in the background while leaving the readout window open throughout the entire scan
    Scans the Attocube and record the fluorescence while applying a fixed MW tone, to map out magnetic field contours
    Based on GalvoScan, but instead of scanning the galvo using DAQ outputs, scans the Attocube using DAQ outouts
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
                   Parameter('mw_power', -45.0, float, 'microwave power in dB'),
                   Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
                   Parameter('freq', 2.87e9, float, 'start frequency of scan in Hz'),
                   Parameter('normalize_scan', False, bool, 'Normalize scan by taking each line twice, with MW off and on')
                   ])
    ]

    _INSTRUMENTS = {'microwave_generator': MicrowaveGenerator, 'NI6259': NI6259, 'NI9263': NI9263, 'NI9402': NI9402, 'ANC300': ANC300, 'PB': B26PulseBlaster}
    _SCRIPTS = {'pulsed_line_scan': PulsedEsrLineScan}
    _ACQ_TYPE = 'grid'

    def setup_scan(self):
        """
        Setup the scan, i.e. identify the instruments and set up sample rate and such
        :return:
        """
        # defines which daqs contain the input and output based on user selection of daq interface
        if self.settings['daq_type'] == 'PCI':
            self.daq_out = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            self.daq_out = self.instruments['NI9263']['instance']

        # checks that requested daqs are actually physically present in the system
        device_list = NI6259.get_connected_devices()
        if not (self.daq_out.settings['device'] in device_list):
            self.log('The requested output daq ' + self.daq_out.settings[
                'device'] + ' is not connected to this computer. Possible daqs are '
                     + str(device_list) + '. Please choose one of these and try again.')
            raise AttributeError

        line_scan_settings = self.scripts['pulsed_galvo_line_scan'].settings
        line_scan_settings['read_out'] = self.settings['read_out']
        line_scan_settings['mw_generator_switching_time'] = self.settings['settle_time']
        line_scan_settings['freq_start'] = self.settings['point_a']['x']
        line_scan_settings['freq_stop'] = self.settings['point_b']['x']
        line_scan_settings['freq_points'] = self.settings['num_points']['x']
        line_scan_settings['randomize'] = False
        if self.settings['RoI_mode'] == 'center':
            line_scan_settings['range_type'] = 'center_range'
        elif self.settings['RoI_mode'] == 'corner':
            line_scan_settings['range_type'] = 'start_stop'

        sequence_duration = self.scripts['pulsed_galvo_line_scan']._create_pulse_sequences(get_duration=True)/1e9
        line_scan_settings['num_averages'] = int(self.settings['time_per_pt']/sequence_duration)
        line_scan_settings['averaging_block_size'] = line_scan_settings['num_averages']
        print('Running %i averages per pixel' % line_scan_settings['num_averages'])



    def is_valid(self, pulse_sequences=None, verbose=True):
        """
        Validates the pulse sequence, also displays the laser/MW duty cycle
        Returns: None
        """

        self.scripts['pulsed_galvo_line_scan'].settings['read_out'] = self.settings['read_out']
        return self.scripts['pulsed_galvo_line_scan'].is_valid()