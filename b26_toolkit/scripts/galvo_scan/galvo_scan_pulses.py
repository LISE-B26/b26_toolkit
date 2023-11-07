import numpy as np
import time

from b26_toolkit.instruments import NI6259, NI9263, NI9402, PiezoController, MicrowaveGenerator, ANC300, B26PulseBlaster, Pulse
from b26_toolkit.scripts.galvo_scan.galvo_scan import GalvoScan
from b26_toolkit.scripts.galvo_scan.galvo_scan import AttoScanANC300
from b26_toolkit.scripts.pulse_sequences.pulsed_esr import PulsedESRFastSingle, PulsedESRSingleBlind
from pylabcontrol.core import Script, Parameter
from b26_toolkit.scripts.set_laser import SetLaser
from b26_toolkit.scripts.pulse_sequences.stroboscopic_readout import StroboscopicReadout
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence

class GalvoScanStrobed(GalvoScan):
    _DEFAULT_SETTINGS = [
        Parameter('point_a',
                  [Parameter('x', 0, float, 'x-coordinate'),
                   Parameter('y', 0, float, 'y-coordinate')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 1.0, float, 'x-coordinate'),
                   Parameter('y', 1.0, float, 'y-coordinate')
                   ]),
        Parameter('RoI_mode', 'center', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                                                   corner: pta and ptb are diagonal corners of rectangle.\n \
                                                                   center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 64, int, 'number of x points to scan'),
                   Parameter('y', 64, int, 'number of y points to scan')
                   ]),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    _SCRIPTS = {'strobe_readout': StroboscopicReadout, 'SetLaser': SetLaser}

    def setup_scan(self):
        """
        setup the scan, i.e. identify the instruments and set up sample rate and such
        """
        self.scripts['SetLaser'].settings['daq_type'] = 'PCI'
        self.data['counts'] = []

    def read_line(self, y_pos):
        set_laser_script = self.scripts['SetLaser']
        set_laser_script.settings['point']['y'] = y_pos

        line_data = np.zeros(len(self.x_array))

        self.data['counts'].append([])

        for i in range(len(self.x_array)):
            print(i)
            if self._abort:
                break

            set_laser_script.settings['point']['x'] = self.x_array[i]
            set_laser_script.run()

            self.scripts['strobe_readout'].run()
            counts = self.scripts['strobe_readout'].data['counts']

            line_data[i] = np.mean(counts)
            self.data['counts'][-1].append(counts)
        return line_data

    def get_galvo_location(self):
        """
        Returns the current position of the galvo. Requires a daq with analog inputs internally routed to the analog
        outputs (ex. NI6259. Note that the cDAQ does not have this capability).
        Returns: list with two floats, which give the x and y position of the galvo mirror
        """
        point = self.scripts['SetLaser'].settings['point']
        return [point['x'], point['y']]

    def set_galvo_location(self, galvo_position):
        """
        sets the current position of the galvo
        galvo_position: list with two floats, which give the x and y position of the galvo mirror
        """
        if galvo_position[0] > 1 or galvo_position[0] < -1 or galvo_position[1] > 1 or galvo_position[1] < -1:
            raise ValueError('The script attempted to set the galvo position to an illegal position outside of +- 1 V')

        self.scripts['SetLaser'].settings['point']['x'] = galvo_position[0]
        self.scripts['SetLaser'].settings['point']['y'] = galvo_position[1]

class GalvoScanPoint(AttoScanANC300):
    """
    Get rid of Pulsed stuff here. This is still useful as GalvoPointScan
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
        Parameter('time_per_pt', .01, [.01, .015, .02, .025, .03, .035, .04, .05, .1],
                  'time in s to measure at each point'),
        Parameter('settle_time', .01, [.001, .005, .01], 'wait time between points to allow galvo to settle, \n'
                                                         'time_per_pt must be divisible by settle_time!'),
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
                   Parameter('power_out', -45.0, float, 'output power (dBm)'),
                   Parameter('freq', 2.82e9, float, 'frequency of MW'),
                   Parameter('normalize_scan', False, bool, 'Normalize scan by taking each line twice, with MW off and on'),
                   Parameter('turn_off_after', False, bool, 'NOT IMPLEMENTED; if true MW output is turned off after the measurement'),
                   Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')])
    ]

    _INSTRUMENTS = {'microwave_generator': MicrowaveGenerator, 'NI6259': NI6259, 'NI9263': NI9263, 'NI9402': NI9402, 'ANC300': ANC300}
    _SCRIPTS = {'PulsedESRFastSingle': PulsedESRFastSingle}
    _ACQ_TYPE = 'point'

    def read_point(self, x, y):
        self.daq_out.set_analog_voltages(
            {self.settings['DAQ_channels']['x_ao_channel']: x,
             self.settings['DAQ_channels']['y_ao_channel']: y})

        self.scripts['PulsedESRFastSingle'].run()
        esr_counts = self.scripts['PulsedESRFastSingle'].data['esr_counts']
        print(esr_counts)
        return esr_counts[0]


class AttoScanANC300Pulsed(AttoScanANC300):
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
        Parameter('time_per_pt', .01, [.01, .015, .02, .025, .03, .035, .04, .05, .1],
                  'time in s to measure at each point'),
        Parameter('settle_time', .01, [.001, .005, .01], 'wait time between points to allow galvo to settle, \n'
                                                         'time_per_pt must be divisible by settle_time!'),
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
                   Parameter('power_out', -45.0, float, 'output power (dBm)'),
                   Parameter('freq', 2.82e9, float, 'frequency of MW'),
                   Parameter('normalize_scan', False, bool, 'Normalize scan by taking each line twice, with MW off and on'),
                   Parameter('turn_off_after', False, bool, 'NOT IMPLEMENTED; if true MW output is turned off after the measurement'),
                   Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')])
    ]

    _INSTRUMENTS = {'microwave_generator': MicrowaveGenerator, 'NI6259': NI6259, 'NI9263': NI9263, 'NI9402': NI9402, 'ANC300': ANC300, 'PB': B26PulseBlaster}
    _SCRIPTS = {'PulsedESRSingleBlind': PulsedESRSingleBlind}
    _ACQ_TYPE = 'line'

    def before_line_scan(self):

        self.scripts['PulsedESRSingleBlind'].run(verbose=False)

    def setup_scan(self):
        sequence_duration = self.scripts['PulsedESRSingleBlind']._create_pulse_sequences(get_duration=True) / 1e9
        line_duration = (self.settings['time_per_pt'] + self.settings['settle_time']) * self.settings['num_points']['x']

        # Repeat sequence to cover entire duration of scanning a single line, with the addition of a hard-coded time (0.3) to account for any delays
        num_sequences = int((line_duration + 0.3) / sequence_duration)
        self.scripts['PulsedESRSingleBlind'].settings['averaging_block_size'] = num_sequences
        self.log('Running %i sequences over line duration of %.2f s' % (num_sequences, line_duration))

        self.scripts['PulsedESRSingleBlind'].settings['num_averages'] = self.scripts['PulsedESRSingleBlind'].settings['averaging_block_size']
        self.scripts['PulsedESRSingleBlind'].settings['freq_start'] = self.settings['microwaves']['freq']
        self.scripts['PulsedESRSingleBlind'].settings['mw_power'] = self.settings['microwaves']['power_out']

        mw_duration = self.scripts['PulsedESRSingleBlind'].settings['tau_mw'] * self.scripts['PulsedESRSingleBlind'].settings['repetitions'] / 1e9
        mw_duty_cycle = mw_duration / sequence_duration

        equiv_cw_power = self.scripts['PulsedESRSingleBlind'].settings['mw_power'] + np.log10(mw_duty_cycle)*10
        self.log('MW duty cycle: %.2e; equiv. CW power: %.1f dBm' % (mw_duty_cycle, equiv_cw_power))
        super(AttoScanANC300Pulsed, self).setup_scan()


    def time_per_pt_actual(self):
        sequence_duration = self.scripts['PulsedESRSingleBlind']._create_pulse_sequences(get_duration=True)
        laser_duration = self.scripts['PulsedESRSingleBlind'].settings['read_out']['nv_reset_time']*self.scripts['PulsedESRSingleBlind'].settings['repetitions']
        duty_cycle = laser_duration/sequence_duration # Proportion of total time when laser is on
        return self.settings['time_per_pt']*duty_cycle

    def before_scan(self):
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
        super(AttoScanANC300Pulsed, self).before_scan()
