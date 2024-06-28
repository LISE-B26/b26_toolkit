from pylabcontrol.core import Parameter
from b26_toolkit.instruments import NI6259, NI9263, NI9402, MicrowaveGenerator, B26PulseBlaster
from b26_toolkit.scripts.galvo_scan.galvo_scan import GalvoScanSafe
from b26_toolkit.scripts.pulse_sequences.laser_pulses import LaserPulses

class GalvoScanPulsed(GalvoScanSafe):
    """
    Uses the DAQ analog outputs to sweep the voltages sent to the galvo mirrors, and collect the photon count at each voltage.
    This results in an image in the current field of view of the objective.
    Unlike GalvoScan, which leaves the laser on during the entire scan, this script sends in laser pulses
    Readout window is left on during the whole time (to minimize DAQ overhead), so for long time_per_pt, dark counts can wash out the fluorescence signal
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
        Parameter('RoI_mode', 'center', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                                               corner: pta and ptb are diagonal corners of rectangle.\n \
                                                               center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 64, int, 'number of x points to scan'),
                   Parameter('y', 64, int, 'number of y points to scan')
                   ]),
        Parameter('time_per_pt', .002, float, 'time in s to measure at each point'),
        Parameter('laser_pulse_duration', 500, float, 'duration of each laser pulse (ns)'),
        Parameter('laser_off_time', 500, float, 'time between each laser pulse (ns'),
        Parameter('settle_time', .0002, float, 'wait time between points to allow galvo to settle'),
        Parameter('safety_threshold', 10, float, 'PulseBlaster will turn off the laser when a line contains a pixel exceeding this threshold (kCt/s); '
                                                 'set to -1 to disable this feature'),
        Parameter('turn_off_laser_after', True, bool, 'turn off laser after scan'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('DAQ_channels',
                  [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
                   Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
                   Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'], 'Daq channel used for counter')
                   ]),
        Parameter('ending_behavior', 'leave_at_corner', ['return_to_start', 'return_to_origin', 'leave_at_corner'], 'return to the corn'),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    _INSTRUMENTS = {'microwave_generator': MicrowaveGenerator, 'NI6259': NI6259, 'NI9263': NI9263, 'NI9402': NI9402, 'PB': B26PulseBlaster}
    _SCRIPTS = {'LaserPulses': LaserPulses}
    _ACQ_TYPE = 'line'

    def before_scan(self):
        """
        Runs something before starting the scan.
        This can be something like moving the laser to a certain position, running findNV, or turning on MW
        Returns:

        """
        self.instruments['PB']['instance'].update({'laser': {'status': False}})

    def before_line_scan(self):
        self.scripts['LaserPulses'].run(verbose=False)

    def setup_scan(self):
        sequence_duration = self.scripts['LaserPulses']._create_pulse_sequences(get_duration=True) / 1e9
        line_duration = (self.settings['time_per_pt'] + self.settings['settle_time']) * self.settings['num_points']['x']

        # Repeat sequence to cover entire duration of scanning a single line, with the addition of a hard-coded time (0.3) to account for any delays
        num_sequences = int((line_duration+.3) / sequence_duration)
        self.scripts['LaserPulses'].settings['averaging_block_size'] = num_sequences
        self.log('Running %i sequences over line duration of %.2f s' % (num_sequences, line_duration))

        self.scripts['LaserPulses'].settings['laser_pulse_duration'] = self.scripts['LaserPulses'].settings['laser_pulse_duration']
        self.scripts['LaserPulses'].settings['laser_off_time'] = self.scripts['LaserPulses'].settings['laser_off_time']
        self.scripts['LaserPulses'].settings['num_averages'] = self.scripts['LaserPulses'].settings['averaging_block_size']
        self.is_valid()
        super(GalvoScanPulsed, self).setup_scan()

    def time_per_pt_actual(self):
        return self.settings['time_per_pt']*self.scripts['LaserPulses'].laser_duties[0]

    def is_valid(self, pulse_sequences=None, verbose=True):
        """
        Validates the pulse sequence, also displays the laser/MW duty cycle
        Returns: None
        """
        self.scripts['LaserPulses'].settings['laser_pulse_duration'] = self.settings['laser_pulse_duration']
        self.scripts['LaserPulses'].settings['laser_off_time'] = self.settings['laser_off_time']
        return self.scripts['LaserPulses'].is_valid()

