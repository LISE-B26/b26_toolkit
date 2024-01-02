import numpy as np
from b26_toolkit.instruments import NI6259, NI9263, NI9402, B26PulseBlaster
from b26_toolkit.scripts.galvo_scan.galvo_scan import GalvoScanSafe
from b26_toolkit.scripts.pulse_sequences.param_sweep.pulsed_galvo_line_scan import PulsedGalvoLineScan
from pylabcontrol.core import Parameter


class GalvoScanStrobe(GalvoScanSafe):
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
        Parameter('read_out', [
            Parameter('meas_time', 340, float, 'measurement time (in ns)'),
            Parameter('nv_reset_time', 800, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 80, int, 'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_readout', 260, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
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

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9263': NI9263, 'NI9402': NI9402, 'PB': B26PulseBlaster}
    _SCRIPTS = {'pulsed_line_scan': PulsedGalvoLineScan}
    _ACQ_TYPE = 'grid'

    def before_scan(self):
        """
        Runs something before starting the scan.
        This can be something like moving the laser to a certain position, running findNV, or turning on MW
        Returns:

        """
        self.instruments['PB']['instance'].update({'laser': {'status': False}})

    def before_line_scan(self):
        # self.scripts['LaserPulses'].run(verbose=False)
        pass

    def setup_scan(self):
        """
        setup the scan, i.e. identify the instruments and set up sample rate and such


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

        line_scan_settings = self.scripts['pulsed_line_scan'].settings
        line_scan_settings['read_out'] = self.settings['read_out']
        line_scan_settings['settle_time'] = self.settings['settle_time']
        line_scan_settings['voltage_start'] = self.settings['point_a']['x']
        line_scan_settings['voltage_stop'] = self.settings['point_b']['x']
        line_scan_settings['voltage_points'] = self.settings['num_points']['x']
        line_scan_settings['randomize'] = False
        if self.settings['RoI_mode'] == 'center':
            line_scan_settings['range_type'] = 'center_range'
        elif self.settings['RoI_mode'] == 'corner':
            line_scan_settings['range_type'] = 'start_stop'

        sequence_duration = self.scripts['pulsed_line_scan']._create_pulse_sequences(get_duration=True)/1e9
        line_scan_settings['num_averages'] = int(self.settings['time_per_pt']/sequence_duration)
        line_scan_settings['averaging_block_size'] = line_scan_settings['num_averages']
        print('Running %i averages per pixel' % line_scan_settings['num_averages'])

    def read_line(self, y_pos):

        self.before_line_scan()
        self.initPt = [self.x_array[0], y_pos]
        self.daq_out.set_analog_voltages(
            {self.settings['DAQ_channels']['x_ao_channel']: self.initPt[0],
             self.settings['DAQ_channels']['y_ao_channel']: self.initPt[1]})

        self.scripts['pulsed_line_scan'].run(verbose=False)
        normalizedData = self.scripts['pulsed_line_scan'].data['count_data']
        normalizedData = np.array(normalizedData)[:, 0]

        if 'safety_threshold' in self.settings:
            safety_threshold = self.settings['safety_threshold']
            if safety_threshold != -1 and any(pixel > safety_threshold for pixel in normalizedData):
                self.safety_threshold_exceeded = True
                self.safety_exceeded_behavior()
        return normalizedData

    def time_per_pt_actual(self):
        # Actual time spent measuring counts at each point. Normally the same as time_per_pt setting, but might be much shorter if using pulse sequences
        return self.settings['time_per_pt']

    def is_valid(self, pulse_sequences=None, verbose=True):
        """
        Validates the pulse sequence, also displays the laser/MW duty cycle
        Returns: None
        """

        self.scripts['pulsed_line_scan'].settings['read_out'] = self.settings['read_out']
        return self.scripts['pulsed_line_scan'].is_valid()

