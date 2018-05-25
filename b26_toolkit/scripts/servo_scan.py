import numpy as np

from pylabcontrol.core import Script, Parameter
from b26_toolkit.instruments import KDC001
from b26_toolkit.scripts import GalvoScan



class ServoScan(Script):
    """
    Sweeps the three axis stage in a given pattern and runs the selected loop_script at each position
    """

    _DEFAULT_SETTINGS = [
        Parameter('loop_script', 'GalvoScan', ['GalvoScan'], 'Script to execute at each servo position'),
        Parameter('channel1',
                  [Parameter('axis', 'z', ['None', 'x', 'y', 'z', 'xy', 'xz', 'yz'], 'axis to scan on'),
                   Parameter('min_value', 0.0, float, 'starting point of scan'),
                   Parameter('max_value', 1.0, float, 'ending point of scan'),
                   Parameter('num_points', 10, int, 'number of points for scan'),
                   Parameter('min_value2', 0.0, float, 'starting point of scan (second dimension)'),
                   Parameter('max_value2', 1.0, float, 'ending point of scan (second dimension)'),
                   Parameter('num_points2', 10, int, 'number of points for scan (second dimension)')
                  ]),
        Parameter('channel2',
                  [Parameter('axis', 'y', ['None', 'x', 'y', 'z', 'xy', 'xz', 'yz'], 'axis to scan on'),
                   Parameter('min_value', 0.0, float, 'starting point of scan'),
                   Parameter('max_value', 1.0, float, 'ending point of scan'),
                   Parameter('num_points', 10, int, 'number of points for scan'),
                   Parameter('min_value2', 0.0, float, 'starting point of scan (second dimension)'),
                   Parameter('max_value2', 1.0, float, 'ending point of scan (second dimension)'),
                   Parameter('num_points2', 10, int, 'number of points for scan (second dimension)')
                   ]),
        Parameter('channel3',
                  [Parameter('axis', 'x', ['None', 'x', 'y', 'z', 'xy', 'xz', 'yz'], 'axis to scan on'),
                   Parameter('min_value', 0.0, float, 'starting point of scan'),
                   Parameter('max_value', 1.0, float, 'ending point of scan'),
                   Parameter('num_points', 10, int, 'number of points for scan'),
                   Parameter('min_value2', 0.0, float, 'starting point of scan (second dimension)'),
                   Parameter('max_value2', 1.0, float, 'ending point of scan (second dimension)'),
                   Parameter('num_points2', 10, int, 'number of points for scan (second dimension)')
                   ])
    ]

    _INSTRUMENTS = {'XServo': KDC001, 'YServo': KDC001, 'ZServo': KDC001}

    _SCRIPTS = {'GalvoScan': GalvoScan}

    def _function(self):
        """
        Runs the servo scan
        """
        #when instruments are initially loaded after script is imported, all will have the same serial number and two
        #will fail to connect. This instead pushes the serial numbers from the script dropdown and forces a reconnect
        #so after this step, all three should be connected properly.
        for key in list(self.instruments.keys()):
            self.instruments[key]['instance'].update({'serial_number': self.instruments[key]['settings']['serial_number']})

        def _get_channel(channel_name):
            """
            Takes a channel name and returns the zero, one, or two servo instrument instances corresponding to the
            movement for that channel
            Args:
                channel_name: 'channel1', 'channel2', or 'channel3'

            Returns: two servo instrument instances (or None) correponding to the servo(s) used to move in the specified
            direction
            """
            if self.settings[channel_name]['axis'] == 'None':
                return None, None
            elif self.settings[channel_name]['axis'] == 'x':
                return self.instruments['XServo']['instance'], None
            elif self.settings[channel_name]['axis'] == 'y':
                return self.instruments['YServo']['instance'], None
            elif self.settings[channel_name]['axis'] == 'z':
                return self.instruments['ZServo']['instance'], None
            elif self.settings[channel_name]['axis'] == 'xy':
                return self.instruments['XServo']['instance'], self.instruments['YServo']['instance']
            elif self.settings[channel_name]['axis'] == 'xz':
                return self.instruments['XServo']['instance'], self.instruments['ZServo']['instance']
            elif self.settings[channel_name]['axis'] == 'yz':
                return self.instruments['YServo']['instance'], self.instruments['ZServo']['instance']

        def _get_channel_positions(channel_name):
            """
            Returns a list of the servo positions for the scan
            Args:
                channel_name: 'channel1', 'channel2', or 'channel3'

            Returns: Two values, either [0], None for an empty channel, [position list], None for a linear scan, or
            [position list], [position list2] for a diagonal scan

            """
            if self.settings[channel_name]['axis'] == 'None':
                return([0], None)
            elif self.settings[channel_name]['axis'] in ['x', 'y', 'z']:
                return(np.linspace(self.settings[channel_name]['min_value'], self.settings[channel_name]['max_value'], self.settings[channel_name]['num_points']), None)
            elif self.settings[channel_name]['axis'] in ['xy', 'xz', 'yz']:
                return(np.linspace(self.settings[channel_name]['min_value'], self.settings[channel_name]['max_value'], self.settings[channel_name]['num_points']),
                      np.linspace(self.settings[channel_name]['min_value2'], self.settings[channel_name]['max_value2'], self.settings[channel_name]['num_points2']))

        channel1_servos = _get_channel('channel1')
        channel2_servos = _get_channel('channel2')
        channel3_servos = _get_channel('channel3')

        channel1_pos, channel1_pos2 = _get_channel_positions('channel1')
        channel2_pos, channel2_pos2 = _get_channel_positions('channel2')
        channel3_pos, channel3_pos2 = _get_channel_positions('channel3')

        for index, pos1 in enumerate(channel1_pos):
            if channel1_servos[0]:
                channel1_servos[0].position = float(pos1)
                if channel1_servos[1]:
                    channel1_servos[1].position = float(channel1_pos2[index])
            for index, pos2 in enumerate(channel2_pos):
                if channel2_servos[0]:
                    channel2_servos[0].position = float(pos2)
                    if channel2_servos[1]:
                        channel2_servos[1].position = float(channel2_pos2[index])
                for index, pos3 in enumerate(channel3_pos):
                    if self._abort:
                        return
                    if channel3_servos[0]:
                        channel3_servos[0].position = float(pos3)
                        if channel3_servos[1]:
                            channel3_servos[1].position = float(channel3_pos2[index])
                        self.settings['tag'] = self.settings['tag'] + '_' + str(pos1) + '_' + str(pos2) + '_' + str(pos3)
                        self.scripts[self.settings['loop_script']].run()

    def plot(self, axes_list):
        if self.scripts[self.settings['loop_script']].is_running:
            self.scripts[self.settings['loop_script']].plot(axes_list)