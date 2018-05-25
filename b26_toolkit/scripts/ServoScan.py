import numpy as np

from pylabcontrol.core import Script, Parameter
from b26_toolkit.instruments import KDC001
from b26_toolkit.scripts import GalvoScan, ESR



class ServoScan(Script):
    """
    GalvoScan uses the apd, daq, and galvo to sweep across voltages while counting photons at each voltage,
    resulting in an image in the current field of view of the objective.
    """

    _DEFAULT_SETTINGS = [
        Parameter('loop_script', 'GalvoScan', ['GalvoScan', 'ESR'], 'Script to execute at each servo position'),
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

    _SCRIPTS = {'GalvoScan': GalvoScan, 'ESR': ESR}

    def _function(self):
        def _get_channel(channel_name):
            if self.settings[channel_name]['axis'] == 'None':
                return None, None
            elif self.settings[channel_name]['axis'] == 'x':
                return self.instruments['XServo'], None
            elif self.settings[channel_name]['axis'] == 'y':
                return self.instruments['YServo'], None
            elif self.settings[channel_name]['axis'] == 'z':
                return self.instruments['ZServo'], None
            elif self.settings[channel_name]['axis'] == 'xy':
                return self.instruments['XServo'], self.instruments['YServo']
            elif self.settings[channel_name]['axis'] == 'xz':
                return self.instruments['XServo'], self.instruments['ZServo']
            elif self.settings[channel_name]['axis'] == 'yz':
                return self.instruments['YServo'], self.instruments['ZServo']

        def _get_channel_positions(channel_name):
            if self.settings[channel_name]['axis'] == 'None':
                return([0], None)
            elif self.settings[channel_name]['axis'] in ['x', 'y', 'z']:
                return([np.linspace(self.settings['min_value'], self.settings['max_value'], self.settings['num_points'])], None)
            elif self.settings[channel_name]['axis'] in ['xy', 'xz', 'yz']:
                return([np.linspace(self.settings['min_value'], self.settings['max_value'], self.settings['num_points'])],
                      [np.linspace(self.settings['min_value2'], self.settings['max_value2'], self.settings['num_points2'])])

        channel1_servos = _get_channel('channel1')
        channel2_servos = _get_channel('channel2')
        channel3_servos = _get_channel('channel3')

        channel1_pos, channel1_pos2 = _get_channel_positions(self.settings['channel_1']['axis'])
        channel2_pos, channel2_pos2 = _get_channel_positions(self.settings['channel_2']['axis'])
        channel3_pos, channel3_pos2 = _get_channel_positions(self.settings['channel_3']['axis'])

        for index, pos1 in enumerate(channel1_pos):
            channel1_servos[0].position = pos1
            if channel1_servos[1]:
                channel1_servos[1].position = channel1_pos2[index]
            for index, pos2 in enumerate(channel2_pos):
                channel2_servos[0].position = pos2
                if channel2_servos[1]:
                    channel2_servos[1].position = channel2_pos2[index]
                for index, pos3 in enumerate(channel3_pos):
                    channel3_servos[0].position = pos3
                    if channel3_servos[1]:
                        channel3_servos[1].position = channel3_pos2[index]
                        script = self.settings['loop_script']
                        script.settings['tag'] = script.settings['tag'] + '_' + str(pos1) + '_' + str(pos2) + '_' + str(pos3)
                        self.scripts[self.settings['loop_script']].run()
