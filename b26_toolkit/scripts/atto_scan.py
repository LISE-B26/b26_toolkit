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

from b26_toolkit.instruments import Attocube
from pylabcontrol.core import Parameter, Script


class AttoStep(Script):
    # COMMENT_ME
    _DEFAULT_SETTINGS = [
        Parameter('axis', 'z', ['x', 'y', 'z'], 'Axis to step on'),
        Parameter('direction', 'Up', ['Up', 'Down'], 'step direction, up or down in voltage (or on physical switch)')
    ]

    _INSTRUMENTS = {'attocube': Attocube}
    _SCRIPTS = {}

    def __init__(self, instruments = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Default script initialization
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, log_function= log_function, data_path = data_path)

    def _function(self):
        """
        Performs a single attocube step with the voltage and frequency, and in the direction, specified in settings
        """
        attocube = self.instruments['attocube']['instance']
        attocube_voltage = self.instruments['attocube']['settings'][self.settings['axis']]['voltage']
        attocube.update({self.settings['axis']: {'voltage': attocube_voltage}})
        attocube_freq = self.instruments['attocube']['settings'][self.settings['axis']]['freq']
        attocube.update({self.settings['axis']: {'freq': attocube_freq}})
        if self.settings['direction'] == 'Up':
            dir = 0
        elif self.settings['direction'] == 'Down':
            dir = 1
        self.instruments['attocube']['instance'].step(self.settings['axis'], dir)

if __name__ == '__main__':
    script, failed, instr = Script.load_and_append({'AttoStep': 'AttoStep'})

    print(script)
    print(failed)
    print(instr)
    # fp = Find_Points(settings={'path': 'Z:/Lab/Cantilever/Measurements/__tmp__', 'tag':'nvs'})
    # fp.run()

    # plt.pcolor(fp.data['image'])
    # print(fp.data['image_gaussian'].shape)
    # plt.pcolor(fp.data['image'])
    # plt.imshow(fp.data['image'], cmap = 'pink', interpolation = 'nearest')
    #
    #
    # for x in fp.data['NV_positions']:
    #     plt.plot(x[0],x[1],'ro')
    #
    # plt.show()

    # plt.figure()
    # plt.imshow(fp.data['image_gaussian'])
    # Axes3D.plot(fp.data['image_gaussian'])
    # plt.show()
    # print(max(fp.data['image']))
    # print(max(fp.data['image_gaussian'].flatten()))
    # print('NV_positions', fp.data['NV_positions'])
