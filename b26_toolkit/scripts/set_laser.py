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
from matplotlib import patches

from b26_toolkit.instruments import NI6259, NI9263, PiezoController, NI9263_02, ANC300
from pylabcontrol.core import Script, Parameter
import time

class SetLaser(Script):
    """
    Sets the voltages sent by the DAQ analog outputs to the galvo mirrors. Typically used to move the galvo such that the laser spot is on a point of interest.
    """

    _DEFAULT_SETTINGS = [
        Parameter('point',
                  [Parameter('x', -0.4, float, 'x-coordinate'),
                   Parameter('y', -0.4, float, 'y-coordinate')
                   ]),
        Parameter('DAQ_channels',
            [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
            Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output')
            ]),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    _INSTRUMENTS = {'NI6259': NI6259, 'NI9263': NI9263}

    _SCRIPTS = {}

    def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_function= log_function, data_path = data_path)

        self._save_fig = False

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        if 'ANC300' in self.instruments:
            # Seems like ANC300 has weird behavior when it gets numbers like 16.666666
            #print('Rounding to position 2 dp')
            pt = (round(self.settings['point']['x'], 2) / self.scale(), round(self.settings['point']['y'], 2) / self.scale())
        else:
            pt = (self.settings['point']['x']/self.scale(), self.settings['point']['y']/self.scale())
            #print('NOT rounding position to 2 dp')
        try:
            self.check_bounds(pt[0], pt[1])
        except AttributeError:
            return

        # daq API only accepts either one point and one channel or multiple points and multiple channels
        pt = np.transpose(np.column_stack((pt[0], pt[1])))
        pt = (np.repeat(pt, 2, axis=1))

        self._setup_daq()

        t_start = time.time()
        task = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], pt)
        self.daq_out.run(task)
        self.daq_out.waitToFinish(task)
        self.daq_out.stop(task)
        self.log('laser set to Vx={:.4f}, Vy={:.4f}'.format(self.settings['point']['x'], self.settings['point']['y']))
        print(time.time()-t_start)


        # if self.settings['daq_type'] == 'PCI':
        #     self.daq_out = self.instruments['NI6259']['instance']
        # elif self.settings['daq_type'] == 'cDAQ':
        #     if 'NI9263' in self.instruments:
        #         self.daq_out = self.instruments['NI9263']['instance']
        #     elif 'NI9263_02' in self.instruments:
        #         print('using ni9263_02')
        #         self.daq_out = self.instruments['NI9263_02']['instance']

    def scale(self):
        return 1

    def check_bounds(self, x, y):
        pass

    def _setup_daq(self):
        # defines which daqs contain the input and output based on user selection of daq interface
        if self.settings['daq_type'] == 'PCI':
            self.daq_out = self.instruments['NI6259']['instance']
        elif self.settings['daq_type'] == 'cDAQ':
            if 'NI9263' in self.instruments:
                self.daq_out = self.instruments['NI9263']['instance']
            elif 'NI9263_02' in self.instruments:
                print('using ni9263_02')
                self.daq_out = self.instruments['NI9263_02']['instance']

    def get_galvo_position(self):
        """
        reads the current position from the x and y channels and returns it
        Returns:

        """
        if self.settings['daq_type'] == 'PCI':
            galvo_position = self.daq_out.get_analog_voltages([
                self.settings['DAQ_channels']['x_ao_channel'],
                self.settings['DAQ_channels']['y_ao_channel']]
            )
        elif self.settings['daq_type'] == 'cDAQ':
            print("WARNING cDAQ doesn't allow to read values")
            galvo_position = []

        return galvo_position
    #must be passed figure with galvo plot on first axis
    def plot(self, figure_list):
        axes_Image = figure_list[0].axes[0]

        # removes patches
        [child.remove() for child in axes_Image.get_children() if isinstance(child, patches.Circle)]
        xticks = axes_Image.get_xticks()
        radius = (xticks[-1]-xticks[0])*0.01


        patch = patches.Circle((self.settings['point']['x'], self.settings['point']['y']), radius, fc='r', alpha = .8)
        axes_Image.add_patch(patch)


class SetLaser_cDAQ(SetLaser):
    """
    Sets the voltages sent by the DAQ analog outputs to the galvo mirrors. Typically used to move the galvo such that the laser spot is on a point of interest.
    """

    _DEFAULT_SETTINGS = [
        Parameter('point',
                  [Parameter('x', 0, float, 'x-coordinate'),
                   Parameter('y', 0, float, 'y-coordinate')
                   ]),
        Parameter('DAQ_channels',
            [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
            Parameter('y_ao_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output')
            ])
    ]

    _INSTRUMENTS = {'daq_out':  NI9263}

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        pt = (self.settings['point']['x'], self.settings['point']['y'])

        # daq API only accepts either one point and one channel or multiple points and multiple channels
        pt = np.transpose(np.column_stack((pt[0],pt[1])))
        pt = (np.repeat(pt, 2, axis=1))

        task = self.instruments['daq_out']['instance'].setup_AO([self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], pt)
        self.instruments['daq_out']['instance'].run(task)
        self.instruments['daq_out']['instance'].waitToFinish(task)
        self.instruments['daq_out']['instance'].stop(task)
        self.log_message()

    def log_message(self):
        """
        Displays a log message after running the script
        """
        self.log('laser set to Vx={:.4}, Vy={:.4}'.format(self.settings['point']['x'], self.settings['point']['y']))


class SetAttoPiezoController(SetLaser):
    """
    Sets the voltages sent by the DAQ analog outputs to a piezo controller, which typically controls an Attocube. Basically same as SetLaser, but this script
    also requires the piezo controller instrument to read out its voltage limits/scaling factor.
    """
    _DEFAULT_SETTINGS = [
        Parameter('point',
                  [Parameter('x', 0., float, 'x-coordinate'),
                   Parameter('y', 0., float, 'y-coordinate')
                   ]),
        Parameter('DAQ_channels',
                  [Parameter('x_ao_channel', 'ao2', ['ao0', 'ao1', 'ao2', 'ao3'],
                             'Daq channel used for x voltage analog output'),
                   Parameter('y_ao_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'],
                             'Daq channel used for y voltage analog output')
                   ]),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    _INSTRUMENTS = {'piezo_controller': PiezoController, 'NI6259':  NI6259, 'NI9263': NI9263}

    def check_bounds(self, x, y):
        if x < 0 or y < 0:
            self.log('Piezo cannot accept negative voltages!')
            raise AttributeError

    def scale(self):
        voltage_limit = int(self.instruments['piezo_controller']['instance'].read_probes('voltage_limit'))
        #print(voltage_limit)
        if voltage_limit == 75:
            scale = 7.5
        elif voltage_limit == 100:
            scale = 10
        elif voltage_limit == 150:
            scale = 15
        return scale

    def log_message(self):
        """
        Displays a log message after running the script
        """
        self.log('Attocube set to Vx={:.4}, Vy={:.4}'.format(self.settings['point']['x'], self.settings['point']['y']))


class SetAtto(SetAttoPiezoController):
    """
    Sets the voltages sent by the DAQ analog outputs to an Attocube controller. Basically same as SetLaser, but this script
    also requires the Attocube controller instrument to configure its outputs (e.g. enable DC input mode to accept DAQ voltages)
     """
    _DEFAULT_SETTINGS = [
        Parameter('point',
                  [Parameter('x', 0., float, 'x-coordinate'),
                   Parameter('y', 0., float, 'y-coordinate')
                   ]),
        Parameter('zero_first', False, bool, 'Zero piezo inputs before moving to setpoint, to avoid hysteresis effects'),
        Parameter('DAQ_channels',
                  [Parameter('x_ao_channel', 'ao2', ['ao0', 'ao1', 'ao2', 'ao3'],
                             'Daq channel used for x voltage analog output'),
                   Parameter('y_ao_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'],
                             'Daq channel used for y voltage analog output')
                   ]),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]
    _INSTRUMENTS = {'NI6259': NI6259, 'NI9263': NI9263, 'ANC300': ANC300}

    def scale(self):
        return 15

    def _function(self):
        def set_ao(x, y):
            pt = (x / self.scale(), y / self.scale())
            try:
                self.check_bounds(pt[0], pt[1])
            except AttributeError:
                return

            # daq API only accepts either one point and one channel or multiple points and multiple channels
            pt = np.transpose(np.column_stack((pt[0], pt[1])))
            pt = (np.repeat(pt, 2, axis=1))

            self._setup_daq()

            task = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], pt)
            self.daq_out.run(task)
            self.daq_out.waitToFinish(task)
            self.daq_out.stop(task)
            self.log('laser set to Vx={:.4}, Vy={:.4}'.format(float(x), float(y)))



        if self.settings['zero_first']:
            set_ao(0, 0)
            self.instruments['ANC300']['instance']._set_mode(1, 'input')
            self.instruments['ANC300']['instance']._set_mode(2, 'input')

        set_ao(self.settings['point']['x'], self.settings['point']['y'])
        self.instruments['ANC300']['instance']._set_mode(1, 'input')
        self.instruments['ANC300']['instance']._set_mode(2, 'input')

        #self.log_message()


class SetLaserInterferometer(SetLaser):
    """
    Simply SetLaser with a different scaling factor for the voltages. Wrote this to set IR laser spot using kinematic mount piezos, but can be much more general.
    """
    _INSTRUMENTS = {'NI9263_02': NI9263_02}


class SetLaserSingleAxis(SetLaser):
    """
    Sets the voltages sent by the DAQ analog outputs to the galvo mirrors. Typically used to move the galvo such that the laser spot is on a point of interest.
    This script sets the voltage of ONLY ONE AXIS. This script is meant as a subscript for something like GalvoGridScan.
    """

    _DEFAULT_SETTINGS = [
        Parameter('point', 0, float, 'voltage (V)'),
        Parameter('axis', 'x', ['x', 'y'], 'galvo axis'),
        Parameter('DAQ_channels',
            [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
            Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output')
            ]),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        pt = self.settings['point']/self.scale()
        #print('NOT rounding position to 2 dp')
        try:
            self.check_bounds(pt)
        except AttributeError:
            return

        # daq API only accepts either one point and one channel or multiple points and multiple channels
        pt = self.settings['point']
        self._setup_daq()
        channel_dict = {'x': self.settings['DAQ_channels']['x_ao_channel'], 'y': self.settings['DAQ_channels']['y_ao_channel']}
        task = self.daq_out.setup_AO([channel_dict[self.settings['axis']]], [pt])
        self.daq_out.run(task)
        self.daq_out.waitToFinish(task)
        self.daq_out.stop(task)
        self.log('laser set to V%s=%.4f' % (self.settings['axis'], pt))

    def plot(self, figure_list):
        pass

    def check_bounds(self, x):
        pass


class SetAttoSingleAxis(SetAtto):
    """
    Sets the voltages sent by the DAQ analog outputs to the Attocube controller.
    This script sets the voltage of ONLY ONE AXIS. This script is meant as a subscript for something like PulsedEsrLineScan.
    """

    _DEFAULT_SETTINGS = [
        Parameter('point', 0, float, 'voltage (V)'),
        Parameter('axis', 'x', ['x', 'y'], 'Attocube axis'),
        Parameter('zero_first', False, bool, 'Zero piezo inputs before moving to setpoint, to avoid hysteresis effects'),
        Parameter('DAQ_channels',
            [Parameter('x_ao_channel', 'ao2', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
            Parameter('y_ao_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output')
            ]),
        Parameter('daq_type', 'cDAQ', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        def set_ao(v):
            pt = v / self.scale()
            try:
                self.check_bounds(pt)
            except AttributeError:
                return

            # daq API only accepts either one point and one channel or multiple points and multiple channels
            self._setup_daq()
            channel_dict = {'x': self.settings['DAQ_channels']['x_ao_channel'], 'y': self.settings['DAQ_channels']['y_ao_channel']}
            task = self.daq_out.setup_AO([channel_dict[self.settings['axis']]], [pt])
            #task = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], [pt])
            self.daq_out.run(task)
            self.daq_out.waitToFinish(task)
            self.daq_out.stop(task)
            self.log('Attocube set to V%s=%.4f' % (self.settings['axis'], v))

        if self.settings['zero_first']:
            set_ao(0)
            self.instruments['ANC300']['instance']._set_mode(1, 'input')
            self.instruments['ANC300']['instance']._set_mode(2, 'input')

        set_ao(self.settings['point'])
        self.instruments['ANC300']['instance']._set_mode(1, 'input')
        self.instruments['ANC300']['instance']._set_mode(2, 'input')

    def check_bounds(self, v):
        if v < 0:
            self.log('Piezo cannot accept negative voltages!')
            raise AttributeError

    def plot(self, figure_list):
        pass


if __name__ == '__main__':
    from pylabcontrol.core import Instrument

    # instruments, instruments_failed = Instrument.load_and_append({'daq':  'NI6259'})

    script, failed, instruments = Script.load_and_append(script_dict={'SetLaser_cDAQ': 'SetLaser_cDAQ'})

    print(script)
    print(failed)
    # print(instruments)