import numpy as np

from b26_toolkit.src.instruments import NI6259, NI9263, NI9402
from b26_toolkit.src.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from PyLabControl.src.core import Script, Parameter
from galvo_scan_generic import GalvoScanGeneric


class GalvoScanDAQ(GalvoScanGeneric):
    """
    GalvoScan uses the apd, daq, and galvo to sweep across voltages while counting photons at each voltage,
    resulting in an image in the current field of view of the objective.
    """

    _DEFAULT_SETTINGS = [
        Parameter('DAQ_channels',
                   [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
                    Parameter('y_ao_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
                    Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1'], 'Daq channel used for counter')
                  ])
    ]

    _INSTRUMENTS = {'daq':  NI6259}

    def __init__(self, instruments, name=None, settings=None, log_function=None, data_path=None):
        '''
        Initializes GalvoScan script for use in gui

        Args:
            instruments: list of instrument objects
            name: name to give to instantiated script object
            settings: dictionary of new settings to pass in to override defaults
            log_function: log function passed from the gui to direct log calls to the gui log
            data_path: path to save data

        '''
        self._DEFAULT_SETTINGS = self._DEFAULT_SETTINGS + GalvoScanGeneric._DEFAULT_SETTINGS
        Script.__init__(self, name, settings=settings, instruments=instruments, log_function=log_function, data_path = data_path)

    def setup_scan(self):

        self.clockAdjust = int((self.settings['time_per_pt'] + self.settings['settle_time']) / self.settings['settle_time'])

        [xVmin, xVmax, yVmax, yVmin] = self.pts_to_extent(self.settings['point_a'],self.settings['point_b'],self.settings['RoI_mode'])


        self.x_array = np.repeat(np.linspace(xVmin, xVmax, self.settings['num_points']['x'], endpoint=True),self.clockAdjust)
        self.y_array = np.linspace(yVmin, yVmax, self.settings['num_points']['y'], endpoint=True)
        sample_rate = float(1) / self.settings['settle_time']
        self.instruments['daq']['instance'].settings['analog_output'][
            self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
        self.instruments['daq']['instance'].settings['analog_output'][
            self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate
        self.instruments['daq']['instance'].settings['digital_input'][
            self.settings['DAQ_channels']['counter_channel']]['sample_rate'] = sample_rate

    def get_galvo_location(self):
        """
        returns the current position of the galvo
        Returns: list with two floats, which give the x and y position of the galvo mirror
        """
        galvo_position = self.instruments['daq']['instance'].get_analog_voltages([
            self.settings['DAQ_channels']['x_ao_channel'],
            self.settings['DAQ_channels']['y_ao_channel']]
        )
        return galvo_position

    def set_galvo_location(self, galvo_position):
        """
        sets the current position of the galvo
        galvo_position: list with two floats, which give the x and y position of the galvo mirror
        """
        if galvo_position[0] > 1 or galvo_position[0] < -1 or galvo_position[1] > 1 or galvo_position[1] < -1:
            raise ValueError('The script attempted to set the galvo position to an illegal position outside of +- 1 V')

        pt = galvo_position
        daq = self.instruments['daq']['instance']
        # daq API only accepts either one point and one channel or multiple points and multiple channels
        pt = np.transpose(np.column_stack((pt[0],pt[1])))
        pt = (np.repeat(pt, 2, axis=1))

        daq.setup_AO([self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], pt)
        daq.AO_run()
        daq.AO_waitToFinish()
        daq.AO_stop()
