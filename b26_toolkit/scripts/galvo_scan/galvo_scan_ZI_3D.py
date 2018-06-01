from pylabcontrol.core import Script, Parameter
from .galvo_scan_ZI import GalvoScanZI
from b26_toolkit.instruments.newport_smc100 import SMC100
import numpy as np

class GalvoScanZI3D(Script):
    _DEFAULT_SETTINGS = [
        Parameter('z_axis_center_position', 6000, float, 'center point of autofocus sweep'),
        Parameter('scan_width', 5, float, 'distance (in V or mm) between the minimum and maximum points of the range'),
        Parameter('num_sweep_points', 10, int, 'number of values to sweep between min and max voltage'),
    ]

    _SCRIPTS = {'galvoscanZI': GalvoScanZI}
    _INSTRUMENTS = {'z_driver': SMC100}

    def __init__(self, scripts, instruments = None, name = None, settings = None, log_function = None, data_path = None):
        """

        """
        Script.__init__(self, name, settings, instruments, scripts, log_function= log_function, data_path = data_path)

    def _function(self):
        min_pos = self.settings['z_axis_center_position'] - self.settings['scan_width']/2.0
        max_pos = self.settings['z_axis_center_position'] + self.settings['scan_width']/2.0

        z_pos_list = np.linspace(min_pos, max_pos, self.settings['num_sweep_points'])
        self.data['sweep_voltages'] = z_pos_list
        for z_pos in z_pos_list:
            print(('zpos', z_pos))
            self._step_piezo(z_pos)
            self.scripts['galvoscanZI'].run()


    def _step_piezo(self, position):
        """
        steps the piezo.  Has to be overwritten specifically for each different hardware realization
        voltage: target piezo voltage
        wait_time: settle time after voltage step
        """
        z_driver = self.instruments['z_driver']['instance']
        # set the voltage on the piezo
        try:
            z_driver.position = float(position)
        except ValueError:
            raise
            self.log('requested value not permitted. Did not set value to {:0.3f} um'.format(position))

    def plot(self, figure_list):
        self.scripts['galvoscanZI'].plot(figure_list)