from pylabcontrol.core import Script, Parameter
from b26_toolkit.instruments import UEyeCamera
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new


class TakeImage(Script):
    """
    This wraps the ueye_camera get_frame method in a Script. This can be used to take images from the gui, or used like
    GalvoScan in Autofocus or similar Scripts.
    """

    _DEFAULT_SETTINGS = [
        Parameter('width', 800, int, 'width in pixels of image'),
        Parameter('height', 600, int, 'height in pixels of image')
    ]

    _INSTRUMENTS = {'Camera': UEyeCamera}

    _SCRIPTS = {}


    def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Sets up script and sends width and height settings to camera
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_function= log_function, data_path = data_path)
        self.instruments['Camera']['instance'].update({'width': self.settings['width'], 'height': self.settings['height']})

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """
        self.data = {'image_data': self.instruments['Camera']['instance'].get_frame(), 'extent': [0,self.settings['width'], self.settings['height'], 0]}

    def _plot(self, axes_list, data = None):
        """
        Plots camera image to first axis
        :param axes_list: list containing axis to plot to as zeroth element
        :param data:
        """
        plot_fluorescence_new(self.data['image_data'], self.data['extent'], axes_list[0])

if __name__ == '__main__':
    script, failed, instruments = Script.load_and_append(script_dict={'TakeImage': 'TakeImage'})

    print(script)
    print(failed)