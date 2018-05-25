# # from https://stackoverflow.com/questions/40563139/ueye-camera-with-python-on-windows
from pylabcontrol.src.core import Parameter, Instrument
import cv2


class UEyeCamera(Instrument):
    """
    This class implements a UEye-compatable camera.
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('width', 800, int, 'width in pixels of image'),
        Parameter('height', 600, int, 'height in pixels of image')
    ])

    _PROBES = {}

    def __init__(self, name=None, settings=None):
        super(UEyeCamera, self).__init__(name, settings)
        self.cam = cv2.VideoCapture(0)
        assert self.cam.isOpened(), "Could not open camera!"
        self.cam.set(cv2.cv.CV_CAP_PROP_FRAME_HEIGHT, self.settings['height'])
        self.cam.set(cv2.cv.CV_CAP_PROP_FRAME_WIDTH, self.settings['width'])

    def update(self, settings):
        """
        Updates internal settings, as well as the pixel width and height on the physical device
        Args:
            settings: A dictionary in the form of settings as seen in default settings
        """
        super(UEyeCamera, self).update(settings)
        for key, value in settings.iteritems():
            if key == 'width':
                self.cam.set(cv2.CV_CAP_PROP_FRAME_WIDTH, self.settings['width'])
            elif key == 'height':
                self.cam.set(cv2.CV_CAP_PROP_FRAME_HEIGHT, self.settings['height'])

    def read_probes(self, key):
        pass

    def get_frame(self):
        """
        Reads and returns a single frame from the camera.

        Returns:
            A 2d numpy array containing the image
        """
        ret, buff = self.cam.read()

        if ret: #if read operation is successful
            return cv2.cvtColor(buff, cv2.COLOR_BGR2GRAY)
        else:
            raise EnvironmentError("Could not successfully take image from camera.")

if __name__ == '__main__':
    a = UEyeCamera()