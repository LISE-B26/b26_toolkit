# # from https://stackoverflow.com/questions/40563139/ueye-camera-with-python-on-windows
from pylabcontrol.core import Parameter, Instrument
import cv2
import numpy as np


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
        # the following assumes only one camera is connected to the computer, since it will connect to camera number 0.
        self.cam = cv2.VideoCapture(0)
        assert self.cam.isOpened(), "Could not open camera!"
        self.cam.set(cv2.CAP_PROP_FRAME_HEIGHT, self.settings['height'])
        self.cam.set(cv2.CAP_PROP_FRAME_WIDTH, self.settings['width'])

    def update(self, settings):
        """
        Updates internal settings, as well as the pixel width and height on the physical device
        Args:
            settings: A dictionary in the form of settings as seen in default settings
        """
        super(UEyeCamera, self).update(settings)
        for key, value in settings.iteritems():
            if key == 'width':
                self.cam.set(cv2.CAP_PROP_FRAME_WIDTH, self.settings['width'])
            elif key == 'height':
                self.cam.set(cv2.CAP_PROP_FRAME_HEIGHT, self.settings['height'])

    def read_probes(self, key):
        pass

    def get_frame(self):
        """
        Reads and returns a single frame from the camera.

        Returns:
            A 2d numpy array containing the image
        """

        is_successful, bgr_image = self.cam.read()

        # try again if we get a completely black image
        if np.count_nonzero(bgr_image) == 0:
            is_successful, bgr_image = self.cam.read()

        if is_successful:
            return cv2.cvtColor(bgr_image, cv2.COLOR_BGR2GRAY)
        else:
            raise EnvironmentError("Could not successfully take image from camera.")


if __name__ == '__main__':
    camera = UEyeCamera()
    frame = camera.get_frame()
    cv2.imshow('Captured Frame', frame)
    cv2.waitKey(0)
    cv2.destroyAllWindows()