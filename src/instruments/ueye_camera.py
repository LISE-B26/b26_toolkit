# # from https://stackoverflow.com/questions/40563139/ueye-camera-with-python-on-windows
#
from PyLabControl.src.core import Parameter, Instrument

import sys
sys.path.append('C:\ProgramData\Anaconda2\Lib\site-packages')

try:
    import cv2
except:
    pass

import atexit

import matplotlib.pyplot as plt


class UEyeCamera(Instrument):
    """
    This class implements a UEye-compatable camera.
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('width', 800, int, 'width in pixels of image'),
        Parameter('height', 600, int, 'height in pixels of image')
    ])

    def __init__(self, name=None, settings=None):

        super(UEyeCamera, self).__init__(name, settings)
        self.cam = cv2.VideoCapture(0)

        # JG commented out the following two lines because hr and wr are not use anywhere!!
        # hr = self.cam.set(cv2.cv.CV_CAP_PROP_FRAME_HEIGHT, self.settings['height'])
        # wr = self.cam.set(cv2.cv.CV_CAP_PROP_FRAME_WIDTH, self.settings['width'])
        atexit.register(self.disconnect) #makes sure that camera is released when class exits

    def update(self, settings):
        '''
        Updates internal settings, as well as the pixel width and height on the physical device
        Args:
            settings: A dictionary in the form of settings as seen in default settings
        '''
        super(UEyeCamera, self).update(settings)
        # for key, value in settings.iteritems():
        #     if key == 'width':
        #         self.cam.set(cv2.CV_CAP_PROP_FRAME_WIDTH, self.settings['width'])
        #     elif key == 'height':
        #         self.cam.set(cv2.CV_CAP_PROP_FRAME_HEIGHT, self.settings['height'])

    # is often somehow getting called without init, or failing to get called at the proper time. Use atexit.register
    # instead
    # def __del__(self):
    #     try:
    #         self.cam.release()
    #     except AttributeError:
    #         #catches issue where del runs without init, so self.cam is never created
    #         pass


    def disconnect(self):
        """
        Releases camera from program control.
        """
        try:
            self.cam.release()
        except AttributeError:
            #catches issue where del runs without init, so self.cam is never created
            pass

    def get_frame(self):
        """
        Reads and returns a single frame from the camera.
        :return: A 2d numpy array containing the image
        """
        ret, buff = self.cam.read()
        if ret: #if read operation is successful
            return cv2.cvtColor(buff, cv2.COLOR_BGR2GRAY)
        else:
            raise EnvironmentError

# a = UEyeCamera()
# frame = a.get_frame()
#
# plt.imshow(frame)
# plt.show()