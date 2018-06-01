# clr is python for .net
import clr # run pip install pythonnet
import sys, os
import time
from pylabcontrol.core.read_write_functions import get_config_value
from pylabcontrol.core import Parameter, Instrument


# import dll for SMC100


dll_path = get_config_value('SMC100_DLL_PATH',
                            os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.txt'))

print(dll_path)

if dll_path:

    sys.path.insert(0, dll_path)
    # Uses python for .net to add dll assembly to namespace
    try: 
        clr.AddReference('Newport.SMC100.CommandInterface')
        import CommandInterfaceSMC100
    except Exception as exception_details:
        print('Could not load SMC100 dll from path specified in configuration file. Check path is correct.')
        print(('exception encountered: ' + str(exception_details)))
        CommandInterfaceSMC100 = None
else:
    print("Could not import SMC100CommandInterface, will not be able to initialize SMC100 object.")
    CommandInterfaceSMC100 = None


########################################################################################################################
## INSTRUCTIONS ON USING THIS DLL AND PYTHON FOR .NET
# Most of the functions in this DLL take three arguments. The first is the controller number (hardcoded as 1 for now),
# the second is either a number to input or a by reference number to hold output, and the third is a string passed by
# reference to record errors. The functions also a return an int, 0 for success and -1 for error. Since passing by
# reference doesn't exist in python, python for .net will return all by reference values as part of an array. For
# example, the function to get position with one return value and two by reference values returns
# [int result, float position, str error]. However, YOU MUST GIVE DUMMY REFERENCE ARGUMENTS. The function to get
# position takes (int, float byref, str byref), and all three must be provided even if the latter two are just constants
# that get thrown out. This is contrary to the example in the provided documentation, but if you only provide (int),
# python won't be able to match arguments and resolve the function.
########################################################################################################################
#status codes
# https://www.newport.com/medias/sys_master/images/images/h8d/h3a/8797263101982/SMC100CC-SMC100PP-User-Manual.pdf
# see page 65
MOVING = '28'
DONE_MOVING = '33'

class SMC100(Instrument):
    """
Class to control the Newport SMC100 stepper motor driver. Class controlled over USB via DLL.
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('port', 'COM3', str, 'serial number written on device'),
        Parameter('position', 25000.0, float, 'servo position (from 0 to 25000 in um)'),
        Parameter('velocity', 1000, float, 'servo velocity (from 0 to 1000 in um/s)'),
        Parameter('height_lower_limit', 9, float, 'lowest position servo can move to (in mm)')
    ])

    def __init__(self, name = None, settings = None):
        """
        Initializes connection to the device.
        Args:
            name: device name
            settings: dictionary containing desired settings for instrument
        """
        if CommandInterfaceSMC100 is None:
            print("SMC100 dll not imported --- cannot initialize SMC100 object without dll")
            raise

        super(SMC100, self).__init__(name, settings)

        self.SMC = CommandInterfaceSMC100.SMC100()
        result = self.SMC.OpenInstrument(self.settings['port'])
        if result == -1:
            print(('Failed to open device on port' + str(self.settings['port'])))
            raise

    def update(self, settings):
        """
        Updates internal settings, as well as the position set on the physical device
        Args:
            settings: A dictionary in the form of settings as seen in default settings
        """
        super(SMC100, self).update(settings)
        for key, value in settings.items():
            if key == 'position':
                self._set_position(value/1e3)
            elif key == 'velocity':
                self._set_velocity(value/1e3)

    @property
    def _PROBES(self):
        return{
            'position': 'motor position in um',
            'velocity': 'motor velocity in um/s'
        }

    def read_probes(self, key):
        """
        requestes value from the instrument and returns it
        Args:
            key: name of requested value

        Returns: reads values from instrument

        """
        assert key in list(self._PROBES.keys())
        assert isinstance(key, str)

        if key == 'position':
            return self._get_position()*1e3
        elif key == 'velocity':
            return self._get_velocity()*1e3

    def __del__(self):
        """
        Cleans up SMC100 connection
        :PostState: SMC100 is disconnected
        """
        self.SMC.CloseInstrument()

    def _get_position(self):
        """
        Reads current position of motor.

        Returns: current position of motor in mm

        """
        i_ref = -1
        s_ref = ''
        result, response, errString = self.SMC.TP(1, i_ref, s_ref)
        if result == -1:
            print(('ERROR: ' + errString))
            raise
        return response

    def _set_position(self, pos):
        """
        Sets motor to desired position.

        Args:
            pos: position to set in mm. Must be less than height_lower_limit in settings.

        """
        if pos < self.settings['height_lower_limit']:
            raise ValueError('cannot set position below height_lower_limit')

        s_ref = ''
        # reenable computer control if keypad was used last
        self._enable_computer_control() #reenable computer control if keypad was used last
        #start movement
        result, errString = self.SMC.PA_Set(1, pos, s_ref)
        if result == -1:
            print(('ERROR: ' + errString))
            raise
        #block until movement is done
        while True:
            result, ErrorCode, StatusCode, errString = self.SMC.TS(1, s_ref, s_ref, s_ref)
            if result == -1:
                print(('ERROR: ' + errString))
                raise
            if StatusCode in [DONE_MOVING, '34']:
                break
            time.sleep(.1)
            # print('motor is moving - gui blocked StatusCode: {:s}'.format(str(StatusCode)))
    def _get_velocity(self):
        i_ref = -1
        s_ref = ''
        result, response, errString = self.SMC.VA_Get(1, i_ref, s_ref)
        if result == -1:
            print(('ERROR: ' + errString))
            raise
        return response

    def _set_velocity(self, vel):
        s_ref = ''
        result, errString = self.SMC.VA_Set(1, vel, s_ref)
        if result == -1:
            print(('ERROR: ' + errString))
            raise

    def _enable_computer_control(self):
        s_ref = ''
        result, errString = self.SMC.MM_Set(1, 1, s_ref)
        if result == -1:
            print(('ERROR: ' + errString))
            raise

if __name__ == "__main__":
    a = SMC100()
# print(a.position)
# a._enable_computer_control()
# a._set_position(15)
# a._set_position(25)
# a._set_position(10)
# print(a._get_position())
# inst, failed = Instrument.load_and_append({'SMC100': SMC100})
# print(inst)
# print(failed)szZzz y
