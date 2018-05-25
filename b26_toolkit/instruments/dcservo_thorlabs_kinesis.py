"""
    This file is part of b26_toolkit, a PyLabControl add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell
    PyLabControl is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    PyLabControl is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with PyLabControl.  If not, see <http://www.gnu.org/licenses/>.
"""
# clr is python for .net
import clr # run pip install pythonnet
import sys, os
from pylabcontrol.core.read_write_functions import get_config_value
# JG: July 27 2016 uncommented folowing line: don't use import *!
# from pylabcontrol.core.instruments import *
from pylabcontrol.core import Parameter, Instrument
import ctypes
import time

dll_path = get_config_value('KINESIS_DLL_PATH', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.txt'))
if dll_path:
    sys.path.insert(0, dll_path)
# makes each dll, corresponding to a namespace, avaliable to python at runtime
    try:
        clr.AddReference('ThorLabs.MotionControl.DeviceManagerCLI')
        clr.AddReference('Thorlabs.MotionControl.TCube.DCServoCLI')
        clr.AddReference('Thorlabs.MotionControl.KCube.DCServoCLI')
        clr.AddReference('System')
        # imports classes from the namespaces. All read as unresolved references because python doesn't know about the dlls
        # until runtime
        # adds .NET stuctures corresponding to primitives
        from System import Decimal, Double, String
        # names will be in red due to the fact that pycharm can't see the thorlabs library until runtime
        from Thorlabs.MotionControl.DeviceManagerCLI import DeviceManagerCLI
        from Thorlabs.MotionControl.TCube.DCServoCLI import TCubeDCServo
        from Thorlabs.MotionControl.KCube.DCServoCLI import KCubeDCServo
    except Exception as exception_details:
        print("Could not load Thorlabs dll's to control Thorlabs servos.")
        print("exception details " + str(exception_details))
        DeviceManagerCLI = None
        TCubeDCServo = None
        KCubeDCServo = None
else:
    print("Could not import Thorlabs dll's --- will not be able to initialize SMC100 object.")
    DeviceManagerCLI = None
    TCubeDCServo = None
    KCubeDCServo = None
class ThorlabsServo(Instrument):
    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 27501971, int, 'serial number written on device'),
        Parameter('position', 0, float, 'servo position (from 0 to 6 in mm)'),
        Parameter('velocity', 0, float, 'servo maximum velocity in mm/s')
    ])
    def __init__(self, name = None, settings = None):
        """
        Initializes and connects to motors. Must define a variable self.Servo that aliases the correct device, ex.
        self.Servo = KCubeDCServo, and must call self._connect at some point.
        Args:
            name:
            settings:
        """
        raise NotImplementedError
    def _connect(self):
        """
        Connects to the servo using parameters defined elsewhere in the code.
        """
        raise NotImplementedError
    def update(self, settings):
        '''
        Updates internal settings, as well as the position and velocity set on the physical device. If the instrument's
        serial number is changed, will disconnect from the current device (if any) and reconnect to the new device.
        Args:
            settings: A dictionary in the form of settings as seen in default settings
        '''
        super(ThorlabsServo, self).update(settings)
        #update will usually trigger in super.__init__, which is generally the first line of self.__init__. At this
        #point, the servo won't be connected, so don't make changes in hardware at this point. After the end of the
        #__init__ when initialization is complete and the servo is connected, we can then have changes update hardware
        if self._settings_initialized:
            for key, value in settings.iteritems():
                if key == 'position':
                    self._move_servo(value)
                elif key == 'velocity':
                    self._set_velocity(value)
                elif key == 'serial_number':
                    #if it exists, disconnect from previous device
                    try:
                        self.device.StopPolling()
                        self.device.Disconnect()
                    #either the old connection was severed or there was no old connection, so now (re)connect
                    finally:
                        self._connect()
    @property
    def _PROBES(self):
        return{
            'position': 'servo position in mm',
            'velocity': 'servo velocity in mm/s'
        }
    def read_probes(self, key):
        assert key in self._PROBES.keys()
        assert isinstance(key, str)
        #query always returns string, need to cast to proper return type
        if key in ['position']:
            return self._get_position()
        elif key in ['velocity']:
            return self._get_velocity()
    @property
    def is_connected(self):
        return self.connected
        # DeviceManagerCLI.BuildDeviceList()
        # return(str(self.settings['serial_number']) in DeviceManagerCLI.GetDeviceList(self.Servo.DevicePrefix))
    # def __del__(self):
    #     '''
    #     Cleans up TDC001 connection
    #     :PostState: TDC001 is disconnected
    #     '''
    #     self.device.StopPolling()
    #     self.device.Disconnect()
    def goto_home(self):
        '''
        Recenters device at the home position. Raises an exception on failure.
        '''
        try:
            self.device.Home(60000)
        except Exception:
            print("Failed to move to position")
            raise
    def _move_servo(self, position, velocity = 0):
        '''
        Move servo to given position with given maximum velocity. Raises an exception on failure.
        :param position: position in mm, ranges from 0-6
        :param velocity: maximum velocity in mm/s, ranges from 0-2.5
        :PostState: servo has moved
        '''
        try:
            if(velocity != 0):
                self._set_velocity(velocity)
            # print("Moving Device to " + str(position))
            print('initialized3', self.device.IsSettingsInitialized())
            self.device.MoveTo(self._Py_Decimal(position), 60000)
        except Exception:
            print("Failed to move to position")
            raise
    def _get_position(self):
        '''
        :return: position of servo
        '''
        if(not self.connected):
            return -1
        else:
            return self._Undo_Decimal(self.device.Position)
    def _set_velocity(self, velocity):
        '''
        :param maximum velocity in mm/s, ranges from 0-2.5
        :PostState: velocity changed in hardware
        '''
        if(velocity != 0):
            velPars = self.device.GetVelocityParams()
            velPars.MaxVelocity = self._Py_Decimal(velocity)
            self.device.SetVelocityParams(velPars)
    def _get_velocity(self):
        '''
        :return: maximum velocity setting
        '''
        if(not self.connected):
            return -1
        else:
            return self._Undo_Decimal(self.device.GetVelocityParams().MaxVelocity)
    def _Py_Decimal(self, value):
        '''
        Casting a python double to System.Decimal results in the Decimal having only integer values, likely due to an
        improper selection of the overloaded Decimal function. Casting it first to System.Double, which always maintains
        precision, then from Double to Decimal, where the proper overloaded function is clear, bypasses this issue
        :param value: a python double
        :return: the input as a System.Decimal
        '''
        return Decimal(Double(value))
    def _Undo_Decimal(self, value):
        '''
        Casting back from System.Decimal to a python float fails due to overloading issues, but one can successfully
        cast back to a string. Thus, we use a two-part cast to return to python numeric types
        :param value: a System.Decimal
        :return: the input as a python float
        '''
        return float(str(value))

class TDC001(ThorlabsServo):
    """
    Class to control the thorlabs TDC001 servo. Note that ALL DLL FUNCTIONS TAKING NUMERIC INPUT REQUIRE A SYSTEM.DECIMAL
    VALUE. Check help doc at C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.DotNet_API for the DLL api.
    The class communicates with the device over USB.
    """
    def __init__(self, name = None, settings = None):
        super(ThorlabsServo, self).__init__(name, settings)
        self.Servo = TCubeDCServo
        self._connect()
    def _connect(self):
        try:
            DeviceManagerCLI.BuildDeviceList()
            serial_number_list = DeviceManagerCLI.GetDeviceList(self.Servo.DevicePrefix)
        except (Exception):
            print("Exception raised by BuildDeviceList")
        if not (str(self.settings['serial_number']) in serial_number_list):
            print(str(self.settings['serial_number']) + " is not a valid serial number")
            raise
        self.device = self.Servo.CreateTCubeDCServo(str(self.settings['serial_number']))
        if self.device is None:
            error_msg = self.settings['serial_number'] + " is not a TCubeDCServo"
            raise AttributeError(error_msg)
        try:
            self.device.Connect(str(self.settings['serial_number']))
        except Exception:
            print('Failed to open device ' + str(self.settings['serial_number']))
            # raise
        if not self.device.IsSettingsInitialized():
            try:
                self.device.WaitForSettingsInitialized(5000)
            except Exception:
                print("Settings failed to initialize")
                raise
        self.device.StartPolling(250)
        motorSettings = self.device.GetMotorConfiguration(str(self.settings['serial_number']))
        currentDeviceSettings = self.device.MotorDeviceSettings

class TLI_DeviceInfo(ctypes.Structure):
    _fields_ = [("typeID", ctypes.c_ulong),
                ("description", (65 * ctypes.c_char)),
                ("serialNo", (9 * ctypes.c_char)),
                ("PID", ctypes.c_ulong),
                ("isKnownType", ctypes.c_bool),
                ("motorType", ctypes.c_int),
                ("isPiezoDevice", ctypes.c_bool),
                ("isLaser", ctypes.c_bool),
                ("isCustomType", ctypes.c_bool),
                ("isRack", ctypes.c_bool),
                ("maxChannels", ctypes.c_short)]

class KDC001(Instrument):
    """
    Class to control the thorlabs KDC001 servo. Note that ALL DLL FUNCTIONS TAKING NUMERIC INPUT REQUIRE A SYSTEM.DECIMAL
    VALUE. Check help doc at C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.DotNet_API for the DLL api.
    The class communicates with the device over USB.
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 27501971, int, 'serial number written on device'),
        Parameter('position', 25.0, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 0.0, float, 'servo maximum velocity (0 to 2.6) [mm/s]')
    ])

    def __init__(self, name=None, settings=None):
        super().__init__(name, settings)
        self.servo_library = ctypes.cdll.LoadLibrary('C:\\Program Files\\Thorlabs\\Kinesis\\Thorlabs.MotionControl.KCube.DCServo.dll')
        self.position_encoder2mm_conversion_factor = 34304
        #TODO: figure out what the conversion factor is for the velocity!
        # self.velocity_encoder2mm_conversion_factor =
        self.manually_set_library_inputs_and_outputs()
        self._connect()

    def manually_set_library_inputs_and_outputs(self):
        """
        Sets the input and output types for each servo library call we make.

        """
        self.servo_library.TLI_BuildDeviceList.restypes = ctypes.c_short

        self.servo_library.TLI_GetDeviceListSize.restypes = ctypes.c_short

        self.servo_library.TLI_GetDeviceListByTypeExt.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_ulong,
                                                                  ctypes.c_int]
        self.servo_library.TLI_GetDeviceListByTypeExt.restypes = ctypes.c_short

        self.servo_library.TLI_GetDeviceInfo.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(TLI_DeviceInfo)]
        self.servo_library.TLI_GetDeviceInfo.restypes = ctypes.c_short

        self.servo_library.CC_StartPolling.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int]
        self.servo_library.CC_StartPolling.restypes = ctypes.c_bool

        self.servo_library.CC_Open.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self.servo_library.CC_Open.restypes = ctypes.c_short

        self.servo_library.CC_Close.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self.servo_library.CC_Close.restypes = ctypes.c_short

        self.servo_library.CC_ClearMessageQueue.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self.servo_library.CC_ClearMessageQueue.restypes = ctypes.c_short

        self.servo_library.CC_WaitForMessage.argtypes = [ctypes.POINTER(ctypes.c_char),
                                                         ctypes.POINTER(ctypes.c_ushort),
                                                         ctypes.POINTER(ctypes.c_ushort),
                                                         ctypes.POINTER(ctypes.c_ulong)]
        self.servo_library.CC_WaitForMessage.restypes = ctypes.c_bool

        self.servo_library.CC_GetVelParams.argtypes = [ctypes.POINTER(ctypes.c_char),
                                                       ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        self.servo_library.CC_GetVelParams.restypes = ctypes.c_short

        self.servo_library.CC_MoveToPosition.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int]
        self.servo_library.CC_MoveToPosition.restypes = ctypes.c_short

        self.servo_library.CC_SetVelParams.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int, ctypes.c_int]
        self.servo_library.CC_SetVelParams.restypes = ctypes.c_short

        self.servo_library.CC_StopPolling.argtypes = [ctypes.POINTER(ctypes.c_char)]

        self.servo_library.CC_Home.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self.servo_library.CC_Home.restypes = ctypes.c_short

        self.servo_library.CC_GetPosition.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self.servo_library.CC_GetPosition.restypes = ctypes.c_int

    def _connect(self, verbose=True):

        # this version of the serial number is useful
        serial_num = ctypes.c_char_p(bytes(str(self.settings['serial_number']), "utf-8"))

        if self.servo_library.TLI_BuildDeviceList() == 0:
            num_devices = self.servo_library.TLI_GetDeviceListSize()
            if verbose:
                print('Number of devices detected: {0}'.format(num_devices))

            # The servo library fills in a byte string of connected device serial numbers, separated by a comma.
            # The length of this string will be the length of the serial number (8 bytes), plus a byte for a comma,
            # for each connected device. We subtract 1 since the last entry does not have a comma.
            # Here, we first pre-allocate the string, then send it to the library function, and then examine it.

            connected_devices_string_length = num_devices*9 -1
            string_of_serial_nums = ctypes.create_string_buffer(connected_devices_string_length)
            self.servo_library.TLI_GetDeviceListByTypeExt(string_of_serial_nums,
                                                          ctypes.c_ulong(connected_devices_string_length),
                                                          ctypes.c_int(27))

            list_of_connected_serial_nums = string_of_serial_nums.raw.decode("utf-8").split(',')
            if str(self.settings['serial_number']) not in list_of_connected_serial_nums:
                error_msg = 'No servo with the given serial number was detected.\n'
                error_msg += 'Connected Devices: {0}\n'.format(list_of_connected_serial_nums)
                error_msg += 'Given Serial Number: {0}'.format(self.settings['serial_number'])
                raise AttributeError(error_msg)

            elif verbose:
                print('Found device with matching serial number.')
                device_info = TLI_DeviceInfo()
                self.servo_library.TLI_GetDeviceInfo(serial_num, ctypes.byref(device_info))
                print("Description: ", device_info.description)
                print("Serial No: ", device_info.serialNo)
                print("USB PID: ", device_info.PID)

            if not self.servo_library.CC_Open(serial_num):
                if verbose:
                    print('Starting to poll')
                milliseconds = ctypes.c_int(200)
                self.servo_library.CC_StartPolling(serial_num, milliseconds)
                time.sleep(3)
                self.servo_library.CC_ClearMessageQueue(serial_num)
                self.servo_library.CC_Home(serial_num)
                if verbose:
                    print('Device is homing.')

                message_type = ctypes.c_ushort(0)
                message_id = ctypes.c_ushort(0)
                message_data = ctypes.c_ulong(0)
                self.servo_library.CC_WaitForMessage(serial_num, message_type, message_id, message_data)
                if verbose:
                    print('printing all messages from K Cube while we wait for finished homing indicator')
                while message_type.value != 2 or message_id.value != 0:
                    if verbose:
                        print('message_type: {0}'.format(message_type))
                        print('message id: {0}'.format(message_id))
                        print('message data: {0}\n'.format(message_data))
                    self.servo_library.CC_WaitForMessage(serial_num, message_type, message_id, message_data)

                if verbose:
                    print('message_type: {0}'.format(message_type))
                    print('message id: {0}'.format(message_id))
                    print('message data: {0}\n'.format(message_data))
                    print('Device finished homing.')

                # if self.settings['velocity'] > 0:
                #     self.servo_library.CC_SetVelParams(serial_num, ctypes.c)

                self.servo_library.CC_ClearMessageQueue(serial_num)
                position = ctypes.c_int(int(self.position_encoder2mm_conversion_factor * self.settings['position']))
                self.servo_library.CC_MoveToPosition(serial_num, position)
                if verbose:
                    print('Device is moving to indicated position ({0} mm).'.format(self.settings['position']))

                self.servo_library.CC_WaitForMessage(serial_num, message_type, message_id, message_data)
                while message_type.value != 2 or message_id.value != 1:
                    print('message_type: {0}'.format(message_type))
                    print('message id: {0}'.format(message_id))
                    print('message data: {0}\n'.format(message_data))
                    self.servo_library.CC_WaitForMessage(serial_num, message_type, message_id, message_data)

                if verbose:
                    position = self.servo_library.CC_GetPosition(serial_num)/self.position_encoder2mm_conversion_factor
                    print('Device now at position {0}'.format(position))

                self.servo_library.CC_StopPolling(serial_num)
                self.servo_library.CC_Close(serial_num)

if __name__ == '__main__':
    #A test function for the device. Tries to connect to the
    a = KDC001()
    # print(a.is_connected)
    # a._move_servo(5)