from pylabcontrol.core import Parameter, Instrument
import ctypes
import time

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
        Parameter('position', 3.0, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 1.0, float, 'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
    ])

    _PROBES = {}

    def __init__(self, name=None, settings=None):

        super(KDC001, self).__init__(name, settings)
        self._servo_library = ctypes.cdll.LoadLibrary('C:\\Program Files\\Thorlabs\\Kinesis\\Thorlabs.MotionControl.KCube.DCServo.dll')

        # conversion from mm to device units at
        # https://www.thorlabs.com/software/apt/APT_Communications_Protocol_Rev_15.pdf, page 16
        self._position_encoder2mm_conversion_factor = 34304
        self._velocity_encoder2mm_conversion_factor = 767367.49
        self._acceleration_encoder2mm_conversion_factor = 261.93
        self._acceleration = 1 * self._acceleration_encoder2mm_conversion_factor # use hard-coded acceleration of 1 mm/s^2 (can be changed to parameter if we want)
        self._manually_set_library_inputs_and_outputs()
        self._connect()

    def _manually_set_library_inputs_and_outputs(self):
        """
        Sets the input and output types for each servo library call we make.

        """
        self._servo_library.TLI_BuildDeviceList.restypes = ctypes.c_short

        self._servo_library.TLI_GetDeviceListSize.restypes = ctypes.c_short

        self._servo_library.TLI_GetDeviceListByTypeExt.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_ulong,
                                                                   ctypes.c_int]
        self._servo_library.TLI_GetDeviceListByTypeExt.restypes = ctypes.c_short

        self._servo_library.TLI_GetDeviceInfo.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(TLI_DeviceInfo)]
        self._servo_library.TLI_GetDeviceInfo.restypes = ctypes.c_short

        self._servo_library.CC_StartPolling.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int]
        self._servo_library.CC_StartPolling.restypes = ctypes.c_bool

        self._servo_library.CC_Open.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_Open.restypes = ctypes.c_short

        self._servo_library.CC_Close.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_Close.restypes = ctypes.c_short

        self._servo_library.CC_ClearMessageQueue.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_ClearMessageQueue.restypes = ctypes.c_short

        self._servo_library.CC_WaitForMessage.argtypes = [ctypes.POINTER(ctypes.c_char),
                                                          ctypes.POINTER(ctypes.c_ushort),
                                                          ctypes.POINTER(ctypes.c_ushort),
                                                          ctypes.POINTER(ctypes.c_ulong)]
        self._servo_library.CC_WaitForMessage.restypes = ctypes.c_bool

        self._servo_library.CC_GetVelParams.argtypes = [ctypes.POINTER(ctypes.c_char),
                                                        ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]

        self._servo_library.CC_GetVelParams.restypes = ctypes.c_short

        self._servo_library.CC_MoveToPosition.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int]
        self._servo_library.CC_MoveToPosition.restypes = ctypes.c_short

        self._servo_library.CC_SetVelParams.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int, ctypes.c_int]
        self._servo_library.CC_SetVelParams.restypes = ctypes.c_short

        self._servo_library.CC_StopPolling.argtypes = [ctypes.POINTER(ctypes.c_char)]

        self._servo_library.CC_Home.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_Home.restypes = ctypes.c_short

        self._servo_library.CC_GetPosition.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_GetPosition.restypes = ctypes.c_int

    def _connect(self, verbose=True):

        """
        makes connection to the device through the serial port. Also opens connection to the device.
        """

        # this version of the serial number is useful
        self._serial_num = ctypes.c_char_p(bytes(str(self.settings['serial_number']), "utf-8"))

        if self._servo_library.TLI_BuildDeviceList() == 0:
            num_devices = self._servo_library.TLI_GetDeviceListSize()
            if verbose:
                print('Number of devices detected: {0}'.format(num_devices))

            # The servo library fills in a byte string of connected device serial numbers, separated by a comma.
            # The length of this string will be the length of the serial number (8 bytes), plus a byte for a comma,
            # for each connected device.
            # Here, we first pre-allocate the string, then send it to the library function, and then examine it.

            connected_devices_string_length = num_devices*9
            string_of_serial_nums = ctypes.create_string_buffer(connected_devices_string_length)
            self._servo_library.TLI_GetDeviceListByTypeExt(string_of_serial_nums,
                                                           ctypes.c_ulong(connected_devices_string_length+1),
                                                           ctypes.c_int(27))

            list_of_connected_serial_nums = string_of_serial_nums.raw.decode("utf-8")

            if str(self.settings['serial_number']) not in list_of_connected_serial_nums:
                error_msg = 'No servo with the given serial number was detected.\n'
                error_msg += 'Given Serial Number: {0}\n'.format(self.settings['serial_number'])
                if list_of_connected_serial_nums:
                    error_msg += 'Connected Devices: {0}\n'.format(list_of_connected_serial_nums)
                else:
                    error_msg += 'Connected Devices: None\n'
                raise AttributeError(error_msg)

            elif verbose:
                print('Found device with matching serial number.')
                device_info = TLI_DeviceInfo()
                self._servo_library.TLI_GetDeviceInfo(self._serial_num, ctypes.byref(device_info))
                print("Description: ", device_info.description)
                print("Serial No: ", device_info.serialNo)
                print("USB PID: ", device_info.PID)

        if self.settings['velocity'] > 0.0: # update the velocity if not set to default
            self.set_velocity()

        # open connection to the device. The other option is to NOT call _open_device, and require the script to handle device opening and closing
        self._open_device()
        self.position = self.get_position()
        # todo(emma): implement self.get_velocity()
        #self.velocity = self.get_velocity()

    def set_position(self, verbose=True):

        """
        sets position of device in mm
        """
        print('executing set_position')
        position = ctypes.c_int(int(self._position_encoder2mm_conversion_factor * self.settings['position']))
        self._servo_library.CC_MoveToPosition(self._serial_num, position) # command sent to device
        message_type = ctypes.c_ushort(0)
        message_id = ctypes.c_ushort(0)
        message_data = ctypes.c_ulong(0)
        if verbose:
            print('Device is moving to indicated position ({0} mm)'.format(self.settings['position']))

        # wait for it to actually move
        self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
        while message_type.value != 2 or message_id.value != 1:
            self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)

        if verbose:
            position = self.get_position()
            print('Device now at position {0} mm'.format(position))

    def home(self, verbose=True):

        """
        Calibrates origin: scans full range and sets stage position to zero
        """

        self._servo_library.CC_Home(self._serial_num) # command sent to device!
        if verbose:
            print('Device is homing.')
        # wait for it to actually home
        message_type = ctypes.c_ushort(0)
        message_id = ctypes.c_ushort(0)
        message_data = ctypes.c_ulong(0)
        self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
        while message_type.value != 2 or message_id.value != 0:
            self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)

        if verbose:
            print('Device finished homing')

    def _open_device(self, verbose=True):

        """
        Opens the serial connection to the device, which lets you send commands to and receive data from the device.
        Also starts polling all of the devices connected to the computer.
        """

        if not self._servo_library.CC_Open(
                self._serial_num):
            if verbose:
                print('Starting to poll')
            milliseconds = ctypes.c_int(200)
            self._servo_library.CC_StartPolling(self._serial_num,
                                                milliseconds)  # continuous polling allows you to keep track of what the instrument is doing in real time
            time.sleep(3)
            self._servo_library.CC_ClearMessageQueue(self._serial_num)

    def _close_device(self, verbose=True):

        """
        Stops polling and closes the serial connection.
        """

        self._servo_library.CC_StopPolling(self._serial_num)
        self._servo_library.CC_Close(self._serial_num)
        if verbose:
            print('Device closed')

    def get_position(self, verbose=True):

        """
        returns position of stage in mm
        """

        position = self._servo_library.CC_GetPosition(
            self._serial_num) / self._position_encoder2mm_conversion_factor
        if verbose:
            print('Position of device is currently {0} mm'.format(position))
        return position

    def get_velocity(self, verbose=True):

        """
        returns velocity (when moving) of stage in mm/s
        """
       # raise NotImplementedError

        ## todo: fix the input arguments of get_velocity ER 20180519
        velocity_pointer_type = ctypes.POINTER(ctypes.c_long) #ctypes.POINTER(ctypes.c_long)#ctypes.LP_c_long()#ctypes.POINTER(ctypes.c_int)
        acceleration_pointer_type = ctypes.POINTER(ctypes.c_long) #ctypes.POINTER(ctypes.c_long)#ctypes.LP_c_long()#ctypes.POINTER(ctypes.c_int)
        velocity_pointer = ctypes.byref(ctypes.c_long())#velocity_pointer_type()
        acceleration_pointer = ctypes.byref(ctypes.c_long())#acceleration_pointer_type()
#        print(len(acceleration_pointer))
        print(acceleration_pointer)
        vel = self._servo_library.CC_GetVelParams(self.serial_num, velocity_pointer, acceleration_pointer)
        if verbose:
            print('Velocity of device is currently {0} mm/s'.format(vel))
        return vel

    def set_velocity(self, verbose=True):

        """
        sets velocity (when moving) of stage in mm/s
        """

        print('executing set_velocity')

        velocity = ctypes.c_int(int(self._velocity_encoder2mm_conversion_factor * self.settings['velocity']))
        acceleration = ctypes.c_int(int(self._acceleration))
        self._servo_library.CC_SetVelParams(self._serial_num, acceleration, velocity)
        if verbose:
            pass
          #  set_vel = self.get_velocity()
          # todo(emma): implement get_velocity()
           # print('Device velocity was set to {0} mm/s'.format(set_vel))

    def update(self, settings, verbose=True):
        print('update called!')
        super(KDC001, self).update(settings)
        print('got to here after super')
        if self._settings_initialized:
            for key, value in settings.items():
                print('updating keys')
                print(key, value)
                if key == 'position':
                    print('key is position')
                    self.set_position()
                elif key == 'velocity':
                    self.set_velocity()

    def read_probes(self, key = None):
        pass

    def __del__(self):
        # when done with device we have to close the connections
        # called when GUI is closed
        self._close_device()

class B26KDC001(KDC001):
    '''

    Same as KDC001 except one more parameters: axis , which will be associated with each serial number

    '''
    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 27501971, int, 'serial number written on device'),
        Parameter('position', 6.0, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 0.0, float,
                  'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used'),
        Parameter('axis', 'x', ['x', 'y', 'z'], 'axis associated with the serial number')
    ])

    _PROBES = {}

    def _set_serial_and_safety_lims(self):
        '''

        Hardcoded safety bounds and serial number values for our setup

        '''
        print('self.settings:')
        print(self.settings)
        if self.settings['axis'] == 'x':
            self.serial_number = 27501971
            self.min_pos = 0.
            self.max_pos = 1.
        elif self.settings['axis'] == 'y':
            self.serial_number = 27501971
            self.min_pos = 0.
            self.max_pos = 2.
        elif self.settings['axis'] == 'z':
            self.serial_number = 27501971
            self.min_pos = 0.
            self.max_pos = 3.
        else:
            raise AttributeError('no such axis!!')

    def update(self, settings):
        # call the update_parameter_list to update the parameter list
        print('self.settings before _set_serial_and_safety_limits:')
        print(self.settings)

        super(B26KDC001, self).update(settings)
        self._set_serial_and_safety_lims()
        print('self.settings after _set_serial_and_safety_limits:')
        print(self.settings)

    def read_probes(self, key):
        pass

    def set_position(self, verbose=True):
        # check if safety limits are met
        print('executing set_position in lower class')
        if self.settings['position'] < self.settings['max_pos'] and self.settings['position'] > self.settings['min_pos']:
            print('calling super set position')
            super(B26KDC001, self).set_position(self, verbose)
        else:
            print('didnt make the cut!')
            raise AttributeError('position is outside safety limits!! Doing nothing.')

if __name__ == '__main__':
    a = B26KDC001()
  #  a.set_velocity()
  #  a.set_position()