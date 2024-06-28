from pylabcontrol.core import Parameter, Instrument
import ctypes
import time
import clr

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
    VALUE. Check help doc at C:/Program Files/Thorlabs/Kinesis/Thorlabs.MotionControl.DotNet_API for the DLL api.
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
        import sys
        sys.path.insert(0, 'C:\\Program Files\\Thorlabs\\Kinesis\\')
        clr.AddReference('ThorLabs.MotionControl.KCube.DCServoCLI')
        self._servo_library = ctypes.cdll.LoadLibrary('C:\\Program Files\\Thorlabs\\Kinesis\\Thorlabs.MotionControl.KCube.DCServo.dll')

        #os.environ['PATH'] = '' + os.pathsep + os.environ['PATH']
       # self._servo_library = ctypes.cdll.LoadLibrary(r"C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.KCube.DCServo.dll")


        # conversion from mm to device units at
        # https://www.thorlabs.com/software/apt/APT_Communications_Protocol_Rev_15.pdf, page 16
        self._position_encoder2mm_conversion_factor = 34304
        self._velocity_encoder2mm_conversion_factor = 767367.49
        self._acceleration_encoder2mm_conversion_factor = 261.93
        self._angle_encoder2deg_conversion_factor = 1638.4 # 12288
        self._acceleration = 4 * self._acceleration_encoder2mm_conversion_factor # use hard-coded acceleration of 1 mm/s^2 (can be changed to parameter if we want)
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
        self._serial_num_int = self.settings['serial_number']
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

        if self.settings['velocity'] > 0.0:  # update the velocity if not set to default
            self.set_velocity()

        # open connection to the device. The other option is to NOT call _open_device, and require the script to handle device opening and closing
        self._open_device()

        # get the current position of the device and set the position setting to this value. That makes sure that the device does not move when you create the instrument.
        self.settings['position'] = self.get_position()
        # todo(emma): implement self.get_velocity()
        #self.velocity = self.get_velocity()

    def set_position(self, verbose=False):

        """
        sets position of device in mm
        """
        position = ctypes.c_int(int(self._position_encoder2mm_conversion_factor * self.settings['position']))
        self._servo_library.CC_MoveToPosition(self._serial_num, position)  # command sent to device
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
            position = self.get_position(verbose=True)
            #print('Device now at position {0} mm'.format(position))

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

    def get_position(self, verbose=False):

        """
        returns position of stage in mm
        """

        position = self._servo_library.CC_GetPosition(
            self._serial_num) / self._position_encoder2mm_conversion_factor
        if verbose:
            print('Position of device is currently {0} mm'.format(position))
            print('Position of device is currently %i d.u.'%(position*self._position_encoder2mm_conversion_factor))
            print('Position of device is currently %.3f deg' % (position * self._position_encoder2mm_conversion_factor/self._angle_encoder2deg_conversion_factor))
        return position

    def get_velocity(self):

        """
        returns velocity (when moving) of stage in mm/s
        """

        ## todo: fix the input arguments of get_velocity ER 20180519
        velocity_pointer = ctypes.pointer(ctypes.c_int())#velocity_pointer_type()
        acceleration_pointer = ctypes.pointer(ctypes.c_int())#acceleration_pointer_type()
        vel = self._servo_library.CC_GetVelParams(self._serial_num, acceleration_pointer, velocity_pointer)

        return velocity_pointer.contents.value/self._velocity_encoder2mm_conversion_factor, acceleration_pointer.contents.value/self._acceleration_encoder2mm_conversion_factor

    def get_trigger_config(self):

        """
        Configure trigger 1&2 type (e.g. trigger at max velocity or at certain steps), and set trigger polarity
        :return:
        """
        """
        Trigger mode:
        0: Trigger disabled
        1: Trigger Input - General purpose logic input
        2: Trigger Input - Move relative using relative move parameters
        3: Trigger Input - Move absolute using absolute move parameters
        4: Trigger Input - Perform a Home action
        10: Trigger Output - General purpose output
        11: Trigger Output - Set when device moving
        12: Trigger Output - Set when at max velocity
        13: Trigger Output - Set when at predefine position steps
        14: Trigger Output - TBD mode
        
        Trigger polarity:
        0: Low
        1: High
        """
        trigger1Mode_pointer = ctypes.pointer(ctypes.c_int16())
        trigger1Polarity_pointer = ctypes.pointer(ctypes.c_int16())
        trigger2Mode_pointer = ctypes.pointer(ctypes.c_int16())
        trigger2Polarity_pointer = ctypes.pointer(ctypes.c_int16())

        self._servo_library.CC_GetTriggerConfigParams(self._serial_num, trigger1Mode_pointer, trigger1Polarity_pointer,
                                                      trigger2Mode_pointer, trigger2Polarity_pointer)

        trigger1Mode = trigger1Mode_pointer.contents.value
        trigger1Polarity = trigger1Polarity_pointer.contents.value
        trigger2Mode = trigger2Mode_pointer.contents.value
        trigger2Polarity = trigger2Polarity_pointer.contents.value

        return trigger1Mode, trigger1Polarity, trigger2Mode, trigger2Polarity

    def set_trigger_config(self, trigger_config):

        """
        Configure trigger 1&2 type (e.g. trigger at max velocity or at certain steps), and set trigger polarity
        :return:
        """
        """
        Trigger mode:
        0: Trigger disabled
        1: Trigger Input - General purpose logic input
        2: Trigger Input - Move relative using relative move parameters
        3: Trigger Input - Move absolute using absolute move parameters
        4: Trigger Input - Perform a Home action
        10: Trigger Output - General purpose output
        11: Trigger Output - Set when device moving
        12: Trigger Output - Set when at max velocity
        13: Trigger Output - Set when at predefine position steps
        14: Trigger Output - TBD mode

        Trigger polarity:
        0: Low
        1: High
        """

        trigger1Mode_pointer = ctypes.c_int16(int(trigger_config['trigger1Mode']))
        trigger1Polarity_pointer = ctypes.c_int16(int(trigger_config['trigger1Polarity']))
        trigger2Mode_pointer = ctypes.c_int16(int(trigger_config['trigger2Mode']))
        trigger2Polarity_pointer = ctypes.c_int16(int(trigger_config['trigger2Polarity']))

        set_trigger_config = self._servo_library.CC_SetTriggerConfigParams(self._serial_num, trigger1Mode_pointer, trigger1Polarity_pointer,
                                                      trigger2Mode_pointer, trigger2Polarity_pointer)

        return set_trigger_config


    def get_trigger_params(self, verbose=True):

        """
        Returns trigger parameters
        """

        triggerStartPositionFwd_pointer = ctypes.pointer(ctypes.c_int32())
        triggerIntervalFwd_pointer = ctypes.pointer(ctypes.c_int32())
        triggerPulseCountFwd_pointer = ctypes.pointer(ctypes.c_int32())
        triggerStartPositionRev_pointer = ctypes.pointer(ctypes.c_int32())
        triggerIntervalRev_pointer = ctypes.pointer(ctypes.c_int32())
        triggerPulseCountRev_pointer = ctypes.pointer(ctypes.c_int32())
        triggerPulseWidth_pointer = ctypes.pointer(ctypes.c_int32())
        cycleCount_pointer = ctypes.pointer(ctypes.c_int32())


        # Returns 0 if no error. Variable not used
        get_trigger_params = self._servo_library.CC_GetTriggerParamsParams(self._serial_num, triggerStartPositionFwd_pointer, triggerIntervalFwd_pointer,
                                                                       triggerPulseCountFwd_pointer, triggerStartPositionRev_pointer,
                                                                       triggerIntervalRev_pointer, triggerPulseCountRev_pointer, triggerPulseWidth_pointer,
                                                                       cycleCount_pointer)

        # All distances in mm, all time values in us
        triggerStartPositionFwd = triggerStartPositionFwd_pointer.contents.value / self._position_encoder2mm_conversion_factor
        triggerIntervalFwd = triggerIntervalFwd_pointer.contents.value / self._position_encoder2mm_conversion_factor
        triggerPulseCountFwd = triggerPulseCountFwd_pointer.contents.value
        triggerStartPositionRev = triggerStartPositionRev_pointer.contents.value / self._position_encoder2mm_conversion_factor
        triggerIntervalRev = triggerIntervalRev_pointer.contents.value / self._position_encoder2mm_conversion_factor
        triggerPulseCountRev = triggerPulseCountRev_pointer.contents.value
        triggerPulseWidth = triggerPulseWidth_pointer.contents.value
        cycleCount = cycleCount_pointer.contents.value

        return triggerStartPositionFwd, triggerIntervalFwd, triggerPulseCountFwd, triggerStartPositionRev, triggerIntervalRev, triggerPulseCountRev, \
               triggerPulseWidth, cycleCount

    def set_trigger_params(self, trigger_params, verbose=True):
        """
        :param trigger_params: Dictionary containing trigger parameters
        :param verbose:
        :return:
        """

        triggerStartPositionFwd_pointer = ctypes.c_int32(int(trigger_params['triggerStartPositionFwd'] * self._position_encoder2mm_conversion_factor))
        triggerIntervalFwd_pointer = ctypes.c_int32(int(trigger_params['triggerIntervalFwd'] * self._position_encoder2mm_conversion_factor))
        triggerPulseCountFwd_pointer = ctypes.c_int32(int(trigger_params['triggerPulseCountFwd']))
        triggerStartPositionRev_pointer = ctypes.c_int32(int(trigger_params['triggerStartPositionRev'] * self._position_encoder2mm_conversion_factor))
        triggerIntervalRev_pointer = ctypes.c_int32(int(trigger_params['triggerIntervalRev'] * self._position_encoder2mm_conversion_factor))
        triggerPulseCountRev_pointer = ctypes.c_int32(int(trigger_params['triggerPulseCountRev']))
        triggerPulseWidth_pointer = ctypes.c_int32(int(trigger_params['triggerPulseWidth']))
        cycleCount_pointer = ctypes.c_int32(int(trigger_params['cycleCount']))

        set_trigger_params = self._servo_library.CC_SetTriggerParamsParams(self._serial_num, triggerStartPositionFwd_pointer, triggerIntervalFwd_pointer,
                                                                       triggerPulseCountFwd_pointer, triggerStartPositionRev_pointer,
                                                                       triggerIntervalRev_pointer, triggerPulseCountRev_pointer, triggerPulseWidth_pointer,
                                                                       cycleCount_pointer)
        return set_trigger_params

    def set_velocity(self):

        """
        sets velocity (when moving) of stage in mm/s
        """
        # maximum velocity of instrument
        max_vel = 2.6
        if self.settings['velocity'] > max_vel:
            self.settings['velocity'] = max_vel

        velocity = ctypes.c_int(int(self._velocity_encoder2mm_conversion_factor * self.settings['velocity']))
        acceleration = ctypes.c_int(int(self._acceleration))
        self._servo_library.CC_SetVelParams(self._serial_num, acceleration, velocity)

    def set_acceleration(self, acc):

        """
        sets acceleration (when moving) of stage in mm/s^2
        """
        # maximum velocity of instrument
        max_acc = 4
        if acc > max_acc:
            acc = max_acc

        velocity = ctypes.c_int(int(self._velocity_encoder2mm_conversion_factor * self.settings['velocity']))
        acceleration = ctypes.c_int(int(self._acceleration_encoder2mm_conversion_factor * acc))
        self._servo_library.CC_SetVelParams(self._serial_num, acceleration, velocity)

    def update(self, settings, verbose=True):
        super(KDC001, self).update(settings)
        if self._settings_initialized:
            for key, value in settings.items():
                if key == 'position':
                    self.set_position(verbose=verbose)
                elif key == 'velocity':
                    self.set_velocity(verbose=verbose)

    def read_probes(self, key = None):
        assert key in list(self._PROBES.keys())
        assert isinstance(key, str)

        if key in ['position']:
            return (self.get_position(True))
        if key in ['serial_number']:
            return (self._serial_num_int)
        elif key in ['velocity']:
            # todo(emma): implement get_velocity
            pass

    def __del__(self):
        # when done with device we have to close the connections
        # called when GUI is closed
        self._close_device()

class B26KDC001x(KDC001):
    '''

    Same as KDC001 except adds safety limits and specifies serial number for the x axis

    '''
    _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
        #Parameter('serial_number', 27002905, int, 'serial number written on device'),
        Parameter('serial_number', 27501971, int, 'serial number written on device'),
        Parameter('position', 3.0, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 0.0, float,
                  'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used'),
    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        self.max_pos = 25.5
        self.min_pos = -3.
        super(B26KDC001x, self).__init__()

    def set_position(self, verbose=False):
        '''

        sets position of the stage x if safety limits are met

        '''

        # check if safety limits are met
        if self.settings['position'] <= self.max_pos and self.settings['position'] >= self.min_pos:
            super(B26KDC001x, self).set_position(verbose=verbose)
        else:
            print('didnt make the safety cut! doing nothing')
          #  raise AttributeError('position is outside safety limits!! Doing nothing.')

class B26KDC001y(KDC001):
    '''

    Same as KDC001 except adds safety limits and specifies serial number for the y axis

    '''
    _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
        Parameter('serial_number', 27001862, int, 'serial number written on device'),
        Parameter('position', 25.0, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 0.0, float,
                  'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        self.max_pos = 25.5
        self.min_pos = -3.
        super(B26KDC001y, self).__init__()

    def set_position(self, verbose=True):
        '''

        sets position of the stage y if safety limits are met

        '''

        # check if safety limits are met
        if self.settings['position'] <= self.max_pos and self.settings['position'] >= self.min_pos:
            super(B26KDC001y, self).set_position(verbose=verbose)
        else:
            print('didnt make the safety cut! doing nothing')
          #  raise AttributeError('position is outside safety limits!! Doing nothing.')


class B26KDC001z(KDC001):
    '''

    Same as KDC001 except adds safety limits and specifies serial number for the y axis

    '''
    _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
        Parameter('serial_number', 27001863, int, 'serial number written on device'),
        Parameter('position', 25.0, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 0.0, float,
                  'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        self.max_pos = 25.5
        self.min_pos = -3.
        super(B26KDC001z, self).__init__()

    def set_position(self, verbose=True):
        '''

        sets position of the stage y if safety limits are met

        '''

        # check if safety limits are met
        if self.settings['position'] <= self.max_pos and self.settings['position'] >= self.min_pos:
            super(B26KDC001z, self).set_position(verbose=verbose)
        else:
            print('didnt make the safety cut! doing nothing')
          #  raise AttributeError('position is outside safety limits!! Doing nothing.')


class LockboxToggleArm(KDC001):
    '''

    Same as KDC001 except adds safety limits and specifies serial number for the y axis

    '''
    _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
        Parameter('serial_number', 27002905, int, 'serial number written on device'),
        Parameter('position', 25, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 0.0, float,
                  'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used'),
        Parameter('lock', False, bool, 'toggle locking lever on front panel')
    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        self.max_pos = 35  # for mal/warm1, MM 20190813
        self.min_pos = 25
        super(LockboxToggleArm, self).__init__()

    def set_position(self, verbose=False):
        '''

        sets position of the stage z if safety limits are met

        '''

        # check if safety limits are met
        if self.settings['position'] <= self.max_pos and self.settings['position'] >= self.min_pos:
            self._set_position(verbose=verbose)
        else:
            print('didnt make the safety cut! doing nothing')
          #  raise AttributeError('position is outside safety limits!! Doing nothing.')

    def _set_position(self, verbose=False):

        """
        sets position of device in mm
        """
        position = ctypes.c_int(int(self._angle_encoder2deg_conversion_factor * self.settings['position']))
        self._servo_library.CC_MoveToPosition(self._serial_num, position)  # command sent to device
        message_type = ctypes.c_ushort(0)
        message_id = ctypes.c_ushort(0)
        message_data = ctypes.c_ulong(0)
        if verbose:
            print('Device is moving to indicated position ({0} deg)'.format(self.settings['position']))

        # wait for it to actually move
        self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
        while message_type.value != 2 or message_id.value != 1:
            self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)

        if verbose:
            position = self.get_position(verbose=True)
            print('Device now at position {0} deg'.format(position))

    def get_position(self, verbose=False):

        """
        returns position of stage in deg
        """

        position = self._servo_library.CC_GetPosition(
            self._serial_num) / self._angle_encoder2deg_conversion_factor
        if verbose:
            print('Position of device is currently %.3f deg' % position)
        return position

    def update(self, settings, verbose=True):
        super(KDC001, self).update(settings)
        if self._settings_initialized:
            for key, value in settings.items():
                if key == 'position':
                    self.set_position(verbose=verbose)
                elif key == 'velocity':
                    self.set_velocity(verbose=verbose)
                elif key == 'lock':
                    if settings['lock']:
                        self.settings['position'] = self.min_pos
                        self.set_position()
                    else:
                        self.settings['position'] = self.max_pos
                        self.set_position()


if __name__ == '__main__':
    a = LockboxToggleArm()
    print(a.get_position(verbose=True))
    # a.settings['position'] = 26
    # a.set_position()
    print(a.get_position(verbose=True))
    a._close_device()
