"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    b26_toolkit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    b26_toolkit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with b26_toolkit.  If not, see <http://www.gnu.org/licenses/>.
"""
import ctypes, serial
import time, os
import warnings
from pylabcontrol.core.read_write_functions import get_config_value
from pylabcontrol.core import Instrument, Parameter

int32 = ctypes.c_long
uInt32 = ctypes.c_ulong
uInt64 = ctypes.c_ulonglong
float64 = ctypes.c_double

#define built-in error codes
NCB_Ok = 0
NCB_Error = -1
NCB_Timeout = 1
NCB_NotConnected = 2
NCB_DriverError = 3
NCB_BootIgnored = 4
NCB_FileNotFound = 5
NCB_InvalidParam = 6
NCB_DeviceLocked = 7
NCB_NotSpecifiedParam = 8

# converts x,y,z to axis number in controller

ANC350_axes = {
    'x': int32(1),
    'y': int32(2),
    'z': int32(0)
}

ANC300_axes = {
    'x': 1,
    'y': 2,
    'z': 3
}

# c struct used as return type for some old_functions
class PositionerInfo(ctypes.Structure):
    _fields_ = [(("id"), ctypes.c_int32), (("locked"), ctypes.c_bool)]

class Attocube(Instrument):
    '''
    Generic class for implementing attocube controllers. Needs to be different since ANC300 communicates through serial
    and ANC350 with DLL.
    '''

    _DEFAULT_SETTINGS = Parameter([
        Parameter('x_on', False, [True, False], 'toggle axis on and off'),
        Parameter('x_voltage', 30, float, 'voltage on x axis'),
        Parameter('x_freq', 100, float, 'x frequency in Hz'),
        Parameter('x_pos', 0., float, 'x position in um'),
        Parameter('y_on', False, [True, False], 'toggle axis on and off'),
        Parameter('y_voltage', 30, float, 'voltage on y axis'),
        Parameter('y_freq', 100, float, 'y frequency in Hz'),
        Parameter('y_pos', 0., float, 'y position in um'),
        Parameter('z_on', False, [True, False], 'toggle axis on and off'),
        Parameter('z_voltage', 30, float, 'voltage on z axis'),
        Parameter('z_freq', 100, float, 'z frequency in Hz'),
        Parameter('z_pos', 0., float, 'z position in um')
    ])

    _AXES = ['x', 'y', 'z']

    def _set_frequency(self, axis, freq):
        '''
        Sets frequency of attocube axis
        :param axis: axis number to set (defined by controller)
        :param freq: frequency to set axis to in Hz (int)
        '''
        raise NotImplementedError

    def _get_frequency(self, axis):
        '''
        Gets frequency of attocube axis
        :param axis: axis number to set (defined by controller)
        '''
        raise NotImplementedError

    def _set_amplitude(self, axis, amplitude):
        '''
        Sets amplitude of attocube axis
        :param axis: axis number to set (defined by controller)
        :param amplitude: amplitude to set axis to in V (float)
        '''
        raise NotImplementedError

    def _get_amplitude(self, axis):
        '''
        Gets amplitude of attocube axis
        :param axis: axis number to set (defined by controller)
        '''
        raise NotImplementedError

    def _cap_measure(self, axis):
        '''
        Measures capacitance for specified axis
        :param axis: axis number to set (defined by controller)
        '''
        raise NotImplementedError

    def step(self, axis, dir):
        '''
        Take single step
        :param axis: axis to take step along (str: x, y, z)
        :param dir: direction to take step in (int: 0 for positive, 1 for negative)
        '''
        raise NotImplementedError

    def multistep(self, axis, num_steps):
        '''
        Take multiple steps
        :param axis: axis to take step along (str: x, y, z)
        :param num_steps: number of steps to take. num_steps < 0 for negative direction (int)
        '''
        raise NotImplementedError

    def _convert_axis(self, axis):
        '''
        Convert axis name to number
        :param axis: axis name (x, y, z) (str)
        :return: axis number (defined by controller)
        '''
        raise NotImplementedError

class ANC300(Attocube):
    '''
    Class to control an attocube using a supplied ANC300 controller. Has been tested controlling a stack of two
    ANPx101res and one ANPz101res, but it should work with any controllers supporting the same low level serial
    commands. The commands can be found in the ANC300 manual. The class communicates with the device over the USB slave
    port.
    '''


    _DEFAULT_SETTINGS = Parameter([
        Parameter('port', 'COM15', str, 'serial port on which to connect'),
        Parameter('baudrate', 9600, int, 'baudrate of connection'),
        Parameter('timeout', 1., float, 'connection timeout in seconds'),
        Parameter('x_voltage', 30, float, 'voltage on x axis'),
        Parameter('x_freq', 100, int, 'x frequency in Hz'),
        Parameter('y_voltage', 30, float, 'voltage on y axis'),
        Parameter('y_freq', 100, int, 'y frequency in Hz'),
        Parameter('z_voltage', 30, float, 'voltage on x axis'),
        Parameter('z_freq', 100, int, 'x frequency in Hz')
    ])
    _WRITE_TIMEOUT = 0.03  # in seconds

    def __init__(self, name=None, settings=None):
        '''
        Connects to attocube controller ANC300 through USB using PySerial. Throws exception if unable to connect.
        baudrate was not specified in manual but 9600 seems to work.
        :param name: name of instruments
        :param settings: settings to update instrument with
        '''
        super().__init__(name, settings)
        self._is_connected = False
        try:
            self._connect(self.settings['port'],
                         self.settings['baudrate'],
                         self.settings['timeout'])
        except Exception as e:
            print(('Attocube not detected. Check connection.', UserWarning))
            raise e

    def update(self, settings):
        '''
        Updates the internal settings, updating voltage or frequency
        Args:
            settings: a dictionary in the same form as settings with the new values
        '''
        super().update(settings)
        for key, value in settings.items():
            split_key = key.split('_')

            # if axis setting
            if split_key[0] in self._AXES:
                # updates axis parameters
                if split_key[1] == 'voltage':
                    self._set_amplitude(self._convert_axis(split_key[0]), value)
                elif split_key[1] == 'freq':
                    self._set_frequency(self._convert_axis(split_key[0]), value)
                else:
                    raise ValueError('No such key')

            # update connection parameters
            elif key == 'port':
                self.ser.port = value
            elif key == 'baudrate':
                self.ser.baudrate = value
            elif key == 'timeout':
                self.ser.timeout = value
            else:
                raise ValueError('No such key')

    def _connect(self, port, baudrate=9600, timeout=1.):
        '''
        Connects to ANC300 through USB serial connection. Throws exception if unable to connect
        :param port: virtual port to connect to (usually COM)
        :param baudrate: baudrate of connection
        :param timeout: connection timeout in seconds
        '''
        self.ser = serial.Serial(port=port,
                                 baudrate=baudrate,
                                 timeout=timeout,
                                 write_timeout=self._WRITE_TIMEOUT)
        self._is_connected = True

    def __del__(self):
        '''
        close connection upon destruction of instance
        '''
        if self._is_connected:
            self.ser.close()

    def _set_frequency(self, axis, freq):
        '''
        Sets frequency of attocube axis
        :param axis: axis number to set (int)
        :param freq: frequency to set axis to in Hz (int)
        '''
        self.ser.write('setf {} {}\n'.format(axis, freq).encode())
        self._get_OK()

    def _get_frequency(self, axis):
        '''
        Gets frequency of attocube axis
        :param axis: axis number to set (int)
        '''
        return int(self._get_param(axis, 'getf', 'frequency'))

    def _set_amplitude(self, axis, amplitude):
        '''
        Sets amplitude of attocube axis
        :param axis: axis number to set (int)
        :param amplitude: amplitude to set axis to in V (float)
        '''
        self.ser.write('setv {} {}\n'.format(axis, amplitude).encode())
        self._get_OK()

    def _get_amplitude(self, axis):
        '''
        Gets amplitude of attocube axis
        :param axis: axis number to set (int)
        '''
        return self._get_param(axis, 'getv', 'voltage')

    def _cap_measure(self, axis):
        '''
        Measures capacitance for specified axis
        :param axis: axis number to set (int)
        '''

        # Start capacitance measurement
        self.ser.write('setm {} cap\n'.format(axis).encode())

        # Wait for capacitance measurement to finish
        self.ser.write('capw {}\n'.format(axis).encode())

        # Get OK for starting and waiting
        self._get_OK()
        self._get_OK()

        # Acquire capacitance value
        return self._get_param(axis, 'getc', 'capacitance')

    def _get_param(self, axis, command, param):
        '''
        Get parameter (frequency, amplitude, etc) from attocube controller
        :param axis: axis number to acquire parameter from (int)
        :param command: command for getting parameter (i.e. getf, getc) (str)
        :param param: name of paramter to be read (str)
        :return: parameter value
        '''
        # Send command to get param
        self.ser.write((command + ' {}\n').format(axis).encode())

        # read line containing command
        self.ser.readline()

        # read line that should contain value
        reply = self.ser.readline().decode()

        # If error or no value given, signal exception
        if param + ' = ' not in reply:
            self._get_OK()
            raise Exception
        elif reply == 'OK\r\n':
            raise Exception

        # get OK from command
        self._get_OK()

        # get value
        return float(reply.split(' ')[2])

    def step(self, axis, dir):
        '''
        Take single step
        :param axis: axis to take step along (str: x, y, z)
        :param dir: direction to take step in (int: 0 for positive, 1 for negative)
        '''
        # If positive direction
        if dir == 0:
            self.multistep(axis, 1)
        # If negative direction
        elif dir == 1:
            self.multistep(axis, -1)
        else:
            print('dir parameter is incorrect value')
            raise ValueError

    def multistep(self, axis, num_steps):
        '''
        Take multiple steps
        :param axis: axis to take step along (str: x, y, z)
        :param num_steps: number of steps to take. num_steps < 0 for negative direction (int)
        '''

        # Do nothing if no steps
        if num_steps == 0:
            return

        # Set to stepping mode and get OK
        axis = self._convert_axis(axis)
        self.ser.write('setm {} stp\n'.format(axis).encode())
        self._get_OK()

        # Perform steps
        if num_steps > 0:
            self.ser.write('stepu {} {}\n'.format(axis, num_steps).encode())
        else:
            self.ser.write('stepd {} {}\n'.format(axis, -num_steps).encode())

        # Wait until stepping is done
        self.ser.write('stepw {}\n'.format(axis).encode())

        # Get OK for stepping and waiting
        self._get_OK()
        self._get_OK()


    def _convert_axis(self, axis):
        '''
        Convert axis name to number
        :param axis: axis name (x, y, z) (str)
        :return: axis number
        '''
        if axis not in self._AXES:
            raise ValueError('No such axis available')

        return ANC300_axes[axis]

    def _get_OK(self):
        '''
        Read acknowledgment line after command.
        '''

        # Read command recently written and then what should be OK
        line1 = self.ser.readline().decode()
        line2 = self.ser.readline().decode()

        # If OK not received, probably error
        if 'OK\r\n' not in [line1, line2]:

            # Read line until OK received, then raise exception
            line1 = self.ser.readline().decode()
            while line1 != 'OK\r\n':
                line1 = self.ser.readline().decode()
            print('Error')
            raise Exception

    @property
    def _PROBES(self):
        return {
            'x_voltage': 'the voltage of the x direction (with respect to the camera)',
            'x_freq': 'the frequency of the x direction (with respect to the camera)',
            'x_cap': 'the capacitance of the piezo in the x direction (with respect to the camera)',
            'y_voltage': 'the voltage of the y direction (with respect to the camera)',
            'y_freq': 'the frequency of the y direction (with respect to the camera)',
            'y_cap': 'the capacitance of the piezo in the y direction (with respect to the camera)',
            'z_voltage': 'the voltage of the z direction (with respect to the camera)',
            'z_freq': 'the frequency of the z direction (with respect to the camera)',
            'z_cap': 'the capacitance of the piezo in the z direction (with respect to the camera)'
        }

    def read_probes(self, key):
        assert key in list(self._PROBES.keys())
        assert isinstance(key, str)


        if key in [el + '_voltage' for el in self._AXES]:
            return self._get_amplitude(self._convert_axis(key.split('_')[0]))
        elif key in [el + '_freq' for el in self._AXES]:
            return self._get_frequency(self._convert_axis(key.split('_')[0]))
        elif key in [el + '_cap' for el in self._AXES]:
            return self._cap_measure(self._convert_axis(key.split('_')[0]))

class ANC350(Attocube):
    '''
    Class to control an attocube using a supplied controller. Has been tested on an
    ANC350 controlling a stack of two ANPx101res and one ANPz101res, but it should
    work with any controllers supporting the same low level dll commands.
    The path to the dll should be set in config.txt in this folder
    Note that we use the 1.5 version of the dll, the 2.0 version cannot be read properly
    and may be written in a non-ctypes compatible language
    The class communicates with the device over USB.
    '''

    def __init__(self, name = None, settings = None):
        # Load DLL and check that attocube is connected to computer. If no DLL, continue to work but throw a warning
        # and all future calls will fail. If a DLL but no instrument connected, throw an error.
        dll_path = get_config_value('ATTOCUBE_DLL_PATH',os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.txt'))
        try:
            self.attocube = ctypes.WinDLL(dll_path)
            # self.attocube = ctypes.WinDLL('C:/Users/Experiment/Downloads/attocube/Software/ANC350_Software_v1.5.15/ANC350_DLL/Win_64Bit/pylabcontrol/anc350v2.dll')
            dll_detected = True
        except WindowsError:
            # make a fake Attocube instrument
            dll_detected = False
            warnings.warn("Attocube DLL not found. If it should be present, check the path.")

        if dll_detected == True:
            try:
                self.pi = PositionerInfo()
                dev_count = self.attocube.PositionerCheck(ctypes.byref(self.pi))
                device_handle = int32()
                self._check_error(self.attocube.PositionerConnect(0, ctypes.byref(device_handle)))
                self._check_error(self.attocube.PositionerClose(device_handle))
            except Exception:
                print(('Attocube not detected. Check connection.', UserWarning))

        super().__init__(name, settings)

    def update(self, settings):
        '''
        Updates the internal settings, as well as turning the attocube channel on or off, updating
        voltage or frequency, or moving to the given position
        Args:
            settings: a dictionary in the same form as settings with the new values
        '''
        super().update(settings)
        for key, value in settings.items():
            split_key = key.split('_')
            if split_key[0] in self._AXES:
                if split_key[1] == 'on':
                    self._toggle_axis(self._convert_axis(key), value)
                elif split_key[1] == 'pos':
                    self.move_absolute(self._convert_axis(key), value)
                elif split_key[1] == 'voltage':
                    self._set_amplitude(self._convert_axis(key), value)
                elif split_key[1] == 'freq':
                    self._set_frequency(self._convert_axis(key), value)
                else:
                    raise ValueError('No such key')
            else:
                raise ValueError('No such key')

    @property
    def is_connected(self):
        '''
        Check if attocube controller is connected
        Returns: True if controller is connected, false otherwise
        '''
        #connecting fails if device is not connected, so this catches that error
        try:
            device_handle = int32()
            self._check_error(self.attocube.PositionerConnect(0, ctypes.byref(device_handle)))
            self._check_error(self.attocube.PositionerClose(device_handle))
            return True
        except Exception:
            return False

    def _toggle_axis(self, axis, on):
        '''
        Turn axis on or off
        :param axis: axis_x, axis_y, or axis_z
        :param on: True or False
        '''
        device_handle = int32()
        self._check_error(self.attocube.PositionerConnect(0, ctypes.byref(device_handle)))
        self._check_error(self.attocube.PositionerSetOutput(device_handle, axis, ctypes.c_bool(on)))
        self._check_error(self.attocube.PositionerClose(device_handle))

    def _set_frequency(self, axis, freq):
        '''
        :param axis: axis_x, axis_y, or axis_z
        :param freq: frequency to set in Hz
        '''
        assert (freq <= 2000)
        device_handle = int32()
        self._check_error(self.attocube.PositionerConnect(0, ctypes.byref(device_handle)))
        self._check_error(self.attocube.PositionerFrequency(device_handle, axis, int32(int(freq))))
        self._check_error(self.attocube.PositionerClose(device_handle))

    def _get_frequency(self, axis):
        '''
        :param axis: axis_x, axis_y, or axis_z
        :return: current frequency of axis in Hz
        '''
        device_handle = int32()
        freq = int32()
        self._check_error(self.attocube.PositionerConnect(0, ctypes.byref(device_handle)))
        self._check_error(self.attocube.PositionerGetFrequency(device_handle, axis, ctypes.byref(freq)))
        self._check_error(self.attocube.PositionerClose(device_handle))
        return freq.value

    def _set_amplitude(self, axis, amplitude):
        '''
        :param axis: axis: axis_x, axis_y, or axis_z
        :param amplitude: amplitude in V
        '''
        assert(amplitude <= 60)
        device_handle = int32()
        amplitude *= 1000
        self._check_error(self.attocube.PositionerConnect(0, ctypes.byref(device_handle)))
        self._check_error(self.attocube.PositionerAmplitude(device_handle, axis, int32(int(amplitude))))
        self._check_error(self.attocube.PositionerClose(device_handle))

    def _get_amplitude(self, axis):
        '''
        :param axis: axis_x, axis_y, or axis_z
        :return: amplitude in V
        '''
        device_handle = int32()
        amplitude = int32()
        self._check_error(self.attocube.PositionerConnect(0, ctypes.byref(device_handle)))
        self._check_error(self.attocube.PositionerGetAmplitude(device_handle, axis, ctypes.byref(amplitude)))
        self._check_error(self.attocube.PositionerClose(device_handle))
        return (amplitude.value / 1000.0)

    def _get_position(self, axis):
        '''
        :param axis: axis_x, axis_y, or axis_z
        :return: position of axis in um
        '''
        device_handle = int32()
        position = int32()
        self._check_error(self.attocube.PositionerConnect(0, ctypes.byref(device_handle)))
        # wait command needed since polling rate of attocube is 20 Hz. Empirically determined that .2 is lowest value
        # that always works. No idea why no other function also needs this wait command
        time.sleep(.2)
        self._check_error(self.attocube.PositionerGetPosition(device_handle, axis, ctypes.byref(position)))
        self._check_error(self.attocube.PositionerClose(device_handle))
        return position.value/1000.0

    def _cap_measure(self, axis):
        '''
        :param axis: axis_x, axis_y, or axis_z
        :return: Capacitance in uF
        '''
        device_handle = int32()
        capacitance = int32()
        self._check_error(self.attocube.PositionerConnect(0, ctypes.byref(device_handle)))
        self._check_error(self.attocube.PositionerCapMeasure(device_handle, axis, ctypes.byref(capacitance)))
        self._check_error(self.attocube.PositionerClose(device_handle))
        return capacitance.value

    def move_absolute(self, axis, position):
        '''
        Precondition: Must set voltage and frequency sufficiently low that ANC's internal feedback will be able to
        settle on the appropriate position (ex. 7V, 100Hz). Otherwise, fluctuates around target position and never stops
        :param axis: axis_x, axis_y, or axis_z
        :param position: position of axis to move to in um
        '''
        device_handle = int32()
        self._check_error(self.attocube.PositionerConnect(0, ctypes.byref(device_handle)))
        self._check_error(self.attocube.PositionerMoveAbsolute(device_handle, self._convert_axis(axis), int32(int(position * 1000.0))))
        self._check_error(self.attocube.PositionerClose(device_handle))

    def step(self, axis, dir):
        '''
        Move a single step on the given axis in the given direction
        Args:
            axis: 'x', 'y', or 'z'
            dir: 0 for forwards, 1 for backwards

        '''
        device_handle = int32()
        axis = self._convert_axis(axis)
        self._check_error(self.attocube.PositionerConnect(0, ctypes.byref(device_handle)))
        self._check_error(self.attocube.PositionerMoveSingleStep(device_handle, axis, int32(dir)))
        self._check_error(self.attocube.PositionerClose(device_handle))

    def multistep(self, axis, num_steps):
        dir = int(num_steps < 0)
        device_handle = int32()
        axis = self._convert_axis(axis)

        self._check_error(self.attocube.PositionerConnect(0, ctypes.byref(device_handle)))
        for _ in range(abs(num_steps)):
            self._check_error(self.attocube.PositionerMoveSingleStep(device_handle, axis, int32(dir)))
        self._check_error(self.attocube.PositionerClose(device_handle))


    def _convert_axis(self, axis):
        if axis not in self._AXES:
            raise ValueError('No such axis available')

        return ANC350_axes[axis]

    @staticmethod
    def _check_error(code):
        '''
        Translates error codes to human readable message
        :param code: input error code (integer 0-8)
        :poststate: message printed to stdout
        '''
        if(code == NCB_Ok):
            return
        elif(code == NCB_BootIgnored):
            print( "Warning: boot ignored\n" )
            raise Exception
        elif(code == NCB_Error):
            print( "Error: unspecific\n" )
            raise Exception
        elif(code == NCB_Timeout):
            print( "Error: comm. timeout\n" )
            raise Exception
        elif(code == NCB_NotConnected):
            print( "Error: not connected\n" )
            raise Exception
        elif(code == NCB_DriverError):
            print( "Error: driver error\n" )
            raise Exception
        elif(code == NCB_FileNotFound):
            print( "Error: file not found\n" )
            raise Exception
        elif(code == NCB_InvalidParam):
            print( "Error: invalid parameter\n" )
            raise Exception
        elif(code == NCB_DeviceLocked):
            print( "Error: device locked\n" )
            raise Exception
        elif(code == NCB_NotSpecifiedParam):
            print( "Error: unspec. parameter\n" )
            raise Exception
        else:
            print( "Error: unknown\n" )
            raise Exception

    @property
    def _PROBES(self):
        return {
            'x_pos': 'the position the x direction (with respect to the camera) in um',
            'x_voltage': 'the voltage of the x direction (with respect to the camera)',
            'x_freq': 'the frequency of the x direction (with respect to the camera)',
            'x_cap': 'the capacitance of the piezo in the x direction (with respect to the camera)',
            'y_pos': 'the position the y direction (with respect to the camera) in um',
            'y_voltage': 'the voltage of the y direction (with respect to the camera)',
            'y_freq': 'the frequency of the y direction (with respect to the camera)',
            'y_cap': 'the capacitance of the piezo in the y direction (with respect to the camera)',
            'z_pos': 'the position the z direction (with respect to the camera) in um',
            'z_voltage': 'the voltage of the z direction (with respect to the camera)',
            'z_freq': 'the frequency of the z direction (with respect to the camera)',
            'z_cap': 'the capacitance of the piezo in the z direction (with respect to the camera)'
        }

    def read_probes(self, key):
        assert key in list(self._PROBES.keys())
        assert isinstance(key, str)

        if key in [el + '_pos' for el in self._AXES]:  # ['x_pos', 'y_pos', 'z_pos']:
            return self._get_position(self._convert_axis(key[0]))
        elif key in [el + '_voltage' for el in self._AXES]:
            return self._get_amplitude(self._convert_axis(key[0]))
        elif key in [el + '_freq' for el in self._AXES]:
            return self._get_frequency(self._convert_axis(key[0]))
        elif key in [el + '_cap' for el in self._AXES]:
            return self._cap_measure(self._convert_axis(key[0]))

if __name__ == '__main__':
    print((os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.txt')))

    try:
        a = ANC300()
        a.multistep('x', 100)
    except Exception:
        print('yike')
    # a.update({'x': {'voltage': 20}})
    print((a, a.is_connected))
