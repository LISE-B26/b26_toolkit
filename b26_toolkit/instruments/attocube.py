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
    _DEFAULT_SETTINGS = Parameter([
        Parameter('x',
                  [
                      Parameter('on', False, [True, False], 'x axis on'),
                      Parameter('pos', 0.0, float, 'x axis position in um'),
                      Parameter('voltage', 30, float, 'voltage on x axis'),
                      Parameter('freq', 100, float, 'x frequency in Hz')
                  ]
                  ),
        Parameter('y',
                  [
                      Parameter('on', False, [True, False], 'y axis on'),
                      Parameter('pos', 0, float, 'y axis position in um'),
                      Parameter('voltage', 30, float, 'voltage on y axis'),
                      Parameter('freq', 100, float, 'y frequency in Hz')
                  ]
                  ),
        Parameter('z',
                  [
                      Parameter('on', False, [True, False], 'z axis on'),
                      Parameter('pos', 0, float, 'x axis position in um'),
                      Parameter('voltage', 30, float, 'voltage on x axis'),
                      Parameter('freq', 100, float, 'x frequency in Hz')
                  ]
                  )
    ])

    _AXES = ['x', 'y', 'z']

    def _toggle_axis(self, axis, on):
        raise NotImplementedError

    def _set_frequency(self, axis, freq):
        raise NotImplementedError

    def _get_frequency(self, axis):
        raise NotImplementedError

    def _set_amplitude(self, axis, amplitude):
        raise NotImplementedError

    def _get_amplitude(self, axis):
        raise NotImplementedError

    def _cap_measure(self, axis):
        raise NotImplementedError

    def step(self, axis, dir):
        raise NotImplementedError

    def multistep(self, axis, num_steps):
        raise NotImplementedError

    def _convert_axis(self, axis):
        raise NotImplementedError

class ANC300(Attocube):

    _DEFAULT_SETTINGS = Parameter([
        Parameter('port', 'COM23', str, 'serial port on which to connect'),
        Parameter('baudrate', 9600, int, 'baudrate of connection'),
        Parameter('timeout', 1., float, 'connection timeout in seconds'),
        Parameter('x_voltage', 30, float, 'voltage on x axis'),
        Parameter('x_freq', 100, float, 'x frequency in Hz'),
        Parameter('y_voltage', 30, float, 'voltage on y axis'),
        Parameter('y_freq', 100, float, 'y frequency in Hz'),
        Parameter('z_voltage', 30, float, 'voltage on x axis'),
        Parameter('z_freq', 100, float, 'x frequency in Hz')
    ])
    _WRITE_TIMEOUT = 0.03 # seconds

    def __init__(self, name=None, settings=None):

        super().__init__(name, settings)
        self._is_connected = False
        try:
            self._connect(self.settings['port'],
                         self.settings['baudrate'],
                         self.settings['timeout'])
            # self.update(self.settings)
        except Exception as e:
            print(('Attocube not detected. Check connection.', UserWarning))
            raise e

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
                if split_key[1] == 'voltage':
                    self._set_amplitude(self._convert_axis(split_key[0]), value)
                elif split_key[1] == 'freq':
                    self._set_frequency(self._convert_axis(split_key[0]), value)
                else:
                    raise ValueError('No such key')

    def _connect(self, port, baudrate=9600, timeout=1.):
        self.ser = serial.Serial(port=port,
                                 baudrate=baudrate,
                                 timeout=timeout,
                                 write_timeout=self._WRITE_TIMEOUT)
        self._is_connected = True

    def __del__(self):
        if self._is_connected:
            self.ser.close()

    def _toggle_axis(self, axis, on):
        #TODO
        return

    def _set_frequency(self, axis, freq):
        self.ser.write('setf {} {}\n'.format(axis, freq).encode())
        self._get_OK()

    def _get_frequency(self, axis):
        self.ser.write('getf {}\n'.format(axis).encode())
        self.ser.readline()
        reply = self.ser.readline().decode()
        if 'frequency = ' not in reply:
            self._get_OK()
            raise Exception
        elif reply == 'OK\r\n':
            raise Exception

        freq = float(reply.split(' ')[2])
        self._get_OK()
        return freq

    def _set_amplitude(self, axis, amplitude):
        self.ser.write('setv {} {}\n'.format(axis, amplitude).encode())
        self._get_OK()

    def _get_amplitude(self, axis):
        self.ser.write('getv {}\n'.format(axis).encode())
        self.ser.readline()
        reply = self.ser.readline().decode()
        if 'voltage = ' not in reply:
            self._get_OK()
            raise Exception
        elif reply == 'OK\r\n':
            raise Exception

        voltage = float(reply.split(' ')[2])
        self._get_OK()
        return voltage

    def _cap_measure(self, axis):
        self.ser.write('getc {}\n'.format(axis).encode())
        # wait
        self.ser.write('capw {}\n'.format(axis).encode())

        self.ser.readline()
        reply = self.ser.readline().decode()
        if 'capacitance = ' not in reply:
            self._get_OK()
            raise Exception
        elif reply == 'OK\r\n':
            raise Exception

        cap = float(reply.split(' ')[2])

        self._get_OK()
        self._get_OK()
        return cap

    def step(self, axis, dir):
        if dir == 0:
            self.multistep(axis, 1)
        elif dir == 1:
            self.multistep(axis, -1)
        else:
            print('dir parameter is incorrect value')
            raise ValueError

    def multistep(self, axis, num_steps):
        axis = self._convert_axis(axis)
        self.ser.write('setm {} stp\n'.format(axis).encode())
        self._get_OK()

        if num_steps > 0:
            self.ser.write('stepu {} {}\n'.format(axis, num_steps).encode())
        else:
            self.ser.write('stepd {} {}\n'.format(axis, -num_steps).encode())


        self.ser.write('stepw {}\n'.format(axis).encode())
        self._get_OK()
        self._get_OK()


    def _convert_axis(self, axis):
        if axis not in self._AXES:
            raise ValueError('No such axis available')

        return ANC300_axes[axis]

    def _get_OK(self):
        line1 = self.ser.readline().decode()
        line2 = self.ser.readline().decode()
        if 'OK\r\n' not in [line1, line2]:
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

# class ANC300XY(ANC300):
#     _DEFAULT_SETTINGS = Parameter([
#         Parameter('port', 'COM3', str, 'serial port on which to connect'),
#         Parameter('baudrate', 9600, int, 'baudrate of connection'),
#         Parameter('timeout', 1., float, 'connection timeout in seconds'),
#         Parameter('x',
#                   [
#                       Parameter('on', False, [True, False], 'x axis on'),
#                       Parameter('voltage', 30, float, 'voltage on x axis'),
#                       Parameter('freq', 1000, float, 'x frequency in Hz')
#                   ]
#                   ),
#         Parameter('y',
#                   [
#                       Parameter('on', False, [True, False], 'y axis on'),
#                       Parameter('voltage', 30, float, 'voltage on y axis'),
#                       Parameter('freq', 1000, float, 'y frequency in Hz')
#                   ]
#                   ),
#     ])
#
#     _AXES = ['x', 'y']
#
#     @property
#     def _PROBES(self):
#         return{
#             'x_voltage': 'the voltage of the x direction (with respect to the camera)',
#             'x_freq': 'the frequency of the x direction (with respect to the camera)',
#             'x_cap': 'the capacitance of the piezo in the x direction (with respect to the camera)',
#             'y_voltage': 'the voltage of the y direction (with respect to the camera)',
#             'y_freq': 'the frequency of the y direction (with respect to the camera)',
#             'y_cap': 'the capacitance of the piezo in the y direction (with respect to the camera)',
#         }

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

    # _DEFAULT_SETTINGS = Parameter([
    #     Parameter('x',
    #               [
    #                   Parameter('on', False, [True, False], 'x axis on'),
    #                   Parameter('pos', 0.0, float, 'x axis position in um'),
    #                   Parameter('voltage', 30, float, 'voltage on x axis'),
    #                   Parameter('freq', 1000, float, 'x frequency in Hz')
    #               ]
    #               ),
    #     Parameter('y',
    #               [
    #                   Parameter('on', False, [True, False], 'y axis on'),
    #                   Parameter('pos', 0, float, 'y axis position in um'),
    #                   Parameter('voltage', 30, float, 'voltage on y axis'),
    #                   Parameter('freq', 1000, float, 'y frequency in Hz')
    #               ]
    #               ),
    #     Parameter('z',
    #               [
    #                   Parameter('on', False, [True, False], 'z axis on'),
    #                   Parameter('pos', 0, float, 'x axis position in um'),
    #                   Parameter('voltage', 30, float, 'voltage on x axis'),
    #                   Parameter('freq', 1000, float, 'x frequency in Hz')
    #               ]
    #               )
    # ])

    # _AXES = ['x', 'y', 'z']

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
        super(Attocube, self).update(settings)
        for key, value in settings.items():
            if isinstance(value, dict) and key in self._AXES:
                for sub_key, sub_value in sorted(value.items()):
                    if sub_key == 'on':
                        self._toggle_axis(self._convert_axis(key), sub_value)
                    elif sub_key == 'pos':
                        self.move_absolute(self._convert_axis(key), sub_value)
                    elif sub_key == 'voltage':
                        self._set_amplitude(self._convert_axis(key), sub_value)
                    elif sub_key == 'freq':
                        self._set_frequency(self._convert_axis(key), sub_value)
                    else:
                        raise ValueError('No such key')
            else:
                raise ValueError('No such key')


    # def update(self, settings):
    #     '''
    #     Updates the internal settings, as well as turning the attocube channel on or off, updating
    #     voltage or frequency, or moving to the given position
    #     Args:
    #         settings: a dictionary in the same form as settings with the new values
    #     '''
    #     super(Attocube, self).update(settings)
    #     for key, value in settings.items():
    #         if isinstance(value, dict) and key in self._AXES:
    #             for sub_key, sub_value in sorted(value.items()):
    #                 if sub_key == 'on':
    #                     self._toggle_axis(self._convert_axis(key), sub_value)
    #                 elif sub_key == 'pos':
    #                     self._move_absolute(self._convert_axis(key), sub_value)
    #                 elif sub_key == 'voltage':
    #                     self._set_amplitude(self._convert_axis(key), sub_value)
    #                 elif sub_key == 'freq':
    #                     self._set_frequency(self._convert_axis(key), sub_value)
    #                 else:
    #                     raise ValueError('No such key')
    #         else:
    #             raise ValueError('No such key')


    # @property
    # def _PROBES(self):
    #     return{
    #         'x_pos': 'the position the x direction (with respect to the camera) in um',
    #         'x_voltage': 'the voltage of the x direction (with respect to the camera)',
    #         'x_freq': 'the frequency of the x direction (with respect to the camera)',
    #         'x_cap': 'the capacitance of the piezo in the x direction (with respect to the camera)',
    #         'y_pos': 'the position the y direction (with respect to the camera) in um',
    #         'y_voltage': 'the voltage of the y direction (with respect to the camera)',
    #         'y_freq': 'the frequency of the y direction (with respect to the camera)',
    #         'y_cap': 'the capacitance of the piezo in the y direction (with respect to the camera)',
    #         'z_pos': 'the position the z direction (with respect to the camera) in um',
    #         'z_voltage': 'the voltage of the z direction (with respect to the camera)',
    #         'z_freq': 'the frequency of the z direction (with respect to the camera)',
    #         'z_cap': 'the capacitance of the piezo in the z direction (with respect to the camera)'
    #     }

    # def read_probes(self, key):
    #     assert key in list(self._PROBES.keys())
    #     assert isinstance(key, str)
    #
    #     if key in [el + '_pos' for el in self._AXES]:#['x_pos', 'y_pos', 'z_pos']:
    #         return self._get_position(self._convert_axis(key[0]))
    #     elif key in [el + '_voltage' for el in self._AXES]:
    #         return self._get_amplitude(self._convert_axis(key[0]))
    #     elif key in [el + '_freq' for el in self._AXES]:
    #         return self._get_frequency(self._convert_axis(key[0]))
    #     elif key in [el + '_cap' for el in self._AXES]:
    #         return self._cap_measure(self._convert_axis(key[0]))

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

# class ANC350XY(ANC350):
#     _DEFAULT_SETTINGS = Parameter([
#         Parameter('x',
#                   [
#                       Parameter('on', False, [True, False], 'x axis on'),
#                       Parameter('pos', 0.0, float, 'x axis position in um'),
#                       Parameter('voltage', 30, float, 'voltage on x axis'),
#                       Parameter('freq', 1000, float, 'x frequency in Hz')
#                   ]
#                   ),
#         Parameter('y',
#                   [
#                       Parameter('on', False, [True, False], 'y axis on'),
#                       Parameter('pos', 0, float, 'y axis position in um'),
#                       Parameter('voltage', 30, float, 'voltage on y axis'),
#                       Parameter('freq', 1000, float, 'y frequency in Hz')
#                   ]
#                   )
#     ])
#
#     _AXES = ['x', 'y']
#
#     @property
#     def _PROBES(self):
#         return{
#             'x_pos': 'the position the x direction (with respect to the camera) in um',
#             'x_voltage': 'the voltage of the x direction (with respect to the camera)',
#             'x_freq': 'the frequency of the x direction (with respect to the camera)',
#             'x_cap': 'the capacitance of the piezo in the x direction (with respect to the camera)',
#             'y_pos': 'the position the y direction (with respect to the camera) in um',
#             'y_voltage': 'the voltage of the y direction (with respect to the camera)',
#             'y_freq': 'the frequency of the y direction (with respect to the camera)',
#             'y_cap': 'the capacitance of the piezo in the y direction (with respect to the camera)',
#         }

if __name__ == '__main__':
    print((os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.txt')))

    try:
        a = ANC300XY()
        a.multistep('x', 100)
    except Exception:
        print('yike')
    # a.update({'x': {'voltage': 20}})
    print((a, a.is_connected))
