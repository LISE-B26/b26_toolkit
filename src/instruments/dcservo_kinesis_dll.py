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
from PyLabControl.src.core.read_write_functions import get_config_value

dll_path = get_config_value('KINESIS_DLL_PATH', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.txt'))
sys.path.insert(0,dll_path)
# JG: July 27 2016 uncommented folowing line: don't use import *!
# from PyLabControl.src.core.instruments import *
from PyLabControl.src.core import Parameter, Instrument

# ctypes DLL load failed: Probably a C++ dll was provided, which is incompatable with ctypes, possibly due to name
# mangling. Instead, we use the .net framework with python for .net to interface with the dll
#ctypes.cdll.LoadLibrary("Thorlabs.MotionControl.TCube.DCServo.dll")

# makes each dll, corresponding to a namespace, avaliable to python at runtime
clr.AddReference('ThorLabs.MotionControl.DeviceManagerCLI')
clr.AddReference('Thorlabs.MotionControl.TCube.DCServoCLI')
clr.AddReference('Thorlabs.MotionControl.KCube.DCServoCLI')
clr.AddReference('System')

# imports classes from the namespaces. All read as unresolved references because python doesn't know about the dlls
# until runtime
from Thorlabs.MotionControl.DeviceManagerCLI import DeviceManagerCLI
from Thorlabs.MotionControl.TCube.DCServoCLI import TCubeDCServo
from Thorlabs.MotionControl.TCube.DCServoCLI import KCubeDCServo
# adds .NET stuctures corresponding to primitives
from System import Decimal, Double

class ThorlabsServo(Instrument):
    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 83832028, int, 'serial number written on device'),
        Parameter('position', 0, float, 'servo position (from 0 to 6 in mm)'),
        Parameter('velocity', 0, float, 'servo maximum velocity in mm/s')
    ])

    def __init__(self, name = None, settings = None):
        raise NotImplementedError

    def update(self, settings):
        '''
        Updates internal settings, as well as the position and velocity set on the physical device
        Args:
            settings: A dictionary in the form of settings as seen in default settings
        '''
        super(ThorlabsServo, self).update(settings)
        for key, value in settings.iteritems():
            if key == 'position':
                self._move_servo(value)
            elif key == 'velocity':
                self._set_velocity(value)

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
        DeviceManagerCLI.BuildDeviceList()
        return(str(self.settings['serial_number']) in DeviceManagerCLI.GetDeviceList(self.Servo.DevicePrefix))

    def __del__(self):
        '''
        Cleans up TDC001 connection
        :PostState: TDC001 is disconnected
        '''
        self.device.StopPolling()
        self.device.Disconnect()

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
            self.device.MoveTo(self._Py_Decimal(position), 60000)
        except Exception:
            print("Failed to move to position")
            raise

    def _get_position(self):
        '''
        :return: position of servo
        '''
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
    '''
    Class to control the thorlabs TDC001 servo. Note that ALL DLL FUNCTIONS TAKING NUMERIC INPUT REQUIRE A SYSTEM.DECIMAL
    VALUE. Check help doc at C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.DotNet_API for the DLL api.
    The class communicates with the device over USB.
    '''

    def __init__(self, name = None, settings = None):
        super(ThorlabsServo, self).__init__(name, settings)
        self.Servo = TCubeDCServo
        try:
            DeviceManagerCLI.BuildDeviceList()
            serial_number_list = DeviceManagerCLI.GetDeviceList(self.Servo.DevicePrefix)
        except (Exception):
            print("Exception raised by BuildDeviceList")
        if not (str(self.settings['serial_number']) in serial_number_list):
            print(str(self.settings['serial_number']) + " is not a valid serial number")
            raise

        self.device = self.Servo.CreateTCubeDCServo(str(self.settings['serial_number']))
        if(self.device == None):
            print(self.settings['serial_number'] + " is not a TCubeDCServo")
            raise

        try:
            self.device.Connect(str(self.settings['serial_number']))
        except Exception:
            print('Failed to open device ' + str(self.settings['serial_number']))
            raise

        if not self.device.IsSettingsInitialized():
            try:
                self.device.WaitForSettingsInitialized(5000)
            except Exception:
                print("Settings failed to initialize")
                raise

        self.device.StartPolling(250)

        motorSettings = self.device.GetMotorConfiguration(str(self.settings['serial_number']))
        currentDeviceSettings = self.device.MotorDeviceSettings

class KDC001(ThorlabsServo):
    '''
    Class to control the thorlabs TDC001 servo. Note that ALL DLL FUNCTIONS TAKING NUMERIC INPUT REQUIRE A SYSTEM.DECIMAL
    VALUE. Check help doc at C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.DotNet_API for the DLL api.
    The class communicates with the device over USB.
    '''

    def __init__(self, name=None, settings=None):
        super(ThorlabsServo, self).__init__(name, settings)
        self.Servo = KCubeDCServo
        try:
            DeviceManagerCLI.BuildDeviceList()
            serial_number_list = DeviceManagerCLI.GetDeviceList(self.Servo.DevicePrefix)
        except (Exception):
            print("Exception raised by BuildDeviceList")
        if not (str(self.settings['serial_number']) in serial_number_list):
            print(str(self.settings['serial_number']) + " is not a valid serial number")
            raise

        self.device = self.Servo.CreateKCubeDCServo(str(self.settings['serial_number']))
        if (self.device == None):
            print(self.settings['serial_number'] + " is not a TCubeDCServo")
            raise

        try:
            self.device.Connect(str(self.settings['serial_number']))
        except Exception:
            print('Failed to open device ' + str(self.settings['serial_number']))
            raise

        if not self.device.IsSettingsInitialized():
            try:
                self.device.WaitForSettingsInitialized(5000)
            except Exception:
                print("Settings failed to initialize")
                raise

        self.device.StartPolling(250)

        motorSettings = self.device.GetMotorConfiguration(str(self.settings['serial_number']))
        currentDeviceSettings = self.device.MotorDeviceSettings


if __name__ == '__main__':
    #A test function for the device. Tries to connect to the
    a = TDC001()
    a.is_connected