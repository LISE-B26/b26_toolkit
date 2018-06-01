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

import os
from pylabcontrol.core.read_write_functions import get_config_value
from ctypes import *
# TODO: find a way to call lib from a folder which doesn't contrain the bitfile of the FPGA (now it has to be place in the same directory as the python file to work

dll_path = get_config_value('NIFPGA_DLL_PATH', os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'instruments\\config.txt')) + 'read_ai_ao/read_ai_ao.dll'
# print('dll_path', dll_path)
dll_path = 'C:\\Users\\Experiment\\PycharmProjects\\b26_toolkit\\pylabcontrol\\labview_fpga_lib\\read_ai_ao/read_ai_ao.dll'
_libfpga = WinDLL(dll_path)


_libfpga.start_fpga.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.start_fpga.restype = None

_libfpga.stop_fpga.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.stop_fpga.restype = None

# read inputs
_libfpga.read_AI0.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AI0.restype = c_int16

_libfpga.read_AI1.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AI1.restype = c_int16

_libfpga.read_AI2.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AI2.restype = c_int16

_libfpga.read_AI3.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AI3.restype = c_int16

_libfpga.read_AI4.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AI4.restype = c_int16

_libfpga.read_AI5.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AI5.restype = c_int16

_libfpga.read_AI6.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AI6.restype = c_int16

_libfpga.read_AI7.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AI7.restype = c_int16


_libfpga.read_AO0.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AO0.restype = c_int16

_libfpga.read_AO1.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AO1.restype = c_int16

_libfpga.read_AO2.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AO2.restype = c_int16

_libfpga.read_AO3.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AO3.restype = c_int16

_libfpga.read_AO4.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AO4.restype = c_int16

_libfpga.read_AO5.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AO5.restype = c_int16

_libfpga.read_AO6.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AO6.restype = c_int16

_libfpga.read_AO7.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_AO7.restype = c_int16

# set outputs
_libfpga.set_AO0.argtypes = [c_int16, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.set_AO0.restype = None

_libfpga.set_AO1.argtypes = [c_int16, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.set_AO1.restype = None

_libfpga.set_AO2.argtypes = [c_int16, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.set_AO2.restype = None

_libfpga.set_AO3.argtypes = [c_int16, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.set_AO3.restype = None

_libfpga.set_AO4.argtypes = [c_int16, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.set_AO4.restype = None

_libfpga.set_AO5.argtypes = [c_int16, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.set_AO5.restype = None

_libfpga.set_AO6.argtypes = [c_int16, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.set_AO6.restype = None

_libfpga.set_AO7.argtypes = [c_int16, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.set_AO7.restype = None


_libfpga.read_DIO0.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_DIO0.restype = c_bool

_libfpga.read_DIO1.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_DIO1.restype = c_bool

_libfpga.read_DIO2.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_DIO3.restype = c_bool

_libfpga.read_DIO3.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.read_DIO3.restype = c_bool


_libfpga.set_DIO4.argtypes = [c_bool, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.set_DIO4.restype = None

_libfpga.set_DIO5.argtypes = [c_bool, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.set_DIO5.restype = None

_libfpga.set_DIO6.argtypes = [c_bool, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.set_DIO6.restype = None

_libfpga.set_DIO7.argtypes = [c_bool, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.set_DIO7.restype = None


def start_fpga(session, status):
    return _libfpga.start_fpga(byref(session), byref(status))

def stop_fpga(session, status):
    return _libfpga.stop_fpga(byref(session), byref(status))


# read Analog inputs
def read_AI0(session, status):
    return _libfpga.read_AI0(byref(session), byref(status))

def read_AI1(session, status):
    return _libfpga.read_AI1(byref(session), byref(status))

def read_AI2(session, status):
    return _libfpga.read_AI2(byref(session), byref(status))

def read_AI3(session, status):
    return _libfpga.read_AI3(byref(session), byref(status))

def read_AI4(session, status):
    print('sadasdsad')
    return _libfpga.read_AI4(byref(session), byref(status))

def read_AI5(session, status):
    return _libfpga.read_AI5(byref(session), byref(status))

def read_AI6(session, status):
    return _libfpga.read_AI6(byref(session), byref(status))

def read_AI7(session, status):
    return _libfpga.read_AI7(byref(session), byref(status))


def read_AO0(session, status):
    return _libfpga.read_AO0(byref(session), byref(status))

def read_AO1(session, status):
    return _libfpga.read_AO1(byref(session), byref(status))

def read_AO2(session, status):
    return _libfpga.read_AO2(byref(session), byref(status))

def read_AO3(session, status):
    return _libfpga.read_AO3(byref(session), byref(status))

def read_AO4(session, status):
    return _libfpga.read_AO4(byref(session), byref(status))

def read_AO5(session, status):
    return _libfpga.read_AO5(byref(session), byref(status))

def read_AO6(session, status):
    return _libfpga.read_AO6(byref(session), byref(status))

def read_AO7(session, status):
    return _libfpga.read_AO7(byref(session), byref(status))

def set_AO0(value, session, status):
    return _libfpga.set_AO0(value, byref(session), byref(status))

def set_AO1(value, session, status):
    return _libfpga.set_AO1(value, byref(session), byref(status))

def set_AO2(value, session, status):
    return _libfpga.set_AO2(value, byref(session), byref(status))

def set_AO3(value, session, status):
    return _libfpga.set_AO3(value, byref(session), byref(status))

def set_AO4(value, session, status):
    return _libfpga.set_AO4(value, byref(session), byref(status))

def set_AO5(value, session, status):
    return _libfpga.set_AO5(value, byref(session), byref(status))

def set_AO6(value, session, status):
    return _libfpga.set_AO6(value, byref(session), byref(status))

def set_AO7(value, session, status):
    return _libfpga.set_AO7(value, byref(session), byref(status))

def read_DIO0(session, status):
    return _libfpga.read_DIO0(byref(session), byref(status))

def read_DIO1(session, status):
    return _libfpga.read_DIO1(byref(session), byref(status))

def read_DIO2(session, status):
    return _libfpga.read_DIO2(byref(session), byref(status))

def read_DIO3(session, status):
    return _libfpga.read_DIO3(byref(session), byref(status))

def set_DIO4(value, session, status):
    return _libfpga.set_DIO4(value, byref(session), byref(status))

def set_DIO5(value, session, status):
    return _libfpga.set_DIO5(value, byref(session), byref(status))

def set_DIO6(value, session, status):
    return _libfpga.set_DIO6(value, byref(session), byref(status))

def set_DIO7(value, session, status):
    return _libfpga.set_DIO7(value, byref(session), byref(status))

class NI7845R(object):
    session = c_uint32()
    status = c_int32()

    def __init__(self):
        """
        object to establish communication with the NI FPGA
        note this has to be implemented for each bitfile (i.e. labview fpga program)
        because the extual implementation of start in the .c code calls a different bitfile
        """
        pass


    def start(self):
        start_fpga(self.session, self.status)
        # print(self.status.value)
        if self.status.value != 0:
            if int(self.status.value) ==  -63101:
                print("ERROR 63101: Bitfile not found")
            else:
                print(('ERROR IN STARTING FPGA  (ERROR CODE: ', self.status.value, ')'))
        return self.status

    def stop(self):
        stop_fpga(self.session, self.status)

    # @property
    # def device_temperature(self):
    #     read_DeviceTemperature(self.session, self.status)



if __name__ == '__main__':

    fpga = NI7845R()
    fpga.start()
    fpga.stop()