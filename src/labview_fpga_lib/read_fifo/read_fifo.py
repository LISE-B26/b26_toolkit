"""
    This file is part of b26_toolkit, a PyLabControl add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    Foobar is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""

"""
Created on Wed Oct 22 17:46:08 2014

@author: Jan Gieseler


Wrapper for c-compiled FPGA_PID_Loop_Simple.vi (complied into FPGA_PID_lib.dll)
here we wrap the c-functions that can are then accessed by the higher level
FPGA_PID_Loop_Simple.py which defines the higher level Python objects that are then used to interact with the NI FPGA

"""

# TODO: reading of analog input gives only positive values (as if cast into unsigned integer, try to figure out why that is

from ctypes import *
import numpy as np
import os
from PyLabControl.src.core.read_write_functions import get_config_value
# =========================================================================
# ======= LOAD DLL ========================================================
# =========================================================================
dll_path = get_config_value('NIFPGA_DLL_PATH', os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'config.txt')) + 'read_fifo/read_fifo.dll'
_libfpga = WinDLL(dll_path)
# =========================================================================
# ======= DEFINE SETTER FUNCTIONS =========================================
# =========================================================================
# name of dictionary entry is the name of the function
# value of the dictionary entry is the data type that is passed to the function

setter_functions = {
    "set_ElementsToWrite": c_int32,
    "set_SamplePeriodsAcq": c_uint32,
    "set_AcquireData": c_bool,
    "set_Stop": c_bool
}

for fun_name in setter_functions:
    setattr( _libfpga, "{:s}.argtypes".format(fun_name), [setter_functions[fun_name], POINTER(c_uint32), POINTER(c_int32)])
    setattr( _libfpga, "{:s}.restype".format(fun_name), None)
    exec("""def {:s}(value, session, status):
        return _libfpga.{:s}(value, byref(session), byref(status))""".format(fun_name, fun_name))

# =========================================================================
# ======= DEFINE GETTER FUNCTIONS =========================================
# =========================================================================
# name of dictionary entry is the name of the function
# value of the dictionary entry is the data type that is returned from the function
getter_functions = {
    "start_fpga": None,
    "stop_fpga": None,
    "read_AI0": c_int16,
    "read_AI1": c_int16,
    "read_ElementsWritten": c_int32,
    "read_SamplePeriodsAcq": c_uint32,
    "read_LoopTicksAcq": c_uint32,
    "read_LoopRateLimitAcq": c_bool,
    "read_TimeOutAcq": c_bool,
    "read_FPGARunning": c_bool,
    "read_DMATimeOut": c_bool,
    "read_AcquireData": c_bool
}

for fun_name in getter_functions:
    setattr( _libfpga, "{:s}.argtypes".format(fun_name), [POINTER(c_uint32), POINTER(c_int32)])
    setattr( _libfpga, "{:s}.restype".format(fun_name), getter_functions[fun_name])
    exec("""def {:s}(session, status):
        return _libfpga.{:s}(byref(session), byref(status))""".format(fun_name, fun_name))


# =========================================================================
# ======= DEFINE FIFO FUNCTIONS =========================================
# =========================================================================
_libfpga.configure_FIFO_AI.argtypes = [c_uint32, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.configure_FIFO_AI.restype = c_uint32
def configure_FIFO_AI(requestedDepth, session, status):
    return _libfpga.configure_FIFO_AI(requestedDepth, byref(session), byref(status))


# _libfpga.set_FifoTimeout.argtypes = [c_int32, POINTER(c_uint32),
#                                      POINTER(c_int32)]
# _libfpga.set_FifoTimeout.restype = None


# start FIFO
_libfpga.start_FIFO_AI.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.start_FIFO_AI.restype = None
def start_FIFO_AI(session, status):
    return _libfpga.start_FIFO_AI(byref(session), byref(status))

# stop FIFO
_libfpga.stop_FIFO_AI.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.stop_FIFO_AI.restype = None
def stop_FIFO_AI(session, status):
    return _libfpga.stop_FIFO_AI(byref(session), byref(status))


# read FIFO
_libfpga.read_FIFO_AI.argtypes = [POINTER(c_uint32), c_int32,
                                  POINTER(c_uint32), POINTER(c_int32),
                                  POINTER(c_int32)]
_libfpga.read_FIFO_AI.restype = None

def read_FIFO_AI(size, session, status):
    AI1 = (c_int16*size)()
    AI2 = (c_int16*size)()
    elements_remaining = c_int32()

    _libfpga.read_FIFO_AI_unpack(AI1, AI2, size, byref(session), byref(status), byref(elements_remaining))


    return {'AI1': np.array(AI1), 'AI2': np.array(AI2), 'elements_remaining': elements_remaining.value}



# _libfpga.read_FIFO_AI_unpack.argtypes = [POINTER(c_int16), POINTER(c_int16),
#                                          c_int32, POINTER(c_uint32),
#                                          POINTER(c_int32), POINTER(c_int32)]
# _libfpga.read_FIFO_AI_unpack.restype = None

#
# def read_FIFO_conv(size, session, status, ticks=56):
#     """Reads a block of elements from the FPGA FIFO and determines the time
#     array corresponding to them.
#     """
#     set_LoopTicks(ticks, session, status)
#
#     [ai1, ai2, elements_remaining] = read_FIFO_AI(size, session, status)
#
#     if elements_remaining == size:
#         print("Warning: FIFO full and elements might get lost.")
#
#     # ai0 = int_to_voltage(array(list(ai0)))
#     # ai1 = int_to_voltage(array(list(ai1)))
#     # ai2 = int_to_voltage(array(list(ai2)))
#
#     # times = cumsum(array(list(ticks))) * 25e-9
#
#     return ai0, ai1, ai2, times





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