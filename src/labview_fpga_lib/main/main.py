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
    along with PyLabControl.  If not, see <http://www.gnu.org/licenses/>.
"""

from ctypes import *
import numpy as np
import os
from PyLabControl.src.core.read_write_functions import get_config_value
# TODO: find a way to call lib from a folder which doesn't contrain the bitfile of the FPGA (now it has to be place in the same directory as the python file to work

# =========================================================================
# ======= LOAD DLL ========================================================
# =========================================================================

try:
    _libfpga = WinDLL('C:/Users/Experiment/PycharmProjects/b26_toolkit/src/labview_fpga_lib/main/main.dll')
except:
    raise ImportError

# =========================================================================
# ======= DEFINE constants =========================================
# =========================================================================
IDLE = 0
READ_IO = 1
GALVO_SCAN = 2



# =========================================================================
# ======= DEFINE SETTER FUNCTIONS =========================================
# =========================================================================
# name of dictionary entry is the name of the function
# value of the dictionary entry is the data type that is passed to the function
setter_functions = {
    "set_Nx": c_int16,
    "set_Vmin_x": c_int16,
    "set_dVmin_x": c_int16,
    "set_Ny": c_int16,
    "set_Vmin_y": c_int16,
    "set_dVmin_y": c_int16,
    "set_scanmode_x": c_uint8,
    "set_scanmode_y": c_uint8,
    "set_detector_mode": c_uint8,
    "set_settle_time": c_uint16,
    "set_meas_per_pt": c_uint16,
    'set_stop_all': c_bool,
    "set_run_mode": c_uint16,
    "set_count_ms": c_uint32,
}
# analog channels
for i in range(8):
    setter_functions.update({'set_AO{:d}'.format(i): c_int16})
# digital channels
for i in range(4):
    setter_functions.update({'set_DIO{:d}'.format(i+4): c_bool})
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
    "reset_fpga": None,
    "read_Nx": c_int32,
    "read_Ny": c_int32,
    'read_run_mode':c_uint16,
    'read_stop_all':c_bool,
    'read_executing_subvi':c_bool,
    'read_meas_per_pt':c_uint16,
    'read_settle_time':c_uint16,
    "read_elements_written_to_dma": c_int32
}

# analog channels
# for i in range(8):
#     getter_functions.update({'read_AI{:d}'.format(i): c_int16})
#     getter_functions.update({'read_AO{:d}'.format(i): c_int16})
# digital channels
for i in range(4):
    getter_functions.update({'read_DIO{:d}'.format(i): c_bool})

for fun_name in getter_functions:
    setattr( _libfpga, "{:s}.argtypes".format(fun_name), [POINTER(c_uint32), POINTER(c_int32)])
    setattr( _libfpga, "{:s}.restype".format(fun_name), getter_functions[fun_name])
    exec ("""def {:s}(session, status):
        return _libfpga.{:s}(byref(session), byref(status))""".format(fun_name, fun_name))


getter_functions = {
}

# analog channels
for i in range(8):
    getter_functions.update({'read_AI{:d}'.format(i): c_int16})
    getter_functions.update({'read_AO{:d}'.format(i): c_int16})

for fun_name in getter_functions:
    setattr( _libfpga, "{:s}.argtypes".format(fun_name), [POINTER(c_uint32), POINTER(c_int32)])
    setattr( _libfpga, "{:s}.restype".format(fun_name), getter_functions[fun_name])
    # # force the analog signals into a c_int16 type (for some reason they seem to be c_uint16
    exec ("""def {:s}(session, status):
        value = _libfpga.{:s}(byref(session), byref(status))
        value = c_int16(value)
        value = value.value
        return value""".format(fun_name, fun_name))
    # exec ("""def {:s}(session, status):
    #     return _libfpga.{:s}(byref(session), byref(status))""".format(fun_name, fun_name))


# =========================================================================
# ======= DEFINE FIFO FUNCTIONS =========================================
# =========================================================================
_libfpga.configure_FIFO.argtypes = [c_uint32, POINTER(c_uint32), POINTER(c_int32)]
_libfpga.configure_FIFO.restype = c_uint32
def configure_FIFO(requestedDepth, session, status):
    return _libfpga.configure_FIFO(requestedDepth, byref(session), byref(status))

# start FIFO
_libfpga.start_FIFO.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.start_FIFO.restype = None
def start_FIFO(session, status):
    return _libfpga.start_FIFO(byref(session), byref(status))

# stop FIFO
_libfpga.stop_FIFO.argtypes = [POINTER(c_uint32), POINTER(c_int32)]
_libfpga.stop_FIFO.restype = None
def stop_FIFO(session, status):
    return _libfpga.stop_FIFO(byref(session), byref(status))


# read FIFO
_libfpga.read_FIFO.argtypes = [POINTER(c_int32), c_uint32,
                               POINTER(c_uint32), POINTER(c_int32),
                               POINTER(c_uint32)]

_libfpga.read_FIFO.restype = None

def read_FIFO(size, session, status):
    Signal = (c_int32*size)()
    elements_remaining = c_uint32()
    print(('ffffff reading fifo size', size))
    _libfpga.read_FIFO(Signal, size, byref(session), byref(status), byref(elements_remaining))

    print('ffffff reading fifo')
    return {'signal': np.array(Signal), 'elements_remaining': elements_remaining.value}

#
#
# # set values
# def set_scan_parameter(session, status,
#                        vmin_x, dv_x, n_x,
#                        vmin_y, dv_y, n_y,
#                        settle_time = 1000,
#                        loop_time = 10000,
#                        scan_mode = 'forward',
#                        forward_y = True,
#                        ):
#
#     _libfpga.settle_time_CountTicks(value, byref(session), byref(status))

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
        print(('fpga started, status = {:s}'.format(str(self.status.value))))
        # reset_fpga(self.session, self.status)
        # print('fpga reset, status = ', self.status.value)
        #
        # start_fpga(self.session, self.status)
        print(('fpga started, status = {:s}'.format(str(self.status.value))))
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