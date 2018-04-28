"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
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

# from ctypes import c_uint32, c_int32
# import time
# import old_lib.FPGAlib as FPGAlib
#
# session = c_uint32()
# status = c_int32()
# print('open')
# print(session, status)
# FPGAlib.start_fpga(session, status)
# time.sleep(1)
# for i in range(10):
#     x = FPGAlib.read_AI0(session, status)
#     print(x)
# print('close')
# FPGAlib.stop_fpga(session, status)
#
# print(session, status)
# #

from src.instruments.labview_fpga_lib.old import old_lib as NI

if __name__ == '__main__':
    fpga = NI.NI7845R()
    fpga.start()

    AI = NI.AnalogInput(0, fpga)
    for i in range(10):
        x = AI.read()
        print(x)
    fpga.stop()

