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

@author: Erik Hebestreit, Jan Gieseler


Wrapper for c-compiled FPGA_read_inputs.vi

"""

from ctypes import c_uint32, c_int32

import src.instruments.labview_fpga_lib.read_ai_ao.FPGAlib as FPGAlib


class NI7845R(object):
    session = c_uint32()
    status = c_int32()

    def __init__(self):
        pass

    def start(self):
        FPGAlib.start_fpga(self.session, self.status)
        return self.status

    def stop(self):
        FPGAlib.stop_fpga(self.session, self.status)

    # @property
    # def device_temperature(self):
    #     FPGAlib.read_DeviceTemperature(self.session, self.status)

class AnalogInput(object):
    _channel_number = None
    _fpga = None

    def __init__(self, channel, fpga):
        self._channel_number = channel
        self._fpga = fpga

    def read(self):
        return getattr(FPGAlib, 'read_AI%0d' % self._channel_number)(self._fpga.session, self._fpga.is_connected)


class AnalogOutput(object):
    _channel_number = None
    _fpga = None

    def __init__(self, channel, fpga):
        self._channel_number = channel
        self._fpga = fpga

    def write(self, value):
        return getattr(FPGAlib, 'set_AO%0d' % self._channel_number) \
            (value, self._fpga.session,
             self._fpga.is_connected)
