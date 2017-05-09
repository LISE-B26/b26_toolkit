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

from gauge_controller import PressureGauge, PumpLinePressureGauge, ChamberPressureGauge
from spectrum_analyzer import SpectrumAnalyzer
from ni_daq_v2 import NI6259, NI9263, NI9402
from piezo_controller import PiezoController
from zurich_instruments import ZIHF2
from pulse_blaster import B26PulseBlaster, Pulse
from maestro import MaestroLightControl
from attocube import Attocube
from microwave_generator import MicrowaveGenerator
from magnet_coils import MagnetCoils
from temperature_controller import TemperatureController
from montana import CryoStation
from labview_fpga import NI7845RMain
# from labview_fpga import FPGA_GalvoScan
from newport_smc100 import SMC100


# from labview_fpga import NI7845RReadWrite, NI7845RPidSimpleLoop, FPGA_GalvoScan, NI7845RReadFifo



# # old with tr/except
# # try:
# #     from Virtual_PI_Controler import PIControler
# # except:
# #     print("./src/instrument/__init__ warning! PIControler did not load")
#
# try:
#     from gauge_controller import PressureGauge
# except:
#     print("./src/instrument/__init__ warning! PressureGauge did not load")
#
# # try:
# #     from spectrum_analyzer import SpectrumAnalyzer
# # except:
# #     print("./src/instrument/__init__ warning! SpectrumAnalyzer did not load")
#
# try:
#     from ni_daq import DAQ
# except:
#     print("./src/instrument/__init__ warning! DAQ did not load")
#
# try:
#     from piezo_controller import PiezoController
# except:
#     print("./src/instrument/__init__ warning! PiezoController did not load")
#
# try:
#     from zurich_instruments import ZIHF2
# except:
#     print("./src/instrument/__init__ warning! ZIHF2 did not load")
#
# try:
#     from pulse_blaster import B26PulseBlaster, Pulse
# except:
#     print("./src/instrument/__init__ warning! PulseBlaster did not load")
#
# try:
#     from maestro import MaestroLightControl
# except:
#     print("./src/instrument/__init__ warning! MaestroLightControl did not load")
#
# try:
#     from attocube import Attocube
# except:
#     print("./src/instrument/__init__ warning! Attocube did not load")
#
#
# try:
#     from microwave_generator import MicrowaveGenerator
# except:
#     print("./src/instrument/__init__ warning! MicrowaveGenerator did not load")
#
# try:
#     from labview_fpga import NI7845RReadWrite
# except:
#     print("./src/instrument/__init__ warning! NI7845RReadAnalogIO did not load")
#
# try:
#     from labview_fpga import NI7845RPidSimpleLoop
# except:
#     print("./src/instrument/__init__ warning! NI7845RPidSimpleLoop did not load")
#
# try:
#     from labview_fpga import FPGA_GalvoScan
# except:
#     print("./src/instrument/__init__ warning! FPGA_GalvoScan did not load")
#
# try:
#     from labview_fpga import NI7845RReadFifo
# except:
#     print("./src/instrument/__init__ warning! NI7845RReadFifo did not load")
#
#
# # skip for now:
# # try:
# #     from spectrum_analyzer import SpectrumAnalyzer
# # except:
# #     print("./src/instrument/__init__ warning! SpectrumAnalyzer did not load")