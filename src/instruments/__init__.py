from gauge_controller import PressureGauge
# from spectrum_analyzer import SpectrumAnalyzer
from ni_daq import DAQ
from piezo_controller import PiezoController
from zurich_instruments import ZIHF2
from pulse_blaster import B26PulseBlaster, Pulse
from maestro import MaestroLightControl

from attocube import Attocube

from microwave_generator import MicrowaveGenerator
# from labview_fpga import NI7845RReadWrite, NI7845RPidSimpleLoop, NI7845RGalvoScan, NI7845RReadFifo



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
#     from labview_fpga import NI7845RGalvoScan
# except:
#     print("./src/instrument/__init__ warning! NI7845RGalvoScan did not load")
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