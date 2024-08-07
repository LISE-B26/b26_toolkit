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
import platform
if platform.release() != '7':
    from .zurich_instruments import ZIHF2, Hf2Li

from .gauge_controller import PressureGauge, PumpLinePressureGauge, ChamberPressureGauge
from .spectrum_analyzer import SpectrumAnalyzer
from .ni_daq import NI6259, NI9263, NI9402, NI9219, NI9263_02, NI9215
from .piezo_controller import PiezoController, PiezoControllerCold, MDT693A
from .pulse_blaster import B26PulseBlaster, Pulse, PulseBlasterHwTrig
from .maestro import MaestroLightControl
from .attocube import ANC300, ANC350
from .microwave_generator import MicrowaveGenerator, RFGenerator, MicrowaveGenerator2
from .magnet_coils import MagnetCoils
from .montana import CryoStation
from .newport_smc100 import SMC100
from .awg import AFG3022C, AFG3022C_02
from .oscilloscope import RigolOscilloscope
from .keysight_oscilloscope import Oscilloscope
from .thorlabs_kcube import KDC001, TLI_DeviceInfo, B26KDC001x, B26KDC001y, B26KDC001z, LockboxToggleArm
from .magnet_coils import MagnetCoils
from .ueye_camera import UEyeCamera
from .optotune_lens import OptotuneLens
from .arduino import ArduinoZero
from .commander import Commander
from .temperature_controller import TemperatureController, LakeShore211
from .red_laser import WlmMonitorSiV
from .moku_lab import MokuLockInAmplifier
from .afg import AFG3021C
