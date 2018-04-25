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

from .gauge_controller import PressureGauge, PumpLinePressureGauge, ChamberPressureGauge
from .spectrum_analyzer import SpectrumAnalyzer
from .ni_daq import NI6259, NI9263, NI9402
from .piezo_controller import PiezoController
from .zurich_instruments import ZIHF2
from .pulse_blaster import B26PulseBlaster, Pulse
from .maestro import MaestroLightControl
from .attocube import Attocube
from .microwave_generator import MicrowaveGenerator
from .magnet_coils import MagnetCoils
from .temperature_controller import TemperatureController
from .montana import CryoStation
from .newport_smc100 import SMC100
from .awg import AWG
from .keysight_oscilloscope import Oscilloscope
from .dcservo_kinesis_dll import KDC001
from .magnet_coils import MagnetCoils
from .ueye_camera import UEyeCamera