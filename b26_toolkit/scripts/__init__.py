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

# # from test_script import ScriptTest
from .galvo_scan.galvo_scan import GalvoScan, GalvoScanTimetrace
from .galvo_scan.galvo_scan_photodiode import GalvoScanPhotodiode
from .set_laser import SetLaser, SetAtto
from .find_nv import FindNV
from .daq_read_counter import Daq_Read_Counter
from .take_image_camera import TakeImage
from .esr import ESR, ESR_tracking, ESR_simple_lowerupper, ESR_simple
from .esr_dithering import ESR_FM_Dither
from .esr_two_freq_continuous import ESRTwoFreqContinuous
from .spec_analyzer_get_spectrum import SpecAnalyzerGetSpectrum
from .zi_sweeper import ZISweeper
from .zi_high_res_sweep import ZISweeperHighResolution
from .atto_scan import AttoStep
# # from .pulse_sequences import XY8_k, T1, Rabi, PDD, XY4, T1SingleInit, PulsedESR, \
# #     HahnEcho, XY4, XYXY, ReadoutStartTimeWithoutMW, ReadoutStartTime, ReadoutDuration, CPMG, \
# #     HahnEchoManyNVs, RabiPowerSweepSingleTau
from .pulse_sequences.rabi import Rabi
from .esr_and_rabi import ESRAndRabi
# from .spec_analyzer_get_spectrum import KeysightGetSpectrum
from .light_control import ApplyLightControlSettings, CameraOn
from .correlate_images import Track_Correlate_Images, Take_And_Correlate_Images
from .autofocus import AutoFocusDAQ, AutoFocusTwoPoints, AutoFocusTwoPointsFR, AutoFocusDaqSMC, AutoFocusCameraSMC, AutoFocusDAQCold, AutoFocusDAQMax
from .record_pressures import RecordPressures
from .set_magnetic_coils import SetMagneticCoils
from .align_magnetic_field_to_NV import AlignFieldToNV
from .Ni_9263_polarization_controller import Ni9263_BalancePolarization
from .stability_with_microwaves import Stability_With_Microwaves
from .read_temperature_lakeshore import ReadTemperatureLakeshore211
from .daq_read_counter_timetrace import Daq_TimeTrace_NI6259, Daq_TimeTrace_NI9402_NI9219
