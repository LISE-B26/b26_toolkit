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

import time

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import time
from b26_toolkit.instruments import NI9263, NI6259
from pylabcontrol.core import Parameter, Script
from b26_toolkit.data_processing.fit_functions import fit_cose_parameter, cose
from b26_toolkit.labview_fpga_lib.labview_helper_functions.labview_conversion import int_to_voltage, voltage_to_int
from b26_toolkit.instruments import MaestroLightControl
import datetime


class Ni9263_BalancePolarization(Script):
    """
Ni9263_BalancePolarization:
script to balance photodetector to zero by adjusting polarization controller voltages
 uses the Ni9263 as output and the NI DAQ as input
    """

    _DEFAULT_SETTINGS = [
        Parameter('channels', [
            Parameter('channel_WP_1', 'NI9263_ao0', ['NI9263_ao0', 'NI9263_ao1', 'NI9263_ao2'], 'analog channel that controls waveplate 1'),
            Parameter('channel_WP_2', 'NI9263_ao1', ['NI9263_ao0', 'NI9263_ao1', 'NI9263_ao2'], 'analog channel that controls waveplate 2'),
            Parameter('channel_WP_3', 'NI9263_ao2', ['NI9263_ao0', 'NI9263_ao1', 'NI9263_ao2'], 'analog channel that controls waveplate 3'),
            # Parameter('channel_OnOff', 'NI6259_do8', ['NI6259_do8'], 'digital channel that turns polarization controller on/off'),
            Parameter('channel_OnOff', 'NI9263_ao3', ['NI9263_ao3'], 'analog channel that turns polarization controller on/off'),
            Parameter('channel_detector', 'NI6259_ai0', ['NI6259_ai0'], 'analog input channel of the detector signal')
        ]),
        Parameter('setpoints', [
            Parameter('V_1', 2.4, float, 'voltage applied to waveplate 1'),
            Parameter('V_2', 4.0, float, 'voltage applied to waveplate 2'),
            Parameter('V_3', 2.4, float, 'voltage applied to waveplate 3'),
            Parameter('save_result_as_setpoint', False, bool, 'uses the current best result as the new setpoint')
        ]),
        Parameter('optimization',[
            Parameter('target', .01, float, 'target max detector signal'),
            Parameter('settle_time', 2., float, 'settle time (s)'),
            Parameter('WP_control', 2, [1, 2, 3], 'control waveplate'),
            Parameter('dV', 0.01, float, 'initial step size of search algorithm'),
            Parameter('slope', 'negative', ['positive', 'negative'], 'is the slope around the zero crossing is positive or negative')
            # Parameter('start with current', True, bool, 'uses the current output as starting point (instead of setpoint) if not zero'),
        ]),
        Parameter('measure_at_zero',[
            Parameter('on', True, bool, 'if true keep measuring after zero is found'),
            Parameter('N', 10, int, 'number of measurement points after zero is found')
        ])
    ]

    _INSTRUMENTS = {'NI6259': NI6259,
                    'NI9263': NI9263
                    }

    _SCRIPTS = {}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        def get_direction(detector_value, slope):
            """
            give the direction we have to move given the detector value and the slope of the zero crossing
            Args:
                detector_value: detector value
                slope: slope of the zero crossing

            Returns: direction in which to step


            we calculate the directin based on the following truth table (i.e. xor table)

                slope   detector    direction
                +       -           +
                +       +           -
                -       +           +
                -       -           -


            """
            direction = (int(np.sign(detector_value)) == 1) ^ (int(slope) == 1)
            # now map True and False and 1 and -1
            direction = 1 if direction else -1

            return direction


        self.data = {'voltage_waveplate': [], 'detector_signal': [], 'det_signal_cont': []}
        wp_control = self.settings['optimization']['WP_control']
        settle_time = self.settings['optimization']['settle_time']
        v_out = float(self.settings['setpoints']['V_{:d}'.format(wp_control)])
        target = self.settings['optimization']['target']
        detector_value = 2* target

        # convert slope to numeric value
        slope = 1 if self.settings['optimization']['slope'] == 'positive' else -1

        # get the channels
        print((self.settings['channels']['channel_OnOff']))
        control_channel = self.settings['channels']['channel_OnOff'].split('NI9263_')[1]
        channel_out = self.settings['channels']['channel_WP_{:d}'.format(wp_control)].split('NI9263_')[1]
        channel_in = self.settings['channels']['channel_detector'].split('NI6259_')[1]



        # NIDAQ_DIO = self.instruments['NI6259']['instance'] # digital channel
        NIDAQ_AI = self.instruments['NI6259']['instance']  # analog input
        NIDAQ_AO = self.instruments['NI9263']['instance'] # analog output

        # turn controller on
        # NIDAQ_DIO.set_digital_output({control_channel: True}) #switched to analog output since digital output is
                                                                #slightly undervolted and threshold is at 5V
        NIDAQ_AO.set_analog_voltages({control_channel: 5.1})

        # fpga_io.update({'read_io':{control_channel: True}})


        # set the setpoints for all three waveplates
        dictator = {}
        for i in [1, 2, 3]:
            # readout of daq output voltages not avaliable on Ni_9263, no analog input
            # if self.settings['optimization']['start with current'] and i == wp_control:
            #     # value = getattr(fpga_io, channel_out)
            #     value = NIDAQ_AO[channel_out]
            #     if value == 0:
            #         value = float(self.settings['setpoints']['V_{:d}'.format(i)])
            #         self.log('current value is zero take setpoint {:0.3f} V as starting point'.format(value))
            #     else:
            #         # value = int_to_voltage(value)
            #         v_out = value
            #         self.log('use current value {:0.3f} V as starting point'.format(value))
            # else:
            #     value = float(self.settings['setpoints']['V_{:d}'.format(i)])
            value = float(self.settings['setpoints']['V_{:d}'.format(i)])

            daq_wp_channel = self.settings['channels']['channel_WP_{:d}'.format(i)].split('NI9263_')[1]
            dictator.update({daq_wp_channel: value})

        NIDAQ_AO.set_analog_voltages(dictator)
        time.sleep(settle_time)

        crossed_zero = False
        while abs(detector_value) > abs(target):
            if self._abort:
                break

            # set output
            # NIDAQ_AO.update({channel_out: float(v_out)})
            NIDAQ_AO.set_analog_voltages({channel_out: float(v_out)})
            # fpga_io.update({'read_io':{channel_out: float(v_out)}})
            # wait for system to settle
            time.sleep(settle_time)
            # read detector
            # detector_value = NIDAQ_DIO[channel_in]
            detector_value = NIDAQ_AI.get_analog_voltages([channel_in])[0]

            self.data['voltage_waveplate'].append(v_out)
            self.data['detector_signal'].append(detector_value)
            self.progress = 50.
            self.updateProgress.emit(self.progress)

            direction = get_direction(detector_value, slope)
            print(('----> out', v_out, voltage_to_int((v_out))))
            print(('----> det', detector_value, slope, direction))
            # calculate the next step
            if len(self.data['voltage_waveplate']) ==1:
                v_step = self.settings['optimization']['dV'] # start with initial step size
            elif len(self.data['voltage_waveplate']) > 1:

                # check for zero crossing
                if self.data['detector_signal'][-2] * self.data['detector_signal'][-1] < 0:
                    self.log('detected zero crossing!')
                    v_step /=2 # decrease the step size since we are closer to zero

            # calculate next output voltage
            v_out += v_step * direction


            if v_out > 5 or v_out <0:
                v_out = 5 if v_out > 5 else 0 # set the output to be within the range
                slope *= -1

            if min(self.data['voltage_waveplate']) == 0 and max(self.data['voltage_waveplate']) ==5:
                self.log('warning! scanned full range without finding zero. abort!')
                self._abort = True


        self.data['setpoint'] = v_out
        print('starting continuous measurement')

        n_opt = len(self.data['voltage_waveplate'])
        n_cont = self.settings['measure_at_zero']['N']
        if self.settings['measure_at_zero']['on']:
            for i in range(n_cont):
                if self._abort:
                    break
                detector_value = NIDAQ_AI.get_analog_voltages([channel_in])
                self.data['det_signal_cont'].append(detector_value)
                self.progress = 100.* (i+n_opt+1) /(n_cont+ n_opt)

                self.updateProgress.emit(self.progress)
                time.sleep(settle_time)

        if self.settings['setpoints']['save_result_as_setpoint']:
            self.settings['setpoints']['V_{:d}'.format(self.settings['optimization']['WP_control'])] = v_out


    def _plot(self, axes_list):

        if self.data != {}:

            # dt = self.settings['optimization']['settle_time']
            # N = len(self.data['voltage_waveplate'])

            # t = np.linspace(0,N * dt, N)
            volt_range = np.array(self.data['voltage_waveplate'])
            signal = np.array(self.data['detector_signal'])

            if len(self.data['det_signal_cont']) == 0:
                axes_list[0].plot(signal, '-o')
                axes_list[0].hold(False)
                axes_list[0].set_ylabel('detector signal (bits)')

                axes_list[1].plot(volt_range, '-o')
                axes_list[1].hold(False)
                axes_list[0].set_ylabel('wp voltage (V)')
            else:
                axes_list[1].plot(signal/float(2**15), '-o', label = 'detector signal / (2^15)')
                axes_list[1].hold(True)
                axes_list[1].plot(volt_range/5., '-o', label = 'wp voltage / 5V')
                axes_list[1].hold(False)
                axes_list[1].set_ylabel('signal')
                axes_list[1].legend(fontsize=8)

                signal_det = np.array(self.data['det_signal_cont'])
                axes_list[0].plot(signal_det, '-o')
                axes_list[0].set_ylabel('detector signal (bit)')
                axes_list[0].hold(False)

    def _update(self, axes_list):
        volt_range = np.array(self.data['voltage_waveplate'])
        signal = np.array(self.data['detector_signal'])
        signal_det = np.array(self.data['det_signal_cont'])

        if len(self.data['det_signal_cont']) == 1:
            self._plot(self, axes_list)
        elif len(self.data['det_signal_cont']) == 0:

            axes_list[0].lines[0].set_ydata(signal)
            axes_list[1].lines[0].set_ydata(volt_range)
        else:
            axes_list[0].lines[0].set_ydata(signal_det)
            # axes_list[1].lines[0].set_ydata(signal/float(2**15))
            # axes_list[1].lines[1].set_ydata(volt_range/5.)


if __name__ == '__main__':
    pass
    # from pylabcontrol.core import Instrument, Parameter
    #
    # daq, failed = Instrument.load_and_append({'daq': NI9263})
    # print(daq)
    #
    # daq['daq'].ao0 = 0



