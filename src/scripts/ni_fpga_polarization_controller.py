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

import time

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import time
from b26_toolkit.src.instruments.labview_fpga import NI7845RReadWrite
from PyLabControl.src.core import Parameter, Script
from PyLabControl.src.instruments import PIControler
from b26_toolkit.src.labview_fpga_lib.labview_helper_functions.labview_conversion import int_to_voltage
import datetime

class FPGA_PolarizationController(Script):
    """
script to balance photodetector to zero by adjusting polarization controller voltages
    """
    _DEFAULT_SETTINGS = [
        Parameter('channels', [
            Parameter('channel_WP_1', 5, range(8), 'analog channel that controls waveplate 1'),
            Parameter('channel_WP_2', 6, range(8), 'analog channel that controls waveplate 2'),
            Parameter('channel_WP_3', 7, range(8), 'analog channel that controls waveplate 3'),
            Parameter('channel_OnOff', 4, [4, 5, 6, 7], 'digital channel that turns polarization controller on/off'),
            Parameter('channel_detector', 0, range(4), 'analog input channel of the detector signal')
            ]),
        Parameter('setpoints',[
        Parameter('V_1', 2.0, float, 'voltage applied to waveplate 1'),
        Parameter('V_2', 1.0, float, 'voltage applied to waveplate 2'),
        Parameter('V_3', 2.0, float, 'voltage applied to waveplate 3')
        ]),
        Parameter('WP_FB', 2, [1,2,3] , 'waveplate that is used in feedback'),
        Parameter('FB_on', True, bool, 'activate feedback'),
        Parameter('stop_condition',[
            Parameter('threshold', 100., float, 'threshold that determines if detector is balanced'),
            Parameter('length', 10, int, 'number of data points that are used to calculate the stop condition'),
            Parameter('exit_mode', 'rms', ['continuous', 'rms', 'min-max'], 'mode how to determine if zero has been reached')
        ])
    ]

    _INSTRUMENTS = {
        'FPGA_IO': NI7845RReadWrite,
        'controler': PIControler
    }

    _SCRIPTS = {

    }

    def __init__(self, instruments, scripts = None, name=None, settings=None, log_function=None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path = data_path)
    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        # allocate data
        self.data = {'WP_volt': [], 'det_signal': []}

        # get instrument references
        fpga_io = self.instruments['FPGA_IO']['instance']
        controler = self.instruments['controler']['instance']

        controler.update(self.instruments['controler']['settings'])

        # get settings
        time_step = self.instruments['controler']['settings']['time_step']
        threshold= self.settings['stop_condition']['threshold']
        exit_mode = self.settings['stop_condition']['exit_mode']
        exit_length = self.settings['stop_condition']['length']

        # turn  polarization controller on
        control_channel = 'DIO{:d}'.format(self.settings['channels']['channel_OnOff'])
        instrument_settings = {control_channel: True}

        # set the voltages
        for c in [1,2,3]:
            channel_out = 'AO{:d}'.format(self.settings['channels']['channel_WP_{:d}'.format(c)])
            working_point = float(self.settings['setpoints']['V_{:d}'.format(c)])
            instrument_settings.update({channel_out:working_point})

        fpga_io.update(instrument_settings)
        self.log('wait {:0.2f} seconds to settle down'.format(time_step))
        time.sleep(time_step)


        channel_out = 'AO{:d}'.format(self.settings['channels']['channel_WP_{:d}'.format(self.settings['WP_FB'])])
        channel_in = 'AI{:d}'.format(self.settings['channels']['channel_detector'])

        working_point = float(self.settings['setpoints']['V_{:d}'.format(self.settings['WP_FB'])])

        exit_value = 2*threshold
        # exit condition is that the exit value falls below the threshold value
        while exit_value >= threshold:
            if self._abort:
                break

            # get the detector value and calculate the control value
            # use this syntax (getattr) because the channel is variable, i.e. we don't know if fpga_io.AI0, fpga_io.AI1 or fpga_io.AI2
            detector_value = getattr(fpga_io, channel_in)
            control_value = controler.controler_output(detector_value)

            control_value = int_to_voltage(control_value)

            print('control_value: ', control_value, exit_value, detector_value)
            if self.settings['FB_on']:
                fpga_io.update({channel_out: control_value+ working_point})

            self.data['WP_volt'].append(control_value )
            self.data['det_signal'].append(detector_value)

            time.sleep(time_step)

            self.progress = 10
            self.updateProgress.emit(self.progress)

            if len(self.data['det_signal']) < exit_length:
                continue
            # calculate if exit condition is met
            exit_signal = self.data['det_signal'][-exit_length:]
            if exit_mode == 'rms':
                exit_value = np.mean(np.abs(exit_signal))
            elif exit_mode == 'min-max':
                exit_value = max(exit_signal) - min(exit_signal)

        self.progress = 90
        self.updateProgress.emit(self.progress)



    def _plot(self, axes_list):
        axes1, axes2 = axes_list
        axes1.hold(False)
        axes2.hold(False)
        dt = self.instruments['controler']['settings']['time_step']
        N = len(self.data['det_signal'])
        t = dt* np.arange(N)
        axes1.plot(t, self.data['WP_volt'][:N], '-o')
        axes2.plot(t, self.data['det_signal'][:N], '-o')

        axes1.set_xlabel('time (s)')
        axes1.set_ylabel('wp signal')
        axes2.set_xlabel('frequency (Hz)')
        axes2.set_ylabel('detector signal')

        # axes2.plot(self.data['WP_volt'][0:-1], np.diff(self.data['det_signal']), '-o')

class FPGA_PolarizationSignalScan(Script):
    """
script to map out detector response as a function of polarization controller voltage WP2
the script scans the voltage of  channel 2 from 0 to 5 volt and records the detector response
this gives a one dimensional dataset
    """

    _DEFAULT_SETTINGS = [
        Parameter('channels', [
            Parameter('channel_WP_1', 5, range(8), 'analog channel that controls waveplate 1'),
            Parameter('channel_WP_2', 6, range(8), 'analog channel that controls waveplate 2'),
            Parameter('channel_WP_3', 7, range(8), 'analog channel that controls waveplate 3'),
            Parameter('channel_OnOff', 4, [4, 5, 6, 7], 'digital channel that turns polarization controller on/off'),
            Parameter('channel_detector', 0, range(4), 'analog input channel of the detector signal')
            ]),
        Parameter('setpoints',[
            Parameter('V_1', 2.4, float, 'voltage applied to waveplate 1'),
            Parameter('V_2', 4.0, float, 'voltage applied to waveplate 2'),
            Parameter('V_3', 2.4, float, 'voltage applied to waveplate 3')
        ]),
        Parameter('scan',[
            Parameter('WP_scan', 2, [1,2,3] , 'waveplate that is used in scan'),
            Parameter('settle_time', 1.0, float, 'settle time in seconds between voltage steps'),
            Parameter('scan_range', 0.5, float, 'voltage scan range around setpoint (output is capped at 0 and 5V)'),
            Parameter('N', 10, int, 'number of points in scan range')
        ]),
        Parameter('stop_at_zero', True, bool, 'if true scan stop once detector signal crosses 0')
    ]

    _INSTRUMENTS = {
        'FPGA_IO': NI7845RReadWrite
    }

    _SCRIPTS = {

    }

    def __init__(self, instruments, scripts = None, name=None, settings=None, log_function=None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path = data_path)


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        # def calc_progress(v2, volt_range):
        #     dV = np.mean(np.diff(volt_range))
        #     progress = (v2 - min(volt_range)) / (max(volt_range) - min(volt_range))
        #     return 100.*progress

        self.data = {}
        fpga_io = self.instruments['FPGA_IO']['instance']
        # fpga_io.update(self.instruments['FPGA_IO']['settings'])

        # turn controller on
        control_channel = 'DIO{:d}'.format(self.settings['channels']['channel_OnOff'])
        instrument_settings = {control_channel: True}
        fpga_io.update(instrument_settings)

        channel_out = 'AO{:d}'.format(self.settings['scan']['WP_scan'])

        channel_in = 'AI{:d}'.format(self.settings['channels']['channel_detector'])

        dictator = {}
        for i in [1,2,3]:
            channel_out = 'AO{:d}'.format(self.settings['channels']['channel_WP_{:d}'.format(i)])
            dictator.update({channel_out: float(self.settings['setpoints']['V_{:d}'.format(i)])})
        fpga_io.update(dictator)

        settle_time = self.settings['scan']['settle_time']

        time.sleep(settle_time)

        self.data = {'voltage_waveplate':[], 'detector_signal':[], 'det_signal_cont':[]}

        Vo = float(self.settings['setpoints']['V_{:d}'.format(self.settings['scan']['WP_scan'])])
        dV = float( self.settings['scan']['scan_range'])
        N = self.settings['scan']['N']


        volt_range = np.linspace(max(Vo - dV / 2, 0), min(Vo + dV / 2, 5), N)
        crossed_zero = False
        for i, v2 in enumerate(volt_range):
            if self._abort:
                break

            fpga_io.update({channel_out: float(v2)})
            time.sleep(settle_time)
            detector_value = getattr(fpga_io, channel_in)

            self.progress = 100.* i / len(volt_range)

            self.data['voltage_waveplate'].append(v2)
            self.data['detector_signal'].append(detector_value)

            self.updateProgress.emit(self.progress)

            if i > 2 and self.data['detector_signal'][-2] * self.data['detector_signal'][-1] < 0:
                crossed_zero = True
                self.log('detected zero crossing!')
            if self.settings['stop_at_zero'] and crossed_zero:
                break


        setpoint = float(self.calc_setpoint())
        self.data['setpoint'] = setpoint
        self.data['crossed_zero'] = crossed_zero

    def calc_setpoint(self):

        def func(x, a, b):
            return a * x + b

        volt_range = self.data['voltage_waveplate']
        signal = self.data['detector_signal']
        popt, pcov = curve_fit(func, volt_range, signal)

        setpoint = -popt[1] / popt[0]

        return setpoint
    def _plot(self, axes_list):

        if self.data != {}:
            axes = axes_list[0]

            volt_range = self.data['voltage_waveplate']
            signal = self.data['detector_signal']

            signal_cont = self.data['det_signal_cont']


            axes.plot(volt_range, signal, '-o')
            axes.hold(True)

            if len(volt_range)>5:

                # popt, pcov = curve_fit(func, volt_range, signal)
                # axes.plot(volt_range, func(volt_range, popt[0], popt[1]), 'k-')
                setpoint = self.calc_setpoint()
                axes.plot([setpoint, setpoint],[min(signal), max(signal)], 'r--')
                axes.set_title('setpoint = {:0.2f}V'.format(setpoint))
                axes.plot([min(volt_range), max(volt_range)], [0, 0], 'k-')

            axes.hold(False)

            if len(signal_cont)>0:
                axes = axes_list[1]

                axes.plot(signal_cont, '-o')
                axes.hold(False)



class FPGA_BalancePolarization(Script):
    """
    script to bring the detector response to zero
    two channels are set to a fixed voltage while the signal of the third channel is varied until the detector response is zero
    """

    _DEFAULT_SETTINGS = [
        Parameter('goto_zero', True, bool,'true output is set to the calculated setpoint such that expected detector signal is zero'),
        Parameter('read_signal', True, bool, 'keeps script running and just measures the signal'),
        Parameter('max_duration', 100, float, 'maximum time the script runs after finding the setpoint (s)'),
        Parameter('scan_range_target', 0.1, float, 'target scan range'),
        Parameter('settle_time', 1.0, float, 'settle time in seconds between voltage steps'),
    ]

    _INSTRUMENTS = {
    }

    _SCRIPTS = {
        'find_setpoint': FPGA_PolarizationSignalScan
    }

    def __init__(self, instruments = None, scripts = None, name=None, settings=None, log_function=None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """

        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path = data_path)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        self.data['det_signal'] = []
        find_setpoint = self.scripts['find_setpoint']


        dt = self.settings['settle_time']
        dV = float(find_setpoint.settings['scan']['scan_range']) / 2
        dV_target = self.settings['scan_range_target']

        WP_scan = find_setpoint.settings['scan']['WP_scan']
        channel_out = 'AO{:d}'.format(WP_scan)
        channel_in = 'AI{:d}'.format(find_setpoint.settings['channels']['channel_detector'])
        Vo = None

        fpga_io = find_setpoint.instruments['FPGA_IO']['instance']

        settle_time = self.settings['settle_time']

        while dV>dV_target:
            find_setpoint.run()

            if find_setpoint.data['crossed_zero'] is False:
                self.log('warning, couldn\'t find zero crossing - exit scan')
                break
            else:
                Vo = float(find_setpoint.data['setpoint'])
                dV /= 2

                find_setpoint.settings.update({
                    'setpoint': {'V_{:d}'.format(WP_scan):Vo},
                    'scan' : {'scan_range': dV}
                })



        if self.settings['goto_zero'] and Vo is not None:
            fpga_io.update({channel_out: Vo})

        start_time = datetime.datetime.now()

        while self.settings['read_signal']:
            if self._abort:
                break
            detector_value = getattr(fpga_io, channel_in)
            print('detector_value', detector_value)
            self.data['det_signal'].append(detector_value)
            time.sleep(settle_time)
            self.progress = 50
            self.updateProgress.emit(self.progress)
            if (datetime.datetime.now()- start_time).total_seconds()>=self.settings['max_duration']:
                break


    def plot(self, axes_list):
        if self._current_subscript_stage['current_subscript'].name == 'find_setpoint':
            self.scripts['find_setpoint'].plot(axes_list)
        else:
            ax = axes_list[1]

            if len(self.data['det_signal'] )>0:
                ax.plot(self.data['det_signal'])

    # def _plot(self, axes_list):
    #
    #     print('asdad', self.scripts['find setpoint'].is_running, self._current_subscript_stage)
    #     if self.scripts['find setpoint'].is_running:
    #         self.scripts['find setpoint']._plot(axes_list)
    #     else:
    #         self.scripts['polarization_control']._plot(axes_list)




class FPGA_PolarizationSignalMap(Script):
    """
script to map out detector response as a function of polarization controller voltages
the script scans the voltage of all three channels from 0 to 5 volt and records the detector response
this gives a three dimensional dataset
    """

    _DEFAULT_SETTINGS = [
        Parameter('WP_1', [
            Parameter('channel', 5, range(8), 'analog channel that controls waveplate 1'),
            Parameter('V_min', 0.0, float, 'minimum voltage of scan'),
            Parameter('V_max', 5.0, float, 'maximum voltage of scan'),
            Parameter('dV', 0.2, float, 'voltage step of scan')
        ]),
        Parameter('WP_2', [
            Parameter('channel', 6, range(8), 'analog channel that controls waveplate 1'),
            Parameter('V_min', 0.0, float, 'minimum voltage of scan'),
            Parameter('V_max', 5.0, float, 'maximum voltage of scan'),
            Parameter('dV', 0.2, float, 'voltage step of scan')
        ]),
        Parameter('WP_3', [
            Parameter('channel', 7, range(8), 'analog channel that controls waveplate 1'),
            Parameter('V_min', 0.0, float, 'minimum voltage of scan'),
            Parameter('V_max', 5.0, float, 'maximum voltage of scan'),
            Parameter('dV', 0.2, float, 'voltage step of scan')
        ]),
        # Parameter('channel_WP_1', 5, range(8), 'analog channel that controls waveplate 1'),
        # Parameter('channel_WP_2', 6, range(8), 'analog channel that controls waveplate 2'),
        # Parameter('channel_WP_3', 7, range(8), 'analog channel that controls waveplate 3'),
        Parameter('channel_OnOff', 4, [4,5,6,7], 'digital channel that turns polarization controller on/off'),
        Parameter('channel_detector', 0, range(4), 'analog input channel of the detector signal')
    ]

    _INSTRUMENTS = {
        'FPGA_IO': NI7845RReadWrite
    }

    _SCRIPTS = {

    }

    def __init__(self, instruments, scripts = None, name=None, settings=None, log_function=None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path = data_path)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """


        def calc_progress(v1, v3, volt_range_1, volt_range_3):
            dV = np.mean(np.diff(volt_range_1))
            progress_v3 = (v3 - min(volt_range_3)) / (max(volt_range_3) - min(volt_range_3))
            progress_v1 = (v1 - min(volt_range_1)) / (max(volt_range_1) - min(volt_range_1))
            progress = progress_v1 + progress_v3 * dV / (max(volt_range_1) - min(volt_range_1))
            return int(100*progress)

        self.data = {}
        fpga_io = self.instruments['FPGA_IO']['instance']
        # fpga_io.update(self.instruments['FPGA_IO']['settings'])

        # turn controller on
        control_channel = 'DIO{:d}'.format(self.settings['channel_OnOff'])
        instrument_settings = {control_channel: True}
        fpga_io.update(instrument_settings)

        # channel_out_1 = 'AO{:d}'.format(self.settings['channel_WP_{:d}'.format(1)])
        # channel_out_2 = 'AO{:d}'.format(self.settings['channel_WP_{:d}'.format(2)])
        # channel_out_3 = 'AO{:d}'.format(self.settings['channel_WP_{:d}'.format(3)])
        channel_out_1 = 'AO{:d}'.format(self.settings['WP_1']['channel'])
        channel_out_2 = 'AO{:d}'.format(self.settings['WP_2']['channel'])
        channel_out_3 = 'AO{:d}'.format(self.settings['WP_3']['channel'])
        channel_in = 'AI{:d}'.format(self.settings['channel_detector'])

        # set the voltages
        volt_range_1 = np.arange(self.settings['WP_1']['V_min'], self.settings['WP_1']['V_max'], self.settings['WP_1']['dV'])
        volt_range_2 = np.arange(self.settings['WP_2']['V_min'], self.settings['WP_2']['V_max'], self.settings['WP_2']['dV'])
        volt_range_3 = np.arange(self.settings['WP_3']['V_min'], self.settings['WP_3']['V_max'], self.settings['WP_3']['dV'])
        for v1 in volt_range_1:
            signal_1 = float(v1)
            for v3 in volt_range_3:
                signal_3 = float(v3)
                self.log('WP1 = {:0.2f}V, WP3 = {:0.2f}V: wait 3 seconds to settle down'.format(signal_1, signal_3))
                fpga_io.update({channel_out_1: signal_1, channel_out_2: 0, channel_out_3: signal_3})
                time.sleep(3)
                data = []
                for v2 in volt_range_2:
                    signal_2 = float(v2)
                    fpga_io.update({channel_out_2: signal_2})
                    time.sleep(1)
                    detector_value = getattr(fpga_io, channel_in)
                    data.append(detector_value)
                self.data.update({'WP1 = {:0.2f}V, WP3 = {:0.2f}V'.format(signal_1, signal_3) : data})
                progress = calc_progress(v1, v3, volt_range_1, volt_range_3)
                self.updateProgress.emit(progress)



    @staticmethod
    def data_to_pandas_panel(data):
        """
        casts the data in dictionary form into a pandas panel (3D dataset)
        Args:
            data: dictionary produced by the script

        Returns: pandas panel (3D dataset)

        """

        def get_voltages_from_tag(s):
            v1 = float(s.split(',')[0].split('WP1 = ')[1].split('V')[0])
            v3 = float(s.split(',')[1].split(' WP3 = ')[1].split('V')[0])
            return v1, v3

        v1_array = sorted(list(set([get_voltages_from_tag(k)[0] for k, value in data.iteritems()])))
        v3_array = sorted(list(set([get_voltages_from_tag(k)[1] for k, value in data.iteritems()])))
        v2_array = sorted(list(set([get_voltages_from_tag(k)[1] for k, value in data.iteritems()])))

        def get_img_df(v3):
            df = pd.DataFrame(
                [data['WP1 = {:0.2f}V, WP3 = {:0.2f}V'.format(v1, v3)][0:len(v2_array)] for v1 in v1_array],
                index=['{:0.2f}V'.format(v1) for v1 in v1_array],
                columns=['{:0.2f}V'.format(v2) for v2 in v2_array])
            return df

        panel = pd.Panel.from_dict({'{:0.2f}V'.format(v3): get_img_df(v3) for v3 in v3_array})
        return panel.transpose(1, 2, 0)

    @staticmethod
    def plot_pol_map(value, axis, data_panel, plot_axes = None):
        """
        plots the polarization map
        Args:
            value: value of voltage for which to plot 2D map
            axis: axis along which to slice data (1,2,3)
            data_panel: pandas data panel as obtained from data_to_pandas_panel()
            plot_axes: axes on which to plot data
        Returns:

        """
        if isinstance(value, float):
            value = '{:0.2f}V'.format(value)

        v1_array, v2_array, v3_array = data_panel.axes

        if axis == 1:
            data = data_panel[value]
            #         extent=[0, len(v2_array),len(v3_array),0]
            xlabel = 'waveplate 2 (V)'
            ylabel = 'waveplate 3 (V)'
            title = 'waveplate 1 = {:s}'.format(value)
            extent = [float(min(v2_array).split('V')[0]), float(max(v2_array).split('V')[0]),
                      float(max(v3_array).split('V')[0]), float(min(v3_array).split('V')[0])]
        elif axis == 2:
            data = data_panel.major_xs(value)
            #         extent=[0, len(v1_array),len(v3_array),0]
            extent = [float(min(v1_array).split('V')[0]), float(max(v1_array).split('V')[0]),
                      float(max(v3_array).split('V')[0]), float(min(v3_array).split('V')[0])]
            xlabel = 'waveplate 1 (V)'
            ylabel = 'waveplate 3 (V)'
            title = 'waveplate 2 = {:s}'.format(value)
        elif axis == 3:
            data = data_panel.minor_xs(value)
            #         extent=[0, len(v1_array),len(v2_array),0]
            extent = [float(min(v1_array).split('V')[0]), float(max(v1_array).split('V')[0]),
                      float(max(v2_array).split('V')[0]), float(min(v2_array).split('V')[0])]
            xlabel = 'waveplate 1 (V)'
            ylabel = 'waveplate 2 (V)'
            title = 'waveplate 3 = {:s}'.format(value)

        if plot_axes is None:
            import matplotlib.pyplot as plt
            plt.figure()
            im = plt.imshow(data, extent=extent)
            cb = plt.colorbar(im)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(title)
        else:
            im = plot_axes.imshow(data, extent=extent)
            # cb = plt.colorbar(im)
            plot_axes.set_xlabel(xlabel)
            plot_axes.set_ylabel(ylabel)
            plot_axes.set_title(title)
        # return data

    @staticmethod
    def plot_scan(data_panel, v1=None, v2=None, v3=None, plot_axes = None):
        """
        plots the detector signal as a function of the axis for which the value is None
        Args:
            data_panel: data set in form of a pandas data panel
            v1: float or text of form {:0.2f} or None
            v2: float or text of form {:0.2f} or None
            v3: float or text of form {:0.2f} or None

        Returns:

        """
        #     need to get two values
        assert np.sum([x == None for x in [v1, v2, v3]]) == 1, 'function requires exactly two not None values'

        if v1 == None:
            perm = [1, 0, 2]
            value1, value2 = v2, v3
            xlabel = 'waveplate 1 (V)'
            title = 'WP2 = {:0.2f}V, WP3 = {:0.2f}V'.format(v2, v3)
        elif v2 == None:
            perm = [0, 1, 2]
            value1, value2 = v1, v3
            xlabel = 'waveplate 2 (V)'
            title = 'WP1 = {:0.2f}V, WP3 = {:0.2f}V'.format(v1, v3)
        elif v3 == None:
            perm = [0, 2, 1]
            value1, value2 = v1, v2
            xlabel = 'waveplate 3 (V)'
            title = 'WP1 = {:0.2f}V, WP2 = {:0.2f}V'.format(v1, v2)

        if isinstance(value1, float):
            value1 = '{:0.2f}V'.format(value1)
        if isinstance(value2, float):
            value2 = '{:0.2f}V'.format(value2)

        print('p', perm, value1, value2)

        scan_data = data_panel.transpose(perm[0], perm[1], perm[2])[value1][value2]

        if plot_axes is None:
            import matplotlib.pyplot as plt
            plt.figure()
            plt.plot([float(x.split('V')[0]) for x in scan_data.index], scan_data.values)
            plt.xlabel(xlabel)
            plt.ylabel('detector signal')
            plt.title(title)
        else:
            plot_axes.plot([float(x.split('V')[0]) for x in scan_data.index], scan_data.values)
            plot_axes.set_xlabel(xlabel)
            plot_axes.set_ylabel('detector signal')
            plot_axes.set_title(title)

    def plot(self, axes1):

        last_key = sorted(self.data.keys())[-1]
        print('last_key', last_key)
        volt_range = np.arange(0, 5, 0.2)
        axes1.plot(volt_range, self.data[last_key], '-o')

    def stop(self):
        self._abort = True









if __name__ == '__main__':
    script = {}
    instr = {}
    script, failed, instr = Script.load_and_append({'pol_control': 'FPGA_PolarizationController'}, script, instr)

    print(script)
    print(failed)
    print(instr)