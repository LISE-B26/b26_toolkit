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
from copy import deepcopy

import numpy as np
import scipy as sp
from PyQt4.QtCore import pyqtSlot

from b26_toolkit.src.instruments import PiezoController
from b26_toolkit.src.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from PyLabControl.src.core import Parameter, Script
from b26_toolkit.src.scripts import GalvoScan
# from src.scripts import GalvoScanNIFpga


class AutoFocusGeneric(Script):
    """
Autofocus: Takes images at different piezo voltages and uses a heuristic to figure out the point at which the objective
            is focused.
    """

    _DEFAULT_SETTINGS = [
        Parameter('save_images', False, bool, 'save image taken at each voltage'),
        Parameter('piezo_center_voltage', 50, float, 'center point of piezo voltage sweep'),
        Parameter('piezo_voltage_width', 5, float, 'V between the minimum and maximum points of the range'),
        Parameter('num_sweep_points', 10, int, 'number of values to sweep between min and max voltage'),
        Parameter('focusing_optimizer', 'standard_deviation',
                  ['mean', 'standard_deviation', 'normalized_standard_deviation'], 'optimization function for focusing'),
        Parameter('wait_time', 0.1, float),
        Parameter('use_current_piezo_voltage', False, bool, 'Overrides piezo_center_voltage and instead uses the current piezo voltage as the center of the range'),
        Parameter('center_on_current_location', False, bool, 'Check to use current galvo location rather than center point in take_image'),
        Parameter('galvo_return_to_initial', False, bool, 'Check to return galvo location to initial value (before calling autofocus)'),
        # Parameter('galvo_position', 'take_image_pta', ['take_image', 'current_location', 'last_run'], 'select galvo location (center point in acquire_image, current location of galvo or location from previous run)')
    ]

    # the take image script depends on the particular hardware, e.g. DAQ or NIFPGA
    _SCRIPTS = {
        'take_image': NotImplemented
    }

    _INSTRUMENTS = {}

    def __init__(self, scripts, instruments = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings, instruments, scripts, log_function= log_function, data_path = data_path)


    @pyqtSlot(int)
    def _receive_signal(self, progress):
        """
        this function takes care of signals emitted by the subscripts
        the default behaviour is that it just reemits the signal
        Args:
            progress: progress of subscript take_image
        """
        # just ignore the signals from the subscript, we just send out our own signal
        pass
        # sender_emitter = self.sender()
        #
        # if self._current_subscript_stage['current_subscript'] is self.scripts['take_image']:
        #     img_index = self._current_subscript_stage['subscript_exec_count']['take_image']
        #     total_img_count = self.settings['num_sweep_points']
        #     progress = 1.* ((img_index -1 + float(progress)/ 100)/ total_img_count)
        # else:
        #     print('WHERE DID THIS SIGNAL COME FROM??? sender', sender_emitter)
        #
        #     current_image = self.scripts['take_image'].data['image_data']
        #     self.data['current_image'] = deepcopy(current_image)
        #     self.data['extent'] = self.scripts['take_image'].data['extent']
        #
        # self.updateProgress.emit(int(100.*progress))

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """


        def calc_focusing_optimizer(image, optimizer):
            """
            calculates a measure for how well the image is focused
            Args:
                optimizer: one of the three strings: mean, standard_deviation, normalized_standard_deviation
            Returns:  measure for how well the image is focused
            """
            if optimizer == 'mean':
                return np.mean(image)
            elif optimizer == 'standard_deviation':
                return np.std(image)
            elif optimizer == 'normalized_standard_deviation':
                return np.std(image) / np.mean(image)

        def autofocus_loop(sweep_voltages):
            """
            this is the actual autofocus loop
            Args:
                sweep_voltages: array of sweep voltages

            Returns:

            """
            # update instrument

            for index, voltage in enumerate(sweep_voltages):

                if self._abort:
                    self.log('Leaving autofocusing loop')
                    break

                # set the voltage on the piezo
                self._step_piezo(voltage, self.settings['wait_time'])

                # take a galvo scan
                self.scripts['take_image'].run()
                self.data['current_image'] = deepcopy(self.scripts['take_image'].data['image_data'])

                # calculate focusing function for this sweep
                self.data['focus_function_result'].append(
                    calc_focusing_optimizer(self.data['current_image'], self.settings['focusing_optimizer']))

                # save image if the user requests it
                if self.settings['save_images']:
                    self.scripts['take_image'].save_image_to_disk(
                        '{:s}\\image_{:03d}.jpg'.format(self.filename_image, index))
                    self.scripts['take_image'].save_data('{:s}\\image_{:03d}.csv'.format(self.filename_image, index),
                                                         'image_data')

                self.progress = 100. * index / len(sweep_voltages)
                self.updateProgress.emit(self.progress if self.progress < 100 else 99)


        if self.settings['save'] or self.settings['save_images']:
            self.filename_image = '{:s}\\image\\'.format(self.filename())
        else:
            self.filename_image = None

        daq_pt = self._get_galvo_location()

        if self.settings['center_on_current_location']:
            self.scripts['take_image'].settings['point_a'].update({'x': daq_pt[0], 'y': daq_pt[1]})

        min_voltage = self.settings['piezo_center_voltage'] - self.settings['piezo_voltage_width']/2.0
        max_voltage = self.settings['piezo_center_voltage'] + self.settings['piezo_voltage_width']/2.0

        sweep_voltages = np.linspace(min_voltage, max_voltage, self.settings['num_sweep_points'])

        self.data['sweep_voltages'] = sweep_voltages
        self.data['focus_function_result'] = []
        self.data['fit_parameters'] = [0, 0, 0, 0]
        self.data['current_image'] = np.zeros([1,1])
        self.data['extent'] = None

        autofocus_loop(sweep_voltages)


        piezo_voltage, self.data['fit_parameters'] = self.fit_focus()
        self._step_piezo(piezo_voltage, self.settings['wait_time'])

        if self.settings['galvo_return_to_initial']:
            self._set_galvo_location(daq_pt)

    def fit_focus(self):
        '''
        fit the data and set piezo to focus spot

        if fails return None otherwise it returns the voltage for the piezo
        '''

        noise_guess = np.min(self.data['focus_function_result'])
        amplitude_guess = np.max(self.data['focus_function_result']) - noise_guess
        center_guess = np.mean(self.data['sweep_voltages'])
        width_guess = 0.8

        p2 = [noise_guess, amplitude_guess, center_guess, width_guess]

        return_voltage = None
        try:
            p2, success = sp.optimize.curve_fit(self.gaussian, self.data['sweep_voltages'],
                                                self.data['focus_function_result'], p0=p2,
                                                bounds=([0, [np.inf, np.inf, 100., 100.]]), max_nfev=2000)

            return_voltage = p2[2]

            self.log('Found fit parameters: ' + str(p2))
        except(ValueError, RuntimeError):
            average_voltage = np.mean(self.data['sweep_voltages'])
            self.log(
                'Could not converge to fit parameters, setting piezo to middle of sweep range, {0} V'.format(
                    average_voltage))
            # z_piezo.update({'voltage': float(average_voltage)})
        finally:
            # even if there is an exception we want to script to continue
            sweep_voltages = self.data['sweep_voltages']
            if not return_voltage is None:
                if return_voltage > sweep_voltages[-1]:
                    return_voltage = float(sweep_voltages[-1])
                    self.log('Best fit found center to be above max sweep range, setting voltage to max, {0} V'.format(return_voltage))
            elif return_voltage < sweep_voltages[0]:
                return_voltage = float(sweep_voltages[0])
                self.log('Best fit found center to be below min sweep range, setting voltage to min, {0} V'.format(return_voltage))

            return return_voltage, p2

    def _get_galvo_location(self):
        """
        returns the current position of the galvo
        Returns: list with two floats, which give the x and y position of the galvo mirror

        """
        raise NotImplementedError

    def _set_galvo_location(self):
        """
        sets the current position of the galvo
        galvo_position: list with two floats, which give the x and y position of the galvo mirror
        """
        raise NotImplementedError

    def _step_piezo(self, voltage, wait_time):
        """
        steps the piezo.  Has to be overwritten specifically for each different hardware realization
        voltage: target piezo voltage
        wait_time: settle time after voltage step
        """
        raise NotImplementedError

    def _plot(self, axes_list, data = None):
        # fit the data and set piezo to focus spot
        if data is None:
            data  = self.data
        axis_focus, axis_image = axes_list

        # if take image is running we take the data from there otherwise we use the scripts own image data
        if self._current_subscript_stage['current_subscript'] is self.scripts['take_image']:

            if 'image_data' in self.scripts['take_image'].data:
                current_image = self.scripts['take_image'].data['image_data']
                extent = self.scripts['take_image'].data['extent']
            else:
                current_image = None
        else:
            current_image = data['current_image']
            extent = data['extent']
        if current_image is not None:
            plot_fluorescence_new(current_image, extent, axis_image)

        if 'focus_function_result' in data:
            focus_data = data['focus_function_result']
            sweep_voltages = data['sweep_voltages']
            if len(focus_data)>0:
                axis_focus.plot(sweep_voltages[0:len(focus_data)],focus_data)
                if not (np.array_equal(data['fit_parameters'], [0,0,0,0])):
                    axis_focus.plot(sweep_voltages[0:len(focus_data)], self.gaussian(sweep_voltages[0:len(focus_data)], *self.data['fit_parameters']), 'k')
                axis_focus.hold(False)


        axis_focus.set_xlabel('Piezo Voltage [V]')

        if self.settings['focusing_optimizer'] == 'mean':
            ylabel = 'Image Mean [kcounts]'
        elif self.settings['focusing_optimizer'] == 'standard_deviation':
            ylabel = 'Image Standard Deviation [kcounts]'
        elif self.settings['focusing_optimizer'] == 'normalized_standard_deviation':
            ylabel = 'Image Normalized Standard Deviation [arb]'
        else:
            ylabel = self.settings['focusing_optimizer']

        axis_focus.set_ylabel(ylabel)
        axis_focus.set_title('Autofocusing Routine')

    def _update_plot(self, axes_list):
        # fit the data and set piezo to focus spot

        axis_focus, axis_image = axes_list

        # if take image is running we take the data from there otherwise we use the scripts own image data
        if self._current_subscript_stage['current_subscript'] is self.scripts['take_image']:
            if 'image_data' in self.scripts['take_image'].data:
                current_image = self.scripts['take_image'].data['image_data']
            else:
                current_image = None
        else:
            current_image = self.data['current_image']

        if current_image is not None:
            update_fluorescence(current_image, axis_image)

        axis_focus, axis_image = axes_list

        update_fluorescence(self.data['current_image'], axis_image)

        focus_data = self.data['focus_function_result']
        sweep_voltages = self.data['sweep_voltages']
        if len(focus_data) > 0:
            axis_focus.plot(sweep_voltages[0:len(focus_data)], focus_data)
            axis_focus.hold(False)

    def gaussian(self, x, noise, amp, center, width):
        return (noise + amp * np.exp(-1.0 * (np.square((x - center)) / (2 * (width ** 2)))))


class AutoFocusNIFPGA(AutoFocusGeneric):
    """
Autofocus: Takes images at different piezo voltages and uses a heuristic to figure out the point at which the objective
            is focused.
    """

    _SCRIPTS = {
        # 'take_image': GalvoScanNIFpga
    }

    def __init__(self, scripts, instruments = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings, instruments, scripts, log_function= log_function, data_path = data_path)

    def _step_piezo(self, voltage, wait_time):
        """
        steps the piezo.  Has to be overwritten specifically for each different hardware realization
        voltage: target piezo voltage
        wait_time: settle time after voltage step
        """
        fpga_instr = self.scripts['take_image'].instruments['NI7845RGalvoScan']['instance']
        # set the voltage on the piezo
        fpga_instr.piezo = float(voltage)
        time.sleep(wait_time)

    def fit_focus(self):
        '''
        fit the data and set piezo to focus spot
        '''
        gaussian = lambda x, noise, amp, center, width: noise + amp * np.exp(
            -1.0 * (np.square((x - center)) / (2 * (width ** 2))))

        noise_guess = np.min(self.data['focus_function_result'])
        amplitude_guess = np.max(self.data['focus_function_result']) - noise_guess
        center_guess = np.mean(self.data['sweep_voltages'])
        width_guess = 0.8

        reasonable_params = [noise_guess, amplitude_guess, center_guess, width_guess]

        try:
            p2, success = sp.optimize.curve_fit(gaussian, self.data['sweep_voltages'],
                                                self.data['focus_function_result'], p0=reasonable_params,
                                                bounds=([0, [np.inf, np.inf, 100., 100.]]), max_nfev=2000)

            self.log('Found fit parameters: ' + str(p2))

            if p2[2] > sweep_voltages[-1]:
                fpga_instr.piezo = float(sweep_voltages[-1])
                self.log(
                    'Best fit found center to be above max sweep range, setting voltage to max, {0} V'.format(
                        sweep_voltages[-1]))
            elif p2[2] < sweep_voltages[0]:
                fpga_instr.piezo = float(sweep_voltages[0])
                self.log(
                    'Best fit found center to be below min sweep range, setting voltage to min, {0} V'.format(
                        sweep_voltages[0]))
            else:
                fpga_instr.piezo = float(p2[2])
        except(ValueError):
            p2 = [0, 0, 0, 0]
            average_voltage = np.mean(self.data['sweep_voltages'])
            self.log(
                'Could not converge to fit parameters, setting piezo to middle of sweep range, {0} V'.format(
                    average_voltage))
            fpga_instr.piezo = float(average_voltage)

        return p2


class AutoFocusDAQ(AutoFocusGeneric):
    """
Autofocus: Takes images at different piezo voltages and uses a heuristic to figure out the point at which the objective
            is focused.
    """

    # _DEFAULT_SETTINGS = []

    _INSTRUMENTS = {
        'z_piezo': PiezoController
    }
    _SCRIPTS = {
        'take_image': GalvoScan
    }

    def __init__(self, scripts, instruments = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings, instruments, scripts, log_function= log_function, data_path = data_path)

    def _step_piezo(self, voltage, wait_time):
        """
        steps the piezo.  Has to be overwritten specifically for each different hardware realization
        voltage: target piezo voltage
        wait_time: settle time after voltage step
        """
        z_piezo = self.instruments['z_piezo']['instance']
        # set the voltage on the piezo
        z_piezo.voltage = float(voltage)
        time.sleep(wait_time)

    def _get_galvo_location(self):
        """
        returns the current position of the galvo
        Returns: list with two floats, which give the x and y position of the galvo mirror
        """
        galvo_position = self.scripts['take_image'].get_galvo_location()
        # galvo_position = self.scripts['take_image'].instruments['daq']['instance'].get_analog_out_voltages([
        #     self.scripts['take_image'].settings['DAQ_channels']['x_ao_channel'],
        #     self.scripts['take_image'].settings['DAQ_channels']['y_ao_channel']]
        # )
        return galvo_position

    def _set_galvo_location(self, galvo_position):
        """
        sets the current position of the galvo
        galvo_position: list with two floats, which give the x and y position of the galvo mirror
        """
        self.scripts['take_image'].set_galvo_location(galvo_position)
        # pt = galvo_position
        # daq = self.scripts['take_image'].instruments['daq']['instance']
        # # daq API only accepts either one point and one channel or multiple points and multiple channels
        # pt = np.transpose(np.column_stack((pt[0],pt[1])))
        # pt = (np.repeat(pt, 2, axis=1))
        #
        # daq.AO_init([self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], pt)
        # daq.AO_run()
        # daq.AO_waitToFinish()
        # daq.AO_stop()

    def _function(self):
        #update piezo settings
        if self.settings['use_current_piezo_voltage']:
            self.settings['piezo_center_voltage'] = self.instruments['z_piezo']['instance'].read_probes('voltage')
        AutoFocusGeneric._function(self)



if __name__ == '__main__':


    # from PyLabControl.src.core.read_write_functions import load_b26_file
    #
    # in_data = load_b26_file('C:\\b26_tmp\\gui_settings.b26')
    #
    # instruments = in_data['instruments']
    # scripts = in_data['scripts']
    # probes = in_data['probes']
    #
    # instruments, failed = Instrument.load_and_append(instruments)
    # scripts, failed, instruments = Script.load_and_append(
    #     script_dict=scripts,
    #     instruments=instruments,
    #     data_path='c:\\')


    # scripts, loaded_failed, instruments = Script.load_and_append({"af": 'AutoFocusNIFPGA'})
    # print('===++++++===========++++++===========++++++========')
    # print(scripts)
    # print('===++++++===========++++++===========++++++========')
    # print(scripts['af'].scripts['take_image'].instruments['NI7845RGalvoScan'])
    # print(type(scripts['af'].scripts['take_image'].instruments['NI7845RGalvoScan']['settings']))
    # print(type(scripts['af'].scripts['take_image'].instruments['NI7845RGalvoScan']['instance']))
    #
    # self.scripts, failed, self.instruments = Script.load_and_append(
    #     script_dict=scripts,
    #     instruments=self.instruments,
    #     log_function=self.log,
    #     data_path=self.gui_settings['data_folder'])



    # ========================= test explicitely loading for DAQ ========================================
    scripts, loaded_failed, instruments = Script.load_and_append({'take_image': 'GalvoScanDAQ'})
    print('===++++++===========++++++===========++++++========')
    print(scripts)
    print('===++++++===========++++++===========++++++========')

    af = AutoFocusNIFPGA(scripts=scripts)
    print(af)

    # ========================= test explicitely loading for DAQ ========================================
    scripts, loaded_failed, instruments = Script.load_and_append({'take_image': 'GalvoScanDAQ'})
    print('===++++++===========++++++===========++++++========')
    print(scripts)
    print('===++++++===========++++++===========++++++========')

    af = AutoFocusNIFPGA(scripts=scripts)
    print(af)
