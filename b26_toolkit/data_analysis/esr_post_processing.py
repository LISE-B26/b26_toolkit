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

import numpy as np
import pandas as pd
import glob, os
import matplotlib.pyplot as plt
from matplotlib import gridspec
import time
import datetime


from pylabcontrol.core import Script
from b26_toolkit.data_processing.esr_signal_processing import fit_esr, find_nv_peaks
from b26_toolkit.plotting.plots_1d import plot_esr
from pylabcontrol.core.helper_functions import datetime_from_str

from b26_toolkit.data_processing.esr_signal_processing import get_lorentzian_fit_starting_values, fit_lorentzian, fit_double_lorentzian


from b26_toolkit.scripts.find_nv import FindNV

from win32com.client import Dispatch
import pythoncom

# stuff for manual correction
import queue
import ipywidgets as widgets
from IPython.display import display, clear_output
import threading
import sys

freq_to_mag = 1. / (2 * 2.8e6)

def process_esrs(source_folder, target_folder):
    #first autofit esrs to get initial estimate
    fit_data_set = autofit_esrs(source_folder)

    if fit_data_set is None:
        print('Failed to find any ESRs to analyze. Check your filepaths.')
        return

    # define list shared between threads
    nv_type_manual = [''] * len(fit_data_set)
    b_field_manual = [np.nan] * len(fit_data_set)



    # define queue for inter-thread communication
    next_queue = queue.Queue()
    current_id_queue = queue.Queue()

    #set up manual correction gui
    widget_list = process_manual_esrs(nv_type_manual, b_field_manual, next_queue, current_id_queue)

    # launch manual_correction function on separate thread, as otherwise buttons will not respond until that
    # function ends, that function is waiting on buttons to continue, and you deadlock
    thread = threading.Thread(target=manual_correction, args=(source_folder, target_folder, fit_data_set, nv_type_manual, b_field_manual, next_queue, current_id_queue, widget_list[4], widget_list[5], widget_list[6], widget_list[7]))
    thread.start()
    # manual_correction(source_folder, target_folder, fit_data_set, nv_type_manual, b_field_manual, next_queue, current_id_queue, widget_list[4], widget_list[5], widget_list[6], widget_list[7])


    # thread.join()
    # visualize_magnetic_fields(source_folder, target_folder)

    # for widget in widget_list:
    #     widget.close()

def get_nv_type(fit_params):
    "based on the fit figure out the nv type: split / single / no_peak"

    if fit_params is None:
        nv_type = 'no_peak'
    else:
        if len(fit_params) == 4:
            if np.abs(fit_params[2] - 2.87e9) > 0.01e9:
                nv_type = 'split'
            else:
                nv_type = 'single'
        elif len(fit_params) == 6:
            nv_type = 'split'
        else:
            raise TypeError

    return nv_type

def get_B_field(nv_type, fit_params):
    """
    calculate the magentic field (in Gauss) from the nv_type and fit_params
    """

    Bfield = None

    if nv_type == 'single':
        Bfield = 0
    elif nv_type == 'split':
        if len(fit_params) == 4:
            Bfield = np.abs(fit_params[2] - 2.87e9) * 2.
        elif len(fit_params) == 6:
            Bfield = np.abs(fit_params[5] - fit_params[4])
        Bfield /= 5.6e6

    return Bfield

def create_shortcut(src_path, dst_path):
    """
    Creates a shortcut in src_path to dst_path.
    Pywin32 can't handle / escaping in dst_path, must be escaped by \\
    Args:
        src_path: location to create shortcut
        dst_path: target of shortcut
    """
    pythoncom.CoInitialize()
    shell = Dispatch("WScript.Shell")
    shortcut = shell.CreateShortCut(src_path)
    shortcut.Targetpath = dst_path
    shortcut.save()

def autofit_esrs(folder):
    """

    fits the esr data, plots them and asks the user for confirmation, the fit data is saved to the folder target_folder with the same structure as folders

    Args:
        folders: source folder with esr data, this folder shoudl contain a subfolder data_subscripts which contains subfolders *esr* with the esr data
        target_folder: target folder where the output data is saved in form of a .csv file

    Returns: fitdataset as a pandas array

    """
    # loop over all the folders in the data_subscripts subfolder and retrieve fitparameters and position of NV
    esr_folders = glob.glob(os.path.join(folder, './data_subscripts/*esr*'))

    fit_data_set = None

    # classify the nvs according to the following categories
    # by default we set this to na (not available) and try to figure it out based on the data and fitquality
    # nv_type = 'na' # split / single / no_peak / no_nv / na

    for i, esr_folder in enumerate(esr_folders):

        # find the NV index
        pt_id = int(os.path.basename(esr_folder).split('pt_')[-1])

        findnv_folder =  sorted(glob.glob(folder + '/data_subscripts/*find_nv*pt_*{:d}'.format(pt_id)))[0]

        # load data
        data = Script.load_data(esr_folder)
        fit_params = fit_esr(data['frequency'], data['data'])
        nv_type = get_nv_type(fit_params)

        # get initial guess for peaks
        freq_peaks, ampl_peaks = find_nv_peaks(data['frequency'], data['data'])

        # get nv positions
        data_pos = Script.load_data(findnv_folder)
        pos = data_pos['maximum_point']
        pos_init = data_pos['initial_point']

        if fit_params is None:
            fit_data_set_single = {}
        else:
            fit_data_set_single = {'fit_{:d}'.format(i): f for i, f in enumerate(fit_params)}

        fit_data_set_single.update({'id': pt_id, 'NV type': nv_type})

        fit_data_set_single.update({'x': pos['x'][0], 'y': pos['y'][0]})
        fit_data_set_single.update({'xo': pos_init['x'][0], 'yo': pos_init['y'][0]})

        fit_data_set_single.update({'B-field (gauss)': get_B_field(nv_type, fit_params)})

        # convert to dataframe
        fit_data_set_single = pd.DataFrame.from_dict({k: [v] for k, v in fit_data_set_single.items()})

        if fit_data_set is None:
            fit_data_set = pd.DataFrame(fit_data_set_single)
        else:
            fit_data_set = fit_data_set.append(fit_data_set_single, ignore_index=True)

    return fit_data_set

def manual_correction(folder, target_folder, fit_data_set, nv_type_manual, b_field_manual, queue, current_id_queue, lower_peak_widget, upper_peak_widget, lower_fit_widget, upper_fit_widget):
    """
    Backend code to display and fit ESRs, then once input has been received from front-end, incorporate the data into the
    current data set

    Args:
        folders: folder containing the data to analyze
        target_folder: target for processed data
        fit_data_set: starting point data set containing automatically fitted esrs
        nv_type_manual: pointer to empty array to populate with nv type manual corrections
        b_field_manual: pointer to empty array to populate with b field manual corrections
        queue: queue for communication with front-end in separate thread
        lower_peak_widget: widget containing frequency value of lower peak
        upper_peak_widget: widget containing frequency value of upper peak

    Poststate: populates fit_data_set with manual corrections
    """

    lower_peak_manual = [np.nan] * len(fit_data_set)
    upper_peak_manual = [np.nan] * len(fit_data_set)

    filepath = os.path.join(target_folder, folder[2:])
    data_filepath = os.path.join(filepath, 'data-manual.csv')
    if os.path.exists(data_filepath):
        prev_data = pd.read_csv(data_filepath)
        if 'manual_peak_1' in list(prev_data.keys()):
            for i in range(0, len(prev_data['manual_B_field'])):
                b_field_manual[i] = prev_data['manual_B_field'][i]
                nv_type_manual[i] = prev_data['manual_nv_type'][i]
            lower_peak_manual = prev_data['manual_peak_1']
            upper_peak_manual = prev_data['manual_peak_2']


    #TODO: Add saving as you go, add ability to start at arbitrary NV, add ability to specify a next NV number, eliminate peak/correct -> only have 'accept fit'

    try:


        print('STARTING')

        fit_data_set_array = fit_data_set.as_matrix()

        w = widgets.HTML("Event information appears here when you click on the figure")
        display(w)

        # loop over all the folders in the data_subscripts subfolder and retrieve fitparameters and position of NV
        esr_folders = glob.glob(os.path.join(folder, '.\\data_subscripts\\*esr*'))

        # create folder to save images to
        # filepath_image = os.path.join(target_folder, os.path.dirname(folder).split('./')[1])
        # image_folder = os.path.join(filepath_image, '{:s}\\images'.format(os.path.basename(folder)))
        image_folder = os.path.join(target_folder, '{:s}\\images'.format(folder[2:]))
        # image_folder = os.path.normpath(
        #     os.path.abspath(os.path.join(os.path.join(target_folder, 'images'), os.path.basename(folders[0]))))
        if not os.path.exists(image_folder):
            os.makedirs(image_folder)
        if not os.path.exists(os.path.join(image_folder, 'bad_data')):
            os.makedirs(os.path.join(image_folder, 'bad_data'))

        f = plt.figure(figsize=(12, 6))

        def onclick(event):
            if event.button == 1:
                if event.key == 'control':
                    lower_fit_widget.value = event.xdata
                else:
                    lower_peak_widget.value = event.xdata
            elif event.button == 3:
                if event.key == 'control':
                    upper_fit_widget.value = event.xdata
                else:
                    upper_peak_widget.value = event.xdata

        cid = f.canvas.mpl_connect('button_press_event', onclick)

        data_array = []
        data_pos_array = []
        for esr_folder in esr_folders:
            print(esr_folder)
            sys.stdout.flush()
            data = Script.load_data(esr_folder)
            data_array.append(data)
            print('looping')
            sys.stdout.flush()


        nv_folders = glob.glob(folder + '\\data_subscripts\\*find_nv*pt_*')
        for nv_folder in nv_folders:
            data_pos_array.append(Script.load_data(nv_folder))

        while True:


        # for i, esr_folder in enumerate(esr_folders):

            i = current_id_queue.queue[0]
            if i >= len(data_array):
                break

            lower_fit_widget.value = 0
            upper_fit_widget.value = 10e9

            lower_peak_widget.value = 2.87e9
            upper_peak_widget.value = 0

            def display_data(pt_id, lower_peak_widget = None, upper_peak_widget = None, display_fit = True):
                # find the NV index
                # pt_id = int(os.path.basename(esr_folder).split('pt_')[-1])

                # findnv_folder = glob.glob(folder + '\\data_subscripts\\*find_nv*pt_*{:d}'.format(pt_id))[0]

                f.clf()
                gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
                ax0 = plt.subplot(gs[0])
                ax1 = plt.subplot(gs[1])
                ax = [ax0, ax1]
                plt.suptitle('NV #{:d}'.format(pt_id), fontsize=16)

                # load data
                data = data_array[i]
                if lower_fit_widget.value == 0 and upper_fit_widget.value == 10e9:
                    freq = data['frequency']
                    ampl = data['data']
                else:
                    freq = data['frequency'][np.logical_and(data['frequency'] > lower_fit_widget.value, data['frequency'] < upper_fit_widget.value)]
                    ampl = data['data'][np.logical_and(data['frequency'] > lower_fit_widget.value, data['frequency'] < upper_fit_widget.value)]
                if lower_peak_widget is None:
                    fit_params = fit_data_set_array[pt_id, 2:8]
                else:
                    lower_peak = lower_peak_widget.value
                    upper_peak = upper_peak_widget.value
                    if upper_peak == 0:
                        start_vals = get_lorentzian_fit_starting_values(freq, ampl)
                        start_vals[2] = lower_peak
                        start_vals[1] = ampl[np.argmin(np.abs(freq - lower_peak))] - start_vals[0]
                        try:
                            fit_params = fit_lorentzian(freq, ampl, starting_params=start_vals,
                                                 bounds=[(0, -np.inf, 0, 0), (np.inf, 0, np.inf, np.inf)])
                        except:
                            # ESR fit failed!
                            fit_params = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

                    else:
                        center_freq = np.mean(freq)
                        start_vals = []
                        start_vals.append(
                            get_lorentzian_fit_starting_values(freq[freq < center_freq], ampl[freq < center_freq]))
                        start_vals.append(
                            get_lorentzian_fit_starting_values(freq[freq > center_freq], ampl[freq > center_freq]))
                        start_vals = [
                            np.mean([start_vals[0][0], start_vals[1][0]]),  # offset
                            np.sum([start_vals[0][3], start_vals[1][3]]),  # FWHM
                            ampl[np.argmin(np.abs(freq-lower_peak))] - start_vals[0][0], ampl[np.argmin(np.abs(freq-upper_peak))]- start_vals[1][0],  # amplitudes
                            lower_peak, upper_peak  # centers
                        ]
                        try:
                            fit_params = fit_double_lorentzian(freq, ampl, starting_params=start_vals, bounds=
                            [(0, 0, -np.inf, -np.inf, min(freq), min(freq)), (np.inf, np.inf, 0, 0, max(freq), max(freq))])
                        except:
                            fit_params = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

                if len(fit_params) == 4 or np.isnan(fit_params[4]):
                    fit_params = fit_params[0:4]

                # get nv positions
                #             data_pos = {'initial_point': [fit_data_set['xo'].values[pt_id]]}
                data_pos = data_pos_array[i]
                #             pos = data_pos['maximum_point']
                #             pos_init = data_pos['initial_point']

                # plot NV image
                FindNV.plot_data([ax[1]], data_pos)

                # plot data and fits
                # print("fit_params: ", fit_params)

                sys.stdout.flush()

                if display_fit:
                    plot_esr(ax[0], data['frequency'], data['data'], fit_params=fit_params)
                else:
                    plot_esr(ax[0], data['frequency'], data['data'], fit_params=[np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])

                plt.tight_layout()
                plt.subplots_adjust(top=0.85) # Makes room at top of plot for figure suptitle

                plt.draw()
                plt.show()

                return fit_params, pt_id

            fit_params, pt_id = display_data(i)
            if len(fit_params) == 6:
                lower_peak_widget.value = fit_params[4]
                upper_peak_widget.value = fit_params[5]
            elif len(fit_params) == 4:
                lower_peak_widget.value = fit_params[2]
                upper_peak_widget.value = 0

            while True:
                if queue.empty():
                    time.sleep(.5)
                else:
                    value = queue.get()
                    if value == -1:
                        fit_params, point_id = display_data(i, lower_peak_widget=lower_peak_widget, upper_peak_widget=upper_peak_widget)
                        if len(fit_params) == 6:
                            lower_peak_widget.value = fit_params[4]
                            upper_peak_widget.value = fit_params[5]
                        elif len(fit_params) == 4:
                            lower_peak_widget.value = fit_params[2]
                            upper_peak_widget.value = 0
                        continue
                    elif value == -2:
                        display_data(display_fit = False)
                        fit_params = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
                        lower_fit_widget.value = 0
                        upper_fit_widget.value = 10e9
                    else:
                        break

            if nv_type_manual[i] == 'split':
                if np.isnan(fit_params[0]):
                    lower_peak_manual[i] = lower_peak_widget.value
                    upper_peak_manual[i] = upper_peak_widget.value
                    b_field_manual[i] = ((upper_peak_widget.value - lower_peak_widget.value) / 5.6e6)
                elif len(fit_params) == 4:
                    lower_peak_manual[i] = fit_params[2]
                    b_field_manual[i] = (np.abs(2.87e9-fit_params[2]) / 2.8e6)
                else:
                    lower_peak_manual[i] = fit_params[4]
                    upper_peak_manual[i] = fit_params[5]
                    b_field_manual[i] = ((fit_params[5] - fit_params[4]) / 5.6e6)
            elif nv_type_manual[i] == 'single':
                if np.isnan(fit_params[0]):
                    lower_peak_manual[i] = lower_peak_widget.value
                    b_field_manual[i] = 0
                else:
                    lower_peak_manual[i] = fit_params[2]
                    b_field_manual[i] = 0

            if nv_type_manual[i] == '':
                if fit_params is None:
                    f.savefig(os.path.join(os.path.join(image_folder, 'bad_data'), 'esr_pt_{:02d}.jpg'.format(pt_id)))
                else:
                    f.savefig(os.path.join(image_folder, 'esr_pt_{:02d}.jpg'.format(pt_id)))
            else:
                if nv_type_manual[i] in ['bad', 'no_split']:
                    f.savefig(os.path.join(os.path.join(image_folder, 'bad_data'), 'esr_pt_{:02d}.jpg'.format(pt_id)))
                else:
                    f.savefig(os.path.join(image_folder, 'esr_pt_{:02d}.jpg'.format(pt_id)))

            if not os.path.exists(filepath):
                os.makedirs(filepath)
            fit_data_set['manual_B_field'] = b_field_manual
            fit_data_set['manual_nv_type'] = nv_type_manual
            fit_data_set['manual_peak_1'] = lower_peak_manual
            fit_data_set['manual_peak_2'] = upper_peak_manual
            fit_data_set.to_csv(data_filepath)

        f.canvas.mpl_disconnect(cid)
        fit_data_set['manual_B_field'] = b_field_manual
        fit_data_set['manual_nv_type'] = nv_type_manual
        fit_data_set['manual_peak_1'] = lower_peak_manual
        fit_data_set['manual_peak_2'] = upper_peak_manual

        # filepath = os.path.join(target_folder, folder[2:])
        # data_filepath = os.path.join(filepath, 'data-manual.csv')
        # filename = '{:s}\\data-manual.csv'.format(os.path.basename(folder))
        # filepath = os.path.join(target_folder, os.path.dirname(folder).split('./')[1])
        # data_filepath = os.path.join(filepath, filename)

        if not os.path.exists(filepath):
            os.makedirs(filepath)

        fit_data_set.to_csv(data_filepath)

        create_shortcut(os.path.abspath(os.path.join(filepath, 'to_data.lnk')), os.path.abspath(folder))
        create_shortcut(os.path.join(os.path.abspath(folder), 'to_processed.lnk'), os.path.abspath(filepath))

        print('DONE!')

    except Exception as e:
        print(e)
        raise
        # exc_type, exc_value, traceback = sys.exc_info()
        # error_widget.value = str(exc_type) + ', ' + exc_value + ', ' + str(traceback)


def process_manual_esrs(nv_type_manual, b_field_manual, next_queue, current_id_queue):
    """
    Sets up the gui for manual correction.
    Args:
        nv_type_manual: pointer to array containing manual corrections for nv_type
        b_field_manual: pointer to array containing manual corrections for b_field
        next_queue:

    Returns:

    """
    # define queue for inter-thread communication
    current_id_queue.put(0)

    # define buttons to be added to display
    # error_box = widgets.Text(description='error message')
    button_correct = widgets.Button(description="accept fit")
    button_bad = widgets.Button(description="bad")
    button_no_peak = widgets.Button(description="no_peak")
    button_refit = widgets.Button(description = "refit")
    button_clear_refit = widgets.Button(description = "clear refit")
    next_NV_index = widgets.IntText(description = 'next NV number')
    next_NV_index.value = 1
    lower_peak = widgets.FloatText(description='lower (single) peak')
    upper_peak = widgets.FloatText(description='upper (second) peak')
    lower_fit = widgets.FloatText(description = 'fitting lower bound')
    upper_fit = widgets.FloatText(description = 'fitting upper bound')
    widget_list = [button_correct, next_NV_index, button_bad, button_no_peak, lower_peak, upper_peak, lower_fit, upper_fit]

    # box = widgets.VBox([widgets.HBox([button_correct, button_bad]),
    #                     widgets.HBox([button_no_peak, button_peak]),
    #                     widgets.HBox([button_refit, button_clear_refit]),
    #                     widgets.HBox([lower_peak, upper_peak]),
    #                     widgets.HBox([lower_fit, upper_fit])])
    box = widgets.HBox([widgets.VBox([button_correct, button_no_peak, button_refit, lower_peak, lower_fit]),
                        widgets.VBox([next_NV_index, button_bad, button_clear_refit, upper_peak, upper_fit])])

    # display all widgets
    # display(error_box)
    # display(button_correct)
    # display(button_bad)
    # display(button_no_peak)
    # display(button_peak)
    # display(button_refit)
    # display(button_clear_refit)
    # display(lower_peak)
    # display(upper_peak)
    # display(lower_fit)
    # display(upper_fit)
    display(box)


    def button_correct_clicked(b):
        current_id = current_id_queue.get()
        if upper_peak.value == 0:
            if np.abs(lower_peak.value - 2.87e9) > 0.05e9:
                print(('current_id', current_id))
                nv_type_manual[current_id] = 'split'
                # b_field_manual[current_id] = (np.abs(lower_peak.value - 2.87e9) * 2.) / (5.6e6)
                lower_peak.value = 0
            else:
                nv_type_manual[current_id] = 'single'
        else:
            nv_type_manual[current_id] = 'split'
            # b_field_manual[current_id] = (upper_peak.value - lower_peak.value) / (5.6e6)
        # current_id_queue.put(current_id + 1)
        current_id_queue.put(next_NV_index.value)
        next_NV_index.value = next_NV_index.value + 1
        next_queue.put(0)

    button_correct.on_click(button_correct_clicked)

    def button_bad_clicked(b):
        current_id = current_id_queue.get()
        nv_type_manual[current_id] = 'bad'
        # current_id_queue.put(current_id + 1)
        current_id_queue.put(next_NV_index.value)
        next_NV_index.value = next_NV_index.value + 1
        next_queue.put(0)

    button_bad.on_click(button_bad_clicked)

    def button_no_split_clicked(b):
        current_id = current_id_queue.get()
        nv_type_manual[current_id] = 'no_split'
        # current_id_queue.put(current_id + 1)
        current_id_queue.put(next_NV_index.value)
        next_NV_index.value = next_NV_index.value + 1
        next_queue.put(0)

    button_no_peak.on_click(button_no_split_clicked)

    def button_refit_clicked(b):
        next_queue.put(-1)

    button_refit.on_click(button_refit_clicked)


    def button_clear_refit_clicked(b):
        next_queue.put(-2)

    button_clear_refit.on_click(button_clear_refit_clicked)

    return widget_list

def get_select_points(folder):
    """
    find the select_points folder taken previous to the data in folder
    """
    time_folder = datetime_from_str(os.path.basename(folder)[0:15])

    select_points_folders = glob.glob(os.path.join(os.path.dirname(folder), '*select_points*'))

    # get all the time from select_points folders
    select_points_folders_time = [datetime_from_str(os.path.basename(f)[0:15]) for f in select_points_folders]
    # keep only the ones that were take before the data in 'folder' was taken
    select_points_folders_time = sorted(
        [t for t in select_points_folders_time if time_folder - t > datetime.timedelta(0)])
    # take the last select_points folder
    select_points_folders_time = select_points_folders_time[-1].strftime('%y%m%d-%H_%M_%S')
    # find the folder
    select_points_folder = glob.glob(
        os.path.join(os.path.dirname(folder), '{:s}*select_points*'.format(select_points_folders_time)))

    return select_points_folder[0]

def visualize_magnetic_fields(src_folder, target_folder, manual=True):
    filepath = os.path.join(target_folder, os.path.dirname(src_folder).split('./')[1])

    # load the fit_data
    if manual:
        #         df = pd.read_csv(os.path.join(target_folder, '{:s}-manual.csv'.format(folder.split('./')[1].replace('/','-'))), index_col = 'id', skipinitialspace=True)
        filename = '{:s}\\data-manual.csv'.format(os.path.basename(src_folder))
    else:
        filename = '{:s}.csv'.format(os.path.basename(src_folder))
    # df = pd.read_csv(os.path.join(target_folder, '{:s}.csv'.format(folder.split('./')[1].replace('/','-'))), index_col = 'id', skipinitialspace=True)
    df = pd.read_csv(os.path.join(filepath, filename), index_col='id', skipinitialspace=True)
    df = df.drop('Unnamed: 0', 1)

    # include manual corrections
    if manual:
        for id, nv_type in enumerate(df['manual_nv_type']):
            if not pd.isnull(nv_type):
                df.set_value(id, 'NV type', nv_type)
                if nv_type == 'split':
                    df.set_value(id, 'B-field (gauss)', df.get_value(id, 'manual_B_field'))

    # load the image data
    select_points_data = Script.load_data(get_select_points(src_folder))
    image_data = select_points_data['image_data']
    #     points_data = select_points_data['nv_locations']
    extent = select_points_data['extent']

    # prepare figure
    f = plt.figure(figsize=(15, 8))
    gs = gridspec.GridSpec(1, 2, height_ratios=[1, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax = [ax0, ax1]

    # plot the map
    ax0.imshow(image_data, extent=extent, interpolation='nearest', cmap='pink')

    for nv_type, color in zip(['split', 'bad', 'no_peak', 'single'], ['g', 'k', 'b', 'r']):
        subset = df[df['NV type'] == nv_type]
        ax0.scatter(subset['xo'], subset['yo'], c=color, label=nv_type)

    ax0.legend(bbox_to_anchor=(2.3, 1))

    subset = df[df['NV type'] == 'split']
    for i, j, n in zip(subset['xo'], subset['yo'], subset.index):
        corr = +0.005  # adds a little correction to put annotation in marker's centrum
        ax0.annotate(str(n), xy=(i + corr, j + corr), fontsize=8, color='w')
    ax0.set_xlabel('x (V)')
    ax0.set_ylabel('y (V)')

    # plot the fields on a 1D plot
    subset = df[df['NV type'] == 'split']
    ax1.plot(subset['xo'], subset['B-field (gauss)'] / freq_to_mag * 1e-6, 'o')
    ax1.set_title('ESR Splittings')
    ax1.set_xlabel('x coordinate (V)')
    ax1.set_ylabel('Splitting (MHz)')
    ax1.set_xlim([extent[0], extent[1]])

    ax2 = ax1.twinx()
    mn, mx = ax1.get_ylim()
    ax2.set_ylim(mn * freq_to_mag * 1e6, mx * freq_to_mag * 1e6)
    ax2.set_ylabel('Magnetic Field Projection (Gauss)')

    for i, j, n in zip(subset['xo'], subset['B-field (gauss)'] / freq_to_mag * 1e-6, subset.index):
        corr = 0.0005  # adds a little correction to put annotation in marker's centrum
        ax1.annotate(str(n), xy=(i + corr, j + corr))

    # f.set_tight_layout(True)

    #     f.savefig(os.path.join(target_folder, '{:s}.jpg'.format(os.path.basename(folder))),
    #                bbox_inches='tight',
    # #                transparent=True,
    #                pad_inches=0)

    f.savefig(os.path.join(filepath, '{:s}.jpg'.format(os.path.basename(src_folder))),
              bbox_inches='tight',
              #                transparent=True,
              pad_inches=0)