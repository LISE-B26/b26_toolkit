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


    This file contains functions to postprocess esr data acquired by looping over several NVs


"""
import numpy as np
import pandas as pd
import glob, os
import matplotlib.pyplot as plt
from matplotlib import gridspec
import time
import datetime


from PyLabControl.src.core import Script
from b26_toolkit.src.data_processing.esr_signal_processing import fit_esr, find_nv_peaks
from b26_toolkit.src.plotting.plots_1d import plot_esr
from PyLabControl.src.core.helper_functions import datetime_from_str


from b26_toolkit.src.scripts.find_nv import FindNV

from win32com.client import Dispatch
import pythoncom

# stuff for manual correction
import Queue
import ipywidgets as widgets
from IPython.display import display, clear_output
import threading

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
    next_queue = Queue.Queue()

    #set up manual correction gui
    widget_list = process_manual_esrs(nv_type_manual, b_field_manual, next_queue)

    # launch manual_correction function on separate thread, as otherwise buttons will not respond until that
    # function ends, that function is waiting on buttons to continue, and you deadlock
    thread = threading.Thread(target=manual_correction, args=(source_folder, target_folder, fit_data_set, nv_type_manual, b_field_manual, next_queue, widget_list[4], widget_list[5]))
    thread.start()

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

        findnv_folder = glob.glob(folder + '/data_subscripts/*find_nv*pt_*{:02d}'.format(pt_id))[0]

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
        fit_data_set_single = pd.DataFrame.from_dict({k: [v] for k, v in fit_data_set_single.iteritems()})

        if fit_data_set is None:
            fit_data_set = pd.DataFrame(fit_data_set_single)
        else:
            fit_data_set = fit_data_set.append(fit_data_set_single, ignore_index=True)

    return fit_data_set

def manual_correction(folder, target_folder, fit_data_set, nv_type_manual, b_field_manual, queue, lower_peak_widget, upper_peak_widget):
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
    try:
        fit_data_set_array = fit_data_set.as_matrix()

        w = widgets.HTML("Event information appears here when you click on the figure")
        display(w)

        # loop over all the folders in the data_subscripts subfolder and retrieve fitparameters and position of NV
        esr_folders = glob.glob(os.path.join(folder, '.\\data_subscripts\\*esr*'))

        # create folder to save images to
        # filepath_image = os.path.join(target_folder, os.path.dirname(folder).split('./')[1])
        # image_folder = os.path.join(filepath_image, '{:s}\\images'.format(os.path.basename(folder)))
        image_folder = os.path.join(target_folder, '{:s}\\images'.format(folder[2:]))
        print('image_folder', image_folder)
        # image_folder = os.path.normpath(
        #     os.path.abspath(os.path.join(os.path.join(target_folder, 'images'), os.path.basename(folders[0]))))
        if not os.path.exists(image_folder):
            os.makedirs(image_folder)
        if not os.path.exists(os.path.join(image_folder, 'bad_data')):
            os.makedirs(os.path.join(image_folder, 'bad_data'))

        f = plt.figure(figsize=(12, 6))

        def onclick(event):
            if event.button == 1:
                lower_peak_widget.value = event.xdata
            elif event.button == 3:
                upper_peak_widget.value = event.xdata

        cid = f.canvas.mpl_connect('button_press_event', onclick)

        for i, esr_folder in enumerate(esr_folders):

            # find the NV index
            pt_id = int(os.path.basename(esr_folder).split('pt_')[-1])

            findnv_folder = glob.glob(folder + '\\data_subscripts\\*find_nv*pt_{:02d}'.format(pt_id))[0]

            f.clf()
            gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
            ax0 = plt.subplot(gs[0])
            ax1 = plt.subplot(gs[1])
            ax = [ax0, ax1]

            # load data
            data = Script.load_data(esr_folder)
            fit_params = fit_data_set_array[pt_id, 2:8]
            if np.isnan(fit_params[4]):
                fit_params = fit_params[0:4]

            # get nv positions
            #             data_pos = {'initial_point': [fit_data_set['xo'].values[pt_id]]}
            data_pos = Script.load_data(findnv_folder)
            #             pos = data_pos['maximum_point']
            #             pos_init = data_pos['initial_point']


            # plot NV image
            FindNV.plot_data([ax[1]], data_pos)

            # plot data and fits
            print("fit_params: ", fit_params)
            plot_esr(ax[0], data['frequency'], data['data'], fit_params=fit_params)

            plt.tight_layout()

            plt.draw()
            plt.show()

            print('showing image')

            while queue.empty():
                time.sleep(.5)

            if nv_type_manual == '':
                if fit_params is None:
                    f.savefig(os.path.join(os.path.join(image_folder, 'bad_data'), 'esr_pt_{:02d}.jpg'.format(pt_id)))
                else:
                    f.savefig(os.path.join(image_folder, 'esr_pt_{:02d}.jpg'.format(pt_id)))
            else:
                if nv_type_manual[i] in ['bad', 'no_split']:
                    f.savefig(os.path.join(os.path.join(image_folder, 'bad_data'), 'esr_pt_{:02d}.jpg'.format(pt_id)))
                else:
                    f.savefig(os.path.join(image_folder, 'esr_pt_{:02d}.jpg'.format(pt_id)))

            queue.get()

        f.canvas.mpl_disconnect(cid)
        fit_data_set['manual_B_field'] = b_field_manual
        fit_data_set['manual_nv_type'] = nv_type_manual

        filename = '{:s}\\data-manual.csv'.format(os.path.basename(folder))
        filepath = os.path.join(target_folder, os.path.dirname(folder).split('./')[1])
        data_filepath = os.path.join(filepath, filename)

        print('to_save path: ', os.path.join(filepath, filename))


        if not os.path.exists(filepath):
            os.makedirs(filepath)

        fit_data_set.to_csv(data_filepath)

        create_shortcut(os.path.abspath(os.path.join(filepath, '{:s}\\to_data.lnk'.format(os.path.basename(folder)))), os.path.abspath(folder))
        create_shortcut(os.path.join(os.path.abspath(folder), 'to_processed.lnk'), os.path.abspath(os.path.join(filepath, '{:s}'.format(os.path.basename(folder)))))


        print('DONE!')
    except Exception as e:
        print('FAILED')
        print e

def process_manual_esrs(nv_type_manual, b_field_manual, next_queue):
    """
    Sets up the gui for manual correction.
    Args:
        nv_type_manual: pointer to array containing manual corrections for nv_type
        b_field_manual: pointer to array containing manual corrections for b_field
        next_queue:

    Returns:

    """
    # define queue for inter-button communication
    current_id_queue = Queue.Queue()
    current_id_queue.put(0)

    # define buttons to be added to display
    button_correct = widgets.Button(description="correct")
    button_bad = widgets.Button(description="bad")
    button_no_peak = widgets.Button(description="no_peak")
    button_peak = widgets.Button(description="peak")
    lower_peak = widgets.FloatText(description='lower (single) peak')
    upper_peak = widgets.FloatText(description='upper peak')
    widget_list = [button_correct, button_bad, button_no_peak, button_peak, lower_peak, upper_peak]

    # display all widgets
    display(button_correct)
    display(button_bad)
    display(button_no_peak)
    display(button_peak)
    display(lower_peak)
    display(upper_peak)


    def button_correct_clicked(b):
        current_id_queue.put(current_id_queue.get() + 1)
        next_queue.put(0)

    button_correct.on_click(button_correct_clicked)

    def button_bad_clicked(b):
        current_id = current_id_queue.get()
        nv_type_manual[current_id] = 'bad'
        current_id_queue.put(current_id + 1)
        next_queue.put(0)

    button_bad.on_click(button_bad_clicked)

    def button_no_split_clicked(b):
        current_id = current_id_queue.get()
        nv_type_manual[current_id] = 'no_split'
        current_id_queue.put(current_id + 1)
        next_queue.put(0)

    button_no_peak.on_click(button_no_split_clicked)

    def button_peak_clicked(b):
        current_id = current_id_queue.get()
        if upper_peak.value == 0:
            if np.abs(lower_peak.value - 2.87e9) > 0.01e9:
                nv_type_manual[current_id] = 'split'
                b_field_manual[current_id] = (np.abs(lower_peak.value - 2.87e9) * 2.) / (5.6e6)
                lower_peak.value = 0
            else:
                nv_type_manual[current_id] = 'single'
        else:
            nv_type_manual[current_id] = 'split'
            b_field_manual[current_id] = (upper_peak.value - lower_peak.value) / (5.6e6)
        lower_peak.value = 2.87e9
        upper_peak.value = 0
        current_id_queue.put(current_id + 1)
        next_queue.put(0)

    button_peak.on_click(button_peak_clicked)

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