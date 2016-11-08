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

from PyLabControl.src.core import Script
from b26_toolkit.src.data_processing.esr_signal_processing import fit_esr, find_nv_peaks
from b26_toolkit.src.plotting.plots_1d import plot_esr

from b26_toolkit.src.scripts.find_nv import FindNV

# stuff for manual correction
import Queue
import ipywidgets as widgets
from IPython.display import display
from IPython.display import clear_output
import threading


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


def process_esrs(folders, target_folder):
    """

    fits the esr data, plots them and asks the user for confirmation, the fit data is saved to the folder target_folder with the same structure as folders



    Args:
        folders: list of source folders with esr data, this folder shoudl contain a subfolder data_subscripts which contains subfolders *esr* with the esr data
        target_folder: target folder where the output data is saved in form of a .csv file

    Returns: fitdataset as a pandas array

    """
    for folder in folders:
        # loop over all the folders in the data_subscripts subfolder and retrieve fitparameters and position of NV
        esr_folders = glob.glob(os.path.join(folder, './data_subscripts/*esr*'))

        fit_data_set = None

        # classify the nvs according to the following categories
        # by default we set this to na (not available) and try to figure it out based on the data and fitquality
        # nv_type = 'na' # split / single / no_peak / no_nv / na

        for i, esr_folder in enumerate(esr_folders):

            # find the NV index
            pt_id = int(os.path.basename(esr_folder).split('pt_')[-1])

            findnv_folder = glob.glob(folder + '/data_subscripts/*find_nv*pt_{:02d}'.format(pt_id))[0]

            # prepare figure
            f = plt.figure(figsize=(12, 6))
            gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
            ax0 = plt.subplot(gs[0])
            ax1 = plt.subplot(gs[1])
            ax = [ax0, ax1]

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

            # plot NV image
            FindNV.plot_data([ax[1]], data_pos)

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

            # plot data and fits
            plot_esr(ax[0], data['frequency'], data['data'], fit_params=fit_params)
            # plot initial fit guesses
            ax[0].hold(True)
            for freq, ampl in zip(freq_peaks, ampl_peaks):
                ax[0].plot([freq, freq], [min(data['data']), max(data['data'])], 'k--')

            plt.tight_layout()
            plt.close('all')
        # save data set
        filename = '{:s}.csv'.format(os.path.basename(folder))
        filepath = os.path.join(target_folder, os.path.dirname(folder).split('./')[1])
        print('saving to {:s}'.format(filename))

        if not os.path.exists(filepath):
            os.makedirs(filepath)

        fit_data_set.to_csv(os.path.join(filepath, filename))

    return fit_data_set


def manual_correction(folders, fit_data_set, queue, lower_peak_widget, upper_peak_widget, widgets, display, target_folder):


    nv_type_manual = [''] * len(fit_data_set)
    b_field_manual = [np.nan] * len(fit_data_set)


    fit_data_set_array = fit_data_set.as_matrix()

    w = widgets.HTML("Event information appears here when you click on the figure")
    display(w)

    for folder in folders:
        # loop over all the folders in the data_subscripts subfolder and retrieve fitparameters and position of NV
        esr_folders = glob.glob(os.path.join(folder, './data_subscripts/*esr*'))

        # create folder to save images to
        image_folder = os.path.normpath(
            os.path.abspath(os.path.join(os.path.join(target_folder, 'images'), os.path.basename(folders[0]))))
        if not os.path.exists(image_folder):
            os.makedirs(image_folder)
        if not os.path.exists(os.path.join(image_folder, 'bad_data')):
            os.makedirs(os.path.join(image_folder, 'bad_data'))

        f = plt.figure(figsize=(12, 6))

        #         gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
        #         ax0 = plt.subplot(gs[0])
        #         ax1 = plt.subplot(gs[1])
        #         ax = [ax0, ax1]

        def onclick(event):
            #             w.value = 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
            #                       event.button, event.x, event.y, event.xdata, event.ydata)
            if event.button == 1:
                lower_peak_widget.value = event.xdata
            elif event.button == 3:
                upper_peak_widget.value = event.xdata

        cid = f.canvas.mpl_connect('button_press_event', onclick)

        for i, esr_folder in enumerate(esr_folders):

            # find the NV index
            pt_id = int(os.path.basename(esr_folder).split('pt_')[-1])

            findnv_folder = glob.glob(folder + '/data_subscripts/*find_nv*pt_{:02d}'.format(pt_id))[0]

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
            nv_type = fit_data_set_array[pt_id, 1]

            # get nv positions
            #             data_pos = {'initial_point': [fit_data_set['xo'].values[pt_id]]}
            data_pos = Script.load_data(findnv_folder)
            #             pos = data_pos['maximum_point']
            #             pos_init = data_pos['initial_point']


            # plot NV image
            FindNV.plot_data([ax[1]], data_pos)

            # plot data and fits
            plot_esr(ax[0], data['frequency'], data['data'], fit_params=fit_params)

            plt.tight_layout()

            plt.show()

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

        filename = '{:s}-manual.csv'.format(os.path.basename(folder))
        filepath = os.path.join(target_folder, os.path.dirname(folder).split('./')[1])

        if not os.path.exists(filepath):
            os.makedirs(filepath)

        fit_data_set.to_csv(os.path.join(filepath, filename))

        print('DONE!')


def process_manual_esrs(folders, fit_data_set):
    # define queues for inter-thread communication
    next_queue = Queue.Queue()
    current_id_queue = Queue.Queue()
    current_id_queue.put(0)

    # define buttons to be added to display
    button_correct = widgets.Button(description="correct")
    button_bad = widgets.Button(description="bad")
    button_no_peak = widgets.Button(description="no_peak")
    button_peak = widgets.Button(description="peak")
    lower_peak = widgets.FloatText(description='lower (single) peak')
    upper_peak = widgets.FloatText(description='upper peak')

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

    # launch manual_correction function on separate thread, as otherwise buttons will not respond until that
    # function ends, that function is waiting on buttons to continue, and you deadlock
    thread = threading.Thread(target=manual_correction,args=(folders, next_queue, lower_peak, upper_peak))
    thread.start()