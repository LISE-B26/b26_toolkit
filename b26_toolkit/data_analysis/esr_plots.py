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

import glob
import pandas as pd
import numpy as np

import ipywidgets as widgets
from IPython.display import display

import matplotlib.pyplot as plt
from matplotlib import gridspec

import os
import datetime

from pylabcontrol.core.helper_functions import datetime_from_str
from pylabcontrol.core.scripts import Script

freq_to_mag = 1. / (2 * 2.8e6)

# THIS IS A TEST JG
print('THIS IS A TEST JG sadsada')
def get_select_points(folder):
    """
    Finds the filepath to the data from the select_points script taken immediately previous to the data in folder.
    Args:
        folder: reference folder from which to find the closest select_points script

    Returns: folder of select_points script data

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

def distance_point_to_line(pt, line_pt1, line_pt2):
    dist = np.abs((line_pt2[1]-line_pt1[1])*pt[0] - (line_pt2[0]-line_pt1[0])*pt[1] + (line_pt2[0] * line_pt1[1]) - (line_pt2[1] * line_pt1[0]))\
           /np.sqrt((line_pt2[1] - line_pt1[1])**2 + (line_pt2[0] - line_pt1[0])**2)
    dist *= np.sign((line_pt1[0] - line_pt2[0]) * (pt[1] - line_pt2[1]) - (line_pt1[1] - line_pt2[1]) * (pt[0] - line_pt2[0]))
    return(dist)

def visualize_magnetic_fields(folder, manual=True, legend_location = 'upper right', bbox_to_anchor = (1,1), scatter_plot_axis = 'x', line_points = None):
    """
    Creates, plots, and saves two plots of NV data. The first is a plot of the fluoresence image, with a tag on each NV
    indicating whether it is split or unsplit. The second is a plot of the magnitude of the splitting and associated
    magnetic field vs the NV x coordinate.

    Args:
        folder: folder of original data
        manual: True if manual processing has been done first, such as by process_esrs in
                b26_toolkit.pylabcontrol.data_analysis.esr_post_processing. False if only the raw data should be used.
        legend_location: location of legend in plot, using standard matplotlib location tags. Use this argument to move
                         the legend if it covers the NVs.

    """
    # load the fit_data
    if manual:
        df = pd.read_csv(os.path.join('./==processed==', '{:s}/data-manual.csv'.format(folder)), index_col='id',
                         skipinitialspace=True)
    else:
        df = pd.read_csv(os.path.join(folder, '{:s}.csv'.format(folder.split('./')[1].replace('/', '-'))),
                         index_col='id', skipinitialspace=True)
    df = df.drop('Unnamed: 0', 1)

    # include manual corrections
    if manual:
        for id, nv_type in enumerate(df['manual_nv_type']):
            if not pd.isnull(nv_type):
                df.set_value(id, 'NV type', nv_type)
                if nv_type == 'split':
                    df.set_value(id, 'B-field (gauss)', df.get_value(id, 'manual_B_field'))

    # load the image data
    select_points_data = Script.load_data(get_select_points(folder))
    image_data = select_points_data['image_data']
    #     points_data = select_points_data['nv_locations']
    extent = select_points_data['extent']

    pt1x = widgets.FloatText(description='upper point x')
    pt1y = widgets.FloatText(description='upper point y')
    pt2x = widgets.FloatText(description='lower point x')
    pt2y = widgets.FloatText(description='lower point y')

    display(pt1x)
    display(pt1y)
    display(pt2x)
    display(pt2y)

    # prepare figure
    f = plt.figure(figsize=(15, 8))
    gs = gridspec.GridSpec(1, 2, height_ratios=[1, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    def onclick(event):
        if event.button == 1:
            pt1x.value = event.xdata
            pt1y.value = event.ydata
        elif event.button == 3:
            pt2x.value = event.xdata
            pt2y.value = event.ydata

    cid = f.canvas.mpl_connect('button_press_event', onclick)

    # plot the map
    ax0.imshow(image_data, extent=extent, interpolation='nearest', cmap='pink')

    for nv_type, color in zip(['split', 'bad', 'no_peak', 'single'], ['g', 'k', 'b', 'r']):
        subset = df[df['NV type'] == nv_type]
        ax0.scatter(subset['xo'], subset['yo'], c=color, label=nv_type)

    ax0.legend(loc=legend_location, bbox_to_anchor = bbox_to_anchor)

    subset = df[df['NV type'] == 'split']
    for i, j, n in zip(subset['xo'], subset['yo'], subset.index):
        corr = +0.005  # adds a little correction to put annotation in marker's centrum
        ax0.annotate(str(n), xy=(i + corr, j + corr), fontsize=4, color='w')
    ax0.set_xlabel('x (V)')
    ax0.set_ylabel('y (V)')

    # plot the fields on a 1D plot
    subset = df[df['NV type'] == 'split']
    if scatter_plot_axis == 'x':
        x_coor = subset['xo']
    elif scatter_plot_axis == 'y':
        x_coor = subset['yo']
    elif scatter_plot_axis == 'line':
        pt1 = line_points[0]
        pt2 = line_points[1]
        x_coor = []
        for x,y in zip(subset['xo'], subset['yo']):
            x_coor.append(distance_point_to_line([x,y], pt1, pt2))
    ax1.plot(x_coor, subset['B-field (gauss)'] / freq_to_mag * 1e-6, 'go')

    for i, j, n in zip(x_coor, subset['B-field (gauss)'] / freq_to_mag * 1e-6, subset.index):
        corr = 0.0005  # adds a little correction to put annotation in marker's centrum
        ax1.annotate(str(n), xy=(i + corr, j + corr))

    subset = df[df['NV type'] == 'single']
    for i, j, n in zip(subset['xo'], subset['yo'], subset.index):
        corr = +0.005  # adds a little correction to put annotation in marker's centrum
        ax0.annotate(str(n), xy=(i + corr, j + corr), fontsize=4, color='w')

    subset = df[df['NV type'] == 'single']
    subset_dict = {'xo': [], 'yo': [], 'B-field (gauss)': [], 'index': []}
    for nv in subset.iterrows():
        if 'fit_4' in nv[1]:
            subset_dict['xo'].append(nv[1][9])
            subset_dict['yo'].append(nv[1][11])
        else:
            subset_dict['xo'].append(nv[1][7])
            subset_dict['yo'].append(nv[1][9])
        subset_dict['index'].append(nv[0])
        if np.isnan(nv[1][6]) and ((nv[1][4]) > 2.92e9 or (nv[1][4]) < 2.82e9):
            subset_dict['B-field (gauss)'].append((nv[1][4] - 2.87e9) * 2 * freq_to_mag)
        else:
            subset_dict['B-field (gauss)'].append(0)

    if scatter_plot_axis == 'x':
        print('HERE')
        x_coor = subset_dict['xo']
    elif scatter_plot_axis == 'y':
        x_coor = subset_dict['yo']
    elif scatter_plot_axis == 'line':
        pt1 = line_points[0]
        pt2 = line_points[1]
        x_coor = []
        for x, y in zip(subset_dict['xo'], subset_dict['yo']):
            x_coor.append(distance_point_to_line([x, y], pt1, pt2))
    ax1.plot(x_coor, np.array(subset_dict['B-field (gauss)']) / freq_to_mag * 1e-6, 'ro')

    for i, j, n in zip(x_coor, np.array(subset_dict['B-field (gauss)']) / freq_to_mag * 1e-6, subset_dict['index']):
        corr = 0.0005  # adds a little correction to put annotation in marker's centrum
        ax1.annotate(str(n), xy=(i + corr, j + corr))

    ax1.set_title('ESR Splittings')
    if scatter_plot_axis == 'x':
        ax1.set_xlabel('x coordinate (V)')
    if scatter_plot_axis == 'y':
        ax1.set_xlabel('y coordinate (V)')
    if scatter_plot_axis == 'line':
        ax1.set_xlabel('Distance to magnet edge (V)')
    ax1.set_ylabel('Splitting (MHz)')
    if not scatter_plot_axis == 'line':
        if scatter_plot_axis == 'y':
            ax1.set_xlim([extent[3], extent[2]])
        else:
            ax1.set_xlim([extent[0], extent[1]])

    plt.axvline(x=0, color = 'black', linestyle = 'dashed', linewidth = 2)

    ax2 = ax1.twinx()
    mn, mx = ax1.get_ylim()
    ax2.set_ylim(mn * freq_to_mag * 1e6, mx * freq_to_mag * 1e6)
    ax2.set_ylabel('Magnetic Field Projection (Gauss)')


    # f.set_tight_layout(True)

    print(('path', os.path.join('./==processed==', folder[2:], 'splittings_plot.jpg'.format(os.path.basename(folder)))))

    f.savefig(os.path.join('./==processed==', folder[2:], 'splittings_plot.jpg'.format(os.path.basename(folder))),
              bbox_inches='tight',
              pad_inches=0)


def visualize_magnetic_fields_comparison(folders, labels, manual=True, legend_location = 'upper right', bbox_to_anchor = (1,1), scatter_plot_axis = 'x', line_points_array = None, transformation_matrix = None):
    """
    Creates, plots, and saves two plots of NV data. The first is a plot of the fluoresence image, with a tag on each NV
    indicating whether it is split or unsplit. The second is a plot of the magnitude of the splitting and associated
    magnetic field vs the NV x coordinate.

    Args:
        folder: folder of original data
        manual: True if manual processing has been done first, such as by process_esrs in
                b26_toolkit.pylabcontrol.data_analysis.esr_post_processing. False if only the raw data should be used.
        legend_location: location of legend in plot, using standard matplotlib location tags. Use this argument to move
                         the legend if it covers the NVs.

    """
    # load the image data
    select_points_data = Script.load_data(get_select_points(folders[0]))
    image_data = select_points_data['image_data']
    #     points_data = select_points_data['nv_locations']
    extent = select_points_data['extent']

    pt1x = widgets.FloatText(description='upper point x')
    pt1y = widgets.FloatText(description='upper point y')
    pt2x = widgets.FloatText(description='lower point x')
    pt2y = widgets.FloatText(description='lower point y')

    display(pt1x)
    display(pt1y)
    display(pt2x)
    display(pt2y)

    # prepare figure
    f = plt.figure(figsize=(15, 8))
    gs = gridspec.GridSpec(1, 1)
    # ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[0])

    def onclick(event):
        if event.button == 1:
            pt1x.value = event.xdata
            pt1y.value = event.ydata
        elif event.button == 3:
            pt2x.value = event.xdata
            pt2y.value = event.ydata

    cid = f.canvas.mpl_connect('button_press_event', onclick)

    # plot the map
    # ax0.imshow(image_data, extent=extent, interpolation='nearest', cmap='pink')

    if line_points_array is None:
        line_points_array = np.repeat(0,len(folders))

    folder_num = 0
    for folder, label, line_points in zip(folders, labels, line_points_array):
        # load the fit_data
        if manual:
            df = pd.read_csv(os.path.join('./==processed==', '{:s}/data-manual.csv'.format(folder)), index_col='id',
                             skipinitialspace=True)
        else:
            df = pd.read_csv(os.path.join(folder, '{:s}.csv'.format(folder.split('./')[1].replace('/', '-'))),
                             index_col='id', skipinitialspace=True)
        df = df.drop('Unnamed: 0', 1)

        # include manual corrections
        if manual:
            for id, nv_type in enumerate(df['manual_nv_type']):
                if not pd.isnull(nv_type):
                    df.set_value(id, 'NV type', nv_type)
                    if nv_type == 'split':
                        df.set_value(id, 'B-field (gauss)', df.get_value(id, 'manual_B_field'))

        # for nv_type, color in zip(['split', 'bad', 'no_peak', 'single'], ['g', 'k', 'b', 'r']):
        #     subset = df[df['NV type'] == nv_type]
        #     ax0.scatter(subset['xo'], subset['yo'], c=color, label=nv_type)
        #
        # ax0.legend(loc=legend_location, bbox_to_anchor = bbox_to_anchor)
        #
        # subset = df[df['NV type'] == 'split']
        # for i, j, n in zip(subset['xo'], subset['yo'], subset.index):
        #     corr = +0.005  # adds a little correction to put annotation in marker's centrum
        #     ax0.annotate(str(n), xy=(i + corr, j + corr), fontsize=8, color='w')
        # ax0.set_xlabel('x (V)')
        # ax0.set_ylabel('y (V)')

        subset = df[df['NV type'] == 'split']

        if transformation_matrix is None or folder_num == 0:
            nv_xs = subset['xo']
            nv_ys = subset['yo']
        else:
            if folder_num >= 0:
                nv_coors = list()
                coordinates = list(zip(subset['xo'], subset['yo']))
                for c in coordinates:
                    c3 = [c[0], c[1], 1]
                    c3_trans = np.dot(transformation_matrix,c3)
                    nv_coors.append([c3_trans[0], c3_trans[1]])
                    nv_xs, nv_ys = list(zip(*nv_coors))

        # plot the fields on a 1D plot
        if scatter_plot_axis == 'x':
            x_coor = nv_xs
        elif scatter_plot_axis == 'y':
            x_coor = nv_ys
        elif scatter_plot_axis == 'line':
            pt1 = line_points[0]
            pt2 = line_points[1]
            x_coor = []
            for x,y in zip(nv_xs, nv_ys):
                x_coor.append(distance_point_to_line([x,y], pt1, pt2))
        prev_plot = ax1.plot(x_coor, subset['B-field (gauss)'] / freq_to_mag * 1e-6, 'o', label = label)

        for i, j, n in zip(x_coor, subset['B-field (gauss)'] / freq_to_mag * 1e-6, subset.index):
            corr = 0.0005  # adds a little correction to put annotation in marker's centrum
            ax1.annotate(str(n), xy=(i + corr, j + corr))

        # subset = df[df['NV type'] == 'single']
        # for i, j, n in zip(subset['xo'], subset['yo'], subset.index):
        #     corr = +0.005  # adds a little correction to put annotation in marker's centrum
        #     ax0.annotate(str(n), xy=(i + corr, j + corr), fontsize=8, color='w')

        subset = df[df['NV type'] == 'single']
        subset_dict = {'xo': [], 'yo': [], 'B-field (gauss)': [], 'index': []}
        for nv in subset.iterrows():
            if 'fit_4' in nv[1]:
                subset_dict['xo'].append(nv[1][9])
                subset_dict['yo'].append(nv[1][11])
            else:
                subset_dict['xo'].append(nv[1][7])
                subset_dict['yo'].append(nv[1][9])
            subset_dict['index'].append(nv[0])
            if np.isnan(nv[1][6]) and ((nv[1][4]) > 2.92e9 or (nv[1][4]) < 2.82e9):
                subset_dict['B-field (gauss)'].append((nv[1][4] - 2.87e9) * 2 * freq_to_mag)
            else:
                subset_dict['B-field (gauss)'].append(0)

        if scatter_plot_axis == 'x':
            x_coor = subset_dict['xo']
        elif scatter_plot_axis == 'y':
            x_coor = subset_dict['yo']
        elif scatter_plot_axis == 'line':
            pt1 = line_points[0]
            pt2 = line_points[1]
            x_coor = []
            for x, y in zip(subset_dict['xo'], subset_dict['yo']):
                x_coor.append(distance_point_to_line([x, y], pt1, pt2))
        ax1.plot(x_coor, np.array(subset_dict['B-field (gauss)']) / freq_to_mag * 1e-6, 'o', color = prev_plot[0].get_color())


        for i, j, n in zip(x_coor, np.array(subset_dict['B-field (gauss)']) / freq_to_mag * 1e-6, subset_dict['index']):
            corr = 0.0005  # adds a little correction to put annotation in marker's centrum
            ax1.annotate(str(n), xy=(i + corr, j + corr))

        folder_num += 1

    ax1.set_title('ESR Splittings')
    if scatter_plot_axis == 'x':
        ax1.set_xlabel('x coordinate (V)')
    if scatter_plot_axis == 'y':
        ax1.set_xlabel('y coordinate (V)')
    if scatter_plot_axis == 'line':
        ax1.set_xlabel('Distance to magnet edge (V)')
    ax1.set_ylabel('Splitting (MHz)')
    if not scatter_plot_axis == 'line':
        if scatter_plot_axis == 'y':
            ax1.set_xlim([extent[3], extent[2]])
        else:
            ax1.set_xlim([extent[0], extent[1]])

    plt.axvline(x=0, color = 'black', linestyle = 'dashed', linewidth = 2)

    ax2 = ax1.twinx()
    mn, mx = ax1.get_ylim()
    ax2.set_ylim(mn * freq_to_mag * 1e6, mx * freq_to_mag * 1e6)
    ax2.set_ylabel('Magnetic Field Projection (Gauss)')

    ax1.legend()

    # f.set_tight_layout(True)

    print(('path', os.path.join('./==processed==', folder[2:], 'splittings_plot_comparison.jpg'.format(os.path.basename(folder)))))

    f.savefig(os.path.join('./==processed==', folder[2:], 'splittings_plot_comparison.jpg'.format(os.path.basename(folder))),
              bbox_inches='tight',
              pad_inches=0)