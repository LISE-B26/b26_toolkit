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


    This file contains functions to plot previously processed ESR data

"""

import glob
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import gridspec

import os
import datetime

from PyLabControl.src.core.helper_functions import datetime_from_str
from PyLabControl.src.core.scripts import Script

freq_to_mag = 1. / (2 * 2.8e6)


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

def visualize_magnetic_fields(folder, manual=True, legend_location = 'upper right'):
    """
    Creates, plots, and saves two plots of NV data. The first is a plot of the fluoresence image, with a tag on each NV
    indicating whether it is split or unsplit. The second is a plot of the magnitude of the splitting and associated
    magnetic field vs the NV x coordinate.

    Args:
        folder: folder of original data
        manual: True if manual processing has been done first, such as by process_esrs in
                b26_toolkit.src.data_analysis.esr_post_processing. False if only the raw data should be used.
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

    ax0.legend(loc=legend_location)

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

    print('path', os.path.join('./==processed==', folder[2:], 'splittings_plot.jpg'.format(os.path.basename(folder))))

    f.savefig(os.path.join('./==processed==', folder[2:], 'splittings_plot.jpg'.format(os.path.basename(folder))),
              bbox_inches='tight',
              pad_inches=0)