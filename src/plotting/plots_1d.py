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

import matplotlib.patches as patches
import numpy as np
from matplotlib.collections import PatchCollection

from b26_toolkit.src.data_processing.fit_functions import lorentzian, double_lorentzian


def plot_psd(freq, psd, axes, clear = True):
    '''
    plots the power spectral density on to the canvas axes
    :param freq: x-values array of length N
    :param psd: y-values array of length N
    :param axes: target axes object
    :return: None
    '''
    if clear is True:
        print('CLEAR')
        axes.clear()

    if np.mean(freq) > 1e6:
        freq /= 1e6
        unit = 'MHz'
    elif np.mean(freq) > 1e3:
        freq /= 1e3
        unit = 'kHz'

    axes.semilogy(freq, psd)

    axes.set_xlabel('frequency ({:s})'.format(unit))


def plot_esr(axes, frequency, counts, fit_params=None, plot_marker_data = 'b', plot_marker_fit = 'r'):
    """
    plots the esr
    Args:
        axes: axes object
        fit_params: array with fitparameters either length 4 for single peak or length 6 for double peak
        frequency: mw frequency (array)
        counts: counts (array)
        plot_marker:  (str) plot marker

    Returns:

    """

    #  ======== plot data =========
    axes.plot(frequency, counts, plot_marker_data)
    axes.hold(True)

    title = 'ESR'
    fit_data = None

    #  ======== plot fit =========
    if fit_params is not None and fit_params[0] != -1:  # check if fit valid
        if len(fit_params) == 4:
            # single peak
            fit_data = lorentzian(frequency, *fit_params)
            title = 'ESR fo = {:0.4e}, wo = {:0.2e}'.format(fit_params[2], fit_params[1])
        elif len(fit_params) == 6:
            # double peak
            fit_data = double_lorentzian(frequency, *fit_params)
            title = 'ESR f1 = {:0.4e} Hz, f2 = {:0.4e} Hz, wo = {:0.2e}'.format(fit_params[4], fit_params[5],
                                                                                fit_params[1])

    if fit_data is not None:
        axes.plot(frequency, fit_data, plot_marker_fit)

    # if fit_data is not None:  # plot esr and fit data
    #     lines = axes.plot(frequency, data, 'b', frequency, fit_data, 'r')
    # else:  # plot just esr data
    #     lines = axes.plot(frequency, data, 'b')

    axes.set_title(title)
    axes.set_xlabel('Frequency (Hz)')
    axes.set_ylabel('Kcounts/s')
    axes.hold(False)
    # return lines


def plot_pulses(axis, pulse_collection, pulse_colors=None):
    """
    creates a visualization of pulses (in pulse_collection) on a matplotlib axis (axis)

    Args:
        axis: The axis for the matplotlib plot
        pulse_collection: a collection of pulses, named tuples (channel_id, start_time, duration)
        pulse_colors: a dictionary of {channel_id:matplotlib_color} that maps channels to colors

    Returns:

    """

    # create a list of unique instruments from the pulses
    instrument_names = sorted(list(set([pulse.channel_id for pulse in pulse_collection])))

    # assign colors for certain specific channels
    if pulse_colors is None:
        pulse_colors = {'laser': '#50FF00', 'microwave_i': 'r', 'apd_readout': 'k'}

    # find the maximum time from the list of pulses
    max_time = max([pulse.start_time + pulse.duration for pulse in pulse_collection])

    # set axis boundaries
    axis.set_ylim(-0.75, len(instrument_names) - .25)
    axis.set_xlim(0, max_time)

    # label y axis with pulse names
    axis.set_yticks(range(len(instrument_names)))
    axis.set_yticklabels(instrument_names)

    # create horizontal lines for each pulse
    for pulse_plot_y_position in range(0, len(instrument_names)):
        axis.axhline(pulse_plot_y_position - .25, 0.0, max_time, color='k')
    axis.tick_params(axis='y', which=u'both', length=0)  # remove tick marks on y axis

    # create a vertical line denoting the end of the pulse sequence loop
    # axis.axvline(max_time, -0.5, len(instrument_names), color='r')

    # create rectangles for the pulses
    patch_list = []
    for pulse in pulse_collection:
        patch_list.append(
            patches.Rectangle((pulse.start_time, instrument_names.index(pulse.channel_id) - .25), pulse.duration, 0.5,
                              fc=pulse_colors.get(pulse.channel_id, 'b')))

    patch_collection = PatchCollection(patch_list, match_original=True)
    axis.add_collection(patch_collection)

    # label the axis
    axis.set_title('Pulse Visualization')
    axis.set_xlabel('time [ns]')
    axis.set_ylabel('pulse destination')


def update_pulse_plot(axis, pulse_collection, pulse_colors=None):
    """
    updates a previously created plot of pulses, removing the previous ones and adding ones corresponding to
    pulse_collection. The new pulse collection must only contain channel_ids already present on the passed axis

    Args:
        axis: The axis for the matplotlib plot
        pulse_collection: a collection of pulses, named tuples (channel_id, start_time, duration)
        pulse_colors: a dictionary of {channel_id:matplotlib_color} that maps channels to colors

    Returns:

    """

    # assign colors for certain specific channels
    if pulse_colors is None:
        pulse_colors = {'laser': '#50FF00', 'microwave_i': 'r', 'apd_readout': 'k'}

    # get a list of unique instruments from the pulses
    instrument_names = [str(label.get_text()) for label in axis.get_yticklabels()]

    # find the maximum time from the list of pulses
    max_time = max([pulse.start_time + pulse.duration for pulse in pulse_collection])

    axis.set_xlim(0, max_time)

    # remove the previous pulses
    [child.remove() for child in axis.get_children() if isinstance(child, PatchCollection)]

    # create rectangles for the pulses
    patch_list = []
    for pulse in pulse_collection:
        patch_list.append(
            patches.Rectangle((pulse.start_time, instrument_names.index(pulse.channel_id) - .25), pulse.duration, 0.5,
                              fc=pulse_colors.get(pulse.channel_id, 'b')))

    patch_collection = PatchCollection(patch_list, match_original=True)
    axis.add_collection(patch_collection)

def plot_counts(axis, data):
    # COMMENT_ME
    axis.plot(data)
    axis.hold(False)

    axis.set_xlabel('time')
    axis.set_ylabel('kCounts/sec')


def plot_1d_simple_timetrace_ns(axis, times, data_list, y_label='kCounts/sec', title=None):
    """
    plots a time trace for a list of data assuming that the times are give in ns
    Args:
        axis: axis object on which to plot
        times: times in ns (list or array of length N)
        data_list: list of data (size MxN)
        y_label: (optional) label for y axis
        title:  (optional) title

    """
    for counts in data_list:
        axis.plot(times, counts)

    axis.hold(False)

    if max(times) < 1e3:
        x_label = 'time (ns)'
    elif max(times) < 1e6:
        x_label = 'time (us)'
        times *= 1e-3
    elif max(times) < 1e9:
        x_label = 'time (ms)'
        times *= 1e-6
    elif max(times) < 1e12:
        x_label = 'time (s)'
        times *= 1e-9

    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    if title:
        axis.set_title(title)
    axis.set_xlim([min(times), max(times)])


def update_1d_simple(axis, times, counts_list):
    """

    Args:
        axis: axes object
        times: JG: THIS IS NOT USED! WHAT IS IT? => add comment, e.g. for future purpose or delete!
        counts_list: list of

    Returns:

    """
    if len(axis.lines) != len(counts_list):
        counts_list = np.transpose(counts_list)

    assert len(axis.lines) == len(counts_list)
    for index, counts in enumerate(counts_list):
        axis.lines[index].set_ydata(counts)
    axis.relim()
    axis.autoscale_view()
