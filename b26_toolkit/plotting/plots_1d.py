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

import matplotlib.patches as patches
import numpy as np
from matplotlib.collections import PatchCollection

from b26_toolkit.data_processing.fit_functions import lorentzian, double_lorentzian

def plot_psd(freq, psd, axes, y_scaling = 'log', x_scaling = 'lin', color='b'):
    '''
    plots the power spectral density on to the canvas axes
    :param freq: x-values array of length N
    :param psd: y-values array of length N
    :param axes: target axes object
    :return: None
    '''
    unit = 'Hz'
    c_unit = 1.0
    if x_scaling == 'lin':
        if np.mean(freq) > 1e6:
            c_unit = 1e-6
            unit = 'MHz'
        elif np.mean(freq) > 1e3:
            c_unit = 1e-3
            unit = 'kHz'


    if y_scaling == 'log' and x_scaling == 'log':
        axes.loglog(c_unit*freq, psd, color=color)
    elif y_scaling == 'log' and x_scaling == 'lin':
        axes.semilogy(c_unit*freq, psd, color=color)
    elif y_scaling == 'lin' and x_scaling == 'log':
        axes.semilogx(c_unit*freq, psd, color=color)
    elif y_scaling == 'lin' and x_scaling == 'lin':
        axes.plot(c_unit*freq, psd, color=color)

    axes.set_xlabel('frequency ({:s})'.format(unit))

    axes.set_xlim([min(c_unit*freq), max(c_unit*freq)])

def plot_esr(axes, frequency, counts, fit_params=None, plot_marker_data = 'b', plot_marker_fit = 'r', linestyle = '-', marker = '.'):
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
    axes.clear() # ER 20181012 - matplotlib axes.hold() removed in update to 3.0.0
    if len(np.array(counts).shape) == 1:
        if np.shape(frequency) == np.shape(counts): # ER 20190129
            axes.plot(frequency, counts, plot_marker_data, linestyle = linestyle, marker = marker)
    elif len(np.array(counts).shape) == 2:
        for count_row in counts:
            if np.shape(frequency) == np.shape(count_row):  # ER 20190129
                axes.plot(frequency, count_row, plot_marker_data, linestyle=linestyle, marker=marker)

    #axes.hold(True) #ER 20181012

    title = 'ESR'
    fit_data = None

    #  ======== plot fit =========
    if fit_params is not None and len(fit_params) and fit_params[0] != -1:  # check if fit valid
        if len(fit_params) == 4:
            # single peak
            fit_data = lorentzian(frequency, *fit_params)
            title = 'ESR fo = {:0.4e}, wo = {:0.2e}, contrast0 = {:.2%}'.format(fit_params[2], fit_params[3],
                                                                                np.abs(fit_params[1]) / fit_params[0])
        elif len(fit_params) == 6:
            # double peak
            fit_data = double_lorentzian(frequency, *fit_params)
            title = 'ESR f1 = {:0.4e} Hz, f2 = {:0.4e} Hz, wo = {:0.2e},\n contrast1 = {:.2%}, contrast2 = {:.2%}'.format(fit_params[4], fit_params[5],
                                                                                fit_params[1],
                                                                                np.abs(fit_params[2]) / fit_params[0], np.abs(fit_params[3]) / fit_params[0])

    if fit_data is not None:
        axes.plot(frequency, fit_data, plot_marker_fit)

    # if fit_data is not None:  # plot esr and fit data
    #     lines = axes.plot(frequency, data, 'b', frequency, fit_data, 'r')
    # else:  # plot just esr data
    #     lines = axes.plot(frequency, data, 'b')

    axes.set_title(title)
    axes.set_xlabel('Frequency (Hz)')
    axes.set_ylabel('Kcounts/s')

def plot_pulsedesr(axes, frequency, fit_params=None, plot_marker_fit = 'r'):
    """
    plots the esr (simplified for pulsedESR)
    Args:
        axes: axes object
        fit_params: array with fitparameters either length 4 for single peak or length 6 for double peak
        frequency: mw frequency (array)
        counts: counts (array)
        plot_marker:  (str) plot marker

    Returns:

    """

    fit_data = None

    #  ======== plot fit =========
    if fit_params is not None and len(fit_params) and fit_params[0] != -1:  # check if fit valid
        if len(fit_params) == 4:
            # single peak
            fit_data = lorentzian(frequency, *fit_params)
            #title = 'ESR fo = {:0.4e}, wo = {:0.2e}, contrast0 = {:.2%}'.format(fit_params[2], fit_params[3],
            #                                                                    np.abs(fit_params[1]) / fit_params[0])
        elif len(fit_params) == 6:
            # double peak
            fit_data = double_lorentzian(frequency, *fit_params)
            print(fit_data)
            #title = 'ESR f1 = {:0.4e} Hz, f2 = {:0.4e} Hz, wo = {:0.2e},\n contrast1 = {:.2%}, contrast2 = {:.2%}'.format(fit_params[4], fit_params[5],
            #                                                                    fit_params[1],
            #                                                                    np.abs(fit_params[2]) / fit_params[0], np.abs(fit_params[3]) / fit_params[0])

    if fit_data is not None:
        axes.plot(frequency, fit_data, plot_marker_fit)

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
    axis.set_yticks(list(range(len(instrument_names))))
    axis.set_yticklabels(instrument_names)

    # remove all previous pulses
    [child.remove() for child in axis.get_children() if isinstance(child, PatchCollection)]

    # create horizontal lines for each pulse
    for pulse_plot_y_position in range(0, len(instrument_names)):
        axis.axhline(pulse_plot_y_position - .25, 0.0, max_time, color='k')
    axis.tick_params(axis='y', which='both', length=0)  # remove tick marks on y axis

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
    axis.set_ylabel('pulse destination')
    axis.set_xlabel('time [ns]')

    xticks = np.array(axis.get_xticks())

    if np.sum(xticks > 1E9) > 0:
        axis.set_xticklabels(xticks/1e9)
        axis.set_xlabel('time [s]')

    elif np.sum(xticks > 1E6) > 0:
        axis.set_xticklabels(xticks/1e6)
        axis.set_xlabel('time [ms]')

    elif np.sum(xticks > 1E3) > 0:
        axis.set_xticklabels(xticks/1e3)
        axis.set_xlabel('time [us]')


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
    instrument_names_old = [str(label.get_text()) for label in axis.get_yticklabels()]
    instrument_names = sorted(list(set([pulse.channel_id for pulse in pulse_collection])))

    # if the number of pulses has changed call plot instead of update
    if set(instrument_names_old) != set(instrument_names):
        plot_pulses(axis, pulse_collection, pulse_colors = pulse_colors)
    else:
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


        xticks = np.array(axis.get_xticks())

        if np.sum(xticks > 1E9) > 0:
            axis.set_xticklabels(xticks/1e9)
            axis.set_xlabel('time [s]')

        elif np.sum(xticks > 1E6) > 0:
            axis.set_xticklabels(xticks/1e6)
            axis.set_xlabel('time [ms]')

        elif np.sum(xticks > 1E3) > 0:
            axis.set_xticklabels(xticks/1e3)
            axis.set_xlabel('time [us]')


def plot_counts(axis, data):
    """
    plots APD timeseries data

    Args:
        axis: the axis to draw the plot
        data (2d array): APD count timeseries data

    Returns: (none)
    """
    if len(np.shape(data)) > 1:
        for datum in data:
            axis.plot(datum, linewidth=1.25)
    else:
        axis.plot(data, linewidth=1.25)
    # axis.hold(False)

    axis.set_xlabel('sample number')
    axis.set_ylabel('[kCounts/s]')

def update_counts(axis, data, data_size=1):
    if data is None:
        return
    if len(np.shape(data)) > 1:
        for i, datum in enumerate(data):
            axis.lines[i].set_ydata(datum)
            axis.lines[i].set_xdata(range(0, len(datum)))

    else:
        axis.lines[0].set_ydata(data)
        axis.lines[0].set_xdata(range(0, len(data)))
    axis.relim()
    axis.autoscale_view()

def plot_voltage(axis, data):
    """
    plots APD timeseries data

    Args:
        axis: the axis to draw the plot
        data (2d array): APD count timeseries data

    Returns: (none)
    """

    axis.plot(data, linewidth=2.0)
    # axis.hold(False)

    axis.set_xlabel('time')
    axis.set_ylabel('voltage (V)')


def plot_temperature(axis, data, sample_rate, update=False):
    """
    plots the temperature

    Args:
        axis:
        data:
        sample_rate: at which data has been acquired
    Returns:

    """

    time = np.arange(len(data))/float(sample_rate)

    label = 'time (s)'

    if max(time)>60:
        time /= 60.
        label = 'time (min)'
    if max(time)>60:
        time /= 60.
        label = 'time (h)'

    if update:
        axis.lines[0].set_ydata(data)
        axis.lines[0].set_xdata(time)

        axis.relim()
        axis.autoscale_view()
    else:
        axis.plot(time, data)

    axis.set_xlabel(label)
    axis.set_ylabel('temperature (K)')

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

    times = 1.*np.array(times) # cast onto numpy in case we got a list, which breaks some of the commands below

    if max(times) < 1e3:
        x_label = 'time [ns]'
    elif max(times) < 1e6:
        x_label = 'time [us]'
        times *= 1e-3
    elif max(times) < 1e9:
        x_label = 'time [ms]'
        times *= 1e-6
    elif max(times) < 1e12:
        x_label = 'time [s]'
        times *= 1e-9

    for counts in data_list:
        axis.plot(times, counts)

    # axis.hold(False)


    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    axis.grid(True, alpha=.15)
    if title:
        axis.set_title(title)
    axis.set_xlim([min(times), max(times)])

def plot_1d_simple_freq(axis, freqs, data_list, y_label='kCounts/sec', title=None):
    """
    plots a time trace for a list of data assuming that the times are give in ns
    Args:
        axis: axis object on which to plot
        times: times in ns (list or array of length N)
        data_list: list of data (size MxN)
        y_label: (optional) label for y axis
        title:  (optional) title

    """

    freqs = 1.*np.array(freqs) # cast onto numpy in case we got a list, which breaks some of the commands below

    if max(freqs) < 1e3:
        x_label = 'freq [Hz]'
    elif max(freqs) < 1e6:
        x_label = 'freq [kHz]'
        freqs *= 1e-3
    elif max(freqs) < 1e9:
        x_label = 'freq [MHz]'
        freqs *= 1e-6
    elif max(freqs) < 1e12:
        x_label = 'freq [GHz]'
        freqs *= 1e-9

    for counts in data_list:
        axis.plot(freqs, counts)

    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    axis.grid(True, alpha=.1)
    if title:
        axis.set_title(title)
    axis.set_xlim([min(freqs), max(freqs)])


def update_1d_simple(axis, times, counts_list):
    """

    Args:
        axis: axes object
        times: JG: THIS IS NOT USED! WHAT IS IT? => add comment, e.g. for future purpose or delete!
        counts_list: list of

    Returns:

    """


    if len(np.shape(counts_list)) == 1:
        counts_list = [counts_list]
    if len(axis.lines) != len(counts_list):
        counts_list = np.transpose(counts_list)


    #print('===>> JG 20181128: number of lines', len(axis.lines),' number of list ', len(counts_list), ' (they should be the same!!)')
    # if len(axis.lines) != len(counts_list):
    #     print('UUUUUU axes.lines:', len(axis.lines), 'len counts:', len(counts_list))

    #assert len(axis.lines) == len(counts_list)
    if len(axis.lines) != len(counts_list): # don't update the plot if the number of lines to plot isn't equal to the number of counts lists
        print('Number of lines to plot is not equal to number of counts lists!!!')
        print('===>> ER 20181201: number of lines', len(axis.lines), ' number of list ', len(counts_list),
              ' (they should be the same!!)')
        #print('counts_list: ', counts_list)
        return

    else:
        for index, counts in enumerate(counts_list):
            axis.lines[index].set_ydata(counts)
        axis.relim()
        axis.autoscale_view()

def plot_counts_vs_pos(axis, data, pos):
    """
    plots magnet position vs. counts for aligning the field with fluorescence

    """
    axis.set_xlabel('position (mm)')
    axis.set_ylabel('kCounts/sec')

    axis.plot(pos, data, linewidth=2.0)

def update_counts_vs_pos(axis, data, pos):
    """

    updates the plot of the position vs counts when scanning the magnet

    """

    if data is None:
        return

    axis.lines[0].set_ydata(data)
    axis.lines[0].set_xdata(pos)
    axis.relim()
    axis.autoscale_view()
