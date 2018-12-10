import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import optimize

import os
#import datetime_from_str

from pylabcontrol.core.helper_functions import datetime_from_str
from pylabcontrol.core.script import Script

freq_to_mag = 1. / (2 * 2.8e6)
V_to_dist = 60 # convert galvo voltages to distances 1 V is about 60um

## auto reload external files, so that we can edit the external .py file and inmediately see the changes here
#%load_ext autoreload
#%autoreload 2

def get_freqs_and_data(ESR_FOLDER, plot_with_norm, esr_fixed=False, get_times=True):
    '''
    Get the frequencies of the ESR sweep as well as ESR contrast.

    ESR_FOLDER: folder with ESR data in it

    plot_with_norm: normalize to the photodiode signal (i.e. reflectted or incident laser power)

    esr_fixed: corresponds to data taken after we changed the name of the laser data in self.data, for when plot_with_norm is on.

    get_times: record and return the time duration of each ESR experiment.

    '''
    data_esr = []
    times = []
    counter = 0  # use a counter to figure out the number of ESR points actually taken, to trunctate the real space coordinate to

    if plot_with_norm:
        for f in sorted(glob.glob('{:s}/data_subscripts/*'.format(ESR_FOLDER))[-1]):
            data = Script.load_data(f)
            if 'esr_data' not in data or data is None: # not an ESR folder (e.g., find_nv instead)
                continue
            if not esr_fixed:
                full_data = np.divide(data['esr_data'], data['full_laser_data'])
            else:
                full_data = np.divide(data['esr_data'], data['laser_data'])
            data_esr.append(np.mean(full_data, axis=0))
            counter = counter + 1

            if get_times:
                split_path = f.split('-')
                time = split_path.pop(-1)
                day = split_path.pop(-1)
                day = day.split('18')[-1] # will need to change 18 to 19 when it becomes 2019
                time = time.split('_')
                time = float(time[0]) + float(time[1])/60 + float(time[2])/3600 + float(day)*24
                times.append(time)
    else:
        for f in sorted(glob.glob('{:s}/data_subscripts/*'.format(ESR_FOLDER))):
            data = Script.load_data(f)
            if 'esr_data' not in data or data is None: # not an ESR folder (e.g., find_nv instead)
                continue
            data = Script.load_data(f)
            data_esr.append(data['data'])
            counter = counter + 1

            if get_times:
                split_path = f.split('-')
                time = split_path.pop(-1)
                day = split_path.pop(-1)
                day = day.split('18')[-1] # will need to change 18 to 19 when it becomes 2019
                time = time.split('_')
                time = float(time[0]) + float(time[1])/60 + float(time[2])/3600 + float(day)*24
                times.append(time)

    f = data['frequency']
    if not get_times:
        return f, data_esr, counter
    else:
        return f, data_esr, counter, times

def get_pts(PTS_FOLDER, counter, flip = False):
    '''
    Get the points in real space of the scan. V_to_dist is the conversion factor from galvo voltage to a real distance

    flip: allows you to reverse the direction of the plot(e.g. towards or away from the iron/bead)

    '''
    data = Script.load_data(PTS_FOLDER)
    r = []
    for pt in data['nv_locations']:
        r.append( np.sqrt((pt[0]-data['nv_locations'][0][0])**2+(pt[1]-data['nv_locations'][0][1])**2))

    r = V_to_dist * np.array(r)
    r = r[0:counter] # truncate # of points to the actual number of data points taken
    if flip:
        r = np.flipud(r)
    return r


def plot_2d_esr(ESR_FOLDER, PTS_FOLDER, plot_with_norm=True, esr_fixed=False, flip=False):
    '''
    Retrieve the ESR contrast, ESR frequency points, and real space points. Then plot it in a 2D plot.

    plot_with_norm: whether or not to divide the ESR data by the laser / photodiode power drift data
    esr_fixed: whether or not the experiment is with the new ESR script, i.e. setting the center frequency for each point.
        This is relevant because we changed the name of some of the data with the new ESR script.
    flip: whether or not the points scan was taken starting far away or close to the bead / iron
    '''
    f, data_esr, counter = get_freqs_and_data(ESR_FOLDER, plot_with_norm, esr_fixed)

    print('number of ESRs found: ', counter)
    r = get_pts(PTS_FOLDER, counter, flip)

    data_esr_norm = []  # normalize the ESR data to filter out slow laser power drifts
    for d in data_esr:
        data_esr_norm.append(d / np.mean(d))

    min_contrast =  np.min(data_esr_norm)
    max_contrast = np.max(data_esr_norm)
    print('min and max of the contrast: ', min_contrast, max_contrast)
    plt.pcolormesh(f, r, data_esr_norm, vmin=min_contrast, vmax=max_contrast)
    plt.axis([f.min(), f.max(), r.min(), r.max()])
    plt.xlim(min(f), max(f))
    plt.xlabel('frequency (Hz)')
    plt.ylabel('position (um)')
    #    plt.ylim([r.min()-1, r.max()+1])
    plt.ylim([r.min(), r.max()])
    # add zlim
    plt.title('ESR around the magnet')
    plt.show()

def plot_2d_esr_vs_time(ESR_FOLDER, plot_with_norm=False, esr_fixed=False, get_times=True):
    '''
    Retrieve the ESR contrast, ESR frequency arrays, as a function of time. Then plot it in a 2D plot.

    plot_with_norm: whether or not to divide the ESR data by the laser / photodiode power drift data
    esr_fixed: whether or not the experiment is with the new ESR script, i.e. setting the center frequency for each point.
        This is relevant because we changed the name of some of the data with the new ESR script.
    get_times: extract time duration of each ESR experiment to plot on the x axis.
    '''


    if not get_times:
        f, data_esr, counter = get_freqs_and_data(ESR_FOLDER, plot_with_norm, esr_fixed, get_times=False)
        times = np.linspace(0, counter-1, counter)
    else:
        f, data_esr, counter, times = get_freqs_and_data(ESR_FOLDER, plot_with_norm, esr_fixed, get_times=True)
        times = np.array(times)
        times = times - times[0]
        times = times[1:]

    print('number of ESRs found: ', counter)
    data_esr_norm = []
    counter = 0
    for d in data_esr:
        if counter == 0:
            counter += 1
            continue
        data_esr_norm.append(d/np.mean(d))
        counter += 1

    min_contrast = np.min(data_esr_norm)
    max_contrast = np.max(data_esr_norm)
    print('size of data_esr_norm: ', np.size(data_esr_norm))
    print('size of times: ', np.size(times))
    print('times: ', times)
    fig = plt.pcolormesh(f, times, data_esr_norm, vmin=1-(1-min_contrast)/2, vmax=(max_contrast-1)/2+1)  # np.max(data_esr_norm))
    plt.axis([f.min(), f.max(), times.min(), times.max()])
    plt.xlim(min(f), max(f))
    plt.xlabel('frequency (Hz)')
    if not get_times:
        plt.ylabel('time (AU)')
    else:
        plt.ylabel('time (hrs)')
    plt.ylim([times.min(), times.max()])
    plt.title('ESR contrast versus time')
    plt.show()

    return fig

def plot_esr_fit_vs_time(ESR_FOLDER):

    fit_freqs = []
    times = []

    for f in sorted(glob.glob('{:s}/data_subscripts/*'.format(ESR_FOLDER))):
        data = Script.load_data(f)
        if 'esr_data' not in data or data is None: # not an ESR folder (e.g., find_nv instead)
            continue
        fit_params = data['fit_params']
        if fit_params is not None and len(fit_params) and fit_params[0] != -1:  # check if fit valid
            if len(fit_params) == 4:
                # single peak
                fit_freqs.append(fit_params[2])
                split_path = f.split('-')
                time = split_path.pop(-1)
                day = split_path.pop(-1)
                day = day.split('18')[-1]  # will need to change 18 to 19 when it becomes 2019
                time = time.split('_')
                time = float(time[0]) + float(time[1]) / 60 + float(time[2]) / 3600 + float(day) * 24
                times.append(time)
            elif len(fit_params) == 6:
                # double peak, don't update the frequency - the fit may be bad
                continue

    times = np.array(times)
    times = times - times[0]
    times = times[1:]
    fit_freqs = np.array(fit_freqs)/1.e6

    plt.plot(times, fit_freqs[1:])
    plt.xlabel('time (hours)')
    plt.ylabel('fit frequency (MHz)')