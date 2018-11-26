import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import optimize

import os
#import datetime_from_str

from pylabcontrol.core.helper_functions import datetime_from_str
from pylabcontrol.core.scripts import Script

freq_to_mag = 1. / (2 * 2.8e6)
V_to_dist = 60 # convert galvo voltages to distances 1 V is about 60um

## auto reload external files, so that we can edit the external .py file and inmediately see the changes here
#%load_ext autoreload
#%autoreload 2

def get_freqs_and_data(ESR_FOLDER, plot_with_norm, esr_fixed=False):
    '''
    Get the frequencies of the ESR sweep as well as ESR contrast
    '''
    data_esr = []
    counter = 0  # use a counter to figure out the number of ESR points actually taken, to trunctate the real space coordinate to

    if plot_with_norm:
        for f in sorted(glob.glob('{:s}/data_subscripts/*'.format(ESR_FOLDER))):
            data = Script.load_data(f)
            if not esr_fixed:
                full_data = np.divide(data['esr_data'], data['full_laser_data'])
            else:
                full_data = np.divide(data['esr_data'], data['laser_data'])
            data_esr.append(np.mean(full_data, axis=0))
            counter = counter + 1
    else:
        for f in sorted(glob.glob('{:s}/data_subscripts/*'.format(ESR_FOLDER))):
            data = Script.load_data(f)
            data_esr.append(data['data'])
            counter = counter + 1

    f = data['frequency']
    return f, data_esr, counter

def get_pts(PTS_FOLDER, counter, flip = False):
    '''
    Get the points in real space of the scan. V_to_dist is the conversion factor from galvo voltage to a real distance
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
    '''
    f, data_esr, counter = get_freqs_and_data(ESR_FOLDER, plot_with_norm, esr_fixed)

    print('number of ESRs found: ', counter)
    r = get_pts(PTS_FOLDER, counter, flip)

    data_esr_norm = []  # normalize the ESR data to filter out slow laser power drifts
    for d in data_esr:
        data_esr_norm.append(d / np.mean(d))
    print('min and max of the contrast: ', np.min(data_esr_norm), np.max(data_esr_norm))
    plt.pcolormesh(f, r, data_esr_norm, vmin=0.98, vmax=1.01)  # np.max(data_esr_norm))
    plt.axis([f.min(), f.max(), r.min(), r.max()])
    plt.xlim(min(f), max(f))
    plt.xlabel('frequency (Hz)')
    plt.ylabel('position (um)')
    #    plt.ylim([r.min()-1, r.max()+1])
    plt.ylim([r.min(), r.max()])
    # add zlim
    plt.title('ESR around the magnet')
    plt.show()
