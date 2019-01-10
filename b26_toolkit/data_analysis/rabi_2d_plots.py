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
from pylabcontrol.data_processing.signal_processing import power_spectral_density

freq_to_mag = 1. / (2 * 2.8e6)
V_to_dist = 60 # convert galvo voltages to distances 1 V is about 60um


def plot_2d_fft_rabi_vs_powers(RABI_FOLDER, xlog=False, ylog = False, freq_range = None):
    '''

    Args:
        RABI_FOLDER: folder with the rabi data to use

    Returns:
        rabi_data: rabi contrast APD output 2 / APD output 1 versus tau
        taus: tau values in the rabi sweep
        powers: power values used in the 2D sweep, in dBm

    '''

    rabi_data = []
    taus = []
    powers = []
    freqs = []
    for f in sorted(glob.glob('{:s}/data_subscripts/*'.format(RABI_FOLDER))):
        data = Script.load_data(f)
        if 'tau' not in data or data is None: # not a Rabi folder (e.g., find_nv instead)
            continue
        cnts = np.transpose(data['counts'])
        cnts1 = cnts[1]
        cnts0 = cnts[0]
        single_rabi_data = cnts1/cnts0
        taus = data['tau']

        single_rabi_data -= np.mean(single_rabi_data)
        freqs, single_rabi_data = power_spectral_density(single_rabi_data, (taus[1]-taus[0])*1e-9, freq_range=freq_range)

        if ylog:
            single_rabi_data = np.log10(single_rabi_data)
        rabi_data.append(single_rabi_data)
        power = float(f.split('_')[-1])
        powers.append(power)

    if xlog:
        freqs += 1
        freqs = np.log10(freqs)
        print(freqs)

    rabi_data = np.array(rabi_data)[np.argsort(powers)]
    powers = np.sort(np.array(powers))
    plt.pcolormesh(freqs, powers, rabi_data)
    plt.xlabel('frequencies (Hz)')
    plt.ylabel('input power (dBm)')

def plot_2d_rabi_vs_powers(RABI_FOLDER):
    '''

    Args:
        RABI_FOLDER: folder with the rabi data to use

    Returns:
        rabi_data: rabi contrast APD output 2 / APD output 1 versus tau
        taus: tau values in the rabi sweep
        powers: power values used in the 2D sweep, in dBm

    '''

    rabi_data = []
    taus = []
    powers = []
    for f in sorted(glob.glob('{:s}/data_subscripts/*'.format(RABI_FOLDER))):
        data = Script.load_data(f)
        if 'tau' not in data or data is None: # not a Rabi folder (e.g., find_nv instead)
            continue
        cnts = np.transpose(data['counts'])
        cnts1 = cnts[1]
        cnts0 = cnts[0]
        single_rabi_data = cnts1/cnts0
        rabi_data.append(single_rabi_data)
        taus = data['tau']
        power = float(f.split('_')[-1])
        powers.append(power)
    rabi_data = np.array(rabi_data)[np.argsort(powers)]
    powers = np.sort(np.array(powers))
    plt.pcolormesh(taus, powers, rabi_data)
    plt.xlabel('tau values (ns)')
    plt.ylabel('input power (dBm)')

    return rabi_data, taus, powers


