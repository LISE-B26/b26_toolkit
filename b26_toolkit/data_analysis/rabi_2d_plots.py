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


def get_rabi_and_powers(RABI_FOLDER):
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

    print('size of rabi_data: ', np.size(rabi_data))
    print('size of powers: ', np.size(powers))
    print('size of tau values: ', np.size(taus))
    print('powers: ', powers)
    plt.pcolormesh(taus, powers, rabi_data)
    plt.xlabel('tau values (ns)')
    plt.ylabel('input power (dBm)')
  #  plt.xlim(min(taus), max(taus))

