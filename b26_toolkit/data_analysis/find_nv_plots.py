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

V_to_dist = 60 # convert galvo voltages to distances 1 V is about 60um


def plot_find_nv_pos(path, convert_to_um=True):
    '''


    Args:
        path: the filename of the
        convert_to_um:

    Returns:

    '''
    findnv_paths = glob.glob(os.path.abspath(path + '/data_subscripts/*/data_subscripts/*find_nv*'))
    findnv_paths += glob.glob(os.path.abspath(path + '/data_subscripts/*find_nv*'))
    x = []
    y = []
    datetimes = []
    for p in findnv_paths:
        data = Script.load_data(p)
        #     print(data['maximum_point']['x'])
        x.append(data['maximum_point']['x'])
        y.append(data['maximum_point']['y'])
        datetimes.append(Script.load_time(os.path.basename(p)))
    t0 = datetimes[0]
    t = [date - t0 for date in datetimes]
    times = [temp_time.seconds/(60*60) for temp_time in t]# convert times to hours
    x = np.array(x)-np.mean(np.array(x))
    y = np.array(y)-np.mean(np.array(y))

    if convert_to_um:
        x = x*V_to_dist
        y = y*V_to_dist

    labeltime = 'time (hrs)'
    labely = 'y position '
    labelx = 'x position'
    if convert_to_um:
        labely += '(um)'
        labelx += '(um)'
    else:
        labely += '(V)'
        labelx += '(V)'

    fig, ax = plt.subplots(1, 3, figsize = (20, 4))
    ax[0].plot(times, x, marker='o', linestyle='')
    ax[0].set_xlabel(labeltime)
    ax[0].set_ylabel(labelx)
    ax[1].plot(times, y, marker='o', linestyle='')
    ax[1].set_xlabel(labeltime)
    ax[1].set_ylabel(labely)
    ax[2].plot(x, y, marker='o', linestyle='')
    ax[2].set_xlabel(labelx)
    ax[2].set_ylabel(labely)

    return x, y