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

# ==========================================================================================================
# ==========================================================================================================
#
#
# ==========================================================================================================
# ==========================================================================================================

import matplotlib.pyplot as plt

from b26_toolkit.plotting.plots_2d import plot_fluorescence
from pylabcontrol.core import Script


def galvo_images(data_path, target_path = None):
    """
    save load data from galvo scans and save images to target directory
    Args:
        data_path: path to image data
        target_path: target path to save images

    Returns:

    """
    if target_path == None:
        target_path = DATA_PATH

    data  = Script.load_data(DATA_PATH)
    number_of_images = len([k for k in list(data.keys()) if len(k.split('image'))>1])

    for c in range(number_of_images):
        k  = 'image_{:d}'.format(c)
        fig = plt.figure()
        ax = plt.subplot(111)
        plot_fluorescence(data[k], extent =[0.02, 0.17, 0.05, -0.10], axes = ax)
        fig.savefig('{:s}/{:s}.png'.format(TARGET_PATH, k))
        fig.close()

if __name__ == '__main__':

    DATA_PATH = 'Z:\\Lab\\Cantilever\\Measurements\\20160524_Focsing\\160524-15_18_51_reflection\\'

    galvo_images(DATA_PATH)
