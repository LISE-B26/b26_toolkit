"""
    This file is part of b26_toolkit, a PyLabControl add-on for experiments in Harvard LISE B26.
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


    This file contains script to visualize leviation data acquired with a video camera

"""
import matplotlib.pyplot as plt
import pims
from pims import pipeline
from skimage.color import rgb2gray
rgb2gray_pipeline = pipeline(rgb2gray)
from scipy.ndimage.filters import gaussian_filter

import numpy as np

def plot_frame(file_path, frames, xy_position = None, gaussian_filter_width=None, xylim = None):
    """

    plots frames of the video

    Args:
        file_path: path to video file
        frames: integer or list of integers of the frames to be plotted
        xy_position: xy position of the magnet. If provided, this position will be plotted on top of  the image
        gaussian_filter_width: if not None apply Gaussian filter
        xylim: xylim to zoom in on a region, if not specified (ie None) show full image


    """

    if not hasattr(frames, '__len__'):
        frames = [frames]

    v = pims.Video(file_path)
    video = rgb2gray_pipeline(v)

    if not gaussian_filter_width is None:
        gaussian_filter_pipeline = pipeline(gaussian_filter)
        video = gaussian_filter_pipeline(video, gaussian_filter_width)


    frame_shape = np.shape(video[frames[0]])

    for frame in frames:
        plt.figure()
        plt.imshow(video[frame], cmap='pink')
        if not xy_position is None:
            plt.plot(xy_position[frame, 0], xy_position[frame, 1], 'xg', markersize = 30, linewidth = 4)
            # plot also the positins obtained with trackpy
            if len ( xy_position[frame]) == 4:
                plt.plot(xy_position[frame, 2], xy_position[frame, 3], 'xr', markersize=30, linewidth = 2)

        if xylim is None:
            xlim = [0, frame_shape[0]]
            ylim = [0, frame_shape[1]]
        else:
            xlim, ylim = xylim

        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.show()


# older stuff
def plot_fft(x, frame_rate, start_frame = 0, window_width = 10000, plottype= 'lin'):
    """

    Args:
        x:
        frame_rate: frame rate in fps of the original video
        start frame: first frame to use for the window
        window_width: number of frames to use in window
        plottype: 'lin' linear plot, 'log' loglog plot, 'semilogy' semilog y plot

    Returns:

    """


    # x,y = zip(*max_coors)

    plt.figure()

    if plottype== 'semilogy':
        plt.semilogy(freqs[1:], np.abs(fft)[1:])
    elif plottype== 'lin':
        plt.plot(freqs[1:],np.abs(fft)[1:])
    elif plottype== 'log':
        plt.loglog(freqs[1:],np.abs(fft)[1:])

    plt.title("Fourier Transform")
    plt.xlabel('Frequency(Hz)')
    plt.ylabel('Amplitude (arb)')
    plt.xlim((1,frame_rate/2))
    plt.show()


def plot_ringdown(filepaths, frequency, window_width, fps, bead_diameter, axis='x', starting_frame=0,
                  save_filename=None):
    """
    Takes a file or sequential files with a set of bead positions (such as created by extract_motion), and computes and plots the ringdown
    filepaths: a list of filepaths to .csv files containing the bead position trace in (x,y) coordinates
    frequency: oscillation frequency of mode
    window_width: width of window centered on frequency over which to integrate
    fps: frame rate in frames per second of the data
    bead_diameter: diameter (in um) of the bead
    axis: either 'x' or 'y', specifies which direction to look at oscillations in (if using the reflection on
        the right wall, this is x for z mode, y for xy modes)
    starting_frame: all frames before this one will not be included
    save_filename: if provided, uses this path to save the plot, otherwise saves in the original filepath .png
    """
    assert (axis in ['x', 'y'])

    # bead density in kg/m^3
    bead_density = 7500
    # bead density in kg
    bead_mass = 4 * np.pi / 3 * bead_density * (bead_diameter * 1e-6 / 2) ** 3
    data_ringdown = []
    # read all filepaths and combines them into one dataframe for processing
    for filepath in filepaths:
        data_ringdown.append(pd.read_csv(filepath, usecols=[1, 2], names=['x', 'y'], skiprows=1))
    data_ringdown = pd.concat(data_ringdown)

    fit_params, freq, amp_vel, fig = fit_Q(data_ringdown[axis][starting_frame:] * pix_to_um, 1. / fps, bead_mass,
                                           frequency_range=[frequency - window_width / 2.0,
                                                            frequency + window_width / 2.0], tmax=45,
                                           nbin=1e4)
    plt.title('Ringdown at ' + str(frequency) + 'Hz, Q = ' + str(int(frequency * (fit_params[1] * 60))), fontsize=16)

    if (save_filename == None):
        save_filename = filepaths[0][-4:] + '.png'

    fig.savefig(save_filename)