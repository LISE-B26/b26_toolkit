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


    This file contains script to process leviation data acquired with a video camera

"""


import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import scipy.optimize as sopt
import glob

import pims
from pims import pipeline
from skimage.color import rgb2gray
rgb2gray_pipeline = pipeline(rgb2gray)
from scipy.ndimage.filters import gaussian_filter

import trackpy as tp
# status bar
from ipywidgets import FloatProgress
from IPython.display import display
import datetime

# following functions are imported just to make them available within the context of b26toolkit - analysis of levitation data
# from PyLabControl.src.data_processing.signal_processing import power_spectral_density

## TMP JG 20180522 - just because the above import doesn't work

def power_spectral_density(x, time_step, frequency_range = None):
    """
    returns the *single sided* power spectral density of the time trace x which is sampled at intervals time_step
    i.e. integral over the the PSD from f_min~0 to f_max is equal to variance over x
    Args:
        x (array):  timetrace
        time_step (float): sampling interval of x
        freq_range (array or tuple): frequency range in the form [f_min, f_max] to return only the spectrum within this range

    Returns:

    """
    N = len(x)
    p = 2 * np.abs(np.fft.rfft(x)) ** 2 / N * time_step

    f = np.fft.rfftfreq(N, time_step)

    if not frequency_range is None:
        assert len(frequency_range) == 2
        assert frequency_range[1] > frequency_range[0]

        bRange = np.all([(f > frequency_range[0]), (f < frequency_range[1])], axis=0)
        f = f[bRange]
        p = p[bRange]

    return f, p

## end TMP JG 20180522


## stuff by Jan including modified functions from Aaron

def get_video_info(filename):
    """
    Args:
        filename: path to video file

    Returns: a dictionary with the video metadata

    """
    v = pims.Video(filename)

    video_info = {
        'frame_rate': v.frame_rate,
        'duration': v.duration
    }

    return video_info

def reencode_video(filepath, filepath_target = None):
    """

    uses ffmpeg to reencode the video and removes the first second. This is usually
    necessary for the 1200fps videos, which are timestamped weirdly and often have the first second
    corrupted. If you turn this off and receive the error "AVError: [Errno 1094995529] Invalid data found
    when processing input", try turning this on. Will double the runtime of the function.

    Args:
        filepath: path to file of original video
        filepath_target: target file (optional) if None same as input with replacing ".avi" by  "_reencode.avi"

    Returns: nothing but writes reencoded file to disk

    """

    if filepath_target is None:
        filepath_target = filepath.replace('.avi', '_reencode.avi')

    print('start time:\t{:s}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    # command string: ffmpeg -i Z:\...\ringdown.avi -s 1 -c copy Z:\...\ringdown_reencode.avi
    # calls ffmpeg, -i specifies input path, -ss 1 cuts first second of video, -c copy copies
    # the input codec and uses it for the output codec, and the last argument is the output file
    # cutting the first second isn't always necessary, but sometimes the videos will not load without it
    cmd_string = "ffmpeg -i " + filepath + " -ss 1 -c copy " + filepath_target
    # performs system (command line) call
    os.system(cmd_string)
    print('end time:\t{:s}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    print('wrote:\n{:s}'.format(filepath_target))

def extract_motion(filepath, target_path=None,  gaussian_filter_width=2, use_trackpy =False, show_progress = True,
                   trackpy_parameters = None, max_frames = None):
    """
    Takes in a bead video, and saves the bead's extracted position in every frame in a csv, and
    (optionally) an uncorrupted version of the video
    filepath: path to an avi file of the bead
    target_path: (optional) if None same as video path

    gaussian_filter_width: width (in pixels) of the gaussian used to smooth the video

    use_trackpy: if True use Trackpy to better localize the brightest point, need to provid trackpy_parameters

    show_progress: if True show progress (use when called from ipython notebook)

    trackpy_parameters: parameters for trackpy (only needed if use_trackpy is True)

    max_frames: max number of frames to analyze, usefull for debugging, if none do entire video

    returns: path to .csv file
    """



    if target_path is None:
        target_path = os.path.dirname(filepath)

    assert os.path.isdir(target_path)

    if use_trackpy:
        assert 'diameter' in trackpy_parameters
        if not 'minmass' in trackpy_parameters:
            trackpy_parameters['minmass'] = None
        trackpy_parameters['missing_frames'] = []

    target_filename = os.path.join(target_path, os.path.basename(filepath).replace('.avi',  '_data_globalmax.csv'))

    v = pims.Video(filepath)
    processed_v = rgb2gray_pipeline(v)

    gaussian_filter_pipeline = pipeline(gaussian_filter)
    filtered = gaussian_filter_pipeline(processed_v, gaussian_filter_width)

    if max_frames is None:
        max_frames=len(filtered)

    if show_progress:
        f = FloatProgress(min=0, max=max_frames)  # instantiate the bar
        display(f)  # display the bar
    start_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print('start time:\t{:s}'.format(start_time))

    image_size = np.shape(filtered[0]) # size of image
    # find the maximum brightness pixel for every frame, corresponding to some consistent point at the bead
    max_coors = []
    for index, image in enumerate(filtered):
        pixel_max = np.argmax(image)

        po = np.unravel_index(pixel_max, image_size)
        po = [po[1], po[0]] # flip the order to get x, y

        if use_trackpy:
            locate_info = tp.locate(image, trackpy_parameters['diameter'], minmass=trackpy_parameters['minmass'])
            pts = locate_info[['x', 'y']].as_matrix()

            # pick the one that is closest to the original one
            pts = pts[np.argmin(np.array([np.linalg.norm(p - np.array(po)) for p in pts]))]
            if len(pts)< 2:
                trackpy_parameters['missing_frames'].append(index)

            po = po + [pts[0], pts[1]] # append point found by findnv to dataset of current frame
        max_coors.append(po)
        # max_coors.append((pixel_max % image_size[0], pixel_max / image_size[1]))
        if index % 1000 == 0:
            if show_progress:
                f.value += 1000


        # tmp for debugging purpuse end early
        if not max_frames is None:
            if index>max_frames:
                break

    df = pd.DataFrame(max_coors)
    df.to_csv(target_filename)
    end_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print('end time:\t{:s}'.format(end_time))
    print('file written to:\n{:s}'.format(target_filename))

    info = {
        'filename_xy_position': target_filename,
        'image_size':image_size,
        'gaussian_filter_width':gaussian_filter_width,
        'N_frames':len(filtered),
        'use_trackpy':use_trackpy,
        'start_time':start_time,
        'end_time':end_time,
        'max_frames':max_frames
    }

    if use_trackpy:
        info.update({
            'trackpy_parameters': trackpy_parameters
        })


    return info


def load_xy_time_trace(filepath, center_at_zero = True):
    """
    Takes in a filepath to a csv containing the bead positions

    Args:
        filepath: filepath to a csv containing the bead positions
        center_at_zero: is true substract mean, else return raw data

    Returns: the x and y position as a Nx2 matrix where N is the number of datapoints and the first column is x and the second y

    """
    #load data
    xy = pd.read_csv(filepath)
    xy = xy.as_matrix()

    xy = xy[:,1:] # drop first column since this is just the index

    if center_at_zero:
        xy = xy - np.mean(xy, axis=0)

    return xy

    """
    Takes in a filepath to a csv containing the bead positions and a frame rate, and plots the fourier transform
    for a given frame window
    filepath: path to a csv containing the bead (x,y) position for each frame, such as the output from
        extract_motion()
    frame_rate: frame rate in fps of the original video
    start frame: first frame to use for the window
    window_width: number of frames to use in window
    """

## stuff from Aaron
def load_data(path):
    """
    loads the data that has been save with Script.save.
    Args:
        path: path to folder saved by Script.save or raw_data folder within
    Returns:
        a dictionary with the data of form
        data = {param_1_name: param_1_data, ...}
    """


    # check that path exists
    if not os.path.exists(path):
        print(path)
        raise AttributeError('Path given does not exist!')

    # windows can't deal with long filenames (>260 chars) so we have to use the prefix '\\\\?\\'
#         if len(path.split('\\\\?\\')) == 1:
#             path = '\\\\?\\' + os.path.abspath(path)


    # if raw_data folder exists, get a list of directories from within it; otherwise, get names of all .csv files in
    # current directory
    data = {}
    # if self.RAW_DATA_DIR in os.listdir(path): #8/26/16 AK: self not defined in static context
    #     data_files = os.listdir(os.path.join(path, self.RAW_DATA_DIR + '/'))
    #     path = os.path.join(path, self.RAW_DATA_DIR + '/')
    #
    # else:
    if 'raw_data' in os.listdir(path):  #temporarily hardcoded
        data_files = os.listdir(os.path.join(path, 'raw_data' + '/'))
        path = os.path.join(path, 'raw_data' + '/')

    else:
        data_files = glob.glob(os.path.join(path, '*.csv'))

    # If no data files were found, raise error
    if not data_files:
        raise AttributeError('Could not find data files in {:s}'.format(path))

    # import data from each csv
    for data_file in data_files:
        # get data name, read the data from the csv, and save it to dictionary
        data_name = data_file.split('-')[-1][0:-4] # JG: why do we strip of the date?
        imported_data_df = pd.read_csv(os.path.join(path, data_file))

        # check if there are real headers, if the headers are digits than we ignore them because then they are just indecies
        # real headers are strings (however, the digits are also of type str! that why we use the isdigit method)
        column_headers = list(imported_data_df.columns.values)
        if sum([int(x.isdigit()) for x in column_headers]) != len(column_headers):
            data[data_name] = {h: imported_data_df[h].as_matrix() for h in column_headers}
        else:
            # note, np.squeeze removes extraneous length-1 dimensions from the returned 'matrix' from the dataframe
            data[data_name] = np.squeeze(imported_data_df.as_matrix())

    return data

def exponential(x, amp, tau, offset):
    """
    Defines an exponential with argument x, amplitude amp, decay time tau, and offset offset
    """
    return amp * np.exp(-x/tau) + offset


def get_freq_and_amp_timetrace(x, frequency_range, dt, nbin=10e3, velocity_type=True):
    """
    returns the time max freq and totel power in freq range for chunks of timetrace

    x: timetrace

    frequency_range: freq range

    nbin: length of timetrace chunk

    velocity_type:  False - returns totel power of position spectral density in freq range for chunks of timetrace
                    True - returns totel power of velocity spectral density in freq range for chunks of timetrace


    """

    x = x[0:int(np.floor(len(x) / nbin) * nbin)]

    x = x.reshape([int(np.floor(len(x) / nbin)), int(nbin)])

    df = 1. / (nbin * dt)

    f_of_t = []
    a_of_t = []
    for y in x:
        f, p = power_spectral_density(y, dt, frequency_range)
        f_of_t.append(f[np.argmax(p)])

        if velocity_type:
            a_of_t.append(np.sum(p * (2 * np.pi * f) ** 2) * df)
        else:
            a_of_t.append(np.sum(p) * df)

    return f_of_t, a_of_t


def fit_Q(x, dt, bead_mass, frequency_range=[60, 80], tmax=10, nbin=1e4, plot=True):
    """
    x: timetrace
    dt: time interval of datapoints in x
    frequency_range: range of freq over which we get the energy and look for mode freq.
    tmax: max time for fit in min
    nbin: number of points of x over which we extract the freq and energy
    """

    def exponential(x, amp, tau, c):
        return amp * np.exp(-x / tau) + 0

    x -= np.mean(x)  # remove offset

    t = dt * np.arange(len(x)) / 60  # time in min

    freq, amp_vel = get_freq_and_amp_timetrace(x, frequency_range=frequency_range, dt=dt, nbin=nbin, velocity_type=True)

    tbin = np.arange(len(freq)) * nbin * dt / 60  # time of binned data in min

    freq = freq[5:]
    amp_vel = amp_vel[5:]
    tbin = tbin[5:]

    energy_vel = np.array(amp_vel) * bead_mass * 1e-12  # factor 1e-12 because x is in um

    energy_vel /= 1.38e-23  # convert to Kelvins

    energy_min_res = (168) ** 2 * bead_mass * 1e-12
    energy_min_res /= 1.38e-23

    imax = np.sum((tbin < tmax) * 1)  # index corresponding to tmax

    Eo = np.mean(energy_vel[0])  # initial energy
    Ef = np.mean(energy_vel[-1])  # final energy

    #     tau0 = tbin[np.sum((energy_vel>Eo/3)*1)] # time after which initial energy is only 1/3 ~ 1/e

    #     po = [Eo,tau0, Ef]

    po = [Eo, 10, 1e9]

    fit_params, _ = sopt.curve_fit(exponential, tbin[:imax], energy_vel[:imax], p0=po)

    fig = plt.figure()
    ax = plt.subplot(111)
    plt.plot(tbin, energy_vel)
    plt.plot(tbin, exponential(tbin, fit_params[0], fit_params[1], fit_params[2]))

    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='minor', labelsize=12)

    plt.xlabel('time (min)', fontsize=14)
    plt.ylabel('energy (K)', fontsize=14)

    ax.yaxis.offsetText.set_fontsize(12)

    return fit_params, freq, amp_vel, fig






def compute_gamma(filepaths, cd_height, bead_diameter, frequency, window_width, fps, axis='x', starting_frame=0,
                  time_bin_size=1e4):
    """
    Computes and plots the slope of the anharmonicity parameter gamma, which describes the relationship between
    resonator amplitude squared and frequency
    filename: filepath to csv containing bead position data pairs (x,y)
    cd height: height (in um) from superconductor to top of bead during cooldown
    bead_diameter: diameter (in um) of the bead
    frequency: oscillation frequency of mode
    window_width: width of window centered on frequency in which to look for the frequency of the mode
    fps: frame rate in frames per second of the data
    axis: either 'x' or 'y', specifies which direction to look at oscillations in (if using the reflection on
        the right wall, this is x for z mode, y for xy modes)
    starting_frame: all frames before this one will not be included. Allows exclusion of any non-linear part
        of the ringdown
    time_bin_size: number of frames to include in each time bin. Decreasing this increases time resolution, but
        if it's too small there are insufficient points for the fourier transform and it looks like steps
    """
    matplotlib.rcParams.update({'font.size': 15})

    z0 = (cd_height - bead_diameter / 2.0)

    data_ringdown = []
    # read all filepaths and combines them into one dataframe for processing
    for filepath in filepaths:
        data_ringdown.append(pd.read_csv(filepath, usecols=[1, 2], names=['x', 'y'], skiprows=1))
    data_ringdown = pd.concat(data_ringdown)

    freqs_1, amp_1 = get_freq_and_amp_timetrace(data_ringdown[axis][starting_frame:] * pix_to_um,
                                                [frequency - window_width / 2.0, frequency + window_width / 2.0],
                                                1. / fps, nbin=time_bin_size, velocity_type=False)
    freqs_1 = np.array(freqs_1[0:])
    amp_1 = np.array(amp_1[0:]) / z0 ** 2

    def lin(x, m, b):
        return m * x + b

    # first fit to find what the zero_amplitude frequency is in order to properly normalize frequencies
    fit_params, covar = sopt.curve_fit(lin, np.array(amp_1) / z0 ** 2, freqs_1)
    w0 = fit_params[1]

    freqs_1 = freqs_1 / w0

    plt.figure()
    plt.plot(amp_1, freqs_1, '.')
    plt.xlabel('Amplitude Squared (Normalized to cooldown height)')
    plt.ylabel('Frequency (Normalized to zero amplitude)')

    # then fit the normalized
    fit_params, pcov = sopt.curve_fit(lin, amp_1, freqs_1)
    plt.plot(amp_1, lin(amp_1, *fit_params))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    err = np.abs(pcov[0][0]) ** 0.5
    # print('Slope: ', fit_params[0])
    plt.title('Amplitude Squared vs Frequency, gamma = ' + str(round(fit_params[0], 2)))
    ax = plt.gca()
    tx = ax.xaxis.get_offset_text()
    tx.set_fontsize(10)
    plt.savefig(filepath[0][:-19] + '_anharmonicity2.png')