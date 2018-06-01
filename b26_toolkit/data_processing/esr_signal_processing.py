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

from copy import deepcopy

import numpy as np
import scipy.signal as signal
from peakutils.peak import indexes # package required https://pypi.python.org/pypi/PeakUtils

from b26_toolkit.data_processing.fit_functions import fit_lorentzian, get_lorentzian_fit_starting_values, fit_double_lorentzian, lorentzian, double_lorentzian


def fit_esr_old(freq, ampl):
    """
    Returns lorentzian fit parameters for a typical NV esr sweep, giving 4 or 6 parameters depending on if 1 or 2
    lorentzian dips are detected.

    Args:
        freq: 1d array of frequencies which were scanned for esr resonance
        ampl: 1d array of amplitudes corresponding to the frequencies (freq)

    Returns:
        fit parameters for either a single or double lorentizian peak fit to the data, as figured out by find_nv_peaks.
        if a single lorentzian was fit, the list [constant_offset, amplitude, center, fwhm] is returned;
        otherwise, the list [constant_offset, amplitude1, amplitude2, center1, center2, fwhm] is returned

    """
    freq_peaks, ampl_peaks = find_nv_peaks(freq, ampl)

    # figure out if one or two peaks
    if freq_peaks[0] == freq_peaks[1]:
        start_vals = get_lorentzian_fit_starting_values(freq, ampl)
    elif freq_peaks[0] == 0:
        print('data too noisy!!')
        # later we can extend this to fit a Lorenzian to each peak, for now we assume it's noisy data..
        fit = None
    else:
        center_freq = np.mean(freq_peaks)
        start_vals = []
        start_vals.append(get_lorentzian_fit_starting_values(freq[freq < center_freq], ampl[freq < center_freq]))
        start_vals.append(get_lorentzian_fit_starting_values(freq[freq > center_freq], ampl[freq > center_freq]))
        # we already have a good estimate for the peak position, so update that
        start_vals[0][2], start_vals[1][2] = freq_peaks
        start_vals = [
            np.mean([start_vals[0][0], start_vals[1][0]]),  # offset
            np.sum([start_vals[0][3], start_vals[1][3]]),  # FWHM
            start_vals[0][1], start_vals[1][1],  # amplitudes
            start_vals[0][2], start_vals[1][2]  # centers
        ]

    try:
        if len(freq_peaks) == 2:
            fit = fit_double_lorentzian(freq, ampl, starting_params=start_vals)
        elif len(freq_peaks) == 1:
            fit = fit_lorentzian(freq, ampl, starting_params=start_vals)
    except:
        # ESR fit failed!
        fit = None

    return fit

def find_nv_peaks(freq, data, width_Hz=0.005e9, initial_threshold = 0.00, steps_size = 0.05, ax=None):
    """
    finds the single peak or double peak frequencies of esr spectrum
    Args:
        freq: frequnency points of esr spectrum
        data: data points of esr spectrum
        width_Hz: expected width of peak
        ax: optional axes object to plot processed data

    Returns:
        freq_max: peak frequency(ies)
        data_max: esr signal at peak frequency(ies)

    """


    freq_to_mag = 1. / (2 * 2.8e6)  # conversion factor from frequency to magentic field in Gauss (1/(2*gyromag ratio))

    def find_peaks_pts(data, width, initial_threshold, steps_size):
        """
        find the maximum points in data
        Args:
            data: processed data that has one or two Gaussian like features
            width: expected width of features

        Returns:
            idx, data[idx]
            index and value of peaks

        """

        threshold = initial_threshold  # initial threshold
        continue_search = True
        # number of peaks of previous attempt
        number_peaks_previous = -1


        while continue_search:
            idx = indexes(np.array(data), thres=threshold, min_dist=width)
            #             print(idx)
            if len(idx) > 2:
                threshold += steps_size
            elif len(idx) == 2:
                # double peak detected, maybe need to check if reasonable here
                continue_search = False
            elif len(idx) == 1:
                # single peak detected
                continue_search = False
            elif len(idx) == 0:
                # peak detection in this iteration but was successful before
                # this means that we should go back and raise the threshold slower
                if number_peaks_previous >0:
                    threshold -= steps_size # go back to previous threshold that was successful
                    steps_size /= 2. # reduce stepsize by factor 2
                else:
                    # search failed
                    continue_search = False

            number_peaks_previous = len(idx)

        return idx, data[idx]

    def check_double_peak(freq_max, peak_max, fo = 2.878):
        """
        checks if double peak is physical of not return the frequency and value of the physical peak
        Args:
            freq_max: vector of length 2 with the frequencies of the two peaks
            peak_max: value of two peaks
            fo: NV center frequency without Zeeman shift (2.878 GHz)

        Returns:
            freq_max, peak_max
            vectors of length one or two that contain the frequencies and values of the peak(s)

        """
        assert len(freq_max) == 2



        # calculate symmetry with respect to expected center freq
        df = np.abs(freq_max - fo)
        asymmetry_f = max(df) / min(df) - 1
        # calculate symmetry with respect to peak height
        asymmetry_p = np.abs(np.diff(peak_max) / np.mean(peak_max))[0]
        #
        # print('freq_o', np.mean(freq_max),
        #       'asymmetry f', asymmetry_f,
        #       'asymmetry p', asymmetry_p)

        if asymmetry_p > 1:
            freq_max = [freq_max[np.argmax(peak_max)]]

        return freq_max

    #         check symmetry with respect to unshifted NV center symmetry


    # get width in pts
    df = np.mean(np.diff(freq))
    width_pts = int(width_Hz / df)

    #
    sig = deepcopy(data)
    sig /= np.mean(sig)
    sig -= 1
    sig *= -1.0

    # smooth signal with filter
    # int(width_pts/2)*2*3+1: to get a window size about three times the width and odd number
    win = signal.gaussian(int(width_pts / 2) * 2 * 3 + 1, std=width_pts)

    sig_filtered = signal.convolve(sig, win, mode='same') / sum(win)
    max_idx, max_pts = find_peaks_pts(sig_filtered, width_pts, initial_threshold, steps_size)

    if len(max_idx) == 2:
        fo = 2.878 # NV center frequency without Zeeman shift (2.878 GHz)
        if min(freq)> fo or max(freq)<fo:

            freq_max = freq[max_idx]
        else:
            freq_max = check_double_peak(freq[max_idx], max_pts, fo)
        if len(freq_max) == 1:
            data_max = [dx for (fx, dx) in zip(freq, data) if fx == freq_max]
            max_pts = [dx for (fx, dx) in zip(freq, sig_filtered) if fx == freq_max]
        else:
            data_max = data[max_idx]
    elif len(max_idx) == 1:
        freq_max = freq[max_idx]
        data_max = data[max_idx]
    else:
        freq_max = [0, 0]
        data_max = [0, 0]
        max_pts = [0, 0]
        print('No peak found!!!')

    # for a single peak we still return two (identical) values, this makes further processing easier
    if len(freq_max) == 1:
        freq_max = [freq_max[0], freq_max[0]]
        max_pts = [max_pts[0], max_pts[0]]
        data_max = [data_max[0],data_max[0]]

    if ax is not None:
        ax.plot(freq, sig)
        ax.plot(freq, sig_filtered)
        if freq_max[1] != 0:
            ax.plot(freq_max, max_pts, 'o')

        if freq_max[0] == freq_max[1] and freq_max[1] != 0:
            ax.set_title('single: {:e}'.format(freq_max[0]))
        else:
            ax.set_title('mag field: {:e} Gauss'.format(freq_to_mag * np.abs(freq_max[0] - freq_max[1])))

    return freq_max, data_max


def fit_esr(freq, ampl, min_counts = .5, contrast_factor = 1.5, strain_filtering=False, verbose = False):
    """
    Returns lorentzian fit parameters for a typical NV esr sweep, giving 4 or 6 parameters depending on if 1 or 2
    lorentzian dips are detected.
    Args:
        freq: 1d array of frequencies which were scanned for esr resonance
        ampl: 1d array of amplitudes corresponding to the frequencies (freq)
        min_counts: minimum counts for an ESR to not be considered noise
        contrast_factor: require that peaks are at least a factor 3 deeper than the noise (calculated as std dev of data-fit )
    Returns:
        fit parameters for either a single or double lorentizian peak fit to the data, as figured out by find_nv_peaks.
        if a single lorentzian was fit, the list [constant_offset, amplitude, center, fwhm] is returned;
        otherwise, the list [constant_offset, amplitude1, amplitude2, center1, center2, fwhm] is returned
    """
    if strain_filtering:
        MAX_STRAIN = 5e7
    else:
        MAX_STRAIN = np.inf
    F0 = 2.878e9
    MIN_WIDTH = 3*np.mean(np.diff(freq)) # set the minumum width to at least 3 times the sample spacing
    MAX_WIDTH = 100e6  # set the max width 100MHz
    freq_peaks, ampl_peaks = find_nv_peaks(freq, ampl)


    if verbose:
        print(('found peaks at ', freq_peaks))

    if verbose:
        print(('minimum peak width:',  MIN_WIDTH))

    # check if scanning full range for two peaks or half range for one peak
    if max(freq) < F0:
        start_vals = get_lorentzian_fit_starting_values(freq, ampl)
        start_vals[2] = freq_peaks[0]
        try:

            if verbose:
                print(('fit single peak with initial values', start_vals))

            fit = fit_lorentzian(freq, ampl, starting_params=start_vals,
                                 bounds=[(0, -np.inf, 0, 0), (np.inf, 0, np.inf, np.inf)])
        except:
            # ESR fit failed!
            fit = None
            return fit
    else:
        # check for double peaks on one side
        if (len(freq_peaks) == 2) and not (freq_peaks[0] == freq_peaks[1]) and (
                    np.abs((np.mean([freq_peaks[0], freq_peaks[1]]) - F0)) > MAX_STRAIN):
            if np.abs(ampl_peaks[0]) > np.abs(ampl_peaks[1]):
                if np.abs(freq_peaks[0] < F0):
                    freq_peaks[1] = F0 + np.abs(freq_peaks[0] - F0)
                else:
                    freq_peaks[1] = F0 - np.abs(freq_peaks[0] - F0)
            else:
                if np.abs(freq_peaks[1] < F0):
                    freq_peaks[0] = F0 + np.abs(freq_peaks[1] - F0)
                else:
                    freq_peaks[0] = F0 - np.abs(freq_peaks[1] - F0)
        # figure out if one or two peaks
        if freq_peaks[0] == freq_peaks[1]:
            start_vals = get_lorentzian_fit_starting_values(freq, ampl)
            start_vals[2] = freq_peaks[0]
            freq_peaks = [freq_peaks[0]]
        elif freq_peaks[0] == 0:
            print('data too noisy!!')
            # later we can extend this to fit a Lorenzian to each peak, for now we assume it's noisy data..
            fit = None
        else:
            center_freq = np.mean(freq_peaks)
            start_vals = []
            start_vals.append(get_lorentzian_fit_starting_values(freq[freq < center_freq], ampl[freq < center_freq]))
            start_vals.append(get_lorentzian_fit_starting_values(freq[freq > center_freq], ampl[freq > center_freq]))
            # we already have a good estimate for the peak position, so update that
            start_vals[0][2], start_vals[1][2] = freq_peaks
            start_vals = [
                np.mean([start_vals[0][0], start_vals[1][0]]),  # offset
                np.sum([start_vals[0][3], start_vals[1][3]]),  # FWHM
                start_vals[0][1], start_vals[1][1],  # amplitudes
                start_vals[0][2], start_vals[1][2]  # centers
            ]
        try:
            if len(freq_peaks) == 2:
                if verbose:
                    print(('fit double peak with initial values', start_vals))

                fit = fit_double_lorentzian(freq, ampl, starting_params=start_vals, bounds=
                [(0, 0, -np.inf, -np.inf, min(freq), min(freq)), (np.inf, np.inf, 0, 0, max(freq), max(freq))])

            elif len(freq_peaks) == 1:
                if verbose:
                    print(('fit single peak with initial values', start_vals))

                fit = fit_lorentzian(freq, ampl, starting_params=start_vals,
                                     bounds=[(0, -np.inf, 0, 0), (np.inf, 0, np.inf, np.inf)])
                # if this detects a single peak far shifted, throw it out
                if np.abs(fit[2] - F0) > MAX_STRAIN:
                    fit = None
                    return fit
        except:
            # ESR fit failed!
            fit = None
            return fit
        # if offset is < 5 kCounts/sec, definitely all noise
        if fit[0] < min_counts:
            fit = None
            return fit

        # EXPERIMENTAL, REMOVE IF PROBLEMATIC
        # check that amplitude of at least one peak is greater than twice standard deviation (above the noise)
        if len(fit) == 6:
            if ((calc_esr_noise(freq, ampl, fit) * contrast_factor > np.abs(fit[2])) and (
                    calc_esr_noise(freq, ampl, fit) * contrast_factor > np.abs(fit[3]))):
                fit = None
            elif ((calc_esr_noise(freq, ampl, fit) * contrast_factor > min(np.abs(fit[2]), np.abs(fit[3]))) and (
                    calc_esr_noise(freq, ampl, fit) * contrast_factor < max(np.abs(fit[2]), np.abs(fit[3])))):
                if np.abs(fit[2]) > np.abs(fit[3]):
                    start_vals = [start_vals[0], start_vals[2], start_vals[4], start_vals[1]]
                else:
                    start_vals = [start_vals[0], start_vals[3], start_vals[5], start_vals[1]]
                try:
                    fit = fit_lorentzian(freq, ampl, starting_params=start_vals,
                                         bounds=[(0, -np.inf, 0, 0), (np.inf, 0, np.inf, np.inf)])
                    if calc_esr_noise(freq, ampl, fit) * contrast_factor > np.abs(fit[1]):
                        fit = None
                except:
                    fit = None

                # if this detects a single peak far shifted, throw it out
                if not fit is None and np.abs(fit[2] - F0) > MAX_STRAIN:
                    fit = None

        elif len(fit) == 4:
            if calc_esr_noise(freq, ampl, fit) * contrast_factor > np.abs(fit[1]):
                fit = None

        # if the width is exactly less than the minimum width, then it found a junk peak
        if not fit is None:
            if (len(fit) == 6 and int(fit[1]) <= MIN_WIDTH):
                fit = None
            elif (len(fit) == 4 and int(fit[3]) <= MIN_WIDTH):
                fit = None
        if not fit is None:
            if (len(fit) == 6 and int(fit[1]) >= MAX_WIDTH):
                fit = None
            elif (len(fit) == 4 and int(fit[3]) >= MAX_WIDTH):
                fit = None



        # #if the width is exactly equal to the minimum width, then it found a junk peak
        # if  not fit is None and (len(fit) == 6 and int(fit[1]) == MIN_WIDTH):
        #     fit = None

    return fit

def calc_esr_noise(freq, amp, fit_params):
    if fit_params is not None and fit_params[0] != -1:  # check if fit valid
        if len(fit_params) == 4:
            # single peak
            fit_data = lorentzian(freq, *fit_params)
        elif len(fit_params) == 6:
            # double peak
            fit_data = double_lorentzian(freq, *fit_params)
#         average_deviation = np.mean(np.abs(fit_data - amp))
        average_deviation = np.std(fit_data - amp)
        return average_deviation
    else:
        return None