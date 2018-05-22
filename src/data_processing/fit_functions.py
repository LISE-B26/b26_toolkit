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

import numpy as np
from scipy import optimize


# ========= Gaussian fit functions =============================
#===============================================================
def fit_gaussian(x_values, y_values, starting_params=None, bounds=None):
    """

    Args:
        x_values: domain of fit function
        y_values: y-values to fit
        starting_params: reasonable guesses for where to start the fitting optimization of the parameters. This is a
        length 4 list of the form [constant_offset, amplitude, center, width].
        bounds: Optionally, include bounds for the parameters in the gaussian fitting, in the following form:
                ([offset_lb, amplitude_lb, center_lb, width_lb],[offset_ub, amplitude_ub, center_ub, width_ub])

    Returns:
        a length-4 list of [fit_parameters] in the form [constant_offset, amplitude, center, width]

    """

    if bounds:
        fit_params = optimize.curve_fit(gaussian, x_values, y_values, p0=starting_params, bounds=bounds, max_nfev=2000)[0]
    else:
        fit_params = optimize.curve_fit(gaussian, x_values, y_values, p0=starting_params)[0]

    # todo: catch if fit is not successful and return all zeros

    return fit_params
def gaussian(x, constant_offset, amplitude, center, width):
    return constant_offset + amplitude * np.exp(-1.0 * (np.square((x - center)) / (2 * (width ** 2))))

def guess_gaussian_parameter(x_values, y_values):
    """
    guesses the parameters for a Gaussian dataset
    Args:
        x_values:
        y_values:

    Returns: estimated fit parameters for Gaussian fit
    """
    # todo: find a smarter algorith, in particular for the width, for now this function is only used in autofocus, but might be good to generalize
    noise_guess = np.min(y_values)
    amplitude_guess = np.max(y_values) - noise_guess
    center_guess = x_values[np.argmax(y_values)]
    width_guess = 0.8

    return [noise_guess, amplitude_guess, center_guess, width_guess]



# ========= Lorenzian fit functions =============================
#===============================================================
def get_lorentzian_fit_starting_values(x_values, y_values, negative_peak=True):
    """
    estimates the parameter for a Lorentzian fit to the data set
    Note that the Lorentzian is assumed to
    Args:
        x_values:
        y_values:
        negative_peak: if peak is negative or positive
    Returns: estimated parameters as a list: [constant_offset, amplitude, center, fwhm]

    """

    constant_offset = np.mean(y_values)
    amplitude = min(y_values) - max(y_values) + 2 * np.std(y_values)  # min - max overestimates the amplitude, that's why we reduce the amplitude by the stddev
    if negative_peak is False:
        amplitude = -amplitude
    center = np.mean(x_values)
    # fwhm = (max(x_values) - min(x_values)) / 5  # assume that the peak is about 1/5th of the size of the image
    fwhm = .003e9
    return [constant_offset, amplitude, center, fwhm]


def fit_lorentzian(x_values, y_values, starting_params=None, bounds=None):
    """
    fits to lorenzian or two lorenzians: future fit to arbitrarily many lorenzians
    Args:
        x_values: domain of fit function
        y_values: y-values to fit
        starting_params: reasonable guesses for where to start the fitting optimization of the parameters. This is a
        length 4 list of the form [constant_offset, amplitude, center, full_width_half_max] or list of list of length 4
        which are the estimates for each peak.
        bounds: Optionally, include bounds for the parameters in the gaussian fitting, in the following form:
                [(offset_lb, amplitude_lb, center_lb, fwhm_lb), (offset_ub, amplitude_ub, center_ub, fwhm_ub)]

    Returns:
        a length-4 list of [fit_parameters] in the form [constant_offset, amplitude, center, fwhm]

    """

    # defines a lorentzian with amplitude, width, center, and offset to use with opt.curve_fit
    if bounds:
        return optimize.curve_fit(lorentzian, x_values, y_values, p0=starting_params, bounds=bounds, max_nfev=2000)[0]
    else:
        return optimize.curve_fit(lorentzian, x_values, y_values, p0=starting_params)[0]


def fit_double_lorentzian(x_values, y_values, starting_params=None, bounds=None):
    """
    fits to lorenzian or two lorenzians: future fit to arbitrarily many lorenzians
    Args:
        x_values: domain of fit function
        y_values: y-values to fit
        starting_params: reasonable guesses for where to start the fitting optimization of the parameters. This is a
        length 6 list of the form [constant_offset, fwhm, amplitude_1, amplitude_2, center_1, center_2]
        bounds: Optionally, include bounds for the parameters in the gaussian fitting, in the following form:
                [(offset_lb, fwhm_lb, amplitude1_lb, amplitude2_lb, center1_lb, center2_lb),
                (offset_ub, fwhm_ub, amplitude1_ub, amplitude2_ub, center1_ub, center2_ub)]

    Returns:
        a length-6 list of [fit_parameters] in the form [constant_offset, fwhm, amplitude_1, amplitude_2, center_1, center_2]

    """

    # defines a lorentzian with amplitude, width, center, and offset to use with opt.curve_fit
    if bounds:
        return optimize.curve_fit(double_lorentzian, x_values, y_values, p0=starting_params, bounds=bounds, max_nfev=2000)[0]
    else:
        return optimize.curve_fit(double_lorentzian, x_values, y_values, p0=starting_params)[0]


def lorentzian(x, constant_offset, amplitude, center, fwhm):
    """
    Lorentzian curve
    Args:
        x:  numpy array with x-coordinates
        constant_offset: float
        amplitude: float
        center: float
        fwhm: float

    Returns:  numpy array with y-values

    """
    return constant_offset + amplitude * np.square(0.5 * fwhm) / (np.square(x - center) + np.square(0.5 * fwhm))


def double_lorentzian(x, constant_offset, fwhm, amplitude_1, amplitude_2, center_1, center_2):
    """
    Two Lorentzian curves
    Args:
        x:
        constant_offset: float
        fwhm: float
        amplitude_1:
        amplitude_2: float
        center_1: float
        center_2: float

    Returns: numpy array with y-values

    """
    return lorentzian(x, constant_offset / 2, amplitude_1, center_1, fwhm) + lorentzian(x, constant_offset / 2,
                                                                                        amplitude_2, center_2, fwhm)

# ========= Cose fit functions =============================
#===============================================================
def get_ampfreqphase_FFT(qx, dt, n0 = 0, f_range = None, return_Spectra = False):
    '''
    returns estimate of amplitdue, frequency and phase from FFT

    [ax, wx, phi] = get_ampfreqphase_FFT(qx, dt,n0 = 0, f_range=None, return_Spectra = False)
    [ax, wx, phi], [Fx, Ax] = get_ampfreqphase_FFT(qx, dt,n0 = 0, f_range=None, return_Spectra = True)
    input:
        qx: time trace  sampled at intervals dt
        dt: sampling interval

    input (optional):
        n0 = t0/dt: index of time zero
        f_range = [f_x, df]: frequency is looked in intervals f_x +-df respectively
        return_Spectra = True/False: returns spectra over range f_range in addition to [phi, ax, fx]

    output:
        dominant angular frequency, amplitude at that frequency and phase
        method: get fourier component of max signals
    '''

    n = len(qx)
    f = np.fft.fftfreq(n, dt)[0:int(n/2)]

    # look for max frequencies only in certain range
    if f_range is None:
        irange_x = np.arange(int(n/2))
    else:
        [f_x, df] = f_range
        imin = np.argwhere(f >= f_x-df)[0,0]
        imax = np.argwhere(f <= f_x+df)[-1,0] + 1
        irange_x = np.arange(imax-imin+1)+imin

    # convert to int (from float)
    irange_x = [int(x) for x in irange_x]

    # Fourier transforms (remove offset, in case there is a large DC)
    Ax = np.fft.fft(qx-np.mean(qx))[irange_x] / n*2
    Fx = f[irange_x]

    # frequency and amplitude x
    i_max_x = np.argmax(np.abs(Ax))
    fx = Fx[i_max_x]
    ax = np.abs(Ax[i_max_x])
    # phase
    phi = np.angle(Ax[i_max_x] * np.exp(-1j *2 * np.pi * fx * n0))

    if return_Spectra == True:
        return [ax, 2*np.pi*fx, phi], [Fx, Ax]
    else:
        return [ax, 2*np.pi*fx, phi]
def A_fun(qx, w, dt):
    '''
    Ak = A_fun(qx, w, fs)
    input:
        xx: input signal vector length N
        w: omega, w = 2*pi*k*fs / M  vector length K
        dt: sampling interval
    output:
        Ak: spectrum at frequencies w
    '''

    N = len(qx)
    j = 1j

    nn = np.array( [ np.arange(N) ] ).transpose()

    Ak = (1./N) * np.dot(np.array([qx]),  np.exp(-j * np.dot( nn, np.array([w]) ) * dt) )

    return Ak[0]

def guess_cose_parameter(t, y):
    """
    guesses the parameters for a cosinus dataset
    Args:
        t: time vector, here we assume that t is sampled evenly
        y: data vector

    Returns: estimated fit parameters for Sine with offset fit
    """
    dt = np.mean(np.diff(t))
    offset = float(max(y) + min(y)) / 2
    [ax, wx, phi] = get_ampfreqphase_FFT(y-offset, dt)
    

    # if the oscillation is less than a peroiod we take the average of the min and max as the offset otherwise we take the mean
    if max(t) < 2*np.pi/wx:
        offset = float(max(y) + min(y)) / 2
    else:
        offset = np.mean(y)
        ax = np.mean([ax, np.std(y)])# take an average off std and FFT amplitude

    return [ax, wx, phi, offset]


def cose(t, a0, w0, phi0, offset):
    """
        cosine function
    """
    return a0*np.cos(w0*t+phi0)+ offset

def fit_cose_parameter(t, y, verbose = False):
    """
    fits the data to a cosine
    Args:
        t:
        y:

    Returns: [ax, wx, phi, offset] = [amplitude, angular frequency, phase and offset]

    """
    [ax, wx, phi, offset] = guess_cose_parameter(t, y)
    if verbose:
        print(('initial estimates [ax, wx, phi, offset]:', [ax, wx, phi, offset]))
    def cost_function_fit(x):
        """
        cost function for fit to sin
        """
        ao, wo, po, offset = x
        #         sig = sine(x, t)
        sig = ao * np.cos(wo * t + po) + offset
        return np.sum((sig - y)**2)

    opt = optimize.minimize(cost_function_fit, [ax, wx, phi, offset])

    if verbose:
        print(('optimization result:', opt))
    [ax, wx, phi, offset] = opt.x

    return [ax, wx, phi, offset]

def cose_with_decay(t, a0, w0, phi0, offset, tau):
    """
        cosine function
    """
    return a0*np.exp(-t/tau)*np.cos(w0*t+phi0)+ offset


def get_decay_data(t, y, wo, verbose=False):
    """
        average the data y over a oscillation period to smoothout oscillations
    returns: averaged data

    """

    period = 2 * np.pi / wo
    dt = np.mean(np.diff(t))
    index_per_interval = int(period / dt)

    number_of_oscillations = int(np.floor(len(y) / index_per_interval))

    if verbose:
        print((
        'initial estimates [index_per_interval, number_of_oscillations]:', [index_per_interval, number_of_oscillations]))

    decay_y = np.array(
        [np.std(y[index_per_interval * i:index_per_interval * (i + 1)]) for i in range(number_of_oscillations)])
    decay_t = np.array(
        [np.mean(t[index_per_interval * i:index_per_interval * (i + 1)]) for i in range(number_of_oscillations)])
    return np.array(decay_t), np.array(decay_y)

def fit_exp_decay(t, y, offset = False, verbose=False):
    """
    fits the data to a decaying exponential, with or without an offset
    Args:
        t: x data
        y: y data
        offset: False if fit should decay to y=0, True otherwise
        verbose: prints results to screen

    Returns: fit parameters, either [ao, tau, offset] if offset is True, or or [ao, tau] if offset is False
            ao: amplitude above offset (or zero if offset is False)
            tau: decay parameter
            offset: asymptotic value as t->INF

    """
    if verbose:
        print(' ======= fitting exponential decay =======')

    init_params = estimate_exp_decay_parameters(t, y, offset)
    if offset:
        [ao, tau, offset] = optimize.curve_fit(exp_offset, t, y, p0=init_params)[0]
    else:
        [ao, tau] = optimize.curve_fit(exp, t, y, p0=init_params)[0]

    if offset:
        if verbose:
            print(('optimization result:', [ao, tau, offset]))
        return [ao, tau, offset]
    else:
        if verbose:
            print(('optimization result:', [ao, tau]))
        return [ao, tau]

def estimate_exp_decay_parameters(t,y,offset):
    '''
    Returns an initial estimate for exponential decay parameters. Meant to be used with optimize.curve_fit.
    Args:
        t: x data
        y: y data
        offset: False if fit should decay to y=0, True otherwise

    Returns: fit parameter estimate, either [ao, tau, offset] if offset is True, or or [ao, tau] if offset is False
            ao: amplitude above offset (or zero if offset is False)
            tau: decay parameter
            offset: asymptotic value as t->INF

    '''
    if offset:
        offset = y[-1]
    else:
        offset = 0
    total_amp = y[0]
    ao = total_amp-offset
    decay = t[np.argmin(np.abs(y - (total_amp+offset)/2))] #finds time at which the value is closest to midway between the max and min
    if offset:
        return [ao, decay, offset]
    else:
        return [ao, decay]

def exp(t, ao, tau):
    '''
    Exponential decay: ao*E^(t/tau)
    '''
    return np.exp(-t / tau) * ao

def exp_offset(t, ao, tau, offset):
    '''
    Exponential decay with offset: ao*E^(t/tau) + offset
    '''
    return np.exp(-t / tau) * ao + offset


def fit_rabi_decay(t, y, variable_phase=False, verbose=False, return_guess = False):
    """
    fit to a cosine with an exponential envelope
    Args:
        t: time in ns
        y: counts in kilocounts
        variable_phase: if true the phase is a free parameter if false the phase is 0 (cosine)
        return_guess: return also the initial guess parameters that are used in the fit

    """

    [ax, wx, phi, offset] = guess_cose_parameter(t, y)

    # estimate the decay constant
    t_decay, y_decay = get_decay_data(t, y, wx, verbose)
    [_, to] = fit_exp_decay(t_decay, y_decay, )

    if variable_phase:
        # added by ER 7.27.17 to make Rabi frequency from fit always positive
        initial_parameter = [ax, abs(wx), phi, offset, to]
    else:
        initial_parameter = [ax, abs(wx), offset, to]
    verbose = 1
    if verbose:
        if variable_phase:
            print(('initial estimates [ax, wx, phi, offset, tau]:', initial_parameter))
        else:
            print(('initial estimates [ax, wx, offset, tau]:', initial_parameter))

    def cost_function_fit(x):
        """
        cost function for fit to exponentially decaying cosine
        """
        if variable_phase:
            ao, wo, po, offset, to = x
            # added by ER 7.27.17 to make Rabi frequency from fit always positive
            wo = abs(wo)
            sig = cose_with_decay(t, ao, wo, po, offset, to)
            # sig = ao * np.exp(-t / to) * np.cos(wo * t + po) + offset
        else:
            ao, wo, offset, to = x
            # added by ER 7.27.17 to make Rabi frequency from fit always positive
            wo = abs(wo)
            sig = cose_with_decay(t, ao, wo, 0, offset, to)
            # sig = ao * np.exp(-t / to) * np.cos(wo * t) + offset
        return np.sum((sig - y) ** 2)

    opt = optimize.minimize(cost_function_fit, initial_parameter)
    #     opt = optimize.minimize(cost_function_fit, [ax, wx, phi, offset, to],
    #                             bounds=[(None, None), (1.1*wx, 2*wx), (None, None), (None, None), (None, None)])

    #
    # [ax, wx, phi, offset, tau] = opt.x

    if verbose:
        print(('optimization result:', opt))
    if return_guess:
        return opt.x, initial_parameter
    else:
        return opt.x
