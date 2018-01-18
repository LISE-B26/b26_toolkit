# Here we collect code that is a bit higher level and helps to understand nv related measurements


import fields as f
import nv_optical_response as nv
import fields_plot as fp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_best_NV_position(p, max_broadening=100, max_off_axis_field=0.01, nv_id=1, n=[0, 0, 1],
                         wo=500e-9, gammaNV=28e9,
                         arrow_length=1, arrow_angle=200, verbose=False, plot=False, plot_prop='shift'):
    """

    calculates the fields and gradients for parameters p
    and finds the position in space with the best conditions for the spin-mechanics experiment

    p: parameters - dictionary with following entries:
        tag: name identifier (string)
         a: radius in um
         Br: surface magnetization in Teslas
         phi_m: polar angle in deg
         theta_m: azimuthal angle in deg
         d_bead_z: distance top of bead to NV plane
         mu_0: vacuum permeability ( T m /A)
         d_bead_z: distance between bead and z plane
         dx: distance between points (in um)
         x_min, x_max, y_min, y_max: plot dimensions (in um)

    max_broadening = 100 # max broadening in MHz
    max_off_axis_field = 0.01 # max off axis field in Teslas

    nv_id =1 # select a NV [1,2,3,4]
    n = [0,0,1] # direction of gradient


    wo = 500e-9 # focal spot size
    gammaNV = 28e9 # (ESR shift is 28 GHz/T)


    arrow_length = 1
    arrow_angle = 200

    verbose = True
    plot: if true plot the NV shift

    """
    s = nv.nNV[nv_id - 1]  # NV orientation

    # =============== calculate the gradients ==============

    df = f.calc_Gradient_single_dipole(p, s, n, verbose=verbose)
    # calculate gradent along x
    data = f.calc_Gradient_single_dipole(p, s, n, verbose=verbose)
    # calculate gradent along y
    data2 = f.calc_Gradient_single_dipole(p, s, n, verbose=verbose)
    # now calculate the avrg gradient in xy
    df['Gxy'] = np.sqrt(data['G'] ** 2 + data2['G'] ** 2)
    # calcualte the broadening
    df['Broadening'] = df['Gxy'] * gammaNV * wo  # the linewidth broadening due to a gradient in the xy-plane in MHz
    # calculate the magetic field
    data = f.calc_B_field_single_dipole(p, verbose=verbose)
    # df['Bx']
    # on-axis field
    df['Bpar'] = np.abs(np.dot(np.array(data[['Bx', 'By', 'Bz']]), np.array([s]).T))
    # off-axis field
    df['Bperp'] = np.linalg.norm(np.cross(np.array(data[['Bx', 'By', 'Bz']]), np.array([s])), axis=1)
    # total field
    df['Bmag'] = np.linalg.norm(np.array(data[['Bx', 'By', 'Bz']]), axis=1)

    # convert fields into NV frame
    BNV = nv.B_fields_in_NV_frame(np.array(data[['Bx', 'By', 'Bz']]), nv_id)
    esr_freq = nv.esr_frequencies(BNV)
    df['fm'] = esr_freq[:, 0]
    df['fp'] = esr_freq[:, 1]

    if verbose:
        print('Calculated fields and gradients at {:d} points'.format(len(df)))

    # =============== find the best position ==============
    # get the points where the xy gradient is less than the specified limit
    x = df.loc[(np.abs(df['Broadening']) < max_broadening)]
    if verbose:
        print('Limited to xy inhomogeneous broadening < {:0.0f} MHz, {:d} points left'.format(max_broadening, len(x)))
    # get the points where the xy gradient is less than the specified limit
    x = x.loc[(np.abs(x['Bperp']) < max_off_axis_field)]
    if verbose:
        print('Limited to off axis field < {:0.0f} mT, {:d} points left'.format(max_off_axis_field * 1e3, len(x)))

    # out of the subset get the point with the highest gradient
    x = x.loc[np.abs(x['G']) == np.max(np.abs(x['G']))]

    # calculate some values for the arrow plot
    xo, yo = np.float(x['x']), np.float(x['y'])
    dx = arrow_length * np.cos(arrow_angle * np.pi / 180)
    dy = arrow_length * np.sin(arrow_angle * np.pi / 180)

    if plot:
        fig = fp.plot_NV_property_map(df, plot_prop)
        plt.arrow(xo + dx, yo + dy, -dx, -dy,
                  head_width=0.3, head_length=0.2, fc='w', ec='w', head_starts_at_zero=False,
                  length_includes_head=True)
        return x, fig
    else:

        return x