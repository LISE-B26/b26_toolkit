# Here we collect code that is a bit higher level and helps to understand nv related measurements


import fields as f
import nv_optical_response as nv
import fields_plot as fp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



def z_rotation(phi):
    """
    rotation matrix for a rotation about the z axis
    phi: rotation angle in degree
    """
    phi/=180/np.pi
    print(phi)
    return np.array([
            [np.cos(phi), -np.sin(phi), 0],
            [np.sin(phi), np.cos(phi), 0],
            [0,0,1]
        ]).T

def get_full_nv_dataset(p, nv_id=1, n=[0, 0, 1], nv_rotation_matrix = None, wo=500e-9, gammaNV=28e9, verbose=False):
    """
    returns a full dataset for Nv number nv_id, based on parameters p

    :param p:
    :param nv_id:
    :param n:
    :param nv_rotation_matrix: matrix that rotates the coordinate system of the nv center, e.g. to account for rotations of the diamond wrt the resonator
    :param wo:
    :param gammaNV:
    :param verbose:
    :return:
    """
    s = nv.nNV[nv_id - 1]  # NV orientation
    if nv_rotation_matrix is not None:
        assert np.shape(nv_rotation_matrix) == (3,3)
        s = np.dot(s, nv_rotation_matrix)
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

    return df


def get_best_NV_position(df, max_broadening=100, max_off_axis_field=0.01, verbose=False):
    """

    calculates the fields and gradients for parameters p
    and finds the position in space with the best conditions for the spin-mechanics experiment

    df: dataframe with all the fields and gradients calculated for Nv with nv_id

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



    return x