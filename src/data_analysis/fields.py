import numpy as np
import pandas as pd
from copy import deepcopy
import time

def calcBfield_oommf(rs, data, info, use_parallel = True, verbose = False):
    '''
    Calculate the magnetic field for a collection of dipoles, where the dipoles have been calculated by OOMMF
    rs: (matrix with dimension m x 3), positions at which field is evaluated in space (in m)
    data: dataframe with columns 'mx', 'my', 'mz', 'x', 'y', 'z' that gives the dipolevector and its location
    info: dictionary with metadata for the dataset, contains 'xstepsize', 'ystepsize', 'zstepsize', which give the spacing of the dipole locations
    use_parallel:  (boolean) if True use parallel execution of code this is not working yet....

    :returns pandas dataframe columns 'Bx', 'By', 'Bz', 'x', 'y', 'z' that gives the fieldvector and its location (=rs)
    '''

    dV = info['xstepsize'] * info['ystepsize'] * info['zstepsize'] * 1e18  # cell volume in um^3

    print(('length of data', len(data)))

    # pick only the data where the magnetization is actually non-zero
    data = data[np.sum(data[['mx', 'my', 'mz']].as_matrix()**2,1)>0]


    DipolePositions = data[['x', 'y', 'z']].as_matrix() * 1e6  # convert from m to um
    m = data[['mx', 'my', 'mz']].as_matrix() * dV  # multiply by the cell volume to get the magnetic dipole moment (1e-6 A um^2 = 1e-18 J/T)


    rs *=1e6# convert from m to um

    B = b_field(rs, DipolePositions, m, use_parallel=use_parallel, verbose=verbose)

    return B

    #
    # print('number of magnetic moments', len(data))
    # print('number of positions', len(rs))
    # if use_parallel:
    #     # try importing the multiprocessing library
    #     from joblib import Parallel, delayed
    #     import multiprocessing
    #     num_cores = multiprocessing.cpu_count()
    #
    #     print('using ', num_cores, ' cores')
    #
    # def process(r):
    #     return b_field_single_pt(r, DipolePositions, m)
    #
    #
    #
    # # convert from um to m
    # data_out = {
    #     'x':deepcopy(rs[:,0]*1e-6),
    #     'y':deepcopy(rs[:,1]*1e-6),
    #     'z':deepcopy(rs[:,2]*1e-6)
    # }
    #
    #
    #
    # if use_parallel:
    #
    #     # Parallel(n_jobs=1)(delayed(sqrt)(i ** 2) for i in range(10))
    #     # B = Parallel(n_jobs=num_cores)(delayed(b_field_single_pt)(r, DipolePositions, m) for r in rs)
    #     B = Parallel(n_jobs=num_cores)(delayed(b_field_single_pt)(r, DipolePositions, m) for r in rs)
    #     # B = Parallel(n_jobs=num_cores)(delayed(process)(r) for r in rs)
    #     B = np.array(B)
    #
    # else:
    #     B = np.array([b_field_single_pt(r, DipolePositions, m) for r in rs])
    #
    # # put data into a dictionary
    # data_out['Bx'] = deepcopy(B[:, 0])
    # data_out['By'] = deepcopy(B[:, 1])
    # data_out['Bz'] = deepcopy(B[:, 2])


    # return data as a pandas dataframe
    # return pd.DataFrame.from_dict(data_out)

def b_field_single_pt(r, DipolePositions, m, mu0 =4 * np.pi * 1e-7):
    """
    calculates the magnetic field at position r
    :param r: vector of length 3 position at which field is evaluates (in um)
    :param DipolePositions: matrix Nx3, of positions of dipoles (in um)
    :param m:  matrix Nx3, components dipole moment at position DipolePositions mx, my, mz (in 1e-18 J/T)
    mu0 = 4 * np.pi * 1e-7  # T m /A
    :return:
    """

    # check that DipolePositions and m have the same shape
    assert np.shape(DipolePositions) == np.shape(m)
    #
    # # check that r is a vector of length 3
    assert len(np.shape(r)) == 1
    assert len(r) == 3

    a = np.ones((np.shape(DipolePositions)[0], 1)) * np.array([r]) - DipolePositions #

    rho = np.sqrt(np.sum(a ** 2, 1))

    # if we request thefield at the location of the dipole the field diverges, thus we exclude this value because we only want to get the fields from all the other dipoles
    zero_value_index = np.argwhere(rho == 0)

    rho = np.array([rho]).T * np.ones((1, 3))

    # calculate the vector product of m and a: m*(r-ri)
    ma = np.array([np.sum(m * a, 1)]).T * np.ones((1, 3))
    B = mu0 / (4 * np.pi) * (3. * a * ma / rho ** 5 - m / rho ** 3)  # magnetic field in Tesla

    # exclude the dipole at the location where we calculate the field
    if len(zero_value_index) > 0:
        B[zero_value_index, :] = np.zeros([len(zero_value_index), 3])

    return np.sum(B, 0)

def b_field(rs, DipolePositions, m, use_parallel = True, verbose = False):
    '''
    calculates the magnetic field at multiple positions r
    :param rs:  matrix Mx3, position at which field is evaluated (in um)
    :param DipolePositions: matrix Nx3, of positions of dipoles (in um)
    :param m:  matrix Nx3, components dipole moment at position DipolePositions mx, my, mz (in 1e-18 J/T)
    use_parallel:  (boolean) if True use parallel execution of code
    :param verbose: if True print information as script is executed
    :returns pandas dataframe columns 'Bx', 'By', 'Bz', 'x', 'y', 'z' that gives the fieldvector and its location (=rs)
    '''



    # check that DipolePositions and m have the same shape
    assert np.shape(DipolePositions) == np.shape(m)
    #
    # # check that r is a vector of length 3
    assert np.shape(rs)[1] == 3

    if verbose:
        print(('number of magnetic moments', len(m)))
        print(('number of positions', len(rs)))
    if use_parallel:
        # try importing the multiprocessing library
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if verbose:
            print(('using ', num_cores, ' cores'))


    if use_parallel:

        B = Parallel(n_jobs=num_cores)(delayed(b_field_single_pt)(r, DipolePositions, m) for r in rs)
        B = np.array(B)

    else:
        B = np.array([b_field_single_pt(r, DipolePositions, m) for r in rs])

    if verbose:
        print(('rs shape', np.shape(rs)))
        # print('rs', rs)

    # put data into a dictionary
    data_out = {
        'x':deepcopy(rs[:,0]),
        'y':deepcopy(rs[:,1]),
        'z':deepcopy(rs[:,2]),
        'Bx':deepcopy(B[:, 0]),
        'By': deepcopy(B[:, 1]),
        'Bz': deepcopy(B[:, 2]),
    }


    # return data as a pandas dataframe
    return pd.DataFrame.from_dict(data_out)

def gradient_single_pt(r, DipolePositions, m, s, n, verbose = False, mu0 =4 * np.pi * 1e-7):
    '''
    calculates the magnetic field at position r
    :param r: vector of length 3 position at which field is evaluates (in um)
    DipolePositions: (matrix with dimension N x 3) m dipole location in space in nm
    m:  (matrix with dimension N x 3) magnetic moment of a dipole in 1e-18J/T
    s: (vector of length 3) spin vector no units
    n: (vector of length 3) projection vector of the gradient, e.g. motion of resonator
    mu0 = 4 * np.pi * 1e-7  # T m /A
    Output in T/um
    '''

    # check that DipolePositions and m have the same shape
    assert np.shape(DipolePositions) == np.shape(m)

    if len(np.shape(m)) == 1:
        DipolePositions = np.array([DipolePositions])
        m = np.array([m])
    assert np.shape(m)[1]==3
    N = len(m)

    #
    # # check that r, s and n are vectors of length 3
    assert len(np.shape(r)) == 1
    assert len(r) == 3
    assert len(np.shape(s)) == 1
    assert len(s) == 3
    assert len(np.shape(n)) == 1
    assert len(n) == 3


    a = np.ones((N,1)) * np.array([r])-DipolePositions
    rho = np.sqrt(np.sum(a**2,1))

    # calculate the vector product of m and a: m*(r-ri)
    ma = np.sum(m * a, 1)
    # calculate the vector product of s and a: s*(r-ri)
    sa = np.sum(np.ones((N,1)) * np.array([s])*a,1)
    # calculate the vector product of n and a: n*(r-ri)
    na = np.sum(np.ones((N,1)) * np.array([n])*a,1)
    # calculate the vector product of s and n
    sn = np.dot(s, n)
    # calculate the vector product of m and n
    mn = np.sum(np.ones((N, 1)) * np.array([n]) * m, 1)
    # calculate the vector product of m and s
    ms = np.sum(np.ones((N, 1)) * np.array([s]) * m, 1)

    if verbose:
        print(('rho', rho))
        print(('a', a))
        print(('ma', ma))
        print(('sa', sa))
        print(('na', na))

    gradB = 3. * mu0 / (4 * np.pi * (rho)** 5) * (
            ma * sn + sa * mn + ms * na
            - 5 * (sa * ma / rho ** 2) * na
    )
    if verbose:
        print(('sn', sn))
        print(('mn', mn))
        print(('ms', ms))

    if verbose:
        print(('grad due to every dipole', gradB))

    return np.sum(gradB,0)

def gradient(rs, DipolePositions, m, s, n, use_parallel=True, verbose=False):
    '''
    Calculate the magnetic field gradient for a collection of dipoles at multiple positions r
    :param rs:  matrix Mx3, position at which field is evaluated (in um)
    :param DipolePositions: matrix Nx3, of positions of dipoles (in um)
    :param m:  matrix Nx3, components dipole moment at position DipolePositions mx, my, mz (in 1e-18 J/T)
    :param s: direction of field component for which the gradient is evaluates (e.g. direction of NV center)
    :param n: direction of gradient (e.g. direction of resonator motion)
    use_parallel:  (boolean) if True use parallel execution of code
    :param verbose: if True print information as script is executed
    :returns pandas dataframe columns 'Gx', 'By', 'Bz', 'x', 'y', 'z' that gives the fieldvector and its location (=rs)

    rs: (matrix with dimension M x 3), positions at which field is evaluated in space (in um)
    M: dataframe with columns 'mx', 'my', 'mz', 'x', 'y', 'z' that gives the dipolevector and its location
    info: dictionary with metadata for the dataset, contains 'xstepsize', 'ystepsize', 'zstepsize', which give the spacing of the dipole locations
    use_parallel:  (boolean) if True use parallel execution of code this is not working yet....

    :returns pandas dataframe columns 'G', 'x', 'y', 'z' that gives the Gradient (T/um) of component s along n and its location (=rs)
    '''

    if verbose:
        print(('number of magnetic moments', len(data)))
        print(('number of positions', len(rs)))
    if use_parallel:
        # try importing the multiprocessing library
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if verbose:
            print(('using ', num_cores, ' cores'))


    if use_parallel:
        G = Parallel(n_jobs=num_cores)(delayed(gradient_single_pt)(r, DipolePositions, m, s, n) for r in rs)
        G = np.array(G)
    else:
        G = np.array([gradient_single_pt(r, DipolePositions, m) for r in rs])

    # put data into a dictionary
    data_out = {
        'x': deepcopy(rs[:, 0]),
        'y': deepcopy(rs[:, 1]),
        'z': deepcopy(rs[:, 2]),
        'G': G
    }




    # return data as a pandas dataframe
    return pd.DataFrame.from_dict(data_out)

def b_field_single_dipole(r, DipolePosition, m, mu0 =4 * np.pi * 1e-7, verbose = False):
    """
    calculates the magnetic field at position r
    :param r: matrix Nx3, positions at which field is evaluates (in um)
    :param DipolePosition: vector of length 3 of positions of dipoles (in um)
    :param m:  vector of length 3 components dipole moment at position DipolePosition mx, my, mz (in 1e-18 J/T)
    mu0 = 4 * np.pi * 1e-7  # T m /A
    :return:
    """

    # check that DipolePositions and m have the same shape
    assert np.shape(DipolePosition) == np.shape(m)
    #
    # # check that m is a vector of length 3
    assert len(np.shape(m)) == 1
    assert len(m) == 3

    a = np.ones((len(r), 1)) * np.array([DipolePosition]) - r  #

    rho = np.sqrt(np.sum(a ** 2, 1))

    rho = np.array([rho]).T * np.ones((1, 3))

    # calculate the vector product of m and a: m*(r-ri)
    ma = np.array([np.sum(m * a, 1)]).T * np.ones((1, 3))

    if verbose:
        print(('rho', rho))
        print(('a', a))
        print(('ma', ma))

    B = mu0 / (4 * np.pi) * (3. * a * ma / rho ** 5 - m / rho ** 3)  # magnetic field in Tesla


    return B

def gradient_single_dipole(r, DipolePosition, m, s, n, verbose = False, mu0 =4 * np.pi * 1e-7):
    '''
    calculates the magnetic gradient field at position r
    :param r: vector of length 3 position at which field is evaluates (in um)
    DipolePosition: (vector of length 3) dipole location in space in nm
    m:  (vector of length 3)  magnetic moment of a dipole in 1e-18J/T
    s: (vector of length 3) spin vector no units
    n: (vector of length 3) projection vector of the gradient, e.g. motion of resonator
    mu0 = 4 * np.pi * 1e-7  # T m /A
    Output in T/um
    '''
    # check that DipolePositions and m have the same shape
    assert np.shape(DipolePosition) == np.shape(m)
    #
    # # check that m is a vector of length 3
    assert len(np.shape(m)) == 1
    assert len(m) == 3

    # # check that m, s and n are vectors of length 3
    assert len(np.shape(m)) == 1
    assert len(m) == 3
    assert len(np.shape(s)) == 1
    assert len(s) == 3
    assert len(np.shape(n)) == 1
    assert len(n) == 3

    # normalize unit vectors
    s/=np.linalg.norm(s)
    n /= np.linalg.norm(n)

    a = r- np.ones((len(r), 1)) * np.array(DipolePosition)  #

    rho = np.sqrt(np.sum(a ** 2, 1))

    # calculate the vector product of m and a: m*(r-ri)
    # ma = np.array([np.sum(m * a, 1)]).T * np.ones((1, 3))
    # ma = np.array([np.sum(m * a, 1)]).T
    ma = np.sum(m * a, 1).T
    # calculate the vector product of s and a: s*(r-ri)
    # sa = np.array([np.sum(s * a, 1)]).T * np.ones((1, 3))
    # sa = np.array([np.sum(s * a, 1)]).T
    sa = np.sum(s * a, 1).T
    # calculate the vector product of n and a: n*(r-ri)
    # na = np.array([np.sum(n * a, 1)]).T * np.ones((1, 3))
    # na = np.array([np.sum(n * a, 1)]).T
    na = np.sum(n * a, 1).T

    if verbose:
        print(('rho', rho))
        print(('a', a))
        print(('ma', ma))
        print(('sa', sa))
        print(('na', na))

    # a = np.ones((N,1)) * np.array([r])-DipolePosition
    # rho = np.sqrt(np.sum(a**2,1))

    # # calculate the vector product of m and a: m*(r-ri)
    # ma = np.sum(m * a, 1)
    # # calculate the vector product of s and a: s*(r-ri)
    # sa = np.sum(np.ones((N,1)) * np.array([s])*a,1)
    # # calculate the vector product of n and a: n*(r-ri)
    # na = np.sum(np.ones((N,1)) * np.array([n])*a,1)
    # calculate the vector product of s and n
    sn = np.dot(s, n)
    # calculate the vector product of m and n
    mn = np.dot(m, n)
    # mn = np.sum(np.ones((N, 1)) * np.array([n]) * m, 1)
    # calculate the vector product of m and s
    ms = np.dot(m, s)
    # ms = np.sum(np.ones((N, 1)) * np.array([s]) * m, 1)

    if verbose:
        print(('sn', sn))
        print(('mn', mn))
        print(('ms', ms))

    gradB = 3. * mu0 / (4 * np.pi * (rho)** 5) * (
            ma * sn + sa * mn + ms * na
            - 5 * (sa * ma / rho ** 2) * na
    )

    gradB = 3. * mu0 / (4 * np.pi * (rho)** 5) * (
            ma * sn + sa * mn + ms * na
            - 5 * (sa * ma / rho ** 2) * na
    )

    return gradB

def calc_B_field_single_dipole(p, verbose = False):
    """

    calculate the magnetic fields along a direction s for field component n for a dipole that is located at the origin

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

     :returns
      pandas dataframe with columns 'x', 'y', 'z' and 'Bx', 'By' and 'Bz'
    """

    r, M = p_to_positions(p)
    DipolePosition = np.zeros(3)  # we assume that the magnet is at 0,0,0

    start = time.time()
    data_out = b_field_single_dipole(r, DipolePosition, M, mu0=4 * np.pi * 1e-7)

    end = time.time()
    if verbose:
        print(('duration: {:0.2f} min'.format((end - start) / 60)))


    # create a pandas dataset
    data_out = pd.DataFrame.from_dict(
        {'x': r[:,0],'y':r[:,1], 'z':r[:,2],
         'Bx': data_out[:,0], 'By': data_out[:,1], 'Bz': data_out[:,2]
        }
    )

    # if filename not None:
    #     # save data to csv
    #     out_file = os.path.join('data/', filename)
    #     pd.DataFrame.from_dict(data_out).to_csv(out_file, index=False)

    return data_out

def calc_Gradient_single_dipole(p, s, n, verbose = False):
    """
    calculate the gradients along a direction n for field component s for a dipole that is located at the origin

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
    s: (vector of length 3) spin vector no units
    n: (vector of length 3) projection vector of the gradient, e.g. motion of resonator
    """

    r, M = p_to_positions(p)
    DipolePositions = np.zeros(3)  # we assume that the magnet is at 0,0,0

    start = time.time()
    data_out = gradient_single_dipole(r, DipolePositions, M, s, n)
    end = time.time()
    if verbose:
        print(('duration: {:0.2f} min'.format((end - start) / 60)))

    # create a pandas dataset
    data_out = pd.DataFrame.from_dict(
        {'x': r[:,0],'y':r[:,1], 'z':r[:,2],
         'G': data_out
        }
    )


    return data_out

def p_to_positions(p):
    """
    calculate the positions r where to calculate the fields and the moment vector if the dipole
    p: parameters - dictionary with following entries:
    a: radius in um
    Br: surface magnetization in Teslas
    phi_m: polar angle in deg
    theta_m: azimuthal angle in deg
    d_bead_z: distance top of bead to NV plane
    mu_0: vacuum permeability ( T m /A)
    d_bead_z: distance between bead and z plane
    dx: distance between points (in um)
    x_min, x_max, y_min, y_max: plot dimensions (in um)

    exclude_ring: (optional if existent describes the radius of the ring to be excluded from the positions)


    :returns positions r (in um) and magnetic moment M in (in 1e-18 J/T)
    """
    # calculate the magnetic moment
    M = 4 * np.pi / 3 *p['a']**3* p['Br'] / p['mu_0']
    phi_m = p['phi_m'] * np.pi / 180
    theta_m = p['theta_m'] * np.pi / 180
    M = M * np.array([
        np.cos(phi_m) * np.sin(theta_m),
        np.sin(phi_m) * np.sin(theta_m),
        np.cos(theta_m)

    ])

    #     xmax, ymax = 3,3 # extend in um
    #     xmin, ymin = None,None
    dx = p['dx']

    xmax = p['xmax']
    if 'xmin' not in p:
        xmin = -xmax
    else:
        xmin = p['xmin']
    if 'ymin' not in p:
        ymin = xmin
    else:
        ymin = p['ymin']
    if 'ymax' not in p:
        ymax = xmax
    else:
        ymin = p['ymin']

    zo = p['d_bead_z'] + p['a']

    # calculate the grid
    x = np.arange(xmin, xmax+dx, dx)
    y = np.arange(ymin, ymax+dx, dx)
    Nx, Ny = len(x), len(y)
    X, Y = np.meshgrid(x, y)

    r = np.array([X.flatten(), Y.flatten(), zo * np.ones(len(X.flatten()))]).T

    # if exclude ring is true the we throw away the positions within a radius of exclude_ring
    if 'exclude_ring' in p:
        r = np.array(r[r[:, 0] ** 2 + r[:, 1] ** 2 >= p['exclude_ring']** 2])

    return r, M

def magnetic_moment(a, Br, muo = 4 * np.pi * 1e-7):
    """
    calculates the magentic moment for a sphere with radius a and surface field Br

    :param a: sphere  radius in um
    :param Br: surface field in Teslas
    :param muo: 4 * np.pi * 1e-7
    :return: magentic moment
    """
    M = 4 * np.pi / 3 * a ** 3 * Br /muo
    return M

def p_to_filename(p):
    """
    create filename from parameters

    if parameters are given in a dictionary "p" call as follows: p_to_filename(*p)

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
    """

    filename = '{:s}_a_{:0.1f}um_phi_m_{:2.0f}deg_theta_m_{:2.0f}deg'.format(p['tag'], p['a'], p['phi_m'], p['theta_m'])

    return filename

def field_component(data, component_name = None, s = None):
    """
    returns the field component defined by component_name
    :param data:
    :param component_name:
    :param s:
    :return:
    """

    Bx, By, Bz = data['Bx'].as_matrix(), data['By'].as_matrix(), data['Bz'].as_matrix()

    if component_name is None or component_name == 'Bfield_mag':
        D = 1e4 * np.sqrt(Bx**2 + By**2 + Bz**2)
        label = '$|\mathbf{B}|$ (Gauss)'
    elif component_name in ('Bfield_proj', 'Bfield_par', 'Bfield_long', 'parallel'):
        D = 1e4 * (Bx * s[0] + By * s[1] + Bz * s[2])
        # label = '$\mathbf{B}\cdot \mathbf{S}$ (Gauss)'
        label = '$B_{\parallel}$ (Gauss)'
    elif component_name in ('Bfield_perp', 'Bfield_trans', 'perpendicular'):
        D = 1e4 * np.sqrt((By * s[2] - Bz * s[1]) ** 2 + (Bz * s[0] - Bx * s[2]) ** 2 + (Bx * s[1] - By * s[0]) ** 2)
        label = '$B_{\perp}$ (Gauss)'

    return D, label

#
# def unit_vector(theta, phi):
#
#     return [np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)]


if __name__ == '__main__':
    p = {
        'tag': 'bead_1',
        'a': 1.4,
        'Br': 0.4,
        'phi_m': 90,
        'theta_m': 90,
        'mu_0': 4 * np.pi * 1e-7,
        'd_bead_z': 0,
        'dx': 0.2,
        'xmax': 10
    }


    s = np.array([1, 1, 1])
    n = np.array([1, 0, 0])

    data = calc_Gradient_single_dipole(p, s, n)

    # import read_write as rw
    # import os
    # import datetime
    #
    # # folder = 'Z:\Lab\Cantilever\tmp_jan\oommf_results_2'
    # folder = ''
    #
    # file_mag = os.path.join(folder , 'random_K1_length_1.5um-Oxs_MinDriver-Magnetization-00000-0001715-omf.tsv')
    # file_H = os.path.join(folder , 'random_K1_length_1.5um-Oxs_CGEvolve-H-00000-0001715-ovf.tsv')
    #
    # print(os.path.exists(file_H))
    #
    # data_mag, info_mag = rw.load_ommf_vect_data(file_mag)
    # data_mag.head()
    #
    # data_H, info_H = rw.load_ommf_vect_data(file_H)
    # data_H.head()
    #
    # zo = -5e-9
    # subdata_H = rw.get_slice(data_H, zo, info_mag, 'z')
    # subdata_mag = rw.get_slice(data_mag, zo, info_mag, 'z')
    #
    # r = subdata_H[['x', 'y', 'z']].as_matrix()
    #
    # r = r[0:10, :]
    # print('shape r', np.shape(r))
    #
    #
    #
    # t1 =datetime.datetime.now()
    # dataB = b_field(r, subdata_mag, info_mag, True)
    # t2 = datetime.datetime.now()
    # # dataB = b_field(r, subdata_mag, info_mag, False)
    # # t3 = datetime.datetime.now()
    #
    # print('excution time parallel', str(t2-t1))
    # # print('excution time not parallel', str(t3 - t2))
