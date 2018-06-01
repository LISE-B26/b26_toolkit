import numpy as np
import scipy.optimize as opt
import itertools
from itertools import permutations
# spin 1 matrices
Sx = np.matrix([
    [0, -1j, -1j],
    [1j, 0, 0],
    [1j, 0, 0]
])/np.sqrt(2)

Sy = np.matrix([
    [0, -1, 1],
    [-1, 0, 0],
    [1, 0, 0]
])/np.sqrt(2)

Sz = np.matrix([
    [0, 0, 0],
    [0, -1, 0],
    [0, 0, 1]
])


# the four NV orientations
nNV = 1./np.sqrt(3)*np.array([
    [-1,1,1],
    [1,-1,1],
    [1,1,-1],
    [-1,-1,-1]
])


# corresponding NV angles in degree
#   theta   - angle between Nv direction and z axis
#   phi     - polar angle
nNV_angles = np.array([
    [54.735610317245, 135],
    [54.735610317245, -45],
    [125.26438968275, 45],
    [125.26438968275, -135]
])


def projetion_matrix(theta, phi):
    """
    projection matrix that projects a vector onto spherical coordinate systems that is defined by the angles th
    theta   - angle between Nv direction and z axis in (deg)
    phi     - polar angle (in deg)

    Returns: a 3x3 matrix

    """
    # convert to radians
    t = theta*np.pi / 180
    p = phi* np.pi / 180

    # spherical coordinate system
    nr = np.array([np.cos(p) * np.sin(t), np.sin(p) * np.sin(t), np.cos(t)])
    ntheta = np.array([np.cos(p) * np.cos(t), np.sin(p) * np.cos(t), -np.sin(t)])
    nphi = np.array([-np.sin(p), np.cos(p), 0])

    P = np.array([ntheta, nphi, nr])
    return P

def esr_frequencies(Bfield, gs=27.969, muB=1, hbar=1, Dgs=2.87):
    """
    :param Bfield (in Tesla): magnetic field with components Bx, By, Bz in the NV frame!: 1D-array of length 3 or 2D-array of dim Nx3
    :return:  matrix that gives the esr transition frequencies from diagonalizing the Hamiltonian with external magnetic field
        2 element array matrix if input B-field is 1D-array
        Nx2 array if input B-field is 2D-array of length Nx3
    """

    # if the input is a 1D array we cast it into a 2D array to work with the rest of the code
    if len(np.shape(Bfield))==1:
        assert len(Bfield)==3
        Bfield = [Bfield]
        input_1D = True
    elif len(np.shape(Bfield))==2:
        assert np.shape(Bfield)[1]==3
        input_1D = False



    def get_esr_freq(B):

        Hgs = hamiltonian_nv_spin1(B, gs=gs, muB=muB, hbar=hbar, D=Dgs)
        ev, Ugs = np.linalg.eigh(Hgs)

        return np.array([ev[1]-ev[0], ev[2]-ev[0]])

    esr = np.array([get_esr_freq(B) for B in Bfield])

    if input_1D:
        esr = esr[0]


    return esr*1e9

def esr_frequencies_ensemble(B_lab, gs=27.969, muB=1, hbar=1, Dgs=2.87):
    """
    calculates the esr freq. for the four NV families for a given magnetic field in the lab frame

    B_lab: magnetic field in the lab frame (N x 3) matrix

    returns the esr frequencies for all 4 NV families as M x N x 2 array, where
        M is the number of magnetic fields
        N = 4 is the number of NV families
    """

    f = []
    # get b field in cartesian coordinates

    for i in range(4):
        # get the off axis and on axis field for NV_i
        BNV = B_fields_in_NV_frame(B_lab,i)
        # calculate the ESR freq. for NV_i
        fo = esr_frequencies(BNV, gs=gs, muB=muB, hbar=hbar, Dgs=Dgs)

        f.append(fo)

    # return np.array(f)

    # rearrange so that we return a M x N x 2 array
    # M is the number of magnetic fields
    # N is the number of NV families
    return np.moveaxis(f, 0, 1)

def hamiltonian_nv_spin1(Bfield, gs=27.969, muB=1, hbar=1, D=2.87):
    """
    The hamiltonian for a spin 1 system, this hamiltonian describes the NV gound state as well as the excited state
    :param Bfield: magnetic field in Tesla with components Bx, By, Bz: 1D-array of length 3 or 2D-array of dim Nx3
    :param gs: gyromagnetic ration (per Tesla)
    :param muB:
    :param hbar:
    :param D: (2.87 for ground state, 1.42 for excited state) in GHz
    :return: the NV ground state hamiltonian
        3x3 matrix if input B-field is 1D-array
        Nx3x3 array if input B-field is 2D-array of length Nx3
    """

    # if the input is a 1D array we cast it into a 2D array to work with the rest of the code
    if len(np.shape(Bfield))==1:
        assert len(Bfield)==3
        Bfield = [Bfield]
        input_1D = True
    elif len(np.shape(Bfield))==2:
        assert np.shape(Bfield)[1]==3
        input_1D = False


    H = [hbar * D * Sz**2+gs*muB*(B[0]*Sx+B[1]*Sy+B[2]*Sz) for B in Bfield]

    if input_1D:
        H = H[0]

    return H

def transition_rate_matrix(Bfield, k12, k13, beta, kr = 63.2, k47= 10.8, k57 = 60.7, k71 = 0.8, k72 = 0.4):
    """
    the transition matrix

    :param Bfield: magnetic field with components Bx, By, Bz: 1D-array of length 3 or 2D-array of dim Nx3
    :param k12:
    :param k13:
    :param beta:
    :param kr:
    :param k47:
    :param k57:
    :param k71:
    :param k72:
    :return: transition rate matrix obtained from diagonalizing Hamiltonian
        7x7 matrix if input B-field is 1D-array
        Nx7x7 array if input B-field is 2D-array of length Nx3
    """

    # if the input is a 1D array we cast it into a 2D array to work with the rest of the code
    if len(np.shape(Bfield))==1:
        assert len(Bfield)==3
        Bfield = [Bfield]
        input_1D = True
    elif len(np.shape(Bfield))==2:
        assert np.shape(Bfield)[1]==3
        input_1D = False

    ko = np.matrix([
        [0, k12, k13, kr, 0, 0, k71],
        [k12, 0, 0, 0, kr, 0, k72],
        [k13, 0, 0, 0, 0, kr, k72],
        [kr*beta, 0, 0, 0, 0, 0, 0],
        [0, kr*beta, 0, 0, 0, 0, 0],
        [0, 0, kr*beta, 0, 0, 0, 0],
        [0, 0, 0, k47, k57, k57, 0]
    ])


    def get_k(B):
        U = np.array(coupling_matrix(B).H) # need to double check this here also double check if mathematica is correct!!

        k = np.zeros([7,7])

        for i in range(7):
            for j in range(7):
                k[i,j] = np.dot(np.dot(np.abs(U[i,:])**2,ko),np.abs(U[j,:])**2)

        return k

    k = [get_k(B) for B in Bfield]

    if input_1D:
        k = k[0]

    return k


def coupling_matrix(Bfield, gs=27.969, muB=1, hbar=1, Dgs=2.87, Des=1.42):
    """
    :param Bfield: magnetic field with components Bx, By, Bz: 1D-array of length 3 or 2D-array of dim Nx3
    :return: unitary matrix that diagonalizes the Hamiltonian with external magnetic field
        7x7 matrix if input B-field is 1D-array
        Nx7x7 array if input B-field is 2D-array of length Nx3
    """

    # if the input is a 1D array we cast it into a 2D array to work with the rest of the code
    if len(np.shape(Bfield))==1:
        assert len(Bfield)==3
        Bfield = [Bfield]
        input_1D = True
    elif len(np.shape(Bfield))==2:
        assert np.shape(Bfield)[1]==3
        input_1D = False



    def get_Uo(B):
        Uo = np.matrix(np.zeros([7,7])+0j)

        Hgs = hamiltonian_nv_spin1(B, gs=gs, muB=muB, hbar=hbar, D=Dgs)
        ev, Ugs = np.linalg.eigh(Hgs)

        Hes = hamiltonian_nv_spin1(B, gs=gs, muB=muB, hbar=hbar, D=Des)
        ev, Ues = np.linalg.eigh(Hes)

        Uo[0:3,0:3] = Ugs

        Uo[3:6, 3:6] = Ues

        Uo[6, 6] = 1
        return Uo

    Uo = [get_Uo(B) for B in Bfield]

    if input_1D:
        Uo = Uo[0]


    return Uo

def populations(transition_rates):
    """
    calculates the population by solving the rate equations for the transition rate matrix k
    :param transition_rates:
        2D array with dim MxM
        3D array with dim NxMxM
    :return:
        population of the M levels
        1D array of length M if input is 2D array with dim MxM
        2D array of length NxM if input is 3D array with dim NxMxM
    """

    # if the input is a 2D array we cast it into a 3D array to work with the rest of the code
    if len(np.shape(transition_rates))==2:
        transition_rates = [transition_rates]
        input_2D = True
    elif len(np.shape(transition_rates))==3:
        input_2D = False

    def get_pop(k):
        k = np.matrix(k)
        a = k - np.diag(np.array(np.sum(k, 0))[0])

        a = np.row_stack([a, np.ones([1, len(k)])])
        b = np.hstack([np.zeros(len(k)), [1]])

        n, residuals, rank, s = np.linalg.lstsq(a,b)

        # some extra information from the fit we could use to validate the result
        # print('residuals', residuals)
        # print('rank', rank)
        # print('s', s)
        # print('n', n)

        return n

    n = [get_pop(k) for k in transition_rates]



    if input_2D:
        n = n[0]


    return n

def get_ko(k12, k13, beta, kr = 63.2, k47= 10.8, k57 = 60.7, k71 = 0.8, k72 = 0.4):
    """
    the transition matrix

    :param B: the magnetic field
    :param k12:
    :param k13:
    :param beta:
    :param kr:
    :param k47:
    :param k57:
    :param k71:
    :param k72:
    :return:
    """


    ko = np.matrix([
        [0, k12, k13, kr, 0, 0, k71],
        [k12, 0, 0, 0, kr, 0, k72],
        [k13, 0, 0, 0, 0, kr, k72],
        [kr*beta, 0, 0, 0, 0, 0, 0],
        [0, kr*beta, 0, 0, 0, 0, 0],
        [0, 0, kr*beta, 0, 0, 0, 0],
        [0, 0, 0, k47, k57, k57, 0]
    ])

    return ko

def photoluminescence_rate(transition_rates, populations):
    """

    :param transition_rates: transition rate matrix obtained from diagonalizing Hamiltonian
        7x7 matrix
        Nx7x7 array
    :param populations: population of the M levels
        1D array of length 7
        2D array of length Nx7
    :return: photoluminescence rate
        scalar if populations is a 1D array
        1D-array if populations is a 2D array
    """


    # if the input is a 1D array we cast it into a 2D array to work with the rest of the code
    if len(np.shape(transition_rates))==2:
        assert np.shape(transition_rates)==np.array([7,7])
        assert np.shape(populations) == np.array([7])
        transition_rates = [transition_rates]
        populations = [populations]
        input_1D = True
    elif len(np.shape(transition_rates))==3:
        assert np.shape(transition_rates)[1:3]==(7,7)
        assert np.shape(populations)[1] == 7
        input_1D = False

    r = [np.sum(np.dot(k[3:6, 0:3], pop[3:6])) for k, pop in zip(transition_rates, populations)]
    # r = np.sum(np.dot(k[3:6,0:3], pop[3:6]))
    if input_1D:
        r = r[0]

    return r


def photoluminescence_contrast(Bfield, k12, k13, beta, kr=63.2, k47=10.8, k57=60.7, k71=0.8, k72=0.4):
    """

    :param Bfield: magnetic field with components Bx, By, Bz: 1D-array of length 3 or 2D-array of dim Nx3
    :param k12:
    :param k13:
    :param beta:
    :return: photoluminescence contrast in percent
    """


    k_no_mw = transition_rate_matrix(Bfield, 0, 0, beta, kr=kr, k47=k47, k57=k57, k71=k71, k72=k72)
    k_mw = transition_rate_matrix(Bfield, k12, k13, beta, kr=kr, k47=k47, k57=k57, k71=k71, k72=k72)

    pop_no_mw = populations(k_no_mw)
    pop_mw = populations(k_mw)

    pl_no_mw = photoluminescence_rate(k_no_mw, pop_no_mw)
    pl_mw = photoluminescence_rate(k_mw, pop_mw)


    if len(np.shape(pl_no_mw))==1:
        pl_no_mw = np.array(pl_no_mw)
        pl_mw = np.array(pl_mw)
        c = np.array(pl_no_mw - pl_mw)/np.array(pl_no_mw) * 100.
    else:
        c = np.array(pl_no_mw - pl_mw) / np.array(pl_no_mw) * 100.

    return c


def B_field_from_esr(fp, fn, D=2.8707e9, gamma=27.969e9, angular_freq=False, verbose=False):
    """
    wp, wn: upper and lower esr frequency
    gamma: Gyromagnetic ratio (in GHz/Tesla)
    D: NV zero field splitting (in GHz)
    angular_freq: frequencies are the angular frequencies (default = False)
    returns:
        the field along and perpendicular to the NV axis
    """

    if angular_freq:
        print('WARNING CHECK CODE TO MAKE SURE THAT ALL FREQ. ARE ANGULAR FREQs')

    wp = fp
    wn = fn

    # check that wp>wn if not flip them
    if wp < wn:
        wp, wn = wn, wp

    Bz = np.sqrt(-(D + wp - 2 * wn) * (D + wn - 2 * wp) * (D + wn + wp)) / (3 * gamma * np.sqrt(3 * D))
    Bp = np.sqrt(-(2 * D - wp - wn) * (2 * D + 2 * wn - wp) * (2 * D - wn + 2 * wp)) / (3 * gamma * np.sqrt(3 * D))

    if np.isnan(Bp):
        Bp = 0
        if verbose:
            print(('is nan Bp', Bp))
            print(('fp/fn', fp, fn))

    if np.isnan(Bz):
        Bz = 0
        if verbose:
            print(('is nan Bz', Bz))
            print(('fp/fn', fp, fn))

    return Bz, Bp


def B_field_from_esr_ensemble(frequencies, angular_freq=False):
    """
    calculates the magnetitc field components from the NV ESR frequencies

    frequencies (n x 2 matrix or vector of length 2xn ): upper and lower esr frequency for n peaks
    angular_freq: frequencies are the angular frequencies (default = False)
    returns:
        the three vector components Bx, By, Bz in the frame of the Diamond
    """

    print('WARNING THIS FUNCTION MIGHT BE OBSOLETE TRY calc_bfields_esr_ensemble_mag')
    if angular_freq:
        print('WARNING CHECK CODE TO MAKE SURE THAT ALL FREQ. ARE ANGULAR FREQs')

    if len(np.shape(frequencies)) == 1:
        frequencies = np.reshape(frequencies, [len(frequencies) / 2, 2])

    assert len(frequencies.T) == 2

    B = []
    for f in frequencies:
        B.append(B_field_from_esr(*f))

    return B


def B_fields_in_NV_frame(B_lab, NV_id):
    """
    returns the magnetic field in the frame of the NV center

    B_lab: magnetic field in lab frame (Bx, By, Bz) components matrix with nx3 elements or a vector of length 3
    NV_id: id of NV (number from 0 to 3)


    returns:
        magnetic field in frame of NV center [Bx, By, Bz]_NV as a nx3 matrix or a vector of length 3
    """
    vector = False  # if B_lab is a vector or a matrix
    # if b is a vector convert to matrix
    if len(np.shape(B_lab)) == 1:
        B_lab = np.array([B_lab])
        vector = True

    # assert that the form of B is as expected
    assert len(np.shape(B_lab)) == 2
    assert np.shape(B_lab)[1] == 3

    B_NV = np.dot(projetion_matrix(*nNV_angles[NV_id]), B_lab.T).T

    if vector:
        B_NV = B_NV[0]

    return B_NV


def B_cart(B_mag, theta, phi):
    """
    calculate the cartesian components of the magenetic field

    B_mag(vector): magnitude of magnetic field
    phi: angle of field vector in x-y plane (in deg)
    theta: angle between the z-axis and xy-plane (in deg)

    return: cartesian components of the magenetic field
    """
    p = phi*np.pi/180
    t = theta*np.pi / 180
    Bx = B_mag * np.cos(p) * np.sin(t)
    By = B_mag * np.sin(p) * np.sin(t)
    Bz = B_mag * np.cos(t)

    return np.array([Bx, By, Bz]).T


def B_spher(Bx, By, Bz):
    """
    takes cartesian coordinates and returns spherical coordiantes
    B_mag, theta, phi (in deg)
    """

    B_mag = np.sqrt(Bx ** 2 + By ** 2 + Bz ** 2)
    B_rho = np.sqrt(Bx ** 2 + By ** 2)
    theta = np.arctan2(B_rho, Bz) * 180 / np.pi
    phi = np.arctan2(By, Bx) * 180 / np.pi

    return B_mag, theta, phi

def calc_bfields_esr_ensemble_mag(frequencies, verbose=False):
    """
    Calculates the magnetic fields from the ESR frequencies of N NV families along the NV axis and perpedicular

    frequencies: the esr frequencies, e.g. the frequencies obtained from a fit to ESR data. (in GHz)
    This is a N x 2 matrix, where N is the number of families and 2 are the two frequencies

    of

    This is a vector of length 2*N, where the ordering is NVa_low, NVa_high, NVb_low, NVb_high, etc

    returns:
        Babs:   the absolute value of the magnetic field
        Bs:     field along the NV axis and perpedicular
    """

    if len(np.shape(frequencies)) == 1:
        # reshape to expected N x 2 format
        frequencies = np.reshape(frequencies, (len(frequencies)/2, 2))

    assert len(frequencies.T) == 2

    if verbose:
        print(' ===== calc_bfields_esr_ensemble mag ==== ')

    number_of_families = len(frequencies)

    # calculate the on axis and off axis field for each family
    Bs = np.array([B_field_from_esr(f[0], f[1]) for f in frequencies])

    # calculate the abolute field
    Babs = np.sqrt(np.sum(Bs ** 2, axis=1))
    if verbose:
        print(('consistnecy check: total field should be the same for all families - std_dev', np.std(Babs)/np.mean(Babs)))
        if np.std(Babs)/np.mean(Babs)<1e-4:
            print('PASSED!!!')
        else:
            print('FAILED!!!')
    Babs = np.mean(Babs)

    return Babs, Bs

def calc_bfields_esr_ensemble_xyz(frequencies, verbose=False):
    """
    Calculates the magnetic fields from the ESR frequencies of N NV families

    frequencies: the esr frequencies, e.g. the frequencies obtained from a fit to ESR data. (in GHz)
    This is a N x 2 matrix, where N is the number of families and 2 are the two frequencies

    """

    def find_similar(X, Babs, eps=1e-4):
        """
        we assume that there are two values for each of the N x rows in x that are very similar and we want to get that value

        X: matrix of dimensions N x M where
            N is 3 for the three cartesian field components and
            M = 4 is the number of methods used to get the cartesian field components
        Babs: abolute field obtained from the avrg over all families, use this as a consistency check
        eps: threshold at which two values are considered identical

        """

        def pick_correct_B_field(index):
            """

            picks the correct field from the global set X and returns it
            Args:
                index: index from truth table

            Returns:

            """
            # calculate the actual index from 0 to M-1, that corresponds to the correct value for the cartesian field
            index = np.arange(M ** 2).reshape(M, M)[iu1][index] % M
            return X[:, index]

        N, M = np.shape(X)

        # indecies to pick upper triangle of all - all matrix
        iu1 = np.triu_indices(M, 1)

        # Truth table: here we calculate the difference between all possible combinations and check if it's smaller than eps
        # the result is a N x (M-1)/2 matrix
        truth_table = np.array([np.array([c - np.abs(x) for c in np.abs(x)])[iu1] < eps for x in X])

        if verbose:
            print('     truth_table:')
            for i, t in enumerate(truth_table.T):
                print((i, t))

            # consistency check: we expect that there is one combination, that work for all and thus
            check = max(np.sum(truth_table, axis=0)) == M - 1
            if check:
                print('consistency checked successfully')
            else:
                print('consistency check failed')

        # majority vote to get the index where most agree
        votes = np.sum(truth_table, axis=0)
        index = np.arange(len(votes))[votes == np.max(votes)]
        if verbose:
            print(('index majority vote: ', index, np.sum(truth_table, axis=0)))
        print(('asdad', len(index)))

        # compare the total magnetic field to the expected magnetic field, pick the one that agrees better

        # first calculate the error:
        err = [np.abs(Babs - np.sqrt(np.sum(pick_correct_B_field(i) ** 2))) for i in index]

        # now take the one with the smaller error
        B_cart = pick_correct_B_field(index[np.argmin(err)])


        # if len(index) == 1:
        #     if verbose:
        #         print('one set of fields found !')
        #     B_cart = pick_correct_B_field(index[0])
        #
        # if len(index) >1:
        #     if verbose:
        #         print('more than one set of fields found, check for amplitude consistency')
        #
        #     # compare the total magnetic field to the expected magnetic field, pick the one that agrees better
        #
        #     # first calculate the error:
        #     err = [np.abs(Babs - np.sqrt(np.sum(pick_correct_B_field(i) ** 2))) for i in index]
        #
        #     # now take the one with the smaller error
        #     B_cart = pick_correct_B_field(index[np.argmin(err)])
        #
        #     if verbose:
        #         print('devitation from expected field: ', np.min(err))

        return B_cart

    number_of_families = len(frequencies)

    Babs, Bs = calc_bfields_esr_ensemble_mag(frequencies, verbose)

    Bxyz = None
    # === get Bx, By, Bz ====
    if number_of_families == 4:

        if verbose:
            print(' Bs (on and off axis field) ==== ')
            print(Bs)

            # print(' |theta| ==== ')
            # print(np.arctan2(Bs[:, 0], Bs[:, 1]) * 180 / np.pi)
            # print(np.arctan2(Bs[:, 1], Bs[:, 0]) * 180 / np.pi)
        # calculate different combinations of the on-axis field,
        # from which we know that they should give us the cartesian components
        Bxyz = []
        # combination 1 to get xyz
        Bxyz.append(np.sqrt(3) / 2 * np.array([Bs[0, 0] - Bs[3, 0], Bs[0, 0] + Bs[2, 0], Bs[0, 0] + Bs[1, 0]]))
        # since we don't know the sign combination 1 (sign inverted) to get xyz
        Bxyz.append(np.sqrt(3) / 2 * np.array([Bs[0, 0] + Bs[3, 0], Bs[0, 0] - Bs[2, 0], Bs[0, 0] - Bs[1, 0]]))
        # combination 2 to get xyz
        Bxyz.append(np.sqrt(3) / 2 * np.array([Bs[1, 0] + Bs[2, 0], Bs[1, 0] - Bs[3, 0], Bs[2, 0] - Bs[3, 0]]))
        # since we don't know the sign combination 2 (sign inverted) to get xyz
        Bxyz.append(np.sqrt(3) / 2 * np.array([Bs[1, 0] - Bs[2, 0], Bs[1, 0] + Bs[3, 0], Bs[2, 0] + Bs[3, 0]]))
        Bxyz = np.array(Bxyz).T
        if verbose:
            print(' Bxyz not selected ==== ')
            print((Bxyz.T))

        # at this point Bxyz is a 4x3 array (4 combinations, 3 cartesian components)
        # from all the combinations, we take the one that occurs twice
        Bxyz = find_similar(Bxyz, Babs)
        if verbose:
            print(' Bxyz selected ==== ')
            print((Bxyz.T))

            # rearrange so that we return a M x N x 2 array
            # M is the number of magnetic fields
            # N is the number of NV families
            #     return np.moveaxis(Bxyz, 0, 1)

    if verbose:
        print(' ===== END calc_bfields_esr_ensemble ==== ')

    return Bxyz


def fit_Hamiltonian(freq, verbose=False, try_permutations_fit = True, try_permutations_xyz = True, try_permutations_sign = True):
    """
    retrieve the field amplitude and angles by fitting to the ensemble Hamiltonian

    freq (2 x n matrix or 2xn vector): upper and lower esr frequency for n peaks

    try_permutations_fit: permutate xyz components and optimize for each to reduce error
    try_permutations_xyz: permutate xyz components and pick the one with the lowest reduce error
    try_permutation_sign: checks for all possible combinations of signs and  picks the one with the lowest reduce error
    returns:
        B (in Teslas), theta (in deg), phi (in deg)
    """

    def fit_Hamiltonian_inital_guess(freq, verbose=False):
        """
        Here we guess the initial condition for fitting to the NV Hamiltonian

        """
        # get on and off axis field in NV frames
        Br = np.array([B_field_from_esr(*f) for f in freq.T])

        # estimate the total amplitude
        Br_mag = np.diag(np.sqrt(np.dot(Br, Br.T)))

        # all families should be give similar total fields
        if np.std(Br_mag) / np.mean(Br_mag) > 1e-2:
            print('Warning the total magnetic field estimated from all families differs by more than 1%!')
            print(('relative error', np.std(Br_mag) / np.mean(Br_mag)))
            print(('absolute error (Teslas)', np.std(Br_mag)))

        if verbose:
            print(('estimated field amplitude', np.mean(Br_mag), np.std(Br_mag)))

        NV_max_index = np.argmax(Br[:, 0])
        if verbose:
            print('find the NV with the largest on axis field:')
            print(NV_max_index)

        # as a starting point we use the angle of this NV center as an initial guess
        theta_init, phi_init = nNV_angles[NV_max_index]
        B_mag_init = np.mean(Br_mag)

        if verbose:
            print(('initial guess:', B_mag_init, theta_init, phi_init))
            print(('err', fit_err_fun([theta_init, phi_init], B_mag_init, freq)))

        return B_mag_init, theta_init, phi_init

    B_mag_init, theta_init, phi_init = fit_Hamiltonian_inital_guess(freq, verbose)

    fit = opt.minimize(fit_err_fun, np.array([theta_init, phi_init]), args=(B_mag_init, freq),
                       bounds=((0, 90), (0, 180)))


    theta_r, phi_r = fit.x
    err = fit.fun

    if verbose:
        print('fit result')
        print((' theta_r, phi_r', theta_r, phi_r))


    if try_permutations_fit:
        if verbose:
            print('======== trying permutations')
        for B_it in list(itertools.permutations(B_cart(B_mag_init, theta_r, phi_r))):
            B_it_s = B_spher(*B_it)

            fit = opt.minimize(fit_err_fun, np.array([B_it_s[1], B_it_s[2]]), args=(B_it_s[0], freq),
                           bounds=((0, 180), (-180, 180)))
            if verbose:
                print(('B_it_s', B_it_s))
                print(('fit.x', fit.x))
                print(('fit.fun', fit.fun))
                print((' theta_r, phi_r', theta_r, phi_r))
                print('---------------------------------------')
            # if error reduced keep new value
            if fit.fun < err:
                if verbose:
                    print('found better result (perm)')
                theta_r, phi_r = fit.x
                err = fit.fun
        if verbose:
            print('======== fit result after permutation_fit')
            print(('======== theta_r, phi_r', theta_r, phi_r))
    if try_permutations_xyz:
        if verbose:
            print('======== trying permutating xyz')
        for B_it in list(itertools.permutations(B_cart(B_mag_init, theta_r, phi_r))):
            B_it_s = B_spher(*B_it)
            err_xyz = fit_err_fun([B_it_s[1], B_it_s[2]], B_it_s[0], freq)
            if err_xyz < err:
                if verbose:
                    print('found better result (xyz)')
                theta_r, phi_r = B_it_s[1], B_it_s[2]
                err = err_xyz
        if verbose:
            print('======== fit result after permutation_xyz')
            print(('======== theta_r, phi_r', theta_r, phi_r))
    if try_permutations_sign:
        if verbose:
            print('======== trying permutating sign')
        # for B_it in [-1, 1]:
        #     B_it_s = B_spher(*B_it)
            err_sign = fit_err_fun([B_it_s[1], 180+B_it_s[2]], B_it_s[0], freq)
            if err_sign < err:
                if verbose:
                    print('found better result (err_sign)')
                theta_r, phi_r = B_it_s[1], B_it_s[2]+180
                err = err_sign
        if verbose:
            print('======== fit result after permutation sign')
            print(('======== theta_r, phi_r', theta_r, phi_r))
    return B_mag_init, theta_r, phi_r


def fit_err_fun(x, *argv):
    """

    theta: angle theta (in deg)
    phi: angle phi (in deg)

    argv[0]: B_mag: magnetic field amplitude
    argv[1]: freq: measured frequencies


    call:
    fit_err_fun([theta_init, phi_init], B_mag_init, freq)

    """

    theta, phi = x
    #
    B_mag = argv[0]
    freq = argv[1]

    B_lab = B_cart(B_mag, theta, phi)
    freq_est = esr_frequencies_ensemble(B_lab)

    err = np.mean(np.abs((freq_est - freq) / freq))

    return err


def get_r_dr(nv_locations, magnet_diam, magnet_center):
    """

    use for line data, nv_locations are along a line

    calculates the distance to the magnet and the distance between points
    nv_locations: matrix of shape N x 2, where N is the number of NV positions
    magnet_center: distance from first point in nv_locations to magnet center (must be same units as nv_location!!)
    magnet_diam: magnet diameter (must be same units as nv_location!!)
    """
    r = np.array([np.sqrt(np.sum((pt - nv_locations[0]) ** 2)) for pt in nv_locations])
    r = r + magnet_center
    # take into account that we are at least a radius above the magnet
    r = np.sqrt(r**2+magnet_diam/2)
    dr = np.mean(np.diff(r))

    return r, dr




def get_theta_dr(nv_locations, method='radius'):
    """

    use for ring data, i.e. nv_locations are on a ring

    calculates the angles theta and the distance between points
    nv_locations: matrix of shape N x 2, where N is the number of NV positions
    method: calculate dr based on differences (diff) or based on the radius (radius)
        string 'diff' or 'radius'


    :returns
        theta = angles on ring
        radius = radius of ring
        dr = distance between points on ring
    """
    mag_pos = np.mean(nv_locations, 0)
    # radius of circle (in um)
    radius = np.mean(np.sqrt(np.sum((nv_locations - mag_pos) ** 2, 1)))

    # we measured on a ring, so we calculate the angle for each measurement
    theta = np.linspace(0, 360, len(nv_locations))

    if method == 'diff':
        # calculate differential distance between points
        x = np.sum(np.diff(nv_locations - mag_pos, 0) ** 2, 1)
        np.mean(x), np.std(x)
        dr = np.mean(x)
    #         print('differential from differences:', dr)
    elif method == 'radius':
        # differential theta (spacing betweent thetas)
        dtheta = np.mean(np.diff(theta)) * np.pi / 180
        # calculate the distance between measurements, needed for calculation f gradients
        dr = radius * dtheta
    #         print('differential from radius:', dr)
    else:
        print("unknown method try 'diff' or 'radius'")

    return theta, dr


def sort_esr_frequencies(freq_data, permutate_all = True, verbose = False):
    """
    sorts the frequencies from the measurement by trying all different permutations and minimizing the error in the total field,
    while also maximizing the frequency overlap from one measurement to the next.
    freq_data: frequency data, array with dimensions Ndata (number of datasets), Nfreq (number of frequencies)
    verbose: print information along the way

    todo: implement also maximum overlap of the slope to keep track of upper and lower NV frequency

    permutate_all:
        TRUE: permutate all combinations of frequencies (7! = 5040 for 8 frequencies )
        FALSE: permutate only half of combinations of frequencies (4! = 24 for 4 frequencies)

    :returns
        freqs_sorted: sorted frequencies
        perm_index: index of permutation that has to be applied to freq_data to get freqs_sorted

    """
    freq_data = np.array(freq_data)

    if len(np.shape(freq_data)) == 2:
        Ndata, Nfreq = np.shape(freq_data)
    perm_index = []
    freqs_sorted = []

    def calc_err(freq):
        """
        calculates the error in B field for a give collection of frequencies
        """
        # calculate the fields for each family
        Bs = calc_bfields_esr_ensemble_mag(freq)[1]
        # calculate the total field
        Babs = np.sqrt(np.sum(Bs ** 2, axis=1))
        # calculate the error
        # err = np.std(Babs) / np.mean(Babs)
        err = np.std(Babs)

        return err

    def get_perm_freq(freq, index):
        """
        returns the permutation of frequencies, that corresponds to permutation index "index"
        note that we only permutate the last half of the array, since we are looking for NV pairs

        freq = vector of frequencies (typically length = 8 for all four NV families)
        index = single integer or a list of integers
        """
        Nfreq = len(freq)

        if permutate_all:
            # # for a single index return the list
            # if len(np.shape(index)) == 0:
            #     return list(permutations(freq))[index]
            # # for a list of indecies loop over all indecies and return a list for each
            # else:
            #     return [list(permutations(freq))[i] for i in index]
            # for a single index return the list
            if len(np.shape(index)) == 0:
                return list([freq[0]]) + list(list(permutations(freq[1:]))[index])
            # for a list of indecies loop over all indecies and return a list for each
            else:
                return [list([freq[0]]) + list(list(permutations(freq[1:]))[i]) for i in index]
        else:
            # for a single index return the list
            if len(np.shape(index)) == 0:
                return list(freq[0:Nfreq / 2]) + list(list(permutations(freq[Nfreq / 2:]))[index])
            # for a list of indecies loop over all indecies and return a list for each
            else:
                return [list(freq[0:Nfreq / 2]) + list(list(permutations(freq[Nfreq / 2:]))[i]) for i in index]

    def calc_err_freq(freq0, freqs_perm):
        """
        calculates the error in the frequencies freq_perm (e.g. N x 8) matrix
        and the frequencies freq_0 (e.g. vector of length 8)
        """
        Nperm = len(freqs_perm)
        err = np.sum((np.array(freqs_perm) - np.ones([Nperm, 1]) * np.array([freq0])) ** 2, 1)
        # normalize err
        err = err / np.sum(np.array(freq0)) ** 2

        if verbose:
            print(' ====== calc_err_freq: err =====')
            print((np.shape(err)))

            # print(np.shape(freqs_perm), np.shape(np.ones([Nperm, 1]) * np.array([freq0])))

            print(['{:0.3e}'.format(err[k] / 1e9) for k in range(len(err))])

        return err


    for j, freq in enumerate(freq_data):
        if verbose:
            print(('>>>>>>>>>>>>> RUN <<<<<<<<<<<<<<<<', j))
        # permutate over all four families to find the match that gives the lowest error
        if permutate_all:
            # errs = [calc_err(np.array(freq_perm))
            #         for freq_perm in list(permutations(freq))]
            # concat the first 1 freq and the 7 permutated freqs
            errs = [calc_err(np.array([freq[0]] + list(freq_perm)))
                    for freq_perm in list(permutations(freq[1:]))]
        else:
            # concat the first 4 freq and the 4 permutated freqs
            errs = [calc_err(np.array(list(freq[0:Nfreq / 2]) + list(freq_perm)))
                    for freq_perm in list(permutations(freq[Nfreq / 2:]))]

        perm_indecies_min = np.where(errs == min(errs))[0]  # permutation indecies that minimize the error

        if verbose:
            print(' ====== perm_indecies_min =====')
            print(perm_indecies_min)

        #         print(perm_indecies_min)
        if len(perm_index) == 0:
            # there mightbe several permutations with the same error, we take the first for the first set of freqs
            freq_index_min = 0
        else:
            # if one of the indeces is the same as in the previous dataset,
            # we take the one that minimizes the change in frequency
            # frequencies corresponding to the permutations with the lowest error
            freqs_perm_min = get_perm_freq(freq, perm_indecies_min)
            # out of those permutations find the one that maximizes the freq. overlap
            freq_index_min = np.argmin(calc_err_freq(freqs_sorted[-1], freqs_perm_min))


        perm_index.append(perm_indecies_min[freq_index_min])
        freqs_sorted.append(get_perm_freq(freq, perm_index[-1]))
        if verbose:
            print(' ====== original frequencies =====')
            print(['{:0.3f}'.format(freq[k] / 1e9) for k in range(len(freq))])
            # print(freq)
            print(' ====== permutated frequencies =====')
            # print(freqs_sorted[-1])
            print(['{:0.3f}'.format(freqs_sorted[-1][k] / 1e9) for k in range(len(freqs_sorted[-1]))])

    return freqs_sorted, perm_index

def connect_esr_frequencies(esr_data, verbose=False):
    """
    order the esr_data such that the frequency overlap from one measurement to the next is maximized
    (assuming that the data is continuous and assuming that the esr data is already in the right pairs, i.e.
    the rows are of the form NV1_low, NV1_high, NV2_low, NV2_high etc.

    run this after sorting with sort_esr_frequencies to ensure continuity of data

    returns: the ordered frequency freq
    :param esr_data: esr data of the form M x 2*N, where M is the number of measurements and N is the number of NV families

    :param verbose: if true output more text
    :return: continuous esr data
    """
    def order_frequencies(freq, freq_last):
        """
        order the freq such that the overlap with the previous set is maximized
        Here we assume that the NVs are of the form
            N x 2, where N is the number of NV families
        or
            vector of length 2*N, where frequencies from the same family are next to each other


        freq: frequencies to be ordered
        freq_last: last frequencies that are used as reference

        returns: the ordered frequency freq

        """

        assert np.shape(freq) == np.shape(freq_last)
        if verbose:
            print(['{:0.3f}'.format(freq[k] / 1e9) for k in range(len(freq))])
            print(['{:0.3f}'.format(freq_last[k] / 1e9) for k in range(len(freq_last))])

        # if vector get into the desired form
        if len(np.shape(freq)) == 1:
            freq = np.reshape(freq, [len(freq) / 2, 2])
            freq_last = np.reshape(freq_last, [len(freq_last) / 2, 2])

        freq = np.sort(freq)
        freq_last = np.sort(freq_last)
        # find permutation family that minimizes frequency jumps
        index_min = np.argmin([np.sum((freq_last - np.array(x)) ** 2) for x in list(permutations(freq))])

        return np.reshape(list(permutations(freq))[index_min], len(freq) * 2)

    # sort the first dataset such that
    freq_last = esr_data[0]
    freq_last = np.reshape(freq_last, [len(freq_last) / 2, 2])
    freq_last = np.sort(freq_last)
    freq_last = list(np.reshape(freq_last, len(freq_last) * 2))

    esr_data_sorted = []
    #     freq_last = esr_data[0]
    for i in range(len(esr_data) - 1):
        freq = esr_data[i + 1]
        if verbose:
            print(('==================', i))
        esr_data_sorted.append(order_frequencies(freq, freq_last))

        freq_last = esr_data_sorted[-1]

    return np.array(esr_data_sorted)


def magnetic_moment_and_Br_from_fit(dp, a, r, mu0=4 * np.pi * 1e-7):
    """
    calculate the magentic moment and magnetic surface field from the fit parameter dp
    a: radius of magnet
    r: distance between NV circle and center of magnet
    """
    V = 4 * np.pi / 3 * a ** 3
    m = 4 * np.pi / mu0 * r ** 3 * dp
    Br = m / V * mu0
    return m, Br



if __name__ == '__main__':

    # solution calculated previously
    ref_contrast = [11.259005142154569, 11.162777024524557, 10.881956479400158, 10.43856472476315, 9.8649445385505405]
    ref_Bmag = [ 0, 10, 20, 30, 40]


    h_ref=[np.matrix([[ 0.00+0.j,  0.00+0.j,  0.00+0.j],
        [ 0.00+0.j,  2.87+0.j,  0.00+0.j],
        [ 0.00+0.j,  0.00+0.j,  2.87+0.j]]), np.matrix([[ 0.00000000+0.j        ,  0.00000000-0.01901095j,
          0.00000000-0.01901095j],
        [ 0.00000000+0.01901095j,  2.86229071+0.j        ,  0.00000000+0.j        ],
        [ 0.00000000+0.01901095j,  0.00000000+0.j        ,  2.87770929+0.j        ]]), np.matrix([[ 0.00000000+0.j        ,  0.00000000-0.03802189j,
          0.00000000-0.03802189j],
        [ 0.00000000+0.03802189j,  2.85458142+0.j        ,  0.00000000+0.j        ],
        [ 0.00000000+0.03802189j,  0.00000000+0.j        ,  2.88541858+0.j        ]]), np.matrix([[ 0.00000000+0.j        ,  0.00000000-0.05703284j,
          0.00000000-0.05703284j],
        [ 0.00000000+0.05703284j,  2.84687213+0.j        ,  0.00000000+0.j        ],
        [ 0.00000000+0.05703284j,  0.00000000+0.j        ,  2.89312787+0.j        ]]), np.matrix([[ 0.00000000+0.j        ,  0.00000000-0.07604378j,
          0.00000000-0.07604378j],
        [ 0.00000000+0.07604378j,  2.83916283+0.j        ,  0.00000000+0.j        ],
        [ 0.00000000+0.07604378j,  0.00000000+0.j        ,  2.90083717+0.j        ]])]

    k12 = 1
    k13 = 0
    beta = 0.3

    B= np.array([[0.0000961262, 0., 0.0000275637]]).T
    Bmag = np.arange(0,50,10)
    B = np.dot(B, np.array([Bmag])).T



    # h = [hamiltonian_nv_spin1(b) for b in B]
    # h = coupling_matrix(B)
    c = photoluminescence_contrast(B, k12, k13, beta)
    # k = transition_rate_matrix(B, k12, k13, beta)
    #
    # p = populations(k)
    #
    # pl = photoluminescence_rate(k, p)

    print((np.shape(c)))
    print(c)

    print((np.allclose(c, ref_contrast)))
