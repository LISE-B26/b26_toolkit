# Jan Gieseler
# December 12 2017
# here we test the functionc from nv_optical_response.py

from . import nv_optical_response as nv
import numpy as np


def test_magnetic_field_single_mag(Bmax=0.1, verbose=True):
    """
    here we test how well we can reconstruct the on axis and off axis field from the ESR freq. from a single NV ESR

    Bmax: maximum random field in Teslas
    """
    import random

    # get random magnetic field
    B = random.random() * Bmax
    # get random angle field
    theta = random.random() *180
    phi = random.random() *180

    error = calc_error_single(B, theta, phi, verbose)

    if verbose:
        print(('error', error))
        print(('B, theta, phi', B, theta , phi))

    return error, B, theta, phi


def calc_error_single(B, theta, phi, verbose=False):
    """
    here we test how well we can reconstruct the on axis and off axis field from the ESR freq. from a single NV ESR

    B, theta, phi: amplitude and angles of magnetic field
    """

    Bx, By, Bz= nv.B_cart(B, theta, phi)
    # Bx, By, Bz = B * np.sin(theta) * np.cos(phi), B * np.sin(theta) * np.sin(phi), B * np.cos(theta)
    fn, fp = nv.esr_frequencies([Bx, By, Bz])

    Bon_r, Boff_r = nv.B_field_from_esr(fp, fn)

    # calculate the error in the off axis field
    error_on = np.abs(Boff_r - np.sqrt(Bx ** 2 + By ** 2)) / B
    # calculate the error in the on axis field
    error_off = np.abs(Bon_r - np.abs(Bz)) / B

    # calculate the total error
    error = np.sqrt(error_on ** 2 + error_off ** 2)
    if verbose:

        if error > 9e-1:
            print(('Bfield (T) / theta (deg) / phi (deg)', B, theta, phi))
            print(('Bxyz', Bx, By, Bz))
            print(('Bon, Boff', Bz, np.sqrt(Bx ** 2 + By ** 2)))
            print(('Bon_r, Boff_r', Bon_r, Boff_r))
            #         print('NV frequencies (GHz)', fn, fp)

            #         print('Bfield (T) x/z', Bx, Bz)
            #         print('reconstructed Bfield (T) x/z', Bx_r, Bz_r)
            print(('error_on (%)', error_on * 100))
            print(('error_off (%)', error_off * 100))
            print(('error (%)', error * 100))
    return error


def test_magnetic_field_mag_ensemble(Bmax=0.1, verbose=True):
    """
    here we test how well we can reconstruct the on axis and off axis field from the ESR freq. from a single NV ESR

    Bmax: maximum random field in Teslas
    """
    import random

    # get random magnetic field
    B = random.random() * Bmax
    # get random angle field
    theta = random.random() *180
    phi = random.random() *180

    error = calc_error_ensemble_mag(B, theta, phi, verbose)

    if verbose:
        print(('error', error))
        print(('B, theta, phi', B, theta , phi))

    return error, B, theta, phi
def calc_error_ensemble_mag(B, theta, phi, verbose=False):
    """
    here we test how well we can reconstruct the total magnetic field amplitude
    from the ESR freq. from a N families of NV ESR

    B, theta, phi: amplitude and angles (in deg )of magnetic field
    """

    if verbose:
        print(('B, theta, phi', B, theta, phi))

    # calculate field in lab frame
    B_lab = nv.B_cart(B, theta, phi)

    # if B_lab is a vector, the we turn it into a matrix because that is what we expect in the following lines
    if len(np.shape(B_lab)) == 1:
        #         print('reshaping')
        B_lab = np.array([B_lab])
        # print('B_lab shape', np.shape(B_lab))
    # get theESR freq. for the four families
    freq_ensemble = nv.esr_frequencies_ensemble(B_lab)

    if verbose:
        print(('freq_ensemble', freq_ensemble))

    for i, f in enumerate(freq_ensemble):
        aa = nv.calc_bfields_esr_ensemble_mag(f, verbose)[0]

    # now calculate the B-field from the ensemble esr freq.
    Babs = np.array([nv.calc_bfields_esr_ensemble_mag(f, verbose)[0] for i, f in enumerate(freq_ensemble)])
    # calculate the abs magnetic field
    Bo = np.sqrt(np.sum(B_lab ** 2, axis=1))
    # calculate the difference in the field components (absolute value)
    err = np.abs(Babs - Bo)
    # err = np.std((Babs - Bo)/Bo)

    if len(err) == 1:
        err = err[0]
    if verbose:

        # if err > 1e-1:
        print(('Bfield (T) / theta (deg) / phi (deg)', B, theta, phi))
        print(('B_r', Babs))
        print(('error (%)', err * 100))


    return err


def calc_error_ensemble_fit(B, theta, phi, verbose=False, try_permutations_fit = True, try_permutations_xyz = True, try_permutations_sign = True):
    """
    here we test how well we can reconstruct the x y and z components of the magnetic field
    from the ESR freq. from a N families of NV ESR by fitting to the Hamiltonian

    B, theta, phi: amplitude and angles of magnetic field


    try_permutations_fit: permutate xyz components and optimize for each to reduce error
    try_permutations_xyz: permutate xyz components and pick the one with the lowest reduce error
    """

    if verbose:
        print(('B, theta, phi', B, theta, phi))

    # calculate field in lab frame
    B_lab = nv.B_cart(B, theta, phi)

    # # if B_lab is a vector, the we turn it into a matrix because that is what we expect in the following lines
    # if len(np.shape(B_lab)) == 1:
    #     #         print('reshaping')
    #     B_lab = np.array([B_lab])
    # get the ESR freq. for the four families
    freq_ensemble = nv.esr_frequencies_ensemble(B_lab)
    # now reconstruct the B-field from the ensemble esr freq.
    B_r, theta_r, phi_r = nv.fit_Hamiltonian(freq_ensemble, verbose = verbose,
                                             try_permutations_fit = try_permutations_fit ,
                                             try_permutations_xyz = try_permutations_xyz,
                                             try_permutations_sign = try_permutations_sign)
    # calculate the error in the esr frequencies
    err_esr = nv.fit_err_fun([theta_r, phi_r], B_r, freq_ensemble)

    # calculate the error in the actual fields
    B_lab_r = nv.B_cart(B_r, theta_r, phi_r)
    err_field = np.mean(np.abs(B_lab_r-B_lab)/B)

    # calculate the error in the field amplitudes
    B_lab_r = nv.B_cart(B_r, theta_r, phi_r)
    err_field_abs = np.mean(np.abs(np.abs(B_lab_r)-np.abs(B_lab))/B)

    # calculate the error in the field amplitudes each
    err_field_x = np.mean(np.abs(B_lab_r[0]-B_lab[0])/B)
    err_field_y = np.mean(np.abs(B_lab_r[1] - B_lab[1]) / B)
    err_field_z = np.mean(np.abs(B_lab_r[2] - B_lab[2]) / B)

    # calculate the error in the field amplitudes each
    err_field_x = np.mean(np.abs(B_lab_r[0]-B_lab[0])/B)
    err_field_y = np.mean(np.abs(B_lab_r[1] - B_lab[1]) / B)
    err_field_z = np.mean(np.abs(B_lab_r[2] - B_lab[2]) / B)

    if verbose:

        # if err > 1e-1:
        print(('Bfield (T) / theta (deg) / phi (deg)', B, theta, phi))
        print(('B recovered', B_r, theta_r, phi_r))
        print(('B xyz', B_lab))
        print(('B xyz (recovered)', B_lab_r))


        print(('error esr (%)', err_esr * 100))
        print(('error fields (%)', err_field * 100))
        print(('error fields- abs (%)', err_field_abs * 100))
        print(('error fields- x (%)', err_field_x * 100))
        print(('error fields- y (%)', err_field_y * 100))
        print(('error fields- z (%)', err_field_z * 100))


    return err_esr, err_field, err_field_abs, err_field_x, err_field_y, err_field_z
def test_magnetic_field_fit(Bmax=0.1, verbose=True, try_permutations_fit = True, try_permutations_xyz = True, try_permutations_sign =True):
    """
    here we test how well we can reconstruct the on axis and off axis field from the ESR freq. from a single NV ESR

    Bmax: maximum random field in Teslas

    try_permutations_fit: permutate xyz components and optimize for each to reduce error
    try_permutations_xyz: permutate xyz components and pick the one with the lowest reduce error
    """
    import random

    # get random magnetic field
    B = random.random() * Bmax
    # get random angle field
    theta = random.random() * 90
    phi = random.random() * 180



    error = calc_error_ensemble_fit(B, theta, phi, verbose = verbose,
                                    try_permutations_fit = try_permutations_fit ,
                                    try_permutations_xyz = try_permutations_xyz,
                                    try_permutations_sign = try_permutations_sign)

    if verbose:
        print(('error', error))
        print(('B, theta, phi', B, theta , phi))

    return error, B, theta, phi

def test_magnetic_field_xyz(Bmax=0.1, verbose=True):
    """
    here we test how well we can reconstruct the on axis and off axis field from the ESR freq. from a single NV ESR

    Bmax: maximum random field in Teslas
    """
    import random

    # get random magnetic field
    B = random.random() * Bmax
    # get random angle field
    theta = random.random()*90
    phi = random.random()*180

    error = calc_error_ensemble_xyz(B, theta, phi, verbose)

    if verbose:
        print(('error', error))
        print(('B, theta, phi', B, theta, phi))

    return error, B, theta, phi
def calc_error_ensemble_xyz(B, theta, phi, verbose=False):
    """
    here we test how well we can reconstruct the x y and z components of the magnetic field
    from the ESR freq. from a N families of NV ESR

    B, theta, phi: amplitude and angles of magnetic field
    """

    if verbose:
        print(('B, theta, phi', B, theta * 180 / np.pi, phi * 180 / np.pi))

    # calculate field in lab frame
    B_lab = B_vec(B, theta, phi)

    # if B_lab is a vector, the we turn it into a matrix because that is what we expect in the following lines
    if len(np.shape(B_lab)) == 1:
        #         print('reshaping')
        B_lab = np.array([B_lab])
    # get the ESR freq. for the four families
    freq_ensemble = B_field_to_esr_ensemble_freq(B_lab)
    # now calculate the B-field from the ensemble esr freq.
    B_r = np.array([calc_bfields_esr_ensemble_xyz(f, verbose) for i, f in enumerate(freq_ensemble)])

    # calculate the abs magnetic field
    Bo = (np.mean(B_r ** 2, axis=1) + np.mean(B_lab ** 2, axis=1)) / 2
    # calculate the difference in the field components (absolute value)
    err = np.mean(np.sqrt((B_r ** 2 - B_lab ** 2)** 2), axis=1)
    # normalize by the total abs magnetic field
    err = np.array([0 if b == 0 else e / b for e, b in zip(err, Bo)])

    if len(err) == 1 and verbose:

        # if err > 1e-1:
        print(('Bfield (T) / theta (deg) / phi (deg)', B, theta * 180 / np.pi, phi * 180 / np.pi))
        print(('Bxyz', B_lab))
        print(('B_r', B_r))
        print(('error (%)', err * 100))
    if len(err) == 1:
        err = err[0]

    return err




if __name__ == '__main__':
    # TESTING THE ABSolute field RECOVERY =================================
    # nv.calc_error_ensemble_mag(0.23973219239761479, 92.8218341608137/180*np.pi, 0.0, verbose = True)
    # import random
    #
    # B = random.random() * 0.1
    # theta = random.random() * 90
    # phi = random.random() * 180
    # calc_error_ensemble_fit(B, theta, phi, verbose=True)


    B, theta, phi = 0.04234303209482937, 89.6993911391464, 20.10693684160114

    err =calc_error_ensemble_fit(B, theta, phi, verbose=True)

    print(err)
    # error, B, theta, phi = test_magnetic_field_fit()


    # print('----', error, B, theta, phi)



