import random
from unittest import TestCase
import numpy as np
import random



from b26_toolkit.src.data_analysis import nv_optical_response as nv
from b26_toolkit.src.data_analysis import nv_optical_response_test as nv_test

class NVOpticalResponse(TestCase):
    def setUp(self):
        # get random magnetic field
        self.B_mag = random.random() * 0.1
        # get random angles
        self.theta = random.random() * 180
        self.phi = random.random() * 180
        # self.B_mag, self.phi, self.theta = 0.1, 0.3, 1.3
        pass

    def test01_B_cart(self):
        B_lab = nv.B_cart(self.B_mag, self.theta, self.phi)
        print(B_lab)

    def test02_B_fields_in_NV_frame(self):
        print('testing B_fields_in_NV_frame')
        B_lab = nv.B_cart(self.B_mag, self.theta, self.phi)
        for i in range(4):
            B_NV = nv.B_fields_in_NV_frame(B_lab, i)
            err = np.dot(B_NV, B_NV)-np.dot(B_lab, B_lab)
            print(i, err)

    def test03_B_field_to_esr_ensemble_freq(self):
        print('testing ensemble_freq')
        B_lab = nv.B_cart(self.B_mag, self.theta, self.phi)
        ensemble_freq = nv.esr_frequencies_ensemble(B_lab)

        print('ensemble_freq', ensemble_freq)

    def test04_calc_error_ensemble_mag(self):
        print('testing calc_error_ensemble_mag')

        err = nv_test.calc_error_ensemble_mag(self.B_mag, self.theta, self.phi)

        print('err', err)

    def test_magnetic_field_mag_ensemble(self):
        print('05 testing magnetic_field_mag')
        errs = np.array([nv_test.test_magnetic_field_mag_ensemble(Bmax=1, verbose=False) for i in range(50)])
        print(errs)

    def test06_recovery(self):
        """
        testing if we recover the right on and off axis fields for each NV family
        Returns:

        """
        success = True
        # B, theta, phi = 0.08447537165517684, 1.259739889393438, 1.6332567784906755
        B, theta, phi = self.B_mag, self.theta, self.phi
        B_cart = nv.B_cart(B, theta, phi)

        freq = nv.esr_frequencies_ensemble(B_cart)

        print('B_cart', B_cart)
        for i in range(4):
            BNV = nv.B_fields_in_NV_frame(B_cart, i)

            f = freq[:,i]
            Br = nv.B_field_from_esr(*f)



            err = np.abs((np.abs(Br[0])- np.abs(BNV[2])))/B
            if np.abs(err)>1e-1:
                success = False
                print('err on axis ', err)


            err = np.abs((np.abs(Br[1]) - np.sqrt(BNV[0]**2+BNV[1]**2)) / B)
            if np.abs(err)>1e-1:
                success = False
                print('err off axis ', err)

            if not success:
                print('NV', i, '=================')
                print('BNV', BNV)
                print('Br', Br)
                break

        if not success:

            raise RuntimeWarning


    def test07_inverse_projection(self):
        success = True
        theta = random.random() * 90
        phi = random.random() * 180

        nv.projetion_matrix(theta, phi)

        theta = random.random() * 90
        phi = random.random() * 180

        P = nv.projetion_matrix(theta, phi)
        if np.abs(np.sum(np.abs(np.dot(P, P.T))) - 3) >1e-8:
            success = False

        if not success:
            raise RuntimeError
        else:
            print('SUCCESS: test07_inverse_projection')

    def test08_calc_error_ensemble_fit(self):
        success = True
        B = random.random() * 0.1
        theta = random.random() * 90
        phi = random.random() * 180

        err = nv_test.calc_error_ensemble_fit(B, theta, phi, verbose=False)

        for e in err:
            if e>1e-2:
                success = False
        if not success:
            print(err)
            raise RuntimeError
        else:
            print('SUCCESS: test08_calc_error_ensemble_fit')

    def test09_spher_cart(self):
        """
        here we check the conversion from cartesian to shperical coordinates
        Returns:

        """
        success = True
        # create random values
        B = random.random() * 0.1
        theta = random.random() * 90
        phi = random.random() * 180
        Bc = nv.B_cart(B, theta, phi)
        # convert back to spherical
        Bs = nv.B_spher(*Bc)
        # check for differences
        if np.abs(Bs[0]-B)/B>1e-4:
            success = False
        if np.abs(Bs[1]-theta)/180>1e-4:
            success = False
        if np.abs(Bs[2]-phi)/180>1e-4:
            success = False

        if not success:
            raise RuntimeError
        else:
            print('conversion spher->cart  passed!')

    def test10_cart_spher(self):
        """
        here we check the conversion from cartesian to shperical coordinates
        Returns:

        """
        success = True
        # create random values
        Bx = random.random() * 0.1
        By = random.random() * 0.1
        Bz = random.random() * 0.1

        Bs = nv.B_spher(Bx, By, Bz)
        # convert back to cartesian
        Bc = nv.B_cart(*Bs)

        # print('errors', np.abs((Bc[0]-Bx)/Bx), np.abs((Bc[1]-By)/By), np.abs((Bc[2]-Bz)/Bz))

        # check for differences
        if np.abs((Bc[0]-Bx)/Bx)>1e-7:
            success = False
        if np.abs((Bc[1]-By)/By)>1e-7:
            success = False
        if np.abs((Bc[2]-Bz)/Bz)>1e-7:
            success = False
        print(not success)
        if not success:
            raise RuntimeError
        else:
            print('conversion cart->spher passed!')