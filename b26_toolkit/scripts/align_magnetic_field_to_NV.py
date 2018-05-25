import glob
import pandas as pd
import numpy as np
import time
import threading
from copy import deepcopy

from pylabcontrol.core import Script, Parameter
from b26_toolkit.scripts import FindNV, SetMagneticCoils, ESR, Take_And_Correlate_Images
from b26_toolkit.data_processing.esr_signal_processing import fit_esr
from b26_toolkit.data_processing.coordinate_conversions import cartesian_to_spherical

SPLITTING_TO_FIELD = 5.6 #in MHz/Gauss

class AlignFieldToNV(Script):
    """
This script determines the orientation of a given NV
    """

    _DEFAULT_SETTINGS = [
        Parameter('field_magnitude', 10.0, float, 'Magnitude of magnetic field to apply for alignment, in Gauss'),
        Parameter('stabilize_time', 10.0, float, 'time to wait after setting magnetic field for system to stabilize, in seconds')
    ]

    _INSTRUMENTS = {}

    _SCRIPTS = {'FindNV': FindNV, 'SetMagneticCoils': SetMagneticCoils, 'ESR': ESR, 'Correlate': Take_And_Correlate_Images}


    def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_function= log_function, data_path = data_path)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        # DOCSTRING HERE
        # vector = [proj_x, proj_y, proj_z]
        # cross_slopes = [proj_xy, proj_xz, proj_yz]
        def determine_signs(vector, cross_slopes):
            """
            This function takes in a vector containing the unsigned direction of a vector as well as an array containing
            the projections onto x=y, x=z, y=z, then returns a signed version of the directional vector
            Args:
                vector: the nv vector in the format [x,y,z]
                cross_slopes: the projections onto the additional axes in the format [xy, xz, yz]

            Returns: vector with each component multiplied by the correct signs

            """
            if cross_slopes[0] > max(vector[0], vector[1]):
                if cross_slopes[1] > max(vector[0], vector[2]):
                    if cross_slopes[2] > max(vector[1], vector[2]):
                        return {'x': vector[0] * 1, 'y': vector[1] * 1, 'z': vector[2] * 1}
                    else:
                        return None
                else:
                    if cross_slopes[2] > max(vector[1], vector[2]):
                        return None
                    else:
                        return {'x': vector[0] * 1, 'y': vector[1] * 1, 'z': vector[2] * -1}
            else:
                if cross_slopes[1] > max(vector[0], vector[2]):
                    if cross_slopes[2] > max(vector[1], vector[2]):
                        return None
                    else:
                        return {'x': vector[0] * 1, 'y': vector[1] * -1, 'z': vector[2] * 1}
                else:
                    if cross_slopes[2] > max(vector[1], vector[2]):
                        return {'x': vector[0] * -1, 'y': vector[1] * 1, 'z': vector[2] * 1}
                    else:
                        return None


        init_pt = deepcopy(self.scripts['FindNV'].settings['initial_point'])

        mag = self.settings['field_magnitude']
        self.scripts['SetMagneticCoils'].settings['magnetic_fields']['coordinate_system'] = 'Cartesian'
        magnetic_field_configurations = [
            {'x/r_field': 0, 'y/theta_field': 0, 'z/phi_field': 0},
            {'x/r_field': mag, 'y/theta_field': 0, 'z/phi_field': 0},
            {'x/r_field': 0, 'y/theta_field': -1*mag, 'z/phi_field': 0},
            {'x/r_field': 0, 'y/theta_field': 0, 'z/phi_field': mag},
            {'x/r_field': mag, 'y/theta_field': -1*mag, 'z/phi_field': 0},
            {'x/r_field': mag, 'y/theta_field': 0, 'z/phi_field': mag},
            {'x/r_field': 0, 'y/theta_field': -1*mag, 'z/phi_field': mag}
        ]
        axes = ['0', 'x', 'y', 'z', 'xy', 'xz', 'yz']

        results = {}
        for field, axis in zip(magnetic_field_configurations, axes):
            self.scripts['SetMagneticCoils'].settings['magnetic_fields'] = field
            self.scripts['SetMagneticCoils'].run()
            time.sleep(self.settings['stabilize_time'])
            self.scripts['Correlate'].run()
            shift = self.scripts['Correlate'].data['shift']
            new_pt = {'x': init_pt['x'] + shift[0], 'y': init_pt['y'] + shift[1]}
            self.scripts['FindNV'].settings['initial_point'].update(new_pt)

            self.scripts['FindNV'].run()
            self.scripts['ESR'].run()

            fit_params = fit_esr(self.scripts['ESR'].data['frequency'], self.scripts['ESR'].data['data'])
            if len(fit_params) == 4: #if unsplit
                results.update({axis: 0})
            else:
                meas_field = (fit_params[5]-fit_params[4])/1e6/SPLITTING_TO_FIELD
                results.update({axis: meas_field})
            print(results)

        nv_axis = determine_signs([results['x'], results['y'], results['z']], [results['xy'], results['xz'], results['yz']])
        switch_polarity = []
        if nv_axis['x'] < 0:
            switch_polarity.append('x')
        if nv_axis['y'] < 0:
            switch_polarity.append('y')
        if nv_axis['z'] < 0:
            switch_polarity.append('z')
        if switch_polarity == []:
            switch_polarity.append('No Switch Required')

        # [_, theta, phi] = cartesian_to_spherical(nv_axis)
        #calculates angle based on the always positive results, so always returns angles in first octant. Thus, after
        #switching polarity, then just plug this angle in
        [_, theta, phi] = cartesian_to_spherical([results['x'], results['y'], results['z']])


        self.data = {'nv_axis_direction': [theta, phi], 'switch_polarity': switch_polarity}

    #must be passed figure with galvo plot on first axis
    def _plot(self, axes_list):
        if self._current_subscript_stage['current_subscript'] == self.scripts['FindNV']:
            self.scripts['FindNV']._plot(axes_list)
        elif self._current_subscript_stage['current_subscript'] == self.scripts['ESR']:
            self.scripts['ESR']._plot(axes_list)
        elif self._current_subscript_stage['current_subscript'] == self.scripts['Correlate']:
            self.scripts['Correlate']._plot(axes_list)

    #must be passed figure with galvo plot on first axis
    def _update_plot(self, axes_list):
        if self._current_subscript_stage['current_subscript'] == self.scripts['FindNV']:
            self.scripts['FindNV']._update_plot(axes_list)
        elif self._current_subscript_stage['current_subscript'] == self.scripts['ESR']:
            self.scripts['ESR']._update_plot(axes_list)
        elif self._current_subscript_stage['current_subscript'] == self.scripts['Correlate']:
            self.scripts['Correlate']._update_plot(axes_list)

if __name__ == '__main__':
    script, failed, instruments = Script.load_and_append(script_dict={'AlignFieldToNV':'AlignFieldToNV'})

    print(script)
    print(failed)
    # print(instruments)