import numpy as np

def get_freq_array(freq_a, freq_b, num_points, range_type='start_stop'):
    '''
    Construct a list of values through which we will sweep a parameter.
    Function is called get_freq_array, but the array does not have to be frequency. Can be e.g. voltage values for a piezo

    Returns:
        freq_values: array of the parameter values to be tested
        freq_range: the range of the parameter values to be tested, i.e. maximum frequency - minumum frequency
    '''

    if range_type == 'start_stop':
        if freq_a > freq_b:
            raise ValueError('end freq. must be larger than start freq when range_type is start_stop. Abort script')

        return np.linspace(freq_a, freq_b, num_points)

    elif range_type == 'center_range':
        if freq_a < 2 * freq_b:
            raise ValueError('end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
        return np.linspace(freq_a - freq_b / 2, freq_a + freq_b / 2, num_points)

    else:
        raise KeyError('unknown range parameter. Abort script')