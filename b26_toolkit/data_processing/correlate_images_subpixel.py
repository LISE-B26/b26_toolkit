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
from PIL import Image as im
from scipy import signal
from scipy.ndimage import shift as image_shift
import trackpy as tp
from skimage.filters import sobel
from skimage import data
from skimage.feature import register_translation
from skimage.feature.register_translation import _upsampled_dft
from b26_toolkit.data_processing import find_image_shift

def compare_galvos(image, image_extent, offset_image):
    """
    Takes two images and finds the necessary translation of the second image to match the first. Unlike correlate_images,
    we upsample the images to estimate any subpixel shifts. Also, the two galvos MUST have the same extents!

    Args:
        reference_image: numpy 2D array of pixel values
        reference_image_bounds: numpy array with 4 elements containing the voltage bounds of the reference image
        shifted_image: numpy 2D array of pixel values
        shifted_image_bounds: numpy array with 4 elements containing the voltage bounds of the shifted image

    Returns: ordered pair (y_shift, x_shift) of pixel values

    """

    pix2vol = pixel_to_voltage_conversion_factor(image.shape, image_extent)
    # pixel precision first
    shift, error, diffphase = register_translation(image, offset_image)
    shift = -shift
    shift_volts = shift * pix2vol

    # subpixel precision
    shift_subpixel, error, diffphase = register_translation(image, offset_image, 100)
    shift_subpixel = -shift_subpixel
    shift_subpixel_volts = shift_subpixel * pix2vol

    print("Detected offset (y, x): {}px or {}V".format(shift,shift_volts))
    print("Detected subpixel offset (y, x): {}px or {}V".format(shift_subpixel,shift_subpixel_volts))

    image_product = np.fft.fft2(image) * np.fft.fft2(offset_image).conj()
    cc_image = np.fft.fftshift(np.fft.ifft2(image_product))

    corrected_galvo = image_shift(offset_image, -shift_subpixel)
    difference = image - corrected_galvo

    return shift_subpixel_volts[1],shift_subpixel_volts[0],difference


def pixel_to_voltage_conversion_factor(image_shape, image_extent):
    # COMMENT_ME
    image_x_len, image_y_len = image_shape
    image_x_min, image_x_max, image_y_max, image_y_min = image_extent
    x_voltage = (image_x_max - image_x_min) / image_x_len
    y_voltage = (image_y_max - image_y_min) / image_y_len
    return x_voltage, y_voltage


