import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import scipy.optimize as sopt
import pims
import scipy.spatial
import scipy.ndimage
from collections import Counter
import glob

import trackpy as tp
from pims import pipeline
from skimage.color import rgb2gray
rgb2gray_pipeline = pipeline(rgb2gray)


from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()



def get_maxcoors(image, image_size):
    pixel_max = np.argmax(image)
    return(pixel_max%image_size, pixel_max/image_size)

def get_maxcoors_array(image_array, image_size):
    maxes  = []
    for image in image_array:
        pixel_max = np.argmax(image)
        maxes.append((pixel_max%image_size, pixel_max/image_size))
    return maxes


if __name__ == '__main__':
    filepath = 'Z:\\Lab\\Lev\\videos\\20171206_LevSample_4\\20180114_reset_bead\\cycle_2\\20180115_ringdown\\zmode_93hz_1.avi'
    v = pims.Video(filepath)
    # processed_v = rgb2gray_pipeline(v)

    # gaussian_filter_pipeline = pipeline(scipy.ndimage.filters.gaussian_filter)
    # filtered = gaussian_filter_pipeline(processed_v, 2)
    # filtered = processed_v

    image_size = 200

    import time
    #48 seconds
    start = time.time()

    with Parallel(n_jobs = num_cores) as parallel:
        processed_v = parallel(delayed(rgb2gray)(image) for image in v[:100000])
        filtered = processed_v

        parallel(delayed(get_maxcoors_array)(image_array, image_size) for image_array in np.array_split(filtered, 1000))
        print((time.time() - start))

    # start = time.time()
    # max_coors = []
    # for index, image in enumerate(filtered[0:10000]):
    #     pixel_max = np.argmax(image)
    #     max_coors.append((pixel_max % image_size, pixel_max / image_size))
    # print(time.time() - start)