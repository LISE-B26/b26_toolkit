from pylabcontrol.core import Script
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse


def save_image(dir):
    data = Script.load_data(dir)
    image = data['image_data']
    filename= 'image.png'
    plt.imsave(dir+"/"+filename, image,cmap="Greys")
    return True

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=dir_path)
    args = parser.parse_args()
    dir = r'{}'.format(args.path)
    save_image(dir)




