import numpy as np

import ipywidgets as widgets
from IPython.display import display

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import patches

from pylabcontrol.core.scripts import Script

from b26_toolkit.plotting.plots_2d import plot_fluorescence_new

import queue
import threading

import scipy.spatial

freq_to_mag = 1. / (2 * 2.8e6)

def pair_nv_images(img_path_old, img_path_new, return_queue):
    done_queue = setup_button()

    # launch manual_correction function on separate thread, as otherwise buttons will not respond until that
    # function ends, that function is waiting on buttons to continue, and you deadlock
    thread = threading.Thread(target=perform_affine_transform, args=(img_path_old, img_path_new, done_queue, return_queue))
    thread.start()

def setup_button():
    # define queues for inter-thread communication
    done_queue = queue.Queue()

    # define buttons to be added to display
    button_done = widgets.Button(description="Done")

    # display all widgets
    display(button_done)

    def button_done_clicked(b):
        done_queue.put('DONE')

    button_done.on_click(button_done_clicked)

    return done_queue

def calc_transform_matrix(pt_list_1, pt_list_2):
    x = []
    for pt in pt_list_1:
        x.append([pt[0], pt[1], 1, 0, 0, 0])
        x.append([0, 0, 0, pt[0], pt[1], 1])

    x_prime = []
    for pt in pt_list_2:
        x_prime.extend(pt)

    return np.dot(np.linalg.inv(x), x_prime)

def perform_affine_transform(img_path_old, img_path_new, done_queue, return_queue):
    import time

    pt11_widget = widgets.HTML("Pt 1: 0,0")
    pt12_widget = widgets.HTML("Pt 2: 0,0")
    pt13_widget = widgets.HTML("Pt 3: 0,0")
    pt21_widget = widgets.HTML("Pt 1: 0,0")
    pt22_widget = widgets.HTML("Pt 2: 0,0")
    pt23_widget = widgets.HTML("Pt 3: 0,0")

    display(pt11_widget)
    display(pt12_widget)
    display(pt13_widget)
    display(pt21_widget)
    display(pt22_widget)
    display(pt23_widget)

    #     f = plt.figure(figsize=(12, 6))

    #     ax0 = plt.subplot(121)
    #     ax1 = plt.subplot(122)

    f = plt.figure(figsize=(16, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    def onclick_old(event):
        if ax0.in_axes(event):
            if event.key == '1':
                pt11_widget.value = 'Pt 1: ' + str(event.xdata) + ',' + str(event.ydata)
            elif event.key == '2':
                pt12_widget.value = 'Pt 2: ' + str(event.xdata) + ',' + str(event.ydata)
            elif event.key == '3':
                pt13_widget.value = 'Pt 3: ' + str(event.xdata) + ',' + str(event.ydata)
        elif ax1.in_axes(event):
            if event.key == '1':
                pt21_widget.value = 'Pt 1: ' + str(event.xdata) + ',' + str(event.ydata)
            elif event.key == '2':
                pt22_widget.value = 'Pt 2: ' + str(event.xdata) + ',' + str(event.ydata)
            elif event.key == '3':
                pt23_widget.value = 'Pt 3: ' + str(event.xdata) + ',' + str(event.ydata)

    cid = f.canvas.mpl_connect('button_press_event', onclick_old)

    data_im_old = Script.load_data(img_path_old)
    image_old = data_im_old['image_data']
    extent_old = data_im_old['extent']
    nv_pts_old = data_im_old['nv_locations']

    plot_fluorescence_new(image_old, extent_old, ax0)
    f.get_axes()[2].remove()

    data_im_new = Script.load_data(img_path_new)
    image_new = data_im_new['image_data']
    extent_new = data_im_new['extent']
    nv_pts_new = data_im_new['nv_locations']

    plot_fluorescence_new(image_new, extent_new, ax1)
    f.get_axes()[2].remove()

    plt.show()

    while done_queue.empty():
        time.sleep(.5)

    pt11 = list(map(float, pt11_widget.value.replace(" ", "").replace(':', ',').split(',')[1:3]))
    pt12 = list(map(float, pt12_widget.value.replace(" ", "").replace(':', ',').split(',')[1:3]))
    pt13 = list(map(float, pt13_widget.value.replace(" ", "").replace(':', ',').split(',')[1:3]))
    pt21 = list(map(float, pt21_widget.value.replace(" ", "").replace(':', ',').split(',')[1:3]))
    pt22 = list(map(float, pt22_widget.value.replace(" ", "").replace(':', ',').split(',')[1:3]))
    pt23 = list(map(float, pt23_widget.value.replace(" ", "").replace(':', ',').split(',')[1:3]))

    Avec = calc_transform_matrix((pt11, pt12, pt13), (pt21, pt22, pt23))
    Amat = [[Avec[0], Avec[1], Avec[2]], [Avec[3], Avec[4], Avec[5]], [0, 0, 1]]
    print(Amat)

    shifted_nv_pts_old = []

    for pt in nv_pts_old:
        circ = patches.Circle((pt[0], pt[1]), .0005, fc='k', ec = 'k')
        ax0.add_patch(circ)

        vec = [pt[0], pt[1], 1]
        vec_prime = np.dot(Amat, vec)

        circ = patches.Circle((vec_prime[0], vec_prime[1]), .0005, fc='k', ec = 'k')
        ax1.add_patch(circ)

        shifted_nv_pts_old.append([vec_prime[0], vec_prime[1]])

    for pt in nv_pts_new:
        circ = patches.Circle((pt[0], pt[1]), .0005, fc='g', ec = 'g')
        ax1.add_patch(circ)

    tree = scipy.spatial.KDTree(shifted_nv_pts_old)
    _, new_to_old_map = tree.query(nv_pts_new, distance_upper_bound=.005)
    #kd tree returns value of len(new_to_old_map) if no match found, change this to -1
    new_to_old_map = [x if x != len(nv_pts_old) else -1 for x in new_to_old_map]

    return_queue.put(new_to_old_map)