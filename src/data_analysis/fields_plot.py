
# here we define some plotting functions to plot the data generated in fields

import numpy as np
import matplotlib.pyplot as plt



def plot_Bfield_mag(data, title='', offset_field=np.zeros(3), ax=None):
    """

    :param data: data that has the magnetic fields as a pandas dataframe with elements 'Bx', 'By', 'Bz' and 'x', 'y', 'z' 'B' is in Teslas
    :param title: optional title for plot
    :param offset_field: optional offset field along z
    :param ax: optional axis object on which data is plotted
    :return: figure object with magnetic field plot
    """

    Nx, Ny = len(np.unique(data['x'])), len(np.unique(data['y']))

    C = np.sqrt((data['Bx'] + offset_field[0]) ** 2 +
                (data['By'] + offset_field[1]) ** 2 +
                (data['Bz'] + offset_field[2]) ** 2)
    C = C.reshape(Ny, Nx)

    X = data['x'].reshape(Ny, Nx)
    Y = data['y'].reshape(Ny, Nx)

    xmin, xmax = np.min(X), np.max(X)
    ymin, ymax = np.min(Y), np.max(Y)

    if ax is None:
        fig = plt.figure(figsize=(15, 4))

        plt.pcolormesh(X, Y, C * 1e4)
        plt.colorbar(label='B (G)')
        plt.title(title)
        plt.xlabel('x ($\mu m$)')
        plt.ylabel('y ($\mu m$)')

        plt.xlim([xmin, xmax])
        plt.ylim([ymin, ymax])
        #     plt.clim([0, 500])

        plt.axes().set_aspect('equal')
    else:
        raise NotImplementedError
        #todo: implemet plotting on given axis object, this is a bit tricky if we want a colorbar
        fig = plt.figure(figsize=(15, 4))

    return fig

def plot_G(data, grad_dir = None, title='', ax=None):
    """

    :param data: data that has the magnetic fields as a pandas dataframe with elements 'G' and 'x', 'y', 'z'.
    If it doesn't contain 'G' the value grad_dir has to be specified. G is in T/um
    :param grad_dir: (optional) string that is the key for the dataset to be plotted
    :param title: optional title for plot
    :param ax: optional axis object on which data is plotted
    :return: figure object with magnetic field plot
    """

    Nx, Ny = len(np.unique(data['x'])), len(np.unique(data['y']))

    if grad_dir is None:
        grad_dir = 'G'
    C = data[grad_dir].as_matrix()
    C = C.reshape(Ny, Nx)

    X = data['x'].reshape(Ny, Nx)
    Y = data['y'].reshape(Ny, Nx)

    xmin, xmax = np.min(X), np.max(X)

    ymin, ymax = np.min(Y), np.max(Y)

    if ax is None:
        fig = plt.figure(figsize=(15, 4))

        plt.pcolormesh(X, Y, C)
        plt.colorbar(label='{:s} (T/$\mu$m)'.format(grad_dir))
        plt.title(title)
        plt.xlabel('x ($\mu m$)')
        plt.ylabel('y ($\mu m$)')

        plt.xlim([xmin, xmax])
        plt.ylim([ymin, ymax])
        #     plt.clim([0, 500])

        plt.axes().set_aspect('equal')
    else:
        raise NotImplementedError
        #todo: implemet plotting on given axis object, this is a bit tricky if we want a colorbar
        fig = plt.figure(figsize=(15, 4))

    return fig