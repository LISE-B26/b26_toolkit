
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


def plot_NV_property_map(data, prop='shift', title='', gammaNV=28e9, wo=500e-9, ax=None):
    """

    :param data: data that has the magnetic fields as a pandas dataframe with elements 'Gxy' and 'x', 'y', 'z' 'Gxy' is in T/um
    :param prop: property to be plotted, use any of the following shift, contrast, Bperp, Bpar, Gxy
    :param title: optional title for plot

    :param ax: optional axis object on which data is plotted
    :param gammaNV = 28e9 # (ESR shift is 28 GHz/T)
    :param wo = 500e-9 # focal spot size

    :return: figure object with magnetic field plot
    """

    Nx, Ny = len(np.unique(data['x'])), len(np.unique(data['y']))

    if prop == 'shift':
        C = data['Gxy'] * gammaNV * wo  # the linewidht broadening due to a gradient in the xy-plane in MHz
        title = 'Broadening (MHz)' if title is '' else ''
    elif prop == 'contrast':
        raise NotImplementedError
        C = data['Contrast']
        title = 'Contrast (%)' if title is '' else ''
    elif prop == 'Bperp':
        C = data['Bperp'] * 1e3
        title = 'off-axis field (mT)' if title is '' else ''
    elif prop == 'Bpar':
        C = data['Bpar'] * 1e3
        title = 'on-axis field (mT)' if title is '' else ''
    elif prop == 'Gxy':
        C = data['Gxy']
        title = 'xy-gradient (T/um)' if title is '' else ''
    elif prop == 'Bmag':
        C = data['Bmag'] * 1e3
        title = 'total field (mT)' if title is '' else ''
    elif prop == 'fp':
        C = data['fp'] * 1e-9
        title = 'ESR transition freq wp/2pi (GHz)' if title is '' else ''
    elif prop == 'fm':
        C = data['fm'] * 1e-9
        title = 'ESR transition freq wm/2pi (GHz)' if title is '' else ''
    else:
        print('plot type not recognized try prop = shift, contrast, Bperp, Bpar, Gxy')
        raise NotImplementedError
    C = C.reshape(Ny, Nx)

    X = data['x'].reshape(Ny, Nx)
    Y = data['y'].reshape(Ny, Nx)

    xmin, xmax = np.min(X), np.max(X)
    ymin, ymax = np.min(Y), np.max(Y)

    if ax is None:
        fig = plt.figure(figsize=(15, 4))

        plt.pcolormesh(X, Y, C)
        plt.colorbar(label=title)
        plt.title(title)
        plt.xlabel('x ($\mu m$)')
        plt.ylabel('y ($\mu m$)')

        plt.xlim([xmin, xmax])
        plt.ylim([ymin, ymax])
        #     plt.clim([0, 500])

        plt.axes().set_aspect('equal')
    else:
        raise NotImplementedError
        # todo: implemet plotting on given axis object, this is a bit tricky if we want a colorbar
        fig = plt.figure(figsize=(15, 4))

    return fig