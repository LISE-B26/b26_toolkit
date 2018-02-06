
# here we define some plotting functions to plot the data generated in fields

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_Bfield_mag(data, title='', offset_field=np.zeros(3), ax=None, cmap_name=None):
    """

    :param data: data that has the magnetic fields as a pandas dataframe with elements 'Bx', 'By', 'Bz' and 'x', 'y', 'z' 'B' is in Teslas
    :param title: optional title for plot
    :param offset_field: optional offset field along z
    :param ax: optional axis object on which data is plotted
    :return: figure object with magnetic field plot
    """

    data['Bmag'] = np.sqrt((data['Bx'] + offset_field[0]) ** 2 +
                (data['By'] + offset_field[1]) ** 2 +
                (data['Bz'] + offset_field[2]) ** 2)

    ax = plot_NV_property_map(data, prop='Bmag', title=title,ax=ax, cmap_name=None)


    return ax


def plot_NV_property_map(data, prop='shift', title='', gammaNV=28e9, wo=500e-9, ax=None, cmap_name = None):
    """

    :param data: data that has the magnetic fields as a pandas dataframe with elements 'Gxy' and 'x', 'y', 'z' 'Gxy' is in T/um
    :param prop: property to be plotted, use any of the following shift, contrast, Bperp, Bpar, Gxy
    :param title: optional title for plot

    :param ax: optional axis object on which data is plotted
    :param gammaNV = 28e9 # (ESR shift is 28 GHz/T)
    :param wo = 500e-9 # focal spot size
    :param cmap_name name of colormap

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
    elif prop == 'G':
        C = data['G']
        title = 'gradient (T/um)' if title is '' else ''
    else:
        print('plot type not recognized try prop = shift, contrast, Bperp, Bpar, Gxy')
        raise NotImplementedError
    C = C.reshape(Ny, Nx)

    X = data['x'].reshape(Ny, Nx)
    Y = data['y'].reshape(Ny, Nx)

    xmin, xmax = np.min(X), np.max(X)
    ymin, ymax = np.min(Y), np.max(Y)

    if ax is None:
        fig, ax = plt.subplots()

    fig = ax.get_figure()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    if cmap_name is not None:
        im = ax.pcolormesh(X, Y, C, cmap=plt.get_cmap(cmap_name))
    else:
        im = ax.pcolormesh(X, Y, C)

    fig.colorbar(im, cax=cax, orientation='vertical')


    ax.set_title(title)
    ax.set_xlabel('x ($\mu m$)')
    ax.set_ylabel('y ($\mu m$)')

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    ax.set_aspect('equal')

    return ax


def plot_arrow(x, ax, arrow_length=1, arrow_angle=200):
    """
    plots arrows of length arrow_length with angle arrow_angle at positions defined in x


    :param x: pandas dataframe with entries x and y that determine the end point of the arrow
    :param ax: axis object on which to plot the arrow
    :param arrow_length:
    :param arrow_angle:
    :return:
    """
    dx = arrow_length * np.cos(arrow_angle * np.pi / 180)
    dy = arrow_length * np.sin(arrow_angle * np.pi / 180)
    # calculate some values for the arrow plot
    for index, row in x.iterrows():
        xo, yo = np.float(row['x']), np.float(row['y'])
        ax.arrow(xo + dx, yo + dy, -dx, -dy,
                 head_width=0.3, head_length=0.2, fc='w', ec='w', head_starts_at_zero=False,
                 length_includes_head=True)
    return ax


def plot_ring(radius, ax):
    """
    plots ring onto axes object ax


    :param radius: radius of ring

    :return: axes object
    """

    x = radius * np.cos(np.linspace(0, 2 * np.pi))
    y = radius * np.sin(np.linspace(0, 2 * np.pi))
    # calculate some values for the arrow plot
    ax.plot(x, y, 'w', lw='2')

    return ax
