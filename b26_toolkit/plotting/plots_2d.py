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
from matplotlib.ticker import FormatStrFormatter, FuncFormatter

# todo: delete plot_fluorescence and refactor plot_fluorescence_new to plot_fluorescence
def plot_fluorescence(image_data, extent, axes_image, implot=None, cbar=None, max_counts=-1, axes_colorbar=None,
                      labels = None):
    """

    Args:
        image_data: 2D - array
        extent: vector of length 4, i.e. [x_min, x_max, y_max, y_min]
        axes: axes object on which to plot
        implot: reference to image plot
        labels: labels for plotting [title, x_label, y_label, cbar_label]
    Returns:

    """
    fig = axes_image.get_figure()
    if labels is None:
        labels = ['Confocal Image', r'V$_x$ [V]', r'V$_y$ [V]', 'kcounts/sec']

    if axes_colorbar is None:
        # try to figure out if there is a axis for the colorbar
        fig = axes_image.get_figure()
        number_of_axes = len(fig.axes)
        for index in range(number_of_axes):
            if fig.axes[index] == axes_image and index < number_of_axes - 1:
                axes_colorbar = fig.axes[index + 1]

    if implot is None:
        if max_counts > 0:
            implot = axes_image.imshow(image_data, cmap='pink', interpolation="nearest", extent=extent, vmax=max_counts)
        else:
            implot = axes_image.imshow(image_data, cmap='pink', interpolation="nearest", extent=extent)

        title, x_label, y_label, cbar_label = labels
        axes_image.set_xlabel(x_label)
        axes_image.set_ylabel(y_label)
        axes_image.set_title(title)

        #axes_image.set_xlabel(r'V$_x$ [V]')
        #axes_image.set_ylabel(r'V$_y$ [V]')
        #axes_image.set_title('Confocal Image')
    else:
        implot.set_data(image_data)

    if not max_counts > 0:
        implot.autoscale()

    if axes_colorbar is None and cbar is None:
        cbar = fig.colorbar(implot, label=cbar_label)
    elif cbar is None:
        cbar = fig.colorbar(implot, cax=axes_colorbar, label=cbar_label)
    else:
        cbar.update_bruteforce(implot)
    # todo: tightlayout warning test it this avoids the warning:
    fig.set_tight_layout(True)
    # fig.tight_layout()

    return implot, cbar


def update_fluorescence(image_data, axes_image, max_counts = -1):
    """
    updates a the data in a fluorescence  plot. This is more efficient than replotting from scratch
    Args:
        image_data: 2D - array
        axes_image: axes object on which to plot
        implot: reference to image plot
    Returns:

    """

    if max_counts >= 0:
        image_data = np.clip(image_data, -1, max_counts)

    implot = axes_image.images[0]
    colorbar = implot.colorbar

    implot.set_data(image_data)

    implot.autoscale()

    colorbar_min = np.min(np.where(image_data > 0, image_data, np.inf))
    if np.isinf(colorbar_min):
        colorbar_min = -1
    implot.set_clim(colorbar_min, None)

    if colorbar is not None:
        if max_counts < 0 or np.max(image_data) < max_counts:
            colorbar_max = np.max(image_data)
        else:
            colorbar_max = max_counts

        colorbar_labels = [np.floor(x) for x in np.linspace(colorbar_min, colorbar_max, 5, endpoint=True)]
        if np.abs(colorbar_max - colorbar_min) > 4:
            colorbar.set_ticks(colorbar_labels)

        colorbar.mappable.set_clim(colorbar_min, colorbar_max)
        colorbar.update_normal(implot)


def plot_fluorescence_new(image_data, extent, axes_image, max_counts=-1, colorbar = None, labels = None, aspect=1):
    """
    plots fluorescence data in a 2D plot
    Args:
        image_data: 2D - array
        extent: vector of length 4, i.e. [x_min, x_max, y_max, y_min]
        axes_image: axes object on which to plot
        max_counts: cap colorbar at this value if negative autoscale
        labels: labels for plotting [title, x_label, y_label, cbar_label]

    Returns:

    """
    if max_counts >= 0:
        image_data = np.clip(image_data, -1, max_counts)

    if labels is None:
        labels = ['Confocal Image', r'V$_x$ [V]', r'V$_y$ [V]', 'kcounts/sec']

    if len(image_data[0]) == 1 or len(image_data) == 1:
        extra_x_extent, extra_y_extent = 0, 0
    else:
        extra_x_extent = (extent[1]-extent[0])/float(2*(len(image_data[0])-1))
        extra_y_extent = (extent[2]-extent[3])/float(2*(len(image_data)-1))
    extent = [extent[0] - extra_x_extent, extent[1] + extra_x_extent, extent[2] + extra_y_extent, extent[3] - extra_y_extent]

    fig = axes_image.get_figure()

    implot = axes_image.imshow(image_data, cmap='inferno', interpolation="nearest", extent=extent, aspect=aspect)

    implot.autoscale()
    colorbar_min = np.min(np.where(image_data > 0, image_data, np.inf))
    if np.isinf(colorbar_min):
        colorbar_min = -1
    implot.set_clim(colorbar_min, None)

    title, x_label, y_label, cbar_label = labels
    axes_image.set_xlabel(x_label)
    axes_image.set_ylabel(y_label)
    axes_image.set_title(title)

    def fmt(x, pos):
        return round(x, 4)

    axes_image.xaxis.set_major_formatter(FuncFormatter(fmt))
    axes_image.yaxis.set_major_formatter(FuncFormatter(fmt))

    axes_image.tick_params(axis="x", which="both", rotation=90)

    if max_counts < 0 or np.max(image_data) < max_counts:
        colorbar_max = np.max(image_data)
    else:
        colorbar_max = max_counts
    colorbar_labels = [np.floor(x) for x in np.linspace(colorbar_min, colorbar_max, 5, endpoint=True)]

    if colorbar is None:
        colorbar = fig.colorbar(implot, label=cbar_label)
    else:
        colorbar = fig.colorbar(implot, cax=colorbar.ax, label=cbar_label)

    if np.abs(colorbar_max - colorbar_min) > 4:
        colorbar.set_ticks(colorbar_labels)

    colorbar.mappable.set_clim(colorbar_min, colorbar_max)
    colorbar.update_normal(implot)

def plot_fluorescence_pos(image_data, extent, axes_image, max_counts = -1, colorbar = None):
    """
    plots fluorescence data in a 2D plot
    Args:
        image_data: 2D - array
        extent: vector of length 4, i.e. [x_min, x_max, y_max, y_min]
        axes_image: axes object on which to plot
        max_counts: cap colorbar at this value if negative autoscale

    Returns:

    """

    extent = [extent[2], extent[3], extent[1], extent[0]]

    if max_counts >= 0:
        image_data = np.clip(image_data, 0, max_counts)

    fig = axes_image.get_figure()

    implot = axes_image.imshow(image_data, cmap='pink', interpolation="nearest", extent=extent)
    axes_image.set_xlabel(r'pos$_{outer}$ [mm]')
    axes_image.set_ylabel(r'pos$_{inner}$ [mm]')
    axes_image.set_title('Position scan: NV fluorescence')

    # explicitly round x_ticks because otherwise they have too much precision (~17 decimal points) when displayed
    # on plot
    axes_image.set_xticklabels([round(xticklabel, 4) for xticklabel in axes_image.get_xticks()], rotation=90)

    if np.min(image_data)<200:
        colorbar_min = 0
    else:
        colorbar_min = np.min(image_data)

    if max_counts < 0:
        colorbar_max = np.max(image_data)
        implot.autoscale()
    else:
        colorbar_max = max_counts
    colorbar_labels = [np.floor(x) for x in np.linspace(colorbar_min, colorbar_max, 5, endpoint=True)]

    if max_counts <= 0:
        implot.autoscale()

    if colorbar is None:
        colorbar = fig.colorbar(implot, label='kcounts/sec')
        colorbar.set_ticks(colorbar_labels)
        colorbar.mappable.set_clim(colorbar_min, colorbar_max)
    else:
        colorbar = fig.colorbar(implot, cax=colorbar.ax, label='kcounts/sec')
        colorbar.set_ticks(colorbar_labels)
        colorbar.mappable.set_clim(colorbar_min, colorbar_max)
