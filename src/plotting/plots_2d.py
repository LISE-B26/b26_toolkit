"""
    This file is part of b26_toolkit, a PyLabControl add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    Foobar is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""

# todo: delete plot_fluorescence and refactor plot_fluorescence_new to plot_fluorescence
def plot_fluorescence(image_data, extent, axes_image, implot=None, cbar=None, max_counts=-1, axes_colorbar=None):
    """

    Args:
        image_data: 2D - array
        extent: vector of length 4, i.e. [x_min, x_max, y_max, y_min]
        axes: axes object on which to plot
        implot: reference to image plot
    Returns:

    """
    fig = axes_image.get_figure()

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
        axes_image.set_xlabel('Vx [V]')
        axes_image.set_ylabel('Vy [V]')
        axes_image.set_title('Confocal Image')
    else:
        implot.set_data(image_data)

    if not max_counts > 0:
        implot.autoscale()

    if axes_colorbar is None and cbar is None:
        cbar = fig.colorbar(implot, label='kcounts/sec')
    elif cbar is None:
        cbar = fig.colorbar(implot, cax=axes_colorbar, label='kcounts/sec')
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

    implot = axes_image.images[0]
    colorbar = implot.colorbar

    implot.set_data(image_data)

    if not max_counts > 0:
        implot.autoscale()

    if not colorbar is None:
        colorbar.update_bruteforce(implot)

def plot_fluorescence_new(image_data, extent, axes_image, max_counts = -1, colorbar = None):
    """
    plots fluorescence data in a 2D plot
    Args:
        image_data: 2D - array
        extent: vector of length 4, i.e. [x_min, x_max, y_max, y_min]
        axes_image: axes object on which to plot
        max_counts: cap colorbar at this value if negative autoscale

    Returns:

    """

    fig = axes_image.get_figure()

    if max_counts > 0:
        implot = axes_image.imshow(image_data, cmap='pink', interpolation="nearest", extent=extent, vmax=max_counts)
    else:
        implot = axes_image.imshow(image_data, cmap='pink', interpolation="nearest", extent=extent)
    axes_image.set_xlabel('Vx [V]')
    axes_image.set_ylabel('Vy [V]')
    axes_image.set_title('Confocal Image')

    if not max_counts > 0:
        implot.autoscale()

    if colorbar is None:
        fig.colorbar(implot, label='kcounts/sec')
    else:
        fig.colorbar(implot, cax=colorbar.ax, label='kcounts/sec')
