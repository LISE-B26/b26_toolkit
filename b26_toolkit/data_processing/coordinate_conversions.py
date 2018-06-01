import numpy as np

def cartesian_to_spherical(vector):
    [x, y, z] = vector
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    theta = np.arccos(z / r) * 180.0 / np.pi
    phi = np.arctan2(y, x) * 180.0 / np.pi

    return ([r, theta, phi])


def spherical_to_cartesian(vector):
    [r, theta, phi] = vector
    theta *= np.pi / 180.0
    phi *= np.pi / 180.0
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return ([x, y, z])