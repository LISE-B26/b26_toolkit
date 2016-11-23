import numpy as np

def cartesian_to_spherical(vector):
    [x, y, z] = vector
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)

    return ([r, theta, phi])


def spherical_to_cartesian(vector):
    [r, theta, phi] = vector
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return ([x, y, z])