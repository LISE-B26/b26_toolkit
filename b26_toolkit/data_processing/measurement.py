import numpy as np

def dist_pt_line(pt, line):
    """
    gives the distance from the point pt to the line line
    pt: coordinates of point vector of length 2
    line: coordinates of line as a 2x2 matrix line = [[x1, y1], [x2, y2]]

    Returns: minimal distance between point and line
    """
    line = np.array(line)
    pt = np.array(pt)

    assert len(pt) == 2
    assert np.shape(line) == (2, 2)

    v1 = pt - line[0]
    v2 = line[1] - line[0]
    dist = abs(np.cross(v1, v2) / np.linalg.norm(v2))
    return dist