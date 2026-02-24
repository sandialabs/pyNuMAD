from scipy.spatial.transform import Rotation
import numpy as np
from numpy import ndarray


def rotation(axis, angle):
    """
    Designed to replace matlab's makehgtform
    """
    r = Rotation.from_euler(axis, angle)
    rmatrix = np.eye(4)
    rmatrix[0:3, 0:3] = r.as_matrix()
    return rmatrix


def translation(xtrans, ytrans, ztrans):
    tmatrix = np.eye(4)
    tmatrix[0:3, 3] = [xtrans, ytrans, ztrans]
    return tmatrix


def rotate2d(xyin: ndarray, angle: float) -> ndarray:
    """Rotate a set of 2-D points by *angle* radians about the origin.

    Parameters
    ----------
    xyin : ndarray
        Nx2 array of (x, y) points.
    angle : float
        Rotation angle in radians.

    Returns
    -------
    ndarray
        Nx2 array of rotated (x, y) points.
    """
    xyout1 = np.cos(angle) * xyin[:, 0] - np.sin(angle) * xyin[:, 1]
    xyout2 = np.sin(angle) * xyin[:, 0] + np.cos(angle) * xyin[:, 1]
    return np.stack((xyout1, xyout2), axis=1)


if __name__ == "__main__":
    print(translation(2, 3, 4))
    print(rotation("z", 1))
