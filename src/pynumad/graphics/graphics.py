import matplotlib.pyplot as plt
import numpy as np
from pynumad.utils.interpolation import interpolator_wrap


def plot_airfoil(airfoil):
    """Plot airfoil"""
    fig, ax = plt.subplots()
    # ax[0].plot(self.x,self.y,'.-')
    ax.plot(airfoil.coordinates[:, 0], airfoil.coordinates[:, 1], ".-", color="black")
    ax.plot(airfoil.c, airfoil.camber, color="red")
    # mtx = self.maxthick * np.array([1,1])
    # kn = find(self.c >= self.maxthick,1)
    # mty = self.camber(kn) + self.thickness(kn) * np.array([0.5,- 0.5])
    # line(mtx,mty,'LineStyle',':','Color','k')
    # else:
    fig.show()
    return ax


def plot_regions(blade):
    fig, ax = plt.subplots()

    n = blade.keypoints.shape[0]
    for kn in range(n):
        blade.hgKeypoints[kn] = ax.plot(
            blade.keypoints[kn, 2, :],
            blade.keypoints[kn, 0, :],
            blade.keypoints[kn, 1, :],
        )
    fig.show()
    return ax


def plot_geometry(self):
    fig, ax = plt.subplots()
    n = self.geometry.shape[2]
    for k in range(n):
        self.hgGeometry[k] = ax.plot(
            self.geometry[:, 2, k], self.geometry[:, 0, k], self.geometry[:, 1, k]
        )
    fig.show()
    return ax


def plot_profile(blade, k):
    fig, ax = plt.subplots()
    ax.plot(blade.profiles[:, 0, k], blade.profiles[:, 1, k], ".-")
    fig.show()
    return ax


def plot_component(component):
    """
    TODO docstring
    """
    fig, ax = plt.subplots()
    cpx, cpy = component.getcp()
    ax.plot(cpx, cpy)
    x = np.linspace(0, 1, 100)
    y = np.round(interpolator_wrap(cpx, cpy, x, "pchip", 0))
    ax.plot(x, y)
    plt.title(component.name)
    fig.show()
    return ax
