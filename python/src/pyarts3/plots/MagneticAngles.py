import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def rot(x, y, ang, clockwise=False):
    """Calculates a rotation matrix to rotate vector x towards vector y."""
    if clockwise:
        if ang >= 0:
            axis = np.cross(y, x)
            ang = -ang
        else:
            axis = np.cross(x, y)
            ang = 2 * np.pi - ang
    else:
        axis = np.cross(x, y)

    if np.linalg.norm(axis) < 1e-9:  # Handle collinear vectors
        return np.diag([1, 1, 1])
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    R = np.diag([1, 1, 1]) + np.sin(ang) * K + (1 - np.cos(ang)) * (K @ K)
    return R


def plot(data: pyarts.arts.zeeman.MagneticAngles,
         *,
         mode='normal',
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         N: int = 50,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plots the magnetic angles in 3D.  The axis should be 3D.

    The N parameter controls the number of points used to draw the angle arcs.

    Parameters
    ----------
    data : ~pyarts3.arts.zeeman.MagneticAngles
        The MagneticAngles object containing the data to plot.
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes], optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    N : int, optional
        Number of points to use for drawing the angle arcs. Default is 50.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """
    fig, ax = default_fig_ax(fig, ax, fig_kwargs={
                             "figsize": (8, 8)}, ax_kwargs={"subplot_kw": {"projection": '3d'}})

    k = data.k.norm()
    e1 = data.e1.norm()

    # Correctly handle vector assignments and normalization
    b_vec = pyarts.arts.Vector3([data.u * 1e6, data.v * 1e6, data.w * 1e6])
    H_vec = b_vec
    B_vec = data.B_projected

    H_norm = H_vec.norm()
    B_norm = B_vec.norm()

    # Plot vectors
    select_flat_ax(ax, 0).quiver(
        0, 0, 0, k[0], k[1], k[2], color="blue", label="k (LOS)", **kwargs)
    select_flat_ax(ax, 0).quiver(
        0, 0, 0, e1[0], e1[1], e1[2], color="green", label="e1 (UP)", **kwargs)
    select_flat_ax(ax, 0).quiver(
        0, 0, 0, H_norm[0], H_norm[1], H_norm[2], color="red", label="H_norm", **kwargs)
    select_flat_ax(ax, 0).quiver(0, 0, 0, B_norm[0], B_norm[1], B_norm[2],
                                 color="orange", label="B_norm (PROJ)", **kwargs)

    rotation_eta = [rot(e1, B_norm, eta, True) for eta in np.linspace(0, data.eta, N)]
    v_eta = np.array([r @ e1 / 4 for r in rotation_eta])
    select_flat_ax(ax, 0).plot(v_eta[:, 0], v_eta[:, 1], v_eta[:, 2],
                               label=f"eta={round(np.rad2deg(data.eta), 1)}°", color="magenta", **kwargs)

    # Plot theta arc (rotation from H_norm to k)
    rotation_theta = [rot(H_norm, k, theta)
                      for theta in np.linspace(0, data.theta, N)]
    v_theta = np.array([r @ H_norm / 2 for r in rotation_theta])
    select_flat_ax(ax, 0).plot(v_theta[:, 0], v_theta[:, 1], v_theta[:, 2],
                               label=f"theta={round(np.rad2deg(data.theta), 1)}°", color="gray", **kwargs)

    return fig, ax
