from pyarts3.arts.zeeman import MagneticAngles
from pyarts3.arts import Vector3
import matplotlib.pyplot as plt
import numpy as np

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


def plot(ang: MagneticAngles,
         *,
         mode='normal',
         fig=None,
         ax=None,
         N: int = 50,
         **kwargs):
    """Plots the magnetic angles in 3D.  The axis should be 3D.

    The N parameter controls the number of points used to draw the angle arcs.

    Parameters
    ----------
    ang : ~pyarts3.arts.zeeman.MagneticAngles
        The MagneticAngles object containing the data to plot.
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    N : int, optional
        Number of points to use for drawing the angle arcs. Default is 50.

    Returns
    -------
    fig : As input
        As input.
    ax : As input
        As input.
    """
    if fig is None:
        fig = plt.figure(figsize=(8, 8))

    if ax is None:
        ax = fig.add_subplot(1, 1, 1, projection='3d')

    k = ang.k.norm()
    e1 = ang.e1.norm()

    # Correctly handle vector assignments and normalization
    b_vec = Vector3([ang.u * 1e6, ang.v * 1e6, ang.w * 1e6])
    H_vec = b_vec
    B_vec = ang.B_projected

    H_norm = H_vec.norm()
    B_norm = B_vec.norm()

    # Setup 3D subplot
    azimuth_angle_deg = np.rad2deg(np.arctan2(ang.sa, ang.ca)) % 360
    ax.set_title(
        f"za={round(np.rad2deg(np.arccos(ang.cz)), 1)}°, aa={round(azimuth_angle_deg, 1)}°")
    ax.set_xlim([-1.1, 1.1])
    ax.set_ylim([-1.1, 1.1])
    ax.set_zlim([-1.1, 1.1])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    # Plot vectors
    ax.quiver(0, 0, 0, k[0], k[1], k[2], color="blue", label="k (LOS)")
    ax.quiver(0, 0, 0, e1[0], e1[1], e1[2], color="green", label="e1 (UP)")
    ax.quiver(0, 0, 0, H_norm[0], H_norm[1], H_norm[2], color="red", label="H_norm")
    ax.quiver(0, 0, 0, B_norm[0], B_norm[1], B_norm[2],
              color="orange", label="B_norm (PROJ)")

    rotation_eta = [rot(e1, B_norm, eta, True) for eta in np.linspace(0, ang.eta, N)]
    v_eta = np.array([r @ e1 / 4 for r in rotation_eta])
    ax.plot(v_eta[:, 0], v_eta[:, 1], v_eta[:, 2],
            label=f"eta={round(np.rad2deg(ang.eta), 1)}°", color="magenta")

    # Plot theta arc (rotation from H_norm to k)
    rotation_theta = [rot(H_norm, k, theta)
                      for theta in np.linspace(0, ang.theta, N)]
    v_theta = np.array([r @ H_norm / 2 for r in rotation_theta])
    ax.plot(v_theta[:, 0], v_theta[:, 1], v_theta[:, 2],
            label=f"theta={round(np.rad2deg(ang.theta), 1)}°", color="gray")
    ax.legend(ncol=1)

    return fig, ax
