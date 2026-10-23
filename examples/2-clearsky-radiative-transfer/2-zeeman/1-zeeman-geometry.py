import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

# %% Get magnetic field at a sample position
ell = pyarts.arts.planets.Earth.ellipsoid
pos = [50e3, 0, 0]
mag = pyarts.arts.igrf(pos, ell, t="2000-03-11 14:39:37")

# %% Setup figure for multiple subplots
N = 3
M = 4
fig = plt.figure(figsize=(M * 8, N * 8))

# %% Store computed angles for later verification
angles = []

# %% loop over different LOS directions
for i in range(N):
    for j in range(M):
        los = [np.linspace(0, 180, N)[i], np.linspace(0, 360, M)[j]]

        # Setup 3D subplot
        ax = fig.add_subplot(N, M, i*M + j + 1, projection='3d')

        # Plot and store angles
        ang = pyarts.arts.zeeman.MagneticAngles(mag, los)
        pyarts.plots.MagneticAngles.plot(ang, fig=fig, ax=ax, N=50)
        ax.set_ylim(-1.2, 1.2)
        ax.set_xlim(-1.2, 1.2)
        ax.set_zlim(-1.2, 1.2)
        angles.append([ang.eta, ang.theta])

if "ARTS_HEADLESS" not in os.environ:
    plt.tight_layout()
    plt.show()

assert np.allclose(angles,
                   np.array([[3.01433146,  1.08543309],
                             [0.91993636,  1.08543309],
                             [-1.17445874,  1.08543309],
                             [3.01433146,  1.08543309],
                             [0.23613776,  0.50047915],
                             [0.98534648,  2.13636506],
                             [-1.05139296,  1.91926261],
                             [0.23613776,  0.50047915],
                             [0.12726119,  2.05615956],
                             [2.22165629,  2.05615956],
                             [-1.96713391,  2.05615956],
                             [0.12726119,  2.05615956]])), "Angles do not match expected values."
