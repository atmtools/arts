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
        angles.append([ang.eta, ang.theta])

fig.tight_layout()

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

assert np.allclose(angles,
                   np.array([[3.02347622,  1.07933684],
                             [0.92908112,  1.07933684],
                             [-1.16531399,  1.07933684],
                             [3.02347622,  1.07933684],
                             [0.21669967,  0.50432248],
                             [0.98174228,  2.12671818],
                             [-1.04334554,  1.92599378],
                             [0.21669967,  0.50432248],
                             [0.11811643,  2.06225581],
                             [2.21251154,  2.06225581],
                             [-1.97627867,  2.06225581],
                             [0.11811643,  2.06225581]])), "Angles do not match expected values."
