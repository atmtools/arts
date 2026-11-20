"""Convergence of the RTE options with the step size.

The same path is computed with ever fewer, longer steps, for three absorption
profiles (constant, decreasing and increasing k).  Relative to its finest
step, the brightness temperature of the reference option ("Linear": constant
or magop) must not fall below that of the linear-in-tau source ("Linear Src")
or of the linear-propagation option ("Linear Prop": linprop or magop_linsrc).
"""

from pyarts3 import Workspace
from pyarts3.arts import PropmatVector, SpectralPropmat, PropagationPathPoint, AtmPoint
import numpy as np
import matplotlib.pyplot as plt
import os

ws = Workspace()
ws.freq_grid = [100e9]
ws.surf_fieldEarth()
ws.surf_field['t'] = 100.0

k_c = 1e-2
PROFILES = {
    "constant k": lambda N: np.full(N, k_c),
    "decreasing k": lambda N: np.linspace(k_c, k_c / 100, N),
    "increasing k": lambda N: np.linspace(k_c / 100, k_c, N),
}

#: The compared options: ("Linear", "Linear Prop")
OPTIONS = [("constant", "linprop"), ("magop", "linprop"), ("constant", "magop_linsrc")]
LINSRC = "lintau"


def convergence(profile, options):
    """Brightness temperature per option, from 2**12 steps of 1 down to 2 steps of 2**11."""
    tb = {option: [] for option in options}
    N = 2**12
    scl = 1
    while N >= 2:
        k = profile(N)
        T = np.linspace(200, 300, N)

        ws.spectral_propmat_path = [SpectralPropmat(PropmatVector([x])) for x in k]
        ws.ray_path = [PropagationPathPoint() for _ in k]
        ws.atm_path = [AtmPoint() for _ in k]
        for i in range(N):
            ws.atm_path[i]['t'] = T[i]
            ws.ray_path[i].pos[0] = i * scl
            ws.ray_path[i].pos_type = 'atm'
            ws.ray_path[i].los_type = 'atm'
        ws.ray_path[-1].los_type = 'surface'
        ws.ray_point = ws.ray_path[-1]

        ws.freq_grid_pathFromPath()
        ws.spectral_rad_bkgAgendasAtEndOfPath()
        ws.spectral_rad_transform_operatorSet(option="Tb")
        for option in options:
            ws.rte_option = option
            ws.spectral_radStepByStepEmissionFromPropmat()
            ws.spectral_radApplyForwardUnit()
            tb[option].append(ws.spectral_rad[0][0] * 1.0)

        N = N // 2
        scl = scl * 2
    return {option: np.array(x) for option, x in tb.items()}


all_options = list(dict.fromkeys([LINSRC, *(o for pair in OPTIONS for o in pair)]))
for title, profile in PROFILES.items():
    tb = convergence(profile, all_options)
    linsrc = tb[LINSRC]
    for lin_option, linprop_option in OPTIONS:
        lin, linprop = tb[lin_option], tb[linprop_option]
        x = 2**-np.linspace(1, 12, len(lin))

        if "ARTS_HEADLESS" not in os.environ:
            plt.semilogx(x, lin / lin[0], '-x', label="Linear")
            plt.semilogx(x, linsrc / linsrc[0], '-x', label="Linear Src")
            plt.semilogx(x, linprop / linprop[0], '-x', label="Linear Prop")

            plt.legend()
            plt.xlabel("Relative step-size [-]")
            plt.ylabel("Relative Brightness Temperature [-]")
            plt.title(f"Ratio, {title} ({lin_option}, {linprop_option})")
            plt.show()

        assert np.all(lin / lin[0] >= linprop / linprop[0])
        assert np.all(lin / lin[0] >= linsrc / linsrc[0])
