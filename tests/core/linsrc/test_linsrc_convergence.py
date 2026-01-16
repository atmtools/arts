from pyarts3 import Workspace
from pyarts3.arts import PropmatVector, PropmatMatrix, PropagationPathPoint, StokvecVector, StokvecMatrix, AtmPoint
import numpy as np
import matplotlib.pyplot as plt
import os

ws = Workspace()

N = 2**12
scl = 1

Ts = 100.0
k_c = 1e-2
linprop = []
linsrc = []
lin = []


while N >= 2:
    k = np.ones((N)) * k_c
    T = np.linspace(200, 300, N)

    ws.freq_grid = [100e9]

    ws.spectral_propmat_path = [PropmatVector([x]) for x in k]

    ws.spectral_propmat_jac_path = [PropmatMatrix(np.zeros((0, 1))) for _ in k]

    ws.spectral_nlte_srcvec_path = [StokvecVector([0.0]) for _ in k]

    ws.spectral_nlte_srcvec_jac_path = [StokvecMatrix(np.zeros((0, 1))) for _ in k]

    ws.ray_path = [PropagationPathPoint() for _ in k]

    ws.surf_fieldEarth()
    ws.surf_field['t'] = Ts

    ws.atm_path = [AtmPoint() for _ in k]
    for i in range(N):
        ws.atm_path[i]['t'] = T[i]
        ws.ray_path[i].pos[0] = i * scl
        ws.ray_path[i].pos_type = 'atm'
        ws.ray_path[i].los_type = 'atm'
    ws.ray_path[-1].los_type = 'surface'
    ws.ray_point = ws.ray_path[-1]

    ws.rte_option = "lintau"
    ws.spectral_tramat_pathFromPath()
    ws.freq_grid_pathFromPath()
    ws.spectral_rad_srcvec_pathFromPropmat()
    ws.spectral_rad_bkgAgendasAtEndOfPath()

    ws.spectral_radStepByStepEmission()
    ws.spectral_rad_transform_operatorSet(option="Tb")
    ws.spectral_radApplyForwardUnit()
    linsrc.append(ws.spectral_rad[0][0] * 1.0)

    ws.rte_option = "constant"
    ws.spectral_tramat_pathFromPath()
    ws.spectral_radStepByStepEmission()
    ws.spectral_rad_transform_operatorSet(option="Tb")
    ws.spectral_radApplyForwardUnit()
    lin.append(ws.spectral_rad[0][0] * 1.0)

    ws.rte_option = "linprop"
    ws.spectral_tramat_pathFromPath()
    ws.freq_grid_pathFromPath()
    ws.spectral_rad_srcvec_pathFromPropmat()
    ws.spectral_rad_bkgAgendasAtEndOfPath()

    ws.spectral_radStepByStepEmission()
    ws.spectral_rad_transform_operatorSet(option="Tb")
    ws.spectral_radApplyForwardUnit()
    linprop.append(ws.spectral_rad[0][0] * 1.0)

    N = N // 2
    scl = scl * 2

linprop = np.array(linprop)
linsrc = np.array(linsrc)
lin = np.array(lin)
x = 2**-np.linspace(1, 12, len(lin))

if "ARTS_HEADLESS" not in os.environ:
    plt.semilogx(x, lin / lin[0], '-x', label="Linear")
    plt.semilogx(x, linsrc / linsrc[0], '-x', label="Linear Src")
    plt.semilogx(x, linprop / linprop[0], '-x', label="Linear Prop")

    plt.legend()
    plt.xlabel("Relative step-size [-]")
    plt.ylabel("Relative Brightness Temperature [-]")
    plt.title("Ratio, constant k")
    plt.show()

assert np.all(lin / lin[0] >= linprop / linprop[0])
assert np.all(lin / lin[0] >= linsrc / linsrc[0])

linprop = []
linsrc = []
lin = []
N = 2**12
scl = 1


while N >= 2:
    k = np.linspace(k_c, k_c / 100, N)
    T = np.linspace(200, 300, N)

    ws.freq_grid = [100e9]

    ws.spectral_propmat_path = [PropmatVector([x]) for x in k]

    ws.spectral_propmat_jac_path = [PropmatMatrix(np.zeros((0, 1))) for _ in range(N)]

    ws.spectral_nlte_srcvec_path = [StokvecVector([0.0]) for _ in range(N)]

    ws.spectral_nlte_srcvec_jac_path = [StokvecMatrix(np.zeros((0, 1))) for _ in range(N)]

    ws.ray_path = [PropagationPathPoint() for _ in range(N)]

    ws.surf_fieldEarth()
    ws.surf_field['t'] = Ts

    ws.atm_path = [AtmPoint() for _ in range(N)]
    for i in range(N):
        ws.atm_path[i]['t'] = T[i]
        ws.ray_path[i].pos[0] = i * scl
        ws.ray_path[i].pos_type = 'atm'
        ws.ray_path[i].los_type = 'atm'
    ws.ray_path[-1].los_type = 'surface'
    ws.ray_point = ws.ray_path[-1]

    ws.rte_option = "lintau"
    ws.spectral_tramat_pathFromPath()
    ws.freq_grid_pathFromPath()
    ws.spectral_rad_srcvec_pathFromPropmat()
    ws.spectral_rad_bkgAgendasAtEndOfPath()

    ws.spectral_radStepByStepEmission()
    ws.spectral_rad_transform_operatorSet(option="Tb")
    ws.spectral_radApplyForwardUnit()
    linsrc.append(ws.spectral_rad[0][0] * 1.0)

    ws.rte_option = "constant"
    ws.spectral_tramat_pathFromPath()
    ws.spectral_radStepByStepEmission()
    ws.spectral_rad_transform_operatorSet(option="Tb")
    ws.spectral_radApplyForwardUnit()
    lin.append(ws.spectral_rad[0][0] * 1.0)

    ws.rte_option = "linprop"
    ws.spectral_tramat_pathFromPath()
    ws.freq_grid_pathFromPath()
    ws.spectral_rad_srcvec_pathFromPropmat()
    ws.spectral_rad_bkgAgendasAtEndOfPath()

    ws.spectral_radStepByStepEmission()
    ws.spectral_rad_transform_operatorSet(option="Tb")
    ws.spectral_radApplyForwardUnit()
    linprop.append(ws.spectral_rad[0][0] * 1.0)

    N = N // 2
    scl = scl * 2

linprop = np.array(linprop)
linsrc = np.array(linsrc)
lin = np.array(lin)
x = 2**-np.linspace(1, 12, len(lin))

if "ARTS_HEADLESS" not in os.environ:
    plt.semilogx(x, lin / lin[0], '-x', label="Linear")
    plt.semilogx(x, linsrc / linsrc[0], '-x', label="Linear Src")
    plt.semilogx(x, linprop / linprop[0], '-x', label="Linear Prop")

    plt.legend()
    plt.xlabel("Relative step-size [-]")
    plt.ylabel("Relative Brightness Temperature [-]")
    plt.title("Ratio, decreasing k")
    plt.show()

assert np.all(lin / lin[0] >= linprop / linprop[0])
assert np.all(lin / lin[0] >= linsrc / linsrc[0])


linprop = []
linsrc = []
lin = []
N = 2**12
scl = 1


while N >= 2:
    k = np.linspace(k_c / 100, k_c, N)
    T = np.linspace(200, 300, N)

    ws.freq_grid = [100e9]

    ws.spectral_propmat_path = [PropmatVector([x]) for x in k]

    ws.spectral_propmat_jac_path = [PropmatMatrix(np.zeros((0, 1))) for _ in k]

    ws.spectral_nlte_srcvec_path = [StokvecVector([0.0]) for _ in k]

    ws.spectral_nlte_srcvec_jac_path = [StokvecMatrix(np.zeros((0, 1))) for _ in k]

    ws.ray_path = [PropagationPathPoint() for _ in k]

    ws.surf_fieldEarth()
    ws.surf_field['t'] = Ts

    ws.atm_path = [AtmPoint() for _ in k]
    for i in range(N):
        ws.atm_path[i]['t'] = T[i]
        ws.ray_path[i].pos[0] = i * scl
        ws.ray_path[i].pos_type = 'atm'
        ws.ray_path[i].los_type = 'atm'
    ws.ray_path[-1].los_type = 'surface'
    ws.ray_point = ws.ray_path[-1]

    
    ws.rte_option = "lintau"
    ws.spectral_tramat_pathFromPath()
    ws.freq_grid_pathFromPath()
    ws.spectral_rad_srcvec_pathFromPropmat()
    ws.spectral_rad_bkgAgendasAtEndOfPath()

    ws.spectral_radStepByStepEmission()
    ws.spectral_rad_transform_operatorSet(option="Tb")
    ws.spectral_radApplyForwardUnit()
    linsrc.append(ws.spectral_rad[0][0] * 1.0)

    ws.rte_option = "constant"
    ws.spectral_tramat_pathFromPath()
    ws.spectral_radStepByStepEmission()
    ws.spectral_rad_transform_operatorSet(option="Tb")
    ws.spectral_radApplyForwardUnit()
    lin.append(ws.spectral_rad[0][0] * 1.0)

    ws.rte_option = "linprop"
    ws.spectral_tramat_pathFromPath()
    ws.freq_grid_pathFromPath()
    ws.spectral_rad_srcvec_pathFromPropmat()
    ws.spectral_rad_bkgAgendasAtEndOfPath()

    ws.spectral_radStepByStepEmission()
    ws.spectral_rad_transform_operatorSet(option="Tb")
    ws.spectral_radApplyForwardUnit()
    linprop.append(ws.spectral_rad[0][0] * 1.0)

    N = N // 2
    scl = scl * 2

linprop = np.array(linprop)
linsrc = np.array(linsrc)
lin = np.array(lin)
x = 2**-np.linspace(1, 12, len(lin))

if "ARTS_HEADLESS" not in os.environ:
    plt.semilogx(x, lin / lin[0], '-x', label="Linear")
    plt.semilogx(x, linsrc / linsrc[0], '-x', label="Linear Src")
    plt.semilogx(x, linprop / linprop[0], '-x', label="Linear Prop")

    plt.legend()
    plt.xlabel("Relative step-size [-]")
    plt.ylabel("Relative Brightness Temperature [-]")
    plt.title("Ratio, increasing k")
    plt.show()

assert np.all(lin / lin[0] >= linprop / linprop[0])
assert np.all(lin / lin[0] >= linsrc / linsrc[0])