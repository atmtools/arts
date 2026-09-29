"""Plot O2 60 GHz ECS temperature Jacobians against central differences.

Run this file to display the plot; it also saves ecs_o2_temperature.png.
Uses the same approximate catalogue AIR/Bath model as ecs_o2.py.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

ws = pyarts.Workspace()
ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
bandkey = "O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"
ws.abs_bandsKeepID(id=bandkey)
band = ws.abs_bands[bandkey]
band.lines = [line for line in band.lines if 5e9 < line.f0 < 120e9]
band.lineshape = "VP_ECS_MAKAROV"
ws.WignerInit()
ws.abs_ecs_dataInit()
ws.abs_ecs_dataAddMakarov2020()
ws.abs_ecs_dataAddMeanAir(vmrs=[1], species=["N2"])
ws.abs_bandsCheckEcs()

ws.atm_pointInit()
ws.atm_point.pressure = 1e5
ws.ray_point.los = [180, 0]
ws.atm_point["O2"] = 0.21
ws.atm_point["N2"] = 0.79
ws.freq_grid = np.linspace(50e9, 70e9, 1001)
frequency = np.asarray(ws.freq_grid) / 1e9


def calculate(temperature, jacobian=False):
    ws.atm_point.temperature = temperature
    ws.jac_targetsInit()
    if jacobian:
        ws.jac_targetsAddTemperature()
    ws.spectral_propmatInit()
    ws.spectral_propmatAddLines(no_negative_absorption=False)
    absorption = np.array(ws.spectral_propmat)[:, 0]
    if jacobian:
        return absorption, np.array(ws.spectral_propmat_jac)[0, :, 0]
    return absorption


fig, ax = plt.subplots(3, 1, sharex=True, figsize=(9, 9), layout="constrained")
# Stay inside partition-function interpolation cells when perturbing T.
step = 0.05  # K
errors = []
for temperature in (220.25, 280.25, 330.25):
    absorption, analytic = calculate(temperature, jacobian=True)
    finite_difference = (
        calculate(temperature + step) - calculate(temperature - step)
    ) / (2 * step)
    error = (analytic - finite_difference) / np.max(np.abs(finite_difference))
    errors.append(np.max(np.abs(error)))
    print(f"T={temperature:g} K: max error / peak derivative = {errors[-1]:.3g}")

    color = ax[0].plot(frequency, absorption, label=f"{temperature:g} K")[0].get_color()
    ax[1].plot(frequency, analytic, color=color)
    ax[1].plot(
        frequency, finite_difference, "o", color=color,
        fillstyle="none", markersize=4, markevery=40,
    )
    ax[2].plot(frequency, 100 * error, color=color)

ax[0].set_title("O₂ 60 GHz ECS — temperature check at 1000 hPa")
ax[0].set_ylabel("Absorption [m⁻¹]")
ax[0].legend()
ax[1].set_title("Lines: temperature Jacobian; circles: central difference (±0.05 K)")
ax[1].set_ylabel("∂ absorption / ∂T [m⁻¹ K⁻¹]")
ax[2].set_ylabel("Difference / peak derivative [%]")
ax[2].set_xlabel("Frequency [GHz]")
ax[2].axhline(0, color="0.5", linewidth=0.7)
for axis in ax:
    axis.grid(alpha=0.25)

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

# Allow the ordinary Voigt dF approximation; preserve the plot on failure.
assert np.all(np.isfinite(errors)) and max(errors) < 1e-3
