"""
Simulate an AMSU-A-like microwave sensor loaded from file.

The sensor in amsu-a.xml holds the 15 channels: their pass bands and their
polarization.  AMSU-A scans cross-track with a rotating reflector, so the
polarization at nadir turns with the scan angle.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pa

# fmt: off
POS        = [817e3, 0.0, 0.0]
HEADING    = -90.0  # Flying west
SCAN_ANGLE = 50.0   # Positive is to the right of the flight direction
LOS        = [180.0 - SCAN_ANGLE, HEADING + 90.0]
ELL        = pa.arts.planets.Earth.ellipsoid
# fmt: on


def extract_channels(spectral_rad, measurement_sensor):
    """ Helper to extract internal simulation grid and Stokes I of the channels """
    out = []

    for ch in measurement_sensor:
        v = np.nonzero(np.asarray(ch.weight_matrix.dense)[0, :, 0])[0]
        out.append([np.asarray(ch.f_grid)[v], np.asarray(spectral_rad)[v, 0]])

    return out


# %% Download data and set up workspace
ws = pa.Workspace()
pa.data.download()

# %% Use built-in spectroscopic data.
ws.abs_speciesSet(species=["O2-PWR98", "H2O-PWR98"])
ws.ReadCatalogData()
ws.spectral_propmat_agendaAuto()

# %% Set up a simple atmosphere and surface.
ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pa.arts.SurfaceKey("t")] = 300.0
ws.atm_fieldRead(toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1)

# %% Use a geometric path.
ws.ray_path_observer_agendaSetGeometric()

# %% Load the sensor and place it in a cross-track scan geometry.
sensor = pa.arts.SensorBuilder.fromxml("amsu-a.xml")
ws.measurement_sensor, ws.measurement_sensor_meta = sensor.cross_track([POS], [HEADING], [SCAN_ANGLE], ELL)

# %% Collect on a single grid to help plotting and speed-up calculations - this sacrifices Jacobian flexibilities
ws.measurement_sensor.collect()

# %% Compute the sensor response in Planck units
ws.spectral_rad_transform_operatorSet(option="Tb")
ws.measurement_vecFromSensor(kernel="High Performance")

# %% Compute the internal frequency grid response for plotting
ws.freq_grid = ws.measurement_sensor[0].f_grid
ws.spectral_rad_observer_agendaExecute(obs_pos=POS, obs_los=LOS)
ws.spectral_radApplyUnitFromSpectralRadiance()

# %% Plotting
ch = extract_channels(ws.spectral_rad, ws.measurement_sensor)
numbers = np.arange(1, len(ch) + 1)
colors = plt.get_cmap('tab20')(range(len(ch)))
fig_internal, ax_internal = plt.subplots(figsize=(6, 4))
for i in range(len(ch)):
    ax_internal.plot(ch[i][0]/1e9,
                     ch[i][1],
                     "o:",
                     markersize=2.5,
                     linewidth=1.1,
                     label=f"Ch {numbers[i]}",
                     color=colors[i],
                     )

ax_internal.set_xlabel("Frequency [GHz]")
ax_internal.set_ylabel("Brightness temperature (internal grid) [K]")
ax_internal.set_title("AMSU-A-like channels")
ax_internal.grid(True, alpha=0.3)
ax_internal.legend(ncol=3, fontsize=8)

fig_channels, ax_channels = plt.subplots(figsize=(6, 4))
ax_channels.plot(numbers, ws.measurement_vec, "o")
ax_channels.set_xticks(numbers)
ax_channels.set_xlabel("AMSU-A channel number")
ax_channels.set_ylabel("Channel brightness [K]")
ax_channels.set_title("AMSU-A-like channelized measurements")
ax_channels.grid(True, alpha=0.3)

fig_internal.tight_layout()
fig_channels.tight_layout()

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

ref = [294.8203165065995, 297.1947667477886, 283.30250470594353, 262.1501094031428, 244.5833750280594, 226.45647530368532, 216.194523783748,
       209.69385112980513, 210.0837314485364, 219.26382030306507, 230.25211121603354, 241.29466862910792, 252.4615035571569, 261.42686128739854,
       292.2250588989928]

assert np.allclose(ws.measurement_vec, ref), \
    f"Mismatch to reference simulations.\nReference:   {ref}\nSimulations: {ws.measurement_vec:B,}"
