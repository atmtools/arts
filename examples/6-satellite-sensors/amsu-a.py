"""
Build an AMSU-A-like microwave sensor from a simple channel sheet.

This example is intentionally script-like and intentionally approximate:

- it keeps only O2-PWR98 and H2O-PWR98 absorption
- it uses flat-top lobes instead of a fully measured backend response
- it treats the listed bandwidths as total bandwidth, so each lobe spans
  ``center +/- BW/2``
- it uses one nadir view instead of a full cross-track scan

The point is to show how a spec-sheet style channel list can be turned into
ARTS sensor elements, and how the LO offsets create the realized RF lobes.
"""

from dataclasses import dataclass
import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pa


F0_O2 = 57.290344e9
SAMPLES_PER_LOBE = 11


@dataclass(frozen=True)
class ChannelSpec:
    number: int
    spec_text: str
    reference_hz: float
    bandwidth_hz: float
    polarization: str
    nedt_k: float
    lo_offsets_hz: tuple[float, ...] = ()

# fmt: off
CHANNELS = [
    ChannelSpec(1, "23.800",                23.800e9,  270e6, "QV", 0.30),
    ChannelSpec(2, "31.400",                31.400e9,  180e6, "QV", 0.30),
    ChannelSpec(3, "50.300",                50.300e9,  180e6, "QV", 0.40),
    ChannelSpec(4, "52.800",                52.800e9,  400e6, "QV", 0.25),
    ChannelSpec(5, "53.596 ± 0.115",        53.596e9,  170e6, "QH", 0.25, (0.115e9,)),
    ChannelSpec(6, "54.400",                54.400e9,  400e6, "QH", 0.25),
    ChannelSpec(7, "54.940",                54.940e9,  400e6, "QV", 0.25),
    ChannelSpec(8, "55.500",                55.500e9,  330e6, "QH", 0.25),
    ChannelSpec(9, "f0 = 57.290344",        F0_O2,     330e6, "QH", 0.25),
    ChannelSpec(10, "f0 ± 0.217",           F0_O2,      78e6, "QH", 0.40, (0.217e9,)),
    ChannelSpec(11, "f0 ± 0.3222 ± 0.048",  F0_O2,       6e6, "QH", 0.40, (0.3222e9, 0.048e9)),
    ChannelSpec(12, "f0 ± 0.3222 ± 0.022",  F0_O2,      16e6, "QH", 0.60, (0.3222e9, 0.022e9)),
    ChannelSpec(13, "f0 ± 0.3222 ± 0.010",  F0_O2,       8e6, "QH", 0.80, (0.3222e9, 0.010e9)),
    ChannelSpec(14, "f0 ± 0.3222 ± 0.0045", F0_O2,       3e6, "QH", 1.20, (0.3222e9, 0.0045e9)),
    ChannelSpec(15, "89.000 ± 1.0",         89.000e9, 1000e6, "QV", 0.50, (1.0e9,)),
]
# fmt: on


def split_centers(reference_hz: float, offsets_hz: tuple[float, ...]) -> np.ndarray:
    centers = np.array([reference_hz], dtype=float)

    for offset in offsets_hz:
        centers = np.concatenate([centers - offset, centers + offset])

    return np.sort(centers)


def realized_channel(spec: ChannelSpec, samples_per_lobe: int = SAMPLES_PER_LOBE):
    half_bw = 0.5 * spec.bandwidth_hz
    local_df = np.linspace(-half_bw, half_bw, samples_per_lobe)
    lobe_centers = split_centers(spec.reference_hz, spec.lo_offsets_hz)

    freq_grid = np.concatenate([center + local_df for center in lobe_centers])
    weights = np.ones(freq_grid.size, dtype=float)

    order = np.argsort(freq_grid)
    freq_grid = freq_grid[order]
    weights = weights[order]
    weights /= weights.sum()

    return {
        "reference_hz": spec.reference_hz,
        "lobe_centers_hz": lobe_centers,
        "freq_grid_hz": freq_grid,
        "weights": weights,
        "dfreq_hz": freq_grid - spec.reference_hz,
    }


def make_raw_sensor(spec: ChannelSpec, channel_data: dict):
    raw = pa.arts.SortedGriddedField1()
    raw.dataname = f"amsu-ch{spec.number}"
    raw.gridnames = ["dfreq"]
    raw.grids = [pa.arts.AscendingGrid(channel_data["dfreq_hz"])]
    raw.data = channel_data["weights"]
    return raw


def compute_channel_spectral_rad(ws, pos, los, freq_grid_hz: np.ndarray) -> np.ndarray:
    ws.freq_grid = pa.arts.AscendingGrid(freq_grid_hz)
    ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)
    ws.spectral_radClearskyEmission()
    ws.spectral_radApplyUnitFromSpectralRadiance()
    return np.asarray(ws.spectral_rad.value, dtype=float)[:, 0].copy()


pa.data.download()

ws = pa.workspace.Workspace()

# Keep the atmosphere intentionally simple.
ws.abs_speciesSet(species=["O2-PWR98", "H2O-PWR98"])
ws.ReadCatalogData()
ws.spectral_propmat_agendaAuto()

ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pa.arts.SurfaceKey("t")] = 290.0
ws.atm_fieldRead(toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1)

ws.spectral_rad_transform_operatorSet(option="Tb")
ws.ray_path_observer_agendaSetGeometric()

# AMSU-A is a cross-track scanner, we take a single nadir view for simplicity.
pos = pa.arts.Vector3([833e3, 0.0, 0.0])
los = pa.arts.Vector2([180.0, 0.0])

realized = {spec.number: realized_channel(spec) for spec in CHANNELS}

# First plot: evaluate spectral radiance directly on each channel's own internal grid.
channel_spectra = {}
manual_measurement = []
for spec in CHANNELS:
    channel_data = realized[spec.number]
    spectral_tb = compute_channel_spectral_rad(
        ws, pos, los, channel_data["freq_grid_hz"])
    channel_spectra[spec.number] = spectral_tb
    manual_measurement.append(np.dot(channel_data["weights"], spectral_tb))

manual_measurement = np.asarray(manual_measurement)

# Second step: build the same channel set as an ARTS measurement sensor.
ws.measurement_sensorInit()
for spec in CHANNELS:
    ws.measurement_sensorAddRawSensor(
        freq_grid=pa.arts.AscendingGrid(np.array([spec.reference_hz])),
        pos=pos,
        los=los,
        raw_sensor_perturbation=make_raw_sensor(spec, realized[spec.number]),
        normalize=1,
    )
ws.measurement_sensor.normalize()
ws.measurement_sensor.collect()

ws.measurement_vecFromSensor(kernel="High Performance")
measurement_vec = np.asarray(ws.measurement_vec)

colors = plt.get_cmap('tab20')(range(15))
fig_internal, ax_internal = plt.subplots(figsize=(6, 4))
for spec in CHANNELS:
    channel_data = realized[spec.number]
    ax_internal.plot(
        channel_data["freq_grid_hz"] / 1e9,
        channel_spectra[spec.number],
        marker="o",
        markersize=2.5,
        linewidth=1.1,
        label=f"Ch {spec.number}",
        color=colors[spec.number - 1],
    )

ax_internal.set_xlabel("Frequency [GHz]")
ax_internal.set_ylabel("Brightness temperature (internal grid) [K]")
ax_internal.set_title("AMSU-A-like channels")
ax_internal.grid(True, alpha=0.3)
ax_internal.legend(ncol=3, fontsize=8)

fig_channels, ax_channels = plt.subplots(figsize=(6, 4))
ax_channels.errorbar(
    [spec.number for spec in CHANNELS],
    measurement_vec,
    yerr=[spec.nedt_k for spec in CHANNELS],
    fmt="o",
    capsize=3,
    linewidth=1.0,
)
ax_channels.plot([spec.number for spec in CHANNELS],
                 measurement_vec, linewidth=1.0, alpha=0.8)
ax_channels.set_xticks([spec.number for spec in CHANNELS])
ax_channels.set_xlabel("AMSU-A channel number")
ax_channels.set_ylabel("Channel brightness [K]")
ax_channels.set_title("AMSU-A-like channelized measurements")
ax_channels.grid(True, alpha=0.3)

for x, spec, y in zip([spec.number for spec in CHANNELS], CHANNELS, measurement_vec):
    ax_channels.annotate(
        spec.spec_text,
        (x, y),
        textcoords="offset points",
        xytext=(0, 7),
        ha="center",
        fontsize=7,
        rotation=90,
    )

fig_internal.tight_layout()
fig_channels.tight_layout()

if "ARTS_HEADLESS" not in os.environ:
    plt.show()


assert np.allclose(ws.measurement_vec, [289.290378, 289.54484543, 284.35054685, 273.08020996,
                                        259.90747748, 241.3598937, 228.07526569, 217.2986938,
                                        207.60407201, 213.70169277, 223.87462958, 235.45802133,
                                        246.81775696, 257.07536243, 289.1540646])
