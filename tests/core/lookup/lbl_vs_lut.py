import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt
from time import time


toa = 100e3


def calc_table(ws):
    ws.abs_lookup_dataCalc(
        lat=0.0,
        lon=0.0,
        alt_grid=np.linspace(0, toa, 101),
        temperature_perturbation=np.linspace(-30, 30, 5),
        water_perturbation=np.logspace(-1, 1, 5),
        water_affected_species=["H2O"],
    )


def profiles_table(ws):
    ws.abs_lookup_dataFromProfiles(
        pressure_profile=np.array(ws.atm_field["p"].data.flatten()),
        temperature_profile=np.array(ws.atm_field["t"].data.flatten()),
        vmr_profiles={"H2O": np.array(ws.atm_field["H2O"].data.flatten())},
        temperature_perturbation=np.linspace(-30, 30, 9),
        water_perturbation=np.logspace(-1, 1, 9),
        water_affected_species=["H2O"],
    )


def wide_table(ws):
    ws.abs_lookup_dataSimpleWide(water_affected_species=["H2O"],
                                 pressure_range=[1e-2, 1100e2])


def compare(name, hitran_s, freq_grid, build_table, temperature_offsets, interp):
    print(name)

    # %% Setup workspace

    ws = pyarts.Workspace()

    ws.abs_speciesSet(species=["H2O-161"])

    ws.ReadCatalogData()
    ws.abs_bandsSetCutoff(option="ByLine", value=750e9)

    ws.abs_bands.keep_hitran_s(hitran_s)

    ws.surf_fieldPlanet(option="Earth")
    ws.surf_field["t"] = 295.0

    ws.atm_fieldRead(
        toa=toa, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
    )
    tdata = pyarts.arts.Tensor3(ws.atm_field["t"].data.data)
    wdata = pyarts.arts.Tensor3(ws.atm_field["H2O"].data.data)

    ws.freq_grid = freq_grid

    t = time()
    build_table(ws)
    print(round(1000 * (time() - t)), "ms to train the LUT")

    ws.spectral_rad_transform_operatorSet(option="Tb")
    ws.spectral_rad_space_agendaSet(option="UniformCosmicBackground")
    ws.spectral_rad_surface_agendaSet(option="Blackbody")

    pos = [100e3, 0, 0]
    los = [180.0, 0.0]
    ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)

    # %% Checks and settings for LBL
    water_ratios = [5e-1, 5]
    fig, axes = plt.subplots(1, len(water_ratios), figsize=(12, 4), sharey=True)
    fig.suptitle(name)
    for ax, water_ratio in zip(axes, water_ratios):
        for temperature_offset in temperature_offsets:
            print(water_ratio, temperature_offset)
            ws.atm_field["t"].data.data = tdata + temperature_offset
            ws.atm_field["H2O"].data.data = wdata * water_ratio

            t = time()
            ws.spectral_propmat_agendaAuto(use_abs_lookup_data=0)
            ws.spectral_radClearskyEmission()
            ws.spectral_radApplyUnitFromSpectralRadiance()
            lbl = ws.spectral_rad[:, 0] * 1.0
            print(round(1000 * (time() - t)), "ms to compute the LBL spectral radiance")

            t = time()
            ws.spectral_propmat_agendaAuto(
                f_interp_order=0,
                p_interp_order=5,
                use_abs_lookup_data=1,
                **interp,
            )
            ws.spectral_radClearskyEmission()
            ws.spectral_radApplyUnitFromSpectralRadiance()
            lut = ws.spectral_rad[:, 0] * 1.0
            print(round(1000 * (time() - t)), "ms to compute the LUT spectral radiance")

            ax.plot(
                pyarts.arts.convert.freq2kaycm(ws.freq_grid),
                lbl - lut,
                label=f"{temperature_offset} K",
            )

            assert np.allclose(lbl, lut, atol=1e-3), (name, water_ratio, temperature_offset)

        ax.legend(title="Temperature offset")
        ax.set_xlabel("Frequency [cm$^{-1}$]")
        ax.set_title(f"Water ratio: {100*water_ratio}%")
    axes[0].set_ylabel("LBL - LUT radiance difference [K]")


compare(
    "Lookup table along a path",
    70,
    pyarts.arts.AscendingGrid.fromxml(
        "planets/Earth/Optimized-Flux-Frequencies/LW-flux-optimized-f_grid.xml"
    ),
    calc_table,
    np.linspace(-20, 20, 5),
    {"t_interp_order": 4, "water_interp_order": 4},
)

compare(
    "Lookup table from profiles",
    70,
    pyarts.arts.convert.kaycm2freq(np.linspace(400, 2500, 101)),
    profiles_table,
    np.linspace(-20, 20, 3),
    {"t_interp_order": 4, "water_interp_order": 4},
)

compare(
    "Wide lookup table",
    90,
    pyarts.arts.convert.kaycm2freq(np.linspace(400, 2500, 21)),
    wide_table,
    np.linspace(-20, 20, 5),
    {"t_interp_order": 7},
)
