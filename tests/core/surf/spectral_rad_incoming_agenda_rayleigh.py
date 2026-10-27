"""Tests for the ClearskyRayleighScattering options of spectral_rad_incoming_agenda.

Verifies:
1. option="ClearskyRayleighScattering" equals the manual chain
   ray_path_observer_agendaExecute + ray_path_suns_pathFromPathObserver +
   spectral_radClearskyRayleighScattering
2. option="ClearskyRayleighScatteringOnly" equals the manual chain with
   spectral_radClearskyRayleighScatteringOnly
3. The scattering-only radiance is strictly below the full radiance (the LTE
   emission source is removed) yet still positive (first-order solar Rayleigh
   scattering remains)

The observer sits at the surface looking straight up, so the agenda provides
the radiance incoming at the surface -- the intended use for surface
scattering calculations.
"""

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts


def setup_workspace():
    """Tropical atmosphere with H2O absorption, Earth surface, and a sun at 45 deg elevation."""
    ws = pyarts.Workspace()

    ws.freq_grid = arts.AscendingGrid(np.linspace(22.0e9, 23.0e9, 3))

    # Absorbing atmosphere so that the full chain carries thermal emission
    ws.abs_speciesSet(species=["H2O-161"])
    ws.ReadCatalogData()
    ws.abs_bandsSelectFrequencyByLine(fmin=20e9, fmax=25e9)
    ws.spectral_propmat_agendaAuto()

    ws.atm_fieldRead(toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1)

    ws.surf_fieldEarth()
    ws.surf_field[arts.SurfaceKey("t")] = 288.0

    # Sun 45 degrees up in the sky
    ws.lat = 45.0
    ws.lon = 0.0
    ws.sunBlackbody()
    ws.suns = [ws.sun]

    # The Rayleigh scattering method requires an empty jac_targets
    ws.jac_targets = arts.JacobianTargets()

    ws.spectral_rad_space_agendaSet(option="SunOrCosmicBackground")
    ws.spectral_rad_surface_agendaSet(option="Blackbody")
    ws.ray_path_observer_agendaSetGeometric()
    ws.spectral_propmat_scat_agendaSet(option="AirSimple")

    # Observer at the surface looking straight up
    ws.obs_pos = [0.0, 0.0, 0.0]
    ws.obs_los = [0.0, 0.0]

    return ws


def manual_chain(ws, meta_method):
    """Run the meta-method chain by hand and return the resulting spectral_rad."""
    ws.ray_path_observer_agendaExecute()
    ws.ray_path_suns_pathFromPathObserver(just_hit=1)
    getattr(ws, meta_method)()
    return np.asarray(ws.spectral_rad).copy()


def agenda_option(ws, option):
    """Run the predefined incoming agenda and return the resulting spectral_rad."""
    ws.spectral_rad_incoming_agendaSet(option=option)
    ws.spectral_rad_incoming_agendaExecute()
    return np.asarray(ws.spectral_rad).copy()


# ============================================================================
# Test 1: ClearskyRayleighScattering option equals the manual chain
# ============================================================================

ws = setup_workspace()
rad_full_manual = manual_chain(ws, "spectral_radClearskyRayleighScattering")
rad_full_agenda = agenda_option(ws, "ClearskyRayleighScattering")

assert rad_full_manual.shape == (3, 4)
assert np.all(np.isfinite(rad_full_manual))
np.testing.assert_allclose(rad_full_agenda, rad_full_manual, rtol=1e-12)
print("Test 1 passed: ClearskyRayleighScattering agenda option equals the manual chain")

# ============================================================================
# Test 2: ClearskyRayleighScatteringOnly option equals the manual chain
# ============================================================================

rad_only_manual = manual_chain(ws, "spectral_radClearskyRayleighScatteringOnly")
rad_only_agenda = agenda_option(ws, "ClearskyRayleighScatteringOnly")

assert rad_only_manual.shape == (3, 4)
assert np.all(np.isfinite(rad_only_manual))
np.testing.assert_allclose(rad_only_agenda, rad_only_manual, rtol=1e-12)
print("Test 2 passed: ClearskyRayleighScatteringOnly agenda option equals the manual chain")

# ============================================================================
# Test 3: Scattering-only radiance is positive but below the full radiance
# ============================================================================

assert np.all(rad_only_agenda[:, 0] > 0.0), f"Scattering-only radiance not positive:\n{rad_only_agenda}"
assert np.all(rad_only_agenda[:, 0] < rad_full_agenda[:, 0]), (
    f"Removing thermal emission did not lower the radiance:\n"
    f"  only = {rad_only_agenda[:, 0]}\n"
    f"  full = {rad_full_agenda[:, 0]}"
)
print("Test 3 passed: scattering-only radiance is positive and below the full radiance")
