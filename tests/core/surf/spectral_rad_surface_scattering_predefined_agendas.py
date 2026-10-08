"""Tests for the predefined SurfaceScatteringModel* options of *spectral_rad_surface_agenda*.

The options chain the *spectral_radSurfaceScattering* methods:

* SurfaceScatteringModel            -> Init + Diffuse + Specular + DiffuseDirect + SpecularDirect
* SurfaceScatteringModelDiffuseOnly -> Init + Diffuse + Specular
* SurfaceScatteringModelDirectOnly  -> Init + DiffuseDirect + SpecularDirect

Verifies:
1. Each option builds with the exact expected method chain
2. Each option executes in a workspace with a blended model, a sun, and a quadrature point
3. radiance(DiffuseOnly) + radiance(DirectOnly) == radiance(SurfaceScatteringModel),
   and the same for spectral_rad_jac -- with the sun placed outside every traced
   direction, because the full option gates sun-containing directions (exclude_suns)
4. With no suns the DirectOnly option returns only surface emission
"""

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts

from spectral_rad_surface_scattering_chain import (
    setup_workspace_base,
    blended_models,
    add_surface_mask,
    make_sun,
    stokes_array,
    jac_array,
)

INIT = "spectral_radSurfaceScatteringInit"
DIFFUSE = "spectral_radSurfaceScatteringDiffuse"
SPECULAR = "spectral_radSurfaceScatteringSpecular"
DIFFUSE_DIRECT = "spectral_radSurfaceScatteringDiffuseDirect"
SPECULAR_DIRECT = "spectral_radSurfaceScatteringSpecularDirect"

EXPECTED_CHAINS = {
    "SurfaceScatteringModel": [INIT, DIFFUSE, SPECULAR, DIFFUSE_DIRECT, SPECULAR_DIRECT],
    "SurfaceScatteringModelDiffuseOnly": [INIT, DIFFUSE, SPECULAR],
    "SurfaceScatteringModelDirectOnly": [INIT, DIFFUSE_DIRECT, SPECULAR_DIRECT],
}


def build_workspace(freq_grid, option, suns, with_jac=False):
    """Workspace with blended surface models and the predefined agenda option set."""
    ws = setup_workspace_base(freq_grid, suns=suns)
    add_surface_mask(ws, "lambertian")
    add_surface_mask(ws, "flat_scalar")
    ws.surface_models = blended_models(freq_grid, 0.5, 0.3)

    if with_jac:
        ws.measurement_sensorInit()
        ws.jac_targetsAddSurface(target="t")
        ws.jac_targetsFinalize()

    ws.spectral_rad_surface_agendaSet(option=option)
    return ws


def run_option(freq_grid, option, suns, with_jac=False):
    """Build a workspace for the option, execute the agenda, return (rad, jac)."""
    ws = build_workspace(freq_grid, option, suns, with_jac)
    ws.spectral_rad_surface_agendaExecute()
    return stokes_array(ws, len(freq_grid)), jac_array(ws)


# ============================================================================
# Test 1: Each option builds with the exact expected method chain
# ============================================================================
def test_options_build():
    freq_grid = [10e9, 100e9, 183e9]

    for option, expected in EXPECTED_CHAINS.items():
        ws = build_workspace(freq_grid, option, suns=[make_sun(0.0, 0.0)])
        # finalize() inserts anonymous set-methods (_angle_cut, _refinement) for gin
        # defaults, and the agenda creator presets gins via @-prefixed named inputs
        # (@exclude_suns); neither is part of the method chain proper
        names = [m.name for m in ws.spectral_rad_surface_agenda.methods if not m.name.startswith(("_", "@"))]
        assert names == expected, f"{option}: wrong method chain:\n  got      = {names}\n  expected = {expected}"

    print("Test 1 passed: all three options build with the expected method chains")


# ============================================================================
# Test 2: Each option executes
# ============================================================================
def test_options_execute():
    freq_grid = [10e9, 100e9, 183e9]
    nf = len(freq_grid)

    for option in EXPECTED_CHAINS:
        rad, jac = run_option(freq_grid, option, suns=[make_sun(0.0, 0.0)])

        assert rad.shape == (nf, 4), f"{option}: spectral_rad shape mismatch: {rad.shape}"
        assert jac.shape[1] == nf, f"{option}: spectral_rad_jac shape mismatch: {jac.shape}"
        assert np.all(np.isfinite(rad)), f"{option}: non-finite radiance:\n{rad}"
        assert np.all(rad[:, 0] > 0.0), f"{option}: no surface emission in radiance:\n{rad}"

    print("Test 2 passed: all three options execute and produce finite emission")


# ============================================================================
# Test 3: DiffuseOnly + DirectOnly == SurfaceScatteringModel
# ============================================================================
def test_options_sum_identity():
    """The full option must equal the sum of the two half options.

    The sun is placed at 45 deg, outside every traced direction (the single
    up-looking quadrature point and the zenith mirror/glint direction).  With
    the sun at zenith the full option runs the non-Direct methods with
    exclude_suns = 1 and drops those directions, so the identity would no
    longer hold -- that difference *is* the sun de-duplication, pinned by
    tests/core/surf/spectral_rad_surface_scattering_sun_double_count.py.
    """
    freq_grid = [10e9, 100e9, 183e9]
    suns = [make_sun(45.0, 0.0)]

    rad_full, jac_full = run_option(freq_grid, "SurfaceScatteringModel", suns, with_jac=True)
    rad_diffuse, jac_diffuse = run_option(freq_grid, "SurfaceScatteringModelDiffuseOnly", suns, with_jac=True)
    rad_direct, jac_direct = run_option(freq_grid, "SurfaceScatteringModelDirectOnly", suns, with_jac=True)

    assert jac_full.shape == jac_diffuse.shape == jac_direct.shape, (
        f"Jacobian shape mismatch: {jac_full.shape}, {jac_diffuse.shape}, {jac_direct.shape}"
    )
    assert jac_full.shape[0] > 0, "Expected a non-empty jacobian with a surface target"

    expected = rad_diffuse + rad_direct
    assert np.allclose(rad_full, expected, rtol=1e-10, atol=0.0), (
        f"Option sum identity violated for spectral_rad:\n  full     = {rad_full}\n  expected = {expected}"
    )

    expected_jac = jac_diffuse + jac_direct
    assert np.allclose(jac_full, expected_jac, rtol=1e-10, atol=0.0), (
        f"Option sum identity violated for spectral_rad_jac:\n  full     = {jac_full}\n  expected = {expected_jac}"
    )

    print("Test 3 passed: DiffuseOnly + DirectOnly equals SurfaceScatteringModel for rad and jac")


# ============================================================================
# Test 4: DirectOnly without suns returns only surface emission
# ============================================================================
def test_direct_only_without_suns():
    freq_grid = [10e9, 100e9, 183e9]

    rad_nosun, _ = run_option(freq_grid, "SurfaceScatteringModelDirectOnly", suns=[])
    rad_sun, _ = run_option(freq_grid, "SurfaceScatteringModelDirectOnly", suns=[make_sun(0.0, 0.0)])

    assert np.all(rad_nosun[:, 0] > 0.0), f"DirectOnly without suns lost the surface emission:\n{rad_nosun}"
    assert np.all(rad_sun[:, 0] > rad_nosun[:, 0]), (
        f"Suns added no scattered beam term to DirectOnly:\n  with suns = {rad_sun}\n  no suns   = {rad_nosun}"
    )

    rad_diffuse_nosun, _ = run_option(freq_grid, "SurfaceScatteringModelDiffuseOnly", suns=[])
    rad_diffuse_sun, _ = run_option(freq_grid, "SurfaceScatteringModelDiffuseOnly", suns=[make_sun(0.0, 0.0)])
    assert np.allclose(rad_diffuse_nosun, rad_diffuse_sun, rtol=1e-12, atol=0.0), (
        f"DiffuseOnly must not depend on suns:\n  with suns = {rad_diffuse_sun}\n  no suns   = {rad_diffuse_nosun}"
    )

    print("Test 4 passed: DirectOnly without suns is pure surface emission, DiffuseOnly is sun-independent")


if __name__ == "__main__":
    test_options_build()
    test_options_execute()
    test_options_sum_identity()
    test_direct_only_without_suns()
