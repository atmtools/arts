"""
Transformed Jacobian targets must survive an XML round trip, in ascii and
binary, and behave exactly like the originals afterwards.

Every transformation struct is written and read: relative, logarithmic,
log-relative, and relative humidity for the atmosphere, the first three for the
surface and the subsurface, and the polynomial fits of sensors and errors.
"""

import os
import re
import tempfile

import pyarts3 as pyarts
import numpy as np

NF = 11
f = np.linspace(-5e6, 5e6, NF) + 118750348044.712

ws = pyarts.workspace.Workspace()

ws.abs_bands = pyarts.arts.AbsorptionBands()
ws.abs_speciesSet(species=["H2O", "O3", "CO2"])
ws.water_equivalent_pressure_operatorMK05()
ws.measurement_sensorsSimpleGaussian(
    freq_grid=f, std=1, pos=[100e3, 0, 0], los=[180.0, 0.0], pol="I"
)


def reset_fields():
    ws.surf_fieldPlanet(option="Earth")
    ws.surf_field[pyarts.arts.SurfaceKey("t")] = 295.0
    ws.subsurf_field = pyarts.arts.SubsurfaceField()
    ws.subsurf_field.bottom_depth = -30.0
    ws.subsurf_field["t"] = pyarts.arts.GeodeticField3(
        grids=[np.linspace(0, -30, 11)[::-1], [0], [0]],
        grid_names=["altitude", "latitude", "longitude"],
        data=np.linspace(250, 300, 11).reshape(11, 1, 1),
    )
    ws.atm_fieldRead(
        toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
    )


reset_fields()


def targets(surface_toggle, subsurface_toggle):
    """Targets with every atmospheric transformation and the given surface ones"""
    ws.jac_targetsInit()
    ws.jac_targetsAddTemperature()
    ws.jac_targetsAddSpeciesVMR(species="H2O")
    ws.jac_targetsAddSpeciesVMR(species="O3")
    ws.jac_targetsAddSpeciesVMR(species="CO2")
    ws.jac_targetsAddSurface(target=pyarts.arts.SurfaceKey("t"))
    ws.jac_targetsAddSubsurface(target=pyarts.arts.SubsurfaceKey("t"))
    ws.jac_targetsAddSensorFrequencyPolyOffset(sensor_elem=0, d=1e3, polyorder=0)
    ws.jac_targetsAddErrorPolyFit(t=f - f[0], sensor_elem=0, polyorder=1)
    ws.jac_targetsFinalize()

    ws.jac_targetsToggleRelativeAtmTarget(key="t")
    ws.jac_targetsToggleRelativeHumidityAtmTarget(key="H2O")
    ws.jac_targetsToggleLogarithmicAtmTarget(key="O3")
    ws.jac_targetsToggleLogRelAtmTarget(key="CO2")
    getattr(ws, f"jac_targetsToggle{surface_toggle}SurfaceTarget")(key=pyarts.arts.SurfaceKey("t"))
    getattr(ws, f"jac_targetsToggle{subsurface_toggle}SubsurfaceTarget")(key=pyarts.arts.SubsurfaceKey("t"))
    return pyarts.arts.JacobianTargets(ws.jac_targets)


def outputs(targets):
    """Results that use every forward, inverse, and Jacobian transformation"""
    reset_fields()
    ws.jac_targets = targets

    ws.model_state_vecFromData()
    x = np.array(ws.model_state_vec)

    ws.measurement_jac = np.ones((NF, x.size))
    ws.measurement_jacAtmosphereTransformation()
    ws.measurement_jacSurfaceTransformation()
    ws.measurement_jacSubsurfaceTransformation()
    ws.measurement_jacSensorTransformation()
    jac = np.array(ws.measurement_jac)

    ws.model_state_vec = x * 1.01
    ws.atm_fieldFromModelState(model_state_targets=targets)
    ws.surf_fieldFromModelState(model_state_targets=targets)
    ws.subsurf_fieldFromModelState(model_state_targets=targets)
    ws.measurement_sensorsFromModelState(model_state_targets=targets)
    ws.model_state_vecFromData()
    x_back = np.array(ws.model_state_vec)

    x_err = x * 1.01
    err_target = targets.error[0]
    x_err[err_target.x_start : err_target.x_start + err_target.x_size] = [0.5, 1e-6]
    ws.model_state_vec = x_err
    ws.measurement_vec_errorFromModelState(model_state_targets=targets)
    err = np.array(ws.measurement_vec_error)
    err_jac = np.array(ws.measurement_jac_error)

    ws.measurement_sensorsSimpleGaussian(
        freq_grid=f, std=1, pos=[100e3, 0, 0], los=[180.0, 0.0], pol="I"
    )

    return {"x": x, "jac": jac, "x_back": x_back, "err": err, "err_jac": err_jac}


STRUCTS = {
    "relfwd", "relinv", "logfwd", "loginv", "logrelfwd", "logrelinv", "rhfwd", "rhinv",
    "polyfit_t", "polyinv_t", "polyfit_sensor_offset_t", "polyinv_sensor_offset_t",
}

written = set()
with tempfile.TemporaryDirectory() as tmp:
    for toggle in ("Relative", "Logarithmic", "LogRel"):
        original = targets(toggle, toggle)
        ref = outputs(original)
        assert np.allclose(ref["x_back"], ref["x"] * 1.01), f"{toggle}: original targets are inconsistent"

        for filetype in ("ascii", "binary"):
            fn = original.savexml(os.path.join(tmp, f"targets_{toggle}_{filetype}.xml"), type=filetype)
            if filetype == "ascii":
                with open(fn, encoding="utf-8") as file:
                    written |= set(re.findall(r'<Variant[^>]*type="(\w+)"', file.read()))

            res = outputs(pyarts.arts.JacobianTargets.fromxml(fn))
            for key in ref:
                assert np.array_equal(res[key], ref[key]), f"{toggle}, {filetype}: {key} differs after reading"

assert STRUCTS <= written, f"Transformation structs never written: {sorted(STRUCTS - written)}"
