"""Separate state mappings and derivative targets through the real agendas."""

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts


def check_kernel(option):
    ws = pyarts.Workspace()
    ws.freq_grid = [1e9, 2e9, 3e9]
    ws.surf_fieldEarth()
    ws.surf_field["t"] = 280.0
    ws.atm_fieldInit(toa=0.0)
    ws.abs_bands = arts.AbsorptionBands()
    ws.measurement_sensorSimple(pos=[0, 0, 0], los=[180, 0])
    ws.spectral_rad_transform_operatorSet(option="1")
    ws.jac_targetsInit()
    ws.jac_targetsAddSurface(target="t")
    ws.jac_targetsAddErrorPolyFit(t=[-1, 0, 1], sensor_elem=0, polyorder=1)
    ws.jac_targetsFinalize()
    ws.model_state_targets = ws.jac_targets
    original = str(ws.jac_targets)
    ws.measurement_inversion_agendaSet(option=option)
    for agenda in (ws.inversion_iterate_agenda, ws.measurement_inversion_agenda):
        for name in ("jac_targets", "model_state_targets"):
            assert f"Copies the global *{name}*" not in agenda.document()

    seen = []

    def forward(local):
        nf = len(local.get("freq_grid"))
        nx = local.get("jac_targets").x_size()
        temperature = float(local.get("surf_field")["t"].data)
        seen.append((nx, temperature))
        if temperature == 7:
            raise RuntimeError("deliberate value-only failure")
        radiance = np.zeros((nf, 4))
        radiance[:, 0] = temperature**2
        derivative = np.zeros((nx, nf, 4))
        if nx:
            derivative[0, :, 0] = 2 * temperature
        local.get("spectral_rad").value = arts.StokvecVector(radiance)
        local.get("spectral_rad_jac").value = arts.StokvecMatrix(derivative)
        local.get("ray_path").clear()
        local.get("ray_path").append(arts.PropagationPathPoint())

    callback = arts.CallbackOperator(
        forward,
        ["freq_grid", "jac_targets", "surf_field"],
        ["spectral_rad", "spectral_rad_jac", "ray_path"],
    )
    agenda = arts.Agenda("spectral_rad_observer_agenda")
    agenda.add(arts.Method("analytical_radiance", callback))
    agenda.finalize(True)
    ws.spectral_rad_observer_agenda = agenda

    def evaluate(with_jacobian, temperature, offset, slope):
        ws.model_state_vec = [temperature, offset, slope]
        seen.clear()
        ws.inversion_iterate_agendaExecute(
            jac_targets=ws.jac_targets if with_jacobian else arts.JacobianTargets()
        )
        assert seen and all(
            item == (3 if with_jacobian else 0, temperature) for item in seen
        )
        expected = temperature**2 + offset + slope * np.array([-1, 0, 1])
        np.testing.assert_allclose(ws.measurement_vec_fit, expected, rtol=0, atol=1e-12)
        assert float(ws.surf_field["t"].data) == temperature
        assert str(ws.jac_targets) == str(ws.model_state_targets) == original
        if with_jacobian:
            np.testing.assert_allclose(
                ws.measurement_jac,
                [
                    [2 * temperature, 1, -1],
                    [2 * temperature, 1, 0],
                    [2 * temperature, 1, 1],
                ],
                rtol=0,
                atol=1e-12,
            )
        else:
            assert np.asarray(ws.measurement_jac).size == 0

    evaluate(True, 3, 2, 0.5)
    evaluate(False, 4, -1, 0.25)
    evaluate(False, 5, 3, -2)
    evaluate(True, 3, 2, 0.5)
    ws.model_state_vec = [7, 1, 1]
    try:
        ws.inversion_iterate_agendaExecute(jac_targets=arts.JacobianTargets())
    except RuntimeError as error:
        assert "deliberate value-only failure" in str(error)
    else:
        raise AssertionError("Failing trial was not evaluated")
    assert str(ws.jac_targets) == str(ws.model_state_targets) == original
    evaluate(True, 3, 2, 0.5)

    # Independently eliminate the offset from the normal equations for identity
    # covariances. The positive root of 3*T**3 - 67*T - 6 is the MAP minimum.
    for name in ("model_state_covmat", "measurement_vec_error_covmat"):
        covariance = arts.CovarianceMatrix()
        covariance.blocks = [
            arts.Block(
                arts.Range(0, 3), arts.Range(0, 3), (0, 0), arts.Matrix(np.eye(3))
            )
        ]
        setattr(ws, name, covariance)
    ws.model_state_vec_apriori = [3, 2, 0.5]
    ws.model_state_vec = []
    ws.measurement_vec_fit = []
    ws.measurement_jac = arts.Matrix()
    ws.measurement_vec = [24, 25, 26]
    # OEM must supply its full mapping explicitly, even if the outer workspace
    # mapping is stale or empty.
    ws.model_state_targets = arts.JacobianTargets()
    seen.clear()
    ws.OEM(
        method="lm", stop_dx=1e-12, max_iter=100, lm_ga_settings=arts.OEMLMSettings()
    )
    assert ws.oem_diagnostics[0] == 0, ws.oem_diagnostics
    temperature = max(np.roots([3, 0, -67, -6]))
    expected = [temperature, 2 - 0.75 * (temperature**2 - 23), 5 / 6]
    np.testing.assert_allclose(ws.model_state_vec, expected, rtol=0, atol=1e-6)
    assert any(nx == 0 and temp != 3 for nx, temp in seen)
    assert any(nx == 3 for nx, _ in seen)
    assert str(ws.jac_targets) == original
    np.testing.assert_allclose(
        float(ws.surf_field["t"].data), temperature, rtol=0, atol=1e-6
    )


for kernel in ("LowMemory", "HighPerformance"):
    check_kernel(kernel)
