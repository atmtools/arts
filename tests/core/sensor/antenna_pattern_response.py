import numpy as np
import pyarts3 as pyarts


def assert_close(name, got, expected, atol=1e-10):
    got = np.asarray(got, dtype=float)
    expected = np.asarray(expected, dtype=float)
    print(f"{name}: got      = {got}")
    print(f"{name}: expected = {expected}")
    np.testing.assert_allclose(got, expected, atol=atol, rtol=0.0)


def test_pencil_pattern():
    pattern = pyarts.arts.SensorAntennaPattern.pencil()
    dzen = pyarts.arts.AscendingGrid(np.array([-1.0, 0.0, 1.0]))
    dazi = pyarts.arts.AscendingGrid(np.array([-1.0, 0.0, 1.0]))
    raw = pattern.raw_sensor(dzen, dazi)

    assert_close(
        "pencil raw sensor",
        np.array(raw.data),
        np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
            ]
        ),
    )


def test_gaussian_pattern_and_fwhm_equivalence():
    std_zenith = 1.0
    std_azimuth = 2.0
    fwhm_factor = 2.0 * np.sqrt(2.0 * np.log(2.0))

    pattern_std = pyarts.arts.SensorAntennaPattern.gaussian(
        std_zenith, std_azimuth
    )
    pattern_fwhm = pyarts.arts.SensorAntennaPattern.gaussian_fwhm(
        fwhm_factor * std_zenith, fwhm_factor * std_azimuth
    )

    assert_close(
        "gaussian point evaluation",
        np.array([pattern_std(1.0, 2.0)]),
        np.array([np.exp(-1.0)]),
    )

    dzen = pyarts.arts.AscendingGrid(np.array([-1.0, 0.0, 1.0]))
    dazi = pyarts.arts.AscendingGrid(np.array([-2.0, 0.0, 2.0]))

    std_raw = pattern_std.raw_sensor(dzen, dazi)
    fwhm_raw = pattern_fwhm.raw_sensor(dzen, dazi)
    normalized_raw = pattern_std.raw_sensor(dzen, dazi, normalize=True)

    expected_raw = np.array(
        [
            [np.exp(-1.0), np.exp(-0.5), np.exp(-1.0)],
            [np.exp(-0.5), 1.0, np.exp(-0.5)],
            [np.exp(-1.0), np.exp(-0.5), np.exp(-1.0)],
        ]
    )
    expected_normalized = expected_raw / expected_raw.sum()

    assert_close("gaussian raw sensor", np.array(std_raw.data), expected_raw)
    assert_close(
        "gaussian vs fwhm raw sensor",
        np.array(fwhm_raw.data),
        expected_raw,
    )
    assert_close(
        "gaussian normalized raw sensor",
        np.array(normalized_raw.data),
        expected_normalized,
    )


def test_lookup_pattern():
    pattern = pyarts.arts.SensorAntennaPattern.lookup(
        pyarts.arts.AscendingGrid(np.array([0.0, 1.0])),
        pyarts.arts.AscendingGrid(np.array([0.0, 1.0])),
        pyarts.arts.Matrix(np.array([[1.0, 2.0], [3.0, 4.0]])),
    )

    assert_close(
        "lookup bilinear center",
        np.array([pattern(0.5, 0.5)]),
        np.array([2.5]),
    )
    assert_close(
        "lookup outside support",
        np.array([pattern(2.0, 2.0)]),
        np.array([0.0]),
    )


def test_lookup_pattern_with_typed_angular_grids():
    pattern = pyarts.arts.SensorAntennaPattern.lookup(
        pyarts.arts.ZenGrid(np.array([0.0])),
        pyarts.arts.AziGrid(np.array([0.0, 120.0, 240.0])),
        pyarts.arts.Matrix(np.array([[1.0, 2.0, 3.0]])),
    )

    assert_close(
        "typed angular-grid lookup exact point",
        np.array([pattern(0.0, 120.0)]),
        np.array([2.0]),
    )
    assert_close(
        "typed angular-grid lookup interpolated midpoint",
        np.array([pattern(0.0, 60.0)]),
        np.array([1.5]),
    )
    assert_close(
        "typed angular-grid lookup outside support",
        np.array([pattern(0.0, 300.0)]),
        np.array([0.0]),
    )


def test_pattern_to_sensor_via_raw_sensor():
    pattern = pyarts.arts.SensorAntennaPattern.gaussian(1.0, 2.0)
    dzen = pyarts.arts.AscendingGrid(np.array([-1.0, 0.0, 1.0]))
    dazi = pyarts.arts.AscendingGrid(np.array([-2.0, 0.0, 2.0]))
    raw = pattern.raw_sensor(dzen, dazi, normalize=True)

    ws = pyarts.Workspace()
    ws.measurement_sensor = pyarts.arts.ArrayOfSensorObsel()
    ws.measurement_sensor_meta = pyarts.arts.ArrayOfSensorMetaInfo()
    ws.measurement_sensorAddRawSensor(
        freq_grid=pyarts.arts.AscendingGrid(np.array([100.0, 101.0])),
        pos=pyarts.arts.Vector3(np.array([600000.0, 0.0, 0.0])),
        los=pyarts.arts.Vector2(np.array([180.0, 0.0])),
        raw_sensor_perturbation=raw,
        normalize=0,
    )

    expected_weight_vector = np.array(raw.data).reshape(-1)

    assert_close(
        "sensor count from raw antenna pattern",
        np.array([len(ws.measurement_sensor)]),
        np.array([2.0]),
    )
    assert_close(
        "sensor meta count from raw antenna pattern",
        np.array([len(ws.measurement_sensor_meta)]),
        np.array([1.0]),
    )
    assert_close(
        "sensor poslos size from raw antenna pattern",
        np.array([len(ws.measurement_sensor[0].poslos)]),
        np.array([9.0]),
    )
    assert_close(
        "first sensor poslos entry",
        np.array(ws.measurement_sensor[0].poslos[0]),
        np.array([600000.0, 0.0, 0.0, 179.0, -2.0]),
    )
    assert_close(
        "sensor weights from raw antenna pattern",
        np.array(ws.measurement_sensor[0].weight_matrix.reduce(along_freq=True)),
        expected_weight_vector,
    )


test_pencil_pattern()
test_gaussian_pattern_and_fwhm_equivalence()
test_lookup_pattern()
test_lookup_pattern_with_typed_angular_grids()
test_pattern_to_sensor_via_raw_sensor()

print("\nAll antenna pattern checks passed.")