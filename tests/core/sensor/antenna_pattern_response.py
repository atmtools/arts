import numpy as np
import pyarts3 as pyarts


def assert_close(name, got, expected, atol=1e-12):
    got = np.asarray(got, dtype=float)
    expected = np.asarray(expected, dtype=float)
    print(f"{name}: got      = {got}")
    print(f"{name}: expected = {expected}")
    np.testing.assert_allclose(got, expected, atol=atol, rtol=0.0)


def test_pencil_beam_binding():
    antenna = pyarts.arts.SensorPencilBeamAntenna()
    response = antenna.response

    assert response.ok()
    assert_close("pencil zenith grid", response.grids[0], np.array([0.0]))
    assert_close("pencil azimuth grid", response.grids[1], np.array([0.0]))
    assert_close("pencil response weight",
                 response.data[0, 0], np.array([1.0, 0.0, 0.0, 0.0]))

    samples = antenna(np.array([45.0, 30.0]))
    assert len(samples) == 1
    weight, los = samples[0]

    assert_close("pencil mapped weight", weight, np.array([1.0, 0.0, 0.0, 0.0]))
    assert_close("pencil mapped los", los, np.array([45.0, 30.0]))


def test_gaussian_binding_and_mapping():
    peak_weight = pyarts.arts.Stokvec([2.0, 1.0, 0.0, 0.0])
    antenna = pyarts.arts.SensorGaussianAntenna(
        pyarts.arts.ZenGrid(np.array([0.0, 1.0])),
        pyarts.arts.AziGrid(np.array([0.0, 180.0])),
        1.0,
        2.0,
        peak_weight,
    )
    response = antenna.response

    assert response.ok()
    assert_close("gaussian peak weight", response.data[0, 0], np.array(peak_weight))

    expected_offset_weight = np.exp(-0.5) * np.array(peak_weight)
    assert_close("gaussian [1, 0] weight", response.data[1, 0], expected_offset_weight)
    assert_close("gaussian [1, 180] weight",
                 response.data[1, 1], expected_offset_weight)

    samples = antenna(np.array([45.0, 30.0]))
    assert len(samples) == 4

    assert_close("gaussian bore los", samples[0][1], np.array([45.0, 30.0]))
    assert_close("gaussian local [1, 0] los", samples[2][1], np.array([46.0, 30.0]))
    assert_close("gaussian local [1, 180] los", samples[3][1], np.array([44.0, 30.0]))
    assert_close("gaussian local [1, 0] mapped weight",
                 samples[2][0], expected_offset_weight)
    assert_close("gaussian local [1, 180] mapped weight",
                 samples[3][0], expected_offset_weight)


test_pencil_beam_binding()
test_gaussian_binding_and_mapping()
