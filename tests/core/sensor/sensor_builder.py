import numpy as np
import pyarts3 as pyarts


def assert_close(name, got, expected, atol=1e-12):
    got = np.asarray(got, dtype=float)
    expected = np.asarray(expected, dtype=float)
    print(f"{name}: got      = {got}")
    print(f"{name}: expected = {expected}")
    np.testing.assert_allclose(got, expected, atol=atol, rtol=0.0)


def test_sensor_builder_from_sequences():
    builder = pyarts.arts.SensorBuilder(
        [
            pyarts.arts.SensorBoxChannel(pyarts.arts.AscendingGrid([100.0, 101.0])),
            pyarts.arts.SensorDiracChannel(200.0),
        ],
        pyarts.arts.SensorPencilBeamAntenna(),
    )

    pos = [
        pyarts.arts.Vector3([600e3, 10.0, 20.0]),
        pyarts.arts.Vector3([601e3, 11.0, 21.0]),
    ]
    los = [
        pyarts.arts.Vector2([20.0, 30.0]),
        pyarts.arts.Vector2([40.0, 50.0]),
    ]

    obsels, meta = builder(pos, los)

    assert len(obsels) == 4
    assert len(meta) == 2
    assert meta[0].count == 2
    assert meta[1].count == 2
    assert_close("meta[0] channel grid", meta[0].data.grids[0], [0.0, 1.0])
    assert_close("first pos", obsels[0].poslos[0][:3], pos[0])
    assert_close("first los", obsels[0].poslos[0][3:], los[0])
    assert_close("second geometry pos", obsels[2].poslos[0][:3], pos[1])
    assert_close("second geometry los", obsels[2].poslos[0][3:], los[1])
    assert_close("first channel grid", obsels[0].f_grid, [100.0, 101.0])
    assert_close("second channel grid", obsels[1].f_grid, [200.0])
    assert_close("first channel weights",
                 obsels[0].weight_matrix.dense[0, 0],
                 [0.5, 0.0, 0.0, 0.0])


def test_sensor_builder_from_poslos_vector_and_length_guard():
    builder = pyarts.arts.SensorBuilder(
        [pyarts.arts.SensorDiracChannel()],
        pyarts.arts.SensorPencilBeamAntenna(),
    )

    poslos = pyarts.arts.SensorPosLosVector()
    poslos.value = np.array(
        [
            [600e3, 10.0, 20.0, 20.0, 30.0],
            [601e3, 11.0, 21.0, 40.0, 50.0],
        ]
    )

    obsels, meta = builder(poslos)

    assert len(obsels) == 2
    assert len(meta) == 2
    assert meta[0].count == 1
    assert meta[1].count == 1
    assert_close("poslos builder first los", obsels[0].poslos[0][3:], [20.0, 30.0])
    assert_close("poslos builder second los", obsels[1].poslos[0][3:], [40.0, 50.0])

    try:
        builder([pyarts.arts.Vector3([600e3, 10.0, 20.0])], [
            pyarts.arts.Vector2([20.0, 30.0]),
            pyarts.arts.Vector2([40.0, 50.0]),
        ])
    except ValueError:
        pass
    else:
        raise AssertionError("SensorBuilder must reject mismatching position/LOS lengths")


test_sensor_builder_from_sequences()
test_sensor_builder_from_poslos_vector_and_length_guard()