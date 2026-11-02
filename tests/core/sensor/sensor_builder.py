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
            pyarts.arts.sensor.BoxChannel(pyarts.arts.AscendingGrid([100.0, 101.0])),
            pyarts.arts.sensor.DiracChannel(200.0),
        ],
        pyarts.arts.sensor.PencilBeamAntenna(),
    )

    pos = [
        pyarts.arts.Vector3([600e3, 10.0, 20.0]),
        pyarts.arts.Vector3([601e3, 11.0, 21.0]),
    ]
    los = [
        pyarts.arts.Vector2([20.0, 30.0]),
        pyarts.arts.Vector2([40.0, 50.0]),
    ]

    obsels, meta = builder(pos, los, pyarts.arts.planets.Earth.ellipsoid)

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
        [pyarts.arts.sensor.DiracChannel()],
        pyarts.arts.sensor.PencilBeamAntenna(),
    )

    poslos = pyarts.arts.SensorPosLosVector()
    poslos.value = np.array(
        [
            [600e3, 10.0, 20.0, 20.0, 30.0],
            [601e3, 11.0, 21.0, 40.0, 50.0],
        ]
    )

    obsels, meta = builder(poslos, pyarts.arts.planets.Earth.ellipsoid)

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
        ],  pyarts.arts.planets.Earth.ellipsoid)
    except ValueError:
        pass
    else:
        raise AssertionError(
            "Sensor builder must reject mismatching position/LOS lengths")


def test_sensor_builder_xml_round_trip():
    sensor = pyarts.arts.sensor
    gridded = sensor.AntennaPatternField()
    gridded.grids = [pyarts.arts.ZenGrid([0.0, 0.5, 1.0]),
                     pyarts.arts.AziGrid([0.0, 90.0, 180.0, 270.0])]
    gridded.data = np.array([[1.0] * 4, [0.6, 0.5, 0.6, 0.5], [0.1] * 4])

    antennas = [
        sensor.PencilBeamAntenna(),
        sensor.GaussianAntenna([0.0, 0.5, 1.0], 0.4, azi_grid_size=4),
        sensor.GaussianAiryAntenna([0.0, 0.1, 0.2], 1.0, azi_grid_size=3),
        sensor.GriddedAntennaPattern(gridded),
    ]

    vertical = sensor.BoxChannel(pyarts.arts.AscendingGrid([89e9, 89.5e9, 90e9]))
    vertical.polarization.state = [1.0, 1.0, 0.0, 0.0]
    channels = [vertical, sensor.DiracChannel(183e9)]

    pos = [pyarts.arts.Vector3([800e3, 10.0, 20.0])]
    los = [pyarts.arts.Vector2([150.0, 30.0])]
    ell = pyarts.arts.planets.Earth.ellipsoid

    for antenna in antennas:
        name = type(antenna).__name__
        builder = pyarts.arts.SensorBuilder(channels, antenna)
        builder.preserve_common_frequency_grid = True
        builder.savexml("sensor_builder.xml")
        read = pyarts.arts.SensorBuilder.fromxml("sensor_builder.xml")

        assert type(read.antenna) is type(builder.antenna), name
        assert read.preserve_common_frequency_grid, name
        assert len(read.channels) == len(channels), name

        expected, _ = builder(pos, los, ell)
        got, _ = read(pos, los, ell)
        assert len(got) == len(expected), name
        for g, e in zip(got, expected):
            assert np.array_equal(g.f_grid, e.f_grid), name
            assert np.array_equal(g.poslos, e.poslos), name
            assert np.array_equal(g.weight_matrix.dense, e.weight_matrix.dense), name


test_sensor_builder_from_sequences()
test_sensor_builder_from_poslos_vector_and_length_guard()
test_sensor_builder_xml_round_trip()
