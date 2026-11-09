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
            pyarts.arts.sensors.BoxChannel(pyarts.arts.AscendingGrid([100.0, 101.0])),
            pyarts.arts.sensors.DiracChannel(200.0),
        ],
        pyarts.arts.sensors.PencilBeamAntenna(),
    )

    pos = [
        pyarts.arts.Vector3([600e3, 10.0, 20.0]),
        pyarts.arts.Vector3([601e3, 11.0, 21.0]),
    ]
    los = [
        pyarts.arts.Vector2([20.0, 30.0]),
        pyarts.arts.Vector2([40.0, 50.0]),
    ]

    sensor = builder(pos, los, pyarts.arts.planets.Earth.ellipsoid)

    # Without a name, the sensor is named after its channels; its name is that of its meta
    assert sensor.name == "sensor-with-2-channels"
    assert sensor.meta.dataname == sensor.name
    assert builder(pos, los, pyarts.arts.planets.Earth.ellipsoid, name="radiometer").name == "radiometer"

    # One sensor: the elements by placement and then channel, as its placement by channel meta
    assert len(sensor.elements) == 4
    assert np.shape(sensor.meta.data) == (2, 2)
    assert_close("meta placement grid", sensor.meta.grids[0], [0.0, 1.0])
    assert_close("meta channel grid", sensor.meta.grids[1], [0.0, 1.0])
    assert_close("frequencies of both channels", sensor.freq_grid, [100.0, 101.0, 200.0])
    assert_close("first pos", sensor.geometry[0][:3], pos[0])
    assert_close("first los", sensor.geometry[0][3:], los[0])
    assert_close("second placement pos", sensor.geometry[1][:3], pos[1])
    assert_close("second placement los", sensor.geometry[1][3:], los[1])
    assert_close("box channel weights at the first placement",
                 sensor.elements[0].dense[0, :2],
                 [[0.5, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 0.0]])
    assert_close("Dirac channel weights at the second placement",
                 sensor.elements[3].dense[1, 2],
                 [1.0, 0.0, 0.0, 0.0])

    # Extracting a radiance of each geometry and frequency sums the weights of each element
    spectral_rad = pyarts.arts.StokvecMatrix(np.arange(2 * 3 * 4, dtype=float).reshape(2, 3, 4))
    i = np.asarray(spectral_rad)[..., 0]
    assert_close("extracted elements", sensor.extract(spectral_rad),
                 [0.5 * (i[0, 0] + i[0, 1]), i[0, 2], 0.5 * (i[1, 0] + i[1, 1]), i[1, 2]])


def test_sensor_builder_from_poslos_vector_and_length_guard():
    builder = pyarts.arts.SensorBuilder(
        [pyarts.arts.sensors.DiracChannel()],
        pyarts.arts.sensors.PencilBeamAntenna(),
    )

    poslos = pyarts.arts.SensorPosLosVector()
    poslos.value = np.array(
        [
            [600e3, 10.0, 20.0, 20.0, 30.0],
            [601e3, 11.0, 21.0, 40.0, 50.0],
        ]
    )

    sensor = builder(poslos, pyarts.arts.planets.Earth.ellipsoid)

    assert len(sensor.elements) == 2
    assert np.shape(sensor.meta.data) == (2, 1)
    assert_close("poslos builder first los", sensor.geometry[0][3:], [20.0, 30.0])
    assert_close("poslos builder second los", sensor.geometry[1][3:], [40.0, 50.0])

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
    sensors = pyarts.arts.sensors
    gridded = sensors.AntennaPatternField()
    gridded.grids = [pyarts.arts.ZenGrid([0.0, 0.5, 1.0]),
                     pyarts.arts.AziGrid([0.0, 90.0, 180.0, 270.0])]
    gridded.data = np.array([[1.0] * 4, [0.6, 0.5, 0.6, 0.5], [0.1] * 4])

    antennas = [
        sensors.PencilBeamAntenna(),
        sensors.GaussianAntenna([0.0, 0.5, 1.0], 0.4, azi_grid_size=4),
        sensors.GaussianAiryAntenna([0.0, 0.1, 0.2], 1.0, azi_grid_size=3),
        sensors.GriddedAntennaPattern(gridded),
        sensors.CameraAntenna.thin_lens(2, 3, 0.01, 0.01, 0.05, 1e3,
                                        beam=sensors.GaussianAntenna([0.0, 0.2], 0.1, azi_grid_size=3)),
    ]

    vertical = sensors.BoxChannel(pyarts.arts.AscendingGrid([89e9, 89.5e9, 90e9]))
    vertical.polarization.state = [1.0, 1.0, 0.0, 0.0]
    channels = [vertical, sensors.DiracChannel(183e9)]

    pos = [pyarts.arts.Vector3([800e3, 10.0, 20.0])]
    los = [pyarts.arts.Vector2([150.0, 30.0])]
    ell = pyarts.arts.planets.Earth.ellipsoid

    for antenna in antennas:
        name = type(antenna).__name__
        builder = pyarts.arts.SensorBuilder(channels, antenna)
        builder.savexml("sensor_builder.xml")
        read = pyarts.arts.SensorBuilder.fromxml("sensor_builder.xml")

        assert type(read.antenna) is type(builder.antenna), name
        assert len(read.channels) == len(channels), name

        expected = builder(pos, los, ell)
        got = read(pos, los, ell)
        assert np.array_equal(got.freq_grid, expected.freq_grid), name
        assert np.array_equal(got.geometry, expected.geometry), name
        assert np.array_equal(got.meta.data, expected.meta.data), name
        assert len(got.elements) == len(expected.elements), name
        for g, e in zip(got.elements, expected.elements):
            assert np.array_equal(g.dense, e.dense), name

        # The sensor itself survives a round trip through a file
        expected.savexml("sensor.xml")
        again = pyarts.arts.Sensor.fromxml("sensor.xml")
        assert np.array_equal(again.freq_grid, expected.freq_grid), name
        assert np.array_equal(again.geometry, expected.geometry), name
        assert type(again.meta) is type(expected.meta), name
        assert np.array_equal(again.meta.grids[0], expected.meta.grids[0]), name
        for g, e in zip(again.elements, expected.elements):
            assert np.array_equal(g.dense, e.dense), name


def test_sensor_builder_cross_track():
    sensors = pyarts.arts.sensors

    # Vertical at nadir and turning with the scan, as behind a rotating reflector
    turning = sensors.DiracChannel(89e9)
    turning.polarization = sensors.Polarization([1.0, 1.0, 0.0, 0.0], scan_rotation=1.0)

    # Horizontal at any pointing
    fixed = sensors.DiracChannel(183e9)
    fixed.polarization = "Ih"

    builder = pyarts.arts.SensorBuilder([turning, fixed], sensors.PencilBeamAntenna())

    pos = [
        pyarts.arts.Vector3([817e3, 10.0, 20.0]),
        pyarts.arts.Vector3([817e3, 11.0, 21.0]),
    ]
    heading = [30.0, -170.0]
    scan = [-20.0, 0.0, 50.0]
    ell = pyarts.arts.planets.Earth.ellipsoid

    sensor = builder.cross_track(pos, heading, scan, ell)

    # [180 - |s|, heading - 90] to the left and [180 - |s|, heading + 90] to the right, azimuth in [-180, 180]
    expected_los = [
        [[160.0, -60.0], [180.0, 120.0], [130.0, 120.0]],
        [[160.0, 100.0], [180.0, -80.0], [130.0, -80.0]],
    ]

    # Every position and scan angle is a placement, each with the two channels
    assert len(sensor.elements) == len(pos) * len(scan) * 2
    assert len(sensor.geometry) == len(pos) * len(scan)
    assert np.shape(sensor.meta.data) == (len(pos) * len(scan), 2)
    assert_close("frequencies", sensor.freq_grid, [89e9, 183e9])

    for ipos in range(len(pos)):
        for iscan, s in enumerate(scan):
            iplace = ipos * len(scan) + iscan
            name = f"position {ipos}, scan angle {s}"

            assert_close(f"{name} pos", sensor.geometry[iplace][:3], pos[ipos])
            assert_close(f"{name} los", sensor.geometry[iplace][3:], expected_los[ipos][iscan])

            # The linear polarization turns by scan_rotation times the signed scan angle
            cos2s = np.cos(np.deg2rad(2.0 * s))
            sin2s = np.sin(np.deg2rad(2.0 * s))
            expected_state = [[1.0, cos2s, sin2s, 0.0], [1.0, -1.0, 0.0, 0.0]]

            for ichan in range(2):
                weights = sensor.elements[iplace * 2 + ichan].dense
                chan = f"{name}, channel {ichan}"
                assert_close(f"{chan} weights", weights[iplace, ichan], expected_state[ichan])
                assert_close(f"{chan} only its placement and frequency",
                             np.count_nonzero(np.asarray(weights).any(axis=-1)), 1)

    try:
        builder(pos, [pyarts.arts.Vector2([160.0, 120.0])] * len(pos), ell)
    except RuntimeError:
        pass
    else:
        raise AssertionError(
            "A channel that turns with the scan must not be placed at a bare line of sight")

    try:
        builder.cross_track(pos, heading[:1], scan, ell)
    except RuntimeError:
        pass
    else:
        raise AssertionError(
            "Cross-track builder must reject mismatching position/heading lengths")


test_sensor_builder_from_sequences()
test_sensor_builder_from_poslos_vector_and_length_guard()
test_sensor_builder_xml_round_trip()
test_sensor_builder_cross_track()
