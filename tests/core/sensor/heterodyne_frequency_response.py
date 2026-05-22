import numpy as np
import pyarts3 as pyarts


def db_to_lin(db):
    return 10.0 ** (db / 10.0)


def ranges_as_lists(ranges):
    return [list(map(float, pair)) for pair in ranges]


def affine_from_ranges(global_range, local_range):
    g0, g1 = map(float, global_range)
    l0, l1 = map(float, local_range)
    slope = (g1 - g0) / (l1 - l0)
    intercept = g0 - slope * l0
    return [intercept, slope]


def assert_close(name, got, expected, atol=1e-12):
    got = np.asarray(got, dtype=float)
    expected = np.asarray(expected, dtype=float)
    print(f"{name}: got      = {got}")
    print(f"{name}: expected = {expected}")
    np.testing.assert_allclose(got, expected, atol=atol, rtol=0.0)


def inspect_case(name, selector, expected_global, expected_local, expected_affine):
    print(f"\n=== {name} ===")

    got_global = ranges_as_lists(selector.global_ranges)
    got_local = ranges_as_lists(selector.local_ranges)
    got_affine = [
        affine_from_ranges(global_range, local_range)
        for global_range, local_range in zip(got_global, got_local)
    ]

    print("Affine form per path: global = intercept + slope * local")
    assert_close(f"{name} global_ranges", got_global, expected_global)
    assert_close(f"{name} local_ranges", got_local, expected_local)
    assert_close(f"{name} affine", got_affine, expected_affine)


def make_triangular_filter():
    filt = pyarts.arts.SortedGriddedField1()
    filt.grids = (pyarts.arts.AscendingGrid(np.array([5.0, 10.0, 15.0])),)
    filt.data = pyarts.arts.Vector(np.array([0.0, 1.0, 0.0]))
    filt.gridnames = ("frequency",)
    filt.dataname = "triangular"
    return filt


def make_asymmetric_sideband_filter():
    filt = pyarts.arts.SortedGriddedField1()
    filt.grids = (
        pyarts.arts.AscendingGrid(np.array([7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0])),
    )
    filt.data = pyarts.arts.Vector(
        np.array(
            [
                db_to_lin(-10.0),
                db_to_lin(-10.0),
                db_to_lin(-10.0),
                db_to_lin(0.0),
                db_to_lin(1.0),
                db_to_lin(5.0),
                db_to_lin(10.0),
            ]
        )
    )
    filt.gridnames = ("frequency",)
    filt.dataname = "asymmetric_sidebands"
    return filt


def test_lowpass_then_lo_then_box_spectrometer():
    selector = pyarts.arts.SensorHeterodyneFrequencyRange()
    print("test_lowpass_then_lo_then_box_spectrometer")
    selector.lowpass(15.0)
    print("Apply lowpass keeping 0-15:", selector)
    selector.mix(16.0)
    print("Apply LO at 16            :", selector)

    inspect_case(
        "lowpass -> LO",
        selector,
        expected_global=[[15.0, 0.0]],
        expected_local=[[1.0, 16.0]],
        expected_affine=[[16.0, -1.0]],
    )

    channel = pyarts.arts.SensorBoxChannel(2.6, 4.6, 5)
    response = selector.channel_response(channel)

    assert_close(
        "lowpass -> LO channel points",
        np.array(response.grids[0]),
        np.array([11.4, 11.9, 12.4, 12.9, 13.4]),
    )
    assert_close(
        "lowpass -> LO channel weights",
        np.array(response.data),
        np.full(5, 0.2),
    )


def test_default_positive_range_survives_mixes():
    selector = pyarts.arts.SensorHeterodyneFrequencyRange()
    selector.mix(10.0)

    assert_close(
        "default positive range after one LO global_ranges",
        ranges_as_lists(selector.global_ranges),
        [[10.0, np.inf], [10.0, 0.0]],
    )
    assert_close(
        "default positive range after one LO local_ranges",
        ranges_as_lists(selector.local_ranges),
        [[0.0, np.inf], [0.0, 10.0]],
    )

    selector.mix(1.0)

    assert_close(
        "default positive range after two LOs global_ranges",
        ranges_as_lists(selector.global_ranges),
        [[11.0, np.inf], [11.0, 10.0], [9.0, 0.0], [9.0, 10.0]],
    )
    assert_close(
        "default positive range after two LOs local_ranges",
        ranges_as_lists(selector.local_ranges),
        [[0.0, np.inf], [0.0, 1.0], [0.0, 9.0], [0.0, 1.0]],
    )

    response = selector.channel_response(pyarts.arts.SensorDiracChannel(0.25))

    assert_close(
        "default positive range dirac points after two LOs",
        np.array(response.grids[0]),
        np.array([8.75, 9.25, 10.75, 11.25]),
    )
    assert_close(
        "default positive range dirac weights after two LOs",
        np.array(response.data),
        np.full(4, 0.25),
    )


def test_bandpass_lo_bandpass_lo_chain():
    selector = pyarts.arts.SensorHeterodyneFrequencyRange()
    selector.bandpass(np.array([5.0, 15.0]))
    selector.mix(3.0)
    selector.bandpass(np.array([3.0, 7.0]))
    selector.mix(6.0)

    inspect_case(
        "bandpass -> LO -> bandpass -> LO",
        selector,
        expected_global=[[9.0, 10.0], [9.0, 6.0]],
        expected_local=[[0.0, 1.0], [0.0, 3.0]],
        expected_affine=[[9.0, 1.0], [9.0, -1.0]],
    )

    channel = pyarts.arts.SensorBoxChannel(0.0, 1.0, 3)
    response = selector.channel_response(channel)

    assert_close(
        "bandpass -> LO -> bandpass -> LO channel points",
        np.array(response.grids[0]),
        np.array([8.0, 8.5, 9.0, 9.5, 10.0]),
    )
    assert_close(
        "bandpass -> LO -> bandpass -> LO channel weights",
        np.array(response.data),
        np.full(5, 1.0 / 6.0),
    )


def test_ideal_split_about_lo():
    selector = pyarts.arts.SensorHeterodyneFrequencyRange()
    selector.bandpass(np.array([5.0, 15.0]))
    selector.mix(10.0)

    inspect_case(
        "ideal split about LO",
        selector,
        expected_global=[[10.0, 15.0], [10.0, 5.0]],
        expected_local=[[0.0, 5.0], [0.0, 5.0]],
        expected_affine=[[10.0, 1.0], [10.0, -1.0]],
    )

    assert_close(
        "ideal split path 0 global response at 12",
        selector.global_response(np.array([12.0]), 0),
        np.array([1.0]),
    )
    assert_close(
        "ideal split path 1 global response at 8",
        selector.global_response(np.array([8.0]), 1),
        np.array([1.0]),
    )


def test_weighted_split_about_lo():
    selector = pyarts.arts.SensorHeterodyneFrequencyRange()
    selector.filter(make_triangular_filter())
    selector.mix(10.0)

    inspect_case(
        "weighted split about LO",
        selector,
        expected_global=[[10.0, 15.0], [10.0, 5.0]],
        expected_local=[[0.0, 5.0], [0.0, 5.0]],
        expected_affine=[[10.0, 1.0], [10.0, -1.0]],
    )

    assert_close(
        "weighted split path 0 local response at 2",
        selector.local_response(np.array([2.0]), 0),
        np.array([0.6]),
    )
    assert_close(
        "weighted split path 1 local response at 2",
        selector.local_response(np.array([2.0]), 1),
        np.array([0.6]),
    )

    channel = pyarts.arts.SensorDiracChannel(2.0)
    response = selector.channel_response(channel)

    assert_close(
        "weighted split channel points",
        np.array(response.grids[0]),
        np.array([8.0, 12.0]),
    )
    assert_close(
        "weighted split channel weights",
        np.array(response.data),
        np.array([0.3, 0.3]),
    )


def test_zero_frequency_is_not_double_counted():
    selector = pyarts.arts.SensorHeterodyneFrequencyRange()
    selector.bandpass(np.array([5.0, 15.0]))
    selector.mix(10.0)

    response = selector.channel_response(pyarts.arts.SensorDiracChannel(0.0))

    assert_close(
        "zero-frequency channel point",
        np.array(response.grids[0]),
        np.array([10.0]),
    )
    assert_close(
        "zero-frequency channel weight",
        np.array(response.data),
        np.array([0.5]),
    )


def test_overlapping_sidebands_with_asymmetric_bandpass():
    selector = pyarts.arts.SensorHeterodyneFrequencyRange()
    selector.filter(make_asymmetric_sideband_filter())
    selector.mix(10.0)

    inspect_case(
        "overlapping sidebands with asymmetric bandpass",
        selector,
        expected_global=[[10.0, 13.0], [10.0, 7.0]],
        expected_local=[[0.0, 3.0], [0.0, 3.0]],
        expected_affine=[[10.0, 1.0], [10.0, -1.0]],
    )

    upper_expected = np.array([db_to_lin(1.0), db_to_lin(5.0), db_to_lin(10.0)])
    lower_expected = np.full(3, db_to_lin(-10.0))

    assert_close(
        "overlapping sidebands upper-path local gains",
        selector.local_response(np.array([1.0, 2.0, 3.0]), 0),
        upper_expected,
    )
    assert_close(
        "overlapping sidebands lower-path local gains",
        selector.local_response(np.array([1.0, 2.0, 3.0]), 1),
        lower_expected,
    )

    responses = selector.channel_responses(
        [
            pyarts.arts.SensorDiracChannel(1.0),
            pyarts.arts.SensorDiracChannel(2.0),
            pyarts.arts.SensorDiracChannel(3.0),
        ]
    )

    expected_points = [
        np.array([9.0, 11.0]),
        np.array([8.0, 12.0]),
        np.array([7.0, 13.0]),
    ]
    expected_weights = [
        np.array([db_to_lin(-10.0) / 2.0, db_to_lin(1.0) / 2.0]),
        np.array([db_to_lin(-10.0) / 2.0, db_to_lin(5.0) / 2.0]),
        np.array([db_to_lin(-10.0) / 2.0, db_to_lin(10.0) / 2.0]),
    ]

    for index, response in enumerate(responses):
        assert_close(
            f"overlapping sidebands channel {index} points",
            np.array(response.grids[0]),
            expected_points[index],
        )
        assert_close(
            f"overlapping sidebands channel {index} weights",
            np.array(response.data),
            expected_weights[index],
        )


test_lowpass_then_lo_then_box_spectrometer()
test_default_positive_range_survives_mixes()
test_bandpass_lo_bandpass_lo_chain()
test_ideal_split_about_lo()
test_weighted_split_about_lo()
test_zero_frequency_is_not_double_counted()
test_overlapping_sidebands_with_asymmetric_bandpass()

print("\nAll heterodyne frequency-response checks passed.")