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


def quadrature(channel):
    """Frequencies and normalized integration weights of a channel"""
    q = channel.quadrature()
    w = np.array(q.data)
    return np.array(q.grids[0]), w / w.sum()


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
    selector = pyarts.arts.sensor.HeterodyneFrequencyRange()
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

    # The box maps to one piece, 16 - [2.6, 4.6], with the trapezoid weights of 5 equidistant points
    channel = pyarts.arts.sensor.BoxChannel(2.6, 4.6, 5)
    response = selector.channel_response(channel)

    assert len(response.response) == 1
    assert_close(
        "lowpass -> LO piece points",
        np.array(response.response[0].grids[0]),
        np.array([11.4, 11.9, 12.4, 12.9, 13.4]),
    )
    assert_close("lowpass -> LO piece response", np.array(response.response[0].data), np.ones(5))

    f, w = quadrature(response)
    assert_close("lowpass -> LO channel points", f, [11.4, 11.9, 12.4, 12.9, 13.4])
    assert_close("lowpass -> LO channel weights", w, [1 / 8, 1 / 4, 1 / 4, 1 / 4, 1 / 8])


def test_default_positive_range_survives_mixes():
    selector = pyarts.arts.sensor.HeterodyneFrequencyRange()
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

    # One delta per path, each with the full weight: the response is the sum over the sidebands
    response = selector.channel_response(pyarts.arts.sensor.DiracChannel(0.25))
    q = response.quadrature()

    assert_close(
        "default positive range dirac points after two LOs",
        np.array(q.grids[0]),
        np.array([8.75, 9.25, 10.75, 11.25]),
    )
    assert_close(
        "default positive range dirac weights after two LOs",
        np.array(q.data),
        np.ones(4),
    )


def test_bandpass_lo_bandpass_lo_chain():
    selector = pyarts.arts.sensor.HeterodyneFrequencyRange()
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

    # Two pieces, 9 + [0, 1] and 9 - [0, 1], each with trapezoid weights [1/4, 1/2, 1/4];
    # where they meet at 9 the weights add, as for one box over [8, 10]
    channel = pyarts.arts.sensor.BoxChannel(0.0, 1.0, 3)
    response = selector.channel_response(channel)

    assert len(response.response) == 2
    f, w = quadrature(response)
    assert_close("bandpass -> LO -> bandpass -> LO channel points", f, [8.0, 8.5, 9.0, 9.5, 10.0])
    assert_close("bandpass -> LO -> bandpass -> LO channel weights", w, [1 / 8, 1 / 4, 1 / 4, 1 / 4, 1 / 8])


def test_ideal_split_about_lo():
    selector = pyarts.arts.sensor.HeterodyneFrequencyRange()
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
    selector = pyarts.arts.sensor.HeterodyneFrequencyRange()
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

    channel = pyarts.arts.sensor.DiracChannel(2.0)
    q = selector.channel_response(channel).quadrature()

    assert_close(
        "weighted split channel points",
        np.array(q.grids[0]),
        np.array([8.0, 12.0]),
    )
    assert_close(
        "weighted split channel weights",
        np.array(q.data),
        np.array([0.6, 0.6]),
    )


def test_zero_frequency_receives_both_sidebands():
    # At zero intermediate frequency both sidebands are the LO itself.  The response is the sum over
    # the sidebands, so the LO has the weight of both: 2, as two deltas at any other frequency have
    selector = pyarts.arts.sensor.HeterodyneFrequencyRange()
    selector.bandpass(np.array([5.0, 15.0]))
    selector.mix(10.0)

    at_zero = selector.channel_response(pyarts.arts.sensor.DiracChannel(0.0)).quadrature()
    off_zero = selector.channel_response(pyarts.arts.sensor.DiracChannel(1.0)).quadrature()

    assert_close("zero-frequency channel point", np.array(at_zero.grids[0]), np.array([10.0]))
    assert_close("zero-frequency channel weight", np.array(at_zero.data), np.array([2.0]))
    assert_close("off-zero channel points", np.array(off_zero.grids[0]), np.array([9.0, 11.0]))
    assert_close("off-zero channel weights", np.array(off_zero.data), np.array([1.0, 1.0]))


def test_overlapping_sidebands_with_asymmetric_bandpass():
    selector = pyarts.arts.sensor.HeterodyneFrequencyRange()
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
            pyarts.arts.sensor.DiracChannel(1.0),
            pyarts.arts.sensor.DiracChannel(2.0),
            pyarts.arts.sensor.DiracChannel(3.0),
        ]
    )

    expected_points = [
        np.array([9.0, 11.0]),
        np.array([8.0, 12.0]),
        np.array([7.0, 13.0]),
    ]
    expected_weights = [
        np.array([db_to_lin(-10.0), db_to_lin(1.0)]),
        np.array([db_to_lin(-10.0), db_to_lin(5.0)]),
        np.array([db_to_lin(-10.0), db_to_lin(10.0)]),
    ]

    for index, response in enumerate(responses):
        q = response.quadrature()
        assert_close(
            f"overlapping sidebands channel {index} points",
            np.array(q.grids[0]),
            expected_points[index],
        )
        assert_close(
            f"overlapping sidebands channel {index} weights",
            np.array(q.data),
            expected_weights[index],
        )


def channel_mean_frequency(channel):
    """The channel's measurement of a radiance equal to the frequency: exact for the trapezoid rule"""
    f, w = quadrature(channel)
    return float(np.sum(w * f))


def test_sampling_density_does_not_weight_the_band():
    # A flat filter with frequencies only near the lower edge adds sampling points there.  The box
    # still responds evenly over [0, 100], so the radiance f is measured as 50 whatever the sampling.
    flat = pyarts.arts.SortedGriddedField1()
    grid = np.concatenate([np.linspace(0.0, 10.0, 11), [100.0]])
    flat.grids = (pyarts.arts.AscendingGrid(grid),)
    flat.data = pyarts.arts.Vector(np.ones(grid.size))
    flat.gridnames = ("frequency",)

    selector = pyarts.arts.sensor.HeterodyneFrequencyRange()
    selector.filter(flat)

    response = selector.channel_response(pyarts.arts.sensor.BoxChannel(0.0, 100.0, 3))
    assert_close("densely filtered box mean frequency", channel_mean_frequency(response), 50.0)


def test_single_sideband_region_is_not_double_weighted():
    # A box over intermediate frequencies [0, 50] with the LO at 30 sees [30, 80] in the upper
    # sideband and [0, 30] in the lower, where intermediate frequencies above 30 have no image.
    # The response is the sum over the sidebands, one over [0, 80], so f is measured as 40.
    selector = pyarts.arts.sensor.HeterodyneFrequencyRange()
    selector.mix(30.0)

    response = selector.channel_response(pyarts.arts.sensor.BoxChannel(0.0, 50.0, 11))
    assert_close("partial lower sideband mean frequency", channel_mean_frequency(response), 40.0)


test_lowpass_then_lo_then_box_spectrometer()
test_default_positive_range_survives_mixes()
test_bandpass_lo_bandpass_lo_chain()
test_ideal_split_about_lo()
test_weighted_split_about_lo()
test_zero_frequency_receives_both_sidebands()
test_overlapping_sidebands_with_asymmetric_bandpass()
test_sampling_density_does_not_weight_the_band()
test_single_sideband_region_is_not_double_weighted()

print("\nAll heterodyne frequency-response checks passed.")
