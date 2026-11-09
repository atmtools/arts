"""
A sensor built from a thin-lens CameraAntenna looks where measurement_sensorsAddCamera looks

Both describe the same pixels on the same CCD, so for every pixel and frequency the observation element of
the builder must sample the same line of sight, at the same frequency, with the same weight, in the same
order, and the image metadata must have the same pixel offsets.
"""

import numpy as np
import pyarts3 as pa

POS = [600e3, 10.0, 20.0]
FREQ = [100e9, 200e9]
ELL = pa.arts.planets.Earth.ellipsoid
N_H, N_W = 3, 4
CCD_H, CCD_W = 0.03, 0.02
FOCAL_LENGTH, FOCUS_DISTANCE = 0.05, 100.0


def nonzero_sample(sensor, i):
    """ The line of sight, frequency and Stokes I weight of element i of a sensor, which has one nonzero weight """
    weights = np.asarray(sensor.elements[i].dense)
    ip, jf = np.argwhere(weights[..., 0] != 0.0)[0]
    return np.asarray(sensor.geometry)[ip, 3:], np.asarray(sensor.freq_grid)[jf], weights[ip, jf, 0]


for los in ([90.0, 0.0], [180.0, 0.0], [130.0, 40.0]):
    ws = pa.Workspace()
    ws.freq_grid = FREQ
    ws.measurement_sensorsInit()
    ws.measurement_sensorsAddCamera(pos=POS,
                                   los=los,
                                   n_h=N_H,
                                   n_w=N_W,
                                   ccd_h=CCD_H,
                                   ccd_w=CCD_W,
                                   focal_length=FOCAL_LENGTH,
                                   aperture_diameter=0.025,
                                   focus_distance=FOCUS_DISTANCE)

    camera = pa.arts.sensors.CameraAntenna.thin_lens(N_H, N_W, CCD_H, CCD_W, FOCAL_LENGTH, FOCUS_DISTANCE)
    builder = pa.arts.SensorBuilder([pa.arts.sensors.DiracChannel(f) for f in FREQ], camera)
    sensor = builder([pa.arts.Vector3(POS)], [pa.arts.Vector2(los)], ELL)
    added = ws.measurement_sensors[0]

    assert len(ws.measurement_sensors) == 1, f"one sensor for {los}"
    assert len(sensor.elements) == len(added.elements) == N_H * N_W * len(FREQ), f"element count for {los}"
    for i in range(len(sensor.elements)):
        got_los, got_f, got_w = nonzero_sample(sensor, i)
        exp_los, exp_f, exp_w = nonzero_sample(added, i)
        np.testing.assert_allclose(got_los, exp_los, atol=1e-10, rtol=0.0, err_msg=f"line of sight {i} at {los}")
        assert got_f == exp_f, f"frequency of element {i} at {los}"
        assert got_w == exp_w, f"weight of element {i} at {los}"

    # The builder's image has the placement as its first axis
    got_image = sensor.meta
    exp_image = added.meta
    for k in (0, 1):
        np.testing.assert_allclose(got_image.grids[k + 1], exp_image.grids[k], atol=1e-12, rtol=0.0,
                                   err_msg=f"pixel offsets {k} at {los}")
    assert np.asarray(got_image.data).shape == (1, *np.asarray(exp_image.data).shape), f"image shape at {los}"

pa.plots.plot(camera)
