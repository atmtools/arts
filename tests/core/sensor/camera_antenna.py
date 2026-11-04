"""
A sensor built from a thin-lens CameraAntenna looks where measurement_sensorAddCamera looks

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


def nonzero_sample(obsel):
    """ The line of sight, frequency and Stokes I weight of an obsel with one nonzero weight """
    weights = np.asarray(obsel.weight_matrix.dense)
    ip, jf = np.argwhere(weights[..., 0] != 0.0)[0]
    return np.asarray(obsel.poslos)[ip, 3:], np.asarray(obsel.f_grid)[jf], weights[ip, jf, 0]


for los in ([90.0, 0.0], [180.0, 0.0], [130.0, 40.0]):
    ws = pa.Workspace()
    ws.freq_grid = FREQ
    ws.measurement_sensorInit()
    ws.measurement_sensorAddCamera(pos=POS,
                                   los=los,
                                   n_h=N_H,
                                   n_w=N_W,
                                   ccd_h=CCD_H,
                                   ccd_w=CCD_W,
                                   focal_length=FOCAL_LENGTH,
                                   aperture_diameter=0.025,
                                   focus_distance=FOCUS_DISTANCE)

    camera = pa.arts.sensor.CameraAntenna.thin_lens(N_H, N_W, CCD_H, CCD_W, FOCAL_LENGTH, FOCUS_DISTANCE)
    builder = pa.arts.SensorBuilder([pa.arts.sensor.DiracChannel(f) for f in FREQ], camera)
    obsels, meta = builder([pa.arts.Vector3(POS)], [pa.arts.Vector2(los)], ELL)

    assert len(obsels) == len(ws.measurement_sensor) == N_H * N_W * len(FREQ), f"element count for {los}"
    for i, (got, expected) in enumerate(zip(obsels, ws.measurement_sensor)):
        got_los, got_f, got_w = nonzero_sample(got)
        exp_los, exp_f, exp_w = nonzero_sample(expected)
        np.testing.assert_allclose(got_los, exp_los, atol=1e-10, rtol=0.0, err_msg=f"line of sight {i} at {los}")
        assert got_f == exp_f, f"frequency of element {i} at {los}"
        assert got_w == exp_w, f"weight of element {i} at {los}"

    assert len(meta) == len(ws.measurement_sensor_meta) == 1, f"one image at {los}"
    got_image = meta[0].data
    exp_image = ws.measurement_sensor_meta[0].data
    for k in (0, 1):
        np.testing.assert_allclose(got_image.grids[k], exp_image.grids[k], atol=1e-12, rtol=0.0,
                                   err_msg=f"pixel offsets {k} at {los}")
    assert np.asarray(got_image.data).shape == np.asarray(exp_image.data).shape, f"image shape at {los}"

pa.plots.plot(camera)
