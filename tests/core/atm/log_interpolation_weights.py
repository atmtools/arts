import numpy as np
import pyarts3 as pyarts

# With log interpolation, the value is exp(sum w ln x), so its derivative by a
# stored x is w v / x.  The flat weights are checked against finite differences.

freq = np.array([100e9, 200e9, 400e9])
alt = np.array([0.0, 1e3, 3e3])
lat = np.array([-10.0, 10.0])
lon = np.array([-20.0, 0.0, 20.0])


SCALE = {"frequency": 2e-12, "altitude": 3e-4, "latitude": 0.01, "longitude": -0.02}


def positive(*x, dims):
    """Exponential in all coordinates, so log interpolation is exact inside the grids"""
    return np.exp(0.5 + sum(SCALE[d] * v for d, v in zip(dims, x)))


def check(make, data, pos, dims):
    d = make(data)
    d.log_interpolation = True
    v = d(*pos)
    assert np.isclose(v, positive(*pos, dims=dims), rtol=1e-12, atol=0)

    ws = d.ws(*pos)
    assert len(ws) > 0
    flat = data.flatten()
    for j, w in ws:
        h = 1e-7 * flat[j]
        p = data.copy()
        p.flat[j] += h
        q = make(p)
        q.log_interpolation = True
        assert np.isclose((q(*pos) - v) / h, w, rtol=1e-5, atol=0), (j, w, (q(*pos) - v) / h)


# The atmosphere field data
A, B, C = np.meshgrid(alt, lat, lon, indexing="ij")
names = ["altitude", "latitude", "longitude"]
check(lambda x: pyarts.arts.AtmData(pyarts.arts.GeodeticField3(data=x, grid_names=names, grids=[alt, lat, lon])),
      positive(A, B, C, dims=names), (500.0, 3.0, -7.0), names)

# The spectral fields data
F, A, B, C = np.meshgrid(freq, alt, lat, lon, indexing="ij")
names = ["frequency", "altitude", "latitude", "longitude"]
for T in [pyarts.arts.AtmSpectralData, pyarts.arts.SubsurfaceSpectralData]:
    check(lambda x: T(pyarts.arts.SpectralGeodeticField4(data=x, grid_names=names, grids=[freq, alt, lat, lon])),
          positive(F, A, B, C, dims=names), (150e9, 500.0, 3.0, -7.0), names)

F, B, C = np.meshgrid(freq, lat, lon, indexing="ij")
sdims = ["frequency", "latitude", "longitude"]
check(lambda x: pyarts.arts.SurfaceSpectralData(
    pyarts.arts.SpectralGeodeticField3(data=x, grid_names=["frequency", "latitude", "longitude"], grids=[freq, lat, lon])),
      positive(F, B, C, dims=sdims), (150e9, 3.0, -7.0), sdims)

# The frequency grid gives the same log-interpolated values as each frequency on its own
d = pyarts.arts.SurfaceSpectralData(
    pyarts.arts.SpectralGeodeticField3(data=positive(F, B, C, dims=sdims), grid_names=["frequency", "latitude", "longitude"],
                                       grids=[freq, lat, lon]))
d.log_interpolation = True
fs = [150e9, 200e9, 333e9]
vals = d(pyarts.arts.AscendingGrid(fs), 3.0, -7.0)
assert all(vals[i] == d(f, 3.0, -7.0) for i, f in enumerate(fs))
