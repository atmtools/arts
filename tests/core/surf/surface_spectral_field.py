import numpy as np
import pyarts3 as pyarts

# A field linear in all coordinates, so linear interpolation is exact inside the grids
freq = np.array([100e9, 200e9, 400e9])
lat = np.array([-10.0, 10.0])
lon = np.array([-20.0, 0.0, 20.0])


def linear(f, b, c):
    return 1.0 + 2e-11 * f + 0.01 * b - 0.02 * c


F, B, C = np.meshgrid(freq, lat, lon, indexing="ij")
gf = pyarts.arts.SpectralGeodeticField3(
    data=linear(F, B, C),
    grid_names=["frequency", "latitude", "longitude"],
    grids=[freq, lat, lon],
)

data = pyarts.arts.SurfaceSpectralData(gf)
assert data.data_type == "SpectralGeodeticField3"

pos = (3.0, -7.0)
fs = np.array([150e9, 200e9, 333e9])

# Values
for f in fs:
    assert np.isclose(data(f, *pos), linear(f, *pos), rtol=1e-14, atol=0)

# The frequency grid gives the same values as each frequency on its own
vals = data(pyarts.arts.AscendingGrid(fs), *pos)
assert all(vals[i] == data(fs[i], *pos) for i in range(len(fs)))

# Weights: frequency is the slowest dimension of the flat data
flat = gf.data.flatten()
nspatial = len(lat) * len(lon)


def dense(pairs, size):
    out = np.zeros(size)
    for i, w in pairs:
        out[i] += w
    return out


# The weights reproduce the values from the flat data, and zero weights are left out
for f in fs:
    ws = data.ws(f, *pos)
    assert len(ws) == (4 if f in freq else 8) and all(w != 0 for _, w in ws)
    assert np.isclose(dense(ws, flat.size) @ flat, data(f, *pos), rtol=1e-14, atol=0)

# On a data frequency, only the stored values of that data frequency contribute
for k, f in enumerate(freq):
    ws = data.ws(f, *pos)
    assert len(ws) == 4 and all(k * nspatial <= x < (k + 1) * nspatial for x, _ in ws)

# Extrapolation in frequency
data.frq_upp = "None"
try:
    data(500e9, *pos)
    raise AssertionError("Expected an error above the frequency grid")
except Exception as e:
    assert "frequency" in str(e)

data.frq_upp = "Zero"
assert data(500e9, *pos) == 0.0
assert data.ws(500e9, *pos) == []

data.frq_upp = "Nearest"
assert data(500e9, *pos) == data(400e9, *pos)
assert data.ws(500e9, *pos) == data.ws(400e9, *pos)

data.frq_upp = "Linear"
assert np.isclose(data(500e9, *pos), linear(500e9, *pos), rtol=1e-14, atol=0)
ws = data.ws(500e9, *pos)
assert np.isclose(dense(ws, flat.size) @ flat, data(500e9, *pos), rtol=1e-14, atol=0)

# Extrapolation in space
data.lat_upp = "Zero"
assert all(data(pyarts.arts.AscendingGrid(fs), 20.0, -7.0) == 0.0)
assert data.ws(280e9, 20.0, -7.0) == []
data.lat_upp = "None"

# A single value
num = pyarts.arts.SurfaceSpectralData(0.01)
assert num(280e9, *pos) == 0.01
assert num.ws(280e9, *pos) == [(0, 1.0)]

# A single data frequency: every frequency has the same weights
one = pyarts.arts.SurfaceSpectralData(
    pyarts.arts.SpectralGeodeticField3(
        data=linear(F, B, C)[1:2],
        grid_names=["frequency", "latitude", "longitude"],
        grids=[freq[1:2], lat, lon],
    )
)
assert len(one.ws(fs[0], *pos)) == 4
assert all(one.ws(f, *pos) == one.ws(fs[0], *pos) for f in fs)
assert all(one(f, *pos) == one(200e9, *pos) for f in fs)

# A function of frequency and position has no stored values
fun = pyarts.arts.SurfaceSpectralData(pyarts.arts.NumericTernaryOperator(lambda f, b, c: linear(f, b, c)))
assert fun(280e9, *pos) == linear(280e9, *pos)
assert fun.ws(280e9, *pos) == []

# The field
field = pyarts.arts.SurfaceSpectralField()
field["depolarization"] = data
field["refractive index"] = 1.7
assert "depolarization" in field and "refractive index" in field
assert "density" not in field
assert len(field.keys()) == 2
assert all(isinstance(k, pyarts.arts.SurfaceSpectralPropertyTag) for k in field.keys())
assert pyarts.arts.SurfaceSpectralPropertyTag("depolarization") in field
assert field["refractive index"](280e9, *pos) == 1.7
try:
    field["density"]
    raise AssertionError("Expected a KeyError")
except KeyError:
    pass

# XML round trip
field.savexml("surface_spectral_field.xml")
loaded = pyarts.arts.SurfaceSpectralField.fromxml("surface_spectral_field.xml")
assert loaded["depolarization"](333e9, *pos) == field["depolarization"](333e9, *pos)
assert loaded["refractive index"](333e9, *pos) == 1.7
