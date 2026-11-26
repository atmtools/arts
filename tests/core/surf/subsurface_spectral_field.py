import numpy as np
import pyarts3 as pyarts

# A field linear in all coordinates, so linear interpolation is exact inside the grids
freq = np.array([100e9, 200e9, 400e9])
alt = np.array([-3.0, -1.0, 0.0])
lat = np.array([-10.0, 10.0])
lon = np.array([-20.0, 0.0, 20.0])


def linear(f, a, b, c):
    return 1.0 + 2e-11 * f + 0.3 * a + 0.01 * b - 0.02 * c


F, A, B, C = np.meshgrid(freq, alt, lat, lon, indexing="ij")
gf = pyarts.arts.SpectralGeodeticField4(
    data=linear(F, A, B, C),
    grid_names=["frequency", "altitude", "latitude", "longitude"],
    grids=[freq, alt, lat, lon],
)

data = pyarts.arts.SubsurfaceSpectralData(gf)
assert data.data_type == "SpectralGeodeticField4"

pos = (-0.5, 3.0, -7.0)
fs = np.array([150e9, 200e9, 333e9])

# Values
for f in fs:
    assert np.isclose(data(f, *pos), linear(f, *pos), rtol=1e-14, atol=0)

# The frequency grid gives the same values as each frequency on its own
vals = data(pyarts.arts.AscendingGrid(fs), *pos)
assert all(vals[i] == data(fs[i], *pos) for i in range(len(fs)))

# Weights: frequency is the slowest dimension of the flat data
flat = gf.data.flatten()
nspatial = len(alt) * len(lat) * len(lon)


def dense(pairs, size):
    out = np.zeros(size)
    for i, w in pairs:
        out[i] += w
    return out


# The weights reproduce the values from the flat data, and zero weights are left out
for f in fs:
    ws = data.ws(f, *pos)
    assert len(ws) == (8 if f in freq else 16) and all(w != 0 for _, w in ws)
    assert np.isclose(dense(ws, flat.size) @ flat, data(f, *pos), rtol=1e-14, atol=0)

# On a data frequency, only the stored values of that data frequency contribute
for k, f in enumerate(freq):
    ws = data.ws(f, *pos)
    assert len(ws) == 8 and all(k * nspatial <= x < (k + 1) * nspatial for x, _ in ws)

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
data.alt_upp = "Zero"
assert all(data(pyarts.arts.AscendingGrid(fs), 1.0, 3.0, -7.0) == 0.0)
assert data.ws(280e9, 1.0, 3.0, -7.0) == []
data.alt_upp = "None"

# A single value
num = pyarts.arts.SubsurfaceSpectralData(0.01)
assert num(280e9, *pos) == 0.01
assert num.ws(280e9, *pos) == [(0, 1.0)]

# A single data frequency: every frequency has the same weights
one = pyarts.arts.SubsurfaceSpectralData(
    pyarts.arts.SpectralGeodeticField4(
        data=linear(F, A, B, C)[1:2],
        grid_names=["frequency", "altitude", "latitude", "longitude"],
        grids=[freq[1:2], alt, lat, lon],
    )
)
assert len(one.ws(fs[0], *pos)) == 8
assert all(one.ws(f, *pos) == one.ws(fs[0], *pos) for f in fs)
assert all(one(f, *pos) == one(200e9, *pos) for f in fs)

# A function of frequency and position has no stored values
fun = pyarts.arts.SubsurfaceSpectralData(
    pyarts.arts.NumericQuaternaryOperator(lambda f, a, b, c: linear(f, a, b, c))
)
assert fun(280e9, *pos) == linear(280e9, *pos)
assert fun.ws(280e9, *pos) == []

# The field
field = pyarts.arts.SubsurfaceSpectralField()
field["loss tangent"] = data
field["refractive index"] = 1.7
assert "loss tangent" in field and "refractive index" in field
assert "density" not in field
assert len(field.keys()) == 2
assert all(isinstance(k, pyarts.arts.SubsurfaceSpectralPropertyTag) for k in field.keys())
assert pyarts.arts.SubsurfaceSpectralPropertyTag("loss tangent") in field
assert field["refractive index"](280e9, *pos) == 1.7
try:
    field["density"]
    raise AssertionError("Expected a KeyError")
except KeyError:
    pass

# The parameters at one position, on a frequency grid
point = field(pyarts.arts.AscendingGrid(fs), *pos)
assert isinstance(point, pyarts.arts.SubsurfaceSpectralPoint)
assert len(point.keys()) == 2 and all(isinstance(k, pyarts.arts.SubsurfaceSpectralPropertyTag) for k in point.keys())
for key in field.keys():
    assert key in point
    assert all(np.asarray(point[key]) == np.asarray(field[key](pyarts.arts.AscendingGrid(fs), *pos)))
try:
    point["density"]
    raise AssertionError("Expected a KeyError")
except KeyError:
    pass
point.savexml("SubsurfaceSpectralPoint.xml")
loaded_point = pyarts.arts.SubsurfaceSpectralPoint.fromxml("SubsurfaceSpectralPoint.xml")
assert all(np.array_equal(loaded_point[k], point[k]) for k in point.keys())

# XML round trip
field.savexml("subsurface_spectral_field.xml")
loaded = pyarts.arts.SubsurfaceSpectralField.fromxml("subsurface_spectral_field.xml")
assert loaded["loss tangent"](333e9, *pos) == field["loss tangent"](333e9, *pos)
assert loaded["refractive index"](333e9, *pos) == 1.7
