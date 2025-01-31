import pyarts

# create empty
x0 = pyarts.arts.FastLinearGasOptics()

# create with data
x1 = pyarts.arts.FastLinearGasOptics(offset=[1], slopeT=[2], slopeP=[3])

# create empty
y0 = pyarts.arts.FastNonLinearGasOptics()

# create with some field
y1 = pyarts.arts.FastNonLinearGasOptics(linear=x0, nonlinear=x1)

# Create general structure
z = pyarts.arts.FastGasOptics()


z["H2O"] = y0
z["O3"] = y1
z["O2"] = x0
z["NO"] = x1
