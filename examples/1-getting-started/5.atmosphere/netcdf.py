import pyarts
import numpy as np

# Use defaul
ws = pyarts.Workspace()
ws.absorption_speciesSet(species=["O2", "CO2"])
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)

# The atmospheric field must be uniformly gridded to be saved to a netcdf file
gridded_atm = ws.atmospheric_field.as_gridded(np.linspace(0, 100e3, 101), [0], [0])

# Save only O2 VMR, temperature and pressure to a netcdf file
gridded_atm.to_xarray(keys=["O2", 't', 'p']).to_netcdf("atmosphere.nc")

# Use helper method to read the netcdf file (if it is in any path seen by ARTS)
atmnc = pyarts.data.xarray_open_dataset("atmosphere.nc")

# Convert the xarray dataset to an atmospheric field
atm = pyarts.data.to_atmospheric_field(atmnc)
