This demonstrates how to set your atmospheric field using netcdf.

An atmospheric field is created using standard workspace method
flow.  The atmospheric field is then gridded to a uniform grid
before being written to a netcdf file via the :class:`xarray.Dataset`
format.

An helper method
