from datetime import datetime
import numpy as np
import xarray

import pyarts.pyarts_cpp as cxx


def xsecrecord_to_xarray(self):
    """Convert XsecRecord to xarray dataset"""
    attrs = {'creation_date': str(datetime.now())}
    coords = {
        "coeffs": ("coeffs", range(4), {
            "unit": "1",
            "long_name": "Coefficients p00, p10, p01, p20"
        }),
        "bands": ("bands", range(len(self.fitcoeffs)), {
            "unit": "1",
            "long_name": "Band indices"
        })
    }

    data_vars = {}
    data_vars["fitminpressures"] = ("bands", self.fitminpressures, {
        "unit": "Pa",
        "long_name": "Mininum pressures from fit"
    })
    data_vars["fitmaxpressures"] = ("bands", self.fitmaxpressures, {
        "unit": "Pa",
        "long_name": "Maximum pressures from fit"
    })
    data_vars["fitmintemperatures"] = ("bands", self.fitmintemperatures, {
        "unit": "K",
        "long_name": "Minimum temperatures from fit"
    })
    data_vars["fitmaxtemperatures"] = ("bands", self.fitmaxtemperatures, {
        "unit": "K",
        "long_name": "Maximum temperatures from fit"
    })

    for i, band in enumerate(self.fitcoeffs):
        coords[f"band{i}_fgrid"] = (f"band{i}_fgrid", band.get_grid(0), {
            "unit": "Hz",
            "long_name": f"Frequency grid for band {i}"
        })
        data_vars[f"band{i}_coeffs"] = ([f"band{i}_fgrid", "coeffs"], band.data, {
            'unit': "m",
            'long_name': f'Fit coefficients for band {i}'
        })

    attrs["species"] = str(self.species)
    attrs["version"] = self.version

    return xarray.Dataset(data_vars, coords, attrs)


def xsecrecord_from_xarray(cls, ds):
    """Set XsecRecord from xarray dataset"""
    xr = cls()
    coeff_grid = cxx.ArrayOfString(["p00", "p10", "p01", "p20"])
    gfs = []
    xr.species = ds.attrs["species"]
    xr.version = ds.attrs["version"]
    xr.fitminpressures = ds["fitminpressures"].values
    xr.fitmaxpressures = ds["fitmaxpressures"].values
    xr.fitmintemperatures = ds["fitmintemperatures"].values
    xr.fitmaxtemperatures = ds["fitmaxtemperatures"].values
    for i in ds["bands"].values:
        g = cxx.GriddedField2()
        for j, gridname in enumerate(["frequency grid [Hz]", "fit coefficients [m]"]):
            g.set_grid_name(j, gridname)
        g.set_grid(0, ds[f"band{i}_fgrid"].values)
        g.set_grid(1, coeff_grid)
        g.data = ds[f"band{i}_coeffs"].values
        gfs.append(g)
    xr.fitcoeffs = cxx.ArrayOfGriddedField2(gfs)

    return xr


def xsecrecord_to_netcdf(self, filename):
    """Save XsecRecord to NetCDF file"""
    self.to_xarray().to_netcdf(filename)


def xsecrecord_from_netcdf(cls, filename):
    """Read XsecRecord from NetCDF file"""
    return cls.from_xarray(xarray.load_dataset(filename))


def xsecrecord__eq__(self, other):
    if isinstance(other, cxx.XsecRecord) and \
            self.version == other.version and \
            str(self.species) == str(other.species) and \
            np.all(self.fitminpressures == other.fitminpressures) and \
            np.all(self.fitmaxpressures == other.fitmaxpressures) and \
            np.all(self.fitmintemperatures == other.fitmintemperatures) and \
            np.all(self.fitmaxtemperatures == other.fitmaxtemperatures) and \
            np.all([np.all(np.array(i.data) == np.array(j.data)) for i, j in zip(self.fitcoeffs, other.fitcoeffs)]):
        return True
    else:
        return False


getattr(cxx, "XsecRecord::details").to_xarray = xsecrecord_to_xarray
getattr(cxx, "XsecRecord::details").from_xarray = xsecrecord_from_xarray
getattr(cxx, "XsecRecord::details").to_netcdf = xsecrecord_to_netcdf
getattr(cxx, "XsecRecord::details").from_netcdf = xsecrecord_from_netcdf
getattr(cxx, "XsecRecord::details").__eq__ = xsecrecord__eq__
