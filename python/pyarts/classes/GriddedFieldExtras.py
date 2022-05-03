import copy
import pyarts.pyarts_cpp as cxx

import xarray
import numpy as np
from scipy import interpolate

def extract_slice(g, s=slice(None), axis=0):
    """Return a new GriddedField containing a slice of the current one.
    Parameters:
        s (slice): Slice.
        axis (int): Axis to slice along.
    Returns:
        GriddedField containing sliced grids and data.
    """
    g.checksize_strict()
    
    gf = copy.deepcopy(g)
    axis_grid = gf.get_grid(axis)
    if isinstance(axis_grid, cxx.Vector):
        gf.set_grid(axis, axis_grid[s])
    else:
        gf.set_grid(axis, type(axis_grid)(list(axis_grid)[s]))
    slices = [slice(None)] * gf.dim
    slices[axis] = s
    gf.data = gf.data[tuple(slices)]

    return gf


# FIXME: I have not managed to forward kwargs reliably from C++ to python  // Richard 2022-04-27
def refine_grid(gin, new_grid, axis=0, type="linear"):
    """Interpolate GriddedField axis to a new grid.
    This function replaces a grid of a GriddField and interpolates all
    data to match the new coordinates. :func:`scipy.interpolate.interp1d`
    is used for interpolation.
    Parameters:
        new_grid (ndarray): The coordinates of the interpolated values.
        axis (int): Specifies the axis of data along which to interpolate.
            Interpolation defaults to the first axis of the GriddedField.
        type (str or function): Rescaling type for function if str or rescaling function
    Returns: GriddedField
    """
    if type == "linear":
        fun = np.array
    elif type == "log10":
        fun = np.log10
    elif type == "log":
        fun = np.log
    else:
        fun = type
        
    g = copy.deepcopy(gin)
    
    if len(g.get_grid(axis)) > 1:
        f = interpolate.interp1d(fun(g.get_grid(axis)), g.data, axis=axis)  # FIXME: no **kwargs...
        g.set_grid(axis, new_grid)
        g.data = f(fun(new_grid))
    else:  # if the intention is to create a useful TensorX
        g.data = g.data.repeat(len(new_grid), axis=axis)
        g.data = f(fun(new_grid))

    g.checksize_strict()

    return g


def to_xarray(g):
    """Convert GriddedField to xarray.DataArray object.
    Convert a GriddedField object into a :func:`xarray.DataArray`
    object.  The dataname is used as the DataArray name.
    Returns:
        xarray.DataArray object corresponding to gridded field
    """

    da = xarray.DataArray(g.data)
    da = da.rename(dict((k, v)
        for (k, v) in zip(da.dims, [g.get_grid_name(i) for i in range(g.dim)])
        if v!=""))
    da = da.assign_coords(
        **{name: coor
            for (name, coor) in zip(da.dims, [np.array(g.get_grid(i)) for i in range(g.dim)])
            if len(coor)>0})
    if g.name: da.attrs['data_name'] = g.name
    return da


def from_xarray(cls, da):
    """Create GriddedField from a xarray.DataArray object.
    The data and its dimensions are returned as a :class:`GriddedField` object.
    The DataArray name is used as name for the gridded field. If the attribute
    `data_name` is present, it is used as `dataname` on the :class:`GriddedField`.
    Parameters:
        da (xarray.DataArray): xarray.DataArray containing the dimensions and data.
    Returns:
        GriddedField object.
    """
    obj = cls()
    for i in range(obj.dim):
        obj.set_grid(i, da[da.dims[i]].values)
        obj.set_grid_name(i, da.dims[i])
    if da.values.ndim != obj.dim:
        raise RuntimeError(f"Dimension mismatch: Expected {obj.dim} got {da.values.ndim}")
    obj.data = da.values
    obj.name = str(da.attrs.get('data_name', 'Data'))
    obj.checksize_strict()
    return obj


getattr(cxx, "GriddedField::details").extract_slice = extract_slice
getattr(cxx, "GriddedField::details").refine_grid = refine_grid
getattr(cxx, "GriddedField::details").from_xarray = from_xarray
getattr(cxx, "GriddedField::details").to_xarray = to_xarray
