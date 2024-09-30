from pyarts.arts import interp as cxx

import numpy as np
from enum import Enum


def interp(data, *lags):
    """
    Interpolates the given data using the specified lags.

    Args:
        data: The data to be interpolated.
        *lags: Variable number of lags to be used for interpolation.

    Returns:
        The interpolated data.

    """
    return cxx.interp(data, *lags)


def right_reinterp(data, *lags):
    """
    Reinterprets the data array using interpolation for each left-most index.

    Parameters:
    data (numpy.ndarray): The input data array.
    *lags (float): Variable number of lag values used for interpolation.

    Returns:
    numpy.ndarray: The reinterpreted data array.

    """
    return np.array([interp(data[i], *lags) for i in range(data.shape[0])])


class grid_type(Enum):
    """Enumeration for different types of interpolation grids."""
    linear = 1
    cyclic_0to360 = 2
    cyclic_pm180 = 3


class regular_grid:
    """
    A class representing a regular grid.

    Args:
        grid (list): The grid points.
        func (function, optional): A function to apply to the input values. Defaults to identity function.
        type (str, optional): The type of the grid. Defaults to "linear".

    Attributes:
        grid (list): The grid points.
        type (GridType): The type of the grid.
        func (function): The function to apply to the input values.

    Methods:
        __len__(): Returns the length of the grid.
        iw(x, **kwargs): Interpolates the input value on the grid.

    """

    def __init__(self, grid, *, func=lambda x: x, type="linear"):
        self.grid = grid
        self.type = grid_type[type]
        self.func = func

    def __len__(self):
        return len(self.grid)

    def iw(self, x, **kwargs):
        """
        Interpolate the given values `x` onto the grid defined by `self.grid`.

        Parameters:
            x (array-like): The values to be interpolated.
            **kwargs: Additional keyword arguments.

        Returns:
            tuple: A tuple containing two elements:
                - A list of two indices representing the interpolated positions.
                - The interpolated values on the grid.

        Raises:
            None.
        """
        x = self.func(x)

        match self.type:
            case grid_type.linear:
                lag = cxx.PolyGrid(x, self.grid)
                ipos = lag.pos * 1
                lag.pos = 0
                return [ipos, ipos + 1], lag
            case grid_type.cyclic_0to360:
                lag = cxx.CyclicGrid0to360(x, self.grid)
                ipos = lag.pos * 1
                if ipos == len(self) - 1:
                    lag.pos = 1
                    return [ipos, 0], lag
                else:
                    lag.pos = 0
                    return [ipos, ipos + 1], lag
            case grid_type.cyclic_pm180:
                lag = cxx.CyclicGridPM180(x, self.grid)
                ipos = lag.pos * 1
                if ipos == len(self) - 1:
                    lag.pos = 1
                    return [ipos, 0], lag
                else:
                    lag.pos = 0
                    return [ipos, ipos + 1], lag


class computed_grid:
    """
    Represents a computed grid of atmospheric data.

    Args:
        lat_grid (Grid): The latitude grid.
        lon_grid (Grid): The longitude grid.
        *data: Variable length arguments representing the data arrays.
        func (function): The function used to compute the grid.
        type (str, optional): The type of grid. Defaults to "linear".
        **func_kwargs: Keyword arguments to be passed to the function.

    Attributes:
        lat_grid (Grid): The latitude grid.
        lon_grid (Grid): The longitude grid.
        data (tuple): The data arrays.
        type (str): The type of grid.
        func (function): The function used to compute the grid.
        func_kwargs (dict): The keyword arguments passed to the function.

    Raises:
        AssertionError: If the altitude shape of the data arrays does not match.
        AssertionError: If the latitude shape of the data arrays does not match.
        AssertionError: If the longitude shape of the data arrays does not match.
    """

    def __init__(
        self, lat_grid, lon_grid, *data, func, type="linear", **func_kwargs
    ):
        self.lat_grid = lat_grid
        self.lon_grid = lon_grid
        self.data = data
        self.type = type
        self.func = func
        self.func_kwargs = func_kwargs

        assert all(
            x.shape[0] == len(self) for x in self.data
        ), "Mismatch altitude shape of data"

        assert all(
            x.shape[1] == len(self.lat_grid) for x in self.data
        ), "Mismatch latitude shape of data"

        assert all(
            x.shape[2] == len(self.lon_grid) for x in self.data
        ), "Mismatch longitude shape of data"

    def __len__(self):
        return self.data[0].shape[0]

    def iw(self, x, **iw_kwargs):
            """
            Interpolate the field to a given point or set of points.

            Parameters:
            - x: The point or set of points to interpolate the field to.
            - iw_kwargs: Additional keyword arguments for the interpolation.

            Returns:
            - See `func:regular_grid.iw`.
            """
            ilat, laglat = self.lat_grid.iw(iw_kwargs["lat"])
            ilon, laglon = self.lon_grid.iw(iw_kwargs["lon"])

            grid = self.func(
                *[
                    right_reinterp(data[:, ilat, ilon], laglat, laglon)
                    for data in self.data
                ],
                **self.func_kwargs,
            )

            return regular_grid(grid, type=self.type).iw(x, **iw_kwargs)


class xarray_computed_grid:
    """
    Represents a computed grid using xarray.

    Args:
        lat_grid (Grid): The latitude grid.
        lon_grid (Grid): The longitude grid.
        *data: Variable-length arguments representing the data arrays.
        func (function): The function used to compute the grid.
        type (str, optional): The type of grid interpolation. Defaults to "linear".
        **func_kwargs: Keyword arguments to be passed to the function.

    Raises:
        AssertionError: If the altitude shape of the data arrays does not match.
        AssertionError: If the latitude shape of the data arrays does not match.
        AssertionError: If the longitude shape of the data arrays does not match.
    """

    def __init__(
        self, lat_grid, lon_grid, *data, func, type="linear", **func_kwargs
    ):
        self.lat_grid = lat_grid
        self.lon_grid = lon_grid
        self.data = data
        self.type = type
        self.func = func
        self.func_kwargs = func_kwargs

        assert all(
            x.shape[0] == len(self) for x in self.data
        ), "Mismatch altitude shape of data"

        assert all(
            x.shape[1] == len(self.lat_grid) for x in self.data
        ), "Mismatch latitude shape of data"

        assert all(
            x.shape[2] == len(self.lon_grid) for x in self.data
        ), "Mismatch longitude shape of data"

    def __len__(self):
        return self.data[0].shape[0]

    def iw(self, x, **iw_kwargs):
            """
            Interpolate the field to a given point or set of points.

            Parameters:
            - x: The point or set of points to interpolate the field to.
            - **iw_kwargs: Additional keyword arguments for the interpolation.

            Returns:
            - See `func:regular_grid.iw`.
            """
            ilat, laglat = self.lat_grid.iw(iw_kwargs["lat"])
            ilon, laglon = self.lon_grid.iw(iw_kwargs["lon"])

            grid = self.func(
                *[
                    right_reinterp(data[:, ilat, ilon].data, laglat, laglon)
                    for data in self.data
                ],
                **self.func_kwargs,
            )

            return regular_grid(grid, type=self.type).iw(x, **iw_kwargs)


class regular_gridded_data:
    """
    A class representing regular gridded data.

    Parameters:
    - data: numpy.ndarray
        The data array with shape (n_alt, n_lat, n_lon).
    - alt_grid:
        The altitude grid.  Must be a grid-like object.
    - lat_grid:
        The latitude grid.  Must be a grid-like object.
    - lon_grid:
        The longitude grid.  Must be a grid-like object.
    - func: callable, optional
        A function to apply to the interpolated data. Default is identity function.

    Attributes:
    - data: numpy.ndarray
        The data array.
    - alt_grid:
        The altitude grid.  Must be a grid-like object.
    - lat_grid:
        The latitude grid.  Must be a grid-like object.
    - lon_grid:
        The longitude grid.  Must be a grid-like object.
    - func: callable
        The function to apply to the interpolated data.

    Methods:
    - __call__(self, alt, lat, lon)
        Interpolates the data at the given altitude, latitude, and longitude.

    """

    def __init__(
        self, data, alt_grid, lat_grid, lon_grid, *, func=lambda x: x
    ):
        self.data = data
        self.alt_grid = alt_grid
        self.lat_grid = lat_grid
        self.lon_grid = lon_grid
        self.func = func

        assert self.data.shape == (len(alt_grid), len(lat_grid), len(lon_grid))

    def __call__(self, alt, lat, lon):
        ialt, lagalt = self.alt_grid.iw(alt, alt=alt, lat=lat, lon=lon)
        ilat, laglat = self.lat_grid.iw(lat, alt=alt, lat=lat, lon=lon)
        ilon, laglon = self.lon_grid.iw(lon, alt=alt, lat=lat, lon=lon)

        return self.func(
            interp(self.data[ialt, ilat, ilon], lagalt, laglat, laglon)
        )


class xarray_regular_gridded_data:
    """
    Represents a regular gridded data using xarray.

    Parameters:
    - data: The gridded data array.
    - alt_grid: The altitude grid.
    - lat_grid: The latitude grid.
    - lon_grid: The longitude grid.
    - func: Optional function to apply to the data.

    Attributes:
    - data: The gridded data array.
    - alt_grid: The altitude grid.
    - lat_grid: The latitude grid.
    - lon_grid: The longitude grid.
    - func: The function applied to the data.

    Methods:
    - __call__(self, alt, lat, lon): Interpolates the data at the given altitude, latitude, and longitude.

    """

    def __init__(
        self, data, alt_grid, lat_grid, lon_grid, *, func=lambda x: x
    ):
        self.data = data
        self.alt_grid = alt_grid
        self.lat_grid = lat_grid
        self.lon_grid = lon_grid
        self.func = func

        assert self.data.shape == (len(alt_grid), len(lat_grid), len(lon_grid))

    def __call__(self, alt, lat, lon):
        ialt, lagalt = self.alt_grid.iw(alt, alt=alt, lat=lat, lon=lon)
        ilat, laglat = self.lat_grid.iw(lat, alt=alt, lat=lat, lon=lon)
        ilon, laglon = self.lon_grid.iw(lon, alt=alt, lat=lat, lon=lon)

        return self.func(
            interp(self.data[ialt, ilat, ilon].data, lagalt, laglat, laglon)
        )


class computed_gridded_data:
    """
    A class representing computed gridded data.

    Parameters:
    *data: Variable-length arguments representing the data sources.
    func: The function to apply to the data sources.
    **func_kwargs: Keyword arguments to pass to the function.

    Methods:
    __call__(self, alt, lat, lon): Calls the function with the given altitude, latitude, and longitude.

    Returns:
    The result of applying the function to the data sources.
    """

    def __init__(self, *data, func, **func_kwargs):
        self.data = data
        self.func = func
        self.func_kwargs = func_kwargs

    def __call__(self, alt, lat, lon):
        return self.func(
            *[data(alt, lat, lon) for data in self.data], **self.func_kwargs
        )


class altitude_stiched_gridded_data:
    """
    Represents altitude-stitched gridded data.

    Parameters:
    - data (list): A list of functions representing data at different altitudes.
    - alts (list): A list of altitudes corresponding to the data.

    Attributes:
    - data (list): A list of functions representing data at different altitudes.
    - alts (list): A list of altitudes corresponding to the data.
    """

    def __init__(self, data, alts):
        self.data = data
        self.alts = alts

        assert len(self.data) == len(self.alts)

    def __call__(self, alt, lat, lon):
        """
        Returns the data value at the given altitude, latitude, and longitude.

        Parameters:
        - alt (float): The altitude.
        - lat (float): The latitude.
        - lon (float): The longitude.

        Returns:
        - The data value at the given altitude, latitude, and longitude.
        """
        for i in range(len(self.data)):
            if alt < self.alts[i]:
                return self.data[i](alt, lat, lon)
        self.data[-1](alt, lat, lon)
