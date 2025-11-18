.. _Sec Subsurface field:

The subsurface
##############

The subsurface of ARTS is always fully 3D with coordinates of altitude, geodetic latitude, and longitude.

The subsurface in ARTS is represented by three key components:
the subsurface field, the subsurface field data, and the subsurface point.

The subsurface field contains all physical data of the entire subsurface.
Each property of the subsurface field (temperature field, density field, etc.)
is stored as a subsurface field data object, representing the full 3D field of that property.
The local state of the subsurface field is represented by a subsurface point.
A subsurface point effectively holds the local state of the subsurface at a specific location.

Key notations:

- :class:`~pyarts3.arts.SubsurfaceField`: The subsurface field. An instance is in this text named: ``subsurf_field``.  An example from the workspace is :attr:`~pyarts3.workspace.Workspace.subsurf_field`.
- :class:`~pyarts3.arts.SubsurfaceData`: The subsurface field data. An instance is in this text named: ``subsurf_data``.
- :class:`~pyarts3.arts.SubsurfacePoint`: The local subsurface state. An instance is in this text named: ``subsurf_point``.

.. note::

  This text will also only deal with the core functionality of each of these three classes.
  For more advanced and supplementary functionality, please see the API documentation of the classes themselves.
  We will explicitly not cover file handling, data conversion, or other such topics in this text.

.. warning::

  As the name of our model suggests, ARTS is primarily an atmospheric radiative transfer simulator.
  The subsurface functionality is relatively new and still very experimental.
  Use at own risk.

The full subsurface field
*************************

A subsurface field is represented by the class :class:`~pyarts3.arts.SubsurfaceField`.

The subsurface field holds a collection of subsurface field data, :class:`~pyarts3.arts.SubsurfaceData`,
that can be accessed and modified through a :class:`dict`-like interface.
A subsurface field has a minimum altitude, :attr:`~pyarts3.arts.SubsurfaceField.bottom_depth`, below which it is considered undefined.
It has no maximum altitude, instead relying on an external ellipsoid and elevation field to define the surface of the planet above it.
The ellipsoid and the elevation field are (often) part of the :class:`~pyarts3.arts.SurfaceField` class
and are described in :ref:`Sec Surface Field`.
The subsurface field can be called as if it were a function taking altitude, geodetic latitude, and longitude 
coordinate arguments to compute or extract the local subsurface state, :class:`~pyarts3.arts.SubsurfacePoint`.

The core operations on ``subsurf_field`` are:

- ``subsurf_field[key]``: Accessing relevant subsurface field data: :class:`~pyarts3.arts.SubsurfaceData`. See `Subsurface field/point data access`_ for more information on what constitutes a valid ``key``
- ``subsurf_field.bottom_depth``: The minimum altitude of the subsurface.
  This is the altitude below which the subsurface field is undefined.
- ``subsurf_field(alt, lat, lon)``: Compute the local subsurface state (:class:`~pyarts3.arts.SubsurfacePoint`) at the given coordinate.
  It is an error to call this with ``alt < subsurf_field.bottom_depth``.

Shorthand graph for ``subsurf_field``:

.. graphviz::

  digraph g {
    bgcolor="#00000000";
    rankdir = "TD";
    ratio = auto;
    node [ color = "#0271BB" fontcolor = "white" style = "filled,rounded" shape = "rectangle" ];
    "Named data" [ label = "subsurf_field" style = "filled" ];
    "Access Operator" [ label = "subsurf_field[key]" shape = "ellipse" ];
    "Data attribute" [ label = "subsurf_field.bottom_depth" shape = "ellipse" ];
    "Call operator" [ label = "subsurf_field(alt, lat, lon)" shape = "ellipse" ];
    "Single type of data" [ label = "SubsurfaceData" ];
    "The altitude that defines the bottom of the subsurface";
    "Point-wise state of the subsurface" [ label = "SubsurfacePoint" ];
    "Named data" -> "Access Operator" [ arrowhead = "none" ];
    "Named data" -> "Data attribute" [ arrowhead = "none" ];
    "Named data" -> "Call operator" [ arrowhead = "none" ];
    "Access Operator" -> "Single type of data";
    "Call operator" -> "Point-wise state of the subsurface";
    "Data attribute" -> "The altitude that defines the bottom of the subsurface";
  }

A single subsurface point
*************************

A subsurface point holds the local state of the subsurface.
This is required for local calculations of radiative transfer properties,
such as absorption, scattering, emission, etc.
A subsurface point is represented by an instance of :class:`~pyarts3.arts.SubsurfacePoint`.

The main use on a subsurface point is to access the local, numerical state of the subsurface.

The core operations on ``subsurf_point`` are:

- ``subsurf_point[key]``: The local state as a :class:`float`. See `Subsurface field/point data access`_ for more information on what constitutes a valid ``key``.
- ``subsurf_point.temperature``: The local :attr:`~pyarts3.arts.SubsurfacePoint.temperature` [K] as a :class:`float`.
- ``subsurf_point.density``: The local :attr:`~pyarts3.arts.SubsurfacePoint.density` [kg/m³] as a :class:`float`.

Shorthand graph for ``subsurf_point``:

.. graphviz::

  digraph g {
    bgcolor="#00000000";
    rankdir = "TD";
    ratio = auto;
    node [ color = "#0271BB" fontcolor = "white" style = "filled,rounded" shape = "rectangle" ];
    "Named data" [ label = "subsurf_point" style = "filled" ];
    "Access Operator" [ label = "subsurf_point[key]" shape = "ellipse" ];
    "temperature" [ label = "subsurf_point.temperature" shape = "ellipse" ];
    "density" [ label = "subsurf_point.density" shape = "ellipse" ];
    "float" [ label = "float" ];
    "Named data" -> "Access Operator" [ arrowhead = "none" ];
    "Named data" -> "temperature" [ arrowhead = "none" ];
    "Named data" -> "density" [ arrowhead = "none" ];
    "Access Operator" -> "float";
    "temperature" -> "float";
    "density" -> "float";
  }

.. note::

  The subsurface point does not know where it is in the subsurface.  This information is only available in the subsurface field.
  Positional data must be retained by the user if it is needed for calculations.

Subsurface field/point data access
**********************************

The access operator ``subsurf_field[key]`` is used to get and set subsurface field data (:class:`~pyarts3.arts.SubsurfaceData`)
in the subsurface field through the use of types of keys.
Likewise, the access operator ``subsurf_point[key]`` is used to get and set data in the subsurface point,
though it deals with pure floating point data.
Each type of key is meant to represent a different type of subsurface data.
The following types of keys are available:

- :class:`~pyarts3.arts.SubsurfaceKey`: Basic subsurface data.
  Defines temperature [K] and density [kg/m³].
- :class:`~pyarts3.arts.SubsurfacePropertyTag`: Custom data that belongs to specific models of the subsurface.

Shorthand graph for ``key`` of different types:

.. graphviz::

  digraph g {
    bgcolor="#00000000";
    rankdir = "TD";
    ratio = auto;
    node [ color = "#0271BB" fontcolor = "white" style = "filled,rounded" shape = "rectangle" ];
    "a0" [ label = "key type" style = "filled" ];
    "b0" [ label = "SubsurfaceKey" shape = "ellipse" ];
    "b1" [ label = "SubsurfacePropertyTag" shape = "ellipse" ];
    "c0" [ label = "T, rho" ];
    "c1" [ label = "Custom data" ];
    a0 -> b0 [ arrowhead = "none" ];
    a0 -> b1 [ arrowhead = "none" ];
    b0 -> c0;
    b1 -> c1;
  }

.. tip::

  Both ``subsurf_field["temperature"]`` and ``subsurf_field[pyarts3.arts.SubsurfaceKey.temperature]`` will give
  the same :class:`~pyarts3.arts.SubsurfaceData` back in python.  This is
  because ``pyarts3.arts.SubsurfaceKey("temperature") == pyarts3.arts.SubsurfaceKey.temperature``.
  The same is also true when accessing ``subsurf_point``, though it gives floating point values.

.. note::

  Using python :class:`str` instead of the correct type may in very rare circumstances cause name-collisions.
  Such name-collisions cannot be checked for. If it happens to you, please use the appropriate key
  type manually to correct the problem.

Subsurface field data
*********************

The subsurface field data is a core component of the subsurface field.
It is stored in an instance of :class:`~pyarts3.arts.SubsurfaceData`.
This type holds the entire subsurface data for a single subsurface property,
such as the full 3D temperature field, the full 3D pressure field, etc.
It also holds the logic for how to interpolate and extrapolate this data to any altitude, geodetic latitude, and longitude point.
As such, subsurface field data can also be called as if it were a function taking altitude, geodetic latitude, and longitude
to return the local floating point state of the subsurface property it holds.

These are the core operations on ``subsurf_data``:

- ``subsurf_data.data``: The core data in variant form.  See `Data types`_ for what it represents.
- ``subsurf_data.alt_upp``: The settings for how to extrapolate above the allowed altitude.
  What is "allowed" is defined by the data type.
- ``subsurf_data.alt_low``: The settings for how to extrapolate below the allowed altitude.
  What is "allowed" is defined by the data type.
- ``subsurf_data.lat_upp``: The settings for how to extrapolate above the allowed geodetic latitude.
  What is "allowed" is defined by the data type.
- ``subsurf_data.lat_low``: The settings for how to extrapolate below the allowed geodetic latitude.
  What is "allowed" is defined by the data type.
- ``subsurf_data.lon_upp``: The settings for how to extrapolate above the allowed longitude.
  What is "allowed" is defined by the data type.
- ``subsurf_data.lon_low``: The settings for how to extrapolate below the allowed longitude.
  What is "allowed" is defined by the data type.
- ``subsurf_data(alt, lat, lon)``: Extract the floating point value of the data at one
  specific altitude, geodetic latitude, and longitude.  Returns a single float.
  Cannot respect the bottom of the subsurface because it is not available to the data.
  Instead, will strictly respect the extrapolation settings.

Shorthand graph:

.. graphviz::

  digraph g {
    bgcolor="#00000000";
    rankdir = "TD";
    ratio = auto;
    node [ color = "#0271BB" fontcolor = "white" style = "filled,rounded" shape = "rectangle" ];
    "Named data" [ label = "subsurf_data" style = "filled" ];
    "Data variant" [ label = "subsurf_data.data" shape = "ellipse" ];
    "Extrapolation settings" [ label = <subsurf_data.alt_upp<BR/>subsurf_data.alt_low<BR/>subsurf_data.lat_upp<BR/>subsurf_data.lat_low<BR/>subsurf_data.lon_upp<BR/>subsurf_data.lon_low> shape = "ellipse" ];
    "Call operator -> float" [ label = "subsurf_data(alt, lat, lon)" shape = "ellipse" ];
    "The variant data" [ label = "The data type" ];
    "Type of extrapolation" [ label = "Extrapolation settings" ];
    "float" [ label = "Point-wise data; a float" ];
    "Named data" -> "Data variant" [ arrowhead = "none" ];
    "Named data" -> "Extrapolation settings" [ arrowhead = "none" ];
    "Named data" -> "Call operator -> float" [ arrowhead = "none" ];
    "Data variant" -> "The variant data";
    "Extrapolation settings" -> "Type of extrapolation";
    "Call operator -> float" -> "float";
  }

.. tip:: 
  
  An :class:`~pyarts3.arts.SubsurfaceData` is implicitly constructible from each of the `Data types`_ described below.
  The extrapolation settings will be set to appropriate defaults when an implicit construction takes place.
  These default settings depend on the type and even available data.

.. note::

  If the extrapolation settings or the data itself cannot be used to extract a value at a point using the call-operator,
  the :class:`~pyarts3.arts.SubsurfaceData` will raise an exception.  This is to ensure that the user is aware of the problem.
  Changing the extrapolation settings will likely fix the immediate problem, but be aware that the consequences of doing so
  might yield numerical differences from what was originally expected.

Extrapolation rules
-------------------

The rules for extrapolation is governed by :class:`~pyarts3.arts.InterpolationExtrapolation`.
Please see its documentation for more information.
Extrapolation happens only outside the grids of the data.
Interpreting the data inside a grid is done on a type-by-type basis.

Data types
----------

Below are the types of data that can be stored in the subsurface data.
Each data type has its own rules for how to interpret, interpolate, and extrapolate the data.

.. tip::

  Different subsurface field data types can be mixed in the same subsurface field.
  There are no restrictions on how many types can be used in the same subsurface field.

Numeric
^^^^^^^

:class:`~pyarts3.arts.Numeric` data simply means that the subsurface contains constant data.
Extrapolation rules are not relevant for this data type as it is constant everywhere.
An example of using :class:`~pyarts3.arts.Numeric` as subsurface field data is given in the following code block.

.. plot::
  :include-source:

  import matplotlib.pyplot as plt
  import numpy as np
  import pyarts3 as pyarts

  subsurf_field = pyarts.arts.SubsurfaceField(bottom_depth=-1)
  subsurf_field["t"] = 295
  subsurf_field["rho"] = 0.977

  fig = plt.figure(figsize=(14, 8))
  fig, subs = pyarts.plots.SubsurfaceField.plot(subsurf_field, alts=np.linspace(-1, 0), fig=fig, keys=["t", "rho"])
  subs.flatten()[0].set_title("Temperature profile")
  subs.flatten()[1].set_title("Density profile")
  subs.flatten()[0].set_ylabel("Depth [m]")
  subs.flatten()[0].set_xlabel("Temperature [K]")
  subs.flatten()[1].set_xlabel("Density [kg/m$^3$]")
  plt.show()

GeodeticField3
^^^^^^^^^^^^^^

If the subsurface data is of the type :class:`~pyarts3.arts.GeodeticField3`,
the data is defined on a grid of altitude, geodetic latitude, and longitude.
It interpolates linearly between the grid points when extracting point-wise data.
For sake of this linear interpolation, longitude is treated as a cyclic coordinate between [-180, 180) - please ensure your grid is defined accordingly.
This data type fully respects the rules of extrapolation outside its grid.
An example of using :class:`~pyarts3.arts.GeodeticField3` as subsurface field data is given in the following code block.

.. tip::

  It is possible to use any number of 1-long grids in a :class:`~pyarts3.arts.GeodeticField3` meant for use as a :class:`~pyarts3.arts.SubsurfaceData`.
  The 1-long grids will by default apply the "nearest" interpolation rule for those grids, potentially reducing the subsurface data
  to a 1D profile if only the altitude is given, or even a constant if all three grids are 1-long.

.. note::

  If the :class:`~pyarts3.arts.GeodeticField3` does not cover the full range of the subsurface, the extrapolation rules will be used to
  extrapolate it.  By default, these rules are set to not allow any extrapolation.  This can be changed by setting the
  extrapolation settings as needed.  See headers `Extrapolation rules`_ and `Subsurface field data`_ for more information.

.. warning::

  Even though the longitude grid is cyclic, only longitude values [-540, 540) are allowed when interpolating
  the field.  This is because we need the interpolation to be very fast and this is only possible for single
  cycles of the longitude.  Most algorithm will produce values [-360, 360] for the longitude, so this should
  in practice not be a problem for normal use-cases.  Please still ensure that the grid is defined properly
  or the interpolation routines will fail.

NumericTernaryOperator
^^^^^^^^^^^^^^^^^^^^^^

This operator (:class:`~pyarts3.arts.NumericTernaryOperator`) represents that the subsurface property is purely
a function of altitude, geodetic latitude, and longitude.  The operator takes three arguments and returns a float.
Extrapolation rules are not relevant for this data type as it is a function.
An example of using :class:`~pyarts3.arts.NumericTernaryOperator` as subsurface field data is given in the following code block.

.. plot::
  :include-source:

  import matplotlib.pyplot as plt
  import numpy as np
  import pyarts3 as pyarts

  subsurf_field = pyarts.arts.SubsurfaceField(bottom_depth=-1)
  subsurf_field["t"] = lambda alt, lat, lon: 295 + 5 * alt * 10
  subsurf_field["rho"] = lambda alt, lat, lon: 0.977 - 0.001 * alt

  fig = plt.figure(figsize=(14, 8))
  fig, subs = pyarts.plots.SubsurfaceField.plot(subsurf_field, alts=np.linspace(-1, 0), fig=fig, keys=["t", "rho"])
  subs.flatten()[0].set_title("Temperature profile")
  subs.flatten()[1].set_title("Density profile")
  subs.flatten()[0].set_ylabel("Depth [m]")
  subs.flatten()[0].set_xlabel("Temperature [K]")
  subs.flatten()[1].set_xlabel("Density [kg/m$^3$]")
  plt.show()

.. tip::

  Any kind of python function-like object can be used as
  a :class:`~pyarts3.arts.NumericTernaryOperator`.  It must simply take three floats and return another float.
  If you want to pass in a custom class all you need is to define ``__call__(self, alt, lat, lon)`` for it.
