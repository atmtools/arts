.. _Sec Surface Field:

The surface and planet
######################

The surface in ARTS is represented by three key components: the surface field,
the surface field data, and the surface point.

The surface field contains all physical properties of the entire surface.
Each property of the surface field (temperature, elevation, planetary shape, etc.) is stored
in a surface field data object, representing the full 2D field of that property.
The local state of the surface at a specific point is represented by a surface
point. A surface point holds the local state of the surface at a specific location,
including the normal of the surface at that point.

Key notations:

- :class:`~pyarts.arts.SurfaceField` - The surface field.  An instance is in this text named: `surf_field`.
- :class:`~pyarts.arts.SurfaceData` - The surface field data.  An instance is in this text named: `surf_data`.
- :class:`~pyarts.arts.SurfacePoint` - The local surface state.  An instance is in this text named: `surf_point`.

.. note::

  This text will also only deal with the core functionality of each of these three classes.
  For more advanced and supplementary functionality, please see the API documentation of the classes themselves.
  We will explicitly not cover file handling, data conversion, or other such topics in this text.

The full surface field
**********************

A surface field is represented by the class :class:`~pyarts.arts.SurfaceField`.
It is presumed that a surface field can be intersected by a ray.

The surface field holds a collection of surface field data, :class:`~pyarts.arts.SurfaceData`,
that can be accessed and modified through a :class:`dict`-like interface.
A surface field has an ellipsoid, :attr:`~pyarts.arts.SurfaceField.ellipsoid`, which defines the basic shape of the
planet it represents.
It may contain an elevation field that modifies the shape of the surface.  It is often accompanied by an
atmospheric field that defines the state of the surface above the surface.  See :ref:`Sec surface Field`
for more details on the atmospheric field.
The surface field can be called as if it were a function taking latitude and longitude 
coordinate arguments to compute or extract the local surface state, :class:`~pyarts.arts.SurfacePoint`.

The core operations on ``surf_field`` are:

- ``surf_field[key]``: Accessing relevant surface field data: :class:`~pyarts.arts.SurfaceData`.
  See `surface field/point data access`_ for more information on what constitutes a valid ``key``
- ``surf_field.ellipsoid``: The ellipsoid shape of the surface in meters (``a``, ``b``).
  This is the default shape of the surface field.
- ``surf_field(lat, lon)``: Compute the local surface state (:class:`~pyarts.arts.SurfacePoint`) at the given coordinate.

Shorthand graph for ``surf_field``:

.. graphviz::

  digraph g {
    bgcolor="#00000000";
    rankdir = "TD";
    ratio = auto;
    node [ color = "#0271BB" fontcolor = "white" style = "filled,rounded" shape = "rectangle" ];
    "Named data" [ label = "surf_field" style = "filled" ];
    "Access Operator" [ label = "surf_field[key]" shape = "ellipse" ];
    "Data attribute" [ label = "surf_field.ellipsoid" shape = "ellipse" ];
    "Call operator" [ label = "surf_field(lat, lon)" shape = "ellipse" ];
    "Single type of data" [ label = "SurfaceData" ];
    "Point-wise state of the surface" [ label = "SurfacePoint" ];
    "Named data" -> "Access Operator" [ arrowhead = "none" ];
    "Named data" -> "Data attribute" [ arrowhead = "none" ];
    "Named data" -> "Call operator" [ arrowhead = "none" ];
    "Access Operator" -> "Single type of data";
    "Call operator" -> "Point-wise state of the surface";
    "Data attribute" -> "The ellipsoid shape that defines the basic shape of the planet";
  }

A single surface point
**********************

A surface point holds the local state of the surface.
This is required for local calculations of radiative transfer properties,
such as reflection, emission, etc.
A surface point is represented by an instance of :class:`~pyarts.arts.SurfacePoint`.

The main use on a surface point is to access the local, numerical state of the surface.

The core operations on ``surf_point`` are:

- ``surf_point[key]``: The local state as a :class:`float`. See `surface field/point data access`_ for more information on what constitutes a valid ``key``.
- ``surf_point.elevation``: The local :attr:`~pyarts.arts.SurfacePoint.elevation` [m] as a :class:`float`.
- ``surf_point.temperature``: The local :attr:`~pyarts.arts.SurfacePoint.temperature` [K] as a :class:`float`.
- ``surf_point.normal``: The (:attr:`~pyarts.arts.SurfacePoint.normal`) to the surface [degrees] as a :class:`~pyarts.arts.Vector2`.
  This gives both the zenith angle and the azimuth angle of a down-looking ray.

Shorthand graph for ``surf_point``:

.. graphviz::

  digraph g {
    bgcolor="#00000000";
    rankdir = "TD";
    ratio = auto;
    node [ color = "#0271BB" fontcolor = "white" style = "filled,rounded" shape = "rectangle" ];
    "Named data" [ label = "surf_point" style = "filled" ];
    "Access Operator" [ label = "surf_point[key]" shape = "ellipse" ];
    "elevation" [ label = "surf_point.elevation" shape = "ellipse" ];
    "temperature" [ label = "surf_point.temperature" shape = "ellipse" ];
    "normal" [ label = "surf_point.normal" shape = "ellipse" ];
    "float" [ label = "float" ];
    "Vector2" [ label = "Vector2" ];
    "Named data" -> "Access Operator" [ arrowhead = "none" ];
    "Named data" -> "temperature" [ arrowhead = "none" ];
    "Named data" -> "elevation" [ arrowhead = "none" ];
    "Named data" -> "normal" [ arrowhead = "none" ];
    "Access Operator" -> "float";
    "elevation" -> "float";
    "temperature" -> "float";
    "normal" -> "Vector2";
  }

.. note::

  The surface point does not know where it is in the surface.  This information is only available in the surface field.
  Positional data must be retained by the user if it is needed for calculations.

Surface field/point data access
*******************************

The access operator ``surf_field[key]`` is used to get and set surface field data (:class:`~pyarts.arts.SurfaceData`)
in the surface field through the use of types of keys.
Likewise, the access operator ``surf_point[key]`` is used to get and set data in the surface point,
though it deals with pure floating point data.
Each type of key is meant to represent a different type of surface data.
The following types of keys are available:

- :class:`~pyarts.arts.SurfaceKey`: Basic surface data.
  Defines temperature [K] and elevation [m] components.
- :class:`~pyarts.arts.SurfacePropertyTag`: A custom surface data type.
  This is used to define custom surface data types.

Shorthand graph for ``key`` of different types:

.. graphviz::

  digraph g {
    bgcolor="#00000000";
    rankdir = "TD";
    ratio = auto;
    node [ color = "#0271BB" fontcolor = "white" style = "filled,rounded" shape = "rectangle" ];
    "a0" [ label = "key type" style = "filled" ];
    "b0" [ label = "SurfaceKey" shape = "ellipse" ];
    "b1" [ label = "SurfacePropertyTag" shape = "ellipse" ];
    "c0" [ label = "Temperature, Elevation" ];
    "c1" [ label = "Custom Data" ];
    a0 -> b0 [ arrowhead = "none" ];
    a0 -> b1 [ arrowhead = "none" ];
    b0 -> c0;
    b1 -> c1;
  }

.. tip::

  Both ``surf_field["temperature"]`` and ``surf_field[pyarts.arts.SurfaceKey.temperature]`` will give
  the same :class:`~pyarts.arts.SurfaceData` back in python.  This is
  because ``pyarts.arts.SurfaceKey("temperature") == pyarts.arts.SurfaceKey.temperature``.
  The same is also true when accessing ``surf_point``, though it gives floating point values.

.. note::

  Using python :class:`str` instead of the correct type may in very rare circumstances cause name-collisions.
  Such name-collisions cannot be checked for. If it happens to you, please use the appropriate key
  type manually to correct the problem.

Surface field data
******************

The surface field data is a core component of the surface field.
It is stored in an instance of :class:`~pyarts.arts.SurfaceData`.
This type holds the entire surface data for a single surface property,
such as the full 2D temperature field, the full 2D elevation field, etc.
It also holds the logic for how to interpolate and extrapolate this data to any latitude and longitude point.
As such, surface field data can also be called as if it were a function taking latitude and longitude
to return the local floating point state of the surface property it holds.

These are the core operations on ``surf_data``:

- ``surf_data.data``: The core data in variant form.  See `Data types`_ for what it represents.
- ``surf_data.lat_upp``: The settings for how to extrapolate above the allowed latitude.
  What is "allowed" is defined by the data type.
- ``surf_data.lat_low``: The settings for how to extrapolate below the allowed latitude.
  What is "allowed" is defined by the data type.
- ``surf_data.lon_upp``: The settings for how to extrapolate above the allowed longitude.
  What is "allowed" is defined by the data type.
- ``surf_data.lon_low``: The settings for how to extrapolate below the allowed longitude.
  What is "allowed" is defined by the data type.
- ``surf_data(lat, lon)``: Extract the floating point value of the data at one
  specific latitude and longitude.  Returns a single float.

Shorthand graph:

.. graphviz::

  digraph g {
    bgcolor="#00000000";
    rankdir = "TD";
    ratio = auto;
    node [ color = "#0271BB" fontcolor = "white" style = "filled,rounded" shape = "rectangle" ];
    "Named data" [ label = "surf_data" style = "filled" ];
    "Data variant" [ label = "surf_data.data" shape = "ellipse" ];
    "Extrapolation settings" [ label = <surf_data.lat_upp<BR/>surf_data.lat_low<BR/>surf_data.lon_upp<BR/>surf_data.lon_low> shape = "ellipse" ];
    "Call operator -> float" [ label = "surf_data(lat, lon)" shape = "ellipse" ];
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
  
  A :class:`~pyarts.arts.SurfaceData` is implicitly constructible from each of the `Data types`_ described below.
  The extrapolation settings will be set to appropriate defaults when an implicit construction takes place.
  These default settings depend on the type and even available data.

.. note::

  If the extrapolation settings or the data itself cannot be used to extract a value at a point using the call-operator,
  the :class:`~pyarts.arts.SurfaceData` will raise an exception.  This is to ensure that the user is aware of the problem.
  Changing the extrapolation settings will likely fix the immediate problem, but be aware that the consequences of doing so
  might yield numerical differences from what was originally expected.

Extrapolation rules
-------------------

The rules for extrapolation is governed by :class:`~pyarts.arts.InterpolationExtrapolation`.
Please see its documentation for more information.
Extrapolation happens only outside the grids of the data.
Interpreting the data inside a grid is done on a type-by-type basis.

Data types
----------

Below are the types of data that can be stored in the surface data.
Each data type has its own rules for how to interpret, interpolate, and extrapolate the data.

.. tip::

  Different surface field data types can be mixed in the same surface field.
  There are no restrictions on how many types can be used in the same surface field.

Numeric
^^^^^^^

:class:`~pyarts.arts.Numeric` data simply means that the surface contains constant data.
Extrapolation rules are not relevant for this data type as it is constant everywhere.
An example of using :class:`~pyarts.arts.Numeric` as surface field data is given in the following code block.

.. code-block:: python

  import pyarts

  surf_field = pyarts.arts.SurfaceField("Earth")
  surf_field["h"] = 0.
  surf_field["t"] = 295.

  print(surf_field(0, 0))

GriddedField2
^^^^^^^^^^^^^

If the surface data is of the type :class:`~pyarts.arts.GriddedField2`,
the data is defined on a grid of latitude and longitude.
It interpolates linearly between the grid points when extracting point-wise data.
For sake of this linear interpolation, longitude is treated as a cyclic coordinate.
This data type fully respects the rules of extrapolation outside its grid.

.. note::

  If the :class:`~pyarts.arts.GriddedField2` does not cover the full range of the surface, the extrapolation rules will be used to
  extrapolate it.  By default, these rules are set to not allow any extrapolation.  This can be changed by setting the
  extrapolation settings as needed.  See headers `Extrapolation rules`_ and `surface field data`_ for more information.

NumericBinaryOperator
^^^^^^^^^^^^^^^^^^^^^

This operator (:class:`~pyarts.arts.NumericBinaryOperator`) represents that the surface property is purely
a function of latitude and longitude.  The operator takes two arguments and returns a float.
Extrapolation rules are not relevant for this data type as it is a function.

.. tip::

  Any kind of python function-like object can be used as
  a :class:`~pyarts.arts.NumericBinaryOperator`.  It must simply take two floats and return another float.
  If you want to pass in a custom class all you need is to define ``__call__(self, lat, lon)`` for it.

.. _HITRAN: https://hitran.org/
