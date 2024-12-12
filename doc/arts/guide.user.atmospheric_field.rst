.. _Sec Atmospheric Field:

The atmosphere
##############

The atmosphere in ARTS is represented by three key components:
the atmospheric field, the atmospheric field data, and the atmospheric point.

The atmospheric field contains all physical data of the entire atmosphere.
Each property of the atmospheric field (temperature field, pressure field, VMR fields, etc.)
is stored as an atmospheric field data object, representing the full 3D field of that property.
The local state of the atmospheric field is represented by an atmospheric point.
An atmospheric point effectively holds the local state of the atmosphere at a specific location.

Key notations:

- :class:`~pyarts.arts.AtmField`: The atmospheric field. An instance is in this text named: ``atm_field``.
- :class:`~pyarts.arts.AtmData`: The atmospheric field data. An instance is in this text named: ``atm_data``.
- :class:`~pyarts.arts.AtmPoint`: The local atmospheric state. An instance is in this text named: ``atm_point``.

.. note::

  This text will also only deal with the core functionality of each of these three classes.
  For more advanced and supplementary functionality, please see the API documentation of the classes themselves.
  We will explicitly not cover file handling, data conversion, or other such topics in this text.

The full atmospheric field
**************************

An atmospheric field is represented by the class :class:`~pyarts.arts.AtmField`.  It is
often presumed in ARTS that it is possible to draw a path through the atmosphere that
either hits the top of the atmosphere or the surface of the planet.

The atmospheric field holds a collection of atmospheric field data, :class:`~pyarts.arts.AtmData`,
that can be accessed and modified through a :class:`dict`-like interface.
An atmospheric field has a maximum altitude, :attr:`~pyarts.arts.AtmField.top_of_atmosphere`, above which it is considered undefined.
It has no minimum altitude, instead relying on an external ellipsoid and elevation field to define the surface of the planet below it.
The ellipsoid and the elevation field are (often) part of the :class:`~pyarts.arts.SurfaceField` class
and are described in :ref:`Sec Surface Field`.
The atmospheric field can be called as if it were a function taking altitude, latitude, and longitude 
coordinate arguments to compute or extract the local atmospheric state, :class:`~pyarts.arts.AtmPoint`.

The core operations on ``atm_field`` are:

- ``atm_field[key]``: Accessing relevant atmospheric field data: :class:`~pyarts.arts.AtmData`. See `Atmospheric field/point data access`_ for more information on what constitutes a valid ``key``
- ``atm_field.top_of_atmosphere``: The maximum altitude of the atmosphere in meters.
  This is the altitude above which the atmospheric field is undefined.
- ``atm_field(alt, lat, lon)``: Compute the local atmospheric state (:class:`~pyarts.arts.AtmPoint`) at the given coordinate.
  It is an error to call this with ``alt > atm_field.top_of_atmosphere``.

Shorthand graph for ``atm_field``:

.. graphviz::

  digraph g {
    bgcolor="#00000000";
    rankdir = "TD";
    ratio = auto;
    node [ color = "#0271BB" fontcolor = "white" style = "filled,rounded" shape = "rectangle" ];
    "Named data" [ label = "atm_field" style = "filled" ];
    "Access Operator" [ label = "atm_field[key]" shape = "ellipse" ];
    "Data attribute" [ label = "atm_field.top_of_atmosphere" shape = "ellipse" ];
    "Call operator" [ label = "atm_field(alt, lat, lon)" shape = "ellipse" ];
    "Single type of data" [ label = "AtmData" ];
    "The altitude that defines the top of the atmosphere";
    "Point-wise state of the atmosphere" [ label = "AtmPoint" ];
    "Named data" -> "Access Operator" [ arrowhead = "none" ];
    "Named data" -> "Data attribute" [ arrowhead = "none" ];
    "Named data" -> "Call operator" [ arrowhead = "none" ];
    "Access Operator" -> "Single type of data";
    "Call operator" -> "Point-wise state of the atmosphere";
    "Data attribute" -> "The altitude that defines the top of the atmosphere";
  }

A single atmospheric point
**************************

An atmospheric point holds the local state of the atmosphere.
This is required for local calculations of radiative transfer properties,
such as absorption, scattering, emission, etc.
An atmospheric point is represented by an instance of :class:`~pyarts.arts.AtmPoint`.

The main use on an atmospheric point is to access the local, numerical state of the atmosphere.

The core operations on ``atm_point`` are:

- ``atm_point[key]``: The local state as a :class:`float`. See `Atmospheric field/point data access`_ for more information on what constitutes a valid ``key``.
- ``atm_point.pressure``: The local :attr:`~pyarts.arts.AtmPoint.pressure` [Pa] as a :class:`float`.
- ``atm_point.temperature``: The local :attr:`~pyarts.arts.AtmPoint.temperature` [K] as a :class:`float`.
- ``atm_point.mag``: The local magnetic field (:attr:`~pyarts.arts.AtmPoint.mag`) as a :class:`~pyarts.arts.Vector3` [T].
- ``atm_point.wind``: The local wind field (:attr:`~pyarts.arts.AtmPoint.wind`) as a :class:`~pyarts.arts.Vector3` [m/s].

Shorthand graph for ``atm_point``:

.. graphviz::

  digraph g {
    bgcolor="#00000000";
    rankdir = "TD";
    ratio = auto;
    node [ color = "#0271BB" fontcolor = "white" style = "filled,rounded" shape = "rectangle" ];
    "Named data" [ label = "atm_point" style = "filled" ];
    "Access Operator" [ label = "atm_point[key]" shape = "ellipse" ];
    "pressure" [ label = "atm_point.pressure" shape = "ellipse" ];
    "temperature" [ label = "atm_point.temperature" shape = "ellipse" ];
    "mag" [ label = "atm_point.mag" shape = "ellipse" ];
    "wind" [ label = "atm_point.wind" shape = "ellipse" ];
    "float" [ label = "float" ];
    "Vector3" [ label = "Vector3" ];
    "Named data" -> "Access Operator" [ arrowhead = "none" ];
    "Named data" -> "temperature" [ arrowhead = "none" ];
    "Named data" -> "pressure" [ arrowhead = "none" ];
    "Named data" -> "mag" [ arrowhead = "none" ];
    "Named data" -> "wind" [ arrowhead = "none" ];
    "Access Operator" -> "float";
    "pressure" -> "float";
    "temperature" -> "float";
    "mag" -> "Vector3";
    "wind" -> "Vector3";
  }

.. note::

  The atmospheric point does not know where it is in the atmosphere.  This information is only available in the atmospheric field.
  Positional data must be retained by the user if it is needed for calculations.

Atmospheric field/point data access
***********************************

The access operator ``atm_field[key]`` is used to get and set atmospheric field data (:class:`~pyarts.arts.AtmData`)
in the atmospheric field through the use of types of keys.
Likewise, the access operator ``atm_point[key]`` is used to get and set data in the atmospheric point,
though it deals with pure floating point data.
Each type of key is meant to represent a different type of atmospheric data.
The following types of keys are available:

- :class:`~pyarts.arts.AtmKey`: Basic atmospheric data.
  Defines temperature [K], pressure [Pa], wind [m/s], and magnetic [T] components.
- :class:`~pyarts.arts.SpeciesEnum`: Content of species.
  This most often means "volume mixing ratio" (VMR) but for historical reasons there are exceptions.
  VMRs need not sum up to 1 for practical reasons.
- :class:`~pyarts.arts.SpeciesIsotope`: Isotopologue ratios.
  These are ratios of different isotopologues of the same species.
  As for VMRs they need not sum up to 1 per species for practical reasons.
  These are defaulted to values extracted from `HITRAN`_,
  complemented by other sources as necessary.
- :class:`~pyarts.arts.QuantumIdentifier`: Non-LTE data.
  These are the state distributions of energy levels of molecules required for non-LTE calculations.
- :class:`~pyarts.arts.ScatteringSpeciesProperty`: Scattering properties of the atmosphere.
  These are properties of the atmosphere that are relevant for scattering calculations.

Shorthand graph for ``key`` of different types:

.. graphviz::

  digraph g {
    bgcolor="#00000000";
    rankdir = "TD";
    ratio = auto;
    node [ color = "#0271BB" fontcolor = "white" style = "filled,rounded" shape = "rectangle" ];
    "a0" [ label = "key type" style = "filled" ];
    "b0" [ label = "AtmKey" shape = "ellipse" ];
    "b1" [ label = "SpeciesEnum" shape = "ellipse" ];
    "b2" [ label = "SpeciesIsotope" shape = "ellipse" ];
    "b3" [ label = "QuantumIdentifier" shape = "ellipse" ];
    "b4" [ label = "ScatteringSpeciesProperty" shape = "ellipse" ];
    "c0" [ label = "T, P, Mag, Wind" ];
    "c1" [ label = "VMR of O2, H2O, N2, CO2, ..." ];
    "c2" [ label = "Isotopologue ratios of O2-66, O2-67, H2O-161, H2O-162, ..." ];
    "c3" [ label = "Energy level distributions" ];
    "c4" [ label = "Scattering properties" ];
    a0 -> b0 [ arrowhead = "none" ];
    a0 -> b1 [ arrowhead = "none" ];
    a0 -> b2 [ arrowhead = "none" ];
    a0 -> b3 [ arrowhead = "none" ];
    a0 -> b4 [ arrowhead = "none" ];
    b0 -> c0;
    b1 -> c1;
    b2 -> c2;
    b3 -> c3;
    b4 -> c4;
  }

.. tip::

  Both ``atm_field["temperature"]`` and ``atm_field[pyarts.arts.AtmKey.temperature]`` will give
  the same :class:`~pyarts.arts.AtmData` back in python.  This is
  because ``pyarts.arts.AtmKey("temperature") == pyarts.arts.AtmKey.temperature``.
  The same is also true when accessing ``atm_point``, though it gives floating point values.

.. note::

  Using python :class:`str` instead of the correct type may in very rare circumstances cause name-collisions.
  Such name-collisions cannot be checked for. If it happens to you, please use the appropriate key
  type manually to correct the problem.

Atmospheric field data
**********************

The atmospheric field data is a core component of the atmospheric field.
It is stored in an instance of :class:`~pyarts.arts.AtmData`.
This type holds the entire atmospheric data for a single atmospheric property,
such as the full 3D temperature field, the full 3D pressure field, etc.
It also holds the logic for how to interpolate and extrapolate this data to any altitude, latitude, and longitude point.
As such, atmospheric field data can also be called as if it were a function taking altitude, latitude, and longitude
to return the local floating point state of the atmospheric property it holds.

These are the core operations on ``atm_data``:

- ``atm_data.data``: The core data in variant form.  See `Data types`_ for what it represents.
- ``atm_data.alt_upp``: The settings for how to extrapolate above the allowed altitude.
  What is "allowed" is defined by the data type.
- ``atm_data.alt_low``: The settings for how to extrapolate below the allowed altitude.
  What is "allowed" is defined by the data type.
- ``atm_data.lat_upp``: The settings for how to extrapolate above the allowed latitude.
  What is "allowed" is defined by the data type.
- ``atm_data.lat_low``: The settings for how to extrapolate below the allowed latitude.
  What is "allowed" is defined by the data type.
- ``atm_data.lon_upp``: The settings for how to extrapolate above the allowed longitude.
  What is "allowed" is defined by the data type.
- ``atm_data.lon_low``: The settings for how to extrapolate below the allowed longitude.
  What is "allowed" is defined by the data type.
- ``atm_data(alt, lat, lon)``: Extract the floating point value of the data at one
  specific altitude, latitude, and longitude.  Returns a single float.
  Cannot respect the top of the atmosphere because it is not available to the data.
  Instead, will strictly respect the extrapolation settings.

Shorthand graph:

.. graphviz::

  digraph g {
    bgcolor="#00000000";
    rankdir = "TD";
    ratio = auto;
    node [ color = "#0271BB" fontcolor = "white" style = "filled,rounded" shape = "rectangle" ];
    "Named data" [ label = "atm_data" style = "filled" ];
    "Data variant" [ label = "atm_data.data" shape = "ellipse" ];
    "Extrapolation settings" [ label = <atm_data.alt_upp<BR/>atm_data.alt_low<BR/>atm_data.lat_upp<BR/>atm_data.lat_low<BR/>atm_data.lon_upp<BR/>atm_data.lon_low> shape = "ellipse" ];
    "Call operator -> float" [ label = "atm_data(alt, lat, lon)" shape = "ellipse" ];
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
  
  An :class:`~pyarts.arts.AtmData` is implicitly constructible from each of the `Data types`_ described below.
  The extrapolation settings will be set to appropriate defaults when an implicit construction takes place.
  These default settings depend on the type and even available data.

.. note::

  If the extrapolation settings or the data itself cannot be used to extract a value at a point using the call-operator,
  the :class:`~pyarts.arts.AtmData` will raise an exception.  This is to ensure that the user is aware of the problem.
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

Below are the types of data that can be stored in the atmospheric data.
Each data type has its own rules for how to interpret, interpolate, and extrapolate the data.

.. tip::

  Different atmospheric field data types can be mixed in the same atmospheric field.
  There are no restrictions on how many types can be used in the same atmospheric field.

Numeric
^^^^^^^

:class:`~pyarts.arts.Numeric` data simply means that the atmosphere contains constant data.
Extrapolation rules are not relevant for this data type as it is constant everywhere.
An example of using :class:`~pyarts.arts.Numeric` as atmospheric field data is given in the following code block.

.. plot::
  :include-source:

  import matplotlib.pyplot as plt
  import numpy as np
  import pyarts

  atm_field = pyarts.arts.AtmField(toa=100e3)
  atm_field["mag_u"] = 50e-6
  atm_field["mag_v"] = 0
  atm_field["mag_w"] = 3.14

  fig = plt.figure(figsize=(14, 8))
  fig, subs = pyarts.plots.AtmField.plot(atm_field, alts=np.linspace(0, 100e3), fig=fig, keys=["mag_u", "mag_v", "mag_w"])
  subs[0].set_title("Magnetic profile u-component")
  subs[1].set_title("Magnetic profile v-component")
  subs[2].set_title("Magnetic profile w-component")
  subs[0].set_ylabel("Altitude [m]")
  [sub.set_xlabel("Field strength [T]") for sub in subs]
  plt.show()

GriddedField3
^^^^^^^^^^^^^

If the atmospheric data is of the type :class:`~pyarts.arts.GriddedField3`,
the data is defined on a grid of altitude, latitude, and longitude.
It interpolates linearly between the grid points when extracting point-wise data.
For sake of this linear interpolation, longitude is treated as a cyclic coordinate.
This data type fully respects the rules of extrapolation outside its grid.
An example of using :class:`~pyarts.arts.GriddedField3` as atmospheric field data is given in the following code block.

.. plot::
  :include-source:

  import matplotlib.pyplot as plt
  import numpy as np
  import pyarts

  atm_field = pyarts.arts.AtmField(toa=100e3)
  atm_field["t"] = pyarts.arts.GriddedField3.fromxml("planets/Earth/afgl/tropical/t.xml")
  atm_field["O2"] = pyarts.arts.GriddedField3.fromxml("planets/Earth/afgl/tropical/O2.xml")
  atm_field["H2O"] = pyarts.arts.GriddedField3.fromxml("planets/Earth/afgl/tropical/H2O.xml")

  fig = plt.figure(figsize=(14, 8))
  fig, subs = pyarts.plots.AtmField.plot(atm_field, alts=np.linspace(0, 100e3), fig=fig, keys=["t", "O2", "H2O"])
  subs[0].set_title("Temperature profile")
  subs[1].set_title("O$_2$ VMR profile")
  subs[2].set_title("H$_2$O VMR profile")
  subs[0].set_ylabel("Altitude [m]")
  subs[0].set_xlabel("Temperature [K]")
  subs[1].set_xlabel("O$_2$ VMR [-]")
  subs[2].set_xlabel("H$_2$O VMR [-]")
  subs[2].set_xscale("log")
  plt.show()

.. tip::

  It is possible to use any number of 1-long grids in a :class:`~pyarts.arts.GriddedField3` meant for use as a :class:`~pyarts.arts.AtmData`.
  The 1-long grids will by default apply the "nearest" interpolation rule for those grids, potentially reducing the atmospheric data
  to a 1D profile if only the altitude is given, or even a constant if all three grids are 1-long.

.. note::

  If the :class:`~pyarts.arts.GriddedField3` does not cover the full range of the atmosphere, the extrapolation rules will be used to
  extrapolate it.  By default, these rules are set to not allow any extrapolation.  This can be changed by setting the
  extrapolation settings as needed.  See headers `Extrapolation rules`_ and `Atmospheric field data`_ for more information.

NumericTernaryOperator
^^^^^^^^^^^^^^^^^^^^^^

This operator (:class:`~pyarts.arts.NumericTernaryOperator`) represents that the atmospheric property is purely
a function of altitude, latitude, and longitude.  The operator takes three arguments and returns a float.
Extrapolation rules are not relevant for this data type as it is a function.
An example of using :class:`~pyarts.arts.NumericTernaryOperator` as atmospheric field data is given in the following code block.

.. plot::
  :include-source:

  import matplotlib.pyplot as plt
  import numpy as np
  import pyarts

  h = pyarts.arts.GriddedField3.fromxml("planets/Earth/afgl/tropical/p.xml").grids[0]
  p = pyarts.arts.GriddedField3.fromxml("planets/Earth/afgl/tropical/p.xml").data.flatten()

  def h2p(alt, *args):
      return np.interp(alt, h, p)

  atm_field = pyarts.arts.AtmField(toa=100e3)
  atm_field["O3"] = lambda alt, lat, lon: 6e-6 if 25e3 < alt < 45e3 else 0
  atm_field["p"] = h2p

  fig = plt.figure(figsize=(14, 8))
  fig, subs = pyarts.plots.AtmField.plot(atm_field, alts=np.linspace(0, 100e3), fig=fig, keys=["O3",'p'])
  subs[0].set_title("Ozone profile using lambda-expression")
  subs[0].legend().remove()
  subs[1].set_title("Pressure profile using python function")
  subs[0].set_ylabel("Altitude [m]")
  subs[0].set_xlabel("O$_3$ VMR [-]")
  subs[1].set_xlabel("Pressure [Pa]")
  subs[1].set_xscale("log")
  plt.show()

.. tip::

  Any kind of python function-like object can be used as
  a :class:`~pyarts.arts.NumericTernaryOperator`.  It must simply take three floats and return another float.
  If you want to pass in a custom class all you need is to define ``__call__(self, alt, lat, lon)`` for it.

.. note::

  Some workspace methods populate parts of the atmospheric field with :class:`~pyarts.arts.NumericTernaryOperator` objects.
  One example is :func:`~pyarts.workspace.Workspace.atmospheric_fieldIGRF`.
  These functions are generally faster than manually created :class:`~pyarts.arts.NumericTernaryOperator` in python.
  They have 3 advantages: 1) C++ is faster than python, 2) there is no python wrapper overhead for the function call,
  and 3) we can know if these methods are safe for parallel execution, so we do not need to engage the python GIL.

.. _HITRAN: https://hitran.org/
