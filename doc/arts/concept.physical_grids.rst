Physical grids
##############

There are several grids that are important in ARTS.
This document will list some of them and assumptions
that are made about them implicitly throughout the codebase.

There are tests available to check these grids when it matters.
Some of these tests are enforced by the type-system itself,
others are (sometimes optional) runtime checks.  Note that
breaking these grid rules will likely lead to silent errors
in the calculations.

- Altitude grid.

  An altitude grid is a vector of altitudes, in meters, that
  defines the vertical grid.  If the grid is used to define
  atmospheric altitude, it is assumed that the elements will
  all be in ascending order.  If the grid is used to define
  the subsurface, it is assumed that the elements will be in
  descending order.  The range is regardless of sorting [-inf, inf].

- Latitude grid.

  The latitude grid is a vector of latitudes, in degrees,
  that defines a horizontal coordinate upon a ellipsoid.
  It is always in ascending order, and its range is [-90, 90].

  Use these methods to check the latitude grid for correctness:

  - :func:`~pyarts3.workspace.Workspace.atmospheric_fieldCheck`
  - :func:`~pyarts3.workspace.Workspace.surface_fieldCheck`
  - :func:`~pyarts3.workspace.Workspace.subsurface_fieldCheck`

  The only way to fix this if the grid is not valid is
  to regrid the data manually.

- Longitude grid.

  The longitude grid is a vector of longitudes, in degrees,
  that defines a horizontal coordinate upon a ellipsoid.
  It is always in ascending order, and its range is [-180, 180).

  Use these methods to check the longitudes grid for correctness:

  - :func:`~pyarts3.workspace.Workspace.atmospheric_fieldCheck`
  - :func:`~pyarts3.workspace.Workspace.surface_fieldCheck`
  - :func:`~pyarts3.workspace.Workspace.subsurface_fieldCheck`

  If these fail because of the longitude grid,
  there are methods to fix the cyclicity of the grid:

  - :func:`~pyarts3.workspace.Workspace.atmospheric_fieldFixCyclicity`
  - :func:`~pyarts3.workspace.Workspace.surface_fieldFixCyclicity`
  - :func:`~pyarts3.workspace.Workspace.subsurface_fieldFixCyclicity`

- Frequency grid.

  The frequency grid is a vector of frequencies, in hertz,
  that defines the spectral grid. It is always in ascending order,
  and range is (0, inf).

- Zenith grid.

  The zenith grid is a vector of zenith angles, in degrees,
  that defines the vertical direction in the atmosphere.
  It is always in ascending order, and its range is [0, 180].

- Azimuth grid.

  The azimuth grid is a vector of azimuth angles, in degrees,
  that defines the horizontal direction in the atmosphere.
  It is always in ascending order, and its range is FIXME:
  What are the range we allow [-180, 180) or [0, 360)?.
