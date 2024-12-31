Gridded data
############

A gridded data object is a multidimensional array of data, where each of the dimensions have a grid associated with them.
The grid is a one-dimensional array of values that represent the coordinates of the data in that dimension.

It is defined in the ``matpack`` namespace as a template class called ``gridded_data_t``.
The template parameters are the data type and the type of each grid.

For instance, a 3D gridded field of numeric data with unsorted grids is defined as: ``matpack::gridded_data_t<Numeric, Vector, Vector, Vector>``.
This is in fact the type :class:`pyarts.arts.GriddedField3`.  It holds the name of of the data as a string, the data itself
as a :class:`pyarts.arts.Tensor3`, the name of the grids as strings, and the grids themselves as :class:`pyarts.arts.Vector`.

The internal workings of the gridded data class is that it holds a :doc:`develop.matpack.data` object whose rank is
the same as the number of grids passed in.

The effective layout of the class is:

.. code-block:: C++

  template <typename T, typename... Grids>
  struct gridded_data_t {
    String                               name;
    data_t<T, sizeof...(Grids)>          data;
    std::array<String, sizeof...(Grids)> grid_names;
    std::tuple<Grids...>                 grids;
  };

Member methods
==============

The available member methods are a subset of those available for the ``data_t`` as well as helper methods to deal with interpolation.  The methods are:

- ``grid<I>`` - get the grid at index ``I``.  Example: ``gridded_data.grid<0>()`` gives the first grid.
- ``gridname<I>`` - get the name of the grid at index ``I``.  Example: ``gridded_data.gridname<0>()`` gives the name of the first grid.
- ``shape`` - the shape of the grids.  Example: ``gridded_data.shape() == {4, 5, 6}``.
- ``ok`` - check if the data is valid.  Example: ``gridded_data.ok()``.  Must be called in user-facing code to ensure that the data is valid.
- ``operator[]`` - the access operator.  Forwards to the ``data_t::operator[]``.  Example: ``gridded_data[3, 4, 5]``.
- ``resize`` - resize the data.  Forwards to the ``data_t::resize``.  Example: ``gridded_data.resize({4, 5, 6})``.  Note that the object is no longer ``ok`` after this call.  The grids must be manually fixed.
- ``lag<...>`` - template Lagrange operator.  The template argument is of the interpolation Lagrange type.  The input is either a grid and the allowed extrapolation, or a new grid point.
  The former returns a list of the Lagrange coefficients and the latter returns a single interpolation point.  The output can be used with a manual call to ``my_interp::reinterp`` or ``my_interp::interp``,
  potentially with ``my_interp::interpweights``.  See the interpolation section below for more information.
- ``reinterp<...>`` - template reinterpolation, wraps calling ``lag`` and then ``my_interp::reinterp``.  The input is the new grids of the data.  The output is a data object of the same rank as ``this->data``.  Example:
  ``gridded_data.reinterp<...>(new_grid1, new_grid2, new_grid3)``.  See the interpolation section below for more information.
- ``interp<...>`` - Wraps calling ``lag`` and then ``::interp``.  The input is the values on the grid to interpolate to.  Example: ``gridded_data.interp<...>(my_alt, my_lat, my_lon)``.  See the interpolation section below for more information.

Interpolation
=============

To understand how to use the interpolation methods of the gridded data class, you should first read and understand how
non-linear interpolation works in ARTS.  This is available at :doc:`develop.matpack.interpolation`.  The section of interest
is the one on higher order interpolation.

Rather than go through the exact details, here are a couple of code examples to show how to use the interpolation methods.

For 1D data
-----------

.. code-block:: C++

  const GriddedField1 f{.data_name  = "TestData-1D",
                        .data       = Vector{1, 2, 3, 4, 5, 6},
                        .grid_names = {"x"},
                        .grids      = Vector{2, 4, 6, 8, 10, 12}};
  std::print(std::cout, "Matpack data:\n{:B,Ns}\n", f);

  std::print(std::cout,
             "Matpack data:\n{:B,Ns}\n",
             f.reinterp<LagrangeInterpolation>({2, 3}, 1));
  std::print(std::cout,
             "Matpack data:\n{:B,Ns}\n",
             f.reinterp<FixedLagrangeInterpolation<1>>({2, 3, 4, 1}));
  std::print(
      std::cout, "Numeric: {}\n", f.interp<LagrangeInterpolation>(20, 1));
  std::print(
      std::cout, "Numeric: {}\n", f.interp<FixedLagrangeInterpolation<1>>(-1));

.. code-block:: bash

  Matpack data:
  {
  x: [2, 4, 6, 8, 10, 12],
  TestData-1D: [1, 2, 3, 4, 5, 6]
  }
  Matpack data:
  [1, 1.5]
  Matpack data:
  [1, 1.5, 2, 0.5]
  Numeric: 10
  Numeric: -0.5

For 2D data
-----------

.. code-block:: C++

  const GriddedField2 g{.data_name  = "TestData-2D",
                        .data       = Vector{1, 2, 3, 4, 5, 6}.reshape(2, 3),
                        .grid_names = {"x", "y"},
                        .grids      = {Vector{1, 2}, Vector{1, 2, 3}}};
  std::print(std::cout, "Matpack data:\n{:B,Ns}\n", g);

  std::print(std::cout,
             "Matpack data:\n{:B,Ns}\n",
             g.reinterp<LagrangeInterpolation, LagrangeInterpolation>(
                 {1.5}, {1.5}, 1, 20));
  std::print(
      std::cout,
      "Matpack data:\n{:B,Ns}\n",
      g.reinterp<FixedLagrangeInterpolation<1>, FixedLagrangeInterpolation<1>>(
          {2, 3, 4, 1}, {-5, 3}, 20));
  std::print(std::cout,
             "Numeric: {}\n",
             g.interp<LagrangeInterpolation, LagrangeInterpolation>(0, 0, 1));
  std::print(
      std::cout,
      "Numeric: {}\n",
      g.interp<FixedLagrangeInterpolation<1>, FixedLagrangeInterpolation<1>>(
          -1, -2));


.. code-block:: bash

  Matpack data:
  {
  x: [1, 2],
  y: [1, 2, 3],
  TestData-2D: [
  [1, 2, 3],
  [4, 5, 6]
  ]
  }
  Matpack data:
  [
  [3]
  ]
  Matpack data:
  [
  [-2, 6],
  [1, 9],
  [4, 12],
  [-5, 3]
  ]
  Numeric: -3
  Numeric: -8

Relevant files
==============

The relevant files for the data holding core matpack types are:

- ``matpack/matpack_mdspan_helpers_gridded_data_t.h`` - the ``gridded_data_t`` class.
