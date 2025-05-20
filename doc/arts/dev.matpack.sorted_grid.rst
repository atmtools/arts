Sorted grids
############

Sorted grids are a wrapper around the :class:`~pyarts.arts.Vector` class
that demands that all elements fulfill a certain condition.  This condition
is a template the class take as input, which is used to ensure that the elements are sorted.
It is implemented as the class ``grid_t`` in the C++ code.

There are only two classes named of the type:

- :class:`~pyarts.arts.AscendingGrid` which is a grid where every subsequent element must be strictly not less or equal to the previous element.
- :class:`~pyarts.arts.DescendingGrid` which is a grid where every subsequent element must be strictly not greater or equal to the previous element.

The usefulness of these classes is simply that you know that their values are sorted.

The member methods that they provide are a limited subset of the :class:`~pyarts.arts.Vector` class that do not modify 
state of the :class:`~pyarts.arts.Vector` in any way.  Constant iteration (``begin`` and ``end``), size-queries (``size`` and ``shape``),
as well as constant element and ranged access (``operator[]``).
To create a sorted grid, simply pass a :class:`~pyarts.arts.Vector` to the constructor or the ``operator=``.

For convenience, ``grid_t`` provides a static helper method that can be used to check if your :class:`~pyarts.arts.Vector` or view thereof is sorted.  This method is called:
``grid_t<>::is_sorted(vec)`` and simply returns a boolean value for the state of the vector as a grid.  ``true`` being sorted and ``false`` being not sorted.

Extending a grid
================

Because we want to ensure that all grids are always sorted properly, the
safest way to extend a grid is to simply copy the values it holds into a
:class:`~pyarts.arts.Vector` and then modify your new :class:`~pyarts.arts.Vector` as needed.
The last step is to move then to move your now modified :class:`~pyarts.arts.Vector` back into the grid.

For cases where you know this will not work well, you can use the ``extend_grid_t`` class.
It takes a reference to the grid and exposes the :class:`~pyarts.arts.Vector` interface
directly.  Upon destroying the ``extend_grid_t`` object, the grid is checked to ensure the sorting is correct.

.. warning::

  Do not have a living instance of a ``extend_grid_t``-object in the same context you intend to consume ``grid_t``.
  The sorting is never guaranteed to be correct until the ``extend_grid_t`` object is destroyed.

Relevant files
==============

The relevant files for the data holding core matpack types are:

- ``matpack/matpack_mdspan_helpers_grid_t.h`` - the ``grid_t`` and ``extend_grid_t`` classes.
