Matpack
=======

Our general purpose multidimensional array C++ library is called Matpack. It is
based on ``std::mdspan`` (or rather the experimental version of it) and provides
an interface to create any kind of contiguous multidimensional array.

The core classes are the three data classes - ``Sparse``, ``data_t`` and ``cdata_t`` - and
the two view classes - ``view_t`` and ``strided_view_t``.  There are
also helper classes for sorted grids, gridded data, and band matrices. 

.. toctree::
   :maxdepth: 2
   
   dev.matpack.Sparse
   dev.matpack.data
   dev.matpack.view
   dev.matpack.interpolation
   dev.matpack.sorted_grid
   dev.matpack.gridded_data
   dev.matpack.band_matrix
   
.. tip::

  It is fairly trivial to map the Matpack classes to ``numpy`` arrays for the python interface.
  A lot of code for this is already in place.  If you need to add a new class to Matpack,
  please see the existing classes and add this interface to your new class as well.  It helps a lot :)
