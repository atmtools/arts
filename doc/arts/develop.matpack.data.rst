Multidimensional data
=====================

The two general purpose multidimensional data classes are ``data_t`` and ``cdata_t``.
The former is for dynamically allocated data, while the latter is for data that has
a known size at compile time.
Both classes are implemented in the ``matpack`` namespace.  They are both
template classes that can define their rank and data type.

All operations available for ``view_t`` are available for ``data_t`` and ``cdata_t``.

``data_t``
----------

The ``data_t`` class is a dynamically allocated multidimensional array.  It is
a template of the data type and of the rank of the array.  The rank of an array
is the number of access dimensions it has.  For example, a 2D array (e.g., a matrix) has a rank
2, and a 3D array has a rank 3.

The class holds a ``std::vector`` of the data type and a ``view_t`` describing its shape.
If the rank of the ``data_t`` is 1, several of the ``std::vector`` operations, e.g., ``emplace_back``, ``pop_back``, etc.,
are available.  Please see ``std::vector`` documentation for more information - and add missing operations as you need them.
For all ranks, the ``data_t`` has the following additional member methods beyond those offered by ``view_t``:

- ``resize`` - change the size of the array.
- ``flatten`` - shorthand for ``reshape({size()})``.
- ``reshape`` - change the shape of the array.  This is only available for rvalue arrays.

.. note::

  Copying and assigning to a ``data_t`` will, unlike for ``view_t``, first call ``resize`` and then copy the data.

Some types made available via the python interface are:

- :class:`~pyarts.arts.Vector` - a 1D array of :class:`~pyarts.arts.Numeric`, defined as ``data_t<Numeric, 1>``.
- :class:`~pyarts.arts.Matrix` - a 2D array of :class:`~pyarts.arts.Numeric`, defined as ``data_t<Numeric, 2>``.
- :class:`~pyarts.arts.Tensor3` - a 3D array of :class:`~pyarts.arts.Numeric`, defined as ``data_t<Numeric, 3>``.
- :class:`~pyarts.arts.Tensor4` - a 4D array of :class:`~pyarts.arts.Numeric`, defined as ``data_t<Numeric, 4>``.
- :class:`~pyarts.arts.Tensor5` - a 5D array of :class:`~pyarts.arts.Numeric`, defined as ``data_t<Numeric, 5>``.
- :class:`~pyarts.arts.Tensor6` - a 6D array of :class:`~pyarts.arts.Numeric`, defined as ``data_t<Numeric, 6>``.
- :class:`~pyarts.arts.Tensor7` - a 7D array of :class:`~pyarts.arts.Numeric`, defined as ``data_t<Numeric, 7>``.
- :class:`~pyarts.arts.StokvecVector` - a 1D array of :class:`~pyarts.arts.Stokvec`, defined as ``data_t<Stokvec, 1>``.
- :class:`~pyarts.arts.StokvecMatrix` - a 2D array of :class:`~pyarts.arts.Stokvec`, defined as ``data_t<Stokvec, 2>``.
- :class:`~pyarts.arts.StokvecTensor3` - a 3D array of :class:`~pyarts.arts.Stokvec`, defined as ``data_t<Stokvec, 3>``.
- :class:`~pyarts.arts.StokvecTensor4` - a 4D array of :class:`~pyarts.arts.Stokvec`, defined as ``data_t<Stokvec, 4>``.
- :class:`~pyarts.arts.StokvecTensor5` - a 5D array of :class:`~pyarts.arts.Stokvec`, defined as ``data_t<Stokvec, 5>``.
- :class:`~pyarts.arts.StokvecTensor6` - a 6D array of :class:`~pyarts.arts.Stokvec`, defined as ``data_t<Stokvec, 6>``.
- ... there are many more ...

.. tip::

  The ``std::vector`` guarantees that the held data is contiguous in memory as long as the data type itself is contiguous.
  This can have large performance benefits for contiguous data, such as floating-point numbers.
  Beware that this guarantee is forfeit if the data type is not itself contiguous in memory,
  so a ``data_t<data_t<Numeric, 1>, 1>``, for example, is not contiguous since ``std::vector`` hold pointers to its data.

``cdata_t``
-----------

The ``cdata_t`` class is a statically allocated multidimensional array.  It is
a template of the data type and of the shape of the array.  For example, a ``4x4``
matrix would be defined as ``cdata_t<Numeric, 4, 4>``, and a ``3x3x3`` tensor
would be defined as ``cdata_t<Numeric, 3, 3, 3>``.

The class holds a ``std::array`` of the data type.  It operates similarly to
``data_t``, but with a fixed size.  The class supports all the same operations
as ``data_t``, but without the ability to resize or reshape the array.
In addition, pure ``+`` and ``-`` operations are available for compatible arrays.
To help with use of the array, the standard tuple-interface is also available,
making it possible to write, for instance, ``auto& [a, b, c] = vec3;``,
for a ``cdata_t<Numeric, 3> vec3{};``.

You will mostly find ``cdata_t`` used for fixed-size arrays solving specific problems
in ARTS.  Often, this means it is inherited from as a base class, so that specialized
operations do not contaminate the general-purpose ``cdata_t``.

Some types made available via the python interface are:

- :class:`~pyarts.arts.Stokvec` - a 1D array of 4 :class:`~pyarts.arts.Numeric`, defined as ``cdata_t<Numeric, 4>``.  Holds spectral radiance in ARTS.
- :class:`~pyarts.arts.Propmat` - a 1D array of 7 :class:`~pyarts.arts.Numeric`, defined as ``cdata_t<Numeric, 7>``.  Holds a sparsified propagation matrix in ARTS.
- :class:`~pyarts.arts.Muelmat` - a 2D array of ``4x4`` :class:`~pyarts.arts.Numeric`, defined as ``cdata_t<Numeric, 4, 4>``.  Holds a Mueller matrix in ARTS.
- :class:`~pyarts.arts.Specmat` - a 2D array of ``4x4`` :class:`~pyarts.arts.Numeric`, defined as ``cdata_t<Complex, 4, 4>``.  Holds a spectral matrix in ARTS.
- :class:`~pyarts.arts.Vector2` - a 1D array of 2 :class:`~pyarts.arts.Numeric`, defined as ``cdata_t<Numeric, 2>``.  Used for 2D vectors in ARTS.  Not to be confused with :class:`~pyarts.arts.Vector`.  An example is the line-of-sight vector (zenith, azimuth).
- :class:`~pyarts.arts.Vector3` - a 1D array of 3 :class:`~pyarts.arts.Numeric`, defined as ``cdata_t<Numeric, 3>``. Used for 3D vectors in ARTS.  Not to be confused with :class:`~pyarts.arts.Vector`.  An example is the position vector (altitude, latitude, longitude).
- ...

.. tip::

  The ``std::array`` guarantees that the held data is contiguous in memory.  This can have
  large performance benefits.  Unlike for ``data_t``, a ``std::array`` is always contiguous,
  so a chain of ``cdata_t<cdata_t<Numeric, N>, M>`` is also contiguous and would basically be
  a ``MxN`` matrix in memory.  This is a very useful property for fixed-size arrays,
  since it means that ``data_t<cdata_t<Numeric, N>, M>`` is also a contiguous array.
  This array-type has runtime size ``M`` and compile-time size ``N``.
