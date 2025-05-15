Band matrices
#############

Band matrices are special matrices that have many zero elements away from the diagonal.
These zero elements are not stored in memory, which can save a lot of space,
and they are not used in calculations, which can save a lot of compute time.  On
the other hand, the algorithms for band matrices are more complex and
slower than for dense matrices, so it is important to know that the matrix
is actually banded to use this type.

It is implemented in the class ``band_matrix`` in the ``matpack`` namespace.  The
implementation is just a wrap of ``Matrix`` with special item access.
The only operation it supports is solving :math:`\mathrm{A}\vec{x} = \vec{b}` for :math:`\vec{x}`,
where :math:`\mathrm{A}` is a band matrix and :math:`\vec{b}` is a vector.
It implements this solution by wrapping LAPACK's ``dgbsv``.
Please see its documentation for more details on the layout of the band matrix.

The band matrix is not available via the python interface but used solely for internal
solvers.

Like :doc:`dev.matpack.Sparse`, this is a very simple class that
has very specialized use cases.  It is not recommended using this class
for general purpose linear algebra, but it can be very useful for specific
problems.

Relevant files
==============

The relevant files for the band matrix are:

- ``matpack/matpack_mdspan_helpers_band_matrix.h`` - the header file for the band matrix.
- ``matpack/matpack_mdspan_helpers_band_matrix.cpp`` - the implementation file for the band matrix.
