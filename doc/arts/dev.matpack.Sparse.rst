Sparse matrices
###############

Sparse matrices are special matrices that have many zero elements.
These zero elements are not stored in memory, which can save a lot of space,
and they are not used in calculations, which can save a lot of compute time.  On
the other hand, the algorithms for sparse matrices are more complex and
slower than for dense matrices, so it is important to know that the matrix
is actually sparse to use this type.

The ``Sparse`` class is implemented in the ``matpack`` package.  The
implementation is just a wrap of ``Eigen::SparseMatrix`` for floating-point data.
The class supports basic operations like scaling, addition, subtraction,
as well as the linear algebra operations of matrix-matrix and matrix-vector multiplication.

This class is accessible through the python interface via :class:`~pyarts.arts.Sparse`.

Like :doc:`dev.matpack.band_matrix`, this is a very simple class that
has very specialized use cases.  It is not recommended using this class
for general purpose linear algebra, but it can be very useful for specific
problems.


Relevant files
==============

The relevant files for the data holding core matpack types are:

- ``matpack/matpack_sparse.h`` - the ``Sparse`` class.
