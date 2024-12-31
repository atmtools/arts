Multidimensional views
======================

There are two view classes in Matpack: ``view_t`` and ``strided_view_t``.
The former is for contiguous data, while the latter is for strided data as the name implies.
From a performance perspective, contiguous data is generally faster to access than strided data - often by a significant margin.

A view is either created as part of on of the :doc:`develop.matpack.data` classes or by accessing a subarray of an existing view.
Accessing a subarray in a contiguous manner will create a ``view_t``, while accessing a subarray in a strided manner will create a ``strided_view_t``.

Both classes are implemented in the ``matpack`` namespace.  They are template classes that
define their rank and a ``const``-qualified data type.
They are implemented with public inheritance from ``std::mdspan``.  ``view_t`` inherits from a ``std::mdspan`` with
``std::layout_right`` layout, while ``strided_view_t`` inherits from a ``std::mdspan`` with
``std::layout_stride`` layout.  Both are defined to have dynamic extents in all dimensions.

.. warning::

  It is very important that you follow the ``const``-qualification of the data type in all
  methods and functions that manually create views.  This is followed by the ``matpack`` classes
  internally, but if you work around this you can easily create bugs that are hard to find.

  Do not manually call the constructors of ``view_t`` or ``strided_view_t`` with data pointers or 
  custom ``std::mdspan`` without confirming that you are following the ``const``-qualification.

The ``view_t`` and ``strided_view_t`` classes are not available via the python interface.

Member methods
--------------

The list below is all the methods available for both ``view_t`` and ``strided_view_t``.
The examples will assume a ``4x5`` matrix of ``view_t<Complex, 2>``, a complex matrix, called ``mat`` and ``mat2``.

These are the common member methods for both classes:

- ``size`` - the total number of elements in the view.  Example: ``mat.size() == 20``.
- ``shape`` - the shape of the view. Example: ``mat.shape() == {4, 5}``.
- ``extent(i)`` - the extent of the view in dimension ``i``.  Example: ``mat.extent(0) == 5``, ``mat.extent(1) == 4``.
- ``stride(i)`` - the stride of the view in dimension ``i``.  Example: ``mat.stride(0) == 1``, ``mat.stride(1) == 5``.
- ``operator[]`` - access elements or sub-views of the view.  Examples:
  - ``mat[3, 4]`` accesses the complex element in row 3 and column 4.
  - ``mat[3]`` accesses the 4th row of the matrix, returning a ``view_t`` of size 5.
  - ``mat[joker, 4]`` accesses the 5th column of the matrix, returning it as a ``strided_view_t`` of size 4.
  - ``mat[3, StridedRange(0, 2, 2)]`` accesses the 1st and 3rd elements of the 4th row of the matrix, returning it as a ``strided_view_t`` of size 2.
  - ``mat[3, Range(0, 2)]`` accesses the 1st and 2nd elements of the 4th row of the matrix, returning it as a ``view_t`` of size 2.
- ``begin`` and ``end``.  These are iterators that can be used in range-based for loops.  They allow iterations over the inner view.
  Example: ``for (auto& vecview : mat) { ... }``, will have ``vecview`` be a ``view_t``.  This is equivalent to orderly access to ``mat[0]``, ``mat[1]``, ``mat[2]``, then ``mat[3]``.
- ``elem_begin`` and ``elem_end``.  There are element-wise iterators.  Combined with ``elemwise_range``, they allow iteration of the data in the view element-by-element.
  Example: ``for (auto& elem : elemwise_range(mat)) { ... }`` will iterate over all elements in the matrix.  This is equivalent to orderly access to ``mat[0, 0]``, ``mat[0, 1]``, ..., ``mat[3, 4]``.
  This type of access is much more efficient for ``view_t`` than for ``strided_view_t``.
- ``elem_at(i)``.  Helper method to access the element at index ``i``.  This is equivalent to ``*(elem_begin() + i)``, but for ``strided_view_t``, the implementation is the other way around (``elem_begin`` is implemented as a function of ``elem_at``).
  Example: ``mat.elem_at(0) == mat[0, 0]``.
- ``operator=`` - copy the data from another view, or set all elements to a single value.  Examples:
  - ``mat = 0.0`` will set all elements to zero.
  - ``mat = mat2`` will copy all elements from ``mat2`` to ``mat``.  The two views must have the same shape.
  The available right-hand side argument types are:
  - Same type as the ``*this`` value.  Will just copy the data (unless ``&RHS==&LHS``).  Must have the same shape.
  - ``CONSTANT`` of the same date type as that of the view.  Will set all elements to this value.
  - Any other ``data_t``, ``cdata_t``, ``view_t``, or ``strided_view_t`` of the same rank whose data type is convertible to the data type of the view.  Will copy the data.  Must have the same shape.
  - Any data type we can access via the helper methods ``mdshape`` and ``mdvalue``.  These are meta-functions that help access things like the Eigen matrix library types.  Will copy the data.  Must have the same shape.
- ``+=``, ``-=``, ``*=``, and ``/=``.  Will perform element-wise addition, subtraction, multiplication, and division, respectively.  The right-hand side must either have the same shape or be a constant.  Examples:  ``mat += 1.0``, ``mat -= mat2``, ``mat *= 2.0``, ``mat /= mat2``.
- ``real`` and ``imag``.  Returns a ``strided_view_t`` of only the real and imaginary values of a view to a complex data type.  Example: ``mat.real()`` will return a ``strided_view_t`` of the real values of ``mat``.
- ``ncols``, ``nrows``, ``npages``, ``nbooks``, ``nshelves``, ``nvitrines``, and ``nlibraries``.
  These are helper methods to access the shape of the view in the right-most 7 dimensions.  Ensure you do not call these on too low-rank views.  Example: ``mat.nrows() == 4``.
- ``front`` and ``back``.  These are helper methods to access the first and last element of the view.  Example: ``mat.front() == mat[0]``, ``mat.back() == mat[3]``.
- ``base_set``.  This overwrites the view with another view.
- ``base_md``.  This returns the base ``std::mdspan`` of the view.  Useful for interfacing with libraries that use ``std::mdspan``.

These methods are only available for ``strided_view_t`` but are not available for ``view_t``:

- ``unsafe_view`` - to convert a ``strided_view_t`` to a ``view_t``.  This is unsafe and should be used with caution.  All cases where this is used should be documented.  This will be removed in the future.

These methods are only available for ``view_t`` but are not available for ``strided_view_t``:

- ``view_as`` - view the data as if it had a different shape.  The size must remain constant.  Example: ``mat.view_as(2, 2, 5)`` will produce a ``view_t<Complex, 3>`` of shape ``2x2x5``.

Accessing data - ``operator[]``
-------------------------------

The access operator - ``operator[]`` - is a template that accepts any integer, ``Joker``, ``Range``, or ``StridedRange`` as its arguments.
They work as follows:

- Integers are used to access a single element or dimension in the view.
- ``Joker`` is used to access a whole dimension.  A value of the type ``Joker`` available all through ARTS is ``joker``.
- ``Range`` is used to access a range of elements in a dimension.  It has a constant stride of 1 to allow for contiguous access even in the view it produces.  It has two arguments: the start index of the new view and the number of elements.
- ``StridedRange`` is used to access a range of elements in a dimension with a custom stride.  It has three arguments: the start index of the new view, the number of elements, and the stride.

Calling the access operator on a ``view_t`` will produce either a reference to the internal element, a ``view_t``, or a ``strided_view_t``.
Calling the access operator on a ``strided_view_t`` will produce either a reference to the internal element or another ``strided_view_t``.
The type of the produced view depends on the type of the arguments to the access operator.  All integers will produce a reference to the internal element.
For a ``view_t``, there are two access patterns that produces a ``view_t``, both which requires that the right-most access arguments are of type ``Joker``.
These are:
- All the left-most access arguments are integers.
- All the left-most access arguments are first integers and then followed by a single ``Range`` argument.
Any other combination will produce a ``strided_view_t``.  And there is no safe way to go back to a ``view_t`` from a ``strided_view_t``,
even if you later access your new ``strided_view_t`` in a manner that would normally produce a ``view_t``.


