Interpolation in ARTS
#####################

Linear interpolation
====================

There are no general single-step interpolation functions in ARTS.
Instead, there is a set of useful utility functions that can be used
to achieve interpolation. Roughly, you can separate these into
functions determining grid position arrays, functions determining
interpolation weight tensors, and functions applying the
interpolation. Doing an interpolation thus requires a chain of
function calls:

1. ``gridpos`` (one for each interpolation dimension)
2. ``interpweights``
3. ``interp``

Currently implemented in ARTS is multilinear interpolation in up to 6
dimensions. (Is the 6D case called hexa-linear interpolation?)  The
necessary functions and their interaction will be explained in this
chapter.

Implementation files
--------------------

Variables and functions related to interpolation are defined in the files:

- ``interpolation.h``
- ``interpolation.cc``
- ``test interpolation.cc``

The first two files contain the declarations and implementation,
the last file some usage examples.

Green and blue interpolation
----------------------------

.. figure:: Figs/interpolation/interpolation_types.svg
  :width: 800

  The two different types of interpolation. Green (dotted):
  Interpolation to a new grid, output has same dimension as input,
  in this case 2D. Blue (dashed): Interpolation to a sequence of
  points, output is always 1D.

There are two different types of interpolation in ARTS:

- *Green Interpolation*: Interpolation of a gridded field to a new grid.
- *Blue Interpolation*: Interpolation of a gridded field to asequence of positions.

The figure above illustrates the different types
for a 2D example.

The first step of an interpolation always consists in determining
where your new points are, relative to the original grid. You can do
this separately for each dimension. The positions have to be stored
somehow, which is described in the next section.

Grid checking functions
-----------------------

Before you do an interpolation, you should check that the new grid is
inside the old grid. (Or only slightly outside.) You can use the
convenience function ``chk_interpolation_grids`` for this
purpose, which resides in file ``check_input.cc``. The
function has the following parameters:

.. code-block:: C++

  const String&     which_interpolation   A string describing the 
                                          interpolation for which 
                                          the grids are intended. 
  ConstVectorView   old_grid              The original grid.
  ConstVectorView   new_grid              The new grid.
  const Numeric&    extpolfac             The extrapolation fraction. 
                                          See gridpos function for 
                                          details. Has a default 
                                          value, which is consistent 
                                          with gridpos.  

There is also a special version for the case that the new grid is just
a scalar. What the function does is check if old and new grid for an
interpolation are ok. If not, it throws a detailed runtime error
message.

The parameter ``extpolfac`` determines how much extrapolation
is tolerated. Its default value is 0.5, which means that we allow
extrapolation as far out as half the spacing of the last two grid
points on that edge of the grid.

The ``chk_interpolation_grids`` function is quite thorough.
It checks not only the grid range, but also the proper sorting,
whether there are duplicate values, etc.. It is not completely cheap
computationally. Its intended use is at the beginning of workspace
methods, when you check the input variables and issue runtime errors
if there are any problems. The runtime error thrown also explains in
quite a lot of detail what is actually wrong with the grids.


Grid positions
--------------

A grid position specifies where an interpolation point is, relative
to the original grid. It consists of three parts, an :class:`pyarts3.arts.Index` giving the
original grid index below the interpolation point, a :class:`pyarts3.arts.Numeric`
giving the fractional distance to the next original grid point, and a
:class:`pyarts3.arts.Numeric` giving 1 minus this number. Of course, the last element is
redundant. However, it is efficient to store this, since it is used
many times over. We store the two numerics in a plain C array of
dimension 2. (No need to use a fancy Array or Vector for this, since
the dimension is fixed.) So the structure ``GridPos`` looks like:

.. code-block:: C++

  struct GridPos  {
    Index   idx;      /*!< Original grid index below
                            interpolation point. */
    Numeric fd[2];    /*!< Fractional distance to next point
                            (0<=fd[0]<=1), fd[1] = 1-fd[0]. */
  };

For example, ``idx=3`` and ``fd=0.5`` means that this interpolation point is
half-way between index 3 and 4 of the original grid.  Note, that
"below" in the first paragraph means "with a lower index". If the
original grid is sorted in descending order, the value at the grid
point below the interpolation point will be numerically higher than
the interpolation point.  In other words, grid positions and
fractional distances are defined relative to the order of the original
grid. Examples:

.. code-block:: C++

  old grid = 2 3
  new grid = 2.25
  idx      = 0
  fd[0]    = 0.25

  old grid = 3 2
  new grid = 2.25
  idx      = 0
  fd[0]    = 0.75

Note that ``fd[0]`` is different in the second case, because the old grid
is sorted in descending order. Note also that ``idx`` is the same in
both cases.

Grid positions for a whole new grid are stored in an ``Array<GridPos>``
(called ``ArrayOfGridPos``). 

Setting up grid position arrays
-------------------------------

There is only one function to set up grid position arrays, namely 
``gridpos``:

.. code-block:: C++

  void gridpos( ArrayOfGridPos& gp,
                ConstVectorView old_grid,
                ConstVectorView new_grid 
                const Numeric&  extpolfac=0.5 );

Some points to remember:

- As usual, the output ``gp`` has to have the right dimension. 
  
- The old grid has to be strictly sorted. It can be in ascending
  or descending order. But there must not be any duplicate values.
  Furthermore, the old grid must contain at least two points.
  
- The new grid does not have to be sorted, but the function will be
  faster if it is sorted or mostly sorted. It is ok if the new grid
  contains only one point.
  
- The beauty is, that this is all it needs to do also interpolation in
  higher dimensions: You just have to call ``gridpos`` for all the
  dimensions that you want to interpolate.
  
- Note also, that for this step you do not need the field itself at
  all!

- If you want to use the returned gp object for something else
  than interpolation, you should know that gridpos guarantees the
  following:

  - For the ascending old grid case: 

    .. code-block:: C++

      old_grid[tgp.idx]<=tng || tgp.idx==0

  - And for the descending old grid case: 

    .. code-block:: C++

      old_grid[tgp.idx]>=tng || tgp.idx==0

- Finally, note that parameter ``extpolfac`` plays the
  same role as explained above.

Interpolation weights
---------------------

As explained in the "Numerical Recipes"
:cite:p:`numerical_recipes_C:97`, 2D bi-linear interpolation means, that
the interpolated value is a weighted average of the original field at
the four corner points of the grid square in which the interpolation
point is located. Taking the corner points in the order indicated in the Figure
below, the interpolated value is given by:

.. math::

  y(t,u)
  &=& (1-t)*(1-u)*y_1 \nonumber \\
  & & \mbox{} + t*(1-u)*y_2 \nonumber \\
  & & \mbox{} + (1-t)*u*y_3 \nonumber \\
  & & \mbox{} + t*u*y_4 \nonumber \\
  &=& w_1*y_1 + w_2*y_2 + w_3*y_3 + w_4*y_4

where :math:`t` and :math:`u` are the fractional distances between the
corner points in the two dimensions, :math:`y_i` are the field values
at the corner points, and :math:`w_i` are the interpolation weights.

.. figure:: Figs/interpolation/interpolation_square.svg
  :width: 400

  The grid square for 2D interpolation. The numbers 1... 4
  mark the corner points, IP is the interpolation point, :math:`t` and :math:`u`
  are the fractional distances in the two dimensions.
  (By the way, I have discovered that this is exactly the result that
  you get if you first interpolate linearly in one dimension, then in
  the other. I was playing around with this a bit, but it is the more
  efficient way to pre-calculate the :math:`w_i` and do all dimensions at once.)

How many interpolation weights one needs for a multilinear
interpolation depends on the dimension of the interpolation: There are
exactly :math:`2^n` interpolation weights for an :math:`n` dimensional
interpolation.  These weights have to be computed for each
interpolation point (each grid point of the new grid, if we do a
"green" type interpolation. Or each point in the sequence, if we do a
"blue" type interpolation).

This means, calculating the interpolation weights is not exactly
cheap, especially if one interpolates simultaneously in many
dimensions. On the other hand, one can save a lot by re-using the
weights.  Therefore, interpolation weights in ARTS are stored in a
tensor which has one more dimension than the output field. The last
dimension is for the weight, so this last dimension has the extent 4
in the 2D case, 8 in the 3D case, and so on (always :math:`2^n`).

In the case of a "blue" type interpolation, the weights are
always stored in a matrix, since the output field is always 1D (a
vector). 

Setting up interpolation weight tensors
---------------------------------------

Interpolation weight tensors can be computed by a family of functions,
which are all called ``interpweights``. Which function is actually
used depends on the dimension of the input and output quantities. For
this step we still do not need the actual fields, just the grid
positions.

Blue interpolation
~~~~~~~~~~~~~~~~~~

In this case the functions are:

.. code-block:: C++

  void interpweights( MatrixView itw,
                      const ArrayOfGridPos& cgp );
  void interpweights( MatrixView itw,
                      const ArrayOfGridPos& rgp,
                      const ArrayOfGridPos& cgp );
  void interpweights( MatrixView itw,
                      const ArrayOfGridPos& pgp,
                      const ArrayOfGridPos& rgp,
                      const ArrayOfGridPos& cgp );
  void interpweights( MatrixView itw,
                      const ArrayOfGridPos& vgp,
                      const ArrayOfGridPos& sgp,
                      const ArrayOfGridPos& bgp,
                      const ArrayOfGridPos& pgp,
                      const ArrayOfGridPos& rgp,
                      const ArrayOfGridPos& cgp );

In all cases, the dimension of ``itw`` must be consistent with the
given grid position arrays and the dimension of the interpolation
(last dimension :math:`2^n`). Because the grid position arrays are
interpreted as defining a sequence of positions they must all have
the same length.

Green interpolation
~~~~~~~~~~~~~~~~~~~

In this case the functions are:

.. code-block:: C++

  void interpweights( Tensor3View itw,
                      const ArrayOfGridPos& rgp,
                      const ArrayOfGridPos& cgp );
  void interpweights( Tensor4View itw,
                      const ArrayOfGridPos& pgp,
                      const ArrayOfGridPos& rgp,
                      const ArrayOfGridPos& cgp );
  void interpweights( Tensor5View itw,
                      const ArrayOfGridPos& bgp,
                      const ArrayOfGridPos& pgp,
                      const ArrayOfGridPos& rgp,
                      const ArrayOfGridPos& cgp );
  void interpweights( Tensor6View itw,
                      const ArrayOfGridPos& sgp,
                      const ArrayOfGridPos& bgp,
                      const ArrayOfGridPos& pgp,
                      const ArrayOfGridPos& rgp,
                      const ArrayOfGridPos& cgp );
  void interpweights( Tensor7View itw,
                      const ArrayOfGridPos& vgp,
                      const ArrayOfGridPos& sgp,
                      const ArrayOfGridPos& bgp,
                      const ArrayOfGridPos& pgp,
                      const ArrayOfGridPos& rgp,
                      const ArrayOfGridPos& cgp );

In this case the grid position arrays are interpreted as defining the
grids for the interpolated field, therefore they can have different
lengths. Of course, ``itw`` must be consistent with the length of
all the grid position arrays, and with the dimension of the
interpolation (last dimension :math:`2^n`).

The actual interpolation
------------------------

For this final step we need the grid positions, the
interpolation weights, and the actual fields. For each interpolated
value, the weights are applied to the appropriate original field values
and the sum is taken (see Equation above). The ``interp`` family of functions
performs this step.

Blue interpolation
~~~~~~~~~~~~~~~~~~

.. code-block:: C++

  void interp( VectorView            ia,
              ConstMatrixView       itw,
              ConstVectorView       a,    
              const ArrayOfGridPos& cgp);
  void interp( VectorView            ia,
              ConstMatrixView       itw,
              ConstMatrixView       a,    
              const ArrayOfGridPos& rgp,
              const ArrayOfGridPos& cgp);
  void interp( VectorView            ia,
              ConstMatrixView       itw,
              ConstTensor3View      a,    
              const ArrayOfGridPos& pgp,
              const ArrayOfGridPos& rgp,
              const ArrayOfGridPos& cgp);
  void interp( VectorView            ia,
              ConstMatrixView       itw,
              ConstTensor4View      a,    
              const ArrayOfGridPos& bgp,
              const ArrayOfGridPos& pgp,
              const ArrayOfGridPos& rgp,
              const ArrayOfGridPos& cgp);
  void interp( VectorView            ia,
              ConstMatrixView       itw,
              ConstTensor5View      a,    
              const ArrayOfGridPos& sgp,
              const ArrayOfGridPos& bgp,
              const ArrayOfGridPos& pgp,
              const ArrayOfGridPos& rgp,
              const ArrayOfGridPos& cgp);
  void interp( VectorView            ia,
              ConstMatrixView       itw,
              ConstTensor6View      a,    
              const ArrayOfGridPos& vgp,
              const ArrayOfGridPos& sgp,
              const ArrayOfGridPos& bgp,
              const ArrayOfGridPos& pgp,
              const ArrayOfGridPos& rgp,
              const ArrayOfGridPos& cgp);

Green interpolation
~~~~~~~~~~~~~~~~~~~

.. code-block:: C++

  void interp( MatrixView            ia,
              ConstTensor3View      itw,
              ConstMatrixView       a,   
              const ArrayOfGridPos& rgp,
              const ArrayOfGridPos& cgp);
  void interp( Tensor3View           ia,
              ConstTensor4View      itw,
              ConstTensor3View      a,   
              const ArrayOfGridPos& pgp,
              const ArrayOfGridPos& rgp,
              const ArrayOfGridPos& cgp);
  void interp( Tensor4View           ia,
              ConstTensor5View      itw,
              ConstTensor4View      a,   
              const ArrayOfGridPos& bgp,
              const ArrayOfGridPos& pgp,
              const ArrayOfGridPos& rgp,
              const ArrayOfGridPos& cgp);
  void interp( Tensor5View           ia,
              ConstTensor6View      itw,
              ConstTensor5View      a,   
              const ArrayOfGridPos& sgp,
              const ArrayOfGridPos& bgp,
              const ArrayOfGridPos& pgp,
              const ArrayOfGridPos& rgp,
              const ArrayOfGridPos& cgp);
  void interp( Tensor6View           ia,
              ConstTensor7View      itw,
              ConstTensor6View      a,   
              const ArrayOfGridPos& vgp,
              const ArrayOfGridPos& sgp,
              const ArrayOfGridPos& bgp,
              const ArrayOfGridPos& pgp,
              const ArrayOfGridPos& rgp,
              const ArrayOfGridPos& cgp);

Examples
--------

A simple example
~~~~~~~~~~~~~~~~

This example is contained in file ``test_interpolation.cc``.

.. code-block:: C++

  void test05()
  {
    cout << "Very simple interpolation case\n";

    Vector og(1,5,+1);            // 1, 2, 3, 4, 5
    Vector ng(2,5,0.25);          // 2.0, 2,25, 2.5, 2.75, 3.0

    cout << "Original grid:\n" << og << "\n";
    cout << "New grid:\n" << ng << "\n";

    // To store the grid positions:
    ArrayOfGridPos gp(ng.nelem());

    gridpos(gp,og,ng);
    cout << "Grid positions:\n" << gp;

    // To store interpolation weights:
    Matrix itw(gp.nelem(),2);
    interpweights(itw,gp);
      
    cout << "Interpolation weights:\n" << itw << "\n";

    // Original field:
    Vector of(og.nelem(),0);
    of[2] = 10;                   // 0, 0, 10, 0, 0

    cout << "Original field:\n" << of << "\n";

    // Interpolated field:
    Vector nf(ng.nelem());

    interp(nf, itw, of, gp);

    cout << "New field:\n" << nf << "\n";
  }

Ok, maybe you think this is not so simple, but a
large part of the code is either setting up the example grids and
fields, or output. And here is how the output looks like:

.. code-block:: C++

  Very simple interpolation case
  Original grid:
    1   2   3   4   5
  New grid:
    2 2.25 2.5 2.75   3
  Grid positions:
    1 0    1
    1 0.25 0.75
    1 0.5  0.5
    1 0.75 0.25
    1 1    0
  Interpolation weights:
    1   0
  0.75 0.25
  0.5 0.5
  0.25 0.75
    0   1
  Original field:
    0   0  10   0   0
  New field:
    0 2.5   5 7.5  10

A more elaborate example
~~~~~~~~~~~~~~~~~~~~~~~~~

What if you want to interpolate only some dimensions of a tensor,
while retaining others? --- You have to make a loop yourself, but it
is very easy. Below is an explicit example for a more complicated
interpolation case. (Green type interpolation of all pages of a
Tensor3.) This example is also contained in file
``test_interpolation.cc``.

.. code-block:: C++

  void test04()
  {
    cout << "Green type interpolation of all "
        << "pages of a Tensor3\n";

    // The original Tensor is called a, the new one n. 

    // 10 pages, 20 rows, 30 columns, all grids are: 1,2,3
    Vector  a_pgrid(1,3,1), a_rgrid(1,3,1), a_cgrid(1,3,1); 
    Tensor3 a( a_pgrid.nelem(),
              a_rgrid.nelem(),
              a_cgrid.nelem() ); 
    a = 0;
    // Put some simple numbers in the middle of each page:
    a(0,1,1) = 10;
    a(1,1,1) = 20;
    a(2,1,1) = 30;

    // New row and column grids:
    // 1, 1.5, 2, 2.5, 3
    Vector  n_rgrid(1,5,.5), n_cgrid(1,5,.5); 
    Tensor3 n( a_pgrid.nelem(),
              n_rgrid.nelem(),
              n_cgrid.nelem() ); 

    // So, n has the same number of pages as a, 
    // but more rows and columns.

    // Get the grid position arrays:
    ArrayOfGridPos n_rgp(n_rgrid.nelem()); // For rows.
    ArrayOfGridPos n_cgp(n_cgrid.nelem()); // For columns.

    gridpos( n_rgp, a_rgrid, n_rgrid );
    gridpos( n_cgp, a_cgrid, n_cgrid );

    // Get the interpolation weights:
    Tensor3 itw( n_rgrid.nelem(), n_cgrid.nelem(), 4 );
    interpweights( itw, n_rgp, n_cgp );

    // Do a "green" interpolation for all pages of a:

    for ( Index i=0; i<a.npages(); ++i )
      {
        // Select the current page of both a and n:
        ConstMatrixView ap = a( i,
                                Range(joker), Range(joker) );
        MatrixView      np = n( i,
                                Range(joker), Range(joker) );

        // Do the interpolation:
        interp( np, itw, ap, n_rgp, n_cgp );

        // Note that this is efficient, because interpolation
        // weights and grid positions are re-used.
      }

    cout << "Original field:\n";
    for ( Index i=0; i<a.npages(); ++i )
        cout << "page " << i << ":\n"
            << a(i,Range(joker),Range(joker)) << "\n";

    cout << "Interpolated field:\n";
    for ( Index i=0; i<n.npages(); ++i )
        cout << "page " << i << ":\n"
            << n(i,Range(joker),Range(joker)) << "\n";
  }

The output is:

.. code-block:: C++

  Green type interpolation of all pages of a Tensor3
  Original field:
  page 0:
    0   0   0
    0  10   0
    0   0   0
  page 1:
    0   0   0
    0  20   0
    0   0   0
  page 2:
    0   0   0
    0  30   0
    0   0   0
  Interpolated field:
  page 0:
    0   0   0   0   0
    0 2.5   5 2.5   0
    0   5  10   5   0
    0 2.5   5 2.5   0
    0   0   0   0   0
  page 1:
    0   0   0   0   0
    0   5  10   5   0
    0  10  20  10   0
    0   5  10   5   0
    0   0   0   0   0
  page 2:
    0   0   0   0   0
    0 7.5  15 7.5   0
    0  15  30  15   0
    0 7.5  15 7.5   0
    0   0   0   0   0


Higher order interpolation
==========================

Everything that was written so far in this chapter referred to linear
interpolation, which uses two neighboring data points in the 1D
case. But ARTS also has a framework for higher order polynomial
interpolation. It is defined in the the file

- ``matpack/lagrange_interp.h``

The higher order interpolation framework uses a template class,
``lagrange_interp::lag_t<>``, which is described below.
The class takes two compile-time parameters, the polynomial order
:math:`O` and  grid transformation type.  The rules for how
to transform the grid are in the file ``matpack/lagrange_interp.h``.
Regardless, the template type holds weights and indices to the
original grid, which are used to compute the interpolated value
and which are used to map into the data field.

Three styles of interpolations are implemented as functions:

- ``interp``: Pure interpolation.  Reduces the data field to the
  a single interpolated value.  Takes as many ``lagrange_interp::lag_t<>``
  as the rank of the data field.  May also use an ``interpweights``
  approach similar to the linear case, but this is not required.
  The signature is either ``T interp(data, lag_t<>, lag_t<>, ...)``
  or ``T interp(data, interpweights, lag_t<>, lag_t<>, ...)``, to
  return the interpolated value ``T``.
- ``reinterp``: Reinterpolates the data to a new grid.  Retains
  the same rank as the data field but changes the grid.  This
  takes as many lists of ``lagrange_interp::lag_t<>`` as the
  rank of the data field.  The lists can be of any type that
  follows rules similar to ``std::vector`` or ``std::array``.
  May use ``reinterpweights`` as a multidimensional
  interpolation weights approach similar to the linear case,
  and it may reuse data.  This yields four possible signatures:
  ``void reinterp(data_t&, data, list<lag_t<>>, list<lag_t<>>, ...)``,
  ``void reinterp(data_t&, data, reinterpweights, list<lag_t<>>, list<lag_t<>>, ...)``,
  ``data_t reinterp(data, list<lag_t<>>, list<lag_t<>>, ...)``, and
  ``data_t reinterp(data, reinterpweights, list<lag_t<>>, list<lag_t<>>, ...)``.
- ``flat_interp``: Interpolates a line through the data field composed
  of all the coordinates of the ``lag_t<>`` objects.  This reduces
  the data_field to a vector of values.  For rank 1 data fields,
  this is equivalent to the ``interp`` function.  The methods take
  as many lists of ``lag_t<>`` as the rank of the data field.  These
  must all have the same length, and this length is the size of the
  output vector.  May use ``flat_interpweights`` as a
  multidimensional interpolation weights approach similar to the
  linear case, and it may reuse data.  This yields four possible
  signatures:
  ``void flat_interp(data_t&, data, list<lag_t<>>, list<lag_t<>>, ...)``,
  ``void flat_interp(data_t&, data, flat_interpweights, list<lag_t<>>, list<lag_t<>>, ...)``,
  ``data_t flat_interp(data, list<lag_t<>>, list<lag_t<>>, ...)``, and
  ``data_t flat_interp(data, flat_interpweights, list<lag_t<>>, list<lag_t<>>, ...)``.

Weights
-------

Lagrange weights are computed as:

.. math::

  l_j(x) = \prod_{\substack{0 \leq m \leq O \\ m \neq j}}
  \frac{u(f(x)) - u(f(x _m))}{u(f(x_j)) - u(f(x_m))}

where :math:`f` is a grid scaling function and :math:`u` deals with cyclicity.
The transformation can be whatever, but most common
is to use the identity function :math:`f(x) = x`, which is thus the default.

Grid cyclicity
--------------

If the grid is cyclic :math:`\left[X_0, X_1\right)`,
and for simplicity of writing these examples :math:`f(x) = x`, then
the algorithmic first thing that happens is that :math:`x` is shifted
to be within this range :math:`\left[X_0, X_1\right)`.
If the grid coordinates closest to :math:`x` after the shift are
:math:`\left[x_j, x_{j+1}, \cdots, x_{j+n}\right]`,
and :math:`x_j \leq x \lt x_{j+n}` then :math:`u(x) = x` for all values.
If :math:`x \lt x_j`, then :math:`u(x) = x + X_0 - X_1`
for all :math:`x_i \geq \left(X_1-X_0\right) / 2`.
If :math:`x \gt x_{j+n}`, then :math:`u(x) = x - X_0 + X_1`
for all :math:`x_i \lt \left(X_1-X_0\right) / 2`.
