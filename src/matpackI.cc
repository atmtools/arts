/**
   \file   matpackI.cc

   Most of the implementation of matpack is in the header file, in
   order to allow inlining. But those methods that also require Array
   are here, in order to avoid that matpackI.h depends on array.h.

   \author Stefan Buehler
   \date   2001-09-15
*/

#include "matpackI.h"
#include "array.h"

/** Assignment operator from Array<Numeric>. This copies the data from
    an Array<Numeric> to this VectorView. Dimensions must agree! 
    Resizing would destroy the selection that we might have done in
    this VectorView by setting its range. 

    Array<Numeric> can be useful to collect things in, because there
    is a .push_back method for it. Then, after collecting we usually
    have to transfer the content to a Vector. With this assignment
    operator that's easy. */
VectorView VectorView::operator=(const Array<Numeric>& v)
{
  //  cout << "Assigning VectorView from Array<Numeric>.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.nelem());

  // Iterators for Array:
  Array<Numeric>::const_iterator i=v.begin();
  const Array<Numeric>::const_iterator e=v.end();
  // Iterator for Vector:
  Iterator1D target = begin();

  for ( ; i!=e ; ++i,++target )
    *target = *i;

  return *this;
}

