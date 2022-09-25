////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
   \file   poly_roots.h

   Contains the code to determine roots of polynomials.

   \author Oliver Lemke
   \date 2000-03-06
*/

// ZZZ: All of poly_roots can be removed after switching to ZZZ?

#ifndef poly_roots_h
#define poly_roots_h

#include "arts.h"
#include "matpack_data.h"

int poly_root_solve(Matrix& roots, Vector& coeffs);

#endif /* poly_roots_h */
