#pragma once


////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
  \file   legendre.h

  Contains the code to calculate Legendre polynomials.

  \author Oliver Lemke
  \date 2003-08-14
  */

#include "arts.h"
#include "grids.h"
#include "matpack_data.h"

namespace GSL::Integration {
bool GaussLegendre(Vector &x, Vector &w, Index n);
}  // namespace GSL::Integration
