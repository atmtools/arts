/*!
  \file   scatproperties.h
  \author Sreerekha T.R. 
  \date   Mon May 13 11:32:27 2002
*/
#ifndef scatproperties_h
#define scatproperties_h
#include <stdexcept>
#include <math.h>
#include "arts.h"
#include "matpackI.h"
#include "make_vector.h"
#include "array.h"
#include "logic.h"
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "mystring.h"
void amp2ext(MatrixView ext,
	     ConstVectorView amp_coeffs,
	     const Numeric& freq);

void amp2pha(Tensor4View phasemat,
	     ConstTensor3View amp_coeffs);

void amp2abs(VectorView abs,
	     ConstMatrixView ext,
	     ConstTensor4View pha,
	     ConstVectorView za_grid,
	     ConstVectorView aa_grid);


Numeric AngIntegrate_trapezoid(ConstMatrixView Integrand,
			       ConstVectorView za_grid,
			       ConstVectorView aa_grid);

#endif    // scatproperties_h
