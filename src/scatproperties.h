/*!
  \file   scatproperties.h
  \author Sreerekha T.R. 
  \date   Mon May 13 11:32:27 2002
*/

#ifndef scatproperties_h
#define scatproperties_h

#include "arts.h"
#include "matpackI.h"
#include "matpackIV.h"


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




#endif    // scatproperties_h
