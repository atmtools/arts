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
void amp2ext(MatrixView ext,ConstVectorView amp_coeffs,const Numeric& freq);
void amp2phamat(Tensor4View phasemat,ConstVectorView amp_coeffs);
void amp2abs(VectorView abs,ConstMatrixView ext,ConstTensor4View pha);
void ext_mat_partCalc(MatrixView ext_mat_part,MatrixView ext_mat_spt,VectorView pnd);
void double_trapez(Numeric &Integral,ConstMatrixView Integrand,Numeric &LowLimit1,Numeric &LowLimit2,Numeric &UpLimit1,Numeric &UpLimit2, Numeric &h1,Numeric &h2);
void single_trapez(Numeric &Integral,VectorView Integrand,Numeric &LowLimit1,Numeric &UpLimit1, Numeric &h);
#endif    // scatproperties_h
