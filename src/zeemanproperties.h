/*!
  \file   zeemanproperties.h
  \author Nikolay Koulev <nkoulev@uni-bremen.de>
  \date   Mon Jan 20 16:32:35 2003
  
  \brief  A header file for the functions in zeemanproperties.cc.
  
  
*/
#ifndef zeemanproperties_h
#define zeemanproperties_h

#include "arts.h"
#include "matpackI.h"
#include "matpackIV.h"


void Zeeman (Vector& f_grid,
	     Tensor3& ext_mat_zee,
	     Matrix& abs_vec_zee,
	     Matrix& xi_mat,
	     Matrix& f_z_mat,
	     Numeric& N_r,
	     Numeric& BN_r,
	     Numeric& AN_r,
	     Numeric& f_c,
	     Numeric& a1,
	     Numeric& a2,
	     Numeric& a3);


#endif 
