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
	     Tensor3& ExtZ,
	     Matrix& absvZ);
 

#endif 
