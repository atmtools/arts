/*!
  \file   oem.h
  \author simon <simonpf@chalmers.se>
  \date   Fri Apr 17 16:17:54 2015

  \brief Optimal estimation method for retrieval.
*/

#ifndef oem_h
#define oem_h

#include "matpackI.h"

// Optimal estimation method for linear models.
void linear_oem( VectorView x,
		 ConstVectorView y,
		 ConstVectorView xa,
		 ConstMatrixView K,
		 ConstMatrixView Se,
		 ConstMatrixView Sa,
                 bool mform);

#endif // oem_h
