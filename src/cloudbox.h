/* Copyright (C) 2002 Claudia Emde <claudia@sat.physik.uni-bremen.de>
                      
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. 
*/

 /*!
     \file   cloudbox.h
     \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
     \date   Thu May  23 14:34:05 2002
     
     \brief  Internal cloudbox functions.
     
   */

#ifndef cloudbox_h
#define cloudbox_h

#include <stdexcept>
#include <math.h>
#include "lin_alg.h"
#include "arts.h"
#include "auto_md.h"
#include "matpackI.h"
#include "make_vector.h"
#include "array.h"
#include "logic.h"
#include "ppath.h"
#include "interpolation.h"
#include "physics_funcs.h"



void sto_vecCalc(VectorView sto_vec,
		 ConstMatrixView ext_mat,
		 ConstVectorView abs_vec,
		 ConstVectorView sca_vec,
		 const Numeric& ds,
		 const Numeric& B,
		 const Index& stokes_dim);

void sto_vec1DCalc(VectorView sto_vec,
		   ConstMatrixView ext_mat,
		   ConstVectorView abs_vec,
		   ConstVectorView sca_vec,
		   const Numeric& ds,
		   const Numeric& B,
		   const Index& stokes_dim);


#endif    // cloudbox_h
 
