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


void i_field_update1D(
		     Tensor6View i_field,
		     ConstTensor6View i_field_old,
		     ConstTensor6View amp_mat,
		     ConstTensor6View sca_field,
		     const ArrayOfIndex cloudbox_limits,
		     ConstVectorView scat_za_grid,
		     ConstVectorView scat_aa_grid,
		     ConstVectorView p_grid,
		     ConstVectorView lat_grid,
		     ConstVectorView lon_grid,
		     ConstTensor3View t_field,
		     ConstTensor3View z_field,
		     ConstMatrixView z_ground,
		     ConstMatrixView r_geoid,
		     const Numeric f,
		     const Index blackbody_ground,
		     const Index stokes_dim
		     );

void rte_scat_vecCalc(VectorView sto_vec,
		      ConstMatrixView ext_mat,
		      ConstVectorView abs_vec,
		      ConstVectorView sca_vec,
		      const Numeric ds,
		      const Numeric B);



#endif    // cloudbox_h
 
