/* Copyright (C)  Claudia Emde <claudia@sat.physik.uni-bremen.de>

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
   USA. */



/*===========================================================================
  === File description
  ===========================================================================*/

/*!
  \file   optproperties.cc
  \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
  \date   Thu Mar  6 11:29:59 2003
  
  \brief  This file contains definitions and functions related to the
          optical properties of particles.
  
  
 */


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "matpackVII.h"
#include "array.h"
#include "math_funcs.h"
#include "messages.h"
#include "logic.h"
#include "interpolation.h"
#include "optproperties.h"


//! Transformation of absorption vector.
/*! 
 
Documentation will be written (CE). 
 
*/
void abs_vecTransform(//Output and Input
                      VectorView abs_vec_lab,
                      //Input
                      const Tensor3View abs_vec_data,
                      const VectorView za_datagrid,
                      const VectorView aa_datagrid,
                      const PType& ptype,
                      const Numeric& za_sca,
                      const Numeric& aa_sca)
{
 const Index stokes_dim = abs_vec_lab.nelem();
    
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                         "must be 1,2,3 or 4");
  }

  switch (ptype){

  case PTYPE_GENERAL:
    cout << "Case PTYPE_GENERAL not yet implemented. \n"; 
    break;
    
  case PTYPE_MACROS_ISO:
    {
      // The first element of the vector corresponds to the absorption 
      // coefficient which is stored in the database, the others are 0.

      abs_vec_lab[0] = abs_vec_data(0,0,0);
      break;
    }
   default:
    cout << "Not all particle type cases are implemented\n";
    
  }  
    
}


//! Transformation of extinction matrix.
/*! 
 
Documentation will be written (CE). 
 
*/
void ext_matTransform(//Output and Input
                      MatrixView ext_mat_lab,
                      //Input
                      const Tensor3View ext_mat_data,
                      const VectorView za_datagrid,
                      const VectorView aa_datagrid,
                      const PType& ptype,
                      const Numeric& za_sca,
                      const Numeric& aa_sca)
{
 const Index stokes_dim = ext_mat_lab.ncols();
    
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                         "must be 1,2,3 or 4");
  }

  switch (ptype){

  case PTYPE_GENERAL:
    cout << "Case PTYPE_GENERAL not yet implemented. \n"; 
    break;
    
  case PTYPE_MACROS_ISO:
    {
      assert (ext_mat_data.ncols() == 1);
      
      // In the case of randomly oriented particles the extinction matrix is 
      // diagonal. The value of each element of the diagonal is the
      // extinction cross section, which is stored in the database.
      
      ext_mat_lab(0,0) = ext_mat_data(0,0,0);
      
      if( stokes_dim == 1 ){
        break;
      }
      
      ext_mat_lab(1,1) = ext_mat_data(0,0,0);
      
      if( stokes_dim == 2 ){
        break;
      }
      
      ext_mat_lab(2,2) = ext_mat_data(0,0,0);
      
      if( stokes_dim == 3 ){
        break;
      }
      
      ext_mat_lab(3,3) = ext_mat_data(0,0,0);
      break;
    }

  default:
    cout << "Not all particle type cases are implemented\n";
    
  }
}  
 

//! Transformation of phase matrix.
/*! 
 
Documentation will be written (CE). 
 
*/
void pha_matTransform(//Output
                      MatrixView pha_mat_lab,
                      //Input
                      const Tensor5View pha_mat_data,
                      const VectorView za_datagrid,
                      const VectorView aa_datagrid,
                      const PType& ptype,
                      const Numeric& za_sca,
                      const Numeric& aa_sca,
                      const Numeric& za_inc,
                      const Numeric& aa_inc)
{
  
  const Index stokes_dim = pha_mat_lab.ncols();
    
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                         "must be 1,2,3 or 4");
  }

  switch (ptype){

  case PTYPE_GENERAL:
    cout << "Case PTYPE_GENERAL not yet implemented. \n"; 
    break;
    
  case PTYPE_MACROS_ISO:
    {
    assert (pha_mat_data.ncols() == 6);

    // The transformation formulas can be found in Mishchenkho (p.90)
    // (book: Scattering, absorption and emission ...)
    
    // Calculation of the scattering angle:
    const Numeric theta = acos(cos(za_sca)*cos(za_inc) + 
                                 sin(za_sca)*sin(za_inc)*cos(aa_sca-aa_inc));
    
    // Interpolate the data on the scattering angle:
    Vector pha_mat_int(6);

    GridPos thet_gp;
    gridpos(thet_gp, za_datagrid, theta);
    
    Vector itw(2);
    interpweights(itw, thet_gp);

    for (Index i = 0; i < 6; i++)
      {
        pha_mat_int[i] = interp(itw, pha_mat_data(joker, 0, 0, 0, i), 
               thet_gp);
      }

   
    // Scattering matrix elements:

    const Numeric F11 = pha_mat_int[0];
    const Numeric F12 = pha_mat_int[1];
    const Numeric F22 = pha_mat_int[2];
    const Numeric F33 = pha_mat_int[3];
    const Numeric F34 = pha_mat_int[4];
    const Numeric F44 = pha_mat_int[5];
    
    // Calculation of the phase matrix elements in the laboratory 
    // coordinate system:
    //
    pha_mat_lab(0,0) = F11;

    if( stokes_dim == 1 ){
      break;
    }

    const Numeric sigma1 = acos( (cos(za_sca) - cos(za_inc)*cos(theta))/
                                 (sin(za_inc)*sin(theta)));
    
    const Numeric sigma2 = acos( (cos(za_inc) - cos(za_sca)*cos(theta))/
                                 (sin(za_sca)*sin(theta)));
    
    const Numeric C1 = cos( 2*sigma1 );
    const Numeric C2 = cos( 2*sigma2 );

    const Numeric S1 = sin( 2*sigma1 );
    const Numeric S2 = sin( 2*sigma2 );

    pha_mat_lab(0,1) = C1 * F12;
    pha_mat_lab(1,0) = C2 * F12;
    pha_mat_lab(1,1) = C1 * C2 * F22 - S1 * S2 * F33;

    if( stokes_dim == 2 ){
      break;
    }


    pha_mat_lab(0,2) = -S1 * F12;
    pha_mat_lab(1,2) = -S1 * C2 * F22 - C1 * S2 * F33;
    pha_mat_lab(2,0) = S2 * F12;
    pha_mat_lab(2,1) = C1 * S2 * F22 + S1 * C2 * F33;
    pha_mat_lab(2,2) = -S1 * S2 * F22 + C1 * C2 * F33;

    if( stokes_dim == 3 ){
      break;
    }

    pha_mat_lab(0,3) = 0;
    pha_mat_lab(1,3) = -S2 * F34;
    pha_mat_lab(2,3) = C2 * F34;
    pha_mat_lab(3,0) = 0;
    pha_mat_lab(3,1) = -S1 * F34;
    pha_mat_lab(3,2) = -C1 * F34;
    pha_mat_lab(3,3) = F44;
    
    break;
    }
  default:
    cout << "Not all particle type cases are implemented\n";
      
      }

} 

