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
#include "xml_io.h"
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;

//! Transformation of absorption vector.
/*! 
  In the single scattering database the data of the absorption vector is 
  stored in different coordinate systems, depending on the type of 
  hydrometeor species.

  See AUG for information about different classifications of 
  the hydrometeor species. 

  Output and Input:
  \param abs_vec_lab Absorption vector in Laboratory frame.
  Input:
  \param abs_vec_data Absorption vector in database.
  \param za_datagrid Zenith angle grid in the database.
  \param aa_datagrid Zenith angle grid in the database.
  \param ptype Clasiification of the hydometeor species.
  \param za_sca Zenith angle of scattered direction.
  \param aa_sca Azimuth angle of scattered direction.
     
  \author Claudia Emde
  \date   2003-05-24 
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
      
      abs_vec_lab = 0;

      abs_vec_lab[0] = abs_vec_data(0,0,0);
      break;
    }
   default:
    cout << "Not all particle type cases are implemented\n";
    
  }  
    
}


//! Transformation of extinction matrix.
/*! 
  In the single scattering database the data of the extinction matrix is 
  stored in different coordinate systems, depending on the type of 
  hydrometeor species.

  See AUG for information about different classifications of 
  the hydrometeor species. 

  Output and Input:
  \param ext_mat_lab Absorption vector in Laboratory frame.
  Input:
  \param ext_mat_data Absorption vector in database.
  \param za_datagrid Zenith angle grid in the database.
  \param aa_datagrid Zenith angle grid in the database.
  \param ptype Clasiification of the hydometeor species.
  \param za_sca Zenith angle of scattered direction.
  \param aa_sca Azimuth angle of scattered direction.

  \author Claudia Emde
  \date   2003-05-24 
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
     
      ext_mat_lab = 0.;
      
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
  In the single scattering database the data of the phase matrix is 
  stored in different coordinate systems, depending on the type of 
  hydrometeor  species.

  See AUG for information about different classifications of 
  the hydrometeor species. 

  Output and Input:
  \param ext_mat_lab Absorption vector in Laboratory frame.
  Input:
  \param ext_mat_data Absorption vector in database.
  \param za_datagrid Zenith angle grid in the database.
  \param aa_datagrid Zenith angle grid in the database.
  \param ptype Clasiification of the hydometeor species.
  \param za_sca Zenith angle of scattered direction.
  \param aa_sca Azimuth angle of scattered direction.
  \param za_inc Zenith angle of incoming direction.
  \param aa_inc Azimuth angle of incoming direction.

  \author Claudia Emde
  \date   2003-05-24
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

  const Numeric za_sca_rad = za_sca * DEG2RAD;
  const Numeric za_inc_rad = za_inc * DEG2RAD;
  const Numeric aa_sca_rad = aa_sca * DEG2RAD;
  const Numeric aa_inc_rad = aa_inc * DEG2RAD;
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
      // Calculate the scattering and interpolate the data on the scattering
      // angle:
      Vector pha_mat_int(6);
      Numeric theta_rad;

      // Interpolation of the data on the scattering angle:
      interpolate_scat_angle(pha_mat_int, theta_rad, pha_mat_data,
                             za_datagrid, 
                             aa_datagrid, za_sca_rad, aa_sca_rad,
                             za_inc_rad, aa_inc_rad);

      // Caclulate the phase matrix in the laboratory frame:
      pha_mat_labCalc(pha_mat_lab, pha_mat_int, za_sca, aa_sca, za_inc, 
                      aa_inc, theta_rad);
      
      break;
    }
    
  default:
    cout << "Not all particle type cases are implemented\n";
    
  }
}


//! Interpolate data on the scattering angle.
/*! 
  This function is used for the transformation of the phase matrix 
  from scattering frame to the laboratory frame for randomly oriented
  scattering media (case PTYPE_MACRO_ISO).

  The scattering angle is calculated from the angles defining
  the directions of the incoming and scattered radiation.
  After that the data (which is stored in the data files as a function
  of the scattering angle) is interpolated on the calculated 
  scattering angle.

  Output:
  \param pha_mat_int Interpolated phase matrix.
  \param theta_rad Scattering angle [rad].
  Input:
  \param pha_mat_data Phase matrix in database.
  \param za_datagrid Zenith angle grid in the database.
  \param aa_datagrid Zenith angle grid in the database.
  \param za_sca_rad Zenith angle of scattered direction [rad].
  \param aa_sca_rad Azimuth angle of scattered direction [rad].
  \param za_inc_rad Zenith angle of incoming direction [rad].
  \param aa_inc_rad Azimuth angle of incoming direction [rad].
     
  \author Claudia Emde
  \date   2003-05-13 
*/
void interpolate_scat_angle(//Output:
                            VectorView pha_mat_int,
                            Numeric& theta_rad,
                            //Input:
                            const Tensor5View pha_mat_data,
                            const VectorView za_datagrid,
                            const VectorView aa_datagrid,
                            const Numeric& za_sca_rad,
                            const Numeric& aa_sca_rad,
                            const Numeric& za_inc_rad,
                            const Numeric& aa_inc_rad)
{
  // cout << "Interpolation on scattering angle" << endl;
  assert (pha_mat_data.ncols() == 6);

  // Calculation of the scattering angle:
  theta_rad = acos(cos(za_sca_rad) * cos(za_inc_rad) + 
                   sin(za_sca_rad) * sin(za_inc_rad) * 
                   cos(aa_sca_rad - aa_inc_rad));
  
  const Numeric theta = RAD2DEG * theta_rad;
  
  // Interpolation of the data on the scattering angle:
  GridPos thet_gp;
  gridpos(thet_gp, za_datagrid, theta);
    
  Vector itw(2);
  interpweights(itw, thet_gp);

  for (Index i = 0; i < 6; i++)
    {
      pha_mat_int[i] = interp(itw, pha_mat_data(joker, 0, 0, 0, i), 
                              thet_gp);
    }
} 



//! Calculate phase matrix in laboratory coordinate system.
/*! 
  Transformation function for the phase matrix for the case of
  randomly oriented particles (case PTYPE_MACRO_ISO).
  
  Some of the formulas can be found in 

  Mishchenkho: "Scattering, Absorption and Emission of Light 
  by Small Particles", Cambridge University Press, 2002
  Capter 4

  The full set of formulas will be documented in AUG.

  Output and Input:
  \param pha_mat_lab Phase matrix in laboratory frame.
  Input: 
  \param pha_mat_int Interpolated phase matrix.
  \param za_sca Zenith angle of scattered direction.
  \param aa_sca Azimuth angle of scattered direction.
  \param za_inc Zenith angle of incoming direction.
  \param aa_inc Azimuth angle of incoming direction.
  \param theta_rad Scattering angle [rad].
  
  \author Claudia Emde
  \date   2003-05-13 
*/
void pha_mat_labCalc(//Output:
                      MatrixView pha_mat_lab,
                      //Input:
                      const VectorView pha_mat_int,
                      const Numeric& za_sca,
                      const Numeric& aa_sca,
                      const Numeric& za_inc,
                      const Numeric& aa_inc,
                      const Numeric& theta_rad)
{
  Numeric za_sca_rad = za_sca * DEG2RAD;
  Numeric za_inc_rad = za_inc * DEG2RAD;
  Numeric aa_sca_rad = aa_sca * DEG2RAD;
  Numeric aa_inc_rad = aa_inc * DEG2RAD;
  const Numeric theta = RAD2DEG * theta_rad;
  const Index stokes_dim = pha_mat_lab.ncols();

  // cout << "Transformation of phase matrix:" <<endl; 
  
  // Scattering matrix elements:
  const Numeric F11 = pha_mat_int[0];
  const Numeric F12 = pha_mat_int[1];
  const Numeric F22 = pha_mat_int[2];
  const Numeric F33 = pha_mat_int[3];
  const Numeric F34 = pha_mat_int[4];
  const Numeric F44 = pha_mat_int[5];
  
  // For stokes_dim = 1, we only need Z11=F11:
  pha_mat_lab(0,0) = F11;
  
  if( stokes_dim > 1 ){
    //
    // Several cases have to be considered:
    //
   if(
        // Forward scattering
        ((theta > -.01) && (theta < .01) ) ||
        // Backward scattering
        ((theta > 179.99) && (theta < 180.01)) 
        )
      {
        pha_mat_lab(0,1) = F12;
        pha_mat_lab(1,0) = F12;
        pha_mat_lab(1,1) = F22;
        
        if( stokes_dim > 2 ){
          pha_mat_lab(0,2) = 0;
          pha_mat_lab(1,2) = 0;
          pha_mat_lab(2,0) = 0;
          pha_mat_lab(2,1) = 0;
          pha_mat_lab(2,2) = F33;
          
          if( stokes_dim > 3 ){
            pha_mat_lab(0,3) = 0;
            pha_mat_lab(1,3) = 0;
            pha_mat_lab(2,3) = F34;
            pha_mat_lab(3,0) = 0;
            pha_mat_lab(3,1) = 0;
            pha_mat_lab(3,2) = -F34;
            pha_mat_lab(3,3) = F44;
          }
        }
      }
   //  else if(// Scattering frame equals laboratory frame
//             (za_inc == 0) || (za_sca == 180) ||
//         // Scattering frame is "mirrored" laboratory frame
//         (za_sca == 0) || (za_inc == 180) ||
//         // "Grosskreis"
//         (aa_sca == aa_inc) || (aa_sca == aa_inc-360) || 
//             (aa_inc == aa_sca-360) )
//       {
//         // FIXME: Wich formulas do we need in these special cases. The values
//         // are very small, so the overall error is small ... but this needs to
//         // be fixed
//         pha_mat_lab = 0;
//       } 
    else if( (aa_sca - aa_inc) > 0 && (aa_sca - aa_inc) < 180 ||  
             (aa_sca - aa_inc) > -360 && (aa_sca - aa_inc) < -180 )
      {
        // In these cases we have to take limiting value
        // (according personal communication with Mishchenko)
        if (za_inc_rad == 0)
          za_inc_rad = 1e-6;
        if (za_inc_rad == PI)
           za_inc_rad = PI - 1e-6;
        if (za_sca_rad == 0)
          za_sca_rad = 1e-6; 
        if (za_sca_rad == PI)
           za_sca_rad = PI - 1e-6;
        
        const Numeric cos_sigma1 =  (cos(za_sca_rad) - cos(za_inc_rad)
                                      * cos(theta_rad))/
                                 (sin(za_inc_rad)*sin(theta_rad));
        const Numeric cos_sigma2 =  (cos(za_inc_rad) - cos(za_sca_rad)
                                      *cos(theta_rad))/
                                     (sin(za_sca_rad)*sin(theta_rad));
        
               
        const Numeric C1 = 2 * cos_sigma1 * cos_sigma1 - 1;
        const Numeric C2 = 2 * cos_sigma2 * cos_sigma2 - 1;
        
        const Numeric S1 = 2 * sqrt(1 - cos_sigma1) * cos_sigma1;
        const Numeric S2 = 2 * sqrt(1 - cos_sigma2) * cos_sigma2;
        
        pha_mat_lab(0,1) = C1 * F12;
        pha_mat_lab(1,0) = C2 * F12;
        pha_mat_lab(1,1) = C1 * C2 * F22 - S1 * S2 * F33;
        
        if( stokes_dim > 2 ){
                            
          pha_mat_lab(0,2) = S1 * F12;
          pha_mat_lab(1,2) = S1 * C2 * F22 + C1 * S2 * F33;
          pha_mat_lab(2,0) = -S2 * F12;
          pha_mat_lab(2,1) = -C1 * S2 * F22 - S1 * C2 * F33;
          pha_mat_lab(2,2) = -S1 * S2 * F22 + C1 * C2 * F33;
          
          if( stokes_dim > 3 ){
            
            pha_mat_lab(0,3) = 0;
            pha_mat_lab(1,3) = S2 * F34;
            pha_mat_lab(2,3) = C2 * F34;
            pha_mat_lab(3,0) = 0;
            pha_mat_lab(3,1) = S1 * F34;
            pha_mat_lab(3,2) = -C1 * F34;
            pha_mat_lab(3,3) = F44;
          }
        }     
  }
  else if ( (aa_sca - aa_inc) > -180 && (aa_sca - aa_inc) < 0 ||
              (aa_sca - aa_inc) > 180 && (aa_sca - aa_inc) < 360 )
      {
        // In these cases we have to take limiting value
        // (according personal communication with Mishchenko)
        if (za_inc_rad == 0)
          za_inc_rad = 1e-6;
        if (za_inc_rad == PI)
           za_inc_rad = PI - 1e-6;
        if (za_sca_rad == 0)
          za_sca_rad = 1e-6; 
        if (za_sca_rad == PI)
           za_sca_rad = PI - 1e-6;
        
        const Numeric cos_sigma1 =  (cos(za_sca_rad) - cos(za_inc_rad)
                                     * cos(theta_rad))/
                                     (sin(za_inc_rad)*sin(theta_rad));
        const Numeric cos_sigma2 =  (cos(za_inc_rad) - cos(za_sca_rad)
                                     *cos(theta_rad))/
          (sin(za_sca_rad)*sin(theta_rad));
        
        const Numeric C1 = 2 * cos_sigma1 * cos_sigma1 - 1;
        const Numeric C2 = 2 * cos_sigma2 * cos_sigma2 - 1;
        
        const Numeric S1 = 2 * sqrt(1 - cos_sigma1) * cos_sigma1;
        const Numeric S2 = 2 * sqrt(1 - cos_sigma2) * cos_sigma2;
        
        pha_mat_lab(0,1) = C1 * F12;
        pha_mat_lab(1,0) = C2 * F12;
        pha_mat_lab(1,1) = C1 * C2 * F22 - S1 * S2 * F33;
        
        if( stokes_dim > 2 ){
          pha_mat_lab(0,2) = S1 * F12;
          pha_mat_lab(1,2) = S1 * C2 * F22 + C1 * S2 * F33;
          pha_mat_lab(2,0) = -S2 * F12;
          pha_mat_lab(2,1) = -C1 * S2 * F22 + S1 * C2 * F33;
          pha_mat_lab(2,2) = -S1 * S2 * F22 + C1 * C2 * F33;
          
          if( stokes_dim > 3 ){
            pha_mat_lab(0,3) = 0;
            pha_mat_lab(1,3) = S2 * F34;
            pha_mat_lab(2,3) = C2 * F34;
            pha_mat_lab(3,0) = 0;
            pha_mat_lab(3,1) = S1 * F34;
            pha_mat_lab(3,2) = -C1 * F34;
            pha_mat_lab(3,3) = F44;
          }
        }
      }
   }
}
     
