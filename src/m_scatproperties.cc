/*!
  \file   m_scatproperties.cc
  \author Sreerekha T.R. <rekha@uni-bremen.de>
  \date   \date   Mon Jun 10 11:19:11 2002
  \brief  This file contains workspace methods for calculating the single scattering properties of ice particles.
*/
/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include "arts.h"
#include "array.h"
#include "check_input.h"
#include "matpackI.h"
#include "matpackVI.h"
#include "matpackVII.h"
#include "scatproperties.h"
#include "logic.h"

//! This method computes the extinction matrix for a single particle type
//  from teh amplitude matrix.
/*! 
  
  \param ext_mat_spt Output and Input: extinction matrix for a single particle type.
  \param amp_mat Input : amplitude matrix for each particle type
  \param f_grid  Input : frequency grid
  \param scat_f_index  Input : frequency index
  \param scat_za_index  Input : local zenith angle
  \param scat_aa_index  Input : local azimuth angle
 */
void ext_mat_sptCalc(Tensor3& ext_mat_spt,
		     const Tensor6& amp_mat,
		     const Index& scat_za_index,
		     const Index& scat_aa_index,
		     const Index& scat_f_index,
		     const Vector& f_grid
		     )
{
  Index npt = ext_mat_spt.npages();
  Index stokes_dim = ext_mat_spt.nrows();
  if (amp_mat.ncols() != 8)
    {
      throw runtime_error("amplitude matrix has 8 columns");
    }
  if (stokes_dim > 4 || stokes_dim <1){
    throw runtime_error("the dimension of stokes vector can be only 1,2,3, or 4");
  }
  //find out frequency
  Numeric freq = f_grid[scat_f_index];
  for (Index i=0; i<npt; ++i)
    {
      ConstVectorView amp_coeffs = amp_mat(i,
				      scat_za_index,
				      scat_aa_index,
				      scat_za_index,
				      scat_aa_index,
				      Range(joker)
				      );
      MatrixView ext = ext_mat_spt(i,
				   Range(joker),
				   Range(joker)
				   );
      amp2ext(ext,amp_coeffs,freq);
    }
}
//! this function calculates phase matrix for a single particle type from 
//  the amplitude matrix
/*! 
  
  \param pha_mat_spt  Output and Input: phase matrix for a single particle type
  \param amp_mat  Input : amplitude matrix
  \param scat_za_index  Input : zenith angle index
  \param scat_aa_index  Input : azimuth angle index
  
*/
void pha_mat_sptCalc(Tensor5& pha_mat_spt,
		     const Tensor6& amp_mat,
		     const Index& scat_za_index,
		     const Index& scat_aa_index
		     )
{
  Index npt = pha_mat_spt.nshelves();
  Index stokes_dim = pha_mat_spt.nrows();
  if (pha_mat_spt.nrows() != stokes_dim || pha_mat_spt.ncols() != stokes_dim ){
    throw runtime_error(" pha_mat_spt : should have 4 rows and 4 columns");
  }
  if (amp_mat.ncols() != 8){
    throw runtime_error("amplitude matrix has 8 columns");
  }
  if (stokes_dim > 4 || stokes_dim <1){
    throw runtime_error("the dimension of stokes vector can be only 1,2,3, or 4");
  }
  for (Index i=0; i<npt; ++i)
    {
      const Tensor3 amp_coeffs = amp_mat(i,
					   scat_za_index,
					   scat_aa_index,
					   Range(joker),
					   Range(joker),
					   Range(joker)
					   );
      amp2pha(pha_mat_spt(i,Range(joker),Range(joker),Range(joker),Range(joker)),
	      amp_coeffs);
    }
}
//! this function calculates the absorption coefficient for a single particle type
// from the extinction matrix and phase matrix.
/*! 
  
  \param abs_vec_spt  Output : absorption vector for a single particle type
  \param ext_mat_spt  Input : extinction matrix for a single particle type
  \param pha_mat_spt  Input : phase matrix for  a single particle type
 
*/
void abs_vec_sptCalc(Matrix& abs_vec_spt,
		     const Tensor3& ext_mat_spt,
		     const Tensor5& pha_mat_spt
		     )
{
 Index npt = abs_vec_spt.nrows();
 Index stokes_dim = abs_vec_spt.ncols();
 if (abs_vec_spt.ncols() != stokes_dim ){
    throw runtime_error(" abs_vec_spt : should have 4 columns");
 }
 if (ext_mat_spt.nrows() != stokes_dim || ext_mat_spt.ncols() != stokes_dim ){
   throw runtime_error(" ext_mat_spt : should have 4 rows and 4 columns");
 }
 if (pha_mat_spt.nrows() != stokes_dim || pha_mat_spt.ncols() != stokes_dim ){
   throw runtime_error(" pha_mat_spt : should have 4 rows and 4 columns");
 }
 if (stokes_dim > 4 || stokes_dim <1){
    throw runtime_error("the dimension of stokes vector can be only 1,2,3, or 4");
 }
 for (Index i=0;i<npt;++i)
   {
     ConstMatrixView ext = ext_mat_spt(i,
				  Range(joker),
				  Range(joker));
     ConstTensor4View pha=pha_mat_spt(i,
				 Range(joker),
				 Range(joker),
				 Range(joker),
				 Range(joker));
     amp2abs(abs_vec_spt(i,Range(joker)),ext,pha);
   }	   
}
