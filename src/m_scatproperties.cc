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
#include "auto_md.h"

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
    throw runtime_error(
			"Amplitude matrix must have 8 columns.");
  
  if (stokes_dim > 4 || stokes_dim <1){
    throw runtime_error(
			"The dimension of stokes vector can be "
			"only 1,2,3, or 4");
  }

  //find out frequency
  Numeric freq = f_grid[scat_f_index];
  
  for (Index i = 0; i < npt; ++i)
    {
      ConstVectorView amp_coeffs = amp_mat(i,
				      scat_za_index,
				      scat_aa_index,
				      scat_za_index,
				      scat_aa_index,
				      Range(joker)
				      );
      
      amp2ext(ext_mat_spt(i,Range(joker),Range(joker)),
	      amp_coeffs,
	      freq);
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

  if (pha_mat_spt.nrows() != stokes_dim || 
      pha_mat_spt.ncols() != stokes_dim){
    throw runtime_error(" The tensor pha_mat_spt should have"
			"4 rows and 4 columns.");
  }
  if (amp_mat.ncols() != 8){
    throw runtime_error("Amplitude matrix must have 8 columns.");
  }
  if (stokes_dim > 4 || stokes_dim <1){
    throw runtime_error("The dimension of stokes vector" 
                         "can be only 1,2,3, or 4");
  }
  for (Index i = 0; i < npt; ++i)
    {
      const Tensor3 amp_coeffs = amp_mat(i,
					 scat_za_index,
					 scat_aa_index,
					 Range(joker),
					 Range(joker),
					 Range(joker)
					 );
      amp2pha(pha_mat_spt(i, Range(joker), Range(joker), 
			  Range(joker), Range(joker)),
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
		     const Tensor5& pha_mat_spt,
		     const Vector& scat_za_grid,
		     const Vector& scat_aa_grid
		     )
{
 Index npt = abs_vec_spt.nrows();
 Index stokes_dim = abs_vec_spt.ncols();
 Index nza = pha_mat_spt.nshelves(); 
 Index naa = pha_mat_spt.nbooks(); 
 assert (is_size(scat_za_grid, nza));
 assert (is_size(scat_aa_grid, naa));

 if (abs_vec_spt.ncols() != stokes_dim ){
    throw runtime_error("The vector abs_vec_spt should have 4 columns");
 }

 if (ext_mat_spt.nrows() != stokes_dim || 
     ext_mat_spt.ncols() != stokes_dim ){
   throw runtime_error("The tensor ext_mat_spt should"
		       " have 4 rows and 4 columns");
 }

 if (pha_mat_spt.nrows() != stokes_dim || 
     pha_mat_spt.ncols() != stokes_dim ){
   throw runtime_error("The tensor pha_mat_spt should"
		       "have 4 rows and 4 columns");
 }

 if (stokes_dim > 4 || stokes_dim <1){
    throw runtime_error("The dimension of stokes vector "
                        "can be only 1,2,3, or 4");
 }

 for (Index i = 0; i < npt; ++i)
   {
     ConstMatrixView ext = ext_mat_spt(i,
				       Range(joker),
				       Range(joker));
     
     ConstTensor4View pha=pha_mat_spt(i,
				      Range(joker),
				      Range(joker),
				      Range(joker),
				      Range(joker));
    
     // ConstVectorView za_grid = scat_za_grid[scat_za_index];
     //ConstVectorView aa_grid = scat_aa_grid[scat_aa_index];
     
     amp2abs(abs_vec_spt(i, Range(joker)),
	     ext,
	     pha,
	     scat_za_grid,
	     scat_aa_grid);
   }	   
}


//! Extinction Coefficient Matrix for the particle 
/*! 
  
  This function sums up the extinction matrices for all particle 
  types weighted with particle number density
  \param ext_mat_part Output : physical extinction coefficient 
  for the particles for given angles. 
  \param ext_mat_spt Input : extinction matrix for the single particle type
  \param pnd_field Input : particle number density givs the local 
  concentration for all particles.
  \param cloudbox_limits Input : the limits of the cloud box.
  \param atmosphere_dim Input : he atmospheric dimensionality (now 1 or 3)
*/
void ext_mat_partCalc(Matrix& ext_mat_part,
		      const Tensor3& ext_mat_spt,
		      const Tensor4& pnd_field,
		      const ArrayOfIndex& cloudbox_limits,
		      const Index& atmosphere_dim
		     )
{
  Index N_pt = ext_mat_spt.npages();
  Index N_i = ext_mat_spt.nrows();
   
  for (Index m = 0; m < N_i; m++)
    {
      for (Index n = 0; n < N_i; n++)
	ext_mat_part(m, n) = 0.0;// Initialisation
    }
  if (atmosphere_dim == 1)
    {
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt; l++)
	{ 
	  
	  // the loop are over the p_grid
	  for (Index i = cloudbox_limits[0]; i < cloudbox_limits[1] ; ++i)
	    {
	      
	      // now the last two loops over the stokes dimension.
	      for (Index m = 0; m < N_i; m++)
		{
		  for (Index n = 0; n < N_i; n++)
		    
		    ext_mat_part(m, n) += 
		      (ext_mat_spt(l, m, n) * pnd_field(l, i, 0, 0));
		}
	    }
	}
    }
	  
  if (atmosphere_dim == 3)
    {

       // this is a loop over the different particle types
		  for (Index l = 0; l < N_pt; l++)
		    { 
		      
      // the next three loops are over the p_grid, lat_grid, lon_grid
      for (Index i = cloudbox_limits[0]; i < cloudbox_limits[1] ; ++i)
	{
	  for (Index j = cloudbox_limits[2]; j < cloudbox_limits[3] ; ++j)
	    {
	      for (Index k = cloudbox_limits[4]; k < cloudbox_limits[5] ; ++k)
		{
		  
		 
		      // now the last two loops over the stokes dimension.
		      for (Index m = 0; m < N_i; m++)
			{
			  for (Index n = 0; n < N_i; n++)
			    
			    ext_mat_part(m, n) += 
			      (ext_mat_spt(l, m, n) * pnd_field(l, i, j, k));
			  
			} 
		    }
		}
	    }
	}
    }
} 

//! Aabsorption Vector for the particle 
/*! 
  
  This function sums up the absorption vectors for all particle 
  types weighted with particle number density
  \param abs_vec_part Output : physical absorption vector 
  for the particles for given angles. 
  \param abs_vec_spt Input : absorption for the single particle type
  \param pnd_field Input : particle number density givs the local 
  concentration for all particles.
  \param cloudbox_limits Input : the limits of the cloud box.
  \param atmosphere_dim Input : he atmospheric dimensionality (now 1 or 3)
*/
void abs_vec_partCalc(Vector& abs_vec_part,
		      const Matrix& abs_vec_spt,
		      const Tensor4& pnd_field,
		      const ArrayOfIndex& cloudbox_limits,
		      const Index& atmosphere_dim
		      )
{
  Index N_pt = abs_vec_spt.nrows();
  Index N_i = abs_vec_spt.ncols();
 
 
  for (Index m = 0; m < N_i; ++m)
    {
      
      abs_vec_part[m] = 0.0;// Initialisation
    }
  if (atmosphere_dim == 1)
    {
      // this is a loop over the different particle types
	  for (Index l = 0; l < N_pt ; ++l)
	    {
      // the first loop over the p_grid
      for (Index i = cloudbox_limits[0]; i < cloudbox_limits[1] ; ++i)
	{
	  
	  
	      
	      // now the loop over the stokes dimension.
	      for (Index m = 0; l < N_i; ++m)
		
		abs_vec_part[m] += 
		  (abs_vec_spt(l, m) * pnd_field(l, i, 0, 0));
	      
	    }
	}
    }
  
  if (atmosphere_dim == 3)
    {
       // this is a loop over the different particle types
		  for (Index l = 0; l < N_pt ; ++l)
		    {
		      
      // the next three loops are over the p_grid, lat_grid, lon_grid
      for (Index i = cloudbox_limits[0]; i < cloudbox_limits[1] ; ++i)
	{
	  for (Index j = cloudbox_limits[2]; j < cloudbox_limits[3] ; ++j)
	    {
	      for (Index k = cloudbox_limits[4]; k < cloudbox_limits[5] ; ++k)
		{
		  
		 
		      // now the loop over the stokes dimension.
		      for (Index m = 0; l < N_i; ++m)
			
			abs_vec_part[m] += 
			  (abs_vec_spt(l, m) * pnd_field(l, i, j, k));
		      
		    }
		}
	    }
	}
    }
} 


//! Phase Matrix for the particle 
/*! 
  
  This function sums up the phase matrices for all particle 
  types weighted with particle number density
  \param pha_mat_part Output : physical phase matrix 
  for the particles for given angles. 
  \param pha_mat_spt Input : phase matrix for the single particle type
  \param pnd_field Input : particle number density givs the local 
  concentration for all particles.
  \param cloudbox_limits Input : the limits of the cloud box.
  \param atmosphere_dim Input : he atmospheric dimensionality (now 1 or 3)
*/

void pha_mat_partCalc(Tensor4& pha_mat_part,
		      const Tensor5& pha_mat_spt,
		      const Tensor4& pnd_field,
		      const ArrayOfIndex& cloudbox_limits,
		      const Index& atmosphere_dim
		     )
{

  Index N_pt = pha_mat_spt.nshelves();
  Index Nza = pha_mat_spt.nbooks();
  Index Naa = pha_mat_spt.npages();
  Index N_i = pha_mat_spt.nrows();

  // Initialisation
  for (Index za_index = 0; za_index < Nza; ++ za_index)
    {
      for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
	{
	  for (Index stokes_index = 0; stokes_index < N_i; ++ stokes_index)
	    {
	      for (Index stokes_index = 0; stokes_index < N_i; 
		   ++ stokes_index)
		pha_mat_part(za_index, aa_index, stokes_index, stokes_index)
		  = 0.0;
	    }
	}
    }
  if (atmosphere_dim == 1)
    {
      // this is a loop over the different particle types
	  for (Index pt_index = 0; pt_index < N_pt; ++ pt_index)
	    {
      // the loop are over the p_grid
      for (Index p_index = cloudbox_limits[0];  p_index < cloudbox_limits[1];
	   ++  p_index)
	{
	 
	      // these are loops over zenith angle and azimuth angle
	      for (Index za_index = 0; za_index < Nza; ++ za_index)
		{
		  for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
		    {
	      
		      // now the last two loops over the stokes dimension.
		      for (Index stokes_index = 0; stokes_index < N_i; 
			   ++  stokes_index)
			{
			  for (Index stokes_index = 0; stokes_index < N_i;
			       ++ stokes_index)
			    
			    pha_mat_part(za_index, aa_index,  
					 stokes_index, stokes_index) += 
			      
			      (pha_mat_spt(pt_index, za_index, aa_index,  
					   stokes_index, stokes_index) * 
			       pnd_field(pt_index, p_index, 0, 0));
			}
		    }
		}
	    }
	}
    }
	  
  if (atmosphere_dim == 3)
    {
      // this is a loop over the different particle types
      for (Index pt_index = 0; pt_index < N_pt; ++ pt_index)
	{
	  // the first three loops are over the p_grid, lat_grid, lon_grid
	  for (Index p_index = cloudbox_limits[0];  
	       p_index < cloudbox_limits[1];
	       ++  p_index)
	    {
	      for (Index lat_index = cloudbox_limits[2];  
		   lat_index < cloudbox_limits[3]; ++  lat_index)
		{
		  for (Index lon_index = cloudbox_limits[4];  
		       lon_index < cloudbox_limits[5]; ++  lon_index)
		    {
		      
		      // these are loops over zenith angle and azimuth angle
		      for (Index za_index = 0; za_index < Nza; ++ za_index)
			{
			  for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
			    {
			      
			      // now the last two loops over the stokes 
			      //dimension.
			      for (Index i = 0;  i < N_i; ++  i)
				{
				  for (Index j = 0; j < N_i; ++ j)
				    {
				      
				      pha_mat_part(za_index, aa_index, 
						   i,j ) += 
					(pha_mat_spt(pt_index, za_index, 
						     aa_index, i, j) * 
					 pnd_field(pt_index, p_index,  
						   lat_index, lon_index));
				      
				    } 
				}
			    }
			}
		    }
		}
	    } 
	}
    }

}


//! Total extinction matrix
/*! 
  
  This method sums up the extinction  matrices for  particle 
  and gas
  \param ext_mat Output : Total extinction matrix  
  \param ext_mat_part Input : extinction matrix for particle
  \param ext_mat_part Input : extinction matrix for gas
*/

void ext_matCalc(Matrix& ext_mat,
		  const Matrix& ext_mat_part,
		  const Matrix& ext_mat_gas
		  )
{
  Index N_i = ext_mat_part.nrows(); 
  //CloneSize(ext_mat_gas, ext_mat_part){}
  ext_mat = ext_mat_part;
  ext_mat += ext_mat_gas;

}

//! Total Absorption Vector
/*! 
  
  This method sums up the absorption vectors for  particle 
  and gas
  \param abs_vec Output : Total absorption vector  
  \param abs_vec_part Input : absorption vector for particle
  \param abs_vec_gas Input : absorption vector for gas
*/

void abs_vecCalc(Vector& abs_vec,
		  const Vector& abs_vec_part,
		  const Vector& abs_vec_gas
		  )
{
  Index N_i = abs_vec_part.nelem(); 
  //CloneSize(abs_vec_gas, abs_vec_part){}
  abs_vec = abs_vec_part;
  abs_vec += abs_vec_gas;

}
