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
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "arts.h"
#include "array.h"
#include "check_input.h"
#include "matpackI.h"
#include "matpackVI.h"
#include "matpackVII.h"
#include "scatproperties.h"
#include "logic.h"
#include "auto_md.h"
#include "cloudbox.h"

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

  if (ext_mat_spt.nrows() != stokes_dim || 
      ext_mat_spt.ncols() != stokes_dim){
 
    throw runtime_error(" The dimension of the tensor ext_mat_spt should "
			"agree to stokes_dim");
  }
  cout << "The stokes dimension is :" << stokes_dim<<"\n";
  cout << "The scat_za_index : " << scat_za_index  << " \n " ;
  cout << "The scat_aa_index : " << scat_aa_index  << " \n " ;
  if (amp_mat.ncols() != 8)
    throw runtime_error(
			"Amplitude matrix must have 8 columns.");
  
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error(
			"The dimension of stokes vector can be "
			"only 1,2,3, or 4");
  }

  //find out frequency
  Numeric freq = f_grid[scat_f_index];

  cout << "The frequency: " << freq  << " \n " ;
  
  for (Index i = 0; i < npt; ++i)
    {
      ConstVectorView amp_coeffs = amp_mat(i,
				      scat_za_index,
				      scat_aa_index,
				      scat_za_index,
				      scat_aa_index,
				      Range(joker)
				      );
     
      cout << "The amplitude vector : " << " \n " 
	   << amp_coeffs << " \n " ;
      amp2ext(ext_mat_spt(i,Range(joker),Range(joker)),
	      amp_coeffs,
	      freq);
    }
  cout <<  "The Extinction Matrix for single particle type: " << " \n " 
       <<ext_mat_spt << " \n " ;
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
    throw runtime_error(" The dimension of the tensor pha_mat_spt should "
			"agree to stokes_dim");
  }
  if (amp_mat.ncols() != 8){
    throw runtime_error("Amplitude matrix must have 8 columns.");
  }
  if ((stokes_dim > 4) || (stokes_dim <1)){
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


//! this function calculates the absorption coefficient for a single particle 
//type from the extinction matrix and phase matrix.
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
  Index nza = pha_mat_spt.nbooks(); 
  Index naa = pha_mat_spt.npages(); 
  
  assert (is_size(scat_za_grid, nza));
  assert (is_size(scat_aa_grid, naa));
  
  if (abs_vec_spt.ncols() != stokes_dim ){
    //FIXME: Dimension should agree to stokes_dim !!!
    throw runtime_error("The dimension of  abs_vec_spt should"
			"agree to stokes_dim");
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
  
  if ((stokes_dim > 4) || (stokes_dim <1)){
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
  \param atmosphere_dim Input : he atmospheric dimensionality (now 1 or 3)
  \param scat_p_index Input : Pressure index for scattering calculations.
  \param scat_lat_index Input : Latitude index for scattering calculations.
  \param scat_lon_index Input : Longitude index for scattering calculations.
  
*/
void ext_mat_partCalc(Matrix& ext_mat_part,
		      const Tensor3& ext_mat_spt,
		      const Tensor4& pnd_field,
		      const Index& atmosphere_dim,
		      const Index& scat_p_index,
		      const Index& scat_lat_index,
		      const Index& scat_lon_index 
		     )
{
  Index N_pt = ext_mat_spt.npages();
  Index stokes_dim = ext_mat_spt.nrows();
  
  // (CE:) Size of ext_mat is not known before...
  ext_mat_part.resize(stokes_dim, stokes_dim);
  
  for (Index m = 0; m < stokes_dim; m++)
    {
      for (Index n = 0; n < stokes_dim; n++)
	ext_mat_part(m, n) = 0.0;// Initialisation
    }
  if (atmosphere_dim == 1)
    {
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt; l++)
	{ 
	  
	  // now the last two loops over the stokes dimension.
	  for (Index m = 0; m < stokes_dim; m++)
	    {
	      for (Index n = 0; n < stokes_dim; n++)
	      
	      ext_mat_part(m, n) += 
		(ext_mat_spt(l, m, n) * pnd_field(l, scat_p_index, 0, 0));
	    }
	}
    }
   cout <<  "The Extinction Matrix : " << " \n " 
       <<ext_mat_part << " \n " ;
  if (atmosphere_dim == 3)
    {
      
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt; l++)
	{ 
	  
	  // now the last two loops over the stokes dimension.
	  for (Index m = 0; m < stokes_dim; m++)
	    {
	      for (Index n = 0; n < stokes_dim; n++)
		
		ext_mat_part(m, n) +=  (ext_mat_spt(l, m, n) * 
					pnd_field(l, scat_p_index, 
						  scat_lat_index, 
						  scat_lon_index));
	      
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
  \param atmosphere_dim Input : he atmospheric dimensionality (now 1 or 3)
  \param scat_p_index Input : Pressure index for scattering calculations.
  \param scat_lat_index Input : Latitude index for scattering calculations.
  \param scat_lon_index Input : Longitude index for scattering calculations.
  
*/
void abs_vec_partCalc(Vector& abs_vec_part,
		      const Matrix& abs_vec_spt,
		      const Tensor4& pnd_field,
		      const Index& atmosphere_dim,
		      const Index& scat_p_index,
		      const Index& scat_lat_index,
		      const Index& scat_lon_index 
		      
		      )
{
  Index N_pt = abs_vec_spt.nrows();
  Index stokes_dim = abs_vec_spt.ncols();
  
  //(CE:) Resize abs_vec_part
  abs_vec_part.resize(stokes_dim);
  
  for (Index m = 0; m < stokes_dim; ++m)
    {
      
      abs_vec_part[m] = 0.0;// Initialisation
    }
  if (atmosphere_dim == 1)
    {
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt ; ++l)
	{
	  // now the loop over the stokes dimension.
          //(CE:) in the middle was l instead of m
	  for (Index m = 0; m < stokes_dim; ++m)
	    
	    abs_vec_part[m] += 
	      (abs_vec_spt(l, m) * pnd_field(l, scat_p_index, 0, 0));
	  
	}
    }
  
  
  if (atmosphere_dim == 3)
    {
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt ; ++l)
	{
	  
	  // now the loop over the stokes dimension.
	  for (Index m = 0; m < stokes_dim; ++m)
	    
	    abs_vec_part[m] += (abs_vec_spt(l, m) *
				pnd_field(l, scat_p_index,
					  scat_lat_index, 
					  scat_lon_index));
	  
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
  \param atmosphere_dim Input : he atmospheric dimensionality (now 1 or 3)
  \param scat_p_index Input : Pressure index for scattering calculations.
  \param scat_lat_index Input : Latitude index for scattering calculations.
  \param scat_lon_index Input : Longitude index for scattering calculations.
*/

void pha_mat_partCalc(Tensor4& pha_mat_part,
		      const Tensor5& pha_mat_spt,
		      const Tensor4& pnd_field,
		      const Index& atmosphere_dim,
		      const Index& scat_p_index,
		      const Index& scat_lat_index,
		      const Index& scat_lon_index 
		      
		     )
{

  Index N_pt = pha_mat_spt.nshelves();
  Index Nza = pha_mat_spt.nbooks();
  Index Naa = pha_mat_spt.npages();
  Index stokes_dim = pha_mat_spt.nrows();

  //(CE:) Resize pha_mat_part:
  pha_mat_part.resize(Nza, Naa, stokes_dim, stokes_dim);

  // Initialisation
  for (Index za_index = 0; za_index < Nza; ++ za_index)
    {
      for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
	{
	  for (Index stokes_index = 0; stokes_index < stokes_dim;
	       ++ stokes_index)
	    {
	      for (Index stokes_index = 0; stokes_index < stokes_dim; 
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
	  	  
	  // these are loops over zenith angle and azimuth angle
	  for (Index za_index = 0; za_index < Nza; ++ za_index)
	    {
	      for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
		{
		  
		  // now the last two loops over the stokes dimension.
		  for (Index stokes_index = 0; stokes_index < stokes_dim; 
		       ++  stokes_index)
		    {
		      for (Index stokes_index = 0; stokes_index < stokes_dim;
			   ++ stokes_index)
			
			pha_mat_part(za_index, aa_index,  
				     stokes_index, stokes_index) += 
			  
			  (pha_mat_spt(pt_index, za_index, aa_index,  
				       stokes_index, stokes_index) * 
			   pnd_field(pt_index,scat_p_index, 0, 0));
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
	  
	  // these are loops over zenith angle and azimuth angle
	  for (Index za_index = 0; za_index < Nza; ++ za_index)
	    {
	      for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
		{
		  
		  // now the last two loops over the stokes dimension.
		  for (Index i = 0;  i < stokes_dim; ++  i)
		    {
		      for (Index j = 0; j < stokes_dim; ++ j)
			{
			  
			  pha_mat_part(za_index, aa_index, i,j ) += 
			    (pha_mat_spt(pt_index, za_index, aa_index, i, j) * 
			     pnd_field(pt_index, scat_p_index,  
				       scat_lat_index, scat_lon_index));
			  
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
  Index stokes_dim = ext_mat_part.nrows(); 
  
  
  //(CE:) Define size of ext_mat:
  ext_mat.resize(stokes_dim, stokes_dim);

  ext_mat = ext_mat_part;
  ext_mat += ext_mat_gas;
  cout<<"Totaal extinction matrix"<<"\n"<< ext_mat<<"\n";
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
  Index stokes_dim = abs_vec_part.nelem(); 
   
  //(CE:) Resize abs_vec
  abs_vec.resize(stokes_dim);
  
  abs_vec = abs_vec_part;
  abs_vec += abs_vec_gas;

}
