/*!
  \file   m_optproperties.cc
  \author Sreerekha T.R. <rekha@uni-bremen.de>
  \date   Mon Jun 10 11:19:11 2002 
  \brief  This file contains workspace methods for calculating the optical 
  properties for the radiative transfer. These are extinction matrix,
  absorption vector and scattering vector. 
  The optical properties for the gases can be calculated with or without 
  Zeeman effect.
*/
/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include <iostream>
#include <stdlib.h>
#include <cmath>
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
#include "interpolation.h"
#include "make_vector.h"
#include "xml_io.h"
//! This method computes the extinction matrix for a single particle type
//  from teh amplitude matrix.
/*! 
\param ext_mat_spt Output and Input: extinction matrix for a single 
particle type.
\param amp_mat Input : amplitude matrix for each particle type
\param scat_za_index  Input : local zenith angle
\param scat_aa_index  Input : local azimuth angle
\param scat_f_index  Input : frequency index
\param f_grid  Input : frequency grid
\param stokes_dim  Input : stokes dimension
*/

void ext_mat_sptCalc(
		     Tensor3& ext_mat_spt,
		     const Tensor6& amp_mat,
		     const Index& scat_za_index,
		     const Index& scat_aa_index,
		     const Index& f_index,
		     const Vector& f_grid,
		     const Index& stokes_dim)
		     
{
  Index npt = ext_mat_spt.npages();
    
  if ((stokes_dim > 4) || (stokes_dim <1)){
    throw runtime_error("The dimension of stokes vector "
                        "can be only 1,2,3, or 4");
  }
  
  if (ext_mat_spt.nrows() != stokes_dim || 
      ext_mat_spt.ncols() != stokes_dim){
 
    throw runtime_error(" The dimension of the tensor ext_mat_spt should "
			"agree to stokes_dim");
  }
    
  if (amp_mat.ncols() != 8)
    throw runtime_error(
			"Amplitude matrix must have 8 columns.");
  
 
 
  //find out frequency
  Numeric freq = f_grid[f_index];

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

//! Calculate extinction matrix (spt) for the convergence test.  

/*
  This method computes the extinction matrix for a single particle type
  from the amplitude matrix. This extinction matrix contains only 
  extinction due to scattering, not the extinction due to absorption. 
  The function can be used only for testing.

\param ext_mat_spt Output and Input: extinction matrix for a
single particle type.
\param pha_mat_spt Input : phase matrix for single particle type
\param f_grid  Input : frequency grid
\param scat_za_grid  Input : zenith angle grid
\param scat_aa_grid  Input : azimuth angle grid
\param stokes_dim  Input : stokes dimension

*/
void ext_mat_sptScat(
		     Tensor3& ext_mat_spt,
                     const Tensor5& pha_mat_spt,
		     const Vector& scat_za_grid,
		     const Vector& scat_aa_grid,
		     const Index& stokes_dim)
		    
{

  Index npt = pha_mat_spt.nshelves();
  Index nza = pha_mat_spt.nbooks(); 
  Index naa = pha_mat_spt.npages(); 
  
  assert (is_size(scat_za_grid, nza));
  assert (is_size(scat_aa_grid, naa));

  if ((stokes_dim > 4) || (stokes_dim <1)){
    throw runtime_error("The dimension of stokes vector "
                        "can be only 1,2,3, or 4");
  }
  
  ext_mat_spt.resize(npt, stokes_dim, stokes_dim);


  if (pha_mat_spt.nrows() != stokes_dim || 
      pha_mat_spt.ncols() != stokes_dim ){
    throw runtime_error("The tensor pha_mat_spt should"
			"have 4 rows and 4 columns");
  }
  
 
  for (Index i = 0; i < npt; ++i)
    {
      ConstTensor4View pha=pha_mat_spt(i,
				       Range(joker),
				       Range(joker),
				       Range(joker),
				       Range(joker));
      
      // ConstVectorView za_grid = scat_za_grid[scat_za_index];
      //ConstVectorView aa_grid = scat_aa_grid[scat_aa_index];
      
      amp2ext_scat(ext_mat_spt(i,Range(joker), Range(joker)),
                   pha,
                   scat_za_grid,
                   scat_aa_grid);
     
    }	   

}



//! this function calculates phase matrix for a single particle type from 
//  the amplitude matrix
/*! 
  
 \param pha_mat_spt  Output and Input: phase matrix for a single particle type
 \param amp_mat  Input : amplitude matrix
 \param scat_za_index  Input : zenith angle index
 \param scat_aa_index  Input : azimuth angle index
 \param stokes_dim  Input : stokes dimension
  
*/

void pha_mat_sptCalc(
		     Tensor5& pha_mat_spt,
		     const Tensor6& amp_mat,
		     const Index& scat_za_index,
		     const Index& scat_aa_index,
		     const Index& stokes_dim)
		     
{
  Index npt = pha_mat_spt.nshelves();
 
   if ((stokes_dim > 4) || (stokes_dim <1)){
    throw runtime_error("The dimension of stokes vector" 
                         "can be only 1,2,3, or 4");
  }
 
  if (pha_mat_spt.nrows() != stokes_dim || 
      pha_mat_spt.ncols() != stokes_dim){
    throw runtime_error(" The dimension of the tensor pha_mat_spt should "
			"agree to stokes_dim");
  }
  if (amp_mat.ncols() != 8){
    throw runtime_error("Amplitude matrix must have 8 columns.");
  }
 
  
  for (Index i = 0; i < npt; ++i)
    {
      ConstTensor3View amp_coeffs = amp_mat(i,
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
  \param scat_za_grid Input : zenith angle grid.
  \param scat_aa_grid Input : azimuth angle grid.
  \param stokes_dim  Input : stokes dimension

*/
void abs_vec_sptCalc(
		     Matrix& abs_vec_spt,
		     const Tensor3& ext_mat_spt,
		     const Tensor5& pha_mat_spt,
		     const Vector& scat_za_grid,
		     const Vector& scat_aa_grid,
		     const Index& stokes_dim)
		     
{
  Index npt = abs_vec_spt.nrows();
  Index nza = pha_mat_spt.nbooks(); 
  Index naa = pha_mat_spt.npages(); 
  
  assert (is_size(scat_za_grid, nza));
  assert (is_size(scat_aa_grid, naa));

  if ((stokes_dim > 4) || (stokes_dim <1)){
    throw runtime_error("The dimension of stokes vector "
                        "can be only 1,2,3, or 4");
  }
  
  if (abs_vec_spt.ncols() != stokes_dim ){
    //Dimension should agree to stokes_dim !!!
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
  \param ext_mat Output and input : physical extinction coefficient 
  for the particles for given angles. 
  \param ext_mat_spt Input : extinction matrix for the single particle type
  \param pnd_field Input : particle number density givs the local 
  concentration for all particles.
  \param atmosphere_dim Input : he atmospheric dimensionality (now 1 or 3)
  \param scat_p_index Input : Pressure index for scattering calculations.
  \param scat_lat_index Input : Latitude index for scattering calculations.
  \param scat_lon_index Input : Longitude index for scattering calculations.
  \param stokes_dim  Input : stokes dimension
  
*/
void ext_matAddPart(
		      Tensor3& ext_mat,
                      const Tensor3& ext_mat_spt,
		      const Tensor4& pnd_field,
		      const Index& atmosphere_dim,
		      const Index& scat_p_index,
		      const Index& scat_lat_index,
		      const Index& scat_lon_index,
		      const Index& stokes_dim
		      ) 
		     
{
  Index N_pt = ext_mat_spt.npages();
 
  Matrix ext_mat_part(stokes_dim, stokes_dim, 0.0);

  
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error(
			"The dimension of stokes vector can be "
			"only 1,2,3, or 4");
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
	       //summation of the product of pnd_field and 
		//ext_mat_spt.
	      ext_mat_part(m, n) += 
		(ext_mat_spt(l, m, n) * pnd_field(l, scat_p_index, 0, 0));
	    }
	}

      //Add particle extinction matrix to *ext_mat*.
      ext_mat(0, Range(joker), Range(joker)) += ext_mat_part;
    }
 
  if (atmosphere_dim == 3)
    {
      
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt; l++)
	{ 
	  
	  // now the last two loops over the stokes dimension.
	  for (Index m = 0; m < stokes_dim; m++)
	    {
	      for (Index n = 0; n < stokes_dim; n++)
		 //summation of the product of pnd_field and 
		//ext_mat_spt.
		ext_mat_part(m, n) +=  (ext_mat_spt(l, m, n) * 
					pnd_field(l, scat_p_index, 
						  scat_lat_index, 
						  scat_lon_index));
	      
	    } 
	}

      //Add particle extinction matrix to *ext_mat*.
      ext_mat += ext_mat_part;

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
  \param stokes_dim  Input : stokes dimension
*/
void abs_vecAddPart(
		      Matrix& abs_vec,
                      const Matrix& abs_vec_spt,
		      const Tensor4& pnd_field,
		      const Index& atmosphere_dim,
		      const Index& scat_p_index,
		      const Index& scat_lat_index,
		      const Index& scat_lon_index,
                      const Index& stokes_dim
                      ) 
		    
{
  Index N_pt = abs_vec_spt.nrows();
  
  Vector abs_vec_part(stokes_dim, 0.0);

  if ((stokes_dim > 4) || (stokes_dim <1)){
    throw runtime_error("The dimension of stokes vector "
                        "can be only 1,2,3, or 4");
  } 
 
  if (atmosphere_dim == 1)
    {
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt ; ++l)
	{
	  // now the loop over the stokes dimension.
          //(CE:) in the middle was l instead of m
	  for (Index m = 0; m < stokes_dim; ++m)
	     //summation of the product of pnd_field and 
	    //abs_vec_spt.
	    abs_vec_part[m] += 
	      (abs_vec_spt(l, m) * pnd_field(l, scat_p_index, 0, 0));
	  
	}
      //Add the particle absorption
      abs_vec(0, Range(joker)) += abs_vec_part;
    }
  
  if (atmosphere_dim == 3)
    {
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt ; ++l)
	{
	  
	  // now the loop over the stokes dimension.
	  for (Index m = 0; m < stokes_dim; ++m)
	     //summation of the product of pnd_field and 
	    //abs_vec_spt.
	    abs_vec_part[m] += (abs_vec_spt(l, m) *
				pnd_field(l, scat_p_index,
					  scat_lat_index, 
					  scat_lon_index));
	  
	}
      //Add the particle absorption
      abs_vec(0,Range(joker)) += abs_vec_part;
    }
} 

// //! Method for creating an absorption vector for gases for test purposes
// /*! 
  
// The aim of this method is just to give a reasonable value for absorption 
// coefficients that can be used to test the radiative transfer calculation.
// This method takes pressure grids and absorption values corresponding to them
// from ARTS.1.0 and interpolates it onto the pressure grid under consideration 
// for radiative transfer calculation.

// \param abs_vec_gas Output : The gaseous absorption vector
// \param p_grid Input : pressure grid on whiich the radiative transfer 
// calculations are done
// \param atmosphere_dim Input : atmospheric dimension(here  considered only 1D)
// \param stokes_dim Input : stokes dimension
// \param scat_p_index Input : The index corresponding to the pressure grid
// */
// void abs_vec_gasExample(Vector& abs_vec,
// 			const Vector& p_grid,
// 			const Index& atmosphere_dim,
// 			const Index& stokes_dim,
// 			const Index& scat_p_index) 
// {
//   abs_vec.resize( stokes_dim );
//   Vector abs_gas( p_grid.nelem() );
//   if (atmosphere_dim == 1){

//     //         MakeVector typical_abs(0,0,0,0,0,0,0,0,0,0);

//   //This is a typical absorption calculated from arts.1.0 (tropical)
  
//      // MakeVector typical_abs(0.0218978044006849,
// //                             0.0114605893154005,
// //                             0.00355854032448724,
// //                             0.00161793243125695,
// //                             0.000612486464168897,
// //                             0.000168382126027889,
// //                             3.35522552862195e-05,
// //                             8.52950996338938e-06,
// //                             4.32022347140681e-06,
// //                             3.12994113724593e-06);


//  //This is a typical absorption calculated from arts.1.0 (mls)
  
//     MakeVector typical_abs(0.0164163638663716,
// 			   0.00768271579907592,
// 			   0.00294668635075111,
// 			   0.00125411825404778,
// 			   0.000570445848162073,
// 			   0.000236462958473072,
// 			   4.40975932116215e-05,
// 			   7.31218846316807e-06,
// 			   3.643089167928e-06,
// 			   3.12497184475723e-06);
    
// //This is a typical absorption calculated from arts.1.0 (mlw)
  
//     /*MakeVector typical_abs(0.00433659795310861,
// 			   0.00267195350007522,
// 			   0.00122676231151285,
// 			   0.000524902357458338,
// 			   0.000142742001685155,
// 			   4.65799636740885e-05,
// 			   1.07249251986862e-05,
// 			   5.85896298848906e-06,
// 			   4.94996692358466e-06,
// 			   4.35302994835574e-06);*/


//     //The pressure grid for the above calculation (mlw)
//     MakeVector typical_abs_pgrid(100000,
// 				 77426.3682681127,
// 				 59948.4250318941,
// 				 46415.8883361278,
// 				 35938.1366380463,
// 				 27825.5940220712,
// 				 21544.3469003188,
// 				 16681.0053720006,
// 				 12915.4966501488,
// 				 10000);
//     //p_grid is the new pressure grid
//     ArrayOfGridPos p_gp(p_grid.nelem());

//     gridpos(p_gp, typical_abs_pgrid, p_grid);

//     Matrix itw(p_grid.nelem(),2);

//     interpweights(itw, p_gp);
//     //interpolating absorption coefficients from typical_abs_pgrid to p_grid
//     interp(abs_gas, itw, typical_abs, p_gp);
      
//     for (Index i = 0; i<stokes_dim; ++i)
//       {
// 	if (i == 0){
// 	  abs_vec[i] = abs_gas[scat_p_index];
// 	}
// 	else{
// 	  abs_vec[i] = 0.0;
// 	}
//       }
//   }
  
// }


  
//! Method for creating the extinction matrix.
/*! 

 This method is also for test purposes and is very simple.  It takes 
absorption coefficients from the method abs_vec_gasExample and put them 
along the diagonal of a 4 X 4 matrix (if stokes_dim = 4) and the 
off-diagonal elements are set to zero.

\param ext_mat Output : The extinction matrix corresponding to gaseous 
species
\param abs_scalar_gas Input : absorption from gaseous species
\param f_index : Input: Frequency index. If all frequencies shall be calculated this variable has to be set to 0. 
*/

void ext_matAddGas(Tensor3& ext_mat,
                   const Vector& abs_scalar_gas,
                   const Index& f_index,
                   const Index& atmosphere_dim)
{
  Index stokes_dim = ext_mat.ncols();
  Matrix ext_mat_gas(stokes_dim, stokes_dim,0.0);
  cout<<"stokes_dim"<<stokes_dim<<endl;
  //FIXME: After the scalar ags function is ready.
  if (atmosphere_dim == 1){
  
    for (Index i =0; i < stokes_dim; ++i)
      {
	for (Index j =0; j < stokes_dim; ++j)
	  {
	    
	    if ( i == j){
	      // Dies ist noch Quatsch!!!
	      ext_mat_gas(i,j) = abs_scalar_gas[f_index];
	    }
	    else{
	      ext_mat_gas(i,j) = 0.0;
	    }
	  }
      }
  
  ext_mat(0, Range(joker), Range(joker)) += ext_mat_gas;
  }
 }

//! Method for creating the absorption vector.
/*! 

to be written ... 

\param abs_vec Output : The extinction matrix corresponding to gaseous 
species
\param abs_scalar_gas Input : absorption from gaseous species
\param f_index : Input: Frequency index. If all frequencies shall be calculated this variable has to be set to 0. 
*/

void abs_vecAddGas(Matrix& abs_vec,
                   const Vector& abs_scalar_gas,
                   const Index& f_index,
                   const Index& atmosphere_dim)
{
  Index stokes_dim = abs_vec.ncols();
  Vector abs_vec_gas(stokes_dim);
  

  if (atmosphere_dim == 1){
    
    abs_vec_gas[0] = abs_scalar_gas[f_index];
  
    for (Index i = 1; i < stokes_dim; ++i)
      {
        abs_vec_gas[i] = 0.0;
	  
      }
    
    abs_vec(0, Range(joker)) += abs_vec_gas;
  }
  
 }


//! Phase Matrix for the particle 
/*! 
  
  This function sums up the phase matrices for all particle 
  types weighted with particle number density
  \param pha_mat Output : physical phase matrix 
  for the particles for given angles. 
  \param pha_mat_spt Input : phase matrix for the single particle type
  \param pnd_field Input : particle number density givs the local 
  concentration for all particles.
  \param atmosphere_dim Input : he atmospheric dimensionality (now 1 or 3)
  \param scat_p_index Input : Pressure index for scattering calculations.
  \param scat_lat_index Input : Latitude index for scattering calculations.
  \param scat_lon_index Input : Longitude index for scattering calculations.
*/

void pha_matCalc(
		 Tensor4& pha_mat,
		 const Tensor5& pha_mat_spt,
		 const Tensor4& pnd_field,
		 const Index& atmosphere_dim,
		 const Index& scat_p_index,
		 const Index& scat_lat_index,
		 const Index& scat_lon_index) 
		      
{

  Index N_pt = pha_mat_spt.nshelves();
  Index Nza = pha_mat_spt.nbooks();
  Index Naa = pha_mat_spt.npages();
  Index stokes_dim = pha_mat_spt.nrows();
 
  pha_mat.resize(Nza, Naa, stokes_dim, stokes_dim);

  // Initialisation
  for (Index za_index = 0; za_index < Nza; ++ za_index)
    {
      for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
	{
	  for (Index stokes_index_1 = 0; stokes_index_1 < stokes_dim;
	       ++ stokes_index_1)
	    {
	      for (Index stokes_index_2 = 0; stokes_index_2 < stokes_dim; 
		   ++ stokes_index_2)
		pha_mat(za_index, aa_index, stokes_index_1, stokes_index_2)

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
		  for (Index stokes_index_1 = 0; stokes_index_1 < stokes_dim; 
		       ++  stokes_index_1)
		    {
		      for (Index stokes_index_2 = 0; stokes_index_2 < stokes_dim;
			   ++ stokes_index_2)
			 //summation of the product of pnd_field and 
			  //pha_mat_spt.
			pha_mat(za_index, aa_index,  
				     stokes_index_1, stokes_index_2) += 
			  
			  (pha_mat_spt(pt_index, za_index, aa_index,  
				       stokes_index_1, stokes_index_2) * 
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
			  //summation of the product of pnd_field and 
			  //pha_mat_spt.
			  pha_mat(za_index, aa_index, i,j ) += 
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





// //! Extinction coefficient matrix  for convergence test.
// /*! 
//   This function sums up the convergence extinction matrices for all particle 
//   types weighted with particle number density. This is exactly the same 
//   function as *ext_mat_partCalc*, only with different input, *ext_mat_conv_spt*
//   instead of ext_mat_spt.

//   \param ext_mat_part Output : physical extinction coefficient 
//   for the particles for given angles. 
//   \param ext_mat_conv_spt Input : extinction matrix for the single 
//   particle type
//   \param pnd_field Input : particle number density givs the local 
//   concentration for all particles.
//   \param atmosphere_dim Input : he atmospheric dimensionality (now 1 or 3)
//   \param scat_p_index Input : Pressure index for scattering calculations.
//   \param scat_lat_index Input : Latitude index for scattering calculations.
//   \param scat_lon_index Input : Longitude index for scattering calculations.
// */
// void ext_mat_partScat(
// 		      Matrix& ext_mat_part,
// 		      const Tensor3& ext_mat_spt,
// 		      const Tensor4& pnd_field,
// 		      const Index& atmosphere_dim,
// 		      const Index& scat_p_index,
// 		      const Index& scat_lat_index,
// 		      const Index& scat_lon_index) 
		     
// {
//   Index N_pt = ext_mat_spt.npages();
//   Index stokes_dim = ext_mat_spt.nrows();
  
//   ext_mat_part.resize(stokes_dim, stokes_dim);
  
//   for (Index m = 0; m < stokes_dim; m++)
//     {
//       for (Index n = 0; n < stokes_dim; n++)
// 	ext_mat_part(m, n) = 0.0;// Initialisation
//     }
 
//   if (atmosphere_dim == 1)
//     {
//       // this is a loop over the different particle types
//       for (Index l = 0; l < N_pt; l++)
// 	{ 
	  
// 	  // now the last two loops over the stokes dimension.
// 	  for (Index m = 0; m < stokes_dim; m++)
// 	    {
// 	      for (Index n = 0; n < stokes_dim; n++)
// 	       //summation of the product of pnd_field and 
// 		//ext_mat_conv_spt.
// 	      ext_mat_part(m, n) += 
// 		(ext_mat_spt(l, m, n) * pnd_field(l, scat_p_index, 0, 0));
// 	    }
// 	}
//     }
//   if (atmosphere_dim == 3)
//     {
      
//       // this is a loop over the different particle types
//       for (Index l = 0; l < N_pt; l++)
// 	{ 
	  
// 	  // now the last two loops over the stokes dimension.
// 	  for (Index m = 0; m < stokes_dim; m++)
// 	    {
// 	      for (Index n = 0; n < stokes_dim; n++)
// 		 //summation of the product of pnd_field and 
// 		//ext_mat_conv_spt.
// 		ext_mat_part(m, n) +=  (ext_mat_spt(l, m, n) * 
// 					pnd_field(l, scat_p_index, 
// 						  scat_lat_index, 
// 						  scat_lon_index));
	      
// 	    } 
// 	}
//     }
// } 


//! Method for getting amp_mat from amp_mat_raw.
/*! 
  This method converts the raw amplitude matrix data namely amp_mat_raw 
  to workspace variable amp_mat which can be directly used for calculating  
  scattering properties of particles.  The data type of amp_mat_raw is 
  an ArrayOfArrayOfTensor6 which contains one gridded field for each 
  particle type the user wants.  One data file
  contains amp_mat_raw for one particle type.  One data file has
  amp_mat_raw for one particle for a set of frequencies, incoming and 
  outgoing zenith angles and azimuth angles. The frequencies, angles, 
  and the amplitude matrix are each a Tensor 6.  The size of amp_mat_raw 
  is amp_mat_raw[Npt][7].  where Npt is the number of particle types.
  amp_mat_raw[Npt][0] gives the frequency tensor [N_f, 1, 1, 1, 1, 1] 
  where N_f gives the number of frequencies considered in that particular
  database file. Similarly,
  amp_mat_raw[Npt][1] gives the outgoing zenith angle tensor [1,Nza,1,1,1,1],
  amp_mat_raw[Npt][2] gives the outgoing azimuth angle tensor [1,1,Naa,1,1,1], 
  amp_mat_raw[Npt][3] gives the incoming zentih angle tensor [1,1,1,Nza,1,1], 
  amp_mat_raw[Npt][4] gives the incoming azimuth angle tensor [1,1,1,1,Naa,1],
  amp_mat_raw[Npt][5] is a dummy tensor6 and 
  amp_mat_raw[Npt][6] gives amplitude matrix which is also a tensor 6 of
  size [N_f, N_za, N_aa, N_za, N_aa, 8]. Here, Nza is the number 
  of zenith angles, Naa is the number of azimuth angles and 8 denotes the 
  amplitude matrix elements.  
  
  In this method, we have to interpolate the raw data calculated on a 
  specific angular and frequency grid onto a grid which is set by the user. 
  Since we decide that frequency should be the outermost loop for our 
  radiative transfer calculation the frequency grid contains just one 
  value specified by the index f_index.  The angles for which the 
  calculations are to be done are specified by scat_za_grid and scat_aa_grid.
  
  The output of this method is amp_mat has to be a 
  Tensor6 with the first dimension being that of the particle type, then the
  angles and finally the amplitude matrix element 8. The size of amp_mat is 
  (Npt, Nza, Naa, Nza, Naa, 8).  Note that the dimension frequency is taken
  out.
  
  \param amp_mat Output : Amplitude matrix which is interpolated on the 
  required grids set by the user.
  \param amp_mat_raw Input : The original data as read from the database.
  \param f_index Input : The frequency index.
  \param f_grid Input : The frequency grid as required for RT calculation
  \param scat_za_grid Input : The zenith angle grid for which scattering 
  properties are to be calculated 
  \param scat_aa_grid Input : The azimuth angle grid for which scattering 
  properties are to be calculated 
*/
void amp_matCalc(Tensor6& amp_mat,
		 const ArrayOfArrayOfTensor6& amp_mat_raw,
		 const Index& f_index,
		 const Vector& f_grid,
		 const Vector& scat_za_grid,
		 const Vector& scat_aa_grid)
{
  Index N_pt = amp_mat_raw.nelem();
  
  Index N_za = scat_za_grid.nelem();
  Index N_aa = scat_aa_grid.nelem();
  Index N_i = amp_mat_raw [ 0 ] [ 6 ].ncols();

  if (N_i != 8)
    throw runtime_error(
			"Amplitude matrix must have 8 columns.");

  amp_mat.resize(N_pt, N_za, N_aa, N_za, N_aa, N_i);  
  /*amp_mat_raw is an ArrayOfArrayOfTensor6 from which we can get the 
    amplitude matrix which is a tensor 6 with size(Nf, Nza, Naa, Nza,Naa,8). 
    This is what is directly read from the database.  As mentioned in 
    the beginning in this method we have to interpolate amplitude matrix 
    from the grid on which it is calculated in the database to the grid set
    by the user. The interpolation has to be done for the frequency grid,
    zenith angle grid and azimuth angle grid. Since this interpolation is
    from a gridded field to a new field, it is called a green interpolation.
    For more insight into the interpolation schemes refer to Chapter 8-
    Interpolation  of AUG.
    The new grids are f_grid(only one element in this vector), scat_za_grid, 
    scat_aa_grid and the old grids are deroved from the ArrayOfArrayOfTensor6 
    amp_amt_raw.*/

  //Loop over the particle types. We can get information about the number of
  //particle types from the input ArrayOfArrayOfTensor6 amp_mat_raw
   for (Index ipt = 0; ipt < N_pt; ++ ipt )
     {
       Index N_i = amp_mat_raw [ ipt ] [ 6 ].ncols();
       //calling the interpolation routines.  
       
       //Define the grid position arrays. 

       // for  frequency : 
       ArrayOfGridPos f_gp(1); 

       // for outgoing zenith angle grids : 
       ArrayOfGridPos za_gp(scat_za_grid.nelem());
       
       // for outgoing azimuth angle grids : 
       ArrayOfGridPos aa_gp(scat_aa_grid.nelem());

       // for incoming zenith angle grids : 
       ArrayOfGridPos za_in_gp(scat_za_grid.nelem());

       // for incoming azimuth angle grids : rows
       ArrayOfGridPos aa_in_gp(scat_aa_grid.nelem());
       
       // Set up Grid position arrays by calling the function gridpos.
       
       /*for frequency :

       f_gp is the ArrayOfGridpos.
       
       original frequency grid as in the data base can be got from the 
       ArrayOfTensor6 amp_mat_raw.  
       amp_mat_raw[ipt][0] ( Range(joker), 0, 0, 0, 0, 0).  
       
       f_grid(Range(f_index),1) is the new grid which is also a vector
       formally but with one element given by the f_index*/
       //cout<<f_index<<"\n";
       gridpos (f_gp,
		amp_mat_raw [ ipt ]  [ 0 ] ( Range(joker), 0, 0, 0, 0, 0),
		f_grid[Range(f_index, 1)]);
       
       //like for frquency we can get the gridpostions for the angular grids.

       gridpos (za_gp, 
		amp_mat_raw [ ipt ] [  1 ] ( 0, Range(joker), 0, 0, 0, 0),
		scat_za_grid);
      
       gridpos (aa_gp,
		amp_mat_raw [ ipt ] [ 2 ] ( 0, 0, Range(joker), 0, 0, 0),
	       scat_aa_grid);

       gridpos (za_in_gp, 
		amp_mat_raw [ ipt ] [ 3 ] ( 0, 0, 0, Range(joker), 0, 0),
	       scat_za_grid);

       gridpos (aa_in_gp,
		amp_mat_raw [ ipt ] [ 4] ( 0, 0, 0, 0, Range(joker), 0),
		scat_aa_grid);
       
       /*The interpolation weights. Since here we interpolate 
	 simultaneously in 5 dimensions we require exactly 32 
	 interpolation weights. (There are 2^n interpolation 
	 weights for an n-dimensional interpolation. ).  Since
	 in ARTS we want to save a lot by re-using the weights
	 the interpolation weights are stored in a tensor whcih
	 has one more dimension than the output field. The last 
	 dimension is for weight, this explains the last dimension
	 of interpolation weight to be 32 in this case.
	 */
       
       Tensor6 itw(1, scat_za_grid.nelem(), scat_aa_grid.nelem(),
		   scat_za_grid.nelem(), scat_aa_grid.nelem(), 32);

       //function for computing interpolation weight tensor. For this
       //step also we need only gridpositions, not the field yet.
       interpweights ( itw, f_gp, za_gp, aa_gp, za_in_gp, aa_in_gp );
      
       // This is a green interpolation for all columns of source_amp.
       // Loop over the column
       for (Index i = 0 ; i < N_i ; ++ i)
	 {
	   /*Select the current column of target.  The first element 
	     of target is the running variable ipt, whcih gives the
	     particle type.*/
	   Tensor5View target_amp = amp_mat(Range(ipt, 1),
					    Range(joker),
					    Range(joker),
					    Range(joker),
					    Range(joker),
					    i);
	   
	   //Select the current column of source.  Here the first
	   //element is frequency
	   ConstTensor5View source_amp = 
	     amp_mat_raw [ ipt ][6 + 7* ipt](Range(joker),
				      Range(joker),
				      Range(joker),
				      Range(joker),
				      Range(joker),
				      i);

	   //Interpolation is done :
	   interp (target_amp, 
		   itw, 
		   source_amp, 
		   f_gp, 
		   za_gp,
		   aa_gp,
		   za_in_gp,
		   aa_in_gp);
	   
	 }//close column index
     }//close particle index loop
   //cout<< amp_mat<<"\n";
}


