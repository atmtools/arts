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
\param ext_mat_spt Output and Input: extinction matrix for a single particle type.
\param amp_mat Input : amplitude matrix for each particle type
\param f_grid  Input : frequency grid
\param scat_f_index  Input : frequency index
\param scat_za_index  Input : local zenith angle
\param scat_aa_index  Input : local azimuth angle

*/

void ext_mat_sptCalc(
		     Tensor3& ext_mat_spt,
		     const Tensor6& amp_mat,
		     const Index& scat_za_index,
		     const Index& scat_aa_index,
		     const Index& scat_f_index,
		     const Vector& f_grid)
		     
{
  Index npt = ext_mat_spt.npages();
  Index stokes_dim = ext_mat_spt.nrows();

  if (ext_mat_spt.nrows() != stokes_dim || 
      ext_mat_spt.ncols() != stokes_dim){
 
    throw runtime_error(" The dimension of the tensor ext_mat_spt should "
			"agree to stokes_dim");
  }
  //cout << "The stokes dimension is :" << stokes_dim<<"\n";
  //cout << "The scat_za_index : " << scat_za_index  << " \n " ;
  //cout << "The scat_aa_index : " << scat_aa_index  << " \n " ;
  // cout << "Number of particle type : " << npt  << " \n " ;
  
  
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

  //cout << "The frequency: " << freq  << " \n " ;
  
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

void pha_mat_sptCalc(
		     Tensor5& pha_mat_spt,
		     const Tensor6& amp_mat,
		     const Index& scat_za_index,
		     const Index& scat_aa_index)
		     
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
  \param scat_za_grid Input : zenith angle grid.
  \param scat_aa_grid Input : azimuth angle grid.
 
*/
void abs_vec_sptCalc(
		     Matrix& abs_vec_spt,
		     const Tensor3& ext_mat_spt,
		     const Tensor5& pha_mat_spt,
		     const Vector& scat_za_grid,
		     const Vector& scat_aa_grid)
		     
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
void ext_mat_partCalc(
		      Matrix& ext_mat_part,
		      const Tensor3& ext_mat_spt,
		      const Tensor4& pnd_field,
		      const Index& atmosphere_dim,
		      const Index& scat_p_index,
		      const Index& scat_lat_index,
		      const Index& scat_lon_index) 
		     
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
	       //summation of the product of pnd_field and 
		//ext_mat_spt.
	      ext_mat_part(m, n) += 
		(ext_mat_spt(l, m, n) * pnd_field(l, scat_p_index, 0, 0));
	    }
	}
    }
  cout <<  "The Extinction Matrix for particle: " << " \n " 
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
		 //summation of the product of pnd_field and 
		//ext_mat_spt.
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
void abs_vec_partCalc(
		      Vector& abs_vec_part,
		      const Matrix& abs_vec_spt,
		      const Tensor4& pnd_field,
		      const Index& atmosphere_dim,
		      const Index& scat_p_index,
		      const Index& scat_lat_index,
		      const Index& scat_lon_index) 
		    
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
	     //summation of the product of pnd_field and 
	    //abs_vec_spt.
	    abs_vec_part[m] += 
	      (abs_vec_spt(l, m) * pnd_field(l, scat_p_index, 0, 0));
	  
	}
    }
  
  cout <<  "The Absorption Vector for particle : " << " \n " 
     <<abs_vec_part << " \n " ;
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
      
    }
} 

//! Method for creating an absorption vector for gases for test purposes
/*! 
  
The aim of this method is just to give a reasonable value for absorption 
coefficients that can be used to test the radiative transfer calculation.
This method takes pressure grids and absorption values corresponding to them
from ARTS.1.0 and interpolates it onto the pressure grid under consideration 
for radiative transfer calculation.

\param abs_vec_gas Output : The gaseous absorption vector
\param p_grid Input : pressure grid on whiich the radiative transfer 
calculations are done
\param atmosphere_dim Input : atmospheric dimension(here  considered only 1D)
\param stokes_dim Input : stokes dimension
\param scat_p_index Input : The index corresponding to the pressure grid
*/
void abs_vec_gasExample(Vector& abs_vec_gas,
			const Vector& p_grid,
			const Index& atmosphere_dim,
			const Index& stokes_dim,
			const Index& scat_p_index) 
{
  abs_vec_gas.resize( stokes_dim );
  Vector abs_gas( p_grid.nelem() );
  if (atmosphere_dim == 1){
    //This is a typical absorption calculated from arts.1.0
    MakeVector typical_abs(0.016414284648894,
			   0.00065204114511011,
			   0.000156049846860233,
			   4.54320063675961e-05,
			   1.52191594739311e-05,
			   5.13136166733503e-06,
			   1.37451108959307e-06,
			   6.70848900165098e-07,
			   4.06285725309355e-07,
			   2.57499613700983e-07);
    //The pressure grid for the above calculation
    MakeVector typical_abs_pgrid(100000,
				 77426.3682681127,
				 59948.4250318941,
				 46415.8883361278,
				 35938.1366380463,
				 27825.5940220712,
				 21544.3469003188,
				 16681.0053720006,
				 12915.4966501488,
				 10000);
    //p_grid is the new pressure grid
    ArrayOfGridPos p_gp(p_grid.nelem());

    gridpos(p_gp, typical_abs_pgrid, p_grid);

    Matrix itw(p_grid.nelem(),2);

    interpweights(itw, p_gp);
    //interpolating absorption coefficients from typical_abs_pgrid to p_grid
    interp(abs_gas, itw, typical_abs, p_gp);
      
    for (Index i = 0; i<stokes_dim; ++i)
      {
	if (i == 0){
	  abs_vec_gas[i] = abs_gas[scat_p_index];
	}
	else{
	  abs_vec_gas[i] = 0.0;
	}
      }
  }
  
}  
//! Method for creating an extinction matrix for gases for test purposes
/*! 

 This method is also for test purposes and is very simple.  It takes 
absorption coefficients from the method abs_vec_gasExample and put them 
along the diagonal of a 4 X 4 matrix (if stokes_dim = 4) and the 
off-diagonal elements are set to zero.

\param ext_mat_gas Output : The extinction matrix corresponding to gaseous 
species
\param abs_vec_gas Input : absorption from gaseous species
\param atmosphere_dim Input : atmospheric dimension
\param stokes_dim Input : stokes dimension
\param scat_p_index Input : the pressure index corresponding to the radiative 
transfer calculation
*/
void ext_mat_gasExample(Matrix& ext_mat_gas,
			const Vector& abs_vec_gas,
			const Index& atmosphere_dim,
			const Index& stokes_dim,
			const Index& scat_p_index)
{
  ext_mat_gas.resize(stokes_dim, stokes_dim);
  if (atmosphere_dim == 1){
  
    for (Index i =0; i < stokes_dim; ++i)
      {
	for (Index j =0; j < stokes_dim; ++j)
	  {
	    
	    if ( i == j){
	      
	      ext_mat_gas(i,j) = abs_vec_gas[i];
	    }
	    else{
	      ext_mat_gas(i,j) = 0.0;
	    }
	  }
      }
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

  //(CE:) Resize pha_mat:
  pha_mat.resize(Nza, Naa, stokes_dim, stokes_dim);

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
		pha_mat(za_index, aa_index, stokes_index, stokes_index)

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
			 //summation of the product of pnd_field and 
			  //pha_mat_spt.
			pha_mat(za_index, aa_index,  
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


//! Total extinction matrix
/*! 
  
  This method sums up the extinction  matrices for  particle 
  and gas
  \param ext_mat Output : Total extinction matrix  
  \param ext_mat_part Input : extinction matrix for particle
  \param ext_mat_gas Input : extinction matrix for gas
*/

void ext_matCalc(
		 Matrix& ext_mat,
		 const Matrix& ext_mat_part,
		 const Matrix& ext_mat_gas)
		 
{
  Index stokes_dim = ext_mat_part.nrows(); 
  
  
  //(CE:) Define size of ext_mat:
  ext_mat.resize(stokes_dim, stokes_dim);
  //Addition of the two matrices, ext_mat_part and ext_mat_gas.
  ext_mat = ext_mat_part;
  ext_mat += ext_mat_gas;
  // cout<<"Totaal extinction matrix"<<"\n"<< ext_mat<<"\n";
}

//! Total Absorption Vector
/*! 
  
  This method sums up the absorption vectors for  particle 
  and gas
  \param abs_vec Output : Total absorption vector  
  \param abs_vec_part Input : absorption vector for particle
  \param abs_vec_gas Input : absorption vector for gas
*/

void abs_vecCalc(
		 Vector& abs_vec,
		 const Vector& abs_vec_part,
		 const Vector& abs_vec_gas)
  
{
  Index stokes_dim = abs_vec_part.nelem(); 
   
  //(CE:) Resize abs_vec
  abs_vec.resize(stokes_dim);
  // addition of abs_vec_part and abs_vec_gas
  abs_vec = abs_vec_part;
  abs_vec += abs_vec_gas;

}

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
  value specified by the index scat_f_index.  The angles for which the 
  calculations are to be done are specified by scat_za_grid and scat_aa_grid.
  
  The output of this method is amp_mat has to be a 
  Tensor6 with the first dimension being that of the particle type, then the
  angles and finally the amplitude matrix element 8. The size of amp_mat is 
  (Npt, Nza, Naa, Nza, Naa, 8).  Note that the dimension frequency is taken
  out.
  
  \param amp_mat Output : Amplitude matrix which is interpolated on the 
  required grids set by the user.
  \param amp_mat_raw Input : The original data as read from the database.
  \param scat_f_index Input : The frequency index.
  \param f_grid Input : The frequency grid as required for RT calculation
  \param scat_za_grid Input : The zenith angle grid for which scattering 
  properties are to be calculated 
  \param scat_aa_grid Input : The azimuth angle grid for which scattering 
  properties are to be calculated 
*/
void amp_matCalc(Tensor6& amp_mat,
		 const ArrayOfArrayOfTensor6& amp_mat_raw,
		 const Index& scat_f_index,
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
       
       f_grid(Range(scat_f_index),1) is the new grid which is also a vector
       formally but with one element given by the scat_f_index*/
       //cout<<scat_f_index<<"\n";
       gridpos (f_gp,
		amp_mat_raw [ ipt ]  [ 0 ] ( Range(joker), 0, 0, 0, 0, 0),
		f_grid[Range(scat_f_index, 1)]);
       
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

//! This method interpolates clear sky field on the cloudbox boundary 
//on all grid points inside the cloud box. 

/*! 
  This method uses a linear 3D interpolation scheme to obtain the 
  radiation field on all grid points inside the cloud box form the 
  clear sky field on the cloud bod boundary.
  
  \param i_field Output : Intensity field
  \param scat_i_p Input : Intensity field on cloudbox boundary 
  (equal pressure surfaces)
  \param scat_i_lat Input : Intensity field on cloudbox boundary 
  (equal latitude surfaces)
  \param scat_i_lon Input : Intensity field on cloudbox boundary
  (equal longitude surfaces)
  \param f_grid Input : frequency grid
  \param scat_f_index Input : the frequency index for scattering calculation
  \param p_grid Input : the pressure grid
  \param lat_grid Input : the latitude grid
  \param lon_grid Input : the longitude grid
  \param cloudbox_limits Input : Limits of the cloud box
  \param atmospere_dim Input : dimension of atmosphere
*/
void i_fieldSet(Tensor6& i_field,
		const Tensor7& scat_i_p,
		const Tensor7& scat_i_lat,
		const Tensor7& scat_i_lon,
		const Vector& f_grid,
		const Index& scat_f_index,
		const Vector& p_grid,
		const Vector& lat_grid,
		const Vector& lon_grid,
		const ArrayOfIndex& cloudbox_limits,
		const Index& atmosphere_dim )
{
  if(atmosphere_dim == 1)
    {
      Index  N_f = scat_i_p.nlibraries();
      if (f_grid.nelem() != N_f){
	
	throw runtime_error(" scat_i_p should have same frequency  "
			    " dimension as f_grid");
      }
     
      if(scat_i_p.nvitrines() != 2){
	throw runtime_error("scat_i_p should have only two elements "
			    "in pressure grid which corresponds "
			    "to the two pressure surfaces");
      }
      
     
      Index N_za = scat_i_p.npages() ;
   
      Index N_aa = scat_i_p.nrows();
   
      Index N_i = scat_i_p.ncols();
         
      //1. interpolation - pressure grid
      
      
      i_field.resize((cloudbox_limits[1]- cloudbox_limits[0])+1, 1, 1,  N_za, N_aa, N_i);
      
      /*the old grid is having only two elements, corresponding to the 
	cloudbox_limits and the new grid have elements corresponding to
	all grid points inside the cloudbox plus the cloud_box_limits*/

      ArrayOfGridPos p_gp((cloudbox_limits[1]- cloudbox_limits[0])+1);
      
      gridpos(p_gp,
	      p_grid[Range(cloudbox_limits[0], 
			   2,
			   (cloudbox_limits[1]- cloudbox_limits[0]))],
	      p_grid[Range(cloudbox_limits[0], 
			   (cloudbox_limits[1]- cloudbox_limits[0])+1)]);
      
      Matrix itw((cloudbox_limits[1]- cloudbox_limits[0])+1, 2);
      interpweights ( itw, p_gp );
      
      for (Index za_index = 0; za_index < N_za ; ++ za_index)
	{
	  for (Index aa_index = 0; aa_index < N_aa ; ++ aa_index)
	    {
	      for (Index i = 0 ; i < N_i ; ++ i)
		{
		  
		  VectorView target_field = i_field(Range(joker),
						    0,
						    0,
						    za_index,
						    aa_index,
						    i);
		  
		  ConstVectorView source_field = scat_i_p(scat_f_index,
							  Range(joker),    
							  0,
							  0,
							  za_index,
							  aa_index,
							  i);
		  
		  interp(target_field,
			 itw,
			 source_field,
			 p_gp);
		}
	      
	    }
	}
      xml_write_to_file("i_field_Set.xml", i_field);
      //cout<<i_field<<"\n";
    }
  if(atmosphere_dim == 3)
    {
      Index  N_f = scat_i_p.nlibraries();
      if (scat_i_lat.nlibraries() != N_f || 
	  scat_i_lon.nlibraries() != N_f){
	
	throw runtime_error(" scat_i_p, scat_i_lat, scat_i_lon should have  "
			    "same frequency dimension");
      }
      Index N_p = scat_i_lat.nvitrines();
      if(scat_i_lon.nvitrines() != N_p ||
	 p_grid.nelem()         != N_p ){
	throw runtime_error("scat_i_lat and scat_i_lon should have  "
			    "same pressure grid dimension as p_grid");
      }
      
      Index N_lat = scat_i_p.nshelves();
      
      if(scat_i_lon.nshelves() != N_lat ||
	 lat_grid.nelem()      != N_lat){
	throw runtime_error("scat_i_p and scat_i_lon should have  "
			    "same latitude grid dimension as lat_grid");
      }
  
      Index N_lon = scat_i_p.nbooks();
      if(scat_i_lat.nbooks() != N_lon ||
	 lon_grid.nelem()    != N_lon ){
	throw runtime_error("scat_i_p and scat_i_lat should have  "
			    "same longitude grid dimension as lon_grid");
      }
      if(scat_i_p.nvitrines() != 2){
	throw runtime_error("scat_i_p should have only two elements "
			    "in pressure grid which corresponds "
			    "to the two pressure surfaces");
      }
      
      if(scat_i_lat.nshelves() != 2){
	throw runtime_error("scat_i_lat should have only two elements "
			    "in latitude grid which corresponds "
			    "to the two latitude surfaces");
	
      }
      if(scat_i_lon.nbooks() != 2){
	throw runtime_error("scat_i_lon should have only two elements "
			    "in longitude grid which corresponds "
			    "to the two longitude surfaces");
	
      }
      Index N_za = scat_i_p.npages() ;
      if (scat_i_lat.npages() != N_za || 
	  scat_i_lon.npages() != N_za){
	
	throw runtime_error(" scat_i_p, scat_i_lat, scat_i_lon should have  "
			    "same dimension for zenith angles");
      }
      Index N_aa = scat_i_p.nrows();
      if (scat_i_lat.nrows() != N_aa || 
	  scat_i_lon.nrows() != N_aa){
	
	throw runtime_error(" scat_i_p, scat_i_lat, scat_i_lon should have  "
			    "same dimension for azimuth angles");
      }
      Index N_i = scat_i_p.ncols();
      if (scat_i_lat.ncols() != N_i || 
	  scat_i_lon.ncols() != N_i){
	
	throw runtime_error(" scat_i_p, scat_i_lat, scat_i_lon should have  "
			    "same value for stokes_dim and can take only"
			    "values 1,2,3 or 4");
      }
      
      //1. interpolation - pressure grid, latitude grid and longitude grid
      
   
      //i_field
      
      i_field.resize((cloudbox_limits[1]- cloudbox_limits[0])+1, 
		     (cloudbox_limits[3]- cloudbox_limits[2])+1,
		     (cloudbox_limits[5]- cloudbox_limits[4])+1,
		     N_za, 
		     N_aa,
		     N_i);
      
      
      ArrayOfGridPos p_gp((cloudbox_limits[1]- cloudbox_limits[0])+1);
      ArrayOfGridPos lat_gp((cloudbox_limits[3]- cloudbox_limits[2])+1);
      ArrayOfGridPos lon_gp((cloudbox_limits[5]- cloudbox_limits[4])+1);

      /*the old grid is having only two elements, corresponding to the 
	cloudbox_limits and the new grid have elements corresponding to
	all grid points inside the cloudbox plus the cloud_box_limits*/
      
      gridpos(p_gp,
	      p_grid[Range(cloudbox_limits[0], 
			   2,
			   (cloudbox_limits[1]- cloudbox_limits[0]))],
	      p_grid[Range(cloudbox_limits[0], 
			   (cloudbox_limits[1]- cloudbox_limits[0])+1)]);
      gridpos(lat_gp,
	      lat_grid[Range(cloudbox_limits[2], 
			     2,
			     (cloudbox_limits[3]- cloudbox_limits[2]))],
	      lat_grid[Range(cloudbox_limits[2], 
			     (cloudbox_limits[3]- cloudbox_limits[2])+1)]);
      gridpos(lon_gp,
	      lon_grid[Range(cloudbox_limits[4], 
			     2,
			     (cloudbox_limits[5]- cloudbox_limits[4]))],
	      lon_grid[Range(cloudbox_limits[4], 
			     (cloudbox_limits[5]- cloudbox_limits[4])+1)]);
      
      //interpolation weights corresponding to pressure, latitude and 
      //longitude grids.

      Matrix itw_p((cloudbox_limits[1]- cloudbox_limits[0])+1, 2);
      Matrix itw_lat((cloudbox_limits[3]- cloudbox_limits[2])+1, 2);
      Matrix itw_lon((cloudbox_limits[5]- cloudbox_limits[4])+1, 2);

      interpweights ( itw_p, p_gp );
      interpweights ( itw_lat, lat_gp );
      interpweights ( itw_lon, lon_gp );

      // interpolation - pressure grid
      for (Index lat_index = cloudbox_limits[2]; 
	   lat_index < cloudbox_limits[3] ; ++ lat_index)
	{
	  for (Index lon_index = cloudbox_limits[4]; 
	       lon_index < cloudbox_limits[5] ; ++ lon_index)
	    {
	      for (Index za_index = 0; za_index < N_za ; ++ za_index)
		{
		  for (Index aa_index = 0; aa_index < N_aa ; ++ aa_index)
		    {
		      for (Index i = 0 ; i < N_i ; ++ i)
			{
			  
			  VectorView target_field = i_field(Range(joker),
							    lat_index,
							    lon_index,
							    za_index,
							    aa_index,
							    i);
			  
			  ConstVectorView source_field = scat_i_p(scat_f_index,
								  Range(joker),    
								  lat_index,
								  lon_index,
								  za_index,
								  aa_index,
								  i);
			  
			  interp(target_field,
				 itw_p,
				 source_field,
				 p_gp);
			}
		    }
		}
	    } 
	}
      //interpolation latitude
      for (Index p_index = cloudbox_limits[0]; 
	   p_index < cloudbox_limits[1] ; ++ p_index)
	{
	  for (Index lon_index = cloudbox_limits[4]; 
	       lon_index < cloudbox_limits[5] ; ++ lon_index)
	    {
	      for (Index za_index = 0; za_index < N_za ; ++ za_index)
		{
		  for (Index aa_index = 0; aa_index < N_aa ; ++ aa_index)
		    {
		      for (Index i = 0 ; i < N_i ; ++ i)
			{
			  
			  VectorView target_field = i_field(p_index,
							    Range(joker),
							    lon_index,
							    za_index,
							    aa_index,
							    i);
			  
			  ConstVectorView source_field = scat_i_p(scat_f_index,
								  p_index,
								  Range(joker),    
								  lon_index,
								  za_index,
								  aa_index,
								  i);
			  
			  interp(target_field,
				 itw_lat,
				 source_field,
				 lat_gp);
			}
		    }
		}
	    } 
	}
      //interpolation -longitude
      for (Index p_index = cloudbox_limits[0]; 
	   p_index < cloudbox_limits[1] ; ++ p_index)
	{
	  for (Index lat_index = cloudbox_limits[2]; 
	       lat_index < cloudbox_limits[3] ; ++ lat_index)
	    {
	      for (Index za_index = 0; za_index < N_za ; ++ za_index)
		{
		  for (Index aa_index = 0; aa_index < N_aa ; ++ aa_index)
		    {
		      for (Index i = 0 ; i < N_i ; ++ i)
			{
			  
			  VectorView target_field = i_field(p_index,
							    lat_index,
							    Range(joker),
							    za_index,
							    aa_index,
							    i);
			  
			  ConstVectorView source_field = scat_i_p(scat_f_index,
								  p_index,    
								  lat_index,
								  Range(joker),
								  za_index,
								  aa_index,
								  i);
			  
			  interp(target_field,
				 itw_lon,
				 source_field,
				 lon_gp);
			}
		    }
		}
	    } 
	}
      //end of interpolation
    }//ends atmosphere_dim = 3
}
      
