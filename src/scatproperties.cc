/*!
  \file   scatproperties.cc
  \author Sreerekha T.R. 
  \date   Fri May 10 13:42:40 2002
  
  \brief  This file has functions for calculating extinction and phase matrices 
  from amplitude matrix 
*/

#include <iostream>
#include "scatproperties.h"  
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric SPEED_OF_LIGHT; 
#define Re_S11 amp_coeffs[0]
#define Im_S11 amp_coeffs[1]
#define Re_S12 amp_coeffs[2]
#define Im_S12 amp_coeffs[3]
#define Re_S21 amp_coeffs[4]
#define Im_S21 amp_coeffs[5]
#define Re_S22 amp_coeffs[6]
#define Im_S22 amp_coeffs[7]
#define Re_S11ij amp_coeffs(i,j,0)
#define Im_S11ij amp_coeffs(i,j,1)
#define Re_S12ij amp_coeffs(i,j,2)
#define Im_S12ij amp_coeffs(i,j,3)
#define Re_S21ij amp_coeffs(i,j,4)
#define Im_S21ij amp_coeffs(i,j,5) 
#define Re_S22ij amp_coeffs(i,j,6) 
#define Im_S22ij amp_coeffs(i,j,7)
 
//! Calculates extinction cross-section matrix for a single particle type for
// given combination of  angles. 
/*! 
  
  This function calculates the extinction matrix from the amplitude matrix
  elements for given angles theta, phi, theta', phi'. It is called by the 
  method ext_mat-sptCalc which selects the right angles of all the variables.
  This function is called for each particle type that is specified.

  This function requires as input the amplitude matrix for the particle whose 
  components are represented as 
  
  Numeric § Re_S11=amp_mat(part_type, theta, phi, theta', phi', 0)
  Numeric § Im_S11=amp_mat(part_type, theta, phi, theta', phi', 1)
  Numeric § Re_S12=amp_mat(part_type, theta, phi, theta', phi', 2)
  Numeric § Im_S12=amp_mat(part_type, theta, phi, theta', phi', 3)
  Numeric § Re_S21=amp_mat(part_type, theta, phi, theta', phi', 4)
  Numeric § Im_S21=amp_mat(part_type, theta, phi, theta', phi', 5)
  Numeric § Re_S22=amp_mat(part_type, theta, phi, theta', phi', 6)
  Numeric § Im_S22=amp_mat(part_type, theta, phi, theta', phi', 7)
  
  The output of this funtion is a Matrix, ext(stokes_dim,stokes_dim),
  where stokes_dim is stokes vector dimension specified by the user.

  The input to this function namely amp_coeffs is a vector, amp_coeffs(8).
  In the method ext_mat_sptCalc, the right angles and the right particle
  are picked to give this vector. Another input to this function is 
  the corresponding frequency.

  \param ext Output :Extinciton Matrix for given particle type and 
  combination of angles
  \param amp_coeffs Input :amplitude matrix for given particle type 
  and combination of angle
  \param freq Input : frequency
*/

void amp2ext(MatrixView ext,
	     ConstVectorView amp_coeffs,
	     const Numeric& freq)
{
  Index stokes_dim = ext.nrows();

  if (stokes_dim > 4 || stokes_dim <1){
    throw runtime_error("the dimension of stokes vector"
                         "can be only 1,2,3, or 4");
  }
  
  assert (is_size(ext, stokes_dim, stokes_dim));
  assert (is_size(amp_coeffs, 8));
  
  const Numeric wavelength = SPEED_OF_LIGHT / freq;
  
  ext(0, 0) = wavelength * (Im_S11 + Im_S22) ; 
  
  if( 1 == stokes_dim ){
    return;
  }

  //only the first element is required if stokes_dim =1

  ext(0, 1) = wavelength * (Im_S11 - Im_S22);
  ext(1, 0) = ext(0, 1);
  ext(1, 1) = ext(0, 0);
  
  if(2 == stokes_dim){
    return;
  }

  // if stokes_dim =2 only these 4 elements need be evaluated.

  ext(0, 2) = - wavelength * (Im_S12 + Im_S21);
  ext(1, 2) = wavelength * (Im_S21 - Im_S12);
  ext(2, 0) = ext(0, 2);
  ext(2, 1) = - ext(1, 2);
  ext(2, 2) = ext(0, 0);
  
  if(3 == stokes_dim){
    return;
  }

  // if stokes_dim =3 only these 9 elements need be evaluated.


  ext(0, 3) = wavelength * (Re_S21 - Re_S12);
  ext(1, 3) = - wavelength*(Re_S12 + Re_S21);
  ext(2, 3) = wavelength * (Re_S22 - Re_S11);
  ext(3, 0) = ext(0, 3);
  ext(3, 1) = - ext(1, 3);
  ext(3, 2) = - ext(2, 3);
  ext(3, 3) = ext(0, 0);
  if(4 == stokes_dim){
    return;
  }

   // if stokes_dim =4 all 16 elements need be evaluated.
  
}

//! Calculates phase  matrix for a single particle type for
// given combination of  outgoing angles. 

/*!
  This function calculatest the phase matrix from the amplitude matrix
  elements for the given scattering (outgoing) angles. This function is
  called by the method pha_mat_sptCalc which selects the right outgoing 
  angles.  The function is called for each particle type that is 
  specified.

  The output of this function is a Tensor4, 
  phasemat(theta', phi', stokes_dim, stokes_dim) where stokes_dim is 
  stokes vector dimension specified by the user. This function requires 
  as input the amplitude matrix which is a Tensor3, 
  amp_coeffs(theta', phi', 8).

  \param phasemat Output: phase matrix for the single particle type
  \param amp_coeffs  Input : amplitude matrix
*/

void amp2pha(Tensor4View phasemat,
	     ConstTensor3View amp_coeffs)
{
  Index stokes_dim = phasemat.nrows();
  if (stokes_dim > 4 || stokes_dim <1){
    throw runtime_error("the dimension of stokes vector" 
			" can be only 1,2,3, or 4");
  }
  Index nza = phasemat.nbooks();
  Index naa = phasemat.npages();
  assert (is_size(amp_coeffs, nza, naa, 8));
  assert (is_size(phasemat, nza, naa, stokes_dim, stokes_dim));
  
  for (Index i = 0; i < nza; ++i)
    {
      for (Index j = 0; j < naa; ++j)
	{
	  phasemat(i, j, 0, 0) = 0.5 * (Re_S11ij * Re_S11ij + 
					Im_S11ij * Im_S11ij + 
					Re_S12ij * Re_S12ij + 
					Im_S12ij * Im_S12ij +
					Re_S21ij * Re_S21ij + 
					Im_S21ij * Im_S21ij +
					Re_S22ij * Re_S22ij + 
					Im_S22ij * Im_S22ij);
	}
    }
  
  if(1 == stokes_dim){
    return;
  }
  
  //only the first element is required if stokes_dim =1
  
  for (Index i = 0; i < nza; ++i)
    {
      for (Index j = 0; j < naa; ++j)
	{
	  phasemat(i, j, 0, 1) = 0.5 * (Re_S11ij * Re_S11ij + 
					Im_S11ij * Im_S11ij - 
					Re_S12ij * Re_S12ij - 
					Im_S12ij * Im_S12ij +
					Re_S21ij * Re_S21ij +
					Im_S21ij * Im_S21ij -
					Re_S22ij * Re_S22ij -
					Im_S22ij * Im_S22ij);

	  phasemat(i, j, 1, 0) = 0.5 * (Re_S11ij * Re_S11ij + 
					Im_S11ij * Im_S11ij +
					Re_S12ij * Re_S12ij + 
					Im_S12ij * Im_S12ij - 
					Re_S21ij * Re_S21ij - 
					Im_S21ij * Im_S21ij - 
					Re_S22ij * Re_S22ij - 
					Im_S22ij * Im_S22ij);

	  phasemat(i, j, 1, 1) = 0.5 * (Re_S11ij * Re_S11ij + 
					Im_S11ij * Im_S11ij - 
					Re_S12ij * Re_S12ij - 
					Im_S12ij * Im_S12ij -
					Re_S21ij * Re_S21ij - 
					Im_S21ij * Im_S21ij + 
					Re_S22ij + Re_S22ij + 
					Im_S22ij * Im_S22ij);
	}
    }
	  
  if(2 == stokes_dim){
    return;
  }

  // if stokes_dim =2 only these 4 elements need be evaluated.

  for (Index i = 0; i < nza; ++i)
    {
      for (Index j = 0; j < naa; ++j)
	{
	  phasemat(i, j, 0, 2) = - (Re_S11ij * Re_S12ij +
				    Im_S11ij * Im_S12ij + 
				    Re_S22ij * Re_S21ij + 
				    Im_S22ij * Im_S21ij);
	  
	  phasemat(i, j, 1, 2) = - (Re_S11ij * Re_S12ij +
				    Im_S11ij * Im_S12ij -
				    Re_S22ij * Re_S21ij -
				    Im_S22ij * Im_S21ij);
	  
	  phasemat(i, j, 2, 0) = - (Re_S11ij * Re_S21ij + 
				    Im_S11ij * Im_S21ij + 
				    Re_S22ij * Re_S12ij + 
				    Im_S22ij * Im_S12ij);
	  
	  phasemat(i, j, 2, 1) = - (Re_S11ij * Re_S21ij + 
				    Im_S11ij * Im_S21ij -
				    Re_S22ij * Re_S12ij - 
				    Im_S22ij * Im_S12ij);
	  
	  phasemat(i, j, 2, 2) =  (Re_S11ij * Re_S22ij + 
				   Im_S11ij * Im_S22ij +
				   Re_S12ij * Re_S21ij + 
				   Im_S12ij * Im_S21ij);
	}
    }
  
  if(3 == stokes_dim){
    return;
  }

   // if stokes_dim =3 only these 9 elements need be evaluated.

  for (Index i = 0; i < nza; ++i)
    {
      for (Index j = 0; j < naa; ++j)
	{
	  phasemat(i, j, 0, 3) = - (Im_S11ij * Re_S12ij -
				    Re_S11ij * Im_S12ij -
				    Im_S22ij * Re_S21ij + 
				    Re_S22ij * Im_S21ij);
	  
	  phasemat(i, j, 1, 3) = - (Im_S11ij * Re_S12ij - 
				    Re_S11ij * Im_S12ij + 
				    Im_S22ij * Re_S21ij - 
				    Re_S22ij * Im_S21ij);
	  
	  phasemat(i, j, 2, 3) = (Im_S11ij * Re_S22ij - 
				  Re_S11ij * Im_S22ij +
				  Im_S21ij * Re_S12ij - 
				  Re_S21ij * Im_S12ij);
	  
	  phasemat(i, j, 3, 0) = - (Im_S21ij * Re_S11ij - 
				    Re_S21ij * Im_S11ij +
				    Im_S22ij * Re_S12ij - 
				    Re_S22ij * Im_S12ij);
	  
	  phasemat(i, j, 3, 1) = - (Im_S21ij * Re_S11ij - 
				    Re_S21ij * Im_S11ij - 
				    Im_S22ij * Re_S12ij + 
				    Re_S22ij * Im_S12ij);
	  
	  phasemat(i, j, 3, 2) = (Im_S22ij * Re_S11ij - 
				  Re_S22ij * Im_S11ij - 
				  Im_S12ij * Re_S21ij + 
				  Re_S12ij * Im_S21ij); 
	  
	  phasemat(i, j, 3, 3) = (Re_S22ij * Re_S11ij + 
				  Im_S22ij * Im_S11ij - 
				  Re_S12ij * Re_S21ij - 
				  Im_S12ij * Im_S21ij); 
	}
    }
  
  if(4 == stokes_dim){
    return;
  }

  // if stokes_dim =4 all 16 elements need be evaluated.
}


//! calculates absorption coefficeint vector for given angles.
/*! 
  This function calculates absorption vector from phase matrix elements
  and extinction matrix elements for given angle combinations.  This 
  function is called by the method abs_vec_sptCalc which selects the 
  right angles.  The function is called for each particle type that is
  specified.

  The output of the function is a vector, abs(stokes_dim) where 
  stokes_dim is the stokes vector dimension specified by the user.
  It require as input the extinction matrix of dimension,
  ext(stokes_dim, stokes_dim) and phase matrix of dimension,
  pha(theta', phi', stokes_dim, stokes_dim).
  
  \param abs Output : absorption coefficient vector 
  \param ext  Input : Extinction Matrix
  \param pha  Input : Phase matrix
*/

void amp2abs(VectorView abs,
	     ConstMatrixView ext,
	     ConstTensor4View pha)
{
  Index stokes_dim = abs.nelem();
  if (stokes_dim > 4 || stokes_dim <1){
    throw runtime_error("the dimension of stokes vector" 
			"can be only 1,2,3, or 4");
  }
  Index nza = pha.nbooks();
  Index naa = pha.npages();
  assert (is_size(pha, nza, naa, stokes_dim, stokes_dim));
  assert (is_size(ext, stokes_dim, stokes_dim));
  assert (is_size(abs, stokes_dim));
  Vector za_grid(nza);
  Vector aa_grid(naa);
  Numeric Z11_integrated, Z21_integrated, Z31_integrated, Z41_integrated ;
  Matrix Z11_mat = pha (Range(joker), Range(joker), 0, 0);
  for (Index i = 0;i < nza; ++i)
    {
      for (Index j = 0; j < naa; ++j)
	{
	  Z11_integrated = AngIntegrate_trapezoid(Z11_mat,
						  za_grid,
						  aa_grid);
	}
    }

  abs[0] = ext(0, 0) - Z11_integrated;
  
  if(1 == stokes_dim){
    return;
  }
   
  //only the first element is required if stokes_dim =1
  
  Matrix Z21_mat = pha (Range(joker), Range(joker), 1, 0);
  for (Index i = 0; i < nza; ++i)
    {
      for (Index j = 0; j < naa; ++j)
	{
	  Z21_integrated = AngIntegrate_trapezoid(Z21_mat,
						  za_grid,
						  aa_grid);
	}
    }

  abs[1] = ext(1, 1)- Z21_integrated;
  
  if(2 == stokes_dim){
    return;
  }

  // if stokes_dim =2 only these 2 elements need be evaluated.

  Matrix Z31_mat = pha (Range(joker), Range(joker), 2, 0);
  for (Index i = 0; i < nza; ++i)
    {
      for (Index j = 0; j < naa; ++j)
	{
	  Z31_integrated = AngIntegrate_trapezoid(Z31_mat,
						  za_grid,
						  aa_grid);
	}
    }
  
  abs[2] = ext(2, 2) - Z31_integrated;
  
  if(3 == stokes_dim){
    return;
  }

   // if stokes_dim =3 only these 3 elements need be evaluated.

  Matrix Z41_mat = pha (Range(joker), Range(joker), 3, 0);
  for (Index i = 0; i < nza; ++i)
    {
      for (Index j = 0; j < naa; ++j)
	{
	  
	  Z41_integrated = AngIntegrate_trapezoid(Z41_mat,
						  za_grid,
						  aa_grid);
	}
    }

  abs[3] = ext(3, 3) - Z41_integrated;

  if(4 == stokes_dim){
    return;
  }

  // if stokes_dim =4 all 4 elements need be evaluated.
}

/** 
 * 
 * 
 * @param Integrand The Matrix to be integrated
 * @param za_grid Input : The zenith angle grid 
 * @param aa_grid Input : The azimuth angle grid 
 * 
 * @return The resulting integral
 */
Numeric AngIntegrate_trapezoid(ConstMatrixView Integrand,
			       ConstVectorView za_grid,
			       ConstVectorView aa_grid)
{
  //is_size(za_grid.nelem(),aa_grid.nelem());
  
  Index n = za_grid.nelem();
  Index m = aa_grid.nelem();
  Vector res1(n);
  assert (is_size(Integrand, n, m));
  for (Index i = 0; i < n ; ++i)
    {
      res1[i] = 0.0;
      Numeric sintheta = sin(za_grid[i] * DEG2RAD);
      
      for (Index j = 0; j < m - 1; ++j)
	{
	  res1[i] +=  0.5 * (Integrand(i, j) + Integrand(i, j + 1)) *
	    (aa_grid[j + 1] - aa_grid[j]) * sintheta;
	}
    }
  
  Numeric res = 0.0;
  for (Index i = 0; i < n - 1; ++i)
    {
      res += 0.5 *  (res1[i] + res1[i + 1]) * 
	(za_grid[i + 1] - za_grid[i]);
    }
  
  cout<<res<<"\n";
  return res;
}

//! Extinction Coefficient Matrix for the particle 
/*! 
  
  This function sums up the extinction matries for all particle types weighted with particle number density
  \param ext_mat_part Output : physical extinction coefficient for the particles for given angles. 
  \param ext_mat_spt Input : extinction matrix for the scattering species
  \param pnd Input : particle number density givs the local concentration for all particles.
*/
//FIXME; change all matrix ext_mat_spt to Tensor3 ext_mat_spt(l,j,k)
void ext_mat_partCalc(MatrixView ext_mat_part,
		      MatrixView ext_mat_spt,
		      VectorView pnd)
{
  Index N_pt, j, k, l, N_i;
  for (j = 0; j < N_i; j++)
    {
      for (k = 0; k < N_i; k++)
	ext_mat_part(j, k) = 0.0;// Initialisation
    }
  for (l = 0; l < N_pt; l++)
    {
      for (j = 0; j < N_i; j++)
	{
	  for (k = 0; k < N_i; k++)
	    ext_mat_spt(j, k) *= pnd[l]; //multiplying with particle number density.
	  ext_mat_part(j, k) += ext_mat_spt(j,k); // summing up for all particle types.
	}
    } 
} 
 
