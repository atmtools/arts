
/*!
  \file   scatproperties.cc
  \author Sreerekha T.R. 
  \date   Fri May 10 13:42:40 2002
  
  \brief  This file has functions for calculating extinction and phase matrices 
  from amplitude matrix 
*/

#include <stdexcept>
#include <iostream>
#include "scatproperties.h"  
#include "math_funcs.h"
#include "logic.h"

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

  //cout<<"Wavelength corresponding to this frequency"<<wavelength<<"\n";
  //  cout<<"Imaginariy part of S11"<<Im_S11<<"\n";
  //cout << "Amplitude matrix " << amp_coeffs  << " \n " ;
  ext(0, 0) = wavelength * (Im_S11 + Im_S22) * 1e-6; 
  
  if( 1 == stokes_dim ){
    return;
  }

  //only the first element is required if stokes_dim =1

  ext(0, 1) = wavelength * (Im_S11 - Im_S22) * 1e-6 ;
  ext(1, 0) = ext(0, 1);
  ext(1, 1) = ext(0, 0);
  
  if(2 == stokes_dim){
    return;
  }
  
  // if stokes_dim =2 only these 4 elements need be evaluated.

  ext(0, 2) = - wavelength * (Im_S12 + Im_S21) * 1e-6 ;
  ext(1, 2) = wavelength * (Im_S21 - Im_S12) * 1e-6 ;
  ext(2, 0) = ext(0, 2);
  ext(2, 1) = - ext(1, 2);
  ext(2, 2) = ext(0, 0);

  if(3 == stokes_dim){
    return;
  }

  // if stokes_dim =3 only these 9 elements need be evaluated.


  ext(0, 3) = wavelength * (Re_S21 - Re_S12) * 1e-6 ;
  ext(1, 3) = - wavelength*(Re_S12 + Re_S21) * 1e-6  ;
  ext(2, 3) = wavelength * (Re_S22 - Re_S11) * 1e-6 ;
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
  //cout<<"\n nza, naa"<<nza<<" "<<naa<<"\n"; 
  //assert (is_size(amp_coeffs, nza, naa, 8));
  assert (is_size(phasemat, nza, naa, stokes_dim, stokes_dim));
  
  for (Index i = 0; i < nza; ++i)
    {
      for (Index j = 0; j < naa; ++j)
	{
	   phasemat(i, j, 0, 0) = 0.5 * 1e-12 * (Re_S11ij * Re_S11ij + 
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
	  phasemat(i, j, 0, 1) = 0.5 * 1e-12  * (Re_S11ij * Re_S11ij + 
						 Im_S11ij * Im_S11ij - 
						 Re_S12ij * Re_S12ij - 
						 Im_S12ij * Im_S12ij +
						 Re_S21ij * Re_S21ij +
						 Im_S21ij * Im_S21ij -
						 Re_S22ij * Re_S22ij -
						 Im_S22ij * Im_S22ij);
	  
	  phasemat(i, j, 1, 0) = 0.5 * 1e-12  * (Re_S11ij * Re_S11ij + 
						 Im_S11ij * Im_S11ij +
						 Re_S12ij * Re_S12ij + 
						 Im_S12ij * Im_S12ij - 
						 Re_S21ij * Re_S21ij - 
						 Im_S21ij * Im_S21ij - 
						 Re_S22ij * Re_S22ij - 
						 Im_S22ij * Im_S22ij);
	  
	  phasemat(i, j, 1, 1) = 0.5 * 1e-12  * (Re_S11ij * Re_S11ij + 
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
	  phasemat(i, j, 0, 2) = -  1e-12 * (Re_S11ij * Re_S12ij +
					     Im_S11ij * Im_S12ij + 
					     Re_S22ij * Re_S21ij + 
					     Im_S22ij * Im_S21ij);
	  
	  phasemat(i, j, 1, 2) = - 1e-12 *  (Re_S11ij * Re_S12ij +
					     Im_S11ij * Im_S12ij -
					     Re_S22ij * Re_S21ij -
					     Im_S22ij * Im_S21ij);
	  
	  phasemat(i, j, 2, 0) = - 1e-12 *  (Re_S11ij * Re_S21ij + 
					     Im_S11ij * Im_S21ij + 
					     Re_S22ij * Re_S12ij + 
					     Im_S22ij * Im_S12ij);
	  
	  phasemat(i, j, 2, 1) = - 1e-12 *  (Re_S11ij * Re_S21ij + 
					     Im_S11ij * Im_S21ij -
					     Re_S22ij * Re_S12ij - 
					     Im_S22ij * Im_S12ij);
	  
	  phasemat(i, j, 2, 2) = 1e-12 *  (Re_S11ij * Re_S22ij + 
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
	  phasemat(i, j, 0, 3) = - 1e-12 * (Im_S11ij * Re_S12ij -
					    Re_S11ij * Im_S12ij -
					    Im_S22ij * Re_S21ij + 
					    Re_S22ij * Im_S21ij);
	  
	  phasemat(i, j, 1, 3) = - 1e-12 * (Im_S11ij * Re_S12ij - 
					    Re_S11ij * Im_S12ij + 
					    Im_S22ij * Re_S21ij - 
					    Re_S22ij * Im_S21ij);
	  
	  phasemat(i, j, 2, 3) = 1e-12 * (Im_S11ij * Re_S22ij - 
					  Re_S11ij * Im_S22ij +
					  Im_S21ij * Re_S12ij - 
					  Re_S21ij * Im_S12ij);
	  
	  phasemat(i, j, 3, 0) = - 1e-12 * (Im_S21ij * Re_S11ij - 
					    Re_S21ij * Im_S11ij +
					    Im_S22ij * Re_S12ij - 
					    Re_S22ij * Im_S12ij);
	  
	  phasemat(i, j, 3, 1) = - 1e-12 * (Im_S21ij * Re_S11ij - 
					    Re_S21ij * Im_S11ij - 
					    Im_S22ij * Re_S12ij + 
					    Re_S22ij * Im_S12ij);
	  
	  phasemat(i, j, 3, 2) = 1e-12 * (Im_S22ij * Re_S11ij - 
					  Re_S22ij * Im_S11ij - 
					  Im_S12ij * Re_S21ij + 
					  Re_S12ij * Im_S21ij); 
	  
	  phasemat(i, j, 3, 3) = 1e-12 * (Re_S22ij * Re_S11ij + 
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

//! Calculates extinction coefficeint for the convergence test.
/*! 
  For the convergence test, an extinction coefficient which only contains
  extinction due to particle scattering is required. This can be 
  calculated by integrating over the phase matrix, which is done in this
  function.

  Output:
  \param ext_conv Extinction coefficient for the convergence test.
  Input:
  \param pha  Phase matrix.
  \param za_grid  Zenith angle grid
  \param aa_grid  Azimuth angle grid
*/
void amp2ext_scat(MatrixView ext_scat,
                  ConstTensor4View pha,
                  ConstVectorView za_grid,
                  ConstVectorView aa_grid)
{
  
  Index stokes_dim = pha.nrows();
  if (stokes_dim > 4 || stokes_dim <1){
    throw runtime_error("the dimension of stokes vector" 
			"can be only 1,2,3, or 4");
  }

  Index nza = pha.nbooks();
  Index naa = pha.npages();
  assert (is_size(pha, nza, naa, stokes_dim, stokes_dim));
  assert (is_size(za_grid,nza));
  assert (is_size(aa_grid,naa));
  //Vector za_grid(nza);
  //Vector aa_grid(naa);

  //  ext_conv.resize(stokes_dim, stokes_dim);
  ext_scat = 0.;

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

  ext_scat(0,0) = Z11_integrated;
    
  if(1 == stokes_dim){
    return;
  }
   
  //only the first element is required if stokes_dim =1
  
  //FIXME: How does the extinction matrix without absorption look for 
  // higher dimensions?

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

  ext_scat(1,0) = Z21_integrated;
  
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
  
  ext_scat(2,0) = Z31_integrated;
  
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

  ext_scat(3,0) = Z41_integrated;

  if(4 == stokes_dim){
    return;
  }

  // if stokes_dim =4 all 4 elements need be evaluated.
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
  \param za_grid  Input : Zenith angle grid
  \param aa_grid  Input : Azimuth angle grid
*/

void amp2abs(VectorView abs,
	     ConstMatrixView ext,
	     ConstTensor4View pha,
	     ConstVectorView za_grid,
	     ConstVectorView aa_grid)
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
  assert (is_size(za_grid,nza));
  assert (is_size(aa_grid,naa));
  //Vector za_grid(nza);
  //Vector aa_grid(naa);
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



