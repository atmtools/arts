/*!
  \file   scatproperties.cc
  \author Sreerekha T.R. 
  \date   Fri May 10 13:42:40 2002
  
  \brief  This file has functions for calculating extinction and phase matrices 
  from amplitude matrix 
*/
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

//! Calculates extinction cross-section matrix for a single particle type for
// given combination of  angles. 
/*! 
  
  This function calculates the extinction matrix from the amplitude matrix elements
  for given angles theta, phi, theta', phi'.
  This function requires as input the amplitude matrix for the particle whose 
  components are represented as 
  
  Numeric § Re_S11=amp_mat(part_type,theta, phi, theta', phi',0)
  Numeric § Im_S11=amp_mat(part_typetheta, phi, theta', phi',1)
  Numeric § Re_S12=amp_mat(part_typetheta, phi, theta', phi',2)
  Numeric § Im_S12=amp_mat(part_typetheta, phi, theta', phi',3)
  Numeric § Re_S21=amp_mat(part_typetheta, phi, theta', phi',4)
  Numeric § Im_S21=amp_mat(part_typetheta, phi, theta', phi',5)
  Numeric § Re_S22=amp_mat(part_typetheta, phi, theta', phi',6)
  Numeric § Im_S22=amp_mat(part_typetheta, phi, theta', phi',7)
  
  ext_mat_spt(part_type,stokes_dim,stokes_dim) is tensor 3.  n_p denotes the particle type and n_i denotes the stokes vector dimension.
  \param ext Output :Extinciton Matrix for given particle type and combination of angles
  \param amp_coeffs Input :amplitude matrix for given particle type and combination of angle
  \param freq Input : frequency
*/
//FIXME; change all matrix ext_mat_spt to Tensor3 ext_mat_spt(l,j,k)
void amp2ext(MatrixView ext,
	     ConstVectorView amp_coeffs,
	     const Numeric& freq)
{
  assert (is_size(ext,4,4));
  assert (is_size(amp_coeffs,8));
  const Numeric wavelength=SPEED_OF_LIGHT/freq;
  ext(0,0) = wavelength * (Im_S11+Im_S22) ;
  ext(1,1) = wavelength*(Im_S11+Im_S22);
  ext(2,2) = wavelength*(Im_S11+Im_S22);
  ext(3,3) = wavelength*(Im_S11+Im_S22);
  ext(0,1) = wavelength*(Im_S11-Im_S22);
  ext(1,0) = wavelength*(Im_S11-Im_S22);
  ext(0,2) = -wavelength*(Im_S12+Im_S21);
  ext(2,0) = -wavelength*(Im_S12+Im_S21);
  ext(0,3) = wavelength*(Re_S21-Re_S12);
  ext(3,0) = wavelength*(Re_S21-Re_S12);
  ext(1,2) = wavelength*(Im_S21-Im_S12);
  ext(2,1) = -wavelength*(Im_S21-Im_S12);
  ext(1,3) = -wavelength*(Re_S12+Re_S21);
  ext(3,1) = wavelength*(Re_S12+Re_S21);
  ext(2,3) = wavelength*(Re_S22-Re_S11);
  ext(3,2) = -wavelength*(Re_S22-Re_S11);
}

//! 
/*! 
  \param phasemat Output: phase matrix for the single particle type
  \param amp_coeffs  Input : amplitude matrix
*/
void amp2phamat(Tensor4View phasemat,ConstVectorView amp_coeffs)
{
  Index nza = phasemat.nbooks();
  Index naa = phasemat.npages();
  for (Index i=0;i<nza;++i)
    {
      for (Index j=0;j<naa;++j)
	//phamat(i,j,0,0) = 0;
	//this has to be computed fully.
	cout<<j;
    }
}
//! The function calculates absorption coefficeint vector for given angles.
/*! 
  
  \param abs Output : absorption coefficient vector 
  \param ext  Input : Extinction Matrix
  \param pha  Input : Phase matrix
*/
void amp2abs(VectorView abs,ConstMatrixView ext,ConstTensor4View pha)
{
  Matrix Z11_mat = pha (Range(joker),Range(joker),0,0);
  Matrix Z21_mat = pha (Range(joker),Range(joker),1,0);
  Matrix Z31_mat = pha (Range(joker),Range(joker),2,0);
  Matrix Z41_mat = pha (Range(joker),Range(joker),3,0);
  Numeric low_limit_za=0.0;
  Numeric up_limit_za=180.0;
  Numeric no_steps_za=180.0;
  Numeric step_za= (up_limit_za-low_limit_za)/no_steps_za;
  Numeric low_limit_aa=0.0;
  Numeric up_limit_aa=360.0;
  Numeric no_steps_aa=360.0;
  Numeric step_aa= (up_limit_aa-low_limit_aa)/no_steps_aa;
  Z11_mat *= DEG2RAD;
  Z21_mat *= DEG2RAD;
  Z31_mat *= DEG2RAD;
  Z41_mat *= DEG2RAD;
  Numeric Z11_integrated;
  double_trapez(Z11_integrated,
		Z11_mat,
		low_limit_aa,
		low_limit_za,
		up_limit_aa,
		up_limit_za,
		step_aa,
		step_za);
  Numeric Z21_integrated;
  double_trapez(Z21_integrated,
		Z21_mat,
		low_limit_aa,
		low_limit_za,
		up_limit_aa,
		up_limit_za,
		step_aa,
		step_za);
  Numeric Z31_integrated;
  double_trapez(Z31_integrated,
		Z31_mat,
		low_limit_aa,
		low_limit_za,
		up_limit_aa,
		up_limit_za,
		step_aa,
		step_za);
  Numeric Z41_integrated;
  double_trapez(Z41_integrated,
		Z41_mat,
		low_limit_aa,
		low_limit_za,
		up_limit_aa,
		up_limit_za,
		step_aa,
		step_za);
  abs[0] = ext(0,0)-Z11_integrated;
  abs[1] = ext(1,1)-Z21_integrated;
  abs[2] = ext(2,2)-Z31_integrated;
  abs[3] = ext(3,3)-Z41_integrated;
}
//! Performs integration for a matrix
/*
  \param Integral Output : gives the value after integration
  \param Integrand Input :the expression to be integrated
  \param LowLimit1 Input :Lowerlimit of the inner integral
  \param LowLimit2 Input : Lower Limit of the outer integral
  \param UpLimit1 Input : Upper limit of the inner integral
  \param UpLimit2 Input : Upper limit of the outer integral
  \param h1 Input :  Step length for inner integral.
  \param h2 Input : Step length for inner integral.
*/
void double_trapez(Numeric &Integral,ConstMatrixView Integrand,Numeric &LowLimit1,Numeric &LowLimit2,Numeric &UpLimit1,Numeric &UpLimit2, Numeric &h1, Numeric &h2)
{
  Integral=0.0;
  //Numeric no_elem_i=(UpLimit2-LowLimit2)/h2 +1;
  Vector Integral1(Integrand.nrows());
  
  //cout<<LowLimit1<<" "<<LowLimit2<<" "<<UpLimit1<<" "<<UpLimit2<<"\n";
  cout<<"Number of columns  in Integral1"<<" "<<Integral1.nelem()<<"\n";
  for (Numeric i=LowLimit2;i<=UpLimit2;i=i+h2)
    {
      Integral1[i]=0.0;
      for (Numeric j=LowLimit1;j<=UpLimit1;j=j+h1)
	{
	  Integral1[i]=Integral1[i]+h1*Integrand(i,j);
	  
	}
      Integral1[i] = Integral1[i] - 0.5 * h1*(Integrand(i,LowLimit1)+Integrand(i,UpLimit1));
      Integral=Integral+h2*Integral1[i];
     
      
    }
  //cout << " §$%&//// "<< " " <<Integral1<<" "<<"\n" ;
  Integral = Integral - 0.5 *h2*(Integral1[LowLimit2]+Integral1[UpLimit2]);
  
  cout << " §$%&//// "<< " " <<Integral<<"\n" ;
  
  
}
//! Performs  integration for a vector
/*! 
  \param Integral Output : gives the value after integration
  \param Integrand  Input : the expression to be integrated
  \param LowLimit1  Input : Lowerlimit of the inner integral
  \param UpLimit1 Input : Upper limit of the inner integral
  \param h Input : Step length for inner integral.
*/
void single_trapez(Numeric &Integral,VectorView Integrand,Numeric &LowLimit1,Numeric &UpLimit1,Numeric &h)
{
    Integral=0.0;
  for (Numeric j=LowLimit1;j<=UpLimit1;j=j+h)
    {
     
      Integral =Integral +h*Integrand[j];
    }
  Integral = Integral - 0.5 *h*(Integrand[LowLimit1]+Integrand[UpLimit1]);
}
//! Extinction Coefficient Matrix for the particle 
/*! 
  
  This function sums up the extinction matries for all particle types weighted with particle number density
  \param ext_mat_part Output : physical extinction coefficient for the particles for given angles. 
  \param ext_mat_spt Input : extinction matrix for the scattering species
  \param pnd Input : particle number density givs the local concentration for all particles.
*/
//FIXME; change all matrix ext_mat_spt to Tensor3 ext_mat_spt(l,j,k)
void ext_mat_partCalc(MatrixView ext_mat_part,MatrixView ext_mat_spt,VectorView pnd)
{
  Index N_pt,j,k,l,N_i;
  for (j=0;j< N_i;j++)
    {
      for (k=0;k<N_i;k++)
	ext_mat_part(j,k)=0.0;// Initialisation
    }
  for (l=0;l< N_pt;l++)
    {
      for (j=0;j< N_i;j++)
	{
	  for (k=0;k<N_i;k++)
	    ext_mat_spt(j,k) *= pnd[l]; //multiplying with particle number density.
	  ext_mat_part(j,k) += ext_mat_spt(j,k); // summing up for all particle types.
	}
    }
}
