/*!
  \file   zeemanproperties.cc
  \author Nikolay Koulev, Oliver Lemke 
  \date   Mon Dec 06 15:54:00 2003
  
  \brief  This file has functions for calculating the propagation tensor in the case of Zeeman efect, only for oxygen at the moment.
 
*/
#include <iostream>
#include <cmath>
#include "arts.h"
#include "matpackIII.h"
#include "xml_io.h"
#include "complex.h"


extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT; 
const Numeric twoPI_c=2*PI/SPEED_OF_LIGHT;
const Complex complex_i = pow (Complex (-1., 0), 0.5);
 
void Zeeman (Vector& f_grid,
	     Tensor3& ext_mat_zee,
	     Matrix& abs_vec_zee,
	     Matrix& xi_mat,
	     Matrix& f_z_mat)
{
  
  Array<Complex> N_s(f_grid.nelem(), 0);
  
  
  //! Magnitude of Earth's magnetic field along LOS.
  Numeric B_field;
  
  //! Dummy value for the magnetic field.
  B_field = 50;

  //! Angle bertween the directionn of Earth's magnetic field and LOS.    
  Numeric phi;

  //! Dummy value for phi.
  phi = 127.83867/57.2957914331333;

  //! Pressure in [kPa]- just provisionary now. Later will bw replaced by p_grid.
  Numeric p;

  //! Test value for the pressure.
  p = 0.0044661608;

  //! Temperature in [K].
  Numeric T;

  //! Test value for the temperature.
  T = 217.359;

  //! Doppler broadened line width.
  Numeric gamma;
  
  //! Test value for gamma.
  gamma = pow(10.,-4.);

  //! Pressure broadened line width.
  Numeric agam;

  //! Test value for gamma.
  agam = pow(10.,-4.);

     
  //! This is the rotational quantum number denotation of a given Zeeman transition. It has a value of N for the N+ transitions and -N for the N- transitions.
  Numeric N_r;

  //! This the cenre frequency of the unsplit line.
  Numeric f_c;

  //! First intensity coeffcient for a Zeeman split line at given centre frequency.
  Numeric a1;

  //! Second intensity coeffcient for a Zeeman split line at given centre frequency.
  Numeric a2;

  //! Wavenumber.
  Numeric k_0;

  Matrix N;
  xml_read_from_file ("zeeman_intensity_coeff.xml", N);
  // for (int i=N.nrows()-1;i<N.nrows();i++)
  //{
      N_r = N(10,0); 
      cout << "N_r " << N_r << endl;
      f_c = N(10,1);
      cout.precision(10);
      cout << "f_c " << f_c << endl;
      a1 = N(10,2);
      a2 = N(10,3);

      f_z_mat.resize(2*abs(N_r)+1,3);
      xi_mat.resize(2*abs(N_r)+1,3);
	  //! This the vector consists of 3 possible values for the difference of the magnetic quantum number M for the Zeeman transitions.
	  
	  for (Index k=0; k<3;k++)
	    { Index DeltaM;
	      DeltaM = k-1;

	      //! This is the vector consisting of all 2N+1 values the magnetic quantum M for a given transition with rotational number N.
      
	      for (int j=0 ; j<(2*abs(N_r)+1);j++)
		{ 
		  Numeric M;
		  M = j-abs(N_r);

		 
		  //! This is a normalizing factor which gives the intensity of the individual components of the Zeeman split relative to the the intensity of the unsplit line.
		  Numeric xi;

		  //! This the Zeeman shift coefficients.
		  Numeric eta;
	      
		  if (N_r>0 && DeltaM==0)
		    {
 		      xi = 3*((abs(N_r)+1)*(abs(N_r)+1)-M*M)/((abs(N_r)+1)*(2*abs(N_r)+1)*(2*abs(N_r)+3));
		      xi_mat(j,k)= xi;
		      eta = M*(abs(N_r)-1)/(abs(N_r)*(abs(N_r)+1));
		      cout << "if1" << endl;
		    }
		  else if (N_r<0 && DeltaM==0)
		    {	
		      xi = 3*(abs(N_r)*abs(N_r)-M*M)/(abs(N_r)*(2*abs(N_r)-1)*(2*abs(N_r)+1));
		      xi_mat(j,k)= xi;
		      eta = M*(abs(N_r)+2)/(abs(N_r)*(abs(N_r)+1));
		      cout << "if2" << endl;
		    }
		  else if (N_r>0 && DeltaM==1)
		    {
		      xi = 3*(abs(N_r)+M+1)*(abs(N_r)+M+2)/(4*(abs(N_r)+1)*(2*abs(N_r)+1)*(2*abs(N_r)+3));
		      xi_mat(j,k)= xi;
		      eta = (M*(abs(N_r)-1)+abs(N_r))/(abs(N_r)*(abs(N_r)+1));
		      cout << "if3" << endl;
		    }
		  else if (N_r>0 && DeltaM==-1)
		    {	
		      xi = 3*(abs(N_r)-M+1)*(abs(N_r)-M+2)/(4*(abs(N_r)+1)*(2*abs(N_r)+1)*(2*abs(N_r)+3));
		      xi_mat(j,k)= xi;
		      eta = (M*(abs(N_r)-1)-abs(N_r))/(abs(N_r)*(abs(N_r)+1));
		      cout << "if4" << endl;
		    }
		  else if (N_r<0 && DeltaM==-1)
		    {
		      xi = 3*(abs(N_r)+M)*(abs(N_r)+M-1)/(4*abs(N_r)*(2*abs(N_r)-1)*(2*abs(N_r)+1));
		      xi_mat(j,k)= xi;
		      eta = (M*(abs(N_r)+2)+(abs(N_r)+1))/(abs(N_r)*(abs(N_r)+1));
		      cout << "if5" << endl;
		    }
		  else if (N_r<0 && DeltaM==1)
		    {
		      xi = 3*(abs(N_r)-M)*(abs(N_r)-M-1)/(4*abs(N_r)*(2*abs(N_r)-1)*(2*abs(N_r)+1));
		      xi_mat(j,k)= xi;
		      eta = (M*(abs(N_r)+2)-(abs(N_r)+1))/(abs(N_r)*(abs(N_r)+1));
		      cout << "if6" << endl;
		    }
		  else 
		    {
		      cout << "Error " << endl;
		    }
				
		  //! Centre frequency of the of the individual components of the Zeeman split.
		  Numeric f_z;
		  f_z = f_c + 28.03*pow(10.,-6.)*eta*B_field;
		  f_z_mat(j,k)=f_z;
		  cout.precision(10);
		  cout << "f_z " << f_z << endl;
		  cout << "eta " << eta << endl;
		  cout << "xi " << xi << endl;

		  //! Line strength.
		  Numeric S;
		  S = a1*p*pow(300./T,3)*exp(a2*(1-300./T));
		  
		  for (Index f_grid_index=0; f_grid_index < f_grid.nelem(); f_grid_index++)
		    {
		      //! Complex argument of the special line shape (CEF), used in the case of Zeeman splitting.
		      Complex zeta;
		      zeta = agam/gamma + complex_i*(f_grid[f_grid_index] - f_z)/gamma;
		      //cout << "zeta " << zeta << endl;
		      
		      
		      //! The special line shape, Complex Error Function.
		      Complex CEF;
		      CEF = (122.60793178 + 214.38238869*zeta + 181.92853309*pow(zeta,2) 
			     + 93.15558046*pow(zeta,3) + 30.18014220*pow(zeta,4) + 5.91262621*pow(zeta,5) 
			     + 0.56418958*pow(zeta,6))/(122.60793178 + 352.73062511*zeta + 457.33447878*pow(zeta,2) 
							+ 348.70391772*pow(zeta,3) + 170.35400182*pow(zeta,4) + 53.99290691*pow(zeta,5) 
			   + 10.47985711*pow(zeta,6) + 1.00000000*pow(zeta,7));
		      
		      //!The complex refractive index. 
		      N_s[f_grid_index] += S*xi*CEF;
		      
		    }
		  
		}
	      
	    }
	  //  }
  
  
  
  Matrix I;
  I.resize(2,2);
  I(0,0) = 1;
  I(0,1) = 0;
  I(1,0) = 0;
  I(1,1) = 1;
  
  //! This the propagation tensor at given frequency.
  Array<Complex> G_s(f_grid.nelem()*2*2, 0);
  
  for (int k=0; k<3;k++)
    {
      Numeric DeltaM;
      DeltaM=k-1;
      
      //! This is a matrix accounting for the contributions of the 3 different polarizations of the components of the Zeeman split due to 3 different values of DeltaM.
      Matrix P;
      P.resize(2,2);
      
      if (DeltaM==0)
	{
	  P(0,0) =  0.5*sin(phi)*sin(phi);
	  P(0,1) = -0.5*sin(phi)*sin(phi);
	  P(1,0) =  0.5*sin(phi)*sin(phi);
	  P(1,1) = -0.5*sin(phi)*sin(phi);
	}
      else if (DeltaM==1)
	{
	  P(0,0) = 0.5*(1+cos(phi))*(1+cos(phi));
	  P(0,1) = 0.5*sin(phi)*sin(phi);
	  P(1,0) = 0.5*sin(phi)*sin(phi);
	  P(1,1) = 0.5*(1-cos(phi))*(1-cos(phi));;
	}
      else if (DeltaM==-1)
	{
	  P(0,0) = 0.5*(1-cos(phi))*(1-cos(phi));
	  P(0,1) = 0.5*sin(phi)*sin(phi);
	  P(1,0) = 0.5*sin(phi)*sin(phi);
	  P(1,1) = 0.5*(1+cos(phi))*(1+cos(phi));
	}

      
      for (Index f_grid_index=0 ;f_grid_index< f_grid.nelem(); 
	   f_grid_index++)
	{
	  
	  k_0 = twoPI_c*f_grid[f_grid_index];

	  G_s[f_grid_index*4+0] += complex_i*k_0*P(0,0)*N_s[f_grid_index];
	  G_s[f_grid_index*4+1] += P(0,1)*complex_i*k_0*N_s[f_grid_index];
	  G_s[f_grid_index*4+2] += P(1,0)*complex_i*k_0*N_s[f_grid_index];
	  G_s[f_grid_index*4+3] += complex_i*k_0*P(1,1)*N_s[f_grid_index];
	}
    }
     

  
  //! This is resize of the tensor consisting of the Zeeman Extinction matrix and the frequency grid.
    ext_mat_zee.resize(f_grid.nelem(),4,4);
  
  //! This  is resize of tthe matrix consisting of the Zeeman absorption vector and the frequency grid.
    abs_vec_zee.resize(f_grid.nelem(),4);

  for (Index f_grid_index=0; f_grid_index<f_grid.nelem();f_grid_index++)
    {
      k_0=twoPI_c*f_grid[f_grid_index];

     G_s[f_grid_index*4+0]  = complex_i*k_0*(1. + G_s[f_grid_index*4+0]);
     G_s[f_grid_index*4+1] *= complex_i*k_0;
     G_s[f_grid_index*4+2] *= complex_i*k_0;
     G_s[f_grid_index*4+3]  = complex_i*k_0*(1. + G_s[f_grid_index*4+3]);
    
     //! First row of the extinction matrix.
     ext_mat_zee(f_grid_index,0,0) = 0.5*(2.*real(G_s[f_grid_index*4+0]) + 2.*real(G_s[f_grid_index*4+3]));
     ext_mat_zee(f_grid_index,0,1) = 0.5*(2.*real(G_s[f_grid_index*4+1]) + 2.*real(G_s[f_grid_index*4+2]));
     ext_mat_zee(f_grid_index,0,2) = 0;
     ext_mat_zee(f_grid_index,0,3) = 0.5*(2.*real(G_s[f_grid_index*4+0]) - 2.*real(G_s[f_grid_index*4+3]));

     //! Second row of the extinction matrix.
     ext_mat_zee(f_grid_index,1,0) = 0.5*(2.*real(G_s[f_grid_index*4+1]) + 2.*real(G_s[f_grid_index*4+2]));
     ext_mat_zee(f_grid_index,1,1) = 0.5*(2.*real(G_s[f_grid_index*4+0]) + 2.*real(G_s[f_grid_index*4+3]));
     ext_mat_zee(f_grid_index,1,2) = 0.5*(2.*imag(G_s[f_grid_index*4+0]) - 2.*imag(G_s[f_grid_index*4+3]));
     ext_mat_zee(f_grid_index,1,3) = 0;

     //! Third row of the extinction matrix.
     ext_mat_zee(f_grid_index,2,0) = 0;
     ext_mat_zee(f_grid_index,2,1) = 0.5*(2.*imag(G_s[f_grid_index*4+0]) - 2.*imag(G_s[f_grid_index*4+3]));
     ext_mat_zee(f_grid_index,2,2) = 0.5*(2.*real(G_s[f_grid_index*4+0]) + 2.*real(G_s[f_grid_index*4+3]));
     ext_mat_zee(f_grid_index,2,3) = 0.5*(2.*imag(G_s[f_grid_index*4+1]) + 2.*imag(G_s[f_grid_index*4+2]));

     //! Fourth row of the extinction matrix.
     ext_mat_zee(f_grid_index,3,0) = 0.5*(2.*real(G_s[f_grid_index*4+0]) - 2.*real(G_s[f_grid_index*4+3])); 
     ext_mat_zee(f_grid_index,3,1) = 0;
     ext_mat_zee(f_grid_index,3,2) = 0.5*(2.*imag(G_s[f_grid_index*4+1]) + 2.*imag(G_s[f_grid_index*4+2]));
     ext_mat_zee(f_grid_index,3,3) = 0.5*(2.*real(G_s[f_grid_index*4+0]) + 2.*real(G_s[f_grid_index*4+3]));  
     
     //! Absorption vector components.
     abs_vec_zee(f_grid_index,0) = 2.*real(G_s[f_grid_index*4+0]) + 2.*real(G_s[f_grid_index*4+3]);
     abs_vec_zee(f_grid_index,1) = 2.*real(G_s[f_grid_index*4+1]) + 2.*real(G_s[f_grid_index*4+2]);
     abs_vec_zee(f_grid_index,2) = 0;
     abs_vec_zee(f_grid_index,3) = 2.*real(G_s[f_grid_index*4+0]) - 2.*real(G_s[f_grid_index*4+3]);
    } 


  
}










