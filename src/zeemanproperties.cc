/*! 
  \file   zeemanproperties.cc
  \author Nikolay Koulev, Oliver Lemke 
  \date   Mon Dec 06 15:54:00 2003
  
  \brief  This file has functions for calculating the propagation tensor in the case of Zeeman effect, only for oxygen at the moment.
 
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
const Complex complex_i = pow(Complex (-1., 0), 0.5);
 
void Zeeman (//Input & Output
	     Vector& f_grid, // Frequency grid.
	     
	     //Output
	     Tensor3& ext_mat_zee, // Tensor3 of the Extinction Matrix for the specified frequency grid in [1/m]. 
	     Matrix& abs_vec_zee,  // Matrix of the Absorption Vector for the specified frequency grid in [1/m].
	     Matrix& xi_mat, // Matrix of the relative intensities of the Zeeman components, dimensionless.
	     Matrix& f_z_mat, // Matrix of the frequencies of the Zeeman components in [Hz].
	     
	     //Input
	     Numeric& N_r, // Rotational quantum number indentification of the Zeeman line.
	     Numeric& BN_r, // Number of the split components for a given polarization.
	     Numeric& AN_r, // Maximum value of the magnetic quantum number for a given polarization.
	     Numeric& f_c, // Central frequency of the unsplit line [Hz].
	     Numeric& a1, // First intensity coeficient of the unsplit line in [Hz/Pa].
	     Numeric& a2, // Second intensity coeficient of the unsplit line, dimensionless.
	     Numeric& a3) // Third intensity coeficient of the unsplit line in [Hz/Pa].
	     
  

{
  
  //----------------------------------------------------  
  // Provisionary hardwired values for certain variables.
  //----------------------------------------------------
  
  // Magnitude of Earth's magnetic field along LOS in [T].
  Numeric B_field;
  
  // Dummy value for the magnetic field in [T].
  B_field = 3.5e-5;

  // Angle between the directions of Earth's magnetic field and LOS.    
  Numeric phi;

  // Dummy value for phi.
  phi = 90./57.2957914331333;
  

  // Pressure in [Pa]- just provisionary now. 
  //Later will bw replaced by p_grid.
  Numeric p;

  // Test value for the pressure in [Pa] at 95 km 'midlatitude-summer'.
  p = 6.200000e-02;

  // Temperature in [K].
  Numeric T;

  // Test value for the temperature in [K] at 95 km for 'midlatitude-summer'.
  T = 1.783000e+02; 
  //-------------------------------------------------

    


    //------------------------------------------------------------
    // Definitions of the line intensity and broadening parameters,
    // wavenumber and frequency shift.
    //--------------------------------

    // Doppler broadening line width in [Hz].
    Numeric gamma;
    
    
    // Pressure broadening line width [Hz].
    Numeric agam;
    
    
    // Wavenumber.
    Numeric k_0;

    // Intensity factor, or relative intensity, of the 
    // individual components of the Zeeman split relative 
    // to the the intensity of the unsplit line.
    Numeric xi = 0.;
    
    
    // This the Zeeman frequency shift of the individual 
    // components of the Zeeman split.
    Numeric eta = 0.;
    //--------------------------------------------------





    //---------------------------------------------------------------
    //  Definitions & resizements of the tensor and array quantities.
    //---------------------------------------------------------------
    
    // Magnetic susceptibility tensor at a given frequency. 
    // It's in the form of an array because complex matrices
    // or tensors cannot be handled by ARTS.
    Array<Complex> Chi(f_grid.nelem()*2*2,0);
    
    // Complex refractive index.
    Array<Complex> N_s(f_grid.nelem(),0);
    
    // Propagation tensor at a given frequency. It's in the form of an array 
    // because complex matrices or tensors cannot be handled by ARTS. 
    Array<Complex> G_s(f_grid.nelem()*2*2);
      
    // Resize of the tensor consisting of the Zeeman 
    // Extinction matrix and the frequency grid.
    ext_mat_zee.resize(f_grid.nelem(),4,4);
    
    // Resize of the matrix consisting of the Z
    //Z eeman absorption vector and the frequency grid.
    abs_vec_zee.resize(f_grid.nelem(),4);      

    //---------------------------------------------------------------
    
    
    
    
    
    // Matrix variable for the content of the zeeman_intensity_coeff.xml.
    Matrix N;
    
    // Reading in content from zeeman_intensity_coeff.xml.
    xml_read_from_file ("zeeman_intensity_coeff.xml", N);
    
    
    
     if (N_r>0)
      { 
	
	// Assigning a value of the total angular momentum N of the lower  
	// level for a N+ Zeeman transition (N -> N+1).
	AN_r= abs(N_r);
	
	// Assigning a value of the number of the split components for each 
	// of the three possible polarizations of a N+ Zeeman transition.
	BN_r= 2*abs(N_r)+1;
	
	//! Resize of the frequency shift matrix 
	// for a N+ Zeeman transition.
	f_z_mat.resize(2*Index(abs(N_r))+1,3);
	
	//! Resize of the relative intensity matrix 
	// for a N+ Zeeman transition.
	xi_mat.resize(2*Index(abs(N_r))+1,3);
	
      }
    
    else if (N_r<0)
      
      {
	// Assigning a value of the total angular momentum N of the lower 
	// level for a N- Zeeman transition (N -> N-1).
	AN_r= abs(N_r)-1;
	
	// Assigning a value of the number of the split components for each
	// of the three possible polarizations of a N- Zeeman transition.
	BN_r= 2*abs(N_r)-1;
	
	
	// Resize of the frequency shift matrix 
	// for a N- Zeeman transition.
	f_z_mat.resize(2*Index(abs(N_r))-1,3);
	
	// Resize of the relative intensity matrix 
	// for a N- Zeeman transition.
	xi_mat.resize(2*Index(abs(N_r))-1,3);
      }
    
    
    
    
    
    
    
    // Starting the loop over the frequency grid.
    for (Index f_grid_index=0; f_grid_index < f_grid.nelem(); f_grid_index++)
      {
	
	//  Wavenumber claculation.
	k_0=twoPI_c*f_grid[f_grid_index];
	
	
	// Starting the loop over the 3 possible values of the change of the 
	// magnetic quantum number M upon a Zeeman transition.	  
	for (Index k=0;k<3;k++)
	  { 
	    Index DeltaM;
	    DeltaM = k-1;
	    
	    
	    
	    // Starting the loop over 2N+1 or 2N-1 values of
	    // the magnetic quantum number M for a given transition 
	    for (int j=0 ; j<BN_r;j++)
	      { 
		Numeric M;
		M = j-AN_r;
		
		
		
		// Pi polarization for a N+ Zeeman transition.
		if (N_r>0 && DeltaM==0) //
		  {
		    // Matrix representation of the relative intensity value 
		    xi_mat(j,k) = 3*((abs(N_r)+1)*(abs(N_r)+1)-M*M)
		      /((abs(N_r)+1)*(2*abs(N_r)+1)*(2*abs(N_r)+3));
		    
		   		   		    
		    // Zeeman frequency shift value
		    eta = M*(abs(N_r)-1)/(abs(N_r)*(abs(N_r)+1));
		   
		
		  }
		
		
		// Pi polarization for a N- Zeeman transition.
		else if (N_r<0 && DeltaM==0)
		  
		  {
		    // Matrix representation of the relative intensity values 
		    xi_mat(j,k)= 3*(abs(N_r)*abs(N_r)-M*M)
		      /(abs(N_r)*(2*abs(N_r)-1)*(2*abs(N_r)+1));
		    
		    
		    // Zeeman frequency shift values
		    eta = M*(abs(N_r)+2)/(abs(N_r)*(abs(N_r)+1));
		  		    
		  }
		
		// Sigma+ polarization for a N+ Zeeman transition.
		else if (N_r>0 && DeltaM==1)
		  
		  {
		    // Matrix representation of the relative intensity
		    xi_mat(j,k) = 3*(abs(N_r)+M+1)*(abs(N_r)+M+2)/
		      (4*(abs(N_r)+1)*(2*abs(N_r)+1)*(2*abs(N_r)+3));
		    
		    // Zeeman frequency shift value 
		    eta = (M*(abs(N_r)-1)+abs(N_r))/(abs(N_r)*(abs(N_r)+1));
						
		  }
		
		// Sigma- polarization for a N+ Zeeman transition.
		else if (N_r>0 && DeltaM==-1)
		  
		  {	
		    // Matrix representation of the relative intensity values
		    xi_mat(j,k) = 3*(abs(N_r)-M+1)*(abs(N_r)-M+2)/
		      (4*(abs(N_r)+1)*(2*abs(N_r)+1)*(2*abs(N_r)+3));
		   		    
		    // Zeeman frequency shift value
		    eta = (M*(abs(N_r)-1)-abs(N_r))/(abs(N_r)*(abs(N_r)+1));
						
		  }
		
		//! Sigma+ polarization for a N- Zeeman transition.
		else if (N_r<0 && DeltaM==1)
		  
		  {
		    // Matrix representation of the relative intensity values
		    xi_mat(j,k) = 3*(abs(N_r)-M)*(abs(N_r)-M-1)/
		      (4*abs(N_r)*(2*abs(N_r)-1)*(2*abs(N_r)+1));
		    
		    // Zeeman frequency shift value 
		    eta = (M*(abs(N_r)+2)+(abs(N_r)+1))/(abs(N_r)*(abs(N_r)+1));
		    
		  }
		
		// Sigma- polarization for a N- Zeeman transition.
		else if (N_r<0 && DeltaM==-1)
		  
		  {
		    // Matrix representation of the relative intensity values
		    xi_mat(j,k)= 3*(abs(N_r)+M)*(abs(N_r)+M-1)/
		      (4*abs(N_r)*(2*abs(N_r)-1)*(2*abs(N_r)+1));
		    
		   
		    
		    // Zeeman frequency shift value 
		    eta = (M*(abs(N_r)+2)-(abs(N_r)+1))/(abs(N_r)*(abs(N_r)+1));
		    
		  }
		
		else 
		  {
		    cout << "Error N_r or DeltaM" << endl;
		  }
		
		
		
		// Center frequency of the of the individual 
		// components of the Zeeman split.
		Numeric f_z;
		f_z = f_c + 28.03*1e+9*eta*B_field;
		f_z_mat(j,k)=f_z;
		
	
		
		// Line strength of the unsplit line in [Hz].
		Numeric S;
		Numeric theta = 300./T; 
		S = a1*1e-3*p*1e3*pow(theta,3)*exp(a2*(1-theta));
		
		// Numeric expressions for both the pressure (agam)
		// and Doppler (gamma) broadening line width 
		agam = a3*1e3* p * pow(theta,0.85);
		gamma = (1.096/sqrt(theta))*f_z*1e3;
		
		
		// Complex argument of the special line shape (CEF), 
		// used in the case of Zeeman splitting.
                Complex zeta;
		zeta = agam/gamma + complex_i*(f_grid[f_grid_index] - f_z)/gamma;
		
		
		
		// The special line shape, Complex Error Function.
		Complex CEF;
		CEF = (1/(sqrt(PI)*gamma))*(122.60793178*pow(zeta,0) 

		     + 214.38238869*pow(zeta,1) + 181.92853309*pow(zeta,2) 

		     + 93.15558046*pow(zeta,3) + 30.18014220*pow(zeta,4)

		     + 5.91262621*pow(zeta,5) + 0.56418958*pow(zeta,6))

		  /(122.60793178*pow(zeta,0) + 352.73062511*pow(zeta,1)
 
                     + 457.33447878*pow(zeta,2) + 348.70391772*pow(zeta,3) 

                     + 170.35400182*pow(zeta,4) + 53.99290691*pow(zeta,5)
 
		     + 10.47985711*pow(zeta,6) + 1.00000000*pow(zeta,7));
	      
		  
	      
	      
		  // The complex refractive index. 
		  N_s[f_grid_index] += S*xi*(-1)*complex_i*CEF;
	      
		 
	      
		  //!With Lorenzian.
		  //N_s[f_grid_index] = S*xi/(f_z - f_grid[f_grid_index] + 
		  //complex_i*agam);
	      
		  //!With Lorenzian - Whiting's approximation.
		  //N_s[f_grid_index] = S*xi/(f_z - f_grid[f_grid_index] + 
		  //complex_i*(agam/2 + pow(agam/2*agam/2 + gamma*gamma,0.5)));
	      
	      
	      
		}
	  
	  
	  
	  
	  
	  
	  
	  
	      // Polarization  matrix accounting for the contributions of 
	      // the 3 different polarizations of the components of the Zeeman 
	      // split due to 3 different values of DeltaM.
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
	  
	  
	     
	  

	      // Calculating the halved value of magnetic susceptibility tensor 
	      // at given frequency. It's in the form of an array because 
	      // complex matrices or tensors cannot be handled by ARTS.		
	      Chi[f_grid_index*4+0] += P(0,0)*N_s[f_grid_index];
	      Chi[f_grid_index*4+1] += P(0,1)*N_s[f_grid_index];
	      Chi[f_grid_index*4+2] += P(1,0)*N_s[f_grid_index];
	      Chi[f_grid_index*4+3] += P(1,1)*N_s[f_grid_index];
	  
	  	  
	    }
      
      
      
      
	  
      
      
	  // Calculating the propagation tensor at a given frequency. 
	  // It's in the form of an array because complex matrices or 
	  // tensors cannot be handled by ARTS.
	  G_s[f_grid_index*4+0] = complex_i*k_0*(1. + Chi[f_grid_index*4+0]);
	  G_s[f_grid_index*4+1] = complex_i*k_0*Chi[f_grid_index*4+1];
	  G_s[f_grid_index*4+2] = complex_i*k_0*Chi[f_grid_index*4+2];
	  G_s[f_grid_index*4+3] = complex_i*k_0*(1. + Chi[f_grid_index*4+3]);
      
      
      
	  // First row of the extinction matrix at a given frequency.
	  ext_mat_zee(f_grid_index,0,0) = 
	    0.5*(2.*real(G_s[f_grid_index*4+0]) + 2.*real(G_s[f_grid_index*4+3]));

	  ext_mat_zee(f_grid_index,0,1) = 
	  - 0.5*(2.*real(G_s[f_grid_index*4+3]) - 2.*real(G_s[f_grid_index*4+0]));

	  ext_mat_zee(f_grid_index,0,2) = 
	    0.5*(2.*real(G_s[f_grid_index*4+1]) + 2.*real(G_s[f_grid_index*4+2]));

	  ext_mat_zee(f_grid_index,0,3) = 0.0;
      
	  // Second row of the extinction matrix at a given frequency.
	  ext_mat_zee(f_grid_index,1,0) = 
	  - 0.5*(2.*real(G_s[f_grid_index*4+3]) - 2.*real(G_s[f_grid_index*4+0]));

	  ext_mat_zee(f_grid_index,1,1) = 0.5*(2.*real(G_s[f_grid_index*4+0]) 
					     + 2.*real(G_s[f_grid_index*4+3]));
	  ext_mat_zee(f_grid_index,1,2) = 0.0;
	  ext_mat_zee(f_grid_index,1,3) = 0.5*(2.*imag(G_s[f_grid_index*4+1]) + 
					       2.*imag(G_s[f_grid_index*4+2]));
      
	  // Third row of the extinction matrix.
	  ext_mat_zee(f_grid_index,2,0) = 0.5*(2.*real(G_s[f_grid_index*4+1]) + 
					       2.*real(G_s[f_grid_index*4+2]));
	  ext_mat_zee(f_grid_index,2,1) = 0.0;
	  ext_mat_zee(f_grid_index,2,2) = 0.5*(2.*real(G_s[f_grid_index*4+0]) + 
					       2.*real(G_s[f_grid_index*4+3]));
	  ext_mat_zee(f_grid_index,2,3) = -0.5*(2.*imag(G_s[f_grid_index*4+3]) - 
						2.*imag(G_s[f_grid_index*4+0]));
      
	  // Fourth row of the extinction matrix at a given frequency.
	  ext_mat_zee(f_grid_index,3,0) = 0.0; 
	  ext_mat_zee(f_grid_index,3,1) = -0.5*(2.*imag(G_s[f_grid_index*4+1]) + 
						2.*imag(G_s[f_grid_index*4+2]));
	  ext_mat_zee(f_grid_index,3,2) = 0.5*(2.*imag(G_s[f_grid_index*4+3]) - 
					       2.*imag(G_s[f_grid_index*4+0]));
	  ext_mat_zee(f_grid_index,3,3) = 0.5*(2.*real(G_s[f_grid_index*4+0]) + 
					       2.*real(G_s[f_grid_index*4+3]));
      
	  // Absorption vector components at a given frequency.
	  abs_vec_zee(f_grid_index,0) = 0.5*(2.*real(G_s[f_grid_index*4+0]) + 
					     2.*real(G_s[f_grid_index*4+3]));
	  abs_vec_zee(f_grid_index,1) = -0.5*(2.*real(G_s[f_grid_index*4+3]) - 
					      2.*real(G_s[f_grid_index*4+0]));
	  abs_vec_zee(f_grid_index,2) = 0.5*(2.*real(G_s[f_grid_index*4+1]) + 
					     2.*real(G_s[f_grid_index*4+2]));
	  abs_vec_zee(f_grid_index,3) = 0.0;
      
      
	}
  
  
}










