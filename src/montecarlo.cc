
/* Copyright (C) 2003 Cory Davis <cory@met.ed.ac.uk>
                            
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */



/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   montecarlo.cc
  \author Cory Davis <cory@met.ed.ac.uk>
  \date   2003-06-19 

  \brief  functions used by ScatteringMonteCarlo
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "montecarlo.h"

#ifdef HAVE_SSTREAM
#include <sstream>
#else
#include "sstream.h"
#endif



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/



//! Cloudbox_ppath_rteCalc

/*!
   This function was written to eliminate some repeated code in 
   ScatteringMonteCarlo. Basically this function performs internal and external 
   ppath calculations, calculates incoming radiances, and provides scattering 
   properties and other derived data that is required in the monte carlo 
   simulation for a given propagation direction.


   Parameters not already described elsewhere:
   
   \param    record_ppathcloud   A flag that determines whether ppaths within the 
                                 cloud box are recorded for later plotting in 
				 MATLAB.
   \param    record_ppath        Similar, but for ppaths determining incoming 
                                 radiances at the cloudbox boundaries.  For 
				 normal use both of these parameterss should be 
                                 set to 0.
   \param    photon_number       used to label ppath files
   \param    scattering_order    used to label ppath files



   \author Cory Davis
   \date   2003-06-20
*/

void Cloudbox_ppath_rteCalc(
			     Ppath&                ppathcloud,
			     Ppath&                ppath,
			     Ppath&                ppath_step,
			     Vector&               rte_pos,
			     Vector&               rte_los,
			     Vector&               cum_l_step,
			     ArrayOfMatrix&        TArray,
			     ArrayOfMatrix&        ext_matArray,
			     ArrayOfVector&        abs_vecArray,
			     Vector&               t_ppath,
			     Vector&               scat_za_grid,
			     Vector&               scat_aa_grid,
			     Tensor3&              ext_mat,
			     Matrix&               abs_vec,
			     Numeric&              rte_pressure,
			     Numeric&              rte_temperature,
			     Vector&               rte_vmr_list,
			     Matrix&               i_rte,
			     GridPos&              rte_gp_p,
			     GridPos&              rte_gp_lat,
			     GridPos&              rte_gp_lon,
			     Matrix&               i_space,
			     Matrix&               ground_emission,
			     Matrix&               ground_los, 
			     Tensor4&              ground_refl_coeffs,
			     Index&                f_index,
			     Index&                scat_za_index,
			     Index&                scat_aa_index,
			     Tensor3&              ext_mat_spt,
			     Matrix&               abs_vec_spt,
			     Matrix&               pnd_ppath,
			     const Agenda&         ppath_step_agenda,
			     const Index&          atmosphere_dim,
			     const Vector&         p_grid,
			     const Vector&         lat_grid,
			     const Vector&         lon_grid,
			     const Tensor3&        z_field,
			     const Matrix&         r_geoid,
			     const Matrix&         z_ground,
			     const ArrayOfIndex&   cloudbox_limits,
			     const Index&          record_ppathcloud,
			     const Index&          record_ppath,
			     const Agenda&         opt_prop_gas_agenda,
			     const Agenda&         spt_calc_agenda,
			     const Agenda&         scalar_gas_absorption_agenda,
			     const Index&          stokes_dim,
			     const Tensor3&        t_field,
			     const Tensor4&        vmr_field,
			     const Agenda&         rte_agenda,
			     const Agenda&         i_space_agenda,
			     const Agenda&         ground_refl_agenda,
			     const Vector&         f_grid,
			     const Index&          photon_number,
			     const Index&          scattering_order,
			     const Tensor4&        pnd_field)

{
  // Assign dummies for variables associated with sensor.
  Vector   mblock_za_grid_dummy(1);
           mblock_za_grid_dummy[0] = 0;
  Vector   mblock_aa_grid_dummy(0), sensor_rot_dummy(0);
  Matrix   sensor_pol_dummy;
  Index    antenna_dim_dummy = 1; 
  Sparse   sensor_response_dummy;

  // Dummy for measurement vector
  Vector   y_dummy(0);

  const Index cloudbox_on_dummy=0;
  Matrix sensor_pos(1,3);
  Matrix sensor_los(1,2);
  Tensor7 scat_i_p_dummy;
  Tensor7 scat_i_lat_dummy;
  Tensor7 scat_i_lon_dummy;
 
  //  cout << "Cloudbox_ppathCalc\n";
  Cloudbox_ppathCalc(ppathcloud,ppath_step,ppath_step_agenda,atmosphere_dim,
		     p_grid,lat_grid,lon_grid,z_field,r_geoid,z_ground,
		     cloudbox_limits, rte_pos,rte_los);
  if (record_ppathcloud)
    {
      //Record ppathcloud.  This is useful for debugging and educational 
      //purposes.  It would be completely daft to leave this on in a 
      //real calculation
      ostringstream filename;
      filename <<"ppathcloud" << photon_number <<"_"<<scattering_order;
      String longfilename;
      filename_xml(longfilename,filename.str());
      xml_write_to_file(longfilename, ppathcloud, FILE_TYPE_ASCII);
    }

  cum_l_stepCalc(cum_l_step,ppathcloud);
  
  
  //Calculate array of transmittance matrices
  TArrayCalc(TArray, ext_matArray, abs_vecArray, t_ppath, scat_za_grid, 
	     scat_aa_grid, ext_mat, abs_vec, rte_pressure, rte_temperature, 
	     rte_vmr_list, scat_za_index, scat_aa_index, ext_mat_spt, 
	     abs_vec_spt, pnd_ppath, ppathcloud, opt_prop_gas_agenda, 
	     spt_calc_agenda, scalar_gas_absorption_agenda, stokes_dim, 
	     p_grid, lat_grid, lon_grid, t_field, vmr_field, atmosphere_dim,
	     pnd_field);
  //Calculate contribution from the boundary of the cloud box
  //changed to dummy_rte_pos to see if rte_pos was causing assertion failure at ppath.cc:1880
  //it appears that this was not the case
  Vector dummy_rte_pos = rte_pos;
  Vector dummy_rte_los = rte_los;
  rte_posShift(dummy_rte_pos,dummy_rte_los,rte_gp_p, rte_gp_lat,
	    rte_gp_lon,ppathcloud, atmosphere_dim);
  sensor_pos(0,joker)=dummy_rte_pos;
  sensor_los(0,joker)=dummy_rte_los;
  //call rte_calc without input checking, sensor stuff, or verbosity
  rte_calc( y_dummy, ppath, ppath_step, i_rte, rte_pos, rte_los, rte_gp_p, 
	    rte_gp_lat,
	    rte_gp_lon,i_space, ground_emission, ground_los, ground_refl_coeffs,
	    ppath_step_agenda, rte_agenda, i_space_agenda, ground_refl_agenda, 
	    atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, t_field, 
	    r_geoid, z_ground, cloudbox_on_dummy, cloudbox_limits, 
	    scat_i_p_dummy,scat_i_lat_dummy, scat_i_lon_dummy, scat_za_grid,
	    scat_aa_grid, sensor_response_dummy, sensor_pos,sensor_los,
            sensor_pol_dummy, sensor_rot_dummy,
            f_grid,stokes_dim, antenna_dim_dummy, mblock_za_grid_dummy,
	    mblock_aa_grid_dummy, false, false, true);
  
  for (Index i = 0;i<stokes_dim;i++){assert(!isnan(i_rte(0,i)));}
  
  if (record_ppath)
    {
      //Record ppathcloud.  This is useful for debugging and educational purposes.  It would
      //be completely daft to leave this on in a real calculation
      ostringstream filename;
      filename <<"ppath" << photon_number <<"_"<<scattering_order;
      String longfilename;
      filename_xml(longfilename,filename.str());
      xml_write_to_file(longfilename, ppath, FILE_TYPE_ASCII);
    }
  
  
  f_index=0;//For some strange reason f_index is set to -1 in RteStandard
}



//! cloudbox_ppath_start_stepping

/*!
   This function was derived from Patrick's 'ppath_start_stepping' function.  
   It has been adapted for use within the cloud box and is intended for a 
   3D atmosphere only

   \param   ppath             Output: A Ppath structure.
   \param   atmosphere_dim    The atmospheric dimensionality.
   \param   p_grid            The pressure grid.
   \param   lat_grid          The latitude grid.
   \param   lon_grid          The longitude grid.
   \param   z_field           The field of geometrical altitudes.
   \param   r_geoid           The geoid radius.
   \param   z_ground          Ground altitude.
   \param   cloudbox_on       Flag to activate the cloud box.
   \param   cloudbox_limits   Index limits of the cloud box.
   \param   rte_pos             The position of the sensor.
   \param   rte_los             The line-of-sight of the sensor.

   \author Cory Davis (derived from ppath_start_stepping (Patrick Eriksson))
   \date   2003-06-19
*/
void cloudbox_ppath_start_stepping(
              Ppath&          ppath,
        const Index&          atmosphere_dim,
        ConstVectorView       p_grid,
        ConstVectorView       lat_grid,
        ConstVectorView       lon_grid,
        ConstTensor3View      z_field,
        ConstMatrixView       r_geoid,
        ConstMatrixView       z_ground,
        ConstVectorView       rte_pos,
        ConstVectorView       rte_los )
{
  // This function contains no checks or asserts as it is only a sub-function
  // to ppathCalc where the input is checked carefully.

  // Allocate the ppath structure
  ppath_init_structure(  ppath, atmosphere_dim, 1 );

  // Number of pressure levels
  const Index np = p_grid.nelem();

  // The different atmospheric dimensionalities are handled seperately

  if( atmosphere_dim == 3 )
  {

      // Is the sensor inside the latitude and longitude ranges of the
      // model atmosphere, and below the top of the atmosphere? If
      // yes, is_inside = true. Store geoid and ground radii, grid
      // position and interpolation weights for later use.
      //
      double   rv_geoid=-1, rv_ground=-1;  // -1 to avoid compiler warnings
      GridPos   gp_lat, gp_lon;
      Vector    itw(4);
      
      
	    gridpos( gp_lat, lat_grid, rte_pos[1] );
          gridpos( gp_lon, lon_grid, rte_pos[2] );
          interpweights( itw, gp_lat, gp_lon );

          rv_geoid  = interp( itw, r_geoid, gp_lat, gp_lon );
          rv_ground = rv_geoid + interp( itw, z_ground, gp_lat, gp_lon );

//          out2 << "  sensor altitude        : " << (rte_pos[0]-rv_geoid)/1e3 
//               << " km\n";
//
      // If downwards, calculate geometrical tangent position. If the tangent
      // point is inside the covered latitude range, calculate also the 
      // geometrical altitude of the tangent point and the top of atmosphere.
      //
//      Array<double>  geom_tan_pos(0);
//      double geom_tan_z=-2, geom_tan_atmtop=-1;  // OK values if the variables

      // Put sensor position and LOS in ppath as first guess
      ppath.pos(0,0) = rte_pos[0];
      ppath.pos(0,1) = rte_pos[1];
      ppath.pos(0,2) = rte_pos[2];
      ppath.los(0,0) = rte_los[0];
      ppath.los(0,1) = rte_los[1];


     // Geometrical altitude
      ppath.z[0] = ppath.pos(0,0) - rv_geoid;

     // Use below the values in ppath (instead of rte_pos and rte_los) as 
     // they can be modified on the way.
     
     // Grid positions
      ppath.gp_lat[0].idx   = gp_lat.idx;
      ppath.gp_lat[0].fd[0] = gp_lat.fd[0];
      ppath.gp_lat[0].fd[1] = gp_lat.fd[1];
      ppath.gp_lon[0].idx   = gp_lon.idx;
      ppath.gp_lon[0].fd[0] = gp_lon.fd[0];
      ppath.gp_lon[0].fd[1] = gp_lon.fd[1];

      // Create a vector with the geometrical altitude of the pressure 
      // surfaces for the sensor latitude and use it to get ppath.gp_p.
      Vector z_grid(np);
      z_at_latlon( z_grid, p_grid, lat_grid, lon_grid, z_field, 
                                                          gp_lat, gp_lon );
      gridpos( ppath.gp_p, z_grid, ppath.z );


      // Handle possible numerical problems for grid positions
      gridpos_check_fd( ppath.gp_p[0] );
      gridpos_check_fd( ppath.gp_lat[0] );
      gridpos_check_fd( ppath.gp_lon[0] );

    }  // End 3D
}


//! cum_l_stepCalc

/*!
   Returns a vector of cumulative pathlengths for a given ppath by adding ppath.l_step

   \param cum_l_step    Output: vector of cumulative pathlengths.
   \param ppath         a Ppath.

   \author Cory Davis
   \date   2003-06-19
*/

void cum_l_stepCalc(
		      Vector& cum_l_step,
		      const Ppath& ppath
		      )
  {
    cum_l_step.resize(ppath.np);
    Numeric cumsum = 0.0;
    cum_l_step[0] = 0.0;
    for (Index i=0; i<ppath.np-1; i++)
      {
	cumsum += ppath.l_step[i];
	cum_l_step[i+1] = cumsum;
      }
  }		


		
//! Red 1D Interpolate.
/*! 
  This is a slight modifiaction of Stefan's code to do 1_D interpolation
  to get a Matrix from an array of Matrices

  The dimension of itw must be consistent with the dimension of the
  interpolation (2^n).

  \param itw  Interpolation weights.
  \param a    The field to interpolate.(ArrayOfMatrix)
  \param tc   The grid position for the column dimension.

  \return Interpolated value.

  \author Cory Davis (modified original code by Stefan Buehler)
  \date   2003-06-19
*/

Matrix interp( ConstVectorView itw,
                ArrayOfMatrix a,    
                const GridPos&  tc )
{
#ifndef NDEBUG
  const Numeric sum_check_epsilon = 1e-6;
#endif
  
assert(is_size(itw,2));       // We need 2 interpolation
                                // weights.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one.
  assert( is_same_within_epsilon( itw.sum(),
                                  1,
                                  sum_check_epsilon ) );
  
  // To store interpolated value:
  Matrix tia(a[0].nrows(),a[0].ncols(),0.0);
  Matrix b(a[0].nrows(),a[0].ncols(),0.0);
  Index iti = 0;
  for ( Index c=0; c<2; ++c )
    {
      b=a[tc.idx+c];
      b*=itw[iti];   
      tia += b;
      ++iti;
    }

  return tia;
}		


//! Red 1D Interpolate.
/*! 
  This is a slight modifiaction of Stefan's code to do 1_D interpolation
  to get a Vector from an array of Vectors

  The dimension of itw must be consistent with the dimension of the
  interpolation (2^n).

  \param itw  Interpolation weights.
  \param a    The field to interpolate. (ArrayOfVector)
  \param tc   The grid position for the column dimension.

  \return Interpolated value.
  \author Cory Davis (modified original code by Stefan Buehler)
  \date   2003-06-19

*/

Vector interp( ConstVectorView itw,
                ArrayOfVector a,    
                const GridPos&  tc )
{
#ifndef NDEBUG
  const Numeric sum_check_epsilon = 1e-6;
#endif
  
  assert(is_size(itw,2));       // We need 2 interpolation
                                // weights.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one.
  assert( is_same_within_epsilon( itw.sum(),
                                  1,
                                  sum_check_epsilon ) );
  
  // To store interpolated value:
  Vector tia(a[0].nelem(),0.0);
  Vector b(a[0].nelem(),0.0);
  Index iti = 0;
  for ( Index c=0; c<2; ++c )
    {
      b=a[tc.idx+c];
      b*=itw[iti];   
      tia += b;
      ++iti;
    }

  return tia;
}		


//! interpTarray

/*!
   Interpolates several arrays calculated by TarrayCalc to give values at a 
   given pathlength

   \param   T             Output: transmittance matrix ( I may have made this term up ).
   \param   K_abs         Output: absorption coefficient vector
   \param   temperature   Output: 
   \param   K             Output: extinction matrix at interpolation point
   \param   rte_pos         Output: position at pathlength along ppath
   \param   rte_los         Output: LOS at pathlength along ppath
   \param   gp            Output: Gridpos of interpolated point
   \param   TArray        array of transmittance matrices
   \param   ext_matArray  array of extinction matrices
   \param   abs_vecArray  array of absorption coefficients
   \param   t_ppath       array of temperatures
   \param   cum_l_step    vector of cumulative pathlengths
   \param   pathlength    pathlength at which to calculate above values
   \param   stokes_dim    length of Stokes vector
   \param   ppath         the Ppath


   \author Cory Davis
   \date   2003-06-19
*/


//interpolates TArray and PPath to give T and rte_pos(los) at a given pathlength
void interpTArray(Matrix& T,
		  Vector& K_abs,
		  Numeric& temperature,
		  MatrixView& K,
		  Vector& rte_pos,
		  Vector& rte_los,
		  VectorView& pnd_vec,
		  ArrayOfGridPos& gp,
		  const ArrayOfMatrix& TArray,
		  const ArrayOfMatrix& ext_matArray,
		  const ArrayOfVector& abs_vecArray,
		  const Vector& t_ppath,
		  const Matrix& pnd_ppath,
		  const Vector& cum_l_step,
		  const Numeric& pathlength,
		  const Index& stokes_dim,
		  const Ppath& ppath
		  )
{
  //Internal Declarations
  Matrix incT(stokes_dim,stokes_dim);
  Matrix opt_depth_mat(stokes_dim,stokes_dim);
  Vector itw(2);
  Numeric delta_s;
  Index N_pt=pnd_vec.nelem();
  
  //interpolate transmittance matrix
  gridpos(gp, cum_l_step, pathlength);
  interpweights(itw,gp[0]);
  K = interp(itw,ext_matArray,gp[0]);
  delta_s = pathlength - cum_l_step[gp[0].idx];
  opt_depth_mat = K;
  opt_depth_mat*=-delta_s;
  matrix_exp(incT,opt_depth_mat,6);
  mult(T,TArray[gp[0].idx],incT);
  
  K_abs = interp(itw, abs_vecArray,gp[0]);
 
  temperature=interp(itw,t_ppath,gp[0]);

  for (Index i=0;i<N_pt;i++)
    {
      pnd_vec[i]=interp(itw,pnd_ppath(i,Range(joker)),gp[0]);
    }

  for (Index i=0; i<2; i++)
    {
      rte_pos[i] = interp(itw,ppath.pos(Range(joker),i),gp[0]);
      rte_los[i] = interp(itw,ppath.los(Range(joker),i),gp[0]);
    }
  rte_pos[2] = interp(itw,ppath.pos(Range(joker),2),gp[0]);
}



//!  montecarloGetIncoming 

/*!
   Gets incoming radiance at the cloud box boundary in a single propagation 
direction, determined by the cloudboz Ppath 'pptahcloud'. 
Used in ScatteringMonteCarlo.  

   \author Cory Davis
   \date   2003-11-28
*/

void montecarloGetIncoming(
			   Matrix&               i_rte,
			   Vector&               rte_pos,
			   Vector&               rte_los,
			   GridPos&              rte_gp_p,
			   GridPos&              rte_gp_lat,
			   GridPos&              rte_gp_lon,
			   Ppath&                ppath,
			   Ppath&                ppath_step,
			   Matrix&               i_space,
			   Matrix&               ground_emission,
			   Matrix&               ground_los, 
			   Tensor4&              ground_refl_coeffs,
			   Vector&               scat_za_grid,
			   Vector&               scat_aa_grid,
			   const Agenda&         ppath_step_agenda,
			   const Agenda&         rte_agenda,
			   const Agenda&         i_space_agenda,
			   const Agenda&         ground_refl_agenda,
			   const Tensor3&        t_field,
			   const Vector&         p_grid,
			   const Vector&         lat_grid,
			   const Vector&         lon_grid,
			   const Tensor3&        z_field,
			   const Matrix&         r_geoid,
			   const Matrix&         z_ground,
			   const ArrayOfIndex&   cloudbox_limits,
			   const Ppath&          ppathcloud,
			   const Index&          atmosphere_dim,
			   const Vector&         f_grid,
			   const Index&          stokes_dim
			   )

{
  //call rte_calc without input checking, sensor stuff, or verbosity
  // Assign dummies for variables associated with sensor.
  Vector   mblock_za_grid_dummy(1);
  mblock_za_grid_dummy[0] = 0;
  Vector   mblock_aa_grid_dummy(0), sensor_rot_dummy(0);
  Matrix   sensor_pol_dummy;
  Index    antenna_dim_dummy = 1; 
  Sparse   sensor_response_dummy;
  
  // Dummy for measurement vector
  Vector   y_dummy(0);
  
  const Index cloudbox_on_dummy=0;
  Matrix sensor_pos(1,3);
  Matrix sensor_los(1,2);
  Tensor7 scat_i_p_dummy;
  Tensor7 scat_i_lat_dummy;
  Tensor7 scat_i_lon_dummy;
  
  Vector dummy_rte_pos = rte_pos;
  Vector dummy_rte_los = rte_los;
  rte_posShift(dummy_rte_pos,dummy_rte_los,rte_gp_p, rte_gp_lat,
	       rte_gp_lon,ppathcloud, atmosphere_dim);
  sensor_pos(0,joker)=dummy_rte_pos;
  sensor_los(0,joker)=dummy_rte_los;
  //call rte_calc without input checking, sensor stuff, or verbosity
  rte_calc( y_dummy, ppath, ppath_step, i_rte, rte_pos, rte_los, rte_gp_p, 
	    rte_gp_lat,
	    rte_gp_lon,i_space, ground_emission, ground_los, ground_refl_coeffs,
	    ppath_step_agenda, rte_agenda, i_space_agenda, ground_refl_agenda, 
	    atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, t_field, 
	    r_geoid, z_ground, cloudbox_on_dummy, cloudbox_limits, 
	    scat_i_p_dummy,scat_i_lat_dummy, scat_i_lon_dummy, scat_za_grid,
	    scat_aa_grid, sensor_response_dummy, sensor_pos,sensor_los,
	    sensor_pol_dummy, sensor_rot_dummy,
	    f_grid,stokes_dim, antenna_dim_dummy, mblock_za_grid_dummy,
	    mblock_aa_grid_dummy, false, false, true);
  
  for (Index i = 0;i<stokes_dim;i++){assert(!isnan(i_rte(0,i)));}
}



//! pha_mat_singleCalc
/*!
 Returns the total phase matrix for given incident and scattered directions
. It requires a vector of particle number densities to be precalculated

 \param Z               Output: phase matrix
 \param za_scat         scattered 
 \param aa_scat         and
 \param za_inc          incident
 \param aa_inc          directions
 \param scat_data_raw   workspace variable
 \param stokes_dim      workspace variable
 \param f_index         workspace variable
 \param f_grid          workspace variable
 \param scat_theta      workspace variable 
 \param scat_theta_gps  workspace variable
 \param scat_theta_itws workspace variable
 \param pnd_vec         vector of particle number densities at the point 
                          in question
 \author Cory Davis
 \date   2003-11-27
*/

void pha_mat_singleCalc(
			MatrixView& Z,			
			Numeric za_scat, 
			Numeric aa_scat, 
			Numeric za_inc, 
			Numeric aa_inc,
			const ArrayOfSingleScatteringData& scat_data_mono,
			const Index&          stokes_dim,
			const VectorView& pnd_vec, 
			const Tensor4& scat_theta,//these three are a bit annoying
			const ArrayOfArrayOfArrayOfArrayOfGridPos& scat_theta_gps,//
			const Tensor5& scat_theta_itws)//
{
  Index N_pt=pnd_vec.nelem();

  assert(aa_inc>=-180 && aa_inc<=180);
  assert(aa_scat>=-180 && aa_scat<=180);



  Vector scat_za_grid(2);
  Vector scat_aa_grid(2);

  scat_za_grid[0]=za_inc;
  scat_za_grid[1]=za_scat;

  scat_aa_grid[0]=aa_inc;
  scat_aa_grid[1]=aa_scat;
 
  Index scat_za_index=1;
  Index scat_aa_index=1;
  Index za_inc_idx=0;
  Index aa_inc_idx=0;
  Matrix pha_mat_lab(stokes_dim, stokes_dim, 0.);


  Z=0.0;
  // this is a loop over the different particle types
  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      pha_matTransform(pha_mat_lab,scat_data_mono[i_pt].pha_mat_data(0,joker,joker,joker,joker,joker), 
		       scat_data_mono[i_pt].za_grid, 
		       scat_data_mono[i_pt].aa_grid,
		       scat_data_mono[i_pt].ptype,scat_za_index, scat_aa_index, 
		       za_inc_idx,aa_inc_idx, scat_za_grid, scat_aa_grid,
		       scat_theta, scat_theta_gps, scat_theta_itws);
      pha_mat_lab*=pnd_vec[i_pt];
      Z+=pha_mat_lab;
      	
    }
}

//! ppathRecordMC
/*!
  interpolates montecarlo_p_from_belowCsca by za_scat, and adds over particle
types

\param p_from_belowCsca  Output: p_from_belowCsca summed over particle types
                                 for a given scatterd zenith angle
\param za_scat           scattered zenith angle
\param montecarlo_p_from_belowCsca workspace variable
\param pnd_vec           Vector of particle number densities at the point in 
                         question
\param za_grid           The zenith angle grid used in scattering data files

 \author Cory Davis
 \date   2003-12-04

*/

Numeric p_from_belowCscaCalc(
			  const Numeric& za_scat,
			  const Tensor3& montecarlo_p_from_belowCsca,
			  const VectorView& pnd_vec,
			  const VectorView& za_grid
			  )

{
  Index N_pt=pnd_vec.nelem();
  Numeric p_from_belowCsca=0;
  GridPos gp;
  Vector itw(2);
  gridpos(gp,za_grid,za_scat);
  interpweights(itw,gp);
  for (Index pt_index=0;pt_index<N_pt;pt_index++)
    {
      p_from_belowCsca+=pnd_vec[pt_index]*
	interp(itw,montecarlo_p_from_belowCsca(pt_index,0,joker),gp);
    }
  return p_from_belowCsca;
}



//! ppathRecordMC
/*!
  Stores ppath data in XML format during Monte Carlo simulations.  
This can be useful for educational/diagnostic purposes.
\param ppath             Ppath object to be sotred
\param name              prefix for filename
\param photon_number     both these parameters are used 
\param scattering_order  in the file name
*/
void ppathRecordMC(
		   const Ppath& ppath,
		   const String name,
		   const Index& photon_number,
		   const Index& scattering_order
		   )

{
  //Record ppathcloud.  This is useful for debugging and educational purposes.  It would
  //be completely daft to leave this on in a real calculation
  ostringstream filename;
  filename << name << photon_number <<"_"<<scattering_order;
  String longfilename;
  filename_xml(longfilename,filename.str());
  xml_write_to_file(longfilename, ppath, FILE_TYPE_ASCII);
}


//! Sample_los

/*!
  Implementation one of two line of sight sampling methods determined by the 
  input Index 'sampling_method'
  
  sampling_method==1:
     Randomly samples incident direction using a probability density function 
   proportional to sin(za)(cos(za)+1). sin(za) because dsolid_angle=sin(za)dzadaa
   , and (cos(za)+1) because upward radiances tend to be greater than downward
   radiances.  NOTE: THIS IS TERRIBLE IN OPTICALLY THICK CASES AND WILL BE 
   REMOVED

  sampling_method==2:
     Randomly samples incident direction using a probability density function 
   proportional to sin(za).  Better, but again, not very good.

  THIS WHOLE FUNCTION MAY BE SOON REDUNDANT TO SAMPLE_LOS_Z WHICH USES A PDF 
  PROPORTIONAL TO Z11SINZA

   \param   rte_los          Output: incident line of sight for subsequent 
                                     ray-tracing.		      
   \param   g_los_csc_theta  Output: probability density for the chosen
                                     direction multiplied by sin(za)
   \param   rng              Rng random number generator instance
   \param   sampling_method  choice of sampling method: 1 or 2.
   \author Cory Davis
   \date   2003-06-19
*/

		
void Sample_los (
		 VectorView& rte_los,
		 Numeric& g_los_csc_theta,
		 Rng& rng,
		 const Index& sampling_method
		 )
{
 rte_los[1] = rng.draw()*360-180;
 if (sampling_method == 1)
   {
     rte_los[0] = acos(1-2*sqrt(rng.draw()));
     g_los_csc_theta = 0.25*(1-cos(rte_los[0]))/PI;
   }
 else if (sampling_method == 2)
   {
     rte_los[0] = acos(1-2*rng.draw());
     g_los_csc_theta = 0.25/PI;
   }
 else
   {
     throw runtime_error( "Invalid value for sampling_method. "
                          "sampling_method\n must be 1,2, or 3." );
   }
 rte_los[0]*=RAD2DEG;
}

//! Sample_los_Z
/*!
 Uses the rejection method to sampled za_inc and aa_inc according to a 
 probability density function proportional to 
 Z11(za_scat,aa_scat,za_inc, aa_inc)sin(za_inc)

 \param  new_rte_los      Output: new line of sight 
 \param  g_los_csc_theta  Output: prob dens. for chosen los
 \param	 Z                Output: phase matrix for given incident
                                           and scattered directions
 \param	 rng              Rng (random number generator) object
 \param  rte_los          last line of sight
 \param  scat_data_raw    workspace variable
 \param  stokes_dim       workspace variable
 \param  f_index          workspace variable
 \param  f_grid           workspace variable
 \param  scat_theta       workspace variable
 \param  scat_theta_gps   workspace variable
 \param  scat_theta_itws  workspace variable
 \param  pnd_vec          vector of particle number densities at the point 
                          in question
 \param  Numeric Csca     K11-Kabs1
 \author Cory Davis
 \date   2003-11-27
*/
void Sample_los_Z (
		   VectorView& new_rte_los,
		   Numeric& g_los_csc_theta,
		   MatrixView& Z,
		   Rng& rng,
		   const VectorView& rte_los,
		   const ArrayOfSingleScatteringData& scat_data_mono,
		   const Index&          stokes_dim,
		   const Tensor4& scat_theta, // CE: Included 
		   const ArrayOfArrayOfArrayOfArrayOfGridPos& scat_theta_gps,
		   const Tensor5& scat_theta_itws,
		   const VectorView& pnd_vec,
		   Numeric Csca
		   )
{
  bool tryagain=true;
  Matrix dummyZ(stokes_dim,stokes_dim);
  //The following is based on the assumption that the maximum value of the 
  //phase matrix for a given scattered direction is for forward scattering
  Numeric aa_scat = (rte_los[1]>=0) ?-180+rte_los[1]:180+rte_los[1];
  pha_mat_singleCalc(dummyZ,180-rte_los[0],aa_scat,180-rte_los[0],
		     aa_scat,scat_data_mono,stokes_dim,pnd_vec,
		     scat_theta,scat_theta_gps,
			 scat_theta_itws);
  Numeric Z11max=dummyZ(0,0);  
  ///////////////////////////////////////////////////////////////////////  
  while(tryagain)
    {
      new_rte_los[1] = rng.draw()*360-180;
      new_rte_los[0] = acos(1-2*rng.draw())*RAD2DEG;
      //Calculate Phase matrix////////////////////////////////
      Numeric aa_inc= (new_rte_los[1]>=0) ?
	-180+new_rte_los[1]:180+new_rte_los[1];
      
      pha_mat_singleCalc(Z,180-rte_los[0],aa_scat,180-new_rte_los[0],
			 aa_inc,scat_data_mono,stokes_dim,pnd_vec,scat_theta,
			 scat_theta_gps, scat_theta_itws);
      
      if (rng.draw()<=Z(0,0)/Z11max)//then new los is accepted
	{
	  tryagain=false;
	}
    }
  g_los_csc_theta =Z(0,0)/Csca;
}



//! Sample_ppathlength

/*!
   Randomly samples path length from an exponential distribution based on the 
   extinction matrix along the line of sight

   \param   pathlength        Output: the pathlength.
   \param   g                 Output: the probability density of the returned 
                              pathlength, or if the pathlength > dist_to_boundary
                              it gives the probability of that occurance 
   \param   rng               Rng random number generator instance
   \param   ext_matArray      An array of extinction matrices along the line of 
                              sight.
   \param   TArray            An array of the evolution operator along the line 
                              of sight
   \param   cum_l_step        An array of cumulative distance along the line of 
                              sight
   \param   method            An index selecting the method used to sample 
                              ppathlength.  method==1 uses an exponential PDF 
                              using the extinction coefficient averaged along the
                              line of sight.
   \author Cory Davis
   \date   2003-10-03
*/
void Sample_ppathlength (
			 Numeric& pathlength, 
			 Numeric& g,
			 // Numeric& K11,
			 Rng& rng,
			 const ArrayOfMatrix& ext_matArray,
			 const ArrayOfMatrix& TArray,
			 const ConstVectorView& cum_l_step,
			 Index method
			 )
{
  Index npoints=cum_l_step.nelem();
  Numeric dist_to_boundary=cum_l_step[npoints-1];	
  Numeric K11;
  assert(ext_matArray.nelem()==npoints);
  assert(TArray.nelem()==npoints);

  if(method==1)
    {
      //Since we already have an array 
      //of extinction matrix elements we could choose a number of ways
      //to sample the pathlength.  To start with we'll try using the mean 
      //extinction matrix along the line of sight.
      Matrix K = ext_matArray[0];
      for (Index i=1;i<npoints;i++){
	K+=ext_matArray[i];
      }
      K/=ext_matArray.nelem();
      K11 = K(0,0);
      pathlength = -log(rng.draw())/K11;
      
      if(pathlength>dist_to_boundary)
	{
	  g=exp(-dist_to_boundary*K11);//probability of leaving the cloudbox
	}
      else
	{
	  g=K11*exp(-pathlength*K11);//probability density at the chosen pathlength
	}
    }
  else
    {
      //assert(false);
      Numeric r = rng.draw();
      //inefficient first effort
      Vector T11vector(npoints);
      for (Index i=0;i<npoints;i++)
	{
	  T11vector[i]=TArray[i](0,0);
	}
      GridPos gp;
      Vector itw(2);
      if(r < T11vector[npoints-1])
	{
	  //Do something appropriate. photon has left the building.
	  pathlength=1e9;
	  g=T11vector[npoints-1];
	}
      else
	{
	  gridpos(gp,T11vector,r);
	  interpweights(itw,gp);
	  Numeric T11;
	  T11=interp(itw,T11vector,gp);
	  K11=0.5*(ext_matArray[gp.idx](0,0)+ext_matArray[gp.idx+1](0,0));
	  pathlength=cum_l_step[gp.idx]+log(T11vector[gp.idx]/T11)/K11;
	  g=K11*T11;
	}
    }	
}		

//! Sample_ppathlengthLOS

/*!
   Similar to Sample_ppathlength, but ensures that the sampled point lies within the cloudbox

   \param   ppathlength       Output: the pathlength.
   \param   g                 Output: the probability density of the returned pathlength.
   \param   rng               Rng random number generator instance
   \param   ext_matArray      An array of extinction matrices along the line of sight.
   \param   dist_to_boundary  the distance from the current position to the far boundary 
                              along the line of sight
   \author Cory Davis
   \date   2003-03-10
*/

void Sample_ppathlengthLOS (
			 Numeric& pathlength, 
			 Numeric& g,
			 Rng& rng,
			 const ArrayOfMatrix& ext_matArray,
			 const ArrayOfMatrix& TArray,
			 const ConstVectorView& cum_l_step,
			 Index method
			 )
{
  Index npoints=cum_l_step.nelem();
  Numeric dist_to_boundary=cum_l_step[npoints-1];	
  Numeric K11;
  assert(ext_matArray.nelem()==npoints);
  assert(TArray.nelem()==npoints);
	
  if(method==1)
    {
      //Since we already have an array 
      //of extinction matrix elements we could choose a number of ways
      //to sample the pathlength.  To start with we'll try using the mean 
      //extinction matrix along the line of sight.
      Matrix K = ext_matArray[0];
      for (Index i=1;i<ext_matArray.nelem();i++){
	K+=ext_matArray[i];
      }
      K/=ext_matArray.nelem();
      K11 = K(0,0);
      //	cout << "K11 = "<< K11 <<"\n";
      Numeric r=rng.draw();
      
      pathlength = -log(exp(-K11*dist_to_boundary)*(1-r)+r)/K11;
      
      g=K11*exp(-pathlength*K11)/(1-exp(-K11*dist_to_boundary));
    }
  else
    {
      Numeric r = rng.draw();
      //inefficient first effort
      Vector T11vector(npoints);
      for (Index i=0;i<npoints;i++)
	{
	  T11vector[i]=TArray[i](0,0);
	}
      GridPos gp;
      Vector itw(2);
      gridpos(gp,T11vector,1-r*(1-T11vector[npoints-1]));
      interpweights(itw,gp);
      Numeric T11;
      T11=interp(itw,T11vector,gp);
      K11=0.5*(ext_matArray[gp.idx](0,0)+ext_matArray[gp.idx+1](0,0));
      pathlength=cum_l_step[gp.idx]+log(T11vector[gp.idx]/T11)/K11;
      g=K11*T11/(1-T11vector[npoints-1]);
    }
  assert(pathlength<=dist_to_boundary);
}		

//!  TArrayCalc
/*! 
 
   This function calculates arrays of useful variables along a propagation
   path within the cloudbox. 


   \author Cory Davis
   \date   2003-07-19
*/


		
void TArrayCalc(
		//output
		ArrayOfMatrix& TArray,
		ArrayOfMatrix& ext_matArray,
		ArrayOfVector& abs_vecArray,
		Vector& t_ppath,
		Vector& scat_za_grid,
		Vector& scat_aa_grid,
		Tensor3& ext_mat,
		Matrix& abs_vec,
		Numeric&   rte_pressure,
		Numeric&   rte_temperature,
		Vector&    rte_vmr_list,
		Index&    scat_za_index,
		Index&    scat_aa_index,
		Tensor3& ext_mat_spt,
		Matrix& abs_vec_spt,
		Matrix&  pnd_ppath,
		//input
		const Ppath& ppath,
		const Agenda& opt_prop_gas_agenda,
		const Agenda& spt_calc_agenda,
		const Agenda& scalar_gas_absorption_agenda,
		const Index& stokes_dim,
		const Vector&    p_grid,
		const Vector&    lat_grid,
		const Vector&    lon_grid,
		const Tensor3&   t_field,
		const Tensor4&   vmr_field,
		const Index&     atmosphere_dim,
		const Tensor4&   pnd_field
 		)
{ 
  const Index np=ppath.np;  
  const Index   ns = vmr_field.nbooks();
  const Index N_pt = pnd_field.nbooks();

  TArray.resize(np);
  ext_matArray.resize(np); 
  abs_vecArray.resize(np);
  pnd_ppath.resize(N_pt,np);
  t_ppath.resize(np);
  Matrix opt_depth_mat(stokes_dim,stokes_dim),incT(stokes_dim,stokes_dim);
  Matrix zeroMatrix(stokes_dim,stokes_dim,0.0);
  Matrix identity(stokes_dim,stokes_dim,0.0);
  //Identity matrix
  for (Index i=0; i<stokes_dim; i++){identity(i,i)=1.0;}
 
  Matrix ext_mat_part(stokes_dim, stokes_dim, 0.0);
  Vector abs_vec_part(stokes_dim, 0.0);

//use propagation angles from ppath to form scat_za_grid, scat_aa_grid
  scat_za_grid.resize(ppath.np);
  scat_aa_grid.resize(ppath.np);//I don't think aa can change along a 
				//propagation path		

  scat_za_grid=ppath.los(Range(joker),0);
  scat_za_grid*=-1.0;
  scat_za_grid+=180;
  scat_aa_grid=ppath.los(Range(joker),1);
  scat_aa_grid+=180;

  // Determine the pressure at each propagation path point
  Vector   p_ppath(np);
  Matrix   itw_p(np,2);
  //
  interpweights( itw_p, ppath.gp_p );      
  itw2p( p_ppath, p_grid, ppath.gp_p, itw_p );
  
  // Determine the atmospheric temperature and species VMR at 
  // each propagation path point
  Matrix   vmr_ppath(ns,np), itw_field;
  //
  interp_atmfield_gp2itw( itw_field, atmosphere_dim, p_grid, lat_grid, 
			  lon_grid, ppath.gp_p, ppath.gp_lat, ppath.gp_lon );
  //
  interp_atmfield_by_itw( t_ppath,  atmosphere_dim, p_grid, lat_grid, 
			  lon_grid, t_field, "t_field", 
			  ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
  // 
  for( Index is=0; is<ns; is++ )
    {
      interp_atmfield_by_itw( vmr_ppath(is, joker), atmosphere_dim,
			      p_grid, lat_grid, lon_grid, 
			      vmr_field(is, joker, joker,  joker), 
			      "vmr_field", ppath.gp_p, ppath.gp_lat, 
			      ppath.gp_lon, itw_field );
    }
  
  //Determine the particle number density for every particle type at 
  // each propagation path point
  

  for( Index ip=0; ip<N_pt; ip++ )
    {
      interp_atmfield_by_itw( pnd_ppath(ip, joker), atmosphere_dim,
			      p_grid, lat_grid, lon_grid, 
			      pnd_field(ip, joker, joker,  joker), 
			      "pnd_field", ppath.gp_p, ppath.gp_lat, 
			      ppath.gp_lon, itw_field );
    }



//Create array of extinction matrices corresponding to each point in ppath
  for (scat_za_index=0; scat_za_index<ppath.np; scat_za_index++)
    { 
      scat_aa_index=scat_za_index;
      rte_pressure    = p_ppath[scat_za_index];
      rte_temperature = t_ppath[scat_za_index];
      rte_vmr_list    = vmr_ppath(joker,scat_za_index);
      scalar_gas_absorption_agenda.execute( true );
      opt_prop_gas_agenda.execute( true );
      ext_mat_part=0.0;
      abs_vec_part=0.0;
      //Make sure scat_aa is between -180 and 180
      if (scat_aa_grid[scat_aa_index]>180){scat_aa_grid[scat_aa_index]-=360;}
      //
      //opt_prop_part_agenda.execute( true );
      //use pnd_ppath and ext_mat_spt to get extmat (and similar for abs_vec
      spt_calc_agenda.execute( true );
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt; l++)
        { 
          // now the last two loops over the stokes dimension.
          for (Index m = 0; m < stokes_dim; m++)
            {
	      abs_vec_part[m] += (abs_vec_spt(l, m) * pnd_ppath(l, scat_za_index));
              for (Index n = 0; n < stokes_dim; n++)
		ext_mat_part(m, n) +=  (ext_mat_spt(l, m, n) * 
                                        pnd_ppath(l, scat_za_index));
              
            } 
        }

      //Add particle extinction matrix to *ext_mat*.
      ext_mat(0, Range(joker), Range(joker)) += ext_mat_part;
      abs_vec(0,Range(joker)) += abs_vec_part;

      ext_matArray[scat_za_index]=ext_mat(0,joker,joker);
      abs_vecArray[scat_za_index]=abs_vec(0,joker);
    }
  //create an array of T matrices corresponding to each position in the
  //the first ppath through the cloudbox each grid cell spanned by ppath points
  opt_depth_mat=ext_matArray[0];
  opt_depth_mat+=ext_matArray[1];
  opt_depth_mat*=-ppath.l_step[0]/2;
  TArray[0]=identity;
  for (Index i=1; i<ppath.np;i++)
    {
      //Get the appropriate extinction matrix for the cell in question
      opt_depth_mat=ext_matArray[i];
      opt_depth_mat+=ext_matArray[i-1];
      opt_depth_mat*=-ppath.l_step[i-1]/2;
      matrix_exp(incT,opt_depth_mat,6);
      TArray[i]=zeroMatrix;
      mult(TArray[i],TArray[i-1],incT);
    }

}

