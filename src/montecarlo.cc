
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



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/



//! Cloudbox_ppath_rteCalc

/*!
   This function was written to eliminate some repeated code in ScatteringMonteCarlo.
   Basically this function performs internal and external ppath calculations, calculates
   incoming radiances, and provides scattering properties and other derived data that is 
  required in the monte carlo simulation for a given propagation direction.


   Parameters not already described elsewhere:
   
   \param    record_ppathcloud   A flag that determines whether ppaths within the cloud box are recorded                               for later plotting in MATLAB.
   \param    record_ppath        Similar, but for ppaths determining incoming radiances at the cloudbo                                 x boundaries.  For normal use both of these parameterss should be 
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
			     const Agenda&         opt_prop_part_agenda,
			     const Agenda&         scalar_gas_absorption_agenda,
			     const Index&          stokes_dim,
			     const Tensor3&        t_field,
			     const Tensor4&        vmr_field,
			     const Agenda&         rte_agenda,
			     const Agenda&         i_space_agenda,
			     const Agenda&         ground_refl_agenda,
			     const Vector&         f_grid,
			     const Index&          photon_number,
			     const Index&          scattering_order)

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
  TArrayCalc(TArray, ext_matArray, abs_vecArray, t_ppath, scat_za_grid, scat_aa_grid, ext_mat, abs_vec,
	     rte_pressure, rte_temperature, rte_vmr_list, scat_za_index, scat_aa_index, 
	     ppathcloud, opt_prop_gas_agenda, 
	     opt_prop_part_agenda, scalar_gas_absorption_agenda, stokes_dim, 
	     p_grid, lat_grid, lon_grid, t_field, vmr_field, atmosphere_dim);
  //Calculate contribution from the boundary of the cloud box
  //changed to dummy_rte_pos to see if rte_pos was causing assertion failure at ppath.cc:1880
  //it appears that this was not the case
  Vector dummy_rte_pos = rte_pos;
  Vector dummy_rte_los = rte_los;  
  shift_rte_pos(dummy_rte_pos,dummy_rte_los,ppathcloud);
  sensor_pos(0,joker)=dummy_rte_pos;
  sensor_los(0,joker)=dummy_rte_los;
  //call rte_calc without input checking, sensor stuff, or verbosity
  rte_calc( y_dummy, ppath, ppath_step, i_rte, rte_pos, rte_los, rte_gp_p, rte_gp_lat,
	    rte_gp_lon,i_space, ground_emission, ground_los, ground_refl_coeffs,
	    ppath_step_agenda, rte_agenda, i_space_agenda, ground_refl_agenda, atmosphere_dim,
	    p_grid, lat_grid, lon_grid, z_field, t_field, r_geoid, z_ground, cloudbox_on_dummy,
	    cloudbox_limits, scat_i_p_dummy,scat_i_lat_dummy, scat_i_lon_dummy, scat_za_grid,
	    scat_aa_grid, sensor_response_dummy, sensor_pos,sensor_los,
            sensor_pol_dummy, sensor_rot_dummy,
            f_grid,stokes_dim,
	    antenna_dim_dummy, mblock_za_grid_dummy,mblock_aa_grid_dummy,false, false, true);
  
  
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
const Numeric sum_check_epsilon = 1e-6;
  
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
const Numeric sum_check_epsilon = 1e-6;
  
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
		  Vector& rte_pos,
		  Vector& rte_los,
		  ArrayOfGridPos& gp,
		  const ArrayOfMatrix& TArray,
		  const ArrayOfMatrix& ext_matArray,
		  const ArrayOfVector& abs_vecArray,
		  const Vector& t_ppath,
		  const Vector& cum_l_step,
		  const Numeric& pathlength,
		  const Index& stokes_dim,
		  const Ppath& ppath
		  )
{
  //Internal Declarations
  Matrix incT(stokes_dim,stokes_dim);
  Matrix K(stokes_dim,stokes_dim);
  Matrix opt_depth_mat(stokes_dim,stokes_dim);
  Vector itw(2);
  Numeric delta_s;
  
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

  for (Index i=0; i<2; i++)
    {
      rte_pos[i] = interp(itw,ppath.pos(Range(joker),i),gp[0]);
      rte_los[i] = interp(itw,ppath.los(Range(joker),i),gp[0]);
    }
  rte_pos[2] = interp(itw,ppath.pos(Range(joker),2),gp[0]);
}


//! Sample_los

/*!
   Randomly samples incident direction using a probability density function 
   proportional to sin(za)(cos(za)+1). sin(za) because dsolid_angle=sin(za)dzadaa
   , and (cos(za)+1) because upward radiances tend to be greater than downward
   radiances.

   \param   rte_los       Output: incident line of sight for subsequent ray-tracing.
   \param   rng         Rng random number generator instance
   \author Cory Davis
   \date   2003-06-19
*/

		
void Sample_los (
		 Vector& rte_los,
		 Rng& rng
		 )
{
  rte_los[0] = 180-RAD2DEG*acos(2*sqrt(rng.draw())-1);
  rte_los[1] = rng.draw()*360-180;
}



//! Sample_ppathlength

/*!
   Randomly samples path length from an exponential distribution based on the 
   extinction matrix along the line of sight
   Note: currently the same sequence of pathlengths will be produced for every 
   arts run.  This is useful during development for replicating errors.  For 
   serious calculations this will be changed.

   \param   ppathlength       Output: the pathlength.
   \param   g                 Output: the probability density of the returned pathlength.
   \param   rng               Rng random number generator instance
   \param   ext_matArray      An array of extinction matrices along the line of sight.
   
   \author Cory Davis
   \date   2003-06-19
*/
void Sample_ppathlength (
			 Numeric& pathlength, 
			 Numeric& g,
			 Rng& rng,
			 const ArrayOfMatrix& ext_matArray
			 )
{		
	//Since we already have an array 
	//of extinction matrix elements we could choos a number of ways
	//to sample the pathlength.  To start with we'll try using the 
	//extinction matrix at the point in question. 
	Matrix K = ext_matArray[0];
        Numeric K11 = K(0,0);
	//	cout << "K11 = "<< K11 <<"\n";
	pathlength = -log(rng.draw())/K11;

        g=K11*exp(-pathlength*K11);	
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
		//input
		const Ppath& ppath,
		const Agenda& opt_prop_gas_agenda,
		const Agenda& opt_prop_part_agenda,
		const Agenda& scalar_gas_absorption_agenda,
		const Index& stokes_dim,
		const Vector&    p_grid,
		const Vector&    lat_grid,
		const Vector&    lon_grid,
		const Tensor3&   t_field,
		const Tensor4&   vmr_field,
		const Index&     atmosphere_dim
 		)
{ 
  const Index np=ppath.np;  
  const Index   ns = vmr_field.nbooks();
  TArray.resize(np);
  ext_matArray.resize(np); 
  abs_vecArray.resize(np);
  t_ppath.resize(np);
Matrix opt_depth_mat(stokes_dim,stokes_dim),incT(stokes_dim,stokes_dim);
  Matrix zeroMatrix(stokes_dim,stokes_dim,0.0);
  Matrix identity(stokes_dim,stokes_dim,0.0);
  //Identity matrix
  for (Index i=0; i<stokes_dim; i++){identity(i,i)=1.0;}
 //use propagation angles from ppath to form scat_za_grid, scat_aa_grid
  scat_za_grid.resize(ppath.np);
  scat_aa_grid.resize(ppath.np);//I don't think aa can change along a propagation path

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
			      "vmr_field", ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
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

      //Make sure scat_aa is between -180 and 180
      if (scat_aa_grid[scat_aa_index]>180){scat_aa_grid[scat_aa_index]-=360;}
      opt_prop_part_agenda.execute( true );
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

