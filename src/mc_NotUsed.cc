  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_z_field_is_1D" ),
       DESCRIPTION
       (
        "Flag for MCGeneral and associated methods.\n"
        "\n"
        "If fields outside the cloud box are 1D, this flag can be set to 1\n"
        "and the calculations will be more rapid.\n"
        "\n"
        "Usage: Set by the user.\n"
        ),
       GROUP( "Index" )));



  md_data_raw.push_back
    ( MdRecord
      ( NAME( "mc_IWP_cloud_opt_pathCalc" ),
        DESCRIPTION
        (
         "Calculates the FOV averaged ice water path and cloud optical path\n"
         "for a given viewing direction\n"
         ),
        AUTHORS( "Cory Davis" ),
        OUT( "mc_IWP", "mc_cloud_opt_path", "mc_IWP_error", 
             "mc_cloud_opt_path_error", "mc_iteration_count" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "mc_antenna", "sensor_pos", "sensor_los", "ppath_step_agenda", 
            "p_grid", "lat_grid", "lon_grid", "refellipsoid", "z_surface", 
            "z_field", "t_field", "vmr_field", "edensity_field", 
            "f_grid", "f_index", "cloudbox_limits", "pnd_field", 
            "scat_data_array_mono", "particle_masses", "mc_seed" ),
        GIN( "max_iter" ),
        GIN_TYPE(    "Index" ),
        GIN_DEFAULT( NODEF ),
        GIN_DESC( "Maximum number of iterations." )
        ));
  
/* Workspace method: Doxygen documentation will be auto-generated */
void mc_IWP_cloud_opt_pathCalc(Workspace& ws,
                               //output
                               Numeric& mc_IWP,
                               Numeric& mc_cloud_opt_path,
                               Numeric& mc_IWP_error,
                               Numeric& mc_cloud_opt_path_error,
                               Index&   mc_iteration_count,
                               //input
                               const MCAntenna&      mc_antenna,
                               const Matrix&         sensor_pos,
                               const Matrix&         sensor_los,
                               const Agenda&         ppath_step_agenda,
                               const Vector&         p_grid,
                               const Vector&         lat_grid, 
                               const Vector&         lon_grid, 
                               const Vector&         refellipsoid, 
                               const Matrix&         z_surface,
                               const Tensor3&        z_field, 
                               const Tensor3&        t_field, 
                               const Tensor4&        vmr_field, 
                               const Tensor3&        edensity_field, 
                               const Vector&         f_grid,
                               const Index&          f_index,
                               const ArrayOfIndex&   cloudbox_limits, 
                               const Tensor4&        pnd_field,
                               const ArrayOfSingleScatteringData& scat_data_array_mono,
                               const Matrix&         particle_masses,
                               const Index&          mc_seed,
                               //Keyword params
                               const Index&          max_iter,
                               const Verbosity&      verbosity)
{
  const Numeric f_mono = f_grid[f_index];
  
  if( particle_masses.ncols() != 1 )
    throw runtime_error(
                 "*particle_masses* is only allowed to have a single column" );

  //if antenna is pencil beam just need a single path integral, otherwise
  //for now do monte carlo integration
  if (mc_antenna.get_type()==ANTENNA_TYPE_PENCIL_BEAM)
    {
      iwp_cloud_opt_pathCalc(ws, mc_IWP,mc_cloud_opt_path,sensor_pos(0,joker),
                             sensor_los(0,joker),
                             ppath_step_agenda,p_grid,lat_grid,lon_grid,
                             refellipsoid,z_surface,z_field,t_field,vmr_field,
                             edensity_field, f_mono, cloudbox_limits,
                             pnd_field,scat_data_array_mono,particle_masses,
                             verbosity);
    }
  else
    {
      Numeric iwp,cloud_opt_path;
      Numeric iwp_squared=0;
      Numeric tau_squared=0;
      mc_iteration_count=0;
      mc_IWP=0;
      mc_cloud_opt_path=0;
      Vector local_rte_los(2);
      Rng rng;
      rng.seed(mc_seed, verbosity);
      //Begin Main Loop
      while (mc_iteration_count<max_iter)
        {
          mc_antenna.draw_los(local_rte_los,rng,sensor_los(0,joker));
          mc_iteration_count+=1;
          iwp_cloud_opt_pathCalc(ws, iwp,cloud_opt_path,sensor_pos(0,joker),
                                 local_rte_los,
                                 ppath_step_agenda,p_grid,lat_grid,lon_grid,
                                 refellipsoid,z_surface,z_field,t_field,
                                 vmr_field,edensity_field, f_mono,
                                 cloudbox_limits,
                                 pnd_field,scat_data_array_mono,particle_masses,
                                 verbosity);
          mc_IWP+=iwp;
          iwp_squared+=iwp*iwp;
          mc_cloud_opt_path+=cloud_opt_path;
          tau_squared+=cloud_opt_path*cloud_opt_path;
        }
      mc_IWP/=(Numeric)mc_iteration_count;
      mc_cloud_opt_path/=(Numeric)mc_iteration_count;
      mc_IWP_error=sqrt((iwp_squared/(Numeric)mc_iteration_count-mc_IWP*mc_IWP)
                        /(Numeric)mc_iteration_count);
      mc_cloud_opt_path_error=sqrt((tau_squared/(Numeric)mc_iteration_count-
                                    mc_cloud_opt_path*mc_cloud_opt_path)
                                   /(Numeric)mc_iteration_count);
      
    }

}





// The version that was used up to 2012-11-20:
//! mcPathTraceGeneral
/*!
Performs the tasks of pathlength sampling,
ray tracing (but now only as far as determined by pathlength
sampling) and calculation of the evolution operator and several 
atmospheric variables at the new point.

\author Cory Davis
\date 2005-2-21

*/
void mcPathTraceGeneral(Workspace&            ws,
                        MatrixView            evol_op,
                        Vector&               abs_vec_mono,
                        Numeric&              temperature,
                        MatrixView            ext_mat_mono,
                        Rng&                  rng,
                        Vector&               rte_pos,
                        Vector&               rte_los,
                        Vector&               pnd_vec,
                        Numeric&              g,
                        Ppath&                ppath_step,
                        Index&                termination_flag,
                        bool&                 inside_cloud,
                        const Agenda&         ppath_step_agenda,
                        const Agenda&         propmat_clearsky_agenda,
                        const Index           stokes_dim,
                        const Numeric&        f_mono,
                        const Vector&         p_grid,
                        const Vector&         lat_grid,
                        const Vector&         lon_grid,
                        const Tensor3&        z_field,
                        const Vector&         refellipsoid,
                        const Matrix&         z_surface,
                        const Tensor3&        t_field,
                        const Tensor4&        vmr_field,
                        const Tensor3&        edensity_field,
                        const ArrayOfIndex&   cloudbox_limits,
                        const Tensor4&        pnd_field,
                        const ArrayOfSingleScatteringData& scat_data_array_mono,
                        const Verbosity& verbosity)
                        // 2011-06-17 GH commented out, unused?
                        // const Index           z_field_is_1D)

{ 
  ArrayOfMatrix evol_opArray(2);
  ArrayOfMatrix ext_matArray(2);
  ArrayOfVector abs_vecArray(2),pnd_vecArray(2);
  Vector tArray(2);
  Vector cum_l_step;
  Matrix T(stokes_dim,stokes_dim);
  Numeric k;
  Numeric ds;
  Index np=0;
  Index   istep = 0;            // Counter for number of steps
  Matrix opt_depth_mat(stokes_dim,stokes_dim),incT(stokes_dim,stokes_dim,0.0);
  Matrix old_evol_op(stokes_dim,stokes_dim);
  //g=0;k=0;ds=0;
  //at the start of the path the evolution operator is the identity matrix
  id_mat(evol_op);
  evol_opArray[1]=evol_op;
  //initialise Ppath with ppath_start_stepping
  //Index rubbish=z_field_is_1D;rubbish+=1; // 2011-06-17 GH commented out, unused?
  ppath_start_stepping( ppath_step, 3, p_grid, lat_grid, 
                        lon_grid, z_field, refellipsoid, z_surface,
                        0, cloudbox_limits, false, 
                        rte_pos, rte_los, verbosity );
  //Use cloudbox_ppath_start_stepping to avoid unnecessary z_at_latlon calls.
  //cloudbox_ppath_start_stepping( ppath_step, 3, p_grid, lat_grid, 
  //                               lon_grid, z_field, refellipsoid, z_surface, rte_pos,
  //                               rte_los ,z_field_is_1D);
  Range p_range(cloudbox_limits[0], 
                cloudbox_limits[1]-cloudbox_limits[0]+1);
  Range lat_range(cloudbox_limits[2], 
                cloudbox_limits[3]-cloudbox_limits[2]+1);
  Range lon_range(cloudbox_limits[4], 
                cloudbox_limits[5]-cloudbox_limits[4]+1);
  termination_flag=0;

  inside_cloud=is_inside_cloudbox( ppath_step, cloudbox_limits,true );
  
  if (inside_cloud)
    {
      cloudy_rt_vars_at_gp(ws,ext_mat_mono,abs_vec_mono,pnd_vec,temperature,
                           propmat_clearsky_agenda,
                           stokes_dim, f_mono, ppath_step.gp_p[0], ppath_step.gp_lat[0],
                           ppath_step.gp_lon[0],p_grid[p_range], 
                           t_field(p_range,lat_range,lon_range), 
                           vmr_field(joker,p_range,lat_range,lon_range),pnd_field,
                           scat_data_array_mono, cloudbox_limits,ppath_step.los(0,joker),
                           verbosity);
    }
  else
    {
      clear_rt_vars_at_gp( ws, ext_mat_mono, abs_vec_mono, temperature, 
           propmat_clearsky_agenda, f_mono,
            ppath_step.gp_p[0], ppath_step.gp_lat[0], ppath_step.gp_lon[0],
            p_grid, t_field, vmr_field );
      pnd_vec=0.0;
    }
  ext_matArray[1]=ext_mat_mono;
  abs_vecArray[1]=abs_vec_mono;
  tArray[1]=temperature;
  pnd_vecArray[1]=pnd_vec;
  //draw random number to determine end point
  Numeric r = rng.draw();
  while ((evol_op(0,0)>r)&(!termination_flag))
    {
      istep++;
      //we consider more than 5000
      // path points to be an indication on that the calcululations have
      // got stuck in an infinite loop.
      if( istep > 5000 )
        {
          throw runtime_error(
                            "5000 path points have been reached. Is this an infinite loop?" );
        }
      evol_opArray[0]=evol_opArray[1];
      ext_matArray[0]=ext_matArray[1];
      abs_vecArray[0]=abs_vecArray[1];
      tArray[0]=tArray[1];
      pnd_vecArray[0]=pnd_vecArray[1];
      //perform single path step using ppath_step_geom_3d
      ppath_step_geom_3d(ppath_step, lat_grid, lon_grid, z_field,
                         refellipsoid, z_surface, -1);

      // For debugging:
      // Print( ppath_step, 0, verbosity );

      np=ppath_step.np;//I think this should always be 2.
      cum_l_stepCalc(cum_l_step, ppath_step);
      //path_step should now have two elements.
      //calculate evol_op
      inside_cloud=is_inside_cloudbox( ppath_step, cloudbox_limits, true );
      if (inside_cloud)
        {
          cloudy_rt_vars_at_gp(ws,ext_mat_mono,abs_vec_mono,pnd_vec,temperature,
                               propmat_clearsky_agenda, stokes_dim, f_mono,
                               ppath_step.gp_p[np-1],ppath_step.gp_lat[np-1],
                               ppath_step.gp_lon[np-1],p_grid[p_range], 
                               t_field(p_range,lat_range,lon_range), 
                               vmr_field(joker,p_range,lat_range,lon_range),pnd_field,
                               scat_data_array_mono, cloudbox_limits,ppath_step.los(np-1,joker),
                               verbosity);
        }
      else
        {
          clear_rt_vars_at_gp(ws, ext_mat_mono,abs_vec_mono,temperature, 
               propmat_clearsky_agenda, f_mono,
               ppath_step.gp_p[np-1], ppath_step.gp_lat[np-1],
               ppath_step.gp_lon[np-1], p_grid, t_field, vmr_field);
          pnd_vec=0.0;
        }
      ext_matArray[1]=ext_mat_mono;
      abs_vecArray[1]=abs_vec_mono;
      tArray[1]=temperature;
      pnd_vecArray[1]=pnd_vec;
      opt_depth_mat=ext_matArray[1];
      opt_depth_mat+=ext_matArray[0];
      opt_depth_mat*=-cum_l_step[np-1]/2;
      incT=0;
      if ( stokes_dim == 1 )
        {
          incT(0,0)=exp(opt_depth_mat(0,0));
        }
      else if ( is_diagonal( opt_depth_mat))
        {
          for ( Index j=0;j<stokes_dim;j++)
            {
              incT(j,j)=exp(opt_depth_mat(j,j));
            }
        }
      else
        {
          //matrix_exp(incT,opt_depth_mat,10);
          matrix_exp_p30(incT,opt_depth_mat);
        }
      mult(evol_op,evol_opArray[0],incT);
      evol_opArray[1]=evol_op;
     
      //check whether hit ground or space
      Index bg = ppath_what_background(ppath_step);
      if ( bg==2 )
        {
          //we have hit the surface
          termination_flag=2;
          
        }
      //else if ( is_gridpos_at_index_i( ppath_step.gp_p[np-1], p_grid.nelem() - 1 ) )
      else if (fractional_gp(ppath_step.gp_p[np-1])>=(Numeric)(p_grid.nelem() - 1)-1e-3)
        {
          //we have left the top of the atmosphere
          termination_flag=1;
          
        }


    }
  if (termination_flag)
    {//we must have reached the cloudbox boundary
      //
      rte_pos = ppath_step.pos(np-1,Range(0,3));
      rte_los = ppath_step.los(np-1,joker);
      g=evol_op(0,0);
    }
  else
    {
      //find position...and evol_op..and everything else required at the new
      //scattering/emission point
      // GH 2011-09-14: 
      //   log(incT(0,0)) = log(exp(opt_depth_mat(0, 0))) = opt_depth_mat(0, 0)
      //   Avoid loss of precision, use opt_depth_mat directly
      //k=-log(incT(0,0))/cum_l_step[np-1];//K=K11 only for diagonal ext_mat
      k=-opt_depth_mat(0,0)/cum_l_step[np-1];//K=K11 only for diagonal ext_mat
      ds=log(evol_opArray[0](0,0)/r)/k;
      g=k*r;
      Vector x(2,0.0);
      //interpolate atmospheric variables required later
      //length 2
      ArrayOfGridPos gp(1);
      x[1]=cum_l_step[np-1];
      Vector itw(2);
  
      gridpos(gp, x, ds);
      assert(gp[0].idx==0);
      interpweights(itw,gp[0]);
      interp(ext_mat_mono, itw, ext_matArray, gp[0]);
      opt_depth_mat = ext_mat_mono;
      opt_depth_mat+=ext_matArray[gp[0].idx];
      opt_depth_mat*=-ds/2;
      if ( stokes_dim == 1 )
        {
          incT(0,0)=exp(opt_depth_mat(0,0));
        }
      else if ( is_diagonal( opt_depth_mat))
        {
          for ( Index i=0;i<stokes_dim;i++)
            {
              incT(i,i)=exp(opt_depth_mat(i,i));
            }
        }
      else
        {
          //matrix_exp(incT,opt_depth_mat,10);
          matrix_exp_p30(incT,opt_depth_mat);
        }
      mult(evol_op,evol_opArray[gp[0].idx],incT);
      interp(abs_vec_mono, itw, abs_vecArray,gp[0]);
      temperature=interp(itw,tArray,gp[0]);
      interp(pnd_vec, itw,pnd_vecArray,gp[0]);
      if (np > 2)
        {
          gridpos(gp,cum_l_step,ds);
          interpweights(itw,gp[0]);
        }
      for (Index i=0; i<2; i++)
        {
          rte_pos[i] = interp(itw,ppath_step.pos(Range(joker),i),gp[0]);
          rte_los[i] = interp(itw,ppath_step.los(Range(joker),i),gp[0]);
        }
      rte_pos[2] = interp(itw,ppath_step.pos(Range(joker),2),gp[0]);
    }
  assert(isfinite(g));
}



  md_data_raw.push_back     
    ( MdRecord
      ( NAME( "MCIPA" ),
        DESCRIPTION
        ( "A specialised 3D reversed Monte Carlo radiative algorithm, that\n"
          "mimics independent pixel appoximation simulations.\n"
          ),
        AUTHORS( "Cory Davis" ),
        OUT( "y", "mc_iteration_count", "mc_error", "mc_points" ),
        GOUT(),
        GOUT_TYPE(),
        GOUT_DESC(),
        IN( "mc_antenna", "f_grid", "f_index", "sensor_pos", "sensor_los", 
            "stokes_dim", "atmosphere_dim", "iy_space_agenda", 
            "surface_rtprop_agenda",  
            "propmat_clearsky_agenda", "ppath_step_agenda", "p_grid", "lat_grid",
            "lon_grid", "z_field", "refellipsoid", "z_surface", "t_field", 
            "vmr_field", "edensity_field", "cloudbox_limits", "pnd_field", 
            "scat_data_array_mono", "mc_seed", "iy_unit",
            "mc_std_err", "mc_max_time", "mc_max_iter", "mc_min_iter",
            "mc_z_field_is_1D" ),
        GIN(),
        GIN_TYPE(),
        GIN_DEFAULT(),
        GIN_DESC()
        ));
/* Workspace method: Doxygen documentation will be auto-generated */
void MCIPA(Workspace&            ws,
           Vector&               y,
           Index&                mc_iteration_count,
           Vector&               mc_error,
           Tensor3&              mc_points,
           const MCAntenna&      mc_antenna,
           const Vector&         f_grid,
           const Index&          f_index,
           const Matrix&         sensor_pos,
           const Matrix&         sensor_los,
           const Index&          stokes_dim,
           const Index&          atmosphere_dim,
           const Agenda&         iy_space_agenda,
           const Agenda&         surface_rtprop_agenda,
           const Agenda&         propmat_clearsky_agenda, 
           const Agenda&         ppath_step_agenda,
           const Vector&         p_grid,
           const Vector&         lat_grid, 
           const Vector&         lon_grid, 
           const Tensor3&        z_field, 
           const Vector&         refellipsoid, 
           const Matrix&         z_surface,
           const Tensor3&        t_field, 
           const Tensor4&        vmr_field, 
           const Tensor3&        edensity_field, 
           const ArrayOfIndex&   cloudbox_limits, 
           const Tensor4&        pnd_field,
           const ArrayOfSingleScatteringData& scat_data_array_mono,
           const Index&          mc_seed,
           const String&         y_unit,
           //Keyword params
           const Numeric&        std_err,
           const Index&          max_time,
           const Index&          max_iter,
           const Index&          min_iter,
           const Index&          z_field_is_1D,
           const Verbosity&      verbosity)
{
  if( min_iter < 100 )
    { throw runtime_error( "*mc_min_iter* must be >= 100." ); }

  //Check keyword input
  if (max_time<0 && max_iter<0 && std_err<0){
    throw runtime_error( "At least one of std_err, max_time, and max_iter must be positive" );
  }

  if (! (atmosphere_dim == 3))
    {
      throw runtime_error(
                   "For montecarlo, atmosphere_dim must be 3.");
    }

  Ppath ppath_step;
  Rng rng;                      //Random Nuimber generator
  time_t start_time=time(NULL);
  Index N_pt=pnd_field.nbooks();//Number of particle types
  Vector pnd_vec(N_pt); //Vector of particle number densities used at each point
  Vector Z11maxvector;//Vector holding the maximum phase function for each 
  bool anyptype30=is_anyptype30(scat_data_array_mono);
  if (anyptype30)
    {
      findZ11max(Z11maxvector,scat_data_array_mono);
    }
  rng.seed(mc_seed, verbosity);
  bool keepgoing,inside_cloud; // flag indicating whether to stop tracing a photons path
  Numeric g,temperature,albedo,g_los_csc_theta;
  Matrix A(stokes_dim,stokes_dim),Q(stokes_dim,stokes_dim),evol_op(stokes_dim,stokes_dim),
    ext_mat_mono(stokes_dim,stokes_dim),q(stokes_dim,stokes_dim),newQ(stokes_dim,stokes_dim),
    Z(stokes_dim,stokes_dim);
  q=0.0;newQ=0.0;
  mc_iteration_count=0;
  Vector vector1(stokes_dim),abs_vec_mono(stokes_dim),I_i(stokes_dim),Isum(stokes_dim),
    Isquaredsum(stokes_dim);
  y.resize(stokes_dim);
  y=0;
  Index termination_flag=0;
  mc_error.resize(stokes_dim);
  //local versions of workspace
  Matrix local_iy(1,stokes_dim),local_surface_emission(1,stokes_dim),local_surface_los;
  Tensor4 local_surface_rmatrix;
  Vector local_rte_pos(2);
  Vector local_rte_los(2);
  Vector new_rte_los(2);
    Index np;
  mc_points.resize(p_grid.nelem(),lat_grid.nelem(),lon_grid.nelem());
  mc_points=0;
  Isum=0.0;Isquaredsum=0.0;
  const Numeric f_mono = f_grid[f_index];
  Numeric std_err_i;
  bool convert_to_rjbt=false;
  if ( y_unit == "RJBT" )
    { 
      std_err_i=f_grid[0]*f_grid[0]*2*BOLTZMAN_CONST/SPEED_OF_LIGHT/SPEED_OF_LIGHT*
        std_err;
      convert_to_rjbt=true;
    }
  else if ( y_unit == "1" )
    {
      std_err_i=std_err;
    }
  else
    {
      ostringstream os;
      os << "Invalid value for *y_unit*:" << y_unit <<".\n" 
         << "This method allows only the options \"RJBT\" and \"1\".";
      throw runtime_error( os.str() );
    }
      
  //Begin Main Loop
  while (true)
    {
      mc_iteration_count+=1;
      keepgoing=true; //flag indicating whether to continue tracing a photon path
      //Sample a FOV direction
      mc_antenna.draw_los(local_rte_los,rng,sensor_los(0,joker));
      id_mat(Q);
      local_rte_pos=sensor_pos(0,joker);
      I_i=0.0;

      //Obtain a reference propagation path to use for obtaining optical properties
      //for the IPA method.
      Ppath ppath;
      ppath_calc( ws, ppath, ppath_step_agenda, 3, 
                  p_grid, lat_grid, lon_grid, t_field, z_field, vmr_field,
                  edensity_field, Vector(1, f_mono), refellipsoid,
                  z_surface, 0, cloudbox_limits, local_rte_pos, local_rte_los, 
                  0, verbosity );
      
      while (keepgoing)
        {
          //modified path tracing routine for independent pixel approximation
          mcPathTraceIPA( ws,
                         evol_op, abs_vec_mono, temperature, ext_mat_mono, rng,
                         local_rte_pos, local_rte_los, pnd_vec, g,
                         termination_flag, inside_cloud, 
                         propmat_clearsky_agenda, 
                         stokes_dim, f_mono, p_grid,
                         lat_grid, lon_grid, z_field, refellipsoid, z_surface,
                         t_field, vmr_field, cloudbox_limits, pnd_field,
                         scat_data_array_mono, z_field_is_1D, ppath,
                         verbosity );
            
          np=ppath.np;
          mc_points(ppath.gp_p[np-1].idx,ppath.gp_lat[np-1].idx,
                    ppath.gp_lon[np-1].idx)+=1;
          if (termination_flag==1)
            {
              iy_space_agendaExecute( ws, local_iy, Vector(1,f_grid[f_index]),
                                      local_rte_pos, local_rte_los,
                                      iy_space_agenda );
              mult(vector1,evol_op,local_iy(0,joker));
              mult(I_i,Q,vector1);
              I_i/=g;
              keepgoing=false; //stop here. New photon.
            }
          else if (termination_flag==2)
            {
              //Choose appropriate lat and lon grid points for IPA 
              //if ppath (the reference ppath) hits the srface just take the last point.
              GridPos latgp;
              GridPos longp;
              if (ppath.background=="surface")
                {
                  latgp=ppath.gp_lat[ppath.np-1];
                  longp=ppath.gp_lon[ppath.np-1];
                }
              else
                {
                  //Use lat and lon at the lowest z (ie. the tangent point)
                  Numeric latt, lont, zmin=9e99;
                  for( Index i=0; i<ppath.np; i++ )
                    {
                      if( ppath.pos(i,0) < zmin )
                        {
                          zmin = ppath.pos(i,0);
                          latt = ppath.pos(i,1);
                          lont = ppath.pos(i,2);
                        }
                    }
                  gridpos(latgp,lat_grid,latt);
                  gridpos(longp,lon_grid,lont);
                }
              //decide whether we have reflection or emission
              surface_rtprop_agendaExecute(ws,
                          local_surface_emission, local_surface_los, 
                          local_surface_rmatrix, Vector(1,f_grid[f_index]),
                          local_rte_pos, local_rte_los, surface_rtprop_agenda );
              //deal with blackbody case
              if (local_surface_los.nrows()==0)
                {
                  mult(vector1,evol_op,local_surface_emission(0,joker));
                  mult(I_i,Q,vector1);
                  I_i/=g;
                  keepgoing=false;
                }
              else
                //decide between reflection and emission
                {
                  Numeric R11=local_surface_rmatrix(0,0,0,0);
                  if (rng.draw()>R11)
                    {
                      //then we have emission
                      //Matrix oneminusR(stokes_dim,stokes_dim);
                      //id_mat(oneminusR);
                      //oneminusR-=local_surface_rmatrix(0,0,joker,joker);
                      //oneminusR/=1-R11;
                      //mult(vector1,oneminusR,local_surface_emission(0,joker));
                      mult(vector1,evol_op,local_surface_emission(0,joker));
                      mult(I_i,Q,vector1);
                      I_i/=g*(1-R11);
                      keepgoing=false;
                    }
                  else
                    {
                      //we have reflection
                      local_rte_los=local_surface_los(0,joker);
                      
                      mult(q,evol_op,local_surface_rmatrix(0,0,joker,joker));
                      mult(newQ,Q,q);
                      Q=newQ;
                      Q/=g*R11;
                    }
                }
            }
          else if (inside_cloud)
            {
              //we have another scattering/emission point 
              //Estimate single scattering albedo
              albedo=1-abs_vec_mono[0]/ext_mat_mono(0,0);
              //cout<<"albedo = "<<albedo<<" ext_mat_mono(0,0) = "<<ext_mat_mono(0,0)<<" abs_vec_mono[0] = "<<abs_vec_mono[0]<<"\n";
              //determine whether photon is emitted or scattered
              if (rng.draw()>albedo)
                {
                  //Calculate emission
                  Numeric planck_value = planck( f_grid[0], temperature );
                  Vector emission=abs_vec_mono;
                  emission*=planck_value;
                  Vector emissioncontri(stokes_dim);
                  mult(emissioncontri,evol_op,emission);
                  emissioncontri/=(g*(1-albedo));//yuck!
                  mult(I_i,Q,emissioncontri);
                  keepgoing=false;
                  //cout << "emission contri" <<  I_i[0] << "\n";
                }
              else
                {
                  //we have a scattering event
                  //Sample new line of sight.
                  
                  Sample_los (new_rte_los,g_los_csc_theta,Z,rng,local_rte_los,
                              scat_data_array_mono,stokes_dim,
                              pnd_vec,anyptype30,Z11maxvector,ext_mat_mono(0,0)-abs_vec_mono[0],temperature,
                              verbosity);
                                           
                  Z/=g*g_los_csc_theta*albedo;
                  
                  mult(q,evol_op,Z);
                  mult(newQ,Q,q);
                  Q=newQ;
                  //scattering_order+=1;
                  local_rte_los=new_rte_los;
                  //if (silent==0){cout <<"mc_iteration_count = "<<mc_iteration_count << 
                  //                 ", scattering_order = " <<scattering_order <<"\n";}
                }
            }
          else
            {
              //Must be clear sky emission point
              //Calculate emission
              Numeric planck_value = planck( f_grid[0], temperature );
              Vector emission=abs_vec_mono;
              emission*=planck_value;
              Vector emissioncontri(stokes_dim);
              mult(emissioncontri,evol_op,emission);
              emissioncontri/=g;
              mult(I_i,Q,emissioncontri);
              keepgoing=false;
              //cout << "emission contri" <<  I_i[0] << "\n";
            }
        }
      Isum += I_i;
      for(Index j=0; j<stokes_dim; j++)
        {
          assert(!std::isnan(I_i[j]));
          Isquaredsum[j] += I_i[j]*I_i[j];
        }
      y=Isum;
      y/=(Numeric)mc_iteration_count;
      for(Index j=0; j<stokes_dim; j++) 
        {
          mc_error[j]=sqrt((Isquaredsum[j]/(Numeric)mc_iteration_count-y[j]*y[j])/(Numeric)mc_iteration_count);
        }
      if (std_err>0 && mc_iteration_count>=min_iter && mc_error[0]<std_err_i){break;}
      if (max_time>0 && (Index)(time(NULL)-start_time)>=max_time){break;}
      if (max_iter>0 && mc_iteration_count>=max_iter){break;}
    }
  if ( convert_to_rjbt )
    {
      for(Index j=0; j<stokes_dim; j++) 
        {
          y[j]=invrayjean(y[j],f_grid[0]);
          mc_error[j]=invrayjean(mc_error[j],f_grid[0]);
        }
    }
}


  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_matExtractManually" ),
        DESCRIPTION
        (
         "A simple function for manually extract a single phase matrix.\n"
         "\n"
         "The function returns the phase matrix for a single particle, for\n"
         "scattering from (za_in,aa_in) to (za_out,aa_out).\n"
         "\n"
         "Only a single particle type is handled and *scat_data_array* must\n"
         "have length 1. The frequency is selected by *f_grid* and *f_index*.\n"
         "The temperature is set by *rtp_temperature*.\n"
         ),
        AUTHORS( "Patrick Eriksson" ),
        OUT( ),
        GOUT( "pha_mat_single" ),
        GOUT_TYPE( "Matrix" ),
        GOUT_DESC( 
            "Phase matrix for a single frequency and combination of angles" ),
        IN( "f_grid", "f_index", "stokes_dim", "scat_data_array", 
            "rtp_temperature" ),
        GIN( "za_out", "aa_out", "za_in", "aa_in" ),
        GIN_TYPE( "Numeric", "Numeric", "Numeric", "Numeric" ), 
        GIN_DEFAULT( NODEF, NODEF, NODEF, NODEF ),
        GIN_DESC( "Outgoing zenith angle", "Outgoing azimuth angle",
                  "Incoming zenith angle", "Incoming azimuth angle" )
        ));

/* Workspace method: Doxygen documentation will be auto-generated */
void pha_matExtractManually(Matrix&         pha_mat,
                            const Vector&   f_grid,
                            const Index&    f_index,
                            const Index&    stokes_dim,
                            const ArrayOfSingleScatteringData&   scat_data_array,
                            const Numeric&  rtp_temperature,
                            const Numeric&  za_out, 
                            const Numeric&  aa_out, 
                            const Numeric&  za_in, 
                            const Numeric&  aa_in,
                            const Verbosity& verbosity) 
{
  if( scat_data_array.nelem() > 1 )
    throw runtime_error( "Only one element in *scat_data_array* is allowed." );
  
  ArrayOfSingleScatteringData   scat_data;
  scat_data_array_monoCalc( scat_data, scat_data_array, f_grid, f_index, verbosity);
  
  Vector pnd( 1, 1.0 );
  
  pha_mat.resize( stokes_dim, stokes_dim );
  
  pha_mat_singleCalc( pha_mat, za_out, aa_out, za_in, aa_in,
                      scat_data, stokes_dim, pnd, rtp_temperature,
                      verbosity );
}


//! cum_l_stepCalc

/*!
   Returns a vector of cumulative pathlengths for a given ppath by adding
   up ppath.lstep

   \param cum_l_step    Output: vector of cumulative pathlengths.
   \param ppath         a Ppath.

   \author Cory Davis
   \date   2003-06-19
*/
void cum_l_stepCalc(
                    Vector& cum_l_step,
                    const Ppath& ppath )
{
  cum_l_step.resize(ppath.np);
  Numeric cumsum = 0.0;
  cum_l_step[0]  = 0.0;

  for (Index i=0; i<ppath.np-1; i++)
    {
      cumsum += ppath.lstep[i];
      cum_l_step[i+1] = cumsum;
    }
}             


//! mcPathTraceIPA
/*!
Performs the tasks of pathlength sampling,
ray tracing (but now only as far as determined by pathlength
sampling) and calculation of the evolution operator and several 
atmospheric variables at the new point.  This is the same as mcPathTraceGeneral
modified for the independent pixel approximation.

\author Cory Davis
\date 2005-2-21

*/

void mcPathTraceIPA(Workspace&            ws,
                    MatrixView            evol_op,
                    Vector&               abs_vec_mono,
                    Numeric&              temperature,
                    MatrixView            ext_mat_mono,
                    Rng&                  rng,
                    Vector&               rte_pos,
                    Vector&               rte_los,
                    Vector&               pnd_vec,
                    Numeric&              g,
                    Index&                termination_flag,
                    bool&                 inside_cloud,
                    const Agenda&         propmat_clearsky_agenda,
                    const Index           stokes_dim,
                    const Numeric&        f_mono,
                    const Vector&         p_grid,
                    const Vector&         lat_grid,
                    const Vector&         lon_grid,
                    const Tensor3&        z_field,
                    const Vector&         refellipsoid,
                    const Matrix&         z_surface,
                    const Tensor3&        t_field,
                    const Tensor4&        vmr_field,
                    const ArrayOfIndex&   cloudbox_limits,
                    const Tensor4&        pnd_field,
                    const ArrayOfSingleScatteringData& scat_data_array_mono,
                    const Index           z_field_is_1D,
                    const Ppath&          ppath,
                    const Verbosity&      verbosity)
{ 

  //Internal declarations
  ArrayOfMatrix evol_opArray(2);
  ArrayOfMatrix ext_matArray(2);
  ArrayOfVector abs_vecArray(2),pnd_vecArray(2);
  GridPos gp_p,gp_lat,gp_lon;
  Vector    itw(4);
  Vector tArray(2);
  Vector cum_l_step;
  Vector z_grid(p_grid.nelem());
  Matrix T(stokes_dim,stokes_dim);
  Numeric k;
  Numeric x,y,z;
  Numeric ds,lstep=0.,dx,dy,dz;
  Index   istep = 0;            // Counter for number of steps
  Matrix opt_depth_mat(stokes_dim,stokes_dim),incT(stokes_dim,stokes_dim,0.0);
  Matrix old_evol_op(stokes_dim,stokes_dim);
  Numeric rv_ellips=0.;
  Numeric rv_surface=0.;
  Numeric alt;
  Numeric gpnum;
  const Index np_ipa=ppath.np;
  Index i_closest=0;
  Numeric gp_diff;
  Vector lon_iw(2),lat_iw(2);
  Vector l_vec(2); 
  GridPos gp;
  Vector itw1d(2);
        
  //Initialisations
  
  if (z_field_is_1D)
    {
      z_grid=z_field(joker,0,0);
      rv_ellips=refellipsoid[0];
      rv_surface=rv_ellips+z_surface(0,0);
    }

  id_mat(evol_op);
  evol_opArray[1]=evol_op;
  //initialise Ppath with ppath_start_stepping
  //
  Ppath ppath_step;
  //
  ppath_start_stepping( ppath_step, 3, p_grid, lat_grid, 
                        lon_grid, z_field, refellipsoid, z_surface,
                        0, cloudbox_limits, false, 
                        rte_pos, rte_los, verbosity );

  gp_p=ppath_step.gp_p[0];
  gp_lat=ppath_step.gp_lat[0];
  gp_lon=ppath_step.gp_lon[0];
  rte_pos=ppath_step.pos(0,joker);
  rte_los=ppath_step.los(0,joker);

  // For simplicity, rte_pos holds the radius (not the altitude) until end of
  // function 
  const Numeric rre = refell2d( refellipsoid, lat_grid, gp_lat ); 
  rte_pos[0] += rre;

  Range p_range(cloudbox_limits[0], 
                cloudbox_limits[1]-cloudbox_limits[0]+1);
  Range lat_range(cloudbox_limits[2], 
                cloudbox_limits[3]-cloudbox_limits[2]+1);
  Range lon_range(cloudbox_limits[4], 
                cloudbox_limits[5]-cloudbox_limits[4]+1);
  termination_flag=0;

  inside_cloud=is_inside_cloudbox( ppath_step, cloudbox_limits,false );
  
  if (inside_cloud)
    {
      cloudy_rt_vars_at_gp(ws, ext_mat_mono,abs_vec_mono,pnd_vec,temperature,
                           propmat_clearsky_agenda, stokes_dim, f_mono,
                           gp_p, gp_lat, gp_lon,p_grid[p_range], 
                           t_field(p_range,lat_range,lon_range), 
                           vmr_field(joker,p_range,lat_range,lon_range),
                           pnd_field,scat_data_array_mono, cloudbox_limits,rte_los,
                           verbosity);
    }
  else
    {
      clear_rt_vars_at_gp( ws, ext_mat_mono,abs_vec_mono,temperature, 
                           propmat_clearsky_agenda, f_mono,
                           gp_p, gp_lat, gp_lon, p_grid, t_field, vmr_field);
      pnd_vec=0.0;
    }
  ext_matArray[1]=ext_mat_mono;
  abs_vecArray[1]=abs_vec_mono;
  tArray[1]=temperature;
  pnd_vecArray[1]=pnd_vec;
  //draw random number to determine end point
  Numeric r = rng.draw();
  //
  Numeric ppc, lat0, lon0, za0, aa0;
  //
  while ((evol_op(0,0)>r)&(!termination_flag))
    {
      //check we are not in an infinite loop (assert for ease of debugging
      istep++;
      assert(istep<=5000);
      
      //these array variables hold values from the two ends of a path step.
      //Here we make the values for the last point of the previous step, the values
      //for the first poinst of the current step.
      evol_opArray[0]=evol_opArray[1];
      ext_matArray[0]=ext_matArray[1];
      abs_vecArray[0]=abs_vecArray[1];
      tArray[0]=tArray[1];
      pnd_vecArray[0]=pnd_vecArray[1];
      
      //Choose sensible path step length. This is done by considering the
      //dimensions of the current grid cell and the LOS direction. The aim is
      //to move only one grid cell.
      Numeric Dy=(lat_grid[gp_lat.idx+1]-lat_grid[gp_lat.idx])*DEG2RAD*
                  refell2r(refellipsoid,lat_grid[gp_lat.idx]);
      Numeric Dx=(lon_grid[gp_lon.idx+1]-lon_grid[gp_lon.idx])*
        DEG2RAD*refell2r(refellipsoid,lat_grid[gp_lat.idx])*
        cos(DEG2RAD*rte_pos[1]);
      Numeric Dz=z_field(gp_p.idx+1,gp_lat.idx,gp_lon.idx)-
        z_field(gp_p.idx,gp_lat.idx,gp_lon.idx);
      Numeric Dxy=sqrt(Dx*Dx+Dy*Dy);
      if (Dxy/abs(tan(rte_los[0]*DEG2RAD))>Dz){dz=Dz;}
      else {dz=Dxy/abs(tan(rte_los[0]*DEG2RAD));}
      Numeric dxy;
      if(dz*abs(tan(rte_los[0]*DEG2RAD))>Dxy){dxy=Dxy;}
      else{dxy=dz*abs(tan(rte_los[0]*DEG2RAD));}
      lstep=sqrt(dxy*dxy+dz*dz);
      //calculate new position
      //current position and direction vector
      poslos2cart( x, y, z, dx, dy, dz, rte_pos[0], rte_pos[1], rte_pos[2], 
                   rte_los[0], rte_los[1] );
      lat0 = rte_pos[1];
      lon0 = rte_pos[2];
      za0  = rte_los[0];
      aa0  = rte_los[1];
      ppc  = rte_pos[0]*sin(DEG2RAD*rte_los[0]);
      //new_position
      x+=lstep*dx;
      y+=lstep*dy;
      z+=lstep*dz;
      //back to spherical coords
      cart2poslos(rte_pos[0],rte_pos[1],rte_pos[2],rte_los[0],rte_los[1],
                  x,y,z,dx,dy,dz,ppc,lat0,lon0,za0,aa0);
      //get new grid_positions
      gridpos( gp_lat, lat_grid, rte_pos[1] );
      gridpos( gp_lon, lon_grid, rte_pos[2] );
      if (!z_field_is_1D)
        {
          z_at_latlon( z_grid, p_grid, lat_grid, lon_grid, z_field, gp_lat, gp_lon );
          interpweights( itw, gp_lat, gp_lon );
          rv_ellips  = refell2d(refellipsoid,lat_grid,gp_lat);
          rv_surface = rv_ellips + interp( itw, z_surface, gp_lat, gp_lon );
        }
      alt = rte_pos[0] - rv_ellips;
      //Have we left the top of the atmosphere?
      if ((alt>=z_grid[p_grid.nelem()-1])&(rte_los[0]<=90))
        {
          termination_flag=1;
          //shift relavent position variables to the appropriate point
          //z is linearly interpolated with respect to lat and lon
          //I don't think this is overly important.  Just change gp_p and rte_pos
          alt=z_grid[p_grid.nelem()-1];
          rte_pos[0]=alt+rv_ellips;
        }
      //Have we hit the surface?
      if( (rte_pos[0] <= rv_surface)&(rte_los[0]>=90) )
        {
          termination_flag=2;
          rte_pos[0]=rv_surface;
          alt = rte_pos[0] - rv_ellips;
        }
      gridpos( gp_p, z_grid, alt );

      //now project the new point back to the original IPA path
      //find lat and lon gridpoints on reference ppath_ipa 
      gpnum=fractional_gp(gp_p);
      //search through ppath_ipa to find the closest grid point
      gp_diff=abs(fractional_gp(ppath.gp_p[0])-gpnum);
      for (Index i=1;i<np_ipa;i++)
        {
          Numeric new_gp_diff=abs(fractional_gp(ppath.gp_p[i])-gpnum);
          if (new_gp_diff<=gp_diff)
            {
              gp_diff=new_gp_diff;
              i_closest=i;
            }
          else
            {
              //getting further away - can stop now
              break;
            }
        }
      gp_lat=ppath.gp_lat[i_closest];
      gp_lon=ppath.gp_lon[i_closest];
      interpweights(lon_iw,gp_lon);
      interpweights(lat_iw,gp_lat);
      //      
      if (!z_field_is_1D)
       {
          interpweights( itw, gp_lat, gp_lon );
          rv_ellips  = refell2d(refellipsoid,lat_grid,gp_lat);
        }
      rte_pos[0]=rv_ellips+interp_atmfield_by_gp(3,z_field,
                                                gp_p,gp_lat,gp_lon);
      rte_pos[1]=interp(lat_iw, lat_grid,gp_lat);
      rte_pos[2]=interp(lon_iw, lon_grid,gp_lon);

      //calculate evolution operator
      inside_cloud=is_gp_inside_cloudbox(gp_p,gp_lat,gp_lon,
                                         cloudbox_limits, false );
      //calculate RT variables at new point
      if (inside_cloud)
        {
          cloudy_rt_vars_at_gp(ws,
                               ext_mat_mono, abs_vec_mono, pnd_vec, temperature,
                               propmat_clearsky_agenda, stokes_dim, 
                               f_mono, gp_p,gp_lat, gp_lon, p_grid[p_range],
                               t_field(p_range,lat_range,lon_range), 
                               vmr_field(joker,p_range,lat_range,lon_range),
                               pnd_field, scat_data_array_mono, cloudbox_limits,rte_los,
                               verbosity);
        }
      else
        {
          clear_rt_vars_at_gp( ws, ext_mat_mono, abs_vec_mono, temperature, 
                      propmat_clearsky_agenda, f_mono,
                      gp_p, gp_lat, gp_lon, p_grid, t_field, vmr_field);
          pnd_vec=0.0;
        }
      //put these variables in the last Array slot
      ext_matArray[1]=ext_mat_mono;
      abs_vecArray[1]=abs_vec_mono;
      tArray[1]=temperature;
      pnd_vecArray[1]=pnd_vec;
      opt_depth_mat=ext_matArray[1];
      //now calculate evolution operator
      opt_depth_mat+=ext_matArray[0];
      opt_depth_mat*=-lstep/2;
      incT=0;
      if ( stokes_dim == 1 )
        {
          incT(0,0)=exp(opt_depth_mat(0,0));
        }
      else if ( is_diagonal( opt_depth_mat))
        {
          for ( Index j=0;j<stokes_dim;j++)
            {
              incT(j,j)=exp(opt_depth_mat(j,j));
            }
        }
      else
        {
          matrix_exp_p30(incT,opt_depth_mat);
        }
      mult(evol_op,evol_opArray[0],incT);
      assert(incT(0,0)<1);
      assert(evol_op(0,0)<1);
      evol_opArray[1]=evol_op;

    }
  if (termination_flag)
    {
      //we must have reached the surface or the TOA
      g=evol_op(0,0);
    }
  else
    {
      //then we have met the monte carlo pathlength criteria somewhere inside 
      //the atmosphere and between the two points corresponding to the *Array 
      //variables.

      // If istep is 0 means we never entered the above while loop and
      // thus x, y, z, ds, dx, dy, dz are uninitialized
      assert(istep);

      //find position corresponding to the pathlength criteria
      k=-log(incT(0,0))/lstep;
      //length of path segment required by critera
      ds=log(evol_opArray[0](0,0)/r)/k;
      g=k*r;
      l_vec=0.0;
      l_vec[1]=lstep;
      //gridpos and interpolations for required step length
      gridpos(gp, l_vec, ds);
      interpweights(itw1d,gp);
      //interpolate optical properties at required point
      interp(ext_mat_mono, itw1d, ext_matArray, gp);
      opt_depth_mat = ext_mat_mono;
      opt_depth_mat+=ext_matArray[gp.idx];
      opt_depth_mat*=-ds/2;
      if ( stokes_dim == 1 )
        {
          incT(0,0)=exp(opt_depth_mat(0,0));
        }
      else if ( is_diagonal( opt_depth_mat))
        {
          for ( Index i=0;i<stokes_dim;i++)
            {
              incT(i,i)=exp(opt_depth_mat(i,i));
            }
        }
      else
        {
          //matrix_exp(incT,opt_depth_mat,10);
          matrix_exp_p30(incT,opt_depth_mat);
        }
      mult(evol_op,evol_opArray[gp.idx],incT);
      interp(abs_vec_mono, itw1d, abs_vecArray,gp);
      temperature=interp(itw1d,tArray, gp);
      interp(pnd_vec, itw1d, pnd_vecArray, gp);
      //move cartesion coordinates back to required point.
      x+=(ds-lstep)*dx;
      y+=(ds-lstep)*dy;
      z+=(ds-lstep)*dz;
      //and convert back to spherical.
      cart2poslos(rte_pos[0],rte_pos[1],rte_pos[2],rte_los[0],rte_los[1],
                  x,y,z,dx,dy,dz,ppc,lat0,lon0,za0,aa0); 
    }

  // Go back to altitude in rte_pos
  rte_pos[0] -= rre; 

}



//! iwp_cloud_opt_pathCalc
/*!

Calculated the ice water path(iwp) and cloud optical path (cloud_opt_path) for 
a given sensor_pos_ and sensor_los


\author Cory Davis
\date 2006-6-15

*/
void iwp_cloud_opt_pathCalc(Workspace& ws,
                            Numeric& iwp,
                            Numeric& cloud_opt_path,
                            //input
                            const Vector&         rte_pos,
                            const Vector&         rte_los,
                            const Agenda&         ppath_step_agenda,
                            const Vector&         p_grid,
                            const Vector&         lat_grid, 
                            const Vector&         lon_grid, 
                            const Vector&         refellipsoid, 
                            const Matrix&         z_surface,
                            const Tensor3&        z_field, 
                            const Tensor3&        t_field, 
                            const Tensor4&        vmr_field, 
                            const Tensor3&        edensity_field, 
                            const Numeric&        f_mono,
                            const ArrayOfIndex&   cloudbox_limits, 
                            const Tensor4&        pnd_field,
                            const ArrayOfSingleScatteringData& scat_data_array_mono,
                            const Matrix&         particle_masses,
                            const Verbosity&      verbosity)
{
  //internal declarations
  Ppath ppath;
  Vector local_rte_pos=rte_pos;
  Vector local_rte_los=rte_los;
  iwp=0;
  cloud_opt_path=0;
  //calculate ppath to cloudbox boundary
  ppath_calc( ws, ppath, ppath_step_agenda, 3, 
              p_grid, lat_grid, lon_grid, t_field, z_field, vmr_field,
              edensity_field, Vector(1, f_mono), refellipsoid,
              z_surface, 1, cloudbox_limits, local_rte_pos, local_rte_los, 0,
              verbosity );
  //if this ppath hit a cloud, now take ppath inside cloud
  if (ppath_what_background(ppath)>2)
    {
      //move to the intersection of the cloudbox boundary and the 
      //propagation path
      local_rte_pos=ppath.pos(ppath.np-1,joker);
      local_rte_los=ppath.los(ppath.np-1,joker);

      Range p_range(cloudbox_limits[0], 
                    cloudbox_limits[1]-cloudbox_limits[0]+1);
      Range lat_range(cloudbox_limits[2], 
                      cloudbox_limits[3]-cloudbox_limits[2]+1);
      Range lon_range(cloudbox_limits[4], 
                      cloudbox_limits[5]-cloudbox_limits[4]+1);

      ppath_calc( ws, ppath, ppath_step_agenda, 3, 
                  p_grid, lat_grid, lon_grid, t_field, z_field, vmr_field,
                  edensity_field, Vector(1, f_mono), refellipsoid,
                  z_surface, 1, cloudbox_limits, local_rte_pos, local_rte_los, 
                  1, verbosity );

      Matrix  pnd_ppath(particle_masses.nrows(),ppath.np);
      Vector t_ppath(ppath.np);
      Vector   p_ppath(ppath.np);
      Matrix   vmr_ppath(vmr_field.nbooks(),ppath.np);
      Matrix ext_mat_part(1, 1, 0.0);
      Vector abs_vec_part(1, 0.0);

      cloud_atm_vars_by_gp(p_ppath,t_ppath,vmr_ppath,pnd_ppath,ppath.gp_p,
                           ppath.gp_lat,ppath.gp_lon,cloudbox_limits,p_grid[p_range],
                           t_field(p_range,lat_range,lon_range),
                           vmr_field(joker,p_range,lat_range,lon_range),
                           pnd_field);

      Vector k_vec(ppath.np);
      Vector iwc_vec(ppath.np);
      //calculate integrands
      for (Index i = 0; i < ppath.np ; ++i)
        {
          opt_propCalc(ext_mat_part,abs_vec_part,ppath.los(i,0),ppath.los(i,1),scat_data_array_mono,
                       1, pnd_ppath(joker, i), t_ppath[i], verbosity);
          k_vec[i]=ext_mat_part(0,0);
          Vector pnd_vec=pnd_ppath(joker, i);
          assert(particle_masses.ncols()==1);
          assert(pnd_vec.nelem()==particle_masses.nrows());
          iwc_vec[i]=pnd_vec*particle_masses(joker,0);//hopefully this is the dot product
        }
      //integrate IWP and optical properties
      for (Index i = 0; i < ppath.np-1 ; ++i)
        {
          //Add to path integrals
          iwp+=0.5*(iwc_vec[i+1]+iwc_vec[i])*ppath.lstep[i];
          cloud_opt_path+=0.5*(k_vec[i+1]+k_vec[i])*ppath.lstep[i];
        }
    }
}



//! interpTarray

/*!
   Interpolates several arrays calculated by TarrayCalc to give values at a 
   given pathlength

   \param[out]  T             transmittance matrix ( I may have made this term up ).
   \param[out]  K_abs         absorption coefficient vector
   \param[out]  temperature   temperature
   \param[out]  K             extinction matrix at interpolation point
   \param[out]  rte_pos       position at pathlength along ppath
   \param[out]  rte_los       LOS at pathlength along ppath
   \param[in]   pnd_vec       pnd vector
   \param[in]   TArray        array of transmittance matrices
   \param[in]   ext_matArray  array of extinction matrices
   \param[in]   abs_vecArray  array of absorption coefficients
   \param[in]   t_ppath       array of temperatures
   \param[in]   pnd_ppath     array of pressures
   \param[in]   cum_l_step    vector of cumulative pathlengths
   \param[in]   pathlength    pathlength at which to calculate above values
   \param[in]   stokes_dim    length of Stokes vector
   \param[in]   ppath         the Ppath


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
  Matrix incT(stokes_dim,stokes_dim,0.0);
  Matrix opt_depth_mat(stokes_dim,stokes_dim);
  Vector itw(2);
  Numeric delta_s;
  Index N_pt=pnd_vec.nelem();
  ArrayOfGridPos gp(1);
                  
  //interpolate transmittance matrix
  gridpos(gp, cum_l_step, pathlength);
  interpweights(itw,gp[0]);
  interp(K, itw,ext_matArray,gp[0]);
  delta_s = pathlength - cum_l_step[gp[0].idx];
  opt_depth_mat = K;
  opt_depth_mat+=ext_matArray[gp[0].idx];
  opt_depth_mat*=-delta_s/2;
  if ( stokes_dim == 1 )
    {
      incT(0,0)=exp(opt_depth_mat(0,0));
    }
  else if ( is_diagonal( opt_depth_mat))
    {
      for ( Index i=0;i<stokes_dim;i++)
        {
          incT(i,i)=exp(opt_depth_mat(i,i));
        }
    }
  else
    {
      //matrix_exp(incT,opt_depth_mat,10);
      matrix_exp_p30(incT,opt_depth_mat);
    }
  mult(T,TArray[gp[0].idx],incT);
  
  interp(K_abs, itw, abs_vecArray,gp[0]);
 
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



void mcPathTraceGeneral(
         Workspace&      ws,
         MatrixView      evol_op,
         Vector&         abs_vec_mono,
         Numeric&        temperature,
         MatrixView      ext_mat_mono,
         Rng&            rng,
         Vector&         rte_pos,
         Vector&         rte_los,
         Vector&         pnd_vec,
         Numeric&        g,
         Ppath&          ppath_step,
         Index&          termination_flag,
         bool&           inside_cloud,
   const Agenda&         ppath_step_agenda,
   const Numeric&        ppath_lraytrace,
   const Agenda&         propmat_clearsky_agenda,
   const Index           stokes_dim,
   const Numeric&        f_mono,
   const Vector&         p_grid,
   const Vector&         lat_grid,
   const Vector&         lon_grid,
   const Tensor3&        z_field,
   const Vector&         refellipsoid,
   const Matrix&         z_surface,
   const Tensor3&        t_field,
   const Tensor4&        vmr_field,
   const ArrayOfIndex&   cloudbox_limits,
   const Tensor4&        pnd_field,
   const ArrayOfSingleScatteringData& scat_data_array_mono,
   const Verbosity&      verbosity )
{ 
  ArrayOfMatrix evol_opArray(2);
  ArrayOfMatrix ext_matArray(2);
  ArrayOfVector abs_vecArray(2);
  ArrayOfVector pnd_vecArray(2);
  Matrix        opt_depth_mat(stokes_dim,stokes_dim);
  Matrix        incT(stokes_dim,stokes_dim,0.0);
  Vector        tArray(2);
  Vector        f_grid(1,f_mono);     // Vector version of f_mono
  Matrix        T(stokes_dim,stokes_dim);
  Numeric       k;
  Numeric       ds, dl=-999;
  Index         istep = 0;  // Counter for number of steps
  Matrix        old_evol_op(stokes_dim,stokes_dim);

  CREATE_OUT0;

  //at the start of the path the evolution operator is the identity matrix
  id_mat(evol_op);
  evol_opArray[1]=evol_op;

  // range defining cloudbox
  Range p_range(   cloudbox_limits[0], 
                   cloudbox_limits[1]-cloudbox_limits[0]+1);
  Range lat_range( cloudbox_limits[2], 
                   cloudbox_limits[3]-cloudbox_limits[2]+1 );
  Range lon_range( cloudbox_limits[4], 
                   cloudbox_limits[5]-cloudbox_limits[4]+1 );

  //initialise Ppath with ppath_start_stepping
  ppath_start_stepping( ppath_step, 3, p_grid, lat_grid, 
                        lon_grid, z_field, refellipsoid, z_surface,
                        0, cloudbox_limits, false, 
                        rte_pos, rte_los, verbosity );

  // Index in ppath_step of end point considered presently
  Index ip = 0;

  // Is point ip inside the cloudbox?
  inside_cloud = is_gp_inside_cloudbox( ppath_step.gp_p[ip], 
                                        ppath_step.gp_lat[ip], 
                                        ppath_step.gp_lon[ip], 
                                        cloudbox_limits, 0, 3 );

  // Determine radiative properties at point
  if( inside_cloud )
    {
      cloudy_rt_vars_at_gp( ws, ext_mat_mono, abs_vec_mono, pnd_vec, 
                            temperature, propmat_clearsky_agenda, 
                            stokes_dim, f_mono, ppath_step.gp_p[0], 
                            ppath_step.gp_lat[0], ppath_step.gp_lon[0],
                            p_grid[p_range], 
                            t_field(p_range,lat_range,lon_range), 
                            vmr_field(joker,p_range,lat_range,lon_range),
                            pnd_field, scat_data_array_mono, cloudbox_limits,
                            ppath_step.los(0,joker), verbosity );
    }
  else
    {
      clear_rt_vars_at_gp( ws, ext_mat_mono, abs_vec_mono, temperature, 
                           propmat_clearsky_agenda, f_mono,
                           ppath_step.gp_p[0], ppath_step.gp_lat[0], 
                           ppath_step.gp_lon[0], p_grid, t_field, vmr_field );
      pnd_vec = 0.0;
    }

  // Move the data to end point containers 
  ext_matArray[1] = ext_mat_mono;
  abs_vecArray[1] = abs_vec_mono;
  tArray[1]       = temperature;
  pnd_vecArray[1] = pnd_vec;

  //draw random number to determine end point
  Numeric r = rng.draw();

  termination_flag=0;
  
  while( (evol_op(0,0)>r) && (!termination_flag) )
    {
      istep++;

      if( istep > 5000 )
        {
          throw runtime_error( "5000 path points have been reached. "
                               "Is this an infinite loop?" );
        }

      evol_opArray[0] = evol_opArray[1];
      ext_matArray[0] = ext_matArray[1];
      abs_vecArray[0] = abs_vecArray[1];
      tArray[0]       = tArray[1];
      pnd_vecArray[0] = pnd_vecArray[1];

      // Shall new ppath_step be calculated?
      if( ip == ppath_step.np-1 ) 
        {
          ppath_step_agendaExecute( ws, ppath_step, ppath_lraytrace, t_field, 
                                    z_field, vmr_field, f_grid, 
                                    ppath_step_agenda );
          //Print( ppath_step, 0, verbosity );
          ip = 1;

          inside_cloud = is_gp_inside_cloudbox( ppath_step.gp_p[ip], 
                                                ppath_step.gp_lat[ip], 
                                                ppath_step.gp_lon[ip], 
                                                cloudbox_limits, 0, 3 );
        }
      else
        { ip++; }

      dl = ppath_step.lstep[ip-1];

      if( inside_cloud )
        {
          cloudy_rt_vars_at_gp( ws, ext_mat_mono, abs_vec_mono, pnd_vec,
                                temperature, propmat_clearsky_agenda, 
                                stokes_dim, f_mono, ppath_step.gp_p[ip],
                                ppath_step.gp_lat[ip], ppath_step.gp_lon[ip],
                                p_grid[p_range], 
                                t_field(p_range,lat_range,lon_range), 
                                vmr_field(joker,p_range,lat_range,lon_range),
                                pnd_field, scat_data_array_mono, cloudbox_limits,
                                ppath_step.los(ip,joker), verbosity );
        }
      else
        {
          clear_rt_vars_at_gp( ws, ext_mat_mono, abs_vec_mono, temperature, 
                               propmat_clearsky_agenda, f_mono,
                               ppath_step.gp_p[ip], ppath_step.gp_lat[ip],
                               ppath_step.gp_lon[ip], p_grid, t_field, 
                               vmr_field );
          pnd_vec = 0.0;
        }

      ext_matArray[1] = ext_mat_mono;
      abs_vecArray[1] = abs_vec_mono;
      tArray[1]       = temperature;
      pnd_vecArray[1] = pnd_vec;
      opt_depth_mat   = ext_matArray[1];
      opt_depth_mat  += ext_matArray[0];
      opt_depth_mat  *= -dl/2;
      incT            = 0;

      if( opt_depth_mat(0,0) < -4 )
        {
          out0 << "WARNING: A MC path step of high optical depth ("
               << abs(opt_depth_mat(0,0)) << ")!\n";
        }

      if( stokes_dim == 1 )
        { incT(0,0) = exp( opt_depth_mat(0,0) ); }
      else if( is_diagonal( opt_depth_mat ) )
        {
          for ( Index j=0;j<stokes_dim;j++)
            { incT(j,j) = exp( opt_depth_mat(j,j) ); }
        }
      else
        { matrix_exp_p30( incT, opt_depth_mat ); }
      mult( evol_op, evol_opArray[0], incT );
      evol_opArray[1] = evol_op;
     
      if( evol_op(0,0)>r )
        {
          // Check whether hit ground or space.
          // path_step_agenda just detects surface intersections, and
          // if TOA is reached requires a special check.
          // But we are already ready if evol_op<=r
          if( ip == ppath_step.np - 1 )
            {
              if( ppath_what_background(ppath_step) )
                { termination_flag = 2; }   //we have hit the surface
              else if( fractional_gp(ppath_step.gp_p[ip]) >= 
                                          (Numeric)(p_grid.nelem() - 1)-1e-3 )
                { termination_flag = 1; }  //we are at TOA
            }
        }
    } // while


  if( termination_flag ) 
    { //we must have reached the cloudbox boundary
      rte_pos = ppath_step.pos(ip,joker);
      rte_los = ppath_step.los(ip,joker);
      g       = evol_op(0,0);
    }
  else
    {
      //find position...and evol_op..and everything else required at the new
      //scattering/emission point
      // GH 2011-09-14: 
      //   log(incT(0,0)) = log(exp(opt_depth_mat(0, 0))) = opt_depth_mat(0, 0)
      //   Avoid loss of precision, use opt_depth_mat directly
      //k=-log(incT(0,0))/cum_l_step[np-1];//K=K11 only for diagonal ext_mat
      k  = -opt_depth_mat(0,0) / dl;
      ds = log( evol_opArray[0](0,0) / r ) / k;
      g  = k*r;
      Vector x(2,0.0);
      //interpolate atmospheric variables required later
      ArrayOfGridPos gp(1);
      x[1] = dl;
      Vector itw(2);
  
      gridpos( gp, x, ds );
      assert( gp[0].idx == 0 );
      interpweights( itw, gp[0] );
      interp( ext_mat_mono, itw, ext_matArray, gp[0] );
      opt_depth_mat  = ext_mat_mono;
      opt_depth_mat += ext_matArray[gp[0].idx];
      opt_depth_mat *= -ds/2;
      if( stokes_dim == 1 )
        { incT(0,0) = exp( opt_depth_mat(0,0) ); }
      else if( is_diagonal( opt_depth_mat) )
        {
          for ( Index i=0;i<stokes_dim;i++)
            { incT(i,i) = exp( opt_depth_mat(i,i) ); }
        }
      else
        { matrix_exp_p30( incT, opt_depth_mat ); }
      mult( evol_op, evol_opArray[gp[0].idx], incT );
      interp( abs_vec_mono, itw, abs_vecArray,gp[0] );
      temperature = interp( itw, tArray, gp[0] );
      interp( pnd_vec, itw, pnd_vecArray, gp[0] );
      for( Index i=0; i<2; i++ )
        {
          rte_pos[i] = interp( itw, ppath_step.pos(Range(ip-1,2),i), gp[0] );
          rte_los[i] = interp( itw, ppath_step.los(Range(ip-1,2),i), gp[0] );
        }
      rte_pos[2] = interp( itw, ppath_step.pos(Range(ip-1,2),2), gp[0] );
    }

  assert(isfinite(g));

  // A dirty trick to avoid copying ppath_step
  const Index np = ip+1;
  ppath_step.np  = np;

}



//! matrix_exp_p30
/*!
When we have p30 particles, and therefore the extinction matrix has a 
block diagonal form with 3 independent elements, the matrix expontential
can be calculated very quickly and exactly using this function.

\param M output matrix
\param A input matrix (must be of the form described above)

\author Cory Davis
\date 2005-3-2
*/
void matrix_exp_p30(MatrixView M,
                    ConstMatrixView A)
{
  Index m=A.nrows();
  assert( A.ncols()==m );
  M=0;
  Numeric a=A(0,0);
  Numeric b=A(0,1);
  M(0,0)=cosh(b);
  //M(1,1)=cosh(b);
  M(1,1)=M(0,0);
  M(0,1)=sinh(b);
  //M(1,0)=sinh(b);
  M(1,0)=M(0,1);
  if ( m>2 )
    {
      Numeric c=A(2,3);
      M(2,2)=cos(c);
      if ( m > 3 )
        {
          M(2,3)=sin(c);
          //M(3,2)=-sin(c);
          M(3,2)=-M(2,3);
          M(3,3)=M(2,2);
          //M(3,3)=cos(c); // Added by GH 2011-06-15 as per e-mail 2011-06-13
        }
    }
  M*=exp(a);    
}



