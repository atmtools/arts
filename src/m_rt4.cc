/* Copyright (C) 2016 Jana Mendrok <jana.mendrok@gmail.com>

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
   USA. 
*/
  
/*!
  \file   m_rt4.cc
  \author Jana Mendrok <jana.mendrok@gmail.com>
  \date   2016-05-24
  
  \brief  Workspace functions related to application of scattering solver RT4.
  
  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include "auto_md.h"
#include <complex.h>
#include "messages.h"
#include "m_xml.h"
#include "rt4.h"


#ifdef ENABLE_RT4
/* Workspace method: Doxygen documentation will be auto-generated */
void RT4Calc( Workspace& ws,
                // WS Output:
                Tensor7& doit_i_field,
                Vector& scat_za_grid,
                Vector& scat_aa_grid,
                Index& f_index,
                ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                // WS Input
                const Index& rt4_is_initialized,
                const Index& atmfields_checked,
                const Index& atmgeom_checked,
                const Index& cloudbox_checked,
                const Index& cloudbox_on,
                const ArrayOfIndex& cloudbox_limits, 
                const Agenda& propmat_clearsky_agenda, 
                const Agenda& opt_prop_part_agenda,
                const Agenda& spt_calc_agenda,
                const Agenda& surface_rtprop_agenda,
                const Tensor4& pnd_field,
                const Tensor3& t_field, 
                const Tensor3& z_field, 
                const Tensor4& vmr_field,
                const Vector& p_grid, 
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const Vector& f_grid,
                const Index& stokes_dim,
                const Index& nstreams,
                const Index& non_iso_inc _U_,
                const String& pfct_method,
                const String& quadtype,
                const Index& add_straight_angles,
                const Index& pfct_aa_grid_size,
                const Numeric& pfct_threshold,
                const Numeric& max_delta_tau,
                const Verbosity& verbosity )
{
  // FIXME: so far surface is implictly assumed at lowest atmospheric level.
  // That should be fixed (using z_surface and allowing other altitudes) at some
  // point.

  const String quad_type = quadtype.toupper();
  Index nhza, nhstreams, nummu;
  check_rt4_input( nhstreams, nhza, nummu,
                   cloudbox_on, rt4_is_initialized, "RT4CalcWithRT4surface",
                   atmfields_checked, atmgeom_checked, cloudbox_checked,
                   nstreams, quad_type, add_straight_angles,
                   pnd_field.ncols(), doit_i_field.npages() );

  // in RT4 mu_values is generally only output. however, we need the values for
  // preparing the single scattering data at these angles. therefore, we
  // calculate them here using RT4's proprietary quadrature methods. They
  // simultaneously provide the quadrature weights, too. We keep them so far,
  // might use them for ensuring proper normalization in the preparation of the
  // single scattering data.
  Vector mu_values(nummu, 0.);
  Vector quad_weights(nummu, 0.);

  get_quad_angles( mu_values, quad_weights, scat_za_grid, scat_aa_grid,
                   quad_type, nhstreams, nhza, nummu );
 
  // Preparing surface setup.
  //

  // Initializing surface related interface-related RT4 interface parameters.
  const Index nf = f_grid.nelem();

  // dummy values for parameters not relevant for this ground_type
  Numeric surface_skin_t;
  Vector ground_albedo(nf, 0.);
  Tensor3 ground_reflec(nf,stokes_dim,stokes_dim, 0.);
  Complex gidef(1,0.);
  ComplexVector ground_index(nf, gidef);

  // parameters that will be updated below
  Tensor5 surf_refl_mat(nf,nummu,stokes_dim,nummu,stokes_dim, 0.);
  Tensor3 surf_emis_vec(nf,nummu,stokes_dim, 0.);

  // for now, surface at lowest atm level. later use z_surface or the like
  // for that.
  const Numeric surf_altitude = z_field(0,0,0);
  //const Numeric surf_altitude = z_surface(0,0);

  surf_optpropCalc( ws, surf_refl_mat, surf_emis_vec,
                    surface_rtprop_agenda,
                    f_grid, scat_za_grid, mu_values, 
                    quad_weights, stokes_dim,
                    surf_altitude );

  run_rt4( ws, doit_i_field,
           f_index, f_grid, p_grid, z_field, t_field, vmr_field, pnd_field,
           scat_data, scat_data_mono,
           propmat_clearsky_agenda, opt_prop_part_agenda, spt_calc_agenda,
           cloudbox_limits, stokes_dim, nummu, nhza,
           "A", surface_skin_t,
           ground_albedo, ground_reflec, ground_index,
           surf_refl_mat, surf_emis_vec,
           quad_type, scat_za_grid, mu_values, quad_weights,
           pfct_method, pfct_aa_grid_size, pfct_threshold,
           max_delta_tau,
           verbosity );

  scat_za_grid_adjust( scat_za_grid, mu_values, nummu );

}


/* Workspace method: Doxygen documentation will be auto-generated */
void RT4CalcWithRT4Surface(
                Workspace& ws,
                // WS Output:
                Tensor7& doit_i_field,
                Vector& scat_za_grid,
                Vector& scat_aa_grid,
                Index& f_index,
                ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                // WS Input
                const Index& rt4_is_initialized,
                const Index& atmfields_checked,
                const Index& atmgeom_checked,
                const Index& cloudbox_checked,
                const Index& cloudbox_on,
                const ArrayOfIndex& cloudbox_limits, 
                const Agenda& propmat_clearsky_agenda, 
                const Agenda& opt_prop_part_agenda,
                const Agenda& spt_calc_agenda,
                const Tensor4& pnd_field,
                const Tensor3& t_field, 
                const Tensor3& z_field, 
                const Tensor4& vmr_field,
                const Vector& p_grid, 
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const Vector& f_grid,
                const Index& stokes_dim,
                //const Vector& scat_za_grid,
                const Numeric& surface_skin_t,
                const Vector& surface_scalar_reflectivity,
                const Tensor3& surface_reflectivity,
                const GriddedField3& surface_complex_refr_index,
                const Index& nstreams,
                const Index& non_iso_inc _U_,
                const String& pfct_method,
                const String& groundtype,
                const String& quadtype,
                const Index& add_straight_angles,
                const Index& pfct_aa_grid_size,
                const Numeric& pfct_threshold,
                const Numeric& max_delta_tau,
                const Verbosity& verbosity )
{
  // FIXME: so far surface is implictly assumed at lowest atmospheric level.
  // That should be fixed (using z_surface and allowing other altitudes) at some
  // point.

  const String quad_type = quadtype.toupper();
  Index nhza, nhstreams, nummu;
  check_rt4_input( nhstreams, nhza, nummu,
                   cloudbox_on, rt4_is_initialized, "RT4CalcWithRT4surface",
                   atmfields_checked, atmgeom_checked, cloudbox_checked,
                   nstreams, quad_type, add_straight_angles,
                   pnd_field.ncols(), doit_i_field.npages() );

  // in RT4 mu_values is generally only output. however, we need the values for
  // preparing the single scattering data at these angles. therefore, we
  // calculate them here using RT4's proprietary quadrature methods. They
  // simultaneously provide the quadrature weights, too. We keep them so far,
  // might use them for ensuring proper normalization in the preparation of the
  // single scattering data.
  Vector mu_values(nummu, 0.);
  Vector quad_weights(nummu, 0.);

  get_quad_angles( mu_values, quad_weights, scat_za_grid, scat_aa_grid,
                   quad_type, nhstreams, nhza, nummu );
 
  // Preparing surface setup.
  //

  // Initializing surface related interface-related RT4 interface parameters.
  const Index nf = f_grid.nelem();
  const String ground_type = groundtype.toupper();

  // dummy values for parameters not relevant for this ground_type
  Tensor5 surf_refl_mat(nf,nummu,stokes_dim,nummu,stokes_dim, 0.);
  Tensor3 surf_emis_vec(nf,nummu,stokes_dim, 0.);

  // parameters that will be updated below
  Vector ground_albedo(nf, 0.);
  Tensor3 ground_reflec(nf,stokes_dim,stokes_dim, 0.);
  Complex gidef(1,0.);
  ComplexVector ground_index(nf, gidef);

  get_rt4surf_props( ground_albedo, ground_reflec, ground_index,
                  f_grid, ground_type, surface_skin_t,
                  surface_scalar_reflectivity, surface_reflectivity,
                  surface_complex_refr_index, stokes_dim );

  run_rt4( ws, doit_i_field,
           f_index, f_grid, p_grid, z_field, t_field, vmr_field, pnd_field,
           scat_data, scat_data_mono,
           propmat_clearsky_agenda, opt_prop_part_agenda, spt_calc_agenda,
           cloudbox_limits, stokes_dim, nummu, nhza,
           ground_type, surface_skin_t,
           ground_albedo, ground_reflec, ground_index,
           surf_refl_mat, surf_emis_vec,
           quad_type, scat_za_grid, mu_values, quad_weights,
           pfct_method, pfct_aa_grid_size, pfct_threshold,
           max_delta_tau,
           verbosity );

  scat_za_grid_adjust( scat_za_grid, mu_values, nummu );

}


#else /* ENABLE_RT4 */

void RT4Calc( Workspace&,
                // WS Output:
                Tensor7&,
                Vector&,
                Vector&,
                Index&,
                ArrayOfArrayOfSingleScatteringData&,
                // WS Input
                const Index&,
                const Index&,
                const Index&,
                const Index&,
                const Index&,
                const ArrayOfIndex&,
                const Agenda&,
                const Agenda&,
                const Agenda&,
                const Agenda&,
                const Tensor4&,
                const Tensor3&, 
                const Tensor3&, 
                const Tensor4&,
                const Vector&, 
                const ArrayOfArrayOfSingleScatteringData&,
                const Vector&,
                const Index&,
                const Index&,
                const Index&,
                const String&,
                const String&,
                const Index&,
                const Index&,
                const Numeric&,
                const Numeric&,
                const Verbosity& )
{
    throw runtime_error ("This version of ARTS was compiled without RT4 support.");
}

void RT4CalcWithRT4Surface( Workspace&,
                // WS Output:
                Tensor7&,
                Vector&,
                Vector&,
                Index&,
                ArrayOfArrayOfSingleScatteringData&,
                // WS Input
                const Index&,
                const Index&,
                const Index&,
                const Index&,
                const Index&,
                const ArrayOfIndex&,
                const Agenda&,
                const Agenda&,
                const Agenda&,
                const Tensor4&,
                const Tensor3&, 
                const Tensor3&, 
                const Tensor4&,
                const Vector&, 
                const ArrayOfArrayOfSingleScatteringData&,
                const Vector&,
                const Index&,
                const Numeric&,
                const Vector&,
                const Tensor3& ,
                const GriddedField3&,
                const Index&,
                const Index&,
                const String&,
                const String&,
                const String&,
                const Index&,
                const Index&,
                const Numeric&,
                const Numeric&,
                const Verbosity& )
{
    throw runtime_error ("This version of ARTS was compiled without RT4 support.");
}

#endif /* ENABLE_RT4 */


/* Workspace method: Doxygen documentation will be auto-generated */
void RT4Init(//WS Output
              Tensor7& doit_i_field,
              Index& rt4_is_initialized,
              // WS Input
              const Index& stokes_dim,
              const Index& atmosphere_dim,
              const Vector& f_grid,
              const Index& cloudbox_on,
              const ArrayOfIndex& cloudbox_limits,
              const ArrayOfArrayOfSingleScatteringData& scat_data,
              const Index& nstreams,
              const String& quad_type,
              const Index& add_straight_angles,
              const Verbosity& verbosity _U_ )
{
  if (!cloudbox_on)
  {
    //CREATE_OUT0;
    //rt4_is_initialized = 0;
    //out0 << "  Cloudbox is off, scattering calculation will be skipped.\n";
    //return;
    throw runtime_error( "Cloudbox is off, no scattering calculations to be"
                         "performed." );
  }
  
  // -------------- Check the input ------------------------------
  
  if( atmosphere_dim != 1   )
    throw runtime_error( "For running RT4, atmospheric dimensionality "
                         "must be 1.\n");

  if (stokes_dim < 0 || stokes_dim > 2)
    throw runtime_error( "For running RT4, the dimension of stokes vector "
                         "must be 1 or 2.\n");

  if( cloudbox_limits[0] != 0   )
    {
      ostringstream os;
      os << "RT4 calculations currently only possible with "
         << "lower cloudbox limit\n"
         << "at 0th atmospheric level "
         << "(assumes surface there, ignoring z_surface).\n";
      throw runtime_error(os.str());
    }

  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim )
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*");

  if ( scat_data.empty() )
    throw runtime_error(
                         "No single scattering data present.\n"
                         "See documentation of WSV *scat_data* for options to "
                         "make single scattering data available.\n"
                         );

  if( quad_type.length()>1 )
    {
      ostringstream os;
      os << "Input parameter *quad_type* not allowed to contain more than a "
         << "single character.\n"
         << "Yours has " << quad_type.length() << ".\n";
      throw runtime_error(os.str());
    }

  // If quad!='L' we need to add extra angle 0 & 180deg for final scat_za_grid
  Index neza;
  if( quad_type=="D" || quad_type=="G" )
    {
      if( add_straight_angles )
        neza=2;
      else
        neza=0;
    }
  else if( quad_type=="L" )
    {
      neza=0;
    }
  else
    {
      ostringstream os;
      os << "Unknown quadrature type.\n";
      throw runtime_error(os.str());
    }

  // RT4 can only completely or azimuthally randomly oriented particles.
  bool no_p10=true;
  for( Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++ )
    for( Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++ )
      if( scat_data[i_ss][i_se].ptype != PTYPE_TOTAL_RND &&
          scat_data[i_ss][i_se].ptype != PTYPE_AZIMUTH_RND )
        no_p10=false;
  if( !no_p10 )
    {
      ostringstream os;
      os << "RT4 can only handle scattering elements of type "
         << PTYPE_TOTAL_RND << " (" << PTypeToString(PTYPE_TOTAL_RND) << ") and\n"
         << PTYPE_AZIMUTH_RND << " (" << PTypeToString(PTYPE_AZIMUTH_RND) << "),\n"
         << "but at least one element of other type (" << PTYPE_GENERAL
         << "=" << PTypeToString(PTYPE_GENERAL) << ") is present.\n";
      throw runtime_error( os.str() );
    }
    
  //------------- end of checks ---------------------------------------
  
  const Index Nf = f_grid.nelem();
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  //const Index Nza = scat_za_grid.nelem();

  // Resize and initialize radiation field in the cloudbox
  //doit_i_field.resize( Nf, Np_cloud, 1, 1, Nza, 1, 1 );
  doit_i_field.resize( Nf, Np_cloud, 1, 1, nstreams+neza, 1, stokes_dim );
  doit_i_field = NAN;
  
  rt4_is_initialized = 1;
}



/* Workspace method: Doxygen documentation will be auto-generated */
#ifdef ENABLE_RT4
void RT4Test( Tensor4& out_rad,
              const String& datapath,
              const Verbosity& verbosity )
{
    rt4_test( out_rad, datapath, verbosity );
}
#else
void RT4Test( Tensor4&,
              const String&,
              const Verbosity& )
{
    throw runtime_error ("This version of ARTS was compiled without RT4 support.");
}
#endif

