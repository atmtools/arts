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
  \file   rt4.cc
  \author Jana Mendrok <jana.mendrok@gmail.com>
  \date   2016-05-24
  
  \brief  Contains functions related to application of scattering solver RT4.
  
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include <complex.h>
#include "m_xml.h"
#include "rt4.h"
#include "rte.h"

using std::ostringstream;
using std::runtime_error;


//! gas_optpropCalc
/*!
  Calculates layer averaged gaseous extinction (gas_extinct). This variable is
  required as input for the RT4 subroutine.

  \param ws                    Current workspace
  \param gas_extinct           Layer averaged gas extinction for all layers
  \param propmat_clearsky_agenda as the WSA
  \param t_field               as the WSV 
  \param vmr_field             as the WSV 
  \param p_grid                as the WSV 
  \param f_mono                frequency (single entry vector)
  
  \author Jana Mendrok
  \date   2016-08-08
*/
void gas_optpropCalc( Workspace& ws,
                      VectorView gas_extinct,
                      const Agenda& propmat_clearsky_agenda,
                      ConstTensor3View t_field, 
                      ConstTensor4View vmr_field,
                      ConstVectorView p_grid,
                      ConstVectorView f_mono
                    )
{
  // Initialization
  gas_extinct=0.;

  const Index Np = p_grid.nelem();

  assert( gas_extinct.nelem() == Np-1);


  // Local variables to be used in agendas
  Numeric rtp_temperature_local; 
  Numeric rtp_pressure_local;
  Tensor4 propmat_clearsky_local;
  Vector rtp_vmr_local(vmr_field.nbooks());

  const Vector  rtp_temperature_nlte_local_dummy(0);

  // Calculate layer averaged gaseous extinction
  propmat_clearsky_local = 0.;
  for (Index i = 0; i < Np-1; i++)
    {
      rtp_pressure_local = 0.5 * (p_grid[i] + p_grid[i+1]);
      rtp_temperature_local = 0.5 * (t_field(i,0,0) + t_field(i+1,0,0));
     
      // Average vmrs
      for (Index j = 0; j < vmr_field.nbooks(); j++)
        rtp_vmr_local[j] = 0.5 * (vmr_field(j, i, 0, 0) +
                                  vmr_field(j, i+1, 0, 0));
   
      const Vector rtp_mag_dummy(3,0);
      const Vector ppath_los_dummy;

      //FIXME: do this right?
      Tensor3 nlte_dummy;
      // This is right since there should be only clearsky partials
      ArrayOfTensor3 partial_dummy;
      ArrayOfMatrix partial_source_dummy,partial_nlte_dummy;
      propmat_clearsky_agendaExecute(ws,
                                     propmat_clearsky_local,
                                     nlte_dummy,
                                     partial_dummy,
                                     partial_source_dummy,
                                     partial_nlte_dummy,
                                     ArrayOfRetrievalQuantity(0),
                                     f_mono,  // monochromatic calculation
                                     rtp_mag_dummy,ppath_los_dummy,
                                     rtp_pressure_local, 
                                     rtp_temperature_local, 
                                     rtp_temperature_nlte_local_dummy,
                                     rtp_vmr_local,
                                     propmat_clearsky_agenda);  

      //Assuming non-polarized light and only one frequency
      gas_extinct[Np-2-i] = propmat_clearsky_local(joker,0,0,0).sum();
    }  
}


//! par_optpropCalc
/*!
  Calculates layer averaged gaseous extinction (gas_extinct). This variable is
  required as input for the RT4 subroutine.

  \param ws                    Current workspace
  \param emis_vector           Layer averaged particle absorption for all particle layers
  \param extinct_matrix        Layer averaged particle extinction for all particle layers
  \param scatter_matrix        Layer averaged scattering matrix (azimuth mode 0) for all particle layers
  \param scatlayers            Cloud-to-fullAtm layer association
  \param spt_calc_agenda       as the WSA
  \param opt_prop_part_agenda  as the WSA
  \param pnd_field             as the WSV
  \param t_field               as the WSV
  \param cloudbox_limits       as the WSV 
  
  \author Jana Mendrok
  \date   2016-08-08
*/
void par_optpropCalc( Workspace& ws,
                      Tensor4View emis_vector,
                      Tensor5View extinct_matrix,
                      //VectorView scatlayers,
                      const Agenda& spt_calc_agenda,
                      const Agenda& opt_prop_part_agenda,
                      ConstTensor4View pnd_field,
                      ConstTensor3View t_field,
                      const ArrayOfIndex& cloudbox_limits,
                      const Index& stokes_dim,
                      const Index& nummu
                    )
{
  // Initialization
  extinct_matrix=0.;
  emis_vector=0.;
  
  const Index N_se = pnd_field.nbooks();
  const Index Np_cloud = pnd_field.npages();

  assert( emis_vector.nbooks() == Np_cloud-1);
  assert( extinct_matrix.nshelves() == Np_cloud-1);

  // Local variables to be used in agendas
  Matrix abs_vec_spt_local(N_se, stokes_dim, 0.);
  Tensor3 ext_mat_spt_local(N_se, stokes_dim, stokes_dim, 0.);
  Matrix abs_vec_local;
  Tensor3 ext_mat_local;
  Numeric rtp_temperature_local;
  Tensor4 ext_vector(Np_cloud, 2*nummu, stokes_dim, stokes_dim, 0.);
  Tensor3 abs_vector(Np_cloud, 2*nummu, stokes_dim, 0.);

  // Calculate ext_mat and abs_vec for all pressure points in cloudbox 
  for (Index scat_p_index_local = 0;
             scat_p_index_local < Np_cloud; 
             scat_p_index_local ++)
    {
      rtp_temperature_local =
        t_field(scat_p_index_local+cloudbox_limits[0], 0, 0);

      for (Index iza=0; iza<2*nummu; iza++)
        {
          //Calculate optical properties for all individual scattering elements:
          spt_calc_agendaExecute(ws,
                                 ext_mat_spt_local, 
                                 abs_vec_spt_local,
                                 scat_p_index_local, 0, 0, //position
                                 rtp_temperature_local,
                                 iza, 0, // angles, only needed for aa=0
                                 spt_calc_agenda);

          opt_prop_part_agendaExecute(ws,
                                      ext_mat_local, abs_vec_local, 
                                      ext_mat_spt_local, 
                                      abs_vec_spt_local,
                                      scat_p_index_local, 0, 0, 
                                      opt_prop_part_agenda);

          ext_vector(scat_p_index_local,iza,joker,joker) =
            ext_mat_local(0,joker,joker);
          abs_vector(scat_p_index_local,iza,joker) = abs_vec_local(0,joker);
        }
    }

  // Calculate layer averaged single scattering albedo and layer optical depth
  for (Index scat_p_index_local = 0;
             scat_p_index_local < Np_cloud-1; 
             scat_p_index_local ++)
    {
/*
      if ( (ext_vector(scat_p_index_local,0,0,0)+
            ext_vector(scat_p_index_local+1,0,0,0)) > 0. )
        {
          scatlayers[Np_cloud-2-cloudbox_limits[0]-scat_p_index_local] =
            float(scat_p_index_local);
*/
          for (Index imu=0; imu<nummu; imu++)
            for (Index ist1=0; ist1<stokes_dim; ist1++)
              {
                for (Index ist2=0; ist2<stokes_dim; ist2++)
                  {
                    extinct_matrix(scat_p_index_local,0,imu,ist1,ist2) = .5 *
                      ( ext_vector(scat_p_index_local,imu,ist1,ist2) +
                        ext_vector(scat_p_index_local+1,imu,ist1,ist2) );
                    extinct_matrix(scat_p_index_local,1,imu,ist1,ist2) = .5 *
                      ( ext_vector(scat_p_index_local,nummu+imu,ist1,ist2) +
                        ext_vector(scat_p_index_local+1,nummu+imu,ist1,ist2) );
                  }
                emis_vector(scat_p_index_local,0,imu,ist1) = .5 *
                  ( abs_vector(scat_p_index_local,imu,ist1) +
                    abs_vector(scat_p_index_local+1,imu,ist1) );
                emis_vector(scat_p_index_local,1,imu,ist1) = .5 *
                  ( abs_vector(scat_p_index_local,nummu+imu,ist1) +
                    abs_vector(scat_p_index_local+1,nummu+imu,ist1) );
              }
//        }
    }  
}



//! Calculate radiation field using RT4
/*! 
  Calculate radiation field using Evans' RT4 model (part of PolRadTran).

  This is a direct interface to the (almost orignal) RT4 FORTRAN code. No checks
  of input are made. Function is only to be called through other
  functions/methods, which have to ensure input consistency.

  \param[out] name            descript
  \param[in]  name            descript [unit]

  \author Jana Mendrok
  \date 2016-05-24
*/
void rt4_test( Tensor4& out_rad,
               const Verbosity& verbosity )
{
    //emissivity.resize(4);
    //reflectivity.resize(4);

    Index nstokes=2;
    Index nummu=8;
    Numeric max_delta_tau=1.0E-6;
    String quad_type="L";
    Numeric ground_temp=300.;
    String ground_type="L";
    Numeric ground_albedo=0.05;
    Complex ground_index;
    Numeric sky_temp=0.;
    Numeric wavelength=880.;
    //Index noutlevels=1;
    //ArrayOfIndex outlevels(1);
    //outlevels[0]=1;

    Vector height, temperatures, gas_extinct;
    Tensor5 sca_data;
    Tensor4 ext_data;
    Tensor3 abs_data;
    ReadXML( height, "height", "testdata/z.xml", "", verbosity );
    ReadXML( temperatures, "temperatures", "testdata/T.xml", "", verbosity );
    ReadXML( gas_extinct, "gas_extinct", "testdata/abs_gas.xml", "", verbosity );
    ReadXML( abs_data, "abs_data", "testdata/abs_par.xml", "", verbosity );
    ReadXML( ext_data, "ext_data", "testdata/ext_par.xml", "", verbosity );
    ReadXML( sca_data, "sca_data", "testdata/sca_par.xml", "", verbosity );
    Index num_layers=height.nelem()-1;
    Index num_scatlayers=3;
    Vector scatlayers(num_layers,0.);
    scatlayers[3]=1.;
    scatlayers[4]=2.;
    scatlayers[5]=3.;

    // the read in sca/ext/abs_data is the complete set (and it's in the wrong
    // order for passing it directly to radtrano). before handing over to
    // fortran, we need to reduce it to the number of stokes elements to be
    // used. we can't use views here as all data needs to be continuous in
    // memory; that is, we have to explicitly copy the data we need.
    Tensor6 scatter_matrix(num_scatlayers,4,nummu,nstokes,nummu,nstokes);
    for( Index ii=0; ii<4; ii++ )
      for( Index ij=0; ij<nummu; ij++ )
        for( Index ik=0; ik<nstokes; ik++ )
          for( Index il=0; il<nummu; il++ )
            for( Index im=0; im<nstokes; im++ )
              for( Index in=0; in<num_scatlayers; in++ )
                scatter_matrix(in,ii,ij,ik,il,im) = sca_data(im,il,ik,ij,ii);
    Tensor5 extinct_matrix(num_scatlayers,2,nummu,nstokes,nstokes);
    for( Index ii=0; ii<2; ii++ )
      for( Index ij=0; ij<nummu; ij++ )
        for( Index ik=0; ik<nstokes; ik++ )
          for( Index il=0; il<nstokes; il++ )
            for( Index im=0; im<num_scatlayers; im++ )
              extinct_matrix(im,ii,ij,ik,il) = ext_data(il,ik,ij,ii);
    Tensor4 emis_vector(num_scatlayers,2,nummu,nstokes);
    for( Index ii=0; ii<2; ii++ )
      for( Index ij=0; ij<nummu; ij++ )
        for( Index ik=0; ik<nstokes; ik++ )
          for( Index il=0; il<num_scatlayers; il++ )
            emis_vector(il,ii,ij,ik) = abs_data(ik,ij,ii);

    // Output variables
    Vector mu_values(nummu);
    Tensor3 up_rad(num_layers+1,nummu,nstokes, 0.);
    Tensor3 down_rad(num_layers+1,nummu,nstokes, 0.);

    radtrano_( nstokes,
               nummu,
               max_delta_tau,
               quad_type.c_str(),
               ground_temp,
               ground_type.c_str(),
               ground_albedo,
               ground_index,
               sky_temp,
               wavelength,
               num_layers,
               height.get_c_array(),
               temperatures.get_c_array(),
               gas_extinct.get_c_array(),
               num_scatlayers,
               scatlayers.get_c_array(),
               extinct_matrix.get_c_array(),
               emis_vector.get_c_array(),
               scatter_matrix.get_c_array(),
               //noutlevels,
               //outlevels.get_c_array(),
               mu_values.get_c_array(),
               up_rad.get_c_array(),
               down_rad.get_c_array()
             );

    //so far, output is in
    //    units W/m^2 um sr
    //    dimensions [numlayers+1,nummu,nstokes]
    //    sorted from high to low (altitudes) and 0 to |1| (mu)
    //WriteXML( "ascii", up_rad, "up_rad.xml", 0, "up_rad", "", "", verbosity );
    //WriteXML( "ascii", down_rad, "down_rad.xml", 0, "down_rad", "", "", verbosity );
    
    //to be able to compare with RT4 reference results, reshape output into
    //RT4-output type table (specifically, resort up_rad such that it runs from
    //zenith welling to horizontal, thus forms a continuous angle grid with
    //down_rad. if later changing up_rad/down_rad sorting such that it is in
    //line with doit_i_field, then this has to be adapted as well...
    out_rad.resize(num_layers+1,2,nummu,nstokes);
    for( Index ii=0; ii<nummu; ii++ )
      out_rad(joker,0,ii,joker) = up_rad(joker,nummu-1-ii,joker);
    //out_rad(joker,0,joker,joker) = up_rad;
    out_rad(joker,1,joker,joker) = down_rad;
    //WriteXML( "ascii", out_rad, "out_rad.xml", 0, "out_rad", "", "", verbosity );
}

