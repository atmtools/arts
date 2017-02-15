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

#include "config.h"

#ifdef ENABLE_RT4

#include <cfloat>
#include <stdexcept>
#include <complex.h>
#include "interpolation.h"
#include "m_xml.h"
#include "physics_funcs.h"
#include "optproperties.h"
#include "rt4.h"
#include "rte.h"

using std::ostringstream;
using std::runtime_error;

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric SPEED_OF_LIGHT;

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

  assert( gas_extinct.nelem() == Np-1 );


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
  Calculates layer averaged particle extinction and absorption (extinct_matrix
  and emis_vector)). These variables are required as input for the RT4 subroutine.

  \param ws                    Current workspace
  \param emis_vector           Layer averaged particle absorption for all particle layers
  \param extinct_matrix        Layer averaged particle extinction for all particle layers
  \param spt_calc_agenda       as the WSA
  \param opt_prop_part_agenda  as the WSA
  \param pnd_field             as the WSV
  \param t_field               as the WSV
  \param cloudbox_limits       as the WSV 
  \param stokes_dim            as the WSV
  \param nummu                 Number of single hemisphere polar angles
  
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

  assert( emis_vector.nbooks() == Np_cloud-1 );
  assert( extinct_matrix.nshelves() == Np_cloud-1 );

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

          //Calculate bulk optical properties:
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

  // Calculate layer averaged extnction and absorption
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
                    extinct_matrix(scat_p_index_local,0,imu,ist2,ist1) = .5 *
                      ( ext_vector(scat_p_index_local,imu,ist1,ist2) +
                        ext_vector(scat_p_index_local+1,imu,ist1,ist2) );
                    extinct_matrix(scat_p_index_local,1,imu,ist2,ist1) = .5 *
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


//! sca_optpropCalc
/*!
  Calculates layer (and azimuthal) averaged phase matrix (scatter_matrix). This
  variable is required as input for the RT4 subroutine.

  \param scatter_matrix        Layer averaged scattering matrix (azimuth mode 0) for all particle layers
  \param emis_vector           Layer averaged particle absorption for all particle layers
  \param extinct_matrix        Layer averaged particle extinction for all particle layers
  \param scat_data_mono        as the WSV
  \param pnd_field             as the WSV
  \param stokes_dim            as the WSV
  \param scat_za_grid          as the WSV
  \param quad_weights          Quadrature weights associated with scat_za_grid 
  \param pfct_method           Method for scattering matrix temperature dependance handling
  \param pfct_aa_grid_size     Number of azimuthal grid points in Fourier series decomposition of randomly oriented particles
  
  \author Jana Mendrok
  \date   2016-08-08
*/
void sca_optpropCalc( //Output
                      Tensor6View scatter_matrix,
                      //Input
                      ConstTensor4View emis_vector,
                      ConstTensor5View extinct_matrix,
                      const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                      ConstTensor4View pnd_field,
                      const Index& stokes_dim,
                      const Vector& scat_za_grid,
                      ConstVectorView quad_weights,
                      const String& pfct_method,
                      const Index& pfct_aa_grid_size,
                      const Numeric& pfct_threshold,
                      const Verbosity& verbosity )
{
  // FIXME: do we have numerical issues, too, here in case of tiny pnd? check
  // with Patrick's Disort-issue case.

  // Check that we do indeed have scat_data_mono here.
  if( scat_data_mono[0][0].f_grid.nelem() > 1 )
  {
      ostringstream os;
      os << "Scattering data seems to be scat_data (several freq points),\n"
         << "but scat_data_mono (1 freq point only) is expected here.";
      throw runtime_error( os.str() );
  }

  // Initialization
  scatter_matrix=0.;
  
  const Index N_se = pnd_field.nbooks();
  const Index Np_cloud = pnd_field.npages();
  const Index nza_rt = scat_za_grid.nelem();

  assert( scatter_matrix.nvitrines() == Np_cloud-1 );

  // Check that total number of scattering elements in scat_data and pnd_field
  // agree.
  if( TotalNumberOfElements(scat_data_mono) != N_se )
  {
      ostringstream os;
      os << "Total number of scattering elements in scat_data ("
         << TotalNumberOfElements(scat_data_mono)
         << ") and pnd_field (" << N_se << ") disagree.";
      throw runtime_error( os.str() );
  }

  if( pfct_aa_grid_size < 2 )
  {
      ostringstream os;
      os << "Azimuth grid size for scatt matrix extraction "
         << "(*pfct_aa_grid_size*) must be >1.\n"
         << "Yours is " << pfct_aa_grid_size << ".\n";
      throw runtime_error( os.str() );
  }
  Vector aa_grid;
  nlinspace(aa_grid, 0, 180, pfct_aa_grid_size);

  Index i_se_flat=0;
  Tensor5 sca_mat(N_se,nza_rt,nza_rt,stokes_dim,stokes_dim, 0.);

  // first we extract Z at one T, integrate the azimuth data at each
  // za_inc/za_sca combi (to get the Fourier series 0.th mode), then interpolate
  // to the mu/mu' combis we need in RT.
  for (Index i_ss = 0; i_ss < scat_data_mono.nelem(); i_ss++)
    {
      for (Index i_se = 0; i_se < scat_data_mono[i_ss].nelem(); i_se++)
        {
          SingleScatteringData ssd=scat_data_mono[i_ss][i_se];
          Index i_pfct;
          if( pfct_method=="low" )
            i_pfct = 0;
          else if( pfct_method=="high" )
            i_pfct = ssd.T_grid.nelem()-1;
          else //if( pfct_method=="median" )
            i_pfct = ssd.T_grid.nelem()/2;

          if (ssd.ptype == PTYPE_TOTAL_RND)
            {
              Matrix pha_mat(stokes_dim,stokes_dim, 0.);
              for (Index iza=0; iza<nza_rt; iza++)
                for (Index sza=0; sza<nza_rt; sza++)
                  {
                    Matrix pha_mat_int(stokes_dim,stokes_dim, 0.);
                    for (Index saa=0; saa<pfct_aa_grid_size; saa++)
                      {
                        pha_matTransform( pha_mat(joker,joker),
                                          ssd.pha_mat_data(0,i_pfct,joker,
                                            joker,joker,joker,joker),
                                          ssd.za_grid, ssd.aa_grid, ssd.ptype,
                                          sza, saa, iza, 0,
                                          scat_za_grid, aa_grid,
                                          verbosity );
                        Numeric daa;
                        if (saa==0)
                          daa = (aa_grid[saa+1]-aa_grid[saa])/360.;
                        else if (saa==pfct_aa_grid_size-1)
                          daa = (aa_grid[saa]-aa_grid[saa-1])/360.;
                        else
                          daa = (aa_grid[saa+1]-aa_grid[saa-1])/360.;
                        for (Index ist1=0; ist1<stokes_dim; ist1++)
                          for (Index ist2=0; ist2<stokes_dim; ist2++)
                            pha_mat_int(ist1,ist2) += pha_mat(ist1,ist2) * daa;
                      }
                    sca_mat(i_se_flat,iza,sza,joker,joker) = pha_mat_int;
                  }
            }
          else if (ssd.ptype == PTYPE_AZIMUTH_RND)
            {
              Index nza_se = ssd.za_grid.nelem();
              Index naa_se = ssd.aa_grid.nelem();
              Tensor4 pha_mat_int(nza_se,nza_se,stokes_dim,stokes_dim, 0.);
              ConstVectorView za_datagrid = ssd.za_grid;
              ConstVectorView aa_datagrid = ssd.aa_grid;
              assert( aa_datagrid[0]==0. );
              assert( aa_datagrid[naa_se-1]==180. );

              // first, extracting the phase matrix at the scatt elements own
              // polar angle grid, deriving their respective azimuthal (Fourier
              // series) 0-mode
              for (Index iza=0; iza<nza_se; iza++)
                for (Index sza=0; sza<nza_se; sza++)
                  {
                    for (Index saa=0; saa<naa_se; saa++)
                      {
                        Numeric daa;
                        if (saa==0)
                          daa = (aa_datagrid[saa+1]-aa_datagrid[saa])/360.;
                        else if (saa==naa_se-1)
                          daa = (aa_datagrid[saa]-aa_datagrid[saa-1])/360.;
                        else
                          daa = (aa_datagrid[saa+1]-aa_datagrid[saa-1])/360.;
                        for (Index ist1=0; ist1<stokes_dim; ist1++)
                          for (Index ist2=0; ist2<stokes_dim; ist2++)
                            pha_mat_int(sza,iza,ist1,ist2) +=
                              ssd.pha_mat_data(0,i_pfct,sza,saa,iza,0,ist1*4+ist2) * daa;
                      }
                  }

              // second, interpolating the extracted azimuthal mode to the RT4
              // solver polar angles
              for (Index iza=0; iza<nza_rt; iza++)
                for (Index sza=0; sza<nza_rt; sza++)
                  {
                    GridPos za_sca_gp;
                    GridPos za_inc_gp;
                    Vector itw(4);
                    Matrix pha_mat_lab(stokes_dim,stokes_dim, 0.);
                    Numeric za_sca = scat_za_grid[sza]; 
                    Numeric za_inc = scat_za_grid[iza]; 
       
                    gridpos(za_inc_gp,za_datagrid,za_inc);
                    gridpos(za_sca_gp,za_datagrid,za_sca);
                    interpweights(itw,za_sca_gp,za_inc_gp);
      
                    for (Index ist1=0; ist1<stokes_dim; ist1++)
                      for (Index ist2=0; ist2<stokes_dim; ist2++)
                      {
                        pha_mat_lab(ist1,ist2) = interp(itw,
                          pha_mat_int(Range(joker),Range(joker),ist1,ist2),
                          za_sca_gp,za_inc_gp);
                        //if (ist1+ist2==1)
                        //  pha_mat_lab(ist1,ist2) *= -1.;
                      }

                    sca_mat(i_se_flat,iza,sza,joker,joker) = pha_mat_lab;
                  }
            }
          else
            {
              ostringstream os;
              os << "Unsuitable particle type encountered.";
              throw runtime_error( os.str() );
            }
          i_se_flat++;
        }
    }

  assert( i_se_flat == N_se );
  // now we sum up the Z(mu,mu') over the scattering elements weighted by the
  // pnd_field data, deriving Z(z,mu,mu') and sorting this into
  // scatter_matrix
  // FIXME: it seems that at least for p20, corresponding upward and downward
  // directions have exactly the same (0,0) elements. Can we use this to reduce
  // calc efforts (for sca and its normalization as well as for ext and abs?)
  Index nummu = nza_rt/2;
  for (Index scat_p_index_local = 0;
             scat_p_index_local < Np_cloud-1; 
             scat_p_index_local ++)
    {
      for (Index i_se = 0; i_se < N_se; i_se++)
        {
          Numeric pnd_mean = 0.5 * ( pnd_field(i_se,scat_p_index_local+1,0,0)+
                                     pnd_field(i_se,scat_p_index_local,0,0) );
          if ( pnd_mean > 0. )
            for (Index iza=0; iza<nummu; iza++)
              if ( (extinct_matrix(scat_p_index_local,0,iza,0,0)+
                    extinct_matrix(scat_p_index_local,1,iza,0,0)) > 0. )
                for (Index sza=0; sza<nummu; sza++)
                  for (Index ist1=0; ist1<stokes_dim; ist1++)
                    for (Index ist2=0; ist2<stokes_dim; ist2++)
                    {
                      scatter_matrix(scat_p_index_local,0,iza,ist2,sza,ist1) +=
                        pnd_mean * sca_mat(i_se,iza,sza,ist1,ist2);
                      scatter_matrix(scat_p_index_local,1,iza,ist2,sza,ist1) +=
                        pnd_mean * sca_mat(i_se,nummu+iza,sza,ist1,ist2);
                      scatter_matrix(scat_p_index_local,2,iza,ist2,sza,ist1) +=
                        pnd_mean * sca_mat(i_se,iza,nummu+sza,ist1,ist2);
                      scatter_matrix(scat_p_index_local,3,iza,ist2,sza,ist1) +=
                        pnd_mean * sca_mat(i_se,nummu+iza,nummu+sza,ist1,ist2);
                    }
        }
//      cout << "cloudbox layer #" << scat_p_index_local << "\n";
      for (Index iza=0; iza<nummu; iza++)
        for (Index ih=0; ih<2; ih++)
          if ( extinct_matrix(scat_p_index_local,ih,iza,0,0) > 0. )
            {
              Numeric sca_mat_integ = 0.;

              Numeric ext_nom = extinct_matrix(scat_p_index_local,ih,iza,0,0);
              Numeric sca_nom = ext_nom-emis_vector(scat_p_index_local,ih,iza,0);
              Numeric w0_nom = sca_nom/ext_nom;
              assert( w0_nom>=0. );

              for (Index sza=0; sza<nummu; sza++)
              {
//                SUM2 = SUM2 + 2.0D0*PI*QUAD_WEIGHTS(J2)*
//     .                 ( SCATTER_MATRIX(1,J2,1,J1, L,TSL)
//     .                 + SCATTER_MATRIX(1,J2,1,J1, L+2,TSL) )
                sca_mat_integ += quad_weights[sza] *
                  ( scatter_matrix(scat_p_index_local,ih,iza,0,sza,0)+
                    scatter_matrix(scat_p_index_local,ih+2,iza,0,sza,0) );
              }

              // compare integrated scatt matrix with ext-abs for respective
              // incident polar angle - consistently with scat_dataCheck, we do
              // this in trms of albedo deviation (since PFCT deviations at
              // small albedos matter less than those at high albedos)
//            SUM1 = EMIS_VECTOR(1,J1,L,TSL)-EXTINCT_MATRIX(1,1,J1,L,TSL)
              Numeric pfct_norm = 2.*PI*sca_mat_integ / sca_nom;
              Numeric w0_act = 2.*PI*sca_mat_integ / ext_nom;
              if (abs(w0_act-w0_nom) > pfct_threshold)
              {
                ostringstream os;
                os << "Bulk scattering matrix normalization deviates significantly\n"
                   << "from expected value (" << 1e2*abs(1.-pfct_norm) << "%, "
                   << "resulting in albedo deviation of " << abs(w0_act-w0_nom)
                   << ").\n"
                   << "Something seems wrong with your scattering data "
                   << "(did you run *scat_dataCheck*?)\n"
                   << "or your RT4 setup (try increasing *nstreams* and in case "
                   << "of randomly oriented particles possibly also "
                   << "pfct_aa_grid_size).";
                throw runtime_error( os.str() );
              }
              if (abs(w0_act-w0_nom) > pfct_threshold*0.1 || abs(1.-pfct_norm) > 1e-2)
              {
                CREATE_OUT2;
                out2 << "Warning: The bulk scattering matrix is not well normalized\n"
                     << "Deviating from expected value by " << 1e2*abs(1.-pfct_norm)
                     << "% (and " << abs(w0_act-w0_nom)
                     << " in terms of scattering albedo).\n";
//                cout << "polar angle #" << iza << "." << ih << "\n";
//                cout << "Scattering matrix deviating from expected value by "
//                     << 1e2*abs(1.-norm) << "%.\n";
              }

              // rescale scattering matrix to expected (0,0) value (and scale all
              // other elements accordingly)
              // FIXME: not fully clear whether applying the same rescaling
              // factor is the correct way to do. check out, e.g., Vasilieva
              // (JQSRT 2006) for better approaches.
              scatter_matrix(scat_p_index_local,ih,iza,joker,joker,joker) /=
                pfct_norm;
              scatter_matrix(scat_p_index_local,ih+2,iza,joker,joker,joker) /=
                pfct_norm;
            }
    }
}


//! surf_optpropCalc
/*!
  Calculates bidirectional surface reflection matrices and emission direction
  dependent surface emission terms as required as input for the RT4 subroutine.

  \param surf_refl_mat         Bidirectional surface reflection matrices on RT4 stream directions.
  \param surf_emis_vec         Directional surface emission vector on RT4 stream directions.
  \param extinct_matrix        Layer averaged particle extinction for all particle layers
  \param scat_data_mono        as the WSV
  \param pnd_field             as the WSV
  \param stokes_dim            as the WSV
  \param scat_za_grid          as the WSV
  \param quad_weights          Quadrature weights associated with scat_za_grid 
  \param pfct_method           Method for scattering matrix temperature dependance handling
  \param pfct_aa_grid_size     Number of azimuthal grid points in Fourier series decomposition of randomly oriented particles
  
  \author Jana Mendrok
  \date   2017-02-09
*/
void surf_optpropCalc( Workspace& ws,
                       //Output
                       Tensor5View surf_refl_mat,
                       Tensor3View surf_emis_vec,
                       //Input
                       const Agenda& surface_rtprop_agenda,
                       ConstVectorView f_grid,
                       ConstVectorView scat_za_grid,
                       ConstVectorView mu_values,
                       ConstVectorView quad_weights,
                       const Index& stokes_dim,
                       const Numeric& surf_alt )
{
  // While proprietary RT4 - from the input/user control side - handles only
  // Lambertian and Fresnel, the Doubling&Adding solver core applies a surface
  // reflection matrix and a surface radiance term. The reflection matrix is
  // dependent on incident and reflection direction, ie forms a discrete
  // representation of the bidirectional reflection function; the radiance term
  // is dependent on outgoing polar angle. That is, the solver core can
  // basically handle ANY kind of surface reflection type.
  // 
  // Here, we replace RT4's proprietary surface reflection and radiance term
  // calculation and use ARTS' surface_rtprop_agenda instead to set these
  // variables up.
  // That is, ARTS-RT4 is set up such that it can handle all surface reflection
  // types that ARTS itself can handle.
  // What is required here, is to derive the reflection matrix over all incident
  // and reflected polar angle directions. The surface_rtprop_agenda handles one
  // reflected direction at a time (rtp_los); it is sufficient here to simply
  // loop over all reflected angles as given by the stream directions. However,
  // the incident directions (surface_los) provided by surface_rtprop_agenda
  // depends on the specific, user-defined setup of the agenda. Out of what the
  // agenda call provides, we have to derive the reflection matrix entries for
  // the incident directions as defined by the RT4 stream directions. To be
  // completely general (handling everything from specular via semi-specular to
  // Lambertian) is a bit tricky. Here we decided to allow the ARTS-standard
  // 0.5-grid-spacing extrapolation and set everything outside that range to
  // zero.
  //
  // We do all frequencies here at once (assuming this is the faster variant as
  // the agenda anyway (can) provide output for full f_grid at once and as we
  // have to apply the same inter/extrapolation to all the frequencies).
  //
  // FIXME: Make sure that normalization (energy conservation in the form of
  // power reflection coefficient conservation) is given.
  // FIXME: Allow surface to be elsewhere than at lowest atm level (this
  // requires changes in the surface setting part and more extensive ones in the
  // atmospheric optical property prep part within the frequency loop further
  // below).

  chk_not_empty( "surface_rtprop_agenda", surface_rtprop_agenda );
  
  const Index nf = f_grid.nelem();
  const Index nummu = scat_za_grid.nelem()/2;
  const String B_unit = "R";

  // Local input of surface_rtprop_agenda.
  Vector rtp_pos(1, surf_alt); //atmosphere_dim is 1

  for (Index rmu=0; rmu<nummu; rmu++)
    {
      // Local output of surface_rtprop_agenda.
      Numeric   surface_skin_t;
      Matrix    surface_los;
      Tensor4   surface_rmatrix;
      Matrix    surface_emission;

      // rtp_los is reflected direction, ie upwelling direction, which is >90deg
      // in ARTS. scat_za_grid is sorted here as downwelling (90->0) in 1st
      // half, upwelling (90->180) in second half. that is, here we have to take
      // the second half grid or, alternatively, use 180deg-za[imu].
      Vector rtp_los(1, scat_za_grid[nummu+rmu]);
      //cout << "Doing reflected dir #" << rmu << " at " << rtp_los[0] << " degs\n";

      surface_rtprop_agendaExecute( ws,
                                    surface_skin_t, surface_emission,
                                    surface_los, surface_rmatrix,
                                    f_grid, rtp_pos, rtp_los,
                                    surface_rtprop_agenda );
      Index nsl = surface_los.nrows();
      //cout << "surf_los has " << surface_los.ncols() << " columns and "
      //     << nsl << " rows.\n";
      assert( surface_los.ncols()==1 || nsl==0 );

      // ARTS' surface_emission is equivalent to RT4's gnd_radiance (here:
      // surf_emis_vec) except for ARTS using Planck in terms of frequency,
      // while RT4 uses it in terms of wavelengths.
      // To derive gnd_radiance, we can either use ARTS surface_emission and
      // rescale it by B_lambda(surface_skin_t)/B_freq(surface_skin_t). Or we
      // can create it from surface_rmatrix and B_lambda(surface_skin_t). The
      // latter is more direct regarding use of B(T), but is requires
      // (re-)deriving diffuse reflectivity stokes components (which should be
      // the sum of the incident polar angle dependent reflectivities. so, that
      // might ensure better consistency. but then we should calculate
      // gnd_radiance only after surface_los to scat_za_grid conversion of the
      // reflection matrix). For now, we use the rescaling approach.
      for (Index f_index=0; f_index<nf; f_index++)
        {
          Numeric freq = f_grid[f_index];
          Numeric B_freq = planck(freq,surface_skin_t);
          Numeric B_lambda;
          Numeric wave = 1e6*SPEED_OF_LIGHT/freq;
          planck_function_( surface_skin_t, B_unit.c_str(), wave, B_lambda );
          Numeric B_ratio = B_lambda / B_freq;
          surf_emis_vec(f_index,rmu,joker) = surface_emission(f_index,joker);
          surf_emis_vec(f_index,rmu,joker) *= B_ratio;
        }
                                           

      // now we have to properly fill the RT4 stream directions in incident
      // angle dimension of the reflection matrix for RT4 with the agenda
      // output, which is given on/limited to surface_los. Only values inside
      // the (default) extrapolation range are filled; the rest is left as 0.
      // we do this for each incident stream separately as we need to check them
      // separately whether they are within allowed inter/extrapolation range
      // (ARTS can't do this for all elements of a vector at once. it will fail
      // for the whole vector if there's a single value outside.).

      // as we are rescaling surface_rmatrix within here, we need to keep its
      // original normalization for later renormalization. Create the container
      // here, fill in below.
      Vector R_arts(f_grid.nelem(),0.);

      if (nsl>1) // non-blackbody, non-specular reflection
        {
          for (Index f_index=0; f_index<nf; f_index++)
            R_arts[f_index] = surface_rmatrix(joker,f_index,0,0).sum();

          // Determine angle range weights in surface_rmatrix and de-scale
          // surface_rmatrix with those.
          Vector surf_int_grid(nsl+1);
          surf_int_grid[0] =
            surface_los(0,0)-0.5*(surface_los(1,0)-surface_los(0,0));
          surf_int_grid[nsl] =
            surface_los(nsl-1,0)+0.5*(surface_los(nsl-1,0)-surface_los(nsl-2,0));
          for (Index imu=1; imu<nsl; imu++)
            surf_int_grid[imu] = 0.5*(surface_los(imu-1,0)+surface_los(imu,0));
          surf_int_grid *= DEG2RAD;
          for (Index imu=0; imu<nsl; imu++)
            {
              //Numeric coslow = cos(2.*surf_int_grid[imu]);
              //Numeric coshigh = cos(2.*surf_int_grid[imu+1]);
              //Numeric w = 0.5*(coslow-coshigh);
              Numeric w = 0.5*(cos(2.*surf_int_grid[imu])-cos(2.*surf_int_grid[imu+1]));
              //cout << "at surf_los[" << imu << "]=" << surface_los(imu,0) << ":\n";
              //cout << "  angle weight derives as w = 0.5*(" << coslow << "-"
              //     << coshigh << ") = " << w << "\n";
              //cout << "  de-scaling with w from rmat="
              //     << surface_rmatrix(imu,0,0,0);
              surface_rmatrix(imu,joker,joker,joker) /= w;
              //cout << " to " << surface_rmatrix(imu,0,0,0) << "\n";
            }


          // Testing: interpolation in cos(za)
          //Vector mu_surf_los(nsl);
          //for (Index imu=0; imu<nsl; imu++)
          //  mu_surf_los[imu] = cos(surface_los(imu,0)*DEG2RAD);

          for (Index imu=0; imu<nummu; imu++)
          {
            //cout << "Doing incident dir #" <<imu << " at " << scat_za_grid[imu] << " degs\n";
            try
              {
                GridPos gp_za;
                gridpos( gp_za, surface_los(joker,0), scat_za_grid[imu] );
                // Testing: interpolation in cos(za)
                //gridpos( gp_za, mu_surf_los, mu_values[imu] );
                Vector itw(2);
                interpweights( itw, gp_za );


                // now apply the interpolation weights. since we're not in
                // python, we can't do the whole tensor at once, but need to
                // loop over the dimensions.
                for (Index f_index=0; f_index<nf; f_index++)
                  for (Index sto1=0; sto1<stokes_dim; sto1++)
                    for (Index sto2=0; sto2<stokes_dim; sto2++)
                      surf_refl_mat(f_index,imu,sto2,rmu,sto1) =
                        interp( itw, surface_rmatrix(joker,f_index,sto1,sto2), gp_za );
                // Apply new angle range weights - as this is for RT4, we apply
                // the actual RT4 angle (aka quadrature) weights:
                //cout << "  rescaling with quad weight w=" << quad_weights[imu]
                //     << " from " << surf_refl_mat(0,imu,0,rmu,0);
                surf_refl_mat(joker,imu,joker,rmu,joker) *=
                  (quad_weights[imu]*mu_values[imu]);
                //cout << " to " << surf_refl_mat(0,imu,0,rmu,0) << "\n";
              }
            catch( runtime_error e ) 
              { 
                // nothing to do here. we just leave the reflection matrix
                // entry at the 0.0 it was initialized with.
              }
          }
        }
      else if (nsl>0) // specular reflection
                      // no interpolation, no angle weight rescaling,
                      // just setting diagonal elements of surf_refl_mat.
        {
          for (Index f_index=0; f_index<nf; f_index++)
            R_arts[f_index] = surface_rmatrix(joker,f_index,0,0).sum();

          // surface_los angle should be identical to
          // 180. - (rtp_los=scat_za_grid[nummu+rmu]=180.-scat_za_grid[rmu])
          // = scat_za_grid[rmu].
          // check that, and if so, sort values into refmat(rmu,rmu).
          assert(is_same_within_epsilon(surface_los(0,0),scat_za_grid[rmu],1e-12));
          for (Index f_index=0; f_index<nf; f_index++)
            for (Index sto1=0; sto1<stokes_dim; sto1++)
              for (Index sto2=0; sto2<stokes_dim; sto2++)
                surf_refl_mat(f_index,rmu,sto2,rmu,sto1) =
                  surface_rmatrix(0,f_index,sto1,sto2);

        }
      //else {} // explicit blackbody
                // all surf_refl_mat elements to remain at 0., ie nothing to do.

      //eventually make sure the scaling of surf_refl_mat is correct
      for (Index f_index=0; f_index<nf; f_index++)
      {
        Numeric R_scale = 1.;
        Numeric R_rt4 = surf_refl_mat(f_index,joker,0,rmu,0).sum();
        if ( R_rt4==0. )
          {
            if ( R_arts[f_index]!=0. )
              {
                ostringstream os;
                os << "Something went wrong.\n"
                   << "At reflected stream #" << rmu
                   << ", power reflection coefficient for RT4\n"
                   << "became 0, although the one from surface_rtprop_agenda is "
                   << R_arts << ".\n";
                throw runtime_error(os.str());
              }
          }
        else
          {
            R_scale = R_arts[f_index]/R_rt4;
            surf_refl_mat(f_index,joker,joker,rmu,joker) *= R_scale;
          }
        //Numeric R_rert4 = surf_refl_mat(f_index,joker,0,rmu,0).sum();
        //cout << "at f#" << f_index << " R_arts=" << R_arts[f_index]
        //     << ", R_rt4=" << R_rt4 << ", R_scale=" << R_scale
        //     << ", rescaled R_rt4=" << R_rert4 << "\n";
      }
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
               const String& datapath,
               const Verbosity& verbosity )
{
    //emissivity.resize(4);
    //reflectivity.resize(4);

    Index nstokes=2;
    Index nummu=8;
    Index nuummu=0;
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
    ReadXML( height, "height", datapath+"z.xml", "", verbosity );
    ReadXML( temperatures, "temperatures", datapath+"T.xml", "", verbosity );
    ReadXML( gas_extinct, "gas_extinct", datapath+"abs_gas.xml", "", verbosity );
    ReadXML( abs_data, "abs_data", datapath+"abs_par.xml", "", verbosity );
    ReadXML( ext_data, "ext_data", datapath+"ext_par.xml", "", verbosity );
    ReadXML( sca_data, "sca_data", datapath+"sca_par.xml", "", verbosity );
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

    // dummy parameters necessary due to modified, flexible surface handling
    Tensor4 surf_refl_mat(nummu,nstokes,nummu,nstokes, 0.);
    Matrix surf_emis_vec(nummu,nstokes, 0.);
    Matrix ground_reflec(nstokes,nstokes, 0.);

    // Output variables
    Vector mu_values(nummu);
    Tensor3 up_rad(num_layers+1,nummu,nstokes, 0.);
    Tensor3 down_rad(num_layers+1,nummu,nstokes, 0.);

    radtrano_( nstokes,
               nummu,
               nuummu,
               max_delta_tau,
               quad_type.c_str(),
               ground_temp,
               ground_type.c_str(),
               ground_albedo,
               ground_index,
               ground_reflec.get_c_array(),
               surf_refl_mat.get_c_array(),
               surf_emis_vec.get_c_array(),
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

#endif /* ENABLE_RT4 */
