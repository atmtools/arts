/*!
  \file   doit.cc
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   Wed Jun 04 11:03:57 2003
  
  \brief  This file contains functions to calculate the radiative transfer
  inside the cloudbox using the DOIT method.
  
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "doit.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <workspace.h>
#include "array.h"
#include "arts_constants.h"
#include "arts_conversions.h"
#include "atm.h"
#include "check_input.h"
#include "cloudbox.h"
#include "debug.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpack_data.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "rtepack.h"
#include "rtepack_multitype.h"
#include "rtepack_rtestep.h"
#include "rtepack_stokes_vector.h"
#include "sorting.h"
#include "special_interp.h"
#include "species_tags.h"
#include "xml_io.h"

inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);

//FIXME function name of 'rte_step_doit_replacement' should be replaced by
//proper name
void rte_step_doit_replacement(  //Output and Input:
    Stokvec& stokes_vec,
    Muelmat& trans_mat,
    //Input
    const Propmat& ext_mat_av,
    const Stokvec& abs_vec_av,
    const Stokvec& sca_vec_av,
    const Numeric& lstep,
    const Numeric& rtp_planck_value,
    const bool& trans_is_precalc) {
  ARTS_ASSERT(rtp_planck_value >= 0);
  ARTS_ASSERT(lstep >= 0);

  if (!trans_is_precalc) {
    trans_mat = rtepack::exp(-lstep * ext_mat_av);
  }

  stokes_vec = rtepack::linear_step(
      trans_mat, stokes_vec,
      inv(ext_mat_av) * (abs_vec_av * rtp_planck_value + sca_vec_av));
}

void cloud_fieldsCalc(const Workspace& ws,
                      // Output and Input:
                      Tensor5View ext_mat_field,
                      Tensor4View abs_vec_field,
                      // Input:
                      const Agenda& spt_calc_agenda,
                      const Index& za_index,
                      const Index& aa_index,
                      const ArrayOfIndex& cloudbox_limits,
                      ConstTensor3View t_field,
                      ConstTensor4View pnd_field) {
  // Input variables are checked in the WSMs i_fieldUpdateSeqXXX, from
  // where this function is called.

  const Index N_se = pnd_field.nbooks();

  ARTS_ASSERT(ext_mat_field.ncols() == ext_mat_field.nrows() &&
         ext_mat_field.ncols() == abs_vec_field.ncols());

  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;

  // If atmosohere_dim == 1
  Index Nlat_cloud = 1;
  Index Nlon_cloud = 1;

  Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
  Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;

  // Initialize ext_mat(_spt), abs_vec(_spt)
  // Resize and initialize variables for storing optical properties
  // of all scattering elements.
  ArrayOfStokvecVector abs_vec_spt_local(N_se, StokvecVector{1, Stokvec{0, 0, 0, 0}});
  ArrayOfPropmatVector ext_mat_spt_local(N_se, PropmatVector{1, Propmat{0, 0, 0, 0, 0, 0, 0}});

  StokvecVector abs_vec_local;
  PropmatVector ext_mat_local;
  Numeric rtp_temperature_local;

  // Calculate ext_mat, abs_vec for all points inside the cloudbox.
  // sca_vec can be obtained from the workspace variable doit_scat_field.
  // As we need the average for each layer, it makes sense to calculte
  // the coefficients once and store them in an array instead of
  // calculating at each point the coefficient of the point above and
  // the point below.
  // To use special interpolation functions for atmospheric fields we
  // use ext_mat_field and abs_vec_field:

  // Loop over all positions inside the cloudbox defined by the
  // cloudbox_limits.
  for (Index scat_p_index_local = 0; scat_p_index_local < Np_cloud;
       scat_p_index_local++) {
    for (Index scat_lat_index_local = 0; scat_lat_index_local < Nlat_cloud;
         scat_lat_index_local++) {
      for (Index scat_lon_index_local = 0; scat_lon_index_local < Nlon_cloud;
           scat_lon_index_local++) {
        rtp_temperature_local =
            t_field(scat_p_index_local + cloudbox_limits[0],
                    scat_lat_index_local + cloudbox_limits[2],
                    scat_lon_index_local + cloudbox_limits[4]);

        //Calculate optical properties for individual scattering elements:
        //( Execute agendas silently. )
        spt_calc_agendaExecute(ws,
                               ext_mat_spt_local,
                               abs_vec_spt_local,
                               scat_p_index_local,
                               scat_lat_index_local,
                               scat_lon_index_local,
                               rtp_temperature_local,
                               za_index,
                               aa_index,
                               spt_calc_agenda);
        /*
// so far missing here (accessed through workspace within agenda):
// - scat_data
// - za_grid, aa_grid
// - f_index
              opt_prop_sptFromScat_data(ext_mat_spt_local, abs_vec_spt_local,
                                        scat_data, 1,
                                        za_grid, aa_grid,
                                        za_index, aa_index,
                                        f_index,
                                        rtp_temperature_local,
                                        pnd_field, 
                                        scat_p_index_local,
                                        scat_lat_index_local,
                                        scat_lon_index_local,
                                        );
*/

        opt_prop_bulkCalc(ext_mat_local,
                          abs_vec_local,
                          ext_mat_spt_local,
                          abs_vec_spt_local,
                          Tensor4{pnd_field},
                          scat_p_index_local,
                          scat_lat_index_local,
                          scat_lon_index_local);

        // Store coefficients in arrays for the whole cloudbox.
        auto mmat = rtepack::to_muelmat(ext_mat_local[0]);
        for (Index i = 0; i < 4; i++) {
          abs_vec_field(scat_p_index_local, scat_lat_index_local,
                        scat_lon_index_local, i) = abs_vec_local[0].data[i];
          for (Index j = 0; j < 4; j++) {
            ext_mat_field(scat_p_index_local, scat_lat_index_local,
                          scat_lon_index_local, i, j) = mmat(i, j);
          }
        }
      }
    }
  }
}

void cloud_RT_surface(const Workspace& ws,
                      //Output
                      Tensor6View cloudbox_field_mono,
                      //Input
                      const Agenda& surface_rtprop_agenda,
                      ConstVectorView f_grid,
                      const Index& f_index,
                      const Ppath& ppath_step,
                      const ArrayOfIndex& cloudbox_limits,
                      ConstVectorView za_grid,
                      const Index& za_index) {

  Matrix iy;

  // Local output of surface_rtprop_agenda.
  SurfacePoint surface_point;
  Matrix surface_emission;
  Matrix surface_los;
  Tensor4 surface_rmatrix;

  //Set rte_pos and rte_los to match the last point in ppath.

  Index np = ppath_step.np;

  Vector rte_pos;  // ppath_step.pos contains two columns for 1D
  rte_pos.resize(ppath_step.dim);
  rte_pos = ppath_step.pos(np - 1, Range(0, ppath_step.dim));

  Vector rte_los;
  rte_los.resize(ppath_step.los.ncols());
  rte_los = ppath_step.los(np - 1, joker);

  //Execute the surface_rtprop_agenda which gives the surface
  //parameters.
  surface_rtprop_agendaExecute(ws,
                               surface_point,
                               surface_emission,
                               surface_los,
                               surface_rmatrix,
                               Vector(1, f_grid[f_index]),
                               rte_pos,
                               rte_los,
                               surface_rtprop_agenda);

  iy = surface_emission;

  Index nlos = surface_los.nrows();

  if (nlos > 0) {
    Vector rtmp(4);  // Reflected Stokes vector for 1 frequency

    for (Index ilos = 0; ilos < nlos; ilos++) {
      // Several things needs to be fixed here. As far as I understand it,
      // this works only for specular cases and if the lower cloudbox limit
      // is exactly at the surface (PE, 120828)

      mult(rtmp,
           surface_rmatrix(ilos, 0, joker, joker),
           cloudbox_field_mono(cloudbox_limits[0],
                             0,
                             0,
                             (za_grid.nelem() - 1 - za_index),
                             0,
                             joker));
      iy(0, joker) += rtmp;
    }
  }
  cloudbox_field_mono(cloudbox_limits[0], 0, 0, za_index, 0, joker) =
      iy(0, joker);
}

void cloudbox_field_ngAcceleration(Tensor6& cloudbox_field_mono,
                                 const ArrayOfTensor6& acceleration_input,
                                 const Index& accelerated) {
  const Index N_p = cloudbox_field_mono.nvitrines();
  const Index N_za = cloudbox_field_mono.npages();

  // Loop over 4 components of Stokes Vector
  for (Index i = 0; i < accelerated; ++i) {
    ConstMatrixView S1 = acceleration_input[0](joker, 0, 0, joker, 0, i);
    ConstMatrixView S2 = acceleration_input[1](joker, 0, 0, joker, 0, i);
    ConstMatrixView S3 = acceleration_input[2](joker, 0, 0, joker, 0, i);
    ConstMatrixView S4 = acceleration_input[3](joker, 0, 0, joker, 0, i);

    ConstMatrixView J = S4;
    Matrix Q1;
    Matrix Q2;
    Matrix Q3;
    Numeric A1 = 0;
    Numeric A2B1 = 0;
    Numeric B2 = 0;
    Numeric C1 = 0;
    Numeric C2 = 0;
    Numeric NGA = 0;
    Numeric NGB = 0;

    // Q1 = -2*S3 + S4 + S2

    Q1 = S3;
    Q1 *= -2;
    Q1 += S4;
    Q1 += S2;

    // Q2 = S4 - S3 - S2 + S1
    Q2 = S4;
    Q2 -= S3;
    Q2 -= S2;
    Q2 += S1;

    //Q3 = S4 - S3
    Q3 = S4;
    Q3 -= S3;

    for (Index p_index = 0; p_index < N_p; ++p_index) {
      for (Index za_index = 0; za_index < N_za; ++za_index) {
        A1 += Q1(p_index, za_index) * Q1(p_index, za_index) *
              J(p_index, za_index);
        A2B1 += Q2(p_index, za_index) * Q1(p_index, za_index) *
                J(p_index, za_index);
        B2 += Q2(p_index, za_index) * Q2(p_index, za_index) *
              J(p_index, za_index);
        C1 += Q1(p_index, za_index) * Q3(p_index, za_index) *
              J(p_index, za_index);
        C2 += Q2(p_index, za_index) * Q3(p_index, za_index) *
              J(p_index, za_index);
      }
    }

    NGA = (C1 * B2 - C2 * A2B1) / (A1 * B2 - A2B1 * A2B1);
    NGB = (C2 * A1 - C1 * A2B1) / (A1 * B2 - A2B1 * A2B1);

    if (!std::isnan(NGB) && !std::isnan(NGA)) {
      // Calculating the accelerated field
      for (Index p_index = 0; p_index < N_p; ++p_index) {
        for (Index za_index = 0; za_index < N_za; ++za_index) {
          Q1(p_index, za_index) = (1 - NGA - NGB) * S4(p_index, za_index) +
                                  NGA * S3(p_index, za_index) +
                                  NGB * S2(p_index, za_index);
        }
      }
      cloudbox_field_mono(joker, 0, 0, joker, 0, i) = Q1;
    }
  }
}

void interp_cloud_coeff1D(  //Output
    Tensor3View ext_mat_int,
    MatrixView abs_vec_int,
    MatrixView sca_vec_int,
    MatrixView cloudbox_field_mono_int,
    VectorView t_int,
    MatrixView vmr_list_int,
    VectorView p_int,
    //Input
    ConstTensor5View ext_mat_field,
    ConstTensor4View abs_vec_field,
    ConstTensor6View doit_scat_field,
    ConstTensor6View cloudbox_field_mono,
    ConstTensor3View t_field,
    ConstTensor4View vmr_field,
    ConstVectorView p_grid,
    const Ppath& ppath_step,
    const ArrayOfIndex& cloudbox_limits,
    ConstVectorView za_grid,
    const Index& scat_za_interp) {
  // Gridpositions inside the cloudbox.
  // The optical properties are stored only inside the
  // cloudbox. For interpolation we use grids
  // inside the cloudbox.
  ArrayOfGridPos cloud_gp_p = ppath_step.gp_p;

  for (Index i = 0; i < ppath_step.np; i++)
    cloud_gp_p[i].idx -= cloudbox_limits[0];

  // Grid index for points at upper limit of cloudbox must be shifted
  const Index n1 = cloudbox_limits[1] - cloudbox_limits[0];
  gridpos_upperend_check(cloud_gp_p[0], n1);
  gridpos_upperend_check(cloud_gp_p[ppath_step.np - 1], n1);

  Matrix itw(cloud_gp_p.size(), 2);
  interpweights(itw, cloud_gp_p);

  // The zenith angles of the propagation path are needed as we have to
  // interpolate the intensity field and the scattered field on the
  // right angles.
  Vector los_grid{ppath_step.los(joker, 0)};

  ArrayOfGridPos gp_za(los_grid.nelem());
  gridpos(gp_za, za_grid, los_grid);

  Matrix itw_p_za(cloud_gp_p.size(), 4);
  interpweights(itw_p_za, cloud_gp_p, gp_za);

  // Calculate the average of the coefficients for the layers
  // to be considered in the
  // radiative transfer calculation.

  for (Index i = 0; i < 4; i++) {
    // Extinction matrix requires a second loop
    // over 4
    for (Index j = 0; j < 4; j++) {
      //
      // Interpolation of ext_mat
      //
      interp(ext_mat_int(i, j, joker),
             itw,
             ext_mat_field(joker, 0, 0, i, j),
             cloud_gp_p);
    }
    // Particle absorption vector:
    //
    // Interpolation of abs_vec
    //  //
    interp(
        abs_vec_int(i, joker), itw, abs_vec_field(joker, 0, 0, i), cloud_gp_p);
    //
    // Scattered field:
    //
    //

    if (scat_za_interp == 0)  // linear interpolation
    {
      interp(sca_vec_int(i, joker),
             itw_p_za,
             doit_scat_field(joker, 0, 0, joker, 0, i),
             cloud_gp_p,
             gp_za);
      interp(cloudbox_field_mono_int(i, joker),
             itw_p_za,
             cloudbox_field_mono(joker, 0, 0, joker, 0, i),
             cloud_gp_p,
             gp_za);
    } else if (scat_za_interp == 1)  //polynomial interpolation
    {
      // These intermediate variables are needed because polynomial
      // interpolation is not implemented as multidimensional
      // interpolation.
      Tensor3 sca_vec_int_za(
          4, ppath_step.np, za_grid.nelem(), 0.);
      Tensor3 cloudbox_field_mono_int_za(
          4, ppath_step.np, za_grid.nelem(), 0.);
      for (Index za = 0; za < za_grid.nelem(); za++) {
        interp(sca_vec_int_za(i, joker, za),
               itw,
               doit_scat_field(joker, 0, 0, za, 0, i),
               cloud_gp_p);
        interp(cloudbox_field_mono_int_za(i, joker, za),
               itw,
               cloudbox_field_mono(joker, 0, 0, za, 0, i),
               cloud_gp_p);
      }
      for (Index ip = 0; ip < ppath_step.np; ip++) {
        sca_vec_int(i, ip) = interp_poly(za_grid,
                                         sca_vec_int_za(i, ip, joker),
                                         los_grid[ip],
                                         gp_za[ip]);
        cloudbox_field_mono_int(i, ip) =
            interp_poly(za_grid,
                        cloudbox_field_mono_int_za(i, ip, joker),
                        los_grid[ip],
                        gp_za[ip]);
      }
    }
  }
  //
  // Planck function
  //
  // Interpolate temperature field
  //
  interp(t_int, itw, t_field(joker, 0, 0), ppath_step.gp_p);
  //
  // The vmr_field is needed for the gaseous absorption
  // calculation.
  //
  const Index N_species = vmr_field.nbooks();
  //
  // Interpolated vmr_list, holds a vmr_list for each point in
  // ppath_step.
  //
  Vector vmr_int(ppath_step.np);

  for (Index i_sp = 0; i_sp < N_species; i_sp++) {
    interp(vmr_int, itw, vmr_field(i_sp, joker, 0, 0), ppath_step.gp_p);
    vmr_list_int(i_sp, joker) = vmr_int;
  }
  //
  // Interpolate pressure
  //
  itw2p(p_int, p_grid, ppath_step.gp_p, itw);
}

void za_gridOpt(  //Output:
    Vector& za_grid_opt,
    Matrix& cloudbox_field_opt,
    // Input
    ConstVectorView za_grid_fine,
    ConstTensor6View cloudbox_field_mono,
    const Numeric& acc,
    const Index& scat_za_interp) {
  Index N_za = za_grid_fine.nelem();

  ARTS_ASSERT(cloudbox_field_mono.npages() == N_za);

  Index N_p = cloudbox_field_mono.nvitrines();

  Vector i_approx_interp(N_za);
  Vector za_reduced(2);

  ArrayOfIndex idx;
  idx.push_back(0);
  idx.push_back(N_za - 1);
  ArrayOfIndex idx_unsorted;

  Numeric max_diff = 100;

  ArrayOfGridPos gp_za(N_za);
  Matrix itw(za_grid_fine.nelem(), 2);

  ArrayOfIndex i_sort;
  Vector diff_vec(N_za);
  Vector max_diff_za(N_p);
  ArrayOfIndex ind_za(N_p);
  Numeric max_diff_p;
  Index ind_p = 0;

  while (max_diff > acc) {
    za_reduced.resize(idx.size());
    cloudbox_field_opt.resize(N_p, idx.size());
    max_diff_za = 0.;
    max_diff_p = 0.;

    // Interpolate reduced intensity field on fine za_grid for
    // all pressure levels
    for (Index i_p = 0; i_p < N_p; i_p++) {
      for (Size i_za_red = 0; i_za_red < idx.size(); i_za_red++) {
        za_reduced[i_za_red] = za_grid_fine[idx[i_za_red]];
        cloudbox_field_opt(i_p, i_za_red) =
            cloudbox_field_mono(i_p, 0, 0, idx[i_za_red], 0, 0);
      }
      // Calculate grid positions
      gridpos(gp_za, za_reduced, za_grid_fine);
      //linear interpolation
      if (scat_za_interp == 0 || idx.size() < 3) {
        interpweights(itw, gp_za);
        interp(i_approx_interp, itw, cloudbox_field_opt(i_p, joker), gp_za);
      } else if (scat_za_interp == 1) {
        for (Index i_za = 0; i_za < N_za; i_za++) {
          i_approx_interp[i_za] = interp_poly(za_reduced,
                                              cloudbox_field_opt(i_p, joker),
                                              za_grid_fine[i_za],
                                              gp_za[i_za]);
        }
      } else
        // Interpolation method not defined
        ARTS_ASSERT(false);

      // Calculate differences between approximated i-vector and
      // exact i_vector for the i_p pressure level
      for (Index i_za = 0; i_za < N_za; i_za++) {
        diff_vec[i_za] = abs(cloudbox_field_mono(i_p, 0, 0, i_za, 0, 0) -
                             i_approx_interp[i_za]);
        if (diff_vec[i_za] > max_diff_za[i_p]) {
          max_diff_za[i_p] = diff_vec[i_za];
          ind_za[i_p] = i_za;
        }
      }
      // Take maximum value of max_diff_za
      if (max_diff_za[i_p] > max_diff_p) {
        max_diff_p = max_diff_za[i_p];
        ind_p = i_p;
      }
    }

    //Transform in %
    max_diff =
        max_diff_p / cloudbox_field_mono(ind_p, 0, 0, ind_za[ind_p], 0, 0) * 100.;

    idx.push_back(ind_za[ind_p]);
    idx_unsorted = idx;

    i_sort.resize(idx_unsorted.size());
    get_sorted_indexes(i_sort, idx_unsorted);

    for (Size i = 0; i < idx_unsorted.size(); i++)
      idx[i] = idx_unsorted[i_sort[i]];

    za_reduced.resize(idx.size());
  }

  za_grid_opt.resize(idx.size());
  cloudbox_field_opt.resize(N_p, idx.size());
  for (Size i = 0; i < idx.size(); i++) {
    za_grid_opt[i] = za_grid_fine[idx[i]];
    cloudbox_field_opt(joker, i) = cloudbox_field_mono(joker, 0, 0, idx[i], 0, 0);
  }
}

void doit_scat_fieldNormalize(const Workspace& ws,
                              Tensor6& doit_scat_field,
                              const Tensor6& cloudbox_field_mono,
                              const ArrayOfIndex& cloudbox_limits,
                              const Agenda& spt_calc_agenda,
                              const Vector& za_grid,
                              const Vector& aa_grid,
                              const Tensor4& pnd_field,
                              const Tensor3& t_field,
                              const Numeric& norm_error_threshold,
                              const Index& norm_debug) {
  ARTS_USER_ERROR ("Only 1D is supported here for now");

  // Number of zenith angles.
  const Index Nza = za_grid.nelem();

  ARTS_USER_ERROR_IF (za_grid[0] != 0. || za_grid[Nza - 1] != 180.,
                      "The range of *za_grid* must [0 180].");

  // Number of azimuth angles.
  const Index Naa = aa_grid.nelem();

  ARTS_USER_ERROR_IF (Naa > 1 && (aa_grid[0] != 0. || aa_grid[Naa - 1] != 360.),
                      "The range of *aa_grid* must [0 360].");

  // To use special interpolation functions for atmospheric fields we
  // use ext_mat_field and abs_vec_field:
  Tensor5 ext_mat_field(cloudbox_limits[1] - cloudbox_limits[0] + 1,
                        1,
                        1,
                        4,
                        4,
                        0.);
  Tensor4 abs_vec_field(
      cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1, 4, 0.);

  const Index Np = doit_scat_field.nvitrines();

  Tensor5 doit_scat_ext_field(doit_scat_field.nvitrines(),
                              doit_scat_field.nshelves(),
                              doit_scat_field.nbooks(),
                              doit_scat_field.npages(),
                              doit_scat_field.nrows(),
                              0.);

  Index aa_index_local = 0;

  // Calculate scattering extinction field
  for (Index za_index_local = 0; za_index_local < Nza;
       za_index_local++) {
    // This function has to be called inside the angular loop, as
    // spt_calc_agenda takes *za_index* and *aa_index*
    // from the workspace.
    cloud_fieldsCalc(ws,
                     ext_mat_field,
                     abs_vec_field,
                     spt_calc_agenda,
                     za_index_local,
                     aa_index_local,
                     cloudbox_limits,
                     t_field,
                     pnd_field);

    for (Index p_index = 0;
         p_index <= (cloudbox_limits[1] - cloudbox_limits[0]);
         p_index++) {
      // For all in p_grid (in cloudbox):
      // I_ext = (ext_mat_field - abs_vec_field) * cloudbox_field_mono
      // equivalent to:
      // I_ext = I * (K11-a1) + Q * (K12 - a2) + U * (K13 - a3) + V * (K14 - a4)
      for (Index i = 0; i < 4; i++) {
        doit_scat_ext_field(p_index, 0, 0, za_index_local, 0) +=
            cloudbox_field_mono(p_index, 0, 0, za_index_local, 0, i) *
            (ext_mat_field(p_index, 0, 0, 0, i) -
             abs_vec_field(p_index, 0, 0, i));
      }
    }
  }

  Numeric corr_max = .0;
  Index corr_max_p_index = -1;

  for (Index p_index = 0; p_index < Np; p_index++) {
    // Calculate scattering integrals
    const Numeric scat_int = AngIntegrate_trapezoid(
        doit_scat_field(p_index, 0, 0, joker, 0, 0), za_grid);

    const Numeric scat_ext_int = AngIntegrate_trapezoid(
        doit_scat_ext_field(p_index, 0, 0, joker, 0), za_grid);

    // Calculate factor between scattered extinction field integral
    // and scattered field integral
    const Numeric corr_factor = scat_ext_int / scat_int;

    // If no scattering is present, the correction factor can become
    // inf or nan. We just don't apply it for those cases.
    if (!std::isnan(corr_factor) && !std::isinf(corr_factor)) {
      if (abs(corr_factor) > abs(corr_max)) {
        corr_max = corr_factor;
        corr_max_p_index = p_index;
      }
      if (norm_debug) {
      }
      ARTS_USER_ERROR_IF (abs(1. - corr_factor) > norm_error_threshold,
          "ERROR: DOIT correction factor exceeds threshold (=",
          norm_error_threshold, "): ", std::setprecision(4),
          1. - corr_factor, " at p_index ", p_index, "\n")
      if (abs(1. - corr_factor) > norm_error_threshold / 2.) {
      }

      // Scale scattered field with correction factor
      doit_scat_field(p_index, 0, 0, joker, 0, joker) *= corr_factor;
    } else if (norm_debug) {
    }
  }

  std::ostringstream os;
  if (corr_max_p_index != -1) {
    os << "  Max. DOIT correction factor in this iteration: " << 1. - corr_max
       << " at p_index " << corr_max_p_index << "\n";
  } else {
    os << "  No DOIT correction performed in this iteration.\n";
  }
}
