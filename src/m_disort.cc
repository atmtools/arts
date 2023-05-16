/*!
  \file   m_disort.cc
  \author Claudia Emde <claudia.emde@dlr.de>
          Manfred Brath <manfred.brath@uni-hamburg.de>
  \date   2006-02-06, 2022-01-22
  
  \brief  This file contains functions to use the multiple scattering 
  program (C)DISORT.
  
  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "array.h"
#include "auto_md.h"
#include "disort.h"
#include "m_general.h"
#include "math_funcs.h"
#include "messages.h"
 #include "geodetic.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void DisortCalc(Workspace& ws,
                    // WS Output:
                    Tensor7& cloudbox_field,
                    ArrayOfMatrix& disort_aux,
                    // WS Input
                    const Index& atmfields_checked,
                    const Index& atmgeom_checked,
                    const Index& scat_data_checked,
                    const Index& cloudbox_checked,
                    const Index& cloudbox_on,
                    const ArrayOfIndex& cloudbox_limits,
                    const Agenda& propmat_clearsky_agenda,
                    const Agenda& gas_scattering_agenda,
                    const Index& atmosphere_dim,
                    const Tensor4& pnd_field,
                    const Tensor3& t_field,
                    const Tensor3& z_field,
                    const Tensor4& vmr_field,
                    const Vector& p_grid,
                    const Vector& lat_true,
                    const Vector& lon_true,
                    const Vector& refellipsoid,
                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                    const ArrayOfSun& suns,
                    const Vector& f_grid,
                    const Vector& za_grid,
                    const Vector& aa_grid,
                    const Index& stokes_dim,
                    const Matrix& z_surface,
                    const Numeric& surface_skin_t,
                    const Vector& surface_scalar_reflectivity,
                    const Index& gas_scattering_do,
                    const Index& suns_do,
                    const ArrayOfString& disort_aux_vars,
                    const Index& nstreams,
                    const Index& Npfct,
                    const Index& only_tro,
                    const Index& cdisort_quiet,
                    const Index& emission,
                    const Index& intensity_correction,
                    const Verbosity& verbosity) {
  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on) {
    CREATE_OUT0;
    out0 << "  Cloudbox is off, DISORT calculation will be skipped.\n";
    return;
  }

  check_disort_input(cloudbox_on,
                     atmfields_checked,
                     atmgeom_checked,
                     cloudbox_checked,
                     scat_data_checked,
                     atmosphere_dim,
                     stokes_dim,
                     cloudbox_limits,
                     scat_data,
                     za_grid,
                     nstreams);

  //Check for number of suns
  ARTS_USER_ERROR_IF(suns.nelem() > 1,
                     "The simulation setup contains ",
                     suns.nelem(),
                     " suns. \n"
                     "Disort can handle only one sun.")

  //Check for aa_grid
  ARTS_USER_ERROR_IF(aa_grid.nelem() == 0,
                     "aa_grid has a size of 0.\n",
                     "aa_grid must have at least a size of one.")

  //allocate Varibale for direct (sun) source
  Vector sun_rte_los;
  Vector sun_pos(3);
  Vector cloudboxtop_pos(3);
  Index sun_on = suns_do;
  Numeric scale_factor;
  Index N_aa=aa_grid.nelem();

  if (sun_on){

    Vector lon_grid{lon_true[0] - 0.1, lon_true[0] + 0.1};
    Vector lat_grid{lat_true[0] - 0.1, lat_true[0] + 0.1};

    //Position of sun
    sun_pos = {suns[0].distance, suns[0].latitude, suns[0].longitude};

    // Position of top of cloudbox
    cloudboxtop_pos = {
        z_field(cloudbox_limits[1], 0, 0), lat_true[0], lon_true[0]};

    // calculate local position of sun at top of cloudbox
    rte_losGeometricFromRtePosToRtePos2(sun_rte_los,
                                        3,
                                        lat_grid,
                                        lon_grid,
                                        refellipsoid,
                                        cloudboxtop_pos,
                                        sun_pos,
                                        verbosity);

    //FIXME: IF we want to be correct and include refraction, we must calculate the
    // local position of sun via ppathFromRtePos2. The question is, is this needed,
    // because DISORT does not handle refraction at all.

    // Check if sun is above horizon, if not switch it off
    if (sun_rte_los[0] >= 90) {
      sun_on = 0;

      //set number of azimuth angle to 1
      N_aa = 1;

      CREATE_OUT0;
      out0 << "Sun is below the horizon\n";
      out0 << "Sun is ignored.\n";
      out0 << "cloudbox_field will have no azimuthal dependency!\n";
    }

    //get the cloudbox top distance to earth center.
    Numeric R_TOA = refell2r(refellipsoid,
                             lat_true[0]) +
                    cloudboxtop_pos[0];

    //get the distance between sun and cloudbox top
    Numeric R_Sun2CloudboxTop;
    distance3D(R_Sun2CloudboxTop,
               R_TOA,
               lat_true[0],
               lon_true[0],
               sun_pos[0],
               sun_pos[1],
               sun_pos[2]);

    // Geometric scaling factor, scales the sun spectral irradiance at the surface
    // of the sun to the spectral irradiance of the sun at cloubbox top.
    scale_factor=suns[0].radius*suns[0].radius/
                   (suns[0].radius*suns[0].radius+R_Sun2CloudboxTop*R_Sun2CloudboxTop);



  } else {
    CREATE_OUT3;
    out3 << "Disort calculation encountered aa_grid size larger than 1 in a case when it\n";
    out3 << "does not use aa_grid. Calculations are performed as if there is no aa_grid.\n";

    // set number of azimuth angle to 1
    N_aa = 1;
  }

  init_ifield(cloudbox_field,
                f_grid,
                cloudbox_limits,
                za_grid.nelem(),
                N_aa,
                stokes_dim);


  Vector albedo(f_grid.nelem(), 0.);
  Numeric btemp;

  get_disortsurf_props(
      albedo, btemp, f_grid, surface_skin_t, surface_scalar_reflectivity);

  run_cdisort(ws,
              cloudbox_field,
              disort_aux,
              f_grid,
              p_grid,
              z_field(joker, 0, 0),
              z_surface(0, 0),
              t_field(joker, 0, 0),
              vmr_field(joker, joker, 0, 0),
              pnd_field(joker, joker, 0, 0),
              scat_data,
              suns,
              propmat_clearsky_agenda,
              gas_scattering_agenda,
              cloudbox_limits,
              btemp,
              albedo,
              za_grid,
              aa_grid,
              sun_rte_los,
              gas_scattering_do,
              sun_on,
              disort_aux_vars,
              scale_factor,
              nstreams,
              Npfct,
              only_tro,
              cdisort_quiet,
              emission,
              intensity_correction,
              verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void DisortCalcWithARTSSurface(Workspace& ws,
                    // WS Output:
                    Tensor7& cloudbox_field,
                    ArrayOfMatrix& disort_aux,
                    // WS Input
                    const Index& atmfields_checked,
                    const Index& atmgeom_checked,
                    const Index& scat_data_checked,
                    const Index& cloudbox_checked,
                    const Index& cloudbox_on,
                    const ArrayOfIndex& cloudbox_limits,
                    const Agenda& propmat_clearsky_agenda,
                    const Agenda& surface_rtprop_agenda,
                    const Agenda& gas_scattering_agenda,
                    const Index& atmosphere_dim,
                    const Tensor4& pnd_field,
                    const Tensor3& t_field,
                    const Tensor3& z_field,
                    const Tensor4& vmr_field,
                    const Vector& p_grid,
                    const Vector& lat_true,
                    const Vector& lon_true,
                    const Vector& refellipsoid,
                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                    const ArrayOfSun& suns,
                    const Vector& f_grid,
                    const Vector& za_grid,
                    const Vector& aa_grid,
                    const Index& stokes_dim,
                    const Matrix& z_surface,
                    const Index& gas_scattering_do,
                    const Index& suns_do,
                    const ArrayOfString& disort_aux_vars,
                    const Index& nstreams,
                    const Index& Npfct,
                    const Index& only_tro,
                    const Index& cdisort_quiet,
                    const Index& emission,
                    const Index& intensity_correction,
                    const Numeric& inc_angle,
                    const Verbosity& verbosity) {

  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on) {
    CREATE_OUT0;
    out0 << "  Cloudbox is off, DISORT calculation will be skipped.\n";
    return;
  }

  check_disort_input(cloudbox_on,
                     atmfields_checked,
                     atmgeom_checked,
                     cloudbox_checked,
                     scat_data_checked,
                     atmosphere_dim,
                     stokes_dim,
                     cloudbox_limits,
                     scat_data,
                     za_grid,
                     nstreams);

  //Check for number of suns
  ARTS_USER_ERROR_IF(suns.nelem() > 1,
                     "The simulation setup contains ",
                     suns.nelem(),
                     " suns. \n"
                     "Disort can handle only one sun.")

  //Check for aa_grid
  ARTS_USER_ERROR_IF(aa_grid.nelem() == 0,
                     "aa_grid has a size of 0.\n",
                     "aa_grid must have at least a size of one.")

  //allocate Varibale for direct (sun) source
  Vector sun_rte_los;
  Vector sun_pos(3);
  Vector cloudboxtop_pos(3);
  Index sun_on = suns_do;
  Numeric scale_factor;
  Index N_aa=aa_grid.nelem();

  if (sun_on){

    Vector lon_grid{lon_true[0] - 0.1, lon_true[0] + 0.1};
    Vector lat_grid{lat_true[0] - 0.1, lat_true[0] + 0.1};

    //Position of sun
    sun_pos = {suns[0].distance, suns[0].latitude, suns[0].longitude};

    // Position of top of cloudbox
    cloudboxtop_pos = {
        z_field(cloudbox_limits[1], 0, 0), lat_true[0], lon_true[0]};

    // calculate local position of sun at top of cloudbox
    rte_losGeometricFromRtePosToRtePos2(sun_rte_los,
                                        3,
                                        lat_grid,
                                        lon_grid,
                                        refellipsoid,
                                        cloudboxtop_pos,
                                        sun_pos,
                                        verbosity);

    //FIXME: IF we want to be correct and include refraction, we must calculate the
    // local position of sun via ppathFromRtePos2. The question is, is this needed,
    // because DISORT does not handle refraction at all.

    // Check if sun is above horizon, if not switch it off
    if (sun_rte_los[0] >= 90) {
      sun_on = 0;

      // set number of azimuth angle to 1
      N_aa = 1;

      CREATE_OUT0;
      out0 << "Sun is below the horizon\n";
      out0 << "Sun is ignored.\n";
      out0 << "cloudbox_field will have no azimuthal dependency!\n";

    }

    //get the cloudbox top distance to earth center.
    Numeric R_TOA = refell2r(refellipsoid,
                             lat_true[0]) +
                    cloudboxtop_pos[0];

    //get the distance between sun and cloudbox top
    Numeric R_Sun2CloudboxTop;
    distance3D(R_Sun2CloudboxTop,
               R_TOA,
               lat_true[0],
               lon_true[0],
               sun_pos[0],
               sun_pos[1],
               sun_pos[2]);

    // Geometric scaling factor, scales the sun spectral irradiance at the surface
    // of the sun to the spectral irradiance of the sun at cloubbox top.
    scale_factor=suns[0].radius*suns[0].radius/
                   (suns[0].radius*suns[0].radius+R_Sun2CloudboxTop*R_Sun2CloudboxTop);



  } else {
    CREATE_OUT3;
    out3 << "Disort calculation encountered aa_grid size larger than 1 in a case when it\n";
    out3 << "does not use aa_grid. Calculations are performed as if there is no aa_grid.\n";

    // set number of azimuth angle to 1
    N_aa = 1;
  }

  init_ifield(cloudbox_field,
                f_grid,
                cloudbox_limits,
                za_grid.nelem(),
                N_aa,
                stokes_dim);


  Vector albedo(f_grid.nelem(), 0.);
  Numeric btemp;

  if (inc_angle<0 || inc_angle>90) {
    surf_albedoCalc(ws,
                    albedo,
                    btemp,
                    surface_rtprop_agenda,
                    f_grid,
                    za_grid,
                    z_surface(0,0),
                    verbosity);
  } else {
    surf_albedoCalcSingleAngle(ws,
                               albedo,
                               btemp,
                               surface_rtprop_agenda,
                               f_grid,
                               z_surface(0,0),
                               inc_angle);
  }

  run_cdisort(ws,
              cloudbox_field,
              disort_aux,
              f_grid,
              p_grid,
              z_field(joker, 0, 0),
              z_surface(0, 0),
              t_field(joker, 0, 0),
              vmr_field(joker, joker, 0, 0),
              pnd_field(joker, joker, 0, 0),
              scat_data,
              suns,
              propmat_clearsky_agenda,
              gas_scattering_agenda,
              cloudbox_limits,
              btemp,
              albedo,
              za_grid,
              aa_grid,
              sun_rte_los,
              gas_scattering_do,
              sun_on,
              disort_aux_vars,
              scale_factor,
              nstreams,
              Npfct,
              only_tro,
              cdisort_quiet,
              emission,
              intensity_correction,
              verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void DisortCalcClearsky(Workspace& ws,
                    // WS Output:
                    Tensor7& spectral_radiance_field,
                    ArrayOfMatrix& disort_aux,
                    // WS Input
                    const Index& atmfields_checked,
                    const Index& atmgeom_checked,
                    const Agenda& propmat_clearsky_agenda,
                    const Agenda& gas_scattering_agenda,
                    const Index& atmosphere_dim,
                    const Tensor3& t_field,
                    const Tensor3& z_field,
                    const Tensor4& vmr_field,
                    const Vector& p_grid,
                    const Vector& lat_true,
                    const Vector& lon_true,
                    const Vector& refellipsoid,
                    const ArrayOfSun& suns,
                    const Vector& f_grid,
                    const Vector& za_grid,
                    const Vector& aa_grid,
                    const Index& stokes_dim,
                    const Matrix& z_surface,
                    const Numeric& surface_skin_t,
                    const Vector& surface_scalar_reflectivity,
                    const Index& gas_scattering_do,
                    const Index& suns_do,
                    const ArrayOfString& disort_aux_vars,
                    const Index& nstreams,
                    const Index& cdisort_quiet,
                    const Index& emission,
                    const Index& intensity_correction,
                    const Verbosity& verbosity) {

  if (atmosphere_dim != 1)
    throw runtime_error(
        "For running DISORT, atmospheric dimensionality "
        "must be 1.\n");

  // Set cloudbox to cover complete atmosphere
  Index cloudbox_on;
  ArrayOfIndex cloudbox_limits;
  const Index cloudbox_checked = 1;
  //
  cloudboxSetFullAtm(cloudbox_on,
                     cloudbox_limits,
                     atmosphere_dim,
                     p_grid,
                     Vector(0),
                     Vector(0),
                     0.,
                     verbosity);

  // Create data matching no particles
  Tensor4 pnd_field;
  ArrayOfTensor4 dpnd_field_dx;
  ArrayOfArrayOfSingleScatteringData scat_data;
  const Index scat_data_checked = 1;
  //
  pnd_fieldZero(pnd_field,
                dpnd_field_dx,
                scat_data,
                atmosphere_dim,
                f_grid,
                cloudbox_limits,
                ArrayOfRetrievalQuantity(0),
                verbosity);

  Matrix optical_depth_dummy;

  DisortCalc(ws,
                 // WS Output:
                 spectral_radiance_field,
              disort_aux,
                 // WS Input
                 atmfields_checked,
                 atmgeom_checked,
                 scat_data_checked,
                 cloudbox_checked,
                 cloudbox_on,
                 cloudbox_limits,
                 propmat_clearsky_agenda,
                 gas_scattering_agenda,
                 atmosphere_dim,
                 pnd_field,
                 t_field,
                 z_field,
                 vmr_field,
                 p_grid,
                 lat_true,
                 lon_true,
                 refellipsoid,
                 scat_data,
                 suns,
                 f_grid,
                 za_grid,
                 aa_grid,
                 stokes_dim,
                 z_surface,
                 surface_skin_t,
                 surface_scalar_reflectivity,
                 gas_scattering_do,
                 suns_do,
                 disort_aux_vars,
                 nstreams,
                 181,
                 cdisort_quiet,
                 0,
                 emission,
                 intensity_correction,
                 verbosity);

}

/* Workspace method: Doxygen documentation will be auto-generated */
void DisortCalcIrradiance(Workspace& ws,
                // WS Output:
                Tensor5& spectral_irradiance_field,
                ArrayOfMatrix& disort_aux,
                // WS Input
                const Index& atmfields_checked,
                const Index& atmgeom_checked,
                const Index& scat_data_checked,
                const Agenda& propmat_clearsky_agenda,
                const Agenda& gas_scattering_agenda,
                const Index& atmosphere_dim,
                const Tensor4& pnd_field,
                const Tensor3& t_field,
                const Tensor3& z_field,
                const Tensor4& vmr_field,
                const Vector& p_grid,
                const Vector& lat_true,
                const Vector& lon_true,
                const Vector& refellipsoid,
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const ArrayOfSun& suns,
                const Vector& f_grid,
                const Index& stokes_dim,
                const Matrix& z_surface,
                const Numeric& surface_skin_t,
                const Vector& surface_scalar_reflectivity,
                const Index& gas_scattering_do,
                const Index& suns_do,
                const ArrayOfString& disort_aux_vars,
                const Index& nstreams,
                const Index& Npfct,
                const Index& only_tro,
                const Index& cdisort_quiet,
                const Index& emission,
                const Index& intensity_correction,
                const Verbosity& verbosity) {

  // Set cloudbox to cover complete atmosphere
  Index cloudbox_on;
  ArrayOfIndex cloudbox_limits;
  cloudboxSetFullAtm(cloudbox_on,
                     cloudbox_limits,
                     atmosphere_dim,
                     p_grid,
                     Vector(0),
                     Vector(0),
                     0.,
                     verbosity);

  const Index Nf = f_grid.nelem();
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;

  check_disort_irradiance_input(atmfields_checked,
                                atmgeom_checked,
                                scat_data_checked,
                                atmosphere_dim,
                                stokes_dim,
                                scat_data,
                                nstreams);

  //Check for number of suns
  ARTS_USER_ERROR_IF(suns.nelem() > 1,
                     "The simulation setup contains ",
                     suns.nelem(),
                     " suns. \n"
                     "Disort can handle only one sun.")

  //allocate Varibale for direct (sun) source
  Vector sun_rte_los;
  Vector sun_pos(3);
  Vector cloudboxtop_pos(3);
  Index sun_on = suns_do;
  Numeric scale_factor;

  spectral_irradiance_field.resize(Nf, Np_cloud, 1, 1, 2);
  spectral_irradiance_field = NAN;

  if (sun_on){

    Vector lon_grid{lon_true[0] - 0.1, lon_true[0] + 0.1};
    Vector lat_grid{lat_true[0] - 0.1, lat_true[0] + 0.1};

    //Position of sun
    sun_pos = {suns[0].distance, suns[0].latitude, suns[0].longitude};

    // Position of top of cloudbox
    cloudboxtop_pos = {
        z_field(cloudbox_limits[1], 0, 0), lat_true[0], lon_true[0]};

    // calculate local position of sun at top of cloudbox
    rte_losGeometricFromRtePosToRtePos2(sun_rte_los,
                                        3,
                                        lat_grid,
                                        lon_grid,
                                        refellipsoid,
                                        cloudboxtop_pos,
                                        sun_pos,
                                        verbosity);

    //FIXME: IF we want to be correct and include refraction, we must calculate the
    // local position of sun via ppathFromRtePos2. The question is, is this needed,
    // because DISORT does not handle refraction at all.

    // Check if sun is above horizon, if not switch it off
    if (sun_rte_los[0] >= 90) {
      sun_on = 0;

      CREATE_OUT0;
      out0 << "Sun is below the horizon\n";
      out0 << "Sun is ignored.\n";
    }

    //get the cloudbox top distance to earth center.
    Numeric R_TOA = refell2r(refellipsoid,
                             lat_true[0]) +
                    cloudboxtop_pos[0];

    //get the distance between sun and cloudbox top
    Numeric R_Sun2CloudboxTop;
    distance3D(R_Sun2CloudboxTop,
               R_TOA,
               lat_true[0],
               lon_true[0],
               sun_pos[0],
               sun_pos[1],
               sun_pos[2]);

    // Geometric scaling factor, scales the sun spectral irradiance at the surface
    // of the sun to the spectral irradiance of the sun at cloubbox top.
    scale_factor=suns[0].radius*suns[0].radius/
                   (suns[0].radius*suns[0].radius+R_Sun2CloudboxTop*R_Sun2CloudboxTop);

  }

  Vector albedo(f_grid.nelem(), 0.);
  Numeric btemp;

  get_disortsurf_props(
      albedo, btemp, f_grid, surface_skin_t, surface_scalar_reflectivity);

  run_cdisort_flux(ws,
                   spectral_irradiance_field,
                   disort_aux,
                   f_grid,
                   p_grid,
                   z_field(joker, 0, 0),
                   z_surface(0, 0),
                   t_field(joker, 0, 0),
                   vmr_field(joker, joker, 0, 0),
                   pnd_field(joker, joker, 0, 0),
                   scat_data,
                   suns,
                   propmat_clearsky_agenda,
                   gas_scattering_agenda,
                   cloudbox_limits,
                   btemp,
                   albedo,
                   sun_rte_los,
                   gas_scattering_do,
                   sun_on,
                   disort_aux_vars,
                   scale_factor,
                   nstreams,
                   Npfct,
                   only_tro,
                   cdisort_quiet,
                   emission,
                   intensity_correction,
                   verbosity);
}