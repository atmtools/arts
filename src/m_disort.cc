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
#include "atm.h"
#include <workspace.h>
#include "debug.h"
#include "disort.h"
#include "m_general.h"
#include "math_funcs.h"
#include "matpack_data.h"
 #include "geodetic.h"
#include "species_tags.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void cloudbox_fieldDisort(const Workspace& ws,
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
                    const Tensor4& pnd_field,
                    const AtmField& atm_field,
                    const SurfaceField& surface_field,
                    const Vector& lat_true,
                    const Vector& lon_true,
                    const ArrayOfArrayOfSpeciesTag& abs_species,
                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                    const ArrayOfSun& suns,
                    const Vector& f_grid,
                    const Vector& za_grid,
                    const Vector& aa_grid,
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
                    const Index& intensity_correction) {
  // FIXME: REQUIRES REGULAR GRIDS
  Vector z_grid, lat_grid, lon_grid;
  Tensor3 t_field, p_field, wind_u_field;
  Tensor4 vmr_field;
  ARTS_USER_ERROR("ERROR")
  //const auto& z_grid = atm_field.grid[0];
//const auto& p_field = atm_field[Atm::Key::p].get<const Tensor3&>();
//const auto& t_field = atm_field[Atm::Key::t].get<const Tensor3&>();
//const auto vmr_field = Atm::extract_specs_content(atm_field, abs_species);

  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on) {
    return;
  }

  check_disort_input(cloudbox_on,
                     atmfields_checked,
                     atmgeom_checked,
                     cloudbox_checked,
                     scat_data_checked,
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
    cloudboxtop_pos = {z_grid[cloudbox_limits[1]], lat_true[0], lon_true[0]};

    // calculate local position of sun at top of cloudbox
    ARTS_ASSERT(false)
    /*
    rte_losGeometricFromRtePosToRtePos2(sun_rte_los,
                                        lat_grid,
                                        lon_grid,
                                        refellipsoid,
                                        cloudboxtop_pos,
                                        sun_pos,
                                        );
*/

    //FIXME: IF we want to be correct and include refraction, we must calculate the
    // local position of sun via ppathFromRtePos2. The question is, is this needed,
    // because DISORT does not handle refraction at all.

    // Check if sun is above horizon, if not switch it off
    if (sun_rte_los[0] >= 90) {
      sun_on = 0;

      //set number of azimuth angle to 1
      N_aa = 1;
    }

    //get the cloudbox top distance to earth center.
    Numeric R_TOA;// = refell2r(refellipsoid, lat_true[0]) + cloudboxtop_pos[0];
ARTS_USER_ERROR("ERROR")

    //get the distance between sun and cloudbox top
    Numeric R_Sun2CloudboxTop;
    ARTS_USER_ERROR("ERROR")
  /*  distance3D(R_Sun2CloudboxTop,
               R_TOA,
               lat_true[0],
               lon_true[0],
               sun_pos[0],
               sun_pos[1],
               sun_pos[2]);*/

    // Geometric scaling factor, scales the sun spectral irradiance at the surface
    // of the sun to the spectral irradiance of the sun at cloubbox top.
    scale_factor=suns[0].radius*suns[0].radius/
                   (suns[0].radius*suns[0].radius+R_Sun2CloudboxTop*R_Sun2CloudboxTop);



  } else {
    // set number of azimuth angle to 1
    N_aa = 1;
  }

  init_ifield(cloudbox_field,
                f_grid,
                cloudbox_limits,
                za_grid.nelem(),
                N_aa);


  Vector albedo(f_grid.nelem(), 0.);
  Numeric btemp;

  get_disortsurf_props(
      albedo, btemp, f_grid, surface_skin_t, surface_scalar_reflectivity);

  run_cdisort(ws,
              cloudbox_field,
              disort_aux,
              f_grid,
              p_field(joker, 0, 0),
              z_grid,
              surface_field.single_value(Surf::Key::h, 0, 0),  // FIXME: MUst know lat and lon
              t_field(joker, 0, 0),
              vmr_field(joker, joker, 0, 0),
              pnd_field(joker, joker, 0, 0),
              abs_species,
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
              intensity_correction);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void cloudbox_fieldDisortWithARTSSurface(const Workspace& ws,
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
                    const Tensor4& pnd_field,
                    const AtmField& atm_field,
                    const SurfaceField& surface_field,
                    const Vector& lat_true,
                    const Vector& lon_true,
                    const ArrayOfArrayOfSpeciesTag& abs_species,
                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                    const ArrayOfSun& suns,
                    const Vector& f_grid,
                    const Vector& za_grid,
                    const Vector& aa_grid,
                    const Index& gas_scattering_do,
                    const Index& suns_do,
                    const ArrayOfString& disort_aux_vars,
                    const Index& nstreams,
                    const Index& Npfct,
                    const Index& only_tro,
                    const Index& cdisort_quiet,
                    const Index& emission,
                    const Index& intensity_correction,
                    const Numeric& inc_angle) {
  // FIXME: REQUIRES REGULAR GRIDS
  Vector z_grid, lat_grid, lon_grid;
  Tensor3 t_field, p_field, wind_u_field;
  Tensor4 vmr_field;
  ARTS_USER_ERROR("ERROR")
  //const auto& z_grid = atm_field.grid[0];
  //const auto& p_field = atm_field[Atm::Key::p].get<const Tensor3&>();
  //const auto& t_field = atm_field[Atm::Key::t].get<const Tensor3&>();
  //const auto vmr_field = Atm::extract_specs_content(atm_field, abs_species);

  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on) {
    return;
  }

  check_disort_input(cloudbox_on,
                     atmfields_checked,
                     atmgeom_checked,
                     cloudbox_checked,
                     scat_data_checked,
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
        z_grid[cloudbox_limits[1]], lat_true[0], lon_true[0]};

    // calculate local position of sun at top of cloudbox
    ARTS_ASSERT(false)
    /*
    rte_losGeometricFromRtePosToRtePos2(sun_rte_los,
                                        lat_grid,
                                        lon_grid,
                                        refellipsoid,
                                        cloudboxtop_pos,
                                        sun_pos,
                                        );
*/
    //FIXME: IF we want to be correct and include refraction, we must calculate the
    // local position of sun via ppathFromRtePos2. The question is, is this needed,
    // because DISORT does not handle refraction at all.

    // Check if sun is above horizon, if not switch it off
    if (sun_rte_los[0] >= 90) {
      sun_on = 0;

      // set number of azimuth angle to 1
      N_aa = 1;

    }

    //get the cloudbox top distance to earth center.
    Numeric R_TOA;// = refell2r(refellipsoid,
                   //          lat_true[0]) +
                   // cloudboxtop_pos[0];
ARTS_USER_ERROR("ERROR")

    //get the distance between sun and cloudbox top
    Numeric R_Sun2CloudboxTop;
//    distance3D(R_Sun2CloudboxTop,
  //             R_TOA,
    //           lat_true[0],
      //         lon_true[0],
        //       sun_pos[0],
          //     sun_pos[1],
            //   sun_pos[2]);
ARTS_USER_ERROR("ERROR")

    // Geometric scaling factor, scales the sun spectral irradiance at the surface
    // of the sun to the spectral irradiance of the sun at cloubbox top.
    scale_factor=suns[0].radius*suns[0].radius/
                   (suns[0].radius*suns[0].radius+R_Sun2CloudboxTop*R_Sun2CloudboxTop);



  } else {
    // set number of azimuth angle to 1
    N_aa = 1;
  }

  init_ifield(cloudbox_field,
                f_grid,
                cloudbox_limits,
                za_grid.nelem(),
                N_aa);


  Vector albedo(f_grid.nelem(), 0.);
  Numeric btemp;

const Numeric surf_alt = surface_field.single_value(Surf::Key::h, lat_true[0], lon_true[0]);

  if (inc_angle<0 || inc_angle>90) {
    surf_albedoCalc(ws,
                    albedo,
                    btemp,
                    surface_rtprop_agenda,
                    f_grid,
                    za_grid,
                    surf_alt);
  } else {
    surf_albedoCalcSingleAngle(ws,
                               albedo,
                               btemp,
                               surface_rtprop_agenda,
                               f_grid,
                               surf_alt,
                               inc_angle);
  }

  run_cdisort(ws,
              cloudbox_field,
              disort_aux,
              f_grid,
              p_field(joker, 0, 0),
              z_grid,
              surf_alt,
              t_field(joker, 0, 0),
              vmr_field(joker, joker, 0, 0),
              pnd_field(joker, joker, 0, 0),
              abs_species,
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
              intensity_correction);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void spectral_radiance_fieldDisortClearsky(const const Workspace& ws,
                    // WS Output:
                    Tensor7& spectral_radiance_field,
                    ArrayOfMatrix& disort_aux,
                    // WS Input
                    const Index& atmfields_checked,
                    const Index& atmgeom_checked,
                    const Agenda& propmat_clearsky_agenda,
                    const Agenda& gas_scattering_agenda,
                    const AtmField& atm_field,
                    const SurfaceField& surface_field,
                    const Vector& lat_true,
                    const Vector& lon_true,
                    const ArrayOfArrayOfSpeciesTag& abs_species,
                    const ArrayOfSun& suns,
                    const Vector& f_grid,
                    const Vector& za_grid,
                    const Vector& aa_grid,
                    const Numeric& surface_skin_t,
                    const Vector& surface_scalar_reflectivity,
                    const Index& gas_scattering_do,
                    const Index& suns_do,
                    const ArrayOfString& disort_aux_vars,
                    const Index& nstreams,
                    const Index& cdisort_quiet,
                    const Index& emission,
                    const Index& intensity_correction) {

  if (3 != 1)
    throw runtime_error(
        "For running DISORT, atmospheric dimensionality "
        "must be 1.\n");

  // Set cloudbox to cover complete atmosphere
  Index cloudbox_on;
  ArrayOfIndex cloudbox_limits;
  const Index cloudbox_checked = 1;
  //
  ARTS_ASSERT(false)
  /*
  cloudboxSetFullAtm(cloudbox_on,
                     cloudbox_limits,
                     atm_field,
                     0.,
                     );
*/

  // Create data matching no particles
  Tensor4 pnd_field;
  ArrayOfTensor4 dpnd_field_dx;
  ArrayOfArrayOfSingleScatteringData scat_data;
  const Index scat_data_checked = 1;
  //
  pnd_fieldZero(pnd_field,
                dpnd_field_dx,
                scat_data,
                f_grid,
                cloudbox_limits,
                ArrayOfRetrievalQuantity(0));

  Matrix optical_depth_dummy;

  cloudbox_fieldDisort(ws,
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
                 pnd_field,
                 atm_field,
                 surface_field,
                 lat_true,
                 lon_true,
                 abs_species,
                 scat_data,
                 suns,
                 f_grid,
                 za_grid,
                 aa_grid,
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
                 intensity_correction);

}

/* Workspace method: Doxygen documentation will be auto-generated */
void spectral_irradiance_fieldDisort(const Workspace& ws,
                // WS Output:
                Tensor5& spectral_irradiance_field,
                ArrayOfMatrix& disort_aux,
                // WS Input
                const Index& atmfields_checked,
                const Index& atmgeom_checked,
                const Index& scat_data_checked,
                const Agenda& propmat_clearsky_agenda,
                const Agenda& gas_scattering_agenda,
                const Tensor4& pnd_field,
                const AtmField& atm_field,
                const SurfaceField& surface_field,
                const Vector& lat_true,
                const Vector& lon_true,
                const ArrayOfArrayOfSpeciesTag& abs_species,
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const ArrayOfSun& suns,
                const Vector& f_grid,
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
                const Index& intensity_correction) {
  // FIXME: REQUIRES REGULAR GRIDS
  Vector z_grid, lat_grid, lon_grid;
  Tensor3 t_field, p_field, wind_u_field;
  Tensor4 vmr_field;
  ARTS_USER_ERROR("ERROR")
//const auto& z_grid = atm_field.grid[0];
//const auto& p_field = atm_field[Atm::Key::p].get<const Tensor3&>();
//const auto& t_field = atm_field[Atm::Key::t].get<const Tensor3&>();
//const auto vmr_field = Atm::extract_specs_content(atm_field, abs_species);

  // Set cloudbox to cover complete atmosphere
  Index cloudbox_on;
  ArrayOfIndex cloudbox_limits;
  ARTS_ASSERT(false)
  /*
  cloudboxSetFullAtm(cloudbox_on,
                     cloudbox_limits,
                     atm_field,
                     0.,
                     );
*/
  const Index Nf = f_grid.nelem();
  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;

  check_disort_irradiance_input(atmfields_checked,
                                atmgeom_checked,
                                scat_data_checked,
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
    cloudboxtop_pos = {z_grid[cloudbox_limits[1]], lat_true[0], lon_true[0]};

    // calculate local position of sun at top of cloudbox
    ARTS_ASSERT(false)
    /*
    rte_losGeometricFromRtePosToRtePos2(sun_rte_los,
                                        lat_grid,
                                        lon_grid,
                                        refellipsoid,
                                        cloudboxtop_pos,
                                        sun_pos,
                                        );
*/
    //FIXME: IF we want to be correct and include refraction, we must calculate the
    // local position of sun via ppathFromRtePos2. The question is, is this needed,
    // because DISORT does not handle refraction at all.

    // Check if sun is above horizon, if not switch it off
    if (sun_rte_los[0] >= 90) {
      sun_on = 0;
    }

    //get the cloudbox top distance to earth center.
    Numeric R_TOA ;//= refell2r(refellipsoid,
                     //        lat_true[0]) +
                    //cloudboxtop_pos[0];
ARTS_USER_ERROR("ERROR")

    //get the distance between sun and cloudbox top
    Numeric R_Sun2CloudboxTop;
//    distance3D(R_Sun2CloudboxTop,
  //             R_TOA,
    //           lat_true[0],
      //         lon_true[0],
        //       sun_pos[0],
          //     sun_pos[1],
            //   sun_pos[2]);
ARTS_USER_ERROR("ERROR")

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
                   p_field(joker, 0, 0),
                   z_grid,
                   surface_field.single_value(Surf::Key::h, 0, 0),  // FIXME: Must know lat and lon
                   t_field(joker, 0, 0),
                   vmr_field(joker, joker, 0, 0),
                   pnd_field(joker, joker, 0, 0),
                   abs_species,
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
                   intensity_correction);
}