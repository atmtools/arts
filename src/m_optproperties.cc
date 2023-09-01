/*!
  \file   m_optproperties.cc
  \author Sreerekha T.R. <rekha@uni-bremen.de>, 
          Claudia Emde <claudia.emde@dlr.de>
          Cory Davies <cory@met.ed.ac.uk>
  \date   Mon Jun 10 11:19:11 2002 
  \brief  This filecontains workspace methods for calculating the optical 
  properties for the radiative transfer. 

  Optical properties are the extinction matrix, absorption vector and
  scattering vector.  The optical properties for the gases can be
  calculated with or without Zeeman effect.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cfloat>
#include <cmath>
#include "array.h"
#include "arts_constants.h"
#include "arts_conversions.h"
#include <workspace.h>
#include "check_input.h"
#include "exceptions.h"
#include "interpolation.h"
#include "interp.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpack_data.h"
#include "montecarlo.h"
#include "optproperties.h"
#include "rtepack.h"
#include "sorting.h"
#include "xml_io.h"

inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);

#define PART_TYPE scat_data[i_ss][i_se].ptype
#define F_DATAGRID scat_data[i_ss][i_se].f_grid
#define T_DATAGRID scat_data[i_ss][i_se].T_grid
#define ZA_DATAGRID scat_data[i_ss][i_se].za_grid
#define AA_DATAGRID scat_data[i_ss][i_se].aa_grid
#define PHA_MAT_DATA scat_data[i_ss][i_se].pha_mat_data
#define EXT_MAT_DATA scat_data[i_ss][i_se].ext_mat_data
#define ABS_VEC_DATA scat_data[i_ss][i_se].abs_vec_data

// If particle number density is below this value,
// no transformations will be performed
#define PND_LIMIT 1e-12

/* Workspace method: Doxygen documentation will be auto-generated */
void pha_mat_sptFromData(  // Output:
    Tensor5& pha_mat_spt,
    // Input:
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Vector& za_grid,
    const Vector& aa_grid,
    const Index& za_index,  // propagation directions
    const Index& aa_index,
    const Index& f_index,
    const Vector& f_grid,
    const Numeric& rtp_temperature,
    const Tensor4& pnd_field,
    const Index& scat_p_index,
    const Index& scat_lat_index,
    const Index& scat_lon_index) {

  // Determine total number of scattering elements
  const Index N_se_total = TotalNumberOfElements(scat_data);
  if (N_se_total != pnd_field.nbooks()) {
    std::ostringstream os;
    os << "Total number of scattering elements in scat_data "
       << "inconsistent with size of pnd_field.";
    throw std::runtime_error(os.str());
  }
  // as pha_mat_spt is typically initiallized from pnd_field, this theoretically
  // checks the same as the std::runtime_error above. Still, we keep it to be on the
  // save side.
  ARTS_ASSERT(pha_mat_spt.nshelves() == N_se_total);

  // Check that we don't have scat_data_mono here. Only checking the first
  // scat element, assuming the other elements have been processed in the same
  // manner. That's save against having mono data, but not against having
  // individual elements produced with only a single frequency. This, however,
  // has been checked by scat_data_raw reading routines (ScatSpecies/Element*Add/Read).
  // Unsafe, however, remain when ReadXML is used directly or if scat_data(_raw) is
  // (partly) produced from scat_data_singleTmatrix.
  if (scat_data[0][0].f_grid.nelem() < 2) {
    std::ostringstream os;
    os << "Scattering data seems to be *scat_data_mono* (1 freq point only),\n"
       << "but frequency interpolable data (*scat_data* with >=2 freq points) "
       << "is expected here.";
    throw std::runtime_error(os.str());
  }

  const Index N_ss = scat_data.nelem();

  // Phase matrix in laboratory coordinate system. Dimensions:
  // [frequency, za_inc, aa_inc, 4, 4]
  Tensor5 pha_mat_data_int;

  Index i_se_flat = 0;
  // Loop over scattering species
  for (Index i_ss = 0; i_ss < N_ss; i_ss++) {
    const Index N_se = scat_data[i_ss].nelem();

    // Loop over the included scattering elements
    for (Index i_se = 0; i_se < N_se; i_se++) {
      // If the particle number density at a specific point in the
      // atmosphere for the i_se scattering element is zero, we don't need
      // to do the transfromation!
      if (pnd_field(i_se_flat, scat_p_index, scat_lat_index, scat_lon_index) >
          PND_LIMIT) {
        // First we have to transform the data from the coordinate system
        // used in the database (depending on the kind of ptype) to the
        // laboratory coordinate system.

        // Frequency and temperature interpolation:

        // Container for data at one frequency and one temperature.
        pha_mat_data_int.resize(PHA_MAT_DATA.nshelves(),
                                PHA_MAT_DATA.nbooks(),
                                PHA_MAT_DATA.npages(),
                                PHA_MAT_DATA.nrows(),
                                PHA_MAT_DATA.ncols());

        // Gridpositions:
        GridPos freq_gp;
        gridpos(freq_gp, F_DATAGRID, f_grid[f_index]);
        GridPos t_gp;
        Vector itw;

        Index ti = -1;

        if (PHA_MAT_DATA.nvitrines() == 1)  // just 1 T_grid element
        {
          ti = 0;
        } else if (rtp_temperature < 0.)  // coding for 'not interpolate, but
                                          // pick one temperature'
        {
          if (rtp_temperature > -10.)  // lowest T-point
          {
            ti = 0;
          } else if (rtp_temperature > -20.)  // highest T-point
          {
            ti = T_DATAGRID.nelem() - 1;
          } else  // median T-point
          {
            ti = T_DATAGRID.nelem() / 2;
          }
        }

        if (ti < 0)  // temperature interpolation
        {
          std::ostringstream os;
          os << "In pha_mat_sptFromData.\n"
             << "The temperature grid of the scattering data does not\n"
             << "cover the atmospheric temperature at cloud location.\n"
             << "The data should include the value T = " << rtp_temperature
             << " K.";
          chk_interpolation_grids(os.str(), T_DATAGRID, rtp_temperature);

          gridpos(t_gp, T_DATAGRID, rtp_temperature);

          // Interpolation weights:
          itw.resize(4);
          interpweights(itw, freq_gp, t_gp);

          for (Index i_za_sca = 0; i_za_sca < PHA_MAT_DATA.nshelves();
               i_za_sca++)
            for (Index i_aa_sca = 0; i_aa_sca < PHA_MAT_DATA.nbooks();
                 i_aa_sca++)
              for (Index i_za_inc = 0; i_za_inc < PHA_MAT_DATA.npages();
                   i_za_inc++)
                for (Index i_aa_inc = 0; i_aa_inc < PHA_MAT_DATA.nrows();
                     i_aa_inc++)
                  for (Index i = 0; i < PHA_MAT_DATA.ncols(); i++)
                    // Interpolation of phase matrix:
                    pha_mat_data_int(
                        i_za_sca, i_aa_sca, i_za_inc, i_aa_inc, i) =
                        interp(itw,
                               PHA_MAT_DATA(joker,
                                            joker,
                                            i_za_sca,
                                            i_aa_sca,
                                            i_za_inc,
                                            i_aa_inc,
                                            i),
                               freq_gp,
                               t_gp);
        } else {
          // Interpolation weights:
          itw.resize(2);
          interpweights(itw, freq_gp);
          for (Index i_za_sca = 0; i_za_sca < PHA_MAT_DATA.nshelves();
               i_za_sca++)
            for (Index i_aa_sca = 0; i_aa_sca < PHA_MAT_DATA.nbooks();
                 i_aa_sca++)
              for (Index i_za_inc = 0; i_za_inc < PHA_MAT_DATA.npages();
                   i_za_inc++)
                for (Index i_aa_inc = 0; i_aa_inc < PHA_MAT_DATA.nrows();
                     i_aa_inc++)
                  for (Index i = 0; i < PHA_MAT_DATA.ncols(); i++)
                    // Interpolation of phase matrix:
                    pha_mat_data_int(
                        i_za_sca, i_aa_sca, i_za_inc, i_aa_inc, i) =
                        interp(itw,
                               PHA_MAT_DATA(joker,
                                            ti,
                                            i_za_sca,
                                            i_aa_sca,
                                            i_za_inc,
                                            i_aa_inc,
                                            i),
                               freq_gp);
        }

        // Do the transformation into the laboratory coordinate system.
        for (Index za_inc_idx = 0; za_inc_idx < za_grid.nelem();
             za_inc_idx++) {
          for (Index aa_inc_idx = 0; aa_inc_idx < aa_grid.nelem();
               aa_inc_idx++) {
            pha_matTransform(
                pha_mat_spt(i_se_flat, za_inc_idx, aa_inc_idx, joker, joker),
                pha_mat_data_int,
                ZA_DATAGRID,
                AA_DATAGRID,
                PART_TYPE,
                za_index,
                aa_index,
                za_inc_idx,
                aa_inc_idx,
                za_grid,
                aa_grid);
          }
        }
      }
      i_se_flat++;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void pha_mat_sptFromDataDOITOpt(  // Output:
    Tensor5& pha_mat_spt,
    // Input:
    const ArrayOfTensor7& pha_mat_sptDOITOpt,
    const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
    const Index& doit_za_grid_size,
    const Vector& aa_grid,
    const Index& za_index,  // propagation directions
    const Index& aa_index,
    const Numeric& rtp_temperature,
    const Tensor4& pnd_field,
    const Index& scat_p_index,
    const Index& scat_lat_index,
    const Index& scat_lon_index) {
  const Index N_se_total = TotalNumberOfElements(scat_data_mono);

  if (N_se_total != pnd_field.nbooks()) {
    std::ostringstream os;
    os << "Total number of scattering elements in scat_data_mono "
       << "inconsistent with size of pnd_field.";
    throw std::runtime_error(os.str());
  }
  // as pha_mat_spt is typically initiallized from pnd_field, this theoretically
  // checks the same as the std::runtime_error above. Still, we keep it to be on the
  // save side.
  ARTS_ASSERT(pha_mat_spt.nshelves() == N_se_total);

  // 3 = 3
  if (pnd_field.ncols() > 1) {
    ARTS_ASSERT(pha_mat_sptDOITOpt.nelem() == N_se_total);
    // Assuming that if the size is o.k. for one scattering element, it will
    // also be o.k. for the other scattering elements.
    ARTS_ASSERT(pha_mat_sptDOITOpt[0].nlibraries() ==
           scat_data_mono[0][0].T_grid.nelem());
    ARTS_ASSERT(pha_mat_sptDOITOpt[0].nvitrines() == doit_za_grid_size);
    ARTS_ASSERT(pha_mat_sptDOITOpt[0].nshelves() == aa_grid.nelem());
    ARTS_ASSERT(pha_mat_sptDOITOpt[0].nbooks() == doit_za_grid_size);
    ARTS_ASSERT(pha_mat_sptDOITOpt[0].npages() == aa_grid.nelem());
  }

  // 3 = 1, only zenith angle required for scattered directions.
  else if (pnd_field.ncols() == 1) {
    //ARTS_ASSERT(is_size(scat_theta, doit_za_grid_size, 1,
    //                doit_za_grid_size, aa_grid.nelem()));

    ARTS_ASSERT(pha_mat_sptDOITOpt.nelem() == TotalNumberOfElements(scat_data_mono));
    // Assuming that if the size is o.k. for one scattering element, it will
    // also be o.k. for the other scattering elements.
    ARTS_ASSERT(pha_mat_sptDOITOpt[0].nlibraries() ==
           scat_data_mono[0][0].T_grid.nelem());
    ARTS_ASSERT(pha_mat_sptDOITOpt[0].nvitrines() == doit_za_grid_size);
    ARTS_ASSERT(pha_mat_sptDOITOpt[0].nshelves() == 1);
    ARTS_ASSERT(pha_mat_sptDOITOpt[0].nbooks() == doit_za_grid_size);
    ARTS_ASSERT(pha_mat_sptDOITOpt[0].npages() == aa_grid.nelem());
  }

  ARTS_ASSERT(doit_za_grid_size > 0);

  // Check that we do indeed have scat_data_mono here. Only checking the first
  // scat element, assuming the other elements have been processed in the same
  // manner. That's save against having scat_data here if that originated from
  // scat_data_raw reading routines (ScatSpecies/Element*Add/Read), it's not safe
  // against data read by ReadXML directly or if scat_data(_raw) has been (partly)
  // produced from scat_data_singleTmatrix. That would be too costly here,
  // though.
  // Also, we can't check here whether data is at the correct frequency since we
  // don't know f_grid and f_index here (we could pass it in, though).
  if (scat_data_mono[0][0].f_grid.nelem() > 1) {
    std::ostringstream os;
    os << "Scattering data seems to be *scat_data* (several freq points),\n"
       << "but *scat_data_mono* (1 freq point only) is expected here.";
    throw std::runtime_error(os.str());
  }

  // Create equidistant zenith angle grid
  Vector za_grid;
  nlinspace(za_grid, 0, 180, doit_za_grid_size);

  const Index N_ss = scat_data_mono.nelem();

  if (4 > 4 || 4 < 1) {
    throw std::runtime_error(
        "The dimension of the stokes vector \n"
        "must be 1,2,3 or 4");
  }

  GridPos T_gp;
  Vector itw(2);

  // Initialisation
  pha_mat_spt = 0.;

  Index i_se_flat = 0;

  // Do the transformation into the laboratory coordinate system.
  for (Index i_ss = 0; i_ss < N_ss; i_ss++) {
    const Index N_se = scat_data_mono[i_ss].nelem();

    for (Index i_se = 0; i_se < N_se; i_se++) {
      // If the particle number density at a specific point in the
      // atmosphere for the i_se scattering element is zero, we don't need
      // to do the transformation!
      if (pnd_field(i_se_flat, scat_p_index, scat_lat_index, scat_lon_index) >
          PND_LIMIT)  //TRS
      {
        Index nT = scat_data_mono[i_ss][i_se].pha_mat_data.nvitrines();
        Index ti = -1;

        if (nT == 1)  // just 1 T_grid element
        {
          ti = 0;
        } else if (rtp_temperature < 0.)  // coding for 'not interpolate, but
                                          // pick one temperature'
        {
          if (rtp_temperature > -10.)  // lowest T-point
          {
            ti = 0;
          } else if (rtp_temperature > -20.)  // highest T-point
          {
            ti = nT - 1;
          } else  // median T-point
          {
            ti = nT / 2;
          }
        } else {
          std::ostringstream os;
          os << "In pha_mat_sptFromDataDOITOpt.\n"
             << "The temperature grid of the scattering data does not\n"
             << "cover the atmospheric temperature at cloud location.\n"
             << "The data should include the value T = " << rtp_temperature
             << " K.";
          chk_interpolation_grids(
              os.str(), scat_data_mono[i_ss][i_se].T_grid, rtp_temperature);

          // Gridpositions:
          gridpos(T_gp, scat_data_mono[i_ss][i_se].T_grid, rtp_temperature);
          // Interpolation weights:
          interpweights(itw, T_gp);
        }

        for (Index za_inc_idx = 0; za_inc_idx < doit_za_grid_size;
             za_inc_idx++) {
          for (Index aa_inc_idx = 0; aa_inc_idx < aa_grid.nelem();
               aa_inc_idx++) {
            if (ti < 0)  // Temperature interpolation
            {
              for (Index i = 0; i < 4; i++) {
                for (Index j = 0; j < 4; j++) {
                  pha_mat_spt(i_se_flat, za_inc_idx, aa_inc_idx, i, j) =
                      interp(itw,
                             pha_mat_sptDOITOpt[i_se_flat](joker,
                                                           za_index,
                                                           aa_index,
                                                           za_inc_idx,
                                                           aa_inc_idx,
                                                           i,
                                                           j),
                             T_gp);
                }
              }
            } else {
              pha_mat_spt(i_se_flat, za_inc_idx, aa_inc_idx, joker, joker) =
                  pha_mat_sptDOITOpt[i_se_flat](ti,
                                                za_index,
                                                aa_index,
                                                za_inc_idx,
                                                aa_inc_idx,
                                                joker,
                                                joker);
            }
          }
        }
      }  // TRS

      i_se_flat++;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void opt_prop_sptFromData(  // Output and Input:
    ArrayOfPropmatVector& ext_mat_spt,
    ArrayOfStokvecVector& abs_vec_spt,
    // Input:
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Vector& za_grid,
    const Vector& aa_grid,
    const Index& za_index,  // propagation directions
    const Index& aa_index,
    const Index& f_index,
    const Vector& f_grid,
    const Numeric& rtp_temperature,
    const Tensor4& pnd_field,
    const Index& scat_p_index,
    const Index& scat_lat_index,
    const Index& scat_lon_index) {
  const Index N_ss = scat_data.nelem();
  const Numeric za_sca = za_grid[za_index];
  const Numeric aa_sca = aa_grid[aa_index];

  DEBUG_ONLY(const Index N_se_total = TotalNumberOfElements(scat_data);
             if (N_ss) {
               ARTS_ASSERT(ext_mat_spt[0].nelem() == N_se_total);
               ARTS_ASSERT(abs_vec_spt[0].nelem() == N_se_total);
             });

  // Check that we don't have scat_data_mono here. Only checking the first
  // scat element, assuming the other elements have been processed in the same
  // manner. That's save against having mono data, but not against having
  // individual elements produced with only a single frequency. This, however,
  // has been checked by scat_data_raw reading routines (ScatSpecies/Element*Add/Read).
  // Unsafe, however, remain when ReadXML is used directly or if scat_data(_raw) is
  // (partly) produced from scat_data_singleTmatrix.
  if (scat_data[0][0].f_grid.nelem() < 2) {
    std::ostringstream os;
    os << "Scattering data seems to be *scat_data_mono* (1 freq point only),\n"
       << "but frequency interpolable data (*scat_data* with >=2 freq points) "
       << "is expected here.";
    throw std::runtime_error(os.str());
  }

  // Phase matrix in laboratory coordinate system. Dimensions:
  // [frequency, za_inc, aa_inc, 4, 4]
  Tensor3 ext_mat_data_int;
  Tensor3 abs_vec_data_int;

  // Initialisation
  for (auto& pm: ext_mat_spt) pm = 0.;
  for (auto& sv: abs_vec_spt) sv = 0.;

  Index i_se_flat = 0;
  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss < N_ss; i_ss++) {
    const Index N_se = scat_data[i_ss].nelem();

    // Loop over the included scattering elements
    for (Index i_se = 0; i_se < N_se; i_se++) {
      // If the particle number density at a specific point in the
      // atmosphere for the i_se scattering element is zero, we don't need
      // to do the transformation

      if (pnd_field(i_se_flat, scat_p_index, scat_lat_index, scat_lon_index) >
          PND_LIMIT) {
        // First we have to transform the data from the coordinate system
        // used in the database (depending on the kind of ptype) to the
        // laboratory coordinate system.

        // Frequency interpolation:

        // The data is interpolated on one frequency.
        //
        // Resize the variables for the interpolated data:
        //
        ext_mat_data_int.resize(
            EXT_MAT_DATA.npages(), EXT_MAT_DATA.nrows(), EXT_MAT_DATA.ncols());
        //
        abs_vec_data_int.resize(
            ABS_VEC_DATA.npages(), ABS_VEC_DATA.nrows(), ABS_VEC_DATA.ncols());

        // Gridpositions:
        GridPos freq_gp;
        gridpos(freq_gp, F_DATAGRID, f_grid[f_index]);
        GridPos t_gp;
        Vector itw;

        if (T_DATAGRID.nelem() > 1) {
          std::ostringstream os;
          os << "In opt_prop_sptFromData.\n"
             << "The temperature grid of the scattering data does not\n"
             << "cover the atmospheric temperature at cloud location.\n"
             << "The data should include the value T = " << rtp_temperature
             << " K.";
          chk_interpolation_grids(os.str(), T_DATAGRID, rtp_temperature);

          gridpos(t_gp, T_DATAGRID, rtp_temperature);

          // Interpolation weights:
          itw.resize(4);
          interpweights(itw, freq_gp, t_gp);

          for (Index i_za_sca = 0; i_za_sca < EXT_MAT_DATA.npages();
               i_za_sca++) {
            for (Index i_aa_sca = 0; i_aa_sca < EXT_MAT_DATA.nrows();
                 i_aa_sca++) {
              //
              // Interpolation of extinction matrix:
              //
              for (Index i = 0; i < EXT_MAT_DATA.ncols(); i++) {
                ext_mat_data_int(i_za_sca, i_aa_sca, i) =
                    interp(itw,
                           EXT_MAT_DATA(joker, joker, i_za_sca, i_aa_sca, i),
                           freq_gp,
                           t_gp);
              }
            }
          }

          for (Index i_za_sca = 0; i_za_sca < ABS_VEC_DATA.npages();
               i_za_sca++) {
            for (Index i_aa_sca = 0; i_aa_sca < ABS_VEC_DATA.nrows();
                 i_aa_sca++) {
              //
              // Interpolation of absorption vector:
              //
              for (Index i = 0; i < ABS_VEC_DATA.ncols(); i++) {
                abs_vec_data_int(i_za_sca, i_aa_sca, i) =
                    interp(itw,
                           ABS_VEC_DATA(joker, joker, i_za_sca, i_aa_sca, i),
                           freq_gp,
                           t_gp);
              }
            }
          }
        } else {
          // Interpolation weights:
          itw.resize(2);
          interpweights(itw, freq_gp);

          for (Index i_za_sca = 0; i_za_sca < EXT_MAT_DATA.npages();
               i_za_sca++) {
            for (Index i_aa_sca = 0; i_aa_sca < EXT_MAT_DATA.nrows();
                 i_aa_sca++) {
              //
              // Interpolation of extinction matrix:
              //
              for (Index i = 0; i < EXT_MAT_DATA.ncols(); i++) {
                ext_mat_data_int(i_za_sca, i_aa_sca, i) =
                    interp(itw,
                           EXT_MAT_DATA(joker, 0, i_za_sca, i_aa_sca, i),
                           freq_gp);
              }
            }
          }

          for (Index i_za_sca = 0; i_za_sca < ABS_VEC_DATA.npages();
               i_za_sca++) {
            for (Index i_aa_sca = 0; i_aa_sca < ABS_VEC_DATA.nrows();
                 i_aa_sca++) {
              //
              // Interpolation of absorption vector:
              //
              for (Index i = 0; i < ABS_VEC_DATA.ncols(); i++) {
                abs_vec_data_int(i_za_sca, i_aa_sca, i) =
                    interp(itw,
                           ABS_VEC_DATA(joker, 0, i_za_sca, i_aa_sca, i),
                           freq_gp);
              }
            }
          }
        }

        //
        // Do the transformation into the laboratory coordinate system.
        //
        // Extinction matrix:
        //
        ext_matTransform(ext_mat_spt[i_se_flat],
                         ext_mat_data_int,
                         ZA_DATAGRID,
                         AA_DATAGRID,
                         PART_TYPE,
                         za_sca,
                         aa_sca);
        //
        // Absorption vector:
        //
        abs_vecTransform(abs_vec_spt[i_se_flat],
                         abs_vec_data_int,
                         ZA_DATAGRID,
                         AA_DATAGRID,
                         PART_TYPE,
                         za_sca,
                         aa_sca);
      }

      i_se_flat++;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void opt_prop_sptFromScat_data(  // Output and Input:
    ArrayOfPropmatVector& ext_mat_spt,
    ArrayOfStokvecVector& abs_vec_spt,
    // Input:
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Index& scat_data_checked,
    const Vector& za_grid,
    const Vector& aa_grid,
    const Index& za_index,  // propagation directions
    const Index& aa_index,
    const Index& f_index,
    const Numeric& rtp_temperature,
    const Tensor4& pnd_field,
    const Index& scat_p_index,
    const Index& scat_lat_index,
    const Index& scat_lon_index) {
  if (scat_data_checked != 1)
    throw std::runtime_error(
        "The scattering data must be flagged to have "
        "passed a consistency check (scat_data_checked=1).");

  const Index N_ss = scat_data.nelem();
  const Numeric za_sca = za_grid[za_index];
  const Numeric aa_sca = aa_grid[aa_index];

  DEBUG_ONLY(const Index N_se_total = TotalNumberOfElements(scat_data);)
  ARTS_ASSERT(ext_mat_spt.nelem() == N_se_total);
  ARTS_ASSERT(abs_vec_spt.nelem() == N_se_total);

  // Phase matrix in laboratory coordinate system. Dimensions:
  // [frequency, za_inc, aa_inc, 4, 4]
  Tensor3 ext_mat_data_int;
  Tensor3 abs_vec_data_int;

  // Initialisation
  for (auto& pm : ext_mat_spt) pm = 0.;
  for (auto& sv : abs_vec_spt) sv = 0.;

  Index this_f_index;

  Index i_se_flat = 0;
  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss < N_ss; i_ss++) {
    const Index N_se = scat_data[i_ss].nelem();

    // Loop over the included scattering elements
    for (Index i_se = 0; i_se < N_se; i_se++) {
      // If the particle number density at a specific point in the
      // atmosphere for the i_se scattering element is zero, we don't need
      // to do the transformation

      if (pnd_field(i_se_flat, scat_p_index, scat_lat_index, scat_lon_index) >
          PND_LIMIT) {
        // First we have to transform the data from the coordinate system
        // used in the database (depending on the kind of ptype) to the
        // laboratory coordinate system.

        // Resize the variables for the interpolated data (1freq, 1T):
        ext_mat_data_int.resize(
            EXT_MAT_DATA.npages(), EXT_MAT_DATA.nrows(), EXT_MAT_DATA.ncols());
        abs_vec_data_int.resize(
            ABS_VEC_DATA.npages(), ABS_VEC_DATA.nrows(), ABS_VEC_DATA.ncols());

        // Gridpositions and interpolation weights;
        GridPos t_gp;
        Vector itw;
        if (EXT_MAT_DATA.nbooks() > 1 || ABS_VEC_DATA.nbooks() > 1) {
          std::ostringstream os;
          os << "In opt_prop_sptFromScat_data.\n"
             << "The temperature grid of the scattering data does not\n"
             << "cover the atmospheric temperature at cloud location.\n"
             << "The data should include the value T = " << rtp_temperature
             << " K.";
          chk_interpolation_grids(os.str(), T_DATAGRID, rtp_temperature);

          gridpos(t_gp, T_DATAGRID, rtp_temperature);

          // Interpolation weights:
          itw.resize(2);
          interpweights(itw, t_gp);
        }

        // Frequency extraction and temperature interpolation

        if (EXT_MAT_DATA.nshelves() == 1)
          this_f_index = 0;
        else
          this_f_index = f_index;

        if (EXT_MAT_DATA.nbooks() > 1) {
          // Interpolation of extinction matrix:
          for (Index i_za_sca = 0; i_za_sca < EXT_MAT_DATA.npages();
               i_za_sca++) {
            for (Index i_aa_sca = 0; i_aa_sca < EXT_MAT_DATA.nrows();
                 i_aa_sca++) {
              for (Index i = 0; i < EXT_MAT_DATA.ncols(); i++) {
                ext_mat_data_int(i_za_sca, i_aa_sca, i) = interp(
                    itw,
                    EXT_MAT_DATA(this_f_index, joker, i_za_sca, i_aa_sca, i),
                    t_gp);
              }
            }
          }
        } else {
          ext_mat_data_int = EXT_MAT_DATA(this_f_index, 0, joker, joker, joker);
          /*
                  for (Index i_za_sca = 0; i_za_sca < EXT_MAT_DATA.npages();
                       i_za_sca++)
                  {
                      for(Index i_aa_sca = 0; i_aa_sca < EXT_MAT_DATA.nrows();
                          i_aa_sca++)
                      {
                          for (Index i = 0; i < EXT_MAT_DATA.ncols(); i++)
                          {
                              ext_mat_data_int(i_za_sca, i_aa_sca, i) =
                                EXT_MAT_DATA(this_f_index, 0,
                                                 i_za_sca, i_aa_sca, i);
                          }
                      }
                  } */
        }

        if (ABS_VEC_DATA.nshelves() == 1)
          this_f_index = 0;
        else
          this_f_index = f_index;

        if (ABS_VEC_DATA.nbooks() > 1) {
          // Interpolation of absorption vector:
          for (Index i_za_sca = 0; i_za_sca < ABS_VEC_DATA.npages();
               i_za_sca++) {
            for (Index i_aa_sca = 0; i_aa_sca < ABS_VEC_DATA.nrows();
                 i_aa_sca++) {
              for (Index i = 0; i < ABS_VEC_DATA.ncols(); i++) {
                abs_vec_data_int(i_za_sca, i_aa_sca, i) = interp(
                    itw,
                    ABS_VEC_DATA(this_f_index, joker, i_za_sca, i_aa_sca, i),
                    t_gp);
              }
            }
          }
        } else {
          abs_vec_data_int = ABS_VEC_DATA(this_f_index, 0, joker, joker, joker);
          /*
                  for (Index i_za_sca = 0; i_za_sca < ABS_VEC_DATA.npages();
                       i_za_sca++)
                  {
                      for(Index i_aa_sca = 0; i_aa_sca < ABS_VEC_DATA.nrows();
                          i_aa_sca++)
                      {
                          for (Index i = 0; i < ABS_VEC_DATA.ncols(); i++)
                          {
                              abs_vec_data_int(i_za_sca, i_aa_sca, i) =
                                ABS_VEC_DATA(this_f_index, 0,
                                                 i_za_sca, i_aa_sca, i);
                          }
                      }
                  } */
        }

        //
        // Do the transformation into the laboratory coordinate system.
        //
        // Extinction matrix:
        ext_matTransform(ext_mat_spt[i_se_flat],
                         ext_mat_data_int,
                         ZA_DATAGRID,
                         AA_DATAGRID,
                         PART_TYPE,
                         za_sca,
                         aa_sca);
        // Absorption vector:
        abs_vecTransform(abs_vec_spt[i_se_flat],
                         abs_vec_data_int,
                         ZA_DATAGRID,
                         AA_DATAGRID,
                         PART_TYPE,
                         za_sca,
                         aa_sca);
      }
      i_se_flat++;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void opt_prop_bulkCalc(  // Output and Input:
    PropmatVector& ext_mat,
    StokvecVector& abs_vec,
    // Input:
    const ArrayOfPropmatVector& ext_mat_spt,
    const ArrayOfStokvecVector& abs_vec_spt,
    const Tensor4& pnd_field,
    const Index& scat_p_index,
    const Index& scat_lat_index,
    const Index& scat_lon_index) {
  Index N_se = abs_vec_spt.nelem();

  if (ext_mat_spt.nelem() not_eq N_se) {
    std::ostringstream os;
    os << "Number of scattering elements in *abs_vec_spt* and *ext_mat_spt*\n"
       << "does not agree.";
    throw std::runtime_error(os.str());
  }

  ext_mat = PropmatVector(1, Propmat{0, 0, 0, 0, 0, 0, 0});
  abs_vec = StokvecVector(1, Stokvec{0, 0, 0, 0});

  PropmatVector ext_mat_part(1, Propmat{0, 0, 0, 0, 0, 0, 0});
  StokvecVector abs_vec_part(1, Stokvec{0, 0, 0, 0});

  // this is the loop over the different scattering elements
  for (Index l = 0; l < N_se; l++) {
    abs_vec_part[0] += pnd_field(l, scat_p_index, scat_lat_index, scat_lon_index) * abs_vec_spt[l][0];
    ext_mat_part[0] += pnd_field(l, scat_p_index, scat_lat_index, scat_lon_index) * ext_mat_spt[l][0];
  }

  abs_vec += abs_vec_part;
  ext_mat += ext_mat_part;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ext_matAddGas(PropmatVector& ext_mat,
                   const PropmatVector& propmat_clearsky) {
  const Index f_dim = ext_mat.nelem();

  ARTS_USER_ERROR_IF (f_dim != propmat_clearsky.nelem(),
        "Frequency dimension of ext_mat and propmat_clearsky\n"
        "are inconsistent in ext_matAddGas.");

  ext_mat += propmat_clearsky;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_vecAddGas(StokvecVector& abs_vec,
                   const PropmatVector& propmat_clearsky) {
  const Index f_dim = abs_vec.nelem();

  ARTS_USER_ERROR_IF (f_dim != propmat_clearsky.nelem(),
        "Frequency dimension of abs_vec and propmat_clearsky\n"
        "are inconsistent in abs_vecAddGas.");

  for (Index iv=0; iv<f_dim; iv++) {
    abs_vec[iv] += absvec(propmat_clearsky[iv]);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void pha_matCalc(Tensor4& pha_mat,
                 const Tensor5& pha_mat_spt,
                 const Tensor4& pnd_field,
                 const Index& scat_p_index,
                 const Index& scat_lat_index,
                 const Index& scat_lon_index) {
  Index N_se = pha_mat_spt.nshelves();
  Index Nza = pha_mat_spt.nbooks();
  Index Naa = pha_mat_spt.npages();

  pha_mat.resize(Nza, Naa, 4, 4);

  // Initialisation
  pha_mat = 0.0;

  Index ilat = 0;
  Index ilon = 0;
  if (3 > 1) ilat = scat_lat_index;
  if (3 > 2) ilon = scat_lon_index;

  if (3 == 1) {
    // For 1d atmospheres, we additinally integrate the phase matrix over the
    // azimuth, because there is no azimuth dependency of the incoming
    // field.

    Numeric grid_step_size_azimuth = 360. / (Numeric)(Naa - 1) * DEG2RAD;

    // this is a loop over the different scattering elements
    for (Index pt_index = 0; pt_index < N_se; ++pt_index)
      // these are loops over zenith angle and azimuth angle
      for (Index za_index = 0; za_index < Nza; ++za_index)
        for (Index aa_index = 0; aa_index < Naa - 1; ++aa_index)
          // now the last two loops over the stokes dimension.
          for (Index stokes_index_1 = 0; stokes_index_1 < 4;
               ++stokes_index_1)
            for (Index stokes_index_2 = 0; stokes_index_2 < 4;
                 ++stokes_index_2)
              //summation of the product of pnd_field and
              //pha_mat_spt.
              pha_mat(za_index, 0, stokes_index_1, stokes_index_2) +=
                  ((pha_mat_spt(pt_index,
                                za_index,
                                aa_index,
                                stokes_index_1,
                                stokes_index_2) +
                    pha_mat_spt(pt_index,
                                za_index,
                                aa_index + 1,
                                stokes_index_1,
                                stokes_index_2)) /
                   2 * grid_step_size_azimuth *
                   pnd_field(pt_index, scat_p_index, ilat, ilon));
  } else {
    // this is a loop over the different scattering elements
    for (Index pt_index = 0; pt_index < N_se; ++pt_index)
      // these are loops over zenith angle and azimuth angle
      for (Index za_index = 0; za_index < Nza; ++za_index)
        for (Index aa_index = 0; aa_index < Naa; ++aa_index)
          // now the last two loops over the stokes dimension.
          for (Index stokes_index_1 = 0; stokes_index_1 < 4;
               ++stokes_index_1)
            for (Index stokes_index_2 = 0; stokes_index_2 < 4;
                 ++stokes_index_2)
              //summation of the product of pnd_field and
              //pha_mat_spt.
              pha_mat(za_index, aa_index, stokes_index_1, stokes_index_2) +=
                  (pha_mat_spt(pt_index,
                               za_index,
                               aa_index,
                               stokes_index_1,
                               stokes_index_2) *
                   pnd_field(pt_index, scat_p_index, ilat, ilon));
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void scat_dataCheck(  //Input:
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const String& check_type,
    const Numeric& threshold) {
  // FIXME:
  // so far, this works for both scat_data and scat_data_raw. Needs to be
  // adapted, though, once we have WSM that can create Z/K/a with different
  // f/T dimensions than scat_data_single.f/T_grid.

  const Index N_ss = scat_data.nelem();

  // 1) any negative values in Z11, K11, or a1? K11>=a1?
  // 2) scat_data containing any NaN?
  // 3) sca_mat norm sufficiently good (int(Z11)~=K11-a1?)

  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss < N_ss; i_ss++) {
    const Index N_se = scat_data[i_ss].nelem();

    // Loop over the included scattering elements
    for (Index i_se = 0; i_se < N_se; i_se++) {
      for (Index f = 0; f < F_DATAGRID.nelem(); f++) {
        for (Index zai = 0; zai < ABS_VEC_DATA.npages(); zai++)
          for (Index aai = 0; aai < ABS_VEC_DATA.nrows(); aai++) {
            for (Index t = 0; t < T_DATAGRID.nelem(); t++) {
              if (EXT_MAT_DATA(f, t, zai, aai, 0) < 0 ||
                  ABS_VEC_DATA(f, t, zai, aai, 0) < 0) {
                std::ostringstream os;
                os << "Scatt. species #" << i_ss << " element #" << i_se
                   << " contains negative K11 or a1 at f#" << f << ", T#" << t
                   << ", za#" << zai << ", aa#" << aai << "\n";
                throw std::runtime_error(os.str());
              }
              if (EXT_MAT_DATA(f, t, zai, aai, 0) <
                  ABS_VEC_DATA(f, t, zai, aai, 0)) {
                std::ostringstream os;
                os << "Scatt. species #" << i_ss << " element #" << i_se
                   << " has K11<a1 at f#" << f << ", T#" << t << ", za#" << zai
                   << ", aa#" << aai << "\n";
                throw std::runtime_error(os.str());
              }
            }
            // Since allowing pha_mat to have a single T entry only (while
            // T_grid, ext_mat, abs_vec have more), we need a separate T loop
            // for pha_mat
            Index nTpha = PHA_MAT_DATA.nvitrines();
            for (Index t = 0; t < nTpha; t++) {
              for (Index zas = 0; zas < PHA_MAT_DATA.nshelves(); zas++)
                for (Index aas = 0; aas < PHA_MAT_DATA.nbooks(); aas++)
                  if (PHA_MAT_DATA(f, t, zas, aas, zai, aai, 0) < 0) {
                    std::ostringstream os;
                    os << "Scatt. species #" << i_ss << " element #" << i_se
                       << " contains negative Z11 at f#" << f << ", T#" << t
                       << " (of " << nTpha << "), za_sca#" << zas << ", aa_sca#"
                       << aas << ", za_inc#" << zai << ", aa_inc#" << aai
                       << "\n";
                    throw std::runtime_error(os.str());
                  }
            }
          }
      }
    }
  }

  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss < N_ss; i_ss++) {
    const Index N_se = scat_data[i_ss].nelem();

    // Loop over the included scattering elements
    for (Index i_se = 0; i_se < N_se; i_se++) {
      for (Index f = 0; f < F_DATAGRID.nelem(); f++) {
        for (Index zai = 0; zai < ABS_VEC_DATA.npages(); zai++)
          for (Index aai = 0; aai < ABS_VEC_DATA.nrows(); aai++) {
            for (Index t = 0; t < T_DATAGRID.nelem(); t++) {
              for (Index st = 0; st < ABS_VEC_DATA.ncols(); st++)
                if (std::isnan(ABS_VEC_DATA(f, t, zai, aai, st))) {
                  std::ostringstream os;
                  os << "Scatt. species #" << i_ss << " element #" << i_se
                     << " contains NaN in abs_vec at f#" << f << ", T#" << t
                     << ", za#" << zai << ", aa#" << aai << ", stokes #" << st
                     << "\n";
                  throw std::runtime_error(os.str());
                }
              for (Index st = 0; st < EXT_MAT_DATA.ncols(); st++)
                if (std::isnan(EXT_MAT_DATA(f, t, zai, aai, st))) {
                  std::ostringstream os;
                  os << "Scatt. species #" << i_ss << " element #" << i_se
                     << " contains NaN in ext_mat at f#" << f << ", T#" << t
                     << ", za#" << zai << ", aa#" << aai << ", stokes #" << st
                     << "\n";
                  throw std::runtime_error(os.str());
                }
            }
            Index nTpha = PHA_MAT_DATA.nvitrines();
            for (Index t = 0; t < nTpha; t++) {
              for (Index zas = 0; zas < PHA_MAT_DATA.nshelves(); zas++)
                for (Index aas = 0; aas < PHA_MAT_DATA.nbooks(); aas++)
                  for (Index st = 0; st < PHA_MAT_DATA.ncols(); st++)
                    if (std::isnan(
                            PHA_MAT_DATA(f, t, zas, aas, zai, aai, st))) {
                      std::ostringstream os;
                      os << "Scatt. species #" << i_ss << " element #" << i_se
                         << " contains NaN in pha_mat at f#" << f << ", T#" << t
                         << " (of " << nTpha << "), za_sca#" << zas
                         << ", aa_sca#" << aas << ", za_inc#" << zai
                         << ", aa_inc#" << aai << ", stokes #"
                         << "\n";
                      throw std::runtime_error(os.str());
                    }
            }
          }
      }
    }
  }

  if (check_type.toupper() == "ALL") {
    // Loop over the included scattering species
    for (Index i_ss = 0; i_ss < N_ss; i_ss++) {
      const Index N_se = scat_data[i_ss].nelem();

      // Loop over the included scattering elements
      for (Index i_se = 0; i_se < N_se; i_se++)
        if (T_DATAGRID.nelem() == PHA_MAT_DATA.nvitrines()) switch (PART_TYPE) {
            case PTYPE_TOTAL_RND: {
              for (Index f = 0; f < F_DATAGRID.nelem(); f++) {
                for (Index t = 0; t < T_DATAGRID.nelem(); t++) {
                  Numeric Csca = AngIntegrate_trapezoid(
                      PHA_MAT_DATA(f, t, joker, 0, 0, 0, 0), ZA_DATAGRID);
                  Numeric Cext_data = EXT_MAT_DATA(f, t, 0, 0, 0);
                  //Numeric Cabs = Cext_data - Csca;
                  Numeric Cabs_data = ABS_VEC_DATA(f, t, 0, 0, 0);
                  Numeric Csca_data = Cext_data - Cabs_data;

                  //if (abs(Csca/Csca_data-1.)*Csca_data/Cext_data > threshold)
                  // below equivalent to the above
                  // (it's actually the (absolute) albedo deviation!)
                  if (abs(Csca - Csca_data) / Cext_data > threshold) {
                    std::ostringstream os;
                    os << "  Deviations in scat_data too large:\n"
                       << "  scat dev [%] " << 1e2 * Csca / Csca_data - 1e2
                       << " at nominal (actual) albedo of "
                       << Csca_data / Cext_data << " (" << Csca / Cext_data
                       << ").\n"
                       << "  Check entry for scattering element " << i_se
                       << " of scattering species " << i_ss << " at " << f
                       << ".frequency and " << t << ".temperature!\n";
                    throw std::runtime_error(os.str());
                  }
                }
              }
              break;
            }

            case PTYPE_AZIMUTH_RND: {
              for (Index f = 0; f < F_DATAGRID.nelem(); f++) {
                for (Index t = 0; t < T_DATAGRID.nelem(); t++) {
                  for (Index iza = 0; iza < ABS_VEC_DATA.npages(); iza++) {
                    Numeric Csca =
                        2 * AngIntegrate_trapezoid(
                                PHA_MAT_DATA(f, t, joker, joker, iza, 0, 0),
                                ZA_DATAGRID,
                                AA_DATAGRID);
                    Numeric Cext_data = EXT_MAT_DATA(f, t, iza, 0, 0);
                    //Numeric Cabs = Cext_data - Csca;
                    Numeric Cabs_data = ABS_VEC_DATA(f, t, iza, 0, 0);
                    Numeric Csca_data = Cext_data - Cabs_data;

                    //if (abs(Csca/Csca_data-1.)*Csca_data/Cext_data > threshold)
                    // below equivalent to the above
                    // (it's actually the (absolute) albedo deviation!)
                    if (abs(Csca - Csca_data) / Cext_data > threshold) {
                      std::ostringstream os;
                      os << "  Deviations in scat_data too large:\n"
                         << "  scat dev [%] " << 1e2 * Csca / Csca_data - 1e2
                         << " at nominal (actual) albedo of "
                         << Csca_data / Cext_data << " (" << Csca / Cext_data
                         << ").\n"
                         << "  Check entry for scattering element " << i_se
                         << " of scattering species " << i_ss << " at " << f
                         << ". frequency, " << t << ". temperature, and " << iza
                         << ". incident polar angle!\n";
                      throw std::runtime_error(os.str());
                    }
                  }
                }
              }
              break;
            }

            default: {
            }
          }
    }
  } else if (check_type.toupper() == "SANE") {
  } else {
    std::ostringstream os;
    os << "Invalid value for argument *check_type*: '" << check_type << "'.\n";
    os << "Valid values are 'all' or 'none'.";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void scat_dataCalc(ArrayOfArrayOfSingleScatteringData& scat_data,
                   const ArrayOfArrayOfSingleScatteringData& scat_data_raw,
                   const Vector& f_grid,
                   const Index& interp_order)
// FIXME: when we allow K, a, Z to be on different f and T grids, their use in
// the scatt solvers needs to be reviewed again and adaptedto this!
{
  //Extrapolation factor:
  //const Numeric extpolfac = 0.5;

  Index nf = f_grid.nelem();

  // Check, whether single scattering data contains the right frequencies:
  // The check was changed to allow extrapolation at the boundaries of the
  // frequency grid.
  const String which_interpolation = "scat_data_raw.f_grid to f_grid";
  for (Index i_ss = 0; i_ss < scat_data_raw.nelem(); i_ss++) {
    for (Index i_se = 0; i_se < scat_data_raw[i_ss].nelem(); i_se++) {
      // Check for the special case that ssd.f_grid f_grid have only one
      // element. If identical, that's  fine. If not, throw error.
      if (scat_data_raw[i_ss][i_se].f_grid.nelem() == 1 && nf == 1)
        if (!is_same_within_epsilon(scat_data_raw[i_ss][i_se].f_grid[0],
                                    f_grid[0],
                                    2 * DBL_EPSILON)) {
          std::ostringstream os;
          os << "There is a problem with the grids for the following "
             << "interpolation:\n"
             << which_interpolation << "\n"
             << "If original grid has only 1 element, the new grid must also have\n"
             << "only a single element and hold the same value as the original grid.";
          throw std::runtime_error(os.str());
        }

      // check with extrapolation
      chk_interpolation_grids(which_interpolation,
                              scat_data_raw[i_ss][i_se].f_grid,
                              f_grid,
                              interp_order);
    }
  }

  //Initialise scat_data
  scat_data.resize(scat_data_raw.nelem());

  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss < scat_data_raw.nelem(); i_ss++) {
    const Index N_se = scat_data_raw[i_ss].nelem();

    //Initialise scat_data
    scat_data[i_ss].resize(N_se);

    // Loop over the included scattering elements
    for (Index i_se = 0; i_se < N_se; i_se++) {
      //Stuff that doesn't need interpolating
      PART_TYPE = scat_data_raw[i_ss][i_se].ptype;
      F_DATAGRID = f_grid;
      T_DATAGRID = scat_data_raw[i_ss][i_se].T_grid;
      ZA_DATAGRID = scat_data_raw[i_ss][i_se].za_grid;
      AA_DATAGRID = scat_data_raw[i_ss][i_se].aa_grid;

      //Sizing SSD data containers
      PHA_MAT_DATA.resize(nf,
                          scat_data_raw[i_ss][i_se].pha_mat_data.nvitrines(),
                          scat_data_raw[i_ss][i_se].pha_mat_data.nshelves(),
                          scat_data_raw[i_ss][i_se].pha_mat_data.nbooks(),
                          scat_data_raw[i_ss][i_se].pha_mat_data.npages(),
                          scat_data_raw[i_ss][i_se].pha_mat_data.nrows(),
                          scat_data_raw[i_ss][i_se].pha_mat_data.ncols());
      EXT_MAT_DATA.resize(nf,
                          scat_data_raw[i_ss][i_se].ext_mat_data.nbooks(),
                          scat_data_raw[i_ss][i_se].ext_mat_data.npages(),
                          scat_data_raw[i_ss][i_se].ext_mat_data.nrows(),
                          scat_data_raw[i_ss][i_se].ext_mat_data.ncols());
      ABS_VEC_DATA.resize(nf,
                          scat_data_raw[i_ss][i_se].abs_vec_data.nbooks(),
                          scat_data_raw[i_ss][i_se].abs_vec_data.npages(),
                          scat_data_raw[i_ss][i_se].abs_vec_data.nrows(),
                          scat_data_raw[i_ss][i_se].abs_vec_data.ncols());

      const bool single_se_fgrid =
          (scat_data_raw[i_ss][i_se].f_grid.nelem() == 1);
      if (!single_se_fgrid) {
        // Gridpositions:
        const auto lag_freq=my_interp::lagrange_interpolation_list<LagrangeInterpolation>(f_grid, scat_data_raw[i_ss][i_se].f_grid, interp_order);
        const auto itw = interpweights(lag_freq);

        //Phase matrix data
        for (Index t_index = 0;
             t_index < scat_data_raw[i_ss][i_se].pha_mat_data.nvitrines();
             t_index++) {
          for (Index i_za_sca = 0;
               i_za_sca < scat_data_raw[i_ss][i_se].pha_mat_data.nshelves();
               i_za_sca++) {
            for (Index i_aa_sca = 0;
                 i_aa_sca < scat_data_raw[i_ss][i_se].pha_mat_data.nbooks();
                 i_aa_sca++) {
              for (Index i_za_inc = 0;
                   i_za_inc < scat_data_raw[i_ss][i_se].pha_mat_data.npages();
                   i_za_inc++) {
                for (Index i_aa_inc = 0;
                     i_aa_inc < scat_data_raw[i_ss][i_se].pha_mat_data.nrows();
                     i_aa_inc++) {
                  for (Index i = 0;
                       i < scat_data_raw[i_ss][i_se].pha_mat_data.ncols();
                       i++) {
                    reinterp(scat_data[i_ss][i_se].pha_mat_data(joker,
                                                                t_index,
                                                                i_za_sca,
                                                                i_aa_sca,
                                                                i_za_inc,
                                                                i_aa_inc,
                                                                i),
                             scat_data_raw[i_ss][i_se].pha_mat_data(joker,
                                                                    t_index,
                                                                    i_za_sca,
                                                                    i_aa_sca,
                                                                    i_za_inc,
                                                                    i_aa_inc,
                                                                    i),
                             itw,
                             lag_freq);
                  }
                }
              }
            }
          }
        }

        //Extinction matrix data
        for (Index t_index = 0;
             t_index < scat_data_raw[i_ss][i_se].ext_mat_data.nbooks();
             t_index++) {
          for (Index i_za_sca = 0;
               i_za_sca < scat_data_raw[i_ss][i_se].ext_mat_data.npages();
               i_za_sca++) {
            for (Index i_aa_sca = 0;
                 i_aa_sca < scat_data_raw[i_ss][i_se].ext_mat_data.nrows();
                 i_aa_sca++) {
              for (Index i = 0;
                   i < scat_data_raw[i_ss][i_se].ext_mat_data.ncols();
                   i++) {
                reinterp(scat_data[i_ss][i_se].ext_mat_data(
                             joker, t_index, i_za_sca, i_aa_sca, i),
                         scat_data_raw[i_ss][i_se].ext_mat_data(
                             joker, t_index, i_za_sca, i_aa_sca, i),
                         itw,
                         lag_freq);
              }
            }
          }
        }

        //Absorption vector data
        for (Index t_index = 0;
             t_index < scat_data_raw[i_ss][i_se].abs_vec_data.nbooks();
             t_index++) {
          for (Index i_za_sca = 0;
               i_za_sca < scat_data_raw[i_ss][i_se].abs_vec_data.npages();
               i_za_sca++) {
            for (Index i_aa_sca = 0;
                 i_aa_sca < scat_data_raw[i_ss][i_se].abs_vec_data.nrows();
                 i_aa_sca++) {
              for (Index i = 0;
                   i < scat_data_raw[i_ss][i_se].abs_vec_data.ncols();
                   i++) {
                reinterp(scat_data[i_ss][i_se].abs_vec_data(
                             joker, t_index, i_za_sca, i_aa_sca, i),
                         scat_data_raw[i_ss][i_se].abs_vec_data(
                             joker, t_index, i_za_sca, i_aa_sca, i),
                         itw,
                         lag_freq);
              }
            }
          }
        }
      } else {
        ARTS_ASSERT(nf == 1);
        // we do only have one f_grid value in old and new data (and they have
        // been confirmed to be the same), hence only need to copy over/reassign
        // the data.
        scat_data[i_ss][i_se].pha_mat_data =
            scat_data_raw[i_ss][i_se].pha_mat_data;
        scat_data[i_ss][i_se].ext_mat_data =
            scat_data_raw[i_ss][i_se].ext_mat_data;
        scat_data[i_ss][i_se].abs_vec_data =
            scat_data_raw[i_ss][i_se].abs_vec_data;
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void scat_dataReduceT(ArrayOfArrayOfSingleScatteringData& scat_data,
                      const Index& i_ss,
                      const Numeric& T,
                      const Index& interp_order,
                      const Index& phamat_only,
                      const Numeric& threshold) {
  // We are directly acting on the scat_data entries, modifying them
  // individually. That is, we don't need to resize these arrays. Only the
  // pha_mat and probably ext_mat and abs_vec Tensors (in the latter case also
  // T_grid!).

  // Check that species i_ss exists at all in scat_data
  const Index nss = scat_data.nelem();
  if (nss <= i_ss) {
    std::ostringstream os;
    os << "Can not T-reduce scattering species #" << i_ss << ".\n"
       << "*scat_data* contains only " << nss << " scattering species.";
    throw std::runtime_error(os.str());
  }

  // Loop over the included scattering elements
  for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++) {
    // At very first check validity of the scatt elements ptype (so far we only
    // handle PTYPE_TOTAL_RND and PTYPE_AZIMUTH_RND).
    if (PART_TYPE != PTYPE_TOTAL_RND and PART_TYPE != PTYPE_AZIMUTH_RND) {
      std::ostringstream os;
      os << "Only ptypes " << PTYPE_TOTAL_RND << " and " << PTYPE_AZIMUTH_RND
         << " can be handled.\n"
         << "Scattering element #" << i_se << " has ptype " << PART_TYPE << ".";
      throw std::runtime_error(os.str());
    }

    // If ssd.T_grid already has only a single point, we do nothing.
    // This is not necessarily expected behaviour. BUT, it is in line with
    // previous use (that if nT==1, then assume ssd constant in T).
    Index nT = T_DATAGRID.nelem();
    if (nT > 1) {
      // Check, that we not have data that has already been T-reduced (in
      // pha_mat only. complete ssd T-reduce should have been sorted away
      // already above).
      if (PHA_MAT_DATA.nvitrines() != nT) {
        std::ostringstream os;
        os << "Single scattering data of scat element #" << i_se
           << " of scat species #" << i_ss << "\n"
           << "seems to have undergone some temperature grid manipulation in\n"
           << "*pha_mat_data* already. That can not be done twice!";
        throw std::runtime_error(os.str());
      }

      // Check that ext_mat and abs_vec have the same temp dimensions as T_grid.
      // This should always be true, if not it's a bug not a user mistake, hence
      // use ARTS_ASSERT.
      ARTS_ASSERT(EXT_MAT_DATA.nbooks() == nT and ABS_VEC_DATA.nbooks() == nT);

      // Check that T_grid is consistent with requested interpolation order
      std::ostringstream ost;
      ost << "Scattering data temperature interpolation for\n"
          << "scat element #" << i_se << " of scat species #" << i_ss << ".";
      chk_interpolation_grids(ost.str(), T_DATAGRID, T, interp_order);

      // Gridpositions:
      const LagrangeInterpolation lag_T(0, T, T_DATAGRID, interp_order);
      const auto itw = interpweights(lag_T);

      //Sizing of temporary SSD data containers
      Tensor7 phamat_tmp(PHA_MAT_DATA.nlibraries(),
                         1,
                         PHA_MAT_DATA.nshelves(),
                         PHA_MAT_DATA.nbooks(),
                         PHA_MAT_DATA.npages(),
                         PHA_MAT_DATA.nrows(),
                         PHA_MAT_DATA.ncols(),
                         0.);
      Tensor5 extmat_tmp(EXT_MAT_DATA.nshelves(),
                         1,
                         EXT_MAT_DATA.npages(),
                         EXT_MAT_DATA.nrows(),
                         EXT_MAT_DATA.ncols(),
                         0.);
      Tensor5 absvec_tmp(ABS_VEC_DATA.nshelves(),
                         1,
                         ABS_VEC_DATA.npages(),
                         ABS_VEC_DATA.nrows(),
                         ABS_VEC_DATA.ncols(),
                         0.);

      // a1) temp interpol of pha mat
      //We have to apply the interpolation separately for each of the pha_mat
      //entries, i.e. loop over all remaining size dimensions
      //We don't apply any transformation here, but interpolate the actual
      //stored ssd (i.e. not the 4x4matrices, but the 7-16 elements separately).
      for (Index i_f = 0; i_f < PHA_MAT_DATA.nlibraries(); i_f++)
        for (Index i_za1 = 0; i_za1 < PHA_MAT_DATA.nshelves(); i_za1++)
          for (Index i_aa1 = 0; i_aa1 < PHA_MAT_DATA.nbooks(); i_aa1++)
            for (Index i_za2 = 0; i_za2 < PHA_MAT_DATA.npages(); i_za2++)
              for (Index i_aa2 = 0; i_aa2 < PHA_MAT_DATA.nrows(); i_aa2++)
                for (Index i_st = 0; i_st < PHA_MAT_DATA.ncols(); i_st++)
                  phamat_tmp(i_f, 0, i_za1, i_aa1, i_za2, i_aa2, i_st) =
                      interp(PHA_MAT_DATA(
                                 i_f, joker, i_za1, i_aa1, i_za2, i_aa2, i_st),
                             itw,
                             lag_T);

      // a2) temp interpol of ext and abs.
      //We do that regardless of whether they should be reduced or not, because
      //we need them also for norm checking / renorming.
      for (Index i_f = 0; i_f < EXT_MAT_DATA.nshelves(); i_f++)
        for (Index i_za = 0; i_za < EXT_MAT_DATA.npages(); i_za++)
          for (Index i_aa = 0; i_aa < EXT_MAT_DATA.nrows(); i_aa++) {
            for (Index i_st = 0; i_st < EXT_MAT_DATA.ncols(); i_st++)
              extmat_tmp(i_f, 0, i_za, i_aa, i_st) =
                  interp(EXT_MAT_DATA(i_f, joker, i_za, i_aa, i_st), itw, lag_T);
            for (Index i_st = 0; i_st < ABS_VEC_DATA.ncols(); i_st++)
              absvec_tmp(i_f, 0, i_za, i_aa, i_st) =
                  interp(ABS_VEC_DATA(i_f, joker, i_za, i_aa, i_st), itw, lag_T);
          }

      // Norm & other consistency checks.
      // All done separately for PTYPE_TOTAL_RND and PTYPE_AZIMUTH_RND (the
      // latter needs to loop over scat_za_inc).
      //
      // b) calculate norm of T-reduced pha mat
      // c) check pha mat norm vs. sca xs from ext-abs at T_interpol
      // d) Ensure that T-reduced data is consistent/representative of all data.
      //    and throw error/disallow reduction if sca xs varying too much.
      // d1) in case of pha_mat only reduction, the scat xs (pha_mat norm) needs
      // to be consistent with the sca xs from over the ext/abs T_grid. This is
      // essentially an energy conservation issue. That is, we should be as
      // strict here as with pha_mat norm deviations in general (however, should
      // we allow the norm at T_grid to deviate by a threshold from the
      // T_interpol norm (that should be a little looser) or from the ext-abs
      // derived expected norm?).
      // d2) in case of all-ssp reduction, the data should still be
      // representative. b)&c) ensure data consistency in itself, making this
      // rather an error on the SSP as such. Hence, we might be a little more
      // loose here.
      // d) the resulting check for d1) and d2) is the same (ext-abs sca xs at
      // T_interpol vs ext-abs sca xs at T_grid), but we use different
      // thresholds.
      //
      // FIXME?
      // Regarding b)&c) should we also calc norm of original-T pha mats? To get
      // a measure how strong they already deviate from expected norm (As we can
      // not be more exact here than what is already in the original data...).
      // On the other hand, a certain accuracy should be guaranteed from
      // scat_dataCheck already.
      // Hence, for now we skip that (but maybe added later when proves
      // necessary).
      //
      // FIXME?
      // Regarding d1), we could alternatively make sure here that norm at
      // T_interpol is good. And later on ignore any deviations between norm and
      // ext-abs sca xs and instead blindly renorm to expected norm (would that
      // be ok? correct norm here, doesn't imply correct norm at whatever scat
      // angle grids the user is applying. for that, we could in place also calc
      // the original-data norm. but that might be expensive (as we can't do
      // that from ext-abs sca xs, because we don't know to which T that refers.
      // that would go away if we'd actually store pha_mat normed to 1 or 4Pi.
      // but that's prob not going to happen. is it? Another option would be to
      // introduce an additional T_grid, eg T_grid_phamat.). which we actually
      // want to avoid :-/

      Numeric this_threshold;
      String errmsg;
      if (phamat_only) {
        this_threshold = threshold;
        errmsg =
            "T-reduced *pha_mat_data* norm (=sca xs) deviates too "
            "much from non-reduced *ext_mat_data* and *abs_vec_data*:";
      } else {
        this_threshold = 2 * threshold;
        errmsg =
            "T-reduced *scat_data* deviates too much from original "
            "*scat_data*:";
      }

      // The norm-check code is copied and slightly adapted from scat_dataCheck.
      // Might be better to make a functon out of this and use in both places
      // for consistency.
      //
      // FIXME: no checks on higher Stokes elements are done. Should there?
      // Which?
      switch (PART_TYPE) {
        case PTYPE_TOTAL_RND: {
          for (Index f = 0; f < F_DATAGRID.nelem(); f++) {
            // b) calculate norm of T-reduced pha mat
            Numeric Csca = AngIntegrate_trapezoid(
                phamat_tmp(f, 0, joker, 0, 0, 0, 0), ZA_DATAGRID);
            Numeric Cext_data = extmat_tmp(f, 0, 0, 0, 0);
            //Numeric Cabs = Cext_data - Csca;
            Numeric Cabs_data = absvec_tmp(f, 0, 0, 0, 0);
            Numeric Csca_data = Cext_data - Cabs_data;

            /*
            cout << "  Coefficients in data: "
                 << "Cext: " << Cext_data << " Cabs: " << Cabs_data
                 << " Csca: " << Csca_data << "\n"
                 << "  Calculated coefficients: "
                 << "Cabs calc: " << Cabs
                 << " Csca calc: " << Csca << "\n"
                 << "  Deviations "
                 << "Cabs: " << 1e2*Cabs/Cabs_data-1e2
                 << "% Csca: " << 1e2*Csca/Csca_data-1e2
                 << "% Alb: " << (Csca-Csca_data)/Cext_data << "\n";
            */

            // c) check pha mat norm vs. sca xs from ext-abs at T_interpol (as
            // albedo dev check)
            if (abs(Csca - Csca_data) / Cext_data > threshold) {
              std::ostringstream os;
              os << "  Deviations in T-reduced scat_data too large:\n"
                 << "  scat dev [%] " << 1e2 * Csca / Csca_data - 1e2
                 << " at nominal (actual) albedo of " << Csca_data / Cext_data
                 << " (" << Csca / Cext_data << ").\n"
                 << "  Problem occurs for scattering element #" << i_se
                 << " at " << f << ".frequency!\n";
              throw std::runtime_error(os.str());
            }
            Numeric norm_dev = (Csca - Csca) / Cext_data;

            // d) Ensure that T-reduced data is consistent/representative of all data.
            // below use theoretical (ext-abs derived) sca xs as reference.
            Csca = Csca_data;
            for (Index t = 0; t < T_DATAGRID.nelem(); t++) {
              Cext_data = EXT_MAT_DATA(f, t, 0, 0, 0);
              Csca_data = Cext_data - ABS_VEC_DATA(f, t, 0, 0, 0);
              Numeric xs_dev = (Csca - Csca_data) / Cext_data;
              if (abs(norm_dev + (Csca - Csca_data) / Cext_data) >
                  this_threshold)
                std::cout << "Accumulated deviation (abs(" << norm_dev << "+"
                     << xs_dev << ")=" << abs(norm_dev + xs_dev)
                     << " exceeding threshold (" << this_threshold << ").\n";
              if (abs(Csca - Csca_data) / Cext_data > this_threshold) {
                std::ostringstream os;
                os << "  " << errmsg << "\n"
                   << "  scat dev [%] " << 1e2 * Csca / Csca_data - 1e2
                   << " at nominal (actual) albedo of " << Csca_data / Cext_data
                   << " (" << Csca / Cext_data << ").\n"
                   << "  Problem occurs for scattering element #" << i_se
                   << " at " << f << ".frequency and " << t
                   << ".temperature!\n";
                throw std::runtime_error(os.str());
              }
            }
          }
          break;
        }

        case PTYPE_AZIMUTH_RND: {
          for (Index f = 0; f < F_DATAGRID.nelem(); f++) {
            for (Index iza = 0; iza < ABS_VEC_DATA.npages(); iza++) {
              // b) calculate norm of T-reduced pha mat
              Numeric Csca = 2 * AngIntegrate_trapezoid(
                                     phamat_tmp(f, 0, joker, joker, iza, 0, 0),
                                     ZA_DATAGRID,
                                     AA_DATAGRID);
              Numeric Cext_data = extmat_tmp(f, 0, iza, 0, 0);
              //Numeric Cabs = Cext_data - Csca;
              Numeric Cabs_data = absvec_tmp(f, 0, iza, 0, 0);
              Numeric Csca_data = Cext_data - Cabs_data;

              /*
              cout << "  Coefficients in data: "
                   << "Cext: " << Cext_data << " Cabs: " << Cabs_data
                   << " Csca: " << Csca_data << "\n"
                   << "  Calculated coefficients: "
                   << "Cabs calc: " << Cabs
                   << " Csca calc: " << Csca << "\n"
                   << "  Deviations "
                   << "Cabs: " << 1e2*Cabs/Cabs_data-1e2
                   << "% Csca: " << 1e2*Csca/Csca_data-1e2
                   << "% Alb: " << (Csca-Csca_data)/Cext_data << "\n";
              */

              // c) check pha mat norm vs. sca xs from ext-abs at T_interpol (as
              // albedo dev check)
              if (abs(Csca - Csca_data) / Cext_data > threshold) {
                std::ostringstream os;
                os << "  Deviations in T-reduced scat_data too large:\n"
                   << "  scat dev [%] " << 1e2 * Csca / Csca_data - 1e2
                   << " at nominal (actual) albedo of " << Csca_data / Cext_data
                   << " (" << Csca / Cext_data << ").\n"
                   << "  Problem occurs for scattering element #" << i_se
                   << " at " << f << ".frequency, and " << iza
                   << ". incident polar angle!\n";
                throw std::runtime_error(os.str());
              }

              // d) Ensure that T-reduced data is consistent/representative of all data.
              // below use theoretical (ext-abs derived) sca xs as reference.
              Csca = Csca_data;
              for (Index t = 0; t < T_DATAGRID.nelem(); t++) {
                Cext_data = EXT_MAT_DATA(f, t, 0, 0, 0);
                Csca_data = Cext_data - ABS_VEC_DATA(f, t, 0, 0, 0);
                if (abs(Csca - Csca_data) / Cext_data > this_threshold) {
                  std::ostringstream os;
                  os << "  " << errmsg << "\n"
                     << "  scat dev [%] " << 1e2 * Csca / Csca_data - 1e2
                     << " at nominal (actual) albedo of "
                     << Csca_data / Cext_data << " (" << Csca / Cext_data
                     << ").\n"
                     << "  Problem occurs for scattering element #" << i_se
                     << " at " << f << ".frequency and " << t
                     << ".temperature, and " << iza
                     << ". incident polar angle!\n";
                  throw std::runtime_error(os.str());
                }
              }
            }
          }
          break;
        }

        default: {
          // other ptype cases already excluded above. i.e. we shouldn't end up
          // here. If we do, that's a bug.
          ARTS_ASSERT(0);
        }
      }

      PHA_MAT_DATA = phamat_tmp;
      //We don't need to reset the scat element's grids!
      //Except for T_grid in the case that we reduce ALL three ssd variables.
      if (!phamat_only) {
        T_DATAGRID.resize(1);
        T_DATAGRID = T;
        EXT_MAT_DATA = extmat_tmp;
        ABS_VEC_DATA = absvec_tmp;
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void scat_data_monoCalc(ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                        const ArrayOfArrayOfSingleScatteringData& scat_data,
                        const Vector& f_grid,
                        const Index& f_index) {
  //Extrapolation factor:
  //const Numeric extpolfac = 0.5;

  // Check, whether single scattering data contains the right frequencies:
  // The check was changed to allow extrapolation at the boundaries of the
  // frequency grid.
  for (Index h = 0; h < scat_data.nelem(); h++) {
    for (Index i = 0; i < scat_data[h].nelem(); i++) {
      // check with extrapolation
      chk_interpolation_grids("scat_data.f_grid to f_grid",
                              scat_data[h][i].f_grid,
                              f_grid[f_index]);
    }
  }

  //Initialise scat_data_mono
  scat_data_mono.resize(scat_data.nelem());

  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++) {
    const Index N_se = scat_data[i_ss].nelem();

    //Initialise scat_data_mono
    scat_data_mono[i_ss].resize(N_se);

    // Loop over the included scattering elements
    for (Index i_se = 0; i_se < N_se; i_se++) {
      // Gridpositions:
      GridPos freq_gp;
      gridpos(freq_gp, F_DATAGRID, f_grid[f_index]);

      // Interpolation weights:
      Vector itw(2);
      interpweights(itw, freq_gp);

      //Stuff that doesn't need interpolating
      scat_data_mono[i_ss][i_se].ptype = PART_TYPE;
      scat_data_mono[i_ss][i_se].f_grid.resize(1);
      scat_data_mono[i_ss][i_se].f_grid = f_grid[f_index];
      scat_data_mono[i_ss][i_se].T_grid = scat_data[i_ss][i_se].T_grid;
      scat_data_mono[i_ss][i_se].za_grid = ZA_DATAGRID;
      scat_data_mono[i_ss][i_se].aa_grid = AA_DATAGRID;

      //Phase matrix data
      scat_data_mono[i_ss][i_se].pha_mat_data.resize(1,
                                                     PHA_MAT_DATA.nvitrines(),
                                                     PHA_MAT_DATA.nshelves(),
                                                     PHA_MAT_DATA.nbooks(),
                                                     PHA_MAT_DATA.npages(),
                                                     PHA_MAT_DATA.nrows(),
                                                     PHA_MAT_DATA.ncols());

      for (Index t_index = 0; t_index < PHA_MAT_DATA.nvitrines(); t_index++) {
        for (Index i_za_sca = 0; i_za_sca < PHA_MAT_DATA.nshelves();
             i_za_sca++) {
          for (Index i_aa_sca = 0; i_aa_sca < PHA_MAT_DATA.nbooks();
               i_aa_sca++) {
            for (Index i_za_inc = 0; i_za_inc < PHA_MAT_DATA.npages();
                 i_za_inc++) {
              for (Index i_aa_inc = 0; i_aa_inc < PHA_MAT_DATA.nrows();
                   i_aa_inc++) {
                for (Index i = 0; i < PHA_MAT_DATA.ncols(); i++) {
                  scat_data_mono[i_ss][i_se].pha_mat_data(
                      0, t_index, i_za_sca, i_aa_sca, i_za_inc, i_aa_inc, i) =
                      interp(itw,
                             PHA_MAT_DATA(joker,
                                          t_index,
                                          i_za_sca,
                                          i_aa_sca,
                                          i_za_inc,
                                          i_aa_inc,
                                          i),
                             freq_gp);
                }
              }
            }
          }
        }
        //Extinction matrix data
        scat_data_mono[i_ss][i_se].ext_mat_data.resize(1,
                                                       T_DATAGRID.nelem(),
                                                       EXT_MAT_DATA.npages(),
                                                       EXT_MAT_DATA.nrows(),
                                                       EXT_MAT_DATA.ncols());
        for (Index i_za_sca = 0; i_za_sca < EXT_MAT_DATA.npages(); i_za_sca++) {
          for (Index i_aa_sca = 0; i_aa_sca < EXT_MAT_DATA.nrows();
               i_aa_sca++) {
            //
            // Interpolation of extinction matrix:
            //
            for (Index i = 0; i < EXT_MAT_DATA.ncols(); i++) {
              scat_data_mono[i_ss][i_se].ext_mat_data(
                  0, t_index, i_za_sca, i_aa_sca, i) =
                  interp(itw,
                         EXT_MAT_DATA(joker, t_index, i_za_sca, i_aa_sca, i),
                         freq_gp);
            }
          }
        }
        //Absorption vector data
        scat_data_mono[i_ss][i_se].abs_vec_data.resize(1,
                                                       T_DATAGRID.nelem(),
                                                       ABS_VEC_DATA.npages(),
                                                       ABS_VEC_DATA.nrows(),
                                                       ABS_VEC_DATA.ncols());
        for (Index i_za_sca = 0; i_za_sca < ABS_VEC_DATA.npages(); i_za_sca++) {
          for (Index i_aa_sca = 0; i_aa_sca < ABS_VEC_DATA.nrows();
               i_aa_sca++) {
            //
            // Interpolation of absorption vector:
            //
            for (Index i = 0; i < ABS_VEC_DATA.ncols(); i++) {
              scat_data_mono[i_ss][i_se].abs_vec_data(
                  0, t_index, i_za_sca, i_aa_sca, i) =
                  interp(itw,
                         ABS_VEC_DATA(joker, t_index, i_za_sca, i_aa_sca, i),
                         freq_gp);
            }
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void scat_data_monoExtract(ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                           const ArrayOfArrayOfSingleScatteringData& scat_data,
                           const Index& f_index) {
  //Initialise scat_data_mono
  scat_data_mono.resize(scat_data.nelem());

  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++) {
    const Index N_se = scat_data[i_ss].nelem();

    //Initialise scat_data_mono
    scat_data_mono[i_ss].resize(N_se);

    // Loop over the included scattering elements
    for (Index i_se = 0; i_se < N_se; i_se++) {
      Index nf = F_DATAGRID.nelem();
      if (nf == 1) {
        scat_data_mono[i_ss][i_se] = scat_data[i_ss][i_se];
      } else {
        //Stuff that doesn't need interpolating
        scat_data_mono[i_ss][i_se].ptype = PART_TYPE;
        scat_data_mono[i_ss][i_se].T_grid = T_DATAGRID;
        scat_data_mono[i_ss][i_se].za_grid = ZA_DATAGRID;
        scat_data_mono[i_ss][i_se].aa_grid = AA_DATAGRID;

        scat_data_mono[i_ss][i_se].f_grid.resize(1);
        scat_data_mono[i_ss][i_se].f_grid = F_DATAGRID[f_index];

        Index this_f_index;

        //Phase matrix data
        /*scat_data_mono[i_ss][i_se].pha_mat_data.resize(1,
                                                 PHA_MAT_DATA.nvitrines(),
                                                 PHA_MAT_DATA.nshelves(),
                                                 PHA_MAT_DATA.nbooks(),
                                                 PHA_MAT_DATA.npages(),
                                                 PHA_MAT_DATA.nrows(),
                                                 PHA_MAT_DATA.ncols());*/
        this_f_index = (PHA_MAT_DATA.nlibraries() == 1) ? 0 : f_index;
        scat_data_mono[i_ss][i_se].pha_mat_data = PHA_MAT_DATA(
            Range(this_f_index, 1), joker, joker, joker, joker, joker, joker);

        //Extinction matrix data
        /*scat_data_mono[i_ss][i_se].ext_mat_data.resize(1, T_DATAGRID.nelem(),
                                                     EXT_MAT_DATA.npages(),
                                                     EXT_MAT_DATA.nrows(),
                                                     EXT_MAT_DATA.ncols());*/
        this_f_index = (EXT_MAT_DATA.nshelves() == 1) ? 0 : f_index;
        scat_data_mono[i_ss][i_se].ext_mat_data =
            EXT_MAT_DATA(Range(this_f_index, 1), joker, joker, joker, joker);

        //Absorption vector data
        /*scat_data_mono[i_ss][i_se].abs_vec_data.resize(1, T_DATAGRID.nelem(),
                                                     ABS_VEC_DATA.npages(),
                                                     ABS_VEC_DATA.nrows(),
                                                     ABS_VEC_DATA.ncols());*/
        this_f_index = (ABS_VEC_DATA.nshelves() == 1) ? 0 : f_index;
        scat_data_mono[i_ss][i_se].abs_vec_data =
            ABS_VEC_DATA(Range(this_f_index, 1), joker, joker, joker, joker);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void opt_prop_sptFromMonoData(  // Output and Input:
    ArrayOfPropmatVector& ext_mat_spt,
    ArrayOfStokvecVector& abs_vec_spt,
    // Input:
    const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
    const Vector& za_grid,
    const Vector& aa_grid,
    const Index& za_index,  // propagation directions
    const Index& aa_index,
    const Numeric& rtp_temperature,
    const Tensor4& pnd_field,
    const Index& scat_p_index,
    const Index& scat_lat_index,
    const Index& scat_lon_index) {
  DEBUG_ONLY(const Index N_se_total = TotalNumberOfElements(scat_data_mono);)
  const Numeric za_sca = za_grid[za_index];
  const Numeric aa_sca = aa_grid[aa_index];

  ARTS_ASSERT(ext_mat_spt.nelem() == N_se_total);
  ARTS_ASSERT(abs_vec_spt.nelem() == N_se_total);

  // Check that we do indeed have scat_data_mono here. Only checking the first
  // scat element, assuming the other elements have been processed in the same
  // manner. That's save against having scat_data here if that originated from
  // scat_data_raw reading routines (ScatSpecies/Element*Add/Read), it's not safe
  // against data read by ReadXML directly or if scat_data(_raw) has been (partly)
  // produced from scat_data_singleTmatrix. That would be too costly here,
  // though.
  // Also, we can't check here whether data is at the correct frequency since we
  // don't know f_grid and f_index here (we could pass it in, though).
  if (scat_data_mono[0][0].f_grid.nelem() > 1) {
    std::ostringstream os;
    os << "Scattering data seems to be *scat_data* (several freq points),\n"
       << "but *scat_data_mono* (1 freq point only) is expected here.";
    throw std::runtime_error(os.str());
  }

  // Initialisation
  for (auto& pm : ext_mat_spt) pm = 0.;
  for (auto& av : abs_vec_spt) av = 0.;

  GridPos t_gp;

  Vector itw(2);

  Index i_se_flat = 0;
  // Loop over the included scattering elements
  for (Index i_ss = 0; i_ss < scat_data_mono.nelem(); i_ss++) {
    // Loop over the included scattering elements
    for (Index i_se = 0; i_se < scat_data_mono[i_ss].nelem(); i_se++) {
      // If the particle number density at a specific point in the
      // atmosphere for the i_se scattering element is zero, we don't need
      // to do the transformation!
      if (pnd_field(i_se_flat, scat_p_index, scat_lat_index, scat_lon_index) >
          PND_LIMIT) {
        // First we have to transform the data from the coordinate system
        // used in the database (depending on the kind of ptype) to the
        // laboratory coordinate system.

        //
        // Do the transformation into the laboratory coordinate system.
        //
        // Extinction matrix:
        //
        Index ext_npages = scat_data_mono[i_ss][i_se].ext_mat_data.npages();
        Index ext_nrows = scat_data_mono[i_ss][i_se].ext_mat_data.nrows();
        Index ext_ncols = scat_data_mono[i_ss][i_se].ext_mat_data.ncols();
        Index abs_npages = scat_data_mono[i_ss][i_se].abs_vec_data.npages();
        Index abs_nrows = scat_data_mono[i_ss][i_se].abs_vec_data.nrows();
        Index abs_ncols = scat_data_mono[i_ss][i_se].abs_vec_data.ncols();

        //Check that scattering data temperature range covers required temperature
        ConstVectorView t_grid = scat_data_mono[i_ss][i_se].T_grid;

        if (t_grid.nelem() > 1) {
          std::ostringstream os;
          os << "In opt_prop_sptFromMonoData.\n"
             << "The temperature grid of the scattering data does not\n"
             << "cover the atmospheric temperature at cloud location.\n"
             << "The data should include the value T = " << rtp_temperature
             << " K.";
          chk_interpolation_grids(os.str(), t_grid, rtp_temperature);

          //interpolate over temperature
          Tensor3 ext_mat_data1temp(ext_npages, ext_nrows, ext_ncols);
          gridpos(t_gp, t_grid, rtp_temperature);
          interpweights(itw, t_gp);
          for (Index i_p = 0; i_p < ext_npages; i_p++) {
            for (Index i_r = 0; i_r < ext_nrows; i_r++) {
              for (Index i_c = 0; i_c < ext_ncols; i_c++) {
                ext_mat_data1temp(i_p, i_r, i_c) =
                    interp(itw,
                           scat_data_mono[i_ss][i_se].ext_mat_data(
                               0, joker, i_p, i_r, i_c),
                           t_gp);
              }
            }
          }
          ext_matTransform(ext_mat_spt[i_se_flat],
                           ext_mat_data1temp,
                           scat_data_mono[i_ss][i_se].za_grid,
                           scat_data_mono[i_ss][i_se].aa_grid,
                           scat_data_mono[i_ss][i_se].ptype,
                           za_sca,
                           aa_sca);
        } else {
          ext_matTransform(ext_mat_spt[i_se_flat],
                           scat_data_mono[i_ss][i_se].ext_mat_data(
                               0, 0, joker, joker, joker),
                           scat_data_mono[i_ss][i_se].za_grid,
                           scat_data_mono[i_ss][i_se].aa_grid,
                           scat_data_mono[i_ss][i_se].ptype,
                           za_sca,
                           aa_sca);
        }
        //
        // Absorption vector:
        //

        if (t_grid.nelem() > 1) {
          Tensor3 abs_vec_data1temp(abs_npages, abs_nrows, abs_ncols);
          //interpolate over temperature
          for (Index i_p = 0; i_p < abs_npages; i_p++) {
            for (Index i_r = 0; i_r < abs_nrows; i_r++) {
              for (Index i_c = 0; i_c < abs_ncols; i_c++) {
                abs_vec_data1temp(i_p, i_r, i_c) =
                    interp(itw,
                           scat_data_mono[i_ss][i_se].abs_vec_data(
                               0, joker, i_p, i_r, i_c),
                           t_gp);
              }
            }
          }
          abs_vecTransform(abs_vec_spt[i_se_flat],
                           abs_vec_data1temp,
                           scat_data_mono[i_ss][i_se].za_grid,
                           scat_data_mono[i_ss][i_se].aa_grid,
                           scat_data_mono[i_ss][i_se].ptype,
                           za_sca,
                           aa_sca);
        } else {
          abs_vecTransform(abs_vec_spt[i_se_flat],
                           scat_data_mono[i_ss][i_se].abs_vec_data(
                               0, 0, joker, joker, joker),
                           scat_data_mono[i_ss][i_se].za_grid,
                           scat_data_mono[i_ss][i_se].aa_grid,
                           scat_data_mono[i_ss][i_se].ptype,
                           za_sca,
                           aa_sca);
        }
      }

      i_se_flat++;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void pha_mat_sptFromMonoData(  // Output:
    Tensor5& pha_mat_spt,
    // Input:
    const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
    const Index& doit_za_grid_size,
    const Vector& aa_grid,
    const Index& za_index,  // propagation directions
    const Index& aa_index,
    const Numeric& rtp_temperature,
    const Tensor4& pnd_field,
    const Index& scat_p_index,
    const Index& scat_lat_index,
    const Index& scat_lon_index) {
  Vector za_grid;
  nlinspace(za_grid, 0, 180, doit_za_grid_size);

  const Index N_se_total = TotalNumberOfElements(scat_data_mono);
  if (N_se_total != pnd_field.nbooks()) {
    std::ostringstream os;
    os << "Total number of scattering elements in *scat_data_mono* "
       << "inconsistent with size of pnd_field.";
    throw std::runtime_error(os.str());
  }
  // as pha_mat_spt is typically initialized from pnd_field, this theoretically
  // checks the same as the std::runtime_error above. Still, we keep it to be on the
  // save side.
  ARTS_ASSERT(pha_mat_spt.nshelves() == N_se_total);

  // Check that we do indeed have scat_data_mono here. Only checking the first
  // scat element, assuming the other elements have been processed in the same
  // manner. That's save against having scat_data here if that originated from
  // scat_data_raw reading routines (ScatSpecies/Element*Add/Read), it's not safe
  // against data read by ReadXML directly or if scat_data(_raw) has been (partly)
  // produced from scat_data_singleTmatrix. That would be too costly here,
  // though.
  // Also, we can't check here whether data is at the correct frequency since we
  // don't know f_grid and f_index here (we could pass it in, though).
  if (scat_data_mono[0][0].f_grid.nelem() > 1) {
    std::ostringstream os;
    os << "Scattering data seems to be *scat_data* (several freq points),\n"
       << "but *scat_data_mono* (1 freq point only) is expected here.";
    throw std::runtime_error(os.str());
  }

  GridPos T_gp = {0, {0, 1}}, Tred_gp;
  Vector itw(2);

  // Initialisation
  pha_mat_spt = 0.;

  Index i_se_flat = 0;
  for (Index i_ss = 0; i_ss < scat_data_mono.nelem(); i_ss++) {
    for (Index i_se = 0; i_se < scat_data_mono[i_ss].nelem(); i_se++) {
      // If the particle number density at a specific point in the
      // atmosphere for scattering element i_se is zero, we don't need to
      // do the transformation!
      if (pnd_field(i_se_flat, scat_p_index, scat_lat_index, scat_lon_index) >
          PND_LIMIT) {
        // Temporary phase matrix which includes all temperatures.
        Index nT = scat_data_mono[i_ss][i_se].pha_mat_data.nvitrines();
        Tensor3 pha_mat_spt_tmp(nT, pha_mat_spt.nrows(), pha_mat_spt.ncols());

        pha_mat_spt_tmp = 0.;

        Index ti = -1;
        if (nT == 1)  // just 1 T_grid element
        {
          ti = 0;
        } else if (rtp_temperature < 0.)  // coding for 'not interpolate, but
                                          // pick one temperature'
        {
          if (rtp_temperature > -10.)  // lowest T-point
          {
            ti = 0;
          } else if (rtp_temperature > -20.)  // highest T-point
          {
            ti = nT - 1;
          } else  // median T-point
          {
            ti = nT / 2;
          }
        } else {
          std::ostringstream os;
          os << "In pha_mat_sptFromMonoData.\n"
             << "The temperature grid of the scattering data does not\n"
             << "cover the atmospheric temperature at cloud location.\n"
             << "The data should include the value T = " << rtp_temperature
             << " K.";
          chk_interpolation_grids(
              os.str(), scat_data_mono[i_ss][i_se].T_grid, rtp_temperature);

          // Gridpositions:
          gridpos(T_gp, scat_data_mono[i_ss][i_se].T_grid, rtp_temperature);
          gridpos_copy(Tred_gp, T_gp);
          Tred_gp.idx = 0;
          // Interpolation weights:
          interpweights(itw, Tred_gp);
        }

        // Do the transformation into the laboratory coordinate system.
        for (Index za_inc_idx = 0; za_inc_idx < doit_za_grid_size;
             za_inc_idx++) {
          for (Index aa_inc_idx = 0; aa_inc_idx < aa_grid.nelem();
               aa_inc_idx++) {
            if (ti < 0)  // Temperature interpolation
            {
              for (Index t_idx = 0; t_idx < 2; t_idx++) {
                pha_matTransform(
                    pha_mat_spt_tmp(t_idx, joker, joker),
                    scat_data_mono[i_ss][i_se].pha_mat_data(
                        0, t_idx + T_gp.idx, joker, joker, joker, joker, joker),
                    scat_data_mono[i_ss][i_se].za_grid,
                    scat_data_mono[i_ss][i_se].aa_grid,
                    scat_data_mono[i_ss][i_se].ptype,
                    za_index,
                    aa_index,
                    za_inc_idx,
                    aa_inc_idx,
                    za_grid,
                    aa_grid);
              }

              for (Index i = 0; i < 4; i++) {
                for (Index j = 0; j < 4; j++) {
                  pha_mat_spt(i_se_flat, za_inc_idx, aa_inc_idx, i, j) =
                      interp(itw, pha_mat_spt_tmp(joker, i, j), Tred_gp);
                }
              }
            } else  // no temperature interpolation required
            {
              pha_matTransform(
                  pha_mat_spt(i_se_flat, za_inc_idx, aa_inc_idx, joker, joker),
                  scat_data_mono[i_ss][i_se].pha_mat_data(
                      0, ti, joker, joker, joker, joker, joker),
                  scat_data_mono[i_ss][i_se].za_grid,
                  scat_data_mono[i_ss][i_se].aa_grid,
                  scat_data_mono[i_ss][i_se].ptype,
                  za_index,
                  aa_index,
                  za_inc_idx,
                  aa_inc_idx,
                  za_grid,
                  aa_grid);
            }
          }
        }
      }

      i_se_flat++;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void pha_mat_sptFromScat_data(  // Output:
    Tensor5& pha_mat_spt,
    // Input:
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Index& scat_data_checked,
    const Vector& za_grid,
    const Vector& aa_grid,
    const Index& za_index,  // propagation directions
    const Index& aa_index,
    const Index& f_index,
    const Numeric& rtp_temperature,
    const Tensor4& pnd_field,
    const Index& scat_p_index,
    const Index& scat_lat_index,
    const Index& scat_lon_index) {
  if (scat_data_checked != 1)
    throw std::runtime_error(
        "The scattering data must be flagged to have "
        "passed a consistency check (scat_data_checked=1).");

  // Determine total number of scattering elements
  const Index N_se_total = TotalNumberOfElements(scat_data);
  if (N_se_total != pnd_field.nbooks()) {
    std::ostringstream os;
    os << "Total number of scattering elements in scat_data "
       << "inconsistent with size of pnd_field.";
    throw std::runtime_error(os.str());
  }
  // as pha_mat_spt is typically initialized from pnd_field, this theoretically
  // checks the same as the std::runtime_error above. Still, we keep it to be on the
  // save side.
  ARTS_ASSERT(pha_mat_spt.nshelves() == N_se_total);

  const Index N_ss = scat_data.nelem();

  // Phase matrix in laboratory coordinate system. Dimensions:
  // [frequency, za_inc, aa_inc, 4, 4]
  Tensor5 pha_mat_data_int;

  Index this_f_index;

  Index i_se_flat = 0;
  // Loop over scattering species
  for (Index i_ss = 0; i_ss < N_ss; i_ss++) {
    const Index N_se = scat_data[i_ss].nelem();

    // Loop over the included scattering elements
    for (Index i_se = 0; i_se < N_se; i_se++) {
      // If the particle number density at a specific point in the
      // atmosphere for the i_se scattering element is zero, we don't need
      // to do the transfromation!
      if (abs(pnd_field(
              i_se_flat, scat_p_index, scat_lat_index, scat_lon_index)) >
          PND_LIMIT) {
        // First we have to transform the data from the coordinate system
        // used in the database (depending on the kind of ptype) to the
        // laboratory coordinate system.

        // Resize the variables for the interpolated data (1freq, 1T):
        pha_mat_data_int.resize(PHA_MAT_DATA.nshelves(),
                                PHA_MAT_DATA.nbooks(),
                                PHA_MAT_DATA.npages(),
                                PHA_MAT_DATA.nrows(),
                                PHA_MAT_DATA.ncols());

        // Frequency extraction and temperature interpolation

        // Gridpositions and interpolation weights;
        GridPos t_gp;
        Vector itw;
        Index this_T_index = -1;
        if (PHA_MAT_DATA.nvitrines() == 1) {
          this_T_index = 0;
        } else if (rtp_temperature < 0.)  // coding for 'not interpolate, but
                                          // pick one temperature'
        {
          if (rtp_temperature > -10.)  // lowest T-point
          {
            this_T_index = 0;
          } else if (rtp_temperature > -20.)  // highest T-point
          {
            this_T_index = PHA_MAT_DATA.nvitrines() - 1;
          } else  // median T-point
          {
            this_T_index = PHA_MAT_DATA.nvitrines() / 2;
          }
        } else {
          std::ostringstream os;
          os << "In pha_mat_sptFromScat_data.\n"
             << "The temperature grid of the scattering data does not\n"
             << "cover the atmospheric temperature at cloud location.\n"
             << "The data should include the value T = " << rtp_temperature
             << " K.";
          chk_interpolation_grids(os.str(), T_DATAGRID, rtp_temperature);

          gridpos(t_gp, T_DATAGRID, rtp_temperature);

          // Interpolation weights:
          itw.resize(2);
          interpweights(itw, t_gp);
        }

        if (PHA_MAT_DATA.nlibraries() == 1)
          this_f_index = 0;
        else
          this_f_index = f_index;

        if (this_T_index < 0) {
          // Interpolation of scattering matrix:
          for (Index i_za_sca = 0; i_za_sca < PHA_MAT_DATA.nshelves();
               i_za_sca++)
            for (Index i_aa_sca = 0; i_aa_sca < PHA_MAT_DATA.nbooks();
                 i_aa_sca++)
              for (Index i_za_inc = 0; i_za_inc < PHA_MAT_DATA.npages();
                   i_za_inc++)
                for (Index i_aa_inc = 0; i_aa_inc < PHA_MAT_DATA.nrows();
                     i_aa_inc++)
                  for (Index i = 0; i < PHA_MAT_DATA.ncols(); i++)
                    pha_mat_data_int(
                        i_za_sca, i_aa_sca, i_za_inc, i_aa_inc, i) =
                        interp(itw,
                               PHA_MAT_DATA(this_f_index,
                                            joker,
                                            i_za_sca,
                                            i_aa_sca,
                                            i_za_inc,
                                            i_aa_inc,
                                            i),
                               t_gp);
        } else {
          pha_mat_data_int = PHA_MAT_DATA(
              this_f_index, this_T_index, joker, joker, joker, joker, joker);
          /*
                  for (Index i_za_sca = 0;
                       i_za_sca < PHA_MAT_DATA.nshelves(); i_za_sca++)
                    for (Index i_aa_sca = 0;
                         i_aa_sca < PHA_MAT_DATA.nbooks(); i_aa_sca++)
                      for (Index i_za_inc = 0;
                           i_za_inc < PHA_MAT_DATA.npages(); i_za_inc++)
                        for (Index i_aa_inc = 0;
                             i_aa_inc < PHA_MAT_DATA.nrows(); i_aa_inc++)
                          for (Index i = 0; i < PHA_MAT_DATA.ncols(); i++)
                            // Interpolation of phase matrix:
                            pha_mat_data_int(i_za_sca, i_aa_sca,
                                             i_za_inc, i_aa_inc, i) =
                                PHA_MAT_DATA(this_f_index, this_T_index,
                                                 i_za_sca, i_aa_sca,
                  */
        }

        // Do the transformation into the laboratory coordinate system.
        for (Index za_inc_idx = 0; za_inc_idx < za_grid.nelem();
             za_inc_idx++) {
          for (Index aa_inc_idx = 0; aa_inc_idx < aa_grid.nelem();
               aa_inc_idx++) {
            pha_matTransform(
                pha_mat_spt(i_se_flat, za_inc_idx, aa_inc_idx, joker, joker),
                pha_mat_data_int,
                ZA_DATAGRID,
                AA_DATAGRID,
                PART_TYPE,
                za_index,
                aa_index,
                za_inc_idx,
                aa_inc_idx,
                za_grid,
                aa_grid);
          }
        }
      }
      i_se_flat++;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesMerge(  //WS Output:
    Tensor4& pnd_field,
    ArrayOfArrayOfSingleScatteringData& scat_data,
    ArrayOfArrayOfScatteringMetaData& scat_meta,
    ArrayOfString& scat_species,
    Index& cloudbox_checked,
    //WS Input:
    const Index& cloudbox_on,
    const ArrayOfIndex& cloudbox_limits,
    const Tensor3& t_field,
    const Tensor3& z_field,
    const Matrix &z_surface) {
  // FIXME:
  // so far, this works for both scat_data and scat_data_raw. Needs to be
  // adapted, though, once we have WSM that can create Z/K/a with different
  // f/T dimensions than scat_data_single.f/T_grid.

  // cloudbox variable state should be ok before entering here
  if (!cloudbox_checked)
    throw std::runtime_error(
        "You must call *cloudbox_checkedCalc* before this method.");
  //however, we modify cloudbox variables. hence force re-checking the new
  //variables by resetting cloudbox_checked to False.
  cloudbox_checked = 0;

  if (3 != 1)
    throw std::runtime_error(
        "Merging scattering elements only works with a 1D atmoshere");

  // Cloudbox on/off?
  if (!cloudbox_on) {
    /* Must initialise pnd_field anyway; but empty */
    pnd_field.resize(0, 0, 0, 0);
    return;
  }

  // ------- setup new pnd_field and scat_data -------------------
  ArrayOfIndex limits(2);
  //pressure
  limits[0] = cloudbox_limits[0];
  limits[1] = cloudbox_limits[1] + 1;

  Tensor4 pnd_field_merged(
      limits[1] - limits[0], limits[1] - limits[0], 1, 1, 0.);

  ArrayOfArrayOfSingleScatteringData scat_data_merged;
  scat_data_merged.resize(1);
  scat_data_merged[0].resize(pnd_field_merged.nbooks());
  ArrayOfArrayOfScatteringMetaData scat_meta_merged;
  scat_meta_merged.resize(1);
  scat_meta_merged[0].resize(pnd_field_merged.nbooks());
  ArrayOfString scat_species_merged;
  scat_species_merged.resize(1);
  scat_species_merged[0] = "mergedfield-mergedpsd";
  for (Index sp = 0; sp < scat_data_merged[0].nelem(); sp++) {
    SingleScatteringData& this_part = scat_data_merged[0][sp];
    this_part.ptype = scat_data[0][0].ptype;
    this_part.description = "Merged scattering elements";
    this_part.f_grid = scat_data[0][0].f_grid;
    this_part.za_grid = scat_data[0][0].za_grid;
    this_part.aa_grid = scat_data[0][0].aa_grid;
    this_part.pha_mat_data.resize(scat_data[0][0].pha_mat_data.nlibraries(),
                                  1,
                                  scat_data[0][0].pha_mat_data.nshelves(),
                                  scat_data[0][0].pha_mat_data.nbooks(),
                                  scat_data[0][0].pha_mat_data.npages(),
                                  scat_data[0][0].pha_mat_data.nrows(),
                                  scat_data[0][0].pha_mat_data.ncols());
    this_part.ext_mat_data.resize(scat_data[0][0].ext_mat_data.nshelves(),
                                  1,
                                  scat_data[0][0].ext_mat_data.npages(),
                                  scat_data[0][0].ext_mat_data.nrows(),
                                  scat_data[0][0].ext_mat_data.ncols());
    this_part.abs_vec_data.resize(scat_data[0][0].abs_vec_data.nshelves(),
                                  1,
                                  scat_data[0][0].abs_vec_data.npages(),
                                  scat_data[0][0].abs_vec_data.nrows(),
                                  scat_data[0][0].abs_vec_data.ncols());
    this_part.pha_mat_data = 0.;
    this_part.ext_mat_data = 0.;
    this_part.abs_vec_data = 0.;
    this_part.T_grid.resize(1);
    this_part.T_grid[0] = t_field(sp, 0, 0);

    ScatteringMetaData& this_meta = scat_meta_merged[0][sp];
    std::ostringstream os;
    os << "Merged scattering element of cloudbox-level #" << sp;
    this_meta.description = os.str();
    this_meta.source = "ARTS internal";
    this_meta.refr_index = "Unknown";
    this_meta.mass = -1.;
    this_meta.diameter_max = -1.;
    this_meta.diameter_volume_equ = -1.;
    this_meta.diameter_area_equ_aerodynamical = -1.;
  }

  // Check that all scattering elements have same ptype and data dimensions
  SingleScatteringData& first_part = scat_data[0][0];
  for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++) {
    for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++) {
      SingleScatteringData& orig_part = scat_data[i_ss][i_se];

      if (orig_part.ptype != first_part.ptype)
        throw std::runtime_error(
            "All scattering elements must have the same type");

      if (orig_part.f_grid.nelem() != first_part.f_grid.nelem())
        throw std::runtime_error(
            "All scattering elements must have the same f_grid");

      if (!is_size(orig_part.pha_mat_data(
                       joker, 0, joker, joker, joker, joker, joker),
                   first_part.pha_mat_data.nlibraries(),
                   first_part.pha_mat_data.nshelves(),
                   first_part.pha_mat_data.nbooks(),
                   first_part.pha_mat_data.npages(),
                   first_part.pha_mat_data.nrows(),
                   first_part.pha_mat_data.ncols()))
        throw std::runtime_error(
            "All scattering elements must have the same pha_mat_data size"
            " (except for temperature).");

      if (!is_size(orig_part.ext_mat_data(joker, 0, joker, joker, joker),
                   first_part.ext_mat_data.nshelves(),
                   first_part.ext_mat_data.npages(),
                   first_part.ext_mat_data.nrows(),
                   first_part.ext_mat_data.ncols()))
        throw std::runtime_error(
            "All scattering elements must have the same ext_mat_data size"
            " (except for temperature).");

      if (!is_size(orig_part.abs_vec_data(joker, 0, joker, joker, joker),
                   first_part.abs_vec_data.nshelves(),
                   first_part.abs_vec_data.npages(),
                   first_part.abs_vec_data.nrows(),
                   first_part.abs_vec_data.ncols()))
        throw std::runtime_error(
            "All scattering elements must have the same abs_vec_data size"
            " (except for temperature).");
    }
  }

  //----- Start pnd_field_merged and scat_data_array_merged calculations -----

  GridPos T_gp;
  Vector itw(2);

  Index nlevels = pnd_field_merged.nbooks();
  // loop over pressure levels in cloudbox
  for (Index i_lv = 0; i_lv < nlevels - 1; i_lv++) {
    pnd_field_merged(i_lv, i_lv, 0, 0) = 1.;

    SingleScatteringData& this_part = scat_data_merged[0][i_lv];
    for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++) {
      for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++) {
        SingleScatteringData& orig_part = scat_data[i_ss][i_se];
        const Index pnd_index = FlattenedIndex(scat_data, i_ss, i_se);

        // If the particle number density at a specific point in the
        // atmosphere for the i_se scattering element is zero, we don't
        // need to do the transformation!
        if (pnd_field(pnd_index, i_lv, 0, 0) > PND_LIMIT)  //TRS
        {
          Numeric temperature = this_part.T_grid[0];
          if (orig_part.T_grid.nelem() > 1) {
            std::ostringstream os;
            os << "The temperature grid of the scattering data "
               << "does not cover the\n"
               << "atmospheric temperature at cloud location. "
               << "The data should\n"
               << "include the value T = " << temperature << " K.\n"
               << "Offending particle is scat_data[" << i_ss << "][" << i_se
               << "]:\n"
               << "Description: " << orig_part.description << "\n";
            chk_interpolation_grids(os.str(), orig_part.T_grid, temperature);

            // Gridpositions:
            gridpos(T_gp, orig_part.T_grid, temperature);
            // Interpolation weights:
            interpweights(itw, T_gp);
          }

          ////////// Extinction matrix and absorption vector

          // Loop over frequencies
          for (Index i_f = 0; i_f < orig_part.pha_mat_data.nlibraries();
               i_f++) {
            // Loop over zenith angles
            for (Index i_za = 0; i_za < orig_part.ext_mat_data.npages();
                 i_za++) {
              // Loop over azimuth angles
              for (Index i_aa = 0; i_aa < orig_part.ext_mat_data.nrows();
                   i_aa++) {
                // Weighted sum of ext_mat_data and abs_vec_data
                if (orig_part.T_grid.nelem() == 1) {
                  Vector v{orig_part.ext_mat_data(i_f, 0, i_za, i_aa, joker)};
                  v *= pnd_field(pnd_index, i_lv, 0, 0);
                  this_part.ext_mat_data(i_f, 0, i_za, 0, joker) += v;

                  v = orig_part.abs_vec_data(i_f, 0, i_za, i_aa, joker);
                  v *= pnd_field(pnd_index, i_lv, 0, 0);
                  this_part.abs_vec_data(i_f, 0, i_za, i_aa, joker) += v;
                } else {
                  for (Index i = 0; i < orig_part.ext_mat_data.ncols(); i++) {
                    // Temperature interpolation
                    this_part.ext_mat_data(i_f, 0, i_za, i_aa, i) +=
                        pnd_field(pnd_index, i_lv, 0, 0) *
                        interp(
                            itw,
                            orig_part.ext_mat_data(i_f, joker, i_za, i_aa, i),
                            T_gp);
                  }
                  for (Index i = 0; i < orig_part.abs_vec_data.ncols(); i++) {
                    // Temperature interpolation
                    this_part.abs_vec_data(i_f, 0, i_za, i_aa, i) +=
                        pnd_field(pnd_index, i_lv, 0, 0) *
                        interp(
                            itw,
                            orig_part.abs_vec_data(i_f, joker, i_za, i_aa, i),
                            T_gp);
                  }
                }
              }
            }

            ////////// Phase matrix

            // Loop over outgoing zenith angles
            for (Index i_za_out = 0;
                 i_za_out < orig_part.pha_mat_data.nshelves();
                 i_za_out++) {
              // Loop over outgoing azimuth angles
              for (Index i_aa_out = 0;
                   i_aa_out < orig_part.pha_mat_data.nbooks();
                   i_aa_out++) {
                // Loop over incoming zenith angles
                for (Index i_za_inc = 0;
                     i_za_inc < orig_part.pha_mat_data.npages();
                     i_za_inc++) {
                  // Loop over incoming azimuth angles
                  for (Index i_aa_inc = 0;
                       i_aa_inc < orig_part.pha_mat_data.nrows();
                       i_aa_inc++) {
                    // Weighted sum of pha_mat_data
                    if (orig_part.T_grid.nelem() == 1) {
                      Vector v{orig_part.pha_mat_data(i_f,
                                                        0,
                                                        i_za_out,
                                                        i_aa_out,
                                                        i_za_inc,
                                                        i_aa_inc,
                                                        joker)};
                      v *= pnd_field(pnd_index, i_lv, 0, 0);
                      this_part.pha_mat_data(i_f,
                                             0,
                                             i_za_out,
                                             i_aa_out,
                                             i_za_inc,
                                             i_aa_inc,
                                             joker) = v;
                    } else {
                      // Temperature interpolation
                      for (Index i = 0; i < orig_part.pha_mat_data.ncols();
                           i++) {
                        this_part.pha_mat_data(i_f,
                                               0,
                                               i_za_out,
                                               i_aa_out,
                                               i_za_inc,
                                               i_aa_inc,
                                               i) +=
                            pnd_field(pnd_index, i_lv, 0, 0) *
                            interp(itw,
                                   orig_part.pha_mat_data(i_f,
                                                          joker,
                                                          i_za_out,
                                                          i_aa_out,
                                                          i_za_inc,
                                                          i_aa_inc,
                                                          i),
                                   T_gp);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Set new pnd_field at lowest altitude to 0 if the cloudbox doesn't touch
  // the ground.
  // The consistency for the original pnd_field has already been ensured by
  // cloudbox_checkedCalc
  if (z_field(cloudbox_limits[0], 0, 0) > z_surface(0, 0))
    pnd_field_merged(0, 0, 0, 0) = 0.;

  pnd_field = pnd_field_merged;
  scat_data = scat_data_merged;
  scat_meta = scat_meta_merged;
  scat_species = scat_species_merged;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ExtractFromMetaSingleScatSpecies(
    //WS Output:
    Vector& meta_param,
    //WS Input:
    const ArrayOfArrayOfScatteringMetaData& scat_meta,
    const String& meta_name,
    const Index &scat_species_index) {
  if (scat_species_index < 0) {
    std::ostringstream os;
    os << "scat_species_index can't be <0!";
    throw std::runtime_error(os.str());
  }

  const Index nss = scat_meta.nelem();

  // check that scat_meta actually has at least scat_species_index elements
  if (!(nss > scat_species_index)) {
    std::ostringstream os;
    os << "Can not extract data for scattering species #" << scat_species_index
       << "\n"
       << "because scat_meta has only " << nss << " elements.";
    throw std::runtime_error(os.str());
  }

  const Index nse = scat_meta[scat_species_index].nelem();
  meta_param.resize(nse);

  for (Index i = 0; i < nse; i++) {
    if (meta_name == "mass")
      meta_param[i] = scat_meta[scat_species_index][i].mass;
    else if (meta_name == "diameter_max")
      meta_param[i] = scat_meta[scat_species_index][i].diameter_max;
    else if (meta_name == "diameter_volume_equ")
      meta_param[i] = scat_meta[scat_species_index][i].diameter_volume_equ;
    else if (meta_name == "diameter_area_equ_aerodynamical")
      meta_param[i] =
          scat_meta[scat_species_index][i].diameter_area_equ_aerodynamical;
    else {
      std::ostringstream os;
      os << "Meta parameter \"" << meta_name << "\"is unknown.";
      throw std::runtime_error(os.str());
    }
  }
}


