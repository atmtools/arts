/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_cloudbox.cc
  \author Patrick Eriksson, Claudia Emde and Sreerekha T. R.
  \date   2002-05-08 

  \brief  Workspace functions related to the definintion of the cloud box.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.*/

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <workspace.h>

#include <stdexcept>

#include "array.h"
#include "arts_constants.h"
#include "arts_conversions.h"
#include "arts_omp.h"
#include "atm.h"
#include "check_input.h"
#include "cloudbox.h"
#include "debug.h"
#include "file.h"
#include "jacobian.h"
#include "math_funcs.h"
#include "optproperties.h"
#include "species_tags.h"
#include "xml_io.h"

inline constexpr Numeric PI = Constant::pi;
inline constexpr Numeric DEG2RAD = Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG = Conversion::rad2deg(1);
inline constexpr Numeric DENSITY_OF_ICE = Constant::density_of_ice_at_0c;

/* Workspace method: Doxygen documentation will be auto-generated */
void particle_fieldCleanup(  //WS Output:
    Tensor4& particle_field_out,
    //WS Input:
    const Tensor4& particle_field_in,
    const Numeric& threshold) {
  if (&particle_field_out != &particle_field_in) {
    particle_field_out = particle_field_in;
  }

  // Check that particle_field contains realistic values
  //(values smaller than the threshold will be set to 0)
  for (Index i = 0; i < particle_field_out.nbooks(); i++) {
    for (Index j = 0; j < particle_field_out.npages(); j++) {
      for (Index k = 0; k < particle_field_out.nrows(); k++) {
        for (Index l = 0; l < particle_field_out.ncols(); l++) {
          if (particle_field_out(i, j, k, l) < threshold) {
            particle_field_out(i, j, k, l) = 0.0;
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesInit(  //WS Output:
    ArrayOfString& scat_species,
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    ArrayOfArrayOfScatteringMetaData& scat_meta,
    Index& scat_data_checked,
    ArrayOfGriddedField3& pnd_field_raw) {
  scat_species.resize(0);
  scat_data_raw.resize(0);
  scat_meta.resize(0);
  pnd_field_raw.resize(0);
  scat_data_checked = 0;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatElementsPndAndScatAdd(  //WS Output:
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    ArrayOfGriddedField3& pnd_field_raw,
    // WS Input (needed for checking the datafiles):
    // Keywords:
    const ArrayOfString& scat_data_files,
    const ArrayOfString& pnd_field_files) {
  //--- Check input ---------------------------------------------------------

  // Atmosphere
  //chk_atm_grids( 3, p_grid, lat_grid, lon_grid );

  //--- Reading the data ---------------------------------------------------

  ARTS_USER_ERROR_IF(
      scat_data_files.size() != pnd_field_files.size(),
      "Number of elements in scat_data and pnd_field filelists is"
      "inconsistent.")

  Index last_species = scat_data_raw.size() - 1;
  if (last_species == -1) {
    scat_data_raw.resize(1);
    last_species = 0;
  }

  // Create empty dummy elements to append to *scat_data_raw* and *pnd_field_raw*.
  SingleScatteringData scat_data_single;
  GriddedField3 pnd_field_data;

  for (Size i = 0; i < scat_data_files.size(); i++) {
    // Append *scat_data_raw* and *pnd_field_raw* with empty Arrays of Tensors.
    scat_data_raw[last_species].push_back(scat_data_single);
    pnd_field_raw.push_back(pnd_field_data);

    xml_read_from_file(
        scat_data_files[i],
        scat_data_raw[last_species][scat_data_raw[last_species].size() - 1]);

    if (pnd_field_files[i].size() < 1) {
    } else {
      xml_read_from_file(pnd_field_files[i],
                         pnd_field_raw[pnd_field_raw.size() - 1]);

      chk_pnd_data(pnd_field_raw[pnd_field_raw.size() - 1], pnd_field_files[i]);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesPndAndScatAdd(  //WS Output:
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    ArrayOfGriddedField3& pnd_field_raw,
    // WS Input(needed for checking the datafiles):
    // Keywords:
    const ArrayOfString& scat_data_files,
    const String& pnd_fieldarray_file) {
  //--- Check input ---------------------------------------------------------

  // Atmosphere
  //chk_atm_grids ( 3, p_grid, lat_grid, lon_grid );

  //--- Reading the data ---------------------------------------------------
  ArrayOfSingleScatteringData arr_ssd;
  arr_ssd.resize(scat_data_files.size());

  for (Size i = 0; i < scat_data_files.size(); i++) {
    xml_read_from_file(scat_data_files[i], arr_ssd[i]);
  }

  // append as new scattering species
  if (scat_data_raw.size() == 0) {
    scat_data_raw.resize(1);
    scat_data_raw[0] = arr_ssd;
  } else
    scat_data_raw.push_back(arr_ssd);

  ArrayOfGriddedField3 pnd_tmp;
  xml_read_from_file(pnd_fieldarray_file, pnd_tmp);

  chk_pnd_raw_data(pnd_tmp, pnd_fieldarray_file);

  // append to pnd_field_raw
  for (Size i = 0; i < pnd_tmp.size(); ++i) pnd_field_raw.push_back(pnd_tmp[i]);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatElementsToabs_speciesAdd(  //WS Output:
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    AtmField& atm_field,
    ArrayOfArrayOfSpeciesTag& abs_species,
    // WS Input (needed for checking the datafiles):
    const Vector& f_grid,
    // Keywords:
    const ArrayOfString& scat_data_files,
    const ArrayOfString& pnd_field_files) {
  //--- Check input ---------------------------------------------------------

  // Atmosphere
  //chk_atm_grids( 3, p_grid, lat_grid, lon_grid );

  // Frequency grid
  ARTS_USER_ERROR_IF(f_grid.empty(), "The frequency grid is empty.");
  chk_if_increasing("f_grid", f_grid);

  //--- Reading the data ---------------------------------------------------

  ARTS_USER_ERROR_IF(
      scat_data_files.size() != pnd_field_files.size(),
      "Number of elements in scat_data and pnd_field filelists is"
      "inconsistent.")

  Index last_species = scat_data_raw.size() - 1;
  if (last_species == -1) {
    scat_data_raw.resize(1);
    last_species = 0;
  }

  // Create empty dummy elements to append to *scat_data_raw* and *pnd_field_raw*.
  SingleScatteringData scat_data_single;
  GriddedField3 pnd_field_data;
  ArrayOfString species(1);
  species[0] = "particles";

  absorption_speciesAdd(abs_species, species);

  for (Size i = 0; i < scat_data_files.size(); i++) {
    // Append *scat_data_raw* and *pnd_field_raw* with empty Arrays of Tensors.
    scat_data_raw[last_species].push_back(scat_data_single);
    atm_field[abs_species.back().Species()] = pnd_field_data;

    xml_read_from_file(
        scat_data_files[i],
        scat_data_raw[last_species][scat_data_raw[last_species].size() - 1]);

    chk_interpolation_grids(
        "scat_data_single.f_grid to f_grid",
        scat_data_raw[last_species][scat_data_raw[last_species].size() - 1]
            .f_grid,
        f_grid);

    if (pnd_field_files[i].size() < 1) {
    } else {
      try {
        xml_read_from_file(
            pnd_field_files[i],
            atm_field[abs_species.back().Species()].get<GriddedField3&>());
      } catch (...) {
        ArrayOfGriddedField3 tmp;
        try {
          xml_read_from_file(pnd_field_files[i], tmp);
          if (tmp.size() == 1) {
            atm_field[abs_species.back().Species()] = tmp[0];
          } else {
            ARTS_USER_ERROR(
                "The file ",
                pnd_field_files[i],
                "\n"
                "is neither GriddedField3 nor a 1-long ArrayOfGriddedField3.\n")
          }
        } catch (...) {
          ARTS_USER_ERROR(
              "The file ",
              pnd_field_files[i],
              " does not exist or\n"
              "its type is neither GriddedField3 nor a 1-long ArrayOfGriddedField3.\n")
        }
      }

      chk_pnd_data(
          atm_field[abs_species.back().Species()].get<const GriddedField3&>(),
          pnd_field_files[i]);
    }
  }
  assert(false);
  // scat_dataCheck(scat_data_raw, "sane", 1e-2);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesScatAndMetaRead(  //WS Output:
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    ArrayOfArrayOfScatteringMetaData& scat_meta,
    // Keywords:
    const ArrayOfString& scat_data_files) {
  //--- Reading the data ---------------------------------------------------
  ArrayOfSingleScatteringData arr_ssd;
  ArrayOfScatteringMetaData arr_smd;

  arr_ssd.resize(scat_data_files.size());
  arr_smd.resize(scat_data_files.size());

  Index meta_naming_conv = 0;

  for (Size i = 0; i < 1 && i < scat_data_files.size(); i++) {
    xml_read_from_file(scat_data_files[i], arr_ssd[i]);

    // make meta data name from scat data name
    ArrayOfString strarr;
    String scat_meta_file;

    if (i == 0) {
      split(strarr, scat_data_files[i], ".xml");
      scat_meta_file = strarr[0] + ".meta.xml";

      try {
        find_xml_file(scat_meta_file);
      } catch (const std::runtime_error&) {
      }

      if (file_exists(scat_meta_file)) {
        xml_read_from_file(scat_meta_file, arr_smd[i]);

        meta_naming_conv = 1;
      } else {
        try {
          split(strarr, scat_data_files[i], "scat_data");
          ARTS_USER_ERROR_IF(
              strarr.size() < 2,
              "Splitting scattering data filename up at 'scat_data' also failed.");
          scat_meta_file = strarr[0] + "scat_meta" + strarr[1];

          xml_read_from_file(scat_meta_file, arr_smd[i]);

          meta_naming_conv = 2;
        } catch (const std::runtime_error& e) {
          ARTS_USER_ERROR(
              "No meta data file following one of the allowed naming "
              "conventions was found.\n"
              "Allowed are "
              "*.meta.xml from *.xml and "
              "*scat_meta* from *scat_data*\n"
              "Scattering meta data file not found: ",
              scat_meta_file,
              "\n",
              e.what())
        }
      }
    }
  }

  ArrayOfString fail_msg;

#pragma omp parallel for if (!arts_omp_in_parallel() &&                       \
                                 scat_data_files.size() > 1)                  \
    num_threads(arts_omp_get_max_threads() > 16 ? 16                          \
                                                : arts_omp_get_max_threads()) \
    shared(scat_data_files, arr_ssd, arr_smd)
  for (Size i = 1; i < scat_data_files.size(); i++) {
    // make meta data name from scat data name
    ArrayOfString strarr;
    String scat_meta_file;
    SingleScatteringData ssd;
    ScatteringMetaData smd;

    try {
      xml_read_from_file(scat_data_files[i], ssd);

      split(strarr, scat_data_files[i], ".xml");
      scat_meta_file = strarr[0] + ".meta.xml";

      if (meta_naming_conv == 1) {
        split(strarr, scat_data_files[i], ".xml");
        scat_meta_file = strarr[0] + ".meta.xml";

        xml_read_from_file(scat_meta_file, smd);
      } else {
        split(strarr, scat_data_files[i], "scat_data");
        scat_meta_file = strarr[0] + "scat_meta" + strarr[1];

        xml_read_from_file(scat_meta_file, smd);
      }
    } catch (const std::exception& e) {
      std::ostringstream os;
      os << "Run-time error reading scattering data : \n" << e.what();
#pragma omp critical(ybatchCalc_push_fail_msg)
      fail_msg.push_back(os.str());
    }

#pragma omp critical(ScatSpeciesScatAndMetaRead_assign_ssd)
    arr_ssd[i] = std::move(ssd);
#pragma omp critical(ScatSpeciesScatAndMetaRead_assign_smd)
    arr_smd[i] = std::move(smd);
  }

  if (fail_msg.size()) {
    std::ostringstream os;
    for (auto& msg : fail_msg) os << msg << '\n';

    ARTS_USER_ERROR("{}", os.str());
  }

  // check if arrays have same size
  chk_scattering_data(arr_ssd, arr_smd);

  // append as new scattering species
  scat_data_raw.push_back(std::move(arr_ssd));
  scat_meta.push_back(std::move(arr_smd));
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatElementsSelect(  //WS Output:
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    ArrayOfArrayOfScatteringMetaData& scat_meta,
    // WS Input:
    const ArrayOfString& scat_species,
    const String& species,
    const String& sizeparam,
    const Numeric& sizemin,
    const Numeric& sizemax,
    const Numeric& tolerance,
    const String& delim) {
  // first check that sizes of scat_species and scat_data_raw/scat_meta agree
  const Size nspecies = scat_species.size();
  ARTS_USER_ERROR_IF(
      nspecies != scat_data_raw.size() || nspecies != scat_meta.size(),
      "Number of scattering species specified by scat_species does\n"
      "not agree with number of scattering species in\n"
      "scat_data_raw or scat_meta:\n"
      "scat_species has ",
      nspecies,
      " entries, while scat_data_raw has ",
      scat_data_raw.size(),
      " and scat_meta has ",
      scat_meta.size(),
      ".")

  // create temporary containers for selected elements
  ArrayOfSingleScatteringData scat_data_raw_tmp;
  ArrayOfScatteringMetaData scat_meta_tmp;

  String partfield_name;
  //find the species to handle: compare 'species' to 'partfield' part of
  //scat_species tags
  Index i_ss = -1;
  for (Size i = 0; i < scat_species.size(); i++) {
    parse_partfield_name(partfield_name, scat_species[i], delim);
    if (partfield_name == species) i_ss = i;
  }
  ARTS_USER_ERROR_IF(i_ss < 0,
                     "Scattering species ",
                     species,
                     " not found among scat_species.")

  // choosing the specified SingleScatteringData and ScatteringMetaData
  if (sizeparam == "diameter_max")
    for (Size i_se = 0; i_se < scat_meta[i_ss].size(); i_se++) {
      // scattering element diameter is extracted from the
      // scattering element's meta data and checked whether it's within size
      // selected range (sizemax < 0 check follows from wildcard usage and
      // means consider all sizes on the upper end)
      if (scat_meta[i_ss][i_se].diameter_max > sizemin - sizemin * tolerance &&
          (sizemax + sizemax * tolerance > scat_meta[i_ss][i_se].diameter_max ||
           sizemax < 0.)) {
        // copy selected scattering element to temp arrays
        scat_data_raw_tmp.push_back(scat_data_raw[i_ss][i_se]);
        scat_meta_tmp.push_back(scat_meta[i_ss][i_se]);
      }
    }
  else if (sizeparam == "diameter_volume_equ")
    for (Size i_se = 0; i_se < scat_meta[i_ss].size(); i_se++) {
      if (scat_meta[i_ss][i_se].diameter_volume_equ >
              sizemin - sizemin * tolerance &&
          (sizemax + sizemax * tolerance >
               scat_meta[i_ss][i_se].diameter_volume_equ ||
           sizemax < 0.)) {
        // copy selected scattering element to temp arrays
        scat_data_raw_tmp.push_back(scat_data_raw[i_ss][i_se]);
        scat_meta_tmp.push_back(scat_meta[i_ss][i_se]);
      }
    }
  else if (sizeparam == "diameter_area_equ_aerodynamical")
    for (Size i_se = 0; i_se < scat_meta[i_ss].size(); i_se++) {
      if (scat_meta[i_ss][i_se].diameter_area_equ_aerodynamical >
              sizemin - sizemin * tolerance &&
          (sizemax + sizemax * tolerance >
               scat_meta[i_ss][i_se].diameter_area_equ_aerodynamical ||
           sizemax < 0.)) {
        // copy selected scattering element to temp arrays
        scat_data_raw_tmp.push_back(scat_data_raw[i_ss][i_se]);
        scat_meta_tmp.push_back(scat_meta[i_ss][i_se]);
      }
    }
  else {
    ARTS_USER_ERROR("Size parameter ", sizeparam, "is unknown.")
  }

  // To use a particle species field without associated scattering element
  // data poses a high risk of accidentially neglecting these species. That's
  // unlikely what the user intends. Hence throw error.
  ARTS_USER_ERROR_IF(
      scat_meta_tmp.size() < 1,
      "For scattering species ",
      species,
      " no scattering "
      "element matching the requested size range found.\n"
      "Check *scat_data_raw* and *scat_meta* input as well as your size limit "
      "selection!")

  scat_meta[i_ss] = std::move(scat_meta_tmp);
  scat_data_raw[i_ss] = std::move(scat_data_raw_tmp);

  // check if array is empty. should never apply (since we checked the re-worked
  // data before and that error should also catch cases that are empty from the
  // beginning).
  ARTS_ASSERT(TotalNumberOfElements(scat_meta));
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesExtendTemperature(  //WS Output:
    ArrayOfArrayOfSingleScatteringData& scat_data_raw,
    // Keywords:
    const ArrayOfString& scat_species,
    const String& species,
    const String& scat_species_delim,
    const Numeric& T_low,
    const Numeric& T_high) {
  const bool do_tl = (T_low >= 0.);
  const bool do_th = (T_high >= 0.);

  if (do_tl || do_th) {
    Index i_ss = -1;
    if (species == "") {
      i_ss = scat_data_raw.size() - 1;
      ARTS_USER_ERROR_IF(
          i_ss == -1,
          "No *scat_data* available. Can not extend temperature range on "
          "inexistent data.")
    } else {
      // first check that sizes of scat_species and scat_data_raw agree
      Size nspecies = scat_species.size();
      ARTS_USER_ERROR_IF(
          nspecies != scat_data_raw.size(),
          "Number of scattering species specified by scat_species does\n"
          "not agree with number of scattering species in *scat_data*:\n"
          "scat_species has ",
          nspecies,
          " entries, while *scat_data* has ",
          scat_data_raw.size(),
          ".")
      String partfield_name;
      //find the species to handle: compare 'species' to 'partfield' part of
      //scat_species tags
      for (Size i = 0; i < scat_species.size(); i++) {
        parse_partfield_name(
            partfield_name, scat_species[i], scat_species_delim);
        if (partfield_name == species) i_ss = i;
      }
      ARTS_USER_ERROR_IF(i_ss < 0,
                         "Scattering species ",
                         species,
                         " not found among scat_species.")
    }

    for (Size i_se = 0; i_se < scat_data_raw[i_ss].size(); i_se++) {
      const SingleScatteringData& ssdo = scat_data_raw[i_ss][i_se];
      const Index nTo = ssdo.T_grid.nelem();
      Index nTn = nTo;
      bool do_htl, do_hth;
      if (nTo > 1) {
        do_htl = (do_tl && (T_low < ssdo.T_grid[0]));
        do_hth = (do_th && (T_high > last(ssdo.T_grid)));
      } else {
        do_htl = false;
        do_hth = false;
      }

      if (do_htl || do_hth) {
        // create new instance of SingleScatteringData
        SingleScatteringData ssdn;
        Index iToff;

        // determine new temperature grid
        if (do_htl) nTn += 1;
        if (do_hth) nTn += 1;
        Vector T_grid_new(nTn);
        if (do_htl) {
          T_grid_new[0] = T_low;
          iToff = 1;
        } else {
          iToff = 0;
        }
        for (Index iT = 0; iT < nTo; iT++)
          T_grid_new[iT + iToff] = scat_data_raw[i_ss][i_se].T_grid[iT];
        if (do_hth) T_grid_new[nTo + iToff] = T_high;
        ssdn.T_grid = std::move(T_grid_new);

        // copy grids and other descriptive data that is to remain identical
        ssdn.ptype = ssdo.ptype;
        std::ostringstream description;
        description << ssdo.description;  // here just copy. we append further
                                          // info below if applicable.
        ssdn.f_grid = ssdo.f_grid;
        ssdn.za_grid = ssdo.za_grid;
        ssdn.aa_grid = ssdo.aa_grid;

        // determine size of current optical property data
        const Index nf = ssdo.f_grid.nelem();
        const Index nzas = ssdo.pha_mat_data.nshelves();
        const Index naas = ssdo.pha_mat_data.nbooks();
        const Index nzai = ssdo.pha_mat_data.npages();
        const Index naai = ssdo.pha_mat_data.nrows();
        const Index nmep = ssdo.pha_mat_data.ncols();
        const Index nmee = ssdo.ext_mat_data.ncols();
        const Index nvea = ssdo.abs_vec_data.ncols();

        // create containers for extended optical property data
        ssdn.pha_mat_data.resize(nf, nTn, nzas, naas, nzai, naai, nmep);
        ssdn.ext_mat_data.resize(nf, nTn, nzai, naai, nmee);
        ssdn.abs_vec_data.resize(nf, nTn, nzai, naai, nvea);

        // copy optical property data
        for (Index iT = 0; iT < nTo; iT++) {
          ssdn.pha_mat_data(
              joker, iT + iToff, joker, joker, joker, joker, joker) =
              ssdo.pha_mat_data(joker, iT, joker, joker, joker, joker, joker);
          ssdn.ext_mat_data(joker, iT + iToff, joker, joker, joker) =
              ssdo.ext_mat_data(joker, iT, joker, joker, joker);
          ssdn.abs_vec_data(joker, iT + iToff, joker, joker, joker) =
              ssdo.abs_vec_data(joker, iT, joker, joker, joker);
        }

        // duplicate optical property data on T-edges if applicable
        if (do_htl) {
          ssdn.pha_mat_data(joker, 0, joker, joker, joker, joker, joker) =
              ssdn.pha_mat_data(joker, 1, joker, joker, joker, joker, joker);
          ssdn.ext_mat_data(joker, 0, joker, joker, joker) =
              ssdn.ext_mat_data(joker, 1, joker, joker, joker);
          ssdn.abs_vec_data(joker, 0, joker, joker, joker) =
              ssdn.abs_vec_data(joker, 1, joker, joker, joker);
          description << "\n"
                      << "Low temperature limit extended by"
                      << " duplicating previous low temperature limit"
                      << " single scattering properties.";
        }
        if (do_hth) {
          ssdn.pha_mat_data(joker, nTn - 1, joker, joker, joker, joker, joker) =
              ssdn.pha_mat_data(
                  joker, nTn - 2, joker, joker, joker, joker, joker);
          ssdn.ext_mat_data(joker, nTn - 1, joker, joker, joker) =
              ssdn.ext_mat_data(joker, nTn - 2, joker, joker, joker);
          ssdn.abs_vec_data(joker, nTn - 1, joker, joker, joker) =
              ssdn.abs_vec_data(joker, nTn - 2, joker, joker, joker);
          description << "\n"
                      << "High temperature limit extended by"
                      << " duplicating previous high temperature limit"
                      << " single scattering properties.";
        }
        ssdn.description = description.str();
        scat_data_raw[i_ss][i_se] = std::move(ssdn);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldExpand1D(Tensor4& pnd_field,
                       const Index& cloudbox_on,
                       const ArrayOfIndex& cloudbox_limits,
                       const Index& nzero) {
  ARTS_USER_ERROR_IF(!cloudbox_on,
                     "No use in calling this method with cloudbox off.");
  ARTS_USER_ERROR_IF(nzero < 1, "The argument *nzero* must be > 0.");

  // Sizes
  const Index npart = pnd_field.nbooks();
  const Index np = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index nlat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
  Index nlon = cloudbox_limits[5] - cloudbox_limits[4] + 1;

  ARTS_USER_ERROR_IF(pnd_field.npages() != np || pnd_field.nrows() != 1 ||
                         pnd_field.ncols() != 1,
                     "The input *pnd_field* is either not 1D or does not "
                     "match pressure size of cloudbox.");

  // Temporary container
  Tensor4 pnd_temp = pnd_field;

  // Resize and fill
  pnd_field.resize(npart, np, nlat, nlon);
  pnd_field = 0;
  //
  for (Index ilon = nzero; ilon < nlon - nzero; ilon++) {
    for (Index ilat = nzero; ilat < nlat - nzero; ilat++) {
      for (Index ip = 0; ip < np; ip++) {
        for (Index is = 0; is < npart; is++) {
          pnd_field(is, ip, ilat, ilon) = pnd_temp(is, ip, 0, 0);
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldZero(  //WS Output:
    Tensor4& pnd_field,
    ArrayOfTensor4& dpnd_field_dx,
    ArrayOfArrayOfSingleScatteringData& scat_data,
    //WS Input:
    const Vector& f_grid,
    const ArrayOfIndex& cloudbox_limits,
    const JacobianTargets& jacobian_targets) {
  ARTS_USER_ERROR_IF(cloudbox_limits.size() != 2 * 3,
                     "*cloudbox_limits* is a vector which contains the"
                     "upper and lower limit of the cloud for all "
                     "atmospheric dimensions. So its dimension must"
                     "be 2 x *3*");

  // Resize pnd_field and set it to 0:
  Index np = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  Index nlat = 1, nlon = 1;
  nlat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
  nlon = cloudbox_limits[5] - cloudbox_limits[4] + 1;

  // no (cloudy) Jacobians with this WSM, hence no setting.
  // but we need to size dpnd_field to be consistent with jacobian_quantities.
  dpnd_field_dx.resize(jacobian_targets.target_count());

  // Do only reset scat_data if it has not been set yet.
  // There's no need otherwise, and it's rather unpractical for testing when
  // doing so: we might want to do some actual calcs with the scat_data later
  // on. So why throw it away?
  const Index N_se = TotalNumberOfElements(scat_data);
  if (N_se > 0) {
    pnd_field.resize(N_se, np, nlat, nlon);
  } else {
    pnd_field.resize(1, np, nlat, nlon);

    //Resize scat_data and set it to 0:
    // Number of scattering elements
    scat_data.resize(1);
    scat_data[0].resize(1);
    scat_data[0][0].ptype = PTYPE_TOTAL_RND;
    scat_data[0][0].description = " ";
    // Grids which contain full ranges which one wants to calculate
    Index nf = f_grid.nelem();
    scat_data[0][0].f_grid.resize(nf);
    scat_data[0][0].f_grid = f_grid;
    Index nT = 1;
    scat_data[0][0].T_grid.resize(nT);
    scat_data[0][0].T_grid = 270.;
    Index nza = 5;
    nlinspace(scat_data[0][0].za_grid, 0, 180, nza);
    // Resize the data arrays
    scat_data[0][0].pha_mat_data.resize(nf, nT, nza, 1, 1, 1, 6);
    scat_data[0][0].pha_mat_data = 0.;
    scat_data[0][0].ext_mat_data.resize(nf, nT, 1, 1, 1);
    scat_data[0][0].ext_mat_data = 0.;
    scat_data[0][0].abs_vec_data.resize(nf, nT, 1, 1, 1);
    scat_data[0][0].abs_vec_data = 0.;
  }

  pnd_field = 0.;
}
