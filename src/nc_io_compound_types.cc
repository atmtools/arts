////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   nc_io_basic_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-09-26

  \brief This file contains functions to handle NetCDF data files.

*/

#include "config.h"

#ifdef ENABLE_NETCDF

#include <cstring>

#include "arts.h"
#include "nc_io.h"
#include "nc_io_types.h"

//=== GasAbsLookup ==========================================================

//! Reads a GasAbsLookup table from a NetCDF file
/*!
 \param[in] ncid    NetCDF file descriptor
 \param[in] gal     GasAbsLookup
 
 \author Oliver Lemke
*/
void nca_read_from_file(const int ncid, GasAbsLookup& gal, const Verbosity&) {
  nca_get_data(ncid, "species", gal.species, true);

  ARTS_USER_ERROR_IF(!gal.species.nelem(),
                     "No species found in lookup table file!");

  nca_get_data(
      ncid, "nonlinear_species", gal.nonlinear_species, true);
  nca_get_data(ncid, "f_grid", gal.f_grid, true);
  nca_get_data(ncid, "p_grid", gal.p_grid, true);
  nca_get_data(ncid, "vmrs_ref", gal.vmrs_ref, true);
  nca_get_data(ncid, "t_ref", gal.t_ref, true);
  nca_get_data(ncid, "t_pert", gal.t_pert, true);
  nca_get_data(ncid, "nls_pert", gal.nls_pert, true);
  nca_get_data(ncid, "xsec", gal.xsec, true);
}

//! Writes a GasAbsLookup table to a NetCDF file
/*!
 \param[in]  ncid    NetCDF file descriptor
 \param[out] gal     GasAbsLookup
 
 \author Oliver Lemke
*/
void nca_write_to_file(const int ncid,
                       const GasAbsLookup& gal,
                       const Verbosity&) {
  int retval;

  int species_strings_varid;
  int species_count_varid;

  ArrayOfIndex species_count(gal.species.nelem());
  Index species_max_strlen = 0;
  char* species_strings = nullptr;

  ARTS_USER_ERROR_IF(!gal.species.nelem(),
                     "Current lookup table contains no species!");

  long species_total_nelems = 0;
  for (Index nspecies = 0; nspecies < gal.species.nelem(); nspecies++) {
    Index nspecies_nelem = gal.species[nspecies].nelem();
    species_total_nelems += nspecies_nelem;
    species_count[nspecies] = nspecies_nelem;

    for (const auto &it : gal.species[nspecies])
      if (it.Name().nelem() > species_max_strlen)
        species_max_strlen = it.Name().nelem();
  }
  species_max_strlen++;

  species_strings = new char[species_total_nelems * species_max_strlen];
  memset(species_strings, 0, species_total_nelems * species_max_strlen);

  Index str_i = 0;
  for (const auto &species : gal.species)
    for (const auto &it2 : species) {
      memccpy(&species_strings[str_i], it2.Name().c_str(), 0,
              species_max_strlen);
      str_i += species_max_strlen;
    }

  species_count_varid =
      nca_def_ArrayOfIndex(ncid, "species_count", species_count);

  std::array<int, 2> species_strings_ncdims;
  nca_def_dim(ncid, "species_strings_nelem", species_total_nelems,
              &species_strings_ncdims[0]);
  nca_def_dim(ncid, "species_strings_length", species_max_strlen,
              &species_strings_ncdims[1]);
  nca_def_var(ncid, "species_strings", NC_CHAR, 2, &species_strings_ncdims[0],
              &species_strings_varid);

  // Define dimensions and variables
  int nonlinear_species_varid =
      nca_def_ArrayOfIndex(ncid, "nonlinear_species", gal.nonlinear_species);
  int f_grid_varid = nca_def_Vector(ncid, "f_grid", gal.f_grid);
  int p_grid_varid = nca_def_Vector(ncid, "p_grid", gal.p_grid);
  int vmrs_ref_varid = nca_def_Matrix(ncid, "vmrs_ref", gal.vmrs_ref);
  int t_ref_varid = nca_def_Vector(ncid, "t_ref", gal.t_ref);
  int t_pert_varid = nca_def_Vector(ncid, "t_pert", gal.t_pert);
  int nls_pert_varid = nca_def_Vector(ncid, "nls_pert", gal.nls_pert);
  int xsec_varid = nca_def_Tensor4(ncid, "xsec", gal.xsec);

  if ((retval = nc_enddef(ncid))) nca_error(retval, "nc_enddef");

  // Write variables
  nca_put_var(ncid, species_count_varid, species_count);
  if (gal.species.nelem()) {
    if ((retval =
             nc_put_var_text(ncid, species_strings_varid, species_strings)))
      nca_error(retval, "nc_put_var");
  }

  delete[] species_strings;

  nca_put_var(
      ncid, nonlinear_species_varid, gal.nonlinear_species);
  nca_put_var(ncid, f_grid_varid, gal.f_grid);
  nca_put_var(ncid, p_grid_varid, gal.p_grid);
  nca_put_var(ncid, vmrs_ref_varid, gal.vmrs_ref);
  nca_put_var(ncid, t_ref_varid, gal.t_ref);
  nca_put_var(ncid, t_pert_varid, gal.t_pert);
  nca_put_var(ncid, nls_pert_varid, gal.nls_pert);
  nca_put_var(ncid, xsec_varid, gal.xsec);
}

////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

#define TMPL_NC_READ_WRITE_FILE_DUMMY(what)                                    \
  void nca_write_to_file(const int, const what &, const Verbosity &) {         \
    ARTS_USER_ERROR("NetCDF support not yet implemented for this type!");      \
  }                                                                            \
  void nca_read_from_file(const int, what &, const Verbosity &) {              \
    ARTS_USER_ERROR("NetCDF support not yet implemented for this type!");      \
  }

TMPL_NC_READ_WRITE_FILE_DUMMY(Agenda)

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE_DUMMY

#endif /* ENABLE_NETCDF */
