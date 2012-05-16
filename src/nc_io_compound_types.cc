/* Copyright (C) 2003-2008 Oliver Lemke <olemke@core-dump.info>

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
   USA. */


////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   nc_io_basic_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-09-26

  \brief This file contains functions to handle NetCDF data files.

*/

#include "config_global.h"

#ifdef ENABLE_NETCDF

#include <cstring>

#include "arts.h"
#include "nc_io.h"
#include "nc_io_types.h"


//=== GasAbsLookup ==========================================================

//! Reads a GasAbsLookup table from a NetCDF file
/*!
 \param ncf     NetCDF file descriptor
 \param gal     GasAbsLookup
 */
void
nc_read_from_file (const int ncid,
                   GasAbsLookup& gal)
{
  nc_get_data_ArrayOfArrayOfSpeciesTag(ncid, "species", gal.species, true);
  nc_get_data_ArrayOfIndex(ncid, "nonlinear_species", gal.nonlinear_species, true);
  nc_get_data_Vector(ncid, "f_grid", gal.f_grid, true);
  nc_get_data_Vector(ncid, "p_grid", gal.p_grid, true);  
  nc_get_data_Matrix(ncid, "vmrs_ref", gal.vmrs_ref, true);
  nc_get_data_Vector(ncid, "t_ref", gal.t_ref, true);
  nc_get_data_Vector(ncid, "t_pert", gal.t_pert, true);
  nc_get_data_Vector(ncid, "nls_pert", gal.nls_pert, true);
  nc_get_data_Tensor4(ncid, "xsec", gal.xsec, true);
}


//! Writes a GasAbsLookup table to a NetCDF file
/*!
 \param ncf     NetCDF file descriptor
 \param gal     GasAbsLookup
 */
void
nc_write_to_file (const int ncid,
                  const GasAbsLookup& gal)
{
  int retval;
  
  int species_strings_ncdims[2], species_strings_varid;
  int species_nelem_ncdims[1], species_nelem_varid;
  long species_total_nelems = 0;
  Index* per_species_nelem = new Index[gal.species.nelem()];
  Index species_max_strlen = 0;
  for(Index nspecies = 0; nspecies < gal.species.nelem(); nspecies++)
  {
    Index nspecies_nelem = gal.species[nspecies].nelem();
    species_total_nelems += nspecies_nelem;
    per_species_nelem[nspecies] = nspecies_nelem;
    
    for (ArrayOfSpeciesTag::const_iterator it = gal.species[nspecies].begin();
         it != gal.species[nspecies].end(); it++)
      if (it->Name().nelem() > species_max_strlen) species_max_strlen = it->Name().nelem();
  }
  species_max_strlen++;
  
  char* species_strings = new char[species_total_nelems*species_max_strlen];
  memset(species_strings, 0, species_total_nelems*species_max_strlen);

  Index str_i = 0;
  for (ArrayOfArrayOfSpeciesTag::const_iterator it1 = gal.species.begin();
       it1 != gal.species.end(); it1++)
    for (ArrayOfSpeciesTag::const_iterator it2 = it1->begin();
         it2 != it1->end(); it2++)
    {
      memccpy(&species_strings[str_i], it2->Name().c_str(), 0, species_max_strlen);
      str_i+=species_max_strlen;
    }
      
  if ((retval = nc_def_dim (ncid, "species_count_nelem", gal.species.nelem(), &species_nelem_ncdims[0])))
    ncerror (retval, "nc_def_dim");
  if ((retval = nc_def_var (ncid, "species_count", NC_LONG, 1, &species_nelem_ncdims[0],
                            &species_nelem_varid)))
    ncerror (retval, "nc_def_var");
  if ((retval = nc_def_dim (ncid, "species_strings_nelem", species_total_nelems, &species_strings_ncdims[0])))
    ncerror (retval, "nc_def_dim");
  if ((retval = nc_def_dim (ncid, "species_strings_length", species_max_strlen, &species_strings_ncdims[1])))
    ncerror (retval, "nc_def_dim");
  if ((retval = nc_def_var (ncid, "species_strings", NC_CHAR, 2, &species_strings_ncdims[0], &species_strings_varid)))
    ncerror (retval, "nc_def_var");
  
  
  
  // Define nonlinear_species dimensions and variable
  int nonlinear_species_ncdims[1], nonlinear_species_varid;
  if (gal.nonlinear_species.nelem())
  {
    if ((retval = nc_def_dim (ncid, "nonlinear_species_nelem", gal.nonlinear_species.nelem(),
                              &nonlinear_species_ncdims[0])))
      ncerror (retval, "nc_def_dim");
    if ((retval = nc_def_var (ncid, "nonlinear_species", NC_LONG, 1,
                              &nonlinear_species_ncdims[0], &nonlinear_species_varid)))
      ncerror (retval, "nc_def_var");
  }
  
  // Define f_grid dimensions and variable
  int f_grid_ncdims[1], f_grid_varid;
  if (gal.f_grid.nelem())
  {
    if ((retval = nc_def_dim (ncid, "f_grid_nelem", gal.f_grid.nelem(), &f_grid_ncdims[0])))
      ncerror (retval, "nc_def_dim");
    if ((retval = nc_def_var (ncid, "f_grid", NC_DOUBLE, 1, &f_grid_ncdims[0], &f_grid_varid)))
      ncerror (retval, "nc_def_var");
  }
  
  // Define p_grid dimensions and variable
  int p_grid_ncdims[1], p_grid_varid;
  if (gal.p_grid.nelem())
  {
    if ((retval = nc_def_dim (ncid, "p_grid_nelem", gal.p_grid.nelem(), &p_grid_ncdims[0])))
      ncerror (retval, "nc_def_dim");
    if ((retval = nc_def_var (ncid, "p_grid", NC_DOUBLE, 1, &p_grid_ncdims[0], &p_grid_varid)))
      ncerror (retval, "nc_def_var");
  }
  
  // Define vmrs_ref dimensions and variable
  int vmrs_ref_ncdims[2], vmrs_ref_varid;
  if (gal.vmrs_ref.nrows() && gal.vmrs_ref.ncols())
  {
    if ((retval = nc_def_dim (ncid, "vmrs_ref_nrows", gal.vmrs_ref.nrows(), &vmrs_ref_ncdims[0])))
      ncerror (retval, "nc_def_dim");
    if ((retval = nc_def_dim (ncid, "vmrs_ref_ncols", gal.vmrs_ref.ncols(), &vmrs_ref_ncdims[1])))
      ncerror (retval, "nc_def_dim");
    if ((retval = nc_def_var (ncid, "vmrs_ref", NC_DOUBLE, 2, &vmrs_ref_ncdims[0], &vmrs_ref_varid)))
      ncerror (retval, "nc_def_var");
  }
  
  // Define t_ref dimensions and variable
  int t_ref_ncdims[1], t_ref_varid;
  if (gal.t_ref.nelem())
  {
    if ((retval = nc_def_dim (ncid, "t_ref_nelem", gal.t_ref.nelem(), &t_ref_ncdims[0])))
      ncerror (retval, "nc_def_dim");
    if ((retval = nc_def_var (ncid, "t_ref", NC_DOUBLE, 1, &t_ref_ncdims[0], &t_ref_varid)))
      ncerror (retval, "nc_def_var");
  }
  
  // Define t_pert dimensions and variable
  int t_pert_ncdims[1], t_pert_varid;
  if (gal.t_pert.nelem())
  {
    if ((retval = nc_def_dim (ncid, "t_pert_nelem", gal.t_pert.nelem(), &t_pert_ncdims[0])))
      ncerror (retval, "nc_def_dim");
    if ((retval = nc_def_var (ncid, "t_pert", NC_DOUBLE, 1, &t_pert_ncdims[0], &t_pert_varid)))
      ncerror (retval, "nc_def_var");
  }
  
  // Define nls_pert dimensions and variable
  int nls_pert_ncdims[1], nls_pert_varid;
  if (gal.nls_pert.nelem())
  {
    if ((retval = nc_def_dim (ncid, "nls_pert_nelem", gal.nls_pert.nelem(), &nls_pert_ncdims[0])))
      ncerror (retval, "nc_def_dim");
    if ((retval = nc_def_var (ncid, "nls_pert", NC_DOUBLE, 1, &nls_pert_ncdims[0], &nls_pert_varid)))
      ncerror (retval, "nc_def_var");
  }
  
  // Define xsec dimensions and variable
  int xsec_ncdims[4], xsec_varid;
  if (gal.xsec.nbooks() && gal.xsec.npages() && gal.xsec.nrows() && gal.xsec.ncols())
  {
    if ((retval = nc_def_dim (ncid, "xsec_nbooks", gal.xsec.nbooks(), &xsec_ncdims[0])))
      ncerror (retval, "nc_def_dim");
    if ((retval = nc_def_dim (ncid, "xsec_npages", gal.xsec.npages(), &xsec_ncdims[1])))
      ncerror (retval, "nc_def_dim");
    if ((retval = nc_def_dim (ncid, "xsec_nrows", gal.xsec.nrows(), &xsec_ncdims[2])))
      ncerror (retval, "nc_def_dim");
    if ((retval = nc_def_dim (ncid, "xsec_ncols", gal.xsec.ncols(), &xsec_ncdims[3])))
      ncerror (retval, "nc_def_dim");
    if ((retval = nc_def_var (ncid, "xsec", NC_DOUBLE, 4, &xsec_ncdims[0], &xsec_varid)))
      ncerror (retval, "nc_def_var");
  }
  
  if ((retval = nc_enddef (ncid))) ncerror (retval, "nc_enddef");
  
  // Write variables
  if (species_total_nelems)
  {
    if ((retval = nc_put_var_long (ncid, species_nelem_varid, per_species_nelem)))
      ncerror (retval, "nc_put_var");
    if ((retval = nc_put_var_text (ncid, species_strings_varid, species_strings)))
      ncerror (retval, "nc_put_var");
  }

  delete [] per_species_nelem;
  delete [] species_strings;
  
  if (gal.nonlinear_species.nelem())
  {
    Index *ind_arr = new Index[10];
    Index i = 0;
    for (ArrayOfIndex::const_iterator it = gal.nonlinear_species.begin();
         it != gal.nonlinear_species.end(); it++, i++)
    {
      ind_arr[i] = *it;
    }
    
    if ((retval = nc_put_var_long (ncid, nonlinear_species_varid, ind_arr)))
      ncerror (retval, "nc_put_var");
    
    delete [] ind_arr;
  }
  
  if (gal.f_grid.nelem())
  {
    if ((retval = nc_put_var_double (ncid, f_grid_varid, gal.f_grid.get_c_array())))
      ncerror (retval, "nc_put_var");
  }
  
  if (gal.p_grid.nelem())
  {
    if ((retval = nc_put_var_double (ncid, p_grid_varid, gal.p_grid.get_c_array())))
      ncerror (retval, "nc_put_var");
  }
  
  if (gal.vmrs_ref.nrows() && gal.vmrs_ref.ncols())
  {
    if ((retval = nc_put_var_double (ncid, vmrs_ref_varid, gal.vmrs_ref.get_c_array())))
      ncerror (retval, "nc_put_var");
  }
  
  if (gal.t_ref.nelem())
  {
    if ((retval = nc_put_var_double (ncid, t_ref_varid, gal.t_ref.get_c_array())))
      ncerror (retval, "nc_put_var");
  }
  
  if (gal.t_pert.nelem())
  {
    if ((retval = nc_put_var_double (ncid, t_pert_varid, gal.t_pert.get_c_array())))
      ncerror (retval, "nc_put_var");
  }
  
  if (gal.nls_pert.nelem())
  {
    if ((retval = nc_put_var_double (ncid, nls_pert_varid, gal.nls_pert.get_c_array())))
      ncerror (retval, "nc_put_var");
  }
  
  if (gal.xsec.nbooks() && gal.xsec.npages() && gal.xsec.nrows() && gal.xsec.ncols())
  {
    if ((retval = nc_put_var_double (ncid, xsec_varid, gal.xsec.get_c_array())))
      ncerror (retval, "nc_put_var");
  }
}


////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

#define TMPL_NC_READ_WRITE_FILE_DUMMY(what) \
  void nc_write_to_file (const int, const what&) \
  { \
    throw runtime_error ("NetCDF support not yet implemented for this type!"); \
  } \
  void nc_read_from_file (const int, what&) \
  { \
    throw runtime_error ("NetCDF support not yet implemented for this type!"); \
  }

//=== Compound Types =======================================================

TMPL_NC_READ_WRITE_FILE_DUMMY( Agenda )
TMPL_NC_READ_WRITE_FILE_DUMMY( GriddedField1 )
TMPL_NC_READ_WRITE_FILE_DUMMY( GriddedField2 )
TMPL_NC_READ_WRITE_FILE_DUMMY( GriddedField3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( GriddedField4 )
TMPL_NC_READ_WRITE_FILE_DUMMY( GridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( IsotopeRecord )
TMPL_NC_READ_WRITE_FILE_DUMMY( MCAntenna )
TMPL_NC_READ_WRITE_FILE_DUMMY( Ppath )
TMPL_NC_READ_WRITE_FILE_DUMMY( RetrievalQuantity )
TMPL_NC_READ_WRITE_FILE_DUMMY( SLIData2 )
TMPL_NC_READ_WRITE_FILE_DUMMY( SingleScatteringData )
TMPL_NC_READ_WRITE_FILE_DUMMY( SpeciesRecord )
TMPL_NC_READ_WRITE_FILE_DUMMY( SpeciesTag )

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE_DUMMY

#endif /* ENABLE_NETCDF */

