/* Copyright (C) 2008 Oliver Lemke <olemke@core-dump.info>

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
  \file   nc_io.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-09-12

  \brief This file contains basic functions to handle NetCDF data files.

*/

#include "arts.h"

#ifdef ENABLE_NETCDF

#include "nc_io.h"
#include "nc_io_types.h"
#include "file.h"
#include "messages.h"
#include "exceptions.h"


////////////////////////////////////////////////////////////////////////////
//   Default file name
////////////////////////////////////////////////////////////////////////////

//! Gives the default filename for the NetCDF formats.
/*!
  The default name is only used if the filename is empty.

  \param filename filename
  \param varname variable name
*/
void
filename_nc(String&  filename, const String&  varname )
{
  if ("" == filename)
    {
      extern const String out_basename;
      filename = out_basename + "." + varname + ".nc";
    }
}



//! Gives the default filename, with file index, for the NetCDF formats.
/*!
  The default name is only used if the filename is empty.

  \param[out] filename   filename
  \param[in]  file_index Index appended to the filename
  \param[in]  varname    variable name
*/
void
filename_nc_with_index(String&  filename,
                       const Index&   file_index,
                       const String&  varname )
{
  if ("" == filename)
    {
      extern const String out_basename;
      ostringstream os;
      os << out_basename << "." << varname << "." << file_index << ".nc";
      filename = os.str();
    }
  else
    {
      ostringstream os;
      os << filename << "." << file_index << ".nc";
      filename = os.str();
    }
}


template<typename T> void
nc_read_from_file(const String& filename,
                  T& type,
                  const Verbosity& verbosity)
{
  CREATE_OUT2
  
  String efilename = expand_path(filename);
  
  out2 << "  Reading " << efilename << '\n';

  int ncid;
  if (nc_open (efilename.c_str(), NC_NOWRITE, &ncid))
    {
      ostringstream os;
      os << "Error reading file: " << efilename << endl;
      throw runtime_error (os.str());
    }

  nc_read_from_file (ncid, type);

  nc_close (ncid);
}


template<typename T> void
nc_write_to_file(const String& filename,
                 const T& type,
                 const Verbosity& verbosity)
{
  CREATE_OUT2
  
  String efilename = expand_path(filename);
  
  out2 << "  Writing " << efilename << '\n';

  int ncid;
  if (nc_create (efilename.c_str(), NC_CLOBBER, &ncid))
    {
      ostringstream os;
      os << "Error writing file: " << efilename << endl;
      throw runtime_error (os.str());
    }

  nc_write_to_file (ncid, type);

  nc_close (ncid);
}


void nc_get_data_int(const int ncid, const String &name, int *data)
{
  int retval, varid;
  if ((retval = nc_inq_varid (ncid, name.c_str(), &varid)))
    ncerror (retval, "nc_inq_varid("+name+")");
  if ((retval = nc_get_var_int (ncid, varid, data)))
    ncerror (retval, "nc_get_var("+name+")");
}


void nc_get_data_long(const int ncid, const String &name, long *data)
{
  int retval, varid;
  if ((retval = nc_inq_varid (ncid, name.c_str(), &varid)))
    ncerror (retval, "nc_inq_varid("+name+")");
  if ((retval = nc_get_var_long (ncid, varid, data)))
    ncerror (retval, "nc_get_var("+name+")");
}


void nc_get_data_double(const int ncid, const String &name, Numeric *data)
{
  int retval, varid;
  if ((retval = nc_inq_varid (ncid, name.c_str(), &varid)))
    ncerror (retval, "nc_inq_varid("+name+")");
  if ((retval = nc_get_var_double (ncid, varid, data)))
    ncerror (retval, "nc_get_var("+name+")");
}


void nc_get_dataa_double(const int ncid, const String &name,
                         size_t start, size_t count, Numeric *data)
{
  int retval, varid;
  if ((retval = nc_inq_varid (ncid, name.c_str(), &varid)))
    ncerror (retval, "nc_inq_varid("+name+")");
  if ((retval = nc_get_vara_double (ncid, varid, &start, &count, data)))
    ncerror (retval, "nc_get_var("+name+")");
}


void nc_get_data_text(const int ncid, const String &name, char *data)
{
  int retval, varid;
  if ((retval = nc_inq_varid (ncid, name.c_str(), &varid)))
    ncerror (retval, "nc_inq_varid("+name+")");
  if ((retval = nc_get_var_text (ncid, varid, data)))
    ncerror (retval, "nc_get_var("+name+")");
}


Index nc_get_dim(const int ncid, const String &name, const bool noerror)
{
  int retval, dimid;
  size_t ndim;
  if ((retval = nc_inq_dimid (ncid, name.c_str(), &dimid)))
    if (!noerror) ncerror (retval, "nc_inq_ndims("+name+")"); else return 0;
  if ((retval = nc_inq_dimlen (ncid, dimid, &ndim)))
    if (!noerror) ncerror (retval, "nc_inq_dimlen("+name+")"); else return 0;

  return (Index)ndim;
}


void nc_get_data_ArrayOfIndex(const int ncid, const String &name, ArrayOfIndex &aoi,
                              const bool noerror)
{
  Index nelem = nc_get_dim (ncid, name+"_nelem", noerror);
  aoi.resize(nelem);
  if (nelem)
  {
    Index *ind_arr = new Index[nelem];
    nc_get_data_long (ncid, name, ind_arr);
    Index i = 0;
    for (ArrayOfIndex::iterator it = aoi.begin(); it != aoi.end(); it++, i++)
    {
      *it = ind_arr[i];
    }
  }
}


void nc_get_data_ArrayOfArrayOfSpeciesTag(const int ncid, const String &name,
                                          ArrayOfArrayOfSpeciesTag &aast,
                                          const bool noerror)
{
  ArrayOfIndex species_count;
  nc_get_data_ArrayOfIndex(ncid, name+"_count", species_count, noerror);
  aast.resize(species_count.nelem());
  if (species_count.nelem())
  {
    Index species_strings_nelem = nc_get_dim (ncid, name+"_strings_nelem", noerror);
    Index species_strings_length = nc_get_dim (ncid, name+"_strings_length", noerror);
    char* species_strings = new char[species_strings_nelem*species_strings_length];
    if (species_count.nelem()) nc_get_data_text(ncid, name+"_strings", species_strings);
    
    Index si = 0;
    for(Index i=0; i < species_count.nelem(); i++)
    {
      aast[i].resize(0);
      for(Index j=0; j < species_count[i]; j++)
      {
        aast[i].push_back(SpeciesTag(&species_strings[si]));
        si += species_strings_length;
      }
    }
    
    delete [] species_strings;
  }
}


void nc_get_data_Vector(const int ncid, const String &name, Vector &v, const bool noerror)
{
  Index nelem = nc_get_dim (ncid, name+"_nelem", noerror);
  v.resize(nelem);
  if (nelem) nc_get_data_double (ncid, name, v.get_c_array());
}


void nc_get_data_Matrix(const int ncid, const String &name, Matrix &m, const bool noerror)
{
  Index nrows = nc_get_dim (ncid, name+"_nrows", noerror);
  Index ncols = nc_get_dim (ncid, name+"_ncols", noerror);
  m.resize(nrows, ncols);
  if (nrows && ncols) nc_get_data_double (ncid, name, m.get_c_array());
}


void nc_get_data_Tensor4(const int ncid, const String &name, Tensor4 &t, const bool noerror)
{
  Index nbooks = nc_get_dim (ncid, name+"_nbooks", noerror);
  Index npages = nc_get_dim (ncid, name+"_npages", noerror);
  Index nrows = nc_get_dim (ncid, name+"_nrows", noerror);
  Index ncols = nc_get_dim (ncid, name+"_ncols", noerror);
  t.resize(nbooks, npages, nrows, ncols);
  if (nbooks && npages && nrows && ncols) nc_get_data_double (ncid, name, t.get_c_array());
}


void ncerror (const int e, const String s)
{
  ostringstream os;
  os << "NetCDF error: " << s << ", " << e;
  throw runtime_error (os.str());
}


// We can't do the instantiation at the beginning of this file, because the
// implementation of nc_write_to_file and nc_read_from_file have to be known.

#include "nc_io_instantiation.h"

#endif /* ENABLE_NETCDF */

