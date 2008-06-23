/* Copyright (C) 2002-2008 Claudia Emde <claudia.emde@dlr.de>

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
  \file   gridded_fields.h
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   Wed Jun 26 17:48:29 2002

  \brief  Reading routines for gridded fields.

  This file contains reading routines for gridded fields. Gridded fields are
  needed to store moredimesional data together with the corresponding grids
  in the same variable. The datatype og a gridded field is always an array
  of tensors.

  For further description see AUG.

*/

#ifndef gridded_fields_h
#define gridded_fields_h

#include "matpackIV.h"
#include "array.h"
#include "mystring.h"

/*! Enumeration containing the possible grid types for gridded fields */
typedef enum {
  GRIDTYPE_NUMERIC,
  GRIDTYPE_STRING
} GridType;

typedef Array<GridType> ArrayOfGridType;

/*! Abstract base class for GriddedFields */
class GField
{
private:
  Index dim;
  String mname;
  Array<GridType> mgridtypes;
  ArrayOfString mgridnames;
  Array<ArrayOfString> mstringgrids;
  ArrayOfVector mnumericgrids;

protected:
  GField() : dim(0),
             mname(),
             mgridtypes(),
             mgridnames(),
             mstringgrids(),
             mnumericgrids()
  { /* Nothing to do here */ };

  GField(const Index d, const String s) : dim(d),
                                          mname(s),
                                          mgridtypes(d, GRIDTYPE_NUMERIC),
                                          mgridnames(d),
                                          mstringgrids(d),
                                          mnumericgrids(d)
  { /* Nothing to do here */ }
                                          

public:
  Index get_dim () const { return dim; }

  const String& get_gridname (Index i) const { return mgridnames[i]; }

  Index get_grid_size (const Index i) const;

  GridType get_gridtype (Index i) const { return mgridtypes[i]; }

  ConstVectorView get_numeric_grid (Index i) const;

  const ArrayOfString& get_string_grid (Index i) const;

  const String& get_name () const { return mname; }

  void set_grid (Index i, const Vector& g);

  void set_grid (Index i, const ArrayOfString& g);

  void set_gridname (Index i, const String& s)
    {
      assert (i < dim);
      mgridnames[i] = s;
    };

  void set_name (const String& s) { mname = s; }

  virtual bool checksize() const { return false; };

  virtual ~GField() { }

  friend ostream& operator<<(ostream& os, const GField& gf);
};


class GField1: public GField, public Vector
{
public:
  GField1() : GField(1, "") {};
  GField1(const String name) : GField(1, name) {};

  virtual bool checksize() const
    {
      return (nelem() == get_grid_size(0));
    }

  friend ostream& operator<<(ostream& os, const GField1& gf);
};


class GField2: public GField, public Matrix
{
public:
  GField2() : GField(2, "") {};
  GField2(const String name) : GField(2, name) {};

  virtual bool checksize() const
    {
      return (ncols() == get_grid_size(1)
              && nrows() == get_grid_size(0));
    }

  friend ostream& operator<<(ostream& os, const GField2& gf);
};


class GField3: public GField, public Tensor3
{
public:
  GField3() : GField(3, "") {};
  GField3(const String name) : GField(3, name) {};

  virtual bool checksize() const
    {
      return (ncols() == get_grid_size(2)
              && nrows() == get_grid_size(1)
              && npages() == get_grid_size(0));
    }

  friend ostream& operator<<(ostream& os, const GField3& gf);
};


class GField4: public GField, public Tensor4
{
public:
  GField4() : GField(4, "") {};
  GField4(const String name) : GField(4, name) {};

  virtual bool checksize() const
    {
      return (ncols() == get_grid_size(3)
              && nrows() == get_grid_size(2)
              && npages() == get_grid_size(1)
              && nbooks() == get_grid_size(0));
    }

  friend ostream& operator<<(ostream& os, const GField4& gf);
};


ostream& operator<<(ostream& os, const GField& gf);
ostream& operator<<(ostream& os, const GField1& gf);
ostream& operator<<(ostream& os, const GField2& gf);
ostream& operator<<(ostream& os, const GField3& gf);
ostream& operator<<(ostream& os, const GField4& gf);



//! Contains a GriddedField3.
/*!
 A GriddedField3 consists of a pressure grid vector, a latitude vector, a
 longitude vector and a Tensor3 for the data itself.

 \author Oliver Lemke
 \date   2003-06-19
 */
typedef struct {
      //! Pressure grid.
      Vector p_grid;
      //! Latitude grid.
      Vector lat_grid;
      //! Longitude grid.
      Vector lon_grid;
      //! Data.
      Tensor3 data;
} GriddedField3;

//! Contains a GriddedField4.
/*!
  A GriddedField4 can be used to store fields of several scalar
  variables on the same p-lat-lon grid. A typical example is to store
  several atmospheric variables from a circulation model, or a
  radiosonde profile.

  A GriddedField4 consists of a field name array, a pressure grid
  vector, a latitude vector, a longitude vector and a Tensor4 for the
  data itself. The dimensions in the data tensor are also sorted in
  the above order.
  
  \author Stefan Buehler
  \date   2007-07-25
 */
typedef struct {
      //! Field names
      ArrayOfString field_names;
      //! Pressure grid.
      Vector p_grid;
      //! Latitude grid.
      Vector lat_grid;
      //! Longitude grid.
      Vector lon_grid;
      //! Data.
      Tensor4 data;
} GriddedField4;

typedef Array< Array<GriddedField3> > ArrayOfArrayOfGriddedField3;
typedef Array<GriddedField3> ArrayOfGriddedField3;
typedef Array<GriddedField4> ArrayOfGriddedField4;

ostream& operator<< (ostream& os, const GriddedField3& gfield3);
ostream& operator<< (ostream& os, const ArrayOfGriddedField3& agfield3);

ostream& operator<< (ostream& os, const GriddedField4& gfield4);
ostream& operator<< (ostream& os, const ArrayOfGriddedField4& agfield4);

#endif
