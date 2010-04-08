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

/*! Abstract base class for gridded fields. */
class GriddedField
{
private:
  Index dim;
  String mname;
  Array<GridType> mgridtypes;
  ArrayOfString mgridnames;
  Array<ArrayOfString> mstringgrids;
  ArrayOfVector mnumericgrids;

protected:
  //! Construct an empty GriddedField
  /*!
    The constructor for GriddedField is protected because it is only used internally
    by the derived classed.
  */
  GriddedField() : dim(0),
             mname(),
             mgridtypes(),
             mgridnames(),
             mstringgrids(),
             mnumericgrids()
  { /* Nothing to do here */ };

  //! Construct a GriddedField
  /*!
    Constructs a GriddedField with the given dimension and name.

    The constructor for GriddedField is protected because it is only used internally
    by the derived classed.

    \param[in] d Dimension.
    \param[in] s Name.
  */
  GriddedField(const Index d, const String s) : dim(d),
                                          mname(s),
                                          mgridtypes(d, GRIDTYPE_NUMERIC),
                                          mgridnames(d),
                                          mstringgrids(d),
                                          mnumericgrids(d)
  { /* Nothing to do here */ }
                                          

public:
  //! Get the dimension of this gridded field.
  /*! \return Dimension. */
  Index get_dim () const { return dim; }

  void copy_grids (const GriddedField& gf);

  //! Get grid name.
  /*!
     Returns the name of the grid with index i.

     \param[in] i Grid index.
     \return      Grid name.
  */
  const String& get_grid_name (Index i) const { return mgridnames[i]; }

  Index get_grid_size (const Index i) const;

  //! Get grid type.
  /*!
     Returns the type of the grid with index i.

     \param[in] i Grid index.
     \return      Grid type.
  */
  GridType get_grid_type (Index i) const { return mgridtypes[i]; }

  ConstVectorView get_numeric_grid (Index i) const;

  VectorView get_numeric_grid (Index i);

  const ArrayOfString& get_string_grid (Index i) const;

  ArrayOfString& get_string_grid (Index i);

  //! Get the name of this gridded field.
  /*! \return Gridded field name. */
  const String& get_name () const { return mname; }

  void set_grid (Index i, const Vector& g);

  void set_grid (Index i, const ArrayOfString& g);

  //! Set grid name.
  /*!
    Sets the name with the given index.

    \param[in] i Grid index.
    \param[in] s Grid name.
  */
  void set_gridname (Index i, const String& s)
    {
      assert (i < dim);
      mgridnames[i] = s;
    };

  //! Set name of this gridded field.
  /*! \param[in] s Gridded field name. */
  void set_name (const String& s) { mname = s; }

  //! Consistency check.
  /*!
    Check if the sizes of the grids match the data dimension.

    This function must be overwritten by the derived classes.

    \return True if sizes match.
  */
  virtual bool checksize() const { return false; };

  //! GriddedField destructor
  virtual ~GriddedField() { }

  friend ostream& operator<<(ostream& os, const GriddedField& gf);
};


class GriddedField1: public GriddedField, public Vector
{
public:
  //! Construct an empty GriddedField1
  GriddedField1() : GriddedField(1, "") {};
  //! Construct an empty GriddedField1 with the given name
  /*! \param[in] s Name. */
  GriddedField1(const String s) : GriddedField(1, s) {};

  //! Consistency check.
  /*!
    Check if the sizes of the grids match the data dimension.

    \return True if sizes match.
  */
  virtual bool checksize() const
    {
      return (nelem() == get_grid_size(0));
    }

  //! Make this GriddedField1 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GriddedField1& gf)
    {
      Vector::resize(gf.get_grid_size(0));
    }

  //! Resize the data vector.
  /*! \see Vector::resize */
  void resize(Index n)
    {
      Vector::resize(n);
    }

  friend ostream& operator<<(ostream& os, const GriddedField1& gf);
};


class GriddedField2: public GriddedField, public Matrix
{
public:
  //! Construct an empty GriddedField2
  GriddedField2() : GriddedField(2, "") {};
  //! Construct an empty GriddedField2 with the given name
  /*! \param[in] s Name. */
  GriddedField2(const String s) : GriddedField(2, s) {};

  //! Consistency check.
  /*!
    Check if the sizes of the grids match the data dimension.

    \return True if sizes match.
  */
  virtual bool checksize() const
    {
      return (ncols() == get_grid_size(1)
              && nrows() == get_grid_size(0));
    }

  //! Make this GriddedField2 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GriddedField2& gf)
    {
      Matrix::resize(gf.get_grid_size(0),
                     gf.get_grid_size(1));
    }

  //! Resize the data matrix.
  /*! \see Matrix::resize */
  void resize(Index r, Index c)
    {
      Matrix::resize(r, c);
    }

  friend ostream& operator<<(ostream& os, const GriddedField2& gf);
};


class GriddedField3: public GriddedField, public Tensor3
{
public:
  //! Construct an empty GriddedField3
  GriddedField3() : GriddedField(3, "") {};
  //! Construct an empty GriddedField3 with the given name
  /*! \param[in] s Name. */
  GriddedField3(const String s) : GriddedField(3, s) {};

  GriddedField3& operator=(Numeric n)
    {
      Tensor3::operator=(n);

      return *this;
    }

  //! Consistency check.
  /*!
    Check if the sizes of the grids match the data dimension.

    \return True if sizes match.
  */
  virtual bool checksize() const
    {
      return (ncols() == get_grid_size(2)
              && nrows() == get_grid_size(1)
              && npages() == get_grid_size(0));
    }

  //! Make this GriddedField3 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GriddedField3& gf)
    {
      Tensor3::resize(gf.get_grid_size(0),
                      gf.get_grid_size(1),
                      gf.get_grid_size(2));
    }

  //! Resize the data tensor.
  /*! \see Tensor3::resize */
  void resize(Index p, Index r, Index c)
    {
      Tensor3::resize(p, r, c);
    }

  friend ostream& operator<<(ostream& os, const GriddedField3& gf);
};


class GriddedField4: public GriddedField, public Tensor4
{
public:
  //! Construct an empty GriddedField4
  GriddedField4() : GriddedField(4, "") {};
  //! Construct an empty GriddedField4 with the given name
  /*! \param[in] s Name. */
  GriddedField4(const String s) : GriddedField(4, s) {};

  //! Consistency check.
  /*!
    Check if the sizes of the grids match the data dimension.

    \return True if sizes match.
  */
  virtual bool checksize() const
    {
      return (ncols() == get_grid_size(3)
              && nrows() == get_grid_size(2)
              && npages() == get_grid_size(1)
              && nbooks() == get_grid_size(0));
    }

  //! Make this GriddedField4 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GriddedField4& gf)
    {
      Tensor4::resize(gf.get_grid_size(0),
                      gf.get_grid_size(1),
                      gf.get_grid_size(2),
                      gf.get_grid_size(3));
    }

  //! Resize the data tensor.
  /*! \see Tensor4::resize */
  void resize(Index b, Index p, Index r, Index c)
    {
      Tensor4::resize(b, p, r, c);
    }

  friend ostream& operator<<(ostream& os, const GriddedField4& gf);
};


/********** Output operators **********/

ostream& operator<<(ostream& os, const GriddedField& gf);
ostream& operator<<(ostream& os, const GriddedField1& gf);
ostream& operator<<(ostream& os, const GriddedField2& gf);
ostream& operator<<(ostream& os, const GriddedField3& gf);
ostream& operator<<(ostream& os, const GriddedField4& gf);

/**************************************/


typedef Array<GriddedField1> ArrayOfGriddedField1;
typedef Array<GriddedField2> ArrayOfGriddedField2;
typedef Array<GriddedField3> ArrayOfGriddedField3;
typedef Array<GriddedField4> ArrayOfGriddedField4;
typedef Array< Array<GriddedField1> > ArrayOfArrayOfGriddedField1;
typedef Array< Array<GriddedField3> > ArrayOfArrayOfGriddedField3;

#endif

