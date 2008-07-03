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
  //! Construct an empty GField
  /*!
    The constructor for GField is protected because it is only used internally
    by the derived classed.
  */
  GField() : dim(0),
             mname(),
             mgridtypes(),
             mgridnames(),
             mstringgrids(),
             mnumericgrids()
  { /* Nothing to do here */ };

  //! Construct a GField
  /*!
    Constructs a GField with the given dimension and name.

    The constructor for GField is protected because it is only used internally
    by the derived classed.

    \param[in] d Dimension.
    \param[in] s Name.
  */
  GField(const Index d, const String s) : dim(d),
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

  void copy_grids (const GField& gf);

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

  //! GField destructor
  virtual ~GField() { }

  friend ostream& operator<<(ostream& os, const GField& gf);
};


class GField1: public GField, public Vector
{
public:
  //! Construct an empty GField1
  GField1() : GField(1, "") {};
  //! Construct an empty GField1 with the given name
  /*! \param[in] s Name. */
  GField1(const String s) : GField(1, s) {};

  //! Consistency check.
  /*!
    Check if the sizes of the grids match the data dimension.

    \return True if sizes match.
  */
  virtual bool checksize() const
    {
      return (nelem() == get_grid_size(0));
    }

  //! Make this GField1 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GField1& gf)
    {
      Vector::resize(gf.get_grid_size(0));
    }

  //! Resize the data vector.
  /*! \see Vector::resize */
  void resize(Index n)
    {
      Vector::resize(n);
    }

  friend ostream& operator<<(ostream& os, const GField1& gf);
};


class GField2: public GField, public Matrix
{
public:
  //! Construct an empty GField2
  GField2() : GField(2, "") {};
  //! Construct an empty GField2 with the given name
  /*! \param[in] s Name. */
  GField2(const String s) : GField(2, s) {};

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

  //! Make this GField2 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GField2& gf)
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

  friend ostream& operator<<(ostream& os, const GField2& gf);
};


class GField3: public GField, public Tensor3
{
public:
  //! Construct an empty GField3
  GField3() : GField(3, "") {};
  //! Construct an empty GField3 with the given name
  /*! \param[in] s Name. */
  GField3(const String s) : GField(3, s) {};

  GField3& operator=(Numeric n)
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

  //! Make this GField3 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GField3& gf)
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

  friend ostream& operator<<(ostream& os, const GField3& gf);
};


class GField4: public GField, public Tensor4
{
public:
  //! Construct an empty GField4
  GField4() : GField(4, "") {};
  //! Construct an empty GField4 with the given name
  /*! \param[in] s Name. */
  GField4(const String s) : GField(4, s) {};

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

  //! Make this GField4 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GField4& gf)
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

  friend ostream& operator<<(ostream& os, const GField4& gf);
};


/********** Output operators **********/

ostream& operator<<(ostream& os, const GField& gf);
ostream& operator<<(ostream& os, const GField1& gf);
ostream& operator<<(ostream& os, const GField2& gf);
ostream& operator<<(ostream& os, const GField3& gf);
ostream& operator<<(ostream& os, const GField4& gf);

/**************************************/


typedef Array<GField1> ArrayOfGField1;
typedef Array<GField2> ArrayOfGField2;
typedef Array<GField3> ArrayOfGField3;
typedef Array<GField4> ArrayOfGField4;
typedef Array< Array<GField1> > ArrayOfArrayOfGField1;
typedef Array< Array<GField3> > ArrayOfArrayOfGField3;

#endif

