/* Copyright (C) 2008-2012 Oliver Lemke <olemke@core-dump.info>

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
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-06-24

  \brief  Implementation of gridded fields.

  This file contains the implementation for gridded fields. Gridded fields are
  needed to store moredimesional data together with the corresponding grids
  in the same variable.

  For further description see ARTS Developer Guide.
*/

#ifndef gridded_fields_h
#define gridded_fields_h

#include <stdexcept>
#include "array.h"
#include "matpackVI.h"
#include "mystring.h"

/*! Enumeration containing the possible grid types for gridded fields */
enum GridType { GRID_TYPE_NUMERIC, GRID_TYPE_STRING };

typedef Array<GridType> ArrayOfGridType;

#define CHECK_ERROR_BOILERPLATE                           \
  "size mismatch between grids and data.\n"               \
  "Note that a grid is allowed to be empty, but in the\n" \
  "data that dimension must have exactly one element.\n"

/*! Abstract base class for gridded fields. */
class GriddedField {
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
  GriddedField()
      : dim(0),
        mname(),
        mgridtypes(),
        mgridnames(),
        mstringgrids(),
        mnumericgrids() { /* Nothing to do here */
  }

  //! Construct a GriddedField
  /*!
    Constructs a GriddedField with the given dimension and name.

    The constructor for GriddedField is protected because it is only used internally
    by the derived classes.

    \param[in] d Dimension.
    \param[in] s Name.
  */
  GriddedField(const Index d, const String& s)
      : dim(d),
        mname(s),
        mgridtypes(d, GRID_TYPE_NUMERIC),
        mgridnames(d),
        mstringgrids(d),
        mnumericgrids(d) { /* Nothing to do here */
  }

 public:
  //! Get the dimension of this gridded field.
  /*! \return Dimension. */
  Index get_dim() const { return dim; }

  void copy_grids(const GriddedField& gf);

  //! Get grid name.
  /*!
     Returns the name of the grid with index i.

     \param[in] i Grid index.
     \return      Grid name.
  */
  const String& get_grid_name(Index i) const { return mgridnames[i]; }

  //! Get the size of a grid.
  /*!
   Returns the size of grid i.

   \param[in]  i  Grid index.
   \return        Grid size.
   */
  Index get_grid_size(Index i) const {
    Index ret = 0;
    assert(i < dim);
    switch (mgridtypes[i]) {
      case GRID_TYPE_NUMERIC:
        ret = mnumericgrids[i].nelem();
        break;
      case GRID_TYPE_STRING:
        ret = mstringgrids[i].nelem();
        break;
    }

    return ret;
  }

  //! Get grid type.
  /*!
     Returns the type of the grid with index i.

     \param[in] i Grid index.
     \return      Grid type.
  */
  GridType get_grid_type(Index i) const { return mgridtypes[i]; }

  ConstVectorView get_numeric_grid(Index i) const;

  VectorView get_numeric_grid(Index i);

  const ArrayOfString& get_string_grid(Index i) const;

  ArrayOfString& get_string_grid(Index i);

  //! Get the name of this gridded field.
  /*! \return Gridded field name. */
  const String& get_name() const { return mname; }

  void set_grid(Index i, const Vector& g);

  void set_grid(Index i, const ArrayOfString& g);

  //! Set grid name.
  /*!
    Sets the name with the given index.

    \param[in] i Grid index.
    \param[in] s Grid name.
  */
  void set_grid_name(Index i, const String& s) {
    assert(i < dim);
    mgridnames[i] = s;
  }

  //! Set name of this gridded field.
  /*! \param[in] s Gridded field name. */
  void set_name(const String& s) { mname = s; }

  //! Consistency check.
  /*!
    Check if the sizes of the grids match the data dimension.

    This function must be overwritten by the derived classes.

    \return True if sizes match.
  */
  virtual bool checksize() const = 0;

  //! Strict consistency check.
  /*!
   Same as GriddedField::checksize but throws runtime_error in case of error.
  */
  virtual void checksize_strict() const = 0;

  //! GriddedField virtual destructor
  virtual ~GriddedField() {}

  friend std::ostream& operator<<(std::ostream& os, const GriddedField& gf);
};

class GriddedField1 : public GriddedField {
 public:
  //! Construct an empty GriddedField1
  GriddedField1() : GriddedField(1, "") {}
  //! Construct an empty GriddedField1 with the given name
  /*! \param[in] s Name. */
  GriddedField1(const String& s) : GriddedField(1, s) {}

  virtual bool checksize() const {
    return (!get_grid_size(0) && data.nelem() == 1) ||
           data.nelem() == get_grid_size(0);
  }

  virtual void checksize_strict() const {
    if (!checksize()) {
      std::ostringstream os;
      os << "GriddedField1 ";
      if (get_name().length()) os << "(" << get_name() << ") ";
      os << CHECK_ERROR_BOILERPLATE;
      os << "Grid";
      if (get_grid_name(0).nelem()) os << " (" << get_grid_name(0) << ")";
      os << " = " << get_grid_size(0) << "\n";
      os << "Data";
      os << " = " << data.nelem();
      throw std::runtime_error(os.str());
    }
  }

  //! Make this GriddedField1 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GriddedField1& gf) { data.resize(gf.get_grid_size(0)); }

  //! Resize the data vector.
  /*! \see Vector::resize */
  void resize(Index n) { data.resize(n); }

  friend std::ostream& operator<<(std::ostream& os, const GriddedField1& gf);

  Vector data;
};

class GriddedField2 : public GriddedField {
 public:
  //! Construct an empty GriddedField2
  GriddedField2() : GriddedField(2, "") {}
  //! Construct an empty GriddedField2 with the given name
  /*! \param[in] s Name. */
  GriddedField2(const String& s) : GriddedField(2, s) {}

  virtual bool checksize() const {
    return ((!get_grid_size(1) && data.ncols() == 1) ||
            data.ncols() == get_grid_size(1)) &&
           ((!get_grid_size(0) && data.nrows() == 1) ||
            data.nrows() == get_grid_size(0));
  }

  virtual void checksize_strict() const {
    if (!checksize()) {
      std::ostringstream os;
      os << "GriddedField2 ";
      if (get_name().length()) os << "(" << get_name() << ") ";
      os << CHECK_ERROR_BOILERPLATE;
      for (Index i = 0; i < 2; i++) {
        os << "Grid " << i;
        if (get_grid_name(i).nelem()) os << " (" << get_grid_name(i) << ")";
        os << " = " << get_grid_size(i) << "\n";
      }
      os << "Data";
      os << " = " << data.nrows() << ", " << data.ncols();
      throw std::runtime_error(os.str());
    }
  }

  //! Make this GriddedField2 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GriddedField2& gf) {
    data.resize(gf.get_grid_size(0), gf.get_grid_size(1));
  }

  //! Resize the data matrix.
  /*! \see Matrix::resize */
  void resize(Index r, Index c) { data.resize(r, c); }

  friend std::ostream& operator<<(std::ostream& os, const GriddedField2& gf);

  Matrix data;
};

class GriddedField3 : public GriddedField {
 public:
  //! Construct an empty GriddedField3
  GriddedField3() : GriddedField(3, "") {}
  //! Construct an empty GriddedField3 with the given name
  /*! \param[in] s Name. */
  GriddedField3(const String& s) : GriddedField(3, s) {}

  GriddedField3& operator=(Numeric n) {
    data = n;

    return *this;
  }

  virtual bool checksize() const {
    return ((!get_grid_size(2) && data.ncols() == 1) ||
            data.ncols() == get_grid_size(2)) &&
           ((!get_grid_size(1) && data.nrows() == 1) ||
            data.nrows() == get_grid_size(1)) &&
           ((!get_grid_size(0) && data.npages() == 1) ||
            data.npages() == get_grid_size(0));
  }

  virtual void checksize_strict() const {
    if (!checksize()) {
      std::ostringstream os;
      os << "GriddedField3 ";
      if (get_name().length()) os << "(" << get_name() << ") ";
      os << CHECK_ERROR_BOILERPLATE;
      for (Index i = 0; i < 3; i++) {
        os << "Grid " << i;
        if (get_grid_name(i).nelem()) os << " (" << get_grid_name(i) << ")";
        os << " = " << get_grid_size(i) << "\n";
      }
      os << "Data";
      os << " = " << data.npages() << ", " << data.nrows() << ", "
         << data.ncols();
      throw std::runtime_error(os.str());
    }
  }

  //! Make this GriddedField3 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GriddedField3& gf) {
    data.resize(gf.get_grid_size(0), gf.get_grid_size(1), gf.get_grid_size(2));
  }

  //! Resize the data tensor.
  /*! \see Tensor3::resize */
  void resize(Index p, Index r, Index c) { data.resize(p, r, c); }

  friend std::ostream& operator<<(std::ostream& os, const GriddedField3& gf);

  Tensor3 data;
};

class GriddedField4 : public GriddedField {
 public:
  //! Construct an empty GriddedField4
  GriddedField4() : GriddedField(4, "") {}
  //! Construct an empty GriddedField4 with the given name
  /*! \param[in] s Name. */
  GriddedField4(const String& s) : GriddedField(4, s) {}

  virtual bool checksize() const {
    return ((!get_grid_size(3) && data.ncols() == 1) ||
            data.ncols() == get_grid_size(3)) &&
           ((!get_grid_size(2) && data.nrows() == 1) ||
            data.nrows() == get_grid_size(2)) &&
           ((!get_grid_size(1) && data.npages() == 1) ||
            data.npages() == get_grid_size(1)) &&
           ((!get_grid_size(0) && data.nbooks() == 1) ||
            data.nbooks() == get_grid_size(0));
  }

  virtual void checksize_strict() const {
    if (!checksize()) {
      std::ostringstream os;
      os << "GriddedField4 ";
      if (get_name().length()) os << "(" << get_name() << ") ";
      os << CHECK_ERROR_BOILERPLATE;
      for (Index i = 0; i < 4; i++) {
        os << "Grid " << i;
        if (get_grid_name(i).nelem()) os << " (" << get_grid_name(i) << ")";
        os << " = " << get_grid_size(i) << "\n";
      }
      os << "Data";
      os << " = " << data.nbooks() << ", " << data.npages() << ", ";
      os << data.nrows() << ", " << data.ncols();
      throw std::runtime_error(os.str());
    }
  }

  //! Make this GriddedField4 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GriddedField4& gf) {
    data.resize(gf.get_grid_size(0),
                gf.get_grid_size(1),
                gf.get_grid_size(2),
                gf.get_grid_size(3));
  }

  //! Resize the data tensor.
  /*! \see Tensor4::resize */
  void resize(Index b, Index p, Index r, Index c) { data.resize(b, p, r, c); }

  friend std::ostream& operator<<(std::ostream& os, const GriddedField4& gf);

  Tensor4 data;
};

class GriddedField5 : public GriddedField {
 public:
  //! Construct an empty GriddedField5
  GriddedField5() : GriddedField(5, "") {}
  //! Construct an empty GriddedField5 with the given name
  /*! \param[in] s Name. */
  GriddedField5(const String& s) : GriddedField(5, s) {}

  virtual bool checksize() const {
    return ((!get_grid_size(4) && data.ncols() == 1) ||
            data.ncols() == get_grid_size(4)) &&
           ((!get_grid_size(3) && data.nrows() == 1) ||
            data.nrows() == get_grid_size(3)) &&
           ((!get_grid_size(2) && data.npages() == 1) ||
            data.npages() == get_grid_size(2)) &&
           ((!get_grid_size(1) && data.nbooks() == 1) ||
            data.nbooks() == get_grid_size(1)) &&
           ((!get_grid_size(0) && data.nshelves() == 1) ||
            data.nshelves() == get_grid_size(0));
  }

  virtual void checksize_strict() const {
    if (!checksize()) {
      std::ostringstream os;
      os << "GriddedField5 ";
      if (get_name().length()) os << "(" << get_name() << ") ";
      os << CHECK_ERROR_BOILERPLATE;
      for (Index i = 0; i < 5; i++) {
        os << "Grid " << i;
        if (get_grid_name(i).nelem()) os << " (" << get_grid_name(i) << ")";
        os << " = " << get_grid_size(i) << "\n";
      }
      os << "Data";
      os << " = " << data.nshelves() << ", " << data.nbooks() << ", ";
      os << data.npages() << ", " << data.nrows() << ", " << data.ncols();
      throw std::runtime_error(os.str());
    }
  }

  //! Make this GriddedField5 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GriddedField5& gf) {
    data.resize(gf.get_grid_size(0),
                gf.get_grid_size(1),
                gf.get_grid_size(2),
                gf.get_grid_size(3),
                gf.get_grid_size(4));
  }

  //! Resize the data tensor.
  /*! \see Tensor5::resize */
  void resize(Index s, Index b, Index p, Index r, Index c) {
    data.resize(s, b, p, r, c);
  }

  friend std::ostream& operator<<(std::ostream& os, const GriddedField5& gf);

  Tensor5 data;
};

class GriddedField6 : public GriddedField {
 public:
  //! Construct an empty GriddedField6
  GriddedField6() : GriddedField(6, "") {}
  //! Construct an empty GriddedField6 with the given name
  /*! \param[in] s Name. */
  GriddedField6(const String& s) : GriddedField(6, s) {}

  virtual bool checksize() const {
    return ((!get_grid_size(5) && data.ncols() == 1) ||
            data.ncols() == get_grid_size(5)) &&
           ((!get_grid_size(4) && data.nrows() == 1) ||
            data.nrows() == get_grid_size(4)) &&
           ((!get_grid_size(3) && data.npages() == 1) ||
            data.npages() == get_grid_size(3)) &&
           ((!get_grid_size(2) && data.nbooks() == 1) ||
            data.nbooks() == get_grid_size(2)) &&
           ((!get_grid_size(1) && data.nshelves() == 1) ||
            data.nshelves() == get_grid_size(1)) &&
           ((!get_grid_size(0) && data.nvitrines() == 1) ||
            data.nvitrines() == get_grid_size(0));
  }

  virtual void checksize_strict() const {
    if (!checksize()) {
      std::ostringstream os;
      os << "GriddedField6 ";
      if (get_name().length()) os << "(" << get_name() << ") ";
      os << CHECK_ERROR_BOILERPLATE;
      for (Index i = 0; i < 5; i++) {
        os << "Grid " << i;
        if (get_grid_name(i).nelem()) os << " (" << get_grid_name(i) << ")";
        os << " = " << get_grid_size(i) << "\n";
      }
      os << "Data";
      os << " = " << data.nvitrines() << data.nshelves() << ", "
         << data.nbooks() << ", ";
      os << data.npages() << ", " << data.nrows() << ", " << data.ncols();
      throw std::runtime_error(os.str());
    }
  }

  //! Make this GriddedField6 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GriddedField6& gf) {
    data.resize(gf.get_grid_size(0),
                gf.get_grid_size(1),
                gf.get_grid_size(2),
                gf.get_grid_size(3),
                gf.get_grid_size(4),
                gf.get_grid_size(5));
  }

  //! Resize the data tensor.
  /*! \see Tensor6::resize */
  void resize(Index v, Index s, Index b, Index p, Index r, Index c) {
    data.resize(v, s, b, p, r, c);
  }

  friend std::ostream& operator<<(std::ostream& os, const GriddedField6& gf);

  Tensor6 data;
};

/********** Output operators **********/

std::ostream& operator<<(std::ostream& os, const GriddedField& gf);
std::ostream& operator<<(std::ostream& os, const GriddedField1& gf);
std::ostream& operator<<(std::ostream& os, const GriddedField2& gf);
std::ostream& operator<<(std::ostream& os, const GriddedField3& gf);
std::ostream& operator<<(std::ostream& os, const GriddedField4& gf);
std::ostream& operator<<(std::ostream& os, const GriddedField5& gf);
std::ostream& operator<<(std::ostream& os, const GriddedField6& gf);

/************ Array types *************/

typedef Array<GriddedField1> ArrayOfGriddedField1;
typedef Array<GriddedField2> ArrayOfGriddedField2;
typedef Array<GriddedField3> ArrayOfGriddedField3;
typedef Array<GriddedField4> ArrayOfGriddedField4;
typedef Array<GriddedField5> ArrayOfGriddedField5;
typedef Array<Array<GriddedField1> > ArrayOfArrayOfGriddedField1;
typedef Array<Array<GriddedField2> > ArrayOfArrayOfGriddedField2;
typedef Array<Array<GriddedField3> > ArrayOfArrayOfGriddedField3;

#undef CHECK_ERROR_BOILERPLATE

#endif
