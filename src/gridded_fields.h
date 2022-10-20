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
#include <utility>

#include "array.h"
#include "matpack_arrays.h"
#include "matpack_data.h"
#include "artstime.h"
#include "mystring.h"

namespace GriddedFieldGrids {
  /** Global constant, Index of the frequency grid in GriddedField1.
    \author Patrick Eriksson
    \date   2008-07-02
*/
inline constexpr Index GFIELD1_F_GRID = 0;

/** Global constant, Index of the pressure grid in GriddedField3.
    \author Oliver Lemke
    \date   2008-06-24
*/
inline constexpr Index GFIELD3_P_GRID = 0;

/** Global constant, Index of the latitude grid in GriddedField3.
    \author Oliver Lemke
    \date   2008-06-24
*/
inline constexpr Index GFIELD3_LAT_GRID = 1;

/** Global constant, Index of the longitude grid in GriddedField3.
    \author Oliver Lemke
    \date   2008-06-24
*/
inline constexpr Index GFIELD3_LON_GRID = 2;

/** Global constant, Index of the field names in GriddedField4.
    \author Oliver Lemke
    \date   2008-06-25
*/
inline constexpr Index GFIELD4_FIELD_NAMES = 0;

/** Global constant, Index of incidence angles in GriddedField4.
    \author Patrick Eriksson
    \date   2008-09-20
*/
inline constexpr Index GFIELD4_IA_GRID = 0;

/** Global constant, Index of the pressure grid in GriddedField4.
    \author Oliver Lemke
    \date   2008-06-25
*/
inline constexpr Index GFIELD4_P_GRID = 1;

/** Global constant, Index of the frequency grid in GriddedField4.
    \author Patrick Eriksson
    \date   2008-07-01
*/
inline constexpr Index GFIELD4_F_GRID = 1;

/** Global constant, Index of the latitude grid in GriddedField4.
    \author Oliver Lemke
    \date   2008-06-25
*/
inline constexpr Index GFIELD4_LAT_GRID = 2;

/** Global constant, Index of the zenith angle grid in GriddedField4.
    \author Patrick Eriksson
    \date   2008-07-01
*/
inline constexpr Index GFIELD4_ZA_GRID = 2;

/** Global constant, Index of the longitude grid in GriddedField4.
    \author Oliver Lemke
    \date   2008-06-25
*/
inline constexpr Index GFIELD4_LON_GRID = 3;

/** Global constant, Index of the azimuth angle grid in GriddedField4.
    \author Patrick Eriksson
    \date   2008-07-01
*/
inline constexpr Index GFIELD4_AA_GRID = 3;
}  // namespace GriddedFieldGrids

/*! Enumeration containing the possible grid types for gridded fields */
enum GridType { GRID_TYPE_NUMERIC, GRID_TYPE_STRING, GRID_TYPE_TIME };

using ArrayOfGridType = Array<GridType>;

#define CHECK_ERROR_BOILERPLATE                           \
  "size mismatch between grids and data.\n"               \
  "Note that a grid is allowed to be empty, but in the\n" \
  "data that dimension must have exactly one element.\n"

template <Index N, typename GriddedFieldType>
String metaErrorData(const GriddedFieldType& gf) {
  std::ostringstream os;
  os << "GriddedField" << N << " ";
  if (gf.get_name().length()) os << "(" << gf.get_name() << ") ";
  os << CHECK_ERROR_BOILERPLATE;
  for (Index i = 0; i < N; i++) {
    os << "Grid";
    if (gf.get_grid_name(i).nelem()) os << " (" << gf.get_grid_name(i) << ")";
    os << " = " << gf.get_grid_size(i) << "\n";
  }
  os << "Data =";
  if constexpr (N > 6)
    os << ' ' <<gf.data.nlibraries();
  if constexpr (N > 5)
    os << ' ' <<gf.data.nvitrines();
  if constexpr (N > 4)
    os << ' ' <<gf.data.nshelves();
  if constexpr (N > 3)
    os << ' ' <<gf.data.nbooks();
  if constexpr (N > 2)
    os << ' ' <<gf.data.npages();
  if constexpr (N > 1)
    os << ' ' <<gf.data.nrows() << ' ' << gf.data.ncols();
  else
    os << ' ' << gf.data.nelem();
  return os.str();
}

/*! Abstract base class for gridded fields. */
class GriddedField {
 private:
  Index dim;
  String mname;
  Array<GridType> mgridtypes;
  ArrayOfString mgridnames;
  Array<ArrayOfString> mstringgrids;
  Array<ArrayOfTime> mtimegrids;
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
        mtimegrids(),
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
  GriddedField(const Index d, String s)
      : dim(d),
        mname(std::move(s)),
        mgridtypes(d, GRID_TYPE_NUMERIC),
        mgridnames(d),
        mstringgrids(d),
        mtimegrids(d),
        mnumericgrids(d) { /* Nothing to do here */
  }

 public:
  //! Defaulted con-/de-structors and assignments so unique pointers work w/o warnings
  virtual ~GriddedField() = default;
  GriddedField(const GriddedField&) = default;
  GriddedField(GriddedField&&) = default;
  GriddedField& operator=(const GriddedField&) = default;
  GriddedField& operator=(GriddedField&&) = default;

  //! Get the dimension of this gridded field.
  /*! \return Dimension. */
  [[nodiscard]] Index get_dim() const { return dim; }

  void copy_grids(const GriddedField& gf);

  //! Get grid name.
  /*!
     Returns the name of the grid with index i.

     \param[in] i Grid index.
     \return      Grid name.
  */
  [[nodiscard]] const String& get_grid_name(Index i) const { return mgridnames[i]; }

  //! Get the size of a grid.
  /*!
   Returns the size of grid i.

   \param[in]  i  Grid index.
   \return        Grid size.
   */
  [[nodiscard]] Index get_grid_size(Index i) const {
    Index ret = 0;
    ARTS_ASSERT(i < dim);
    switch (mgridtypes[i]) {
      case GRID_TYPE_NUMERIC:
        ret = mnumericgrids[i].nelem();
        break;
      case GRID_TYPE_STRING:
        ret = mstringgrids[i].nelem();
        break;
      case GRID_TYPE_TIME:
        ret = mtimegrids[i].nelem();
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
  [[nodiscard]] GridType get_grid_type(Index i) const { return mgridtypes[i]; }

  [[nodiscard]] const Vector& get_numeric_grid(Index i) const;

  Vector& get_numeric_grid(Index i);

  [[nodiscard]] const ArrayOfString& get_string_grid(Index i) const;

  ArrayOfString& get_string_grid(Index i);

  [[nodiscard]] const ArrayOfTime& get_time_grid(Index i) const;

  ArrayOfTime& get_time_grid(Index i);

  //! Get the name of this gridded field.
  /*! \return Gridded field name. */
  [[nodiscard]] const String& get_name() const { return mname; }

  void set_grid(Index i, const Vector& g);

  void set_grid(Index i, const ArrayOfString& g);

  void set_grid(Index i, const ArrayOfTime& g);

  //! Set grid name.
  /*!
    Sets the name with the given index.

    \param[in] i Grid index.
    \param[in] s Grid name.
  */
  void set_grid_name(Index i, const String& s) {
    ARTS_ASSERT(i < dim);
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
  [[nodiscard]] virtual bool checksize() const = 0;

  //! Strict consistency check.
  /*!
   Same as GriddedField::checksize but throws runtime_error in case of error.
  */
  virtual void checksize_strict() const = 0;

  friend std::ostream& operator<<(std::ostream& os, const GriddedField& gf);
};

class GriddedField1 final : public GriddedField {
 public:
  //! Construct an empty GriddedField1
  GriddedField1() : GriddedField(1, "") {}
  //! Construct an empty GriddedField1 with the given name
  /*! \param[in] s Name. */
  GriddedField1(const String& s) : GriddedField(1, s) {}

  [[nodiscard]] bool checksize() const final {
    return (!get_grid_size(0) && data.nelem() == 1) ||
           data.nelem() == get_grid_size(0);
  }

  void checksize_strict() const final {
    ARTS_USER_ERROR_IF (!checksize(), metaErrorData<1, GriddedField1>(*this))
  }

  //! Make this GriddedField1 the same size as the given one.
  /*! \param[in] gf Source gridded field. */
  void resize(const GriddedField1& gf) { data.resize(gf.get_grid_size(0)); }

  //! Resize the data vector.
  /*! \see Vector::resize */
  void resize(Index n) { data.resize(n); }

  friend std::ostream& operator<<(std::ostream& os, const GriddedField1& gf);
  
  friend String metaErrorData<1, GriddedField1>(const GriddedField1& gf);

  Vector data;
};

class GriddedField2 final : public GriddedField {
 public:
  //! Construct an empty GriddedField2
  GriddedField2() : GriddedField(2, "") {}
  //! Construct an empty GriddedField2 with the given name
  /*! \param[in] s Name. */
  GriddedField2(const String& s) : GriddedField(2, s) {}

  [[nodiscard]] bool checksize() const final {
    return ((!get_grid_size(1) && data.ncols() == 1) ||
            data.ncols() == get_grid_size(1)) &&
           ((!get_grid_size(0) && data.nrows() == 1) ||
            data.nrows() == get_grid_size(0));
  }

  void checksize_strict() const final {
    ARTS_USER_ERROR_IF (!checksize(), metaErrorData<2, GriddedField2>(*this))
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
  
  friend String metaErrorData<2, GriddedField2>(const GriddedField2& gf);

  Matrix data;
};

class GriddedField3 final : public GriddedField {
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

  [[nodiscard]] bool checksize() const final {
    return ((!get_grid_size(2) && data.ncols() == 1) ||
            data.ncols() == get_grid_size(2)) &&
           ((!get_grid_size(1) && data.nrows() == 1) ||
            data.nrows() == get_grid_size(1)) &&
           ((!get_grid_size(0) && data.npages() == 1) ||
            data.npages() == get_grid_size(0));
  }

  void checksize_strict() const final {
    ARTS_USER_ERROR_IF (!checksize(), metaErrorData<3, GriddedField3>(*this))
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
  
  friend String metaErrorData<3, GriddedField3>(const GriddedField3& gf);

  Tensor3 data;
};

class GriddedField4 final : public GriddedField {
 public:
  //! Construct an empty GriddedField4
  GriddedField4() : GriddedField(4, "") {}
  //! Construct an empty GriddedField4 with the given name
  /*! \param[in] s Name. */
  GriddedField4(const String& s) : GriddedField(4, s) {}

  [[nodiscard]] bool checksize() const final {
    return ((!get_grid_size(3) && data.ncols() == 1) ||
            data.ncols() == get_grid_size(3)) &&
           ((!get_grid_size(2) && data.nrows() == 1) ||
            data.nrows() == get_grid_size(2)) &&
           ((!get_grid_size(1) && data.npages() == 1) ||
            data.npages() == get_grid_size(1)) &&
           ((!get_grid_size(0) && data.nbooks() == 1) ||
            data.nbooks() == get_grid_size(0));
  }

  void checksize_strict() const final {
    ARTS_USER_ERROR_IF (!checksize(), metaErrorData<4, GriddedField4>(*this))
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
  
  friend String metaErrorData<4, GriddedField4>(const GriddedField4& gf);

  Tensor4 data;
};

class GriddedField5 final : public GriddedField {
 public:
  //! Construct an empty GriddedField5
  GriddedField5() : GriddedField(5, "") {}
  //! Construct an empty GriddedField5 with the given name
  /*! \param[in] s Name. */
  GriddedField5(const String& s) : GriddedField(5, s) {}

  [[nodiscard]] bool checksize() const final {
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

  void checksize_strict() const final {
    ARTS_USER_ERROR_IF (!checksize(), metaErrorData<5, GriddedField5>(*this))
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
  
  friend String metaErrorData<5, GriddedField5>(const GriddedField5& gf);

  Tensor5 data;
};

class GriddedField6 final : public GriddedField {
 public:
  //! Construct an empty GriddedField6
  GriddedField6() : GriddedField(6, "") {}
  //! Construct an empty GriddedField6 with the given name
  /*! \param[in] s Name. */
  GriddedField6(const String& s) : GriddedField(6, s) {}

  [[nodiscard]] bool checksize() const final {
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

  void checksize_strict() const final {
    ARTS_USER_ERROR_IF (!checksize(), metaErrorData<6, GriddedField6>(*this))
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
  
  friend String metaErrorData<6, GriddedField6>(const GriddedField6& gf);

  Tensor6 data;
};

/************ Array types *************/

using ArrayOfGriddedField1 = Array<GriddedField1>;
using ArrayOfGriddedField2 = Array<GriddedField2>;
using ArrayOfGriddedField3 = Array<GriddedField3>;
using ArrayOfGriddedField4 = Array<GriddedField4>;
using ArrayOfGriddedField5 = Array<GriddedField5>;
using ArrayOfArrayOfGriddedField1 = Array<Array<GriddedField1>>;
using ArrayOfArrayOfGriddedField2 = Array<Array<GriddedField2>>;
using ArrayOfArrayOfGriddedField3 = Array<Array<GriddedField3>>;

#undef CHECK_ERROR_BOILERPLATE

#endif
