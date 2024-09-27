////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   xml_io_compound_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains basic functions to handle XML data files.
*/

#include <workspace.h>

#include <algorithm>
#include <concepts>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>

#include "atm.h"
#include "cloudbox.h"
#include "configtypes.h"
#include "covariance_matrix.h"
#include "debug.h"
#include "double_imanip.h"
#include "gridded_data.h"
#include "isotopologues.h"
#include "jacobian.h"
#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"
#include "matpack_data.h"
#include "mystring.h"
#include "obsel.h"
#include "path_point.h"
#include "rtepack.h"
#include "sorted_grid.h"
#include "surf.h"
#include "xml_io.h"
#include "xml_io_base.h"
#include "xml_io_general_types.h"

////////////////////////////////////////////////////////////////////////////
//   Overloaded functions for reading/writing data from/to XML stream
////////////////////////////////////////////////////////////////////////////

//=== CIARecord ================================================

//! Reads CIARecord from XML input stream
/*!
  \param is_xml   XML Input stream
  \param irecord  CIARecord return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(std::istream& is_xml,
                          CIARecord& cr,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  String name;
  String molecule1;
  String molecule2;
  SpeciesEnum species1;
  SpeciesEnum species2;

  tag.read_from_stream(is_xml);
  tag.check_name("CIARecord");
  tag.get_attribute_value("molecule1", molecule1);
  tag.get_attribute_value("molecule2", molecule2);

  species1 = to<SpeciesEnum>(molecule1);
  species2 = to<SpeciesEnum>(molecule2);

  if (not good_enum(species1)) {
    std::ostringstream os;
    os << "Unknown species (1st molecule) in CIARecord: " << molecule1;
    throw std::runtime_error(os.str());
  }
  if (not good_enum(species2)) {
    std::ostringstream os;
    os << "Unknown species (2nd molecule) in CIARecord: " << molecule2;
    throw std::runtime_error(os.str());
  }

  cr.SetSpecies(species1, species2);

  xml_read_from_stream(is_xml, cr.mdata, pbifs);

  tag.read_from_stream(is_xml);
  tag.check_name("/CIARecord");
}

//! Writes CIARecord to XML output stream
/*!
  \param os_xml   XML Output stream
  \param irecord  CIARecord
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const CIARecord& cr,
                         bofstream* pbofs,
                         const String& name [[maybe_unused]]) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("CIARecord");
  open_tag.add_attribute("molecule1", cr.MoleculeName(0));
  open_tag.add_attribute("molecule2", cr.MoleculeName(1));
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, cr.Data(), pbofs, "");

  close_tag.set_name("/CIARecord");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== BlockMatrix =========================================================

//! Reads BlockMatrix from XML input stream
/*!
  \param is_xml    XML Input stream
  \param covmat    BlockMatrix
  \param pbifs     Pointer to binary file stream. NULL for ASCII output.
*/
void xml_read_from_stream(std::istream& is_xml,
                          BlockMatrix& bm,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  tag.read_from_stream(is_xml);
  tag.check_name("BlockMatrix");

  Index not_null;
  tag.get_attribute_value("not_null", not_null);

  if (not_null) {
    String type;
    tag.get_attribute_value("type", type);

    if (type == "Matrix") {
      Matrix m;
      xml_read_from_stream(is_xml, m, pbifs);
      bm = std::make_shared<Matrix>(std::move(m));
    } else {
      Sparse s;
      xml_read_from_stream(is_xml, s, pbifs);
      bm = std::make_shared<Sparse>(std::move(s));
    }
  } else {
    bm = BlockMatrix{};
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/BlockMatrix");
}

//! Write BlockMatrix to XML output stream
/*!
  \param os_xml    XML output stream
  \param covmat    BlockMatrix
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Unused
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const BlockMatrix& bm,
                         bofstream* pbofs,
                         const String& name [[maybe_unused]]) {
  ArtsXMLTag covmat_tag;
  ArtsXMLTag close_tag;

  covmat_tag.set_name("BlockMatrix");
  covmat_tag.add_attribute("not_null", Index(bm.not_null()));

  if (bm.not_null()) {
    const String t = bm.is_dense() ? String{"Matrix"} : String{"Sparse"};
    covmat_tag.add_attribute("type", t);
    covmat_tag.write_to_stream(os_xml);
    os_xml << '\n';

    if (bm.is_dense()) {
      xml_write_to_stream(os_xml, bm.dense(), pbofs, "");
    } else {
      xml_write_to_stream(os_xml, bm.sparse(), pbofs, "");
    }
    os_xml << '\n';
  } else {
    covmat_tag.write_to_stream(os_xml);
    os_xml << '\n';
  }

  close_tag.set_name("/BlockMatrix");
  close_tag.write_to_stream(os_xml);
}

//=== CovarianceMatrix =========================================================

//! Reads CovarianceMatrix from XML input stream
/*!
  \param is_xml    XML Input stream
  \param covmat    CovarianceMatrix
  \param pbifs     Pointer to binary file stream. NULL for ASCII output.
*/
void xml_read_from_stream(std::istream& is_xml,
                          CovarianceMatrix& covmat,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  String name, type;
  Index n_blocks, row_start, row_extent, column_start, column_extent, row_index,
      column_index, is_inverse;

  tag.read_from_stream(is_xml);
  tag.check_name("CovarianceMatrix");
  tag.get_attribute_value("n_blocks", n_blocks);

  covmat = CovarianceMatrix();
  for (Index i = 0; i < n_blocks; i++) {
    tag.read_from_stream(is_xml);
    tag.check_name("Block");

    tag.get_attribute_value("row_index", row_index);
    tag.get_attribute_value("column_index", column_index);
    tag.get_attribute_value("row_start", row_start);
    tag.get_attribute_value("row_extent", row_extent);
    tag.get_attribute_value("column_start", column_start);
    tag.get_attribute_value("column_extent", column_extent);
    tag.get_attribute_value("type", type);
    tag.get_attribute_value("is_inverse", is_inverse);

    Range row_range(row_start, row_extent);
    Range column_range(column_start, column_extent);
    if (type == "Matrix") {
      std::shared_ptr<Matrix> M =
          std::make_shared<Matrix>(row_extent, column_extent);
      xml_read_from_stream(is_xml, *M, pbifs);
      if (!is_inverse) {
        covmat.correlations_.emplace_back(
            row_range,
            column_range,
            std::make_pair(row_index, column_index),
            M);
      } else {
        covmat.inverses_.emplace_back(row_range,
                                      column_range,
                                      std::make_pair(row_index, column_index),
                                      M);
      }
    } else if (type == "Sparse") {
      std::shared_ptr<Sparse> M =
          std::make_shared<Sparse>(row_extent, column_extent);
      xml_read_from_stream(is_xml, *M, pbifs);
      if (!is_inverse) {
        covmat.correlations_.emplace_back(
            row_range,
            column_range,
            std::make_pair(row_index, column_index),
            M);
      } else {
        covmat.inverses_.emplace_back(row_range,
                                      column_range,
                                      std::make_pair(row_index, column_index),
                                      M);
      }
    }
    tag.read_from_stream(is_xml);
    tag.check_name("/Block");
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/CovarianceMatrix");
}

//! Write CovarianceMatrix to XML output stream
/*!
  \param os_xml    XML output stream
  \param covmat    CovarianceMatrix
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Unused
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const CovarianceMatrix& covmat,
                         bofstream* pbofs,
                         const String& name [[maybe_unused]]) {
  ArtsXMLTag covmat_tag;
  ArtsXMLTag close_tag;

  covmat_tag.set_name("CovarianceMatrix");
  covmat_tag.add_attribute(
      "n_blocks", Index(covmat.correlations_.size() + covmat.inverses_.size()));
  covmat_tag.write_to_stream(os_xml);
  os_xml << '\n';
  for (const Block& c : covmat.correlations_) {
    ArtsXMLTag block_tag;
    block_tag.set_name("Block");

    Index i, j;
    std::tie(i, j) = c.get_indices();
    block_tag.add_attribute("row_index", i);
    block_tag.add_attribute("column_index", j);

    Range row_range    = c.get_row_range();
    Range column_range = c.get_column_range();
    block_tag.add_attribute("row_start", row_range.offset);
    block_tag.add_attribute("row_extent", row_range.extent);
    block_tag.add_attribute("column_start", column_range.offset);
    block_tag.add_attribute("column_extent", column_range.extent);
    block_tag.add_attribute("is_inverse", Index(0));
    if (c.is_dense()) {
      block_tag.add_attribute("type", "Matrix");
      block_tag.write_to_stream(os_xml);
      os_xml << '\n';
      xml_write_to_stream(os_xml, c.get_dense(), pbofs, name);
    } else {
      block_tag.add_attribute("type", "Sparse");
      block_tag.write_to_stream(os_xml);
      os_xml << '\n';
      xml_write_to_stream(os_xml, c.get_sparse(), pbofs, name);
    }
    close_tag.set_name("/Block");
    close_tag.write_to_stream(os_xml);
    os_xml << '\n';
  }
  for (const Block& c : covmat.inverses_) {
    ArtsXMLTag block_tag;
    block_tag.set_name("Block");

    Index i, j;
    std::tie(i, j) = c.get_indices();
    block_tag.add_attribute("row_index", i);
    block_tag.add_attribute("column_index", j);

    Range row_range    = c.get_row_range();
    Range column_range = c.get_column_range();
    block_tag.add_attribute("row_start", row_range.offset);
    block_tag.add_attribute("row_extent", row_range.extent);
    block_tag.add_attribute("column_start", column_range.offset);
    block_tag.add_attribute("column_extent", column_range.extent);
    block_tag.add_attribute("is_inverse", Index(1));
    if (c.is_dense()) {
      block_tag.add_attribute("type", "Matrix");
      block_tag.write_to_stream(os_xml);
      os_xml << '\n';
      xml_write_to_stream(os_xml, c.get_dense(), pbofs, name);
    } else {
      block_tag.add_attribute("type", "Sparse");
      block_tag.write_to_stream(os_xml);
      os_xml << '\n';
      xml_write_to_stream(os_xml, c.get_sparse(), pbofs, name);
    }
    close_tag.set_name("/Block");
    close_tag.write_to_stream(os_xml);
    os_xml << '\n';
  }
  os_xml << '\n';
  close_tag.set_name("/CovarianceMatrix");
  close_tag.write_to_stream(os_xml);
}

//=== GasAbsLookup ===========================================================

//! Reads GasAbsLookup from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gal     GasAbsLookup return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(std::istream& is_xml,
                          GasAbsLookup& gal,
                          bifstream* pbifs) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("GasAbsLookup");

  xml_read_from_stream(is_xml, gal.species, pbifs);
  xml_read_from_stream(is_xml, gal.nonlinear_species, pbifs);
  xml_read_from_stream(is_xml, gal.f_grid, pbifs);
  xml_read_from_stream(is_xml, gal.p_grid, pbifs);
  xml_read_from_stream(is_xml, gal.vmrs_ref, pbifs);
  xml_read_from_stream(is_xml, gal.t_ref, pbifs);
  xml_read_from_stream(is_xml, gal.t_pert, pbifs);
  xml_read_from_stream(is_xml, gal.nls_pert, pbifs);
  xml_read_from_stream(is_xml, gal.xsec, pbifs);

  tag.read_from_stream(is_xml);
  tag.check_name("/GasAbsLookup");
}

//! Writes GasAbsLookup to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gal     GasAbsLookup
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const GasAbsLookup& gal,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("GasAbsLookup");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, gal.species, pbofs, "");
  xml_write_to_stream(os_xml, gal.nonlinear_species, pbofs, "NonlinearSpecies");
  xml_write_to_stream(os_xml, gal.f_grid, pbofs, "FrequencyGrid");
  xml_write_to_stream(os_xml, gal.p_grid, pbofs, "PressureGrid");
  xml_write_to_stream(os_xml, gal.vmrs_ref, pbofs, "ReferenceVmrProfiles");
  xml_write_to_stream(os_xml, gal.t_ref, pbofs, "ReferenceTemperatureProfile");
  xml_write_to_stream(os_xml, gal.t_pert, pbofs, "TemperaturePerturbations");
  xml_write_to_stream(
      os_xml, gal.nls_pert, pbofs, "NonlinearSpeciesVmrPerturbations");
  xml_write_to_stream(os_xml, gal.xsec, pbofs, "AbsorptionCrossSections");

  close_tag.set_name("/GasAbsLookup");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== ComplexMatrix ==========================================================

//! Reads ComplexMatrix from XML input stream
/*!
  \param is_xml  XML Input stream
  \param matrix  ComplexMatrix return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(std::istream& is_xml,
                          ComplexMatrix& matrix,
                          bifstream* pbifs) {
  XMLTag tag;
  Index nrows, ncols;

  tag.read_from_stream(is_xml);
  tag.check_name("ComplexMatrix");

  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  matrix.resize(nrows, ncols);

  if (pbifs) {
    pbifs->readDoubleArray(reinterpret_cast<Numeric*>(matrix.data_handle()),
                           nrows * ncols);
  } else {
    for (Index r = 0; r < nrows; r++) {
      for (Index c = 0; c < ncols; c++) {
        is_xml >> double_imanip() >> matrix.real()(r, c) >> matrix.imag()(r, c);
        if (is_xml.fail()) {
          std::ostringstream os;
          os << " near "
             << "\n  Row   : " << r << "\n  Column: " << c;
          xml_data_parse_error(tag, os.str());
        }
      }
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/ComplexMatrix");
}

//! Writes ComplexMatrix to XML output stream
/*!
  \param os_xml  XML Output stream
  \param matrix  ComplexMatrix
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const ComplexMatrix& matrix,
                         bofstream* pbofs,
                         const String& name) {
  XMLTag open_tag;
  XMLTag close_tag;

  open_tag.set_name("ComplexMatrix");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("nrows", matrix.nrows());
  open_tag.add_attribute("ncols", matrix.ncols());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  // Write the elements:
  for (Index r = 0; r < matrix.nrows(); ++r) {
    if (pbofs)
      *pbofs << matrix(r, 0);
    else
      os_xml << matrix(r, 0).real() << ' ' << matrix(r, 0).imag();

    for (Index c = 1; c < matrix.ncols(); ++c) {
      if (pbofs)
        *pbofs << matrix(r, c);
      else
        os_xml << " " << matrix(r, c).real() << ' ' << matrix(r, c).imag();
    }

    if (!pbofs) os_xml << '\n';
  }

  close_tag.set_name("/ComplexMatrix");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== GriddedField ===========================================================

template <Size N, typename T, typename... Grids, Size M = sizeof...(Grids)>
void xml_read_from_stream_recursive(std::istream& is_xml,
                                    matpack::gridded_data<T, Grids...>& gfield,
                                    bifstream* pbifs,
                                    XMLTag& tag) {
  if constexpr (M > N) {
    xml_read_from_stream(is_xml, gfield.template grid<N>(), pbifs);
    xml_read_from_stream_recursive<N + 1>(is_xml, gfield, pbifs, tag);
  }
}

//! Reads the grids for gridded fields from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gfield  GriddedField return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
template <typename T, typename... Grids>
void xml_read_from_stream_gf(std::istream& is_xml,
                             matpack::gridded_data<T, Grids...>& gfield,
                             bifstream* pbifs,
                             Index version) {
  using GF = matpack::gridded_data<T, Grids...>;

  XMLTag tag;

  if (version == 0) {
    const auto reader = [&](auto& grid, String& name) {
      using U = std::decay_t<decltype(grid)>;

      tag.get_attribute_value("name", name);

      if constexpr (std::same_as<U, Vector>) {
        ARTS_USER_ERROR_IF(tag.get_name() != "Vector",
                           "Must be Vector, is {}",
                           tag.get_name());
        xml_parse_from_stream(is_xml, grid, pbifs, tag);
        tag.read_from_stream(is_xml);
        tag.check_name("/Vector");
      } else if constexpr (std::same_as<U, ArrayOfString>) {
        ARTS_USER_ERROR_IF(
            tag.get_name() != "Array", "Must be Array, is {}", tag.get_name());
        String s;
        tag.get_attribute_value("type", s);
        ARTS_USER_ERROR_IF(
            s != "String", "Must be Array<String>, is Array<{}>", s);
        xml_parse_from_stream(is_xml, grid, pbifs, tag);
        tag.read_from_stream(is_xml);
        tag.check_name("/Array");
      } else {
        ARTS_USER_ERROR("Unknown grid type: {}", tag.get_name());
      }
    };

    if constexpr (constexpr Size N = 0; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 1; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 2; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 3; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 4; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 5; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 6; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 7; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 8; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 9; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    xml_read_from_stream(is_xml, gfield.data, pbifs);
  } else if (version == 1) {
    ArrayOfString gridnames;
    xml_read_from_stream(is_xml, gridnames, pbifs);
    ARTS_USER_ERROR_IF(gridnames.size() != gfield.grid_names.size(),
                       "Bad number of grid names:\n{}",
                       gridnames);
    std::ranges::move(gridnames, gfield.grid_names.begin());

    xml_read_from_stream_recursive<0>(is_xml, gfield, pbifs, tag);

    xml_read_from_stream(is_xml, gfield.data, pbifs);
  } else {
    ARTS_USER_ERROR("Unknown version: {}", version)
  }
}

#define GF_READ(GF)                                                        \
  void xml_read_from_stream(                                               \
      std::istream& is_xml, GF& gfield, bifstream* pbifs) {                \
    ArtsXMLTag tag;                                                        \
                                                                           \
    tag.read_from_stream(is_xml);                                          \
    tag.check_name(#GF);                                                   \
                                                                           \
    Index version = 0;                                                     \
    if (tag.has_attribute("version"))                                      \
      tag.get_attribute_value("version", version);                         \
    tag.get_attribute_value("name", gfield.data_name);                     \
                                                                           \
    xml_read_from_stream_gf(is_xml, gfield, pbifs, version);               \
                                                                           \
    tag.read_from_stream(is_xml);                                          \
    tag.check_name("/" #GF);                                               \
                                                                           \
    ARTS_USER_ERROR_IF(not gfield.ok(), "Bad gridded field:\n{}", gfield); \
  }

template <Size N, typename T, typename... Grids, Size M = sizeof...(Grids)>
void xml_write_to_stream_recursive(
    std::ostream& os_xml,
    const matpack::gridded_data<T, Grids...>& gfield,
    bofstream* pbofs) {
  if constexpr (M > N) {
    xml_write_to_stream(os_xml,
                        gfield.template grid<N>(),
                        pbofs,
                        gfield.template gridname<N>());
    xml_write_to_stream_recursive<N + 1>(os_xml, gfield, pbofs);
  }
}

//! Writes the grids for gridded fields to an XML input stream
/*!
  \param os_xml  XML output stream
  \param gfield  GriddedField with the grids
  \param pbofs   Pointer to binary output stream. NULL in case of ASCII file.
*/
template <typename T, typename... Grids>
void xml_write_to_stream_gf(std::ostream& os_xml,
                            const matpack::gridded_data<T, Grids...>& gfield,
                            bofstream* pbofs,
                            const String& /* name */) {
  xml_write_to_stream(
      os_xml,
      ArrayOfString{gfield.grid_names.begin(), gfield.grid_names.end()},
      pbofs,
      "GridNames");

  xml_write_to_stream_recursive<0>(os_xml, gfield, pbofs);

  xml_write_to_stream(os_xml, gfield.data, pbofs, "Data");
}

#define GF_WRITE(GF)                                       \
  void xml_write_to_stream(std::ostream& os_xml,           \
                           const GF& gfield,               \
                           bofstream* pbofs,               \
                           const String&) {                \
    ArtsXMLTag open_tag;                                   \
    ArtsXMLTag close_tag;                                  \
                                                           \
    open_tag.set_name(#GF);                                \
    open_tag.add_attribute("name", gfield.data_name);      \
    open_tag.add_attribute("version", Index{1});           \
                                                           \
    open_tag.write_to_stream(os_xml);                      \
    os_xml << '\n';                                        \
                                                           \
    xml_write_to_stream_gf(os_xml, gfield, pbofs, "Data"); \
                                                           \
    close_tag.set_name("/" #GF);                           \
    close_tag.write_to_stream(os_xml);                     \
    os_xml << '\n';                                        \
  }

#define GF_IO(GF) \
  GF_READ(GF)     \
  GF_WRITE(GF)

GF_IO(GriddedField1)
GF_IO(GriddedField2)
GF_IO(GriddedField3)
GF_IO(GriddedField4)
GF_IO(GriddedField5)
GF_IO(GriddedField6)
GF_IO(NamedGriddedField2)
GF_IO(NamedGriddedField3)
GF_IO(GriddedField1Named)
GF_IO(ComplexGriddedField2)
GF_IO(StokvecGriddedField6)

#undef GF_IO
#undef GF_READ
#undef GF_WRITE

//=== GridPos =====================================================

//! Reads GridPos from XML input stream
/*!
  \param is_xml  XML Input stream
  \param gpos    GridPos return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(std::istream& is_xml,
                          GridPos& gpos,
                          bifstream* pbifs) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("GridPos");

  xml_read_from_stream(is_xml, gpos.idx, pbifs);
  xml_read_from_stream(is_xml, gpos.fd[0], pbifs);
  xml_read_from_stream(is_xml, gpos.fd[1], pbifs);

  tag.read_from_stream(is_xml);
  tag.check_name("/GridPos");
}

//! Writes GridPos to XML output stream
/*!
  \param os_xml  XML Output stream
  \param gpos    GridPos
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const GridPos& gpos,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("GridPos");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(
      os_xml, gpos.idx, pbofs, "OriginalGridIndexBelowInterpolationPoint");
  xml_write_to_stream(
      os_xml, gpos.fd[0], pbofs, "FractionalDistanceToNextPoint_1");
  xml_write_to_stream(
      os_xml, gpos.fd[1], pbofs, "FractionalDistanceToNextPoint_2");

  close_tag.set_name("/GridPos");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== QuantumIdentifier =========================================

//! Reads QuantumIdentifier from XML input stream
/*!
  \param is_xml  XML Input stream
  \param qi      QuantumIdentifier return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(std::istream& is_xml,
                          QuantumIdentifier& qi,
                          bifstream* pbifs [[maybe_unused]]) try {
  static_assert(QuantumIdentifier::version == 1);

  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("QuantumIdentifier");

  Index version;
  if (tag.has_attribute("version")) {
    tag.get_attribute_value("version", version);
  } else {
    version = 0;
  }
  ARTS_USER_ERROR_IF(
      QuantumIdentifier::version < version,
      "The version of this quantum identifier is too new.  You need to upgrade ARTS to use it.")

  try {
    String qi_str;
    parse_xml_tag_content_as_string(is_xml, qi_str);
    qi = QuantumIdentifier(qi_str, version);
  } catch (const std::runtime_error& e) {
    std::ostringstream os;
    os << "Error reading QuantumIdentifier: "
       << "\n"
       << e.what();
    throw std::runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/QuantumIdentifier");
} catch (const std::exception& e) {
  ARTS_USER_ERROR("Error reading QuantumIdentifier: {}", e.what())
}

//! Writes QuantumIdentifier to XML output stream
/*!
  \param os_xml  XML Output stream
  \param qi      QuantumIdentifier
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const QuantumIdentifier& qi,
                         bofstream* pbofs [[maybe_unused]],
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  static_assert(QuantumIdentifier::version == 1);

  open_tag.set_name("QuantumIdentifier");
  open_tag.add_attribute("version", QuantumIdentifier::version);
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  os_xml << qi;

  close_tag.set_name("/QuantumIdentifier");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== SingleScatteringData ======================================

//! Reads SingleScatteringData from XML input stream
/*!
  \param is_xml  XML Input stream
  \param ssdata  SingleScatteringData return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(std::istream& is_xml,
                          SingleScatteringData& ssdata,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  String version;

  tag.read_from_stream(is_xml);
  tag.check_name("SingleScatteringData");
  tag.get_attribute_value("version", version);

  if (version == "3") {
    String ptype_string;
    xml_read_from_stream(is_xml, ptype_string, pbifs);
    ssdata.ptype = PTypeFromString(ptype_string);
  } else if (version == "2") {
    String ptype_string;
    xml_read_from_stream(is_xml, ptype_string, pbifs);
    ssdata.ptype = PType2FromString(ptype_string);
  } else {
    Index ptype;
    xml_read_from_stream(is_xml, ptype, pbifs);
    if (ptype != PTYPE_GENERAL && ptype != PTYPE_TOTAL_RND &&
        ptype != PTYPE_AZIMUTH_RND) {
      std::ostringstream os;
      os << "Ptype value (" << ptype << ") is wrong."
         << "It must be \n"
         << PTYPE_TOTAL_RND << " - totally randomly oriented particles,\n"
         << PTYPE_AZIMUTH_RND
         << " - azimuthally randomly oriented particles, or\n"
         << PTYPE_GENERAL << " - arbitrary oriented particles.\n";
      throw std::runtime_error(os.str());
    }
    ssdata.ptype = PType(ptype);
  }
  xml_read_from_stream(is_xml, ssdata.description, pbifs);
  xml_read_from_stream(is_xml, ssdata.f_grid, pbifs);
  xml_read_from_stream(is_xml, ssdata.T_grid, pbifs);
  xml_read_from_stream(is_xml, ssdata.za_grid, pbifs);
  /* Verify that we have a good coverage for the za grid */
  if ((ssdata.za_grid[0] > 1) ||
      ssdata.za_grid[ssdata.za_grid.size() - 1] < 179) {
    std::ostringstream os;
    os << "Missing data in xml-stream. Expected za_grid: [0, 180]. "
       << "Found za_grid: [" << ssdata.za_grid[0] << ", "
       << ssdata.za_grid[ssdata.za_grid.size() - 1] << "]";
    throw std::runtime_error(os.str());
  }
  xml_read_from_stream(is_xml, ssdata.aa_grid, pbifs);

  xml_read_from_stream(is_xml, ssdata.pha_mat_data, pbifs);
  if (ssdata.pha_mat_data.nlibraries() != ssdata.f_grid.size()) {
    throw std::runtime_error(
        "Number of frequencies in f_grid and pha_mat_data "
        "not matching!!!");
  }

  xml_read_from_stream(is_xml, ssdata.ext_mat_data, pbifs);
  xml_read_from_stream(is_xml, ssdata.abs_vec_data, pbifs);

  tag.read_from_stream(is_xml);
  tag.check_name("/SingleScatteringData");

  if (version != "3" && ssdata.ptype == PTYPE_AZIMUTH_RND) {
    ConvertAzimuthallyRandomSingleScatteringData(ssdata);
  }

  chk_scat_data(ssdata);
}

//! Writes SingleScatteringData to XML output stream
/*!
  \param os_xml  XML Output stream
  \param ssdata  SingleScatteringData
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const SingleScatteringData& ssdata,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("SingleScatteringData");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("version", "3");
  open_tag.write_to_stream(os_xml);

  os_xml << '\n';
  xml_write_to_stream(os_xml, PTypeToString(ssdata.ptype), pbofs, "");
  xml_write_to_stream(os_xml, ssdata.description, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.f_grid, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.T_grid, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.za_grid, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.aa_grid, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.pha_mat_data, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.ext_mat_data, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.abs_vec_data, pbofs, "");

  close_tag.set_name("/SingleScatteringData");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== ScatteringMetaData ======================================

//! Reads ScatteringMetaData from XML input stream
/*!
  \param is_xml  XML Input stream
  \param smdata  ScatteringMetaData return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(std::istream& is_xml,
                          ScatteringMetaData& smdata,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  String version;

  tag.read_from_stream(is_xml);
  tag.check_name("ScatteringMetaData");
  tag.get_attribute_value("version", version);

  if (version != "3") {
    std::ostringstream os;
    os << "Only ScatteringMetaData version 3 can be handled. "
       << "Versions 1 and 2 are obsolete.";
    throw std::runtime_error(os.str());
  }

  xml_read_from_stream(is_xml, smdata.description, pbifs);
  xml_read_from_stream(is_xml, smdata.source, pbifs);
  xml_read_from_stream(is_xml, smdata.refr_index, pbifs);
  xml_read_from_stream(is_xml, smdata.mass, pbifs);
  xml_read_from_stream(is_xml, smdata.diameter_max, pbifs);
  xml_read_from_stream(is_xml, smdata.diameter_volume_equ, pbifs);
  xml_read_from_stream(is_xml, smdata.diameter_area_equ_aerodynamical, pbifs);

  tag.read_from_stream(is_xml);
  tag.check_name("/ScatteringMetaData");
}

//! Writes ScatteringMetaData to XML output stream
/*!
  \param os_xml  XML Output stream
  \param smdata  ScatteringMetaData
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const ScatteringMetaData& smdata,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("ScatteringMetaData");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("version", "3");
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, smdata.description, pbofs, "");
  xml_write_to_stream(os_xml, smdata.source, pbofs, "");
  xml_write_to_stream(os_xml, smdata.refr_index, pbofs, "");
  xml_write_to_stream(os_xml, smdata.mass, pbofs, "");
  xml_write_to_stream(os_xml, smdata.diameter_max, pbofs, "");
  xml_write_to_stream(os_xml, smdata.diameter_volume_equ, pbofs, "");
  xml_write_to_stream(
      os_xml, smdata.diameter_area_equ_aerodynamical, pbofs, "");

  close_tag.set_name("/ScatteringMetaData");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== SLIData2 =====================================================
//! Reads SLIData2 from XML input stream
/*!
  \param is_xml   XML Input stream
  \param slidata  SLIData return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/

void xml_read_from_stream(std::istream& is_xml,
                          SLIData2& slidata,
                          bifstream* pbifs) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("SLIData2");

  xml_read_from_stream(is_xml, slidata.x1a, pbifs);
  xml_read_from_stream(is_xml, slidata.x2a, pbifs);
  xml_read_from_stream(is_xml, slidata.ya, pbifs);

  tag.read_from_stream(is_xml);
  tag.check_name("/SLIData2");
}

void xml_write_to_stream(std::ostream& os_xml,
                         const SLIData2& slidata,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("SLIData2");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, slidata.x1a, pbofs, "");
  xml_write_to_stream(os_xml, slidata.x2a, pbofs, "");
  xml_write_to_stream(os_xml, slidata.ya, pbofs, "");

  close_tag.set_name("/SLIData2");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== SpeciesIsotopologueRatios ===========================================

//! Reads SpeciesIsotopologueRatios from XML input stream
/*!
  \param is_xml   XML Input stream
  \param sap      SpeciesIsotopologueRatios return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(std::istream& is_xml,
                          SpeciesIsotopologueRatios& iso_rat,
                          bifstream* pbifs) {
  ARTS_USER_ERROR_IF(pbifs,
                     "No support for binary IO for (SpeciesIsotopologueRatios)")

  iso_rat = SpeciesIsotopologueRatios{};

  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("SpeciesIsotopologueRatios");

  Index nelem;
  tag.get_attribute_value("nelem", nelem);

  String name;
  Numeric val;
  for (Index n = 0; n < nelem; n++) {
    is_xml >> name >> double_imanip() >> val;

    if (is_xml.fail()) {
      std::ostringstream os;
      os << " near "
         << "\n  Element: " << n;
      xml_data_parse_error(tag, os.str());
    }

    const Index i =
        Species::find_species_index(Species::Tag(name).Isotopologue());
    ARTS_USER_ERROR_IF(
        i < 0 or i >= iso_rat.maxsize,
        "Species: {}"
        " cannot be understood as a species by your compiled version of ARTS",
        name)

    iso_rat.data[i] = val;
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/SpeciesIsotopologueRatios");
}

//! Writes SpeciesIsotopologueRatios to XML output stream
/*!
  \param os_xml   XML Output stream
  \param sap      SpeciesIsotopologueRatios
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const SpeciesIsotopologueRatios& iso_rat,
                         bofstream* pbofs,
                         const String& name)

{
  ARTS_USER_ERROR_IF(pbofs,
                     "No support for binary IO for (SpeciesIsotopologueRatios)")

  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("SpeciesIsotopologueRatios");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("nelem", iso_rat.maxsize);

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  os_xml << iso_rat << '\n';

  close_tag.set_name("/SpeciesIsotopologueRatios");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== SpeciesTag ================================================

//! Reads SpeciesTag from XML input stream
/*!
  \param is_xml  XML Input stream
  \param stag    SpeciesTag return value
*/
/*  param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
                 Ignored because SpeciesTag is always stored in ASCII format.*/
void xml_read_from_stream(std::istream& is_xml,
                          SpeciesTag& stag,
                          bifstream* /* pbifs */) {
  ArtsXMLTag tag;
  std::stringbuf strbuf;
  char dummy;

  tag.read_from_stream(is_xml);
  tag.check_name("SpeciesTag");

  // Skip whitespaces
  bool string_starts_with_quotes = true;
  do {
    is_xml >> dummy;
    switch (dummy) {
      case ' ':
      case '\"':
      case '\n':
      case '\r':
      case '\t':
        break;
      default:
        string_starts_with_quotes = false;
    }
  } while (is_xml.good() && dummy != '"' && string_starts_with_quotes);

  // Throw exception if first char after whitespaces is not a quote
  if (!string_starts_with_quotes) {
    xml_parse_error("SpeciesTag must begin with \"");
  }

  is_xml.get(strbuf, '"');
  if (is_xml.fail()) {
    xml_parse_error("SpeciesTag must end with \"");
  }

  stag = SpeciesTag(strbuf.str());

  // Ignore quote
  is_xml >> dummy;

  tag.read_from_stream(is_xml);
  tag.check_name("/SpeciesTag");
}

//! Writes SpeciesTag to XML output stream
/*!
  \param os_xml  XML Output stream
  \param stag    SpeciesTag
  \param name    Optional name attribute
*/
/*  param pbofs   Pointer to binary file stream. NULL for ASCII output.
                 Ignore because SpeciesTag is always stored in ASCII format. */
void xml_write_to_stream(std::ostream& os_xml,
                         const SpeciesTag& stag,
                         bofstream* /* pbofs */,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("SpeciesTag");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  os_xml << '"' << stag.Name() << '"' << '\n';

  close_tag.set_name("/SpeciesTag");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== Sun =====================================================

//! Reads Sun from XML input stream
/*!
  \param is_xml  XML Input stream
  \param sun    Sun return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(std::istream& is_xml, Sun& sun, bifstream* pbifs) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("Sun");

  xml_read_from_stream(is_xml, sun.description, pbifs);
  xml_read_from_stream(is_xml, sun.spectrum, pbifs);
  xml_read_from_stream(is_xml, sun.radius, pbifs);
  xml_read_from_stream(is_xml, sun.distance, pbifs);
  xml_read_from_stream(is_xml, sun.latitude, pbifs);
  xml_read_from_stream(is_xml, sun.longitude, pbifs);

  tag.read_from_stream(is_xml);
  tag.check_name("/Sun");
}

//! Writes Sun to XML output stream
/*!
  \param os_xml  XML Output stream
  \param sun    Sun
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const Sun& sun,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Sun");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, sun.description, pbofs, "StarType");
  xml_write_to_stream(os_xml, sun.spectrum, pbofs, "StarSpectrum");
  xml_write_to_stream(os_xml, sun.radius, pbofs, "StarRadius");
  xml_write_to_stream(os_xml, sun.distance, pbofs, "StarDistance");
  xml_write_to_stream(os_xml, sun.latitude, pbofs, "StarLatitude");
  xml_write_to_stream(os_xml, sun.longitude, pbofs, "StarLongitude");

  close_tag.set_name("/Sun");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== TelsemAtlas ======================================================

//! Reads TelsemAtlas from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param pm         TelsemAtlas return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          TelsemAtlas& ta,
                          bifstream* pbifs) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("TelsemAtlas");

  xml_read_from_stream(is_xml, ta.ndat, pbifs);
  xml_read_from_stream(is_xml, ta.nchan, pbifs);
  xml_read_from_stream(is_xml, ta.name, pbifs);
  xml_read_from_stream(is_xml, ta.month, pbifs);
  xml_read_from_stream(is_xml, ta.dlat, pbifs);
  xml_read_from_stream(is_xml, ta.emis, pbifs);
  xml_read_from_stream(is_xml, ta.correl, pbifs);
  xml_read_from_stream(is_xml, ta.emis_err, pbifs);
  xml_read_from_stream(is_xml, ta.classes1, pbifs);
  xml_read_from_stream(is_xml, ta.classes2, pbifs);
  xml_read_from_stream(is_xml, ta.cellnums, pbifs);
  ta.telsem_calc_correspondence();
  tag.read_from_stream(is_xml);
  tag.check_name("/TelsemAtlas");
}

//! Writes TelsemAtlas to XML output stream
/*!
 * \param os_xml     XML Output stream
 * \param pm         TelsemAtlas
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const TelsemAtlas& ta,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("TelsemAtlas");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';
  xml_write_to_stream(os_xml, ta.ndat, pbofs, "ndat");
  xml_write_to_stream(os_xml, ta.nchan, pbofs, "nchan");
  xml_write_to_stream(os_xml, ta.name, pbofs, "name");
  xml_write_to_stream(os_xml, ta.month, pbofs, "month");
  xml_write_to_stream(os_xml, ta.dlat, pbofs, "dlat");
  xml_write_to_stream(os_xml, ta.emis, pbofs, "emis");
  xml_write_to_stream(os_xml, ta.correl, pbofs, "correl");
  xml_write_to_stream(os_xml, ta.emis_err, pbofs, "emis_err");
  xml_write_to_stream(os_xml, ta.classes1, pbofs, "class1");
  xml_write_to_stream(os_xml, ta.classes2, pbofs, "class2");
  xml_write_to_stream(os_xml, ta.cellnums, pbofs, "cellnum");
  close_tag.set_name("/TelsemAtlas");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== XsecRecord ======================================================

//! Reads XsecData from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param xd         XsecData return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          XsecRecord& xd,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  Index version;

  tag.read_from_stream(is_xml);
  tag.check_name("XsecRecord");
  tag.get_attribute_value("version", version);
  xd.SetVersion(version);

  String species_name;
  xml_read_from_stream(is_xml, species_name, pbifs);

  const SpeciesEnum species = to<SpeciesEnum>(species_name);
  if (not good_enum(species)) {
    std::ostringstream os;
    os << "  Unknown species in XsecRecord: " << species_name;
    throw std::runtime_error(os.str());
  }
  xd.SetSpecies(species);

  xml_read_from_stream(is_xml, xd.FitMinPressures(), pbifs);
  xml_read_from_stream(is_xml, xd.FitMaxPressures(), pbifs);
  xml_read_from_stream(is_xml, xd.FitMinTemperatures(), pbifs);
  xml_read_from_stream(is_xml, xd.FitMaxTemperatures(), pbifs);
  xml_read_from_stream(is_xml, xd.FitCoeffs(), pbifs);

  for (const auto& fitcoeffs : xd.FitCoeffs()) {
    const Index ncoeff = fitcoeffs.data.ncols();
    ARTS_USER_ERROR_IF(ncoeff != 4,
                       "Wrong number of coefficients, expected 4, found {}",
                       ncoeff)
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/XsecRecord");
}

//! Writes XsecData to XML output stream
/*!
 * \param os_xml     XML Output stream
 * \param xd         XsecData
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const XsecRecord& xd,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("XsecRecord");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("version", xd.Version());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';
  xml_write_to_stream(os_xml, xd.SpeciesName(), pbofs, "species");

  xml_write_to_stream(
      os_xml, xd.FitMinPressures(), pbofs, "Mininum pressures from fit");
  xml_write_to_stream(
      os_xml, xd.FitMaxPressures(), pbofs, "Maximum pressures from fit");
  xml_write_to_stream(
      os_xml, xd.FitMinTemperatures(), pbofs, "Mininum temperatures from fit");
  xml_write_to_stream(
      os_xml, xd.FitMaxTemperatures(), pbofs, "Maximum temperatures from fit");
  xml_write_to_stream(os_xml, xd.FitCoeffs(), pbofs, "Fit coefficients");

  close_tag.set_name("/XsecRecord");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== PredefinedModelData =========================================
/*!
 * \param is_xml     XML Input stream
 * \param pmd        PredefinedModelData return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          PredefinedModelData& pmd,
                          bifstream* pbifs) {
  ARTS_USER_ERROR_IF(pbifs, "No binary data")

  pmd = PredefinedModelData{};  // overwrite

  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("PredefinedModelData");

  Index nelem;
  open_tag.get_attribute_value("nelem", nelem);

  for (Index i = 0; i < nelem; i++) {
    ArtsXMLTag internal_open_tag;
    internal_open_tag.read_from_stream(is_xml);
    internal_open_tag.check_name("Data");

    // Get key
    String key_str;
    internal_open_tag.get_attribute_value("key", key_str);
    SpeciesIsotope isot{key_str};
    ARTS_USER_ERROR_IF(not isot.OK(), "Cannot understand key: \"{}\"", key_str)

    String sizes_str;
    internal_open_tag.get_attribute_value("sizes", sizes_str);

    Index sizes_len;
    internal_open_tag.get_attribute_value("sizes_nelem", sizes_len);

    String data_type;
    internal_open_tag.get_attribute_value("type", data_type);
    auto data = Absorption::PredefinedModel::model_data(data_type);

    std::vector<std::size_t> sizes(sizes_len);
    std::istringstream values(sizes_str);
    for (auto& sz : sizes) values >> sz;

    std::visit(
        [&](auto& v) {
          v.resize(sizes);
          is_xml >> v;
          pmd.data[isot] = v;
        },
        data);

    ArtsXMLTag internal_close_tag;
    internal_close_tag.read_from_stream(is_xml);
    internal_close_tag.check_name("/Data");
  }

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/PredefinedModelData");
}

/*!
 * \param os_xml     XML Output stream
 * \param rvb        PredefinedModelData
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const PredefinedModelData& pmd,
                         bofstream* pbofs,
                         const String& name) {
  ARTS_USER_ERROR_IF(pbofs, "No binary data")

  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("PredefinedModelData");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("nelem", Index(pmd.data.size()));
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (auto& [key, data] : pmd.data) {
    ArtsXMLTag internal_open_tag;
    internal_open_tag.set_name("Data");
    internal_open_tag.add_attribute("key", key.FullName());
    internal_open_tag.add_attribute(
        "type", String{Absorption::PredefinedModel::model_name(data)});

    const std::vector<Size> sizes =
        std::visit([&](auto& v) { return v.sizes(); }, data);
    internal_open_tag.add_attribute("sizes_nelem", Index(sizes.size()));

    String sizes_str = "";
    for (std::size_t i = 0; i < sizes.size(); i++) {
      if (i > 0) sizes_str += " ";
      sizes_str += var_string(sizes[i]);
    }
    internal_open_tag.add_attribute("sizes", sizes_str);

    internal_open_tag.write_to_stream(os_xml);
    os_xml << '\n';

    // Set values
    std::visit([&](auto& v) { os_xml << v; }, data);
    os_xml << '\n';

    ArtsXMLTag internal_close_tag;
    internal_close_tag.set_name("/Data");
    internal_close_tag.write_to_stream(os_xml);
    os_xml << '\n';
  }

  close_tag.set_name("/PredefinedModelData");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== AtmField =========================================
void xml_read_from_stream_helper(std::istream& is_xml,
                                 AtmKeyVal& key_val,
                                 Atm::Data& data,
                                 bifstream* pbifs) {
  data = Atm::Data{};  // overwrite

  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("AtmData");

  String keytype, key, type;
  open_tag.get_attribute_value("keytype", keytype);
  open_tag.get_attribute_value("key", key);
  open_tag.get_attribute_value("type", type);

  const auto get_extrapol = [&open_tag](auto& k) {
    String x;
    open_tag.get_attribute_value(k, x);
    return to<InterpolationExtrapolation>(x);
  };
  data.alt_low = get_extrapol("alt_low");
  data.alt_upp = get_extrapol("alt_upp");
  data.lat_low = get_extrapol("lat_low");
  data.lat_upp = get_extrapol("lat_upp");
  data.lon_low = get_extrapol("lon_low");
  data.lon_upp = get_extrapol("lon_upp");

  if (keytype == "AtmKey")
    key_val = to<AtmKey>(key);
  else if (keytype == "Species")
    key_val = to<SpeciesEnum>(key);
  else if (keytype == "IsotopeRecord")
    key_val = Species::Isotopologues[Species::find_species_index(key)];
  else if (keytype == "QuantumIdentifier")
    key_val = QuantumIdentifier(key);
  else
    ARTS_USER_ERROR("Cannot understand the keytype: \"{}\"", keytype)

  if (type == "GriddedField3")
    data.data = GriddedField3{};
  else if (type == "Numeric")
    data.data = Numeric{};
  else if (type == "FunctionalData")
    data.data = Atm::FunctionalData{
        Atm::FunctionalDataAlwaysThrow{"Cannot restore functional data"}};
  else
    ARTS_USER_ERROR("Cannot understand the data type: \"{}\"", type)

  std::visit(
      [&](auto& v) {
        if constexpr (std::same_as<std::remove_cvref_t<decltype(v)>,
                                   Atm::FunctionalData>) {
          String x;
          xml_read_from_stream(is_xml, x, pbifs);
          v = Atm::FunctionalData{Atm::FunctionalDataAlwaysThrow{x}};
        } else {
          xml_read_from_stream(is_xml, v, pbifs);
        }
      },
      data.data);

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/AtmData");
}

/*!
 * \param is_xml     XML Input stream
 * \param atm        AtmField return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          AtmField& atm,
                          bifstream* pbifs) {
  atm = AtmField{};  // overwrite

  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("AtmField");

  Index n;
  open_tag.get_attribute_value("nelem", n);

  open_tag.get_attribute_value("toa", atm.top_of_atmosphere);

  for (Index i = 0; i < n; i++) {
    AtmKeyVal key;
    Atm::Data data;
    xml_read_from_stream_helper(is_xml, key, data, pbifs);
    atm[key] = std::move(data);
  }

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/AtmField");
}

void xml_write_to_stream_helper(std::ostream& os_xml,
                                const AtmKeyVal& key,
                                const Atm::Data& data,
                                bofstream* pbofs) {
  ArtsXMLTag open_data_tag;
  open_data_tag.set_name("AtmData");

  std::visit(
      [&](auto& key_val) {
        if constexpr (Atm::isAtmKey<decltype(key_val)>)
          open_data_tag.add_attribute("keytype", "AtmKey");
        else if constexpr (Atm::isSpecies<decltype(key_val)>)
          open_data_tag.add_attribute("keytype", "Species");
        else if constexpr (Atm::isSpeciesIsotope<decltype(key_val)>)
          open_data_tag.add_attribute("keytype", "IsotopeRecord");
        else if constexpr (Atm::isQuantumIdentifier<decltype(key_val)>)
          open_data_tag.add_attribute("keytype", "QuantumIdentifier");
        else
          ARTS_ASSERT(false,
                      "New key type is not yet handled by writing routine!")

        if constexpr (Atm::isSpecies<decltype(key_val)>)
          open_data_tag.add_attribute("key", String{toString<1>(key_val)});
        else
          open_data_tag.add_attribute("key", var_string(key_val));

        open_data_tag.add_attribute("type", data.data_type());
        open_data_tag.add_attribute("alt_low", String{toString(data.alt_low)});
        open_data_tag.add_attribute("alt_upp", String{toString(data.alt_upp)});
        open_data_tag.add_attribute("lat_low", String{toString(data.lat_low)});
        open_data_tag.add_attribute("lat_upp", String{toString(data.lat_upp)});
        open_data_tag.add_attribute("lon_low", String{toString(data.lon_low)});
        open_data_tag.add_attribute("lon_upp", String{toString(data.lon_upp)});
        open_data_tag.write_to_stream(os_xml);
        os_xml << '\n';

        std::visit(
            [&](auto&& data_type [[maybe_unused]]) {
              if constexpr (std::same_as<
                                std::remove_cvref_t<decltype(data_type)>,
                                Atm::FunctionalData>)
                xml_write_to_stream(
                    os_xml,
                    String{var_string(
                        "Data for ",
                        key_val,
                        " read from file as functional must be set explicitly")},
                    pbofs,
                    "Functional Data Error");
              else
                xml_write_to_stream(os_xml, data_type, pbofs, "Data");
            },
            data.data);
      },
      key);

  ArtsXMLTag close_data_tag;
  close_data_tag.set_name("/AtmData");
  close_data_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

/*!
 * \param os_xml     XML Output stream
 * \param atm        AtmField
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const AtmField& atm,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("AtmField");
  if (name.length()) open_tag.add_attribute("name", name);

  //! List of all KEY values
  const auto keys = atm.keys();

  xml_set_stream_precision(os_xml);

  open_tag.add_attribute("nelem", static_cast<Index>(keys.size()));
  open_tag.add_attribute("toa", atm.top_of_atmosphere);
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (auto& key : keys) {
    xml_write_to_stream_helper(os_xml, key, atm[key], pbofs);
  }

  close_tag.set_name("/AtmField");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== AtmPoint =========================================
/*!
 * \param is_xml     XML Input stream
 * \param atm        AtmPoint return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          AtmPoint& atm,
                          bifstream* pbifs) {
  atm = AtmPoint{};  // overwrite

  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("AtmPoint");

  Index nspec, nnlte, nother, nisot;
  open_tag.get_attribute_value("nspec", nspec);
  open_tag.get_attribute_value("nnlte", nnlte);
  open_tag.get_attribute_value("nother", nother);
  open_tag.get_attribute_value("nisot", nisot);

  const Index n = nspec + nnlte + nother + nisot;
  for (Index i = 0; i < n; i++) {
    Numeric v;
    String k;
    xml_read_from_stream(is_xml, k, pbifs);
    xml_read_from_stream(is_xml, v, pbifs);

    if (nother > 0) {
      atm[to<AtmKey>(k)] = v;
      nother--;
    } else if (nspec > 0) {
      atm[to<SpeciesEnum>(k)] = v;
      nspec--;
    } else if (nnlte > 0) {
      atm[QuantumIdentifier(k)] = v;
      nnlte--;
    } else if (nisot > 0) {
      atm[SpeciesIsotope(k)] = v;
      nisot--;
    } else {
      ARTS_ASSERT(false)
    }
  }
  ARTS_ASSERT((nspec + nnlte + nother) == 0)

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/AtmPoint");
}

/*!
 * \param os_xml     XML Output stream
 * \param atm        AtmPoint
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const AtmPoint& atm,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("AtmPoint");
  if (name.length()) open_tag.add_attribute("name", name);

  //! List of all KEY values
  auto keys          = atm.keys();
  const Index nspec  = atm.nspec();
  const Index nnlte  = atm.nnlte();
  const Index nother = atm.nother();
  const Index nisot  = atm.nisot();

  open_tag.add_attribute("nspec", nspec);
  open_tag.add_attribute("nnlte", nnlte);
  open_tag.add_attribute("nother", nother);
  open_tag.add_attribute("nisot", nisot);
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (auto& key : keys) {
    std::visit(
        [&](auto&& key_val) {
          xml_write_to_stream(
              os_xml, String{var_string(key_val)}, pbofs, "Data Key");
          xml_write_to_stream(os_xml, atm[key_val], pbofs, "Data");
        },
        key);
  }

  close_tag.set_name("/AtmPoint");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== SurfaceField =========================================
void xml_read_from_stream_helper(std::istream& is_xml,
                                 SurfaceKeyVal& key_val,
                                 Surf::Data& data,
                                 bifstream* pbifs) {
  data = Surf::Data{};  // overwrite

  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("SurfaceData");

  String keytype, key, type;
  open_tag.get_attribute_value("keytype", keytype);
  open_tag.get_attribute_value("key", key);
  open_tag.get_attribute_value("type", type);

  const auto get_extrapol = [&open_tag](auto& k) {
    String x;
    open_tag.get_attribute_value(k, x);
    return to<InterpolationExtrapolation>(x);
  };
  data.lat_low = get_extrapol("lat_low");
  data.lat_upp = get_extrapol("lat_upp");
  data.lon_low = get_extrapol("lon_low");
  data.lon_upp = get_extrapol("lon_upp");

  if (keytype == "SurfaceKey")
    key_val = to<SurfaceKey>(key);
  else if (keytype == "SurfaceTypeTag")
    key_val = SurfaceTypeTag{key};
  else
    ARTS_USER_ERROR("Cannot understand the keytype: \"{}\"", keytype)

  if (type == "GriddedField2")
    data.data = GriddedField2{};
  else if (type == "Numeric")
    data.data = Numeric{};
  else if (type == "FunctionalData")
    data.data = Surf::FunctionalData{
        Surf::FunctionalDataAlwaysThrow{"Cannot restore functional data"}};
  else
    ARTS_USER_ERROR("Cannot understand the data type: \"{}\"", type)

  std::visit(
      [&](auto& v) {
        if constexpr (std::same_as<std::remove_cvref_t<decltype(v)>,
                                   Surf::FunctionalData>) {
          String x;
          xml_read_from_stream(is_xml, x, pbifs);
          v = Surf::FunctionalData{Surf::FunctionalDataAlwaysThrow{x}};
        } else {
          xml_read_from_stream(is_xml, v, pbifs);
        }
      },
      data.data);

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/SurfaceData");
}

/*!
 * \param is_xml     XML Input stream
 * \param surf       SurfaceField return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          SurfaceField& surf,
                          bifstream* pbifs) {
  surf = SurfaceField{};  // overwrite

  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("SurfaceField");

  Index n;
  open_tag.get_attribute_value("nelem", n);

  for (Index i = 0; i < n; i++) {
    SurfaceKeyVal key;
    Surf::Data data;
    xml_read_from_stream_helper(is_xml, key, data, pbifs);
    surf[key] = std::move(data);
  }

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/SurfaceField");
}

void xml_write_to_stream_helper(std::ostream& os_xml,
                                const SurfaceKeyVal& key,
                                const Surf::Data& data,
                                bofstream* pbofs) {
  ArtsXMLTag open_data_tag;
  open_data_tag.set_name("SurfaceData");

  std::visit(
      [&](auto& key_val) {
        if constexpr (Surf::isSurfaceKey<decltype(key_val)>)
          open_data_tag.add_attribute("keytype", "AtmKey");
        else if constexpr (Surf::isSurfaceTypeTag<decltype(key_val)>)
          open_data_tag.add_attribute("keytype", "SurfaceTypeTag");
        else
          ARTS_ASSERT(false,
                      "New key type is not yet handled by writing routine!")

        open_data_tag.add_attribute("key", var_string(key_val));
        open_data_tag.add_attribute("type", data.data_type());
        open_data_tag.add_attribute("lat_low", String{toString(data.lat_low)});
        open_data_tag.add_attribute("lat_upp", String{toString(data.lat_upp)});
        open_data_tag.add_attribute("lon_low", String{toString(data.lon_low)});
        open_data_tag.add_attribute("lon_upp", String{toString(data.lon_upp)});
        open_data_tag.write_to_stream(os_xml);
        os_xml << '\n';

        std::visit(
            [&](auto&& data_type [[maybe_unused]]) {
              if constexpr (std::same_as<
                                std::remove_cvref_t<decltype(data_type)>,
                                Surf::FunctionalData>)
                xml_write_to_stream(
                    os_xml,
                    String{var_string(
                        "Data for ",
                        key_val,
                        " read from file as functional must be set explicitly")},
                    pbofs,
                    "Functional Data Error");
              else
                xml_write_to_stream(os_xml, data_type, pbofs, "Data");
            },
            data.data);
      },
      key);

  ArtsXMLTag close_data_tag;
  close_data_tag.set_name("/SurfaceData");
  close_data_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

/*!
 * \param os_xml     XML Output stream
 * \param surf       SurfaceField
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const SurfaceField& surf,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("SurfaceField");
  if (name.length()) open_tag.add_attribute("name", name);

  //! List of all KEY values
  const auto keys = surf.keys();

  xml_set_stream_precision(os_xml);

  open_tag.add_attribute("nelem", static_cast<Index>(keys.size()));
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (auto& key : keys) {
    xml_write_to_stream_helper(os_xml, key, surf[key], pbofs);
  }

  close_tag.set_name("/SurfaceField");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== SurfacePoint =========================================
/*!
 * \param is_xml     XML Input stream
 * \param surf       SurfacePoint return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          SurfacePoint& surf,
                          bifstream* pbifs) {
  surf = SurfacePoint{};  // overwrite

  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("SurfacePoint");

  Index ntype, nother;
  open_tag.get_attribute_value("ntype", ntype);
  open_tag.get_attribute_value("nother", nother);
  open_tag.get_attribute_value("zenith", surf.normal[0]);
  open_tag.get_attribute_value("azimuth", surf.normal[1]);

  const Index n = ntype + nother;
  for (Index i = 0; i < n; i++) {
    Numeric v;
    String k;
    xml_read_from_stream(is_xml, k, pbifs);
    xml_read_from_stream(is_xml, v, pbifs);

    if (nother > 0) {
      surf[to<SurfaceKey>(k)] = v;
      nother--;
    } else if (ntype > 0) {
      surf[SurfaceTypeTag{k}] = v;
      ntype--;
    } else {
      ARTS_ASSERT(false)
    }
  }
  ARTS_ASSERT((ntype + nother) == 0)

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/SurfacePoint");
}

/*!
 * \param os_xml     XML Output stream
 * \param surf       SurfacePoint
 * \param pbofs      Pointer to binary file stream. NULL for ASCII output.
 * \param name       Optional name attribute
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const SurfacePoint& surf,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("SurfacePoint");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("zenith", surf.normal[0]);
  open_tag.add_attribute("azimuth", surf.normal[1]);

  //! List of all KEY values
  auto keys          = surf.keys();
  const Index ntype  = surf.ntype();
  const Index nother = surf.nother();

  open_tag.add_attribute("ntype", ntype);
  open_tag.add_attribute("nother", nother);
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (auto& key : keys) {
    std::visit(
        [&](auto&& key_val) {
          xml_write_to_stream(
              os_xml, String{var_string(key_val)}, pbofs, "Data Key");
          xml_write_to_stream(os_xml, surf[key_val], pbofs, "Data");
        },
        key);
  }

  close_tag.set_name("/SurfacePoint");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== Wsv ================================================

void xml_read_from_stream(std::istream& is_xml, Wsv& wsv, bifstream* pbifs) {
  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("Wsv");

  String type;
  open_tag.get_attribute_value("type", type);

  wsv = Wsv::from_named_type(type);
  std::visit([&](auto& x) { xml_read_from_stream(is_xml, *x, pbifs); },
             wsv.value());

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/Wsv");
}

void xml_write_to_stream(std::ostream& os_xml,
                         const Wsv& wsv,
                         bofstream* pbofs,
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Wsv");
  open_tag.add_attribute("type", String{wsv.type_name()});
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  std::visit([&](auto& x) { xml_write_to_stream(os_xml, *x, pbofs, ""); },
             wsv.value());

  close_tag.set_name("/Wsv");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== Method ================================================

void xml_read_from_stream(std::istream& is_xml, Method& m, bifstream* pbifs) {
  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("Method");

  String name;
  open_tag.get_attribute_value("name", name);

  String type;
  open_tag.get_attribute_value("type", type);

  if (type == "call") {
    ArrayOfString args_;
    xml_read_from_stream(is_xml, args_, pbifs);

    const std::vector<std::string> args{args_.begin(), args_.end()};
    m = Method{name, args, std::unordered_map<std::string, std::string>{}};
  } else if (type == "value") {
    Index overwrite;
    open_tag.get_attribute_value("overwrite", overwrite);

    Wsv wsv;
    xml_read_from_stream(is_xml, wsv, pbifs);

    m = Method{name, wsv, static_cast<bool>(overwrite)};
  } else {
    ARTS_USER_ERROR("Bad type: \"{}\"", type)
  }

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/Method");
}

void xml_write_to_stream(std::ostream& os_xml,
                         const Method& m,
                         bofstream* pbofs,
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Method");
  open_tag.add_attribute("name", m.get_name());

  const bool is_value = static_cast<bool>(m.get_setval());

  if (is_value) {
    open_tag.add_attribute("type", "value");
    open_tag.add_attribute("overwrite", static_cast<Index>(m.overwrite()));
    open_tag.write_to_stream(os_xml);
    os_xml << '\n';

    xml_write_to_stream(os_xml, m.get_setval().value(), pbofs, "");
  } else {
    open_tag.add_attribute("type", "call");
    open_tag.write_to_stream(os_xml);
    os_xml << '\n';

    ArrayOfString args;
    for (const auto& a : m.get_outs()) {
      args.emplace_back(a);
    }
    for (const auto& a : m.get_ins()) {
      if (std::ranges::any_of(args, Cmp::eq(a))) continue;
      args.emplace_back(a);
    }

    xml_write_to_stream(os_xml, args, pbofs, "");
  }

  close_tag.set_name("/Method");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== Agenda ================================================

void xml_read_from_stream(std::istream& is_xml, Agenda& a, bifstream* pbifs) {
  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("Agenda");

  a = [&open_tag]() {
    String x;
    open_tag.get_attribute_value("name", x);
    return Agenda{x};
  }();

  const Index size = [&open_tag]() {
    Index x;
    open_tag.get_attribute_value("size", x);
    return x;
  }();

  for (Index i = 0; i < size; i++) {
    Method m;
    xml_read_from_stream(is_xml, m, pbifs);
    a.add(m);
  }

  a.finalize(true);

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/Agenda");
}

void xml_write_to_stream(std::ostream& os_xml,
                         const Agenda& a,
                         bofstream* pbofs,
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Agenda");
  open_tag.add_attribute("name", a.get_name());
  open_tag.add_attribute("size", static_cast<Index>(a.get_methods().size()));
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (const auto& m : a.get_methods()) {
    xml_write_to_stream(os_xml, m, pbofs, "");
  }

  close_tag.set_name("/Agenda");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== Any ================================================

void xml_read_from_stream(std::istream& is_xml, Any& a, bifstream*) {
  a = Any{};

  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("Any");

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/Any");
}

void xml_write_to_stream(std::ostream& os_xml,
                         const Any&,
                         bofstream*,
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Any");
  open_tag.write_to_stream(os_xml);
  close_tag.set_name("/Any");
  close_tag.write_to_stream(os_xml);
}

//=== MCAntenna ================================================

void xml_read_from_stream(std::istream& is_xml,
                          MCAntenna& m,
                          bifstream* pbifs) {
  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("MCAntenna");

  Index n;
  xml_read_from_stream(is_xml, n, pbifs);
  m.atype = static_cast<AntennaType>(n);
  xml_read_from_stream(is_xml, m.sigma_aa, pbifs);
  xml_read_from_stream(is_xml, m.sigma_za, pbifs);
  xml_read_from_stream(is_xml, m.aa_grid, pbifs);
  xml_read_from_stream(is_xml, m.za_grid, pbifs);
  xml_read_from_stream(is_xml, m.G_lookup, pbifs);

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/MCAntenna");
}

void xml_write_to_stream(std::ostream& os_xml,
                         const MCAntenna& m,
                         bofstream* pbofs,
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("MCAntenna");
  open_tag.write_to_stream(os_xml);

  const auto n = static_cast<Index>(m.atype);
  xml_write_to_stream(os_xml, n, pbofs, "");
  xml_write_to_stream(os_xml, m.sigma_aa, pbofs, "");
  xml_write_to_stream(os_xml, m.sigma_za, pbofs, "");
  xml_write_to_stream(os_xml, m.aa_grid, pbofs, "");
  xml_write_to_stream(os_xml, m.za_grid, pbofs, "");
  xml_write_to_stream(os_xml, m.G_lookup, pbofs, "");

  close_tag.set_name("/MCAntenna");
  close_tag.write_to_stream(os_xml);
}

//=== TessemNN ================================================

void xml_read_from_stream(std::istream& is_xml, TessemNN& t, bifstream* pbifs) {
  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("TessemNN");

  xml_read_from_stream(is_xml, t.nb_inputs, pbifs);
  xml_read_from_stream(is_xml, t.nb_outputs, pbifs);
  xml_read_from_stream(is_xml, t.nb_cache, pbifs);
  xml_read_from_stream(is_xml, t.b1, pbifs);
  xml_read_from_stream(is_xml, t.b2, pbifs);
  xml_read_from_stream(is_xml, t.w1, pbifs);
  xml_read_from_stream(is_xml, t.w2, pbifs);
  xml_read_from_stream(is_xml, t.x_min, pbifs);
  xml_read_from_stream(is_xml, t.x_max, pbifs);
  xml_read_from_stream(is_xml, t.y_min, pbifs);
  xml_read_from_stream(is_xml, t.y_max, pbifs);

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/TessemNN");
}

void xml_write_to_stream(std::ostream& os_xml,
                         const TessemNN& t,
                         bofstream* pbofs,
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("TessemNN");
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, t.nb_inputs, pbofs, "");
  xml_write_to_stream(os_xml, t.nb_outputs, pbofs, "");
  xml_write_to_stream(os_xml, t.nb_cache, pbofs, "");
  xml_write_to_stream(os_xml, t.b1, pbofs, "");
  xml_write_to_stream(os_xml, t.b2, pbofs, "");
  xml_write_to_stream(os_xml, t.w1, pbofs, "");
  xml_write_to_stream(os_xml, t.w2, pbofs, "");
  xml_write_to_stream(os_xml, t.x_min, pbofs, "");
  xml_write_to_stream(os_xml, t.x_max, pbofs, "");
  xml_write_to_stream(os_xml, t.y_min, pbofs, "");
  xml_write_to_stream(os_xml, t.y_max, pbofs, "");

  close_tag.set_name("/TessemNN");
  close_tag.write_to_stream(os_xml);
}

////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

// NOTE: These cannot be fixed by adding the missing functions, because
//       they cannot be completely restructured by the data they own.
//       If you need them to be, consider a redesign of them

//=== JacobianTargets =========================================

void xml_read_from_stream(std::istream&,
                          JacobianTargets&,
                          bifstream* /* pbifs */) {
  ARTS_USER_ERROR("Method not implemented!");
}

void xml_write_to_stream(std::ostream&,
                         const JacobianTargets&,
                         bofstream* /* pbofs */,
                         const String& /* name */) {
  ARTS_USER_ERROR("Method not implemented!");
}

//=== SpectralRadianceOperator =========================================

void xml_read_from_stream(std::istream&,
                          SpectralRadianceOperator&,
                          bifstream* /* pbifs */) {
  ARTS_USER_ERROR("Method not implemented!");
}

void xml_write_to_stream(std::ostream&,
                         const SpectralRadianceOperator&,
                         bofstream* /* pbofs */,
                         const String& /* name */) {
  ARTS_USER_ERROR("Method not implemented!");
}

//=== AbsorptionBand =========================================

void xml_read_from_stream(std::istream& is_xml,
                          lbl::band_data& data,
                          bifstream* pbifs) try {
  ARTS_USER_ERROR_IF(pbifs not_eq nullptr, "No binary data")

  String tag;
  Index nelem;

  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("AbsorptionBandData");

  open_tag.get_attribute_value("lineshape", tag);
  data.lineshape = to<LineByLineLineshape>(tag);

  open_tag.get_attribute_value("cutoff_type", tag);
  data.cutoff = to<LineByLineCutoffType>(tag);

  open_tag.get_attribute_value("cutoff_value", data.cutoff_value);

  open_tag.get_attribute_value("nelem", nelem);
  data.lines.resize(0);
  data.lines.reserve(nelem);

  for (Index j = 0; j < nelem; j++) {
    is_xml >> data.lines.emplace_back();
  }

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/AbsorptionBandData");
}
ARTS_METHOD_ERROR_CATCH

void xml_write_to_stream(std::ostream& os_xml,
                         const lbl::band_data& data,
                         bofstream* pbofs,
                         const String& name) {
  ARTS_USER_ERROR_IF(pbofs not_eq nullptr, "No binary data")

  ArtsXMLTag open_tag;
  open_tag.set_name("AbsorptionBandData");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("lineshape", String{toString(data.lineshape)});
  open_tag.add_attribute("cutoff_type", String{toString(data.cutoff)});
  open_tag.add_attribute("cutoff_value", data.cutoff_value);
  open_tag.add_attribute("nelem", static_cast<Index>(data.lines.size()));
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (auto& line : data) {
    os_xml << line << '\n';
  }

  ArtsXMLTag close_tag;
  close_tag.set_name("/AbsorptionBandData");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

void xml_read_from_stream(std::istream& is_xml,
                          AbsorptionBand& data,
                          bifstream* pbifs) try {
  ArtsXMLTag open_tag;
  open_tag.read_from_stream(is_xml);
  open_tag.check_name("AbsorptionBand");

  xml_read_from_stream(is_xml, data.key, pbifs);
  xml_read_from_stream(is_xml, data.data, pbifs);

  ArtsXMLTag close_tag;
  close_tag.read_from_stream(is_xml);
  close_tag.check_name("/AbsorptionBand");
}
ARTS_METHOD_ERROR_CATCH

void xml_write_to_stream(std::ostream& os_xml,
                         const AbsorptionBand& data,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  open_tag.set_name("AbsorptionBand");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, data.key, pbofs, "");
  xml_write_to_stream(os_xml, data.data, pbofs, "");

  ArtsXMLTag close_tag;
  close_tag.set_name("/AbsorptionBand");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== LinemixingEcsData =========================================
struct meta_data {
  String name;
  String value{};

  void add(ArtsXMLTag& tag, const meta_data& m) {
    tag.add_attribute(m.name, m.value);
  }

  void get(ArtsXMLTag& tag, meta_data& m) {
    tag.get_attribute_value(m.name, m.value);
  }

  meta_data(String n) : name(std::move(n)) {}
  meta_data(String n, String v) : name(std::move(n)), value(std::move(v)) {}
  meta_data(String n, auto v) : name(std::move(n)), value(var_string(v)) {}
};

struct Empty {};

template <bool read>
struct tag {
  static constexpr bool write = not read;
  using T = std::conditional_t<read, std::istream, std::ostream>;

  T& stream;
  String name;
  std::conditional_t<read, std::unordered_map<String, String>, Empty> data{};

  tag(T& s, String n) : stream(s), name(std::move(n)) {
    if constexpr (write) {
      ArtsXMLTag open_tag;
      open_tag.set_name(name);
      open_tag.write_to_stream(stream);
      stream << '\n';
    } else if (read) {
      ArtsXMLTag open_tag;
      open_tag.read_from_stream(stream);
      open_tag.check_name(name);
    }
  }

  template <typename... Ts>
  tag(T& s, String n, Ts&&... args) : stream(s), name(std::move(n)) {
    if constexpr (write) {
      ArtsXMLTag open_tag;
      open_tag.set_name(name);
      (open_tag.add_attribute(args.name, args.value), ...);
      open_tag.write_to_stream(stream);
      stream << '\n';
    } else if (read) {
      ArtsXMLTag open_tag;
      open_tag.read_from_stream(stream);
      open_tag.check_name(name);
      (open_tag.get_attribute_value(args, data[args]), ...);
    }
  }

  void close() {
    if constexpr (std::same_as<T, std::ostream>) {
      ArtsXMLTag close_tag;
      close_tag.set_name("/" + name);
      close_tag.write_to_stream(stream);
      stream << '\n';
    } else if constexpr (std::same_as<T, std::istream>) {
      ArtsXMLTag close_tag;
      close_tag.read_from_stream(stream);
      close_tag.check_name("/" + name);
    }
  }

  template <typename T = String>
  [[nodiscard]] T get(const String& key) const
    requires(read)
  {
    auto& x = data.at(key);
    if constexpr (std::same_as<T, Index>) {
      return std::stoi(x);
    } else if constexpr (requires { to<T>(x); }) {
      auto s = to<T>(x);
      return s;
    } else {
      return T{x};
    }
  }
};

template <typename... Ts>
tag(std::ostream& s, String n, Ts&&...) -> tag<false>;

template <typename... Ts>
tag(std::istream& s, String n, Ts&&...) -> tag<true>;

void xml_read_from_stream(std::istream& is_xml,
                          LinemixingEcsData& ecsd,
                          bifstream* pbifs) {
  ARTS_USER_ERROR_IF(pbifs not_eq nullptr, "No binary data")

  tag rtag{is_xml, "LinemixingEcsData", "nelem"};
  const auto nisot = rtag.get<Index>("nelem");

  // Get values
  ecsd.clear();
  ecsd.reserve(nisot);
  for (Index j = 0; j < nisot; j++) {
    tag isot_tag{is_xml, "isotdata", "nelem", "isot"};
    const auto nspec   = isot_tag.get<Index>("nelem");
    const auto isotkey = isot_tag.get<SpeciesIsotope>("isot");

    auto& spec_data_map = ecsd[isotkey];
    for (Index i = 0; i < nspec; i++) {
      tag spec_tag{is_xml, "specdata", "spec"};
      const auto speckey = spec_tag.get<SpeciesEnum>("spec");

      auto& spec_data = spec_data_map[speckey];
      is_xml >> spec_data.beta >> spec_data.collisional_distance >>
          spec_data.lambda >> spec_data.scaling;
      spec_tag.close();
    }
    isot_tag.close();
  }
  rtag.close();
}

void xml_write_to_stream(std::ostream& os_xml,
                         const LinemixingEcsData& r,
                         bofstream* pbofs,
                         const String&) {
  ARTS_USER_ERROR_IF(pbofs not_eq nullptr, "No binary data")

  tag rtag{os_xml,
           "LinemixingEcsData",
           meta_data{"nelem", static_cast<Index>(r.size())}};

  // Set values
  for (auto& [key, value] : r) {
    tag isot_tag{os_xml,
                 "isotdata",
                 meta_data{"nelem", static_cast<Index>(value.size())},
                 meta_data{"isot", key.FullName()}};

    for (auto& [ikey, ivalue] : value) {
      tag spec_tag{os_xml, "specdata", meta_data{"spec", toString<1>(ikey)}};

      os_xml << ivalue.beta << ' ' << ivalue.collisional_distance << ' '
             << ivalue.lambda << ' ' << ivalue.scaling;
      spec_tag.close();
    }
    isot_tag.close();
  }
  rtag.close();
}

//! SpeciesIsotope
#include <iostream>
void xml_read_from_stream(std::istream& is_xml, SpeciesIsotope& s, bifstream*) {
  tag stag{is_xml, "SpeciesIsotope", "isot"};
  s = stag.get<SpeciesIsotope>("isot");
  stag.close();
}

void xml_write_to_stream(std::ostream& os_xml,
                         const SpeciesIsotope& s,
                         bofstream*,
                         const String&) {
  tag stag{os_xml, "SpeciesIsotope", meta_data{"isot", s.FullName()}};
  stag.close();
}

//! ArrayOfPropagationPathPoint

void xml_read_from_stream(std::istream& is_xml,
                          PropagationPathPoint& ppp,
                          bifstream* pbifs) {
  tag stag{is_xml, "PropagationPathPoint", "pos_t", "los_t"};

  ppp.pos_type = stag.get<PathPositionType>("pos_t");
  ppp.los_type = stag.get<PathPositionType>("los_t");

  if (pbifs == nullptr) {
    is_xml >> double_imanip{} >> ppp.pos[0] >> ppp.pos[1] >> ppp.pos[2] >>
        ppp.los[0] >> ppp.los[1] >> ppp.nreal >> ppp.ngroup;
  } else {
    pbifs->readDoubleArray(&ppp.pos[0], 7);
  }

  stag.close();
}

void xml_write_to_stream(std::ostream& os_xml,
                         const PropagationPathPoint& ppp,
                         bofstream* pbofs,
                         const String&) {
  tag stag{os_xml,
           "PropagationPathPoint",
           meta_data{"pos_t", String{toString(ppp.pos_type)}},
           meta_data{"los_t", String{toString(ppp.los_type)}}};

  if (pbofs == nullptr) {
    os_xml << ppp.pos[0] << ' ' << ppp.pos[1] << ' ' << ppp.pos[2] << ' '
           << ppp.los[0] << ' ' << ppp.los[1] << ' ' << ppp.nreal << ' '
           << ppp.ngroup;
  } else {
    *pbofs << ppp.pos[0] << ppp.pos[1] << ppp.pos[2] << ppp.los[0] << ppp.los[1]
           << ppp.nreal << ppp.ngroup;
  }

  stag.close();
}

//! AscendingGrid

void xml_read_from_stream(std::istream& is_xml,
                          AscendingGrid& g,
                          bifstream* pbifs) try {
  tag stag{is_xml, "AscendingGrid"};

  Vector x;
  xml_read_from_stream(is_xml, x, pbifs);
  g = std::move(x);

  stag.close();
}
ARTS_METHOD_ERROR_CATCH

void xml_write_to_stream(std::ostream& os_xml,
                         const AscendingGrid& g,
                         bofstream* pbofs,
                         const String&) {
  tag stag{os_xml, "AscendingGrid"};
  xml_write_to_stream(os_xml, static_cast<const Vector&>(g), pbofs, "");
  stag.close();
}

//! SensorObsel

void xml_read_from_stream(std::istream& is_xml,
                          SensorObsel& g,
                          bifstream* pbifs) {
  tag stag{is_xml, "SensorObsel"};

  AscendingGrid f_grid;
  SensorPosLosVector poslos_grid;
  StokvecMatrix weight_matrix;

  xml_read_from_stream(is_xml, f_grid, pbifs);
  xml_read_from_stream(is_xml, poslos_grid, pbifs);
  xml_read_from_stream(is_xml, weight_matrix, pbifs);

  g = SensorObsel{f_grid, poslos_grid, weight_matrix};

  stag.close();
}

void xml_write_to_stream(std::ostream& os_xml,
                         const SensorObsel& g,
                         bofstream* pbofs,
                         const String&) {
  tag stag{os_xml, "SensorObsel"};

  xml_write_to_stream(os_xml, g.f_grid(), pbofs, "f_grid");
  xml_write_to_stream(os_xml, g.poslos_grid(), pbofs, "poslos");
  xml_write_to_stream(os_xml, g.weight_matrix(), pbofs, "weight_matrix");

  stag.close();
}

//! SensorPosLos

void xml_read_from_stream(std::istream& is_xml,
                          SensorPosLos& g,
                          bifstream* pbifs) {
  tag stag{is_xml, "SensorPosLos"};

  if (pbifs) {
    pbifs->readDoubleArray(g.pos.data.data(), 3);
    pbifs->readDoubleArray(g.los.data.data(), 2);
  } else {
    is_xml >> double_imanip() >> g.pos[0] >> g.pos[1] >> g.pos[2] >> g.los[0] >>
        g.los[1];
  }

  stag.close();
}

void xml_write_to_stream(std::ostream& os_xml,
                         const SensorPosLos& g,
                         bofstream* pbofs,
                         const String&) {
  tag stag{os_xml, "SensorPosLos"};

  if (pbofs) {
    *pbofs << g.pos[0] << g.pos[1] << g.pos[2] << g.los[0] << g.los[1];
  } else {
    os_xml << g.pos[0] << ' ' << g.pos[1] << ' ' << g.pos[2] << ' ' << g.los[0]
           << ' ' << g.los[1];
  }

  stag.close();
}

//! SensorPosLosVector

void xml_read_from_stream(std::istream& is_xml,
                          SensorPosLosVector& g,
                          bifstream* pbifs) {
  tag stag{is_xml, "SensorPosLosVector", "nelem"};

  const auto n = stag.get<Index>("nelem");
  g.resize(n);
  for (auto& i : g) {
    xml_read_from_stream(is_xml, i, pbifs);
  }

  stag.close();
}

void xml_write_to_stream(std::ostream& os_xml,
                         const SensorPosLosVector& g,
                         bofstream* pbofs,
                         const String&) {
  tag stag{os_xml,
           "SensorPosLosVector",
           meta_data{"nelem", static_cast<Index>(g.size())}};

  for (const auto& i : g) {
    xml_write_to_stream(os_xml, i, pbofs, "poslos");
  }

  stag.close();
}

//! DisortBDRF

void xml_read_from_stream(std::istream&, DisortBDRF&, bifstream*) {
  ARTS_USER_ERROR("Method not implemented!");
}

void xml_write_to_stream(std::ostream&,
                         const DisortBDRF&,
                         bofstream*,
                         const String&) {
  ARTS_USER_ERROR("Method not implemented!");
}

//! MatrixOfDisortBDRF

void xml_read_from_stream(std::istream& is_xml,
                          MatrixOfDisortBDRF& x,
                          bifstream*) {
  x.resize(0, 0);

  ArtsXMLTag tag;
  tag.read_from_stream(is_xml);
  tag.check_name("MatrixOfDisortBDRF");

  Index n;
  tag.get_attribute_value("ncols", n);
  ARTS_USER_ERROR_IF(n != 0, "Only for zero size as workaround")
  tag.get_attribute_value("nrows", n);
  ARTS_USER_ERROR_IF(n != 0, "Only for zero size as workaround")

  tag.read_from_stream(is_xml);
  tag.check_name("/MatrixOfDisortBDRF");
}

void xml_write_to_stream(std::ostream& os_xml,
                         const MatrixOfDisortBDRF& x,
                         bofstream*,
                         const String& name) {
  if (x.size() == 0) {
    os_xml << std::format(
        R"(<MatrixOfDisortBDRF ncols="{}" nrows="{}" name="{}"> </MatrixOfDisortBDRF>)",
        x.ncols(),
        x.nrows(),
        name);
  } else {
    ARTS_USER_ERROR("Cannot store BDRF data to XML")
  }
}

//=== AtmKeyVal =========================================================

//! Reads AtmKeyVal from XML input stream
/*!
  \param is_xml    XML Input stream
  \param covmat    AtmKeyVal
  \param pbifs     Pointer to binary file stream. NULL for ASCII output.
*/
void xml_read_from_stream(std::istream& is_xml,
                          AtmKeyVal& atmt,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  tag.read_from_stream(is_xml);
  tag.check_name("AtmKeyVal");

  String type;
  tag.get_attribute_value("type", type);

  if (type == "AtmKey"s) {
    atmt = AtmKey{};
  } else if (type == "SpeciesEnum"s) {
    atmt = SpeciesEnum{};
  } else if (type == "SpeciesIsotope"s) {
    atmt = SpeciesIsotope{};
  } else if (type == "QuantumIdentifier"s) {
    atmt = QuantumIdentifier{};
  } else if (type == "ScatteringSpeciesProperty") {
    atmt = ScatteringSpeciesProperty{};
  } else {
    ARTS_USER_ERROR(R"(Cannot understand type: "{}")", type);
  }

  std::visit(
      [&](auto& data_type) { xml_read_from_stream(is_xml, data_type, pbifs); },
      atmt);

  tag.read_from_stream(is_xml);
  tag.check_name("/AtmKeyVal");
}

//! Write AtmKeyVal to XML output stream
/*!
  \param os_xml    XML output stream
  \param covmat    AtmKeyVal
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Unused
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const AtmKeyVal& atmt,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("AtmKeyVal");
  open_tag.add_attribute(
      "type",
      std::visit(
          []<typename T>(const T&) {
            if constexpr (std::same_as<T, AtmKey>)
              return "AtmKey"s;
            else if constexpr (std::same_as<T, SpeciesEnum>)
              return "SpeciesEnum"s;
            else if constexpr (std::same_as<T, SpeciesIsotope>)
              return "SpeciesIsotope"s;
            else if constexpr (std::same_as<T, QuantumIdentifier>)
              return "QuantumIdentifier"s;
            else
              return "UnknownType"s;
          },
          atmt));
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  std::visit(
      [&](const auto& data_type) {
        xml_write_to_stream(os_xml, data_type, pbofs, name);
      },
      atmt);

  close_tag.set_name("/AtmKeyVal");
  close_tag.write_to_stream(os_xml);
}

//=== SurfaceKeyVal =========================================================

//! Reads SurfaceKeyVal from XML input stream
/*!
  \param is_xml    XML Input stream
  \param covmat    SurfaceKeyVal
  \param pbifs     Pointer to binary file stream. NULL for ASCII output.
*/
void xml_read_from_stream(std::istream& is_xml,
                          SurfaceKeyVal& surft,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  tag.read_from_stream(is_xml);
  tag.check_name("SurfaceKeyVal");

  String type;
  tag.get_attribute_value("type", type);

  if (type == "SurfaceKey"s) {
    surft = SurfaceKey{};
  } else if (type == "SurfaceTypeTag"s) {
    surft = SurfaceTypeTag{};
  } else if (type == "SurfacePropertyTag"s) {
    surft = SurfacePropertyTag{};
  } else {
    ARTS_USER_ERROR(R"(Cannot understand type: "{}")", type);
  }

  std::visit(
      [&](auto& data_type) { xml_read_from_stream(is_xml, data_type, pbifs); },
      surft);

  tag.read_from_stream(is_xml);
  tag.check_name("/SurfaceKeyVal");
}

//! Write SurfaceKeyVal to XML output stream
/*!
  \param os_xml    XML output stream
  \param covmat    SurfaceKeyVal
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Unused
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const SurfaceKeyVal& surft,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("SurfaceKeyVal");
  open_tag.add_attribute(
      "type",
      std::visit(
          []<typename T>(const T&) {
            if constexpr (std::same_as<T, SurfaceKey>)
              return "SurfaceKey"s;
            else if constexpr (std::same_as<T, SurfaceTypeTag>)
              return "SurfaceTypeTag"s;
            else if constexpr (std::same_as<T, SurfacePropertyTag>)
              return "SurfacePropertyTag"s;
            else
              return "UnknownType"s;
          },
          surft));
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  std::visit(
      [&](const auto& data_type) {
        xml_write_to_stream(os_xml, data_type, pbofs, name);
      },
      surft);

  close_tag.set_name("/SurfaceKeyVal");
  close_tag.write_to_stream(os_xml);
}

//=== LblLineKey =========================================================

//! Reads LblLineKey from XML input stream
/*!
  \param is_xml    XML Input stream
  \param covmat    LblLineKey
  \param pbifs     Pointer to binary file stream. NULL for ASCII output.
*/
void xml_read_from_stream(std::istream& is_xml,
                          LblLineKey& lblt,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  tag.read_from_stream(is_xml);
  tag.check_name("LblLineKey");

  String var;
  tag.get_attribute_value("ls_var", var);
  if (var == "UnknownVariable"s) {
    lblt.ls_var = static_cast<LineShapeModelVariable>(-1);
  } else {
    lblt.ls_var = to<LineShapeModelVariable>(var);
  }

  tag.add_attribute("ls_coeff", var);
  if (var == "UnknownVariable"s) {
    lblt.ls_coeff = static_cast<LineShapeModelCoefficient>(-1);
  } else {
    lblt.ls_coeff = to<LineShapeModelCoefficient>(var);
  }

  tag.add_attribute("var", var);
  if (var == "UnknownVariable"s) {
    lblt.var = static_cast<LineByLineVariable>(-1);
  } else {
    lblt.var = to<LineByLineVariable>(var);
  }

  xml_read_from_stream(is_xml, lblt.band, pbifs);
  Index line;
  xml_read_from_stream(is_xml, line, pbifs);
  lblt.line = static_cast<Size>(line);
  Index spec;
  xml_read_from_stream(is_xml, spec, pbifs);
  lblt.spec = static_cast<Size>(spec);

  tag.read_from_stream(is_xml);
  tag.check_name("/LblLineKey");
}

//! Write LblLineKey to XML output stream
/*!
  \param os_xml    XML output stream
  \param covmat    LblLineKey
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Unused
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const LblLineKey& lblt,
                         bofstream* pbofs,
                         const String& name [[maybe_unused]]) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  const auto goodenum = [](auto x) -> String {
    return good_enum(x) ? String{toString(x)} : "UnknownVariable"s;
  };

  open_tag.set_name("LblLineKey");
  open_tag.add_attribute("ls_var", goodenum(lblt.ls_var));
  open_tag.add_attribute("ls_coeff", goodenum(lblt.ls_coeff));
  open_tag.add_attribute("var", goodenum(lblt.var));
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  const auto sz = [](Size x) -> Index {
    return (x < std::numeric_limits<Index>::max()) ? static_cast<Index>(x)
                                                   : Index{-1};
  };

  xml_write_to_stream(os_xml, lblt.band, pbofs, "");
  xml_write_to_stream(os_xml, sz(lblt.line), pbofs, "");
  xml_write_to_stream(os_xml, sz(lblt.spec), pbofs, "");

  close_tag.set_name("/LblLineKey");
  close_tag.write_to_stream(os_xml);
}

//=== JacobianTargetType =========================================================

//! Reads JacobianTargetType from XML input stream
/*!
  \param is_xml    XML Input stream
  \param covmat    JacobianTargetType
  \param pbifs     Pointer to binary file stream. NULL for ASCII output.
*/
void xml_read_from_stream(std::istream& is_xml,
                          JacobianTargetType& jtt,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  tag.read_from_stream(is_xml);
  tag.check_name("JacobianTargetType");

  String type;
  tag.get_attribute_value("type", type);
  if (type == "AtmKeyVal"s) {
    jtt.target = AtmKeyVal{};
  } else if (type == "SurfaceKeyVal"s) {
    jtt.target = SurfaceKeyVal{};
  } else if (type == "LblLineKey"s) {
    jtt.target = LblLineKey{};
  } else {
    ARTS_USER_ERROR(R"(Cannot understand type: "{}")", type);
  }

  std::visit(
      [&](auto& data_type) { xml_read_from_stream(is_xml, data_type, pbifs); },
      jtt.target);

  tag.read_from_stream(is_xml);
  tag.check_name("/JacobianTargetType");
}

//! Write JacobianTargetType to XML output stream
/*!
  \param os_xml    XML output stream
  \param covmat    JacobianTargetType
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Unused
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const JacobianTargetType& jtt,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("JacobianTargetType");
  open_tag.add_attribute("type", jtt.type());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  std::visit(
      [&](const auto& data_type) {
        xml_write_to_stream(os_xml, data_type, pbofs, name);
      },
      jtt.target);

  close_tag.set_name("/JacobianTargetType");
  close_tag.write_to_stream(os_xml);
}

//=== JacobianTargetsDiagonalCovarianceMatrixMap =========================================================

//! Reads JacobianTargetsDiagonalCovarianceMatrixMap from XML input stream
/*!
  \param is_xml    XML Input stream
  \param covmat    JacobianTargetsDiagonalCovarianceMatrixMap
  \param pbifs     Pointer to binary file stream. NULL for ASCII output.
*/
void xml_read_from_stream(std::istream& is_xml,
                          JacobianTargetsDiagonalCovarianceMatrixMap& jtmap,
                          bifstream* pbifs) {
  jtmap.map.clear();

  ArtsXMLTag tag;
  tag.read_from_stream(is_xml);
  tag.check_name("JacobianTargetsDiagonalCovarianceMatrixMap");

  Index size;
  tag.get_attribute_value("size", size);

  for (Index i = 0; i < size; i++) {
    JacobianTargetType key;
    BlockMatrix matrix;
    BlockMatrix inverse;

    xml_read_from_stream(is_xml, key, pbifs);
    xml_read_from_stream(is_xml, matrix, pbifs);
    xml_read_from_stream(is_xml, inverse, pbifs);

    jtmap.set(key, matrix, inverse);
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/JacobianTargetsDiagonalCovarianceMatrixMap");
}

//! Write JacobianTargetsDiagonalCovarianceMatrixMap to XML output stream
/*!
  \param os_xml    XML output stream
  \param covmat    JacobianTargetsDiagonalCovarianceMatrixMap
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Unused
*/
void xml_write_to_stream(
    std::ostream& os_xml,
    const JacobianTargetsDiagonalCovarianceMatrixMap& jtmap,
    bofstream* pbofs,
    const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("JacobianTargetsDiagonalCovarianceMatrixMap");
  open_tag.add_attribute("size", static_cast<Index>(jtmap.map.size()));
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (const auto& [key, value] : jtmap) {
    xml_write_to_stream(os_xml, key, pbofs, "key");
    xml_write_to_stream(os_xml, value.first, pbofs, "matrix");
    xml_write_to_stream(os_xml, value.second, pbofs, "inverse");
  }

  close_tag.set_name("/JacobianTargetsDiagonalCovarianceMatrixMap");
  close_tag.write_to_stream(os_xml);
}

//=== DisortSettings =========================================================

//! Reads DisortSettings from XML input stream
/*!
  \param is_xml    XML Input stream
  \param v         DisortSettings
  \param pbifs     Pointer to binary file stream. NULL for ASCII output.
*/
void xml_read_from_stream(std::istream& is_xml,
                          DisortSettings& v,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  tag.read_from_stream(is_xml);
  tag.check_name("DisortSettings");

  xml_read_from_stream(is_xml, v.quadrature_dimension, pbifs);
  xml_read_from_stream(is_xml, v.legendre_polynomial_dimension, pbifs);
  xml_read_from_stream(is_xml, v.fourier_mode_dimension, pbifs);
  xml_read_from_stream(is_xml, v.nfreq, pbifs);
  xml_read_from_stream(is_xml, v.nlay, pbifs);
  xml_read_from_stream(is_xml, v.solar_azimuth_angle, pbifs);
  xml_read_from_stream(is_xml, v.solar_zenith_angle, pbifs);
  xml_read_from_stream(is_xml, v.solar_source, pbifs);
  xml_read_from_stream(
      is_xml, v.bidirectional_reflectance_distribution_functions, pbifs);
  xml_read_from_stream(is_xml, v.optical_thicknesses, pbifs);
  xml_read_from_stream(is_xml, v.single_scattering_albedo, pbifs);
  xml_read_from_stream(is_xml, v.fractional_scattering, pbifs);
  xml_read_from_stream(is_xml, v.source_polynomial, pbifs);
  xml_read_from_stream(is_xml, v.legendre_coefficients, pbifs);
  xml_read_from_stream(is_xml, v.positive_boundary_condition, pbifs);
  xml_read_from_stream(is_xml, v.negative_boundary_condition, pbifs);

  tag.read_from_stream(is_xml);
  tag.check_name("/DisortSettings");
}

//! Write DisortSettings to XML output stream
/*!
  \param os_xml    XML output stream
  \param covmat    DisortSettings
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Unused
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const DisortSettings& v,
                         bofstream* pbofs,
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("DisortSettings");
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(
      os_xml, v.quadrature_dimension, pbofs, "quadrature_dimension");
  xml_write_to_stream(os_xml,
                      v.legendre_polynomial_dimension,
                      pbofs,
                      "legendre_polynomial_dimension");
  xml_write_to_stream(
      os_xml, v.fourier_mode_dimension, pbofs, "fourier_mode_dimension");
  xml_write_to_stream(os_xml, v.nfreq, pbofs, "nfreq");
  xml_write_to_stream(os_xml, v.nlay, pbofs, "nlay");
  xml_write_to_stream(os_xml, v.solar_azimuth_angle, pbofs, "solaz");
  xml_write_to_stream(os_xml, v.solar_zenith_angle, pbofs, "solza");
  xml_write_to_stream(os_xml, v.solar_source, pbofs, "solsrc");
  xml_write_to_stream(os_xml,
                      v.bidirectional_reflectance_distribution_functions,
                      pbofs,
                      "BRDF");
  xml_write_to_stream(os_xml, v.optical_thicknesses, pbofs, "Tau");
  xml_write_to_stream(os_xml, v.single_scattering_albedo, pbofs, "albedo");
  xml_write_to_stream(
      os_xml, v.fractional_scattering, pbofs, "fractional_scattering");
  xml_write_to_stream(os_xml, v.source_polynomial, pbofs, "source_polynomial");
  xml_write_to_stream(
      os_xml, v.legendre_coefficients, pbofs, "legendre_coefficients");
  xml_write_to_stream(os_xml,
                      v.positive_boundary_condition,
                      pbofs,
                      "positive_boundary_condition");
  xml_write_to_stream(os_xml,
                      v.negative_boundary_condition,
                      pbofs,
                      "negative_boundary_condition");

  close_tag.set_name("/DisortSettings");
  close_tag.write_to_stream(os_xml);
}
