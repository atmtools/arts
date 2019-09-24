/* Copyright (C) 2003-2012 Oliver Lemke <olemke@core-dump.info>

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
  \file   xml_io_basic_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains basic functions to handle XML data files.

*/

#include "arts.h"
#include "file.h"
#include "xml_io.h"
#include "xml_io_private.h"
#include "xml_io_types.h"

////////////////////////////////////////////////////////////////////////////
//   Overloaded functions for reading/writing data from/to XML stream
////////////////////////////////////////////////////////////////////////////

//=== Index ==================================================================

//! Reads Index from XML input stream
/*!
  \param is_xml  XML Input stream
  \param index   Index return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          Index& index,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("Index");

  if (pbifs) {
    *pbifs >> index;
    if (pbifs->fail()) {
      xml_data_parse_error(tag, "");
    }
  } else {
    is_xml >> index;
    if (is_xml.fail()) {
      xml_data_parse_error(tag, "");
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Index");
}

//! Writes Index to XML output stream
/*!
  \param os_xml  XML Output stream
  \param index   Index value
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const Index& index,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Index");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);

  if (pbofs)
    *pbofs << index;
  else
    os_xml << index;

  close_tag.set_name("/Index");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== Matrix ==========================================================

//! Reads Matrix from XML input stream
/*!
  \param is_xml  XML Input stream
  \param matrix  Matrix return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          Matrix& matrix,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nrows, ncols;

  tag.read_from_stream(is_xml);
  tag.check_name("Matrix");

  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  matrix.resize(nrows, ncols);

  for (Index r = 0; r < nrows; r++) {
    for (Index c = 0; c < ncols; c++) {
      if (pbifs) {
        *pbifs >> matrix(r, c);
        if (pbifs->fail()) {
          ostringstream os;
          os << " near "
             << "\n  Row   : " << r << "\n  Column: " << c;
          xml_data_parse_error(tag, os.str());
        }
      } else {
        is_xml >> double_imanip() >> matrix(r, c);
        if (is_xml.fail()) {
          ostringstream os;
          os << " near "
             << "\n  Row   : " << r << "\n  Column: " << c;
          xml_data_parse_error(tag, os.str());
        }
      }
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Matrix");
}

//! Writes Matrix to XML output stream
/*!
  \param os_xml  XML Output stream
  \param matrix  Matrix
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const Matrix& matrix,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Matrix");
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
      os_xml << matrix(r, 0);

    for (Index c = 1; c < matrix.ncols(); ++c) {
      if (pbofs)
        *pbofs << matrix(r, c);
      else
        os_xml << " " << matrix(r, c);
    }

    if (!pbofs) os_xml << '\n';
  }

  close_tag.set_name("/Matrix");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== Numeric =========================================================

//! Reads Numeric from XML input stream
/*!
  \param is_xml   XML Input stream
  \param numeric  Numeric return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          Numeric& numeric,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("Numeric");

  if (pbifs) {
    *pbifs >> numeric;
    if (pbifs->fail()) {
      xml_data_parse_error(tag, "");
    }
  } else {
    is_xml >> double_imanip() >> numeric;
    if (is_xml.fail()) {
      xml_data_parse_error(tag, "");
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Numeric");
}

//! Writes Numeric to XML output stream
/*!
  \param os_xml   XML Output stream
  \param numeric  Numeric value
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const Numeric& numeric,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Numeric");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);

  xml_set_stream_precision(os_xml);

  if (pbofs)
    *pbofs << numeric;
  else
    os_xml << numeric;

  close_tag.set_name("/Numeric");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== Rational =========================================================

//! Reads Rational from XML input stream
/*!
 * \param is_xml   XML Input stream
 * \param rational  Rational return value
 * \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          Rational& rational,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("Rational");

  if (pbifs) {
    *pbifs >> rational;
    if (pbifs->fail()) {
      xml_data_parse_error(tag, "");
    }
  } else {
    is_xml >> rational;
    if (is_xml.fail()) {
      xml_data_parse_error(tag, "");
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Rational");
}

//! Writes Rational to XML output stream
/*!
 * \param os_xml   XML Output stream
 * \param rational Rational value
 * \param pbofs    Pointer to binary file stream. NULL for ASCII output.
 * \param name     Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const Rational& rational,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Rational");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);

  if (pbofs)
    *pbofs << rational;
  else
    os_xml << rational;

  close_tag.set_name("/Rational");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== Sparse ====================================================

//! Reads Sparse from XML input stream
/*!
  \param is_xml  XML Input stream
  \param sparse  Sparse return value
  \param pbifs   Pointer to binary input stream, NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          Sparse& sparse,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nrows, ncols, nnz;

  tag.read_from_stream(is_xml);
  tag.check_name("Sparse");

  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  sparse.resize(nrows, ncols);

  tag.read_from_stream(is_xml);
  tag.check_name("RowIndex");
  tag.get_attribute_value("nelem", nnz);

  ArrayOfIndex rowind(nnz), colind(nnz);
  Vector data(nnz);

  for (Index i = 0; i < nnz; i++) {
    if (pbifs) {
      *pbifs >> rowind[i];
      if (pbifs->fail()) {
        ostringstream os;
        os << " near "
           << "\n  Row index: " << i;
        xml_data_parse_error(tag, os.str());
      }
    } else {
      is_xml >> rowind[i];
      if (is_xml.fail()) {
        ostringstream os;
        os << " near "
           << "\n  Row index: " << i;
        xml_data_parse_error(tag, os.str());
      }
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/RowIndex");

  tag.read_from_stream(is_xml);
  tag.check_name("ColIndex");

  for (Index i = 0; i < nnz; i++) {
    if (pbifs) {
      *pbifs >> colind[i];
      if (pbifs->fail()) {
        ostringstream os;
        os << " near "
           << "\n  Column index: " << i;
        xml_data_parse_error(tag, os.str());
      }
    } else {
      is_xml >> colind[i];
      if (is_xml.fail()) {
        ostringstream os;
        os << " near "
           << "\n  Column index: " << i;
        xml_data_parse_error(tag, os.str());
      }
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/ColIndex");

  tag.read_from_stream(is_xml);
  tag.check_name("SparseData");

  for (Index i = 0; i < nnz; i++) {
    if (pbifs) {
      *pbifs >> data[i];
      if (pbifs->fail()) {
        ostringstream os;
        os << " near "
           << "\n  Data element: " << i;
        xml_data_parse_error(tag, os.str());
      }
    } else {
      is_xml >> double_imanip() >> data[i];
      if (is_xml.fail()) {
        ostringstream os;
        os << " near "
           << "\n  Data element: " << i;
        xml_data_parse_error(tag, os.str());
      }
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/SparseData");

  tag.read_from_stream(is_xml);
  tag.check_name("/Sparse");

  sparse.insert_elements(nnz, rowind, colind, data);
}

//! Writes Sparse to XML output stream
/*!
  \param os_xml  XML Output stream
  \param sparse  Sparse
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const Sparse& sparse,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag sparse_tag(verbosity);
  ArtsXMLTag row_tag(verbosity);
  ArtsXMLTag col_tag(verbosity);
  ArtsXMLTag data_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  sparse_tag.set_name("Sparse");
  if (name.length()) sparse_tag.add_attribute("name", name);
  sparse_tag.add_attribute("nrows", sparse.nrows());
  sparse_tag.add_attribute("ncols", sparse.ncols());
  //sparse_tag.add_attribute ("nnz", sparse.nnz());
  row_tag.set_name("RowIndex");
  row_tag.add_attribute("nelem", sparse.nnz());
  col_tag.set_name("ColIndex");
  col_tag.add_attribute("nelem", sparse.nnz());
  data_tag.set_name("SparseData");
  data_tag.add_attribute("nelem", sparse.nnz());

  sparse_tag.write_to_stream(os_xml);
  os_xml << '\n';

  row_tag.write_to_stream(os_xml);
  os_xml << '\n';

  ArrayOfIndex rowind(sparse.nnz()), colind(sparse.nnz());
  Vector data(sparse.nnz());
  sparse.list_elements(data, rowind, colind);

  // Write row indices.

  for (Index i = 0; i < sparse.nnz(); i++) {
    if (pbofs)
      //FIXME: It should be the longer lines
      *pbofs << rowind[i];
    else
      os_xml << rowind[i] << '\n';
  }

  close_tag.set_name("/RowIndex");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';

  col_tag.write_to_stream(os_xml);
  os_xml << '\n';

  // Write column indices.

  for (Index i = 0; i < sparse.nnz(); i++) {
    if (pbofs)
      //FIXME: It should be the longer lines
      *pbofs << colind[i];
    else
      os_xml << colind[i] << '\n';
  }

  close_tag.set_name("/ColIndex");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';

  data_tag.write_to_stream(os_xml);
  os_xml << '\n';
  xml_set_stream_precision(os_xml);

  // Write data.

  for (Index i = 0; i < sparse.nnz(); i++) {
    if (pbofs)
      *pbofs << data[i];
    else
      os_xml << data[i] << ' ';
  }
  os_xml << '\n';
  close_tag.set_name("/SparseData");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';

  close_tag.set_name("/Sparse");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== String ===========================================================

//! Reads String from XML input stream
/*!
  \param is_xml  XML Input stream
  \param str     String return value
*/
/*  param pbifs  Pointer to binary input stream. NULL in case of ASCII file.
                 Ignored because strings are always stored in ASCII format. */
void xml_read_from_stream(istream& is_xml,
                          String& str,
                          bifstream* /* pbifs */,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  char dummy;

  tag.read_from_stream(is_xml);
  tag.check_name("String");

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
    xml_parse_error("String must begin with \"");
  }

  //catch case where string is empty. CPD 29/8/05
  dummy = (char)is_xml.peek();
  if (dummy == '"') {
    str = "";
  } else {
    stringbuf strbuf;

    is_xml.get(strbuf, '"');
    if (is_xml.fail()) {
      xml_parse_error("String must end with \"");
    }
    str = strbuf.str();
  }

  // Ignore quote
  is_xml >> dummy;

  tag.read_from_stream(is_xml);
  tag.check_name("/String");
}

//! Writes String to XML output stream
/*!
  \param os_xml  XML Output stream
  \param str     String value
  \param name    Optional name attribute
*/
/*  param pbofs   Pointer to binary file stream. NULL for ASCII output.
                 Ignored because strings are always in ASCII format. */
void xml_write_to_stream(ostream& os_xml,
                         const String& str,
                         bofstream* /* pbofs */,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("String");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);

  os_xml << '\"' << str << '\"';

  close_tag.set_name("/String");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== Tensor3 ================================================================

//! Reads Tensor3 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param tensor  Tensor return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          Tensor3& tensor,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index npages, nrows, ncols;

  tag.read_from_stream(is_xml);
  tag.check_name("Tensor3");

  tag.get_attribute_value("npages", npages);
  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  tensor.resize(npages, nrows, ncols);

  for (Index p = 0; p < npages; p++) {
    for (Index r = 0; r < nrows; r++) {
      for (Index c = 0; c < ncols; c++) {
        if (pbifs) {
          *pbifs >> tensor(p, r, c);
          if (pbifs->fail()) {
            ostringstream os;
            os << " near "
               << "\n  Page  : " << p << "\n  Row   : " << r
               << "\n  Column: " << c;
            xml_data_parse_error(tag, os.str());
          }
        } else {
          is_xml >> double_imanip() >> tensor(p, r, c);
          if (is_xml.fail()) {
            ostringstream os;
            os << " near "
               << "\n  Page  : " << p << "\n  Row   : " << r
               << "\n  Column: " << c;
            xml_data_parse_error(tag, os.str());
          }
        }
      }
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Tensor3");
}

//! Writes Tensor3 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param tensor  Tensor
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const Tensor3& tensor,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Tensor3");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("npages", tensor.npages());
  open_tag.add_attribute("nrows", tensor.nrows());
  open_tag.add_attribute("ncols", tensor.ncols());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  // Write the elements:
  for (Index p = 0; p < tensor.npages(); ++p) {
    for (Index r = 0; r < tensor.nrows(); ++r) {
      if (pbofs)
        *pbofs << tensor(p, r, 0);
      else
        os_xml << tensor(p, r, 0);
      for (Index c = 1; c < tensor.ncols(); ++c) {
        if (pbofs)
          *pbofs << tensor(p, r, c);
        else
          os_xml << " " << tensor(p, r, c);
      }
      if (!pbofs) os_xml << '\n';
    }
  }

  close_tag.set_name("/Tensor3");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== Tensor4 =========================================================

//! Reads Tensor4 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param tensor  Tensor return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          Tensor4& tensor,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nbooks, npages, nrows, ncols;

  tag.read_from_stream(is_xml);
  tag.check_name("Tensor4");

  tag.get_attribute_value("nbooks", nbooks);
  tag.get_attribute_value("npages", npages);
  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  tensor.resize(nbooks, npages, nrows, ncols);

  for (Index b = 0; b < nbooks; b++) {
    for (Index p = 0; p < npages; p++) {
      for (Index r = 0; r < nrows; r++) {
        for (Index c = 0; c < ncols; c++) {
          if (pbifs) {
            *pbifs >> tensor(b, p, r, c);
            if (pbifs->fail()) {
              ostringstream os;
              os << " near "
                 << "\n  Book  : " << b << "\n  Page  : " << p
                 << "\n  Row   : " << r << "\n  Column: " << c;
              xml_data_parse_error(tag, os.str());
            }
          } else {
            is_xml >> double_imanip() >> tensor(b, p, r, c);
            if (is_xml.fail()) {
              ostringstream os;
              os << " near "
                 << "\n  Book  : " << b << "\n  Page  : " << p
                 << "\n  Row   : " << r << "\n  Column: " << c;
              xml_data_parse_error(tag, os.str());
            }
          }
        }
      }
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Tensor4");
}

//! Writes Tensor4 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param tensor  Tensor
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const Tensor4& tensor,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Tensor4");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("nbooks", tensor.nbooks());
  open_tag.add_attribute("npages", tensor.npages());
  open_tag.add_attribute("nrows", tensor.nrows());
  open_tag.add_attribute("ncols", tensor.ncols());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  // Write the elements:
  for (Index b = 0; b < tensor.nbooks(); ++b) {
    for (Index p = 0; p < tensor.npages(); ++p) {
      for (Index r = 0; r < tensor.nrows(); ++r) {
        if (pbofs)
          *pbofs << tensor(b, p, r, 0);
        else
          os_xml << tensor(b, p, r, 0);
        for (Index c = 1; c < tensor.ncols(); ++c) {
          if (pbofs)
            *pbofs << tensor(b, p, r, c);
          else
            os_xml << " " << tensor(b, p, r, c);
        }
        if (!pbofs) os_xml << '\n';
      }
    }
  }

  close_tag.set_name("/Tensor4");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== Tensor5 =========================================================

//! Reads Tensor5 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param tensor  Tensor return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          Tensor5& tensor,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nshelves, nbooks, npages, nrows, ncols;

  tag.read_from_stream(is_xml);
  tag.check_name("Tensor5");

  tag.get_attribute_value("nshelves", nshelves);
  tag.get_attribute_value("nbooks", nbooks);
  tag.get_attribute_value("npages", npages);
  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  tensor.resize(nshelves, nbooks, npages, nrows, ncols);

  for (Index s = 0; s < nshelves; s++) {
    for (Index b = 0; b < nbooks; b++) {
      for (Index p = 0; p < npages; p++) {
        for (Index r = 0; r < nrows; r++) {
          for (Index c = 0; c < ncols; c++) {
            if (pbifs) {
              *pbifs >> tensor(s, b, p, r, c);
              if (pbifs->fail()) {
                ostringstream os;
                os << " near "
                   << "\n  Shelf : " << s << "\n  Book  : " << b
                   << "\n  Page  : " << p << "\n  Row   : " << r
                   << "\n  Column: " << c;
                xml_data_parse_error(tag, os.str());
              }
            } else {
              is_xml >> double_imanip() >> tensor(s, b, p, r, c);
              if (is_xml.fail()) {
                ostringstream os;
                os << " near "
                   << "\n  Shelf : " << s << "\n  Book  : " << b
                   << "\n  Page  : " << p << "\n  Row   : " << r
                   << "\n  Column: " << c;
                xml_data_parse_error(tag, os.str());
              }
            }
          }
        }
      }
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Tensor5");
}

//! Writes Tensor5 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param tensor  Tensor
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const Tensor5& tensor,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Tensor5");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("nshelves", tensor.nshelves());
  open_tag.add_attribute("nbooks", tensor.nbooks());
  open_tag.add_attribute("npages", tensor.npages());
  open_tag.add_attribute("nrows", tensor.nrows());
  open_tag.add_attribute("ncols", tensor.ncols());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  // Write the elements:
  for (Index s = 0; s < tensor.nshelves(); ++s) {
    for (Index b = 0; b < tensor.nbooks(); ++b) {
      for (Index p = 0; p < tensor.npages(); ++p) {
        for (Index r = 0; r < tensor.nrows(); ++r) {
          if (pbofs)
            *pbofs << tensor(s, b, p, r, 0);
          else
            os_xml << tensor(s, b, p, r, 0);
          for (Index c = 1; c < tensor.ncols(); ++c) {
            if (pbofs)
              *pbofs << tensor(s, b, p, r, c);
            else
              os_xml << " " << tensor(s, b, p, r, c);
          }
          if (!pbofs) os_xml << '\n';
        }
      }
    }
  }

  close_tag.set_name("/Tensor5");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== Tensor6 =========================================================

//! Reads Tensor6 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param tensor  Tensor return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          Tensor6& tensor,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nvitrines, nshelves, nbooks, npages, nrows, ncols;

  tag.read_from_stream(is_xml);
  tag.check_name("Tensor6");

  tag.get_attribute_value("nvitrines", nvitrines);
  tag.get_attribute_value("nshelves", nshelves);
  tag.get_attribute_value("nbooks", nbooks);
  tag.get_attribute_value("npages", npages);
  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  tensor.resize(nvitrines, nshelves, nbooks, npages, nrows, ncols);

  for (Index v = 0; v < nvitrines; v++) {
    for (Index s = 0; s < nshelves; s++) {
      for (Index b = 0; b < nbooks; b++) {
        for (Index p = 0; p < npages; p++) {
          for (Index r = 0; r < nrows; r++) {
            for (Index c = 0; c < ncols; c++) {
              if (pbifs) {
                *pbifs >> tensor(v, s, b, p, r, c);
                if (pbifs->fail()) {
                  ostringstream os;
                  os << " near "
                     << "\n  Vitrine: " << v << "\n  Shelf  : " << s
                     << "\n  Book   : " << b << "\n  Page   : " << p
                     << "\n  Row    : " << r << "\n  Column : " << c;
                  xml_data_parse_error(tag, os.str());
                }
              } else {
                is_xml >> double_imanip() >> tensor(v, s, b, p, r, c);
                if (is_xml.fail()) {
                  ostringstream os;
                  os << " near "
                     << "\n  Vitrine: " << v << "\n  Shelf  : " << s
                     << "\n  Book   : " << b << "\n  Page   : " << p
                     << "\n  Row    : " << r << "\n  Column : " << c;
                  xml_data_parse_error(tag, os.str());
                }
              }
            }
          }
        }
      }
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Tensor6");
}

//! Writes Tensor6 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param tensor  Tensor
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const Tensor6& tensor,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Tensor6");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("nvitrines", tensor.nvitrines());
  open_tag.add_attribute("nshelves", tensor.nshelves());
  open_tag.add_attribute("nbooks", tensor.nbooks());
  open_tag.add_attribute("npages", tensor.npages());
  open_tag.add_attribute("nrows", tensor.nrows());
  open_tag.add_attribute("ncols", tensor.ncols());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  // Write the elements:
  for (Index v = 0; v < tensor.nvitrines(); ++v) {
    for (Index s = 0; s < tensor.nshelves(); ++s) {
      for (Index b = 0; b < tensor.nbooks(); ++b) {
        for (Index p = 0; p < tensor.npages(); ++p) {
          for (Index r = 0; r < tensor.nrows(); ++r) {
            if (pbofs)
              *pbofs << tensor(v, s, b, p, r, 0);
            else
              os_xml << tensor(v, s, b, p, r, 0);
            for (Index c = 1; c < tensor.ncols(); ++c) {
              if (pbofs)
                *pbofs << tensor(v, s, b, p, r, c);
              else
                os_xml << " " << tensor(v, s, b, p, r, c);
            }
            if (!pbofs) os_xml << '\n';
          }
        }
      }
    }
  }

  close_tag.set_name("/Tensor6");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== Tensor7 =========================================================

//! Reads Tensor7 from XML input stream
/*!
  \param is_xml  XML Input stream
  \param tensor  Tensor return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          Tensor7& tensor,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nlibraries, nvitrines, nshelves, nbooks, npages, nrows, ncols;

  tag.read_from_stream(is_xml);
  tag.check_name("Tensor7");

  tag.get_attribute_value("nlibraries", nlibraries);
  tag.get_attribute_value("nvitrines", nvitrines);
  tag.get_attribute_value("nshelves", nshelves);
  tag.get_attribute_value("nbooks", nbooks);
  tag.get_attribute_value("npages", npages);
  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  tensor.resize(nlibraries, nvitrines, nshelves, nbooks, npages, nrows, ncols);

  for (Index l = 0; l < nlibraries; l++) {
    for (Index v = 0; v < nvitrines; v++) {
      for (Index s = 0; s < nshelves; s++) {
        for (Index b = 0; b < nbooks; b++) {
          for (Index p = 0; p < npages; p++) {
            for (Index r = 0; r < nrows; r++) {
              for (Index c = 0; c < ncols; c++) {
                if (pbifs) {
                  *pbifs >> tensor(l, v, s, b, p, r, c);
                  if (pbifs->fail()) {
                    ostringstream os;
                    os << " near "
                       << "\n  Library: " << l << "\n  Vitrine: " << v
                       << "\n  Shelf  : " << s << "\n  Book   : " << b
                       << "\n  Page   : " << p << "\n  Row    : " << r
                       << "\n  Column : " << c;
                    xml_data_parse_error(tag, os.str());
                  }
                } else {
                  is_xml >> double_imanip() >> tensor(l, v, s, b, p, r, c);
                  if (is_xml.fail()) {
                    ostringstream os;
                    os << " near "
                       << "\n  Library: " << l << "\n  Vitrine: " << v
                       << "\n  Shelf  : " << s << "\n  Book   : " << b
                       << "\n  Page   : " << p << "\n  Row    : " << r
                       << "\n  Column : " << c;
                    xml_data_parse_error(tag, os.str());
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Tensor7");
}

//! Writes Tensor7 to XML output stream
/*!
  \param os_xml  XML Output stream
  \param tensor  Tensor
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const Tensor7& tensor,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Tensor7");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("nlibraries", tensor.nlibraries());
  open_tag.add_attribute("nvitrines", tensor.nvitrines());
  open_tag.add_attribute("nshelves", tensor.nshelves());
  open_tag.add_attribute("nbooks", tensor.nbooks());
  open_tag.add_attribute("npages", tensor.npages());
  open_tag.add_attribute("nrows", tensor.nrows());
  open_tag.add_attribute("ncols", tensor.ncols());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  // Write the elements:
  for (Index l = 0; l < tensor.nlibraries(); ++l) {
    for (Index v = 0; v < tensor.nvitrines(); ++v) {
      for (Index s = 0; s < tensor.nshelves(); ++s) {
        for (Index b = 0; b < tensor.nbooks(); ++b) {
          for (Index p = 0; p < tensor.npages(); ++p) {
            for (Index r = 0; r < tensor.nrows(); ++r) {
              if (pbofs)
                *pbofs << tensor(l, v, s, b, p, r, 0);
              else
                os_xml << tensor(l, v, s, b, p, r, 0);
              for (Index c = 1; c < tensor.ncols(); ++c) {
                if (pbofs)
                  *pbofs << tensor(l, v, s, b, p, r, c);
                else
                  os_xml << " " << tensor(l, v, s, b, p, r, c);
              }
              if (!pbofs) os_xml << '\n';
            }
          }
        }
      }
    }
  }

  close_tag.set_name("/Tensor7");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== Vector ==========================================================

//! Parses Vector from XML input stream
/*!
  \param is_xml  XML Input stream
  \param vector  Vector return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
  \param tag     XML tag object
*/
void xml_parse_from_stream(istream& is_xml,
                           Vector& vector,
                           bifstream* pbifs,
                           ArtsXMLTag& tag,
                           const Verbosity&) {
  Index nelem;

  tag.get_attribute_value("nelem", nelem);
  vector.resize(nelem);

  for (Index n = 0; n < nelem; n++) {
    if (pbifs) {
      *pbifs >> vector[n];
      if (pbifs->fail()) {
        ostringstream os;
        os << " near "
           << "\n  Element: " << n;
        xml_data_parse_error(tag, os.str());
      }
    } else {
      is_xml >> double_imanip() >> vector[n];
      if (is_xml.fail()) {
        ostringstream os;
        os << " near "
           << "\n  Element: " << n;
        xml_data_parse_error(tag, os.str());
      }
    }
  }
}

//! Reads Vector from XML input stream
/*!
  \param is_xml  XML Input stream
  \param vector  Vector return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          Vector& vector,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("Vector");

  xml_parse_from_stream(is_xml, vector, pbifs, tag, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/Vector");
}

//! Writes Vector to XML output stream
/*!
  \param os_xml  XML Output stream
  \param vector  Vector
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const Vector& vector,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);
  Index n = vector.nelem();
  ostringstream v;

  // Convert nelem to string
  v << n;

  open_tag.set_name("Vector");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("nelem", v.str());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (Index i = 0; i < n; ++i)
    if (pbofs)
      *pbofs << vector[i];
    else
      os_xml << vector[i] << '\n';

  close_tag.set_name("/Vector");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== TransmissionMatrix ================================================================

//! Reads TransmissionMatrix from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param tm      TransmissionMatrix return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          TransmissionMatrix& tm,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index stokes_dim, nf;

  tag.read_from_stream(is_xml);
  tag.check_name("TransmissionMatrix");

  tag.get_attribute_value("Stokes", stokes_dim);
  tag.get_attribute_value("Freqs", nf);
  tm = TransmissionMatrix(nf, stokes_dim);
  if (pbifs) {
    *pbifs >> tm;
    if (pbifs->fail()) {
      ostringstream os;
      os << "TransmissionMatrix has wrong dimensions";
      xml_data_parse_error(tag, os.str());
    }
  } else {
    is_xml >> tm;
    if (is_xml.fail()) {
      ostringstream os;
      os << "TransmissionMatrix has wrong dimensions";
      xml_data_parse_error(tag, os.str());
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/TransmissionMatrix");
}

//! Writes RadiationVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param tm      TransmissionMatrix
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const TransmissionMatrix& tm,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("TransmissionMatrix");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("Stokes", tm.StokesDim());
  open_tag.add_attribute("Freqs", tm.Frequencies());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  if (pbofs)
    *pbofs << tm;
  else
    os_xml << tm;

  close_tag.set_name("/TransmissionMatrix");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== RadiationVector ================================================================

//! Reads RadiationVector from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param rv      RadiationVector return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          RadiationVector& rv,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index stokes_dim, nf;

  tag.read_from_stream(is_xml);
  tag.check_name("RadiationVector");

  tag.get_attribute_value("Stokes", stokes_dim);
  tag.get_attribute_value("Freqs", nf);
  rv = RadiationVector(nf, stokes_dim);
  if (pbifs) {
    *pbifs >> rv;
    if (pbifs->fail()) {
      ostringstream os;
      os << "RadiationVector has wrong dimensions";
      xml_data_parse_error(tag, os.str());
    }
  } else {
    is_xml >> rv;
    if (is_xml.fail()) {
      ostringstream os;
      os << "RadiationVector has wrong dimensions";
      xml_data_parse_error(tag, os.str());
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/RadiationVector");
}

//! Writes RadiationVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param rv      RadiationVector
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const RadiationVector& rv,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("RadiationVector");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("Stokes", rv.StokesDim());
  open_tag.add_attribute("Freqs", rv.Frequencies());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  if (pbofs)
    *pbofs << rv;
  else
    os_xml << rv;

  close_tag.set_name("/RadiationVector");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== AbsorptionLines ================================================================

//! Reads AbsorptionLines from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param al      AbsorptionLines return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          AbsorptionLines& al,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  
  tag.read_from_stream(is_xml);
  tag.check_name("AbsorptionLines");
  
  // Number of lines
  Index nlines;
  tag.get_attribute_value("nlines", nlines);
  
  // Species of the lines
  SpeciesTag spec;
  tag.get_attribute_value("species", spec);
  
  // Cutoff type
  String s_cutoff;
  tag.get_attribute_value("cutofftype", s_cutoff);
  const Absorption::CutoffType cutoff = Absorption::string2cutofftype(s_cutoff);
  
  // Mirroring type
  String s_mirroring;
  tag.get_attribute_value("mirroringtype", s_mirroring);
  const Absorption::MirroringType mirroring = Absorption::string2mirroringtype(s_mirroring);
  
  // Line population type
  String s_population;
  tag.get_attribute_value("populationtype", s_population);
  const Absorption::PopulationType population = Absorption::string2populationtype(s_population);
  
  // Normalization type
  String s_normalization;
  tag.get_attribute_value("normalizationtype", s_normalization);
  const Absorption::NormalizationType normalization = Absorption::string2normalizationtype(s_normalization);
  
  // Normalization type
  String s_lineshapetype;
  tag.get_attribute_value("lineshapetype", s_lineshapetype);
  const LineShape::Type lineshapetype = LineShape::string2shapetype(s_lineshapetype);
  
  /** Reference temperature for all parameters of the lines */
  Numeric T0;
  tag.get_attribute_value("T0", T0);
  
  /** cutoff frequency */
  Numeric cutofffreq;
  tag.get_attribute_value("cutofffreq", cutofffreq);
  
  /** linemixing limit */
  Numeric linemixinglimit;
  tag.get_attribute_value("linemixinglimit", linemixinglimit);
  
  /** List of local quantum numbers, these must be defined */
  std::vector<QuantumNumberType> localquanta;
  tag.get_attribute_value("localquanta", localquanta);
  
  /** Catalog ID */  
  QuantumNumbers upperglobalquanta;
  tag.get_attribute_value("upperglobalquanta", upperglobalquanta);
  QuantumNumbers lowerglobalquanta;
  tag.get_attribute_value("lowerglobalquanta", lowerglobalquanta);
  
  /** A list of broadening species */
  ArrayOfSpeciesTag broadeningspecies;
  bool selfbroadening;
  bool bathbroadening;
  tag.get_attribute_value("broadeningspecies", broadeningspecies, selfbroadening, bathbroadening);
  
  String temperaturemodes;
  tag.get_attribute_value("temperaturemodes", temperaturemodes);
  auto metamodel = LineShape::MetaData2ModelShape(temperaturemodes);

  al = AbsorptionLines(selfbroadening, bathbroadening,
                       nlines, cutoff, mirroring,
                       population, normalization,
                       lineshapetype, T0, cutofffreq,
                       linemixinglimit, QuantumIdentifier(spec.Species(),
                       spec.Isotopologue(), upperglobalquanta, lowerglobalquanta),
                       localquanta, broadeningspecies, metamodel);
  
  if (pbifs) {
//     *pbifs >> al;
//     if (pbifs->fail()) {
      ostringstream os;
      os << "AbsorptionLines cannot do binary IO yet";
      xml_data_parse_error(tag, os.str());
//     }
  } else {
    is_xml >> al;
    if (is_xml.fail()) {
      ostringstream os;
      os << "AbsorptionLines has wrong dimensions";
      xml_data_parse_error(tag, os.str());
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/AbsorptionLines");
}

//! Writes AbsorptionLines to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param al      AbsorptionLines
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(ostream& os_xml,
                         const AbsorptionLines& al,
                         bofstream* pbofs,
                         const String&,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_comment_tag(verbosity);
  ArtsXMLTag close_comment_tag(verbosity);
  open_comment_tag.set_name("comment");
  open_comment_tag.write_to_stream(os_xml);
  os_xml << al.MetaData();
  close_comment_tag.set_name("/comment");
  close_comment_tag.write_to_stream(os_xml);
  os_xml << '\n';

  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("AbsorptionLines");
  open_tag.add_attribute("nlines", al.NumLines());
  open_tag.add_attribute("species", al.SpeciesName());
  open_tag.add_attribute("cutofftype", Absorption::cutofftype2string(al.Cutoff()));
  open_tag.add_attribute("mirroringtype", Absorption::mirroringtype2string(al.Mirroring()));
  open_tag.add_attribute("populationtype", Absorption::populationtype2string(al.Population()));
  open_tag.add_attribute("normalizationtype", Absorption::normalizationtype2string(al.Normalization()));
  open_tag.add_attribute("lineshapetype", LineShape::shapetype2string(al.LineShapeType()));
  open_tag.add_attribute("T0", al.T0());
  open_tag.add_attribute("cutofffreq", al.CutoffFreqValue());
  open_tag.add_attribute("linemixinglimit", al.LinemixingLimit());
  open_tag.add_attribute("localquanta", al.LocalQuanta());
  open_tag.add_attribute("upperglobalquanta", al.UpperQuantumNumbers());
  open_tag.add_attribute("lowerglobalquanta", al.LowerQuantumNumbers());
  open_tag.add_attribute("broadeningspecies", al.BroadeningSpecies(), al.Self(), al.Bath());
  open_tag.add_attribute("temperaturemodes", al.LineShapeMetaData());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  if (pbofs) {
    *pbofs << al;
    throw std::runtime_error("No binary AbsorptionLines IO yet");
  }
  else
    os_xml << al;

  close_tag.set_name("/AbsorptionLines");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

// FIXME: These should be implemented, sooner or later...

void xml_read_from_stream(istream&,
                          Timer&,
                          bifstream* /* pbifs */,
                          const Verbosity&) {
  throw runtime_error("Method not implemented!");
}

void xml_write_to_stream(ostream&,
                         const Timer&,
                         bofstream* /* pbofs */,
                         const String& /* name */,
                         const Verbosity&) {
  throw runtime_error("Method not implemented!");
}
