#include "xml_io_stream_matpack_mdspan.h"

void old_xml_io_parse(std::istream& is,
                      Vector& vector,
                      bifstream* pbifs,
                      XMLTag& tag) {
  Index nelem;

  tag.get_attribute_value("nelem", nelem);
  vector.resize(nelem);

  if (pbifs) {
    pbifs->readDoubleArray(vector.data_handle(), vector.size());
  } else {
    for (Index n = 0; n < nelem; n++) {
      is >> double_imanip() >> vector[n];
      if (is.fail()) {
        std::ostringstream os;
        os << " near "
           << "\n  Element: " << n;
        xml_data_parse_error(tag, os.str());
      }
    }
  }
}

void old_xml_io_read(XMLTag& tag,
                     std::istream& is,
                     Vector& vector,
                     bifstream* pbifs) {
  tag.check_name("Vector"sv);

  old_xml_io_parse(is, vector, pbifs, tag);

  tag.read_from_stream(is);
  tag.check_name("/Vector"sv);
}

void old_xml_io_read(XMLTag& tag,
                     std::istream& is,
                     Matrix& matrix,
                     bifstream* pbifs) {
  Index nrows, ncols;

  tag.check_name("Matrix"sv);

  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  matrix.resize(nrows, ncols);

  if (pbifs) {
    pbifs->readDoubleArray(matrix.data_handle(), nrows * ncols);
  } else {
    for (Index r = 0; r < nrows; r++) {
      for (Index c = 0; c < ncols; c++) {
        is >> double_imanip() >> matrix[r, c];
        if (is.fail()) {
          std::ostringstream os;
          os << " near "
             << "\n  Row   : " << r << "\n  Column: " << c;
          xml_data_parse_error(tag, os.str());
        }
      }
    }
  }

  tag.read_from_stream(is);
  tag.check_name("/Matrix"sv);
}

void old_xml_io_read(XMLTag& tag,
                     std::istream& is,
                     Tensor3& tensor,
                     bifstream* pbifs) {
  Index npages, nrows, ncols;

  tag.check_name("Tensor3"sv);

  tag.get_attribute_value("npages", npages);
  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  tensor.resize(npages, nrows, ncols);

  if (pbifs) {
    pbifs->readDoubleArray(tensor.data_handle(), npages * nrows * ncols);
  } else {
    for (Index p = 0; p < npages; p++) {
      for (Index r = 0; r < nrows; r++) {
        for (Index c = 0; c < ncols; c++) {
          is >> double_imanip() >> tensor[p, r, c];
          if (is.fail()) {
            std::ostringstream os;
            os << " near "
               << "\n  Page  : " << p << "\n  Row   : " << r
               << "\n  Column: " << c;
            xml_data_parse_error(tag, os.str());
          }
        }
      }
    }
  }

  tag.read_from_stream(is);
  tag.check_name("/Tensor3"sv);
}

void old_xml_io_read(XMLTag& tag,
                     std::istream& is,
                     Tensor4& tensor,
                     bifstream* pbifs) {
  Index nbooks, npages, nrows, ncols;

  tag.check_name("Tensor4"sv);

  tag.get_attribute_value("nbooks", nbooks);
  tag.get_attribute_value("npages", npages);
  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  tensor = Tensor4(nbooks, npages, nrows, ncols);

  if (pbifs) {
    pbifs->readDoubleArray(tensor.data_handle(),
                           nbooks * npages * nrows * ncols);
  } else {
    for (Index b = 0; b < nbooks; b++) {
      for (Index p = 0; p < npages; p++) {
        for (Index r = 0; r < nrows; r++) {
          for (Index c = 0; c < ncols; c++) {
            is >> double_imanip() >> tensor[b, p, r, c];
            if (is.fail()) {
              std::ostringstream os;
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

  tag.read_from_stream(is);
  tag.check_name("/Tensor4"sv);
}

void old_xml_io_read(XMLTag& tag,
                     std::istream& is,
                     Tensor5& tensor,
                     bifstream* pbifs) {
  Index nshelves, nbooks, npages, nrows, ncols;

  tag.check_name("Tensor5"sv);

  tag.get_attribute_value("nshelves", nshelves);
  tag.get_attribute_value("nbooks", nbooks);
  tag.get_attribute_value("npages", npages);
  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  tensor.resize(nshelves, nbooks, npages, nrows, ncols);

  if (pbifs) {
    pbifs->readDoubleArray(tensor.data_handle(),
                           nshelves * nbooks * npages * nrows * ncols);
  } else {
    for (Index s = 0; s < nshelves; s++) {
      for (Index b = 0; b < nbooks; b++) {
        for (Index p = 0; p < npages; p++) {
          for (Index r = 0; r < nrows; r++) {
            for (Index c = 0; c < ncols; c++) {
              is >> double_imanip() >> tensor[s, b, p, r, c];
              if (is.fail()) {
                std::ostringstream os;
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

  tag.read_from_stream(is);
  tag.check_name("/Tensor5"sv);
}

void old_xml_io_read(XMLTag& tag,
                     std::istream& is,
                     Tensor6& tensor,
                     bifstream* pbifs) {
  Index nvitrines, nshelves, nbooks, npages, nrows, ncols;

  tag.check_name("Tensor6"sv);

  tag.get_attribute_value("nvitrines", nvitrines);
  tag.get_attribute_value("nshelves", nshelves);
  tag.get_attribute_value("nbooks", nbooks);
  tag.get_attribute_value("npages", npages);
  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  tensor.resize(nvitrines, nshelves, nbooks, npages, nrows, ncols);

  if (pbifs) {
    pbifs->readDoubleArray(
        tensor.data_handle(),
        nvitrines * nshelves * nbooks * npages * nrows * ncols);
  } else {
    for (Index v = 0; v < nvitrines; v++) {
      for (Index s = 0; s < nshelves; s++) {
        for (Index b = 0; b < nbooks; b++) {
          for (Index p = 0; p < npages; p++) {
            for (Index r = 0; r < nrows; r++) {
              for (Index c = 0; c < ncols; c++) {
                is >> double_imanip() >> tensor[v, s, b, p, r, c];
                if (is.fail()) {
                  std::ostringstream os;
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

  tag.read_from_stream(is);
  tag.check_name("/Tensor6"sv);
}

void old_xml_io_read(XMLTag& tag,
                     std::istream& is,
                     Tensor7& tensor,
                     bifstream* pbifs) {
  Index nlibraries, nvitrines, nshelves, nbooks, npages, nrows, ncols;

  tag.check_name("Tensor7"sv);

  tag.get_attribute_value("nlibraries", nlibraries);
  tag.get_attribute_value("nvitrines", nvitrines);
  tag.get_attribute_value("nshelves", nshelves);
  tag.get_attribute_value("nbooks", nbooks);
  tag.get_attribute_value("npages", npages);
  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  tensor.resize(nlibraries, nvitrines, nshelves, nbooks, npages, nrows, ncols);

  if (pbifs) {
    pbifs->readDoubleArray(
        tensor.data_handle(),
        nlibraries * nvitrines * nshelves * nbooks * npages * nrows * ncols);
  } else {
    for (Index l = 0; l < nlibraries; l++) {
      for (Index v = 0; v < nvitrines; v++) {
        for (Index s = 0; s < nshelves; s++) {
          for (Index b = 0; b < nbooks; b++) {
            for (Index p = 0; p < npages; p++) {
              for (Index r = 0; r < nrows; r++) {
                for (Index c = 0; c < ncols; c++) {
                  is >> double_imanip() >> tensor[l, v, s, b, p, r, c];
                  if (is.fail()) {
                    std::ostringstream os;
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

  tag.read_from_stream(is);
  tag.check_name("/Tensor7"sv);
}

void xml_parse_from_stream(std::istream& is,
                           ComplexVector& vector,
                           bifstream* pbifs,
                           XMLTag& tag) {
  Index nelem;

  tag.get_attribute_value("nelem", nelem);
  vector.resize(nelem);

  if (pbifs) {
    pbifs->readDoubleArray(reinterpret_cast<double*>(vector.data_handle()),
                           2 * vector.size());
  } else {
    for (Index n = 0; n < nelem; n++) {
      is >> double_imanip() >> real_val(vector[n]) >> imag_val(vector[n]);
      if (is.fail()) {
        std::ostringstream os;
        os << " near "
           << "\n  Element: " << n;
        xml_data_parse_error(tag, os.str());
      }
    }
  }
}

void old_xml_io_read(XMLTag& tag,
                     std::istream& is,
                     ComplexVector& vector,
                     bifstream* pbifs) {
  tag.check_name("ComplexVector"sv);

  xml_parse_from_stream(is, vector, pbifs, tag);

  tag.read_from_stream(is);
  tag.check_name("/ComplexVector"sv);
}
