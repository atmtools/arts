#include "rtepack_spectral_matrix.h"

void xml_io_stream<rtepack::specmat>::write(std::ostream &os,
                                            const rtepack::specmat &x,
                                            bofstream *pbofs,
                                            std::string_view name) {
  xml_io_stream<rtepack::cmat44>::write(os, x, pbofs, name);
}

void xml_io_stream<rtepack::specmat>::read(std::istream &is,
                                           rtepack::specmat &x,
                                           bifstream *pbifs) {
  xml_io_stream<rtepack::cmat44>::read(is, x, pbifs);
}

void xml_io_stream<rtepack::specmat>::put(std::span<const rtepack::specmat> x,
                                          bofstream *pbofs) {
  xml_io_stream<rtepack::cmat44>::put(
      std::span{reinterpret_cast<const rtepack::cmat44 *>(x.data()), x.size()},
      pbofs);
}

void xml_io_stream<rtepack::specmat>::get(std::span<rtepack::specmat> x,
                                          bifstream *pbifs) {
  xml_io_stream<rtepack::cmat44>::get(
      std::span{reinterpret_cast<rtepack::cmat44 *>(x.data()), x.size()},
      pbifs);
}

void xml_io_stream<rtepack::specmat>::parse(std::span<rtepack::specmat> x,
                                            std::istream &is) {
  xml_io_stream<rtepack::cmat44>::parse(
      std::span{reinterpret_cast<rtepack::cmat44 *>(x.data()), x.size()}, is);
}
