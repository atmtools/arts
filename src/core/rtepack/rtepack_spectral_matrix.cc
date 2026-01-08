#include "rtepack_spectral_matrix.h"

#include <Faddeeva.hh>

namespace rtepack {
specmat dawson(const specmat &A) {
  using Faddeeva::Dawson;

  return specmat{Dawson(A[0, 0]),
                 Dawson(A[0, 1]),
                 Dawson(A[0, 2]),
                 Dawson(A[0, 3]),
                 Dawson(A[1, 0]),
                 Dawson(A[1, 1]),
                 Dawson(A[1, 2]),
                 Dawson(A[1, 3]),
                 Dawson(A[2, 0]),
                 Dawson(A[2, 1]),
                 Dawson(A[2, 2]),
                 Dawson(A[2, 3]),
                 Dawson(A[3, 0]),
                 Dawson(A[3, 1]),
                 Dawson(A[3, 2]),
                 Dawson(A[3, 3])};
}

specmat erfcx(const specmat &m) {
  using Faddeeva::erfcx;
  return specmat{erfcx(m[0, 0]),
                 erfcx(m[0, 1]),
                 erfcx(m[0, 2]),
                 erfcx(m[0, 3]),
                 erfcx(m[1, 0]),
                 erfcx(m[1, 1]),
                 erfcx(m[1, 2]),
                 erfcx(m[1, 3]),
                 erfcx(m[2, 0]),
                 erfcx(m[2, 1]),
                 erfcx(m[2, 2]),
                 erfcx(m[2, 3]),
                 erfcx(m[3, 0]),
                 erfcx(m[3, 1]),
                 erfcx(m[3, 2]),
                 erfcx(m[3, 3])};
}
}  // namespace rtepack

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
