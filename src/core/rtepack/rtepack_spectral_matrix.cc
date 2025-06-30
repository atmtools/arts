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

void xml_io_stream<rtepack::specmat>::put(const rtepack::specmat *const x,
                                          bofstream *pbofs,
                                          Size n) {
  xml_io_stream<rtepack::cmat44>::put(x, pbofs, n);
}

void xml_io_stream<rtepack::specmat>::get(rtepack::specmat *x,
                                          bifstream *pbifs,
                                          Size n) {
  xml_io_stream<rtepack::cmat44>::get(x, pbifs, n);
}
