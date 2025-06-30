#include "rtepack_propagation_matrix.h"

namespace rtepack {
propmat_vector operator*(Numeric x, const propmat_vector_const_view &y) {
  propmat_vector z(y.size());
  for (Size i = 0; i < y.size(); ++i) {
    z[i] = x * y[i];
  }
  return z;
}
}  // namespace rtepack

void xml_io_stream<rtepack::propmat>::write(std::ostream &os,
                                            const rtepack::propmat &x,
                                            bofstream *pbofs,
                                            std::string_view name) {
  xml_io_stream<rtepack::vec7>::write(os, x, pbofs, name);
}

void xml_io_stream<rtepack::propmat>::read(std::istream &is,
                                           rtepack::propmat &x,
                                           bifstream *pbifs) {
  xml_io_stream<rtepack::vec7>::read(is, x, pbifs);
}

void xml_io_stream<rtepack::propmat>::put(const rtepack::propmat *const x,
                                          bofstream *pbofs,
                                          Size n) {
  xml_io_stream<rtepack::vec7>::put(x, pbofs, n);
}

void xml_io_stream<rtepack::propmat>::get(rtepack::propmat *x,
                                          bifstream *pbifs,
                                          Size n) {
  xml_io_stream<rtepack::vec7>::get(x, pbifs, n);
}
