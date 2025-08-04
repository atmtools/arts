#include "rtepack_propagation_matrix.h"

namespace rtepack {
bool propmat::is_polarized() const {
  return std::ranges::any_of(begin() + 1, end(), Cmp::ne(0.0));
};

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

void xml_io_stream<rtepack::propmat>::put(std::span<const rtepack::propmat> x,
                                          bofstream *pbofs) {
  xml_io_stream<rtepack::vec7>::put(
      std::span{reinterpret_cast<const rtepack::vec7 *>(x.data()), x.size()},
      pbofs);
}

void xml_io_stream<rtepack::propmat>::get(std::span<rtepack::propmat> x,
                                          bifstream *pbifs) {
  xml_io_stream<rtepack::vec7>::get(
      std::span{reinterpret_cast<rtepack::vec7 *>(x.data()), x.size()}, pbifs);
}

void xml_io_stream<rtepack::propmat>::parse(std::span<rtepack::propmat> x,
                                            std::istream &is) {
  xml_io_stream<rtepack::vec7>::parse(
      std::span{reinterpret_cast<rtepack::vec7 *>(x.data()), x.size()}, is);
}
