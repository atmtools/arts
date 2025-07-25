#include "rtepack_mueller_matrix.h"

#include <configtypes.h>

namespace rtepack {
void forward_cumulative_transmission(Array<muelmat_vector> &Pi,
                                     const Array<muelmat_vector> &T) {
  const Size N = T.size();

  Pi.resize(N);

  if (N == 0) return;

  const Index nv = T.front().size();
  for (auto &p : Pi) p.resize(nv);

  Pi.front() = T.front();

  for (Size i = 1; i < N; i++) {
    for (auto &&[Pi1, Pi0, T] :
         std::ranges::views::zip(Pi[i], Pi[i - 1], T[i])) {
      Pi1 = Pi0 * T;
    }
  }
}

Array<muelmat_vector> forward_cumulative_transmission(
    const Array<muelmat_vector> &T) {
  Array<muelmat_vector> out;
  forward_cumulative_transmission(out, T);
  return out;
}
}  // namespace rtepack

void xml_io_stream<rtepack::muelmat>::write(std::ostream &os,
                                            const rtepack::muelmat &x,
                                            bofstream *pbofs,
                                            std::string_view name) {
  xml_io_stream<rtepack::mat44>::write(os, x, pbofs, name);
}

void xml_io_stream<rtepack::muelmat>::read(std::istream &is,
                                           rtepack::muelmat &x,
                                           bifstream *pbifs) {
  xml_io_stream<rtepack::mat44>::read(is, x, pbifs);
}

void xml_io_stream<rtepack::muelmat>::put(std::span<const rtepack::muelmat> x,
                                          bofstream *pbofs) {
  xml_io_stream<rtepack::mat44>::put(
      std::span{reinterpret_cast<const rtepack::mat44 *>(x.data()), x.size()},
      pbofs);
}

void xml_io_stream<rtepack::muelmat>::get(std::span<rtepack::muelmat> x,
                                          bifstream *pbifs) {
  xml_io_stream<rtepack::mat44>::get(
      std::span{reinterpret_cast<rtepack::mat44 *>(x.data()), x.size()}, pbifs);
}

void xml_io_stream<rtepack::muelmat>::parse(std::span<rtepack::muelmat> x,
                                            std::istream &is) {
  xml_io_stream<rtepack::mat44>::parse(
      std::span{reinterpret_cast<rtepack::mat44 *>(x.data()), x.size()}, is);
}
