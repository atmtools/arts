#pragma once

#include "rtepack_mueller_matrix.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_spectral_matrix.h"

namespace rtepack {
struct TransmittanceMatrix {
  TransmittanceOption option{};

  muelmat_matrix T{};
  muelmat_matrix L{};
  muelmat_matrix P{};

  muelmat_tensor4 dT{};
  muelmat_tensor4 dL{};

  void init(const std::span<const propmat_vector> &K,
            const std::span<const propmat_matrix> &dK,
            const ConstVectorView &r,
            const ConstTensor3View &dr,
            const TransmittanceOption);

  void init(const std::span<const propmat> &K,
            const std::span<const propmat_vector> &dK,
            const ConstVectorView &r,
            const ConstTensor3View &dr,
            const TransmittanceOption);

  void check(Size np, Size nq, Size nf, const std::string_view caller) const;

  //! Should return nf, np, nq [unchecked]
  [[nodiscard]] std::array<Size, 3> shape() const noexcept;

 private:
  void constant(const std::span<const propmat_vector> &K,
                const std::span<const propmat_matrix> &dK,
                const ConstVectorView &r,
                const ConstTensor3View &dr);

  void linsrc(const std::span<const propmat_vector> &K,
              const std::span<const propmat_matrix> &dK,
              const ConstVectorView &r,
              const ConstTensor3View &dr);

  void linprop(const std::span<const propmat_vector> &K,
               const std::span<const propmat_matrix> &dK,
               const ConstVectorView &r,
               const ConstTensor3View &dr);

  void constant(const std::span<const propmat> &K,
                const std::span<const propmat_vector> &dK,
                const ConstVectorView &r,
                const ConstTensor3View &dr);

  void linsrc(const std::span<const propmat> &K,
              const std::span<const propmat_vector> &dK,
              const ConstVectorView &r,
              const ConstTensor3View &dr);

  void linprop(const std::span<const propmat> &K,
               const std::span<const propmat_vector> &dK,
               const ConstVectorView &r,
               const ConstTensor3View &dr);
};

struct tran {
  Numeric a, exp_a;
  Numeric b, c, d, u, v, w;
  Numeric b2, c2, d2, u2, v2, w2;
  Numeric B, C, S;
  Numeric x2, y2, x, y, cy, sy, cx, sx;
  Numeric ix, iy, inv_x2y2;
  Numeric C0, C1, C2, C3;
  bool polarized, x_zero, y_zero, both_zero, either_zero;

  constexpr tran() = default;

  tran(const propmat &k1, const propmat &k2, const Numeric r);

  [[nodiscard]] muelmat operator()() const noexcept;
  [[nodiscard]] muelmat expm1() const noexcept;
  [[nodiscard]] muelmat linsrc() const noexcept;
  [[nodiscard]] muelmat linsrc_linprop(const muelmat &t,
                                       const propmat &k1,
                                       const propmat &k2,
                                       const Numeric r) const noexcept;

  [[nodiscard]] muelmat deriv(const muelmat &t,
                              const propmat &k1,
                              const propmat &k2,
                              const propmat &dk,
                              const Numeric r,
                              const Numeric dr) const;
  [[nodiscard]] muelmat linsrc_deriv(const propmat &dk,
                                     const Numeric r,
                                     const Numeric dr) const;
  [[nodiscard]] muelmat linsrc_linprop_deriv(const muelmat &lambda,
                                             const muelmat &t,
                                             const propmat &k1,
                                             const propmat &k2,
                                             const propmat &dk_in,
                                             const muelmat &dt,
                                             const Numeric r,
                                             const Numeric dr,
                                             bool k1_deriv) const;
};

muelmat exp(propmat k, Numeric r = 1.0);

propmat logK(const muelmat &m);

specmat sqrt(const propmat &pm);
}  // namespace rtepack

using rtepack::TransmittanceMatrix;

template <>
struct xml_io_stream_name<TransmittanceMatrix> {
  static constexpr std::string_view name = "TransmittanceMatrix";
};

template <>
struct xml_io_stream_aggregate<TransmittanceMatrix> {
  static constexpr bool value = true;
};

template <>
struct std::formatter<TransmittanceMatrix> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const TransmittanceMatrix &v,
                              FmtContext &ctx) const {
    const std::string_view sep = tags.sep();
    return tags.format(
        ctx, v.option, sep, v.T, sep, v.L, sep, v.P, sep, v.dT, sep, v.dL);
  }
};
