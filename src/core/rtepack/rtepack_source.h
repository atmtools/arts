#pragma once

#include "rtepack_propagation_matrix.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
struct SourceVector {
  stokvec_matrix J{};
  stokvec_tensor3 dJ{};

  //! Spectral frequency initialization with NLTE and without Scattering
  void init(const std::span<const propmat_vector> &K,
            const std::span<const propmat_matrix> &dK,
            const std::span<const stokvec_vector> &nlte,
            const std::span<const stokvec_matrix> &dnlte,
            const std::span<const AscendingGrid> &freq_grid,
            const std::span<const Numeric> &ts,
            const Size &it);

  //! Single frequency initialization with NLTE and without Scattering
  void init(const std::span<const propmat> &K,
            const std::span<const propmat_vector> &dK,
            const std::span<const stokvec> &nlte,
            const std::span<const stokvec_vector> &dnlte,
            const std::span<const Numeric> &freq,
            const std::span<const Numeric> &ts,
            const Size &it);

  //! Always call to ensure the object is valid
  void check(Size np, Size nq, Size nf, const std::string_view caller) const;

  //! Should return nf, np, nq [unchecked]
  [[nodiscard]] std::array<Size, 3> shape() const noexcept;
};

void level_nlte(stokvec_vector_view J,
                stokvec_matrix_view dJ,
                const propmat_vector_const_view &K,
                const stokvec_vector_const_view &S,
                const propmat_matrix_const_view &dK,
                const stokvec_matrix_const_view &dS,
                const ConstVectorView &f,
                const Numeric &t,
                const Index &it);
}  // namespace rtepack

using rtepack::SourceVector;

template <>
struct xml_io_stream_name<SourceVector> {
  static constexpr std::string_view name = "SourceVector";
};

template <>
struct xml_io_stream_aggregate<SourceVector> {
  static constexpr bool value = true;
};

template <>
struct std::formatter<SourceVector> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SourceVector &v, FmtContext &ctx) const {
    return tags.format(ctx, v.J, tags.sep(), v.dJ);
  }
};
