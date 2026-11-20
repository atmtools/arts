#pragma once

#include <array.h>
#include <xml.h>

#include "rtepack_propagation_matrix.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
/** The propagation of one point at several frequencies
 *
 * The propagation matrix K and the non-LTE source nlte have one value per
 * frequency.  Their derivatives dK and dnlte have one row per Jacobian target
 * and one column per frequency.
 *
 * nlte_on tells if there is a non-LTE source.  While it is false, nlte and
 * dnlte are all zeros and are not used.  Whoever writes to nlte or dnlte
 * must set it.
 */
struct SpectralPropmat {
  propmat_vector K{};
  propmat_matrix dK{};
  stokvec_vector nlte{};
  stokvec_matrix dnlte{};
  bool           nlte_on{false};

  //! Sizes all to nf frequencies and nq Jacobian targets and sets them to zero, without a non-LTE source
  void init(Size nf, Size nq);

  //! Whether all have the same number of frequencies and the derivatives the same number of Jacobian targets
  [[nodiscard]] bool ok() const;

  //! The number of frequencies and of Jacobian targets [unchecked]
  [[nodiscard]] std::array<Index, 2> shape() const noexcept;

  //! The number of sizes of the shape
  [[nodiscard]] static constexpr Size rank() { return 2; }
};
}  // namespace rtepack

using rtepack::SpectralPropmat;
using ArrayOfSpectralPropmat        = Array<SpectralPropmat>;
using ArrayOfArrayOfSpectralPropmat = Array<ArrayOfSpectralPropmat>;

template <> struct xml_io_stream_name<SpectralPropmat> {
  static constexpr std::string_view name = "SpectralPropmat";
};

template <> struct xml_io_stream_aggregate<SpectralPropmat> {
  static constexpr bool value = true;
};

template <> struct std::formatter<SpectralPropmat> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const SpectralPropmat &v, FmtContext &ctx) const {
    return tags.format(ctx, v.K, tags.sep(), v.dK, tags.sep(), v.nlte, tags.sep(), v.dnlte, tags.sep(), v.nlte_on);
  }
};
