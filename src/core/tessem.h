/*!
  \file   tessem.h

  \brief  This file contains functions that are adapted from TESSEM
  code which is used to calculate surface emissivity.
*/

#ifndef tessem_h
#define tessem_h

#include <fstream>

#include "matpack_data.h"

struct TessemNN {
  Index nb_inputs;
  Index nb_outputs;
  Index nb_cache;
  Vector b1;
  Vector b2;
  Matrix w1;
  Matrix w2;
  Vector x_min;
  Vector x_max;
  Vector y_min;
  Vector y_max;

  friend std::ostream& operator<<(std::ostream& os, const TessemNN&) {
    return os;
  }
};

void tessem_read_ascii(std::ifstream& is, TessemNN& net);

void tessem_prop_nn(VectorView ny, const TessemNN& net, ConstVectorView nx);

template <>
struct std::formatter<TessemNN> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const TessemNN& v, FmtContext& ctx) const {
    const std::string_view sep = tags.sep(true);

    tags.add_if_bracket(ctx, '[');
    tags.format(ctx,
                v.nb_inputs,
                sep,
                v.nb_outputs,
                sep,
                v.nb_cache,
                sep,
                v.b1,
                sep,
                v.b2,
                sep,
                v.w1,
                sep,
                v.w2,
                sep,
                v.x_min,
                sep,
                v.x_max,
                sep,
                v.y_min,
                sep,
                v.y_max);
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

#endif /* tessem_h */
