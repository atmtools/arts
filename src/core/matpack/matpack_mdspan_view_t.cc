#include "matpack_mdspan_view_t.h"

namespace{
template <typename T>
std::string to_string_impl(const matpack::view_t<const T, 2>& md,
                           format_tags& tags,
                           const std::span<const Size> nl) {
  using std::ranges::views::take, std::ranges::views::drop;

  std::string out;
  out.reserve(md.size() * 10);

  const Size N = nl.size();
  const Size d = tags.bracket ? tags.depth : 0;

  const auto row_sep = [nl](Size i) -> Size {
    for (Size j = 0; j < nl.size(); j++) {
      if ((i % nl[j]) == 0) return j;
    }

    return nl.size();
  };

  const auto add_sep = [&out, sep = tags.sep()]() {
    for (auto c : sep) out.push_back(c);
  };

  const auto add_if_bracket = [&out, test = tags.bracket](char c) {
    if (test) out.push_back(c);
  };

  Size i          = 0;
  const auto view = elemwise_range(md);
  const Size sz   = md.size();

  for (Size j = 0; j < d; j++) add_if_bracket('[');

  if (tags.short_str and sz > 8) {
    out += tags.vformat(view[0]);
    add_sep();
    out += tags.vformat(view[1]);
    add_sep();
    out += tags.vformat(view[2]);
    add_sep();
    out.push_back('.');
    out.push_back('.');
    out.push_back('.');
    add_sep();
    out += tags.vformat(view[sz - 3]);
    add_sep();
    out += tags.vformat(view[sz - 2]);
    add_sep();
    out += tags.vformat(view[sz - 1]);
  } else {
    if (sz > 0) out += tags.vformat(view.front());

    for (auto&& e : view | drop(1)) {
      ++i;

      if (const Size spaces = row_sep(i); spaces < N) {
        for (Size j = spaces; j < N; j++) add_if_bracket(']');
        add_sep();
        out.push_back('\n');
        for (Size j = 0; j < spaces; j++) add_if_bracket(' ');
        for (Size j = spaces; j < N; j++) add_if_bracket('[');
      } else {
        add_sep();
      }

      out += tags.vformat(e);
    }
  }

  for (Size i = 0; i < d; i++) add_if_bracket(']');

  return out;
}
}  // namespace

std::string to_string(const matpack::view_t<const Numeric, 2>& x,
                      format_tags& tags,
                      const std::span<const Size> nl) {
  return to_string_impl(x, tags, nl);
}

std::string to_string(const matpack::view_t<const Complex, 2>& x,
                      format_tags& tags,
                      const std::span<const Size> nl) {
  return to_string_impl(x, tags, nl);
}

std::string to_string(const matpack::view_t<Numeric, 2>& x,
                      format_tags& tags,
                      const std::span<const Size> nl) {
  return to_string(matpack::view_t<const Numeric, 2>{x}, tags, nl);
}

std::string to_string(const matpack::view_t<Complex, 2>& x,
                      format_tags& tags,
                      const std::span<const Size> nl) {
  return to_string(matpack::view_t<const Complex, 2>{x}, tags, nl);
}
