#include <matpack.h>
#include <xml.h>

#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <tuple>

namespace {
//! Argument type only used here, so the registration below cannot clash with any other
struct Arg {};

using Func = std::function<Numeric(Arg)>;

Numeric one(Arg) { return 1.0; }
Numeric two(Arg) { return 2.0; }
Numeric three(Arg) { return 3.0; }

constexpr Numeric (*one_ptr)(Arg) = &one;
constexpr Numeric (*two_ptr)(Arg) = &two;
}  // namespace

template <> struct xml_io_stream_functional<Func> {
  using func_t    = Numeric (*)(Arg);
  using structs_t = void;
  static constexpr std::array<const func_t*, 2> funcs{&one_ptr, &two_ptr};
};

namespace {
void require(bool condition, std::string_view message) {
  if (not condition) throw std::runtime_error(std::string(message));
}

Func round_trip(const Func& f) {
  std::stringstream ss;
  xml_write_to_stream(ss, f, nullptr);

  Func out;
  xml_read_from_stream(ss, out, nullptr);
  return out;
}
}  // namespace

int main() try {
  require(round_trip(Func{&one})(Arg{}) == 1.0, "Registered function one did not round trip");
  require(round_trip(Func{&two})(Arg{}) == 2.0, "Registered function two did not round trip");
  require(not round_trip(Func{}), "Empty function did not round trip");

  try {
    std::ignore = round_trip(Func{&three});
  } catch (const std::exception&) { return 0; }

  throw std::runtime_error("Unregistered function was written");
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
