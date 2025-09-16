#include "functional_numeric_ternary.h"

#include "functional_atm.h"
#include "functional_atm_field.h"
#include "functional_atm_mag_field.h"

namespace ternary {
namespace {
template <typename T>
concept jacable = requires(T a) {
  { a.x() } -> std::same_as<VectorView>;
} and requires(const T a) {
  { a.x() } -> std::same_as<ConstVectorView>;
} and requires(const T a) {
  {
    a.w(Numeric{}, Numeric{}, Numeric{})
  } -> std::same_as<std::vector<std::pair<Index, Numeric>>>;
};

using var_t = xml_io_stream_functional<NumericTernary>::structs_t;

template <Size I = 0>
ConstVectorView get_x_tmpl(const NumericTernary& f) {
  if constexpr (I < std::variant_size_v<var_t>) {
    using T = std::variant_alternative_t<I, var_t>;

    if constexpr (jacable<T>) {
      if (const T* ptr = f.template target<T>(); ptr != nullptr)
        return ptr->x();
    }

    return get_x_tmpl<I + 1>(f);
  } else {
    return ConstVectorView{};
  }
}

template <Size I = 0>
VectorView get_x_tmpl(NumericTernary& f) {
  if constexpr (I < std::variant_size_v<var_t>) {
    using T = std::variant_alternative_t<I, var_t>;

    if constexpr (jacable<T>) {
      if (T* ptr = f.template target<T>(); ptr != nullptr) return ptr->x();
    }

    return get_x_tmpl<I + 1>(f);
  } else {
    return VectorView{};
  }
}

template <Size I = 0>
std::vector<std::pair<Index, Numeric>> get_w_tmpl(const NumericTernary& f,
                                                  Numeric alt [[maybe_unused]],
                                                  Numeric lat [[maybe_unused]],
                                                  Numeric lon
                                                  [[maybe_unused]]) {
  if constexpr (I < std::variant_size_v<var_t>) {
    using T = std::variant_alternative_t<I, var_t>;

    if constexpr (jacable<T>) {
      if (const T* ptr = f.template target<T>(); ptr != nullptr)
        return ptr->w(alt, lat, lon);
    }

    return get_w_tmpl<I + 1>(f, alt, lat, lon);
  } else {
    return {};
  }
}

static_assert(jacable<Atm::SchmidthLegendre>);
static_assert(jacable<Atm::MagnitudeField>);
static_assert(jacable<Atm::External>);
}  // namespace

ConstVectorView get_xc(const NumericTernary& f) { return get_x_tmpl(f); }

VectorView get_xm(NumericTernary& f) { return get_x_tmpl(f); }

std::vector<std::pair<Index, Numeric>> get_w(const NumericTernary& f,
                                             Numeric alt,
                                             Numeric lat,
                                             Numeric lon) {
  return get_w_tmpl(f, alt, lat, lon);
}

bool is_field_function(const NumericTernary& f) {
  if (const auto& ptr = f.target<Atm::External>(); ptr != nullptr)
    return ptr->is_field_function;

  return f.target<Atm::SchmidthLegendre>() != nullptr ||
         f.target<Atm::MagnitudeField>() != nullptr;
}
}  // namespace ternary
