#pragma once

#include <absorptionlines.h>
#include <matpack_concepts.h>

namespace fwd::lbl {
namespace band {
template <typename T>
concept sizeable = requires(const T& t) {
  { t.size() } -> std::convertible_to<std::size_t>;
};

template <typename T>
concept constructable = std::constructible_from<std::remove_cvref_t<T>> and
                        std::constructible_from<std::remove_cvref_t<T>,
                                                Numeric,
                                                Numeric,
                                                SpeciesIsotopologueRatios,
                                                ArrayOfArrayOfSpeciesTag,
                                                Vector,
                                                ArrayOfArrayOfAbsorptionLines>;
}  // namespace band

template <typename T>
concept num_callable = requires(T t, const Numeric& f) {
  { t.at(f) } -> std::convertible_to<Complex>;
};

template <typename T>
concept vec_callable = requires(T t, const Vector& f) {
  { t.at(f) } -> std::convertible_to<ComplexVector>;
} and requires(T t, ExhaustiveComplexVectorView out, const Vector& f) {
  { t.at(out, f) };
};

template <typename T>
concept singleable = num_callable<T>;

template <typename T>
concept bandable = band::sizeable<T> and band::constructable<T> and
                   vec_callable<T> and num_callable<T>;

template <typename T>
concept list_singleable = requires(T t) {
  { t.size() } -> std::convertible_to<std::size_t>;
  { t[0] } -> singleable;
  { t.begin() } -> std::random_access_iterator;
  { t.end() } -> std::random_access_iterator;
};

template <typename T>
concept cutoffable = requires(T t) {
  { T::limit } -> std::convertible_to<Numeric>;
};

constexpr bool has_cutoff(cutoffable auto&& t) {
  return t.limit < std::numeric_limits<Numeric>::max();
}
}  // namespace fwd::lbl
