#pragma once

#include <concepts>
#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

#include "matpack.h"

/** The core concepts of fields
 *
 * A field holds data by key.  The data of a key is one parameter as a
 * function of some coordinates, e.g., altitude, latitude and longitude.
 * However it is stored (on a grid, as a single value, or as a function), the
 * data gives
 *
 * - at(x...): the parameter at the coordinates;
 * - flat_weight(weights, x...): the derivative of that parameter with respect
 *   to the stored values, as the flat index into flat_view() and the weight of
 *   each stored value that contributes (none for a function).  weights is
 *   cleared first, so a caller reuses its allocation over many calls;
 * - flat_view(): the stored values as one vector, which a retrieval reads and
 *   writes;
 * - ok(): whether the stored data is consistent.
 *
 * The data of a spectral field has frequency as an extra, first, coordinate,
 * so its flat weights are those of a frequency and a position.  Given a
 * frequency grid in place of the frequency, at(freq_grid, x...) gives the
 * parameter as one value per frequency.
 *
 * The concepts take the types of the coordinates, so the atmospheric and
 * subsurface fields are fields of Numeric, Numeric, Numeric (altitude,
 * latitude, longitude) and the surface field of Numeric, Numeric (latitude,
 * longitude).  Their spectral variants have the same coordinates besides
 * frequency.
 */

//! The data of a field of the coordinates Coords..., see field_concepts.h
template <typename T, typename... Coords>
concept field_data = requires(const std::remove_cvref_t<T>&           c,
                              std::remove_cvref_t<T>&                 m,
                              std::vector<std::pair<Index, Numeric>>& weights,
                              const Coords&... x) {
  { c.at(x...) } -> std::same_as<Numeric>;
  { c.flat_weight(weights, x...) } -> std::same_as<void>;
  { c.flat_view() } -> std::same_as<ConstVectorView>;
  { m.flat_view() } -> std::same_as<VectorView>;
  { c.ok() } -> std::same_as<bool>;
};

//! The data of a spectral field of the coordinates Coords... besides frequency, see field_concepts.h
template <typename T, typename... Coords>
concept spectral_field_data =
    field_data<T, Numeric, Coords...> and
    requires(const std::remove_cvref_t<T>& c, const AscendingGrid& freq_grid, const Coords&... x) {
      { c.at(freq_grid, x...) } -> std::same_as<Vector>;
    };

//! A field of the coordinates Coords... with data by Key
template <typename T, typename Key, typename... Coords>
concept field = requires(const T& c, T& m, const Key& key) {
  { c[key] } -> field_data<Coords...>;
  { m[key] } -> field_data<Coords...>;
  { c.contains(key) } -> std::same_as<bool>;
  { c.keys() } -> std::ranges::range;
};

//! A spectral field of the coordinates Coords... besides frequency with data by Key
template <typename T, typename Key, typename... Coords>
concept spectral_field = requires(const T& c, T& m, const Key& key) {
  { c[key] } -> spectral_field_data<Coords...>;
  { m[key] } -> spectral_field_data<Coords...>;
  { c.contains(key) } -> std::same_as<bool>;
  { c.keys() } -> std::ranges::range;
};
