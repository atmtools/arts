#pragma once

/*! An empty class, used for template parameters when something is not needed

Please add more methods as required
*/
struct Empty {
template <typename ...T>
constexpr Empty(T&&...) noexcept {}

template <typename ...T>
constexpr void emplace_back(T&&...) const noexcept {}

template <typename ...T>
constexpr void reserve(T&&...) const noexcept {}
};