#pragma once

#include <functional>
#include <ostream>

#include "debug.h"

#include <matpack.h>

#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-template-friend"
#endif

template <typename R, typename ... Args>
struct CustomOperator {
  using func_t = std::function<R(Args...)>;
  func_t f;

  friend std::ostream& operator<<(std::ostream& os, const CustomOperator&);

  R operator()(Args... args) const { 
    if (f) return f(args...);
    ARTS_USER_ERROR("CustomOperator not set");
  }
};

#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

using NumericUnaryOperator = CustomOperator<Numeric, Numeric>;

using NumericTernaryOperator = CustomOperator<Numeric, Numeric, Numeric, Numeric>;
