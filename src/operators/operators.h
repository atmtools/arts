#pragma once

#include <functional>
#include <ostream>
#include "debug.h"

#include <matpack.h>

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

using NumericUnaryOperator = CustomOperator<Numeric, Numeric>;
