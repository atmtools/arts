#include "matpack_mdspan.h"
#include <iostream>

#include <cstdlib>
#include <cxxabi.h>
#include <memory>
#include <typeinfo>
std::string demangle(const char *name) {
  int status = -4; // some arbitrary value to eliminate the compiler warning

  // enable c++11 by passing the flag -std=c++11 to g++
  std::unique_ptr<char, void (*)(void *)> res{
      abi::__cxa_demangle(name, NULL, NULL, &status), std::free};

  return (status == 0) ? res.get() : name;
}

template <class T> std::string type(const T &t) {
  return demangle(typeid(t).name());
}

int main() {
  using namespace matpack::md;

  auto x = TMPVector(4);
  const auto y = TMPMatrix(4, 3, 4);
  auto yc = TMPMatrix(4, 3, 4);
  auto z = TMPTensor3(4, 3, 2, 4);

  z(1, 2, 1) += 1;
  z(1, 1, 1) += 1;

  for (auto v : z) {
    v(0, 0) += 3;
    for (auto t : v) {
      t[1] += 5;
      for (auto &s : t)
        s *= 2;
    }
  }
  for (auto &v : x)
    v += 2;

  x = std::move(z).flatten();
  std::cout << "SHOULD BE SAME (TOP IS STORED CONST CHAR *):\n14 18 8 18 8 18 "
               "14 18 8 20 8 20 14 18 8 18 8 "
               "18 14 18 8 18 8 18\n";
  std::cout << x << '\n';
  std::cout << z << '\n';
  z = std::move(std::move(x)).reshape(2, 3, 4);
  std::cout << x << '\n';
  std::cout << z << '\n';
  z = std::move(std::move(z).flatten()).reshape(4, 3, 2);
  std::cout << x << '\n';
  std::cout << z << '\n';
  yc = std::move(std::move(z).flatten()).reshape(4, 6);
  std::cout << x << '\n';
  std::cout << yc << '\n';

  z = simple_view<double, 3, true>{std::move(yc).reshape(4,3,2)};
  std::cout << z << '\n';
}
