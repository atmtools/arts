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

int main () {
using namespace matpack::md;

  auto x = TMPVector(4);
  const auto y = TMPMatrix(4, 3, 4);
  auto z = TMPTensor3(4, 3, 2, 4);

  z(1,2,1) += 1;
  z(1,1,1) += 1;

  std::cout << "x\n";
  std::cout << x << '\n';
  std::cout << "y\n";
  std::cout << y << '\n';
  std::cout << "z\n";
  std::cout << z << '\n';

 std::cout << type((z(joker, joker, joker))) << '\n';
 std::cout << type((z(joker, joker, 1))) << '\n';
 std::cout << type(z(1, 1, 1)) << '\n';

 std::cout << z(1, 2, joker) << '\n';

 for (auto v: z) v(0, 0) += 5;
  std::cout << "z\n";
  std::cout << z << '\n';

 // This really should not be working ///
 for (auto v: y) v[0] += 3;
  std::cout << "y\n";
  std::cout << y << '\n';
 std::cout << type(y.begin()) << '\n';

  for (auto& v: x) v += 4;
  std::cout << "x\n";
  std::cout << x << '\n';

  std::cout << type(z) << '\n';
  for (auto a: z) {
    std::cout << type(a) << '\n';
    for (auto b: a) {
      std::cout << type(b) << '\n';
      for (auto& c: b) {
        std::cout << type(c) << '\n';
      }
    }
  }
}

