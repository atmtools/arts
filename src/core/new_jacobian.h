#pragma once

#include <atm.h>
#include <matpack.h>

#include <functional>
#include <numeric>
#include <variant>
#include <vector>

namespace Jacobian {
struct AtmTarget {
  AtmKeyVal type;
  Size start;
  Size size;
  std::function<void(ExhaustiveVectorView, const AtmField&, const ConstVectorView&)> set{
      [](ExhaustiveVectorView xnew,
         const AtmField&,
         const ConstVectorView& xold) { xnew = xold; }};

  friend std::ostream& operator<<(std::ostream& os, const AtmTarget& target) {
    return os << "Atmosphere key value: " << target.type << ", starting at "
              << target.start << " of size " << target.size;
  }

  void update(AtmField& atm, const Vector& x) const {
    ARTS_USER_ERROR_IF(static_cast<Size>(x.size()) < (start + size),
                       "Got too small vector.")

    ARTS_USER_ERROR_IF(not atm.contains(type),
                       "Atmosphere does not contain key value ",
                       type,
                       '.')

    auto xnew = atm[type].flat_view();
    auto xold_d = x[Range(start, size)];
    ARTS_USER_ERROR_IF(
        xnew.size() not_eq xold_d.size(),
        "Problem with sizes.  \n"
        "Did you change your atmosphere since you set the jacobian targets?")

    set(xnew, atm, xold_d);
  }
};

struct Targets {
  std::vector<AtmTarget> atm;

  [[nodiscard]] Size size() const {
    const auto sz = [](const auto& x) { return x.size; };
    return std::transform_reduce(
        atm.begin(), atm.end(), Size{0}, std::plus<>{}, sz);
  }
};
}  // namespace Jacobian

using JacobianTargets = Jacobian::Targets;
