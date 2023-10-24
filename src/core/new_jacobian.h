#pragma once

#include <atm.h>
#include <matpack.h>

#include <functional>
#include <limits>
#include <numeric>
#include <variant>
#include <vector>
#include "surf.h"

namespace Jacobian {
struct AtmTarget {
  AtmKeyVal type;
  Size x_start{std::numeric_limits<Size>::max()};
  Size x_size{std::numeric_limits<Size>::max()};
  Size target_pos{std::numeric_limits<Size>::max()};
  std::function<void(
      ExhaustiveVectorView, const AtmField&, const ExhaustiveConstVectorView&)>
      set{[](ExhaustiveVectorView xnew,
             const AtmField&,
             const ConstVectorView& xold) { xnew = xold; }};
  std::function<void(
      ExhaustiveVectorView, const AtmField&, const ExhaustiveConstVectorView&)>
      unset{[](VectorView xnew,
             const AtmField&,
             const ExhaustiveConstVectorView& xold) { xnew = xold; }};

  friend std::ostream& operator<<(std::ostream& os, const AtmTarget& target) {
    return os << "Atmosphere key value: " << target.type << ", starting at "
              << target.x_start << " of size " << target.x_size;
  }

  void update(AtmField& atm, const Vector& x) const {
    ARTS_USER_ERROR_IF(static_cast<Size>(x.size()) < (x_start + x_size),
                       "Got too small vector.")

    ARTS_USER_ERROR_IF(not atm.contains(type),
                       "Atmosphere does not contain key value ",
                       type,
                       '.')

    auto xnew = atm[type].flat_view();
    auto xold_d = x.slice(x_start, x_size);
    ARTS_USER_ERROR_IF(
        xnew.size() not_eq xold_d.size(),
        "Problem with sizes.  \n"
        "Did you change your atmosphere since you set the jacobian targets?")

    set(xnew, atm, xold_d);
  }

  void update(Vector& x, const AtmField& atm) const {
    ARTS_USER_ERROR_IF(static_cast<Size>(x.size()) < (x_start + x_size),
                       "Got too small vector.")

    ARTS_USER_ERROR_IF(not atm.contains(type),
                       "Atmosphere does not contain key value ",
                       type,
                       '.')

    auto xnew = x.slice(x_start, x_size);
    auto xold_d = atm[type].flat_view();
    ARTS_USER_ERROR_IF(
        xnew.size() not_eq xold_d.size(),
        "Problem with sizes.  \n"
        "Did you change your atmosphere since you set the jacobian targets?")

    unset(xnew, atm, xold_d);
  }
};


struct SurfaceTarget {
  SurfaceKeyVal type;
  Size x_start{std::numeric_limits<Size>::max()};
  Size x_size{std::numeric_limits<Size>::max()};
  Size target_pos{std::numeric_limits<Size>::max()};
  std::function<void(
      ExhaustiveVectorView, const SurfaceField&, const ExhaustiveConstVectorView&)>
      set{[](ExhaustiveVectorView xnew,
             const SurfaceField&,
             const ConstVectorView& xold) { xnew = xold; }};
  std::function<void(
      ExhaustiveVectorView, const SurfaceField&, const ExhaustiveConstVectorView&)>
      unset{[](VectorView xnew,
             const SurfaceField&,
             const ExhaustiveConstVectorView& xold) { xnew = xold; }};

  friend std::ostream& operator<<(std::ostream& os, const SurfaceTarget& target) {
    return os << "Surface key value: " << target.type << ", starting at "
              << target.x_start << " of size " << target.x_size;
  }

  void update(SurfaceField& surf, const Vector& x) const {
    ARTS_USER_ERROR_IF(static_cast<Size>(x.size()) < (x_start + x_size),
                       "Got too small vector.")

    ARTS_USER_ERROR_IF(not surf.contains(type),
                       "Surface does not contain key value ",
                       type,
                       '.')

    auto xnew = surf[type].flat_view();
    auto xold_d = x.slice(x_start, x_size);
    ARTS_USER_ERROR_IF(
        xnew.size() not_eq xold_d.size(),
        "Problem with sizes.\n"
        "Did you change your surface since you set the jacobian targets?")

    set(xnew, surf, xold_d);
  }

  void update(Vector& x, const SurfaceField& surf) const {
    ARTS_USER_ERROR_IF(static_cast<Size>(x.size()) < (x_start + x_size),
                       "Got too small vector.")

    ARTS_USER_ERROR_IF(not surf.contains(type),
                       "Surface does not contain key value ",
                       type,
                       '.')

    auto xnew = x.slice(x_start, x_size);
    auto xold_d = surf[type].flat_view();
    ARTS_USER_ERROR_IF(
        xnew.size() not_eq xold_d.size(),
        "Problem with sizes.\n"
        "Did you change your surface since you set the jacobian targets?")

    unset(xnew, surf, xold_d);
  }
};

struct Targets {
  std::vector<AtmTarget> atm;
  std::vector<SurfaceTarget> surf;

  //! The length of the x-vector that this would produce
  [[nodiscard]] Size x_size() const {
    const auto sz = [](const auto& x) { return x.x_size; };
    return std::transform_reduce(
        atm.begin(), atm.end(), Size{0}, std::plus<>{}, sz) + std::transform_reduce(
        surf.begin(), surf.end(), Size{0}, std::plus<>{}, sz);
  }

  //! The number of targets
  [[nodiscard]] Size target_count() const {
    return atm.size() + surf.size();
  }

  void throwing_check(Size xsize) const {
    const auto t_size = target_count();

    ARTS_USER_ERROR_IF(xsize not_eq x_size(),
                       "The size of the x-vector does not match the size of the targets.")
    for (auto& a : atm) {
      ARTS_USER_ERROR_IF(a.x_start + a.x_size >= xsize,
                         "The target ",
                         a,
                         " is out of bounds of the x-vector.")
      ARTS_USER_ERROR_IF(t_size >= a.target_pos,
                         "The target ",
                         a,
                         " is out of bounds of the target vector.")
    }

    for (auto& a : surf) {
      ARTS_USER_ERROR_IF(a.x_start + a.x_size >= xsize,
                         "The target ",
                         a,
                         " is out of bounds of the x-vector.")
      ARTS_USER_ERROR_IF(t_size >= a.target_pos,
                         "The target ",
                         a,
                         " is out of bounds of the target vector.")
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const Targets& targets) {
    os << "Jacobian targets:\n";
    for (const auto& t : targets.atm) {
      os << "  " << t << '\n';
    }
    for (const auto& t : targets.surf) {
      os << "  " << t << '\n';
    }
    return os;
  }
};
}  // namespace Jacobian

using JacobianTargets = Jacobian::Targets;
