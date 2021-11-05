#include "partfun.h"

#include <algorithm>
#include <cmath>

#include "lineshape.h"
#include "physics_funcs.h"

#include <Faddeeva/Faddeeva.hh>

using Constant::inv_pi;
using Constant::inv_sqrt_pi;
using Constant::pi;
using Constant::pow2;
using Constant::pow3;
using Constant::pow4;
using Constant::sqrt_ln_2;
using Constant::sqrt_pi;

constexpr Numeric ln_16 =
    2.772588722239781237668928485832706272302000537441021016482720037973574487879;

namespace LineShape {
Complex Doppler::operator()(Numeric f) noexcept {
  x = (f - mF0) * invGD;
  F = invGD * inv_sqrt_pi * std::exp(-pow2(x));
  return F;
}

Complex Voigt::operator()(Numeric f) noexcept {
  real_val(z) = invGD * (f - mF0);
  F = inv_sqrt_pi * invGD * Faddeeva::w(z);
  dF = 2 * invGD * (Complex(0, inv_pi * invGD) - z * F);
  return F;
}

SpeedDependentVoigt::SpeedDependentVoigt(Numeric F0_noshift, const Output &ls,
                                         Numeric GD_div_F0, Numeric dZ) noexcept
    : mF0(F0_noshift + dZ + ls.D0 - 1.5 * ls.D2),
      invGD(sqrt_ln_2 / nonstd::abs(GD_div_F0 * mF0)),
      invc2(1.0 / Complex(ls.G2, ls.D2)), dx(Complex(ls.G0 - 1.5 * ls.G2, mF0)),
      x(dx * invc2), sqrty(invc2 / (2 * invGD)),
      calcs(init(Complex(ls.G2, ls.D2))) {
  calc();
}

Complex SpeedDependentVoigt::dFdf() const noexcept {
  switch (calcs) {
  case CalcType::Full:
    return invGD * invc2 * (dw1 - dw2) / (2 * sqrt_pi * sq);
  case CalcType::Voigt:
    return dw1 * pow2(invGD) * inv_sqrt_pi;
  case CalcType::LowXandHighY:
    return dw1 * pow2(invGD) * inv_sqrt_pi -
           dw2 * invGD * invc2 / (2 * sqrt_pi * sq);
  case CalcType::LowYandLowX:
    return pow2(invc2) * (-dw1 * sq + 1i * w1) / (sqrt_pi * sq);
  case CalcType::LowYandHighX:
    return 1i * pow2(invc2) * (x - 3) / (pi * pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdF0() const noexcept {
  switch (calcs) {
  case CalcType::Full:
    return (4 * pow2(invGD) * (-w1 + w2) * sq +
            1i * invc2 *
                (dw1 * (Complex(0, 2 * mF0 * pow2(invGD)) - 2 * invGD * sq +
                        invc2) -
                 dw2 * (Complex(0, 2 * mF0 * pow2(invGD)) + 2 * invGD * sq +
                        invc2))) /
           (4 * sqrt_pi * invGD * mF0 * sq);
  case CalcType::Voigt:
    return -invGD * (Complex(0, invGD) * dw1 * (dx - Complex(0, mF0)) + w1) /
           (sqrt_pi * mF0);
  case CalcType::LowXandHighY:
    return (4 * pow2(invGD) * (-w1 + w2) * sq -
            1i * (4 * dw1 * pow3(invGD) * (dx - Complex(0, mF0)) * sq +
                  dw2 * invc2 *
                      (Complex(0, 2 * mF0 * pow2(invGD)) + 2 * invGD * sq +
                       invc2))) /
           (4 * sqrt_pi * invGD * mF0 * sq);
  case CalcType::LowYandLowX:
    return pow2(invc2) * (dw1 * sq - 1i * w1) / (sqrt_pi * sq);
  case CalcType::LowYandHighX:
    return 1i * pow2(invc2) * (3 - x) / (pi * pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdD0(Numeric dD0dD0) const noexcept {
  switch (calcs) {
  case CalcType::Full:
    return -dD0dD0 *
           (4 * pow2(invGD) * (w1 - w2) * sq +
            1i * invc2 *
                (-dw1 * (Complex(0, 2 * mF0 * pow2(invGD)) - 2 * invGD * sq +
                         invc2) +
                 dw2 * (Complex(0, 2 * mF0 * pow2(invGD)) + 2 * invGD * sq +
                        invc2))) /
           (4 * sqrt_pi * invGD * mF0 * sq);
  case CalcType::Voigt:
    return -dD0dD0 * invGD *
           (Complex(0, invGD) * dw1 * (dx - Complex(0, mF0)) + w1) /
           (sqrt_pi * mF0);
  case CalcType::LowXandHighY:
    return -dD0dD0 *
           (4 * pow2(invGD) * (w1 - w2) * sq +
            1i * (4 * dw1 * pow3(invGD) * (dx - Complex(0, mF0)) * sq +
                  dw2 * invc2 *
                      (Complex(0, 2 * mF0 * pow2(invGD)) + 2 * invGD * sq +
                       invc2))) /
           (4 * sqrt_pi * invGD * mF0 * sq);
  case CalcType::LowYandLowX:
    return dD0dD0 * pow2(invc2) * (dw1 * sq - 1i * w1) / (sqrt_pi * sq);
  case CalcType::LowYandHighX:
    return -Complex(0, dD0dD0) * pow2(invc2) * (x - 3) / (pi * pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdG0(Numeric dG0dG0) const noexcept {
  switch (calcs) {
  case CalcType::Full:
    return Complex(0, dG0dG0) * invGD * invc2 * (dw1 - dw2) /
           (2 * sqrt_pi * sq);
  case CalcType::Voigt:
    return Complex(0, dG0dG0) * dw1 * pow2(invGD) * inv_sqrt_pi;
  case CalcType::LowXandHighY:
    return Complex(0, dG0dG0) * invGD * (2 * dw1 * invGD * sq - dw2 * invc2) /
           (2 * sqrt_pi * sq);
  case CalcType::LowYandLowX:
    return -dG0dG0 * pow2(invc2) * (1i * dw1 * sq + w1) / (sqrt_pi * sq);
  case CalcType::LowYandHighX:
    return -dG0dG0 * pow2(invc2) * (x - 3) / (pi * pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdD2(Numeric dD2dD2) const noexcept {
  switch (calcs) {
  case CalcType::Full:
    return dD2dD2 *
           (12 * pow2(invGD) * (w1 - w2) * sq +
            1i * invc2 *
                (dw1 * (-Complex(0, 2 * mF0 * pow2(invGD)) *
                            (2 * dx * invc2 + 3) +
                        4 * Complex(0, invGD) * invc2 * mF0 * sq +
                        6 * invGD * sq - Complex(0, 2 * mF0) * pow2(invc2) -
                        3 * invc2) +
                 dw2 *
                     (Complex(0, 2 * mF0 * pow2(invGD)) * (2 * dx * invc2 + 3) +
                      4 * Complex(0, invGD) * invc2 * mF0 * sq +
                      6 * invGD * sq + Complex(0, 2 * mF0) * pow2(invc2) +
                      3 * invc2))) /
           (8 * sqrt_pi * invGD * mF0 * sq);
  case CalcType::Voigt:
    return 3 * dD2dD2 * invGD *
           (Complex(0, invGD) * dw1 * (dx - Complex(0, mF0)) + w1) /
           (2 * sqrt_pi * mF0);
  case CalcType::LowXandHighY:
    return dD2dD2 *
           (12 * pow2(invGD) * (w1 - w2) * sq +
            1i *
                (12 * dw1 * pow3(invGD) * (dx - Complex(0, mF0)) * sq +
                 dw2 * invc2 *
                     (Complex(0, 2 * mF0 * pow2(invGD)) * (2 * dx * invc2 + 3) +
                      4 * Complex(0, invGD) * invc2 * mF0 * sq +
                      6 * invGD * sq + Complex(0, 2 * mF0) * pow2(invc2) +
                      3 * invc2))) /
           (8 * sqrt_pi * invGD * mF0 * sq);
  case CalcType::LowYandLowX:
    return dD2dD2 * pow2(invc2) *
           (4 * 1i * sq * (sqrt_pi * w1 * sq - 1) -
            sqrt_pi * (dw1 * sq - 1i * w1) * (2 * dx * invc2 + 3)) /
           (2 * pi * sq);
  case CalcType::LowYandHighX:
    return Complex(0, dD2dD2) * pow2(invc2) *
           (-x * (2 * x - 3) + (x - 3) * (2 * dx * invc2 + 3)) /
           (2 * pi * pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdG2(Numeric dG2dG2) const noexcept {
  switch (calcs) {
  case CalcType::Full:
    return Complex(0, dG2dG2) * invc2 *
           (dw1 * (-pow2(invGD) * (2 * dx * invc2 + 3) +
                   2 * invGD * invc2 * sq - pow2(invc2)) +
            dw2 * (pow2(invGD) * (2 * dx * invc2 + 3) + 2 * invGD * invc2 * sq +
                   pow2(invc2))) /
           (4 * sqrt_pi * invGD * sq);
  case CalcType::Voigt:
    return -3 * Complex(0, dG2dG2) * dw1 * pow2(invGD) / (2 * sqrt_pi);
  case CalcType::LowXandHighY:
    return Complex(0, dG2dG2) *
           (-6 * dw1 * pow3(invGD) * sq +
            dw2 * invc2 *
                (pow2(invGD) * (2 * dx * invc2 + 3) + 2 * invGD * invc2 * sq +
                 pow2(invc2))) /
           (4 * sqrt_pi * invGD * sq);
  case CalcType::LowYandLowX:
    return dG2dG2 * pow2(invc2) *
           (4 * sq * (sqrt_pi * w1 * sq - 1) +
            sqrt_pi * (2 * dx * invc2 + 3) * (1i * dw1 * sq + w1)) /
           (2 * pi * sq);
  case CalcType::LowYandHighX:
    return dG2dG2 * pow2(invc2) *
           (-x * (2 * x - 3) + (x - 3) * (2 * dx * invc2 + 3)) /
           (2 * pi * pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdH(Numeric dZ) const noexcept {
  switch (calcs) {
  case CalcType::Full:
    return -dZ *
           (4 * pow2(invGD) * (w1 - w2) * sq +
            1i * invc2 *
                (-dw1 * (Complex(0, 2 * mF0 * pow2(invGD)) - 2 * invGD * sq +
                         invc2) +
                 dw2 * (Complex(0, 2 * mF0 * pow2(invGD)) + 2 * invGD * sq +
                        invc2))) /
           (4 * sqrt_pi * invGD * mF0 * sq);
  case CalcType::Voigt:
    return -dZ * invGD *
           (Complex(0, invGD) * dw1 * (dx - Complex(0, mF0)) + w1) /
           (sqrt_pi * mF0);
  case CalcType::LowXandHighY:
    return -dZ *
           (4 * pow2(invGD) * (w1 - w2) * sq +
            1i * (4 * dw1 * pow3(invGD) * (dx - Complex(0, mF0)) * sq +
                  dw2 * invc2 *
                      (Complex(0, 2 * mF0 * pow2(invGD)) + 2 * invGD * sq +
                       invc2))) /
           (4 * sqrt_pi * invGD * mF0 * sq);
  case CalcType::LowYandLowX:
    return dZ * pow2(invc2) * (dw1 * sq - 1i * w1) / (sqrt_pi * sq);
  case CalcType::LowYandHighX:
    return -Complex(0, dZ) * pow2(invc2) * (x - 3) / (pi * pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdVMR(const Output &d) const noexcept {
  switch (calcs) {
  case CalcType::Full:
    return (-4 * pow2(invGD) * (2 * d.D0 - 3 * d.D2) * (w1 - w2) * sq +
            1i * invc2 *
                (dw1 * (-2 * pow2(invGD) * mF0 *
                            (Complex(3 * d.G2 - 2 * d.G0, 3 * d.D2 - 2 * d.D0) +
                             2 * dx * invc2 * Complex(d.G2, d.D2)) +
                        4 * invGD * invc2 * mF0 * sq * Complex(d.G2, d.D2) -
                        2 * invGD * (2 * d.D0 - 3 * d.D2) * sq -
                        2 * pow2(invc2) * mF0 * Complex(d.G2, d.D2) +
                        invc2 * (2 * d.D0 - 3 * d.D2)) -
                 dw2 * (-2 * pow2(invGD) * mF0 *
                            (Complex(3 * d.G2 - 2 * d.G0, 3 * d.D2 - 2 * d.D0) +
                             2 * dx * invc2 * Complex(d.G2, d.D2)) -
                        4 * invGD * invc2 * mF0 * sq * Complex(d.G2, d.D2) +
                        2 * invGD * (2 * d.D0 - 3 * d.D2) * sq -
                        2 * pow2(invc2) * mF0 * Complex(d.G2, d.D2) +
                        invc2 * (2 * d.D0 - 3 * d.D2)))) /
           (8 * sqrt_pi * invGD * mF0 * sq);
  case CalcType::Voigt:
    return -invGD *
           (Complex(0, invGD) * dw1 *
                (dx * (2 * d.D0 - 3 * d.D2) -
                 mF0 * Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2)) +
            w1 * (2 * d.D0 - 3 * d.D2)) /
           (2 * sqrt_pi * mF0);
  case CalcType::LowXandHighY:
    return -(4 * pow2(invGD) * (2 * d.D0 - 3 * d.D2) * (w1 - w2) * sq +
             1i * (4 * dw1 * pow3(invGD) * sq *
                       (dx * (2 * d.D0 - 3 * d.D2) -
                        mF0 *
                            Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2)) +
                   dw2 * invc2 *
                       (2 * pow2(invGD) * mF0 *
                            (Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) -
                             2 * dx * invc2 * Complex(d.G2, d.D2)) -
                        4 * invGD * invc2 * mF0 * sq * Complex(d.G2, d.D2) +
                        2 * invGD * (2 * d.D0 - 3 * d.D2) * sq -
                        2 * pow2(invc2) * mF0 * Complex(d.G2, d.D2) +
                        invc2 * (2 * d.D0 - 3 * d.D2)))) /
           (8 * sqrt_pi * invGD * mF0 * sq);
  case CalcType::LowYandLowX:
    return pow2(invc2) *
           (4 * sq * Complex(d.G2, d.D2) * (sqrt_pi * w1 * sq - 1) -
            sqrt_pi * (1i * dw1 * sq + w1) *
                (Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) -
                 2 * dx * invc2 * Complex(d.G2, d.D2))) /
           (2 * pi * sq);
  case CalcType::LowYandHighX:
    return -pow2(invc2) *
           (x * (2 * x - 3) * Complex(d.G2, d.D2) +
            (x - 3) * (Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) -
                       2 * dx * invc2 * Complex(d.G2, d.D2))) /
           (2 * pi * pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdT(const Output &d, Numeric T) const noexcept {
  switch (calcs) {
  case CalcType::Full:
    return (-pow2(invGD) * (w1 - w2) * sq *
                (T * (2 * d.D0 - 3 * d.D2) * sqrt_ln_2 - pow3(invGD) * mF0) *
                ln_16 +
            1i * invc2 *
                (dw1 *
                     (T * invGD * invc2 * mF0 * sq * Complex(d.G2, d.D2) *
                          ln_16 -
                      2 * invGD * sq *
                          (T * (2 * d.D0 - 3 * d.D2) * sqrt_ln_2 -
                           pow3(invGD) * mF0) *
                          sqrt_ln_2 -
                      (2 * T * pow2(invGD) * mF0 *
                           (Complex(3 * d.G2 - 2 * d.G0, 3 * d.D2 - 2 * d.D0) +
                            2 * dx * invc2 * Complex(d.G2, d.D2)) *
                           sqrt_ln_2 +
                       2 * T * pow2(invc2) * mF0 * Complex(d.G2, d.D2) *
                           sqrt_ln_2 +
                       invc2 * (T * (-2 * d.D0 + 3 * d.D2) * sqrt_ln_2 +
                                pow3(invGD) * mF0)) *
                          sqrt_ln_2) +
                 dw2 *
                     (T * invGD * invc2 * mF0 * sq * Complex(d.G2, d.D2) *
                          ln_16 +
                      2 * invGD * sq *
                          (T * (-2 * d.D0 + 3 * d.D2) * sqrt_ln_2 +
                           pow3(invGD) * mF0) *
                          sqrt_ln_2 +
                      (-2 * T * pow2(invGD) * mF0 *
                           (Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) -
                            2 * dx * invc2 * Complex(d.G2, d.D2)) *
                           sqrt_ln_2 +
                       2 * T * pow2(invc2) * mF0 * Complex(d.G2, d.D2) *
                           sqrt_ln_2 +
                       invc2 * (T * (-2 * d.D0 + 3 * d.D2) * sqrt_ln_2 +
                                pow3(invGD) * mF0)) *
                          sqrt_ln_2)) *
                sqrt_ln_2) /
           (8 * sqrt_pi * T * invGD * mF0 * sq * pow3(sqrt_ln_2));
  case CalcType::Voigt:
    return -invGD *
           (-Complex(0, invGD) * dw1 *
                (T * mF0 * Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) *
                     sqrt_ln_2 -
                 dx * (T * (2 * d.D0 - 3 * d.D2) * sqrt_ln_2 -
                       pow3(invGD) * mF0)) +
            w1 * (T * (2 * d.D0 - 3 * d.D2) * sqrt_ln_2 - pow3(invGD) * mF0)) /
           (2 * sqrt_pi * T * mF0 * sqrt_ln_2);
  case CalcType::LowXandHighY:
    return (-pow2(invGD) * (w1 - w2) * sq *
                (T * (2 * d.D0 - 3 * d.D2) * sqrt_ln_2 - pow3(invGD) * mF0) *
                ln_16 +
            1i * (dw1 * pow3(invGD) * sq *
                      (T * mF0 *
                           Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) *
                           sqrt_ln_2 -
                       dx * (T * (2 * d.D0 - 3 * d.D2) * sqrt_ln_2 -
                             pow3(invGD) * mF0)) *
                      ln_16 +
                  dw2 * invc2 *
                      (T * invGD * invc2 * mF0 * sq * Complex(d.G2, d.D2) *
                           ln_16 +
                       2 * invGD * sq *
                           (T * (-2 * d.D0 + 3 * d.D2) * sqrt_ln_2 +
                            pow3(invGD) * mF0) *
                           sqrt_ln_2 +
                       (-2 * T * pow2(invGD) * mF0 *
                            (Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) -
                             2 * dx * invc2 * Complex(d.G2, d.D2)) *
                            sqrt_ln_2 +
                        2 * T * pow2(invc2) * mF0 * Complex(d.G2, d.D2) *
                            sqrt_ln_2 +
                        invc2 * (T * (-2 * d.D0 + 3 * d.D2) * sqrt_ln_2 +
                                 pow3(invGD) * mF0)) *
                           sqrt_ln_2) *
                      sqrt_ln_2)) /
           (8 * sqrt_pi * T * invGD * mF0 * sq * pow3(sqrt_ln_2));
  case CalcType::LowYandLowX:
    return pow2(invc2) *
           (4 * sq * Complex(d.G2, d.D2) * (sqrt_pi * w1 * sq - 1) -
            sqrt_pi * (1i * dw1 * sq + w1) *
                (Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) -
                 2 * dx * invc2 * Complex(d.G2, d.D2))) /
           (2 * pi * sq);
  case CalcType::LowYandHighX:
    return -pow2(invc2) *
           (x * (2 * x - 3) * Complex(d.G2, d.D2) +
            (x - 3) * (Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) -
                       2 * dx * invc2 * Complex(d.G2, d.D2))) /
           (2 * pi * pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::operator()(Numeric f) noexcept {
  imag_val(dx) = mF0 - f;
  x = dx * invc2;
  update_calcs();
  calc();
  return F;
}

constexpr Numeric abs_squared(Complex z) noexcept {
  return pow2(z.real()) + pow2(z.imag());
}

SpeedDependentVoigt::CalcType
SpeedDependentVoigt::init(const Complex c2) const noexcept {
  if (abs_squared(c2) == 0)
    return CalcType::Voigt;
  if (abs_squared(x) <= 9e-16 * abs_squared(sqrty * sqrty))
    return CalcType::LowXandHighY;
  if ((abs_squared(sqrty * sqrty) <= 1.e-30 * abs_squared(x)) and abs_squared(std::sqrt(x)) <= 16.e6)
    return CalcType::LowYandLowX; // Weird case, untested
  if ((abs_squared(sqrty * sqrty) <= 1.e-30 * abs_squared(x)))
    return CalcType::LowYandHighX;
  return CalcType::Full;
}

void SpeedDependentVoigt::update_calcs() noexcept {
  if (calcs not_eq CalcType::Voigt)
    calcs = init(Complex(1, 1));
}

void SpeedDependentVoigt::calc() noexcept {
  switch (calcs) {
  case CalcType::Full:
    sq = std::sqrt(x + sqrty * sqrty);
    w1 = Faddeeva::w(1i * (sq - sqrty));
    w2 = Faddeeva::w(1i * (sq + sqrty));
    F = inv_sqrt_pi * invGD * (w1 - w2);
    dw1 = 2i * (inv_sqrt_pi - (sq - sqrty) * w1);
    dw2 = 2i * (inv_sqrt_pi - (sq + sqrty) * w2);
    break;
  case CalcType::Voigt:
    w1 = Faddeeva::w(1i * dx * invGD);
    F = inv_sqrt_pi * invGD * w1;
    dw1 = 2i * (inv_sqrt_pi - dx * invGD * w1);
    break;
  case CalcType::LowXandHighY:
    sq = std::sqrt(x + sqrty * sqrty);
    w1 = Faddeeva::w(1i * dx * invGD);
    w2 = Faddeeva::w(1i * (sq + sqrty));
    F = inv_sqrt_pi * invGD * (w1 - w2);
    dw1 = 2i * (inv_sqrt_pi - dx * invGD * w1);
    dw2 = 2i * (inv_sqrt_pi - (sq + sqrty) * w2);
    break;
  case CalcType::LowYandLowX:
    sq = std::sqrt(x);
    w1 = Faddeeva::w(1i * sq);
    F = 2 * inv_pi * invc2 * (1 - sqrt_pi * sq * w1);
    dw1 = 2i * (inv_sqrt_pi - sq * w1);
    break;
  case CalcType::LowYandHighX:
    F = inv_pi * invc2 * (1 / x - 1.5 / pow2(x));
    break;
  }
}

HartmannTran::HartmannTran(Numeric F0_noshift, const Output &ls,
                           Numeric GD_div_F0, Numeric dZ) noexcept
    : G0(ls.G0), D0(ls.D0), G2(ls.G2), D2(ls.D2), FVC(ls.FVC), ETA(ls.ETA),
      mF0(F0_noshift + dZ + (1 - ls.ETA) * (ls.D0 - 1.5 * ls.D2)),
      invGD(sqrt_ln_2 / nonstd::abs(GD_div_F0 * mF0)),
      deltax(ls.FVC + (1 - ls.ETA) * (ls.G0 - 3 * ls.G2 / 2), mF0),
      sqrty(1 / (2 * (1 - ls.ETA) * Complex(ls.G2, ls.D2) * invGD)) {
  calc();
}

Complex HartmannTran::dFdf() const noexcept {
  constexpr Complex ddeltax = -1i;
  Complex dx = -ddeltax / ((ETA - 1) * Complex(G2, D2));
  Complex dsqrtxy = dx / (2 * sqrtxy);

  switch (calcs) {
  case CalcType::Full: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dsqrtxy;
    Complex dA = Complex(0, sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
    Complex dB =
        sqrt_pi *
        ((pow2(z1) - 1) * 1i * dw1 * dz1 - (pow2(z2) - 1) * 1i * dw2 * dz2 +
         2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) /
        (2 * sqrty * (ETA - 1) * Complex(G2, D2));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tLowZ: {
    Complex dz1 = invGD * ddeltax;
    Complex dA = Complex(0, sqrt_pi * invGD) * dw1 * dz1;
    Complex dB =
        -invGD *
        (sqrt_pi * ((pow2(z1) - 1) * 1i * dw1 * dz1 + 2 * w1 * z1 * dz1) - dz1);
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tHighZ: {
    Complex dz1 = invGD * ddeltax;
    Complex dA = Complex(0, sqrt_pi * invGD) * dw1 * dz1;
    Complex dB = Complex(0, sqrt_pi * invGD) * dw1 * dz1 -
                 invGD * dz1 / (2 * pow2(z1)) +
                 9 * invGD * dz1 / (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowXandHighY: {
    Complex dz1 = invGD * ddeltax;
    Complex dz2 = dsqrtxy;
    Complex dA = Complex(0, sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
    Complex dB = Complex(0, sqrt_pi * invGD) * dw1 * dz1 -
                 invGD * dz1 / (2 * pow2(z1)) +
                 9 * invGD * dz1 / (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandLowX: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dx / (2 * sqrtx);
    Complex dA = 2 * sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) /
                 ((ETA - 1) * Complex(G2, D2));
    Complex dB =
        -(2 * sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) *
              (2 * pow2(sqrty) + x - 1) +
          2 * sqrt_pi * w1 * dz1 + Complex(0, 2 * sqrt_pi) * z1 * dw1 * dz1 +
          2 * (sqrt_pi * w2 * z2 - 1) * dx) /
        ((ETA - 1) * Complex(G2, D2));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandHighX: {
    Complex dz1 = dsqrtxy;
    Complex dA = (x - 3) * dx / ((ETA - 1) * Complex(G2, D2) * pow3(x));
    Complex dB =
        (-2 * sqrt_pi * (w1 * dz1 + z1 * 1i * dw1 * dz1) * pow3(x) -
         (x - 3) * (2 * pow2(sqrty) + x - 1) * dx + (2 * x - 3) * x * dx / 2) /
        ((ETA - 1) * Complex(G2, D2) * pow3(x));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  }
  return {};
}

Complex HartmannTran::dFdF0() const noexcept {
  Numeric dGD = (1 / (invGD * mF0));
  Numeric dinvGD = -dGD * pow2(invGD);
  Complex dsqrty = dinvGD / (2 * (ETA - 1) * Complex(G2, D2) * pow2(invGD));
  constexpr Complex ddeltax = 1i;
  Complex dx = -ddeltax / ((ETA - 1) * Complex(G2, D2));
  Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;

  switch (calcs) {
  case CalcType::Full: {
    Complex dz1 = dsqrtxy - dsqrty;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB =
        sqrt_pi *
        ((-(pow2(z1) - 1) * w1 + (pow2(z2) - 1) * w2) * dsqrty +
         ((pow2(z1) - 1) * 1i * dw1 * dz1 - (pow2(z2) - 1) * 1i * dw2 * dz2 +
          2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) *
             sqrty) /
        (2 * (ETA - 1) * Complex(G2, D2) * pow2(sqrty));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tLowZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        -(sqrt_pi * ((pow2(z1) - 1) * 1i * dw1 * dz1 + 2 * w1 * z1 * dz1) -
          dz1) *
            invGD -
        (sqrt_pi * (pow2(z1) - 1) * w1 - z1) * dinvGD;
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tHighZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowXandHighY: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandLowX: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dx / (2 * sqrtx);
    Complex dA = 2 * sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) /
                 ((ETA - 1) * Complex(G2, D2));
    Complex dB =
        -(2 * sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) *
              (2 * pow2(sqrty) + x - 1) +
          2 * sqrt_pi * w1 * dz1 + Complex(0, 2 * sqrt_pi) * z1 * dw1 * dz1 +
          2 * (4 * sqrty * dsqrty + dx) * (sqrt_pi * w2 * z2 - 1)) /
        ((ETA - 1) * Complex(G2, D2));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandHighX: {
    Complex dz1 = dsqrtxy;
    Complex dA = (x - 3) * dx / ((ETA - 1) * Complex(G2, D2) * pow3(x));
    Complex dB = (-2 * sqrt_pi * (w1 * dz1 + z1 * 1i * dw1 * dz1) * pow3(x) +
                  (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x / 2 -
                  (x - 3) * (2 * pow2(sqrty) + x - 1) * dx) /
                 ((ETA - 1) * Complex(G2, D2) * pow3(x));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  }
  return {};
}

Complex HartmannTran::dFdD0(Numeric dD0) const noexcept {
  Numeric dmF0 = (1 - ETA) * dD0;
  Numeric dGD = (dmF0 / (invGD * mF0));
  Numeric dinvGD = -dGD * pow2(invGD);
  Complex dsqrty = dinvGD / (2 * (ETA - 1) * Complex(G2, D2) * pow2(invGD));
  Complex ddeltax = Complex(0, 1 - ETA) * dD0;
  Complex dx = -ddeltax / ((ETA - 1) * Complex(G2, D2));
  Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;

  switch (calcs) {
  case CalcType::Full: {
    Complex dz1 = dsqrtxy - dsqrty;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB =
        sqrt_pi *
        ((-(pow2(z1) - 1) * w1 + (pow2(z2) - 1) * w2) * dsqrty +
         ((pow2(z1) - 1) * 1i * dw1 * dz1 - (pow2(z2) - 1) * 1i * dw2 * dz2 +
          2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) *
             sqrty) /
        (2 * (ETA - 1) * Complex(G2, D2) * pow2(sqrty));
    Complex dK = ETA * Complex(G2, D2) * dB + Complex(0, ETA * dD0) * A +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tLowZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        -(sqrt_pi * ((pow2(z1) - 1) * 1i * dw1 * dz1 + 2 * w1 * z1 * dz1) -
          dz1) *
            invGD -
        (sqrt_pi * (pow2(z1) - 1) * w1 - z1) * dinvGD;
    Complex dK = ETA * Complex(G2, D2) * dB + Complex(0, ETA * dD0) * A +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tHighZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB + Complex(0, ETA * dD0) * A +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowXandHighY: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB + Complex(0, ETA * dD0) * A +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandLowX: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dx / (2 * sqrtx);
    Complex dA = 2 * sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) /
                 ((ETA - 1) * Complex(G2, D2));
    Complex dB =
        -(2 * sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) *
              (2 * pow2(sqrty) + x - 1) +
          2 * sqrt_pi * w1 * dz1 + Complex(0, 2 * sqrt_pi) * z1 * dw1 * dz1 +
          2 * (4 * sqrty * dsqrty + dx) * (sqrt_pi * w2 * z2 - 1)) /
        ((ETA - 1) * Complex(G2, D2));
    Complex dK = ETA * Complex(G2, D2) * dB + Complex(0, ETA * dD0) * A +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandHighX: {
    Complex dz1 = dsqrtxy;
    Complex dA = (x - 3) * dx / ((ETA - 1) * Complex(G2, D2) * pow3(x));
    Complex dB = (-2 * sqrt_pi * (w1 * dz1 + z1 * 1i * dw1 * dz1) * pow3(x) +
                  (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x / 2 -
                  (x - 3) * (2 * pow2(sqrty) + x - 1) * dx) /
                 ((ETA - 1) * Complex(G2, D2) * pow3(x));
    Complex dK = ETA * Complex(G2, D2) * dB + Complex(0, ETA * dD0) * A +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  }
  return {};
}

Complex HartmannTran::dFdG0(Numeric dG0) const noexcept {
  Numeric ddeltax = (1 - ETA) * dG0;
  Complex dx = -ddeltax / ((ETA - 1) * Complex(G2, D2));
  Complex dsqrtxy = dx / (2 * sqrtxy);

  switch (calcs) {
  case CalcType::Full: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dsqrtxy;
    Complex dA = Complex(0, sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
    Complex dB =
        sqrt_pi *
        ((pow2(z1) - 1) * 1i * dw1 * dz1 - (pow2(z2) - 1) * 1i * dw2 * dz2 +
         2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) /
        (2 * sqrty * (ETA - 1) * Complex(G2, D2));
    Complex dK = ETA * Complex(G2, D2) * dB + ETA * A * dG0 +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tLowZ: {
    Complex dz1 = invGD * ddeltax;
    Complex dA = Complex(0, sqrt_pi * invGD) * dw1 * dz1;
    Complex dB =
        -invGD *
        (sqrt_pi * ((pow2(z1) - 1) * 1i * dw1 * dz1 + 2 * w1 * z1 * dz1) - dz1);
    Complex dK = ETA * Complex(G2, D2) * dB + ETA * A * dG0 +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tHighZ: {
    Complex dz1 = invGD * ddeltax;
    Complex dA = Complex(0, sqrt_pi * invGD) * dw1 * dz1;
    Complex dB = Complex(0, sqrt_pi * invGD) * dw1 * dz1 -
                 invGD * dz1 / (2 * pow2(z1)) +
                 9 * invGD * dz1 / (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB + ETA * A * dG0 +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowXandHighY: {
    Complex dz1 = invGD * ddeltax;
    Complex dz2 = dsqrtxy;
    Complex dA = Complex(0, sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
    Complex dB = Complex(0, sqrt_pi * invGD) * dw1 * dz1 -
                 invGD * dz1 / (2 * pow2(z1)) +
                 9 * invGD * dz1 / (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB + ETA * A * dG0 +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandLowX: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dx / (2 * sqrtx);
    Complex dA = 2 * sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) /
                 ((ETA - 1) * Complex(G2, D2));
    Complex dB =
        -(2 * sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) *
              (2 * pow2(sqrty) + x - 1) +
          2 * sqrt_pi * w1 * dz1 + Complex(0, 2 * sqrt_pi) * z1 * dw1 * dz1 +
          2 * (sqrt_pi * w2 * z2 - 1) * dx) /
        ((ETA - 1) * Complex(G2, D2));
    Complex dK = ETA * Complex(G2, D2) * dB + ETA * A * dG0 +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandHighX: {
    Complex dz1 = dsqrtxy;
    Complex dA = (x - 3) * dx / ((ETA - 1) * Complex(G2, D2) * pow3(x));
    Complex dB =
        (-2 * sqrt_pi * (w1 * dz1 + z1 * 1i * dw1 * dz1) * pow3(x) -
         (x - 3) * (2 * pow2(sqrty) + x - 1) * dx + (2 * x - 3) * x * dx / 2) /
        ((ETA - 1) * Complex(G2, D2) * pow3(x));
    Complex dK = ETA * Complex(G2, D2) * dB + ETA * A * dG0 +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  }
  return {};
}

Complex HartmannTran::dFdD2(Numeric dD2) const noexcept {
  Numeric dmF0 = -3 * (1 - ETA) * dD2 / 2;
  Numeric dGD = (dmF0 / (invGD * mF0));
  Numeric dinvGD = -dGD * pow2(invGD);
  Complex dsqrty = (Complex(G2, D2) * dinvGD + Complex(0, invGD) * dD2) /
                   (2 * (ETA - 1) * pow2(Complex(G2, D2)) * pow2(invGD));
  Complex ddeltax = 1.5 * Complex(0, ETA - 1) * dD2;
  Complex dx = (-Complex(G2, D2) * ddeltax + Complex(0, dD2) * deltax) /
               ((ETA - 1) * pow2(Complex(G2, D2)));
  Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;

  switch (calcs) {
  case CalcType::Full: {
    Complex dz1 = dsqrtxy - dsqrty;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB = (sqrt_pi * Complex(G2, D2) *
                      ((-(pow2(z1) - 1) * w1 + (pow2(z2) - 1) * w2) * dsqrty +
                       ((pow2(z1) - 1) * 1i * dw1 * dz1 -
                        (pow2(z2) - 1) * 1i * dw2 * dz2 + 2 * w1 * z1 * dz1 -
                        2 * w2 * z2 * dz2) *
                           sqrty) -
                  1i *
                      (sqrt_pi * (pow2(z1) - 1) * w1 -
                       sqrt_pi * (pow2(z2) - 1) * w2 + 2 * sqrty) *
                      sqrty * dD2) /
                 (2 * (ETA - 1) * pow2(Complex(G2, D2)) * pow2(sqrty));
    Complex dK = ETA * Complex(G2, D2) * dB - Complex(0, 1.5 * dD2 * ETA) * A +
                 Complex(0, dD2 * ETA) * B +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tLowZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        -(sqrt_pi * ((pow2(z1) - 1) * 1i * dw1 * dz1 + 2 * w1 * z1 * dz1) -
          dz1) *
            invGD -
        (sqrt_pi * (pow2(z1) - 1) * w1 - z1) * dinvGD;
    Complex dK = ETA * Complex(G2, D2) * dB - Complex(0, 1.5 * dD2 * ETA) * A +
                 Complex(0, dD2 * ETA) * B +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tHighZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB - Complex(0, 1.5 * dD2 * ETA) * A +
                 Complex(0, dD2 * ETA) * B +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowXandHighY: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB - Complex(0, 1.5 * dD2 * ETA) * A +
                 Complex(0, dD2 * ETA) * B +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandLowX: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dx / (2 * sqrtx);
    Complex dA = 2 *
                 (sqrt_pi * Complex(G2, D2) * (w2 * dz2 + z2 * 1i * dw2 * dz2) -
                  1i * (sqrt_pi * w2 * z2 - 1) * dD2) /
                 ((ETA - 1) * pow2(Complex(G2, D2)));
    Complex dB =
        (-2 * Complex(G2, D2) *
             (sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) *
                  (2 * pow2(sqrty) + x - 1) +
              sqrt_pi * w1 * dz1 + Complex(0, sqrt_pi) * z1 * dw1 * dz1 +
              (4 * sqrty * dsqrty + dx) * (sqrt_pi * w2 * z2 - 1)) +
         1i *
             (2 * sqrt_pi * w1 * z1 +
              2 * (sqrt_pi * w2 * z2 - 1) * (2 * pow2(sqrty) + x - 1) - 1) *
             dD2) /
        ((ETA - 1) * pow2(Complex(G2, D2)));
    Complex dK = ETA * Complex(G2, D2) * dB - Complex(0, 1.5 * dD2 * ETA) * A +
                 Complex(0, dD2 * ETA) * B +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandHighX: {
    Complex dz1 = dsqrtxy;
    Complex dA =
        (Complex(G2, D2) * (x - 3) * dx + 1i * (2 * x - 3) * x * dD2 / 2) /
        ((ETA - 1) * pow2(Complex(G2, D2)) * pow3(x));
    Complex dB =
        (Complex(G2, D2) *
             (-4 * sqrt_pi * (w1 * dz1 + z1 * 1i * dw1 * dz1) * pow3(x) +
              (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x -
              2 * (x - 3) * (2 * pow2(sqrty) + x - 1) * dx) -
         1i *
             (2 * (-2 * sqrt_pi * w1 * z1 + 1) * pow2(x) +
              (2 * x - 3) * (2 * pow2(sqrty) + x - 1)) *
             x * dD2) /
        (2 * (ETA - 1) * pow2(Complex(G2, D2)) * pow3(x));
    Complex dK = ETA * Complex(G2, D2) * dB - Complex(0, 1.5 * dD2 * ETA) * A +
                 Complex(0, dD2 * ETA) * B +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  }
  return {};
}

Complex HartmannTran::dFdG2(Numeric dG2) const noexcept {
  Complex dsqrty = dG2 / (2 * invGD * (ETA - 1) * pow2(Complex(G2, D2)));
  Numeric ddeltax = 3 * (ETA - 1) * dG2 / 2;
  Complex dx = (-Complex(G2, D2) * ddeltax + deltax * dG2) /
               ((ETA - 1) * pow2(Complex(G2, D2)));
  Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;

  switch (calcs) {
  case CalcType::Full: {
    Complex dz1 = dsqrtxy - dsqrty;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = Complex(0, sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
    Complex dB = (sqrt_pi * Complex(G2, D2) *
                      ((-(pow2(z1) - 1) * w1 + (pow2(z2) - 1) * w2) * dsqrty +
                       ((pow2(z1) - 1) * 1i * dw1 * dz1 -
                        (pow2(z2) - 1) * 1i * dw2 * dz2 + 2 * w1 * z1 * dz1 -
                        2 * w2 * z2 * dz2) *
                           sqrty) -
                  (sqrt_pi * (pow2(z1) - 1) * w1 -
                   sqrt_pi * (pow2(z2) - 1) * w2 + 2 * sqrty) *
                      sqrty * dG2) /
                 (2 * (ETA - 1) * pow2(Complex(G2, D2)) * pow2(sqrty));
    Complex dK = ETA * Complex(G2, D2) * dB - 3 * ETA * A * dG2 / 2 +
                 ETA * B * dG2 +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tLowZ: {
    Complex dz1 = invGD * ddeltax;
    Complex dA = Complex(0, sqrt_pi * invGD) * dw1 * dz1;
    Complex dB =
        -invGD *
        (sqrt_pi * ((pow2(z1) - 1) * 1i * dw1 * dz1 + 2 * w1 * z1 * dz1) - dz1);
    Complex dK = ETA * Complex(G2, D2) * dB - 3 * ETA * A * dG2 / 2 +
                 ETA * B * dG2 +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tHighZ: {
    Complex dz1 = invGD * ddeltax;
    Complex dA = Complex(0, sqrt_pi * invGD) * dw1 * dz1;
    Complex dB = Complex(0, sqrt_pi * invGD) * dw1 * dz1 -
                 invGD * dz1 / (2 * pow2(z1)) +
                 9 * invGD * dz1 / (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB - 3 * ETA * A * dG2 / 2 +
                 ETA * B * dG2 +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowXandHighY: {
    Complex dz1 = invGD * ddeltax;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = Complex(0, sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
    Complex dB = Complex(0, sqrt_pi * invGD) * dw1 * dz1 -
                 invGD * dz1 / (2 * pow2(z1)) +
                 9 * invGD * dz1 / (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB - 3 * ETA * A * dG2 / 2 +
                 ETA * B * dG2 +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandLowX: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dx / (2 * sqrtx);
    Complex dA = 2 *
                 (sqrt_pi * Complex(G2, D2) * (w2 * dz2 + z2 * 1i * dw2 * dz2) -
                  (sqrt_pi * w2 * z2 - 1) * dG2) /
                 ((ETA - 1) * pow2(Complex(G2, D2)));
    Complex dB =
        (-2 * Complex(G2, D2) *
             (sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) *
                  (2 * pow2(sqrty) + x - 1) +
              sqrt_pi * w1 * dz1 + Complex(0, sqrt_pi) * z1 * dw1 * dz1 +
              (4 * sqrty * dsqrty + dx) * (sqrt_pi * w2 * z2 - 1)) +
         (2 * sqrt_pi * w1 * z1 +
          2 * (sqrt_pi * w2 * z2 - 1) * (2 * pow2(sqrty) + x - 1) - 1) *
             dG2) /
        ((ETA - 1) * pow2(Complex(G2, D2)));
    Complex dK = ETA * Complex(G2, D2) * dB - 3 * ETA * A * dG2 / 2 +
                 ETA * B * dG2 +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandHighX: {
    Complex dz1 = dsqrtxy;
    Complex dA = (Complex(G2, D2) * (x - 3) * dx + (2 * x - 3) * x * dG2 / 2) /
                 ((ETA - 1) * pow2(Complex(G2, D2)) * pow3(x));
    Complex dB =
        (Complex(G2, D2) *
             (-4 * sqrt_pi * (w1 * dz1 + z1 * 1i * dw1 * dz1) * pow3(x) +
              (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x -
              2 * (x - 3) * (2 * pow2(sqrty) + x - 1) * dx) -
         (2 * (-2 * sqrt_pi * w1 * z1 + 1) * pow2(x) +
          (2 * x - 3) * (2 * pow2(sqrty) + x - 1)) *
             x * dG2) /
        (2 * (ETA - 1) * pow2(Complex(G2, D2)) * pow3(x));
    Complex dK = ETA * Complex(G2, D2) * dB - 3 * ETA * A * dG2 / 2 +
                 ETA * B * dG2 +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  }
  return {};
}

Complex HartmannTran::dFdFVC(Numeric dFVC) const noexcept {
  Numeric ddeltax = dFVC;
  Complex dx = -ddeltax / ((ETA - 1) * Complex(G2, D2));
  Complex dsqrtxy = dx / (2 * sqrtxy);

  switch (calcs) {
  case CalcType::Full: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dsqrtxy;
    Complex dA = Complex(0, sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
    Complex dB =
        sqrt_pi *
        ((pow2(z1) - 1) * 1i * dw1 * dz1 - (pow2(z2) - 1) * 1i * dw2 * dz2 +
         2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) /
        (2 * sqrty * (ETA - 1) * Complex(G2, D2));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA -
                 A * dFVC;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tLowZ: {
    Complex dz1 = invGD * ddeltax;
    Complex dA = Complex(0, sqrt_pi * invGD) * dw1 * dz1;
    Complex dB =
        -invGD *
        (sqrt_pi * ((pow2(z1) - 1) * 1i * dw1 * dz1 + 2 * w1 * z1 * dz1) - dz1);
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA -
                 A * dFVC;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tHighZ: {
    Complex dz1 = invGD * ddeltax;
    Complex dA = Complex(0, sqrt_pi * invGD) * dw1 * dz1;
    Complex dB = Complex(0, sqrt_pi * invGD) * dw1 * dz1 -
                 invGD * dz1 / (2 * pow2(z1)) +
                 9 * invGD * dz1 / (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA -
                 A * dFVC;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowXandHighY: {
    Complex dz1 = invGD * ddeltax;
    Complex dz2 = dsqrtxy;
    Complex dA = Complex(0, sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
    Complex dB = Complex(0, sqrt_pi * invGD) * dw1 * dz1 -
                 invGD * dz1 / (2 * pow2(z1)) +
                 9 * invGD * dz1 / (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA -
                 A * dFVC;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandLowX: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dx / (2 * sqrtx);
    Complex dA = 2 * sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) /
                 ((ETA - 1) * Complex(G2, D2));
    Complex dB =
        -(2 * sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) *
              (2 * pow2(sqrty) + x - 1) +
          2 * sqrt_pi * w1 * dz1 + Complex(0, 2 * sqrt_pi) * z1 * dw1 * dz1 +
          2 * (sqrt_pi * w2 * z2 - 1) * dx) /
        ((ETA - 1) * Complex(G2, D2));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA -
                 A * dFVC;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandHighX: {
    Complex dz1 = dsqrtxy;
    Complex dA = (x - 3) * dx / ((ETA - 1) * Complex(G2, D2) * pow3(x));
    Complex dB =
        (-2 * sqrt_pi * (w1 * dz1 + z1 * 1i * dw1 * dz1) * pow3(x) -
         (x - 3) * (2 * pow2(sqrty) + x - 1) * dx + (2 * x - 3) * x * dx / 2) /
        ((ETA - 1) * Complex(G2, D2) * pow3(x));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA -
                 A * dFVC;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  }
  return {};
}

Complex HartmannTran::dFdETA(Numeric dETA) const noexcept {
  Numeric dmF0 = -(D0 - 3 * D2 / 2) * dETA;
  Numeric dGD = (dmF0 / (invGD * mF0));
  Numeric dinvGD = -dGD * pow2(invGD);
  Complex dsqrty = ((ETA - 1) * dinvGD + invGD * dETA) /
                   (2 * Complex(G2, D2) * pow2(ETA - 1) * pow2(invGD));
  Complex ddeltax = -dETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2);
  Complex dx = (-(ETA - 1) * ddeltax + deltax * dETA) /
               (Complex(G2, D2) * pow2(ETA - 1));
  Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;

  switch (calcs) {
  case CalcType::Full: {
    Complex dz1 = dsqrtxy - dsqrty;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB = (sqrt_pi *
                      ((-(pow2(z1) - 1) * w1 + (pow2(z2) - 1) * w2) * dsqrty +
                       ((pow2(z1) - 1) * 1i * dw1 * dz1 -
                        (pow2(z2) - 1) * 1i * dw2 * dz2 + 2 * w1 * z1 * dz1 -
                        2 * w2 * z2 * dz2) *
                           sqrty) *
                      (ETA - 1) -
                  (sqrt_pi * (pow2(z1) - 1) * w1 -
                   sqrt_pi * (pow2(z2) - 1) * w2 + 2 * sqrty) *
                      sqrty * dETA) /
                 (2 * Complex(G2, D2) * pow2(ETA - 1) * pow2(sqrty));
    Complex dK = (-FVC + Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA) * dA +
                 Complex(G2, D2) * B * dETA + Complex(G2, D2) * ETA * dB -
                 Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * A * dETA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tLowZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        -(sqrt_pi * ((pow2(z1) - 1) * 1i * dw1 * dz1 + 2 * w1 * z1 * dz1) -
          dz1) *
            invGD -
        (sqrt_pi * (pow2(z1) - 1) * w1 - z1) * dinvGD;
    Complex dK = (-FVC + Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA) * dA +
                 Complex(G2, D2) * B * dETA + Complex(G2, D2) * ETA * dB -
                 Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * A * dETA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tHighZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = (-FVC + Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA) * dA +
                 Complex(G2, D2) * B * dETA + Complex(G2, D2) * ETA * dB -
                 Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * A * dETA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowXandHighY: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = (-FVC + Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA) * dA +
                 Complex(G2, D2) * B * dETA + Complex(G2, D2) * ETA * dB -
                 Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * A * dETA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandLowX: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dx / (2 * sqrtx);
    Complex dA = 2 *
                 (sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) * (ETA - 1) -
                  (sqrt_pi * w2 * z2 - 1) * dETA) /
                 (Complex(G2, D2) * pow2(ETA - 1));
    Complex dB =
        (-2 * (ETA - 1) *
             (sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) *
                  (2 * pow2(sqrty) + x - 1) +
              sqrt_pi * w1 * dz1 + Complex(0, sqrt_pi) * z1 * dw1 * dz1 +
              (4 * sqrty * dsqrty + dx) * (sqrt_pi * w2 * z2 - 1)) +
         (2 * sqrt_pi * w1 * z1 +
          2 * (sqrt_pi * w2 * z2 - 1) * (2 * pow2(sqrty) + x - 1) - 1) *
             dETA) /
        (Complex(G2, D2) * pow2(ETA - 1));
    Complex dK = (-FVC + Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA) * dA +
                 Complex(G2, D2) * B * dETA + Complex(G2, D2) * ETA * dB -
                 Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * A * dETA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandHighX: {
    Complex dz1 = dsqrtxy;
    Complex dA = ((ETA - 1) * (x - 3) * dx + (2 * x - 3) * x * dETA / 2) /
                 (Complex(G2, D2) * pow2(ETA - 1) * pow3(x));
    Complex dB = (-(2 * (-2 * sqrt_pi * w1 * z1 + 1) * pow2(x) +
                    (2 * x - 3) * (2 * pow2(sqrty) + x - 1)) *
                      x * dETA +
                  (ETA - 1) * (-4 * sqrt_pi * (w1 * dz1 + z1 * 1i * dw1 * dz1) *
                                   pow3(x) +
                               (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x -
                               2 * (x - 3) * (2 * pow2(sqrty) + x - 1) * dx)) /
                 (2 * Complex(G2, D2) * pow2(ETA - 1) * pow3(x));
    Complex dK = (-FVC + Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA) * dA +
                 Complex(G2, D2) * B * dETA + Complex(G2, D2) * ETA * dB -
                 Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * A * dETA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  }
  return {};
}

Complex HartmannTran::dFdH(Numeric dZ) const noexcept {
  Numeric dmF0 = dZ;
  Numeric dGD = (dmF0 / (invGD * mF0));
  Numeric dinvGD = -dGD * pow2(invGD);
  Complex dsqrty = dinvGD / (2 * (ETA - 1) * Complex(G2, D2) * pow2(invGD));
  Complex ddeltax = Complex(0, dZ);
  Complex dx = -ddeltax / ((ETA - 1) * Complex(G2, D2));
  Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;

  switch (calcs) {
  case CalcType::Full: {
    Complex dz1 = dsqrtxy - dsqrty;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB =
        sqrt_pi *
        ((-(pow2(z1) - 1) * w1 + (pow2(z2) - 1) * w2) * dsqrty +
         ((pow2(z1) - 1) * 1i * dw1 * dz1 - (pow2(z2) - 1) * 1i * dw2 * dz2 +
          2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) *
             sqrty) /
        (2 * (ETA - 1) * Complex(G2, D2) * pow2(sqrty));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tLowZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        -(sqrt_pi * ((pow2(z1) - 1) * 1i * dw1 * dz1 + 2 * w1 * z1 * dz1) -
          dz1) *
            invGD -
        (sqrt_pi * (pow2(z1) - 1) * w1 - z1) * dinvGD;
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tHighZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowXandHighY: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandLowX: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dx / (2 * sqrtx);
    Complex dA = 2 * sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) /
                 ((ETA - 1) * Complex(G2, D2));
    Complex dB =
        -(2 * sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) *
              (2 * pow2(sqrty) + x - 1) +
          2 * sqrt_pi * w1 * dz1 + Complex(0, 2 * sqrt_pi) * z1 * dw1 * dz1 +
          2 * (4 * sqrty * dsqrty + dx) * (sqrt_pi * w2 * z2 - 1)) /
        ((ETA - 1) * Complex(G2, D2));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandHighX: {
    Complex dz1 = dsqrtxy;
    Complex dA = (x - 3) * dx / ((ETA - 1) * Complex(G2, D2) * pow3(x));
    Complex dB = (-2 * sqrt_pi * (w1 * dz1 + z1 * 1i * dw1 * dz1) * pow3(x) +
                  (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x / 2 -
                  (x - 3) * (2 * pow2(sqrty) + x - 1) * dx) /
                 ((ETA - 1) * Complex(G2, D2) * pow3(x));
    Complex dK = ETA * Complex(G2, D2) * dB +
                 (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  }
  return {};
}

Complex HartmannTran::dFdVMR(const Output &d) const noexcept {
  Numeric dmF0 = (1 - ETA) * (d.D0 - 3 * d.D2 / 2) - (D0 - 3 * D2 / 2) * d.ETA;
  Numeric dGD = (dmF0 / (invGD * mF0));
  Numeric dinvGD = -dGD * pow2(invGD);
  Complex dsqrty =
      (Complex(G2, D2) * (ETA - 1) * dinvGD + Complex(G2, D2) * invGD * d.ETA +
       Complex(d.G2, d.D2) * (ETA - 1) * invGD) /
      (2 * pow2(Complex(G2, D2)) * pow2(ETA - 1) * pow2(invGD));
  Complex ddeltax = -(ETA - 1) * Complex(d.G0 - 1.5 * d.G2, d.D0 - 1.5 * d.D2) -
                    Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * d.ETA + d.FVC;
  Complex dx = (-Complex(G2, D2) * (ETA - 1) * ddeltax +
                Complex(G2, D2) * deltax * d.ETA +
                Complex(d.G2, d.D2) * (ETA - 1) * deltax) /
               (pow2(Complex(G2, D2)) * pow2(ETA - 1));
  Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;

  switch (calcs) {
  case CalcType::Full: {
    Complex dz1 = dsqrtxy - dsqrty;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB = (sqrt_pi * Complex(G2, D2) *
                      ((-(pow2(z1) - 1) * w1 + (pow2(z2) - 1) * w2) * dsqrty +
                       ((pow2(z1) - 1) * 1i * dw1 * dz1 -
                        (pow2(z2) - 1) * 1i * dw2 * dz2 + 2 * w1 * z1 * dz1 -
                        2 * w2 * z2 * dz2) *
                           sqrty) *
                      (ETA - 1) -
                  Complex(G2, D2) *
                      (sqrt_pi * (pow2(z1) - 1) * w1 -
                       sqrt_pi * (pow2(z2) - 1) * w2 + 2 * sqrty) *
                      sqrty * d.ETA -
                  Complex(d.G2, d.D2) * (ETA - 1) *
                      (sqrt_pi * (pow2(z1) - 1) * w1 -
                       sqrt_pi * (pow2(z2) - 1) * w2 + 2 * sqrty) *
                      sqrty) /
                 (2 * pow2(Complex(G2, D2)) * pow2(ETA - 1) * pow2(sqrty));
    Complex dK = Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                 Complex(d.G2, d.D2) * B * ETA +
                 (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                 (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                  Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                     A;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tLowZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        -(sqrt_pi * ((pow2(z1) - 1) * 1i * dw1 * dz1 + 2 * w1 * z1 * dz1) -
          dz1) *
            invGD -
        (sqrt_pi * (pow2(z1) - 1) * w1 - z1) * dinvGD;
    Complex dK = Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                 Complex(d.G2, d.D2) * B * ETA +
                 (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                 (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                  Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                     A;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tHighZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                 Complex(d.G2, d.D2) * B * ETA +
                 (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                 (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                  Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                     A;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowXandHighY: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                 Complex(d.G2, d.D2) * B * ETA +
                 (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                 (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                  Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                     A;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandLowX: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dx / (2 * sqrtx);
    Complex dA = 2 *
                 (sqrt_pi * Complex(G2, D2) * (w2 * dz2 + z2 * 1i * dw2 * dz2) *
                      (ETA - 1) -
                  Complex(G2, D2) * (sqrt_pi * w2 * z2 - 1) * d.ETA -
                  Complex(d.G2, d.D2) * (sqrt_pi * w2 * z2 - 1) * (ETA - 1)) /
                 (pow2(Complex(G2, D2)) * pow2(ETA - 1));
    Complex dB =
        (-2 * Complex(G2, D2) * (ETA - 1) *
             (sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) *
                  (2 * pow2(sqrty) + x - 1) +
              sqrt_pi * w1 * dz1 + Complex(0, sqrt_pi) * z1 * dw1 * dz1 +
              (4 * sqrty * dsqrty + dx) * (sqrt_pi * w2 * z2 - 1)) +
         Complex(G2, D2) *
             (2 * sqrt_pi * w1 * z1 +
              2 * (sqrt_pi * w2 * z2 - 1) * (2 * pow2(sqrty) + x - 1) - 1) *
             d.ETA +
         Complex(d.G2, d.D2) * (ETA - 1) *
             (2 * sqrt_pi * w1 * z1 +
              2 * (sqrt_pi * w2 * z2 - 1) * (2 * pow2(sqrty) + x - 1) - 1)) /
        (pow2(Complex(G2, D2)) * pow2(ETA - 1));
    Complex dK = Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                 Complex(d.G2, d.D2) * B * ETA +
                 (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                 (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                  Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                     A;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandHighX: {
    Complex dz1 = dsqrtxy;
    Complex dA = (2 * Complex(G2, D2) * (ETA - 1) * (x - 3) * dx +
                  Complex(G2, D2) * (2 * x - 3) * x * d.ETA +
                  Complex(d.G2, d.D2) * (ETA - 1) * (2 * x - 3) * x) /
                 (2 * pow2(Complex(G2, D2)) * pow2(ETA - 1) * pow3(x));
    Complex dB =
        (-Complex(G2, D2) *
             (2 * (-2 * sqrt_pi * w1 * z1 + 1) * pow2(x) +
              (2 * x - 3) * (2 * pow2(sqrty) + x - 1)) *
             x * d.ETA +
         Complex(G2, D2) * (ETA - 1) *
             (-4 * sqrt_pi * (w1 * dz1 + z1 * 1i * dw1 * dz1) * pow3(x) +
              (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x -
              2 * (x - 3) * (2 * pow2(sqrty) + x - 1) * dx) -
         Complex(d.G2, d.D2) *
             (2 * (-2 * sqrt_pi * w1 * z1 + 1) * pow2(x) +
              (2 * x - 3) * (2 * pow2(sqrty) + x - 1)) *
             (ETA - 1) * x) /
        (2 * pow2(Complex(G2, D2)) * pow2(ETA - 1) * pow3(x));
    Complex dK = Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                 Complex(d.G2, d.D2) * B * ETA +
                 (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                 (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                  Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                     A;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  }
  return {};
}

Complex HartmannTran::dFdT(const Output &d, Numeric T) const noexcept {
  Numeric dmF0 = (1 - ETA) * (d.D0 - 3 * d.D2 / 2) - (D0 - 3 * D2 / 2) * d.ETA;
  Numeric dGD = (dmF0 / (invGD * mF0)) - invGD * invGD / (2 * T * sqrt_ln_2);
  Numeric dinvGD = -dGD * pow2(invGD);
  Complex dsqrty =
      (Complex(G2, D2) * (ETA - 1) * dinvGD + Complex(G2, D2) * invGD * d.ETA +
       Complex(d.G2, d.D2) * (ETA - 1) * invGD) /
      (2 * pow2(Complex(G2, D2)) * pow2(ETA - 1) * pow2(invGD));
  Complex ddeltax = -(ETA - 1) * Complex(d.G0 - 1.5 * d.G2, d.D0 - 1.5 * d.D2) -
                    Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * d.ETA + d.FVC;
  Complex dx = (-Complex(G2, D2) * (ETA - 1) * ddeltax +
                Complex(G2, D2) * deltax * d.ETA +
                Complex(d.G2, d.D2) * (ETA - 1) * deltax) /
               (pow2(Complex(G2, D2)) * pow2(ETA - 1));
  Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;

  switch (calcs) {
  case CalcType::Full: {
    Complex dz1 = dsqrtxy - dsqrty;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB = (sqrt_pi * Complex(G2, D2) *
                      ((-(pow2(z1) - 1) * w1 + (pow2(z2) - 1) * w2) * dsqrty +
                       ((pow2(z1) - 1) * 1i * dw1 * dz1 -
                        (pow2(z2) - 1) * 1i * dw2 * dz2 + 2 * w1 * z1 * dz1 -
                        2 * w2 * z2 * dz2) *
                           sqrty) *
                      (ETA - 1) -
                  Complex(G2, D2) *
                      (sqrt_pi * (pow2(z1) - 1) * w1 -
                       sqrt_pi * (pow2(z2) - 1) * w2 + 2 * sqrty) *
                      sqrty * d.ETA -
                  Complex(d.G2, d.D2) * (ETA - 1) *
                      (sqrt_pi * (pow2(z1) - 1) * w1 -
                       sqrt_pi * (pow2(z2) - 1) * w2 + 2 * sqrty) *
                      sqrty) /
                 (2 * pow2(Complex(G2, D2)) * pow2(ETA - 1) * pow2(sqrty));
    Complex dK = Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                 Complex(d.G2, d.D2) * B * ETA +
                 (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                 (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                  Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                     A;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tLowZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        -(sqrt_pi * ((pow2(z1) - 1) * 1i * dw1 * dz1 + 2 * w1 * z1 * dz1) -
          dz1) *
            invGD -
        (sqrt_pi * (pow2(z1) - 1) * w1 - z1) * dinvGD;
    Complex dK = Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                 Complex(d.G2, d.D2) * B * ETA +
                 (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                 (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                  Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                     A;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::Noc2tHighZ: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dA = sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                 Complex(d.G2, d.D2) * B * ETA +
                 (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                 (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                  Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                     A;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowXandHighY: {
    Complex dz1 = deltax * dinvGD + invGD * ddeltax;
    Complex dz2 = dsqrtxy + dsqrty;
    Complex dA = sqrt_pi * ((w1 - w2) * dinvGD +
                            (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
    Complex dB =
        ((4 * sqrt_pi * w1 * pow3(z1) + 2 * pow2(z1) - 3) * z1 * dinvGD +
         (Complex(0, 4 * sqrt_pi) * pow4(z1) * dw1 * dz1 - 2 * pow2(z1) * dz1 +
          9 * dz1) *
             invGD) /
        (4 * pow4(z1));
    Complex dK = Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                 Complex(d.G2, d.D2) * B * ETA +
                 (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                 (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                  Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                     A;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandLowX: {
    Complex dz1 = dsqrtxy;
    Complex dz2 = dx / (2 * sqrtx);
    Complex dA = 2 *
                 (sqrt_pi * Complex(G2, D2) * (w2 * dz2 + z2 * 1i * dw2 * dz2) *
                      (ETA - 1) -
                  Complex(G2, D2) * (sqrt_pi * w2 * z2 - 1) * d.ETA -
                  Complex(d.G2, d.D2) * (sqrt_pi * w2 * z2 - 1) * (ETA - 1)) /
                 (pow2(Complex(G2, D2)) * pow2(ETA - 1));
    Complex dB =
        (-2 * Complex(G2, D2) * (ETA - 1) *
             (sqrt_pi * (w2 * dz2 + z2 * 1i * dw2 * dz2) *
                  (2 * pow2(sqrty) + x - 1) +
              sqrt_pi * w1 * dz1 + Complex(0, sqrt_pi) * z1 * dw1 * dz1 +
              (4 * sqrty * dsqrty + dx) * (sqrt_pi * w2 * z2 - 1)) +
         Complex(G2, D2) *
             (2 * sqrt_pi * w1 * z1 +
              2 * (sqrt_pi * w2 * z2 - 1) * (2 * pow2(sqrty) + x - 1) - 1) *
             d.ETA +
         Complex(d.G2, d.D2) * (ETA - 1) *
             (2 * sqrt_pi * w1 * z1 +
              2 * (sqrt_pi * w2 * z2 - 1) * (2 * pow2(sqrty) + x - 1) - 1)) /
        (pow2(Complex(G2, D2)) * pow2(ETA - 1));
    Complex dK = Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                 Complex(d.G2, d.D2) * B * ETA +
                 (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                 (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                  Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                     A;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  case CalcType::LowYandHighX: {
    Complex dz1 = dsqrtxy;
    Complex dA = (2 * Complex(G2, D2) * (ETA - 1) * (x - 3) * dx +
                  Complex(G2, D2) * (2 * x - 3) * x * d.ETA +
                  Complex(d.G2, d.D2) * (ETA - 1) * (2 * x - 3) * x) /
                 (2 * pow2(Complex(G2, D2)) * pow2(ETA - 1) * pow3(x));
    Complex dB =
        (-Complex(G2, D2) *
             (2 * (-2 * sqrt_pi * w1 * z1 + 1) * pow2(x) +
              (2 * x - 3) * (2 * pow2(sqrty) + x - 1)) *
             x * d.ETA +
         Complex(G2, D2) * (ETA - 1) *
             (-4 * sqrt_pi * (w1 * dz1 + z1 * 1i * dw1 * dz1) * pow3(x) +
              (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x -
              2 * (x - 3) * (2 * pow2(sqrty) + x - 1) * dx) -
         Complex(d.G2, d.D2) *
             (2 * (-2 * sqrt_pi * w1 * z1 + 1) * pow2(x) +
              (2 * x - 3) * (2 * pow2(sqrty) + x - 1)) *
             (ETA - 1) * x) /
        (2 * pow2(Complex(G2, D2)) * pow2(ETA - 1) * pow3(x));
    Complex dK = Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                 Complex(d.G2, d.D2) * B * ETA +
                 (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                 (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                  Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                     A;
    return inv_pi * (-A * dK + K * dA) / pow2(K);
  }
  }
  return {};
}

Complex HartmannTran::operator()(Numeric f) noexcept {
  imag_val(deltax) = mF0 - f;
  x = deltax / ((1 - ETA) * Complex(G2, D2));
  sqrtxy = std::sqrt(x + sqrty * sqrty);
  update_calcs();
  calc();
  return F;
}

HartmannTran::CalcType HartmannTran::init(const Complex c2t) const noexcept {
  if (abs_squared(c2t) == 0)
    return CalcType::Noc2tHighZ; // nb. Value of high/low changes elsewhere
  if (abs_squared(x) <= 9e-16 * abs_squared(sqrty * sqrty))
    return CalcType::LowXandHighY;
  if ((abs_squared(sqrty * sqrty) <= 1.e-30 * abs_squared(x)) and
           abs_squared(std::sqrt(x)) <= 16.e6)
    return CalcType::LowYandLowX; // Weird case, untested
  if ((abs_squared(sqrty * sqrty) <= 1.e-30 * abs_squared(x)))
    return CalcType::LowYandHighX;
  return CalcType::Full;
}

void HartmannTran::update_calcs() noexcept {
  calcs = init((1 - ETA) * Complex(G2, D2));
}

void HartmannTran::calc() noexcept {
  switch (calcs) {
  case CalcType::Full:
    z1 = sqrtxy - sqrty;
    z2 = sqrtxy + sqrty;
    w1 = Faddeeva::w(1i * z1);
    w2 = Faddeeva::w(1i * z2);
    A = sqrt_pi * invGD * (w1 - w2);
    B = (-1 + sqrt_pi / (2 * sqrty) * (1 - pow2(z1)) * w1 -
         sqrt_pi / (2 * sqrty) * (1 - pow2(z2)) * w2) /
        ((1 - ETA) * Complex(G2, D2));
    break;
  case CalcType::Noc2tLowZ:
  case CalcType::Noc2tHighZ:
    z1 = deltax * invGD;
    w1 = Faddeeva::w(1i * z1);
    A = sqrt_pi * invGD * w1;
    if (abs_squared(z1) < 16e6) {
      calcs = CalcType::Noc2tLowZ;
      B = sqrt_pi * invGD * ((1 - pow2(z1)) * w1 + z1 / sqrt_pi);
    } else {
      calcs = CalcType::Noc2tHighZ;
      B = invGD * (sqrt_pi * w1 + 1 / z1 / 2 - 3 / pow3(z1) / 4);
    }
    break;
  case CalcType::LowXandHighY:
    z1 = deltax * invGD;
    z2 = sqrtxy + sqrty;
    w1 = Faddeeva::w(1i * z1);
    w2 = Faddeeva::w(1i * z2);
    A = sqrt_pi * invGD * (w1 - w2);
    B = invGD * (sqrt_pi * w1 + 1 / z1 / 2 - 3 / pow3(z1) / 4);
    break;
  case CalcType::LowYandLowX:
    sqrtx = std::sqrt(x);
    z1 = sqrtxy;
    z2 = sqrtx;
    w1 = Faddeeva::w(1i * z1);
    w2 = Faddeeva::w(1i * z2);
    A = (2 * sqrt_pi / ((1 - ETA) * Complex(G2, D2))) * (inv_sqrt_pi - z2 * w2);
    B = (1 / ((1 - ETA) * Complex(G2, D2))) *
        (-1 +
         2 * sqrt_pi * (1 - x - 2 * sqrty * sqrty) * (1 / sqrt_pi - z2 * w2) +
         2 * sqrt_pi * z1 * w1);
    break;
  case CalcType::LowYandHighX:
    z1 = sqrtxy;
    w1 = Faddeeva::w(1i * z1);
    A = (1 / ((1 - ETA) * Complex(G2, D2))) * (1 / x - 3 / pow2(x) / 2);
    B = (1 / ((1 - ETA) * Complex(G2, D2))) *
        (-1 + (1 - x - 2 * sqrty * sqrty) * (1 / x - 3 / pow2(x) / 2) +
         2 * sqrt_pi * z1 * w1);
    break;
  }

  dw1 = 2i * (inv_sqrt_pi - z1 * w1);
  dw2 = 2i * (inv_sqrt_pi - z2 * w2);
  K = 1 - (FVC - ETA * (Complex(G0, D0) - 3 * Complex(G2, D2) / 2)) * A +
      ETA * Complex(G2, D2) * B;
  F = inv_pi * A / K;
}

VanVleckHuber::VanVleckHuber(Numeric F0, Numeric T) noexcept
    : c1(Constant::h / (2.0 * Constant::k * T)), tanh_c1f0(std::tanh(c1 * F0)),
      inv_denom(1.0 / (F0 * tanh_c1f0)) {}

Numeric VanVleckHuber::dNdT(Numeric T, Numeric f) const noexcept {
  const Numeric part1 =
      (1 - pow2(tanh_c1f0)) * c1 * f * tanh_c1f / (T * pow2(tanh_c1f0));
  const Numeric part2 = (pow2(tanh_c1f) - 1) * c1 * pow2(f) * inv_denom / T;
  return part1 + part2;
}

Numeric VanVleckHuber::dNdf(Numeric f) const noexcept {
  const Numeric part1 = tanh_c1f * inv_denom;
  const Numeric part2 = c1 * f * inv_denom / (1 - pow2(tanh_c1f));
  return part1 + part2;
}

Numeric VanVleckHuber::dNdF0() const noexcept {
  const Numeric part1 = -N * tanh_c1f0 * inv_denom;
  const Numeric part2 = c1 * N * (pow2(tanh_c1f0) - 1) / tanh_c1f0;
  return part1 + part2;
}

Numeric VanVleckHuber::operator()(Numeric f) noexcept {
  tanh_c1f = std::tanh(c1 * f);
  N = f * tanh_c1f * inv_denom;
  return N;
}

RosenkranzQuadratic::RosenkranzQuadratic(Numeric F0, Numeric T) noexcept
    : fac((Constant::h) / (2.0 * Constant::k * T) /
          std::sinh((Constant::h * F0) / (2.0 * Constant::k * T)) * (1.0 / F0)),
      dfacdT((fac / T) *
             (fac * F0 * F0 *
                  std::cosh((Constant::h * F0) / (2.0 * Constant::k * T)) -
              1)),
      dfacdF0(-fac * (1 / F0 + fac * F0 *
                                   std::cosh((Constant::h * F0) /
                                             (2.0 * Constant::k * T)))) {}

Numeric RosenkranzQuadratic::dNdT(Numeric, Numeric f) const noexcept {
  return f * f * dfacdT;
}

Numeric RosenkranzQuadratic::dNdf(Numeric f) const noexcept {
  return 2.0 * f * fac;
}

Numeric RosenkranzQuadratic::dNdF0() const noexcept {
  return N * dfacdF0 / fac;
}

Numeric RosenkranzQuadratic::operator()(Numeric f) noexcept {
  N = fac * f * f;
  return N;
}

Numeric SimpleFrequencyScaling::dNdT(Numeric t_ [[maybe_unused]], Numeric f) const ARTS_NOEXCEPT {
  ARTS_ASSERT(t_ == T, "Temperature is stored internally in SimpleFrequencyScaling\n"
  "but for interface reasons this function nevertheless takes temprature as an input\n"
  "For some reason, the two temperatures disagree, so regardless, you have encountered\n"
  "a path of the code that should never be encountered.  The two temperatures are: ",
  T, " and ", t_, " K")
  
  return 
  - N * Constant::h * F0 * expF0 / (Constant::k * t_ * t_ * expm1F0)
  + Constant::h * f * f * std::exp(- (Constant::h * f) / (Constant::k * t_)) /
  (F0 * Constant::k * t_ * t_ * expm1F0);
}

Numeric SimpleFrequencyScaling::dNdf(Numeric f) const noexcept {
  return N / f - N * Constant::h *
  std::exp(- (Constant::h * f) / (Constant::k * T)) / 
  (Constant::k * T * std::expm1(- (Constant::h * f) / (Constant::k * T)));
}

Numeric SimpleFrequencyScaling::operator()(Numeric f) noexcept {
  N = (f / F0) * (std::expm1(- Constant::h * f / (Constant::k * T)) / expm1F0);
  return N;
}

LocalThermodynamicEquilibrium::LocalThermodynamicEquilibrium(
    Numeric I0, Numeric T0, Numeric T, Numeric F0, Numeric E0, Numeric QT,
    Numeric QT0, Numeric dQTdT, Numeric r, Numeric drdSELFVMR, const Numeric drdT) noexcept
    : LocalThermodynamicEquilibrium(
          I0, r, drdSELFVMR, drdT, QT0, QT, dQTdT, boltzman_ratio(T, T0, E0),
          dboltzman_ratio_dT_div_boltzmann_ratio(T, E0),
          stimulated_relative_emission(F0, T0, T),
          dstimulated_relative_emission_dT(F0, T0, T),
          dstimulated_relative_emission_dF0(F0, T0, T)) {
          }
    
struct FullNonLocalThermodynamicEquilibriumInitialization {
  Numeric k, dkdF0, dkdr1, dkdr2,
          e, dedF0, dedr2,
          B, dBdT, dBdF0;
          
  FullNonLocalThermodynamicEquilibriumInitialization(Numeric F0, Numeric A21, Numeric T, Numeric r1, Numeric r2, Numeric c2, Numeric c3, Numeric x) noexcept :
  k(c3 * (r1 * x - r2) * (A21 / c2)),
  dkdF0(- 2.0 * k / F0),
  dkdr1(c3 * x * (A21 / c2)),
  dkdr2(- c3 * (A21 / c2)),
  e(c3 * r2 * A21),
  dedF0(e / F0),
  dedr2(c3 * A21),
  B(2 * Constant::h / Constant::pow2(Constant::c) * Constant::pow3(F0) / std::expm1((Constant::h / Constant::k * F0) / T)),
  dBdT(Constant::pow2(B) * Constant::pow2(Constant::c) * std::exp((Constant::h / Constant::k * F0) / T) / (2*Constant::k*Constant::pow2(F0*T))),
  dBdF0(3 * B / F0 - Constant::pow2(B) * Constant::pow2(Constant::c) * std::exp((Constant::h / Constant::k * F0) / T) / (2*Constant::k*T*Constant::pow3(F0)))
  {}
  
  constexpr FullNonLocalThermodynamicEquilibrium operator()(
    Numeric r, Numeric drdSELFVMR, Numeric drdT
  ) && noexcept {
    return FullNonLocalThermodynamicEquilibrium(r, drdSELFVMR, drdT,
                                                k, dkdF0, dkdr1, dkdr2,
                                                e, dedF0, dedr2,
                                                B, dBdT, dBdF0);
  }
};
    
FullNonLocalThermodynamicEquilibrium::FullNonLocalThermodynamicEquilibrium(
  Numeric F0, Numeric A21, Numeric T,
  Numeric g1, Numeric g2, Numeric r1,
  Numeric r2, Numeric r, Numeric drdSELFVMR, Numeric drdT) noexcept :
  FullNonLocalThermodynamicEquilibrium(FullNonLocalThermodynamicEquilibriumInitialization(F0, A21, T, r1, r2, c0 * F0 * F0 * F0, c1 * F0, g2 / g1)(r, drdSELFVMR, drdT)) {}

struct VibrationalTemperaturesNonLocalThermodynamicEquilibriumInitializer {
  Numeric K1, dK1dT, 
          K2, dK2dT, dK2dF0, 
          K3, dK3dT, dK3dF0, dK3dTl, dK3dTu,
          K4, dK4dT, dK4dTu,
          B, dBdT, dBdF0;
  
  VibrationalTemperaturesNonLocalThermodynamicEquilibriumInitializer(Numeric T, Numeric T0, Numeric F0, Numeric E0, Numeric Tl, Numeric Evl, Numeric Tu, Numeric Evu, Numeric gamma, Numeric gamma_ref, Numeric r_low, Numeric r_upp
  ) noexcept :
  K1(boltzman_ratio(T, T0, E0)),
  dK1dT(dboltzman_ratio_dT(K1, T, E0)),
  K2(stimulated_relative_emission(gamma, gamma_ref)),
  dK2dT(dstimulated_relative_emission_dT(gamma, gamma_ref, F0, T)),
  dK2dF0(dstimulated_relative_emission_dF0(gamma, gamma_ref, T, T0)),
  K3(absorption_nlte_ratio(gamma, r_upp, r_low)),
  dK3dT(dabsorption_nlte_rate_dT(gamma, T, F0, Evl, Evu, K4, r_low)),
  dK3dF0(dabsorption_nlte_rate_dF0(gamma, T, K4, r_low)),
  dK3dTl(dabsorption_nlte_rate_dTl(gamma, T, Tl, Evl, r_low)),
  dK3dTu(dabsorption_nlte_rate_dTu(gamma, T, Tu, Evu, K4)),
  K4(boltzman_ratio(Tu, T, Evu)),
  dK4dT(dboltzman_ratio_dT(K4, T, Evu)),
  dK4dTu(dboltzman_ratio_dT(K4, Tu, Evu)),
  B(2 * Constant::h / Constant::pow2(Constant::c) * Constant::pow3(F0) / std::expm1((Constant::h / Constant::k * F0) / T)),
  dBdT(Constant::pow2(B) * Constant::pow2(Constant::c) * std::exp((Constant::h / Constant::k * F0) / T) / (2*Constant::k*Constant::pow2(F0*T))),
  dBdF0(3 * B / F0 - Constant::pow2(B) * Constant::pow2(Constant::c) * std::exp((Constant::h / Constant::k * F0) / T) / (2*Constant::k*T*Constant::pow3(F0)))
  {}
  
  constexpr VibrationalTemperaturesNonLocalThermodynamicEquilibrium operator()(Numeric I0,
                                                                               Numeric QT0, Numeric QT, Numeric dQTdT,
                                                                               Numeric r, Numeric drdSELFVMR, Numeric drdT) && noexcept {
    return VibrationalTemperaturesNonLocalThermodynamicEquilibrium(I0,
                                                                   QT0, QT, dQTdT,
                                                                   r, drdSELFVMR, drdT,
                                                                   K1, dK1dT, 
                                                                   K2, dK2dT, dK2dF0, 
                                                                   K3, dK3dT, dK3dF0, dK3dTl, dK3dTu,
                                                                   K4, dK4dT, dK4dTu,
                                                                   B, dBdT, dBdF0);
  }
};

VibrationalTemperaturesNonLocalThermodynamicEquilibrium::
    VibrationalTemperaturesNonLocalThermodynamicEquilibrium(
        Numeric I0, Numeric T0, Numeric T, Numeric Tl, Numeric Tu, Numeric F0,
        Numeric E0, Numeric Evl, Numeric Evu, Numeric QT, Numeric QT0,
        Numeric dQTdT, Numeric r, Numeric drdSELFVMR, Numeric drdT) noexcept
    :
    VibrationalTemperaturesNonLocalThermodynamicEquilibrium(
      VibrationalTemperaturesNonLocalThermodynamicEquilibriumInitializer(T, T0, F0, E0, Tl, Evl, Tu, Evu,
                                                                         stimulated_emission(T, F0),
                                                                         stimulated_emission(T0, F0),
                                                                         boltzman_ratio(Tl, T, Evl),
                                                                         boltzman_ratio(Tu, T, Evu)
      )(I0, QT0, QT, dQTdT, r, drdSELFVMR, drdT)
    ) {}

#define CutInternalDerivativesImpl(X, Y)                                      \
  else if (deriv == Jacobian::Line::Shape##X##Y) {                            \
    const Numeric d = value.n;                                                \
    const Complex dFm =                                                       \
      std::conj(ls_mirr.dFd##X(d)  - ls_mirr_cut.dFd##X(d));                  \
    const Complex dFls = ls.dFd##X(d) - ls_cut.dFd##X(d) + dFm;               \
    com.dF[com.jac_pos(iv, ij)] += S * LM * dFls;                             \
    if (do_nlte) {                                                            \
      com.dN[com.jac_pos(iv, ij)] += DS * LM * dFls;                          \
    }                                                                         \
  }

#define CutInternalDerivatives(X)                                             \
CutInternalDerivativesImpl(X, X0) CutInternalDerivativesImpl(X, X1)           \
  CutInternalDerivativesImpl(X, X2) CutInternalDerivativesImpl(X, X3)

#define InternalDerivativesImpl(X, Y)                                         \
  else if (deriv == Jacobian::Line::Shape##X##Y) {                            \
    const Numeric d = value.n;                                                \
    const Complex dFm = std::conj(ls_mirr.dFd##X(d));                         \
    const Complex dFls = ls.dFd##X(d) + dFm;                                  \
    com.dF[com.jac_pos(iv, ij)] += S * LM * dFls;                             \
    if (do_nlte) {                                                            \
      com.dN[com.jac_pos(iv, ij)] += DS * LM * dFls;                          \
    }                                                                         \
  }

#define InternalDerivatives(X)                                                \
  InternalDerivativesImpl(X, X0) InternalDerivativesImpl(X, X1)               \
    InternalDerivativesImpl(X, X2) InternalDerivativesImpl(X, X3)

#define InternalDerivativesSetupImpl(X, Y)                                    \
  else if (deriv == Jacobian::Line::Shape##X##Y) {                            \
    const Index pos =                                                         \
      band.BroadeningSpeciesPosition(deriv.Target().LineSpecies());           \
    if (pos >= 0) {                                                           \
    derivs[ij].value.n = band.Line(i).LineShape().d##X##_d##Y(                \
        T, band.T0(), P, pos, vmrs);                                          \
    } else {                                                                  \
      derivs[ij].value.n = 0;                                                 \
    }                                                                         \
  }

#define InternalDerivativesSetup(X)                                           \
  InternalDerivativesSetupImpl(X, X0) InternalDerivativesSetupImpl(X, X1)     \
    InternalDerivativesSetupImpl(X, X2) InternalDerivativesSetupImpl(X, X3)

#define InternalDerivativesGImpl(X)                                           \
  else if (deriv == Jacobian::Line::ShapeG##X) {                              \
    const Numeric dLM = value.n;                                              \
    com.dF[com.jac_pos(iv, ij)] += dLM * S * Fls;                             \
    if (do_nlte) {                                                            \
      com.dN[com.jac_pos(iv, ij)] += dLM * DS * Fls;                          \
    } \
  }

#define InternalDerivativesG                                                  \
  InternalDerivativesGImpl(X0) InternalDerivativesGImpl(X1)                   \
    InternalDerivativesGImpl(X2) InternalDerivativesGImpl(X3)

#define InternalDerivativesYImpl(X)                                           \
  else if (deriv == Jacobian::Line::ShapeY##X) {                              \
    const Complex dLM = Complex(0, -value.n);                                 \
    com.dF[com.jac_pos(iv, ij)] += dLM * S * Fls;                             \
    if (do_nlte) {                                                            \
      com.dN[com.jac_pos(iv, ij)] += dLM * DS * Fls;                          \
    }                                                                         \
  }

#define InternalDerivativesY                                                  \
  InternalDerivativesYImpl(X0) InternalDerivativesYImpl(X1)                   \
    InternalDerivativesYImpl(X2) InternalDerivativesYImpl(X3)

//! Struct to keep the cutoff limited range values
struct CutoffRange {
  Index start, size;
};

/** Gets the start and size of a range such that
 *
 * \f[ fl \leq f_grid[i] \leq fu \f]
 *
 * for i such that all fl <= f_grid[Range(out.first, out.second, 1)] <= fu
 *
 * @param[in] fl Lower frequency limit
 * @param[in] fu Upper frequency limit
 * @param[in] f_grid As WSV, must be sorted
 * @return out so that the Range above can be formed
 */
CutoffRange limited_range(const Numeric fl, const Numeric fu, const Vector &f_grid) ARTS_NOEXCEPT {
  ARTS_ASSERT(fu > fl);
  const Numeric * it0 = f_grid.get_c_array();
  const Numeric * itn = it0 + f_grid.size();
  const Numeric * itl = std::lower_bound(it0, itn, fl);
  return CutoffRange{std::distance(it0, itl), std::distance(itl, std::upper_bound(itl, itn, fu))};
}

//! Struct to keep the cutoff limited range values and the sparse limits
struct SparseLimitRange {
  Index start, size;
  Index start_low, size_low;
  Index start_upp, size_upp;
};

SparseLimitRange linear_sparse_limited_range(const Numeric flc,
                                             const Numeric fuc,
                                             const Numeric fls,
                                             const Numeric fus,
                                             const Vector& f_grid,
                                             const Vector& sparse_f_grid) ARTS_NOEXCEPT {
  ARTS_ASSERT(fls > flc);
  ARTS_ASSERT(fus > fls);
  ARTS_ASSERT(fuc > fus);
  
  const Index nvs = sparse_f_grid.size();
  
  // Find bounds in sparse
  const Numeric * it0s = sparse_f_grid.get_c_array();
  const Numeric * itns = it0s + nvs;
  const Numeric * itlc = std::lower_bound(it0s, itns, std::nextafter(flc, fuc));  // lower cutoff
  const Numeric * ituc = std::upper_bound(itlc, itns, std::nextafter(fuc, flc));  // upper cutoff
  const Numeric * itls = std::upper_bound(itlc, ituc, fls);  // lower sparse
  const Numeric * itus = std::lower_bound(itls, ituc, fus);  // upper sparse
  
  /* Start and size of sparse adjusted to the 2-grid
   * 
   * The interface to the dense grid is altered slightly so that
   * a complete set-of-2 quadratic interopolation points becomes
   * a guarantee.  Their start/end points ignore this restriction
   * but have been defined to not contain anything extra
   */
  Index beg_lr = std::distance(it0s, itlc); /*while (beg_lr % 2) --beg_lr;*/
  Index end_lr = std::distance(it0s, itls); while (end_lr % 2) --end_lr;
  Index beg_ur = std::distance(it0s, itus); while (beg_ur % 2) ++beg_ur;
  Index end_ur = std::distance(it0s, ituc); /*while (end_ur % 2) ++end_ur;*/
  
  // Find new limits
  const Numeric fl = (end_lr <= 0 or end_lr >= nvs) ? flc : sparse_f_grid[end_lr];
  const Numeric fu = (beg_ur <= 0 or beg_ur >= nvs) ? fuc : sparse_f_grid[beg_ur];
  
  // Find bounds in dense
  const Numeric * it0 = f_grid.get_c_array();
  const Numeric * itn = it0 + f_grid.size();
  const Numeric * itl = std::lower_bound(it0, itn, fl);  // include boundary
  const Numeric * itu = std::upper_bound(itl, itn, std::nextafter(fu, fl));  // dismiss boundary
  
  return {
    std::distance(it0, itl), std::distance(itl, itu),
    beg_lr, end_lr - beg_lr,
    beg_ur, end_ur - beg_ur
  };
}

SparseLimitRange quad_sparse_limited_range(const Numeric flc,
                                           const Numeric fuc,
                                           const Numeric fls,
                                           const Numeric fus,
                                           const Vector& f_grid,
                                           const Vector& sparse_f_grid) ARTS_NOEXCEPT {
  ARTS_ASSERT(fls > flc);
  ARTS_ASSERT(fus > fls);
  ARTS_ASSERT(fuc > fus);
  
  const Index nvs = sparse_f_grid.size();
  const Index nv = f_grid.size();
  
  // Find bounds in sparse
  const Numeric * const it0s = sparse_f_grid.get_c_array();
  const Numeric * const itns = it0s + nvs;
  const Numeric * itlc = std::lower_bound(it0s, itns, std::nextafter(flc, fuc));  // lower cutoff
  const Numeric * ituc = std::upper_bound(itlc, itns, std::nextafter(fuc, flc));  // upper cutoff
  const Numeric * itls = std::upper_bound(itlc, ituc, std::nextafter(fls, flc));  // lower sparse
  const Numeric * itus = std::lower_bound(itls, ituc, std::nextafter(fus, fuc));  // upper sparse
  
  /* Start and size of sparse adjusted to the 3-grid
   * 
   * The interface to the dense grid is altered slightly so that
   * a complete set-of-3 quadratic interpolation points becomes
   * a guarantee.  Their end points are modified to contain
   * only set-of-3.  If the range is not modified correctly,
   * you risk having negative interpolation out of your range.
   * 
   * To avoid negative absorption entirely, the range between
   * the true cutoff and the outer-most sparse grid point must
   * have at least two values.  If they have only one value,
   * the quadratic interpolation will lead to negative values
   * in the outer range.  So there's a small range that has to
   * be ignored
   */
  while (std::distance(it0s, itls) % 3) --itls;
  while (std::distance(itlc, itls) % 3 == 1) ++itlc;  // skip some cutoff
  while (std::distance(it0s, itus) % 3) ++itus;
  while (std::distance(itus, ituc) % 3 == 1) --ituc;  // skip some cutoff
  
  // Find bounds in dense
  const Numeric * const it0 = f_grid.get_c_array();
  const Numeric * const itn = it0 + nv;
  const Numeric * itl;
  const Numeric * itu;
  
  if (itls not_eq itns) {
    itl = std::lower_bound(it0, itn, *itls);  // include boundary
  } else {
    itl = std::lower_bound(it0, itn, flc);  // include boundary
  }
  
  if (itus not_eq itns and itl not_eq itn) {
    itu = std::upper_bound(itl, itn, std::nextafter(*itus, *itl));  // dismiss boundary
  } else {
    itu = std::lower_bound(itl, itn, std::nextafter(fuc, flc));  // include boundary
  }
  
  return {
    std::distance(it0, itl), std::distance(itl, itu),
    std::distance(it0s, itlc), std::distance(itlc, itls),
    std::distance(it0s, itus), std::distance(itus, ituc)
  };
}

/** Data struct for keeping derivative keys and values
 *
 * If ever more types of values are necessary, please
 * add the type to the Values union
 */
struct Derivatives {
  Absorption::QuantumIdentifierLineTarget target;
  union Values {
    Output o;
    Numeric n;
    constexpr Values() noexcept : n() {}
  } value;
  Index jac_pos{};
  const RetrievalQuantity * deriv{nullptr};
  Derivatives() noexcept : target(), value() {}
};

//! Helper to keep function signature clean
using ArrayOfDerivatives = Array<Derivatives>;

//! Helper function to find the last relevant derivative
Index active_nelem(const ArrayOfDerivatives& derivs) noexcept {
  return std::distance(derivs.cbegin(),
                      std::find_if(derivs.cbegin(), derivs.cend(),
                                    [](auto& deriv) {
                                      return deriv.deriv == nullptr;
                                    }));
}

struct ComputeValues {
  Complex * const F;
  Complex * const dF;
  Complex * const N;
  Complex * const dN;
  const Numeric * const f;
  const Index size;
  const ArrayOfDerivatives& derivs;
  const Index jac_size;
  const Index max_jac_size;
  const bool do_nlte;
  
  [[nodiscard]] Index jac_pos(Index iv, Index ij) const noexcept {return jac_size * iv + ij;}
  
  ComputeValues(ComplexVector &F_,
                ComplexMatrix &dF_,
                ComplexVector &N_,
                ComplexMatrix &dN_,
                const Vector &f_grid,
                const Index start,
                const Index nv,
                const ArrayOfDerivatives& derivs_,
                const bool do_nlte_) noexcept :
  F(F_.get_c_array()+start),
  dF(dF_.get_c_array()+start*derivs_.size()),
  N(N_.get_c_array()+start),
  dN(dN_.get_c_array()+start*derivs_.size()),
  f(f_grid.get_c_array()+start),
  size(nv), derivs(derivs_), jac_size(derivs_.size()), max_jac_size(active_nelem(derivs)), do_nlte(do_nlte_) {
  }
  
  ComputeValues(Complex &F_, std::vector<Complex> &dF_,
                Complex &N_, std::vector<Complex> &dN_,
                const Numeric &f_lim,
                const ArrayOfDerivatives& derivs_, const bool do_nlte_) noexcept :
  F(&F_), dF(dF_.data()), N(&N_), dN(dN_.data()), f(&f_lim), size(1),
  derivs(derivs_), jac_size(derivs.nelem()), max_jac_size(active_nelem(derivs)), do_nlte(do_nlte_) {
  }
  
  ComputeValues& operator-=(const ComputeValues& cut) ARTS_NOEXCEPT {
    ARTS_ASSERT(cut.size == 1, "Not a cutoff limit")
    ARTS_ASSERT(cut.jac_size == jac_size, "Not from the same Jacobian type")
    
    for (Index iv=0; iv<size; iv++) {
      F[iv] -= *cut.F;
      for (Index ij=0; ij<max_jac_size; ij++) {
        dF[jac_pos(iv, derivs[ij].jac_pos)] -= cut.dF[derivs[ij].jac_pos];
      }
    }
    if (do_nlte) {
      for (Index iv=0; iv<size; iv++) {
        N[iv] -= *cut.N;
        for (Index ij=0; ij<max_jac_size; ij++) {
          dN[jac_pos(iv, derivs[ij].jac_pos)] -= cut.dN[derivs[ij].jac_pos];
        }
      }
    }
    return *this;
  }
};

/** Cutoff frequency loop of the line shape call
 *
 * This simply adds to the four output vectors/matrices for
 * each frequency grid.
 *
 * The equation solved is the inner part of
 *
 * \f[
 * F(f) = \sum_i S_{lm} S_n(f) S_z S_i \left[F_i\left(f\right) +
 * F^M_i\left(f\right) - F_i\left(f_c\right) - F^M_i\left(f_c\right)\right], \f]
 *
 * where \f$ f \f$ is the frequency, \f$ i \f$ is the line index,
 * \f$ S_{lm} \f$ is the line mixing scaling,
 * \f$ S_n(f) \f$ is the line normalization at \f$ f \f$, \f$ S_z \f$
 * is the Zeeman line strength scaling, \f$ S_i \f$ is the line strength,
 * \f$ F_i\left(f\right) \f$ is the line shape at \f$ f \f$,
 * \f$ F^M_i\left(f\right) \f$ is the mirrored line shape at \f$ f \f$,
 * and \f$ f_c \f$ is the cutoff frequency.
 *
 * Each of \f$ S_n \f$, \f$ S_i \f$, \f$ F_i\left(f\right) \f$, and
 * \f$ F^M_i\left(f\right) \f$ are setup to work as functional classes.
 * They each define a set of possible derivatives and these possible
 * derivatives are used to setup the chain rules required to computed
 * \f$ \partial F(f) / \partial x \f$.  Also, \f$ S_i \f$ is setup so
 * to return the NLTE ratios required for computing
 *
 * \f[
 * N(f) = \sum_i S_{lm} S_n(f) S_z N_i \left[F_i\left(f\right) +
 * F^M_i\left(f\right) - F_i\left(f_c\right) - F^M_i\left(f_c\right)\right], \f]
 *
 * with the same rules applied to \f$ \partial N(f) / \partial x \f$ as on
 * \f$ \partial F(f) / \partial x \f$.
 *
 * Note that even though not directly implied by this notation, most of the
 * variables used in these two equations depends on some line parameters.  Also
 * note that the idea is to add full cross-section at once so as to limit the
 * number of allocations in nested runs, so the calculations must happen in
 * the inner most loop.
 *
 * This function is not possible to run on multiple cores.  Such parallelisms
 * must happen at a much higher level.
 *
 * @param[in] jacobian_quantities As WSV
 * @param[in] T The atmospheric temperature
 * @param[in] do_nlte Flag for whether or not NLTE will be computed
 * @param[in] LM The line mixing scaling. \f$ S_{lm} \f$
 * @param[in] ls_str The line strength calculator. \f$ S_i \f$
 * @param[in,out] ls_norm The normalization calculator. \f$ S_n \f$
 * @param[in] derivs A list of pre-computed derivative values and keys
 * @param[in,out] F The cross-section. \f$ F \f$
 * @param[in,out] dF The cross-section's derivatives. \f$ \partial F / \partial
 * x \f$
 * @param[in,out] N The cross-section ratio of the NLTE source. \f$ N \f$
 * @param[in,out] dN The cross-section ratio of the NLTE source's derivatives.
 * \f$ \partial N / \partial x \f$
 * @param[in] f_grid The frequency grid. \f$ \left[ f_0, \cdots, f_n \right] \f$
 * @param[in] dfdH The derivative of the change of frequency w.r.t. magnetic
 * field strength
 * @param[in] Sz The relative Zeeman strength. \f$ S_z \f$
 * @param[in] ls The line shape calculator. \f$ F_i \f$
 * @param[in] ls_mirr The mirrored line shape calculator. \f$ F^M_i \f$
 */
void cutoff_frequency_loop(ComputeValues &com,
                           Calculator &ls,
                           Calculator &ls_mirr,
                           Normalizer &ls_norm,
                           const Calculator &ls_cut,
                           const Calculator &ls_mirr_cut,
                           const IntensityCalculator &ls_str,
                           const ArrayOfDerivatives &derivs, 
                           const Complex LM,
                           const Numeric &T, 
                           const Numeric &dfdH,
                           const Numeric &Sz, 
                           const Species::Species self_species) ARTS_NOEXCEPT {
  const Index nv = com.size;
  const bool do_nlte = com.do_nlte;
  
  const Numeric Si = ls_str.S();
  const Numeric DNi = ls_str.N();

  for (Index iv = 0; iv < nv; iv++) {
    const Numeric f = com.f[iv];

    const Numeric Sn = ls_norm(f);
    const Numeric S = Sz * Sn * Si;
    const Numeric DS = Sz * Sn * DNi;
    const Complex Fm = std::conj(ls_mirr(f) - ls_mirr_cut.F());
    const Complex Fls = ls(f) - ls_cut.F() + Fm;
    com.F[iv] += S * LM * Fls;
    if (do_nlte) {
      com.N[iv] += DS * LM * Fls;
    }
    
    for (const auto& [lt, value, ij, deriv_ptr]: derivs) {
      if (not deriv_ptr) break;
      const auto& deriv = *deriv_ptr;

      if (deriv == Jacobian::Atm::Temperature) {
        const auto &dXdT = value.o;
        const Numeric dSn = ls_norm.dNdT(T, f);
        const Numeric dSi = ls_str.dSdT();
        const Complex dLM(dXdT.G, -dXdT.Y);
        const Complex dFm = std::conj(ls_mirr.dFdT(mirroredOutput(dXdT), T) - ls_mirr_cut.dFdT(mirroredOutput(dXdT), T));
        const Complex dFls = ls.dFdT(dXdT, T) - ls_cut.dFdT(dXdT, T) + dFm;
        com.dF[com.jac_pos(iv, ij)] +=
            S * (LM * dFls + dLM * Fls) + Sz * (dSn * Si + Sn * dSi) * LM * Fls;
        if (do_nlte) {
          const Numeric dDSi = ls_str.dNdT();
          com.dN[com.jac_pos(iv, ij)] += DS * (LM * dFls + dLM * Fls) +
                        Sz * (dSn * Si + Sn * dDSi) * LM * Fls;
        }
      } else if (deriv.is_wind()) {
        const Complex dFm = std::conj(ls_mirr.dFdf() - ls_mirr_cut.dFdf());
        const Complex dFls = ls.dFdf() - ls_cut.dFdf() + dFm;
        const Numeric dS = Sz * ls_norm.dNdf(f) * Si;
        com.dF[com.jac_pos(iv, ij)] += S * LM * dFls + dS * LM * Fls;
        if (do_nlte) {
          const Numeric dDS = Sz * ls_norm.dNdf(f) * DNi;
          com.dN[com.jac_pos(iv, ij)] += DS * LM * dFls + dDS * LM * Fls;
        }
      } else if (deriv.is_mag()) {
        const Complex dFm = std::conj(ls_mirr.dFdH(-dfdH) - ls_mirr_cut.dFdH(-dfdH));
        const Complex dFls = ls.dFdH(dfdH) - ls_cut.dFdH(dfdH) + dFm;
        com.dF[com.jac_pos(iv, ij)] += S * LM * dFls;
        if (do_nlte) {
          com.dN[com.jac_pos(iv, ij)] += DS * LM * dFls;
        }
      } else if (deriv == Jacobian::Special::ArrayOfSpeciesTagVMR) {
        com.dF[com.jac_pos(iv, ij)] += Sz * Sn * ls_str.dSdSELFVMR() * LM * Fls;
        if (do_nlte) {
          com.dN[com.jac_pos(iv, ij)] += Sz * Sn * ls_str.dNdSELFVMR() * LM * Fls;
        }
      } else if (deriv.Target().needQuantumIdentity()) {
        if (deriv == Jacobian::Line::VMR) {
          const auto &dXdVMR = value.o;
          const Complex dFm = std::conj(ls_mirr.dFdVMR(mirroredOutput(dXdVMR)) - ls_mirr_cut.dFdVMR(mirroredOutput(dXdVMR)));
          const Complex dLM(dXdVMR.G, -dXdVMR.Y);
          const Complex dFls = ls.dFdVMR(dXdVMR) - ls_cut.dFdVMR(dXdVMR) + dFm;
          com.dF[com.jac_pos(iv, ij)] += S * LM * dFls + dLM * S * Fls;
          if (self_species == deriv.QuantumIdentity().Species()) {
            com.dF[com.jac_pos(iv, ij)] += ls_str.dSdSELFVMR() * LM * Fls;
          }
          if (do_nlte) {
            com.dN[com.jac_pos(iv, ij)] += DS * LM * dFls + dLM * DS * Fls;
            if (self_species == deriv.QuantumIdentity().Species()) {
              com.dF[com.jac_pos(iv, ij)] += ls_str.dNdSELFVMR() * LM * Fls;
            }
          }
        } else {
          if (lt == Absorption::QuantumIdentifierLineTargetType::Line) {
            if (deriv == Jacobian::Line::Center) {
              const Numeric d = ls_norm.dNdF0() * Si + ls_str.dSdF0() * Sn;
              const Complex dFm = std::conj(ls_mirr.dFdF0() - ls_mirr_cut.dFdF0());
              const Complex dFls = ls.dFdF0() - ls_cut.dFdF0() + dFm;
              const Numeric dS = Sz * d;
              com.dF[com.jac_pos(iv, ij)] += S * LM * dFls + dS * LM * Fls;
              if (do_nlte) {
                const Numeric Dd = ls_norm.dNdF0() * DNi + ls_str.dNdF0() * Sn;
                const Numeric DdS = Sz * Dd;
                com.dN[com.jac_pos(iv, ij)] += DS * LM * dFls + DdS * LM * Fls;
              }
            } else if (deriv == Jacobian::Line::Strength) {
              const Numeric dS = Sz * Sn * ls_str.dSdI0();
              com.dF[com.jac_pos(iv, ij)] += dS * LM * Fls;
              if (do_nlte) {
                const Numeric DdS = Sz * Sn * ls_str.dNdI0();
                com.dN[com.jac_pos(iv, ij)] += DdS * LM * Fls;
              }
            }
            // All line shape derivatives
            CutInternalDerivatives(G0)
            CutInternalDerivatives(D0)
            CutInternalDerivatives(G2)
            CutInternalDerivatives(D2)
            CutInternalDerivatives(ETA)
            CutInternalDerivatives(FVC)
            InternalDerivativesY
            InternalDerivativesG
            CutInternalDerivatives(DV)
          } else if (lt == Absorption::QuantumIdentifierLineTargetType::Level) {
            if (deriv == Jacobian::Line::NLTE) {
              if (lt.upper) {
                const Numeric dS = Sz * Sn * ls_str.dSdNLTEu();
                com.dF[com.jac_pos(iv, ij)] += dS * LM * Fls;
                
                if (do_nlte) {
                  const Numeric DdS = Sz * Sn * ls_str.dNdNLTEu();
                  com.dN[com.jac_pos(iv, ij)] += DdS * LM * Fls;
                }
              }

              if (lt.lower) {
                const Numeric dS = Sz * Sn * ls_str.dSdNLTEl();
                com.dF[com.jac_pos(iv, ij)] += dS * LM * Fls;
                
                if (do_nlte) {
                  const Numeric DdS = Sz * Sn * ls_str.dNdNLTEl();
                  com.dN[com.jac_pos(iv, ij)] += DdS * LM * Fls;
                }
              }
            }
          }
        }
      }
    }
  }
}

/** Frequency loop of the line shape call
 *
 * This simply adds to the four output vectors/matrices for
 * each frequency grid.
 *
 * The equation solved is the inner part of
 *
 * \f[
 * F(f) = \sum_i S_{lm} S_n(f) S_z S_i \left[F_i\left(f\right) +
 * F^M_i\left(f\right)\right], \f]
 *
 * where \f$ f \f$ is the frequency, \f$ i \f$ is the line index,
 * \f$ S_{lm} \f$ is the line mixing scaling,
 * \f$ S_n(f) \f$ is the line normalization at \f$ f \f$, \f$ S_z \f$
 * is the Zeeman line strength scaling, \f$ S_i \f$ is the line strength,
 * \f$ F_i\left(f\right) \f$ is the line shape at \f$ f \f$, and
 * \f$ F^M_i\left(f\right) \f$ is the mirrored line shape at \f$ f \f$.
 *
 * Each of \f$ S_n \f$, \f$ S_i \f$, \f$ F_i\left(f\right) \f$, and
 * \f$ F^M_i\left(f\right) \f$ are setup to work as functional classes.
 * They each define a set of possible derivatives and these possible
 * derivatives are used to setup the chain rules required to computed
 * \f$ \partial F(f) / \partial x \f$.  Also, \f$ S_i \f$ is setup so
 * to return the NLTE ratios required for computing
 *
 * \f[
 * N(f) = \sum_i S_{lm} S_n(f) S_z N_i \left[F_i\left(f\right) +
 * F^M_i\left(f\right)\right], \f]
 *
 * with the same rules applied to \f$ \partial N(f) / \partial x \f$ as on
 * \f$ \partial F(f) / \partial x \f$.
 *
 * Note that even though not directly implied by this notation, most of the
 * variables used in these two equations depends on some line parameters.  Also
 * note that the idea is to add full cross-section at once so as to limit the
 * number of allocations in nested runs, so the calculations must happen in
 * the inner most loop.
 *
 * This function is not possible to run on multiple cores.  Such parallelisms
 * must happen at a much higher level.
 *
 * @param[in] jacobian_quantities As WSV
 * @param[in] T The atmospheric temperature
 * @param[in] do_nlte Flag for whether or not NLTE will be computed
 * @param[in] LM The line mixing scaling. \f$ S_{lm} \f$
 * @param[in] ls_str The line strength calculator. \f$ S_i \f$
 * @param[in,out] ls_norm The normalization calculator. \f$ S_n \f$
 * @param[in] derivs A list of pre-computed derivative values and keys
 * @param[in,out] F The cross-section. \f$ F \f$
 * @param[in,out] dF The cross-section's derivatives. \f$ \partial F / \partial
 * x \f$
 * @param[in,out] N The cross-section ratio of the NLTE source. \f$ N \f$
 * @param[in,out] dN The cross-section ratio of the NLTE source's derivatives.
 * \f$ \partial N / \partial x \f$
 * @param[in] f_grid The frequency grid. \f$ \left[ f_0, \cdots, f_n \right] \f$
 * @param[in] dfdH The derivative of the change of frequency w.r.t. magnetic
 * field strength
 * @param[in] Sz The relative Zeeman strength. \f$ S_z \f$
 * @param[in] ls The line shape calculator. \f$ F_i \f$
 * @param[in] ls_mirr The mirrored line shape calculator. \f$ F^M_i \f$
 */
void frequency_loop(ComputeValues &com,
                    Calculator &ls,
                    Calculator &ls_mirr,
                    Normalizer &ls_norm,
                    const IntensityCalculator &ls_str, 
                    const ArrayOfDerivatives &derivs, 
                    const Complex LM,
                    const Numeric &T, 
                    const Numeric &dfdH,
                    const Numeric &Sz, 
                    const Species::Species self_species) ARTS_NOEXCEPT {
  const Index nv = com.size;
  const bool do_nlte = com.do_nlte;
  
  const Numeric Si = ls_str.S();
  const Numeric DNi = ls_str.N();

  for (Index iv = 0; iv < nv; iv++) {
    const Numeric f = com.f[iv];

    const Numeric Sn = ls_norm(f);
    const Numeric S = Sz * Sn * Si;
    const Numeric DS = Sz * Sn * DNi;
    const Complex Fm = std::conj(ls_mirr(f));
    const Complex Fls = ls(f) + Fm;
    com.F[iv] += S * LM * Fls;
    if (do_nlte) {
      com.N[iv] += DS * LM * Fls;
    }
    
    for (const auto& [lt, value, ij, deriv_ptr]: derivs) {
      if (not deriv_ptr) break;
      const auto& deriv = *deriv_ptr;

      if (deriv == Jacobian::Atm::Temperature) {
        const auto &dXdT = value.o;
        const Numeric dSn = ls_norm.dNdT(T, f);
        const Numeric dSi = ls_str.dSdT();
        const Complex dLM(dXdT.G, -dXdT.Y);
        const Complex dFm = std::conj(ls_mirr.dFdT(mirroredOutput(dXdT), T));
        const Complex dFls = ls.dFdT(dXdT, T) + dFm;
        com.dF[com.jac_pos(iv, ij)] +=
            S * (LM * dFls + dLM * Fls) + Sz * (dSn * Si + Sn * dSi) * LM * Fls;
        if (do_nlte) {
          const Numeric dDSi = ls_str.dNdT();
          com.dN[com.jac_pos(iv, ij)] += DS * (LM * dFls + dLM * Fls) +
                        Sz * (dSn * Si + Sn * dDSi) * LM * Fls;
        }
      } else if (deriv.is_wind()) {
        const Complex dFm = std::conj(ls_mirr.dFdf());
        const Complex dFls = ls.dFdf() + dFm;
        const Numeric dS = Sz * ls_norm.dNdf(f) * Si;
        com.dF[com.jac_pos(iv, ij)] += S * LM * dFls + dS * LM * Fls;
        if (do_nlte) {
          const Numeric dDS = Sz * ls_norm.dNdf(f) * DNi;
          com.dN[com.jac_pos(iv, ij)] += DS * LM * dFls + dDS * LM * Fls;
        }
      } else if (deriv.is_mag()) {
        const Complex dFm = std::conj(ls_mirr.dFdH(-dfdH));
        const Complex dFls = ls.dFdH(dfdH) + dFm;
        com.dF[com.jac_pos(iv, ij)] += S * LM * dFls;
        if (do_nlte) {
          com.dN[com.jac_pos(iv, ij)] += DS * LM * dFls;
        }
      } else if (deriv == Jacobian::Special::ArrayOfSpeciesTagVMR) {
        com.dF[com.jac_pos(iv, ij)] += Sz * Sn * ls_str.dSdSELFVMR() * LM * Fls;
        if (do_nlte) {
              com.dN[com.jac_pos(iv, ij)] += Sz * Sn * ls_str.dNdSELFVMR() * LM * Fls;
        }
      } else if (deriv.Target().needQuantumIdentity()) {
        if (deriv == Jacobian::Line::VMR) {
          const auto &dXdVMR = value.o;
          const Complex dFm = std::conj(ls_mirr.dFdVMR(mirroredOutput(dXdVMR)));
          const Complex dLM(dXdVMR.G, -dXdVMR.Y);
          const Complex dFls = ls.dFdVMR(dXdVMR) + dFm;
          com.dF[com.jac_pos(iv, ij)] += S * LM * dFls + dLM * S * Fls;
          if (self_species == deriv.QuantumIdentity().Species()) {
            com.dF[com.jac_pos(iv, ij)] += ls_str.dSdSELFVMR() * LM * Fls;
          }
          if (do_nlte) {
            com.dN[com.jac_pos(iv, ij)] += DS * LM * dFls + dLM * DS * Fls;
            if (self_species == deriv.QuantumIdentity().Species()) {
              com.dN[com.jac_pos(iv, ij)] += ls_str.dNdSELFVMR() * LM * Fls;
            }
          }
        } else {
          if (lt == Absorption::QuantumIdentifierLineTargetType::Line) {
            if (deriv == Jacobian::Line::Center) {
              const Numeric d = ls_norm.dNdF0() * Si + ls_str.dSdF0() * Sn;
              const Complex dFm = std::conj(ls_mirr.dFdF0());
              const Complex dFls = ls.dFdF0() + dFm;
              const Numeric dS = Sz * d;
              com.dF[com.jac_pos(iv, ij)] += S * LM * dFls + dS * LM * Fls;
              if (do_nlte) {
                const Numeric Dd = ls_norm.dNdF0() * DNi + ls_str.dNdF0() * Sn;
                const Numeric DdS = Sz * Dd;
                com.dN[com.jac_pos(iv, ij)] += DS * LM * dFls + DdS * LM * Fls;
              }
            } else if (deriv == Jacobian::Line::Strength) {
              const Numeric dS = Sz * Sn * ls_str.dSdI0();
              com.dF[com.jac_pos(iv, ij)] += dS * LM * Fls;
              if (do_nlte) {
                const Numeric DdS = Sz * Sn * ls_str.dNdI0();
                com.dN[com.jac_pos(iv, ij)] += DdS * LM * Fls;
              }
            }
            // All line shape derivatives
            InternalDerivatives(G0)
            InternalDerivatives(D0)
            InternalDerivatives(G2)
            InternalDerivatives(D2)
            InternalDerivatives(ETA)
            InternalDerivatives(FVC)
            InternalDerivativesY
            InternalDerivativesG
            InternalDerivatives(DV)
          } else if (lt == Absorption::QuantumIdentifierLineTargetType::Level) {
            if (deriv == Jacobian::Line::NLTE) {
              if (lt.upper) {
                const Numeric dS = Sz * Sn * ls_str.dSdNLTEu();
                com.dF[com.jac_pos(iv, ij)] += dS * LM * Fls;
                
                if (do_nlte) {
                  const Numeric DdS = Sz * Sn * ls_str.dNdNLTEu();
                  com.dN[com.jac_pos(iv, ij)] += DdS * LM * Fls;
                }
              }

              if (lt.lower) {
                const Numeric dS = Sz * Sn * ls_str.dSdNLTEl();
                com.dF[com.jac_pos(iv, ij)] += dS * LM * Fls;
                
                if (do_nlte) {
                  const Numeric DdS = Sz * Sn * ls_str.dNdNLTEl();
                  com.dN[com.jac_pos(iv, ij)] += DdS * LM * Fls;
                }
              }
            }
          }
        }
      }
    }
  }
}

/** Cutoff considerations of the line shape
 *
 * This function takes care of setting up cutoff and line shape considerations
 * for the frequency loop function it wraps.  Internally, the cutoff is
 * calculated from the band information and a view of the correct data is sent
 * on.
 *
 * The Zeeman effect is also considered internally if applicable (or ignored
 * otherwise)
 *
 * This function is not possible to run on multiple cores.  Such parallelisms
 * must happen at a much higher level.
 *
 * @param[in,out] F The cross-section. \f$ F \f$
 * @param[in,out] dF The cross-section's derivatives. \f$ \partial F / \partial
 * x \f$
 * @param[in,out] N The cross-section ratio of the NLTE source. \f$ N \f$
 * @param[in,out] dN The cross-section ratio of the NLTE source's derivatives.
 * \f$ \partial N / \partial x \f$
 * @param[in] f_grid The frequency grid. \f$ \left[ f_0, \cdots, f_n \right] \f$
 * @param[in] band The absorption band
 * @param[in] jacobian_quantities As WSV
 * @param[in] T The atmospheric temperature
 * @param[in] do_nlte Flag for whether or not NLTE will be computed
 * @param[in] H The magnetic field magnitude
 * @param[in] do_zeeman Flag for whether this is part of some Zeeman
 * calculations
 * @param[in] zeeman_polarization The type of Zeeman polarization to consider
 * (if any)
 * @param[in] f_mean The mean frequency of the absorption band
 * @param[in] DC The Doppler broadening constant of the band
 * @param[in] i The line index
 * @param[in] X The line shape model parameters of the atmosphere
 * @param[in] LM The line mixing scaling. \f$ S_{lm} \f$
 * @param[in] ls_str The line strength calculator. \f$ S_i \f$
 * @param[in] ls_norm The normalization calculator. \f$ S_n \f$
 * @param[in] derivs A list of pre-computed derivative values and keys
 */
void cutoff_loop_sparse_linear(ComputeData &com,
                               ComputeData &sparse_com,
                               Normalizer ls_norm,
                               const IntensityCalculator ls_str,
                               const AbsorptionLines &band,
                               const ArrayOfDerivatives &derivs,
                               const Output X,
                               const Numeric &T,
                               const Numeric &H,
                               const Numeric &sparse_lim,
                               const Numeric &DC,
                               const Index i,
                               const bool do_zeeman,
                               const Zeeman::Polarization zeeman_polarization) ARTS_NOEXCEPT {
  // Basic settings
  const bool do_nlte = ls_str.do_nlte();
  const bool do_cutoff = band.Cutoff() not_eq Absorption::CutoffType::None;
  const Numeric fu = band.CutoffFreq(i, X.D0);
  const Numeric fl = band.CutoffFreqMinus(i, X.D0);
  const Numeric fus = band.F0(i) + sparse_lim;
  const Numeric fls = band.F0(i) - sparse_lim;

  // Find sparse and dense ranges
  const auto [dense_start, dense_size,
              sparse_low_start, sparse_low_size,
              sparse_upp_start, sparse_upp_size] = linear_sparse_limited_range(fl, fu, fls, fus,
                                                                               com.f_grid,
                                                                               sparse_com.f_grid);
  if ((dense_size + sparse_low_size + sparse_upp_size) == 0) return;
  
  // Get the compute data view
  ComputeValues comval(com.F, com.dF, com.N, com.dN, com.f_grid, dense_start, dense_size, derivs, do_nlte);
  
  // Get views of the sparse data
  ComputeValues sparse_low_range(sparse_com.F, sparse_com.dF, sparse_com.N, sparse_com.dN, sparse_com.f_grid, sparse_low_start, sparse_low_size, derivs, do_nlte);
  ComputeValues sparse_upp_range(sparse_com.F, sparse_com.dF, sparse_com.N, sparse_com.dN, sparse_com.f_grid, sparse_upp_start, sparse_upp_size, derivs, do_nlte);
  
  const Index nz = do_zeeman ? band.ZeemanCount(i, zeeman_polarization) : 1;
  for (Index iz = 0; iz < nz; iz++) {
    const Numeric dfdH = do_zeeman ? band.ZeemanSplitting(i, zeeman_polarization, iz) : 0;
    const Numeric Sz = do_zeeman ? band.ZeemanStrength(i, zeeman_polarization, iz) : 1;
    const Complex LM = Complex(1 + X.G, -X.Y);
    Calculator ls(band.LineShapeType(), band.F0(i), X, DC, dfdH * H, band.Mirroring() == Absorption::MirroringType::Manual);
    Calculator ls_mirr(band.Mirroring(), band.LineShapeType(), band.F0(i), X, DC, dfdH * H);
    
    if (do_cutoff) {
      // Initialize and set the cutoff values
      Calculator ls_cut = ls;
      Calculator ls_mirr_cut = ls_mirr;
      ls_cut(fu);
      ls_mirr_cut(fu);
      
      cutoff_frequency_loop(comval, ls, ls_mirr, ls_norm, ls_cut, ls_mirr_cut, ls_str,
                     derivs, LM, T, dfdH, Sz, band.Species());
      
      if (sparse_low_size) {
        cutoff_frequency_loop(sparse_low_range, ls, ls_mirr, ls_norm, ls_cut, ls_mirr_cut, ls_str,
                       derivs, LM, T, dfdH, Sz, band.Species());
      }
      
      if (sparse_upp_size) {
        cutoff_frequency_loop(sparse_upp_range, ls, ls_mirr, ls_norm, ls_cut, ls_mirr_cut, ls_str,
                       derivs, LM, T, dfdH, Sz, band.Species());
      }
      
    } else {
      frequency_loop(comval, ls, ls_mirr, ls_norm, ls_str,
                    derivs, LM, T, dfdH, Sz, band.Species());
      
      if (sparse_low_size) {
        frequency_loop(sparse_low_range, ls, ls_mirr, ls_norm, ls_str,
                      derivs, LM, T, dfdH, Sz, band.Species());
      }
      
      if (sparse_upp_size) {
        frequency_loop(sparse_upp_range, ls, ls_mirr, ls_norm, ls_str,
                      derivs, LM, T, dfdH, Sz, band.Species());
      }
    }
  }
}

void cutoff_loop_sparse_triple(ComputeData &com,
                               ComputeData &sparse_com,
                               Normalizer ls_norm,
                               const IntensityCalculator ls_str,
                               const AbsorptionLines &band,
                               const ArrayOfDerivatives &derivs,
                               const Output X,
                               const Numeric &T,
                               const Numeric &H,
                               const Numeric &sparse_lim,
                               const Numeric &DC,
                               const Index i,
                               const bool do_zeeman,
                               const Zeeman::Polarization zeeman_polarization) ARTS_NOEXCEPT {
  // Basic settings
  const bool do_nlte = ls_str.do_nlte();
  const bool do_cutoff = band.Cutoff() not_eq Absorption::CutoffType::None;
  const Numeric fu = band.CutoffFreq(i, X.D0);
  const Numeric fl = band.CutoffFreqMinus(i, X.D0);
  const Numeric fus = band.F0(i) + sparse_lim;
  const Numeric fls = band.F0(i) - sparse_lim;

  // Find sparse and dense ranges
  const auto [dense_start, dense_size,
              sparse_low_start, sparse_low_size,
              sparse_upp_start, sparse_upp_size] = quad_sparse_limited_range(fl, fu, fls, fus,
                                                                             com.f_grid,
                                                                             sparse_com.f_grid);
  if ((dense_size + sparse_low_size + sparse_upp_size) == 0) return;
  
  // Get the compute data view
  ComputeValues comval(com.F, com.dF, com.N, com.dN, com.f_grid, dense_start, dense_size, derivs, do_nlte);
  
  // Get views of the sparse data
  ComputeValues sparse_low_range(sparse_com.F, sparse_com.dF, sparse_com.N, sparse_com.dN, sparse_com.f_grid, sparse_low_start, sparse_low_size, derivs, do_nlte);
  ComputeValues sparse_upp_range(sparse_com.F, sparse_com.dF, sparse_com.N, sparse_com.dN, sparse_com.f_grid, sparse_upp_start, sparse_upp_size, derivs, do_nlte);
  
  const Index nz = do_zeeman ? band.ZeemanCount(i, zeeman_polarization) : 1;
  for (Index iz = 0; iz < nz; iz++) {
    const Numeric dfdH = do_zeeman ? band.ZeemanSplitting(i, zeeman_polarization, iz) : 0;
    const Numeric Sz = do_zeeman ? band.ZeemanStrength(i, zeeman_polarization, iz) : 1;
    const Complex LM = Complex(1 + X.G, -X.Y);
    Calculator ls(band.LineShapeType(), band.F0(i), X, DC, dfdH * H, band.Mirroring() == Absorption::MirroringType::Manual);
    Calculator ls_mirr(band.Mirroring(), band.LineShapeType(), band.F0(i), X, DC, dfdH * H);
    
    if (do_cutoff) {
      // Initialize and set the cutoff values
      Calculator ls_cut = ls;
      Calculator ls_mirr_cut = ls_mirr;
      ls_cut(fu);
      ls_mirr_cut(fu);
      
      cutoff_frequency_loop(comval, ls, ls_mirr, ls_norm, ls_cut, ls_mirr_cut, ls_str,
                            derivs, LM, T, dfdH, Sz, band.Species());
      
      if (sparse_low_size) {
        cutoff_frequency_loop(sparse_low_range, ls, ls_mirr, ls_norm, ls_cut, ls_mirr_cut, ls_str,
                              derivs, LM, T, dfdH, Sz, band.Species());
      }
      
      if (sparse_upp_size) {
        cutoff_frequency_loop(sparse_upp_range, ls, ls_mirr, ls_norm, ls_cut, ls_mirr_cut, ls_str,
                              derivs, LM, T, dfdH, Sz, band.Species());
      }
      
    } else {
      frequency_loop(comval, ls, ls_mirr, ls_norm, ls_str,
                     derivs, LM, T, dfdH, Sz, band.Species());
      
      if (sparse_low_size) {
        frequency_loop(sparse_low_range, ls, ls_mirr, ls_norm, ls_str,
                       derivs, LM, T, dfdH, Sz, band.Species());
      }
      
      if (sparse_upp_size) {
        frequency_loop(sparse_upp_range, ls, ls_mirr, ls_norm, ls_str,
                       derivs, LM, T, dfdH, Sz, band.Species());
      }
    }
  }
}

/** Cutoff considerations of the line shape
 *
 * This function takes care of setting up cutoff and line shape considerations
 * for the frequency loop function it wraps.  Internally, the cutoff is
 * calculated from the band information and a view of the correct data is sent
 * on.
 *
 * The Zeeman effect is also considered internally if applicable (or ignored
 * otherwise)
 *
 * This function is not possible to run on multiple cores.  Such parallelisms
 * must happen at a much higher level.
 *
 * @param[in,out] F The cross-section. \f$ F \f$
 * @param[in,out] dF The cross-section's derivatives. \f$ \partial F / \partial
 * x \f$
 * @param[in,out] N The cross-section ratio of the NLTE source. \f$ N \f$
 * @param[in,out] dN The cross-section ratio of the NLTE source's derivatives.
 * \f$ \partial N / \partial x \f$
 * @param[in] f_grid The frequency grid. \f$ \left[ f_0, \cdots, f_n \right] \f$
 * @param[in] band The absorption band
 * @param[in] jacobian_quantities As WSV
 * @param[in] T The atmospheric temperature
 * @param[in] do_nlte Flag for whether or not NLTE will be computed
 * @param[in] H The magnetic field magnitude
 * @param[in] do_zeeman Flag for whether this is part of some Zeeman
 * calculations
 * @param[in] zeeman_polarization The type of Zeeman polarization to consider
 * (if any)
 * @param[in] f_mean The mean frequency of the absorption band
 * @param[in] DC The Doppler broadening constant of the band
 * @param[in] i The line index
 * @param[in] X The line shape model parameters of the atmosphere
 * @param[in] LM The line mixing scaling. \f$ S_{lm} \f$
 * @param[in] ls_str The line strength calculator. \f$ S_i \f$
 * @param[in] ls_norm The normalization calculator. \f$ S_n \f$
 * @param[in] derivs A list of pre-computed derivative values and keys
 */
void cutoff_loop(ComputeData &com,
                 Normalizer ls_norm,
                 const IntensityCalculator ls_str,
                 const AbsorptionLines &band,
                 const ArrayOfDerivatives &derivs,
                 const Output X,
                 const Numeric &T,
                 const Numeric &H,
                 const Numeric &DC,
                 const Index i,
                 const bool do_zeeman,
                 const Zeeman::Polarization zeeman_polarization) ARTS_NOEXCEPT {
  // Basic settings
  const bool do_nlte = ls_str.do_nlte();
  const bool do_cutoff = band.Cutoff() not_eq Absorption::CutoffType::None;
  const Numeric fu = band.CutoffFreq(i, X.D0);
  const Numeric fl = band.CutoffFreqMinus(i, X.D0);

  // Only for the cutoff-range
  const auto [cutstart, cutsize] = limited_range(fl, fu, com.f_grid);
  if (not cutsize) return;
  
  // Get the compute data view
  ComputeValues comval(com.F, com.dF, com.N, com.dN, com.f_grid, cutstart, cutsize, derivs, do_nlte);
  
  const Index nz = do_zeeman ? band.ZeemanCount(i, zeeman_polarization) : 1;
  for (Index iz = 0; iz < nz; iz++) {
    const Numeric dfdH = do_zeeman ? band.ZeemanSplitting(i, zeeman_polarization, iz) : 0;
    const Numeric Sz = do_zeeman ? band.ZeemanStrength(i, zeeman_polarization, iz) : 1;
    const Complex LM = Complex(1 + X.G, -X.Y);
    Calculator ls(band.LineShapeType(), band.F0(i), X, DC, dfdH * H, band.Mirroring() == Absorption::MirroringType::Manual);
    Calculator ls_mirr(band.Mirroring(), band.LineShapeType(), band.F0(i), X, DC, dfdH * H);
    
    if (do_cutoff) {
      // Initialize and set the cutoff values
      Calculator ls_cut = ls;
      Calculator ls_mirr_cut = ls_mirr;
      ls_cut(fu);
      ls_mirr_cut(fu);
      
      cutoff_frequency_loop(comval, ls, ls_mirr, ls_norm, ls_cut, ls_mirr_cut, ls_str,
                            derivs, LM, T, dfdH, Sz, band.Species());
      
    } else {
      frequency_loop(comval, ls, ls_mirr, ls_norm, ls_str,
                     derivs, LM, T, dfdH, Sz, band.Species());
    }
  }
}

/** Loop all the lines of the band
 *
 * This function is not possible to run on multiple cores.  Such parallelisms
 * must happen at a much higher level.
 *
 * @param[in,out] F The cross-section. \f$ F \f$
 * @param[in,out] dF The cross-section's derivatives. \f$ \partial F / \partial
 * x \f$
 * @param[in,out] N The cross-section ratio of the NLTE source. \f$ N \f$
 * @param[in,out] dN The cross-section ratio of the NLTE source's derivatives.
 * \f$ \partial N / \partial x \f$
 * @param[in] f_grid The frequency grid. \f$ \left[ f_0, \cdots, f_n \right] \f$
 * @param[in] band The absorption band
 * @param[in] jacobian_quantities As WSV
 * @param[in] nlte A map of NLTE data
 * @param[in] vmrs The band volume mixing ratio
 * @param[in] isot_ratio The isotopic ratio of the isotopologue
 * @param[in] P The atmospheric pressure
 * @param[in] T The atmospheric temperature
 * @param[in] do_nlte Flag for whether or not NLTE will be computed
 * @param[in] H The magnetic field magnitude
 * @param[in] do_zeeman Flag for whether this is part of some Zeeman
 * calculations
 * @param[in] zeeman_polarization The type of Zeeman polarization to consider
 * (if any)
 * @param[in] f_mean The mean frequency of the absorption band
 * @param[in] QT The partition function at the temperature
 * @param[in] QT0 The partition function at the reference temperature
 * @param[in] dQTdT The derivative of the partition function at the temperature
 * wrt temperature
 */
void line_loop(ComputeData &com,
               ComputeData &sparse_com,
               const AbsorptionLines &band,
               const ArrayOfRetrievalQuantity &jacobian_quantities,
               const EnergyLevelMap &nlte,
               const Vector &vmrs,
               const ArrayOfSpeciesTag& self_tag,
               const Numeric &P,
               const Numeric &T,
               const Numeric &H,
               const Numeric &sparse_lim,
               const Numeric QT,
               const Numeric QT0,
               const Numeric dQTdT,
               const Numeric r,
               const Numeric drdSELFVMR,
               const Numeric drdT,
               const bool do_zeeman,
               const Zeeman::Polarization zeeman_polarization,
               const Options::LblSpeedup speedup_type) ARTS_NOEXCEPT {
  const Index nj = jacobian_quantities.nelem();
  const Index nl = band.NumLines();
  
  // Derivatives are allocated ahead of all loops
  ArrayOfDerivatives derivs(nj);

  // Doppler constant
  const Numeric DC = band.DopplerConstant(T);

  for (Index i = 0; i < nl; i++) {
    // Pre-compute the derivatives
    for (Index ij = 0; ij < nj; ij++) {
      const auto& deriv = jacobian_quantities[ij];
      derivs[ij].jac_pos = -1;
      derivs[ij].deriv = nullptr;
      
      if (not deriv.propmattype()) continue;
      derivs[ij].jac_pos = ij;
      derivs[ij].deriv = &deriv;

      if (deriv == Jacobian::Atm::Temperature) {
        derivs[ij].value.o = band.ShapeParameters_dT(i, T, P, vmrs);
      } else if (deriv == Jacobian::Special::ArrayOfSpeciesTagVMR) {
        if (not (deriv == self_tag)) {  // Remove if its not good
          derivs[ij].jac_pos = -1;
          derivs[ij].deriv = nullptr;
        }
      } else if (deriv.Target().needQuantumIdentity()) {
        if (deriv == Jacobian::Line::VMR) {
          derivs[ij].value.o = band.ShapeParameters_dVMR(i, T, P, deriv.QuantumIdentity());
        } else {
          auto &lt =
              derivs[ij].target = {deriv.Target().QuantumIdentity(), band, i};
          if (lt == Absorption::QuantumIdentifierLineTargetType::Line) {
            if constexpr (false) {/*skip so the rest can be a else-if block*/}
            // All line shape derivatives
            InternalDerivativesSetup(G0)
            InternalDerivativesSetup(D0)
            InternalDerivativesSetup(G2)
            InternalDerivativesSetup(D2)
            InternalDerivativesSetup(ETA)
            InternalDerivativesSetup(FVC)
            InternalDerivativesSetup(Y)
            InternalDerivativesSetup(G)
            InternalDerivativesSetup(DV)
          }
        }
      }
    }
    std::remove_if(derivs.begin(), derivs.end(), [](Derivatives& dd) { return dd.deriv == nullptr; });

    // Call cut off loop with or without sparsity
    switch (speedup_type) {
      case Options::LblSpeedup::None:
        cutoff_loop(com,
                    Normalizer(band.Normalization(), band.F0(i), T),
                    IntensityCalculator(T, QT, QT0, dQTdT, r, drdSELFVMR, drdT, nlte, band, i),
                    band, derivs, band.ShapeParameters(i, T, P, vmrs), T, H,
                    DC, i, do_zeeman, zeeman_polarization);
        break;
      case Options::LblSpeedup::QuadraticIndependent:
        cutoff_loop_sparse_triple(com, sparse_com,
                    Normalizer(band.Normalization(), band.F0(i), T),
                    IntensityCalculator(T, QT, QT0, dQTdT, r, drdSELFVMR, drdT, nlte, band, i),
                    band, derivs, band.ShapeParameters(i, T, P, vmrs), T, H, sparse_lim,
                    DC, i, do_zeeman, zeeman_polarization);
        break;
      case Options::LblSpeedup::LinearIndependent:
        cutoff_loop_sparse_linear(com, sparse_com,
                                  Normalizer(band.Normalization(), band.F0(i), T),
                                  IntensityCalculator(T, QT, QT0, dQTdT, r, drdSELFVMR, drdT, nlte, band, i),
                                  band, derivs, band.ShapeParameters(i, T, P, vmrs), T, H, sparse_lim,
                                  DC, i, do_zeeman, zeeman_polarization);
        break;
      case Options::LblSpeedup::FINAL: { /* Leave last */ }
    }
  }
}

/** Compute the line shape in its entirety
 *
 * For a single line the line shape is
 *
 * \f[ F_i = S_{z_i}  S_{n_i}  S_i  LM_i  F_i( \cdots ), \f]
 *
 * where \f$ S_{z_i} \f$ is the Zeeman scaling, \f$ S_{n_i} \f$ is the
 * normalization scaling, \f$ S_i \f$ is the line strength scaling, \f$ LM_i \f$
 * is the line mixing scaling, and \f$ F_i( \cdots )\f$ is the shape.
 *
 * and the derivatives are
 *
 * \f[
 *  \frac{\partial F_l}{\partial t} = S_{z_i} \left(
 *  \frac{\partial S_{n_i}}{\partial t}  S_i  LM_i  F_i( \cdots ) +
 *  S_{n_i}  \frac{\partial S_i}{\partial t}  LM_i  F_i( \cdots ) +
 *  S_{n_i}  S_i  \frac{\partial LM_i}{\partial t} F_i( \cdots ) +
 *  S_{n_i}  S_i  LM_i \frac{\partial F_i( \cdots )}{\partial t} \right),
 * \f]
 *
 * where \f$ t \f$ is some arbitrary variable.
 */
void compute(ComputeData &com,
             ComputeData &sparse_com,
             const AbsorptionLines &band,
             const ArrayOfRetrievalQuantity &jacobian_quantities,
             const EnergyLevelMap &nlte,
             const Vector &vmrs,
             const ArrayOfSpeciesTag& self_tag,
             const Numeric& self_vmr, const Numeric &isot_ratio, const Numeric &P, const Numeric &T, const Numeric &H,
             const Numeric &sparse_lim,
             const bool do_zeeman, const Zeeman::Polarization zeeman_polarization,
             const Options::LblSpeedup speedup_type) ARTS_NOEXCEPT {
  [[maybe_unused]] const Index nj = jacobian_quantities.nelem();
  const Index nl = band.NumLines();
  const Index nv = com.f_grid.nelem();

  // Tests that must be true while calling this function
  ARTS_ASSERT(H >= 0, "Only for positive H.  You provided: ", H)
  ARTS_ASSERT(P > 0, "Only for abs positive P.  You provided: ", P)
  ARTS_ASSERT(T > 0, "Only for abs positive T.  You provided: ", T)
  ARTS_ASSERT(band.OK(), "Band is poorly constructed.  You need to use "
                         "a detailed debugger to find out why.")
  ARTS_ASSERT(com.F.size() == nv, "F is wrong size.  Size is (", com.F.size(),
              ") but should be: (", nv, ')')
  ARTS_ASSERT(not com.do_nlte or com.N.size() == nv, "N is wrong size.  Size is (",
              com.N.size(), ") but should be (", nv, ')')
  ARTS_ASSERT(nj == 0 or (com.dF.nrows() == nv and com.dF.ncols() == nj),
              "dF is wrong size.  Size is (", com.dF.nrows(), " x ", com.dF.ncols(),
              ") but should be: (", nv, " x ", nj, ")")
  ARTS_ASSERT(nj == 0 or not com.do_nlte or (com.dN.nrows() == nv and com.dN.ncols() == nj),
              "dN is wrong size.  Size is (", com.dN.nrows(), " x ", com.dN.ncols(),
              ") but should be: (", nv, " x ", nj, ")")
  ARTS_ASSERT((sparse_lim > 0 and sparse_com.f_grid.size() > 1) or (sparse_lim == 0), "Sparse limit is either 0, or the sparse frequency grid has to have upper and lower values")

  // Early return test
  if (nv == 0 or nl == 0 or (Absorption::relaxationtype_relmat(band.Population()) and band.DoLineMixing(P))) {
    return; // No line-by-line computations required/wanted
  }
  
  const Numeric dnumdensdVMR = isot_ratio * number_density(P, T);
  line_loop(com, sparse_com, band, jacobian_quantities, nlte, vmrs, self_tag, P, T, H, sparse_lim,
            single_partition_function(T, band.Isotopologue()),
            single_partition_function(band.T0(), band.Isotopologue()),
            dsingle_partition_function_dT(T, band.Isotopologue()),
            self_vmr * dnumdensdVMR, dnumdensdVMR, self_vmr * isot_ratio * dnumber_density_dt(P, T),
            do_zeeman, zeeman_polarization, speedup_type);
}

#undef InternalDerivatives
#undef InternalDerivativesG
#undef InternalDerivativesY


Index sparse_f_grid_red(const Vector& f_grid, const Numeric& sparse_df) noexcept {
  if (f_grid.nelem())
    return f_grid.nelem() / Index(1 + std::abs(f_grid[f_grid.nelem() - 1] - f_grid[0]) / sparse_df);  
  return 0;
}

Vector linear_sparse_f_grid(const Vector& f_grid, const Numeric& sparse_df) ARTS_NOEXCEPT {
  const Index n = sparse_f_grid_red(f_grid, sparse_df);
  const Index nv = f_grid.nelem();
  
  if (nv and n) {
    std::vector<Numeric> sparse_f_grid;
    for (Index iv=0; iv<nv-n; iv+=n) {
      sparse_f_grid.emplace_back(f_grid[iv]);
      sparse_f_grid.emplace_back(f_grid[iv + n]);
    }
    
    const Numeric f0 = sparse_f_grid.back();
    if (f0 not_eq f_grid[nv-1]) {
      sparse_f_grid.emplace_back(f0);
      sparse_f_grid.emplace_back(f_grid[nv-1]);
    }
    
    return sparse_f_grid;
  }
  return Vector(0);
}

bool good_linear_sparse_f_grid(const Vector& f_grid_dense, const Vector& f_grid_sparse) noexcept {
  const Index nf_sparse=f_grid_sparse.nelem();
  const Index nf_dense=f_grid_dense.nelem();

  if (nf_sparse == 1)
    return false;

  if(nf_sparse and nf_dense)
    return f_grid_dense[0] >= f_grid_sparse[0] and f_grid_dense[nf_dense-1] <= f_grid_sparse[nf_sparse-1];
  
  return true;
}

Vector triple_sparse_f_grid(const Vector& f_grid, const Numeric& sparse_df) noexcept {
  const Index n = sparse_f_grid_red(f_grid, sparse_df);
  const Index nv = f_grid.nelem();
  
  if (nv and n > 2) {
    std::vector<Numeric> sparse_f_grid;
    for (Index iv=0; iv<nv-n; iv+=n) {
      sparse_f_grid.emplace_back(f_grid[iv]);
      sparse_f_grid.emplace_back(f_grid[iv] + 0.5 * (f_grid[iv+n] - f_grid[iv]));
      sparse_f_grid.emplace_back(f_grid[iv+n]);
    }
    
    const Numeric f0 = sparse_f_grid.back();
    if (f0 not_eq f_grid[nv-1]) {
      sparse_f_grid.emplace_back(f0);
      sparse_f_grid.emplace_back(f0 + 0.5 * (f_grid[nv-1] - f0));
      sparse_f_grid.emplace_back(f_grid[nv-1]);
    }
    
    return sparse_f_grid;
  }
  return Vector(0);
}

void ComputeData::interp_add_even(const ComputeData& sparse) ARTS_NOEXCEPT {
  const Index nv = f_grid.nelem();
  const Index sparse_nv = sparse.f_grid.nelem();
  const Index nj = dF.ncols();
  
  ARTS_ASSERT(do_nlte == sparse.do_nlte, "Must have the same NLTE status")
  ARTS_ASSERT(sparse_nv > 1, "Must have at least two sparse grid-points")
  ARTS_ASSERT(nv == 0 or (f_grid[0] == sparse.f_grid[0] and f_grid[nv - 1] >= sparse.f_grid[sparse_nv - 1]),
              "If there are any dense frequency points, then the sparse frequency points must fully cover them")
  ARTS_ASSERT(not (sparse_nv % 2), "Must be multiple of to")
  
  Index sparse_iv=0;
  Numeric f1 = sparse.f_grid[sparse_iv + 1];
  Numeric f0 = sparse.f_grid[sparse_iv];
  Numeric inv = 1.0 / (f1 - f0);
  for (Index iv=0; iv<nv; iv++) {
    if (sparse_iv < (sparse_nv - 2) and f1 == f_grid[iv]) {
      sparse_iv += 2;
      f1 = sparse.f_grid[sparse_iv + 1];
      f0 = sparse.f_grid[sparse_iv];
      inv = 1.0 / (f1 - f0);
    }
    
    const Numeric xm0 = f_grid[iv] - f0;
    const Numeric xm1 = f_grid[iv] - f1;
    const Numeric l0 = - xm1 * inv;
    const Numeric l1 = xm0 * inv;
    
    F[iv] += l0 * sparse.F[sparse_iv] + l1 * sparse.F[sparse_iv + 1];
    for (Index ij=0; ij<nj; ij++) {
      dF(iv, ij) += l0 * sparse.dF(sparse_iv, ij) + l1 * sparse.dF(sparse_iv + 1, ij);
    }
    if (do_nlte) {
      N[iv] += l0 * sparse.N[sparse_iv] + l1 * sparse.N[sparse_iv + 1];
      for (Index ij=0; ij<nj; ij++) {
        dN(iv, ij) += l0 * sparse.dN(sparse_iv, ij) + l1 * sparse.dN(sparse_iv + 1, ij);
      }
    }
  }
}

void ComputeData::interp_add_triplequad(const ComputeData& sparse) ARTS_NOEXCEPT {
  const Index nv = f_grid.nelem();
  const Index sparse_nv = sparse.f_grid.nelem();
  const Index nj = dF.ncols();
  
  ARTS_ASSERT(do_nlte == sparse.do_nlte, "Must have the same NLTE status")
  ARTS_ASSERT(sparse_nv > 2, "Must have at least three sparse grid-points")
  ARTS_ASSERT(nv == 0 or (f_grid[0] == sparse.f_grid[0] and f_grid[nv - 1] >= sparse.f_grid[sparse_nv - 1]),
              "If there are any dense frequency points, then the sparse frequency points must fully cover them")
  ARTS_ASSERT(not (sparse_nv % 3), "Must be multiple of three")
  
  Index sparse_iv=0;
  Numeric f2 = sparse.f_grid[sparse_iv + 2];
  Numeric f1 = sparse.f_grid[sparse_iv + 1];
  Numeric f0 = sparse.f_grid[sparse_iv];
  Numeric inv = 1.0 / Constant::pow2(f1 - f0);
  for (Index iv = 0; iv < nv; iv++) {
    if (sparse_iv < (sparse_nv - 3) and f2 == f_grid[iv]) {
      sparse_iv += 3;
      f2 = sparse.f_grid[sparse_iv + 2];
      f1 = sparse.f_grid[sparse_iv + 1];
      f0 = sparse.f_grid[sparse_iv];
      inv = 1.0 / Constant::pow2(f1 - f0);
    }
    ARTS_ASSERT (f_grid[iv] >= f0 and (f_grid[iv] < f2 or (f2 == f_grid[iv] and sparse_iv == sparse_nv - 3)),
                 "Out of range frequency grid.  Must be caught earlier.\n"
                 "The sparse range is from: ", f0, " to ", f2, " with ", f1, " as the half-way grid point.\n"
                 "The dense frequency is ", f_grid[iv], " and the indices are: sparse_iv=", sparse_iv, "; iv=", iv)
    
    const Numeric xm0 = f_grid[iv] - f0;
    const Numeric xm1 = f_grid[iv] - f1;
    const Numeric xm2 = f_grid[iv] - f2;
    const Numeric l0 = 0.5 * xm1 * xm2 * inv;  // --
    const Numeric l1 = - xm0 * xm2 * inv;      // +-
    const Numeric l2 = 0.5 * xm0 * xm1 * inv;  // ++
    
    F[iv] += l0 * sparse.F[sparse_iv] + l1 * sparse.F[sparse_iv + 1] + l2 * sparse.F[sparse_iv + 2];
    for (Index ij=0; ij<nj; ij++) {
      dF(iv, ij) += l0 * sparse.dF(sparse_iv, ij) + l1 * sparse.dF(sparse_iv + 1, ij) + l2 * sparse.dF(sparse_iv + 2, ij);
    }
    if (do_nlte) {
      N[iv] += l0 * sparse.N[sparse_iv] + l1 * sparse.N[sparse_iv + 1] + l2 * sparse.N[sparse_iv + 2];
      for (Index ij=0; ij<nj; ij++) {
        dN(iv, ij) += l0 * sparse.dN(sparse_iv, ij) + l1 * sparse.dN(sparse_iv + 1, ij) + l2 * sparse.dN(sparse_iv + 2, ij);
      }
    }
  }
}
} // namespace LineShape
