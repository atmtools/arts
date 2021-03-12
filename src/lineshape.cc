#include <memory_resource>

#include "lineshape.h"

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
Doppler::Doppler(Numeric F0_noshift, Numeric DC, Numeric dZ) noexcept
    : mF0(F0_noshift + dZ), invGD(sqrt_ln_2 / std::abs(DC * mF0)) {}

Complex Doppler::dFdT(const Output &, Numeric T) const noexcept {
  return F * (2 * pow2(x) - 1) / (2 * T);
}

Complex Doppler::dFdf() const noexcept { return -2 * invGD * F * x; }

Complex Doppler::dFdF0() const noexcept {
  return F * (2 * x * (invGD * mF0 + x) - 1) / mF0;
}

Complex Doppler::dFdH(Numeric dZ) const noexcept {
  return dZ * (F * (2 * x * (invGD * mF0 + x) - 1) / mF0);
}

Complex Doppler::operator()(Numeric f) noexcept {
  x = (f - mF0) * invGD;
  F = invGD * inv_sqrt_pi * std::exp(-pow2(x));
  return F;
}

Voigt::Voigt(Numeric F0_noshift, const Output &ls, Numeric DC,
             Numeric dZ) noexcept
    : mF0(F0_noshift + dZ + ls.D0 + ls.DV), invGD(sqrt_ln_2 / std::abs(DC * mF0)),
      z(invGD * Complex(-mF0, ls.G0)) {}

Complex Voigt::dFdf() const noexcept { return dF; }

Complex Voigt::dFdF0() const noexcept { return -dF; }

Complex Voigt::dFdDV(Numeric d) const noexcept { return -d * dF; }

Complex Voigt::dFdD0(Numeric d) const noexcept { return -d * dF; }

Complex Voigt::dFdG0(Numeric d) const noexcept { return Complex(0, d) * dF; }

Complex Voigt::dFdH(Numeric dZ) const noexcept { return -dZ * dF; }

Complex Voigt::dFdVMR(const Output &d) const noexcept {
  return Complex(-d.D0 - d.DV, d.G0) * dF;
}

Complex Voigt::dFdT(const Output &d, Numeric T) const noexcept {
  return -(F * invGD + dF * z) * (2 * T * (d.D0 + d.DV) + mF0) /
             (2 * T * invGD * mF0) +
         Complex(-d.D0 - d.DV, d.G0) * dF;
}

Complex Voigt::operator()(Numeric f) noexcept {
  reinterpret_cast<Numeric(&)[2]>(z)[0] = invGD * (f - mF0);
  F = inv_sqrt_pi * invGD * Faddeeva::w(z);
  dF = 2 * invGD * (Complex(0, inv_pi * invGD) - z * F);
  return F;
}

SpeedDependentVoigt::SpeedDependentVoigt(Numeric F0_noshift, const Output &ls,
                                         Numeric GD_div_F0, Numeric dZ) noexcept
    : mF0(F0_noshift + dZ + ls.D0 - 1.5 * ls.D2),
      invGD(sqrt_ln_2 / std::abs(GD_div_F0 * mF0)), invc2(1.0 / Complex(ls.G2, ls.D2)),
      dx(Complex(ls.G0 - 1.5 * ls.G2, mF0)), x(dx * invc2),
      sqrty(invc2 / (2 * invGD)), calcs(init(Complex(ls.G2, ls.D2))) {
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
  reinterpret_cast<Numeric(&)[2]>(dx)[1] = mF0 - f;
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
  else if (abs_squared(x) <= 9e-16 * abs_squared(sqrty * sqrty))
    return CalcType::LowXandHighY;
  else if ((abs_squared(sqrty * sqrty) <= 1.e-30 * abs_squared(x)) and
           abs_squared(std::sqrt(x)) <= 16.e6)
    return CalcType::LowYandLowX; // Weird case, untested
  else if ((abs_squared(sqrty * sqrty) <= 1.e-30 * abs_squared(x)))
    return CalcType::LowYandHighX;
  else
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
      invGD(sqrt_ln_2 / std::abs(GD_div_F0 * mF0)),
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
  reinterpret_cast<Numeric(&)[2]>(deltax)[1] = mF0 - f;
  x = deltax / ((1 - ETA) * Complex(G2, D2));
  sqrtxy = std::sqrt(x + sqrty * sqrty);
  update_calcs();
  calc();
  return F;
}

HartmannTran::CalcType HartmannTran::init(const Complex c2t) const noexcept {
  if (abs_squared(c2t) == 0)
    return CalcType::Noc2tHighZ; // nb. Value of high/low changes elsewhere
  else if (abs_squared(x) <= 9e-16 * abs_squared(sqrty * sqrty))
    return CalcType::LowXandHighY;
  else if ((abs_squared(sqrty * sqrty) <= 1.e-30 * abs_squared(x)) and
           abs_squared(std::sqrt(x)) <= 16.e6)
    return CalcType::LowYandLowX; // Weird case, untested
  else if ((abs_squared(sqrty * sqrty) <= 1.e-30 * abs_squared(x)))
    return CalcType::LowYandHighX;
  else
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

LocalThermodynamicEquilibrium::LocalThermodynamicEquilibrium(
    Numeric I0, Numeric T0, Numeric T, Numeric F0, Numeric E0, Numeric QT,
    Numeric QT0, Numeric dQTdT, Numeric r) noexcept
    : LocalThermodynamicEquilibrium(
          I0, r, QT0, QT, dQTdT, boltzman_ratio(T, T0, E0),
          dboltzman_ratio_dT_div_boltzmann_ratio(T, E0),
          stimulated_relative_emission(F0, T0, T),
          dstimulated_relative_emission_dT(F0, T0, T),
          dstimulated_relative_emission_dF0(F0, T0, T)) {}

FullNonLocalThermodynamicEquilibrium::FullNonLocalThermodynamicEquilibrium(
    Numeric F0, Numeric A21, Numeric T, Numeric g1, Numeric g2, Numeric r1,
    Numeric r2, Numeric r) noexcept
    : dNdTval(
          ((c1 * F0) * r2 * A21) * Conversion::hz2joule(F0) *
          (std::exp(Conversion::hz2joule(F0) / Conversion::kelvin2joule(T))) /
          ((c0 * F0 * F0 * F0) * Conversion::kelvin2joule(T) * T)),

      dSdF0val(c1 * (r1 * (g2 / g1) - r2) * (A21 / (c0 * F0 * F0 * F0)) -
               3.0 *
                   ((c1 * F0) * (r1 * (g2 / g1) - r2) *
                    (A21 / (c0 * F0 * F0 * F0))) /
                   F0),
      dNdF0val(
          ((c1 * F0) * r2 * A21) *
              (Constant::h *
                   (std::exp(Conversion::hz2joule(F0) /
                             Conversion::kelvin2joule(T))) /
                   ((c0 * F0 * F0 * F0) * Conversion::kelvin2joule(T)) -
               3.0 *
                   ((c0 * F0 * F0 * F0) /
                    std::expm1(Conversion::hz2joule(F0) /
                               Conversion::kelvin2joule(T))) /
                   F0) +
          c1 * r2 * A21 /
              ((c0 * F0 * F0 * F0) / std::expm1(Conversion::hz2joule(F0) /
                                                Conversion::kelvin2joule(T)))),

      dSdr2(-(c1 * F0) * A21 / (c0 * F0 * F0 * F0)),
      dNdr2((c1 * F0) * A21 /
            ((c0 * F0 * F0 * F0) / std::expm1(Conversion::hz2joule(F0) /
                                              Conversion::kelvin2joule(T)))),
      dSdr1((c1 * F0) * (g2 / g1) * A21 / (c0 * F0 * F0 * F0)),

      S(r * ((c1 * F0) * (r1 * (g2 / g1) - r2) * (A21 / (c0 * F0 * F0 * F0)))),
      N(r *
        (((c1 * F0) * r2 * A21) /
             ((c0 * F0 * F0 * F0) / std::expm1(Conversion::hz2joule(F0) /
                                               Conversion::kelvin2joule(T))) -
         ((c1 * F0) * (r1 * (g2 / g1) - r2) * (A21 / (c0 * F0 * F0 * F0))))) {}

VibrationalTemperaturesNonLocalThermodynamicEquilibrium::
    VibrationalTemperaturesNonLocalThermodynamicEquilibrium(
        Numeric I0, Numeric T0, Numeric T, Numeric Tl, Numeric Tu, Numeric F0,
        Numeric E0, Numeric Evl, Numeric Evu, Numeric QT, Numeric QT0,
        Numeric dQTdT, Numeric r) noexcept
    :

      dSdI0val(r * boltzman_ratio(T, T0, E0) *
               stimulated_relative_emission(stimulated_emission(T, F0),
                                            stimulated_emission(T0, F0)) *
               absorption_nlte_ratio(stimulated_emission(T, F0),
                                     boltzman_ratio(Tu, T, Evu),
                                     boltzman_ratio(Tl, T, Evl)) *
               QT0 / QT),
      dNdI0val(r * boltzman_ratio(T, T0, E0) *
               stimulated_relative_emission(stimulated_emission(T, F0),
                                            stimulated_emission(T0, F0)) *
               boltzman_ratio(Tu, T, Evu) * QT0 / QT),

      dSdTval(I0 * r * boltzman_ratio(T, T0, E0) *
                  dstimulated_relative_emission_dT(stimulated_emission(T, F0),
                                                   stimulated_emission(T0, F0),
                                                   F0, T) *
                  absorption_nlte_ratio(stimulated_emission(T, F0),
                                        boltzman_ratio(Tu, T, Evu),
                                        boltzman_ratio(Tl, T, Evl)) *
                  QT0 / QT +
              I0 * r * dboltzman_ratio_dT(boltzman_ratio(T, T0, E0), T, E0) *
                  stimulated_relative_emission(stimulated_emission(T, F0),
                                               stimulated_emission(T0, F0)) *
                  absorption_nlte_ratio(stimulated_emission(T, F0),
                                        boltzman_ratio(Tu, T, Evu),
                                        boltzman_ratio(Tl, T, Evl)) *
                  QT0 / QT +
              I0 * r * boltzman_ratio(T, T0, E0) *
                  stimulated_relative_emission(stimulated_emission(T, F0),
                                               stimulated_emission(T0, F0)) *
                  dabsorption_nlte_rate_dT(stimulated_emission(T, F0), T, F0,
                                           Evl, Evu, boltzman_ratio(Tu, T, Evu),
                                           boltzman_ratio(Tl, T, Evl)) *
                  QT0 / QT -
              I0 * dSdI0val * dQTdT / QT),
      dNdTval(I0 * r * boltzman_ratio(T, T0, E0) *
                  dstimulated_relative_emission_dT(stimulated_emission(T, F0),
                                                   stimulated_emission(T0, F0),
                                                   F0, T) *
                  boltzman_ratio(Tu, T, Evu) * QT0 / QT +
              I0 * r * dboltzman_ratio_dT(boltzman_ratio(T, T0, E0), T, E0) *
                  stimulated_relative_emission(stimulated_emission(T, F0),
                                               stimulated_emission(T0, F0)) *
                  boltzman_ratio(Tu, T, Evu) * QT0 / QT +
              I0 * r * boltzman_ratio(T, T0, E0) *
                  stimulated_relative_emission(stimulated_emission(T, F0),
                                               stimulated_emission(T0, F0)) *
                  dboltzman_ratio_dT(boltzman_ratio(Tu, T, Evu), T, Evu) * QT0 /
                  QT -
              I0 * dSdI0val * dQTdT / QT),

      dSdF0val(I0 * r * boltzman_ratio(T, T0, E0) *
                   dstimulated_relative_emission_dF0(
                       stimulated_emission(T, F0), stimulated_emission(T0, F0),
                       T, T0) *
                   absorption_nlte_ratio(stimulated_emission(T, F0),
                                         boltzman_ratio(Tu, T, Evu),
                                         boltzman_ratio(Tl, T, Evl)) *
                   QT0 / QT +
               I0 * r * boltzman_ratio(T, T0, E0) *
                   stimulated_relative_emission(stimulated_emission(T, F0),
                                                stimulated_emission(T0, F0)) *
                   dabsorption_nlte_rate_dF0(stimulated_emission(T, F0), T,
                                             boltzman_ratio(Tu, T, Evu),
                                             boltzman_ratio(Tl, T, Evl)) *
                   QT0 / QT),
      dNdF0val(I0 * r * boltzman_ratio(T, T0, E0) *
               dstimulated_relative_emission_dF0(stimulated_emission(T, F0),
                                                 stimulated_emission(T0, F0), T,
                                                 T0) *
               boltzman_ratio(Tu, T, Evu) * QT0 / QT),

      dSdTl(I0 * r * boltzman_ratio(T, T0, E0) *
            stimulated_relative_emission(stimulated_emission(T, F0),
                                         stimulated_emission(T0, F0)) *
            dabsorption_nlte_rate_dTl(stimulated_emission(T, F0), T, Tl, Evl,
                                      boltzman_ratio(Tl, T, Evl)) *
            QT0 / QT),
      dSdTu(I0 * r * boltzman_ratio(T, T0, E0) *
            stimulated_relative_emission(stimulated_emission(T, F0),
                                         stimulated_emission(T0, F0)) *
            dabsorption_nlte_rate_dTu(stimulated_emission(T, F0), T, Tu, Evu,
                                      boltzman_ratio(Tu, T, Evu)) *
            QT0 / QT),
      dNdTu(I0 * r * boltzman_ratio(T, T0, E0) *
            stimulated_relative_emission(stimulated_emission(T, F0),
                                         stimulated_emission(T0, F0)) *
            dboltzman_ratio_dT(boltzman_ratio(Tu, T, Evu), Tu, Evu) * QT0 / QT),

      S(I0 * dSdI0val), N(I0 * dNdI0val - I0 * dSdI0val) {}

Calculator line_shape_selection(const Type type, const Numeric F0,
                                const Output &X, const Numeric DC,
                                const Numeric DZ) {
  switch (type) {
  case Type::DP:
    return Doppler(F0, DC, DZ);
  case Type::LP:
    return Lorentz(F0, X);
  case Type::VP:
    return Voigt(F0, X, DC, DZ);
  case Type::SDVP:
    return SpeedDependentVoigt(F0, X, DC, DZ);
  case Type::HTP:
    return HartmannTran(F0, X, DC, DZ);
  case Type::FINAL: { /*leave last*/
  }
  }

  return Noshape{};
}

Calculator mirror_line_shape_selection(const Absorption::MirroringType mirror,
                                       const Type type, const Numeric F0,
                                       const Output &X, const Numeric DC,
                                       const Numeric DZ) {
  switch (mirror) {
  case Absorption::MirroringType::Lorentz:
    return Lorentz(-F0, mirroredOutput(X));
  case Absorption::MirroringType::SameAsLineShape:
    return line_shape_selection(type, -F0, mirroredOutput(X), -DC, -DZ);
  case Absorption::MirroringType::Manual:
    return Noshape{};
  case Absorption::MirroringType::None:
    return Noshape{};
  case Absorption::MirroringType::FINAL: { /*leave last*/
  }
  }

  return Noshape{};
}

Normalizer normalizer_selection(const Absorption::NormalizationType type,
                                const Numeric F0, const Numeric T) {
  switch (type) {
  case Absorption::NormalizationType::None:
    return Nonorm{};
  case Absorption::NormalizationType::RQ:
    return RosenkranzQuadratic(F0, T);
  case Absorption::NormalizationType::VVH:
    return VanVleckHuber(F0, T);
  case Absorption::NormalizationType::VVW:
    return VanVleckWeisskopf(F0);
  case Absorption::NormalizationType::FINAL: { /*leave last*/
  }
  }

  return Nonorm{};
}

IntensityCalculator linestrength_selection(const Numeric T, const Numeric QT,
                                           const Numeric QT0,
                                           const Numeric dQTdT, const Numeric r,
                                           const EnergyLevelMap &nlte,
                                           const Absorption::Lines &band,
                                           const Index line_index) {
  const auto &line = band.Line(line_index);
  switch (band.Population()) {
  case Absorption::PopulationType::ByHITRANFullRelmat:
  case Absorption::PopulationType::ByHITRANRosenkranzRelmat:
  case Absorption::PopulationType::ByMakarovFullRelmat:
  case Absorption::PopulationType::LTE:
    return LocalThermodynamicEquilibrium(line.I0(), band.T0(), T, line.F0(),
                                         line.E0(), QT, QT0, dQTdT, r);
  case Absorption::PopulationType::NLTE: {
    const auto [r_low, r_upp] = nlte.get_ratio_params(band, line_index);
    return FullNonLocalThermodynamicEquilibrium(
        line.F0(), line.A(), T, line.g_low(), line.g_upp(), r_low, r_upp, r);
  }
  case Absorption::PopulationType::VibTemps: {
    const auto [E_low, E_upp, T_low, T_upp] =
        nlte.get_vibtemp_params(band, line_index, T);
    return VibrationalTemperaturesNonLocalThermodynamicEquilibrium(
        line.I0(), band.T0(), T, T_low, T_upp, line.F0(), line.E0(), E_low,
        E_upp, QT, QT0, dQTdT, r);
  }
  case Absorption::PopulationType::FINAL: { /*leave last*/
  }
  }

  return Nostrength{};
}

#define InternalDerivativesImpl(X, Y)                                          \
  else if (deriv == Jacobian::Line::Shape##X##Y) {                             \
    const Numeric d = line_derivs[iline_derivs++];                             \
    const Complex dFm = std::visit(                                            \
        [d](auto &&LS) { return std::conj(LS.dFd##X(d)); },                    \
        ls_mirr);                                                              \
    const Complex dFls =                                                       \
        std::visit([d](auto &&LS) { return LS.dFd##X(d); },                    \
                   ls) +                                                       \
        dFm;                                                                   \
    dF(iv, ij) += S * LM * dFls;                                               \
  }

#define InternalDerivatives(X)                                                 \
  InternalDerivativesImpl(X, X0) InternalDerivativesImpl(X, X1)                \
      InternalDerivativesImpl(X, X2) InternalDerivativesImpl(X, X3)


#define InternalDerivativesSetupImpl(X, Y)                                     \
  else if (deriv == Jacobian::Line::Shape##X##Y) {                             \
      line_derivs.emplace_back(band.Line(i).LineShape().d##X##_d##Y(           \
          T, band.T0(), P, deriv.Target().Position(), vmrs));                  \
  }
  
#define InternalDerivativesSetup(X)                                            \
  InternalDerivativesSetupImpl(X, X0) InternalDerivativesSetupImpl(X, X1)      \
  InternalDerivativesSetupImpl(X, X2) InternalDerivativesSetupImpl(X, X3)

#define InternalDerivativesGImpl(X)                                            \
  else if (deriv == Jacobian::Line::ShapeG##X) {                               \
    const Numeric dLM = line_derivs[iline_derivs++];                           \
    dF(iv, ij) += dLM * S * Fls;                                               \
  }

#define InternalDerivativesG                                                   \
  InternalDerivativesGImpl(X0) InternalDerivativesGImpl(X1)                    \
      InternalDerivativesGImpl(X2) InternalDerivativesGImpl(X3)

#define InternalDerivativesYImpl(X)                                            \
  else if (deriv == Jacobian::Line::ShapeY##X) {                               \
    const Complex dLM = Complex(0, -line_derivs[iline_derivs++]);             \
    dF(iv, ij) += dLM * S * Fls;                                               \
  }

#define InternalDerivativesY                                                   \
  InternalDerivativesYImpl(X0) InternalDerivativesYImpl(X1)                    \
      InternalDerivativesYImpl(X2) InternalDerivativesYImpl(X3)

inline
void frequency_loop(ComplexVector &F, ComplexMatrix &dF,
                    ComplexVector &N [[maybe_unused]],
                    ComplexMatrix &dN [[maybe_unused]],
                    const std::pmr::vector<Numeric> &line_derivs,
                    const std::pmr::vector<Output> &lsmp_derivs, Calculator &ls,
                    Calculator &ls_mirr, Normalizer &ls_norm,
                    IntensityCalculator &ls_str, const Vector &f_grid,
                    const AbsorptionLines &band,
                    const ArrayOfRetrievalQuantity &jacobian_quantities,
                    const Complex &LM, const Numeric &Si,
                    const Numeric &DNi, const Numeric &Sz,
                    const Numeric &T, const Numeric &fu, const Numeric &fl,
                    const Numeric &dfdH, const Index &i,
                    const Index &nj, const Index &nv, bool do_nlte) ARTS_NOEXCEPT {
  for (Index iv = 0; iv < nv; iv++) {
    const Numeric f = f_grid[iv];
    if (f > fu or f < fl) {
      // NOTE: f > fu means that fu might be computed twice for cutoff
      continue;
    }

    const Numeric Sn = std::visit([f](auto &&LSN) { return LSN(f); }, ls_norm);
    const Numeric S = Sz * Sn * Si;
    [[maybe_unused]] const Numeric DS = Sz * Sn * DNi;
    const Complex Fm =
        std::visit([f](auto &&LS) { return std::conj(LS(f)); }, ls_mirr);
    const Complex Fls = std::visit([f](auto &&LS) { return LS(f); }, ls) + Fm;
    F[iv] += S * LM * Fls;
    if (do_nlte) {
      N[iv] += DS * LM * Fls;
    }

    Index ilsmp_derivs=0, iline_derivs=0;
    for (Index ij = 0; ij < nj; ij++) {
      if (not propmattype_index(jacobian_quantities, ij)) {
        continue;
      }
      
      const auto &deriv = jacobian_quantities[ij];

      if (deriv == Jacobian::Atm::Temperature) {
        const auto& dXdT = lsmp_derivs[ilsmp_derivs++];
        const Numeric dSn = std::visit([T, f](auto &&LSN) { return LSN.dNdT(T, f); }, ls_norm);
        const Numeric dSi = std::visit([](auto &&LSN) { return LSN.dSdT(); }, ls_str);
        const Complex dLM(dXdT.G, -dXdT.Y);
        const Complex dFm = std::visit([dXdT, T](auto &&LS) { return std::conj(LS.dFdT(mirroredOutput(dXdT), T)); }, ls_mirr);
        const Complex dFls = std::visit([dXdT, T](auto &&LS) { return LS.dFdT(dXdT, T); }, ls) + dFm;
        dF(iv, ij) += S * (LM * dFls + dLM * Fls) + Sz * (dSn * Si + Sn * dSi) * LM * Fls;
        if (do_nlte) {
          const Numeric dDSi = std::visit([](auto &&LSN) { return LSN.dNdT(); }, ls_str);
          dN(iv, ij) += DS * (LM * dFls + dLM * Fls) + Sz * (dSn * Si + Sn * dDSi) * LM * Fls;
        }
      } else if (is_wind_parameter(deriv)) {
        const Complex dFm =
            std::visit([](auto &&LS) { return std::conj(LS.dFdf()); }, ls_mirr);
        const Complex dFls =
            std::visit([](auto &&LS) { return LS.dFdf(); }, ls) + dFm;
        dF(iv, ij) += S * LM * dFls;
      } else if (is_magnetic_parameter(deriv)) {
        const Complex dFm = std::visit(
            [dfdH](auto &&LS) { return std::conj(LS.dFdH(-dfdH)); }, ls_mirr);
        const Complex dFls =
            std::visit([dfdH](auto &&LS) { return LS.dFdH(dfdH); }, ls) + dFm;
        dF(iv, ij) += S * LM * dFls;
      } else if (deriv.Target().needQuantumIdentity()) {
        if (deriv == Jacobian::Line::VMR) {
          const auto& dXdVMR = lsmp_derivs[ilsmp_derivs++];
          const Complex dFm = std::visit(
              [dXdVMR](auto &&LS) {
                return std::conj(LS.dFdVMR(mirroredOutput(dXdVMR)));
              },
              ls_mirr);
          const Complex dLM(dXdVMR.G, -dXdVMR.Y);
          const Complex dFls =
              std::visit([dXdVMR](
                             auto &&LS) { return LS.dFdVMR(dXdVMR); },
                         ls) +
              dFm;
          dF(iv, ij) += S * LM * dFls + dLM * S * Fls;
        } else {
          const Absorption::QuantumIdentifierLineTarget lt(
              deriv.Target().QuantumIdentity(), band, i);
          if (lt == Absorption::QuantumIdentifierLineTargetType::Line) {
            if (deriv == Jacobian::Line::Center) {
              const Numeric d =
                  std::visit([](auto &&LS) { return LS.dNdF0(); }, ls_norm) * Si +
                  std::visit([](auto &&LS) { return LS.dSdF0(); }, ls_str) *  Sn;
              const Complex dFm = std::visit(
                  [](auto &&LS) { return std::conj(LS.dFdF0()); }, ls_mirr);
              const Complex dFls =
                  std::visit([](auto &&LS) { return LS.dFdF0(); }, ls) + dFm;
              const Numeric dS = Sz * d;
              dF(iv, ij) += S * LM * dFls + dS * LM * Fls;
            } else if (deriv == Jacobian::Line::Strength) {
              const Numeric dS = Sz * Sn * line_derivs[iline_derivs++];
              dF(iv, ij) += dS * LM * Fls;
            }
            // All line shape derivatives
            InternalDerivatives(G0) InternalDerivatives(D0)
                InternalDerivatives(G2) InternalDerivatives(D2)
                    InternalDerivatives(ETA) InternalDerivatives(FVC)
                        InternalDerivativesY InternalDerivativesG
                            InternalDerivatives(DV)
          } else if (lt == Absorption::QuantumIdentifierLineTargetType::Level) {
            if (deriv == Jacobian::Line::NLTE) {
              if (lt.upper) {
                const Numeric dS = Sz * Sn * line_derivs[iline_derivs++];;
                dF(iv, ij) += dS * LM * Fls;
              }

              if (lt.lower) {
                const Numeric dS = Sz * Sn * line_derivs[iline_derivs++];;
                dF(iv, ij) += dS * LM * Fls;
              }
            }
          }
        }
      }
    }
  }
}

inline
void zeeman_loop(ComplexVector &F, ComplexMatrix &dF, ComplexVector &N,
                 ComplexMatrix &dN,
                 const std::pmr::vector<Numeric> &line_derivs,
                 const std::pmr::vector<Output> &lsmp_derivs, Normalizer &ls_norm,
                 IntensityCalculator &ls_str, const Vector &f_grid,
                 const AbsorptionLines &band,
                 const ArrayOfRetrievalQuantity &jacobian_quantities,
                 const Output &X, const Complex &LM,
                 const Numeric &Si, const Numeric &DSi, const Numeric &f_mean,
                 const Numeric &T, const Numeric &DC,
                 const Numeric &H, const Index &i, const Index &nz,
                 const Index &nj, const Index &nv,
                 bool do_nlte, bool do_zeeman,
                 const Zeeman::Polarization &zeeman_polarization
                 [[maybe_unused]]) ARTS_NOEXCEPT {
  const Numeric fu = band.CutoffFreq(i), fl = band.CutoffFreqMinus(i, f_mean);

  for (Index iz = 0; iz < nz; iz++) {
    const Numeric dfdH =
        do_zeeman ? band.ZeemanSplitting(i, zeeman_polarization, iz) : 0;
    const Numeric Sz =
        do_zeeman ? band.ZeemanStrength(i, zeeman_polarization, iz) : 1;

    // Line shape
    Calculator ls =
        line_shape_selection(band.LineShapeType(), band.F0(i), X, DC, dfdH * H);

    // Line mirroring around F0 = 0
    Calculator ls_mirr = mirror_line_shape_selection(
        band.Mirroring(), band.LineShapeType(), band.F0(i), X, DC, dfdH * H);

    frequency_loop(F, dF, N, dN, line_derivs, lsmp_derivs, ls, ls_mirr,
                   ls_norm, ls_str, f_grid, band, jacobian_quantities,
                   LM, Si, DSi, Sz, T, fu, fl, dfdH, i,
                   nj, nv, do_nlte);
  }

  // Deal with cutoff if necessary (on stack unless there's a derivative)
  // Added limitation of function: f_grid must be sorted from low-to-high
  if (nv and f_grid[0] not_eq fu and
      fu not_eq std::numeric_limits<Numeric>::max()) {
    Numeric f = fu;
    Complex Fc(0, 0);
    Complex Nc(0, 0);
    Vector f_grid_cut(&f, Range(0, 1));
    ComplexVector F_cut(&Fc, Range(0, 1));
    ComplexVector N_cut(&Nc, Range(0, 1));
    ComplexMatrix dF_cut(1, nj, 0);
    ComplexMatrix dN_cut(1, nj, 0);

    // Compute the cutoff
    for (Index iz = 0; iz < nz; iz++) {
      const Numeric dfdH =
          do_zeeman ? band.ZeemanSplitting(i, zeeman_polarization, iz) : 0;
      const Numeric Sz =
          do_zeeman ? band.ZeemanStrength(i, zeeman_polarization, iz) : 1;

      // Line shape
      Calculator ls = line_shape_selection(band.LineShapeType(), band.F0(i), X,
                                           DC, dfdH * H);

      // Line mirroring around F0 = 0
      Calculator ls_mirr = mirror_line_shape_selection(
          band.Mirroring(), band.LineShapeType(), band.F0(i), X, DC, dfdH * H);

      // Call for just 1 freq...
      frequency_loop(F_cut, dF_cut, N_cut, dN_cut, line_derivs,
                     lsmp_derivs, ls, ls_mirr, ls_norm, ls_str,
                     f_grid_cut, band, jacobian_quantities, LM,
                     Si, DSi, Sz, T, fu, fl, dfdH, i, nj, nv, do_nlte);
    }

    // Remove the cutoff
    for (Index iv = 0; iv < nv; iv++) {
      if (f_grid[iv] > fu or f_grid[iv] < fl)
        continue;
      F[iv] -= Fc;
      N[iv] -= Nc;
      
      for (Index ij = 0; ij < nj; ij++) {
        if (not propmattype_index(jacobian_quantities, ij))
          continue;
        dF(iv, ij) -= dF_cut(0, ij);
        dN(iv, ij) -= dN_cut(0, ij);
      }
    }
  }
}

inline
void line_loop(ComplexVector &F, ComplexMatrix &dF, ComplexVector &N,
               ComplexMatrix &dN, const Vector &f_grid,
               const AbsorptionLines &band,
               const ArrayOfRetrievalQuantity &jacobian_quantities,
               const EnergyLevelMap &nlte, const Vector &vmrs,
               const Numeric &f_mean, const Numeric &isot_ratio,
               const Numeric &P, const Numeric &T, const Numeric &DC,
               const Numeric &QT, const Numeric &QT0, const Numeric &dQTdT,
               const Numeric &H, const Index &nj, const Index &nl,
               const Index &nv, bool do_nlte, bool do_zeeman,
               const Zeeman::Polarization &zeeman_polarization) ARTS_NOEXCEPT {
  // Programming trickery (sets pre-allocated size limits)
  constexpr std::size_t preallocated_count_lsmp = 4;  // Temperature and two VMR
  constexpr std::size_t preallocated_count_line = 16;
  static_assert(preallocated_count_lsmp and preallocated_count_line, "Must be above 0");
  using Preallocator = std::pmr::monotonic_buffer_resource;
  
  for (Index i = 0; i < nl; i++) {
    const Index nz = do_zeeman ? band.ZeemanCount(i, zeeman_polarization) : 1;
    
    const Output X = band.ShapeParameters(i, T, P, vmrs);
    const Complex LM(1 + X.G, -X.Y);
    
    // Line strength value
    IntensityCalculator ls_str =
    linestrength_selection(T, QT, QT0, dQTdT, isot_ratio, nlte, band, i);
    const Numeric Si = std::visit([](auto &&S) { return S.S; }, ls_str);
    const Numeric DSi = std::visit([](auto &&S) { return S.N; }, ls_str);
    
    // Line shape normalization
    Normalizer ls_norm =
    normalizer_selection(band.Normalization(), band.F0(i), T);
    
    // Derivatives
    std::array<Output, preallocated_count_lsmp> lsmp_mem;
    std::array<Numeric, preallocated_count_line> line_mem;
    Preallocator lsmp_alloc(lsmp_mem.data(), sizeof(lsmp_mem));
    Preallocator line_alloc(line_mem.data(), sizeof(line_mem));
    std::pmr::vector<Numeric> line_derivs(&line_alloc);
    std::pmr::vector<Output> lsmp_derivs(&lsmp_alloc);
    
    // Pre-compute the derivatives
    for (Index ij = 0; ij < nj; ij++) {
      if (not propmattype_index(jacobian_quantities, ij)) {
        continue;
      }
      
      const auto &deriv = jacobian_quantities[ij];
      
      if (deriv == Jacobian::Atm::Temperature) {
        lsmp_derivs.emplace_back(band.ShapeParameters_dT(i, T, P, vmrs));
      } else if (deriv.Target().needQuantumIdentity()) {
        if (deriv == Jacobian::Line::VMR) {
          lsmp_derivs.emplace_back(band.ShapeParameters_dVMR(i, T, P, deriv.QuantumIdentity()));
        } else {
          const Absorption::QuantumIdentifierLineTarget lt(deriv.Target().QuantumIdentity(), band, i);
          if (lt == Absorption::QuantumIdentifierLineTargetType::Line) {
            if (deriv == Jacobian::Line::Strength) {
              line_derivs.emplace_back(std::visit([](auto &&LS) { return LS.dSdI0(); }, ls_str));
            }
            // All line shape derivatives
            InternalDerivativesSetup(G0) InternalDerivativesSetup(D0)
            InternalDerivativesSetup(G2) InternalDerivativesSetup(D2)
            InternalDerivativesSetup(ETA) InternalDerivativesSetup(FVC)
            InternalDerivativesSetup(Y) InternalDerivativesSetup(G)
            InternalDerivativesSetup(DV)
          } else if (lt == Absorption::QuantumIdentifierLineTargetType::Level) {
            if (deriv == Jacobian::Line::NLTE) {
              if (lt.upper) {
                line_derivs.emplace_back(std::visit([](auto &&LS) { return LS.dSdNLTEu(); }, ls_str));
              }
              
              if (lt.lower) {
                line_derivs.emplace_back(std::visit([](auto &&LS) { return LS.dSdNLTEl(); }, ls_str));
              }
            }
          }
        }
      }
    }

    zeeman_loop(
        F, dF, N, dN, line_derivs, lsmp_derivs, ls_norm, ls_str, f_grid, band,
        jacobian_quantities, X, LM, Si, DSi, f_mean, T, DC, H, i, nz,
        nj, nv, do_nlte, do_zeeman, zeeman_polarization);
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
void compute(ComplexVector &F, ComplexMatrix &dF, ComplexVector &N,
             ComplexMatrix &dN, const Vector &f_grid,
             const AbsorptionLines &band,
             const ArrayOfRetrievalQuantity &jacobian_quantities,
             const EnergyLevelMap &nlte,
             const SpeciesAuxData::AuxType &partfun_type,
             const ArrayOfGriddedField1 &partfun_data, const Vector &vmrs,
             const Numeric &isot_ratio, const Numeric &P, const Numeric &T,
             const bool do_nlte, const Numeric &H, const bool do_zeeman,
             const Zeeman::Polarization zeeman_polarization) ARTS_NOEXCEPT {

  // Get sizes
  const Index nj = jacobian_quantities.nelem();
  const Index nl = band.NumLines();
  const Index nv = f_grid.nelem();

  const Numeric DC = band.DopplerConstant(T);

  // Tests that must be true while calling this function
  ARTS_ASSERT(H >= 0, "Only for positive H.  You provided: ", H)
  ARTS_ASSERT(P > 0, "Only for abs positive P.  You provided: ", P)
  ARTS_ASSERT(T > 0, "Only for abs positive T.  You provided: ", T)
  ARTS_ASSERT(band.OK(), "Band is poorly constructed.  You need to use "
                         "a detailed debugger to find out why.")
  ARTS_ASSERT(F.size() == nv,
              "F is wrong size.  Size is (", F.size(), ") but should be: (", nv, ')')
  ARTS_ASSERT(not do_nlte or N.size() == nv,
              "N is wrong size.  Size is (", N.size(), ") but should be (", nv, ')')
  ARTS_ASSERT(nj == 0 or (dF.nrows() == nv and dF.ncols() == nj),
              "dF is wrong size.  Size is (", dF.nrows(), " x ", dF.ncols(), ") but should be: (", nv, " x ", nj, ")")
  ARTS_ASSERT(nj == 0 or not do_nlte or (dN.nrows() == nv and dN.ncols() == nj),
              "dN is wrong size.  Size is (", dN.nrows(), " x ", dN.ncols(), ") but should be: (", nv, " x ", nj, ")")
  
  // Implementation error
  ARTS_ASSERT(do_nlte and nj, "Still no support for Non-LTE derivatives")

  // Early return test
  if (nl == 0 or (Absorption::relaxationtype_relmat(band.Population()) and
                  band.DoLineMixing(P))) {
    return; // No line-by-line computations required/wanted
  }

  // Mean averaged frequency
  const Numeric f_mean = band.F_mean();

  // Partition function at T
  const Numeric QT = single_partition_function(T, partfun_type, partfun_data);

  // Partition function derivative wrt T at T
  const Numeric dQTdT =
      dsingle_partition_function_dT(T, partfun_type, partfun_data);

  // Partition function at T0
  const Numeric QT0 =
      single_partition_function(band.T0(), partfun_type, partfun_data);
      
  line_loop(F, dF, N, dN, f_grid,
            band, jacobian_quantities, nlte, vmrs, f_mean,
            isot_ratio, P, T, DC, QT, QT0, dQTdT, H, nj, nl,
            nv, do_nlte, do_zeeman,zeeman_polarization);
}

#undef InternalDerivatives
#undef InternalDerivativesG
#undef InternalDerivativesY
} // namespace LineShape
