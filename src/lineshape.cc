#include <Faddeeva/Faddeeva.hh>

#include "lineshape.h"

namespace LineShape {
Doppler::Doppler (Numeric F0_noshift, Numeric DC, Numeric dZ) noexcept
    : mF0(F0_noshift + dZ),
    invGD(Constant::sqrt_ln_2 / (DC * mF0)),
    x(-mF0 * invGD),
    F(invGD * Constant::inv_sqrt_pi * std::exp(-Constant::pow2(x)), 0) {}

Complex Doppler::dFdT(Output, Numeric T) const noexcept {
  return F * (2 * Constant::pow2(x) - 1) / (2 * T);
}

Complex Doppler::dFdf() const noexcept { return -2 * invGD * F * x; }

Complex Doppler::dFdF0() const noexcept { return F*(2*x*(invGD*mF0 + x) - 1)/mF0; }

Complex Doppler::dFdH(Numeric dZ) const noexcept { return dZ * dFdF0(); }

Complex Doppler::operator()(Numeric f) noexcept {
  x = (f - mF0) * invGD;
  F = invGD * Constant::inv_sqrt_pi * std::exp(-Constant::pow2(x));
  return F;
}

Voigt::Voigt(Numeric F0_noshift, Output ls, Numeric DC, Numeric dZ) noexcept
      : mF0(F0_noshift + dZ + ls.D0 + ls.DV),
        invGD(Constant::sqrt_ln_2 / (DC * mF0)),
        z(invGD * Complex(-mF0, ls.G0)),
        F(Constant::inv_sqrt_pi * invGD * Faddeeva::w(z)),
        dF(2 * (Complex(0, invGD * Constant::inv_pi) - z * F)) {}

Complex Voigt::dFdf() const noexcept { return dF; }

Complex Voigt::dFdF0() const noexcept { return -F / mF0 - dF; }

Complex Voigt::dFdDV(Numeric d) const noexcept { return d * dFdF0(); }

Complex Voigt::dFdD0(Numeric d) const noexcept { return d * dFdF0(); }

Complex Voigt::dFdG0(Numeric d) const noexcept { return Complex(0, d) * dF; }

Complex Voigt::dFdH(Numeric dZ) const noexcept { return -dZ * dF; }

Complex Voigt::dFdVMR(Output d) const noexcept { return Complex(-d.D0 - d.DV, d.G0) * dF; }

Complex Voigt::dFdT(Output d, Numeric T) const noexcept {
  return -(F * invGD + dF * z) * (2 * T * (d.D0 + d.DV) + mF0) /
              (2 * T * invGD * mF0) +
          dFdVMR(d);
}

Complex Voigt::operator()(Numeric f) noexcept {
  reinterpret_cast<Numeric(&)[2]>(z)[0] = invGD * (f - mF0);
  F = Constant::inv_sqrt_pi * invGD * Faddeeva::w(z);
  dF = 2 * invGD * (Complex(0, Constant::inv_pi * invGD) - z * F);
  return F;
}

constexpr Numeric ln_16 = 2.772588722239781237668928485832706272302000537441021016482720037973574487879;

SpeedDependentVoigt::SpeedDependentVoigt(Numeric F0_noshift, Output ls, Numeric GD_div_F0, Numeric dZ) noexcept
    : mF0(F0_noshift + dZ + ls.D0 - 1.5 * ls.D2),
      invGD(Constant::sqrt_ln_2 / (GD_div_F0 * mF0)),
      invc2(1.0 / Complex(ls.G2, ls.D2)),
      dx(Complex(ls.G0 - 1.5 * ls.G2, mF0)),
      x(dx * invc2),
      sqrty(invc2 / (2 * invGD)),
      calcs(init(Complex(ls.G2, ls.D2))) {
  calc();
}

Complex SpeedDependentVoigt::dFdf() const noexcept {
  switch (calcs) {
    case CalcType::Full:
      return invGD * invc2 * (dw1 - dw2) / (2 * Constant::sqrt_pi * sq);
    case CalcType::Voigt:
      return dw1 * Constant::pow2(invGD) * Constant::inv_sqrt_pi;
    case CalcType::LowXandHighY:
      return dw1 * Constant::pow2(invGD) * Constant::inv_sqrt_pi -
      dw2 * invGD * invc2 / (2 * Constant::sqrt_pi * sq);
    case CalcType::LowYandLowX:
      return Constant::pow2(invc2) * (-dw1 * sq + Complex(0, 1) * w1) /
      (Constant::sqrt_pi * sq);
    case CalcType::LowYandHighX:
      return Complex(0, 1) * Constant::pow2(invc2) * (x - 3) /
      (Constant::pi * Constant::pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdF0() const noexcept {
  switch (calcs) {
    case CalcType::Full:
      return (4 * Constant::pow2(invGD) * (-w1 + w2) * sq +
      Complex(0, 1) * invc2 *
      (dw1 * (Complex(0, 2 * mF0 * Constant::pow2(invGD)) -
      2 * invGD * sq + invc2) -
      dw2 * (Complex(0, 2 * mF0 * Constant::pow2(invGD)) +
      2 * invGD * sq + invc2))) /
      (4 * Constant::sqrt_pi * invGD * mF0 * sq);
    case CalcType::Voigt:
      return -invGD *
      (Complex(0, invGD) * dw1 * (dx - Complex(0, mF0)) + w1) /
      (Constant::sqrt_pi * mF0);
    case CalcType::LowXandHighY:
      return (4 * Constant::pow2(invGD) * (-w1 + w2) * sq -
      Complex(0, 1) *
      (4 * dw1 * Constant::pow3(invGD) * (dx - Complex(0, mF0)) *
      sq +
      dw2 * invc2 *
      (Complex(0, 2 * mF0 * Constant::pow2(invGD)) +
      2 * invGD * sq + invc2))) /
      (4 * Constant::sqrt_pi * invGD * mF0 * sq);
    case CalcType::LowYandLowX:
      return Constant::pow2(invc2) * (dw1 * sq - Complex(0, 1) * w1) /
      (Constant::sqrt_pi * sq);
    case CalcType::LowYandHighX:
      return Complex(0, 1) * Constant::pow2(invc2) * (3 - x) /
      (Constant::pi * Constant::pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdD0(Numeric dD0dD0) const noexcept {
  switch (calcs) {
    case CalcType::Full:
      return -dD0dD0 *
      (4 * Constant::pow2(invGD) * (w1 - w2) * sq +
      Complex(0, 1) * invc2 *
      (-dw1 * (Complex(0, 2 * mF0 * Constant::pow2(invGD)) -
      2 * invGD * sq + invc2) +
      dw2 * (Complex(0, 2 * mF0 * Constant::pow2(invGD)) +
      2 * invGD * sq + invc2))) /
      (4 * Constant::sqrt_pi * invGD * mF0 * sq);
    case CalcType::Voigt:
      return -dD0dD0 * invGD *
      (Complex(0, invGD) * dw1 * (dx - Complex(0, mF0)) + w1) /
      (Constant::sqrt_pi * mF0);
    case CalcType::LowXandHighY:
      return -dD0dD0 *
      (4 * Constant::pow2(invGD) * (w1 - w2) * sq +
      Complex(0, 1) *
      (4 * dw1 * Constant::pow3(invGD) * (dx - Complex(0, mF0)) *
      sq +
      dw2 * invc2 *
      (Complex(0, 2 * mF0 * Constant::pow2(invGD)) +
      2 * invGD * sq + invc2))) /
      (4 * Constant::sqrt_pi * invGD * mF0 * sq);
    case CalcType::LowYandLowX:
      return dD0dD0 * Constant::pow2(invc2) *
      (dw1 * sq - Complex(0, 1) * w1) / (Constant::sqrt_pi * sq);
    case CalcType::LowYandHighX:
      return -Complex(0, dD0dD0) * Constant::pow2(invc2) * (x - 3) /
      (Constant::pi * Constant::pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdG0(Numeric dG0dG0) const noexcept {
  switch (calcs) {
    case CalcType::Full:
      return Complex(0, dG0dG0) * invGD * invc2 * (dw1 - dw2) /
      (2 * Constant::sqrt_pi * sq);
    case CalcType::Voigt:
      return Complex(0, dG0dG0) * dw1 * Constant::pow2(invGD) *
      Constant::inv_sqrt_pi;
    case CalcType::LowXandHighY:
      return Complex(0, dG0dG0) * invGD *
      (2 * dw1 * invGD * sq - dw2 * invc2) /
      (2 * Constant::sqrt_pi * sq);
    case CalcType::LowYandLowX:
      return -dG0dG0 * Constant::pow2(invc2) *
      (Complex(0, 1) * dw1 * sq + w1) / (Constant::sqrt_pi * sq);
    case CalcType::LowYandHighX:
      return -dG0dG0 * Constant::pow2(invc2) * (x - 3) /
      (Constant::pi * Constant::pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdD2(Numeric dD2dD2) const noexcept {
  switch (calcs) {
    case CalcType::Full:
      return dD2dD2 *
      (12 * Constant::pow2(invGD) * (w1 - w2) * sq +
      Complex(0, 1) * invc2 *
      (dw1 * (-Complex(0, 2 * mF0 * Constant::pow2(invGD)) *
      (2 * dx * invc2 + 3) +
      4 * Complex(0, invGD) * invc2 * mF0 * sq +
      6 * invGD * sq -
      Complex(0, 2 * mF0) * Constant::pow2(invc2) -
      3 * invc2) +
      dw2 * (Complex(0, 2 * mF0 * Constant::pow2(invGD)) *
      (2 * dx * invc2 + 3) +
      4 * Complex(0, invGD) * invc2 * mF0 * sq +
      6 * invGD * sq +
      Complex(0, 2 * mF0) * Constant::pow2(invc2) +
      3 * invc2))) /
      (8 * Constant::sqrt_pi * invGD * mF0 * sq);
    case CalcType::Voigt:
      return 3 * dD2dD2 * invGD *
      (Complex(0, invGD) * dw1 * (dx - Complex(0, mF0)) + w1) /
      (2 * Constant::sqrt_pi * mF0);
    case CalcType::LowXandHighY:
      return dD2dD2 *
      (12 * Constant::pow2(invGD) * (w1 - w2) * sq +
      Complex(0, 1) *
      (12 * dw1 * Constant::pow3(invGD) * (dx - Complex(0, mF0)) *
      sq +
      dw2 * invc2 *
      (Complex(0, 2 * mF0 * Constant::pow2(invGD)) *
      (2 * dx * invc2 + 3) +
      4 * Complex(0, invGD) * invc2 * mF0 * sq +
      6 * invGD * sq +
      Complex(0, 2 * mF0) * Constant::pow2(invc2) +
      3 * invc2))) /
      (8 * Constant::sqrt_pi * invGD * mF0 * sq);
    case CalcType::LowYandLowX:
      return dD2dD2 * Constant::pow2(invc2) *
      (4 * Complex(0, 1) * sq * (Constant::sqrt_pi * w1 * sq - 1) -
      Constant::sqrt_pi * (dw1 * sq - Complex(0, 1) * w1) *
      (2 * dx * invc2 + 3)) /
      (2 * Constant::pi * sq);
    case CalcType::LowYandHighX:
      return Complex(0, dD2dD2) * Constant::pow2(invc2) *
      (-x * (2 * x - 3) + (x - 3) * (2 * dx * invc2 + 3)) /
      (2 * Constant::pi * Constant::pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdG2(Numeric dG2dG2) const noexcept {
  switch (calcs) {
    case CalcType::Full:
      return Complex(0, dG2dG2) * invc2 *
      (dw1 * (-Constant::pow2(invGD) * (2 * dx * invc2 + 3) +
      2 * invGD * invc2 * sq - Constant::pow2(invc2)) +
      dw2 * (Constant::pow2(invGD) * (2 * dx * invc2 + 3) +
      2 * invGD * invc2 * sq + Constant::pow2(invc2))) /
      (4 * Constant::sqrt_pi * invGD * sq);
    case CalcType::Voigt:
      return -3 * Complex(0, dG2dG2) * dw1 * Constant::pow2(invGD) /
      (2 * Constant::sqrt_pi);
    case CalcType::LowXandHighY:
      return Complex(0, dG2dG2) *
      (-6 * dw1 * Constant::pow3(invGD) * sq +
      dw2 * invc2 *
      (Constant::pow2(invGD) * (2 * dx * invc2 + 3) +
      2 * invGD * invc2 * sq + Constant::pow2(invc2))) /
      (4 * Constant::sqrt_pi * invGD * sq);
    case CalcType::LowYandLowX:
      return dG2dG2 * Constant::pow2(invc2) *
      (4 * sq * (Constant::sqrt_pi * w1 * sq - 1) +
      Constant::sqrt_pi * (2 * dx * invc2 + 3) *
      (Complex(0, 1) * dw1 * sq + w1)) /
      (2 * Constant::pi * sq);
    case CalcType::LowYandHighX:
      return dG2dG2 * Constant::pow2(invc2) *
      (-x * (2 * x - 3) + (x - 3) * (2 * dx * invc2 + 3)) /
      (2 * Constant::pi * Constant::pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdH(Numeric dZ) const noexcept {
  switch (calcs) {
    case CalcType::Full:
      return -dZ *
      (4 * Constant::pow2(invGD) * (w1 - w2) * sq +
      Complex(0, 1) * invc2 *
      (-dw1 * (Complex(0, 2 * mF0 * Constant::pow2(invGD)) -
      2 * invGD * sq + invc2) +
      dw2 * (Complex(0, 2 * mF0 * Constant::pow2(invGD)) +
      2 * invGD * sq + invc2))) /
      (4 * Constant::sqrt_pi * invGD * mF0 * sq);
    case CalcType::Voigt:
      return -dZ * invGD *
      (Complex(0, invGD) * dw1 * (dx - Complex(0, mF0)) + w1) /
      (Constant::sqrt_pi * mF0);
    case CalcType::LowXandHighY:
      return -dZ *
      (4 * Constant::pow2(invGD) * (w1 - w2) * sq +
      Complex(0, 1) *
      (4 * dw1 * Constant::pow3(invGD) * (dx - Complex(0, mF0)) *
      sq +
      dw2 * invc2 *
      (Complex(0, 2 * mF0 * Constant::pow2(invGD)) +
      2 * invGD * sq + invc2))) /
      (4 * Constant::sqrt_pi * invGD * mF0 * sq);
    case CalcType::LowYandLowX:
      return dZ * Constant::pow2(invc2) * (dw1 * sq - Complex(0, 1) * w1) /
      (Constant::sqrt_pi * sq);
    case CalcType::LowYandHighX:
      return -Complex(0, dZ) * Constant::pow2(invc2) * (x - 3) /
      (Constant::pi * Constant::pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdVMR(Output d) const noexcept {
  switch (calcs) {
    case CalcType::Full:
      return (-4 * Constant::pow2(invGD) * (2 * d.D0 - 3 * d.D2) * (w1 - w2) *
      sq +
      Complex(0, 1) * invc2 *
      (dw1 * (-2 * Constant::pow2(invGD) * mF0 *
      (Complex(3 * d.G2 - 2 * d.G0,
               3 * d.D2 - 2 * d.D0) +
      2 * dx * invc2 * Complex(d.G2, d.D2)) +
      4 * invGD * invc2 * mF0 * sq *
      Complex(d.G2, d.D2) -
      2 * invGD * (2 * d.D0 - 3 * d.D2) * sq -
      2 * Constant::pow2(invc2) * mF0 *
      Complex(d.G2, d.D2) +
      invc2 * (2 * d.D0 - 3 * d.D2)) -
      dw2 * (-2 * Constant::pow2(invGD) * mF0 *
      (Complex(3 * d.G2 - 2 * d.G0,
               3 * d.D2 - 2 * d.D0) +
      2 * dx * invc2 * Complex(d.G2, d.D2)) -
      4 * invGD * invc2 * mF0 * sq *
      Complex(d.G2, d.D2) +
      2 * invGD * (2 * d.D0 - 3 * d.D2) * sq -
      2 * Constant::pow2(invc2) * mF0 *
      Complex(d.G2, d.D2) +
      invc2 * (2 * d.D0 - 3 * d.D2)))) /
      (8 * Constant::sqrt_pi * invGD * mF0 * sq);
    case CalcType::Voigt:
      return -invGD *
      (Complex(0, invGD) * dw1 *
      (dx * (2 * d.D0 - 3 * d.D2) -
      mF0 * Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2)) +
      w1 * (2 * d.D0 - 3 * d.D2)) /
      (2 * Constant::sqrt_pi * mF0);
    case CalcType::LowXandHighY:
      return -(4 * Constant::pow2(invGD) * (2 * d.D0 - 3 * d.D2) * (w1 - w2) *
      sq +
      Complex(0, 1) *
      (4 * dw1 * Constant::pow3(invGD) * sq *
      (dx * (2 * d.D0 - 3 * d.D2) -
      mF0 * Complex(2 * d.G0 - 3 * d.G2,
                    2 * d.D0 - 3 * d.D2)) +
      dw2 * invc2 *
      (2 * Constant::pow2(invGD) * mF0 *
      (Complex(2 * d.G0 - 3 * d.G2,
               2 * d.D0 - 3 * d.D2) -
               2 * dx * invc2 * Complex(d.G2, d.D2)) -
               4 * invGD * invc2 * mF0 * sq *
               Complex(d.G2, d.D2) +
               2 * invGD * (2 * d.D0 - 3 * d.D2) * sq -
               2 * Constant::pow2(invc2) * mF0 *
               Complex(d.G2, d.D2) +
               invc2 * (2 * d.D0 - 3 * d.D2)))) /
               (8 * Constant::sqrt_pi * invGD * mF0 * sq);
    case CalcType::LowYandLowX:
      return Constant::pow2(invc2) *
      (4 * sq * Complex(d.G2, d.D2) *
      (Constant::sqrt_pi * w1 * sq - 1) -
      Constant::sqrt_pi * (Complex(0, 1) * dw1 * sq + w1) *
      (Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) -
      2 * dx * invc2 * Complex(d.G2, d.D2))) /
      (2 * Constant::pi * sq);
    case CalcType::LowYandHighX:
      return -Constant::pow2(invc2) *
      (x * (2 * x - 3) * Complex(d.G2, d.D2) +
      (x - 3) * (Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) -
      2 * dx * invc2 * Complex(d.G2, d.D2))) /
      (2 * Constant::pi * Constant::pow3(x));
  }
  return {};
}

Complex SpeedDependentVoigt::dFdT(Output d, Numeric T) const noexcept {
  switch (calcs) {
    case CalcType::Full:
      return (-Constant::pow2(invGD) * (w1 - w2) * sq *
      (T * (2 * d.D0 - 3 * d.D2) * Constant::sqrt_ln_2 -
      Constant::pow3(invGD) * mF0) *
      ln_16 +
      Complex(0, 1) * invc2 *
      (dw1 *
      (T * invGD * invc2 * mF0 * sq *
      Complex(d.G2, d.D2) * ln_16 -
      2 * invGD * sq *
      (T * (2 * d.D0 - 3 * d.D2) * Constant::sqrt_ln_2 -
      Constant::pow3(invGD) * mF0) *
      Constant::sqrt_ln_2 -
      (2 * T * Constant::pow2(invGD) * mF0 *
      (Complex(3 * d.G2 - 2 * d.G0,
               3 * d.D2 - 2 * d.D0) +
      2 * dx * invc2 * Complex(d.G2, d.D2)) *
      Constant::sqrt_ln_2 +
      2 * T * Constant::pow2(invc2) * mF0 *
      Complex(d.G2, d.D2) * Constant::sqrt_ln_2 +
      invc2 * (T * (-2 * d.D0 + 3 * d.D2) *
      Constant::sqrt_ln_2 +
      Constant::pow3(invGD) * mF0)) *
      Constant::sqrt_ln_2) +
      dw2 * (T * invGD * invc2 * mF0 * sq *
      Complex(d.G2, d.D2) * ln_16 +
      2 * invGD * sq *
      (T * (-2 * d.D0 + 3 * d.D2) *
      Constant::sqrt_ln_2 +
      Constant::pow3(invGD) * mF0) *
      Constant::sqrt_ln_2 +
      (-2 * T * Constant::pow2(invGD) * mF0 *
      (Complex(2 * d.G0 - 3 * d.G2,
               2 * d.D0 - 3 * d.D2) -
               2 * dx * invc2 * Complex(d.G2, d.D2)) *
               Constant::sqrt_ln_2 +
               2 * T * Constant::pow2(invc2) * mF0 *
               Complex(d.G2, d.D2) * Constant::sqrt_ln_2 +
               invc2 * (T * (-2 * d.D0 + 3 * d.D2) *
               Constant::sqrt_ln_2 +
               Constant::pow3(invGD) * mF0)) *
               Constant::sqrt_ln_2)) *
               Constant::sqrt_ln_2) /
               (8 * Constant::sqrt_pi * T * invGD * mF0 * sq *
               Constant::pow3(Constant::sqrt_ln_2));
    case CalcType::Voigt:
      return -invGD *
      (-Complex(0, invGD) * dw1 *
      (T * mF0 *
      Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) *
      Constant::sqrt_ln_2 -
      dx * (T * (2 * d.D0 - 3 * d.D2) * Constant::sqrt_ln_2 -
      Constant::pow3(invGD) * mF0)) +
      w1 * (T * (2 * d.D0 - 3 * d.D2) * Constant::sqrt_ln_2 -
      Constant::pow3(invGD) * mF0)) /
      (2 * Constant::sqrt_pi * T * mF0 * Constant::sqrt_ln_2);
    case CalcType::LowXandHighY:
      return (-Constant::pow2(invGD) * (w1 - w2) * sq *
      (T * (2 * d.D0 - 3 * d.D2) * Constant::sqrt_ln_2 -
      Constant::pow3(invGD) * mF0) *
      ln_16 +
      Complex(0, 1) *
      (dw1 * Constant::pow3(invGD) * sq *
      (T * mF0 *
      Complex(2 * d.G0 - 3 * d.G2,
              2 * d.D0 - 3 * d.D2) *
              Constant::sqrt_ln_2 -
              dx *
              (T * (2 * d.D0 - 3 * d.D2) * Constant::sqrt_ln_2 -
              Constant::pow3(invGD) * mF0)) *
              ln_16 +
              dw2 * invc2 *
              (T * invGD * invc2 * mF0 * sq *
              Complex(d.G2, d.D2) * ln_16 +
              2 * invGD * sq *
              (T * (-2 * d.D0 + 3 * d.D2) *
              Constant::sqrt_ln_2 +
              Constant::pow3(invGD) * mF0) *
              Constant::sqrt_ln_2 +
              (-2 * T * Constant::pow2(invGD) * mF0 *
              (Complex(2 * d.G0 - 3 * d.G2,
                       2 * d.D0 - 3 * d.D2) -
                       2 * dx * invc2 * Complex(d.G2, d.D2)) *
                       Constant::sqrt_ln_2 +
                       2 * T * Constant::pow2(invc2) * mF0 *
                       Complex(d.G2, d.D2) * Constant::sqrt_ln_2 +
                       invc2 * (T * (-2 * d.D0 + 3 * d.D2) *
                       Constant::sqrt_ln_2 +
                       Constant::pow3(invGD) * mF0)) *
                       Constant::sqrt_ln_2) *
                       Constant::sqrt_ln_2)) /
                       (8 * Constant::sqrt_pi * T * invGD * mF0 * sq *
                       Constant::pow3(Constant::sqrt_ln_2));
    case CalcType::LowYandLowX:
      return Constant::pow2(invc2) *
      (4 * sq * Complex(d.G2, d.D2) *
      (Constant::sqrt_pi * w1 * sq - 1) -
      Constant::sqrt_pi * (Complex(0, 1) * dw1 * sq + w1) *
      (Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) -
      2 * dx * invc2 * Complex(d.G2, d.D2))) /
      (2 * Constant::pi * sq);
    case CalcType::LowYandHighX:
      return -Constant::pow2(invc2) *
      (x * (2 * x - 3) * Complex(d.G2, d.D2) +
      (x - 3) * (Complex(2 * d.G0 - 3 * d.G2, 2 * d.D0 - 3 * d.D2) -
      2 * dx * invc2 * Complex(d.G2, d.D2))) /
      (2 * Constant::pi * Constant::pow3(x));
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

constexpr Numeric abs_squared(Complex z) noexcept { return Constant::pow2(z.real()) + Constant::pow2(z.imag()); }

SpeedDependentVoigt::CalcType SpeedDependentVoigt::init(const Complex c2) const noexcept {
  if (abs_squared(c2) == 0)
    return CalcType::Voigt;
  else if (abs_squared(x) <= 9e-16 * abs_squared(sqrty * sqrty))
    return CalcType::LowXandHighY;
  else if ((abs_squared(sqrty * sqrty) <= 1.e-30 * abs_squared(x)) and
    abs_squared(std::sqrt(x)) <= 16.e6)
    return CalcType::LowYandLowX;  // Weird case, untested
    else if ((abs_squared(sqrty * sqrty) <= 1.e-30 * abs_squared(x)))
      return CalcType::LowYandHighX;
    else
      return CalcType::Full;
}

void SpeedDependentVoigt::update_calcs() noexcept {
  if (calcs not_eq CalcType::Voigt) calcs = init(Complex(1, 1));
}

void SpeedDependentVoigt::calc() noexcept {
  switch (calcs) {
    case CalcType::Full:
      sq = std::sqrt(x + sqrty * sqrty);
      w1 = Faddeeva::w(Complex(0, 1) * (sq - sqrty));
      w2 = Faddeeva::w(Complex(0, 1) * (sq + sqrty));
      F = Constant::inv_sqrt_pi * invGD * (w1 - w2);
      dw1 = Complex(0, 2) * (Constant::inv_sqrt_pi - (sq - sqrty) * w1);
      dw2 = Complex(0, 2) * (Constant::inv_sqrt_pi - (sq + sqrty) * w2);
      break;
    case CalcType::Voigt:
      w1 = Faddeeva::w(Complex(0, 1) * dx * invGD);
      F = Constant::inv_sqrt_pi * invGD * w1;
      dw1 = Complex(0, 2) * (Constant::inv_sqrt_pi - dx * invGD * w1);
      break;
    case CalcType::LowXandHighY:
      sq = std::sqrt(x + sqrty * sqrty);
      w1 = Faddeeva::w(Complex(0, 1) * dx * invGD);
      w2 = Faddeeva::w(Complex(0, 1) * (sq + sqrty));
      F = Constant::inv_sqrt_pi * invGD * (w1 - w2);
      dw1 = Complex(0, 2) * (Constant::inv_sqrt_pi - dx * invGD * w1);
      dw2 = Complex(0, 2) * (Constant::inv_sqrt_pi - (sq + sqrty) * w2);
      break;
    case CalcType::LowYandLowX:
      sq = std::sqrt(x);
      w1 = Faddeeva::w(Complex(0, 1) * sq);
      F = 2 * Constant::inv_pi * invc2 * (1 - Constant::sqrt_pi * sq * w1);
      dw1 = Complex(0, 2) * (Constant::inv_sqrt_pi - sq * w1);
      break;
    case CalcType::LowYandHighX:
      F = Constant::inv_pi * invc2 * (1 / x - 1.5 / Constant::pow2(x));
      break;
  }
}

HartmannTran::HartmannTran(Numeric F0_noshift, Output ls, Numeric GD_div_F0, Numeric dZ) noexcept
             : G0(ls.G0),
             D0(ls.D0),
             G2(ls.G2),
             D2(ls.D2),
             FVC(ls.FVC),
             ETA(ls.ETA),
             mF0(F0_noshift + dZ + (1 - ls.ETA) * (ls.D0 - 1.5 * ls.D2)),
             invGD(Constant::sqrt_ln_2 / (GD_div_F0 * mF0)),
             deltax(ls.FVC + (1 - ls.ETA) * (ls.G0 - 3 * ls.G2 / 2), mF0),
             sqrty(1 / (2 * (1 - ls.ETA) * Complex(ls.G2, ls.D2) * invGD)) {
               calc();
             }
             
             Complex HartmannTran::dFdf() const noexcept {
               constexpr Complex ddeltax = Complex(0, -1);
               Complex dx = -ddeltax / ((ETA - 1) * Complex(G2, D2));
               Complex dsqrtxy = dx / (2 * sqrtxy);
               
               switch (calcs) {
                 case CalcType::Full: {
                   Complex dz1 = dsqrtxy;
                   Complex dz2 = dsqrtxy;
                   Complex dA =
                   Complex(0, Constant::sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
                   Complex dB = Constant::sqrt_pi *
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 -
                   (Constant::pow2(z2) - 1) * Complex(0, 1) * dw2 * dz2 +
                   2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) /
                   (2 * sqrty * (ETA - 1) * Complex(G2, D2));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tLowZ: {
                   Complex dz1 = invGD * ddeltax;
                   Complex dA = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1;
                   Complex dB =
                   -invGD * (Constant::sqrt_pi * ((Constant::pow2(z1) - 1) *
                   Complex(0, 1) * dw1 * dz1 +
                   2 * w1 * z1 * dz1) -
                   dz1);
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tHighZ: {
                   Complex dz1 = invGD * ddeltax;
                   Complex dA = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1;
                   Complex dB = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1 -
                   invGD * dz1 / (2 * Constant::pow2(z1)) +
                   9 * invGD * dz1 / (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowXandHighY: {
                   Complex dz1 = invGD * ddeltax;
                   Complex dz2 = dsqrtxy;
                   Complex dA =
                   Complex(0, Constant::sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
                   Complex dB = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1 -
                   invGD * dz1 / (2 * Constant::pow2(z1)) +
                   9 * invGD * dz1 / (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandLowX: {
                   Complex dz1 = dsqrtxy;
                   Complex dz2 = dx / (2 * sqrtx);
                   Complex dA = 2 * Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) /
                   ((ETA - 1) * Complex(G2, D2));
                   Complex dB = -(2 * Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) *
                   (2 * Constant::pow2(sqrty) + x - 1) +
                   2 * Constant::sqrt_pi * w1 * dz1 +
                   Complex(0, 2 * Constant::sqrt_pi) * z1 * dw1 * dz1 +
                   2 * (Constant::sqrt_pi * w2 * z2 - 1) * dx) /
                   ((ETA - 1) * Complex(G2, D2));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandHighX: {
                   Complex dz1 = dsqrtxy;
                   Complex dA =
                   (x - 3) * dx / ((ETA - 1) * Complex(G2, D2) * Constant::pow3(x));
                   Complex dB = (-2 * Constant::sqrt_pi *
                   (w1 * dz1 + z1 * Complex(0, 1) * dw1 * dz1) *
                   Constant::pow3(x) -
                   (x - 3) * (2 * Constant::pow2(sqrty) + x - 1) * dx +
                   (2 * x - 3) * x * dx / 2) /
                   ((ETA - 1) * Complex(G2, D2) * Constant::pow3(x));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
               }
               return {};
             }
             
             Complex HartmannTran::dFdF0() const noexcept {
               Numeric dGD = (1 / (invGD * mF0));
               Numeric dinvGD = -dGD * Constant::pow2(invGD);
               Complex dsqrty =
               dinvGD / (2 * (ETA - 1) * Complex(G2, D2) * Constant::pow2(invGD));
               constexpr Complex ddeltax = Complex(0, 1);
               Complex dx = -ddeltax / ((ETA - 1) * Complex(G2, D2));
               Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;
               
               switch (calcs) {
                 case CalcType::Full: {
                   Complex dz1 = dsqrtxy - dsqrty;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB =
                   Constant::sqrt_pi *
                   ((-(Constant::pow2(z1) - 1) * w1 + (Constant::pow2(z2) - 1) * w2) *
                   dsqrty +
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 -
                   (Constant::pow2(z2) - 1) * Complex(0, 1) * dw2 * dz2 +
                   2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) *
                   sqrty) /
                   (2 * (ETA - 1) * Complex(G2, D2) * Constant::pow2(sqrty));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tLowZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB =
                   -(Constant::sqrt_pi *
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 +
                   2 * w1 * z1 * dz1) -
                   dz1) *
                   invGD -
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 - z1) * dinvGD;
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tHighZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowXandHighY: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandLowX: {
                   Complex dz1 = dsqrtxy;
                   Complex dz2 = dx / (2 * sqrtx);
                   Complex dA = 2 * Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) /
                   ((ETA - 1) * Complex(G2, D2));
                   Complex dB = -(2 * Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) *
                   (2 * Constant::pow2(sqrty) + x - 1) +
                   2 * Constant::sqrt_pi * w1 * dz1 +
                   Complex(0, 2 * Constant::sqrt_pi) * z1 * dw1 * dz1 +
                   2 * (4 * sqrty * dsqrty + dx) *
                   (Constant::sqrt_pi * w2 * z2 - 1)) /
                   ((ETA - 1) * Complex(G2, D2));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandHighX: {
                   Complex dz1 = dsqrtxy;
                   Complex dA =
                   (x - 3) * dx / ((ETA - 1) * Complex(G2, D2) * Constant::pow3(x));
                   Complex dB = (-2 * Constant::sqrt_pi *
                   (w1 * dz1 + z1 * Complex(0, 1) * dw1 * dz1) *
                   Constant::pow3(x) +
                   (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x / 2 -
                   (x - 3) * (2 * Constant::pow2(sqrty) + x - 1) * dx) /
                   ((ETA - 1) * Complex(G2, D2) * Constant::pow3(x));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
               }
               return {};
             }
             
             Complex HartmannTran::dFdD0(Numeric dD0) const noexcept {
               Numeric dmF0 = (1 - ETA) * dD0;
               Numeric dGD = (dmF0 / (invGD * mF0));
               Numeric dinvGD = -dGD * Constant::pow2(invGD);
               Complex dsqrty =
               dinvGD / (2 * (ETA - 1) * Complex(G2, D2) * Constant::pow2(invGD));
               Complex ddeltax = Complex(0, 1 - ETA) * dD0;
               Complex dx = -ddeltax / ((ETA - 1) * Complex(G2, D2));
               Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;
               
               switch (calcs) {
                 case CalcType::Full: {
                   Complex dz1 = dsqrtxy - dsqrty;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB =
                   Constant::sqrt_pi *
                   ((-(Constant::pow2(z1) - 1) * w1 + (Constant::pow2(z2) - 1) * w2) *
                   dsqrty +
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 -
                   (Constant::pow2(z2) - 1) * Complex(0, 1) * dw2 * dz2 +
                   2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) *
                   sqrty) /
                   (2 * (ETA - 1) * Complex(G2, D2) * Constant::pow2(sqrty));
                   Complex dK = ETA * Complex(G2, D2) * dB + Complex(0, ETA * dD0) * A +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tLowZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB =
                   -(Constant::sqrt_pi *
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 +
                   2 * w1 * z1 * dz1) -
                   dz1) *
                   invGD -
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 - z1) * dinvGD;
                   Complex dK = ETA * Complex(G2, D2) * dB + Complex(0, ETA * dD0) * A +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tHighZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB + Complex(0, ETA * dD0) * A +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowXandHighY: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB + Complex(0, ETA * dD0) * A +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandLowX: {
                   Complex dz1 = dsqrtxy;
                   Complex dz2 = dx / (2 * sqrtx);
                   Complex dA = 2 * Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) /
                   ((ETA - 1) * Complex(G2, D2));
                   Complex dB = -(2 * Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) *
                   (2 * Constant::pow2(sqrty) + x - 1) +
                   2 * Constant::sqrt_pi * w1 * dz1 +
                   Complex(0, 2 * Constant::sqrt_pi) * z1 * dw1 * dz1 +
                   2 * (4 * sqrty * dsqrty + dx) *
                   (Constant::sqrt_pi * w2 * z2 - 1)) /
                   ((ETA - 1) * Complex(G2, D2));
                   Complex dK = ETA * Complex(G2, D2) * dB + Complex(0, ETA * dD0) * A +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandHighX: {
                   Complex dz1 = dsqrtxy;
                   Complex dA =
                   (x - 3) * dx / ((ETA - 1) * Complex(G2, D2) * Constant::pow3(x));
                   Complex dB = (-2 * Constant::sqrt_pi *
                   (w1 * dz1 + z1 * Complex(0, 1) * dw1 * dz1) *
                   Constant::pow3(x) +
                   (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x / 2 -
                   (x - 3) * (2 * Constant::pow2(sqrty) + x - 1) * dx) /
                   ((ETA - 1) * Complex(G2, D2) * Constant::pow3(x));
                   Complex dK = ETA * Complex(G2, D2) * dB + Complex(0, ETA * dD0) * A +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
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
                   Complex dA =
                   Complex(0, Constant::sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
                   Complex dB = Constant::sqrt_pi *
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 -
                   (Constant::pow2(z2) - 1) * Complex(0, 1) * dw2 * dz2 +
                   2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) /
                   (2 * sqrty * (ETA - 1) * Complex(G2, D2));
                   Complex dK = ETA * Complex(G2, D2) * dB + ETA * A * dG0 +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tLowZ: {
                   Complex dz1 = invGD * ddeltax;
                   Complex dA = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1;
                   Complex dB =
                   -invGD * (Constant::sqrt_pi * ((Constant::pow2(z1) - 1) *
                   Complex(0, 1) * dw1 * dz1 +
                   2 * w1 * z1 * dz1) -
                   dz1);
                   Complex dK = ETA * Complex(G2, D2) * dB + ETA * A * dG0 +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tHighZ: {
                   Complex dz1 = invGD * ddeltax;
                   Complex dA = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1;
                   Complex dB = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1 -
                   invGD * dz1 / (2 * Constant::pow2(z1)) +
                   9 * invGD * dz1 / (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB + ETA * A * dG0 +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowXandHighY: {
                   Complex dz1 = invGD * ddeltax;
                   Complex dz2 = dsqrtxy;
                   Complex dA =
                   Complex(0, Constant::sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
                   Complex dB = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1 -
                   invGD * dz1 / (2 * Constant::pow2(z1)) +
                   9 * invGD * dz1 / (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB + ETA * A * dG0 +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandLowX: {
                   Complex dz1 = dsqrtxy;
                   Complex dz2 = dx / (2 * sqrtx);
                   Complex dA = 2 * Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) /
                   ((ETA - 1) * Complex(G2, D2));
                   Complex dB = -(2 * Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) *
                   (2 * Constant::pow2(sqrty) + x - 1) +
                   2 * Constant::sqrt_pi * w1 * dz1 +
                   Complex(0, 2 * Constant::sqrt_pi) * z1 * dw1 * dz1 +
                   2 * (Constant::sqrt_pi * w2 * z2 - 1) * dx) /
                   ((ETA - 1) * Complex(G2, D2));
                   Complex dK = ETA * Complex(G2, D2) * dB + ETA * A * dG0 +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandHighX: {
                   Complex dz1 = dsqrtxy;
                   Complex dA =
                   (x - 3) * dx / ((ETA - 1) * Complex(G2, D2) * Constant::pow3(x));
                   Complex dB = (-2 * Constant::sqrt_pi *
                   (w1 * dz1 + z1 * Complex(0, 1) * dw1 * dz1) *
                   Constant::pow3(x) -
                   (x - 3) * (2 * Constant::pow2(sqrty) + x - 1) * dx +
                   (2 * x - 3) * x * dx / 2) /
                   ((ETA - 1) * Complex(G2, D2) * Constant::pow3(x));
                   Complex dK = ETA * Complex(G2, D2) * dB + ETA * A * dG0 +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
               }
               return {};
             }
             
             Complex HartmannTran::dFdD2(Numeric dD2) const noexcept {
               Numeric dmF0 = -3 * (1 - ETA) * dD2 / 2;
               Numeric dGD = (dmF0 / (invGD * mF0));
               Numeric dinvGD = -dGD * Constant::pow2(invGD);
               Complex dsqrty = (Complex(G2, D2) * dinvGD + Complex(0, invGD) * dD2) /
               (2 * (ETA - 1) * Constant::pow2(Complex(G2, D2)) *
               Constant::pow2(invGD));
               Complex ddeltax = 1.5 * Complex(0, ETA - 1) * dD2;
               Complex dx = (-Complex(G2, D2) * ddeltax + Complex(0, dD2) * deltax) /
               ((ETA - 1) * Constant::pow2(Complex(G2, D2)));
               Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;
               
               switch (calcs) {
                 case CalcType::Full: {
                   Complex dz1 = dsqrtxy - dsqrty;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB =
                   (Constant::sqrt_pi * Complex(G2, D2) *
                   ((-(Constant::pow2(z1) - 1) * w1 +
                   (Constant::pow2(z2) - 1) * w2) *
                   dsqrty +
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 -
                   (Constant::pow2(z2) - 1) * Complex(0, 1) * dw2 * dz2 +
                   2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) *
                   sqrty) -
                   Complex(0, 1) *
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 -
                   Constant::sqrt_pi * (Constant::pow2(z2) - 1) * w2 +
                   2 * sqrty) *
                   sqrty * dD2) /
                   (2 * (ETA - 1) * Constant::pow2(Complex(G2, D2)) *
                   Constant::pow2(sqrty));
                   Complex dK = ETA * Complex(G2, D2) * dB -
                   Complex(0, 1.5 * dD2 * ETA) * A +
                   Complex(0, dD2 * ETA) * B +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tLowZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB =
                   -(Constant::sqrt_pi *
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 +
                   2 * w1 * z1 * dz1) -
                   dz1) *
                   invGD -
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 - z1) * dinvGD;
                   Complex dK = ETA * Complex(G2, D2) * dB -
                   Complex(0, 1.5 * dD2 * ETA) * A +
                   Complex(0, dD2 * ETA) * B +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tHighZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB -
                   Complex(0, 1.5 * dD2 * ETA) * A +
                   Complex(0, dD2 * ETA) * B +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowXandHighY: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB -
                   Complex(0, 1.5 * dD2 * ETA) * A +
                   Complex(0, dD2 * ETA) * B +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandLowX: {
                   Complex dz1 = dsqrtxy;
                   Complex dz2 = dx / (2 * sqrtx);
                   Complex dA = 2 *
                   (Constant::sqrt_pi * Complex(G2, D2) *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) -
                   Complex(0, 1) * (Constant::sqrt_pi * w2 * z2 - 1) * dD2) /
                   ((ETA - 1) * Constant::pow2(Complex(G2, D2)));
                   Complex dB = (-2 * Complex(G2, D2) *
                   (Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) *
                   (2 * Constant::pow2(sqrty) + x - 1) +
                   Constant::sqrt_pi * w1 * dz1 +
                   Complex(0, Constant::sqrt_pi) * z1 * dw1 * dz1 +
                   (4 * sqrty * dsqrty + dx) *
                   (Constant::sqrt_pi * w2 * z2 - 1)) +
                   Complex(0, 1) *
                   (2 * Constant::sqrt_pi * w1 * z1 +
                   2 * (Constant::sqrt_pi * w2 * z2 - 1) *
                   (2 * Constant::pow2(sqrty) + x - 1) -
                   1) *
                   dD2) /
                   ((ETA - 1) * Constant::pow2(Complex(G2, D2)));
                   Complex dK = ETA * Complex(G2, D2) * dB -
                   Complex(0, 1.5 * dD2 * ETA) * A +
                   Complex(0, dD2 * ETA) * B +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandHighX: {
                   Complex dz1 = dsqrtxy;
                   Complex dA =
                   (Complex(G2, D2) * (x - 3) * dx +
                   Complex(0, 1) * (2 * x - 3) * x * dD2 / 2) /
                   ((ETA - 1) * Constant::pow2(Complex(G2, D2)) * Constant::pow3(x));
                   Complex dB =
                   (Complex(G2, D2) *
                   (-4 * Constant::sqrt_pi *
                   (w1 * dz1 + z1 * Complex(0, 1) * dw1 * dz1) *
                   Constant::pow3(x) +
                   (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x -
                   2 * (x - 3) * (2 * Constant::pow2(sqrty) + x - 1) * dx) -
                   Complex(0, 1) *
                   (2 * (-2 * Constant::sqrt_pi * w1 * z1 + 1) *
                   Constant::pow2(x) +
                   (2 * x - 3) * (2 * Constant::pow2(sqrty) + x - 1)) *
                   x * dD2) /
                   (2 * (ETA - 1) * Constant::pow2(Complex(G2, D2)) *
                   Constant::pow3(x));
                   Complex dK = ETA * Complex(G2, D2) * dB -
                   Complex(0, 1.5 * dD2 * ETA) * A +
                   Complex(0, dD2 * ETA) * B +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
               }
               return {};
             }
             
             Complex HartmannTran::dFdG2(Numeric dG2) const noexcept {
               Complex dsqrty =
               dG2 / (2 * invGD * (ETA - 1) * Constant::pow2(Complex(G2, D2)));
               Numeric ddeltax = 3 * (ETA - 1) * dG2 / 2;
               Complex dx = (-Complex(G2, D2) * ddeltax + deltax * dG2) /
               ((ETA - 1) * Constant::pow2(Complex(G2, D2)));
               Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;
               
               switch (calcs) {
                 case CalcType::Full: {
                   Complex dz1 = dsqrtxy - dsqrty;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Complex(0, Constant::sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
                   Complex dB =
                   (Constant::sqrt_pi * Complex(G2, D2) *
                   ((-(Constant::pow2(z1) - 1) * w1 +
                   (Constant::pow2(z2) - 1) * w2) *
                   dsqrty +
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 -
                   (Constant::pow2(z2) - 1) * Complex(0, 1) * dw2 * dz2 +
                   2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) *
                   sqrty) -
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 -
                   Constant::sqrt_pi * (Constant::pow2(z2) - 1) * w2 + 2 * sqrty) *
                   sqrty * dG2) /
                   (2 * (ETA - 1) * Constant::pow2(Complex(G2, D2)) *
                   Constant::pow2(sqrty));
                   Complex dK = ETA * Complex(G2, D2) * dB - 3 * ETA * A * dG2 / 2 +
                   ETA * B * dG2 +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tLowZ: {
                   Complex dz1 = invGD * ddeltax;
                   Complex dA = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1;
                   Complex dB =
                   -invGD * (Constant::sqrt_pi * ((Constant::pow2(z1) - 1) *
                   Complex(0, 1) * dw1 * dz1 +
                   2 * w1 * z1 * dz1) -
                   dz1);
                   Complex dK = ETA * Complex(G2, D2) * dB - 3 * ETA * A * dG2 / 2 +
                   ETA * B * dG2 +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tHighZ: {
                   Complex dz1 = invGD * ddeltax;
                   Complex dA = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1;
                   Complex dB = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1 -
                   invGD * dz1 / (2 * Constant::pow2(z1)) +
                   9 * invGD * dz1 / (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB - 3 * ETA * A * dG2 / 2 +
                   ETA * B * dG2 +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowXandHighY: {
                   Complex dz1 = invGD * ddeltax;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Complex(0, Constant::sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
                   Complex dB = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1 -
                   invGD * dz1 / (2 * Constant::pow2(z1)) +
                   9 * invGD * dz1 / (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB - 3 * ETA * A * dG2 / 2 +
                   ETA * B * dG2 +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandLowX: {
                   Complex dz1 = dsqrtxy;
                   Complex dz2 = dx / (2 * sqrtx);
                   Complex dA = 2 *
                   (Constant::sqrt_pi * Complex(G2, D2) *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) -
                   (Constant::sqrt_pi * w2 * z2 - 1) * dG2) /
                   ((ETA - 1) * Constant::pow2(Complex(G2, D2)));
                   Complex dB = (-2 * Complex(G2, D2) *
                   (Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) *
                   (2 * Constant::pow2(sqrty) + x - 1) +
                   Constant::sqrt_pi * w1 * dz1 +
                   Complex(0, Constant::sqrt_pi) * z1 * dw1 * dz1 +
                   (4 * sqrty * dsqrty + dx) *
                   (Constant::sqrt_pi * w2 * z2 - 1)) +
                   (2 * Constant::sqrt_pi * w1 * z1 +
                   2 * (Constant::sqrt_pi * w2 * z2 - 1) *
                   (2 * Constant::pow2(sqrty) + x - 1) -
                   1) *
                   dG2) /
                   ((ETA - 1) * Constant::pow2(Complex(G2, D2)));
                   Complex dK = ETA * Complex(G2, D2) * dB - 3 * ETA * A * dG2 / 2 +
                   ETA * B * dG2 +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandHighX: {
                   Complex dz1 = dsqrtxy;
                   Complex dA =
                   (Complex(G2, D2) * (x - 3) * dx + (2 * x - 3) * x * dG2 / 2) /
                   ((ETA - 1) * Constant::pow2(Complex(G2, D2)) * Constant::pow3(x));
                   Complex dB =
                   (Complex(G2, D2) *
                   (-4 * Constant::sqrt_pi *
                   (w1 * dz1 + z1 * Complex(0, 1) * dw1 * dz1) *
                   Constant::pow3(x) +
                   (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x -
                   2 * (x - 3) * (2 * Constant::pow2(sqrty) + x - 1) * dx) -
                   (2 * (-2 * Constant::sqrt_pi * w1 * z1 + 1) * Constant::pow2(x) +
                   (2 * x - 3) * (2 * Constant::pow2(sqrty) + x - 1)) *
                   x * dG2) /
                   (2 * (ETA - 1) * Constant::pow2(Complex(G2, D2)) *
                   Constant::pow3(x));
                   Complex dK = ETA * Complex(G2, D2) * dB - 3 * ETA * A * dG2 / 2 +
                   ETA * B * dG2 +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
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
                   Complex dA =
                   Complex(0, Constant::sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
                   Complex dB = Constant::sqrt_pi *
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 -
                   (Constant::pow2(z2) - 1) * Complex(0, 1) * dw2 * dz2 +
                   2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) /
                   (2 * sqrty * (ETA - 1) * Complex(G2, D2));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA -
                   A * dFVC;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tLowZ: {
                   Complex dz1 = invGD * ddeltax;
                   Complex dA = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1;
                   Complex dB =
                   -invGD * (Constant::sqrt_pi * ((Constant::pow2(z1) - 1) *
                   Complex(0, 1) * dw1 * dz1 +
                   2 * w1 * z1 * dz1) -
                   dz1);
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA -
                   A * dFVC;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tHighZ: {
                   Complex dz1 = invGD * ddeltax;
                   Complex dA = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1;
                   Complex dB = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1 -
                   invGD * dz1 / (2 * Constant::pow2(z1)) +
                   9 * invGD * dz1 / (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA -
                   A * dFVC;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowXandHighY: {
                   Complex dz1 = invGD * ddeltax;
                   Complex dz2 = dsqrtxy;
                   Complex dA =
                   Complex(0, Constant::sqrt_pi * invGD) * (dw1 * dz1 - dw2 * dz2);
                   Complex dB = Complex(0, Constant::sqrt_pi * invGD) * dw1 * dz1 -
                   invGD * dz1 / (2 * Constant::pow2(z1)) +
                   9 * invGD * dz1 / (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA -
                   A * dFVC;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandLowX: {
                   Complex dz1 = dsqrtxy;
                   Complex dz2 = dx / (2 * sqrtx);
                   Complex dA = 2 * Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) /
                   ((ETA - 1) * Complex(G2, D2));
                   Complex dB = -(2 * Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) *
                   (2 * Constant::pow2(sqrty) + x - 1) +
                   2 * Constant::sqrt_pi * w1 * dz1 +
                   Complex(0, 2 * Constant::sqrt_pi) * z1 * dw1 * dz1 +
                   2 * (Constant::sqrt_pi * w2 * z2 - 1) * dx) /
                   ((ETA - 1) * Complex(G2, D2));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA -
                   A * dFVC;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandHighX: {
                   Complex dz1 = dsqrtxy;
                   Complex dA =
                   (x - 3) * dx / ((ETA - 1) * Complex(G2, D2) * Constant::pow3(x));
                   Complex dB = (-2 * Constant::sqrt_pi *
                   (w1 * dz1 + z1 * Complex(0, 1) * dw1 * dz1) *
                   Constant::pow3(x) -
                   (x - 3) * (2 * Constant::pow2(sqrty) + x - 1) * dx +
                   (2 * x - 3) * x * dx / 2) /
                   ((ETA - 1) * Complex(G2, D2) * Constant::pow3(x));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA -
                   A * dFVC;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
               }
               return {};
             }
             
             Complex HartmannTran::dFdETA(Numeric dETA) const noexcept {
               Numeric dmF0 = -(D0 - 3 * D2 / 2) * dETA;
               Numeric dGD = (dmF0 / (invGD * mF0));
               Numeric dinvGD = -dGD * Constant::pow2(invGD);
               Complex dsqrty =
               ((ETA - 1) * dinvGD + invGD * dETA) /
               (2 * Complex(G2, D2) * Constant::pow2(ETA - 1) * Constant::pow2(invGD));
               Complex ddeltax = -dETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2);
               Complex dx = (-(ETA - 1) * ddeltax + deltax * dETA) /
               (Complex(G2, D2) * Constant::pow2(ETA - 1));
               Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;
               
               switch (calcs) {
                 case CalcType::Full: {
                   Complex dz1 = dsqrtxy - dsqrty;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB =
                   (Constant::sqrt_pi *
                   ((-(Constant::pow2(z1) - 1) * w1 +
                   (Constant::pow2(z2) - 1) * w2) *
                   dsqrty +
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 -
                   (Constant::pow2(z2) - 1) * Complex(0, 1) * dw2 * dz2 +
                   2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) *
                   sqrty) *
                   (ETA - 1) -
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 -
                   Constant::sqrt_pi * (Constant::pow2(z2) - 1) * w2 + 2 * sqrty) *
                   sqrty * dETA) /
                   (2 * Complex(G2, D2) * Constant::pow2(ETA - 1) *
                   Constant::pow2(sqrty));
                   Complex dK = (-FVC + Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA) * dA +
                   Complex(G2, D2) * B * dETA + Complex(G2, D2) * ETA * dB -
                   Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * A * dETA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tLowZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB =
                   -(Constant::sqrt_pi *
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 +
                   2 * w1 * z1 * dz1) -
                   dz1) *
                   invGD -
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 - z1) * dinvGD;
                   Complex dK = (-FVC + Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA) * dA +
                   Complex(G2, D2) * B * dETA + Complex(G2, D2) * ETA * dB -
                   Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * A * dETA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tHighZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK = (-FVC + Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA) * dA +
                   Complex(G2, D2) * B * dETA + Complex(G2, D2) * ETA * dB -
                   Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * A * dETA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowXandHighY: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK = (-FVC + Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA) * dA +
                   Complex(G2, D2) * B * dETA + Complex(G2, D2) * ETA * dB -
                   Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * A * dETA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandLowX: {
                   Complex dz1 = dsqrtxy;
                   Complex dz2 = dx / (2 * sqrtx);
                   Complex dA =
                   2 *
                   (Constant::sqrt_pi * (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) *
                   (ETA - 1) -
                   (Constant::sqrt_pi * w2 * z2 - 1) * dETA) /
                   (Complex(G2, D2) * Constant::pow2(ETA - 1));
                   Complex dB = (-2 * (ETA - 1) *
                   (Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) *
                   (2 * Constant::pow2(sqrty) + x - 1) +
                   Constant::sqrt_pi * w1 * dz1 +
                   Complex(0, Constant::sqrt_pi) * z1 * dw1 * dz1 +
                   (4 * sqrty * dsqrty + dx) *
                   (Constant::sqrt_pi * w2 * z2 - 1)) +
                   (2 * Constant::sqrt_pi * w1 * z1 +
                   2 * (Constant::sqrt_pi * w2 * z2 - 1) *
                   (2 * Constant::pow2(sqrty) + x - 1) -
                   1) *
                   dETA) /
                   (Complex(G2, D2) * Constant::pow2(ETA - 1));
                   Complex dK = (-FVC + Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA) * dA +
                   Complex(G2, D2) * B * dETA + Complex(G2, D2) * ETA * dB -
                   Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * A * dETA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandHighX: {
                   Complex dz1 = dsqrtxy;
                   Complex dA =
                   ((ETA - 1) * (x - 3) * dx + (2 * x - 3) * x * dETA / 2) /
                   (Complex(G2, D2) * Constant::pow2(ETA - 1) * Constant::pow3(x));
                   Complex dB =
                   (-(2 * (-2 * Constant::sqrt_pi * w1 * z1 + 1) * Constant::pow2(x) +
                   (2 * x - 3) * (2 * Constant::pow2(sqrty) + x - 1)) *
                   x * dETA +
                   (ETA - 1) *
                   (-4 * Constant::sqrt_pi *
                   (w1 * dz1 + z1 * Complex(0, 1) * dw1 * dz1) *
                   Constant::pow3(x) +
                   (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x -
                   2 * (x - 3) * (2 * Constant::pow2(sqrty) + x - 1) * dx)) /
                   (2 * Complex(G2, D2) * Constant::pow2(ETA - 1) * Constant::pow3(x));
                   Complex dK = (-FVC + Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA) * dA +
                   Complex(G2, D2) * B * dETA + Complex(G2, D2) * ETA * dB -
                   Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * A * dETA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
               }
               return {};
             }
             
             Complex HartmannTran::dFdH(Numeric dZ) const noexcept {
               Numeric dmF0 = dZ;
               Numeric dGD = (dmF0 / (invGD * mF0));
               Numeric dinvGD = -dGD * Constant::pow2(invGD);
               Complex dsqrty =
               dinvGD / (2 * (ETA - 1) * Complex(G2, D2) * Constant::pow2(invGD));
               Complex ddeltax = Complex(0, dZ);
               Complex dx = -ddeltax / ((ETA - 1) * Complex(G2, D2));
               Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;
               
               switch (calcs) {
                 case CalcType::Full: {
                   Complex dz1 = dsqrtxy - dsqrty;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB =
                   Constant::sqrt_pi *
                   ((-(Constant::pow2(z1) - 1) * w1 + (Constant::pow2(z2) - 1) * w2) *
                   dsqrty +
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 -
                   (Constant::pow2(z2) - 1) * Complex(0, 1) * dw2 * dz2 +
                   2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) *
                   sqrty) /
                   (2 * (ETA - 1) * Complex(G2, D2) * Constant::pow2(sqrty));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tLowZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB =
                   -(Constant::sqrt_pi *
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 +
                   2 * w1 * z1 * dz1) -
                   dz1) *
                   invGD -
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 - z1) * dinvGD;
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tHighZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowXandHighY: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandLowX: {
                   Complex dz1 = dsqrtxy;
                   Complex dz2 = dx / (2 * sqrtx);
                   Complex dA = 2 * Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) /
                   ((ETA - 1) * Complex(G2, D2));
                   Complex dB = -(2 * Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) *
                   (2 * Constant::pow2(sqrty) + x - 1) +
                   2 * Constant::sqrt_pi * w1 * dz1 +
                   Complex(0, 2 * Constant::sqrt_pi) * z1 * dw1 * dz1 +
                   2 * (4 * sqrty * dsqrty + dx) *
                   (Constant::sqrt_pi * w2 * z2 - 1)) /
                   ((ETA - 1) * Complex(G2, D2));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandHighX: {
                   Complex dz1 = dsqrtxy;
                   Complex dA =
                   (x - 3) * dx / ((ETA - 1) * Complex(G2, D2) * Constant::pow3(x));
                   Complex dB = (-2 * Constant::sqrt_pi *
                   (w1 * dz1 + z1 * Complex(0, 1) * dw1 * dz1) *
                   Constant::pow3(x) +
                   (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x / 2 -
                   (x - 3) * (2 * Constant::pow2(sqrty) + x - 1) * dx) /
                   ((ETA - 1) * Complex(G2, D2) * Constant::pow3(x));
                   Complex dK = ETA * Complex(G2, D2) * dB +
                   (ETA * Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) - FVC) * dA;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
               }
               return {};
             }
             
             Complex HartmannTran::dFdVMR(Output d) const noexcept {
               Numeric dmF0 = (1 - ETA) * (d.D0 - 3 * d.D2 / 2) - (D0 - 3 * D2 / 2) * d.ETA;
               Numeric dGD = (dmF0 / (invGD * mF0));
               Numeric dinvGD = -dGD * Constant::pow2(invGD);
               Complex dsqrty = (Complex(G2, D2) * (ETA - 1) * dinvGD +
               Complex(G2, D2) * invGD * d.ETA +
               Complex(d.G2, d.D2) * (ETA - 1) * invGD) /
               (2 * Constant::pow2(Complex(G2, D2)) *
               Constant::pow2(ETA - 1) * Constant::pow2(invGD));
               Complex ddeltax =
               -(ETA - 1) * Complex(d.G0 - 1.5 * d.G2, d.D0 - 1.5 * d.D2) -
               Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * d.ETA + d.FVC;
               Complex dx = (-Complex(G2, D2) * (ETA - 1) * ddeltax +
               Complex(G2, D2) * deltax * d.ETA +
               Complex(d.G2, d.D2) * (ETA - 1) * deltax) /
               (Constant::pow2(Complex(G2, D2)) * Constant::pow2(ETA - 1));
               Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;
               
               switch (calcs) {
                 case CalcType::Full: {
                   Complex dz1 = dsqrtxy - dsqrty;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB =
                   (Constant::sqrt_pi * Complex(G2, D2) *
                   ((-(Constant::pow2(z1) - 1) * w1 +
                   (Constant::pow2(z2) - 1) * w2) *
                   dsqrty +
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 -
                   (Constant::pow2(z2) - 1) * Complex(0, 1) * dw2 * dz2 +
                   2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) *
                   sqrty) *
                   (ETA - 1) -
                   Complex(G2, D2) *
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 -
                   Constant::sqrt_pi * (Constant::pow2(z2) - 1) * w2 +
                   2 * sqrty) *
                   sqrty * d.ETA -
                   Complex(d.G2, d.D2) * (ETA - 1) *
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 -
                   Constant::sqrt_pi * (Constant::pow2(z2) - 1) * w2 +
                   2 * sqrty) *
                   sqrty) /
                   (2 * Constant::pow2(Complex(G2, D2)) * Constant::pow2(ETA - 1) *
                   Constant::pow2(sqrty));
                   Complex dK =
                   Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                   Complex(d.G2, d.D2) * B * ETA +
                   (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                   (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                   Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                   A;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tLowZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB =
                   -(Constant::sqrt_pi *
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 +
                   2 * w1 * z1 * dz1) -
                   dz1) *
                   invGD -
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 - z1) * dinvGD;
                   Complex dK =
                   Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                   Complex(d.G2, d.D2) * B * ETA +
                   (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                   (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                   Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                   A;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tHighZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK =
                   Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                   Complex(d.G2, d.D2) * B * ETA +
                   (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                   (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                   Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                   A;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowXandHighY: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK =
                   Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                   Complex(d.G2, d.D2) * B * ETA +
                   (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                   (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                   Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                   A;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandLowX: {
                   Complex dz1 = dsqrtxy;
                   Complex dz2 = dx / (2 * sqrtx);
                   Complex dA =
                   2 *
                   (Constant::sqrt_pi * Complex(G2, D2) *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) * (ETA - 1) -
                   Complex(G2, D2) * (Constant::sqrt_pi * w2 * z2 - 1) * d.ETA -
                   Complex(d.G2, d.D2) * (Constant::sqrt_pi * w2 * z2 - 1) *
                   (ETA - 1)) /
                   (Constant::pow2(Complex(G2, D2)) * Constant::pow2(ETA - 1));
                   Complex dB =
                   (-2 * Complex(G2, D2) * (ETA - 1) *
                   (Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) *
                   (2 * Constant::pow2(sqrty) + x - 1) +
                   Constant::sqrt_pi * w1 * dz1 +
                   Complex(0, Constant::sqrt_pi) * z1 * dw1 * dz1 +
                   (4 * sqrty * dsqrty + dx) *
                   (Constant::sqrt_pi * w2 * z2 - 1)) +
                   Complex(G2, D2) *
                   (2 * Constant::sqrt_pi * w1 * z1 +
                   2 * (Constant::sqrt_pi * w2 * z2 - 1) *
                   (2 * Constant::pow2(sqrty) + x - 1) -
                   1) *
                   d.ETA +
                   Complex(d.G2, d.D2) * (ETA - 1) *
                   (2 * Constant::sqrt_pi * w1 * z1 +
                   2 * (Constant::sqrt_pi * w2 * z2 - 1) *
                   (2 * Constant::pow2(sqrty) + x - 1) -
                   1)) /
                   (Constant::pow2(Complex(G2, D2)) * Constant::pow2(ETA - 1));
                   Complex dK =
                   Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                   Complex(d.G2, d.D2) * B * ETA +
                   (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                   (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                   Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                   A;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandHighX: {
                   Complex dz1 = dsqrtxy;
                   Complex dA = (2 * Complex(G2, D2) * (ETA - 1) * (x - 3) * dx +
                   Complex(G2, D2) * (2 * x - 3) * x * d.ETA +
                   Complex(d.G2, d.D2) * (ETA - 1) * (2 * x - 3) * x) /
                   (2 * Constant::pow2(Complex(G2, D2)) *
                   Constant::pow2(ETA - 1) * Constant::pow3(x));
                   Complex dB =
                   (-Complex(G2, D2) *
                   (2 * (-2 * Constant::sqrt_pi * w1 * z1 + 1) *
                   Constant::pow2(x) +
                   (2 * x - 3) * (2 * Constant::pow2(sqrty) + x - 1)) *
                   x * d.ETA +
                   Complex(G2, D2) * (ETA - 1) *
                   (-4 * Constant::sqrt_pi *
                   (w1 * dz1 + z1 * Complex(0, 1) * dw1 * dz1) *
                   Constant::pow3(x) +
                   (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x -
                   2 * (x - 3) * (2 * Constant::pow2(sqrty) + x - 1) * dx) -
                   Complex(d.G2, d.D2) *
                   (2 * (-2 * Constant::sqrt_pi * w1 * z1 + 1) *
                   Constant::pow2(x) +
                   (2 * x - 3) * (2 * Constant::pow2(sqrty) + x - 1)) *
                   (ETA - 1) * x) /
                   (2 * Constant::pow2(Complex(G2, D2)) * Constant::pow2(ETA - 1) *
                   Constant::pow3(x));
                   Complex dK =
                   Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                   Complex(d.G2, d.D2) * B * ETA +
                   (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                   (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                   Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                   A;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
               }
               return {};
             }
             
             Complex HartmannTran::dFdT(Output d, Numeric T) const noexcept {
               Numeric dmF0 = (1 - ETA) * (d.D0 - 3 * d.D2 / 2) - (D0 - 3 * D2 / 2) * d.ETA;
               Numeric dGD =
               (dmF0 / (invGD * mF0)) - invGD * invGD / (2 * T * Constant::sqrt_ln_2);
               Numeric dinvGD = -dGD * Constant::pow2(invGD);
               Complex dsqrty = (Complex(G2, D2) * (ETA - 1) * dinvGD +
               Complex(G2, D2) * invGD * d.ETA +
               Complex(d.G2, d.D2) * (ETA - 1) * invGD) /
               (2 * Constant::pow2(Complex(G2, D2)) *
               Constant::pow2(ETA - 1) * Constant::pow2(invGD));
               Complex ddeltax =
               -(ETA - 1) * Complex(d.G0 - 1.5 * d.G2, d.D0 - 1.5 * d.D2) -
               Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * d.ETA + d.FVC;
               Complex dx = (-Complex(G2, D2) * (ETA - 1) * ddeltax +
               Complex(G2, D2) * deltax * d.ETA +
               Complex(d.G2, d.D2) * (ETA - 1) * deltax) /
               (Constant::pow2(Complex(G2, D2)) * Constant::pow2(ETA - 1));
               Complex dsqrtxy = (sqrty * dsqrty + dx / 2) / sqrtxy;
               
               switch (calcs) {
                 case CalcType::Full: {
                   Complex dz1 = dsqrtxy - dsqrty;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB =
                   (Constant::sqrt_pi * Complex(G2, D2) *
                   ((-(Constant::pow2(z1) - 1) * w1 +
                   (Constant::pow2(z2) - 1) * w2) *
                   dsqrty +
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 -
                   (Constant::pow2(z2) - 1) * Complex(0, 1) * dw2 * dz2 +
                   2 * w1 * z1 * dz1 - 2 * w2 * z2 * dz2) *
                   sqrty) *
                   (ETA - 1) -
                   Complex(G2, D2) *
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 -
                   Constant::sqrt_pi * (Constant::pow2(z2) - 1) * w2 +
                   2 * sqrty) *
                   sqrty * d.ETA -
                   Complex(d.G2, d.D2) * (ETA - 1) *
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 -
                   Constant::sqrt_pi * (Constant::pow2(z2) - 1) * w2 +
                   2 * sqrty) *
                   sqrty) /
                   (2 * Constant::pow2(Complex(G2, D2)) * Constant::pow2(ETA - 1) *
                   Constant::pow2(sqrty));
                   Complex dK =
                   Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                   Complex(d.G2, d.D2) * B * ETA +
                   (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                   (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                   Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                   A;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tLowZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB =
                   -(Constant::sqrt_pi *
                   ((Constant::pow2(z1) - 1) * Complex(0, 1) * dw1 * dz1 +
                   2 * w1 * z1 * dz1) -
                   dz1) *
                   invGD -
                   (Constant::sqrt_pi * (Constant::pow2(z1) - 1) * w1 - z1) * dinvGD;
                   Complex dK =
                   Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                   Complex(d.G2, d.D2) * B * ETA +
                   (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                   (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                   Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                   A;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::Noc2tHighZ: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dA =
                   Constant::sqrt_pi * (Complex(0, invGD) * dw1 * dz1 + w1 * dinvGD);
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK =
                   Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                   Complex(d.G2, d.D2) * B * ETA +
                   (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                   (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                   Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                   A;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowXandHighY: {
                   Complex dz1 = deltax * dinvGD + invGD * ddeltax;
                   Complex dz2 = dsqrtxy + dsqrty;
                   Complex dA =
                   Constant::sqrt_pi * ((w1 - w2) * dinvGD +
                   (Complex(0, invGD) * (dw1 * dz1 - dw2 * dz2)));
                   Complex dB = ((4 * Constant::sqrt_pi * w1 * Constant::pow3(z1) +
                   2 * Constant::pow2(z1) - 3) *
                   z1 * dinvGD +
                   (Complex(0, 4 * Constant::sqrt_pi) * Constant::pow4(z1) *
                   dw1 * dz1 -
                   2 * Constant::pow2(z1) * dz1 + 9 * dz1) *
                   invGD) /
                   (4 * Constant::pow4(z1));
                   Complex dK =
                   Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                   Complex(d.G2, d.D2) * B * ETA +
                   (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                   (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                   Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                   A;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandLowX: {
                   Complex dz1 = dsqrtxy;
                   Complex dz2 = dx / (2 * sqrtx);
                   Complex dA =
                   2 *
                   (Constant::sqrt_pi * Complex(G2, D2) *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) * (ETA - 1) -
                   Complex(G2, D2) * (Constant::sqrt_pi * w2 * z2 - 1) * d.ETA -
                   Complex(d.G2, d.D2) * (Constant::sqrt_pi * w2 * z2 - 1) *
                   (ETA - 1)) /
                   (Constant::pow2(Complex(G2, D2)) * Constant::pow2(ETA - 1));
                   Complex dB =
                   (-2 * Complex(G2, D2) * (ETA - 1) *
                   (Constant::sqrt_pi *
                   (w2 * dz2 + z2 * Complex(0, 1) * dw2 * dz2) *
                   (2 * Constant::pow2(sqrty) + x - 1) +
                   Constant::sqrt_pi * w1 * dz1 +
                   Complex(0, Constant::sqrt_pi) * z1 * dw1 * dz1 +
                   (4 * sqrty * dsqrty + dx) *
                   (Constant::sqrt_pi * w2 * z2 - 1)) +
                   Complex(G2, D2) *
                   (2 * Constant::sqrt_pi * w1 * z1 +
                   2 * (Constant::sqrt_pi * w2 * z2 - 1) *
                   (2 * Constant::pow2(sqrty) + x - 1) -
                   1) *
                   d.ETA +
                   Complex(d.G2, d.D2) * (ETA - 1) *
                   (2 * Constant::sqrt_pi * w1 * z1 +
                   2 * (Constant::sqrt_pi * w2 * z2 - 1) *
                   (2 * Constant::pow2(sqrty) + x - 1) -
                   1)) /
                   (Constant::pow2(Complex(G2, D2)) * Constant::pow2(ETA - 1));
                   Complex dK =
                   Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                   Complex(d.G2, d.D2) * B * ETA +
                   (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                   (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                   Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                   A;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
                 }
                 case CalcType::LowYandHighX: {
                   Complex dz1 = dsqrtxy;
                   Complex dA = (2 * Complex(G2, D2) * (ETA - 1) * (x - 3) * dx +
                   Complex(G2, D2) * (2 * x - 3) * x * d.ETA +
                   Complex(d.G2, d.D2) * (ETA - 1) * (2 * x - 3) * x) /
                   (2 * Constant::pow2(Complex(G2, D2)) *
                   Constant::pow2(ETA - 1) * Constant::pow3(x));
                   Complex dB =
                   (-Complex(G2, D2) *
                   (2 * (-2 * Constant::sqrt_pi * w1 * z1 + 1) *
                   Constant::pow2(x) +
                   (2 * x - 3) * (2 * Constant::pow2(sqrty) + x - 1)) *
                   x * d.ETA +
                   Complex(G2, D2) * (ETA - 1) *
                   (-4 * Constant::sqrt_pi *
                   (w1 * dz1 + z1 * Complex(0, 1) * dw1 * dz1) *
                   Constant::pow3(x) +
                   (4 * sqrty * dsqrty + dx) * (2 * x - 3) * x -
                   2 * (x - 3) * (2 * Constant::pow2(sqrty) + x - 1) * dx) -
                   Complex(d.G2, d.D2) *
                   (2 * (-2 * Constant::sqrt_pi * w1 * z1 + 1) *
                   Constant::pow2(x) +
                   (2 * x - 3) * (2 * Constant::pow2(sqrty) + x - 1)) *
                   (ETA - 1) * x) /
                   (2 * Constant::pow2(Complex(G2, D2)) * Constant::pow2(ETA - 1) *
                   Constant::pow3(x));
                   Complex dK =
                   Complex(G2, D2) * B * d.ETA + Complex(G2, D2) * ETA * dB +
                   Complex(d.G2, d.D2) * B * ETA +
                   (Complex(G0 - 1.5 * G2, D0 - 1.5 * D2) * ETA - FVC) * dA +
                   (-Complex(1.5 * G2 - G0, 1.5 * D2 - D0) * d.ETA -
                   Complex(1.5 * d.G2 - d.G0, 1.5 * d.D2 - d.D0) * ETA - d.FVC) *
                   A;
                   return Constant::inv_pi * (-A * dK + K * dA) / Constant::pow2(K);
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
    return CalcType::Noc2tHighZ;  // nb. Value of high/low changes elsewhere
    else if (abs_squared(x) <= 9e-16 * abs_squared(sqrty * sqrty))
      return CalcType::LowXandHighY;
    else if ((abs_squared(sqrty * sqrty) <= 1.e-30 * abs_squared(x)) and
      abs_squared(std::sqrt(x)) <= 16.e6)
      return CalcType::LowYandLowX;  // Weird case, untested
      else if ((abs_squared(sqrty * sqrty) <= 1.e-30 * abs_squared(x)))
        return CalcType::LowYandHighX;
      else
        return CalcType::Full;
}

void HartmannTran::update_calcs() noexcept { calcs = init((1 - ETA) * Complex(G2, D2)); }

void HartmannTran::calc() noexcept {
  switch (calcs) {
    case CalcType::Full:
      z1 = sqrtxy - sqrty;
      z2 = sqrtxy + sqrty;
      w1 = Faddeeva::w(Complex(0, 1) * z1);
      w2 = Faddeeva::w(Complex(0, 1) * z2);
      A = Constant::sqrt_pi * invGD * (w1 - w2);
      B = (-1 +
      Constant::sqrt_pi / (2 * sqrty) * (1 - Constant::pow2(z1)) * w1 -
      Constant::sqrt_pi / (2 * sqrty) * (1 - Constant::pow2(z2)) * w2) /
      ((1 - ETA) * Complex(G2, D2));
      break;
    case CalcType::Noc2tLowZ:
    case CalcType::Noc2tHighZ:
      z1 = deltax * invGD;
      w1 = Faddeeva::w(Complex(0, 1) * z1);
      A = Constant::sqrt_pi * invGD * w1;
      if (abs_squared(z1) < 16e6) {
        calcs = CalcType::Noc2tLowZ;
        B = Constant::sqrt_pi * invGD *
        ((1 - Constant::pow2(z1)) * w1 + z1 / Constant::sqrt_pi);
      } else {
        calcs = CalcType::Noc2tHighZ;
        B = invGD * (Constant::sqrt_pi * w1 + 1 / z1 / 2 -
        3 / Constant::pow3(z1) / 4);
      }
      break;
    case CalcType::LowXandHighY:
      z1 = deltax * invGD;
      z2 = sqrtxy + sqrty;
      w1 = Faddeeva::w(Complex(0, 1) * z1);
      w2 = Faddeeva::w(Complex(0, 1) * z2);
      A = Constant::sqrt_pi * invGD * (w1 - w2);
      B = invGD *
      (Constant::sqrt_pi * w1 + 1 / z1 / 2 - 3 / Constant::pow3(z1) / 4);
      break;
    case CalcType::LowYandLowX:
      sqrtx = std::sqrt(x);
      z1 = sqrtxy;
      z2 = sqrtx;
      w1 = Faddeeva::w(Complex(0, 1) * z1);
      w2 = Faddeeva::w(Complex(0, 1) * z2);
      A = (2 * Constant::sqrt_pi / ((1 - ETA) * Complex(G2, D2))) *
      (Constant::inv_sqrt_pi - z2 * w2);
      B = (1 / ((1 - ETA) * Complex(G2, D2))) *
      (-1 +
      2 * Constant::sqrt_pi * (1 - x - 2 * sqrty * sqrty) *
      (1 / Constant::sqrt_pi - z2 * w2) +
      2 * Constant::sqrt_pi * z1 * w1);
      break;
    case CalcType::LowYandHighX:
      z1 = sqrtxy;
      w1 = Faddeeva::w(Complex(0, 1) * z1);
      A = (1 / ((1 - ETA) * Complex(G2, D2))) *
      (1 / x - 3 / Constant::pow2(x) / 2);
      B = (1 / ((1 - ETA) * Complex(G2, D2))) *
      (-1 +
      (1 - x - 2 * sqrty * sqrty) * (1 / x - 3 / Constant::pow2(x) / 2) +
      2 * Constant::sqrt_pi * z1 * w1);
      break;
  }
  
  dw1 = Complex(0, 2) * (Constant::inv_sqrt_pi - z1 * w1);
  dw2 = Complex(0, 2) * (Constant::inv_sqrt_pi - z2 * w2);
  K = 1 - (FVC - ETA * (Complex(G0, D0) - 3 * Complex(G2, D2) / 2)) * A +
  ETA * Complex(G2, D2) * B;
  F = Constant::inv_pi * A / K;
}

Calculator line_shape_selection(
  const Type type,
  const Numeric F0,
  const Output X,
  const Numeric DC,
  const Numeric DZ) {
  switch (type) {
    case Type::DP:   return Doppler            (F0,    DC, DZ);
    case Type::LP:   return Lorentz            (F0, X        );
    case Type::VP:   return Voigt              (F0, X, DC, DZ);
    case Type::SDVP: return SpeedDependentVoigt(F0, X, DC, DZ);
    case Type::HTP:  return HartmannTran       (F0, X, DC, DZ);
    case Type::FINAL: { /*leave last*/ }
  }
  
  return Noshape{};
}

#define InternalDerivatives(X)                                                                  \
else if (deriv == Jacobian::Line::Shape##X## X0) {                                              \
  if (iz == 0 and iv == 0)                                                                                  \
    line_derivs[ij] =                                                                           \
      band.Line(i).LineShape().d##X## _dX0(T, band.T0(), P, deriv.Target().Position(), vmrs);   \
  dF[ij][iv] += LM * Sz * std::visit([d##X=line_derivs.at(ij)](auto&& LS){return LS.dFd##X(d##X);}, lsm); \
} else if (deriv == Jacobian::Line::Shape##X## X1) {                                            \
  if (iz == 0 and iv == 0)                                                                                  \
    line_derivs[ij] =                                                                           \
      band.Line(i).LineShape().d##X## _dX1(T, band.T0(), P, deriv.Target().Position(), vmrs);   \
  dF[ij][iv] += LM * Sz * std::visit([d##X=line_derivs.at(ij)](auto&& LS){return LS.dFd##X(d##X);}, lsm); \
} else if (deriv == Jacobian::Line::Shape##X## X2) {                                            \
  if (iz == 0 and iv == 0)                                                                                  \
    line_derivs[ij] =                                                                           \
      band.Line(i).LineShape().d##X## _dX2(T, band.T0(), P, deriv.Target().Position(), vmrs);   \
  dF[ij][iv] += LM * Sz * std::visit([d##X=line_derivs.at(ij)](auto&& LS){return LS.dFd##X(d##X);}, lsm); \
} else if (deriv == Jacobian::Line::Shape##X## X3) {                                            \
  if (iz == 0 and iv == 0)                                                                                  \
    line_derivs[ij] =                                                                           \
      band.Line(i).LineShape().d##X## _dX3(T, band.T0(), P, deriv.Target().Position(), vmrs);   \
  dF[ij][iv] += LM * Sz * std::visit([d##X=line_derivs.at(ij)](auto&& LS){return LS.dFd##X(d##X);}, lsm); \
}

#define InternalDerivativesG                                                                  \
else if (deriv == Jacobian::Line::ShapeGX0) {                                              \
  if (iz == 0 and iv == 0)                                                                                  \
    line_derivs[ij] =                                                                           \
    band.Line(i).LineShape().dG_dX0(T, band.T0(), P, deriv.Target().Position(), vmrs);   \
    const Numeric dLM = line_derivs.at(ij); \
    dF[ij][iv] += LM * Sz * std::visit([dG=line_derivs.at(ij)](auto&& LS){return LS.dFdG(dG);}, lsm) + dLM * Sz * std::visit([](auto&& LS){return LS.F;}, lsm); \
} else if (deriv == Jacobian::Line::ShapeGX1) {                                            \
  if (iz == 0 and iv == 0)                                                                                  \
    line_derivs[ij] =                                                                           \
    band.Line(i).LineShape().dG_dX1(T, band.T0(), P, deriv.Target().Position(), vmrs);   \
    const Numeric dLM = line_derivs.at(ij); \
    dF[ij][iv] += LM * Sz * std::visit([dG=line_derivs.at(ij)](auto&& LS){return LS.dFdG(dG);}, lsm) + dLM * Sz * std::visit([](auto&& LS){return LS.F;}, lsm); \
} else if (deriv == Jacobian::Line::ShapeGX2) {                                            \
  if (iz == 0 and iv == 0)                                                                                  \
    line_derivs[ij] =                                                                           \
    band.Line(i).LineShape().dG_dX2(T, band.T0(), P, deriv.Target().Position(), vmrs);   \
    const Numeric dLM = line_derivs.at(ij); \
    dF[ij][iv] += LM * Sz * std::visit([dG=line_derivs.at(ij)](auto&& LS){return LS.dFdG(dG);}, lsm) + dLM * Sz * std::visit([](auto&& LS){return LS.F;}, lsm); \
} else if (deriv == Jacobian::Line::ShapeGX3) {                                            \
  if (iz == 0 and iv == 0)                                                                                  \
    line_derivs[ij] =                                                                           \
    band.Line(i).LineShape().dG_dX3(T, band.T0(), P, deriv.Target().Position(), vmrs);   \
    const Numeric dLM = line_derivs.at(ij); \
    dF[ij][iv] += LM * Sz * std::visit([dG=line_derivs.at(ij)](auto&& LS){return LS.dFdG(dG);}, lsm) + dLM * Sz * std::visit([](auto&& LS){return LS.F;}, lsm); \
}

#define InternalDerivativesY                                                                 \
else if (deriv == Jacobian::Line::ShapeYX0) {                                              \
  if (iz == 0 and iv == 0)                                                                                  \
    line_derivs[ij] =                                                                           \
    band.Line(i).LineShape().dY_dX0(T, band.T0(), P, deriv.Target().Position(), vmrs);   \
    const Complex dLM = Complex(0, -line_derivs.at(ij)); \
    dF[ij][iv] += LM * Sz * std::visit([dY=line_derivs.at(ij)](auto&& LS){return LS.dFdY(dY);}, lsm) + dLM * Sz * std::visit([](auto&& LS){return LS.F;}, lsm); \
} else if (deriv == Jacobian::Line::ShapeYX1) {                                            \
  if (iz == 0 and iv == 0)                                                                                  \
    line_derivs[ij] =                                                                           \
    band.Line(i).LineShape().dY_dX1(T, band.T0(), P, deriv.Target().Position(), vmrs);   \
    const Complex dLM = Complex(0, -line_derivs.at(ij)); \
    dF[ij][iv] += LM * Sz * std::visit([dY=line_derivs.at(ij)](auto&& LS){return LS.dFdY(dY);}, lsm) + dLM * Sz * std::visit([](auto&& LS){return LS.F;}, lsm); \
} else if (deriv == Jacobian::Line::ShapeYX2) {                                            \
  if (iz == 0 and iv == 0)                                                                                  \
    line_derivs[ij] =                                                                           \
    band.Line(i).LineShape().dY_dX2(T, band.T0(), P, deriv.Target().Position(), vmrs);   \
    const Complex dLM = Complex(0, -line_derivs.at(ij)); \
    dF[ij][iv] += LM * Sz * std::visit([dY=line_derivs.at(ij)](auto&& LS){return LS.dFdY(dY);}, lsm) + dLM * Sz * std::visit([](auto&& LS){return LS.F;}, lsm); \
} else if (deriv == Jacobian::Line::ShapeYX3) {                                            \
  if (iz == 0 and iv == 0)                                                                                  \
    line_derivs[ij] =                                                                           \
    band.Line(i).LineShape().dY_dX3(T, band.T0(), P, deriv.Target().Position(), vmrs);   \
    const Complex dLM = Complex(0, -line_derivs.at(ij)); \
    dF[ij][iv] += LM * Sz * std::visit([dY=line_derivs.at(ij)](auto&& LS){return LS.dFdY(dY);}, lsm) + dLM * Sz * std::visit([](auto&& LS){return LS.F;}, lsm); \
}

void compute(ComplexVector& F,
             ArrayOfComplexVector& dF,
             const ConstVectorView& f_grid,
             const AbsorptionLines& band,
             const ArrayOfRetrievalQuantity& jacobian_quantities,
             const Vector& vmrs,
             const Numeric& P,
             const Numeric& T,
             const Numeric& H,
             const bool zeeman,
             const Zeeman::Polarization zeeman_polarization) ARTS_NOEXCEPT {
  // Single line parameter derivatives
  thread_local std::map<Index, Numeric> line_derivs;  // Special storage duration!
  // Multiple line shape model parameter derivatives 
  thread_local std::map<Index, Output > lsmp_derivs;  // Special storage duration!
  
  const Index nj = jacobian_quantities.nelem();
  const Index nl = band.NumLines();
  const Index nv = f_grid.nelem();
  
  const Numeric DC = band.DopplerConstant(T);
  
  // Tests that must be true while calling this function
  ARTS_ASSERT(H >= 0, "Only for positive H.  You provided: ", H)
  ARTS_ASSERT(P >= 0, "Only for positive P.  You provided: ", P)
  ARTS_ASSERT(T > 0, "Only for abs positive T.  You provided: ", T)
  ARTS_ASSERT(band.OK(), "Band is poorly constructed.  You need to use detailed debugger to find out why.")
  ARTS_ASSERT(Index(F.size()) >= nv, "Must have smaller f_grid size than main output.\n"
    "You have: ", nv, "-sized f_grid and ", F.size(), "-sized main output")
  ARTS_ASSERT(Index(dF.size()) >= nj, "Must have smaller jacobian_quantities size than derivative output.\n"
    "You have: ", nj, "-sized jacobian_quantities and ", dF.size(), "-sized derivative output")
  
  // Reset all to zero:
  F = 0;
  dF = F;
  
  const bool do_temperature = do_temperature_jacobian(jacobian_quantities);
  
  if (nl == 0 or (Absorption::relaxationtype_relmat(band.Population()) and band.DoLineMixing(P))) {
    return;  // No line-by-line computations required/wanted
  }
  
  for (Index i=0; i<nl; i++) {
    const Index nz = zeeman ? band.ZeemanCount(i, zeeman_polarization) : 1;
    
    const Output X = band.ShapeParameters(i, T, P, vmrs);
    const Complex LM(1 + X.G, - X.Y);
    if (i == 0) std::cout << LM << '\n';
    
    for (Index iz=0; iz<nz; iz++) {
      const Numeric dfdH = zeeman ? band.ZeemanSplitting(i, zeeman_polarization, iz) : 0;
      const Numeric Sz = zeeman ? band.ZeemanStrength(i, zeeman_polarization, iz) : 1;
      Calculator lsm = line_shape_selection(band.LineShapeType(), band.F0(i), X, DC, dfdH*H);
      
      for (Index iv=0; iv<nv; iv++) {
        F[iv] += LM * Sz * std::visit([f=f_grid[iv]](auto&& LS){return LS(f);}, lsm);
        for (Index ij=0; ij<nj; ij++) {
          if (not propmattype_index(jacobian_quantities, ij)) continue;
          const auto& deriv = jacobian_quantities[ij];
          ARTS_ASSERT(dF[ij].size() == nv, "Must have same freq count as main output for derivatives\n"
            "Should have: ", nv, "Has instead: ", dF[ij].size(), "\nFor derivative type: ", deriv.Target())
          
          if (deriv == Jacobian::Atm::Temperature) {
            if (iz == 0 and iv == 0) lsmp_derivs[ij] = do_temperature ? band.ShapeParameters_dT(i, T, P, vmrs) : Output{};
            const Complex dLM(lsmp_derivs.at(ij).G, - lsmp_derivs.at(ij).Y);
            dF[ij][iv] += Sz *(LM * std::visit([dXdT=lsmp_derivs.at(ij), T](auto&& LS){return LS.dFdT(dXdT, T);}, lsm) +
                              dLM * std::visit([](auto&& LS){return LS.F;}, lsm));
          } else if (is_wind_parameter(deriv)) {
            dF[ij][iv] += LM * Sz * std::visit([](auto&& LS){return LS.dFdf();}, lsm);
          } else if (is_magnetic_parameter(deriv)) {
            dF[ij][iv] += LM * Sz * std::visit([dfdH](auto&& LS){return LS.dFdH(dfdH);}, lsm);
          } else if (deriv.Target().needQuantumIdentity()) {
            if (deriv == Jacobian::Line::VMR) {
              if (iz == 0 and iv == 0) lsmp_derivs[ij] = band.ShapeParameters_dVMR(i, T, P, deriv.QuantumIdentity());
              const Complex dLM(lsmp_derivs.at(ij).G, - lsmp_derivs.at(ij).Y);
              dF[ij][iv] += LM * Sz * std::visit([dXdVMR=lsmp_derivs.at(ij)](auto&& LS){return LS.dFdVMR(dXdVMR);}, lsm) +
                           dLM * Sz * std::visit([](auto&& LS){return LS.F;}, lsm);
            } else {
              const Absorption::QuantumIdentifierLineTarget lt(deriv.Target().QuantumIdentity(), band, i);
              if (lt != Absorption::QuantumIdentifierLineTargetType::Line) continue;
              
              if (deriv == Jacobian::Line::Center) {
                dF[ij][iv] += LM * Sz * std::visit([](auto&& LS){return LS.dFdF0();}, lsm);
              }
              InternalDerivatives(G0)
              InternalDerivatives(D0)
              InternalDerivatives(G2)
              InternalDerivatives(D2)
              InternalDerivatives(ETA)
              InternalDerivatives(FVC)
              InternalDerivativesY  // Special case...
              InternalDerivativesG  // Special case...
              InternalDerivatives(DV)
            }
          }
        }
      }
    }
  }
}

#undef InternalDerivatives
#undef InternalDerivativesG
#undef InternalDerivativesY
}
