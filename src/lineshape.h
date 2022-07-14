#ifndef lineshapes_h
#define lineshapes_h

#include <variant>

#include "constants.h"
#include "energylevelmap.h"
#include "linescaling.h"
#include "nonstd.h"

namespace LineShape {
struct Noshape {
  static constexpr Complex F = Complex(0, 0);

  static constexpr Complex dFdT(const Output &, Numeric) noexcept { return 0; }
  static constexpr Complex dFdf() noexcept { return 0; }
  static constexpr Complex dFdF0() noexcept { return 0; }
  static constexpr Complex dFdH(Numeric) noexcept { return 0; }
  static constexpr Complex dFdFVC(Numeric) noexcept { return 0; }
  static constexpr Complex dFdETA(Numeric) noexcept { return 0; }
  static constexpr Complex dFdVMR(const Output &) noexcept { return 0; }
  static constexpr Complex dFdDV(Numeric) noexcept { return 0; }
  static constexpr Complex dFdD0(Numeric) noexcept { return 0; }
  static constexpr Complex dFdG0(Numeric) noexcept { return 0; }
  static constexpr Complex dFdD2(Numeric) noexcept { return 0; }
  static constexpr Complex dFdG2(Numeric) noexcept { return 0; }

  constexpr Complex operator()(Numeric) const noexcept { return F; }
};  // Noshape

struct Doppler {
  Complex F;
  Numeric x;

  Numeric mF0;
  Numeric invGD;

  constexpr Doppler(Numeric F0_noshift, Numeric DC, Numeric dZ) noexcept
      : F(), x(), mF0(F0_noshift + dZ), invGD(1.0 / nonstd::abs(DC * mF0)) {}

  [[nodiscard]] constexpr Complex dFdT(const Output &,
                                       Numeric T) const noexcept {
    return F * (2 * Constant::pow2(x) - 1) / (2 * T);
  }
  [[nodiscard]] constexpr Complex dFdf() const noexcept {
    return -2 * invGD * F * x;
  }
  [[nodiscard]] constexpr Complex dFdF0() const noexcept {
    return F * (2 * x * (invGD * mF0 + x) - 1) / mF0;
  }
  [[nodiscard]] constexpr Complex dFdH(Numeric dZ) const noexcept {
    return dZ * (F * (2 * x * (invGD * mF0 + x) - 1) / mF0);
  }
  static constexpr Complex dFdFVC(Numeric) noexcept { return 0; }
  static constexpr Complex dFdETA(Numeric) noexcept { return 0; }
  static constexpr Complex dFdVMR(const Output &) noexcept { return 0; }
  static constexpr Complex dFdDV(Numeric) noexcept { return 0; }
  static constexpr Complex dFdD0(Numeric) noexcept { return 0; }
  static constexpr Complex dFdG0(Numeric) noexcept { return 0; }
  static constexpr Complex dFdD2(Numeric) noexcept { return 0; }
  static constexpr Complex dFdG2(Numeric) noexcept { return 0; }

  Complex operator()(Numeric f) noexcept;
};  // Doppler

struct Lorentz {
  Complex F;
  Complex dF;

  Numeric mF0;
  Numeric G0;

  constexpr Lorentz(Numeric F0_noshift, const Output &ls) noexcept
      : mF0(F0_noshift + ls.D0 + ls.DV), G0(ls.G0) {}

  [[nodiscard]] constexpr Complex dFdVMR(const Output &d) const noexcept {
    return Complex(d.G0, d.D0 + d.DV) * dF;
  }
  [[nodiscard]] constexpr Complex dFdT(const Output &d,
                                       Numeric) const noexcept {
    return dFdVMR(d);
  }
  [[nodiscard]] [[nodiscard]] constexpr Complex dFdf() const noexcept {
    return Complex(0, -1) * dF;
  }
  [[nodiscard]] constexpr Complex dFdF0() const noexcept {
    return Complex(0, 1) * dF;
  }
  [[nodiscard]] constexpr Complex dFdDV(Numeric d) const noexcept {
    return d * dFdF0();
  }
  [[nodiscard]] constexpr Complex dFdD0(Numeric d) const noexcept {
    return d * dFdF0();
  }
  static constexpr Complex dFdH(Numeric) noexcept { return 0; }
  [[nodiscard]] constexpr Complex dFdG0(Numeric d) const noexcept {
    return d * dF;
  }
  static constexpr Complex dFdD2(Numeric) noexcept { return 0; }
  static constexpr Complex dFdG2(Numeric) { return 0; }
  static constexpr Complex dFdETA(Numeric) noexcept { return 0; }
  static constexpr Complex dFdFVC(Numeric) noexcept { return 0; }

  constexpr Complex operator()(Numeric f) noexcept {
    F = Constant::inv_pi / Complex(G0, mF0 - f);
    dF = -Constant::pi * Constant::pow2(F);
    return F;
  }
};  // Lorentz

struct Voigt {
  Complex F;
  Complex dF;

  Numeric mF0;
  Numeric invGD;
  Complex z;

  constexpr Voigt(Numeric F0_noshift,
                  const Output &ls,
                  Numeric DC,
                  Numeric dZ) noexcept
      : F(),
        dF(),
        mF0(F0_noshift + dZ + ls.D0 + ls.DV),
        invGD(1.0 / nonstd::abs(DC * mF0)),
        z(invGD * Complex(-mF0, ls.G0)) {}

  [[nodiscard]] constexpr Complex dFdf() const noexcept { return dF; }
  [[nodiscard]] constexpr Complex dFdF0() const noexcept { return -dF; }
  [[nodiscard]] constexpr Complex dFdDV(Numeric d) const noexcept {
    return -d * dF;
  }
  [[nodiscard]] constexpr Complex dFdD0(Numeric d) const noexcept {
    return -d * dF;
  }
  [[nodiscard]] constexpr Complex dFdG0(Numeric d) const noexcept {
    return Complex(0, d) * dF;
  }
  [[nodiscard]] constexpr Complex dFdH(Numeric dZ) const noexcept {
    return -dZ * dF;
  }
  [[nodiscard]] constexpr Complex dFdVMR(const Output &d) const noexcept {
    return Complex(-d.D0 - d.DV, d.G0) * dF;
  }
  [[nodiscard]] constexpr Complex dFdT(const Output &d,
                                       Numeric T) const noexcept {
    return -(F * invGD + dF * z) * (2 * T * (d.D0 + d.DV) + mF0) /
               (2 * T * invGD * mF0) +
           Complex(-d.D0 - d.DV, d.G0) * dF;
  }
  static constexpr Complex dFdETA(Numeric) noexcept { return 0; }
  static constexpr Complex dFdFVC(Numeric) noexcept { return 0; }
  static constexpr Complex dFdD2(Numeric) noexcept { return 0; }
  static constexpr Complex dFdG2(Numeric) { return 0; }

  Complex operator()(Numeric f) noexcept;

  [[nodiscard]] bool OK() const noexcept { return invGD > 0; }
};  // Voigt

struct SpeedDependentVoigt {
  enum struct CalcType : char {
    Voigt,
    LowXandHighY,
    LowYandLowX,
    LowYandHighX,
    Full
  };

  Complex F;

  Numeric mF0;
  Numeric invGD;
  Complex invc2;
  Complex dx;
  Complex x;
  Complex sqrty;
  CalcType calcs;
  Complex sq;
  Complex w1;
  Complex w2;
  Complex dw1;
  Complex dw2;

  SpeedDependentVoigt(Numeric F0_noshift,
                      const Output &ls,
                      Numeric GD_div_F0,
                      Numeric dZ) noexcept;

  [[nodiscard]] Complex dFdf() const noexcept;
  [[nodiscard]] Complex dFdF0() const noexcept;
  [[nodiscard]] Complex dFdD0(Numeric dD0dD0) const noexcept;
  [[nodiscard]] Complex dFdG0(Numeric dG0dG0) const noexcept;
  [[nodiscard]] Complex dFdD2(Numeric dD2dD2) const noexcept;
  [[nodiscard]] Complex dFdG2(Numeric dG2dG2) const noexcept;
  [[nodiscard]] Complex dFdH(Numeric dZ) const noexcept;
  [[nodiscard]] Complex dFdVMR(const Output &d) const noexcept;
  [[nodiscard]] Complex dFdT(const Output &d, Numeric T) const noexcept;
  static constexpr Complex dFdETA(Numeric) noexcept { return 0; }
  static constexpr Complex dFdFVC(Numeric) noexcept { return 0; }
  static constexpr Complex dFdDV(Numeric) noexcept { return 0; }

  Complex operator()(Numeric f) noexcept;

  [[nodiscard]] CalcType init(const Complex c2) const noexcept;
  void update_calcs() noexcept;
  void calc() noexcept;
};  // SpeedDependentVoigt

struct HartmannTran {
  enum struct CalcType : char {
    Noc2tLowZ,
    Noc2tHighZ,
    LowXandHighY,
    LowYandLowX,
    LowYandHighX,
    Full
  };

  Complex F;

  Numeric G0;
  Numeric D0;
  Numeric G2;
  Numeric D2;
  Numeric FVC;
  Numeric ETA;
  Numeric mF0;
  Numeric invGD;
  Complex deltax;
  Complex sqrty;

  CalcType calcs;

  Complex x;
  Complex sqrtxy;
  Complex sqrtx;
  Complex z1;
  Complex z2;
  Complex w1;
  Complex w2;
  Complex A;
  Complex B;
  Complex K;
  Complex dw1;
  Complex dw2;

  HartmannTran(Numeric F0_noshift,
               const Output &ls,
               Numeric GD_div_F0,
               Numeric dZ) noexcept;

  [[nodiscard]] Complex dFdf() const noexcept;
  [[nodiscard]] Complex dFdF0() const noexcept;
  [[nodiscard]] Complex dFdD0(Numeric dD0) const noexcept;
  [[nodiscard]] Complex dFdG0(Numeric dG0) const noexcept;
  [[nodiscard]] Complex dFdD2(Numeric dD2) const noexcept;
  [[nodiscard]] Complex dFdG2(Numeric dG2) const noexcept;
  [[nodiscard]] Complex dFdFVC(Numeric dFVC) const noexcept;
  [[nodiscard]] Complex dFdETA(Numeric dETA) const noexcept;
  [[nodiscard]] Complex dFdH(Numeric dZ) const noexcept;
  [[nodiscard]] Complex dFdVMR(const Output &d) const noexcept;
  [[nodiscard]] Complex dFdT(const Output &d, Numeric T) const noexcept;
  static constexpr Complex dFdDV(Numeric) noexcept { return 0; }

  Complex operator()(Numeric f) noexcept;

  [[nodiscard]] CalcType init(const Complex c2t) const noexcept;
  void update_calcs() noexcept;
  void calc() noexcept;
};  // HartmannTran

struct Nonorm {
  static constexpr Numeric N = 1.0;

  constexpr Nonorm() noexcept = default;

  static constexpr Numeric dNdT(Numeric, Numeric) noexcept { return 0; }
  static constexpr Numeric dNdf(Numeric) noexcept { return 0; }
  [[nodiscard]] constexpr Numeric dNdF0() const noexcept { return 0; }

  constexpr Numeric operator()(Numeric) noexcept { return N; }
};  // Nonorm

struct VanVleckHuber {
  Numeric N;

  Numeric c1;
  Numeric tanh_c1f0;
  Numeric inv_denom;
  Numeric tanh_c1f;

  VanVleckHuber(Numeric F0, Numeric T) noexcept;

  [[nodiscard]] Numeric dNdT(Numeric T, Numeric f) const noexcept;
  [[nodiscard]] Numeric dNdf(Numeric f) const noexcept;
  [[nodiscard]] Numeric dNdF0() const noexcept;

  Numeric operator()(Numeric f) noexcept;
};  // VanVleckHuber

struct VanVleckWeisskopf {
  Numeric N;

  Numeric invF0;

  constexpr VanVleckWeisskopf(Numeric F0) noexcept : N(1), invF0(1.0 / F0) {}

  static constexpr Numeric dNdT(Numeric, Numeric) noexcept { return 0; }
  [[nodiscard]] constexpr Numeric dNdf(Numeric f) const noexcept {
    return 2.0 * f * Constant::pow2(invF0);
  }
  [[nodiscard]] constexpr Numeric dNdF0() const noexcept {
    return -2.0 * N * invF0;
  }

  constexpr Numeric operator()(Numeric f) noexcept {
    N = Constant::pow2(f * invF0);
    return N;
  }
};  // VanVleckWeisskopf

struct RosenkranzQuadratic {
  Numeric N;

  Numeric fac;
  Numeric dfacdT;
  Numeric dfacdF0;

  RosenkranzQuadratic(Numeric F0, Numeric T) noexcept;

  [[nodiscard]] Numeric dNdT(Numeric, Numeric f) const noexcept;
  [[nodiscard]] Numeric dNdf(Numeric f) const noexcept;
  [[nodiscard]] Numeric dNdF0() const noexcept;

  Numeric operator()(Numeric f) noexcept;
};  // RosenkranzQuadratic

struct SimpleFrequencyScaling {
  Numeric N;

  Numeric T;
  Numeric F0;
  Numeric expF0;
  Numeric expm1F0;

  constexpr SimpleFrequencyScaling(Numeric exp,
                                   Numeric expm1,
                                   Numeric F0_,
                                   Numeric T_)
      : N(1.0), T(T_), F0(F0_), expF0(exp), expm1F0(expm1) {}

  SimpleFrequencyScaling(Numeric F0_, Numeric T_) noexcept
      : SimpleFrequencyScaling(
            std::exp(-(Constant::h * F0_) / (Constant::k * T_)),
            std::expm1(-(Constant::h * F0_) / (Constant::k * T_)),
            F0_,
            T_) {}

  [[nodiscard]] Numeric dNdT(Numeric t_ [[maybe_unused]],
                             Numeric f) const ARTS_NOEXCEPT;
  [[nodiscard]] Numeric dNdf(Numeric f) const noexcept;
  [[nodiscard]] constexpr Numeric dNdF0() const noexcept {
    return -N / F0 + N * Constant::h * expF0 / (Constant::k * T * expm1F0);
  }

  Numeric operator()(Numeric f) noexcept;
};  // SimpleFrequencyScaling

struct Nostrength {
  static constexpr Numeric S = 1.0;
  static constexpr Numeric N = 0.0;

  static constexpr Numeric dSdT() noexcept { return 0; }
  static constexpr Numeric dSdI0() noexcept { return 0; }
  static constexpr Numeric dSdF0() noexcept { return 0; }
  static constexpr Numeric dSdNLTEu() noexcept { return 0; }
  static constexpr Numeric dSdNLTEl() noexcept { return 0; }
  static constexpr Numeric dSdSELFVMR() noexcept { return 0; }

  static constexpr Numeric dNdT() noexcept { return 0; }
  static constexpr Numeric dNdI0() noexcept { return 0; }
  static constexpr Numeric dNdF0() noexcept { return 0; }
  static constexpr Numeric dNdNLTEu() noexcept { return 0; }
  static constexpr Numeric dNdNLTEl() noexcept { return 0; }
  static constexpr Numeric dNdSELFVMR() noexcept { return 0; }

  static constexpr bool do_nlte() noexcept { return false; }
};  // Nostrength

struct LocalThermodynamicEquilibrium {
  Numeric S;
  static constexpr Numeric N = 0.0;

  Numeric dSdI0val;
  Numeric dSdTval;
  Numeric dSdF0val;
  Numeric dSdSELFVMRval;

  constexpr LocalThermodynamicEquilibrium(Numeric I0,
                                          Numeric r,
                                          Numeric drdSELFVMR,
                                          Numeric drdT,
                                          Numeric QT0,
                                          Numeric QT,
                                          Numeric dQTdT,
                                          Numeric br,
                                          Numeric dbr_dT_rat,
                                          Numeric stim,
                                          Numeric dstim_dT,
                                          Numeric dstim_dF0) noexcept
      : S(),
        dSdI0val(r * br * stim * QT0 / QT),
        dSdTval(I0 * (r * br * dstim_dT * QT0 / QT +
                      dSdI0val * (dbr_dT_rat - dQTdT / QT) +
                      drdT * br * stim * QT0 / QT)),
        dSdF0val(r * I0 * br * dstim_dF0 * QT0 / QT),
        dSdSELFVMRval(drdSELFVMR * I0 * br * stim * QT0 / QT) {
    S = I0 * dSdI0val;
  }

  LocalThermodynamicEquilibrium(Numeric I0,
                                Numeric T0,
                                Numeric T,
                                Numeric F0,
                                Numeric E0,
                                Numeric QT,
                                Numeric QT0,
                                Numeric dQTdT,
                                Numeric r,
                                Numeric drdSELFVMR,
                                Numeric drdT) noexcept;

  [[nodiscard]] constexpr Numeric dSdT() const noexcept { return dSdTval; }
  [[nodiscard]] constexpr Numeric dSdI0() const noexcept { return dSdI0val; }
  [[nodiscard]] constexpr Numeric dSdF0() const noexcept { return dSdF0val; }
  static constexpr Numeric dSdNLTEu() noexcept { return 0; }
  static constexpr Numeric dSdNLTEl() noexcept { return 0; }
  [[nodiscard]] constexpr Numeric dSdSELFVMR() const noexcept {
    return dSdSELFVMRval;
  }

  static constexpr Numeric dNdT() noexcept { return 0; }
  static constexpr Numeric dNdI0() noexcept { return 0; }
  static constexpr Numeric dNdF0() noexcept { return 0; }
  static constexpr Numeric dNdNLTEu() noexcept { return 0; }
  static constexpr Numeric dNdNLTEl() noexcept { return 0; }
  static constexpr Numeric dNdSELFVMR() noexcept { return 0; }

  static constexpr bool do_nlte() noexcept { return false; }
};  // LocalThermodynamicEquilibrium

struct FullNonLocalThermodynamicEquilibrium {
  Numeric S;
  Numeric N;

  static constexpr Numeric c0 = 2.0 * Constant::h / Constant::pow2(Constant::c);
  static constexpr Numeric c1 = Constant::h / (4 * Constant::pi);

  Numeric dSdTval;
  Numeric dNdTval;

  Numeric dSdF0val;
  Numeric dNdF0val;

  Numeric dSdr1;
  Numeric dSdr2;
  Numeric dNdr2;

  Numeric dSdSELFVMRval;
  Numeric dNdSELFVMRval;

  constexpr FullNonLocalThermodynamicEquilibrium(Numeric r,
                                                 Numeric drdSELFVMR,
                                                 Numeric drdt,
                                                 Numeric k,
                                                 Numeric dkdF0,
                                                 Numeric dkdr1,
                                                 Numeric dkdr2,
                                                 Numeric e,
                                                 Numeric dedF0,
                                                 Numeric dedr2,
                                                 Numeric B,
                                                 Numeric dBdT,
                                                 Numeric dBdF0) noexcept
      : S(),
        N(),
        dSdTval(drdt * k),
        dNdTval(drdt * (e - k * B) - r * k * dBdT),
        dSdF0val(r * dkdF0),
        dNdF0val(r * (dedF0 - dkdF0 * B - k * dBdF0)),
        dSdr1(r * dkdr1),
        dSdr2(r * dkdr2),
        dNdr2(r * (dedr2 - dkdr2 * B)),
        dSdSELFVMRval(drdSELFVMR * k),
        dNdSELFVMRval(drdSELFVMR * (e - k * B)) {
    S = r * k;
    N = r * (e - k * B);
  }

  FullNonLocalThermodynamicEquilibrium(Numeric F0,
                                       Numeric A21,
                                       Numeric T,
                                       Numeric g1,
                                       Numeric g2,
                                       Numeric r1,
                                       Numeric r2,
                                       Numeric r,
                                       Numeric drdSELFVMR,
                                       Numeric drdT) noexcept;

  [[nodiscard]] constexpr Numeric dSdT() const noexcept { return 0; }
  static constexpr Numeric dSdI0() noexcept { return 0; }
  [[nodiscard]] constexpr Numeric dSdF0() const noexcept { return dSdF0val; }
  [[nodiscard]] constexpr Numeric dSdNLTEu() const noexcept { return dSdr1; }
  [[nodiscard]] constexpr Numeric dSdNLTEl() const noexcept { return dSdr2; }
  [[nodiscard]] constexpr Numeric dSdSELFVMR() const noexcept {
    return dSdSELFVMRval;
  }

  [[nodiscard]] constexpr Numeric dNdT() const noexcept { return dNdTval; }
  static constexpr Numeric dNdI0() noexcept { return 0; }
  [[nodiscard]] constexpr Numeric dNdF0() const noexcept { return dNdF0val; }
  [[nodiscard]] constexpr Numeric dNdNLTEu() const noexcept { return -dSdr1; }
  [[nodiscard]] constexpr Numeric dNdNLTEl() const noexcept { return dNdr2; }
  [[nodiscard]] constexpr Numeric dNdSELFVMR() const noexcept {
    return dNdSELFVMRval;
  }

  static constexpr bool do_nlte() noexcept { return true; }
};  // FullNonLocalThermodynamicEquilibrium

struct VibrationalTemperaturesNonLocalThermodynamicEquilibrium {
  Numeric S;
  Numeric N;

  Numeric dSdI0val;
  Numeric dNdI0val;

  Numeric dSdTval;
  Numeric dNdTval;

  Numeric dSdF0val;
  Numeric dNdF0val;

  Numeric dSdTl;
  Numeric dSdTu;
  Numeric dNdTu;

  Numeric dSdSELFVMRval;
  Numeric dNdSELFVMRval;

  constexpr VibrationalTemperaturesNonLocalThermodynamicEquilibrium(
      Numeric I0,
      Numeric QT0,
      Numeric QT,
      Numeric dQTdT,
      Numeric r,
      Numeric drdSELFVMR,
      Numeric drdT,
      Numeric K1,
      Numeric dK1dT,
      Numeric K2,
      Numeric dK2dT,
      Numeric dK2dF0,
      Numeric K3,
      Numeric dK3dT,
      Numeric dK3dF0,
      Numeric dK3dTl,
      Numeric dK3dTu,
      Numeric K4,
      Numeric dK4dT,
      Numeric dK4dTu,
      Numeric B,
      Numeric dBdT,
      Numeric dBdF0) noexcept
      : S(),
        N(),
        dSdI0val(r * QT0 / QT * K1 * K2 * K3),
        dNdI0val(B * r * QT0 / QT * K1 * K2 * (K4 - K3)),
        dSdTval(I0 * (drdT * QT0 / QT * K1 * K2 * K3 -
                      r * dQTdT * QT0 / Constant::pow2(QT) * K1 * K2 * K3 +
                      r * QT0 / QT * dK1dT * K2 * K3 +
                      r * QT0 / QT * K1 * dK2dT * K3 +
                      r * QT0 / QT * K1 * K2 * dK3dT)),
        dNdTval(I0 * (dBdT * r * QT0 / QT * K1 * K2 * (K4 - K3) +
                      B * drdT * QT0 / QT * K1 * K2 * (K4 - K3) -
                      B * r * dQTdT * QT0 / Constant::pow2(QT) * K1 * K2 *
                          (K4 - K3) +
                      B * r * QT0 / QT * dK1dT * K2 * (K4 - K3) +
                      B * r * QT0 / QT * K1 * dK2dT * (K4 - K3) +
                      B * r * QT0 / QT * K1 * K2 * (dK4dT - dK3dT))),
        dSdF0val(I0 * (r * QT0 / QT * K1 * dK2dF0 * K3 +
                       r * QT0 / QT * K1 * K2 * dK3dF0)),
        dNdF0val(I0 * (dBdF0 * r * QT0 / QT * K1 * K2 * (K4 - K3) +
                       B * r * QT0 / QT * K1 * dK2dF0 * (K4 - K3) -
                       B * r * QT0 / QT * K1 * K2 * dK3dF0)),
        dSdTl(I0 * r * QT0 / QT * K1 * K2 * dK3dTl),
        dSdTu(I0 * r * QT0 / QT * K1 * K2 * dK3dTu),
        dNdTu(B * r * QT0 / QT * K1 * K2 * (dK4dTu - dK3dTu)),
        dSdSELFVMRval(I0 * drdSELFVMR * QT0 / QT * K1 * K2 * K3),
        dNdSELFVMRval(I0 * B * drdSELFVMR * QT0 / QT * K1 * K2 * (K4 - K3)) {
    S = I0 * dSdI0val;
    N = I0 * dNdI0val;
  }

  VibrationalTemperaturesNonLocalThermodynamicEquilibrium(
      Numeric I0,
      Numeric T0,
      Numeric T,
      Numeric Tl,
      Numeric Tu,
      Numeric F0,
      Numeric E0,
      Numeric Evl,
      Numeric Evu,
      Numeric QT,
      Numeric QT0,
      Numeric dQTdT,
      Numeric r,
      Numeric drdSELFVMR,
      Numeric drdT) noexcept;

  [[nodiscard]] constexpr Numeric dSdT() const noexcept { return dSdTval; }
  [[nodiscard]] constexpr Numeric dSdI0() const noexcept { return dSdI0val; }
  [[nodiscard]] constexpr Numeric dSdF0() const noexcept { return dSdF0val; }
  [[nodiscard]] constexpr Numeric dSdNLTEu() const noexcept { return dSdTl; }
  [[nodiscard]] constexpr Numeric dSdNLTEl() const noexcept { return dSdTu; }
  [[nodiscard]] constexpr Numeric dSdSELFVMR() const noexcept {
    return dSdSELFVMRval;
  }

  [[nodiscard]] constexpr Numeric dNdT() const noexcept { return dNdTval; }
  [[nodiscard]] constexpr Numeric dNdI0() const noexcept { return dNdI0val; }
  [[nodiscard]] constexpr Numeric dNdF0() const noexcept { return dNdF0val; }
  [[nodiscard]] constexpr Numeric dNdNLTEu() const noexcept { return -dSdTl; }
  [[nodiscard]] constexpr Numeric dNdNLTEl() const noexcept { return dNdTu; }
  [[nodiscard]] constexpr Numeric dNdSELFVMR() const noexcept {
    return dNdSELFVMRval;
  }

  static constexpr bool do_nlte() noexcept { return true; }
};  // VibrationalTemperaturesNonLocalThermodynamicEquilibrium

//! Line shape calculator.
class Calculator {
  using Variant = std::variant<Noshape,
                               Doppler,
                               Lorentz,
                               Voigt,
                               SpeedDependentVoigt,
                               HartmannTran>;
  Variant ls;

 public:
  [[nodiscard]] Complex dFdT(const Output &dXdT, Numeric T) const noexcept;

  [[nodiscard]] Complex dFdf() const noexcept;

  [[nodiscard]] Complex dFdF0() const noexcept;

  [[nodiscard]] Complex dFdH(Numeric dfdH) const noexcept;

  [[nodiscard]] Complex dFdVMR(const Output &dXdVMR) const noexcept;

  [[nodiscard]] Complex dFdFVC(Numeric d) const noexcept;

  [[nodiscard]] Complex dFdETA(Numeric d) const noexcept;

  [[nodiscard]] Complex dFdDV(Numeric d) const noexcept;

  [[nodiscard]] Complex dFdD0(Numeric d) const noexcept;

  [[nodiscard]] Complex dFdG0(Numeric d) const noexcept;

  [[nodiscard]] Complex dFdD2(Numeric d) const noexcept;

  [[nodiscard]] Complex dFdG2(Numeric d) const noexcept;

  [[nodiscard]] Complex F() const noexcept;

  //! Call operator on frequency.  Must call this before any of the derivatives
  Complex operator()(Numeric f) noexcept;

  Calculator(const Type type,
             const Numeric F0,
             const Output &X,
             const Numeric DC,
             const Numeric DZ,
             bool manually_mirrored) noexcept;

  Calculator(const Absorption::MirroringType mirror,
             const Type type,
             const Numeric F0,
             const Output &X,
             const Numeric DC,
             const Numeric DZ);
};  // Calculator

class Normalizer {
  using Variant = std::variant<Nonorm,
                               VanVleckHuber,
                               VanVleckWeisskopf,
                               RosenkranzQuadratic,
                               SimpleFrequencyScaling>;
  Variant ls_norm;

 public:
  [[nodiscard]] Numeric dNdT(Numeric T, Numeric f) const noexcept;

  [[nodiscard]] Numeric dNdf(Numeric f) const noexcept;

  [[nodiscard]] Numeric dNdF0() const noexcept;

  [[nodiscard]] Numeric operator()(Numeric f) noexcept;

  Normalizer(const Absorption::NormalizationType type,
             const Numeric F0,
             const Numeric T) noexcept;
};  // Normalizer

/** Class encapsulating all supported types of intensity calculations of individual absorption lines */
class IntensityCalculator {
  using Variant =
      std::variant<Nostrength,
                   LocalThermodynamicEquilibrium,
                   FullNonLocalThermodynamicEquilibrium,
                   VibrationalTemperaturesNonLocalThermodynamicEquilibrium>;
  Variant ls_str;

 public:
  /** The line strength absorption */
  [[nodiscard]] Numeric S() const noexcept;

  /** The line strength absorption derivative wrt temperature */
  [[nodiscard]] Numeric dSdT() const noexcept;

  /** The line strength absorption derivative wrt the reference line strength */
  [[nodiscard]] Numeric dSdI0() const noexcept;

  /** The line strength absorption derivative wrt the line center */
  [[nodiscard]] Numeric dSdF0() const noexcept;

  /** The line strength absorption derivative wrt either the upper state number density distribution or its vibration temperature */
  [[nodiscard]] Numeric dSdNLTEu() const noexcept;

  /** The line strength absorption derivative wrt either the lower state number density distribution or its vibration temperature */
  [[nodiscard]] Numeric dSdNLTEl() const noexcept;

  /** The line strength derivative wrt the VMR of the band's species */
  [[nodiscard]] Numeric dSdSELFVMR() const noexcept;

  /** The line strength source offset */
  [[nodiscard]] Numeric N() const noexcept;

  /** The line strength source offset derivative wrt temperature */
  [[nodiscard]] Numeric dNdT() const noexcept;

  /** The line strength source offset derivative wrt the reference line strength */
  [[nodiscard]] Numeric dNdI0() const noexcept;

  /** The line strength source offset derivative wrt the line center */
  [[nodiscard]] Numeric dNdF0() const noexcept;

  /** The line strength source offset derivative wrt either the upper state number density distribution or its vibration temperature */
  [[nodiscard]] Numeric dNdNLTEu() const noexcept;

  /** The line strength source offset derivative wrt either the lower state number density distribution or its vibration temperature */
  [[nodiscard]] Numeric dNdNLTEl() const noexcept;

  /** The line source offset derivative wrt the VMR of the band's species */
  [[nodiscard]] Numeric dNdSELFVMR() const noexcept;

  /** Whether or not NLTE is possible with the selected intensity variant */
  [[nodiscard]] constexpr bool do_nlte() const noexcept {
    return std::visit([](auto &&S) { return S.do_nlte(); }, ls_str);
  }

  IntensityCalculator(const Numeric T,
                      const Numeric QT,
                      const Numeric QT0,
                      const Numeric dQTdT,
                      const Numeric r,
                      const Numeric drdSELFVMR,
                      const Numeric drdT,
                      const EnergyLevelMap &nlte,
                      const Absorption::Lines &band,
                      const Index line_index) noexcept;
};  // IntensityCalculator

/** Main computational data for the line shape and strength calculations */
struct ComputeData {
  ComplexVector F, N;
  ComplexMatrix dF, dN;
  const Vector &f_grid;
  const bool do_nlte;

  ComputeData(const Vector &f,
              const ArrayOfRetrievalQuantity &jacobian_quantities,
              const bool nlte) noexcept
      : F(f.nelem(), 0),
        N(nlte ? f.nelem() : 0, 0),
        dF(f.nelem(), jacobian_quantities.nelem(), 0),
        dN(nlte ? f.nelem() : 0, nlte ? jacobian_quantities.nelem() : 0, 0),
        f_grid(f),
        do_nlte(nlte) {}

  void reset() noexcept {
    F = 0;
    N = 0;
    dF = 0;
    dN = 0;
  }

  /** Add a sparse grid to this grid via linear interpolation
    *
    * @param[in] sparse The sparsely gridded data
  */
  void interp_add_even(const ComputeData &sparse) ARTS_NOEXCEPT;

  /** Add a sparse grid to this grid via square interpolation
    *
    * @param[in] sparse The sparsely gridded data
  */
  void interp_add_triplequad(const ComputeData &sparse) ARTS_NOEXCEPT;

  /** All four fields are set to zero at i if F[i].real() < 0 */
  void enforce_positive_absorption() noexcept {
    const Index nf = f_grid.nelem();
    for (Index i = 0; i < nf; i++) {
      if (F[i].real() < 0) {
        F[i] = 0;
        dF(i, joker) = 0;
        if (do_nlte) {
          N[i] = 0;
          dN(i, joker) = 0;
        }
      }
    }
  }

  /** Adds two identical compute data fields, the size must be identical */
  ComputeData &operator+=(const ComputeData &other) {
    F += other.F;
    N += other.N;
    dF += other.dF;
    dN += other.dN;
    return *this;
  }
};

/** Compute the absorption of an absorption band
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
 *
 * @param[inout] com Main computations variable.  Should be initialized and may have been used before.
 * @param[inout] sparse_com Sparse computations variable.  Should be initialized and may have been used before.
 * @param[in] band The absorption band
 * @param[in] jacobian_quantities As WSV
 * @param[in] rtp_nlte As WSV
 * @param[in] vmrs The volume mixing ratios of the band's line shape model
 * @param[in] self_tag The species tag from abs_species this band belongs to (only used for derivatives)
 * @param[in] self_vmr The volume mixing of the band's species.
 * @param[in] isot_ratio The sotopologue ratio of the band's species
 * @param[in] rtp_pressure As WSV
 * @param[in] rtp_temperature As WSV
 * @param[in] H The magnetic field strength in Teslas
 * @param[in] sparse_lim The frequency separating the sparse and dense frequency grid calculations
 * @param[in] zeeman_polarization Type of Zeeman polarization
 * @param[in] speedup_type Type of sparse grid interactions
 * @param[in] robust If true, a band with line mixing parameters guarantees non-negative output by allocating its own com and sparse_com for local calculations
 */
void compute(ComputeData &com,
             ComputeData &sparse_com,
             const AbsorptionLines &band,
             const ArrayOfRetrievalQuantity &jacobian_quantities,
             const EnergyLevelMap &rtp_nlte,
             const Vector &vmrs,
             const ArrayOfSpeciesTag &self_tag,
             const Numeric &self_vmr,
             const Numeric &isot_ratio,
             const Numeric &rtp_pressure,
             const Numeric &rtp_temperature,
             const Numeric &H,
             const Numeric &sparse_lim,
             const Zeeman::Polarization zeeman_polarization,
             const Options::LblSpeedup speedup_type,
             const bool robust) ARTS_NOEXCEPT;

Vector linear_sparse_f_grid(const Vector &f_grid,
                            const Numeric &sparse_df) ARTS_NOEXCEPT;

bool good_linear_sparse_f_grid(const Vector &f_grid_dense,
                               const Vector &f_grid_sparse) noexcept;

Vector triple_sparse_f_grid(const Vector &f_grid,
                            const Numeric &sparse_df) noexcept;

}  // namespace LineShape

#endif  // lineshapes_h
