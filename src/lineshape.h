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
}; // Noshape

struct Doppler {
  Complex F;
  Numeric x;
  
  Numeric mF0;
  Numeric invGD;

  constexpr Doppler(Numeric F0_noshift, Numeric DC, Numeric dZ) noexcept
  : F(), x(), mF0(F0_noshift + dZ), invGD(1.0 / nonstd::abs(DC * mF0)) {}

  constexpr Complex dFdT(const Output &, Numeric T) const noexcept {
    return F * (2 * Constant::pow2(x) - 1) / (2 * T);
  }
  constexpr Complex dFdf() const noexcept { return -2 * invGD * F * x; }
  constexpr Complex dFdF0() const noexcept {
    return F * (2 * x * (invGD * mF0 + x) - 1) / mF0;
  }
  constexpr Complex dFdH(Numeric dZ) const noexcept {
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
}; // Doppler

struct Lorentz {
  Complex F;
  Complex dF;
  
  Numeric mF0;
  Numeric G0;

  constexpr Lorentz(Numeric F0_noshift, const Output &ls) noexcept
      : mF0(F0_noshift + ls.D0 + ls.DV), G0(ls.G0) {}

  constexpr Complex dFdVMR(const Output &d) const noexcept {
    return Complex(d.G0, d.D0 + d.DV) * dF;
  }
  constexpr Complex dFdT(const Output &d, Numeric) const noexcept {
    return dFdVMR(d);
  }
  constexpr Complex dFdf() const noexcept { return Complex(0, -1) * dF; }
  constexpr Complex dFdF0() const noexcept { return Complex(0, 1) * dF; }
  constexpr Complex dFdDV(Numeric d) const noexcept { return d * dFdF0(); }
  constexpr Complex dFdD0(Numeric d) const noexcept { return d * dFdF0(); }
  static constexpr Complex dFdH(Numeric) noexcept { return 0; }
  constexpr Complex dFdG0(Numeric d) const noexcept { return d * dF; }
  static constexpr Complex dFdD2(Numeric) noexcept { return 0; }
  static constexpr Complex dFdG2(Numeric) { return 0; }
  static constexpr Complex dFdETA(Numeric) noexcept { return 0; }
  static constexpr Complex dFdFVC(Numeric) noexcept { return 0; }

  constexpr Complex operator()(Numeric f) noexcept {
    F = Constant::inv_pi / Complex(G0, mF0 - f);
    dF = -Constant::pi * Constant::pow2(F);
    return F;
  }
}; // Lorentz

struct Voigt {
  Complex F;
  Complex dF;
  
  Numeric mF0;
  Numeric invGD;
  Complex z;

  constexpr Voigt(Numeric F0_noshift, const Output &ls, Numeric DC, Numeric dZ) noexcept
  : F(), dF(), mF0(F0_noshift + dZ + ls.D0 + ls.DV),
  invGD(1.0 / nonstd::abs(DC * mF0)), z(invGD * Complex(-mF0, ls.G0)) {}

  constexpr Complex dFdf() const noexcept { return dF; }
  constexpr Complex dFdF0() const noexcept { return -dF; }
  constexpr Complex dFdDV(Numeric d) const noexcept { return -d * dF; }
  constexpr Complex dFdD0(Numeric d) const noexcept { return -d * dF; }
  constexpr Complex dFdG0(Numeric d) const noexcept { return Complex(0, d) * dF; }
  constexpr Complex dFdH(Numeric dZ) const noexcept { return -dZ * dF; }
  constexpr Complex dFdVMR(const Output &d) const noexcept {
    return Complex(-d.D0 - d.DV, d.G0) * dF;
  }
  constexpr Complex dFdT(const Output &d, Numeric T) const noexcept {
    return -(F * invGD + dF * z) * (2 * T * (d.D0 + d.DV) + mF0) /
    (2 * T * invGD * mF0) +
    Complex(-d.D0 - d.DV, d.G0) * dF;
  }
  static constexpr Complex dFdETA(Numeric) noexcept { return 0; }
  static constexpr Complex dFdFVC(Numeric) noexcept { return 0; }
  static constexpr Complex dFdD2(Numeric) noexcept { return 0; }
  static constexpr Complex dFdG2(Numeric) { return 0; }

  Complex operator()(Numeric f) noexcept;

  bool OK() const noexcept { return invGD > 0; }
}; // Voigt

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

  SpeedDependentVoigt(Numeric F0_noshift, const Output &ls, Numeric GD_div_F0,
                      Numeric dZ) noexcept;

  Complex dFdf() const noexcept;
  Complex dFdF0() const noexcept;
  Complex dFdD0(Numeric dD0dD0) const noexcept;
  Complex dFdG0(Numeric dG0dG0) const noexcept;
  Complex dFdD2(Numeric dD2dD2) const noexcept;
  Complex dFdG2(Numeric dG2dG2) const noexcept;
  Complex dFdH(Numeric dZ) const noexcept;
  Complex dFdVMR(const Output &d) const noexcept;
  Complex dFdT(const Output &d, Numeric T) const noexcept;
  static constexpr Complex dFdETA(Numeric) noexcept { return 0; }
  static constexpr Complex dFdFVC(Numeric) noexcept { return 0; }
  static constexpr Complex dFdDV(Numeric) noexcept { return 0; }

  Complex operator()(Numeric f) noexcept;

  CalcType init(const Complex c2) const noexcept;
  void update_calcs() noexcept;
  void calc() noexcept;
}; // SpeedDependentVoigt

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

  HartmannTran(Numeric F0_noshift, const Output &ls, Numeric GD_div_F0,
               Numeric dZ) noexcept;

  Complex dFdf() const noexcept;
  Complex dFdF0() const noexcept;
  Complex dFdD0(Numeric dD0) const noexcept;
  Complex dFdG0(Numeric dG0) const noexcept;
  Complex dFdD2(Numeric dD2) const noexcept;
  Complex dFdG2(Numeric dG2) const noexcept;
  Complex dFdFVC(Numeric dFVC) const noexcept;
  Complex dFdETA(Numeric dETA) const noexcept;
  Complex dFdH(Numeric dZ) const noexcept;
  Complex dFdVMR(const Output &d) const noexcept;
  Complex dFdT(const Output &d, Numeric T) const noexcept;
  static constexpr Complex dFdDV(Numeric) noexcept { return 0; }

  Complex operator()(Numeric f) noexcept;

  CalcType init(const Complex c2t) const noexcept;
  void update_calcs() noexcept;
  void calc() noexcept;
}; // HartmannTran

struct Nonorm {
  static constexpr Numeric N = 1.0;

  constexpr Nonorm() noexcept {}

  static constexpr Numeric dNdT(Numeric, Numeric) noexcept { return 0; }
  static constexpr Numeric dNdf(Numeric) noexcept { return 0; }
  constexpr Numeric dNdF0() const noexcept { return 0; }

  constexpr Numeric operator()(Numeric) noexcept { return N; }
}; // Nonorm

struct VanVleckHuber {
  Numeric N;
  
  Numeric c1;
  Numeric tanh_c1f0;
  Numeric inv_denom;
  Numeric tanh_c1f;

  VanVleckHuber(Numeric F0, Numeric T) noexcept;

  Numeric dNdT(Numeric T, Numeric f) const noexcept;
  Numeric dNdf(Numeric f) const noexcept;
  Numeric dNdF0() const noexcept;

  Numeric operator()(Numeric f) noexcept;
}; // VanVleckHuber

struct VanVleckWeisskopf {
  Numeric N;
  
  Numeric invF0;

  constexpr VanVleckWeisskopf(Numeric F0) noexcept : N(1), invF0(1.0 / F0) {}

  static constexpr Numeric dNdT(Numeric, Numeric) noexcept { return 0; }
  constexpr Numeric dNdf(Numeric f) const noexcept {
    return 2.0 * f * Constant::pow2(invF0);
  }
  constexpr Numeric dNdF0() const noexcept { return -2.0 * N * invF0; }

  constexpr Numeric operator()(Numeric f) noexcept {
    N = Constant::pow2(f * invF0);
    return N;
  }
}; // VanVleckWeisskopf

struct RosenkranzQuadratic {
  Numeric N;
  
  Numeric fac;
  Numeric dfacdT;
  Numeric dfacdF0;

  RosenkranzQuadratic(Numeric F0, Numeric T) noexcept;

  Numeric dNdT(Numeric, Numeric f) const noexcept;
  Numeric dNdf(Numeric f) const noexcept;
  Numeric dNdF0() const noexcept;

  Numeric operator()(Numeric f) noexcept;
}; // RosenkranzQuadratic

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
  
  static constexpr bool do_nlte() noexcept {return false;}
}; // Nostrength

struct LocalThermodynamicEquilibrium {
  Numeric S;
  static constexpr Numeric N = 0.0;
  
  Numeric dSdI0val;
  Numeric dSdTval;
  Numeric dSdF0val;
  Numeric dSdSELFVMRval;

  constexpr LocalThermodynamicEquilibrium(Numeric I0, Numeric r, Numeric drdSELFVMR, Numeric drdT, Numeric QT0,
                                          Numeric QT, Numeric dQTdT, Numeric br,
                                          Numeric dbr_dT_rat, Numeric stim,
                                          Numeric dstim_dT, Numeric dstim_dF0) noexcept
      : S(), dSdI0val(r * br * stim * QT0 / QT),
        dSdTval(I0 * (r * br * dstim_dT * QT0 / QT + dSdI0val * (dbr_dT_rat - dQTdT / QT) + drdT * br * stim * QT0 / QT)),
        dSdF0val(r * I0 * br * dstim_dF0 * QT0 / QT),
        dSdSELFVMRval(drdSELFVMR * I0 * br * stim * QT0 / QT)
  {
    S = I0 * dSdI0val;
  }

  LocalThermodynamicEquilibrium(Numeric I0, Numeric T0, Numeric T, Numeric F0,
                                Numeric E0, Numeric QT, Numeric QT0,
                                Numeric dQTdT, Numeric r, Numeric drdSELFVMR, Numeric drdT) noexcept;

  constexpr Numeric dSdT() const noexcept { return dSdTval; }
  constexpr Numeric dSdI0() const noexcept { return dSdI0val; }
  constexpr Numeric dSdF0() const noexcept { return dSdF0val; }
  static constexpr Numeric dSdNLTEu() noexcept { return 0; }
  static constexpr Numeric dSdNLTEl() noexcept { return 0; }
  constexpr Numeric dSdSELFVMR() const noexcept { return dSdSELFVMRval; }

  static constexpr Numeric dNdT() noexcept { return 0; }
  static constexpr Numeric dNdI0() noexcept { return 0; }
  static constexpr Numeric dNdF0() noexcept { return 0; }
  static constexpr Numeric dNdNLTEu() noexcept { return 0; }
  static constexpr Numeric dNdNLTEl() noexcept { return 0; }
  static constexpr Numeric dNdSELFVMR() noexcept { return 0; }
  
  static constexpr bool do_nlte() noexcept {return false;}
}; // LocalThermodynamicEquilibrium

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
  
  constexpr FullNonLocalThermodynamicEquilibrium(
    Numeric r, Numeric drdSELFVMR, Numeric drdt,
    Numeric k, Numeric dkdF0, Numeric dkdr1, Numeric dkdr2,
    Numeric e, Numeric dedF0, Numeric dedr2,
    Numeric B, Numeric dBdT, Numeric dBdF0) noexcept :
  S(), N(),
  dSdTval(drdt * k),
  dNdTval(drdt * (e - k * B) - r * k * dBdT),
  dSdF0val(r * dkdF0),
  dNdF0val(r * (dedF0 - dkdF0 * B - k * dBdF0)),
  dSdr1(r * dkdr1),
  dSdr2(r * dkdr2),
  dNdr2(r * (dedr2 - dkdr2 * B)),
  dSdSELFVMRval(drdSELFVMR * k),
  dNdSELFVMRval(drdSELFVMR * (e - k * B))
  {
    S = r * k;
    N = r * (e - k * B);
  }
  
  FullNonLocalThermodynamicEquilibrium(Numeric F0, Numeric A21, Numeric T,
                                       Numeric g1, Numeric g2, Numeric r1,
                                       Numeric r2, Numeric r, Numeric drdSELFVMR, Numeric drdT) noexcept;
                                         

  constexpr Numeric dSdT() const noexcept { return 0; }
  static constexpr Numeric dSdI0() noexcept { return 0; }
  constexpr Numeric dSdF0() const noexcept { return dSdF0val; }
  constexpr Numeric dSdNLTEu() const noexcept { return dSdr1; }
  constexpr Numeric dSdNLTEl() const noexcept { return dSdr2; }
  constexpr Numeric dSdSELFVMR() const noexcept { return dSdSELFVMRval; }

  constexpr Numeric dNdT() const noexcept { return dNdTval; }
  static constexpr Numeric dNdI0() noexcept { return 0; }
  constexpr Numeric dNdF0() const noexcept { return dNdF0val; }
  constexpr Numeric dNdNLTEu() const noexcept { return -dSdr1; }
  constexpr Numeric dNdNLTEl() const noexcept { return dNdr2; }
  constexpr Numeric dNdSELFVMR() const noexcept { return dNdSELFVMRval; }
  
  static constexpr bool do_nlte() noexcept {return true;}
}; // FullNonLocalThermodynamicEquilibrium

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

  constexpr VibrationalTemperaturesNonLocalThermodynamicEquilibrium(Numeric I0,
                                                                    Numeric QT0, Numeric QT, Numeric dQTdT,
                                                                    Numeric r, Numeric drdSELFVMR, Numeric drdT,
                                                                    Numeric K1, Numeric dK1dT, 
                                                                    Numeric K2, Numeric dK2dT, Numeric dK2dF0, 
                                                                    Numeric K3, Numeric dK3dT, Numeric dK3dF0, Numeric dK3dTl, Numeric dK3dTu,
                                                                    Numeric K4, Numeric dK4dT, Numeric dK4dTu,
                                                                    Numeric B, Numeric dBdT, Numeric dBdF0) noexcept :
  S(), N(),
  dSdI0val(    r * QT0 / QT * K1 * K2 * K3),
  dNdI0val(B * r * QT0 / QT * K1 * K2 * (K4 - K3)),
  dSdTval(I0 * (drdT * QT0 / QT * K1 * K2 * K3 -
                r * dQTdT * QT0 / Constant::pow2(QT) * K1 * K2 * K3 + 
                r * QT0 / QT * dK1dT * K2 * K3 + 
                r * QT0 / QT * K1 * dK2dT * K3 + 
                r * QT0 / QT * K1 * K2 * dK3dT)),
  dNdTval(I0 * (dBdT * r * QT0 / QT * K1 * K2 * (K4 - K3) +
                B * drdT * QT0 / QT * K1 * K2 * (K4 - K3) -
                B * r * dQTdT * QT0 / Constant::pow2(QT) * K1 * K2 * (K4 - K3) +
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
  dNdSELFVMRval(I0 * B * drdSELFVMR * QT0 / QT * K1 * K2 * (K4 - K3))
  {
    S = I0 * dSdI0val;
    N = I0 * dNdI0val;
  }

  VibrationalTemperaturesNonLocalThermodynamicEquilibrium(
      Numeric I0, Numeric T0, Numeric T, Numeric Tl, Numeric Tu, Numeric F0,
      Numeric E0, Numeric Evl, Numeric Evu, Numeric QT, Numeric QT0,
      Numeric dQTdT, Numeric r, Numeric drdSELFVMR, Numeric drdT) noexcept;

      constexpr Numeric dSdT() const noexcept { return dSdTval; }
  constexpr Numeric dSdI0() const noexcept { return dSdI0val; }
  constexpr Numeric dSdF0() const noexcept { return dSdF0val; }
  constexpr Numeric dSdNLTEu() const noexcept { return dSdTl; }
  constexpr Numeric dSdNLTEl() const noexcept { return dSdTu; }
  constexpr Numeric dSdSELFVMR() const noexcept { return dSdSELFVMRval; }

  constexpr Numeric dNdT() const noexcept { return dNdTval; }
  constexpr Numeric dNdI0() const noexcept { return dNdI0val; }
  constexpr Numeric dNdF0() const noexcept { return dNdF0val; }
  constexpr Numeric dNdNLTEu() const noexcept { return -dSdTl; }
  constexpr Numeric dNdNLTEl() const noexcept { return dNdTu; }
  constexpr Numeric dNdSELFVMR() const noexcept { return dNdSELFVMRval; }
  
  static constexpr bool do_nlte() noexcept {return true;}
}; // VibrationalTemperaturesNonLocalThermodynamicEquilibrium

typedef std::variant<Noshape, Doppler, Lorentz, Voigt, SpeedDependentVoigt,
                     HartmannTran>
    Calculator;

typedef std::variant<Nonorm, VanVleckHuber, VanVleckWeisskopf,
                     RosenkranzQuadratic>
    Normalizer;

typedef std::variant<Nostrength, LocalThermodynamicEquilibrium,
                     FullNonLocalThermodynamicEquilibrium,
                     VibrationalTemperaturesNonLocalThermodynamicEquilibrium>
    IntensityCalculator;

struct ComputeData {
  ComplexVector F, N;
  ComplexMatrix dF, dN;
  const Vector & f_grid;
  const bool do_nlte;
  
  ComputeData(const Vector& f, const ArrayOfRetrievalQuantity &jacobian_quantities, const bool nlte) noexcept : 
  F(f.nelem(), 0),
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
  
  void interp_add_even(const ComputeData& sparse) ARTS_NOEXCEPT {
    const Index nv = f_grid.nelem();
    const Index nj = dF.ncols();
    
    ARTS_ASSERT(do_nlte == sparse.do_nlte, "Must have the same NLTE status")
    ARTS_ASSERT(sparse.f_grid.nelem() > 1, "Must have at least two sparse grid-points")
    ARTS_ASSERT(nv == 0 or (f_grid[0] >= sparse.f_grid[0] and f_grid[nv - 1] <= sparse.f_grid[sparse.f_grid.nelem() - 1]),
                "If there are any dense frequency points, then the sparse frequency points must fully cover them")
    
    Index sparse_iv=1;
    const Numeric invdf = 1.0 / (sparse.f_grid[sparse_iv] - sparse.f_grid[sparse_iv-1]);
    for (Index iv=0; iv<nv; iv++) {
      while (sparse.f_grid[sparse_iv] < f_grid[iv]) sparse_iv++;  // sparse_iv cannot exceed sparse_nv by assert
      
      // Interpolation weight
      const Numeric x = (f_grid[iv] - sparse.f_grid[sparse_iv-1]) * invdf;
      
      F[iv] += sparse.F[sparse_iv-1] + x * (sparse.F[sparse_iv] - sparse.F[sparse_iv-1]);
      for (Index ij=0; ij<nj; ij++) {
        dF(iv, ij) += sparse.dF(sparse_iv-1, ij) + x * (sparse.dF(sparse_iv, ij) - sparse.dF(sparse_iv-1, ij));
      }
      if (do_nlte) {
        N[iv] += sparse.N[sparse_iv-1] + x * (sparse.N[sparse_iv] - sparse.N[sparse_iv-1]);
        for (Index ij=0; ij<nj; ij++) {
          dN(iv, ij) += sparse.dN(sparse_iv-1, ij) + x * (sparse.dN(sparse_iv, ij) - sparse.dN(sparse_iv-1, ij));
        }
      }
    }
  }
  
  void interp_add_triplequad(const ComputeData& sparse) ARTS_NOEXCEPT {
    const Index nv = f_grid.nelem();
    const Index sparse_nv = sparse.f_grid.nelem();
    const Index nj = dF.ncols();
    
    ARTS_ASSERT(do_nlte == sparse.do_nlte, "Must have the same NLTE status")
    ARTS_ASSERT(sparse_nv > 2, "Must have at least three sparse grid-points")
    ARTS_ASSERT(nv == 0 or (f_grid[0] == sparse.f_grid[0] and f_grid[nv - 1] >= sparse.f_grid[sparse_nv - 1]),
                "If there are any dense frequency points, then the sparse frequency points must fully cover them")
    ARTS_ASSERT(not (sparse_nv % 3), "Must be multiple of three")
    
    Index sparse_iv=0;
    const Numeric invdf2 = 1.0 / Constant::pow2(sparse.f_grid[sparse_iv+1] - sparse.f_grid[sparse_iv]);
    const Numeric invdf2_last = 1.0 / Constant::pow2(sparse.f_grid[sparse_nv - 1] - sparse.f_grid[sparse_nv - 2]);
    for (Index iv = 0; iv < nv; iv++) {
      
      if ((sparse_iv not_eq sparse_nv - 3) and 
          (sparse.f_grid[sparse_iv + 2] == f_grid[iv])) {
        sparse_iv += 3;
      }
      
      if (sparse_iv == sparse_nv - 3) {
        
        // Interpolation weights
        // l0 is ((f - f1) / (f0 - f1)) * ((f - f2) / (f0 - f2))
        // f0 - f2 == 2 * (f0 - f1)
        // l0 is ((f - f1) / (f0 - f1)) * 0.5 * ((f - f2) / (f0 - f1))
        // 1.0 / (f0 - f1) == - 1 / df
        // l0 = 0.5 * ((f - f1) / -df) * ((f - f2) / -df)
        // l0 = 0.5 * (f - f1) * (f - f2) / df^2
        const Numeric xm0 = f_grid[iv] - sparse.f_grid[sparse_nv - 3];
        const Numeric xm1 = f_grid[iv] - sparse.f_grid[sparse_nv - 2];
        const Numeric xm2 = f_grid[iv] - sparse.f_grid[sparse_nv - 1];
        const Numeric l0 = 0.5 * xm1 * xm2 * invdf2_last;  // --
        const Numeric l1 = - xm0 * xm2 * invdf2_last;      // +-
        const Numeric l2 = 0.5 * xm0 * xm1 * invdf2_last;  // ++
        
        F[iv] += l0 * sparse.F[sparse_nv - 3] + l1 * sparse.F[sparse_nv - 2] + l2 * sparse.F[sparse_nv - 1];
        for (Index ij=0; ij<nj; ij++) {
          dF(iv, ij) += l0 * sparse.dF(sparse_nv - 3, ij) + l1 * sparse.dF(sparse_nv - 2, ij) + l2 * sparse.dF(sparse_nv - 1, ij);
        }
        if (do_nlte) {
          N[iv] += l0 * sparse.N[sparse_nv - 3] + l1 * sparse.N[sparse_nv - 2] + l2 * sparse.N[sparse_nv - 1];
          for (Index ij=0; ij<nj; ij++) {
            dN(iv, ij) += l0 * sparse.dN(sparse_nv - 3, ij) + l1 * sparse.dN(sparse_nv - 2, ij) + l2 * sparse.dN(sparse_nv - 1, ij);
          }
        }
      } else {
        
        // Interpolation weight
        const Numeric xm0 = f_grid[iv] - sparse.f_grid[sparse_iv + 0];
        const Numeric xm1 = f_grid[iv] - sparse.f_grid[sparse_iv + 1];
        const Numeric xm2 = f_grid[iv] - sparse.f_grid[sparse_iv + 2];
        const Numeric l0 = 0.5 * xm1 * xm2 * invdf2;  // --
        const Numeric l1 = - xm0 * xm2 * invdf2;      // +-
        const Numeric l2 = 0.5 * xm0 * xm1 * invdf2;  // ++
        
        F[iv] += l0 * sparse.F[sparse_iv + 0] + l1 * sparse.F[sparse_iv + 1] + l2 * sparse.F[sparse_iv + 2];
        for (Index ij=0; ij<nj; ij++) {
          dF(iv, ij) += l0 * sparse.dF(sparse_iv + 0, ij) + l1 * sparse.dF(sparse_iv + 1, ij) + l2 * sparse.dF(sparse_iv + 2, ij);
        }
        if (do_nlte) {
          N[iv] += l0 * sparse.N[sparse_iv + 0] + l1 * sparse.N[sparse_iv + 1] + l2 * sparse.N[sparse_iv + 2];
          for (Index ij=0; ij<nj; ij++) {
            dN(iv, ij) += l0 * sparse.dN(sparse_iv + 0, ij) + l1 * sparse.dN(sparse_iv + 1, ij) + l2 * sparse.dN(sparse_iv + 2, ij);
          }
        }
      }
    }
  }
};

void compute(ComputeData &com,
             ComputeData &sparse_com,
             const AbsorptionLines &band,
             const ArrayOfRetrievalQuantity &jacobian_quantities,
             const EnergyLevelMap &nlte,
             const SpeciesAuxData::AuxType &partfun_type,
             const ArrayOfGriddedField1 &partfun_data, const Vector &vmrs,
             const Numeric &self_vmr, const Numeric &isot_ratio, const Numeric &P, const Numeric &T, const Numeric &H, const Numeric &sparse_lim,
             const bool do_zeeman, const Zeeman::Polarization zeeman_polarization, const Options::LblSpeedup speedup_type) ARTS_NOEXCEPT;

} // namespace LineShape

#endif // lineshapes_h
