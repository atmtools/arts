#ifndef lineshapes_h
#define lineshapes_h

#include <variant>

#include "constants.h"
#include "energylevelmap.h"
#include "linescaling.h"

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

class Doppler {
  Numeric mF0;
  Numeric invGD;
  Numeric x;

public:
  Complex F;

  Doppler(Numeric F0_noshift, Numeric DC, Numeric dZ) noexcept;

  Complex dFdT(const Output &, Numeric T) const noexcept;
  Complex dFdf() const noexcept;
  Complex dFdF0() const noexcept;
  Complex dFdH(Numeric dZ) const noexcept;
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

class Lorentz {
  Numeric mF0;
  Numeric G0;

public:
  Complex F;

private:
  Complex dF;

public:
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

class Voigt {
  Numeric mF0;
  Numeric invGD;
  Complex z;

public:
  Complex F;

private:
  Complex dF;

public:
  Voigt(Numeric F0_noshift, const Output &ls, Numeric DC, Numeric dZ) noexcept;

  Complex dFdf() const noexcept;
  Complex dFdF0() const noexcept;
  Complex dFdDV(Numeric d) const noexcept;
  Complex dFdD0(Numeric d) const noexcept;
  Complex dFdG0(Numeric d) const noexcept;
  Complex dFdH(Numeric dZ) const noexcept;
  Complex dFdVMR(const Output &d) const noexcept;
  Complex dFdT(const Output &d, Numeric T) const noexcept;
  static constexpr Complex dFdETA(Numeric) noexcept { return 0; }
  static constexpr Complex dFdFVC(Numeric) noexcept { return 0; }
  static constexpr Complex dFdD2(Numeric) noexcept { return 0; }
  static constexpr Complex dFdG2(Numeric) { return 0; }

  Complex operator()(Numeric f) noexcept;

  bool OK() const noexcept { return invGD > 0; }
}; // Voigt

class SpeedDependentVoigt {
  enum class CalcType : char {
    Voigt,
    LowXandHighY,
    LowYandLowX,
    LowYandHighX,
    Full
  };

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

public:
  Complex F;

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

private:
  CalcType init(const Complex c2) const noexcept;
  void update_calcs() noexcept;
  void calc() noexcept;
}; // SpeedDependentVoigt

class HartmannTran {
  enum class CalcType : char {
    Noc2tLowZ,
    Noc2tHighZ,
    LowXandHighY,
    LowYandLowX,
    LowYandHighX,
    Full
  };

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

public:
  Complex F;

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

private:
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

class VanVleckHuber {
  Numeric c1;
  Numeric tanh_c1f0;
  Numeric inv_denom;
  Numeric tanh_c1f;

public:
  Numeric N;

  VanVleckHuber(Numeric F0, Numeric T) noexcept;

  Numeric dNdT(Numeric T, Numeric f) const noexcept;
  Numeric dNdf(Numeric f) const noexcept;
  Numeric dNdF0() const noexcept;

  Numeric operator()(Numeric f) noexcept;
}; // VanVleckHuber

class VanVleckWeisskopf {
  Numeric invF0;

public:
  Numeric N;

  constexpr VanVleckWeisskopf(Numeric F0) noexcept : invF0(1.0 / F0), N(1) {}

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

class RosenkranzQuadratic {
  Numeric fac;
  Numeric dfacdT;
  Numeric dfacdF0;

public:
  Numeric N;

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
}; // Nostrength

class LocalThermodynamicEquilibrium {
  Numeric dSdI0val;
  Numeric dSdTval;
  Numeric dSdF0val;
  Numeric dSdSELFVMRval;

public:
  Numeric S;
  static constexpr Numeric N = 0.0;

  constexpr LocalThermodynamicEquilibrium(Numeric I0, Numeric r, Numeric drdSELFVMR, Numeric drdT, Numeric QT0,
                                          Numeric QT, Numeric dQTdT, Numeric br,
                                          Numeric dbr_dF0_rat, Numeric stim,
                                          Numeric dstim_dT, Numeric dstim_dF0) noexcept
      : dSdI0val(r * br * stim * QT0 / QT),
        dSdTval(I0 * (r * br * dstim_dT * QT0 / QT + dSdI0val * (dbr_dF0_rat - dQTdT / QT) + drdT * br * stim * QT0 / QT)),
        dSdF0val(r * I0 * br * dstim_dF0 * QT0 / QT), dSdSELFVMRval(drdSELFVMR * I0 * br * stim * QT0 / QT), S(I0 * dSdI0val) {}

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
}; // LocalThermodynamicEquilibrium

class FullNonLocalThermodynamicEquilibrium {
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
  
public:
  Numeric S;
  Numeric N;
  
  constexpr FullNonLocalThermodynamicEquilibrium(
    Numeric r, Numeric drdSELFVMR, Numeric drdt,
    Numeric k, Numeric dkdF0, Numeric dkdr1, Numeric dkdr2,
    Numeric e, Numeric dedF0, Numeric dedr2,
    Numeric B, Numeric dBdT, Numeric dBdF0) noexcept :
  dSdTval(drdt * k),
  dNdTval(drdt * (e - k * B) - r * k * dBdT),
  dSdF0val(r * dkdF0),
  dNdF0val(r * (dedF0 - dkdF0 * B - k * dBdF0)),
  dSdr1(r * dkdr1),
  dSdr2(r * dkdr2),
  dNdr2(r * (dedr2 - dkdr2 * B)),
  dSdSELFVMRval(drdSELFVMR * k),
  dNdSELFVMRval(drdSELFVMR * (e - k * B)),
  S(r * k),
  N(r * (e - k * B))
  {}
  
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
}; // FullNonLocalThermodynamicEquilibrium

class VibrationalTemperaturesNonLocalThermodynamicEquilibrium {
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

public:
  Numeric S;
  Numeric N;
  
  constexpr VibrationalTemperaturesNonLocalThermodynamicEquilibrium(Numeric I0,
                                                                    Numeric QT0, Numeric QT, Numeric dQTdT,
                                                                    Numeric r, Numeric drdSELFVMR, Numeric drdT,
                                                                    Numeric K1, Numeric dK1dT, 
                                                                    Numeric K2, Numeric dK2dT, Numeric dK2dF0, 
                                                                    Numeric K3, Numeric dK3dT, Numeric dK3dF0, Numeric dK3dTl, Numeric dK3dTu,
                                                                    Numeric K4, Numeric dK4dT, Numeric dK4dTu,
                                                                    Numeric B, Numeric dBdT, Numeric dBdF0) noexcept :
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
  dNdSELFVMRval(I0 * B * drdSELFVMR * QT0 / QT * K1 * K2 * (K4 - K3)),
  S(I0 * dSdI0val),
  N(I0 * dNdI0val)
  {}

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

void compute(ComplexVector &F, ComplexMatrix &dF, ComplexVector &N,
             ComplexMatrix &dN, const Vector &f_grid,
             const AbsorptionLines &band,
             const ArrayOfRetrievalQuantity &jacobian_quantities,
             const EnergyLevelMap &nlte,
             const SpeciesAuxData::AuxType &partfun_type,
             const ArrayOfGriddedField1 &partfun_data, const Vector &vmrs,
             const Numeric &isot_ratio, const Numeric &P, const Numeric &T,
             const bool do_nlte = false, const Numeric &H = 0,
             const bool do_zeeman = false,
             const Zeeman::Polarization zeeman_polarization =
                 Zeeman::Polarization::Pi, const Numeric& self_vmr=std::numeric_limits<Numeric>::quiet_NaN(), const bool do_numden=false) ARTS_NOEXCEPT;

} // namespace LineShape

#endif // lineshapes_h
