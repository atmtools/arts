#ifndef lineshapes_h
#define lineshapes_h

#include <variant>

#include "energylevelmap.h"

namespace LineShape {
struct Noshape {
  static constexpr Complex F = Complex(0, 0);

  static constexpr Complex dFdT(Output, Numeric) noexcept { return 0; }
  static constexpr Complex dFdf() noexcept { return 0; }
  static constexpr Complex dFdF0() noexcept { return 0; }
  static constexpr Complex dFdH(Numeric) noexcept { return 0; }
  static constexpr Complex dFdFVC(Numeric) noexcept { return 0; }
  static constexpr Complex dFdETA(Numeric) noexcept { return 0; }
  static constexpr Complex dFdVMR(Output) noexcept { return 0; }
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

  Complex dFdT(Output, Numeric T) const noexcept;
  Complex dFdf() const noexcept;
  Complex dFdF0() const noexcept;
  Complex dFdH(Numeric dZ) const noexcept;
  static constexpr Complex dFdFVC(Numeric) noexcept { return 0; }
  static constexpr Complex dFdETA(Numeric) noexcept { return 0; }
  static constexpr Complex dFdVMR(Output) noexcept { return 0; }
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
  constexpr Lorentz(Numeric F0_noshift, Output ls) noexcept
      : mF0(F0_noshift + ls.D0 + ls.DV), G0(ls.G0) {}

  constexpr Complex dFdVMR(Output d) const noexcept {
    return Complex(d.G0, d.D0 + d.DV) * dF;
  }
  constexpr Complex dFdT(Output d, Numeric) const noexcept { return dFdVMR(d); }
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
  Voigt(Numeric F0_noshift, Output ls, Numeric DC, Numeric dZ) noexcept;

  Complex dFdf() const noexcept;
  Complex dFdF0() const noexcept;
  Complex dFdDV(Numeric d) const noexcept;
  Complex dFdD0(Numeric d) const noexcept;
  Complex dFdG0(Numeric d) const noexcept;
  Complex dFdH(Numeric dZ) const noexcept;
  Complex dFdVMR(Output d) const noexcept;
  Complex dFdT(Output d, Numeric T) const noexcept;
  static constexpr Complex dFdETA(Numeric) noexcept { return 0; }
  static constexpr Complex dFdFVC(Numeric) noexcept { return 0; }
  static constexpr Complex dFdD2(Numeric) noexcept { return 0; }
  static constexpr Complex dFdG2(Numeric) { return 0; }

  Complex operator()(Numeric f) noexcept;
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

  SpeedDependentVoigt(Numeric F0_noshift, Output ls, Numeric GD_div_F0,
                      Numeric dZ) noexcept;

  Complex dFdf() const noexcept;
  Complex dFdF0() const noexcept;
  Complex dFdD0(Numeric dD0dD0) const noexcept;
  Complex dFdG0(Numeric dG0dG0) const noexcept;
  Complex dFdD2(Numeric dD2dD2) const noexcept;
  Complex dFdG2(Numeric dG2dG2) const noexcept;
  Complex dFdH(Numeric dZ) const noexcept;
  Complex dFdVMR(Output d) const noexcept;
  Complex dFdT(Output d, Numeric T) const noexcept;
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

  HartmannTran(Numeric F0_noshift, Output ls, Numeric GD_div_F0,
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
  Complex dFdVMR(Output d) const noexcept;
  Complex dFdT(Output d, Numeric T) const noexcept;
  static constexpr Complex dFdDV(Numeric) noexcept { return 0; }

  Complex operator()(Numeric f) noexcept;

private:
  CalcType init(const Complex c2t) const noexcept;
  void update_calcs() noexcept;
  void calc() noexcept;
}; // HartmannTran

typedef std::variant<Noshape, Doppler, Lorentz, Voigt, SpeedDependentVoigt,
                     HartmannTran>
    Calculator;

void compute(ComplexVector &F, ArrayOfComplexVector &dF,
             const ConstVectorView &f_grid, const AbsorptionLines &band,
             const ArrayOfRetrievalQuantity &derivatives_data,
             const Vector &vmrs, const Numeric &P, const Numeric &T,
             const Numeric &H, const bool zeeman,
             const Zeeman::Polarization zeeman_polarization) ARTS_NOEXCEPT;

} // namespace LineShape

#endif // lineshapes_h
