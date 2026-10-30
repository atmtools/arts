/**
 * @file   physics_funcs.h
 * @author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
 * @date   2002-05-08
 *
 * @brief This file contains declerations of functions of physical character.
*/

#ifndef physics_h
#define physics_h

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <arts_constants.h>
#include <arts_constexpr_math.h>
#include <arts_conversions.h>
#include <matpack.h>

/*===========================================================================
  === Functions in physics_funcs.h
  ===========================================================================*/

Numeric barometric_heightformula(const Numeric& p, const Numeric& dh);

Numeric dinvplanckdI(const Numeric& i, const Numeric& f);

void fresnel(Complex& Rv, Complex& Rh, const Complex& n1, const Complex& n2, const Numeric& theta);
std::pair<Complex, Complex> fresnel(const Complex& n1, const Complex& n2, const Numeric& theta);

Numeric invplanck(const Numeric& i, const Numeric& f);

Numeric invrayjean(const Numeric& i, const Numeric& f);

/** number_density
 * 
 * Calculates the atmospheric number density.
 *
 * @param[in]  p  Pressure.
 * @param[in]  t  Temperature.
 *
 * @return     Number density.
 *
 * @author Patrick Eriksson
 * @date   2000-04-08
 */
constexpr Numeric number_density(Numeric p, Numeric t) noexcept { return p / (Constant::k * t); }

/** dnumber_density_dT
 * 
 * Calculates the atmospheric number density derivative with temperature.
 *
 * @param[in]  p  Pressure.
 * @param[in]  t  Temperature.
 *
 * @return     Number density.
 *
 * @author Richard Larsson
 * @date   2015-09-22
 */
constexpr Numeric dnumber_density_dt(Numeric p, Numeric t) noexcept { return -p / (Constant::k * Math::pow2(t)); }

Numeric planck(const Numeric& f, const Numeric& t);

void planck(StridedVectorView b, const ConstVectorView& f, const Numeric& t);

Vector planck(const StridedConstVectorView& f, const Numeric& t);

Numeric dplanck_dt(const Numeric& f, const Numeric& t);

void dplanck_dt(VectorView b, const ConstVectorView& f, const Numeric& t);

Vector dplanck_dt(const ConstVectorView& f, const Numeric& t);

Numeric dplanck_df(const Numeric& f, const Numeric& t);

Vector dplanck_df(const ConstVectorView& f, const Numeric& t);

Numeric rayjean(const Numeric& f, const Numeric& t);

/** Refractive index of water and steam in the visible and near infrared.
 *
 * Implements the formulation of Harvey et al. (1998), "Revised
 * Formulation for the Refractive Index of Water and Steam as a Function of
 * Wavelength, Temperature and Density".  The formulation only supplies the
 * real part of the refractive index.
 *
 * @param frequency Frequency [Hz].
 * @param temperature Temperature [K].
 * @param density Water mass density [kg/m3].
 * @param check_validity Enforce the published validity range.
 * @return Real refractive index.
 */
Numeric refractive_index_water_visible_nir_harvey98(Numeric frequency,
                                                    Numeric temperature,
                                                    Numeric density,
                                                    bool    check_validity = true);

/** Complex permittivity of liquid water in the microwave.
 *
 * Implements the double-Debye model of Liebe, Hufford and Manabe (1991), "A
 * model for the complex permittivity of water at frequencies below 1 THz",
 * as ARTS2's complex_n_water_liebe93 and Atmlab's eps_water_liebe93.m.  The
 * relaxation-frequency coefficient is 146, not the 146.4 of epswater93.m.
 * The imaginary part is positive.
 *
 * @param frequency Frequency [Hz].
 * @param temperature Temperature [K].
 * @param check_validity Enforce 10 GHz <= frequency <= 1000 GHz and
 * 233.15 K <= temperature <= 373.15 K.
 * @return Complex relative permittivity.
 */
Complex permittivity_water_liebe93(Numeric frequency, Numeric temperature, bool check_validity = true);

#endif  // physics_h
