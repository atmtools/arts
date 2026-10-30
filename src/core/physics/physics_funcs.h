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

/** Dispersion from the imaginary part of a complex absorption coefficient
 *
 * The real part of the complex absorption coefficient, k, is the absorption
 * coefficient, 4 pi f Im(n) / c.  Its imaginary part gives the refractivity,
 * n - 1 = -c Im(k) / (4 pi f), where n is the complex refractive index.
 *
 * Use this also for derivatives of Im(k) with regards to anything but frequency.
 *
 * @param[in] f Frequency [Hz]
 * @param[in] k_imag Imaginary part of the complex absorption coefficient [1/m]
 * @return The refractivity, n - 1 [-]
 */
constexpr Numeric dispersion_from_imag_absorption(Numeric f, Numeric k_imag) noexcept {
  return -Constant::c * k_imag / (4.0 * Constant::pi * f);
}

/** Frequency derivative of dispersion_from_imag_absorption
 *
 * @param[in] f Frequency [Hz]
 * @param[in] k_imag Imaginary part of the complex absorption coefficient [1/m]
 * @param[in] dk_imag_df Frequency derivative of k_imag [1/(m Hz)]
 * @return The frequency derivative of the refractivity, n - 1 [1/Hz]
 */
constexpr Numeric ddispersion_from_imag_absorption_df(Numeric f, Numeric k_imag, Numeric dk_imag_df) noexcept {
  return -Constant::c * (dk_imag_df - k_imag / f) / (4.0 * Constant::pi * f);
}

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

#endif  // physics_h
