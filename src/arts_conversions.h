/* Copyright (C) 2019
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/*!
 * \file   arts_conversions.h
 * \brief  Common ARTS conversions
 *
 * Where possible these are going to
 * 
 * \author Richard Larsson
 * \date   2019-04-01
 */

#ifndef CONVERSIONS_IN_ARTS_H
#define CONVERSIONS_IN_ARTS_H

#include <cmath>

#include "arts_constants.h"

/** Namespace containing several practical unit conversions, physical and mathematical **/
namespace Conversion {
using namespace Constant;

/** Converts degrees to radians */
constexpr auto deg2rad(auto x) noexcept { return x * one_degree_in_radians; }

/** Converts radians to degrees */
constexpr auto rad2deg(auto x) noexcept { return x / one_degree_in_radians; }

/** Returns the cosine of the deg2rad of the input */
auto cosd(auto x) noexcept { return std::cos(deg2rad(x)); }

/** Returns the sine of the deg2rad of the input */
auto sind(auto x) noexcept { return std::sin(deg2rad(x)); }

/** Returns the tangent of the deg2rad of the input */
auto tand(auto x) noexcept { return std::tan(deg2rad(x)); }

/** Returns rad2deg of the arc-cosine of the input */
auto acosd(auto x) noexcept { return rad2deg(std::acos(x)); }

/** Returns rad2deg of the arc-sine of the input */
auto asind(auto x) noexcept { return rad2deg(std::asin(x)); }

/** Returns rad2deg of the arc-tangent of the input */
auto atand(auto x) noexcept { return rad2deg(std::atan(x)); }

/** Returns rad2deg of the arc-tangent of inputs #T1/#T2  */
auto atan2d(auto y, auto x) noexcept { return rad2deg(std::atan2(y, x)); }

/** Conversion from Kayser wavenumber to Hz */
constexpr auto kaycm2freq(auto x) noexcept { return x * (100 * c); }

/** Conversion from Hz to Kayser wavenumber */
constexpr auto freq2kaycm(auto x) noexcept { return x / (100 * c); }

/** Conversion from Angular wavenumber to Hz */
constexpr auto angcm2freq(auto x) noexcept {
  return x * kaycm2freq(inv_two_pi);
}

/** Conversion from Hz to Angular wavenumber */
constexpr auto freq2angcm(auto x) noexcept {
  return x / kaycm2freq(inv_two_pi);
}

/** Conversion from Angular Hz to Hz */
constexpr auto angfreq2freq(auto x) noexcept { return x * inv_two_pi; }

/** Conversion from Hz to Angular Hz */
constexpr auto freq2angfreq(auto x) noexcept { return x * two_pi; }

/** Conversion from wavelength to Hz */
constexpr auto wavelen2freq(auto x) noexcept { return c / x; }

/** Conversion from Hz to wavelength */
constexpr auto freq2wavelen(auto x) noexcept { return c / x; }

/** Conversion from wavelength to Hz */
constexpr auto hz2ghz(auto x) noexcept { return x * 1e-9; }

/** Conversion from Hz to wavelength */
constexpr auto ghz2hz(auto x) noexcept { return x * 1e9; }

/** Conversion from Atm to Pa */
constexpr auto atm2pa(auto x) noexcept { return x * 101'325.0; }

/** Conversion from Pa to Atm */
constexpr auto pa2atm(auto x) noexcept { return x / 101'325.0; }

/** Conversion from bar to Pa */
constexpr auto bar2pa(auto x) noexcept { return x * 1e5; }

/** Conversion from Pa to bar */
constexpr auto pa2bar(auto x) noexcept { return x * 1e-5; }

/** Conversion from hPa to Pa */
constexpr auto hpa2pa(auto x) noexcept { return x * 1e2; }

/** Conversion from Pa to hPa */
constexpr auto pa2hpa(auto x) noexcept { return x * 1e-2; }

/** Conversion from hPa to bar */
constexpr auto hpa2bar(auto x) noexcept { return x * 1e-3; }

/** Conversion from bar to hPa */
constexpr auto bar2hpa(auto x) noexcept{ return x * 1e3; }

/** Conversion from Torr to Pa */
constexpr auto torr2pa(auto x) noexcept { return x * atm2pa(1.0 / 760.0); }

/** Conversion from Pa to Torr */
constexpr auto pa2torr(auto x) noexcept { return x / atm2pa(1.0 / 760.0); }

/** Conversion from MHz/Torr to Hz/Pa */
constexpr auto mhz_per_torr2hz_per_pa(auto x) noexcept {
  return x * pa2torr(1e6);
}

/** Conversion from C to K */
constexpr auto celsius2kelvin(auto x) noexcept { return x + 273.15; }

/** Conversion from K to C */
constexpr auto kelvin2celsius(auto x) noexcept { return x - 273.15; }

/** Conversion from cm-1 per molecule per cm^2 to Hz per molecule per m^2 **/
constexpr auto kaycm_per_cmsquared2hz_per_msquared(auto x) noexcept {
  return x * kaycm2freq(1e-4);
}

/** Conversion from Hz per molecule per m^2 to cm-1 per molecule per cm^2 **/
constexpr auto hz_per_msquared2kaycm_per_cmsquared(auto x) noexcept {
  return x * freq2kaycm(1e4);
}

/** Conversion from cm-1 per atmosphere to Hz per Pascal **/
constexpr auto kaycm_per_atm2hz_per_pa(auto x) noexcept {
  return x * kaycm2freq(pa2atm(1));
}

/** Conversion from Hz per Pascal to cm-1 per atmosphere **/
constexpr auto hz_per_pa2kaycm_per_atm(auto x) noexcept {
  return x * freq2kaycm(atm2pa(1));
}

/** Conversion from cm-1 to Joule **/
constexpr auto kaycm2joule(auto x) noexcept { return x * kaycm2freq(h); }

/** Conversion from MHz to Joule **/
constexpr auto hz2joule(auto x) noexcept { return x * h; }

/** Conversion from MHz to Joule **/
constexpr auto mhz2joule(auto x) noexcept { return hz2joule(x) * 1e6; }

/** Conversion from Kelvin to Joule **/
constexpr auto kelvin2joule(auto x) noexcept { return x * k; }

/** Conversion from Hz to Joule **/
constexpr auto joule2hz(auto x) noexcept { return x / h; }

/** Conversion from Joule to cm-1 **/
constexpr auto joule2kaycm(auto x) noexcept { return x / kaycm2freq(h); }

/** Conversion from Å to meter **/
constexpr auto angstrom2meter(auto x) noexcept { return x * 1e-10; }

/** Conversion from meter to Å **/
constexpr auto meter2angstrom(auto x) noexcept { return x * 1e10; }
};  // namespace Conversion

#endif
