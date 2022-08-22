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
template <class T>
constexpr auto deg2rad(T x) noexcept -> decltype(x * one_degree_in_radians) {
  return x * one_degree_in_radians;
}

/** Converts radians to degrees */
template <class T>
constexpr auto rad2deg(T x) noexcept -> decltype(x / one_degree_in_radians) {
  return x / one_degree_in_radians;
}

/** Returns the cosine of the deg2rad of the input */
template <class T>
auto cosd(T x) noexcept -> decltype(std::cos(deg2rad(x))) {
  return std::cos(deg2rad(x));
}

/** Returns the sine of the deg2rad of the input */
template <class T>
auto sind(T x) noexcept -> decltype(std::sin(deg2rad(x))) {
  return std::sin(deg2rad(x));
}

/** Returns the tangent of the deg2rad of the input */
template <class T>
auto tand(T x) noexcept -> decltype(std::tan(deg2rad(x))) {
  return std::tan(deg2rad(x));
}

/** Returns rad2deg of the arc-cosine of the input */
template <class T>
auto acosd(T x) noexcept -> decltype(rad2deg(std::acos(x))) {
  return rad2deg(std::acos(x));
}

/** Returns rad2deg of the arc-sine of the input */
template <class T>
auto asind(T x) noexcept -> decltype(rad2deg(std::asin(x))) {
  return rad2deg(std::asin(x));
}

/** Returns rad2deg of the arc-tangent of the input */
template <class T>
auto atand(T x) noexcept -> decltype(rad2deg(std::atan(x))) {
  return rad2deg(std::atan(x));
}

/** Returns rad2deg of the arc-tangent of inputs #T1/#T2  */
template <class T1, class T2>
auto atan2d(T1 y, T2 x) noexcept -> decltype(rad2deg(std::atan2(y, x))) {
  return rad2deg(std::atan2(y, x));
}

/** Conversion from Kayser wavenumber to Hz */
template <class T>
constexpr auto kaycm2freq(T x) noexcept -> decltype(x * (100 * c)) {
  return x * (100 * c);
}

/** Conversion from Hz to Kayser wavenumber */
template <class T>
constexpr auto freq2kaycm(T x) noexcept -> decltype(x / (100 * c)) {
  return x / (100 * c);
}

/** Conversion from Angular wavenumber to Hz */
template <class T>
constexpr auto angcm2freq(T x) noexcept -> decltype(kaycm2freq(inv_two_pi)) {
  return x * kaycm2freq(inv_two_pi);
}

/** Conversion from Hz to Angular wavenumber */
template <class T>
constexpr auto freq2angcm(T x) noexcept
    -> decltype(x / kaycm2freq(inv_two_pi)) {
  return x / kaycm2freq(inv_two_pi);
}

/** Conversion from Angular Hz to Hz */
template <class T>
constexpr auto angfreq2freq(T x) noexcept -> decltype(x * inv_two_pi) {
  return x * inv_two_pi;
}

/** Conversion from Hz to Angular Hz */
template <class T>
constexpr auto freq2angfreq(T x) noexcept -> decltype(x * two_pi) {
  return x * two_pi;
}

/** Conversion from wavelength to Hz */
template <class T>
constexpr auto wavelen2freq(T x) noexcept -> decltype(c / x) {
  return c / x;
}

/** Conversion from Hz to wavelength */
template <class T>
constexpr auto freq2wavelen(T x) noexcept -> decltype(c / x) {
  return c / x;
}

/** Conversion from wavelength to Hz */
template <class T>
constexpr auto hz2ghz(T x) noexcept -> decltype(x * 1e-9) {
  return x * 1e-9;
}

/** Conversion from Hz to wavelength */
template <class T>
constexpr auto ghz2hz(T x) noexcept -> decltype(x * 1e9) {
  return x * 1e9;
}

/** Conversion from Atm to Pa */
template <class T>
constexpr auto atm2pa(T x) noexcept -> decltype(x * 101'325.0) {
  return x * 101'325.0;
}

/** Conversion from Pa to Atm */
template <class T>
constexpr auto pa2atm(T x) noexcept -> decltype(x / 101'325.0) {
  return x / 101'325.0;
}

/** Conversion from bar to Pa */
template <class T>
constexpr auto bar2pa(T x) noexcept -> decltype(x * 1e5) {
  return x * 1e5;
}

/** Conversion from Pa to bar */
template <class T>
constexpr auto pa2bar(T x) noexcept -> decltype(x * 1e-5) {
  return x * 1e-5;
}

/** Conversion from Torr to Pa */
template <class T>
constexpr auto torr2pa(T x) noexcept -> decltype(x * atm2pa(1.0 / 760.0)) {
  return x * atm2pa(1.0 / 760.0);
}

/** Conversion from Pa to Torr */
template <class T>
constexpr auto pa2torr(T x) noexcept -> decltype(x / atm2pa(1.0 / 760.0)) {
  return x / atm2pa(1.0 / 760.0);
}

/** Conversion from MHz/Torr to Hz/Pa */
template <class T>
constexpr auto mhz_per_torr2hz_per_pa(T x) noexcept
    -> decltype(x * pa2torr(1e6)) {
  return x * pa2torr(1e6);
}

/** Conversion from C to K */
template <class T>
constexpr auto celsius2kelvin(T x) noexcept -> decltype(x + 273.15) {
  return x + 273.15;
}

/** Conversion from K to C */
template <class T>
constexpr auto kelvin2celsius(T x) noexcept -> decltype(x - 273.15) {
  return x - 273.15;
}

/** Conversion from cm-1 per molecule per cm^2 to Hz per molecule per m^2 **/
template <class T>
constexpr auto kaycm_per_cmsquared2hz_per_msquared(T x) noexcept
    -> decltype(x * kaycm2freq(1e-4)) {
  return x * kaycm2freq(1e-4);
}

/** Conversion from Hz per molecule per m^2 to cm-1 per molecule per cm^2 **/
template <class T>
constexpr auto hz_per_msquared2kaycm_per_cmsquared(T x) noexcept
    -> decltype(x * freq2kaycm(1e4)) {
  return x * freq2kaycm(1e4);
}

/** Conversion from cm-1 per atmosphere to Hz per Pascal **/
template <class T>
constexpr auto kaycm_per_atm2hz_per_pa(T x) noexcept
    -> decltype(x * kaycm2freq(pa2atm(1))) {
  return x * kaycm2freq(pa2atm(1));
}

/** Conversion from Hz per Pascal to cm-1 per atmosphere **/
template <class T>
constexpr auto hz_per_pa2kaycm_per_atm(T x) noexcept
    -> decltype(x * freq2kaycm(atm2pa(1))) {
  return x * freq2kaycm(atm2pa(1));
}

/** Conversion from cm-1 to Joule **/
template <class T>
constexpr auto kaycm2joule(T x) noexcept -> decltype(x * kaycm2freq(h)) {
  return x * kaycm2freq(h);
}

/** Conversion from MHz to Joule **/
template <class T>
constexpr auto hz2joule(T x) noexcept -> decltype(x * h) {
  return x * h;
}

/** Conversion from MHz to Joule **/
template <class T>
constexpr auto mhz2joule(T x) noexcept -> decltype(hz2joule(x) * 1e6) {
  return hz2joule(x) * 1e6;
}

/** Conversion from Kelvin to Joule **/
template <class T>
constexpr auto kelvin2joule(T x) noexcept -> decltype(x * k) {
  return x * k;
}

/** Conversion from Hz to Joule **/
template <class T>
constexpr auto joule2hz(T x) noexcept -> decltype(x / h) {
  return x / h;
}

/** Conversion from Joule to cm-1 **/
template <class T>
constexpr auto joule2kaycm(T x) noexcept -> decltype(x / kaycm2freq(h)) {
  return x / kaycm2freq(h);
}

/** Conversion from Å to meter **/
template <class T>
constexpr auto angstrom2meter(T x) noexcept -> decltype(x * 1e-10) {
  return x * 1e-10;
}

/** Conversion from meter to Å **/
template <class T>
constexpr auto meter2angstrom(T x) noexcept -> decltype(x * 1e10) {
  return x * 1e10;
}
};  // namespace Conversion

#endif
