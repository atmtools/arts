/* Copyright (C) 2020
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
 * @file   predefined_absorption_models.cc
 * @author Richard Larsson
 * @date   2020-01-29
 * 
 * @brief  Full absorption models of various kinds
 */

#include "predefined_absorption_models.h"
#include "lin_alg.h"
#include "linescaling.h"
#include "wigner_functions.h"


constexpr auto necs2020 = 38;
constexpr auto nlines_mpm2020 = necs2020 + 6;


constexpr LineShape::SingleSpeciesModel init_mpm2020_slsm(Numeric g00,
                                                          Numeric y0,
                                                          Numeric y1,
                                                          Numeric g0,
                                                          Numeric g1,
                                                          Numeric dv0,
                                                          Numeric dv1,
                                                          Numeric x)
{
  return LineShape::SingleSpeciesModel(
  {LineShape::TemperatureModel::T1, g00, x, NAN, NAN},
  {LineShape::TemperatureModel::None, NAN, NAN, NAN, NAN},
  {LineShape::TemperatureModel::None, NAN, NAN, NAN, NAN},
  {LineShape::TemperatureModel::None, NAN, NAN, NAN, NAN},
  {LineShape::TemperatureModel::None, NAN, NAN, NAN, NAN},
  {LineShape::TemperatureModel::None, NAN, NAN, NAN, NAN},
  {LineShape::TemperatureModel::T4, y0, y1, x, NAN},
  {LineShape::TemperatureModel::T4, g0, g1, 2*x, NAN},
  {LineShape::TemperatureModel::T4, dv0, dv1, 2*x, NAN});
}


constexpr std::array<LineShape::SingleSpeciesModel, nlines_mpm2020> init_mpm2020_lsm()
{
  // Pressure broadening [1/Pa] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> g00 =
    {1.685E+4, 1.703E+4, 1.513E+4, 1.495E+4, 1.433E+4, 
      1.408E+4, 1.353E+4, 1.353E+4, 1.303E+4, 1.319E+4, 
      1.262E+4, 1.265E+4, 1.238E+4, 1.217E+4, 1.207E+4, 
      1.207E+4, 1.137E+4, 1.137E+4, 1.101E+4, 1.101E+4, 
      1.037E+4, 1.038E+4, 9.96E+3, 9.96E+3, 9.55E+3, 
      9.55E+3, 9.06E+3, 9.06E+3, 8.58E+3, 8.58E+3, 
      8.11E+3, 8.11E+3, 7.64E+3, 7.64E+3, 7.17E+3, 
      7.17E+3, 6.69E+3, 6.69E+3, 1.64E+4, 1.64E+4, 
      1.60E+4, 1.60E+4, 1.62E+4, 1.47E+4,};
  
  // First order line mixing first coefficient [1/Pa] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> y0 = 
    {-4.1E-7, 0.00000277, -0.00000372, 0.00000559, -0.00000573, 
      0.00000618, -0.00000366, 0.00000278, -8.9E-7, -2.1E-7, 
      6.0E-7, -0.00000152, 0.00000216, -0.00000293, 0.00000373, 
      -0.00000436, 0.00000491, -0.00000542, 0.00000571, -0.00000613, 
      0.00000636, -0.00000670, 0.00000690, -0.00000718, 0.00000740, 
      -0.00000763, 0.00000788, -0.00000807, 0.00000834, -0.00000849, 
      0.00000876, -0.00000887, 0.00000915, -0.00000922, 0.00000950, 
      -0.00000955, 0.00000987, -0.00000988, 0.00000, 0.00000, 
      0.00000, 0.00000, 0.00000, 0.00000,};
  
  // First order line mixing second coefficient [1/Pa] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> y1 =
    {0.00000, 0.00000124, -2E-8, 8E-8, 4.5E-7, 
      -9.3E-7, 0.00000264, -0.00000351, 0.00000359, -0.00000416, 
      0.00000326, -0.00000353, 0.00000484, -0.00000503, 0.00000579, 
      -0.00000590, 0.00000616, -0.00000619, 0.00000611, -0.00000609, 
      0.00000574, -0.00000568, 0.00000574, -0.00000566, 0.0000060, 
      -0.0000059, 0.0000063, -0.0000062, 0.0000064, -0.0000063, 
      0.0000065, -0.0000064, 0.0000065, -0.0000064, 0.0000065, 
      -0.0000064, 0.0000064, -0.0000062, 0.00000, 0.00000, 
      0.00000, 0.00000, 0.00000, 0.00000, };
  
  // Second order line mixing strength adjustment first coefficient [1/Pa^2] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> g0 =
    {-6.95E-14, -9.0E-12, -1.03E-11, -2.39E-11, -1.72E-11, 
      -1.71E-11, 2.8E-12, 1.50E-11, 1.32E-11, 1.70E-11, 
      8.7E-12, 6.9E-12, 8.3E-12, 6.7E-12, 7E-13, 
      1.6E-12, -2.1E-12, -6.6E-12, -9.5E-12, -1.15E-11, 
      -1.18E-11, -1.40E-11, -1.73E-11, -1.86E-11, -2.17E-11, 
      -2.27E-11, -2.34E-11, -2.42E-11, -2.66E-11, -2.72E-11, 
      -3.01E-11, -3.04E-11, -3.34E-11, -3.33E-11, -3.61E-11, 
      -3.58E-11, -3.48E-11, -3.44E-11, 0E-10, 0E-10, 
      0E-10, 0E-10, 0E-10, 0E-10,};
  
  // Second order line mixing strength adjustment second coefficient [1/Pa^2] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> g1 =
    {0E-10, -4.5E-12, 7E-13, 3.3E-12, 8.1E-12, 
      1.62E-11, 1.79E-11, 2.25E-11, 5.4E-12, 3E-13, 
      4E-14, -4.7E-12, -3.4E-12, -7.1E-12, -1.80E-11, 
      -2.10E-11, -2.85E-11, -3.23E-11, -3.63E-11, -3.80E-11, 
      -3.78E-11, -3.87E-11, -3.92E-11, -3.94E-11, -4.24E-11, 
      -4.22E-11, -4.65E-11, -4.6E-11, -5.1E-11, -5.0E-11, 
      -5.5E-11, -5.4E-11, -5.8E-11, -5.6E-11, -6.2E-11, 
      -5.9E-11, -6.8E-11, -6.5E-11, 0E-10, 0E-10, 
      0E-10, 0E-10, 0E-10, 0E-10, };
  
  // Second order line mixing frequency adjustment first coefficient [Hz/Pa^2] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> dv0 =
    {-0.000028, 0.000597, -0.00195, 0.0032, -0.00475, 
      0.00541, -0.00232, 0.00154, 0.00007, -0.00084, 
      -0.00025, -0.00014, -0.00004, -0.00020, 0.0005, 
      -0.00066, 0.00072, -0.0008, 0.00064, -0.00070, 
      0.00056, -0.00060, 0.00047, -0.00049, 0.00040, 
      -0.00041, 0.00036, -0.00037, 0.00033, -0.00034, 
      0.00032, -0.00032, 0.00030, -0.00030, 0.00028, 
      -0.00029, 0.00029, -0.00029, 0.0, 0.0, 
      0.0, 0.0, 0.0, 0.0, };
  
  // Second order line mixing frequency adjustment second coefficient [Hz/Pa^2] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> dv1 =
    {-0.000039, 0.0009, -0.0012, 0.0016, -0.0027, 
      0.0029, 0.0006, -0.0015, 0.0010, -0.0014, 
      -0.0013, 0.0013, 0.0004, -0.0005, 0.0010, 
      -0.0010, 0.0010, -0.0011, 0.0008, -0.0009, 
      0.0003, -0.0003, 0.00009, -0.00009, 0.00017, 
      -0.00016, 0.00024, -0.00023, 0.00024, -0.00024, 
      0.00024, -0.00020, 0.00017, -0.00016, 0.00013, 
      -0.00012, 0.00005, -0.00004, 0.0, 0.0, 
      0.0, 0.0, 0.0, 0.0, };
  
  // Temperature scaling exponent [-]
  constexpr Numeric x = 0.754;
  
  return {init_mpm2020_slsm(g00[0], y0[0], y1[0], g0[0], g1[0], dv0[0], dv1[0], x),
          init_mpm2020_slsm(g00[1], y0[1], y1[1], g0[1], g1[1], dv0[1], dv1[1], x),
          init_mpm2020_slsm(g00[2], y0[2], y1[2], g0[2], g1[2], dv0[2], dv1[2], x),
          init_mpm2020_slsm(g00[3], y0[3], y1[3], g0[3], g1[3], dv0[3], dv1[3], x),
          init_mpm2020_slsm(g00[4], y0[4], y1[4], g0[4], g1[4], dv0[4], dv1[4], x),
          init_mpm2020_slsm(g00[5], y0[5], y1[5], g0[5], g1[5], dv0[5], dv1[5], x),
          init_mpm2020_slsm(g00[6], y0[6], y1[6], g0[6], g1[6], dv0[6], dv1[6], x),
          init_mpm2020_slsm(g00[7], y0[7], y1[7], g0[7], g1[7], dv0[7], dv1[7], x),
          init_mpm2020_slsm(g00[8], y0[8], y1[8], g0[8], g1[8], dv0[8], dv1[8], x),
          init_mpm2020_slsm(g00[9], y0[9], y1[9], g0[9], g1[9], dv0[9], dv1[9], x),
          init_mpm2020_slsm(g00[10], y0[10], y1[10], g0[10], g1[10], dv0[10], dv1[10], x),
          init_mpm2020_slsm(g00[11], y0[11], y1[11], g0[11], g1[11], dv0[11], dv1[11], x),
          init_mpm2020_slsm(g00[12], y0[12], y1[12], g0[12], g1[12], dv0[12], dv1[12], x),
          init_mpm2020_slsm(g00[13], y0[13], y1[13], g0[13], g1[13], dv0[13], dv1[13], x),
          init_mpm2020_slsm(g00[14], y0[14], y1[14], g0[14], g1[14], dv0[14], dv1[14], x),
          init_mpm2020_slsm(g00[15], y0[15], y1[15], g0[15], g1[15], dv0[15], dv1[15], x),
          init_mpm2020_slsm(g00[16], y0[16], y1[16], g0[16], g1[16], dv0[16], dv1[16], x),
          init_mpm2020_slsm(g00[17], y0[17], y1[17], g0[17], g1[17], dv0[17], dv1[17], x),
          init_mpm2020_slsm(g00[18], y0[18], y1[18], g0[18], g1[18], dv0[18], dv1[18], x),
          init_mpm2020_slsm(g00[19], y0[19], y1[19], g0[19], g1[19], dv0[19], dv1[19], x),
          init_mpm2020_slsm(g00[20], y0[20], y1[20], g0[20], g1[20], dv0[20], dv1[20], x),
          init_mpm2020_slsm(g00[21], y0[21], y1[21], g0[21], g1[21], dv0[21], dv1[21], x),
          init_mpm2020_slsm(g00[22], y0[22], y1[22], g0[22], g1[22], dv0[22], dv1[22], x),
          init_mpm2020_slsm(g00[23], y0[23], y1[23], g0[23], g1[23], dv0[23], dv1[23], x),
          init_mpm2020_slsm(g00[24], y0[24], y1[24], g0[24], g1[24], dv0[24], dv1[24], x),
          init_mpm2020_slsm(g00[25], y0[25], y1[25], g0[25], g1[25], dv0[25], dv1[25], x),
          init_mpm2020_slsm(g00[26], y0[26], y1[26], g0[26], g1[26], dv0[26], dv1[26], x),
          init_mpm2020_slsm(g00[27], y0[27], y1[27], g0[27], g1[27], dv0[27], dv1[27], x),
          init_mpm2020_slsm(g00[28], y0[28], y1[28], g0[28], g1[28], dv0[28], dv1[28], x),
          init_mpm2020_slsm(g00[29], y0[29], y1[29], g0[29], g1[29], dv0[29], dv1[29], x),
          init_mpm2020_slsm(g00[30], y0[30], y1[30], g0[30], g1[30], dv0[30], dv1[30], x),
          init_mpm2020_slsm(g00[31], y0[31], y1[31], g0[31], g1[31], dv0[31], dv1[31], x),
          init_mpm2020_slsm(g00[32], y0[32], y1[32], g0[32], g1[32], dv0[32], dv1[32], x),
          init_mpm2020_slsm(g00[33], y0[33], y1[33], g0[33], g1[33], dv0[33], dv1[33], x),
          init_mpm2020_slsm(g00[34], y0[34], y1[34], g0[34], g1[34], dv0[34], dv1[34], x),
          init_mpm2020_slsm(g00[35], y0[35], y1[35], g0[35], g1[35], dv0[35], dv1[35], x),
          init_mpm2020_slsm(g00[36], y0[36], y1[36], g0[36], g1[36], dv0[36], dv1[36], x),
          init_mpm2020_slsm(g00[37], y0[37], y1[37], g0[37], g1[37], dv0[37], dv1[37], x),
          init_mpm2020_slsm(g00[38], y0[38], y1[38], g0[38], g1[38], dv0[38], dv1[38], x),
          init_mpm2020_slsm(g00[39], y0[39], y1[39], g0[39], g1[39], dv0[39], dv1[39], x),
          init_mpm2020_slsm(g00[40], y0[40], y1[40], g0[40], g1[40], dv0[40], dv1[40], x),
          init_mpm2020_slsm(g00[41], y0[41], y1[41], g0[41], g1[41], dv0[41], dv1[41], x),
          init_mpm2020_slsm(g00[42], y0[42], y1[42], g0[42], g1[42], dv0[42], dv1[42], x),
          init_mpm2020_slsm(g00[43], y0[43], y1[43], g0[43], g1[43], dv0[43], dv1[43], x),};
}


constexpr QuantumIdentifier init_mpm2020_qid(Index species, Index isot, Rational Jup, Rational Jlow, Rational Nup, Rational Nlow)
{
  return QuantumIdentifier(species, isot, QuantumNumbers(Jup, Nup, 0), QuantumNumbers(Jlow, Nlow, 0));
}


constexpr std::array<QuantumIdentifier, nlines_mpm2020> init_mpm2020_qids(const Index& species, const Index& isot)
{
  // N of upper level
  constexpr std::array<Index, nlines_mpm2020> Np = {
    1, 1, 3, 3, 5, 5, 7, 7, 9, 9,
    11, 11, 13, 13, 15, 15, 17, 17,
    19, 19,  21, 21, 23, 23, 25, 25,
    27, 27, 29, 29, 31, 31, 33, 33,
    35, 35, 37, 37, 1, 1, 1, 3, 3, 3};
  
  // N of lower level
  constexpr std::array<Index, nlines_mpm2020> Npp = {
    1, 1, 3, 3, 5, 5, 7, 7, 9, 9,
    11, 11, 13, 13, 15, 15, 17, 17,
    19, 19,  21, 21, 23, 23, 25, 25,
    27, 27, 29, 29, 31, 31, 33, 33,
    35, 35, 37, 37, 3, 3, 3, 5, 5, 5};
  
  // J of upper level
  constexpr std::array<Index, nlines_mpm2020> Jp = {
    1, 1, 3, 3, 5, 5, 7, 7, 9, 9,
    11, 11, 13, 13, 15, 15, 17, 17,
    19, 19,  21, 21, 23, 23, 25, 25,
    27, 27, 29, 29, 31, 31, 33, 33,
    35, 35, 37, 37, 1, 2, 2, 3, 4, 4};
  
  // J of lower level
  constexpr std::array<Index, nlines_mpm2020> Jpp = {
    0, 2, 2, 4, 4, 6, 6, 8, 8, 10,
    10, 12, 12, 14, 14, 16, 16, 18,
    18, 20,  20, 22, 22, 24, 24, 26,
    26, 28, 28, 30, 30, 32, 32, 34,
    34, 36, 36, 38, 2, 2, 3, 4, 4, 5};
  
  return {QuantumIdentifier(species, isot, QuantumNumbers(Jp[0], Np[0], 0), QuantumNumbers(Jpp[0], Npp[0], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[1], Np[1], 0), QuantumNumbers(Jpp[1], Npp[1], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[2], Np[2], 0), QuantumNumbers(Jpp[2], Npp[2], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[3], Np[3], 0), QuantumNumbers(Jpp[3], Npp[3], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[4], Np[4], 0), QuantumNumbers(Jpp[4], Npp[4], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[5], Np[5], 0), QuantumNumbers(Jpp[5], Npp[5], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[6], Np[6], 0), QuantumNumbers(Jpp[6], Npp[6], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[7], Np[7], 0), QuantumNumbers(Jpp[7], Npp[7], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[8], Np[8], 0), QuantumNumbers(Jpp[8], Npp[8], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[9], Np[9], 0), QuantumNumbers(Jpp[9], Npp[9], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[10], Np[10], 0), QuantumNumbers(Jpp[10], Npp[10], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[11], Np[11], 0), QuantumNumbers(Jpp[11], Npp[11], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[12], Np[12], 0), QuantumNumbers(Jpp[12], Npp[12], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[13], Np[13], 0), QuantumNumbers(Jpp[13], Npp[13], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[14], Np[14], 0), QuantumNumbers(Jpp[14], Npp[14], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[15], Np[15], 0), QuantumNumbers(Jpp[15], Npp[15], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[16], Np[16], 0), QuantumNumbers(Jpp[16], Npp[16], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[17], Np[17], 0), QuantumNumbers(Jpp[17], Npp[17], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[18], Np[18], 0), QuantumNumbers(Jpp[18], Npp[18], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[19], Np[19], 0), QuantumNumbers(Jpp[19], Npp[19], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[20], Np[20], 0), QuantumNumbers(Jpp[20], Npp[20], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[21], Np[21], 0), QuantumNumbers(Jpp[21], Npp[21], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[22], Np[22], 0), QuantumNumbers(Jpp[22], Npp[22], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[23], Np[23], 0), QuantumNumbers(Jpp[23], Npp[23], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[24], Np[24], 0), QuantumNumbers(Jpp[24], Npp[24], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[25], Np[25], 0), QuantumNumbers(Jpp[25], Npp[25], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[26], Np[26], 0), QuantumNumbers(Jpp[26], Npp[26], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[27], Np[27], 0), QuantumNumbers(Jpp[27], Npp[27], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[28], Np[28], 0), QuantumNumbers(Jpp[28], Npp[28], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[29], Np[29], 0), QuantumNumbers(Jpp[29], Npp[29], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[30], Np[30], 0), QuantumNumbers(Jpp[30], Npp[30], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[31], Np[31], 0), QuantumNumbers(Jpp[31], Npp[31], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[32], Np[32], 0), QuantumNumbers(Jpp[32], Npp[32], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[33], Np[33], 0), QuantumNumbers(Jpp[33], Npp[33], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[34], Np[34], 0), QuantumNumbers(Jpp[34], Npp[34], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[35], Np[35], 0), QuantumNumbers(Jpp[35], Npp[35], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[36], Np[36], 0), QuantumNumbers(Jpp[36], Npp[36], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[37], Np[37], 0), QuantumNumbers(Jpp[37], Npp[37], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[38], Np[38], 0), QuantumNumbers(Jpp[38], Npp[38], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[39], Np[39], 0), QuantumNumbers(Jpp[39], Npp[39], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[40], Np[40], 0), QuantumNumbers(Jpp[40], Npp[40], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[41], Np[41], 0), QuantumNumbers(Jpp[41], Npp[41], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[42], Np[42], 0), QuantumNumbers(Jpp[42], Npp[42], 0)),
          QuantumIdentifier(species, isot, QuantumNumbers(Jp[43], Np[43], 0), QuantumNumbers(Jpp[43], Npp[43], 0)), };
}


void Absorption::PredefinedModel::makarov2020_o2_lines_mpm(Matrix& xsec,
                                                           ArrayOfMatrix& dxsec,
                                                           const Vector& f,
                                                           const Vector& p,
                                                           const Vector& t,
                                                           const Vector& water_vmr,
                                                           const ArrayOfRetrievalQuantity& jacs,
                                                           const ArrayOfIndex& jacs_pos)
{
  using Constant::pi;
  using Constant::sqrt_pi;
  using Constant::inv_sqrt_pi;
  using Constant::pow2;
  using Constant::pow3;
  
  // Central frequency [Hz]
  constexpr std::array<Numeric, nlines_mpm2020> f0 = {
    1.18750334E+11, 5.6264774E+10, 6.2486253E+10, 5.8446588E+10, 6.0306056E+10, 
    5.9590983E+10, 5.9164204E+10, 6.0434778E+10, 5.8323877E+10, 6.1150562E+10, 
    5.7612486E+10, 6.1800158E+10, 5.6968211E+10, 6.2411220E+10, 5.6363399E+10, 
    6.2997984E+10, 5.5783815E+10, 6.3568526E+10, 5.5221384E+10, 6.4127775E+10, 
    5.4671180E+10, 6.4678910E+10, 5.4130025E+10, 6.5224078E+10, 5.3595775E+10, 
    6.5764779E+10, 5.3066934E+10, 6.6302096E+10, 5.2542418E+10, 6.6836834E+10, 
    5.2021429E+10, 6.7369601E+10, 5.1503360E+10, 6.7900868E+10, 5.0987745E+10, 
    6.8431006E+10, 5.0474214E+10, 6.8960312E+10, 3.68498246E+11, 4.24763020E+11, 
    4.87249273E+11, 7.15392902E+11, 7.73839490E+11, 8.34145546E+11, };
  
  // Intensity [1 / Pa] (rounded to 10 digits because at most 9 digits exist in f0)
  constexpr std::array<Numeric, nlines_mpm2020> intens = {
    1.591521878E-21, 1.941172240E-21, 4.834543970E-21, 4.959264029E-21, 7.010386457E-21, 
    7.051673348E-21, 8.085012578E-21, 8.108262250E-21, 8.145673278E-21, 8.149757320E-21, 
    7.396406085E-21, 7.401923754E-21, 6.162286575E-21, 6.168475265E-21, 4.749226167E-21, 
    4.754435107E-21, 3.405982896E-21, 3.408455562E-21, 2.282498656E-21, 2.283934341E-21, 
    1.432692459E-21, 1.433513473E-21, 8.439995690E-22, 8.443521837E-22, 4.672706507E-22, 
    4.676049313E-22, 2.435008301E-22, 2.437304596E-22, 1.195038747E-22, 1.196873412E-22, 
    5.532759045E-23, 5.537261239E-23, 2.416832398E-23, 2.418989865E-23, 9.969285671E-24, 
    9.977543709E-24, 3.882541154E-24, 3.888101811E-24, 3.676253816E-23, 3.017524005E-22, 
    9.792882227E-23, 2.756166168E-23, 1.486462215E-22, 4.411918954E-23, };
  
  // Temperature intensity modifier
  constexpr std::array<Numeric, nlines_mpm2020> a2 = {
    0.01, 0.014, 0.083, 0.083, 0.207, 0.207, 0.387, 0.386,
    0.621, 0.621, 0.910, 0.910, 1.255, 1.255, 1.654, 1.654,
    2.109, 2.108, 2.618, 2.617, 3.182, 3.181, 3.800, 3.800,
    4.474, 4.473, 5.201, 5.200, 5.983, 5.982, 6.819, 6.818,
    7.709, 7.708, 8.653, 8.652, 9.651, 9.650, 0.048, 0.044,
    0.049, 0.145, 0.141, 0.145};
  
  // Line shape model in SI units
  constexpr auto lsm = init_mpm2020_lsm();
  
  // Reference temperature [K]
  constexpr Numeric t0 = 300;
  
  // FIXME: Can be constexpr when SpeciesTag is constexpr 
  auto species = SpeciesTag("O2-66");
  const std::array<QuantumIdentifier, nlines_mpm2020> qids = 
    init_mpm2020_qids(species.Species(), species.Isotopologue());
  
  // Model setting
  const bool do_temp_deriv = do_temperature_jacobian(jacs);
  
  // Per pressure level
  #pragma omp parallel for if (not arts_omp_in_parallel() and p.nelem() >= arts_omp_get_max_threads()) schedule(guided) 
  for (Index ip=0; ip<p.nelem(); ip++) {
    const Numeric theta = t0 / t[ip];
    const Numeric theta_m1 = theta - 1;
    const Numeric theta_3 = pow3(theta);
    const Numeric GD_div_F0 = Linefunctions::DopplerConstant(t[ip], species.SpeciesMass());
    
    for (Index i=0; i<nlines_mpm2020; i++) {
      const Numeric invGD = 1 / (GD_div_F0 * f0[i]);
      const Numeric fac = sqrt_pi * invGD;
      const Numeric ST = theta_3 * p[ip] * intens[i] * std::exp(-a2[i] * theta_m1);
      const Numeric G0 = (1 + 0.1*water_vmr[ip]) * p[ip] * lsm[i].compute(t[ip], t0, LineShape::Variable::G0);
      const Numeric Y = p[ip] * lsm[i].compute(t[ip], t0, LineShape::Variable::Y);
      const Numeric G = pow2( p[ip]) * lsm[i].compute(t[ip], t0, LineShape::Variable::G);
      const Numeric DV = pow2(p[ip]) * lsm[i].compute(t[ip], t0, LineShape::Variable::DV);
      
      const Numeric dinvGD_dT = do_temp_deriv ? - invGD * Linefunctions::dDopplerConstant_dT(t[ip], GD_div_F0) : 0;
      const Numeric dST_dT = do_temp_deriv ? (a2[i]*t0 - 3*t[ip]) / pow2(t[ip]) * ST : 0;
      const Numeric dG0_dT = do_temp_deriv ? (1 + 0.1*water_vmr[ip]) * p[ip] * lsm[i].compute_dT(t[ip], t0, LineShape::Variable::G0) : 0;
      const Numeric dY_dT = do_temp_deriv ? p[ip] * lsm[i].compute_dT(t[ip], t0, LineShape::Variable::Y) : 0;
      const Numeric dG_dT = do_temp_deriv ? pow2(p[ip]) * lsm[i].compute_dT(t[ip], t0, LineShape::Variable::G) : 0;
      const Numeric dDV_dT = do_temp_deriv ? pow2(p[ip]) * lsm[i].compute_dT(t[ip], t0, LineShape::Variable::DV) : 0;
      
      for (Index j=0; j<f.nelem(); j++) {
        const Complex z = Complex(f0[i] + DV - f[j], G0) * invGD;
        const Complex Fv = fac * Faddeeva::w(z);
        const Complex Flm = 1 / Complex(G0, f[j] + f0[i] + DV);
        
        const Complex abs = std::real(
          /* around line center */
          Complex(1 + G, Y) * Fv +
          /* mirrored line far from line center */
          Complex(1 + G, -Y) * Flm);
        
        xsec(j, ip) += ST * pow2(f[j]) * abs.real();
        
        if (jacs_pos.nelem()) {
          const Complex dw = 2 * (Complex(0, fac * inv_sqrt_pi) - z * Fv);
          const Complex dm = - pi * pow2(Flm);
          
          for (Index iq=0; iq<jacs_pos.nelem(); iq++) {
            const auto& deriv = jacs[jacs_pos[iq]];
            
            if (deriv == JacPropMatType::Temperature) {
              const Complex dFv = dw * (invGD * Complex(dDV_dT, dG0_dT) - dinvGD_dT) + Fv * dinvGD_dT;
              const Complex dFlm = dm * Complex(dG0_dT, dDV_dT);
              dxsec[iq](j, ip) += pow2(f[j]) * (ST * std::real(
                /* around line center */
                Complex(1 + G, Y) * dFv + Complex(dG_dT, dY_dT) * Fv +
                /* mirrored line far from line center */
                Complex(1 + G, -Y) * dFlm + Complex(G, -dY_dT) * Flm) + abs.real() * dST_dT);
            } else if (is_frequency_parameter(deriv)) {
              const Complex dFv = - dw * invGD;
              const Complex dFlm = Complex(0, 1) * dm;
              dxsec[iq](j, ip) += ST * (pow2(f[j]) * std::real(
                /* around line center */
                Complex(1 + G, Y) * dFv +
                /* mirrored line far from line center */
                Complex(1 + G, -Y) * dFlm) + 2 * abs.real() * f[j]);
            } else if (deriv.QuantumIdentity().In(qids[i])) {
              if (deriv == JacPropMatType::LineShapeG0X0) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  /* around line center */
                  Complex(1 + G, Y) * Complex(0, 1) * dw * invGD +
                  /* mirrored line far from line center */
                  Complex(1 + G, -Y) * dm) * 
                    lsm[i].compute_dX0(t[ip], t0, LineShape::Variable::G0);
              } else if (deriv == JacPropMatType::LineShapeG0X1) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  /* around line center */
                  Complex(1 + G, Y) * Complex(0, 1) * dw * invGD +
                  /* mirrored line far from line center */
                  Complex(1 + G, -Y) * dm) * 
                    lsm[i].compute_dX1(t[ip], t0, LineShape::Variable::DV);
              } else if (deriv == JacPropMatType::LineShapeDVX0) {
                const Complex dFv = dw * invGD;
                const Complex dFlm = Complex(0, 1) * dm;
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  /* around line center */
                  Complex(1 + G, Y) * dFv +
                  /* mirrored line far from line center */
                  Complex(1 + G, -Y) * dFlm) * 
                    lsm[i].compute_dX0(t[ip], t0, LineShape::Variable::DV);;
              } else if (deriv == JacPropMatType::LineShapeDVX1) {
                const Complex dFv = dw * invGD;
                const Complex dFlm = Complex(0, 1) * dm;
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  /* around line center */
                  Complex(1 + G, Y) * dFv +
                  /* mirrored line far from line center */
                  Complex(1 + G, -Y) * dFlm) * 
                    lsm[i].compute_dX1(t[ip], t0, LineShape::Variable::DV);;
              } else if (deriv == JacPropMatType::LineShapeDVX2) {
                const Complex dFv = dw * invGD;
                const Complex dFlm = Complex(0, 1) * dm;
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  /* around line center */
                  Complex(1 + G, Y) * dFv +
                  /* mirrored line far from line center */
                  Complex(1 + G, -Y) * dFlm) * 
                    lsm[i].compute_dX2(t[ip], t0, LineShape::Variable::DV);
              } else if (deriv == JacPropMatType::LineShapeGX0) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  /* around line center */
                  Fv +
                  /* mirrored line far from line center */
                  Flm) * 
                    lsm[i].compute_dX0(t[ip], t0, LineShape::Variable::G);
              } else if (deriv == JacPropMatType::LineShapeYX0) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  /* around line center */
                  Fv +
                  /* mirrored line far from line center */
                  Flm) * 
                    lsm[i].compute_dX0(t[ip], t0, LineShape::Variable::Y);
              } else if (deriv == JacPropMatType::LineShapeGX1) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  /* around line center */
                  Fv -
                  /* mirrored line far from line center */
                  Flm) * 
                    lsm[i].compute_dX1(t[ip], t0, LineShape::Variable::G);
              } else if (deriv == JacPropMatType::LineShapeYX1) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  /* around line center */
                  Fv +
                  /* mirrored line far from line center */
                  Flm) * 
                    lsm[i].compute_dX1(t[ip], t0, LineShape::Variable::Y);
              } else if (deriv == JacPropMatType::LineShapeGX2) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  /* around line center */
                  Fv -
                  /* mirrored line far from line center */
                  Flm) * 
                    lsm[i].compute_dX2(t[ip], t0, LineShape::Variable::G);
              } else if (deriv == JacPropMatType::LineShapeYX2) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  /* around line center */
                  Fv -
                  /* mirrored line far from line center */
                  Flm) * 
                    lsm[i].compute_dX2(t[ip], t0, LineShape::Variable::Y);
              } else if (deriv == JacPropMatType::LineCenter) {
                const Complex dFv = Fv / f0[i] - dw * invGD + dw * z / f0[i];
                const Complex dFlm = Complex(0, 1) * dm;
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  /* around line center */
                  Complex(1 + G, Y) * dFv +
                  /* mirrored line far from line center */
                  Complex(1 + G, -Y) * dFlm);
              } else if (deriv == JacPropMatType::LineStrength) {
                dxsec[iq](j, ip) += theta_3 * p[ip] * std::exp(-a2[i] * theta_m1) * pow2(f[j]) * abs.real();
              } else if (deriv == JacPropMatType::LineSpecialParameter1) {
                dxsec[iq](j, ip) += -theta_m1 * ST * pow2(f[j]) * abs.real();
              }
            }
          }
        }
      }
    }
  }
}


void normalize_relaxation_matrix(Eigen::Ref<Eigen::MatrixXcd> W,
                                 const Eigen::Ref<Eigen::ArrayXd> rho,
                                 const Eigen::Ref<Eigen::ArrayXd> d)
{
  // FIXME:  
  // rho[i] / rho[j] or rho[j] / rho[i] in the next for-loop.  [See {Prop. 1}]
  // The current order is the same as for the second-to-next for-loop.  [See {Prop. 2}]
  // If the first loop changes, the renormalization in the second loop
  // must change...
  
  // Population density balanced parameters
  for (Index i=0; i<necs2020; i++) {
    for (Index j=0; j<i; j++) {
      if (i not_eq j) {
        W(i, j) = W(j, i) * (rho[i] / rho[j]);  // {Prop. 1}  --- division-order should change?
      }
    }
  }
  
  // Balance by the reduced dipole
  for (Index i=0; i<necs2020; i++) {
    auto norm = - std::imag(d[i] * W(i, i)) / std::imag(W.row(i) * d.abs().matrix() - std::abs(d[i]) * W(i, i));  // {Prop. 2}  --- dot-product axis should to column?
    for (Index j=0; j<necs2020; j++) {
      if (i not_eq j) {
        W(i, j) *= norm;  // {Prop. 2}  --- If yes, then change W(i, j) to W(j, i)
      }
    }
  }
}


void Absorption::PredefinedModel::makarov2020_o2_lines_ecs(ComplexVector& I, const Vector& f, Numeric P, Numeric T, Numeric water_vmr)
{
  using Constant::h;
  using Constant::k;
  using Constant::pi;
  using Constant::pow3;
  
  constexpr auto qids = init_mpm2020_qids(0, 0);
  constexpr auto lsm = init_mpm2020_lsm();
  constexpr auto T0 = 300;
  
  // Central frequency [Hz]
  constexpr std::array<Numeric, nlines_mpm2020> f0 = {
    1.18750334E+11, 5.6264774E+10, 6.2486253E+10, 5.8446588E+10, 6.0306056E+10, 
    5.9590983E+10, 5.9164204E+10, 6.0434778E+10, 5.8323877E+10, 6.1150562E+10, 
    5.7612486E+10, 6.1800158E+10, 5.6968211E+10, 6.2411220E+10, 5.6363399E+10, 
    6.2997984E+10, 5.5783815E+10, 6.3568526E+10, 5.5221384E+10, 6.4127775E+10, 
    5.4671180E+10, 6.4678910E+10, 5.4130025E+10, 6.5224078E+10, 5.3595775E+10, 
    6.5764779E+10, 5.3066934E+10, 6.6302096E+10, 5.2542418E+10, 6.6836834E+10, 
    5.2021429E+10, 6.7369601E+10, 5.1503360E+10, 6.7900868E+10, 5.0987745E+10, 
    6.8431006E+10, 5.0474214E+10, 6.8960312E+10, 3.68498246E+11, 4.24763020E+11, 
    4.87249273E+11, 7.15392902E+11, 7.73839490E+11, 8.34145546E+11, };
    
  // Intensity [1 / Pa] (rounded to 10 digits because at most 9 digits exist in f0)
  constexpr std::array<Numeric, nlines_mpm2020> intens = {
    1.591521878E-21, 1.941172240E-21, 4.834543970E-21, 4.959264029E-21, 7.010386457E-21, 
    7.051673348E-21, 8.085012578E-21, 8.108262250E-21, 8.145673278E-21, 8.149757320E-21, 
    7.396406085E-21, 7.401923754E-21, 6.162286575E-21, 6.168475265E-21, 4.749226167E-21, 
    4.754435107E-21, 3.405982896E-21, 3.408455562E-21, 2.282498656E-21, 2.283934341E-21, 
    1.432692459E-21, 1.433513473E-21, 8.439995690E-22, 8.443521837E-22, 4.672706507E-22, 
    4.676049313E-22, 2.435008301E-22, 2.437304596E-22, 1.195038747E-22, 1.196873412E-22, 
    5.532759045E-23, 5.537261239E-23, 2.416832398E-23, 2.418989865E-23, 9.969285671E-24, 
    9.977543709E-24, 3.882541154E-24, 3.888101811E-24, 3.676253816E-23, 3.017524005E-22, 
    9.792882227E-23, 2.756166168E-23, 1.486462215E-22, 4.411918954E-23, };
    
  // Temperature intensity modifier
  constexpr std::array<Numeric, nlines_mpm2020> a2 = {
    0.01, 0.014, 0.083, 0.083, 0.207, 0.207, 0.387, 0.386,
    0.621, 0.621, 0.910, 0.910, 1.255, 1.255, 1.654, 1.654,
    2.109, 2.108, 2.618, 2.617, 3.182, 3.181, 3.800, 3.800,
    4.474, 4.473, 5.201, 5.200, 5.983, 5.982, 6.819, 6.818,
    7.709, 7.708, 8.653, 8.652, 9.651, 9.650, 0.048, 0.044,
    0.049, 0.145, 0.141, 0.145};
  
  // Computed arrays values
  Eigen::Array<Numeric, necs2020, 1> d;
  Eigen::Array<Numeric, necs2020, 1> rho;
  Eigen::Matrix<Complex, necs2020, necs2020> W;
  
  // Diagonal and upper triangular values
  for (Index i=0; i<necs2020; i++) {
    auto Ji = qids[i].UpperQuantumNumber(QuantumNumberType::J);
    auto Jf = qids[i].LowerQuantumNumber(QuantumNumberType::J);
    auto Ni = qids[i].UpperQuantumNumber(QuantumNumberType::N);
    auto Nf = qids[i].LowerQuantumNumber(QuantumNumberType::N);
    for (Index j=i; j<necs2020; j++) {
      auto Ji_p = qids[j].UpperQuantumNumber(QuantumNumberType::J);
      auto Jf_p = qids[j].LowerQuantumNumber(QuantumNumberType::J);
      auto Ni_p = qids[j].UpperQuantumNumber(QuantumNumberType::N);
      auto Nf_p = qids[j].LowerQuantumNumber(QuantumNumberType::N);
      if (i not_eq j) {
        W(i, j) = 1i * o2_ecs_wigner_symbol_tran(Ji, Jf, Ni, Nf, 1, 1, Ji_p, Jf_p, Ni_p, Nf_p, 1, T);
      } else {
        W(i, j) = 1i * (1 + 0.1*water_vmr) * P * lsm[i].compute(T, T0, LineShape::Variable::G0);
      }
    }
  }
  
  std::transform(qids.cbegin(), qids.cbegin()+necs2020, d.data(), [](auto& qns){return 
    o2_makarov2013_reduced_dipole (qns.UpperQuantumNumber(QuantumNumberType::J),
                                   qns.LowerQuantumNumber(QuantumNumberType::J),
                                   qns.UpperQuantumNumber(QuantumNumberType::N));});
  
  std::transform(qids.cbegin(), qids.cbegin()+necs2020, rho.data(), [T](auto& qns){return 
    boltzman_factor(T, o2_ecs_erot_jn_same(qns.LowerQuantumNumber(QuantumNumberType::J)));});
  
  normalize_relaxation_matrix(W, rho, d);
  
  // Rescale dipole to be proportional to line strength
  const Numeric theta = T0 / T;
  const Numeric theta_m1 = theta - 1;
  const Numeric theta_3 = pow3(theta);
  for (Index i=0; i<necs2020; i++) {
    const Numeric ST = pi * theta_3 * P * intens[i] * std::exp(-a2[i] * theta_m1);
    const Numeric sgn = std::abs(d[i]) / d[i];
    d[i] = sgn * f0[i] * std::sqrt(ST / rho[i]);
  }
  
  /*if (not full)*/
  for (Index i=0; i<necs2020; i++) {
    W(i, i) += f0[i];
  }
  
  const Eigen::ComplexEigenSolver<Eigen::Matrix<Complex, necs2020, necs2020>> eV(W, true);
  auto& D = eV.eigenvalues(); 
  auto& V = eV.eigenvectors(); 
  auto& Vinv = W = eV.eigenvectors().inverse();  // Reuse W memory but with different &name
  
  Eigen::Array<Complex, necs2020, 1> B; B *= 0;
  for (Index m=0; m<necs2020; m++) {
    for (Index i=0; i<necs2020; i++) {
      for (Index j=0; j<necs2020; j++) {
        B[m] += rho[i] * d[i] * d[j] * V(j, m) * Vinv(m, i);
      }
    }
  }

//  // Lorentz profile!
//  const Numeric div = Constant::inv_pi / rho.cwiseProduct(d.cwiseAbs2()).sum();
//  for (Index i=0; i<f.nelem(); i++) {
//    I[i] = div * (B.array() / (f[i] - D.array())).sum();
  
  auto GD0 = Linefunctions::DopplerConstant(T, SpeciesTag("O2-66").SpeciesMass());
  for (Index i=0; i<necs2020; i++) {
    for (Index iv=0; iv<f.nelem(); iv++) {
      I[iv] +=  (Constant::inv_sqrt_pi / (GD0 * D[i].real())) * Faddeeva::w((f[iv] - std::conj(D[i])) / (GD0 * D[i].real())) * B[i];
    }
  }
}
