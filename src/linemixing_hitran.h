/* Copyright (C) 2020
 * Richard Larsson <larsson@mps.mpg.de>
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

/**
 * @file linemixing_hitran.h
 * @author Richard Larsson
 * @date 2020-06-23
 * 
 * @brief Namespace and functions to deal with HITRAN linemixing
 */

#ifndef LINEMIXING_HITRAN_H
#define LINEMIXING_HITRAN_H

#include "absorptionlines.h"
#include "complex.h"
#include "constants.h"
#include "linescaling.h"
#include "matpackIV.h"
#include "mystring.h"


struct HitranRelaxationMatrixData {
  Tensor4 W0pp, B0pp;
  Tensor4 W0rp, B0rp;
  Tensor4 W0qp, B0qp;
  
  Tensor4 W0pr, B0pr;
  Tensor4 W0rr, B0rr;
  Tensor4 W0qr, B0qr;
  
  Tensor4 W0pq, B0pq;
  Tensor4 W0rq, B0rq;
  Tensor4 W0qq, B0qq;
  HitranRelaxationMatrixData() {}
  friend std::ostream& operator<<(std::ostream& os, const HitranRelaxationMatrixData& hit) {
    return os << hit.W0pp << '\n' << hit.B0pp << '\n'
              << hit.W0rp << '\n' << hit.B0rp << '\n'
              << hit.W0qp << '\n' << hit.B0qp << '\n'
              << hit.W0pr << '\n' << hit.B0pr << '\n'
              << hit.W0rr << '\n' << hit.B0rr << '\n'
              << hit.W0qr << '\n' << hit.B0qr << '\n'
              << hit.W0pq << '\n' << hit.B0pq << '\n'
              << hit.W0rq << '\n' << hit.B0rq << '\n'
              << hit.W0qq << '\n' << hit.B0qq << '\n';
  }
};


namespace lm_hitran_2017 {
enum class calctype {
  FullVP,
  FullRosenkranz,
  FullW,
  SDVP,
  SDRosenkranz,
  SDW,
  NoneVP,
  NoneRosenkranz,
  NoneW
};

Vector compute(const Numeric p,
               const Numeric t,
               const Numeric xco2,
               const Numeric xh2o,
               const ConstVectorView invcm_grid,
               const Numeric stotmax,
               const calctype type=calctype::FullW);

/** Compute the absorptionlines
 * 
 * @param[in] hitran Hitran data for the relaxation matrix calculations
 * @param[in] bands List of absorption bands
 * @param[in] P Pressure in Pascal
 * @param[in] T Temperatures in Kelvin
 * @param[in] vmrs VMR ratio.  Must be 3-long containing {air, h2o, co2} vmrs
 * @param[in] f_grid Frequency grid in Hz
 * @param[in] partition_functions As WSV
 * @return The absorption, a vector same length as f_grid
 */
Vector compute(const HitranRelaxationMatrixData& hitran,
               const ArrayOfAbsorptionLines& bands,
               const Numeric P,
               const Numeric T,
               const ConstVectorView vmrs,
               const ConstVectorView f_grid);

/** Class that controls ReadFromLineMixingStream output */
enum class ModeOfLineMixing {
  VP,  // Sets LineShape::VP, will not use LineMixing code; Sets ByLTE mode
  VP_Y,  // Sets LineShape::VP, will use LineMixing code with pressure > linemixinglimit;  Sets ByRosenkranzRelmatLTE mode
  SDVP,  // Sets LineShape::SDVP, will not use LineMixing code; Sets ByLTE mode
  SDVP_Y,  // Sets LineShape::SDVP, will use LineMixing code with pressure > linemixinglimit;  Sets ByHITRANRosenkranzRelmat mode
  FullW, // Sets LineShape::Lorentz, will use LineMixing code with pressure > linemixinglimit;  Sets ByHITRANFullRelmat mode
  VP_W  // Sets LineShape::Voigt, will use LineMixing code with pressure > linemixinglimit;  Sets ByHITRANFullRelmat mode
};

constexpr bool typeVP(ModeOfLineMixing x)
{
  return x == ModeOfLineMixing::VP or x == ModeOfLineMixing::VP_Y or x == ModeOfLineMixing::FullW or x == ModeOfLineMixing::VP_W;
}

constexpr bool typeLP(ModeOfLineMixing x)
{
  return x == ModeOfLineMixing::FullW;
}

constexpr bool typeFull(ModeOfLineMixing x)
{
  return x == ModeOfLineMixing::FullW or x == ModeOfLineMixing::VP_W;
}

/** Read from HITRAN online line mixing file
 * 
 * @param[out] hitran Hitran data for the relaxation matrix calculations
 * @param[out] bands List of absorption bands
 * @param[in] basedir The base directory of the HITRAN line mixing files
 * @param[in] linemixinglimit The pressure limit for using line mixing instead of pure Voigt
 * @param[in] fmin Minimum frequency
 * @param[in] fmax Maximum frequency
 * @param[in] stot Minimum band-strength
 * @param[in] mode The type of calculations
 */ 
void read(HitranRelaxationMatrixData& hitran, ArrayOfAbsorptionLines& bands, const String& basedir, const Numeric linemixinglimit, const Numeric fmin, const Numeric fmax, const Numeric stot, const ModeOfLineMixing mode);

void hitran_lm_eigenvalue_adaptation(AbsorptionLines& band,
                                     const Vector& temperatures,
                                     const HitranRelaxationMatrixData& hitran,
                                     const Numeric P0,
                                     const Index ord);
};  // lm_hitran_2017

#endif  // LINEMIXING_HITRAN_H
