#ifndef LINEMIXING_HITRAN_H
#define LINEMIXING_HITRAN_H

#include "absorptionlines.h"
#include "complex.h"
#include "constants.h"
#include "linescaling.h"
#include "matpackIV.h"
#include "mystring.h"

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

struct HitranRelaxationMatrixData {
  ArrayOfAbsorptionLines bands;
  
  Tensor4 W0pp, B0pp;
  Tensor4 W0rp, B0rp;
  Tensor4 W0qp, B0qp;
  
  Tensor4 W0pr, B0pr;
  Tensor4 W0rr, B0rr;
  Tensor4 W0qr, B0qr;
  
  Tensor4 W0pq, B0pq;
  Tensor4 W0rq, B0rq;
  Tensor4 W0qq, B0qq;
  
  ArrayOfVector dip0;
  ArrayOfVector pop0;
};

Vector compute(const HitranRelaxationMatrixData& hitran,
               const Numeric P,
               const Numeric T,
               const ConstVectorView vmrs,
               const ConstVectorView f_grid,
               const SpeciesAuxData& partition_functions);

/** Class that controls ReadFromLineMixingStream output */
enum class ModeOfLineMixing {
  VP,  // Sets LineShape::VP, will not use LineMixing code; Sets ByLTE mode
  VP_Y,  // Sets LineShape::VP, will use LineMixing code with pressure > linemixinglimit;  Sets ByRosenkranzRelmatLTE mode
  SDVP,  // Sets LineShape::SDVP, will not use LineMixing code; Sets ByLTE mode
  SDVP_Y,  // Sets LineShape::SDVP, will use LineMixing code with pressure > linemixinglimit;  Sets ByHITRANRosenkranzRelmat mode
  FullW  // Sets LineShape::Lorentz, will use LineMixing code with pressure > linemixinglimit;  Sets ByHITRANFullRelmat mode
};

constexpr bool typeVP(ModeOfLineMixing x)
{
  return x == ModeOfLineMixing::VP or x == ModeOfLineMixing::VP_Y or x == ModeOfLineMixing::FullW;
}

/** Read from HITRAN online line mixing file
 * 
 * @param[in] basedir The base directory of the HITRAN line mixing files
 * @param[in] linemixinglimit The pressure limit for using line mixing instead of pure Voigt
 * @param[in] fmin Minimum frequency
 * @param[in] fmax Maximum frequency
 * @param[in] mode The type of calculations
 * @return HitranRelaxationMatrixData
 */ 
HitranRelaxationMatrixData read(const String& basedir, const Numeric linemixinglimit, const Numeric fmin, const Numeric fmax, const Numeric stot, const ModeOfLineMixing mode);

};  // lm_hitran_2017

#endif  // LINEMIXING_HITRAN_H
