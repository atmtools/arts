#ifndef LINEMIXING_HITRAN_H
#define LINEMIXING_HITRAN_H

#include "complex.h"
#include "constants.h"
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

};  // lm_hitran_2017

#endif  // LINEMIXING_HITRAN_H
