#ifndef LINEMIXING_HITRAN_H
#define LINEMIXING_HITRAN_H

#include "complex.h"
#include "constants.h"
#include "matpackIV.h"
#include "mystring.h"

namespace lm_hitran_2017 {
enum class calctype {
  VP,
  VPALT,
  VPPURE,
  SDVP,
  Rosenkranz,
  Full
};

Vector compute(const Numeric p,
               const Numeric t,
               const Numeric xco2,
               const Numeric xh2o,
               const Numeric sigmin,
               const Numeric sigmax,
               const Numeric stotmax,
               const Numeric dsig,
               const calctype type=calctype::Full);

};  // lm_hitran_2017

#endif  // LINEMIXING_HITRAN_H
