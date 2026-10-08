#include "radintg3.h"

#include <debug.h>

#include <array>

namespace polradtran::rt3 {
void initialize(Numeric          delta_z,
                ConstVectorView  mu_values,
                Numeric          extinction,
                Numeric          albedo,
                ConstTensor5View phase_function,
                Tensor5View      reflect,
                Tensor5View      trans) {
  const Index nummu   = mu_values.size();
  const Index nstokes = reflect.extent(4);
  ARTS_USER_ERROR_IF(nummu < 1 or nstokes < 1 or
                         reflect.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         trans.shape() != reflect.shape() or
                         phase_function.shape() != (std::array<Index, 5>{4, nummu, nstokes, nummu, nstokes}),
                     "INITIALIZE with {} mu_values needs phase_function [4, nummu, nstokes, nummu, nstokes] and "
                     "reflect and trans [2, nummu, nstokes, nummu, nstokes]; got {:B,}, {:B,} and {:B,}",
                     nummu,
                     phase_function.shape(),
                     reflect.shape(),
                     trans.shape());

  // The rows of angle j1, [l, joker, joker, j1, joker]:
  // REFLECT = f ALBEDO PHASE_FUNCTION(2 or 3) and
  // TRANS = DIAG - f (DIAG - ALBEDO PHASE_FUNCTION(1 or 4)).  The diagonal
  // of TRANS is kept in that form, 1 minus a small number rounded once: it
  // carries the thin layer's extinction, which the doubling amplifies.
  for (Index j1 = 0; j1 < nummu; j1++) {
    const Numeric f = delta_z / mu_values[j1] * extinction;
    for (Index l = 0; l < 2; l++) {
      reflect[l, joker, joker, j1, joker]  = phase_function[l + 1, joker, joker, j1, joker];
      reflect[l, joker, joker, j1, joker] *= f * albedo;
      trans[l, joker, joker, j1, joker]    = phase_function[3 * l, joker, joker, j1, joker];
      trans[l, joker, joker, j1, joker]   *= f * albedo;
      for (Index i = 0; i < nstokes; i++)
        trans[l, j1, i, j1, i] = 1.0 - f * (1.0 - albedo * phase_function[3 * l, j1, i, j1, i]);
    }
  }
}

void initial_source(Numeric          delta_z,
                    ConstVectorView  mu_values,
                    Numeric          extinction,
                    ConstTensor3View source_vector,
                    Tensor3View      source) {
  const Index nummu   = mu_values.size();
  const Index nstokes = source.extent(2);
  ARTS_USER_ERROR_IF(nummu < 1 or nstokes < 1 or source.shape() != (std::array<Index, 3>{2, nummu, nstokes}) or
                         source_vector.shape() != source.shape(),
                     "INITIAL_SOURCE with {} mu_values needs source_vector and source [2, nummu, nstokes]; got "
                     "{:B,} and {:B,}",
                     nummu,
                     source_vector.shape(),
                     source.shape());

  for (Index j = 0; j < nummu; j++) {
    const Numeric tmp        = delta_z / mu_values[j];
    source[joker, j, joker]  = source_vector[joker, j, joker];
    source[joker, j, joker] *= tmp * extinction;
  }
}
}  // namespace polradtran::rt3
