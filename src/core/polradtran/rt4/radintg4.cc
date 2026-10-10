#include "radintg4.h"

#include <arts_constants.h>
#include <debug.h>

#include <array>

namespace polradtran::rt4 {
void initialize(Numeric          delta_z,
                ConstVectorView  mu_values,
                ConstVectorView  quad_weights,
                Numeric          gas_extinct,
                ConstTensor4View extinct_matrix,
                ConstTensor5View scatter_matrix,
                Tensor5View      reflect,
                Tensor5View      trans) {
  const Index nummu   = mu_values.extent(0);
  const Index nstokes = reflect.extent(4);
  ARTS_USER_ERROR_IF(quad_weights.extent(0) != nummu or
                         extinct_matrix.shape() != (std::array<Index, 4>{2, nummu, nstokes, nstokes}) or
                         scatter_matrix.shape() != (std::array<Index, 5>{4, nummu, nstokes, nummu, nstokes}) or
                         reflect.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         trans.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}),
                     "INITIALIZE with {} mu_values needs quad_weights [nummu], extinct_matrix [2, nummu, nstokes, "
                     "nstokes], scatter_matrix [4, nummu, nstokes, nummu, nstokes] and reflect and trans [2, nummu, "
                     "nstokes, nummu, nstokes]; got {}, {:B,}, {:B,}, {:B,} and {:B,}",
                     nummu,
                     quad_weights.extent(0),
                     extinct_matrix.shape(),
                     scatter_matrix.shape(),
                     reflect.shape(),
                     trans.shape());

  const Numeric c = Constant::two_pi;

  for (Index i2 = 0; i2 < nstokes; i2++) {
    for (Index j2 = 0; j2 < nummu; j2++) {
      const Numeric tmp = delta_z / mu_values[j2];
      for (Index i1 = 0; i1 < nstokes; i1++) {
        Numeric gext = 0.0;
        if (i1 == i2) gext = gas_extinct;
        for (Index j1 = 0; j1 < nummu; j1++) {
          reflect[0, j1, i1, j2, i2] = c * tmp * quad_weights[j1] * scatter_matrix[1, j1, i1, j2, i2];
          reflect[1, j1, i1, j2, i2] = c * tmp * quad_weights[j1] * scatter_matrix[2, j1, i1, j2, i2];
          Numeric diag               = 0.0;
          if (i1 == i2 and j1 == j2) diag = 1.0;
          Numeric ext = 0.0;
          if (j1 == j2) ext = extinct_matrix[0, j2, i1, i2] + gext;
          trans[0, j1, i1, j2, i2] = diag - tmp * (ext - c * quad_weights[j1] * scatter_matrix[0, j1, i1, j2, i2]);
          if (j1 == j2) ext = extinct_matrix[1, j2, i1, i2] + gext;
          trans[1, j1, i1, j2, i2] = diag - tmp * (ext - c * quad_weights[j1] * scatter_matrix[3, j1, i1, j2, i2]);
        }
      }
    }
  }
}

void initial_source(Numeric          delta_z,
                    ConstVectorView  mu_values,
                    Numeric          planck,
                    ConstTensor3View emis_vector,
                    Numeric          gas_extinct,
                    Tensor3View      source) {
  const Index nummu   = mu_values.extent(0);
  const Index nstokes = source.extent(2);
  ARTS_USER_ERROR_IF(emis_vector.shape() != (std::array<Index, 3>{2, nummu, nstokes}) or
                         source.shape() != (std::array<Index, 3>{2, nummu, nstokes}),
                     "INITIAL_SOURCE with {} mu_values needs emis_vector and source [2, nummu, nstokes]; got "
                     "{:B,} and {:B,}",
                     nummu,
                     emis_vector.shape(),
                     source.shape());

  source = 0.0;

  for (Index i = 0; i < nstokes; i++) {
    for (Index j = 0; j < nummu; j++) {
      Numeric ext = 0.0;
      if (i == 0) ext = gas_extinct;
      const Numeric tmp = planck * delta_z / mu_values[j];
      source[0, j, i]   = tmp * (emis_vector[0, j, i] + ext);
      source[1, j, i]   = tmp * (emis_vector[1, j, i] + ext);
    }
  }
}
}  // namespace polradtran::rt4
