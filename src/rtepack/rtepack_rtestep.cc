#include "rtepack_rtestep.h"
#include "rtepack_mueller_matrix.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
void two_level_linear_emission_step(
    stokvec_vector_view I, stokvec_matrix_view dI1, stokvec_matrix_view dI2,
    const stokvec_vector_const_view &J1, const stokvec_vector_const_view &J2,
    const stokvec_matrix_const_view &dJ1, const stokvec_matrix_const_view &dJ2,
    const muelmat_vector_const_view &T, const muelmat_vector_const_view &PiT,
    const muelmat_matrix_const_view &dT1,
    const muelmat_matrix_const_view &dT2) {
  const Index N = I.nelem();
  const Index M = dI1.nrows();
  ARTS_ASSERT(N == dI1.ncols() and N == dI2.ncols() and N == J1.nelem() and
              N == J2.nelem() and N == dJ1.ncols() and N == dJ2.ncols() and
              N == T.nelem() and N == PiT.nelem() and N == dT1.ncols() and
              N == dT2.ncols())
  ARTS_ASSERT(M == dI2.nrows() and M == dJ1.nrows() and M == dJ2.nrows() and
              M == dT1.nrows() and M == dT2.nrows())

  for (Index i = 0; i < N; i++) {
    const auto src = avg(J1[i], J2[i]);
    I[i] = I[i] - src; //! ImJ

    for (Index j = 0; j < M; j++) {
      dI1(j, i) += PiT[i] * dT1(j, i) * I[i] + dJ1(j, i) - T[i] * dJ1(j, i);
      dI2(j, i) += PiT[i] * dT2(j, i) * I[i] + dJ2(j, i) - T[i] * dJ2(j, i);
    }

    I[i] = T[i] * I[i] + src;
  }
}

void two_level_linear_emission_step(stokvec_vector_view I,
                                    const stokvec_vector_const_view &J1,
                                    const stokvec_vector_const_view &J2,
                                    const muelmat_vector_const_view &T) {
  const Index N = I.nelem();
  ARTS_ASSERT(N == J1.nelem() and N == J2.nelem() and N == T.nelem())

  for (Index i = 0; i < N; i++) {
    const auto src = avg(J1[i], J2[i]);
    I[i] = T[i] * (I[i] - src) + src;
  }
}

void two_level_linear_transmission_step(stokvec_vector_view I,
                                        stokvec_matrix_view dI1,
                                        stokvec_matrix_view dI2,
                                        const muelmat_vector_const_view &T,
                                        const muelmat_vector_const_view &PiT,
                                        const muelmat_matrix_const_view &dT1,
                                        const muelmat_matrix_const_view &dT2) {
  const Index N = I.nelem();
  const Index M = dI1.nrows();
  ARTS_ASSERT(N == dI1.ncols() and N == dI2.ncols() and N == T.nelem() and
              N == PiT.nelem() and N == dT1.ncols() and N == dT2.ncols())
  ARTS_ASSERT(M == dI2.nrows() and M == dT1.nrows() and M == dT2.nrows())

  for (Index i = 0; i < N; i++) {
    for (Index j = 0; j < M; j++) {
      dI1(j, i) += PiT[i] * dT1(j, i) * I[i];
      dI2(j, i) += PiT[i] * dT2(j, i) * I[i];
    }
    I[i] = T[i] * I[i];
  }
}
} // namespace rtepack