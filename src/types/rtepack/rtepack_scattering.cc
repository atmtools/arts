#include "rtepack_scattering.h"
#include "rtepack_mueller_matrix.h"
#include "rtepack_multitype.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
Array<muelmat_vector> bulk_backscatter(const ConstTensor5View &Pe,
                                       const ConstMatrixView &pnd) {
  ARTS_ASSERT(Pe.ncols() == 4 and Pe.nrows() == 4)
  const Index nv = Pe.npages();
  const Index np = Pe.nbooks();
  const Index ne = Pe.nshelves();

  Array<muelmat_vector> aotm(np, muelmat_vector(nv));
  if (ne == 0)
    return aotm;

  for (Index ip = 0; ip < np; ip++) {
    for (Index iv = 0; iv < nv; iv++) {
      aotm[ip][iv] = pnd(0, ip) * muelmat(Pe(0, ip, iv, joker, joker));
    }
  }

  for (Index ie = 1; ie < ne; ie++) {
    for (Index ip = 0; ip < np; ip++) {
      for (Index iv = 0; iv < nv; iv++) {
        aotm[ip][iv] += pnd(ie, ip) * muelmat(Pe(ie, ip, iv, joker, joker));
      }
    }
  }
  return aotm;
}

Array<muelmat_matrix>
bulk_backscatter_derivative(const ConstTensor5View &Pe,
                            const ArrayOfMatrix &dpnd_dx) {
  ARTS_ASSERT(Pe.ncols() == 4 and Pe.nrows() == 4)

  const Index nv = Pe.npages();
  const Index np = Pe.nbooks();
  const Index ne = Pe.nshelves();
  const Size nq = dpnd_dx.size();

  Array<muelmat_matrix> aoaotm(np, muelmat_matrix(nq, nv));

  for (Index ip = 0; ip < np; ip++) {
    for (Size iq = 0; iq < nq; iq++) {
      aoaotm[ip][iq] = 0.0;
      for (Index iv = 0; iv < nv; iv++) {
        for (Index ie = 0; ie < ne; ie++) {
          aoaotm[ip](iq, iv) +=
              dpnd_dx[iq](ie, ip) * muelmat(Pe(ie, ip, iv, joker, joker));
        }
      }
    }
  }
  return aoaotm;
}

void setBackscatterTransmission(stokvec_vector& out,
                                                 const stokvec_vector& I0,
                                                 const muelmat_vector& Tr,
                                                 const muelmat_vector& Tf,
                                                 const muelmat_vector& Z) {
  for (Index i = 0; i < out.nelem(); i++)
    out[i] = Tr[i] * Z[i] * Tf[i] * I0[i];
}

void setBackscatterTransmissionDerivative(
  stokvec_matrix& out,
    const stokvec_vector& I0,
    const muelmat_vector& Tr,
    const muelmat_vector& Tf,
    const muelmat_matrix& dZ) {
  for (Index j = 0; j < dZ.nrows(); j++)
  for (Index i = 0; i < dZ.ncols(); i++)
    out(j, i) += Tr[i] * dZ(j, i) * Tf[i] * I0[i];
}

void bulk_backscatter_commutative_transmission_rte(
    Array<stokvec_vector> &I,
    Array<Array<stokvec_matrix>> &dI,
    const stokvec_vector &I_incoming,
    const Array<muelmat_vector> &T,
    const Array<muelmat_vector> &PiTf,
    const Array<muelmat_vector> &PiTr,
    const Array<muelmat_vector> &Z,
    const Array<muelmat_matrix> &dT1,
    const Array<muelmat_matrix> &dT2,
    const Array<muelmat_matrix> &dZ) {
  const Size np = dT1.size();
  const Index nv = np ? dT1.front().ncols() : 0;
  const Index nq = np ? dT1.front().nrows() : 0;

  // For all transmission, the I-vector is the same
  for (Size ip = 0; ip < np; ip++)
    setBackscatterTransmission(I[ip], I_incoming, PiTr[ip], PiTf[ip], Z[ip]);

  for (Size ip = 0; ip < np; ip++) {
    setBackscatterTransmissionDerivative(dI[ip][ip], I_incoming, PiTr[ip],
                                         PiTf[ip], dZ[ip]);
  }

  for (Size ip = 0; ip < np; ip++) {
    for (Size j = ip; j < np; j++) {
      for (Index iq = 0; iq < nq; iq++) {
        for (Index iv = 0; iv < nv; iv++) {
          dI[ip][j](iq, iv) +=
              T[ip][iv] *(
              (dT1[ip](iq, iv) + dT2[ip](iq, iv)) * I[j][iv]);

          if (j < np - 1 and j > ip)
            dI[ip][j](iq, iv) +=
                inv(T[ip + 1][iv]) * (
                (dT1[ip](iq, iv) + dT2[ip](iq, iv)) * I[j][iv]);
        }
      }
    }
  }
}
} // namespace rtepack