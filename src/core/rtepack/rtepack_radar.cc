#include "rtepack_radar.h"

#include <debug.h>

#include <utility>

#include "rtepack_multitype.h"
#include "rtepack_transmission.h"

namespace rtepack {
namespace {
constexpr bool uses_magnus(const TransmittanceOption option) {
  switch (option) {
    case TransmittanceOption::constant:
    case TransmittanceOption::linsrc:
    case TransmittanceOption::linprop:      return false;
    case TransmittanceOption::magop:
    case TransmittanceOption::magop_linsrc: return true;
  }
  std::unreachable();
}

//! Transmission over a layer from start to end, and its derivatives with regards to either end
struct layer_transmission {
  propmat start;
  propmat end;
  Numeric r;
  bool    magnus;
  tran    state;
  muelmat T;

  layer_transmission(const propmat &k_start, const propmat &k_end, const Numeric r_, const TransmittanceOption option)
      : start(k_start),
        end(k_end),
        r(r_),
        magnus(uses_magnus(option)),
        state(magnus ? tran{k_start, k_end, r_, MagnusOperator::magnus} : tran{k_start, k_end, r_}),
        T(state()) {}

  [[nodiscard]] muelmat d_start(const propmat &dk) const {
    return magnus ? state.magnus_deriv(T, start, end, dk, r, 0.0, true) : state.deriv(T, -0.5 * r * dk);
  }

  [[nodiscard]] muelmat d_end(const propmat &dk) const {
    return magnus ? state.magnus_deriv(T, start, end, dk, r, 0.0, false) : state.deriv(T, -0.5 * r * dk);
  }
};
}  // namespace

void radar_single_scattering(stokvec_matrix_view                    I,
                             stokvec_tensor4_view                   dI,
                             const stokvec                         &transmitted,
                             const std::span<const propmat_vector> &K_out,
                             const std::span<const propmat_matrix> &dK_out,
                             const std::span<const propmat_vector> &K_ret,
                             const std::span<const propmat_matrix> &dK_ret,
                             const std::span<const muelmat_vector> &Z,
                             const std::span<const muelmat_matrix> &dZ,
                             const ConstVectorView                 &r,
                             const std::span<const radar_sample>   &samples,
                             const TransmittanceOption              option) {
  const Size np    = K_out.size();
  const Size ngate = I.nrows();
  const Size nf    = I.ncols();
  const Size nq    = dI.nrows();

  if (static_cast<Size>(dI.nbooks()) != ngate or static_cast<Size>(dI.npages()) != np or
      static_cast<Size>(dI.ncols()) != nf)
    arts::user_error("Radar return shapes do not agree: I: {:B,}, dI: {:B,}, np: {}", I.shape(), dI.shape(), np);
  if (K_ret.size() != np or Z.size() != np or dK_out.size() != np or dK_ret.size() != np or dZ.size() != np or
      r.size() + 1 != std::max<Size>(np, 1))
    arts::user_error(
        "Radar level sizes do not agree: K_out: {}, K_ret: {}, Z: {}, dK_out: {}, dK_ret: {}, dZ: {}, r: {}",
        np,
        K_ret.size(),
        Z.size(),
        dK_out.size(),
        dK_ret.size(),
        dZ.size(),
        r.size());
  for (const auto &sample : samples)
    if (sample.gate >= ngate or sample.layer == 0 or sample.layer >= np)
      arts::user_error("Radar sample in layer {} for gate {} is outside {} levels and {} gates",
                       sample.layer,
                       sample.gate,
                       np,
                       ngate);

  I  = stokvec{};
  dI = stokvec{};
  if (samples.empty()) return;

  // Tf * s and Tr up to the near level of the current layer, and their
  // derivatives with regards to every level.
  stokvec_vector  Tfs(nf, transmitted);
  muelmat_vector  Tr(nf, muelmat::id());
  stokvec_tensor3 dTfs(np, nq, nf, stokvec{});
  muelmat_tensor3 dTr(np, nq, nf, muelmat{});

  Size is = 0;
  for (Size i = 1; i < np; i++) {
    for (; is < samples.size() and samples[is].layer == i; is++) {
      const auto &[layer, x, w, gate] = samples[is];
      const Numeric rx                = x * r[i - 1];

      for (Size iv = 0; iv < nf; iv++) {
        const propmat ko = (1.0 - x) * K_out[i - 1][iv] + x * K_out[i][iv];
        const propmat kr = (1.0 - x) * K_ret[i - 1][iv] + x * K_ret[i][iv];
        const muelmat Zx = (1.0 - x) * Z[i - 1][iv] + x * Z[i][iv];

        const layer_transmission out{K_out[i - 1][iv], ko, rx, option};
        const layer_transmission ret{kr, K_ret[i - 1][iv], rx, option};

        const stokvec tfs  = out.T * Tfs[iv];
        const stokvec zs   = Zx * tfs;
        const muelmat Trx  = Tr[iv] * ret.T;
        I[gate, iv]       += w * (Trx * zs);

        if (nq == 0) continue;

        // Through the cumulative transmissions of all previous layers
        const muelmat M = Trx * Zx * out.T;
        const stokvec u = ret.T * zs;
        for (Size j = 0; j < i; j++)
          for (Size iq = 0; iq < nq; iq++) dI[gate, j, iq, iv] += w * (M * dTfs[j, iq, iv] + dTr[j, iq, iv] * u);

        // Through the partial layer, whose inner end is interpolated
        for (Size iq = 0; iq < nq; iq++) {
          const propmat &dko0 = dK_out[i - 1][iq, iv];
          const propmat &dko1 = dK_out[i][iq, iv];
          const propmat &dkr0 = dK_ret[i - 1][iq, iv];
          const propmat &dkr1 = dK_ret[i][iq, iv];

          const muelmat dout0 = out.d_start(dko0) + out.d_end((1.0 - x) * dko0);
          const muelmat dout1 = out.d_end(x * dko1);
          const muelmat dret0 = ret.d_start((1.0 - x) * dkr0) + ret.d_end(dkr0);
          const muelmat dret1 = ret.d_start(x * dkr1);

          dI[gate, i - 1, iq, iv] += w * (Tr[iv] * (dret0 * zs) + Trx * ((1.0 - x) * dZ[i - 1][iq, iv] * tfs) +
                                          Trx * (Zx * (dout0 * Tfs[iv])));
          dI[gate, i, iq, iv] +=
              w * (Tr[iv] * (dret1 * zs) + Trx * (x * dZ[i][iq, iv] * tfs) + Trx * (Zx * (dout1 * Tfs[iv])));
        }
      }
    }

    if (is == samples.size()) break;

    for (Size iv = 0; iv < nf; iv++) {
      const layer_transmission out{K_out[i - 1][iv], K_out[i][iv], r[i - 1], option};
      const layer_transmission ret{K_ret[i][iv], K_ret[i - 1][iv], r[i - 1], option};

      for (Size j = 0; j < i; j++) {
        for (Size iq = 0; iq < nq; iq++) {
          dTfs[j, iq, iv] = out.T * dTfs[j, iq, iv];
          dTr[j, iq, iv]  = dTr[j, iq, iv] * ret.T;
        }
      }

      for (Size iq = 0; iq < nq; iq++) {
        dTfs[i - 1, iq, iv] += out.d_start(dK_out[i - 1][iq, iv]) * Tfs[iv];
        dTfs[i, iq, iv]     += out.d_end(dK_out[i][iq, iv]) * Tfs[iv];
        dTr[i - 1, iq, iv]  += Tr[iv] * ret.d_end(dK_ret[i - 1][iq, iv]);
        dTr[i, iq, iv]      += Tr[iv] * ret.d_start(dK_ret[i][iq, iv]);
      }

      Tfs[iv] = out.T * Tfs[iv];
      Tr[iv]  = Tr[iv] * ret.T;
    }
  }

  if (is != samples.size()) arts::user_error("Radar samples must be sorted by layer");
}
}  // namespace rtepack
