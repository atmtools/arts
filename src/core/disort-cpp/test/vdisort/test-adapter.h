#pragma once

#include <debug.h>
#include <vdisort.h>

namespace vdisort_test {
inline rtepack::stokvec as_stokvec(const ConstVectorView& value) {
  ARTS_USER_ERROR_IF(
      value.size() != vdisort::stokes_dimension, "A Stokes vector must have four components, got {}", value.size());
  return {value[0], value[1], value[2], value[3]};
}

inline rtepack::muelmat as_muelmat(const ConstMatrixView& value) {
  ARTS_USER_ERROR_IF(value.nrows() != vdisort::stokes_dimension or value.ncols() != vdisort::stokes_dimension,
                     "A Mueller matrix must have shape [4, 4], got {:B,}",
                     value.shape());
  rtepack::muelmat out{0.0};
  for (Index i = 0; i < vdisort::stokes_dimension; ++i)
    for (Index j = 0; j < vdisort::stokes_dimension; ++j) out[i, j] = value[i, j];
  return out;
}

inline vdisort::phase_matrix_data as_phase_data(const Tensor7& value) {
  const auto [nalpha, nfourier, nlayers, nout, nin, nstokes_out, nstokes_in] = value.shape();
  ARTS_USER_ERROR_IF(nstokes_out != vdisort::stokes_dimension or nstokes_in != vdisort::stokes_dimension,
                     "The last two phase-matrix dimensions must both be 4, got {:B,}",
                     value.shape());
  vdisort::phase_matrix_data out(nalpha, nfourier, nlayers, nout, nin);
  for (Index alpha = 0; alpha < nalpha; ++alpha)
    for (Index m = 0; m < nfourier; ++m)
      for (Index layer = 0; layer < nlayers; ++layer)
        for (Index i = 0; i < nout; ++i)
          for (Index j = 0; j < nin; ++j) out[alpha, m, layer, i, j] = as_muelmat(value[alpha, m, layer, i, j]);
  return out;
}

inline vdisort::beam_phase_matrix_data as_beam_phase_data(const Tensor6& value) {
  if (value.size() == 0) return {};
  const auto [nalpha, nfourier, nlayers, nout, nstokes_out, nstokes_in] = value.shape();
  ARTS_USER_ERROR_IF(nstokes_out != vdisort::stokes_dimension or nstokes_in != vdisort::stokes_dimension,
                     "The last two beam phase-matrix dimensions must both be 4, got {:B,}",
                     value.shape());
  vdisort::beam_phase_matrix_data out(nalpha, nfourier, nlayers, nout);
  for (Index alpha = 0; alpha < nalpha; ++alpha)
    for (Index m = 0; m < nfourier; ++m)
      for (Index layer = 0; layer < nlayers; ++layer)
        for (Index i = 0; i < nout; ++i) out[alpha, m, layer, i] = as_muelmat(value[alpha, m, layer, i]);
  return out;
}

inline rtepack::stokvec_tensor3 as_boundary_data(const Tensor4& value) {
  const auto [nalpha, nfourier, nstream, nstokes] = value.shape();
  ARTS_USER_ERROR_IF(
      nstokes != vdisort::stokes_dimension, "The boundary Stokes dimension must be 4, got {:B,}", value.shape());
  rtepack::stokvec_tensor3 out(nalpha, nfourier, nstream);
  for (Index alpha = 0; alpha < nalpha; ++alpha)
    for (Index m = 0; m < nfourier; ++m)
      for (Index i = 0; i < nstream; ++i) out[alpha, m, i] = as_stokvec(value[alpha, m, i]);
  return out;
}

inline rtepack::stokvec_matrix as_source_data(const Tensor3& value) {
  const auto [nlayers, ncoeffs, nstokes] = value.shape();
  ARTS_USER_ERROR_IF(
      nstokes != vdisort::stokes_dimension, "The source Stokes dimension must be 4, got {:B,}", value.shape());
  rtepack::stokvec_matrix out(nlayers, ncoeffs);
  for (Index layer = 0; layer < nlayers; ++layer)
    for (Index coefficient = 0; coefficient < ncoeffs; ++coefficient)
      out[layer, coefficient] = as_stokvec(value[layer, coefficient]);
  return out;
}

inline vdisort::main_data make_solver(Index                      nquad,
                                      Index                      nfourier,
                                      AscendingGrid              tau,
                                      Vector                     omega,
                                      Tensor7                    phase,
                                      Tensor4                    boundary_up,
                                      Tensor4                    boundary_down,
                                      Tensor3                    source,
                                      std::vector<vdisort::BDRF> brdf,
                                      Numeric                    mu0,
                                      Vector                     beam,
                                      Numeric                    phi0,
                                      Tensor6                    beam_phase    = {},
                                      Vector                     source_scale  = {},
                                      Vector                     source_offset = {}) try {
  return vdisort::main_data(nquad,
                            nfourier,
                            std::move(tau),
                            std::move(omega),
                            as_phase_data(phase),
                            as_boundary_data(boundary_up),
                            as_boundary_data(boundary_down),
                            as_source_data(source),
                            std::move(brdf),
                            mu0,
                            as_stokvec(beam),
                            phi0,
                            as_beam_phase_data(beam_phase),
                            std::move(source_scale),
                            std::move(source_offset));
} catch (std::exception& e) { throw std::runtime_error(std::format("Error in make_solver:\n{}", e.what())); }
}  // namespace vdisort_test
