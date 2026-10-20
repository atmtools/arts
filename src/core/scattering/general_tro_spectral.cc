#include "general_tro_spectral.h"

#include <debug.h>

#include <memory>

ScatteringTroSpectralVector ScatteringGeneralSpectralTRO::get_bulk_scattering_properties_tro_spectral(
    const AtmPoint& atm_point, const Vector& f_grid, Index degree) const {
  ARTS_USER_ERROR_IF(degree < 0, "The Legendre degree must be >= 0, got {}", degree)

  // The user's function must give exactly what was asked
  auto        out = f(atm_point, f_grid, degree);
  const Index nf  = static_cast<Index>(f_grid.size());
  ARTS_USER_ERROR_IF(not out.phase_matrix.has_value(),
                     "The user's scattering function gives no phase matrix; it must give its Legendre series")
  ARTS_USER_ERROR_IF(out.phase_matrix->nrows() != nf or out.phase_matrix->ncols() != degree + 1,
                     "The user's scattering function gives a Legendre series of shape [{}, {}], but [{}, {}] "
                     "([frequencies, degree + 1]) was asked for",
                     out.phase_matrix->nrows(),
                     out.phase_matrix->ncols(),
                     nf,
                     degree + 1)
  ARTS_USER_ERROR_IF(
      static_cast<Index>(out.extinction_matrix.size()) != nf or static_cast<Index>(out.absorption_vector.size()) != nf,
      "The user's scattering function gives {} extinction matrices and {} absorption vectors, but "
      "there are {} frequencies",
      out.extinction_matrix.size(),
      out.absorption_vector.size(),
      nf)
  return out;
}

scattering::BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Fourier>
ScatteringGeneralSpectralTRO::get_bulk_scattering_properties_aro_fourier(const AtmPoint& atm_point,
                                                                         const Vector&   f_grid,
                                                                         const Vector&   za_inc_grid,
                                                                         const Vector&   za_scat_grid,
                                                                         Index           max_mode) const {
  ARTS_USER_ERROR_IF(degree < 0,
                     "The azimuthal Fourier modes of a user's scattering function need the degree of its Legendre "
                     "series; set degree (>= 0) on the ScatteringGeneralSpectralTRO, it is {}",
                     degree)

  const auto bulk = get_bulk_scattering_properties_tro_spectral(atm_point, f_grid, degree)
                        .to_general(std::make_shared<Vector>(f_grid));
  const auto za_inc  = std::make_shared<const Vector>(za_inc_grid);
  const auto za_scat = std::make_shared<const scattering::ZenithAngleGrid>(
      scattering::IrregularZenithAngleGrid(za_scat_grid));
  return {.phase_matrix      = bulk.phase_matrix->to_lab_frame_fourier_modes(za_inc, za_scat, max_mode),
          .extinction_matrix = bulk.extinction_matrix.to_gridded().to_lab_frame(za_inc).to_fourier(),
          .absorption_vector = bulk.absorption_vector.to_gridded().to_lab_frame(za_inc).to_fourier()};
}

ScatteringTroSpectralVector& ScatteringTroSpectralVector::operator+=(const ScatteringTroSpectralVector& other) {
  if (phase_matrix.has_value() and other.phase_matrix.has_value()) {
    phase_matrix.value() += other.phase_matrix.value();
  }

  extinction_matrix += other.extinction_matrix;
  absorption_vector += other.absorption_vector;

  return *this;
}

scattering::PhaseMatrixData<scattering::Format::TRO, scattering::Representation::Spectral>
ScatteringTroSpectralVector::to_general(const CompactPlanarMuelmatMatrix& phase_matrix,
                                        const std::shared_ptr<Vector>&    f_grid_ptr) {
  ARTS_USER_ERROR_IF(not f_grid_ptr, "f grid must be provided")

  ARTS_USER_ERROR_IF(static_cast<Size>(phase_matrix.nrows()) != f_grid_ptr->size(),
                     "Phase matrix and f grid must have the same size")

  const auto t_grid = std::make_shared<Vector>(Vector{0.0});

  const Index l = phase_matrix.ncols() - 1;
  ARTS_USER_ERROR_IF(l < 0, "Legendre degree must be non-negative, is {}", l);

  scattering::PhaseMatrixData<scattering::Format::TRO, scattering::Representation::Spectral> pm{t_grid, f_grid_ptr, l};

  pm[0, joker, joker] = phase_matrix;

  return pm;
}

scattering::ExtinctionMatrixData<scattering::Format::TRO, scattering::Representation::Spectral>
ScatteringTroSpectralVector::to_general(const CompactPlanarPropmatVector& extinction_matrix,
                                        const std::shared_ptr<Vector>&    f_grid_ptr) {
  ARTS_USER_ERROR_IF(not f_grid_ptr, "f grid must be provided")

  ARTS_USER_ERROR_IF(extinction_matrix.size() != f_grid_ptr->size(),
                     "Extinction matrix and f grid must have the same size")

  const auto t_grid = std::make_shared<Vector>(Vector{0.0});

  scattering::ExtinctionMatrixData<scattering::Format::TRO, scattering::Representation::Spectral> emd{t_grid,
                                                                                                      f_grid_ptr};

  emd[0, joker] = extinction_matrix;

  return emd;
}

scattering::AbsorptionVectorData<scattering::Format::TRO, scattering::Representation::Spectral>
ScatteringTroSpectralVector::to_general(const CompactPlanarStokvecVector& absorption_vector,
                                        const std::shared_ptr<Vector>&    f_grid_ptr) {
  ARTS_USER_ERROR_IF(not f_grid_ptr, "f grid must be provided")

  ARTS_USER_ERROR_IF(absorption_vector.size() != f_grid_ptr->size(),
                     "Absorption vector and f grid must have the same size")

  const auto t_grid = std::make_shared<Vector>(Vector{0.0});

  scattering::AbsorptionVectorData<scattering::Format::TRO, scattering::Representation::Spectral> av{t_grid,
                                                                                                     f_grid_ptr};

  av[0, joker] = absorption_vector;

  return av;
}

ScatteringTroSpectralVector::general_t ScatteringTroSpectralVector::to_general(const std::shared_ptr<Vector>& f) const {
  ARTS_USER_ERROR_IF(not phase_matrix.has_value(), "No phase matrix")

  if (f) {
    return {.phase_matrix      = to_general(phase_matrix.value(), f),
            .extinction_matrix = to_general(extinction_matrix, f),
            .absorption_vector = to_general(absorption_vector, f)};
  }

  const Index nf = phase_matrix.value().nrows();
  const auto  fs = std::make_shared<Vector>(nlinspace(1, static_cast<Numeric>(nf), nf));
  return {.phase_matrix      = to_general(phase_matrix.value(), fs),
          .extinction_matrix = to_general(extinction_matrix, fs),
          .absorption_vector = to_general(absorption_vector, fs)};
}

ScatteringTroSpectralVector::gridded_t ScatteringTroSpectralVector::to_lab_frame(
    std::shared_ptr<const Vector>                      za_inc_grid,
    std::shared_ptr<const Vector>                      delta_aa_grid,
    std::shared_ptr<const scattering::ZenithAngleGrid> za_scat_grid_new) const {
  return to_general().to_lab_frame(std::move(za_inc_grid), std::move(delta_aa_grid), std::move(za_scat_grid_new));
}

void xml_io_stream<ScatteringGeneralSpectralTRO>::write(std::ostream&                       os,
                                                        const ScatteringGeneralSpectralTRO& x,
                                                        bofstream*                          pbofs,
                                                        std::string_view                    name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.f, pbofs);
  xml_write_to_stream(os, x.degree, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<ScatteringGeneralSpectralTRO>::read(std::istream&                 is,
                                                       ScatteringGeneralSpectralTRO& x,
                                                       bifstream*                    pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.f, pbifs);
  xml_read_from_stream(is, x.degree, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
