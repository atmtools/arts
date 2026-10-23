#pragma once

#include <optional>
#include <ranges>
#include <utility>

#include "absorption_vector.h"
#include "enumsSizeParameter.h"
#include "extinction_matrix.h"
#include "mie.h"
#include "optproperties.h"
#include "phase_matrix.h"
#include "psd.h"

namespace scattering {

struct ParticleProperties {
  std::string name             = "";
  std::string source           = "";
  std::string refractive_index = "";
  double      mass             = 0.0;
  double      d_veq            = 0.0;
  double      d_max            = 0.0;
};

/** Single scattering data.
 *
 * The SingleScatteringData class is a container that holds single scattering
 * data of a single scattering particle. Optionally, it also holds additional
 * properties of the corresponding particle.
 *
 * The single scattering data comprises data for the phase and extinction
 * matrix, the absorption vector and back- and forward scattering matrices
 * for a range of temperatures, frequencies and directions. Since the phase
 * matrix data stands for the largest part of the scattering data it may
 * be empty in cases where phase matrix data is not needed.
 *
 */
template <Format format, Representation repr> struct SingleScatteringData {
 public:
  static SingleScatteringData<Format::TRO, Representation::Gridded> liquid_sphere(const StridedVectorView &t_grid,
                                                                                  const StridedVectorView &f_grid,
                                                                                  Numeric                  diameter,
                                                                                  const ZenithAngleGrid   &za_grid) {
    ComplexMatrix refractive_index(t_grid.size(), f_grid.size());
    for (Index it = 0; it < t_grid.ncols(); ++it)
      for (Index jf = 0; jf < f_grid.ncols(); ++jf)
        refractive_index[it, jf] = refr_index_water_ellison07(f_grid[jf], t_grid[it]);
    auto result                         = sphere(t_grid, f_grid, diameter, za_grid, refractive_index, 1e3);
    result.properties->refractive_index = "Ellison (2007)";
    return result;
  }

  static SingleScatteringData<Format::TRO, Representation::Gridded> sphere(const StridedVectorView &t_grid,
                                                                           const StridedVectorView &f_grid,
                                                                           Numeric                  diameter,
                                                                           const ZenithAngleGrid   &za_grid,
                                                                           const ComplexMatrix     &refractive_index,
                                                                           Numeric                  density) {
    if (refractive_index.nrows() != t_grid.ncols() || refractive_index.ncols() != f_grid.ncols())
      arts::user_error("Refractive-index shape must match temperature/frequency grids.");
    if (!std::isfinite(diameter) || diameter <= 0 || !std::isfinite(density) || density <= 0)
      arts::user_error("Diameter and density must be finite and positive.");
    for (auto t : t_grid) {
      if (!std::isfinite(t) || t <= 0) arts::user_error("Temperatures must be finite and positive.");
    }
    for (auto f : f_grid) {
      if (!std::isfinite(f) || f <= 0) arts::user_error("Frequencies must be finite and positive.");
    }
    for (auto row : refractive_index)
      for (auto m : row) {
        if (!std::isfinite(m.real()) || m.real() <= 0 || !std::isfinite(m.imag()) || m.imag() < 0)
          arts::user_error("Refractive index must have positive real and nonnegative imaginary parts.");
      }
    auto t_grid_ptr  = std::make_shared<Vector>(t_grid);
    auto f_grid_ptr  = std::make_shared<Vector>(f_grid);
    auto za_grid_ptr = std::make_shared<ZenithAngleGrid>(za_grid);

    PhaseMatrixData<Format::TRO, Representation::Gridded>      phase_matrix(t_grid_ptr, f_grid_ptr, za_grid_ptr);
    ExtinctionMatrixData<Format::TRO, Representation::Gridded> extinction_matrix(t_grid_ptr, f_grid_ptr);
    AbsorptionVectorData<Format::TRO, Representation::Gridded> absorption_vector(t_grid_ptr, f_grid_ptr);
    BackscatterMatrixData<Format::TRO>                         backscatter_matrix(t_grid_ptr, f_grid_ptr);
    ForwardscatterMatrixData<Format::TRO>                      forwardscatter_matrix(t_grid_ptr, f_grid_ptr);

    // Evaluate exact endpoints even when the requested angular grid omits them.
    const auto &requested_angles = grid_vector(*za_grid_ptr);
    Vector      angles(requested_angles.size() + 2);
    for (Index i = 0; i < requested_angles.ncols(); ++i) angles[i] = requested_angles[i];
    angles[requested_angles.size()]     = 0;
    angles[requested_angles.size() + 1] = 180;
    for (size_t temp_ind = 0; temp_ind < t_grid_ptr->size(); ++temp_ind) {
      for (size_t freq_ind = 0; freq_ind < f_grid_ptr->size(); ++freq_ind) {
        Numeric freq = f_grid_ptr->operator[](freq_ind);
        auto    sphere =
            MieSphere(Constant::speed_of_light / freq, diameter / 2.0, refractive_index[temp_ind, freq_ind], angles);
        auto optical = sphere.get_scattering_matrix_compact();
        for (Index ia = 0; ia < requested_angles.ncols(); ++ia) phase_matrix[temp_ind, freq_ind, ia] = optical[ia];
        forwardscatter_matrix[temp_ind, freq_ind] = optical[requested_angles.size()];
        backscatter_matrix[temp_ind, freq_ind]    = optical[requested_angles.size() + 1];
        extinction_matrix[temp_ind, freq_ind].A() = sphere.get_extinction_coeff();
        absorption_vector[temp_ind, freq_ind].I() = sphere.get_absorption_coeff();
      }
    }

    auto pprops = ParticleProperties{.name             = "Mie Sphere",
                                     .source           = "ARTS Mie solver",
                                     .refractive_index = "User-supplied temperature/frequency grid",
                                     .mass             = density * Constant::pi / 6 * std::pow(diameter, 3),
                                     .d_veq            = diameter,
                                     .d_max            = diameter};

    return SingleScatteringData(
        pprops, phase_matrix, extinction_matrix, absorption_vector, backscatter_matrix, forwardscatter_matrix);
  }

  static SingleScatteringData<Format::TRO, Representation::Gridded> from_legacy_tro(::SingleScatteringData ssd,
                                                                                    ::ScatteringMetaData   smd) {
    if (ssd.ptype != PType::PTYPE_TOTAL_RND)
      arts::user_error(
          "Converting legacy ARTS single scattering data to TRO format requires"
          " the input data to have PType::PTYPR_TOTAL_RND.");

    auto t_grid       = std::make_shared<Vector>(ssd.T_grid);
    auto f_grid       = std::make_shared<Vector>(ssd.f_grid);
    auto za_scat_grid = std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(ssd.za_grid));

    PhaseMatrixData<Format::TRO, Representation::Gridded>      phase_matrix(t_grid, f_grid, za_scat_grid);
    ExtinctionMatrixData<Format::TRO, Representation::Gridded> extinction_matrix(t_grid, f_grid);
    AbsorptionVectorData<Format::TRO, Representation::Gridded> absorption_vector(t_grid, f_grid);

    for (Size i_t = 0; i_t < t_grid->size(); ++i_t) {
      for (Size i_f = 0; i_f < f_grid->size(); ++i_f) {
        for (Index i_za_scat = 0; i_za_scat < grid_size(*za_scat_grid); ++i_za_scat) {
          phase_matrix[i_t, i_f, i_za_scat] =
              rtepack::compact_planar_muelmat{ssd.pha_mat_data[i_f, i_t, i_za_scat, 0, 0, 0, Range{0, 6}]};
        }
        absorption_vector[i_t, i_f].I() = ssd.abs_vec_data[i_f, i_t, 0, 0, 0];
        extinction_matrix[i_t, i_f].A() = ssd.ext_mat_data[i_f, i_t, 0, 0, 0];
      }
    }

    //auto backscatter_matrix = BackscatterMatrixData<Format::TRO>{t_grid, f_grid,};// = phase_matrix.extract_backscatter_matrix();
    auto backscatter_matrix = phase_matrix.extract_backscatter_matrix();
    auto forwardscatter_matrix =
        ForwardscatterMatrixData<Format::TRO>{t_grid, f_grid};  // = phase_matrix.extract_forwardscatter_matrix();

    auto properties = ParticleProperties{.name             = smd.description,
                                         .source           = smd.source,
                                         .refractive_index = smd.refr_index,
                                         .mass             = smd.mass,
                                         .d_veq            = smd.diameter_volume_equ,
                                         .d_max            = smd.diameter_max};

    return SingleScatteringData<Format::TRO, Representation::Gridded>(
        properties, phase_matrix, extinction_matrix, absorption_vector, backscatter_matrix, forwardscatter_matrix);
  }

  static SingleScatteringData<Format::ARO, Representation::Gridded> from_legacy_aro(::SingleScatteringData ssd,
                                                                                    ::ScatteringMetaData   smd) {
    if (ssd.ptype != PType::PTYPE_AZIMUTH_RND)
      arts::user_error("Converting legacy scattering data to ARO format requires PType::PTYPE_AZIMUTH_RND.");
    if (ssd.pha_mat_data.extent(5) != 1)
      arts::user_error("Legacy ARO conversion requires a singleton incident-azimuth dimension.");

    auto t_grid      = std::make_shared<Vector>(ssd.T_grid);
    auto f_grid      = std::make_shared<Vector>(ssd.f_grid);
    auto za_inc_grid = std::make_shared<Vector>(ssd.za_grid);
    if (ssd.aa_grid.empty() || ssd.aa_grid.front() != 0.0 || ssd.aa_grid.back() != 180.0)
      arts::user_error("Legacy ARO azimuth grid must span 0 to 180 degrees.");
    Vector signed_delta_aa(2 * ssd.aa_grid.size() - 1);
    for (Size i = 1; i < ssd.aa_grid.size(); ++i) { signed_delta_aa[ssd.aa_grid.size() - 1 - i] = -ssd.aa_grid[i]; }
    for (Size i = 0; i < ssd.aa_grid.size(); ++i) { signed_delta_aa[ssd.aa_grid.size() - 1 + i] = ssd.aa_grid[i]; }
    auto delta_aa_grid = std::make_shared<Vector>(std::move(signed_delta_aa));
    auto za_scat_grid  = std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(ssd.za_grid));

    PhaseMatrixData<Format::ARO, Representation::Gridded> phase_matrix(
        t_grid, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);
    ExtinctionMatrixData<Format::ARO, Representation::Gridded> extinction_matrix(t_grid, f_grid, za_inc_grid);
    AbsorptionVectorData<Format::ARO, Representation::Gridded> absorption_vector(t_grid, f_grid, za_inc_grid);

    for (Size i_t = 0; i_t < t_grid->size(); ++i_t) {
      for (Size i_f = 0; i_f < f_grid->size(); ++i_f) {
        for (Size i_za_inc = 0; i_za_inc < za_inc_grid->size(); ++i_za_inc) {
          for (Size i_delta_aa = 0; i_delta_aa < delta_aa_grid->size(); ++i_delta_aa) {
            const bool negative = (*delta_aa_grid)[i_delta_aa] < 0.0;
            const Size source_aa =
                negative ? ssd.aa_grid.size() - 1 - i_delta_aa : i_delta_aa - (ssd.aa_grid.size() - 1);
            for (Size i_za_scat = 0; i_za_scat < static_cast<Size>(grid_size(*za_scat_grid)); ++i_za_scat) {
              for (Size i_s = 0; i_s < 16; ++i_s) {
                const bool odd =
                    i_s == 2 || i_s == 3 || i_s == 6 || i_s == 7 || i_s == 8 || i_s == 9 || i_s == 12 || i_s == 13;
                phase_matrix[i_t, i_f, i_za_inc, i_delta_aa, i_za_scat].data[i_s] =
                    (negative && odd ? -1.0 : 1.0) * ssd.pha_mat_data[i_f, i_t, i_za_scat, source_aa, i_za_inc, 0, i_s];
              }
            }
          }
          extinction_matrix[i_t, i_f, i_za_inc] =
              rtepack::compact_planar_propmat{ssd.ext_mat_data[i_f, i_t, i_za_inc, 0, Range{0, 3}]};
          absorption_vector[i_t, i_f, i_za_inc] =
              rtepack::compact_planar_stokvec{ssd.abs_vec_data[i_f, i_t, i_za_inc, 0, Range{0, 2}]};
        }
      }
    }

    auto properties            = ParticleProperties{.name             = smd.description,
                                                    .source           = smd.source,
                                                    .refractive_index = smd.refr_index,
                                                    .mass             = smd.mass,
                                                    .d_veq            = smd.diameter_volume_equ,
                                                    .d_max            = smd.diameter_max};
    auto backscatter_matrix    = phase_matrix.extract_backscatter_matrix();
    auto forwardscatter_matrix = phase_matrix.extract_forwardscatter_matrix();
    return SingleScatteringData<Format::ARO, Representation::Gridded>(
        properties, phase_matrix, extinction_matrix, absorption_vector, backscatter_matrix, forwardscatter_matrix);
  }

  /** Create SingleScatteringData container without particle propreties.
   *
   * @param phase_matrix_ The phase matrix data.
   * @param extinction_matrix_ The extinction matrix data.
   * @param absorption_vector_ The absorption vector data.
   * @param backscatter_matrix_ The backscatter matrix.
   * @param forwardscatter_matrix_ The forwardscatter matrix.
   */
  SingleScatteringData(PhaseMatrixData<format, repr>      phase_matrix_,
                       ExtinctionMatrixData<format, repr> extinction_matrix_,
                       AbsorptionVectorData<format, repr> absorption_vector_,
                       BackscatterMatrixData<format>      backscatter_matrix_,
                       ForwardscatterMatrixData<format>   forwardscatter_matrix_)
      : phase_matrix(phase_matrix_),
        extinction_matrix(extinction_matrix_),
        absorption_vector(absorption_vector_),
        backscatter_matrix(backscatter_matrix_),
        forwardscatter_matrix(forwardscatter_matrix_) {}

  /** Create SingleScatteringDat container with particle propreties.
   *
   * @param properties_ The properties of the particle.
   * @param phase_matrix_ The phase matrix data.
   * @param extinction_matrix_ The extinction matrix data.
   * @param absorption_vector_ The absorption vector data.
   * @param backscatter_matrix_ The backscatter matrix.
   * @param forwardscatter_matrix_ The forwardscatter matrix.
   */
  SingleScatteringData(ParticleProperties                 properties_,
                       PhaseMatrixData<format, repr>      phase_matrix_,
                       ExtinctionMatrixData<format, repr> extinction_matrix_,
                       AbsorptionVectorData<format, repr> absorption_vector_,
                       BackscatterMatrixData<format>      backscatter_matrix_,
                       ForwardscatterMatrixData<format>   forwardscatter_matrix_)
      : properties(properties_),
        phase_matrix(phase_matrix_),
        extinction_matrix(extinction_matrix_),
        absorption_vector(absorption_vector_),
        backscatter_matrix(backscatter_matrix_),
        forwardscatter_matrix(forwardscatter_matrix_) {}

  /** Create SingleScatteringDat container with particle propreties.
   *
   * @param properties_ The properties of the particle.
   * @param phase_matrix_ The phase matrix data.
   * @param extinction_matrix_ The extinction matrix data.
   * @param absorption_vector_ The absorption vector data.
   * @param backscatter_matrix_ The backscatter matrix.
   * @param forwardscatter_matrix_ The forwardscatter matrix.
   */
  SingleScatteringData(std::optional<ParticleProperties>            properties_,
                       std::optional<PhaseMatrixData<format, repr>> phase_matrix_,
                       ExtinctionMatrixData<format, repr>           extinction_matrix_,
                       AbsorptionVectorData<format, repr>           absorption_vector_,
                       BackscatterMatrixData<format>                backscatter_matrix_,
                       ForwardscatterMatrixData<format>             forwardscatter_matrix_)
      : properties(std::move(properties_)),
        phase_matrix(std::move(phase_matrix_)),
        extinction_matrix(std::move(extinction_matrix_)),
        absorption_vector(std::move(absorption_vector_)),
        backscatter_matrix(std::move(backscatter_matrix_)),
        forwardscatter_matrix(std::move(forwardscatter_matrix_)) {}

  SingleScatteringData()                                        = default;
  SingleScatteringData(const SingleScatteringData &)            = default;
  SingleScatteringData(SingleScatteringData &&)                 = default;
  SingleScatteringData &operator=(const SingleScatteringData &) = default;
  SingleScatteringData &operator=(SingleScatteringData &&)      = default;

  static constexpr Format get_format() noexcept { return format; }

  static constexpr Representation get_representation() noexcept { return repr; }

  [[nodiscard]] std::optional<Numeric> get_mass() const {
    auto extract_size = [](const ParticleProperties &part_props) { return part_props.mass; };
    return properties.transform(extract_size);
  }

  [[nodiscard]] std::optional<Numeric> get_size(SizeParameter param) const {
    auto extract_size = [&param](const ParticleProperties &part_props) {
      switch (param) {
        using enum SizeParameter;
        case Mass: return part_props.mass;
        case DMax: return part_props.d_max;
        case DVeq: return part_props.d_veq;
      }
      std::unreachable();
    };
    return properties.transform(extract_size);
  }

  [[nodiscard]] SingleScatteringData regrid(const ScatteringDataGrids &grids, const RegridWeights &weights) const {
    return SingleScatteringData(properties,
                                phase_matrix.transform([&grids](const auto &pm) { return pm.regrid(grids); }),
                                extinction_matrix.regrid(grids, weights),
                                absorption_vector.regrid(grids, weights),
                                backscatter_matrix.regrid(grids, weights),
                                forwardscatter_matrix.regrid(grids, weights));
  }

  [[nodiscard]] SingleScatteringData regrid(const ScatteringDataGrids &grids) const {
    return SingleScatteringData(properties,
                                phase_matrix.transform([&grids](const auto &pm) { return pm.regrid(grids); }),
                                extinction_matrix.regrid(grids),
                                absorption_vector.regrid(grids),
                                backscatter_matrix.regrid(grids),
                                forwardscatter_matrix.regrid(grids));
  }

  /** The spectral form: for TRO data the Legendre series to degree l (m must be 0), for ARO data the SHT of degree l
   *  and order m.
   *
   * See PhaseMatrixData::to_spectral of the format.  The extinction matrix,
   * absorption vector and back- and forward-scatter matrices are unchanged.
   */
  [[nodiscard]] SingleScatteringData<format, Representation::Spectral> to_spectral(Index l, Index m = 0) const {
    if ((format == Format::TRO) && (m > 0))
      arts::user_error("Order of SHT representation must be 0 for scattering data in TRO format");
    auto new_phase_matrix = phase_matrix.transform([&l, &m](const auto &pm) { return pm.to_spectral(l, m); });
    return SingleScatteringData<format, Representation::Spectral>(properties,
                                                                  new_phase_matrix,
                                                                  extinction_matrix.to_spectral(),
                                                                  absorption_vector.to_spectral(),
                                                                  backscatter_matrix,
                                                                  forwardscatter_matrix);
  }

  /** The azimuthal Fourier modes to m = max_mode of the laboratory-frame phase matrix of a TRO Legendre series
   *
   * The temperatures and frequencies are those of grids, interpolated; the
   * incidence and scattering zenith angles are those of grids, exactly.
   */
  [[nodiscard]] SingleScatteringData<Format::ARO, Representation::Fourier> to_lab_frame_fourier_modes(
      const ScatteringDataGrids &grids, Index max_mode) const
      requires(format == Format::TRO and repr == Representation::Spectral) {
    if (not grids.za_inc_grid or not grids.za_scat_grid)
      arts::user_error("Laboratory-frame Fourier modes need incidence and scattering zenith-angle grids");
    const ScatteringDataGrids tf_grids(grids.t_grid, grids.f_grid);
    auto                      new_pm = phase_matrix.transform([&](const auto &pm) {
      return pm.regrid(tf_grids).to_lab_frame_fourier_modes(grids.za_inc_grid, grids.za_scat_grid, max_mode);
    });
    auto                      new_em = extinction_matrix.regrid(tf_grids).to_lab_frame(grids.za_inc_grid).to_fourier();
    auto                      new_av = absorption_vector.regrid(tf_grids).to_lab_frame(grids.za_inc_grid).to_fourier();
    auto new_bsm = BackscatterMatrixData<Format::ARO>(backscatter_matrix.regrid(tf_grids), grids.za_inc_grid);
    auto new_fsm = ForwardscatterMatrixData<Format::ARO>(forwardscatter_matrix.regrid(tf_grids), grids.za_inc_grid);
    return SingleScatteringData<Format::ARO, Representation::Fourier>(
        properties, new_pm, new_em, new_av, new_bsm, new_fsm);
  }

  /** The azimuthal Fourier modes to m = max_mode of gridded ARO data, on their own zenith grids
   *
   * See PhaseMatrixData<ARO, Gridded>::to_fourier.
   */
  [[nodiscard]] SingleScatteringData<Format::ARO, Representation::Fourier> to_fourier(Index max_mode) const
      requires(format == Format::ARO and repr == Representation::Gridded) {
    return SingleScatteringData<Format::ARO, Representation::Fourier>(
        properties,
        phase_matrix.transform([max_mode](const auto &pm) { return pm.to_fourier(max_mode); }),
        extinction_matrix.to_fourier(),
        absorption_vector.to_fourier(),
        backscatter_matrix,
        forwardscatter_matrix);
  }

  /** The azimuthal Fourier modes to m = max_mode of SHT ARO data at the scattering zenith angles of a grid
   *
   * See PhaseMatrixData<ARO, Spectral>::to_fourier.
   */
  [[nodiscard]] SingleScatteringData<Format::ARO, Representation::Fourier> to_fourier(
      std::shared_ptr<const ZenithAngleGrid> za_scat_grid, Index max_mode) const
      requires(format == Format::ARO and repr == Representation::Spectral) {
    return SingleScatteringData<Format::ARO, Representation::Fourier>(
        properties,
        phase_matrix.transform([&](const auto &pm) { return pm.to_fourier(za_scat_grid, max_mode); }),
        extinction_matrix.to_fourier(),
        absorption_vector.to_fourier(),
        backscatter_matrix,
        forwardscatter_matrix);
  }

  [[nodiscard]] SingleScatteringData<format, Representation::Gridded> to_gridded() const {
    auto new_phase_matrix = phase_matrix.transform([](const auto &pm) { return pm.to_gridded(); });
    return SingleScatteringData<format, Representation::Gridded>(properties,
                                                                 new_phase_matrix,
                                                                 extinction_matrix.to_gridded(),
                                                                 absorption_vector.to_gridded(),
                                                                 backscatter_matrix,
                                                                 forwardscatter_matrix);
  }

  /** The series at the scattering angles of a grid, exactly
   *
   * @param za_scat_grid The scattering angles
   */
  [[nodiscard]] SingleScatteringData<format, Representation::Gridded> to_gridded(
      std::shared_ptr<const ZenithAngleGrid> za_scat_grid) const
      requires(format == Format::TRO and repr == Representation::Spectral) {
    auto new_phase_matrix = phase_matrix.transform([&](const auto &pm) { return pm.to_gridded(za_scat_grid); });
    return SingleScatteringData<format, Representation::Gridded>(properties,
                                                                 new_phase_matrix,
                                                                 extinction_matrix.to_gridded(),
                                                                 absorption_vector.to_gridded(),
                                                                 backscatter_matrix,
                                                                 forwardscatter_matrix);
  }

  [[nodiscard]] SingleScatteringData<Format::ARO, Representation::Gridded> to_lab_frame(
      const ScatteringDataGrids &grids) const {
    if constexpr (format == Format::ARO) { return regrid(grids); }
    // Interpolate the compact TRO quantities before expanding them.  The
    // laboratory-frame grids are then already exact, avoiding both redundant
    // work and degenerate singleton-grid interpolation at zenith/nadir.
    const ScatteringDataGrids tf_grids(grids.t_grid, grids.f_grid);
    auto                      new_pm  = phase_matrix.transform([&](const auto &pm) {
      if constexpr (repr == Representation::Spectral) {
        // The series at the exact scattering angle of every direction pair
        return pm.regrid(tf_grids).to_lab_frame(grids.za_inc_grid, grids.aa_scat_grid, grids.za_scat_grid);
      } else {
        // TRO's zenith grid is a scattering-angle grid, not the laboratory
        // output-zenith grid.  Preserve it while interpolating T/f.
        const ScatteringDataGrids tro_grids(grids.t_grid, grids.f_grid, pm.get_za_scat_grid());
        return pm.regrid(tro_grids).to_lab_frame(grids.za_inc_grid, grids.aa_scat_grid, grids.za_scat_grid);
      }
    });
    auto                      new_em  = extinction_matrix.regrid(tf_grids).to_gridded().to_lab_frame(grids.za_inc_grid);
    auto                      new_av  = absorption_vector.regrid(tf_grids).to_gridded().to_lab_frame(grids.za_inc_grid);
    auto                      new_bsm = BackscatterMatrixData<Format::ARO>(backscatter_matrix, grids.za_inc_grid);
    auto                      new_fsm = ForwardscatterMatrixData<Format::ARO>(forwardscatter_matrix, grids.za_inc_grid);
    return SingleScatteringData<Format::ARO, Representation::Gridded>(
        properties, new_pm, new_em, new_av, new_bsm, new_fsm);
  }

  std::optional<ParticleProperties>            properties;
  std::optional<PhaseMatrixData<format, repr>> phase_matrix;
  ExtinctionMatrixData<format, repr>           extinction_matrix;
  AbsorptionVectorData<format, repr>           absorption_vector;
  BackscatterMatrixData<format>                backscatter_matrix;
  ForwardscatterMatrixData<format>             forwardscatter_matrix;
};

/** The Legendre series to degree of gridded TRO data, and how well it represents them.
 *
 * The report's normalisation_error is (2 pi int F11 dcos(Theta) - (K11 - a1))
 * / K11 of the series, which a solver that takes its albedo from K11 - a1
 * needs to be small.  A free function: MSVC cannot instantiate the pair of a
 * class template inside that class.
 */
inline std::pair<SingleScatteringData<Format::TRO, Representation::Spectral>, LegendreReport> to_spectral_with_report(
    const SingleScatteringData<Format::TRO, Representation::Gridded> &ssd, Index degree) {
  if (not ssd.phase_matrix) arts::user_error("Scattering data without a phase matrix have no Legendre series");
  auto       spectral = ssd.to_spectral(degree);
  auto       report   = ssd.phase_matrix->legendre_report(*spectral.phase_matrix);
  const auto integral = spectral.phase_matrix->integrate_phase_matrix();
  for (Index i_t = 0; i_t < integral.extent(0); ++i_t) {
    for (Index i_f = 0; i_f < integral.extent(1); ++i_f) {
      const Numeric k11 = ssd.extinction_matrix[i_t, i_f].A(), a1 = ssd.absorption_vector[i_t, i_f].I();
      report.normalisation_error[i_t, i_f] = (integral[i_t, i_f].F11() - (k11 - a1)) / k11;
    }
  }
  return {std::move(spectral), std::move(report)};
}

template <Format format, Representation repr> class ArrayOfSingleScatteringData
    : public std::vector<SingleScatteringData<format, repr>> {
 public:
  //ArrayOfSingleScatteringData regrid(ScatteringDataGrids grids) {
  //  return ArrayOfSingleScatteringData{};
  //}

  //  ArrayOfSingleScatteringData<format, Representation::Gridded> to_gridded() {
  //    return ArrayOfSingleScatteringData{};
  //  }
  //
  //  SingleScatteringData<format, repr> calculate_bulk_properties(Vector pnd,
  //                                                                                   Numeric temperature,
  //                                                                                   bool include_phase_matrix=true)
  //  {
  //    return SingleScatteringData<format, repr>{};
  //  }

 private:
};

template <Format format, Representation repr>
std::ostream &operator<<(std::ostream &out, SingleScatteringData<format, repr>) {
  out << "SingleScatteringData" << '\n';
  out << "\t Format:          " << format << '\n';
  out << "\t Representation:  " << repr << '\n';
  return out;
}

}  // namespace scattering
