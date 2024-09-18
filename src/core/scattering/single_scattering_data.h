#pragma once

#include <optional>

#include "optproperties.h"
#include "scattering/absorption_vector.h"
#include "scattering/extinction_matrix.h"
#include "scattering/phase_matrix.h"

namespace scattering {

struct ParticleProperties {
  std::string name = "";
  std::string source = "";
  std::string refractive_index = "";
  double mass = 0.0;
  double d_veq = 0.0;
  double d_max = 0.0;
};

/** Single scattering data.
 *
 * The SingleScatteringData class is a container that holds single scattering
 * data of scattering particles. Optionally, it also holds additional
 * properties of the corresponding particle.
 *
 * The single scattering data comprises data for the phase and extinction
 * matrix, the absorption vector and back- and forward scattering matrices
 * for a range of temperatures, frequencies and directions. Since the phase
 * matrix data stands for the largest part of the scattering data it may
 * be empty in cases where phase matrix data is not needed.
 *
 */
template <std::floating_point Scalar,
          Format format,
          Representation repr,
          Index stokes_dim>
struct SingleScatteringData {
 public:
  static SingleScatteringData<Numeric, Format::TRO, Representation::Gridded, 4>
  from_legacy_tro(::SingleScatteringData ssd, ::ScatteringMetaData smd) {
    ARTS_USER_ERROR_IF(
        ssd.ptype != PType::PTYPE_TOTAL_RND,
        "Converting legacy ARTS single scattering data to TRO format requires"
        " the input data to have PType::PTYPR_TOTAL_RND.");

    auto t_grid = std::make_shared<Vector>(ssd.T_grid);
    auto f_grid = std::make_shared<Vector>(ssd.f_grid);
    auto za_scat_grid = std::make_shared<IrregularLatitudeGrid>(ssd.za_grid);

    PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded, 4>
        phase_matrix(t_grid, f_grid, za_scat_grid);
    ExtinctionMatrixData<Numeric, Format::TRO, Representation::Gridded, 4>
        extinction_matrix(t_grid, f_grid);
    AbsorptionVectorData<Numeric, Format::TRO, Representation::Gridded, 4>
        absorption_vector(t_grid, f_grid);

    for (Index i_t = 0; i_t < t_grid->size(); ++i_t) {
      for (Index i_f = 0; i_f < f_grid->size(); ++i_f) {
        for (Index i_za_scat = 0; i_za_scat < za_scat_grid->size();
             ++i_za_scat) {
          for (Index i_s = 0; i_s < phase_matrix.n_stokes_coeffs; ++i_s) {
            phase_matrix(i_t, i_f, i_za_scat, i_s) =
                ssd.pha_mat_data(i_f, i_t, i_za_scat, 0, 0, 0, i_s);
          }
        }
        for (Index i_s = 0; i_s < absorption_vector.n_stokes_coeffs; ++i_s) {
          absorption_vector(i_t, i_f, i_s) =
              ssd.abs_vec_data(i_f, i_t, 0, 0, i_s);
        }
        for (Index i_s = 0; i_s < extinction_matrix.n_stokes_coeffs; ++i_s) {
          extinction_matrix(i_t, i_f, i_s) =
              ssd.ext_mat_data(i_f, i_t, 0, 0, i_s);
        }
      }
    }

    auto backscatter_matrix = phase_matrix.extract_backscatter_matrix();
    auto forwardscatter_matrix = phase_matrix.extract_forwardscatter_matrix();

    auto properties = ParticleProperties{
      smd.description, smd.source, smd.refr_index, smd.mass, smd.diameter_volume_equ, smd.diameter_max
    };

    return SingleScatteringData<Numeric, Format::TRO, Representation::Gridded, 4>(properties,
                                                                                  phase_matrix,
                                                                                  extinction_matrix,
                                                                                  absorption_vector,
                                                                                  backscatter_matrix,
                                                                                  forwardscatter_matrix);
  }

  /** Create SingleScatteringData container without particle propreties.
   *
   * @param phase_matrix_ The phase matrix data.
   * @param extinction_matrix_ The extinction matrix data.
   * @param absorption_vector_ The absorption vector data.
   * @param backscatter_matrix_ The backscatter matrix.
   * @param forwardscatter_matrix_ The forwardscatter matrix.
   */
  SingleScatteringData(
      PhaseMatrixData<Scalar, format, repr, stokes_dim> phase_matrix_,
      ExtinctionMatrixData<Scalar, format, repr, stokes_dim> extinction_matrix_,
      AbsorptionVectorData<Scalar, format, repr, stokes_dim> absorption_vector_,
      BackscatterMatrixData<Scalar, format, stokes_dim> backscatter_matrix_,
      ForwardscatterMatrixData<Scalar, format, stokes_dim>
          forwardscatter_matrix_)
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
  SingleScatteringData(
      ParticleProperties properties_,
      PhaseMatrixData<Scalar, format, repr, stokes_dim> phase_matrix_,
      ExtinctionMatrixData<Scalar, format, repr, stokes_dim> extinction_matrix_,
      AbsorptionVectorData<Scalar, format, repr, stokes_dim> absorption_vector_,
      BackscatterMatrixData<Scalar, format, stokes_dim> backscatter_matrix_,
      ForwardscatterMatrixData<Scalar, format, stokes_dim> forwardscatter_matrix_
                       )
      : properties(properties_),
        phase_matrix(phase_matrix_),
        extinction_matrix(extinction_matrix_),
        absorption_vector(absorption_vector_),
        backscatter_matrix(backscatter_matrix_),
        forwardscatter_matrix(forwardscatter_matrix_) {}

  std::optional<ParticleProperties> properties;
  std::optional<PhaseMatrixData<Scalar, format, repr, stokes_dim>> phase_matrix;
  ExtinctionMatrixData<Scalar, format, repr, stokes_dim> extinction_matrix;
  AbsorptionVectorData<Scalar, format, repr, stokes_dim> absorption_vector;
  BackscatterMatrixData<Scalar, format, stokes_dim> backscatter_matrix;
  ForwardscatterMatrixData<Scalar, format, stokes_dim> forwardscatter_matrix;
};

template <std::floating_point Scalar,
          Format format,
          Representation repr,
          Index stokes_dim>
class ArrayOfSingleScatteringData
    : public std::vector<
          SingleScatteringData<Scalar, format, repr, stokes_dim>> {
 public:
  //ArrayOfSingleScatteringData regrid(ScatteringDataGrids grids) {
  //  return ArrayOfSingleScatteringData{};
  //}

  //  ArrayOfSingleScatteringData<Scalar, format, Representation::Gridded, stokes_dim> to_gridded() {
  //    return ArrayOfSingleScatteringData{};
  //  }
  //
  //  SingleScatteringData<Scalar, format, repr, stokes_dim> calculate_bulk_properties(Vector pnd,
  //                                                                                   Numeric temperature,
  //                                                                                   bool include_phase_matrix=true)
  //  {
  //    return SingleScatteringData<Scalar, format, repr, stokes_dim>{};
  //  }

 private:
};

template <std::floating_point Scalar,
          Format format,
          Representation repr,
          Index stokes_dim>
std::ostream &operator<<(
    std::ostream &out,
    SingleScatteringData<Scalar, format, repr, stokes_dim> ssd) {
  out << "SingleScatteringData" << std::endl;
  out << "\t Format:          " << format << std::endl;
  out << "\t Representation:  " << repr << std::endl;
  out << "\t Max. stokes dim: " << stokes_dim << std::endl;
  return out;
}

}  // namespace scattering
