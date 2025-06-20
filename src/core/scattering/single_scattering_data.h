#pragma once

#include <optional>
#include <ranges>

#include "mie.h"
#include "optproperties.h"
#include "scattering/absorption_vector.h"
#include "scattering/extinction_matrix.h"
#include "scattering/phase_matrix.h"
#include "scattering/psd.h"


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
template <std::floating_point Scalar,
          Format format,
          Representation repr>
struct SingleScatteringData {
 public:

  static SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>
  liquid_sphere(const StridedVectorView &t_grid,
                const StridedVectorView &f_grid,
                Numeric diameter,
                const ZenithAngleGrid &za_grid) {

    auto t_grid_ptr = std::make_shared<Vector>(t_grid);
    auto f_grid_ptr = std::make_shared<Vector>(f_grid);
    auto za_grid_ptr = std::make_shared<ZenithAngleGrid>(za_grid);

    PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded> phase_matrix(t_grid_ptr,
                                                                                f_grid_ptr,
                                                                                za_grid_ptr);
    ExtinctionMatrixData<Numeric, Format::TRO, Representation::Gridded> extinction_matrix(t_grid_ptr, f_grid_ptr);
    AbsorptionVectorData<Numeric, Format::TRO, Representation::Gridded> absorption_vector(t_grid_ptr, f_grid_ptr);
    BackscatterMatrixData<Numeric, Format::TRO> backscatter_matrix(t_grid_ptr, f_grid_ptr);
    ForwardscatterMatrixData<Numeric, Format::TRO> forwardscatter_matrix(t_grid_ptr, f_grid_ptr);

    for (size_t temp_ind = 0; temp_ind < t_grid_ptr->size(); ++temp_ind) {
      Numeric temp = t_grid_ptr->operator[](temp_ind);
      for (size_t freq_ind = 0; freq_ind < f_grid_ptr->size(); ++freq_ind) {
        Numeric freq = f_grid_ptr->operator[](freq_ind);
        auto sphere = MieSphere<Scalar>::Liquid(freq, temp, diameter / 2.0, grid_vector(*za_grid_ptr));
        phase_matrix[temp_ind, freq_ind] = sphere.get_scattering_matrix_compact();
        extinction_matrix[temp_ind, freq_ind] = sphere.get_extinction_coeff();
        absorption_vector[temp_ind, freq_ind] = sphere.get_absorption_coeff();
      }
    }

    auto pprops = ParticleProperties{
      "Mie Sphere",
      "ARTS Mie solver",
      "Ellison (2007)",
      1e3 * 4.0 * Constant::pi * std::pow(diameter / 2.0, 3),
      diameter,
      diameter
    };

    return SingleScatteringData(pprops,
                                phase_matrix,
                                extinction_matrix,
                                absorption_vector,
                                backscatter_matrix,
                                forwardscatter_matrix);

  }


  static SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>
  from_legacy_tro(::SingleScatteringData ssd, ::ScatteringMetaData smd) {
    ARTS_USER_ERROR_IF(
        ssd.ptype != PType::PTYPE_TOTAL_RND,
        "Converting legacy ARTS single scattering data to TRO format requires"
        " the input data to have PType::PTYPR_TOTAL_RND.");

    auto t_grid = std::make_shared<Vector>(ssd.T_grid);
    auto f_grid = std::make_shared<Vector>(ssd.f_grid);
    auto za_scat_grid = std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(ssd.za_grid));

    PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded>
        phase_matrix(t_grid, f_grid, za_scat_grid);
    ExtinctionMatrixData<Numeric, Format::TRO, Representation::Gridded>
        extinction_matrix(t_grid, f_grid);
    AbsorptionVectorData<Numeric, Format::TRO, Representation::Gridded>
        absorption_vector(t_grid, f_grid);

    for (Size i_t = 0; i_t < t_grid->size(); ++i_t) {
      for (Size i_f = 0; i_f < f_grid->size(); ++i_f) {
        for (Index i_za_scat = 0; i_za_scat < grid_size(*za_scat_grid);
             ++i_za_scat) {
          for (Index i_s = 0; i_s < phase_matrix.n_stokes_coeffs; ++i_s) {
            phase_matrix[i_t, i_f, i_za_scat, i_s] =
                ssd.pha_mat_data[i_f, i_t, i_za_scat, 0, 0, 0, i_s];
          }
        }
        for (Index i_s = 0; i_s < absorption_vector.n_stokes_coeffs; ++i_s) {
          absorption_vector[i_t, i_f, i_s] =
              ssd.abs_vec_data[i_f, i_t, 0, 0, i_s];
        }
        for (Index i_s = 0; i_s < extinction_matrix.n_stokes_coeffs; ++i_s) {
          extinction_matrix[i_t, i_f, i_s] =
              ssd.ext_mat_data[i_f, i_t, 0, 0, i_s];
        }
      }
    }

    //auto backscatter_matrix = BackscatterMatrixData<Numeric, Format::TRO>{t_grid, f_grid,};// = phase_matrix.extract_backscatter_matrix();
    auto backscatter_matrix = phase_matrix.extract_backscatter_matrix();
    auto forwardscatter_matrix = ForwardscatterMatrixData<Numeric, Format::TRO>{t_grid, f_grid};// = phase_matrix.extract_forwardscatter_matrix();

    auto properties = ParticleProperties{
      smd.description, smd.source, smd.refr_index, smd.mass, smd.diameter_volume_equ, smd.diameter_max
    };

    return SingleScatteringData<Numeric, Format::TRO, Representation::Gridded>(properties,
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
      PhaseMatrixData<Scalar, format, repr> phase_matrix_,
      ExtinctionMatrixData<Scalar, format, repr> extinction_matrix_,
      AbsorptionVectorData<Scalar, format, repr> absorption_vector_,
      BackscatterMatrixData<Scalar, format> backscatter_matrix_,
      ForwardscatterMatrixData<Scalar, format>
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
      PhaseMatrixData<Scalar, format, repr> phase_matrix_,
      ExtinctionMatrixData<Scalar, format, repr> extinction_matrix_,
      AbsorptionVectorData<Scalar, format, repr> absorption_vector_,
      BackscatterMatrixData<Scalar, format> backscatter_matrix_,
      ForwardscatterMatrixData<Scalar, format> forwardscatter_matrix_
                       )
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
  SingleScatteringData(
      std::optional<ParticleProperties> properties_,
      std::optional<PhaseMatrixData<Scalar, format, repr>> phase_matrix_,
      ExtinctionMatrixData<Scalar, format, repr> extinction_matrix_,
      AbsorptionVectorData<Scalar, format, repr> absorption_vector_,
      BackscatterMatrixData<Scalar, format> backscatter_matrix_,
      ForwardscatterMatrixData<Scalar, format> forwardscatter_matrix_
                       )
      : properties(properties_),
        phase_matrix(phase_matrix_),
        extinction_matrix(extinction_matrix_),
        absorption_vector(absorption_vector_),
        backscatter_matrix(backscatter_matrix_),
        forwardscatter_matrix(forwardscatter_matrix_) {}

  SingleScatteringData(const SingleScatteringData &) = default;

  constexpr Format get_format() const noexcept {
    return format;
  }

   constexpr Representation get_representation() const noexcept {
    return repr;
  }

  std::optional<Numeric> get_mass() const {
    auto extract_size = [](const ParticleProperties &part_props) {return part_props.mass;};
    return properties.transform(extract_size);
  }

  std::optional<Numeric> get_size(SizeParameter param) const {auto extract_size = [&param](const ParticleProperties &part_props) {
      if (param ==  SizeParameter::Mass) {
        return part_props.mass;
      } else if (param == SizeParameter::DMax) {
        return part_props.d_max;
      } else if (param == SizeParameter::DVeq) {
        return part_props.d_veq;
      }
      ARTS_USER_ERROR("Encountered unsupported size parameter.");
    };
    return properties.transform(extract_size);
  }

  SingleScatteringData regrid(const ScatteringDataGrids& grids,
                              const RegridWeights& weights) const {
    return SingleScatteringData(properties,
                                phase_matrix.transform([&grids](const auto& pm) {return pm.regrid(grids);}),
                                extinction_matrix.regrid(grids, weights),
                                absorption_vector.regrid(grids, weights),
                                backscatter_matrix.regrid(grids, weights),
                                forwardscatter_matrix.regrid(grids, weights));
  }

  SingleScatteringData regrid(const ScatteringDataGrids& grids) const {
    return SingleScatteringData(properties,
                                phase_matrix.transform([&grids](const auto& pm) {return pm.regrid(grids);}),
                                extinction_matrix.regrid(grids),
                                absorption_vector.regrid(grids),
                                backscatter_matrix.regrid(grids),
                                forwardscatter_matrix.regrid(grids));
  }

  SingleScatteringData<Numeric, format, Representation::Spectral> to_spectral(Index l, Index m=0) const {
    ARTS_USER_ERROR_IF((format == Format::TRO) && (m > 0),
                       "Order of SHT representation must be 0 for scattering data in TRO format");
    auto new_phase_matrix = phase_matrix.transform([&l, &m](const auto &pm) {return pm.to_spectral(l, m);});
    return SingleScatteringData<Numeric, format, Representation::Spectral>(properties,
                                                                           new_phase_matrix,
                                                                           extinction_matrix.to_spectral(),
                                                                           absorption_vector.to_spectral(),
                                                                           backscatter_matrix,
                                                                           forwardscatter_matrix);
  }

  SingleScatteringData<Numeric, format, Representation::Gridded> to_gridded() const {
    auto new_phase_matrix = phase_matrix.transform([](const auto &pm) {return pm.to_gridded();});
    return SingleScatteringData<Numeric, format, Representation::Gridded>(properties,
                                                                          new_phase_matrix,
                                                                          extinction_matrix.to_gridded(),
                                                                          absorption_vector.to_gridded(),
                                                                          backscatter_matrix,
                                                                          forwardscatter_matrix);
  }

  SingleScatteringData<Numeric, Format::ARO, Representation::Gridded> to_lab_frame(const ScatteringDataGrids& grids) const {
    if constexpr (format == Format::ARO) {
      return regrid(grids);
    }
    auto new_pm = phase_matrix.transform([&grids](const auto &pm) {
      return pm.to_lab_frame(grids.za_inc_grid, grids.aa_scat_grid, grids.za_scat_grid).regrid(grids);}
      );
    auto new_em = extinction_matrix.to_lab_frame(grids.za_inc_grid);
    auto new_av = absorption_vector.to_lab_frame(grids.za_inc_grid);
    auto new_bsm = BackscatterMatrixData<Numeric, Format::ARO>(backscatter_matrix, grids.za_inc_grid);
    auto new_fsm = ForwardscatterMatrixData<Numeric, Format::ARO>(forwardscatter_matrix, grids.za_inc_grid);
    return SingleScatteringData<Numeric, Format::ARO, Representation::Gridded>(properties,
                                                                               new_pm,
                                                                               new_em.regrid(grids),
                                                                               new_av.regrid(grids),
                                                                               new_bsm,
                                                                               new_fsm);
  }

  std::optional<ParticleProperties> properties;
  std::optional<PhaseMatrixData<Scalar, format, repr>> phase_matrix;
  ExtinctionMatrixData<Scalar, format, repr> extinction_matrix;
  AbsorptionVectorData<Scalar, format, repr> absorption_vector;
  BackscatterMatrixData<Scalar, format> backscatter_matrix;
  ForwardscatterMatrixData<Scalar, format> forwardscatter_matrix;

};

template <std::floating_point Scalar,
          Format format,
          Representation repr>
class ArrayOfSingleScatteringData
    : public std::vector<
          SingleScatteringData<Scalar, format, repr>> {
 public:
  //ArrayOfSingleScatteringData regrid(ScatteringDataGrids grids) {
  //  return ArrayOfSingleScatteringData{};
  //}

  //  ArrayOfSingleScatteringData<Scalar, format, Representation::Gridded> to_gridded() {
  //    return ArrayOfSingleScatteringData{};
  //  }
  //
  //  SingleScatteringData<Scalar, format, repr> calculate_bulk_properties(Vector pnd,
  //                                                                                   Numeric temperature,
  //                                                                                   bool include_phase_matrix=true)
  //  {
  //    return SingleScatteringData<Scalar, format, repr>{};
  //  }

 private:
};

template <std::floating_point Scalar,
          Format format,
          Representation repr>
std::ostream &operator<<(
    std::ostream &out,
    SingleScatteringData<Scalar, format, repr>) {
  out << "SingleScatteringData" << '\n';
  out << "\t Format:          " << format << '\n';
  out << "\t Representation:  " << repr << '\n';
  return out;
}

}  // namespace scattering
