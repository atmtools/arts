#include <core/scattering/particle_habit.h>
#include <core/scattering/single_scattering_data.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <python_interface.h>

#include "hpy_arts.h"
#include "py_macros.h"

NB_MAKE_OPAQUE(scattering::ZenithAngleGrid);

namespace Python {

template <typename Scalar>
void bind_phase_matrix_data_tro_gridded(py::module_& m,
                                        const std::string& class_name) {
  using PMD = scattering::PhaseMatrixData<Scalar,
                                          scattering::Format::TRO,
                                          scattering::Representation::Gridded>;

  py::class_<PMD, matpack::matpack_data<Scalar, 4>>(m, class_name.c_str())
      .def(py::init<std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>,
                    std::shared_ptr<const scattering::ZenithAngleGrid>>(),
           py::arg("t_grid"),
           py::arg("f_grid"),
           py::arg("za_scat_grid"))

      // Bind methods, such as `extract_backscatter_matrix` and `extract_forwardscatter_matrix`
      .def("extract_backscatter_matrix", &PMD::extract_backscatter_matrix)
      .def("extract_forwardscatter_matrix", &PMD::extract_forwardscatter_matrix)

      // Bind grid accessors
      .def("get_t_grid", &PMD::get_t_grid)
      .def("get_f_grid", &PMD::get_f_grid)
      .def("get_za_scat_grid", &PMD::get_za_scat_grid)

      .def("to_spectral", [](const PMD& obj) { obj.to_spectral(); })

      // Bind other member functions
      .def("integrate_phase_matrix", &PMD::integrate_phase_matrix)

      // Bind the extraction of stokes coefficients
      .def(
          "extract_stokes_coeffs",
          &PMD::
              extract_stokes_coeffs)  // You can specify the new stokes_dim here

      // Bind regrid method
      .def("regrid", &PMD::regrid);
}

template <typename Scalar>
void bind_phase_matrix_data_tro_spectral(py::module_& m,
                                         const std::string& class_name) {
  using PMD = scattering::PhaseMatrixData<Scalar,
                                          scattering::Format::TRO,
                                          scattering::Representation::Spectral>;
  py::class_<PMD, matpack::matpack_data<std::complex<Scalar>, 4>>(
      m, class_name.c_str())
      .def(py::init<>())
      .def(py::init<std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>,
                    std::shared_ptr<scattering::SHT>>(),
           py::arg("t_grid"),
           py::arg("f_grid"),
           py::arg("sht"))

      .def("extract_backscatter_matrix", &PMD::extract_backscatter_matrix)
      .def("extract_forwardscatter_matrix", &PMD::extract_forwardscatter_matrix)

      .def("get_t_grid", &PMD::get_t_grid)
      .def("get_f_grid", &PMD::get_f_grid)

      .def("to_gridded", [](const PMD& obj) { return obj.to_gridded(); })
      .def("integrate_phase_matrix", &PMD::integrate_phase_matrix)
      .def(
          "extract_stokes_coeffs",
          &PMD::
              extract_stokes_coeffs)  // You can specify the new stokes_dim here
      .def("regrid", &PMD::regrid);
}

template <typename Scalar, scattering::Representation repr>
void bind_absorption_vector_data_tro(py::module_& m, const std::string& name) {
  using AVD =
      scattering::AbsorptionVectorData<Scalar, scattering::Format::TRO, repr>;
  py::class_<AVD, matpack::matpack_data<Scalar, 3>>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>>(),
           "t_grid"_a,
           "f_grid"_a)
      .def("get_coeff_vector_view", &AVD::get_coeff_vector_view);
}

template <typename Scalar, scattering::Representation repr>
void bind_absorption_vector_data_aro(py::module_& m, const std::string& name) {
  using AVD =
      scattering::AbsorptionVectorData<Scalar, scattering::Format::ARO, repr>;
  py::class_<AVD, matpack::matpack_data<Scalar, 4>>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>>(),
           "t_grid"_a,
           "f_grid"_a,
           "za_inc_grid"_a)
      .def("get_coeff_vector_view", &AVD::get_coeff_vector_view)
      .def("regrid", &AVD::regrid, "grids"_a, "weights"_a);
}

template <typename Scalar, scattering::Representation repr>
void bind_extinction_matrix_data_tro(py::module_& m, const std::string& name) {
  using EMD =
      scattering::ExtinctionMatrixData<Scalar, scattering::Format::TRO, repr>;
  py::class_<EMD, matpack::matpack_data<Scalar, 3>>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>>(),
           "t_grid"_a,
           "f_grid"_a)
      .def("get_coeff_vector_view", &EMD::get_coeff_vector_view)
      .def("regrid", &EMD::regrid, "grids"_a, "weights"_a);
}

template <typename Scalar, scattering::Representation repr>
void bind_extinction_matrix_data_aro(py::module_& m, const std::string& name) {
  using EMD =
      scattering::ExtinctionMatrixData<Scalar, scattering::Format::ARO, repr>;
  py::class_<EMD, matpack::matpack_data<Scalar, 4>>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>>(),
           "t_grid"_a,
           "f_grid"_a,
           "za_inc_grid"_a)
      .def("get_coeff_vector_view", &EMD::get_coeff_vector_view)
      .def("regrid", &EMD::regrid, "grids"_a, "weights"_a);
}

template <typename Scalar,
          scattering::Format format,
          scattering::Representation repr>
void bind_single_scattering_data(py::module_& m, const std::string& name) {
  using SSDClass = scattering::SingleScatteringData<Scalar, format, repr>;

  py::class_<SSDClass>(m, name.c_str())
      .def(py::init<scattering::PhaseMatrixData<Scalar, format, repr>,
                    scattering::ExtinctionMatrixData<Scalar, format, repr>,
                    scattering::AbsorptionVectorData<Scalar, format, repr>,
                    scattering::BackscatterMatrixData<Scalar, format>,
                    scattering::ForwardscatterMatrixData<Scalar, format>>(),
           "phase_matrix"_a,
           "extinction_matrix"_a,
           "absorption_vector"_a,
           "backscatter_matrix"_a,
           "forwardscatter_matrix"_a)
      .def(py::init<scattering::ParticleProperties,
                    scattering::PhaseMatrixData<Scalar, format, repr>,
                    scattering::ExtinctionMatrixData<Scalar, format, repr>,
                    scattering::AbsorptionVectorData<Scalar, format, repr>,
                    scattering::BackscatterMatrixData<Scalar, format>,
                    scattering::ForwardscatterMatrixData<Scalar, format>>(),
           "properties"_a,
           "phase_matrix"_a,
           "extinction_matrix"_a,
           "absorption_vector"_a,
           "backscatter_matrix"_a,
           "forwardscatter_matrix"_a)
      .def_rw("properties", &SSDClass::properties)
      .def_rw("phase_matrix", &SSDClass::phase_matrix)
      .def_rw("extinction_matrix", &SSDClass::extinction_matrix)
      .def_rw("absorption_vector", &SSDClass::absorption_vector)
      .def_rw("backscatter_matrix", &SSDClass::backscatter_matrix)
      .def_rw("forwardscatter_matrix", &SSDClass::forwardscatter_matrix)
      .def_static(
          "from_legacy_tro", &SSDClass::from_legacy_tro, "ssd"_a, "smd"_a)
      .def("__repr__", [](const SSDClass& ssd) {
        std::ostringstream oss;
        oss << ssd;
        return oss.str();
      });
}

template <scattering::Format format, scattering::Representation repr>
void bind_bulk_scattering_properties(py::module_& m, const std::string& name) {
  py::class_<scattering::BulkScatteringProperties<format, repr>>(m,
                                                                 name.c_str())
      .def_rw("phase_matrix",
              &scattering::BulkScatteringProperties<format, repr>::phase_matrix)
      .def_rw("extinction_matrix",
              &scattering::BulkScatteringProperties<format,
                                                    repr>::extinction_matrix)
      .def_rw("absorption_vector",
              &scattering::BulkScatteringProperties<format,
                                                    repr>::absorption_vector);
}

void py_scattering_species(py::module_& m) try {  //
  // ScatSpeciesProperty
  //

  py::class_<ScatteringSpeciesProperty> ssp(m, "ScatteringSpeciesProperty");
  workspace_group_interface(ssp);
  ssp.def(py::init<std::string, ParticulateProperty>(), "Constructor")
      .def_rw("species_name", &ScatteringSpeciesProperty::species_name)
      .def_rw("pproperty", &ScatteringSpeciesProperty::pproperty);

  //
  // Modified gamma PSD
  //

  using BulkScatteringPropertiesTROSpectral =
      std::variant<scattering::BulkScatteringProperties<
          scattering::Format::TRO,
          scattering::Representation::Spectral>>;
  using BulkScatteringPropertiesTROGridded =
      std::variant<scattering::BulkScatteringProperties<
          scattering::Format::TRO,
          scattering::Representation::Gridded>>;
  using BulkScatteringPropertiesAROSpectral =
      std::variant<scattering::BulkScatteringProperties<
          scattering::Format::ARO,
          scattering::Representation::Spectral>>;
  using BulkScatteringPropertiesAROGridded =
      std::variant<scattering::BulkScatteringProperties<
          scattering::Format::ARO,
          scattering::Representation::Gridded>>;

  py::class_<MGDSingleMoment>(m, "MGDSingleMoment");
  py::class_<ScatteringHabit>(m, "ScatteringHabit");
  py::class_<HenyeyGreensteinScatterer>(m, "HenyeyGreensteinScatterer")
      .def(py::init<ExtSSACallback, Numeric>(), "func"_a, "g"_a)
      .def(py::init<>())
      .def(py::init<ScatteringSpeciesProperty,
                    ScatteringSpeciesProperty,
                    Numeric>())
      .def("get_bulk_scattering_properties_tro_spectral",
           [](const HenyeyGreensteinScatterer& hg,
              const AtmPoint& atm_point,
              const Vector& f_grid,
              Index l) {
             return BulkScatteringPropertiesTROSpectral{
                 hg.get_bulk_scattering_properties_tro_spectral(
                     atm_point, f_grid, l)};
           })
      .def("get_bulk_scattering_properties_tro_gridded",
           [](const HenyeyGreensteinScatterer& hg,
              const AtmPoint& atm_point,
              const Vector& f_grid,
              std::shared_ptr<scattering::ZenithAngleGrid> za_grid) {
             return BulkScatteringPropertiesTROGridded{
                 hg.get_bulk_scattering_properties_tro_gridded(
                     atm_point, f_grid, za_grid)};
           });

  py::class_<scattering::IrregularZenithAngleGrid>(m,
                                                   "IrregularZenithAngleGrid")
      .def(py::init<Vector>());
  py::class_<scattering::GaussLegendreGrid>(m, "GaussLegendreGrid")
      .def(py::init<Index>());
  py::class_<scattering::DoubleGaussGrid>(m, "DoubleGaussGrid")
      .def(py::init<Index>());
  py::class_<scattering::LobattoGrid>(m, "LobattoGrid").def(py::init<Index>());
  py::class_<scattering::FejerGrid>(m, "FejerGrid").def(py::init<Index>());

  py::class_<scattering::ZenithAngleGrid>(m, "ZenithAngleGrid")
      .def(py::init<scattering::IrregularZenithAngleGrid>())
      .def(py::init<scattering::GaussLegendreGrid>())
      .def(py::init<scattering::DoubleGaussGrid>())
      .def(py::init<scattering::LobattoGrid>())
      .def(py::init<scattering::FejerGrid>());

  py::class_<ArrayOfScatteringSpecies> aoss(m, "ArrayOfScatteringSpecies");
  aoss.def(py::init<>())
      .def(py::init_implicit<std::vector<ScatteringSpecies>>())
      .def("add", &ArrayOfScatteringSpecies::add)
      .def("get_bulk_scattering_properties_tro_spectral",
           [](const ArrayOfScatteringSpecies& aoss,
              const AtmPoint& atm_point,
              const Vector& f_grid,
              Index l) {
             return BulkScatteringPropertiesTROSpectral{
                 aoss.get_bulk_scattering_properties_tro_spectral(
                     atm_point, f_grid, l)};
           })
      .def("get_bulk_scattering_properties_tro_gridded",
           [](const ArrayOfScatteringSpecies& aoss,
              const AtmPoint& atm_point,
              const Vector& f_grid,
              std::shared_ptr<scattering::ZenithAngleGrid> za_grid) {
             return BulkScatteringPropertiesTROGridded{
                 aoss.get_bulk_scattering_properties_tro_gridded(
                     atm_point, f_grid, za_grid)};
           })
      .def("get_bulk_scattering_properties_aro_gridded",
           [](const ArrayOfScatteringSpecies& aoss,
              const AtmPoint& atm_point,
              const Vector& f_grid,
              const Vector& za_inc_grid,
              const Vector& delta_aa_grid,
              std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) {
             return BulkScatteringPropertiesAROGridded{
                 aoss.get_bulk_scattering_properties_aro_gridded(atm_point,
                                                                 f_grid,
                                                                 za_inc_grid,
                                                                 delta_aa_grid,
                                                                 za_scat_grid)};
           })
      .def("get_bulk_scattering_properties_aro_spectral",
           [](const ArrayOfScatteringSpecies& aoss,
              const AtmPoint& atm_point,
              const Vector& f_grid,
              const Vector& za_inc_grid,
              Index l,
              Index m) {
             return BulkScatteringPropertiesAROSpectral{
                 aoss.get_bulk_scattering_properties_aro_spectral(
                     atm_point, f_grid, za_inc_grid, l, m)};
           });

  workspace_group_interface(aoss);

  bind_phase_matrix_data_tro_gridded<double>(m, "PhaseMatrixDataTROGridded4");
  bind_phase_matrix_data_tro_spectral<double>(m, "PhaseMatrixDataTROSpectral4");

  bind_absorption_vector_data_tro<double, scattering::Representation::Gridded>(
      m, "AbsorptionVectorDataGriddedTRO4");
  bind_absorption_vector_data_tro<double, scattering::Representation::Spectral>(
      m, "AbsorptionVectorDataSpectralTRO4");
  bind_absorption_vector_data_aro<double, scattering::Representation::Gridded>(
      m, "AbsorptionVectorDataGriddedARO4");
  bind_absorption_vector_data_aro<double, scattering::Representation::Spectral>(
      m, "AbsorptionVectorDataSpectralARO4");

  bind_extinction_matrix_data_tro<double, scattering::Representation::Gridded>(
      m, "ExtinctionMatrixDataGriddedTRO4");
  bind_extinction_matrix_data_tro<double, scattering::Representation::Spectral>(
      m, "ExtinctionMatrixDataSpectralTRO4");
  bind_extinction_matrix_data_aro<double, scattering::Representation::Gridded>(
      m, "ExtinctionMatrixDataGriddedARO4");
  bind_extinction_matrix_data_aro<double, scattering::Representation::Spectral>(
      m, "ExtinctionMatrixDataSpectralARO4");

  py::class_<scattering::ParticleProperties>(m, "ParticleProperties")
      .def(py::init<>())
      .def_rw("name", &scattering::ParticleProperties::name)
      .def_rw("source", &scattering::ParticleProperties::source)
      .def_rw("refractive_index",
              &scattering::ParticleProperties::refractive_index)
      .def_rw("mass", &scattering::ParticleProperties::mass)
      .def_rw("d_veq", &scattering::ParticleProperties::d_veq)
      .def_rw("d_max", &scattering::ParticleProperties::d_max);

  bind_single_scattering_data<double,
                              scattering::Format::TRO,
                              scattering::Representation::Gridded>(
      m, "SingleScatteringDataTROGridded4");
  bind_single_scattering_data<double,
                              scattering::Format::ARO,
                              scattering::Representation::Gridded>(
      m, "SingleScatteringDataAROGridded1");
  bind_single_scattering_data<double,
                              scattering::Format::TRO,
                              scattering::Representation::Spectral>(
      m, "SingleScatteringDataTROSpectral4");
  bind_single_scattering_data<double,
                              scattering::Format::ARO,
                              scattering::Representation::Spectral>(
      m, "SingleScatteringDataAROSpectral1");

  bind_bulk_scattering_properties<scattering::Format::TRO,
                                  scattering::Representation::Gridded>(
      m, "BulkScatteringPropertiesTROGridded4");
  bind_bulk_scattering_properties<scattering::Format::TRO,
                                  scattering::Representation::Spectral>(
      m, "BulkScatteringPropertiesTROSpectral4");

  py::class_<ParticleHabit>(m, "ParticleHabit")
      .def_static(
          "from_legacy_tro", &ParticleHabit::from_legacy_tro, "ssd"_a, "smd"_a);

} catch (std::exception& e) {
  throw std::runtime_error(var_string(
      "DEV ERROR:\nCannot initialize scattering species:\n", e.what()));
};
}  // namespace Python
