#include <core/scattering/particle_habit.h>
#include <core/scattering/single_scattering_data.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <python_interface.h>

#include "hpy_arts.h"
#include "py_macros.h"

NB_MAKE_OPAQUE(scattering::ZenithAngleGrid);

namespace Python {

template <typename Scalar, Index stokes_dim>
void bind_phase_matrix_data_tro_gridded(py::module_& m,
                                        const std::string& class_name) {
  using PMD = scattering::PhaseMatrixData<Scalar,
                                          scattering::Format::TRO,
                                          scattering::Representation::Gridded,
                                          stokes_dim>;

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
      .def("extract_stokes_coeffs",
           &PMD::template extract_stokes_coeffs<
               3>)  // You can specify the new stokes_dim here

      // Bind regrid method
      .def("regrid", &PMD::regrid);
}

template <typename Scalar, Index stokes_dim>
void bind_phase_matrix_data_tro_spectral(py::module_& m,
                                         const std::string& class_name) {
  using PMD = scattering::PhaseMatrixData<Scalar,
                                          scattering::Format::TRO,
                                          scattering::Representation::Spectral,
                                          stokes_dim>;
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
      .def("extract_stokes_coeffs",
           &PMD::template extract_stokes_coeffs<
               3>)  // You can specify the new stokes_dim here
      .def("regrid", &PMD::regrid);
}

template <typename Scalar, scattering::Representation repr, Index stokes_dim>
void bind_absorption_vector_data_tro(py::module_& m, const std::string& name) {
  using AVD = scattering::
      AbsorptionVectorData<Scalar, scattering::Format::TRO, repr, stokes_dim>;
  py::class_<AVD, matpack::matpack_data<Scalar, 3>>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>>(),
           "t_grid"_a,
           "f_grid"_a)
      .def("get_coeff_vector_view", &AVD::get_coeff_vector_view);
}

template <typename Scalar, scattering::Representation repr, Index stokes_dim>
void bind_absorption_vector_data_aro(py::module_& m, const std::string& name) {
  using AVD = scattering::
      AbsorptionVectorData<Scalar, scattering::Format::ARO, repr, stokes_dim>;
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

template <typename Scalar, scattering::Representation repr, Index stokes_dim>
void bind_extinction_matrix_data_tro(py::module_& m, const std::string& name) {
  using EMD = scattering::
      ExtinctionMatrixData<Scalar, scattering::Format::TRO, repr, stokes_dim>;
  py::class_<EMD, matpack::matpack_data<Scalar, 3>>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>>(),
           "t_grid"_a,
           "f_grid"_a)
      .def("get_coeff_vector_view", &EMD::get_coeff_vector_view)
      .def("regrid", &EMD::regrid, "grids"_a, "weights"_a);
}

template <typename Scalar, scattering::Representation repr, Index stokes_dim>
void bind_extinction_matrix_data_aro(py::module_& m, const std::string& name) {
  using EMD = scattering::
      ExtinctionMatrixData<Scalar, scattering::Format::ARO, repr, stokes_dim>;
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
          scattering::Representation repr,
          Index stokes_dim>
void bind_single_scattering_data(py::module_& m, const std::string& name) {
  using SSDClass =
      scattering::SingleScatteringData<Scalar, format, repr, stokes_dim>;

  py::class_<SSDClass>(m, name.c_str())
      .def(py::init<
               scattering::PhaseMatrixData<Scalar, format, repr, stokes_dim>,
               scattering::
                   ExtinctionMatrixData<Scalar, format, repr, stokes_dim>,
               scattering::
                   AbsorptionVectorData<Scalar, format, repr, stokes_dim>,
               scattering::BackscatterMatrixData<Scalar, format, stokes_dim>,
               scattering::
                   ForwardscatterMatrixData<Scalar, format, stokes_dim>>(),
           "phase_matrix"_a,
           "extinction_matrix"_a,
           "absorption_vector"_a,
           "backscatter_matrix"_a,
           "forwardscatter_matrix"_a)
      .def(py::init<
               scattering::ParticleProperties,
               scattering::PhaseMatrixData<Scalar, format, repr, stokes_dim>,
               scattering::
                   ExtinctionMatrixData<Scalar, format, repr, stokes_dim>,
               scattering::
                   AbsorptionVectorData<Scalar, format, repr, stokes_dim>,
               scattering::BackscatterMatrixData<Scalar, format, stokes_dim>,
               scattering::
                   ForwardscatterMatrixData<Scalar, format, stokes_dim>>(),
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

template <scattering::Format format,
          scattering::Representation repr,
          Index stokes_dim>
void bind_bulk_scattering_properties(py::module_& m, const std::string& name) {
  py::class_<scattering::BulkScatteringProperties<format, repr, stokes_dim>>(
      m, name.c_str())
      .def_rw("phase_matrix",
              &scattering::BulkScatteringProperties<format, repr, stokes_dim>::
                  phase_matrix)
      .def_rw("extinction_matrix",
              &scattering::BulkScatteringProperties<format, repr, stokes_dim>::
                  extinction_matrix)
      .def_rw("absorption_vector",
              &scattering::BulkScatteringProperties<format, repr, stokes_dim>::
                  absorption_vector);
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

  using BulkScatteringPropertiesTROSpectral = std::variant<
      scattering::BulkScatteringProperties<scattering::Format::TRO,
                                           scattering::Representation::Spectral,
                                           1>,
      scattering::BulkScatteringProperties<scattering::Format::TRO,
                                           scattering::Representation::Spectral,
                                           2>,
      scattering::BulkScatteringProperties<scattering::Format::TRO,
                                           scattering::Representation::Spectral,
                                           3>,
      scattering::BulkScatteringProperties<scattering::Format::TRO,
                                           scattering::Representation::Spectral,
                                           4>>;
  using BulkScatteringPropertiesTROGridded = std::variant<
      scattering::BulkScatteringProperties<scattering::Format::TRO,
                                           scattering::Representation::Gridded,
                                           1>,
      scattering::BulkScatteringProperties<scattering::Format::TRO,
                                           scattering::Representation::Gridded,
                                           2>,
      scattering::BulkScatteringProperties<scattering::Format::TRO,
                                           scattering::Representation::Gridded,
                                           3>,
      scattering::BulkScatteringProperties<scattering::Format::TRO,
                                           scattering::Representation::Gridded,
                                           4>>;
  using BulkScatteringPropertiesAROSpectral = std::variant<
      scattering::BulkScatteringProperties<scattering::Format::ARO,
                                           scattering::Representation::Spectral,
                                           1>,
      scattering::BulkScatteringProperties<scattering::Format::ARO,
                                           scattering::Representation::Spectral,
                                           2>,
      scattering::BulkScatteringProperties<scattering::Format::ARO,
                                           scattering::Representation::Spectral,
                                           3>,
      scattering::BulkScatteringProperties<scattering::Format::ARO,
                                           scattering::Representation::Spectral,
                                           4>>;
  using BulkScatteringPropertiesAROGridded = std::variant<
      scattering::BulkScatteringProperties<scattering::Format::ARO,
                                           scattering::Representation::Gridded,
                                           1>,
      scattering::BulkScatteringProperties<scattering::Format::ARO,
                                           scattering::Representation::Gridded,
                                           2>,
      scattering::BulkScatteringProperties<scattering::Format::ARO,
                                           scattering::Representation::Gridded,
                                           3>,
      scattering::BulkScatteringProperties<scattering::Format::ARO,
                                           scattering::Representation::Gridded,
                                           4>>;

  py::class_<MGDSingleMoment>(m, "MGDSingleMoment");
  py::class_<ScatteringHabit>(m, "ScatteringHabit");
  py::class_<HenyeyGreensteinScatterer>(m, "HenyeyGreensteinScatterer")
      .def(py::init<>())
      .def(py::init<ScatteringSpeciesProperty,
                    ScatteringSpeciesProperty,
                    Numeric>())
      .def("get_bulk_scattering_properties_tro_spectral",
           [](const HenyeyGreensteinScatterer& hg,
              const AtmPoint& atm_point,
              const Vector& f_grid,
              Index l,
              const Index stokes_dim) {
             if (stokes_dim == 1)
               return BulkScatteringPropertiesTROSpectral{
                   hg.get_bulk_scattering_properties_tro_spectral<1>(
                       atm_point, f_grid, l)};
             if (stokes_dim == 2)
               return BulkScatteringPropertiesTROSpectral{
                   hg.get_bulk_scattering_properties_tro_spectral<2>(
                       atm_point, f_grid, l)};
             if (stokes_dim == 3)
               return BulkScatteringPropertiesTROSpectral{
                   hg.get_bulk_scattering_properties_tro_spectral<3>(
                       atm_point, f_grid, l)};
             if (stokes_dim == 4)
               return BulkScatteringPropertiesTROSpectral{
                   hg.get_bulk_scattering_properties_tro_spectral<4>(
                       atm_point, f_grid, l)};
             throw std::runtime_error(
                 "Stokes dim must be one of 1, 2, 3, or 4.");
           })
      .def("get_bulk_scattering_properties_tro_gridded",
           [](const HenyeyGreensteinScatterer& hg,
              const AtmPoint& atm_point,
              const Vector& f_grid,
              std::shared_ptr<scattering::ZenithAngleGrid> za_grid,
              const Index stokes_dim) {
             if (stokes_dim == 1)
               return BulkScatteringPropertiesTROGridded{
                   hg.get_bulk_scattering_properties_tro_gridded<1>(
                       atm_point, f_grid, za_grid)};
             if (stokes_dim == 2)
               return BulkScatteringPropertiesTROGridded{
                   hg.get_bulk_scattering_properties_tro_gridded<2>(
                       atm_point, f_grid, za_grid)};
             if (stokes_dim == 3)
               return BulkScatteringPropertiesTROGridded{
                   hg.get_bulk_scattering_properties_tro_gridded<3>(
                       atm_point, f_grid, za_grid)};
             if (stokes_dim == 4)
               return BulkScatteringPropertiesTROGridded{
                   hg.get_bulk_scattering_properties_tro_gridded<4>(
                       atm_point, f_grid, za_grid)};
             throw std::runtime_error(
                 "Stokes dim must be one of 1, 2, 3, or 4.");
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
      .def("add", &ArrayOfScatteringSpecies::add)
      .def("get_bulk_scattering_properties_tro_spectral",
           [](const ArrayOfScatteringSpecies& aoss,
              const AtmPoint& atm_point,
              const Vector& f_grid,
              Index l,
              const Index stokes_dim) {
             if (stokes_dim == 1)
               return BulkScatteringPropertiesTROSpectral{
                   aoss.get_bulk_scattering_properties_tro_spectral<1>(
                       atm_point, f_grid, l)};
             if (stokes_dim == 2)
               return BulkScatteringPropertiesTROSpectral{
                   aoss.get_bulk_scattering_properties_tro_spectral<2>(
                       atm_point, f_grid, l)};
             if (stokes_dim == 3)
               return BulkScatteringPropertiesTROSpectral{
                   aoss.get_bulk_scattering_properties_tro_spectral<3>(
                       atm_point, f_grid, l)};
             if (stokes_dim == 4)
               return BulkScatteringPropertiesTROSpectral{
                   aoss.get_bulk_scattering_properties_tro_spectral<4>(
                       atm_point, f_grid, l)};
             throw std::runtime_error(
                 "Stokes dim must be one of 1, 2, 3, or 4.");
           })
      .def("get_bulk_scattering_properties_tro_gridded",
           [](const ArrayOfScatteringSpecies& aoss,
              const AtmPoint& atm_point,
              const Vector& f_grid,
              std::shared_ptr<scattering::ZenithAngleGrid> za_grid,
              const Index stokes_dim) {
             if (stokes_dim == 1)
               return BulkScatteringPropertiesTROGridded{
                   aoss.get_bulk_scattering_properties_tro_gridded<1>(
                       atm_point, f_grid, za_grid)};
             if (stokes_dim == 2)
               return BulkScatteringPropertiesTROGridded{
                   aoss.get_bulk_scattering_properties_tro_gridded<2>(
                       atm_point, f_grid, za_grid)};
             if (stokes_dim == 3)
               return BulkScatteringPropertiesTROGridded{
                   aoss.get_bulk_scattering_properties_tro_gridded<3>(
                       atm_point, f_grid, za_grid)};
             if (stokes_dim == 4)
               return BulkScatteringPropertiesTROGridded{
                   aoss.get_bulk_scattering_properties_tro_gridded<4>(
                       atm_point, f_grid, za_grid)};
             throw std::runtime_error(
                 "Stokes dim must be one of 1, 2, 3, or 4.");
           })
      .def("get_bulk_scattering_properties_aro_gridded",
           [](const ArrayOfScatteringSpecies& aoss,
              const AtmPoint& atm_point,
              const Vector& f_grid,
              const Vector& za_inc_grid,
              const Vector& delta_aa_grid,
              std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid,
              const Index stokes_dim) {
             if (stokes_dim == 1)
               return BulkScatteringPropertiesAROGridded{
                   aoss.get_bulk_scattering_properties_aro_gridded<1>(
                       atm_point,
                       f_grid,
                       za_inc_grid,
                       delta_aa_grid,
                       za_scat_grid)};
             if (stokes_dim == 2)
               return BulkScatteringPropertiesAROGridded{
                   aoss.get_bulk_scattering_properties_aro_gridded<2>(
                       atm_point,
                       f_grid,
                       za_inc_grid,
                       delta_aa_grid,
                       za_scat_grid)};
             if (stokes_dim == 3)
               return BulkScatteringPropertiesAROGridded{
                   aoss.get_bulk_scattering_properties_aro_gridded<3>(
                       atm_point,
                       f_grid,
                       za_inc_grid,
                       delta_aa_grid,
                       za_scat_grid)};
             if (stokes_dim == 4)
               return BulkScatteringPropertiesAROGridded{
                   aoss.get_bulk_scattering_properties_aro_gridded<4>(
                       atm_point,
                       f_grid,
                       za_inc_grid,
                       delta_aa_grid,
                       za_scat_grid)};
             throw std::runtime_error(
                 "Stokes dim must be one of 1, 2, 3, or 4.");
           })
      .def("get_bulk_scattering_properties_aro_spectral",
           [](const ArrayOfScatteringSpecies& aoss,
              const AtmPoint& atm_point,
              const Vector& f_grid,
              const Vector& za_inc_grid,
              Index l,
              Index m,
              const Index stokes_dim) {
             if (stokes_dim == 1)
               return BulkScatteringPropertiesAROSpectral{
                   aoss.get_bulk_scattering_properties_aro_spectral<1>(
                       atm_point, f_grid, za_inc_grid, l, m)};
             if (stokes_dim == 2)
               return BulkScatteringPropertiesAROSpectral{
                   aoss.get_bulk_scattering_properties_aro_spectral<2>(
                       atm_point, f_grid, za_inc_grid, l, m)};
             if (stokes_dim == 3)
               return BulkScatteringPropertiesAROSpectral{
                   aoss.get_bulk_scattering_properties_aro_spectral<3>(
                       atm_point, f_grid, za_inc_grid, l, m)};
             if (stokes_dim == 4)
               return BulkScatteringPropertiesAROSpectral{
                   aoss.get_bulk_scattering_properties_aro_spectral<4>(
                       atm_point, f_grid, za_inc_grid, l, m)};
             std::runtime_error("Stokes dim must be one of 1, 2, 3, or 4.");
             return BulkScatteringPropertiesAROSpectral{
                 aoss.get_bulk_scattering_properties_aro_spectral<1>(
                     atm_point, f_grid, za_inc_grid, l, m)};
           });

  workspace_group_interface(aoss);

  bind_phase_matrix_data_tro_gridded<double, 1>(m,
                                                "PhaseMatrixDataTROGridded1");
  bind_phase_matrix_data_tro_gridded<double, 2>(m,
                                                "PhaseMatrixDataTROGridded2");
  bind_phase_matrix_data_tro_gridded<double, 3>(m,
                                                "PhaseMatrixDataTROGridded3");
  bind_phase_matrix_data_tro_gridded<double, 4>(m,
                                                "PhaseMatrixDataTROGridded4");
  bind_phase_matrix_data_tro_spectral<double, 1>(m,
                                                 "PhaseMatrixDataTROSpectral1");
  bind_phase_matrix_data_tro_spectral<double, 2>(m,
                                                 "PhaseMatrixDataTROSpectral2");
  bind_phase_matrix_data_tro_spectral<double, 3>(m,
                                                 "PhaseMatrixDataTROSpectral3");
  bind_phase_matrix_data_tro_spectral<double, 4>(m,
                                                 "PhaseMatrixDataTROSpectral4");

  bind_absorption_vector_data_tro<double,
                                  scattering::Representation::Gridded,
                                  1>(m, "AbsorptionVectorDataGriddedTRO1");
  bind_absorption_vector_data_tro<double,
                                  scattering::Representation::Gridded,
                                  2>(m, "AbsorptionVectorDataGriddedTRO2");
  bind_absorption_vector_data_tro<double,
                                  scattering::Representation::Gridded,
                                  3>(m, "AbsorptionVectorDataGriddedTRO3");
  bind_absorption_vector_data_tro<double,
                                  scattering::Representation::Gridded,
                                  4>(m, "AbsorptionVectorDataGriddedTRO4");
  bind_absorption_vector_data_tro<double,
                                  scattering::Representation::Spectral,
                                  1>(m, "AbsorptionVectorDataSpectralTRO1");
  bind_absorption_vector_data_tro<double,
                                  scattering::Representation::Spectral,
                                  2>(m, "AbsorptionVectorDataSpectralTRO2");
  bind_absorption_vector_data_tro<double,
                                  scattering::Representation::Spectral,
                                  3>(m, "AbsorptionVectorDataSpectralTRO3");
  bind_absorption_vector_data_tro<double,
                                  scattering::Representation::Spectral,
                                  4>(m, "AbsorptionVectorDataSpectralTRO4");
  bind_absorption_vector_data_aro<double,
                                  scattering::Representation::Gridded,
                                  1>(m, "AbsorptionVectorDataGriddedARO1");
  bind_absorption_vector_data_aro<double,
                                  scattering::Representation::Gridded,
                                  2>(m, "AbsorptionVectorDataGriddedARO2");
  bind_absorption_vector_data_aro<double,
                                  scattering::Representation::Gridded,
                                  3>(m, "AbsorptionVectorDataGriddedARO3");
  bind_absorption_vector_data_aro<double,
                                  scattering::Representation::Gridded,
                                  4>(m, "AbsorptionVectorDataGriddedARO4");
  bind_absorption_vector_data_aro<double,
                                  scattering::Representation::Spectral,
                                  1>(m, "AbsorptionVectorDataSpectralARO1");
  bind_absorption_vector_data_aro<double,
                                  scattering::Representation::Spectral,
                                  2>(m, "AbsorptionVectorDataSpectralARO2");
  bind_absorption_vector_data_aro<double,
                                  scattering::Representation::Spectral,
                                  3>(m, "AbsorptionVectorDataSpectralARO3");
  bind_absorption_vector_data_aro<double,
                                  scattering::Representation::Spectral,
                                  4>(m, "AbsorptionVectorDataSpectralARO4");

  bind_extinction_matrix_data_tro<double,
                                  scattering::Representation::Gridded,
                                  1>(m, "ExtinctionMatrixDataGriddedTRO1");
  bind_extinction_matrix_data_tro<double,
                                  scattering::Representation::Gridded,
                                  2>(m, "ExtinctionMatrixDataGriddedTRO2");
  bind_extinction_matrix_data_tro<double,
                                  scattering::Representation::Gridded,
                                  3>(m, "ExtinctionMatrixDataGriddedTRO3");
  bind_extinction_matrix_data_tro<double,
                                  scattering::Representation::Gridded,
                                  4>(m, "ExtinctionMatrixDataGriddedTRO4");
  bind_extinction_matrix_data_tro<double,
                                  scattering::Representation::Spectral,
                                  1>(m, "ExtinctionMatrixDataSpectralTRO1");
  bind_extinction_matrix_data_tro<double,
                                  scattering::Representation::Spectral,
                                  2>(m, "ExtinctionMatrixDataSpectralTRO2");
  bind_extinction_matrix_data_tro<double,
                                  scattering::Representation::Spectral,
                                  3>(m, "ExtinctionMatrixDataSpectralTRO3");
  bind_extinction_matrix_data_tro<double,
                                  scattering::Representation::Spectral,
                                  4>(m, "ExtinctionMatrixDataSpectralTRO4");
  bind_extinction_matrix_data_aro<double,
                                  scattering::Representation::Gridded,
                                  1>(m, "ExtinctionMatrixDataGriddedARO1");
  bind_extinction_matrix_data_aro<double,
                                  scattering::Representation::Gridded,
                                  2>(m, "ExtinctionMatrixDataGriddedARO2");
  bind_extinction_matrix_data_aro<double,
                                  scattering::Representation::Gridded,
                                  3>(m, "ExtinctionMatrixDataGriddedARO3");
  bind_extinction_matrix_data_aro<double,
                                  scattering::Representation::Gridded,
                                  4>(m, "ExtinctionMatrixDataGriddedARO4");
  bind_extinction_matrix_data_aro<double,
                                  scattering::Representation::Spectral,
                                  1>(m, "ExtinctionMatrixDataSpectralARO1");
  bind_extinction_matrix_data_aro<double,
                                  scattering::Representation::Spectral,
                                  2>(m, "ExtinctionMatrixDataSpectralARO2");
  bind_extinction_matrix_data_aro<double,
                                  scattering::Representation::Spectral,
                                  3>(m, "ExtinctionMatrixDataSpectralARO3");
  bind_extinction_matrix_data_aro<double,
                                  scattering::Representation::Spectral,
                                  4>(m, "ExtinctionMatrixDataSpectralARO4");

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
                              scattering::Representation::Gridded,
                              1>(m, "SingleScatteringDataTROGridded1");
  bind_single_scattering_data<double,
                              scattering::Format::TRO,
                              scattering::Representation::Gridded,
                              2>(m, "SingleScatteringDataTROGridded2");
  bind_single_scattering_data<double,
                              scattering::Format::TRO,
                              scattering::Representation::Gridded,
                              3>(m, "SingleScatteringDataTROGridded3");
  bind_single_scattering_data<double,
                              scattering::Format::TRO,
                              scattering::Representation::Gridded,
                              4>(m, "SingleScatteringDataTROGridded4");
  bind_single_scattering_data<double,
                              scattering::Format::ARO,
                              scattering::Representation::Gridded,
                              1>(m, "SingleScatteringDataAROGridded1");
  bind_single_scattering_data<double,
                              scattering::Format::ARO,
                              scattering::Representation::Gridded,
                              2>(m, "SingleScatteringDataAROGridded1");
  bind_single_scattering_data<double,
                              scattering::Format::ARO,
                              scattering::Representation::Gridded,
                              3>(m, "SingleScatteringDataAROGridded1");
  bind_single_scattering_data<double,
                              scattering::Format::ARO,
                              scattering::Representation::Gridded,
                              4>(m, "SingleScatteringDataAROGridded1");
  bind_single_scattering_data<double,
                              scattering::Format::TRO,
                              scattering::Representation::Spectral,
                              1>(m, "SingleScatteringDataTROSpectral1");
  bind_single_scattering_data<double,
                              scattering::Format::TRO,
                              scattering::Representation::Spectral,
                              2>(m, "SingleScatteringDataTROSpectral2");
  bind_single_scattering_data<double,
                              scattering::Format::TRO,
                              scattering::Representation::Spectral,
                              3>(m, "SingleScatteringDataTROSpectral3");
  bind_single_scattering_data<double,
                              scattering::Format::TRO,
                              scattering::Representation::Spectral,
                              4>(m, "SingleScatteringDataTROSpectral4");
  bind_single_scattering_data<double,
                              scattering::Format::ARO,
                              scattering::Representation::Spectral,
                              1>(m, "SingleScatteringDataAROSpectral1");
  bind_single_scattering_data<double,
                              scattering::Format::ARO,
                              scattering::Representation::Spectral,
                              2>(m, "SingleScatteringDataAROSpectral1");
  bind_single_scattering_data<double,
                              scattering::Format::ARO,
                              scattering::Representation::Spectral,
                              3>(m, "SingleScatteringDataAROSpectral1");
  bind_single_scattering_data<double,
                              scattering::Format::ARO,
                              scattering::Representation::Spectral,
                              4>(m, "SingleScatteringDataAROSpectral1");

  bind_bulk_scattering_properties<scattering::Format::TRO,
                                  scattering::Representation::Gridded,
                                  1>(m, "BulkScatteringPropertiesAROGridded1");
  bind_bulk_scattering_properties<scattering::Format::TRO,
                                  scattering::Representation::Spectral,
                                  1>(m, "BulkScatteringPropertiesAROSpectral1");
  bind_bulk_scattering_properties<scattering::Format::TRO,
                                  scattering::Representation::Gridded,
                                  2>(m, "BulkScatteringPropertiesAROGridded2");
  bind_bulk_scattering_properties<scattering::Format::TRO,
                                  scattering::Representation::Spectral,
                                  2>(m, "BulkScatteringPropertiesAROSpectral2");
  bind_bulk_scattering_properties<scattering::Format::TRO,
                                  scattering::Representation::Gridded,
                                  3>(m, "BulkScatteringPropertiesAROGridded3");
  bind_bulk_scattering_properties<scattering::Format::TRO,
                                  scattering::Representation::Spectral,
                                  3>(m, "BulkScatteringPropertiesAROSpectral3");
  bind_bulk_scattering_properties<scattering::Format::TRO,
                                  scattering::Representation::Gridded,
                                  4>(m, "BulkScatteringPropertiesAROGridded4");
  bind_bulk_scattering_properties<scattering::Format::TRO,
                                  scattering::Representation::Spectral,
                                  4>(m, "BulkScatteringPropertiesAROSpectral4");

  py::class_<ParticleHabit>(m, "ParticleHabit")
      .def_static(
          "from_legacy_tro", &ParticleHabit::from_legacy_tro, "ssd"_a, "smd"_a);

} catch (std::exception& e) {
  throw std::runtime_error(var_string(
      "DEV ERROR:\nCannot initialize scattering species:\n", e.what()));
};
}  // namespace Python
