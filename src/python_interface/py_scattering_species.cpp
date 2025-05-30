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
#include <stdexcept>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "py_macros.h"


namespace Python {

template <typename Scalar>
[[nodiscard]]
auto bind_phase_matrix_data_tro_gridded(py::module_& m,
                                        const std::string& class_name) {
  using PMD = scattering::PhaseMatrixData<Scalar,
                                          scattering::Format::TRO,
                                          scattering::Representation::Gridded>;

  py::class_<PMD, matpack::data_t<Scalar, 4>> s(m, class_name.c_str());
  s.def(
       "__init__",
       [](PMD* self,
          std::shared_ptr<const Vector> t_grid,
          std::shared_ptr<const Vector> f_grid,
          scattering::ZenithAngleGrid za_grid) {
         new (self) PMD{
             t_grid,
             f_grid,
             std::make_shared<const scattering::ZenithAngleGrid>(std::move(za_grid))};
       },
       py::arg("t_grid"),
       py::arg("f_grid"),
       py::arg("za_scat_grid"))
      // Bind methods, such as `extract_backscatter_matrix` and `extract_forwardscatter_matrix`
      .def("extract_backscatter_matrix",
           &PMD::extract_backscatter_matrix,
           "Extract backscatter matrix")
      .def("extract_forwardscatter_matrix",
           &PMD::extract_forwardscatter_matrix,
           "Extract forward scatter matrix")

      .def("get_t_grid", &PMD::get_t_grid, "Get temperature grid")
      .def("get_f_grid", &PMD::get_f_grid, "Get frequency grid")
      .def(
          "get_za_scat_grid",
          [](const PMD& pmd) {
            if (pmd.get_za_scat_grid()) return *(pmd.get_za_scat_grid());
            throw std::runtime_error("PMD Zenith angle grid not initialized");
          },
          "Get scattering zenith angle grid")

      .def(
          "to_spectral",
          [](const PMD& obj) {return obj.to_spectral(); },
          "Convert to spectral")

      // Bind other member functions
      .def("integrate_phase_matrix",
           &PMD::integrate_phase_matrix,
           "Integrate phase matrix")

      // Bind the extraction of stokes coefficients
      .def("extract_stokes_coeffs",
           &PMD::extract_stokes_coeffs,
           "Extract stokes coefficients from phase matrix");

      // Bind regrid method
      //.def("regrid", &PMD::regrid, "Regrid phase matrix");
  return s;
}

template <typename Scalar>
[[nodiscard]]
auto bind_phase_matrix_data_tro_spectral(py::module_& m,
                                         const std::string& class_name) {
  using PMD = scattering::PhaseMatrixData<Scalar,
                                          scattering::Format::TRO,
                                          scattering::Representation::Spectral>;
  py::class_<PMD, matpack::data_t<std::complex<Scalar>, 4>> s(
      m, class_name.c_str());
  s.def(py::init<>())
      .def(py::init<std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>,
                    std::shared_ptr<scattering::SHT>>(),
           py::arg("t_grid"),
           py::arg("f_grid"),
           py::arg("sht"))

      .def("extract_backscatter_matrix",
           &PMD::extract_backscatter_matrix,
           "Extract backscatter matrix")
      .def("extract_forwardscatter_matrix",
           &PMD::extract_forwardscatter_matrix,
           "Extract forward scatter matrix")

      .def("get_t_grid", &PMD::get_t_grid, "Get temperature grid")
      .def("get_f_grid", &PMD::get_f_grid, "Get frequency grid")

      .def(
          "to_gridded",
          [](const PMD& obj) { return obj.to_gridded(); },
          "Convert to gridded")
      .def("integrate_phase_matrix",
           &PMD::integrate_phase_matrix,
           "Integrate phase matrix")
      .def("extract_stokes_coeffs",
           &PMD::extract_stokes_coeffs,
           "Extract stokes coefficients from phase matrix");
      //.def("regrid", &PMD::regrid, "Regrid phase matrix");
  return s;
}

template <typename Scalar, scattering::Representation repr>
[[nodiscard]]
auto bind_absorption_vector_data_tro(py::module_& m, const std::string& name) {
  using AVD =
      scattering::AbsorptionVectorData<Scalar, scattering::Format::TRO, repr>;
  py::class_<AVD, matpack::data_t<Scalar, 3>> s(m, name.c_str());
  s.def(py::init<>())
      .def(py::init<std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>>(),
           "t_grid"_a,
           "f_grid"_a)
      .def("get_coeff_vector_view",
           &AVD::get_coeff_vector_view,
           "Get coefficient vector view");
  return s;
}

template <typename Scalar, scattering::Representation repr>
[[nodiscard]]
auto bind_absorption_vector_data_aro(py::module_& m, const std::string& name) {
  using AVD =
      scattering::AbsorptionVectorData<Scalar, scattering::Format::ARO, repr>;
  py::class_<AVD, matpack::data_t<Scalar, 4>> s(m, name.c_str());
  s.def(py::init<>())
      .def(py::init<std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>>(),
           "t_grid"_a,
           "f_grid"_a,
           "za_inc_grid"_a)
      .def("get_coeff_vector_view",
           &AVD::get_coeff_vector_view,
           "Get coefficient vector view");
      //.def("regrid",
      //     &AVD::regrid,
      //     "grids"_a,
      //     "weights"_a,
      //     "Regrid absorption vector");
  return s;
}

template <typename Scalar, scattering::Representation repr>
[[nodiscard]]
auto bind_extinction_matrix_data_tro(py::module_& m, const std::string& name) {
  using EMD =
      scattering::ExtinctionMatrixData<Scalar, scattering::Format::TRO, repr>;
  py::class_<EMD, matpack::data_t<Scalar, 3>> s(m, name.c_str());
  s.def(py::init<>())
      .def(py::init<std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>>(),
           "t_grid"_a,
           "f_grid"_a)
      .def("get_coeff_vector_view",
           &EMD::get_coeff_vector_view,
           "Get coefficient vector view");
      //.def("regrid",
      //     &EMD::regrid,
      //     "grids"_a,
      //     "weights"_a,
      //     "Regrid extinction matrix");
  return s;
}

template <typename Scalar, scattering::Representation repr>
[[nodiscard]]
auto bind_extinction_matrix_data_aro(py::module_& m, const std::string& name) {
  using EMD =
      scattering::ExtinctionMatrixData<Scalar, scattering::Format::ARO, repr>;
  py::class_<EMD, matpack::data_t<Scalar, 4>> s(m, name.c_str());
  s.def(py::init<>())
      .def(py::init<std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>,
                    std::shared_ptr<const Vector>>(),
           "t_grid"_a,
           "f_grid"_a,
           "za_inc_grid"_a)
      .def("get_coeff_vector_view",
           &EMD::get_coeff_vector_view,
           "Get coefficient vector view");
      //.def("regrid",
      //     &EMD::regrid,
      //     "grids"_a,
      //     "weights"_a,
      //     "Regrid extinction matrix");
  return s;
}

template <typename Scalar,
          scattering::Format format,
          scattering::Representation repr>
[[nodiscard]]
auto bind_single_scattering_data(py::module_& m, const std::string& name) {
  using SSDClass = scattering::SingleScatteringData<Scalar, format, repr>;

  py::class_<SSDClass> s(m, name.c_str());
  s.def(py::init<scattering::PhaseMatrixData<Scalar, format, repr>,
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
      .def_rw("properties", &SSDClass::properties, "Particle properties")
      .def_rw("phase_matrix", &SSDClass::phase_matrix, "Phase matrix")
      .def_rw("extinction_matrix",
              &SSDClass::extinction_matrix,
              "Extinction matrix")
      .def_rw("absorption_vector",
              &SSDClass::absorption_vector,
              "Absorption vector")
      .def_rw("backscatter_matrix",
              &SSDClass::backscatter_matrix,
              "Back scatter matrix")
      .def_rw("forwardscatter_matrix",
              &SSDClass::forwardscatter_matrix,
              "Forward scatter matrix")
      .def_static("from_legacy_tro",
                  &SSDClass::from_legacy_tro,
                  "ssd"_a,
                  "smd"_a,
                  "Create from legacy TRO")
      .def("__repr__", [](const SSDClass& ssd) {
        std::ostringstream oss;
        oss << ssd;
        return oss.str();
      });
  return s;
}

template <scattering::Format format, scattering::Representation repr>
[[nodiscard]]
auto bind_bulk_scattering_properties(py::module_& m, const std::string& name) {
  py::class_<scattering::BulkScatteringProperties<format, repr>> s(
      m, name.c_str());
  s.def_rw("phase_matrix",
           &scattering::BulkScatteringProperties<format, repr>::phase_matrix,
           "Phase matrix")
      .def_rw("extinction_matrix",
              &scattering::BulkScatteringProperties<format,
                                                    repr>::extinction_matrix,
              "Extinction matrix")
      .def_rw("absorption_vector",
              &scattering::BulkScatteringProperties<format,
                                                    repr>::absorption_vector,
              "Absorption vector");
  return s;
}

void py_scattering_species(py::module_& m) try {
  py::class_<ScatteringTroSpectralVector> stsv(m,
                                               "ScatteringTroSpectralVector");
  stsv.def_rw("phase_matrix",
              &ScatteringTroSpectralVector::phase_matrix,
              "Phase matrix")
      .def_rw("extinction_matrix",
              &ScatteringTroSpectralVector::extinction_matrix,
              "Extinction matrix")
      .def_rw("absorption_vector",
              &ScatteringTroSpectralVector::absorption_vector,
              "Absorption vector");
  stsv.def_static(
      "to_gridded",
      [](const SpecmatMatrix& phase_matrix, const std::shared_ptr<Vector>& f) {
        if (f) {
          return ScatteringTroSpectralVector::to_general(phase_matrix, f)
              .to_gridded();
        }

        const Index nf = phase_matrix.nrows();
        const auto fs  = std::make_shared<Vector>(
            nlinspace(1, static_cast<Numeric>(nf), nf));
        return ScatteringTroSpectralVector::to_general(phase_matrix, fs)
            .to_gridded();
      },
      "spectral_pm"_a,
      "f"_a = std::shared_ptr<Vector>{nullptr},
      "Convert spectral to gridded");
  str_interface(stsv);
  stsv.doc() = "Scattering TRO spectral vector";

  py::class_<ScatteringGeneralSpectralTRO> sgstro(
      m, "ScatteringGeneralSpectralTRO");
  sgstro.def(py::init<>());
  sgstro.def(py::init_implicit<ScatteringGeneralSpectralTROFunc>());
  sgstro.def(py::init_implicit<ScatteringGeneralSpectralTROFunc::func_t>());
  sgstro.def_rw("f", &ScatteringGeneralSpectralTRO::f, "Frequency grid");
  sgstro.doc() = "Scattering general spectral TRO";
  str_interface(stsv);

  //
  // ScatSpeciesProperty
  //

  py::class_<ScatteringSpeciesProperty> ssp(m, "ScatteringSpeciesProperty");
  workspace_group_interface(ssp);
  ssp.def(py::init<std::string, ParticulateProperty>(), "Constructor")
      .def_rw("species_name",
              &ScatteringSpeciesProperty::species_name,
              "Species name")
      .def_rw("pproperty",
              &ScatteringSpeciesProperty::pproperty,
              "Particulate property");


//   using BulkScatteringPropertiesTROSpectral =
//       std::variant<scattering::BulkScatteringProperties<
//           scattering::Format::TRO,
//           scattering::Representation::Spectral>>;
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

  //
  // Modified gamma PSD
  //

  py::class_<HenyeyGreensteinScatterer>(m, "HenyeyGreensteinScatterer")
      .def(py::init<ExtSSACallback, Numeric>(), "func"_a, "g"_a)
      .def(py::init<>())
      .def(py::init<ScatteringSpeciesProperty,
                    ScatteringSpeciesProperty,
                    Numeric>())
      .def(
          "get_bulk_scattering_properties_tro_spectral",
          [](const HenyeyGreensteinScatterer& hg,
             const AtmPoint& atm_point,
             const Vector& f_grid,
             Index l) {
            return hg.get_bulk_scattering_properties_tro_spectral(
                atm_point, f_grid, l);
          },
          "atm_point"_a,
          "f_grid"_a,
          "l"_a,
          "Get bulk scattering properties")
      .def(
          "get_bulk_scattering_properties_tro_gridded",
          [](const HenyeyGreensteinScatterer& hg,
             const AtmPoint& atm_point,
             const Vector& f_grid,
             scattering::ZenithAngleGrid za_grid) {
            return BulkScatteringPropertiesTROGridded{
                hg.get_bulk_scattering_properties_tro_gridded(
                    atm_point, f_grid, std::make_shared<scattering::ZenithAngleGrid>(std::move(za_grid)))};
          },
          "atm_point"_a,
          "f_grid"_a,
          "za_grid"_a,
          "Get bulk scattering properties")
      .doc() = "Henyey-Greenstein scatterer";

  py::class_<scattering::IrregularZenithAngleGrid> irr_grid(m,
                                                   "IrregularZenithAngleGrid");
      irr_grid.def(py::init<Vector>())
      .def_rw("value", &scattering::IrregularZenithAngleGrid::angles, "Zenith angle grid")
      .doc() = "Irregular zenith angle grid";
      common_ndarray(irr_grid);

  py::class_<scattering::GaussLegendreGrid> gauss_grid(m, "GaussLegendreGrid");
      gauss_grid.def(py::init<Index>())
      .def_rw("value", &scattering::GaussLegendreGrid::angles, "Zenith angle grid for Legendre calculations")
      .doc() = "Gaussian Legendre grid";
      common_ndarray(gauss_grid);

  py::class_<scattering::DoubleGaussGrid> double_gauss_grid(m, "DoubleGaussGrid");
      double_gauss_grid.def(py::init<Index>())
      .def_rw("value", &scattering::DoubleGaussGrid::angles, "Zenith angle grid for Double Gauss calculations")
      .doc() = "Double Gaussian grid";
      common_ndarray(double_gauss_grid);

  py::class_<scattering::LobattoGrid> lobatto_grid(m, "LobattoGrid");
      lobatto_grid.def(py::init<Index>())
      .def_rw("value", &scattering::LobattoGrid::angles, "Zenith angle grid for Lobatto calculations")
      .doc() = "Lobatto grid";
      common_ndarray(lobatto_grid);

  py::class_<scattering::FejerGrid> fejer_grid(m, "FejerGrid");
      fejer_grid.def(py::init<Index>())
      .def_rw("value", &scattering::FejerGrid::angles, "Zenith angle grid for Fejer calculations")
      .doc() = "Fejer grid";
      common_ndarray(fejer_grid);

  py::class_<ArrayOfScatteringSpecies> aoss(m, "ArrayOfScatteringSpecies");
  aoss.def(py::init<>())
      .def(py::init_implicit<std::vector<ScatteringSpecies>>())
      .def("add", &ArrayOfScatteringSpecies::add, "Add scattering species")
      .def(
          "get_bulk_scattering_properties_tro_spectral",
          [](const ArrayOfScatteringSpecies& aoss,
             const AtmPoint& atm_point,
             const Vector& f_grid,
             Index l) {
            return aoss.get_bulk_scattering_properties_tro_spectral(
                atm_point, f_grid, l);
          },
          "atm_point"_a,
          "f_grid"_a,
          "l"_a,
          "Get bulk scattering properties")
      .def(
          "get_bulk_scattering_properties_tro_gridded",
          [](const ArrayOfScatteringSpecies& aoss,
             const AtmPoint& atm_point,
             const Vector& f_grid,
             scattering::ZenithAngleGrid za_grid) {
            return BulkScatteringPropertiesTROGridded{
                aoss.get_bulk_scattering_properties_tro_gridded(
                    atm_point, f_grid, std::make_shared<scattering::ZenithAngleGrid>(std::move(za_grid)))};
          },
          "atm_point"_a,
          "f_grid"_a,
          "za_grid"_a,
          "Get bulk scattering properties")
      .def(
          "get_bulk_scattering_properties_aro_gridded",
          [](const ArrayOfScatteringSpecies& aoss,
             const AtmPoint& atm_point,
             const Vector& f_grid,
             const Vector& za_inc_grid,
             const Vector& delta_aa_grid,
             scattering::ZenithAngleGrid za_scat_grid) {
            return BulkScatteringPropertiesAROGridded{
                aoss.get_bulk_scattering_properties_aro_gridded(atm_point,
                                                                f_grid,
                                                                za_inc_grid,
                                                                delta_aa_grid,
                                                                std::make_shared<scattering::ZenithAngleGrid>(std::move(za_scat_grid)))};
          },
          "atm_point"_a,
          "f_grid"_a,
          "za_inc_grid"_a,
          "delta_aa_grid"_a,
          "za_scat_grid"_a,
          "Get bulk scattering properties")
      .def(
          "get_bulk_scattering_properties_aro_spectral",
          [](const ArrayOfScatteringSpecies& aoss,
             const AtmPoint& atm_point,
             const Vector& f_grid,
             const Vector& za_inc_grid,
             Index l,
             Index m) {
            return BulkScatteringPropertiesAROSpectral{
                aoss.get_bulk_scattering_properties_aro_spectral(
                    atm_point, f_grid, za_inc_grid, l, m)};
          },
          "atm_point"_a,
          "f_grid"_a,
          "za_inc_grid"_a,
          "l"_a,
          "m"_a,
          "Get bulk scattering properties");

  workspace_group_interface(aoss);

  bind_phase_matrix_data_tro_gridded<double>(m, "PhaseMatrixDataTROGridded4")
      .doc() = "Phase matrix data";
  bind_phase_matrix_data_tro_spectral<double>(m, "PhaseMatrixDataTROSpectral4")
      .doc() = "Phase matrix data";

  bind_absorption_vector_data_tro<double, scattering::Representation::Gridded>(
      m, "AbsorptionVectorDataGriddedTRO4")
      .doc() = "Absorption vector data";
  bind_absorption_vector_data_tro<double, scattering::Representation::Spectral>(
      m, "AbsorptionVectorDataSpectralTRO4")
      .doc() = "Absorption vector data";
  bind_absorption_vector_data_aro<double, scattering::Representation::Gridded>(
      m, "AbsorptionVectorDataGriddedARO4")
      .doc() = "Absorption vector data";
  bind_absorption_vector_data_aro<double, scattering::Representation::Spectral>(
      m, "AbsorptionVectorDataSpectralARO4")
      .doc() = "Absorption vector data";

  bind_extinction_matrix_data_tro<double, scattering::Representation::Gridded>(
      m, "ExtinctionMatrixDataGriddedTRO4")
      .doc() = "Extinction matrix data";
  bind_extinction_matrix_data_tro<double, scattering::Representation::Spectral>(
      m, "ExtinctionMatrixDataSpectralTRO4")
      .doc() = "Extinction matrix data";
  bind_extinction_matrix_data_aro<double, scattering::Representation::Gridded>(
      m, "ExtinctionMatrixDataGriddedARO4")
      .doc() = "Extinction matrix data";
  bind_extinction_matrix_data_aro<double, scattering::Representation::Spectral>(
      m, "ExtinctionMatrixDataSpectralARO4")
      .doc() = "Extinction matrix data";

  py::class_<scattering::ParticleProperties>(m, "ParticleProperties")
      .def(py::init<>())
      .def_rw("name", &scattering::ParticleProperties::name, "Name")
      .def_rw("source", &scattering::ParticleProperties::source, "Source")
      .def_rw("refractive_index",
              &scattering::ParticleProperties::refractive_index,
              "Refractive index")
      .def_rw("mass", &scattering::ParticleProperties::mass, "Mass")
      .def_rw("d_veq", &scattering::ParticleProperties::d_veq, "Diameter")
      .def_rw("d_max", &scattering::ParticleProperties::d_max, "Max diameter")
      .doc() = "Particle properties";


  bind_single_scattering_data<double,
                              scattering::Format::TRO,
                              scattering::Representation::Gridded>(
      m, "SingleScatteringDataTROGridded4")
      .doc() = "Single scattering data";

  bind_single_scattering_data<double,
                              scattering::Format::ARO,
                              scattering::Representation::Gridded>(
      m, "SingleScatteringDataAROGridded1")
      .doc() = "Single scattering data";
  bind_single_scattering_data<double,
                              scattering::Format::TRO,
                              scattering::Representation::Spectral>(
      m, "SingleScatteringDataTROSpectral4")
      .doc() = "Single scattering data";
  bind_single_scattering_data<double,
                              scattering::Format::ARO,
                              scattering::Representation::Spectral>(
      m, "SingleScatteringDataAROSpectral1")
      .doc() = "Single scattering data";

  bind_bulk_scattering_properties<scattering::Format::TRO,
                                  scattering::Representation::Gridded>(
      m, "BulkScatteringPropertiesTROGridded4")
      .doc() = "Bulk scattering properties";
  bind_bulk_scattering_properties<scattering::Format::TRO,
                                  scattering::Representation::Spectral>(
      m, "BulkScatteringPropertiesTROSpectral4")
      .doc() = "Bulk scattering properties";



  py::class_<ParticleHabit>(m, "ParticleHabit")
      .def_static("from_legacy_tro",
                  &ParticleHabit::from_legacy_tro,
                  "ssd"_a,
                  "smd"_a,
                  "Create from legacy TRO")
      .def("to_tro_spectral",
           &ParticleHabit::to_tro_spectral,
           "t_grid"_a,
           "f_grid"_a,
           "l"_a,
           "Convert scattering data to TRO spectral format")
      .def("__getitem__",
           [](ParticleHabit &habit, Index ind) {return habit[ind % habit.size()];},
           py::rv_policy::reference_internal)
      .def("__len__", [](ParticleHabit &habit) {return habit.size();})
      .def("size", [](ParticleHabit &habit) {return habit.size();})
      .def("get_sizes", [](ParticleHabit &habit, const SizeParameter &param) {return habit.get_sizes(param);})
      .def("get_size_mass_info", [](ParticleHabit &habit, const SizeParameter &param) {return habit.get_size_mass_info(param);})
      .doc() = "Particle habit";

  py::class_<ScatteringHabit>(m, "ScatteringHabit")
    .def(py::init<const ParticleHabit &, const PSD &, Numeric, Numeric>(),
         py::arg("particle_habit"),
         py::arg("psd"),
         py::arg("mass_size_rel_a") = -1.0,
         py::arg("mass_size_rel_b") = -1.0)
    .def("get_bulk_scattering_properties_tro_spectral",
         &ScatteringHabit::get_bulk_scattering_properties_tro_spectral,
         "point"_a,
         "f_grid"_a,
         "f_tol"_a)
    .doc() = "A scattering habit combines a particle habit with a PSD so that it can be used as a scattering species.";

} catch (std::exception& e) {
  throw std::runtime_error(std::format(
      "DEV ERROR:\nCannot initialize scattering species:\n{}", e.what()));
};
}  // namespace Python
