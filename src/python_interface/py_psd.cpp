#include <core/scattering/particle_habit.h>
#include <core/scattering/properties.h>
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

namespace Python {

void py_psd(py::module_& m) try {
  py::class_<scattering::PSDData>(m, "PSDData")
      .def_ro("values", &scattering::PSDData::values, "PSD values\n\n.. :class:`~pyarts3.arts.Vector`")
      .def(
          "derivative",
          [](const scattering::PSDData& data, const ScatteringSpeciesProperty& property) -> const Vector& {
            return data.derivatives.at(property);
          },
          "property"_a,
          "Return the PSD derivative for an atmospheric property.",
          py::rv_policy::reference_internal)
      .def_prop_ro(
          "derivative_properties",
          [](const scattering::PSDData& data) {
            std::vector<ScatteringSpeciesProperty> result;
            result.reserve(data.derivatives.size());
            for (const auto& [property, derivative] : data.derivatives) {
              static_cast<void>(derivative);
              result.push_back(property);
            }
            return result;
          },
          "Atmospheric properties with available PSD derivatives\n\n.. :class:`list[ScatteringSpeciesProperty]`")
      .doc() = "Particle-size distribution values and derivatives by atmospheric property.";

  py::class_<scattering::MonodispersePSD>(m, "MonodispersePSD")
      .def(py::init<ScatteringSpeciesProperty, Numeric, Numeric>(),
           "number_density"_a,
           "t_min"_a = 0.0,
           "t_max"_a = 350.0)
      .def("evaluate", &scattering::MonodispersePSD::evaluate, "Evaluate PSD at an atmospheric point.")
      .def("evaluate_with_derivatives",
           &scattering::MonodispersePSD::evaluate_with_derivatives,
           "Evaluate the PSD and its derivatives at an atmospheric point.")
      .doc() = "A monodisperse population scaled by an atmospheric number-density field.";

  py::enum_<scattering::MGDTwoMomentType>(m, "MGDTwoMomentType")
      .value("MassMeanSize", scattering::MGDTwoMomentType::MassMeanSize)
      .value("MassNumberDensity", scattering::MGDTwoMomentType::MassNumberDensity)
      .export_values();

  //
  // Modified Gamma Single Moment
  //

  py::class_<scattering::MGDSingleMoment>(m, "MGDSingleMoment")
      .def(py::init<ScatteringSpeciesProperty, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, bool>(),
           "properties"_a,
           "n_alpha"_a,
           "n_b"_a,
           "mu"_a,
           "gamma"_a,
           "t_min"_a,
           "t_max"_a,
           "picky"_a)
      .def(py::init<ScatteringSpeciesProperty, std::string, Numeric, Numeric, bool>(),
           "properties"_a,
           "name"_a,
           "t_min"_a,
           "t_max"_a,
           "picky"_a)
      .def("evaluate", &scattering::MGDSingleMoment::evaluate, "Evaluate PSD at given point.")
      .def("evaluate_with_derivatives",
           &scattering::MGDSingleMoment::evaluate_with_derivatives,
           "Evaluate the PSD and its derivatives at the given point.")
      .doc() = "Modified gamma PSD single moment";

  py::class_<scattering::MGDMass>(m, "MGDMass")
      .def(py::init<ScatteringSpeciesProperty, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, bool>(),
           "mass"_a,
           "n0"_a,
           "mu"_a,
           "lambda_"_a,
           "gamma"_a,
           "t_min"_a,
           "t_max"_a,
           "picky"_a = false)
      .def("evaluate", &scattering::MGDMass::evaluate, "Evaluate PSD at given point.")
      .def("evaluate_with_derivatives",
           &scattering::MGDMass::evaluate_with_derivatives,
           "Evaluate the PSD and its derivatives at the given point.")
      .doc() = "Modified-gamma PSD constrained by mass density; exactly one of ``n0`` and ``lambda_`` must be NaN.";

  py::class_<scattering::MGDTwoMoment>(m, "MGDTwoMoment")
      .def(py::init<ScatteringSpeciesProperty,
                    ScatteringSpeciesProperty,
                    scattering::MGDTwoMomentType,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    bool>(),
           "mass"_a,
           "second_moment"_a,
           "type"_a,
           "mu"_a,
           "gamma"_a,
           "t_min"_a,
           "t_max"_a,
           "picky"_a = false)
      .def("evaluate", &scattering::MGDTwoMoment::evaluate, "Evaluate PSD at given point.")
      .def("evaluate_with_derivatives",
           &scattering::MGDTwoMoment::evaluate_with_derivatives,
           "Evaluate the PSD and its derivatives at the given point.")
      .doc() = "Modified-gamma PSD constrained by mass density and a second moment.";

  py::class_<scattering::DelanoeEtAl14>(m, "DelanoeEtAl14")
      .def(py::init<ScatteringSpeciesProperty, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, bool>(),
           "mass"_a,
           "rho"_a    = 916.7,
           "alpha"_a  = -0.237,
           "beta"_a   = 1.839,
           "t_min"_a  = 0.0,
           "t_max"_a  = 400.0,
           "dm_min"_a = -1.0,
           "picky"_a  = false)
      .def(py::init<ScatteringSpeciesProperty,
                    ScatteringSpeciesProperty,
                    ScatteringSpeciesProperty,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    bool>(),
           "mass"_a,
           "intercept_parameter"_a,
           "mean_size"_a,
           "rho"_a    = 916.7,
           "alpha"_a  = -0.237,
           "beta"_a   = 1.839,
           "t_min"_a  = 0.0,
           "t_max"_a  = 400.0,
           "dm_min"_a = -1.0,
           "picky"_a  = false)
      .def("evaluate", &scattering::DelanoeEtAl14::evaluate, "Evaluate PSD at given point.")
      .def("evaluate_with_derivatives",
           &scattering::DelanoeEtAl14::evaluate_with_derivatives,
           "Evaluate the PSD and its derivatives at the given point.")
      .doc() = "Delanoe et al. (2014) normalized ice PSD.";

  py::class_<scattering::FieldEtAl07>(m, "FieldEtAl07")
      .def(py::init<ScatteringSpeciesProperty, std::string, Numeric, Numeric, Numeric, Numeric, bool>(),
           "mass"_a,
           "regime"_a,
           "t_min"_a     = 0.0,
           "t_max"_a     = 290.0,
           "t_min_psd"_a = 200.0,
           "t_max_psd"_a = 273.15,
           "picky"_a     = false)
      .def("evaluate", &scattering::FieldEtAl07::evaluate, "Evaluate PSD at given point.")
      .def("evaluate_with_derivatives",
           &scattering::FieldEtAl07::evaluate_with_derivatives,
           "Evaluate the PSD and its derivatives at the given point.")
      .doc() = "Field et al. (2007) tropical or midlatitude ice PSD.";

  py::class_<scattering::McFarquharHeymsfield97>(m, "McFarquharHeymsfield97")
      .def(py::init<ScatteringSpeciesProperty, Numeric, Numeric, Numeric, Numeric, bool>(),
           "mass"_a,
           "t_min"_a     = 0.0,
           "t_max"_a     = 280.0,
           "t_min_psd"_a = 180.0,
           "t_max_psd"_a = 273.15,
           "picky"_a     = false)
      .def("evaluate", &scattering::McFarquharHeymsfield97::evaluate, "Evaluate PSD at given point.")
      .def("evaluate_with_derivatives",
           &scattering::McFarquharHeymsfield97::evaluate_with_derivatives,
           "Evaluate the PSD and its derivatives at the given point.")
      .doc() = "McFarquhar and Heymsfield (1997) cloud-ice PSD.";

  //
  // BinnedPSD
  //

  py::class_<scattering::BinnedPSD>(m, "BinnedPSD")
      .def(py::init<SizeParameter, Vector, Vector, Numeric, Numeric>(),
           "size_parameter"_a,
           "bins"_a,
           "counts"_a,
           "t_min"_a,
           "t_max"_a)
      .def("evaluate", &scattering::BinnedPSD::evaluate, "Evaluate PSD at given point.")
      .def("evaluate_with_derivatives",
           &scattering::BinnedPSD::evaluate_with_derivatives,
           "Evaluate the PSD and its derivatives at the given point.")
      .doc() =
      "Binned PSD returning a fixed particle concentration defined over a sequence of size bins with particle number zeros outside of size bins and temperature range.";

} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize scattering species:\n{}", e.what()));
};

}  // namespace Python
