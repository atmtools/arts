#include <core/scattering/particle_habit.h>
#include <core/scattering/single_scattering_data.h>
#include <core/scattering/properties.h>
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

void py_psd(py::module_& m) try {

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
         "picky"_a
         )
    .def(py::init<ScatteringSpeciesProperty, std::string, Numeric, Numeric, bool>(),
         "properties"_a,
         "name"_a,
         "t_min"_a,
         "t_max"_a,
         "picky"_a
         )
    .def("evaluate", &scattering::MGDSingleMoment::evaluate, "Evaluate PSD at given point.")
    .doc() = "Modified gamma PSD single moment";

  //
  // BinnedPSD
  //

  py::class_<scattering::BinnedPSD>(m, "BinnedPSD")
    .def(py::init<SizeParameter, Vector, Vector, Numeric, Numeric>(),
         "size_parameter"_a,
         "bins"_a,
         "counts"_a,
         "t_min"_a,
         "t_max"_a
         )
    .def("evaluate", &scattering::BinnedPSD::evaluate, "Evaluate PSD at given point.")
    .doc() = "Binned PSD returning a fixed particle concentration defined over a sequence of size bins with particle number zeros outside of size bins and temperature range.";


 } catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize scattering species:\n{}", e.what()));
 };

}
