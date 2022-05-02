#include <constants.h>
#include <pybind11/pybind11.h>

namespace Python {
namespace py = pybind11;

struct ConstantDummy {};
struct ConversionDummy {};

void py_constants(py::module_& m) {
  py::class_<ConstantDummy>(m, "constant")
      .def_readonly_static("pi", &Constant::pi);

  py::class_<ConversionDummy>(m, "conversion")
      .def_static("freq2kaycm",
                  [](double x) { return Conversion::freq2kaycm(x); })
      .def_static("kaycm2freq",
                  [](double x) { return Conversion::kaycm2freq(x); });
}
}  // namespace Python