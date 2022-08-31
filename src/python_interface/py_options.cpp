#include <arts_options.h>

#include "py_auto_interface.h"
#include "python_interface.h"

#define DeclareOption(opt_namespace, opt_localname)  {\
    using namespace opt_namespace;\
    auto my_enum_ = py::class_<opt_localname>(\
        opt, #opt_localname);\
    my_enum_.def(py::init([](const std::string& s) {\
      return to##opt_localname##OrThrow(s);\
    }));\
    my_enum_.def_property(\
        "value",\
        [](const opt_localname& x) { return toString(x); },\
        [](opt_localname& x, const std::string& s) {\
          x = to##opt_localname##OrThrow(s);\
        });\
    my_enum_.def("__str__", [](const opt_localname& x) {\
      return toString(x);\
    });\
    my_enum_.doc() =\
        "Options for "\
        #opt_localname;\
\
    for (auto& x : enumtyps::opt_localname##Types) {\
      my_enum_.def_property_readonly_static(\
          std::string(toString(x)).c_str(),\
          [x](const py::object& /* self */) { return toString(x); });\
    }\
    py::implicitly_convertible<std::string, LineShape::TemperatureModel>();\
  }

namespace Python {
void py_options(py::module_& m) {
  auto opt = m.def_submodule("options");
  opt.doc() = "Various named options of Arts";

  DeclareOption(Options, iy_main_agendaDefaultOptions)
  DeclareOption(Options, iy_loop_freqs_agendaDefaultOptions)
  DeclareOption(Options, iy_space_agendaDefaultOptions)
  DeclareOption(Options, iy_surface_agendaDefaultOptions)
  DeclareOption(Options, iy_cloudbox_agendaDefaultOptions)
  DeclareOption(Options, ppath_agendaDefaultOptions)
  DeclareOption(Options, ppath_step_agendaDefaultOptions)
  DeclareOption(Options, refr_index_air_agendaDefaultOptions)
  DeclareOption(Options, water_p_eq_agendaDefaultOptions)
  DeclareOption(Options, gas_scattering_agendaDefaultOptions)
  DeclareOption(Options, surface_rtprop_agendaDefaultOptions)
  DeclareOption(Options, g0_agendaDefaultOptions)
}
}  // namespace Python
