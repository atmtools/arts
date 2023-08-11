#include <python_interface.h>

#include "py_macros.h"


namespace Python {
void py_star(py::module_& m) try {
  artsclass<Sun>(m, "Sun")
      .def(py::init([]() { return std::make_shared<Sun>(); }), "Empty sun")
      .PythonInterfaceCopyValue(Sun)
//      .PythonInterfaceWorkspaceVariableConversion(Sun)
      .PythonInterfaceBasicRepresentation(Sun)
      .PythonInterfaceFileIO(Sun)
      .def_readwrite("description", &Sun::description, ":class:`~pyarts.arts.String` Sun description")
      .def_readwrite("spectrum", &Sun::spectrum, ":class:`~pyarts.arts.Matrix` Sun spectrum, monochrmatic radiance spectrum at the surface of the sun")
      .def_readwrite("radius", &Sun::radius, ":class:`float` Sun radius")
      .def_readwrite("distance", &Sun::distance, ":class:`float` Sun distance")
      .def_readwrite("latitude", &Sun::latitude, ":class:`float` Sun latitude")
      .def_readwrite("longitude", &Sun::longitude, ":class:`float` Sun longitude")
      .def(py::pickle(
          [](const Sun& self) {
            return py::make_tuple(self.description,
                                  self.spectrum,
                                  self.radius,
                                  self.distance,
                                  self.latitude,
                                  self.longitude);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 6, "Invalid state!")
            return std::make_shared<Sun>(Sun{t[0].cast<String>(),
                             t[1].cast<Matrix>(),
                             t[2].cast<Numeric>(),
                             t[3].cast<Numeric>(),
                             t[4].cast<Numeric>(),
                             t[5].cast<Numeric>()});
          })).doc()=R"--(A single sun.
          
Each sun is described by a struct with its spectrum, radius
distance from center of planet to center of sun,
temperature (if possible), latitude in the sky of the planet,
longitude in the sky of the planet and the type )--";

  PythonInterfaceWorkspaceArray(Sun);
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize star\n", e.what()));
}
}  // namespace Python