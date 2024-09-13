#include <py_auto_interface.h>

#include "py_macros.h"


namespace Python {
void py_star(py::module_& m) {
  py::class_<Sun>(m, "Sun")
      .def(py::init([]() { return std::make_unique<Sun>(); }), "Empty sun")
      .PythonInterfaceCopyValue(Sun)
//      .PythonInterfaceWorkspaceVariableConversion(Sun)
      .PythonInterfaceBasicRepresentation(Sun)
      .PythonInterfaceFileIO(Sun)
      .def_readwrite("description", &Sun::description, ":class:`~pyarts.arts.String` Sun description")
      .def_readwrite("spectrum", &Sun::spectrum, ":class:`~pyarts.arts.Matrix` Sun spectrum, spectral irradiance at the position of the sun")
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
            return std::make_unique<Sun>(Sun{t[0].cast<String>(),
                             t[1].cast<Matrix>(),
                             t[2].cast<Numeric>(),
                             t[3].cast<Numeric>(),
                             t[4].cast<Numeric>(),
                             t[5].cast<Numeric>()});
          })).doc()=R"--(A single sun.
          
Each sun is described by a struct with its spectral irradiance, radius
distance from center of planet to center of sun, latitude in the sky of the planet,
and longitude in the sky of the planet.  Optionally, a free-form text description of the sun may be available. )--";

  PythonInterfaceWorkspaceArray(Sun);
}
}  // namespace Python