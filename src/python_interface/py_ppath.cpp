#include <py_auto_interface.h>

#include "py_macros.h"

#include <ppath_struct.h>

namespace Python {
void py_ppath(py::module_& m) {
  py::class_<GridPos>(m, "GridPos")
      .def(py::init([]() { return std::make_unique<GridPos>(); }), "Empty position")
      .PythonInterfaceCopyValue(GridPos)
      .PythonInterfaceWorkspaceVariableConversion(GridPos)
      .def(py::init([](Index i, std::array<Numeric, 2> x) { return std::make_unique<GridPos>(GridPos{i, x}); }), "Complete position")
      .PythonInterfaceFileIO(GridPos)
      .PythonInterfaceBasicRepresentation(GridPos)
      .def_readwrite("idx", &GridPos::idx, ":class:`~pyarts.arts.Index` Index in a grid")
      .def_readwrite("fd", &GridPos::fd, ":class:`list` The two weights")
      .def(py::pickle(
          [](const GridPos& self) {
            return py::make_tuple(self.idx, self.fd);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 2, "Invalid state!")
            return std::make_unique<GridPos>(GridPos{t[0].cast<Index>(), t[1].cast<std::array<Numeric, 2>>()});
          }))
      .PythonInterfaceWorkspaceDocumentation(GridPos);

  py::class_<ArrayOfGridPos>(m, "ArrayOfGridPos")
      .PythonInterfaceArrayDefault(GridPos)
      .PythonInterfaceBasicRepresentation(ArrayOfGridPos).doc() = "List of :class:`~pyarts.arts.GridPos`";
  py::implicitly_convertible<std::vector<GridPos>, ArrayOfGridPos>();

  py::class_<Ppath>(m, "Ppath")
      .def(py::init([]() { return std::make_unique<Ppath>(); }), "Empty path")
      .PythonInterfaceCopyValue(Ppath)
      .PythonInterfaceWorkspaceVariableConversion(Ppath)
      .def("__repr__", [](Ppath&) { return "Ppath"; })
      .PythonInterfaceFileIO(Ppath)
      .def_readwrite("dim", &Ppath::dim, ":class:`int` Atmospheric dimension")
      .def_readwrite("np", &Ppath::np, ":class:`int` Number of path points")
      .def_readwrite("constant", &Ppath::constant, ":class:`float` The propagation path constant (only used for 1D)")
      .def_readwrite("background", &Ppath::background, ":class:`~pyarts.arts.String` Radiative background")
      .def_readwrite("start_pos", &Ppath::start_pos, ":class:`~pyarts.arts.Vector` Start position")
      .def_readwrite("start_los", &Ppath::start_los, ":class:`~pyarts.arts.Vector` Start line-of-sight")
      .def_readwrite("start_lstep", &Ppath::start_lstep, ":class:`float` Length between sensor and atmospheric boundary")
      .def_readwrite("pos", &Ppath::pos, ":class:`~pyarts.arts.Matrix` The distance between start pos and the last position in pos")
      .def_readwrite("los", &Ppath::los, ":class:`~pyarts.arts.Matrix` Line-of-sight at each ppath point")
      .def_readwrite("r", &Ppath::r, ":class:`~pyarts.arts.Vector` Radius of each ppath point")
      .def_readwrite("lstep", &Ppath::lstep, ":class:`~pyarts.arts.Vector` The length between ppath points")
      .def_readwrite("end_pos", &Ppath::end_pos, ":class:`~pyarts.arts.Vector` End position")
      .def_readwrite("end_los", &Ppath::end_los, ":class:`~pyarts.arts.Vector` End line-of-sight")
      .def_readwrite("end_lstep", &Ppath::end_lstep, ":class:`float` The distance between end pos and the first position in pos")
      .def_readwrite("nreal", &Ppath::nreal, ":class:`~pyarts.arts.Vector` The real part of the refractive index at each path position")
      .def_readwrite("ngroup", &Ppath::ngroup, ":class:`~pyarts.arts.Vector` The group index of refraction")
      .def_readwrite("gp_p", &Ppath::gp_p, ":class:`~pyarts.arts.ArrayOfGridPos` Index position with respect to the pressure grid")
      .def_readwrite("gp_lat", &Ppath::gp_lat, ":class:`~pyarts.arts.ArrayOfGridPos` Index position with respect to the latitude grid")
      .def_readwrite("gp_lon", &Ppath::gp_lon, ":class:`~pyarts.arts.ArrayOfGridPos` Index position with respect to the longitude grid")
      .def(py::pickle(
          [](const Ppath& self) {
            return py::make_tuple(self.dim,
                                  self.np,
                                  self.constant,
                                  self.background,
                                  self.start_pos,
                                  self.start_los,
                                  self.start_lstep,
                                  self.pos,
                                  self.los,
                                  self.r,
                                  self.lstep,
                                  self.end_pos,
                                  self.end_los,
                                  self.end_lstep,
                                  self.nreal,
                                  self.ngroup,
                                  self.gp_p,
                                  self.gp_lat,
                                  self.gp_lon);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 19, "Invalid state!")
            return std::make_unique<Ppath>(Ppath{t[0].cast<Index>(),
                                                 t[1].cast<Index>(),
                                                 t[2].cast<Numeric>(),
                                                 t[3].cast<String>(),
                                                 t[4].cast<Vector>(),
                                                 t[5].cast<Vector>(),
                                                 t[6].cast<Numeric>(),
                                                 t[7].cast<Matrix>(),
                                                 t[8].cast<Matrix>(),
                                                 t[9].cast<Vector>(),
                                                 t[10].cast<Vector>(),
                                                 t[11].cast<Vector>(),
                                                 t[12].cast<Vector>(),
                                                 t[13].cast<Numeric>(),
                                                 t[14].cast<Vector>(),
                                                 t[15].cast<Vector>(),
                                                 t[16].cast<ArrayOfGridPos>(),
                                                 t[17].cast<ArrayOfGridPos>(),
                                                 t[18].cast<ArrayOfGridPos>()});
          }))
      .PythonInterfaceWorkspaceDocumentation(Ppath);

  PythonInterfaceWorkspaceArray(Ppath);
}
}  // namespace Python