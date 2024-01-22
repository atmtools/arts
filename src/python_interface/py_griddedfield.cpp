#include <pybind11/cast.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <python_interface.h>

#include <functional>
#include <optional>
#include <stdexcept>
#include <variant>
#include <vector>

#include "array.h"
#include "debug.h"
#include "details.h"
#include "gridded_data.h"
#include "matpack.h"
#include "mystring.h"
#include "py_macros.h"

namespace Python {

template <typename T, typename... Grids>
auto artsgf(py::module_& m, const char* name) {
  using GF = matpack::gridded_data<T, Grids...>;

  auto gf =
      artsclass<GF>(m, name)
          .def(py::init([]() { return std::make_shared<GF>(); }), "Empty field")
          .def(py::init<GF>(), "Copy field")
          .def(py::init<typename GF::data_t,
                        typename GF::grids_t,
                        std::string,
                        std::array<std::string, GF::dim>>(),
               "Full field",
               py::arg("data"),
               py::arg("grids"),
               py::arg("dataname") = "",
               py::arg("gridnames") = std::array<std::string, GF::dim>{});

  gf.def_readwrite("data", &GF::data, "Data field");
  gf.def_readwrite("dataname", &GF::data_name, "Data name");
  gf.def_readwrite("grids", &GF::grids, "Grid fields");
  gf.def_readwrite("gridnames", &GF::grid_names, "Grid name");
  gf.def_property_readonly_static(
      "dim",
      [](const py::object&) { return GF::dim; },
      "Dimension of the field");
  gf.def_property_readonly(
      "shape",
      [](const GF& gd) { return py::tuple(py::cast(gd.shape())); },
      "Shape of the grids field");
  gf.def("check", &GF::check, "Check the field");
  gf.def("__getitem__", [](py::object& x, py::object& y) {
    return x.attr("data").attr("__getitem__")(y);
  });
  gf.def("__setitem__", [](py::object& x, py::object& y, py::object& z) {
    return x.attr("data").attr("__setitem__")(y, z);
  });
  gf.def("__repr__", [](const GF& gd) { return var_string(gd); });
  gf.def("__str__", [](const GF& gd) { return var_string(gd); });

  gf.def(py::pickle(
      [](const py::object& self) {
        return py::make_tuple(self.attr("data"),
                              self.attr("grids"),
                              self.attr("dataname"),
                              self.attr("gridnames"));
      },
      [](const py::tuple& t) {
        ARTS_USER_ERROR_IF(t.size() != 4, "Invalid state!")
        auto gf = std::make_shared<GF>();
        gf->data = t[0].cast<typename GF::data_t>();
        gf->grids = t[1].cast<typename GF::grids_t>();
        gf->data_name = t[2].cast<std::string>();
        gf->grid_names = t[3].cast<std::array<std::string, GF::dim>>();
        return gf;
      }));

  gf.def("to_dict", [](const py::object& gd) {
    py::dict out;

    out["dims"] = py::list{};
    out["coords"] = py::dict{};
    for (Size i = 0; i < GF::dim; ++i) {
      auto gridname = gd.attr("gridnames").attr("__getitem__")(i);
      auto grid = gd.attr("grids").attr("__getitem__")(i);

      const auto dim = py::cast<py::bool_>(gridname)
                           ? gridname
                           : py::str(var_string("dim", i));
      out["dims"].attr("append")(dim);
      out["coords"][dim] = py::dict{};
      out["coords"][dim]["data"] = grid;
      out["coords"][dim]["dims"] = dim;
    }

    out["name"] = gd.attr("dataname");
    out["data"] = gd.attr("data");

    return out;
  });

  gf.def("to_xarray", [](const py::object& gf) {
    py::module_ xarray = py::module_::import("xarray");
    return xarray.attr("DataArray").attr("from_dict")(gf.attr("to_dict")());
  });

  gf.def_static("from_dict", [](const py::dict& d) {
    if (not d.contains("coords"))
      throw std::invalid_argument(
          "cannot convert dict without the key 'coords''");
    if (not d.contains("dims"))
      throw std::invalid_argument(
          "cannot convert dict without the key 'dims''");
    if (not d.contains("data"))
      throw std::invalid_argument(
          "cannot convert dict without the key 'data''");

    auto gd = std::make_shared<GF>();

    gd->data = d["data"].cast<typename GF::data_t>();
    if (d.contains("name")) {
      gd->data_name = d["name"].cast<std::string>();
    }

    gd->grid_names = d["dims"].cast<std::array<std::string, GF::dim>>();

    py::list coords{};
    for (auto& dim : gd->grid_names) {
      coords.append(d["coords"][py::str(dim)]["data"]);
    }
    gd->grids = coords.cast<typename GF::grids_t>();

    return gd;
  });

  gf.def_static("from_xarray", [](const py::object& xarray) {
    return py::type::of<GF>().attr("from_dict")(xarray.attr("to_dict")());
  });

  gf.def("__eq__", [](const GF& a, const GF& b) {
    return a == b;
  });

  return gf;
}

using VectorOrArrayOfString = std::variant<Vector, ArrayOfString>;

void py_griddedfield(py::module_& m) try {
  artsgf<Numeric, Vector>(m, "GriddedField1")
      .PythonInterfaceFileIO(GriddedField1)
      .PythonInterfaceWorkspaceDocumentation(GriddedField1);

  artsgf<Numeric, Vector, Vector>(m, "GriddedField2")
      .PythonInterfaceFileIO(GriddedField2)
      .PythonInterfaceWorkspaceDocumentation(GriddedField2);

  artsgf<Numeric, Vector, Vector, Vector>(m, "GriddedField3")
      .PythonInterfaceFileIO(GriddedField3)
      .PythonInterfaceWorkspaceDocumentation(GriddedField3);

  artsgf<Numeric, Vector, Vector, Vector, Vector>(m, "GriddedField4")
      .PythonInterfaceFileIO(GriddedField4)
      .PythonInterfaceWorkspaceDocumentation(GriddedField4);

  artsgf<Numeric, Vector, Vector, Vector, Vector, Vector>(m, "GriddedField5")
      .PythonInterfaceFileIO(GriddedField5)
      .PythonInterfaceWorkspaceDocumentation(GriddedField5);

  artsgf<Numeric, Vector, Vector, Vector, Vector, Vector, Vector>(
      m, "GriddedField6")
      .PythonInterfaceFileIO(GriddedField6)
      .PythonInterfaceWorkspaceDocumentation(GriddedField6);

  artsgf<Numeric, ArrayOfString, Vector, Vector>(m, "NamedGriddedField2")
      .PythonInterfaceFileIO(NamedGriddedField2)
      .PythonInterfaceWorkspaceDocumentation(NamedGriddedField2);
  artsgf<Numeric, ArrayOfString, Vector, Vector, Vector>(m,
                                                         "NamedGriddedField3")
      .PythonInterfaceFileIO(NamedGriddedField3)
      .PythonInterfaceWorkspaceDocumentation(NamedGriddedField3);

  artsgf<Numeric, Vector, ArrayOfString>(m, "GriddedField1Named")
      .PythonInterfaceFileIO(GriddedField1Named)
      .PythonInterfaceWorkspaceDocumentation(GriddedField1Named);

  artsgf<Complex, Vector, Vector>(m, "ComplexGriddedField2")
      .PythonInterfaceFileIO(ComplexGriddedField2)
      .PythonInterfaceWorkspaceDocumentation(ComplexGriddedField2);

  artsarray<ArrayOfGriddedField1>(m, "ArrayOfGriddedField1")
      .PythonInterfaceFileIO(ArrayOfGriddedField1)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfGriddedField1);

  artsarray<ArrayOfGriddedField2>(m, "ArrayOfGriddedField2")
      .PythonInterfaceFileIO(ArrayOfGriddedField2)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfGriddedField2);

  artsarray<ArrayOfGriddedField3>(m, "ArrayOfGriddedField3")
      .PythonInterfaceFileIO(ArrayOfGriddedField3)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfGriddedField3);

  artsarray<ArrayOfGriddedField4>(m, "ArrayOfGriddedField4")
      .PythonInterfaceFileIO(ArrayOfGriddedField4)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfGriddedField4);

  artsarray<ArrayOfArrayOfGriddedField1>(m, "ArrayOfArrayOfGriddedField1")
      .PythonInterfaceFileIO(ArrayOfArrayOfGriddedField1)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfArrayOfGriddedField1);

  artsarray<ArrayOfArrayOfGriddedField2>(m, "ArrayOfArrayOfGriddedField2")
      .PythonInterfaceFileIO(ArrayOfArrayOfGriddedField2)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfArrayOfGriddedField2);

  artsarray<ArrayOfArrayOfGriddedField3>(m, "ArrayOfArrayOfGriddedField3")
      .PythonInterfaceFileIO(ArrayOfArrayOfGriddedField3)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfArrayOfGriddedField3);

  artsarray<ArrayOfNamedGriddedField2>(m, "ArrayOfNamedGriddedField2")
      .PythonInterfaceFileIO(ArrayOfNamedGriddedField2)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfNamedGriddedField2);

  artsarray<ArrayOfGriddedField1Named>(m, "ArrayOfGriddedField1Named")
      .PythonInterfaceFileIO(ArrayOfGriddedField1Named)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfGriddedField1Named);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize gridded field\n", e.what()));
}
}  // namespace Python