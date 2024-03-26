#include <python_interface.h>

#include <stdexcept>
#include <variant>

#include "debug.h"
#include "gridded_data.h"
#include "mystring.h"
#include "py_macros.h"
#include "sorted_grid.h"

namespace Python {

template <typename T, typename... Grids>
auto& fix_artsgf(artsclass<matpack::gridded_data<T, Grids...>>& gf) {
  using GF = matpack::gridded_data<T, Grids...>;

  gf.def(py::init<std::string,
                  typename GF::data_t,
                  std::array<std::string, GF::dim>,
                  typename GF::grids_t>(),
         "Full field",
         py::arg("dataname") = "",
         py::arg("data"),
         py::arg("gridnames") = std::array<std::string, GF::dim>{},
         py::arg("grids"));

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
  gf.def("ok", &GF::ok, "Check the field");
  gf.def("__getitem__", [](py::object& x, py::object& y) {
    return x.attr("data").attr("__getitem__")(y);
  });
  gf.def("__setitem__", [](py::object& x, py::object& y, py::object& z) {
    return x.attr("data").attr("__setitem__")(y, z);
  });

  gf.def(py::pickle(
      [](const py::object& self) {
        return py::make_tuple(self.attr("data"),
                              self.attr("grids"),
                              self.attr("dataname"),
                              self.attr("gridnames"));
      },
      [](const py::tuple& t) {
        ARTS_USER_ERROR_IF(t.size() != 4, "Invalid state!")
        auto gd = std::make_shared<GF>();
        gd->data = t[0].cast<typename GF::data_t>();
        gd->grids = t[1].cast<typename GF::grids_t>();
        gd->data_name = t[2].cast<std::string>();
        gd->grid_names = t[3].cast<std::array<std::string, GF::dim>>();
        return gd;
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

  gf.def("to_xarray", [](const py::object& gd) {
    py::module_ xarray = py::module_::import("xarray");
    return xarray.attr("DataArray").attr("from_dict")(gd.attr("to_dict")());
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

  gf.def("__eq__", [](const GF& a, const GF& b) { return a == b; });

  return gf;
}

template <typename T, typename... Grids>
auto artsgf(py::module_& m, const char* name) {
  using GF = matpack::gridded_data<T, Grids...>;

  auto gf =
      artsclass<GF>(m, name)
          .def(py::init([]() { return std::make_shared<GF>(); }), "Empty field")
          .def(py::init<GF>(), "Copy field");

  fix_artsgf(gf);

  gf.def("__repr__", [](const GF& gd) { return var_string(gd); });
  gf.def("__str__", [](const GF& gd) { return var_string(gd); });

  return gf;
}

using VectorOrArrayOfString = std::variant<Vector, ArrayOfString>;

void py_griddedfield(py::module_& m) try {
  fix_artsgf(py_staticGriddedField1(m));
  fix_artsgf(py_staticGriddedField2(m));
  fix_artsgf(py_staticGriddedField3(m));
  fix_artsgf(py_staticGriddedField4(m));
  fix_artsgf(py_staticGriddedField5(m));
  fix_artsgf(py_staticGriddedField6(m));

  fix_artsgf(py_staticNamedGriddedField2(m));
  fix_artsgf(py_staticNamedGriddedField3(m));

  fix_artsgf(py_staticGriddedField1Named(m));

  fix_artsgf(py_staticComplexGriddedField2(m));

  fix_artsgf(py_staticStokvecGriddedField6(m));

  py_staticArrayOfGriddedField1(m);
  py_staticArrayOfGriddedField2(m);
  py_staticArrayOfGriddedField3(m);
  py_staticArrayOfGriddedField4(m);

  py_staticArrayOfArrayOfGriddedField1(m);
  py_staticArrayOfArrayOfGriddedField2(m);
  py_staticArrayOfArrayOfGriddedField3(m);

  py_staticArrayOfNamedGriddedField2(m);

  py_staticArrayOfGriddedField1Named(m);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize gridded field\n", e.what()));
}
}  // namespace Python