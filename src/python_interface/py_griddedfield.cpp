#include <python_interface.h>

#include <functional>
#include <optional>
#include <stdexcept>
#include <variant>
#include <vector>

#include "array.h"
#include "debug.h"
#include "details.h"
#include "gridded_fields.h"
#include "mystring.h"
#include "py_macros.h"

namespace Python {

namespace details {
struct GriddedField {
  inline static auto extract_slice{three_args};
  inline static auto refine_grid{five_args};
  inline static auto to_xarray{one_arg};
  inline static auto to_dict{one_arg};
  inline static auto from_xarray{two_args};
};
}  // namespace details

using VectorOrArrayOfString = std::variant<Vector, ArrayOfString>;

void py_griddedfield(py::module_& m) try {
  m.add_object("_cleanupGriddedField", py::capsule([]() {
                 details::GriddedField::extract_slice = details::three_args;
                 details::GriddedField::refine_grid = details::five_args;
                 details::GriddedField::to_xarray = details::one_arg;
                 details::GriddedField::to_dict = details::one_arg;
                 details::GriddedField::from_xarray = details::two_args;
               }));

  artsclass<details::GriddedField>(m, "_detailsGriddedField")
      .def_readwrite_static("extract_slice",
                            &details::GriddedField::extract_slice)
      .def_readwrite_static("refine_grid", &details::GriddedField::refine_grid)
      .def_readwrite_static("to_xarray", &details::GriddedField::to_xarray)
      .def_readwrite_static("to_dict", &details::GriddedField::to_dict)
      .def_readwrite_static("from_xarray", &details::GriddedField::from_xarray);

  artsclass<GriddedField>(m, "_GriddedField")
      .def_property_readonly(
          "dim",
          [](GriddedField& gf) { return gf.get_dim(); },
          ":class:`~pyarts.arts.Index` Dimensionality of field")
      .def_property("name",
                    &GriddedField::get_name,
                    &GriddedField::set_name,
                    ":class:`~pyarts.arts.String` Name of field")
      .def("get_grid",
           [](GriddedField& g, Index i) -> VectorOrArrayOfString {
             if (i < 0 or i > g.get_dim()) throw std::out_of_range("Bad axis");
             if (g.get_grid_type(i) == GRID_TYPE_NUMERIC)
               return g.get_numeric_grid(i);
             return g.get_string_grid(i);
           }, "Get the grid at position i\n\nReturn\n------\n    (:class:`~pyarts.arts.ArrayOfString` or :class:`~pyarts.arts.Vector`) Grid values", py::arg("i"))
      .def("set_grid",
           [](GriddedField& g, Index i, const VectorOrArrayOfString& y) {
             if (i < 0 or i > g.get_dim()) throw std::out_of_range("Bad axis");
             std::visit([&](auto&& z) { g.set_grid(i, z); }, y);
           }, "Set the grid at position i", py::arg("i"), py::arg("new_grid"))
      .def("get_grid_name",
           [](GriddedField& g, Index i) {
             if (i < 0 or i > g.get_dim()) throw std::out_of_range("Bad axis");
             return g.get_grid_name(i);
           }, "Get the grid name at position i\n\nReturn\n------\n    (:class:`~pyarts.arts.String`) Grid name", py::arg("i"))
      .def("set_grid_name",
           [](GriddedField& g, Index i, const String& s) {
             if (i < 0 or i > g.get_dim()) throw std::out_of_range("Bad axis");
             g.set_grid_name(i, s);
           }, "Set the grid name at position i", py::arg("i"), py::arg("new_grid_name"))
      .def_property(
          "grids",
          [](GriddedField& g) {
            const Index n = g.get_dim();
            std::vector<VectorOrArrayOfString> o(n);
            for (Index i = 0; i < n; i++) {
              if (g.get_grid_type(i) == GRID_TYPE_NUMERIC) {
                o[i] = g.get_numeric_grid(i);
              } else {
                o[i] = g.get_string_grid(i);
              }
            }
            return o;
          },
          [](GriddedField& g, const std::vector<VectorOrArrayOfString>& y) {
            const Index n = g.get_dim();
            ARTS_USER_ERROR_IF(std::size_t(n) not_eq y.size(), "Bad size");
            for (Index i = 0; i < n; i++) {
              std::visit([&](auto&& z) { g.set_grid(i, z); }, y[i]);
            }
          }, ":class:`list` Grids of the field")
      .def_property(
          "gridnames",
          [](GriddedField& g) {
            const Index n = g.get_dim();
            std::vector<String> o(n);
            for (Index i = 0; i < n; i++) o[i] = g.get_grid_name(i);
            return o;
          },
          [](GriddedField& g, const ArrayOfString& s) {
            const Index n = g.get_dim();
            ARTS_USER_ERROR_IF(std::size_t(n) not_eq s.size(), "Bad size");
            for (Index i = 0; i < n; i++) {
              g.set_grid_name(i, s[i]);
            }
          }, ":class:`~pyarts.arts.String` Grid names of the field")
      .def("checksize_strict", &GriddedField::checksize_strict, "Check the grid and throw if bad")
      .def("checksize", &GriddedField::checksize, "Check the grid and return its state.\n\nReturn\n------\n    (bool) Good state when ``True``")
      .def(
          "extract_slice",
          [](py::object& g, py::object& s, py::object& i) {
            return details::GriddedField::extract_slice(g, s, i);
          },
          py::arg("slice"),
          py::arg("axis") = 0)
      .def(
          "to_xarray",
          [](py::object& g) { return details::GriddedField::to_xarray(g); })
      .def(
          "to_dict",
          [](py::object& g) { return details::GriddedField::to_dict(g); })
      .def(
          "refine_grid",
          [](py::object& g,
             py::object& s,
             py::object& i,
             py::object& t,
             const py::kwargs& kw) {
            py::object in = py::dict{};
            for (auto& a : kw) in[a.first] = a.second;
            return details::GriddedField::refine_grid(g, s, i, t, in);
          },
          py::arg("new_grid"),
          py::arg("axis") = py::int_(0),
          py::arg("type") = py::str("linear"))
      .doc() = "Base class for gridding fields of data";

  artsclass<GriddedField1, GriddedField>(m, "GriddedField1")
      .def(py::init([]() { return std::make_shared<GriddedField1>(); }), "Empty field")
      .def(py::init([](const String& s) { return std::make_shared<GriddedField1>(s); }), "Empty named field")
      .def(py::init([](const std::array<VectorOrArrayOfString, 1>& grids,
                       const Vector& data,
                       const std::optional<std::array<String, 1>>& gridnames,
                       const String& name) {
             constexpr std::size_t n = 1;
             auto out = std::make_shared<GriddedField1>(name);
             for (std::size_t i = 0; i < n; i++)
               std::visit([&](auto&& g) { out->set_grid(i, g); }, grids[i]);
             if (gridnames.has_value())
               for (std::size_t i = 0; i < n; i++)
                 out->set_grid_name(i, gridnames->at(i));
             out->data = data;
             return out;
           }),
           py::arg("grids"),
           py::arg("data"),
           py::arg("gridnames") = std::nullopt,
           py::arg_v("name", "", "None"), py::doc("Complete field"))
      .PythonInterfaceCopyValue(GriddedField1)
      .PythonInterfaceWorkspaceVariableConversion(GriddedField1)
      .PythonInterfaceFileIO(GriddedField1)
      .PythonInterfaceBasicRepresentation(GriddedField1)
      .PythonInterfaceGriddedField(GriddedField1)
      .PythonInterfaceWorkspaceDocumentation(GriddedField1);

  artsclass<GriddedField2, GriddedField>(m, "GriddedField2")
      .def(py::init([]() { return std::make_shared<GriddedField2>(); }), "Empty field")
      .def(py::init([](const String& s) { return std::make_shared<GriddedField2>(s); }), "Empty named field")
      .def(py::init([](const std::array<VectorOrArrayOfString, 2>& grids,
                       const Matrix& data,
                       const std::optional<std::array<String, 2>>& gridnames,
                       const String& name) {
             constexpr std::size_t n = 2;
             auto out = std::make_shared<GriddedField2>(name);
             for (std::size_t i = 0; i < n; i++)
               std::visit([&](auto&& g) { out->set_grid(i, g); }, grids[i]);
             if (gridnames.has_value())
               for (std::size_t i = 0; i < n; i++)
                 out->set_grid_name(i, gridnames->at(i));
             out->data = data;
             return out;
           }),
           py::arg("grids"),
           py::arg("data"),
           py::arg("gridnames") = std::nullopt,
           py::arg_v("name", "", "None"), "Complete field")
      .PythonInterfaceCopyValue(GriddedField2)
      .PythonInterfaceWorkspaceVariableConversion(GriddedField2)
      .PythonInterfaceFileIO(GriddedField2)
      .PythonInterfaceBasicRepresentation(GriddedField2)
      .PythonInterfaceGriddedField(GriddedField2)
      .PythonInterfaceWorkspaceDocumentation(GriddedField2);

  artsclass<GriddedField3, GriddedField>(m, "GriddedField3")
      .def(py::init([]() { return std::make_shared<GriddedField3>(); }), "Empty field")
      .def(py::init([](const String& s) { return std::make_shared<GriddedField3>(s); }), "Empty named field")
      .def(py::init([](const std::array<VectorOrArrayOfString, 3>& grids,
                       const Tensor3& data,
                       const std::optional<std::array<String, 3>>& gridnames,
                       const String& name) {
             constexpr std::size_t n = 3;
             auto out = std::make_shared<GriddedField3>(name);
             for (std::size_t i = 0; i < n; i++)
               std::visit([&](auto&& g) { out->set_grid(i, g); }, grids[i]);
             if (gridnames.has_value())
               for (std::size_t i = 0; i < n; i++)
                 out->set_grid_name(i, gridnames->at(i));
             out->data = data;
             return out;
           }),
           py::arg("grids"),
           py::arg("data"),
           py::arg("gridnames") = std::nullopt,
           py::arg_v("name", "", "None"), "Complete field")
      .PythonInterfaceCopyValue(GriddedField3)
      .PythonInterfaceWorkspaceVariableConversion(GriddedField3)
      .PythonInterfaceFileIO(GriddedField3)
      .PythonInterfaceBasicRepresentation(GriddedField3)
      .PythonInterfaceGriddedField(GriddedField3)
      .PythonInterfaceWorkspaceDocumentation(GriddedField3);

  artsclass<GriddedField4, GriddedField>(m, "GriddedField4")
      .def(py::init([]() { return std::make_shared<GriddedField4>(); }), "Empty field")
      .def(py::init([](const String& s) { return std::make_shared<GriddedField4>(s); }), "Empty named field")
      .def(py::init([](const std::array<VectorOrArrayOfString, 4>& grids,
                       const Tensor4& data,
                       const std::optional<std::array<String, 4>>& gridnames,
                       const String& name) {
             constexpr std::size_t n = 4;
             auto out = std::make_shared<GriddedField4>(name);
             for (std::size_t i = 0; i < n; i++)
               std::visit([&](auto&& g) { out->set_grid(i, g); }, grids[i]);
             if (gridnames.has_value())
               for (std::size_t i = 0; i < n; i++)
                 out->set_grid_name(i, gridnames->at(i));
             out->data = data;
             return out;
           }),
           py::arg("grids"),
           py::arg("data"),
           py::arg("gridnames") = std::nullopt,
           py::arg_v("name", "", "None"), "Complete field")
      .PythonInterfaceCopyValue(GriddedField4)
      .PythonInterfaceWorkspaceVariableConversion(GriddedField4)
      .PythonInterfaceFileIO(GriddedField4)
      .PythonInterfaceBasicRepresentation(GriddedField4)
      .PythonInterfaceGriddedField(GriddedField4)
      .PythonInterfaceWorkspaceDocumentation(GriddedField4);

  artsclass<GriddedField5, GriddedField>(m, "GriddedField5")
      .def(py::init([]() { return std::make_shared<GriddedField5>(); }), "Empty field")
      .def(py::init([](const String& s) { return std::make_shared<GriddedField5>(s); }), "Empty named field")
      .def(py::init([](const std::array<VectorOrArrayOfString, 5>& grids,
                       const Tensor5& data,
                       const std::optional<std::array<String, 5>>& gridnames,
                       const String& name) {
             constexpr std::size_t n = 5;
             auto out = std::make_shared<GriddedField5>(name);;
             for (std::size_t i = 0; i < n; i++)
               std::visit([&](auto&& g) { out->set_grid(i, g); }, grids[i]);
             if (gridnames.has_value())
               for (std::size_t i = 0; i < n; i++)
                 out->set_grid_name(i, gridnames->at(i));
             out->data = data;
             return out;
           }),
           py::arg("grids"),
           py::arg("data"),
           py::arg("gridnames") = std::nullopt,
           py::arg_v("name", "", "None"), "Complete field")
      .PythonInterfaceCopyValue(GriddedField5)
      .PythonInterfaceWorkspaceVariableConversion(GriddedField5)
      .PythonInterfaceFileIO(GriddedField5)
      .PythonInterfaceBasicRepresentation(GriddedField5)
      .PythonInterfaceGriddedField(GriddedField5)
      .PythonInterfaceWorkspaceDocumentation(GriddedField5);

  artsclass<GriddedField6, GriddedField>(m, "GriddedField6")
      .def(py::init([]() { return std::make_shared<GriddedField6>(); }), "Empty field")
      .def(py::init([](const String& s) { return std::make_shared<GriddedField6>(s); }), "Empty named field")
      .def(py::init([](const std::array<VectorOrArrayOfString, 6>& grids,
                       const Tensor6& data,
                       const std::optional<std::array<String, 6>>& gridnames,
                       const String& name) {
             constexpr std::size_t n = 6;
             auto out = std::make_shared<GriddedField6>(name);
             for (std::size_t i = 0; i < n; i++)
               std::visit([&](auto&& g) { out->set_grid(i, g); }, grids[i]);
             if (gridnames.has_value())
               for (std::size_t i = 0; i < n; i++)
                 out->set_grid_name(i, gridnames->at(i));
             out->data = data;
             return out;
           }),
           py::arg("grids"),
           py::arg("data"),
           py::arg("gridnames") = std::nullopt,
           py::arg_v("name", "", "None"), "Complete field")
      .PythonInterfaceCopyValue(GriddedField6)
      .PythonInterfaceWorkspaceVariableConversion(GriddedField6)
      .PythonInterfaceFileIO(GriddedField6)
      .PythonInterfaceBasicRepresentation(GriddedField6)
      .PythonInterfaceGriddedField(GriddedField6)
      .PythonInterfaceWorkspaceDocumentation(GriddedField6);

  PythonInterfaceWorkspaceArray(GriddedField1);
  PythonInterfaceWorkspaceArray(GriddedField2);
  PythonInterfaceWorkspaceArray(GriddedField3);
  PythonInterfaceWorkspaceArray(GriddedField4);

  PythonInterfaceWorkspaceArray(ArrayOfGriddedField1)
      .def(py::init([](const std::vector<std::vector<GriddedField1>>& x) {
        ArrayOfArrayOfGriddedField1 y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }), "Initialize from lists");
  py::implicitly_convertible<std::vector<std::vector<GriddedField1>>,
                             ArrayOfArrayOfGriddedField1>();

  PythonInterfaceWorkspaceArray(ArrayOfGriddedField2)
      .def(py::init([](const std::vector<std::vector<GriddedField2>>& x) {
        ArrayOfArrayOfGriddedField2 y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }), "Initialize from lists");
  py::implicitly_convertible<std::vector<std::vector<GriddedField2>>,
                             ArrayOfArrayOfGriddedField2>();

  PythonInterfaceWorkspaceArray(ArrayOfGriddedField3)
      .def(py::init([](const std::vector<std::vector<GriddedField3>>& x) {
        ArrayOfArrayOfGriddedField3 y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }), "Initialize from lists");
  py::implicitly_convertible<std::vector<std::vector<GriddedField3>>,
                             ArrayOfArrayOfGriddedField3>();
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize gridded field\n", e.what()));
}
}  // namespace Python