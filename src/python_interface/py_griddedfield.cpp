#include <py_auto_interface.h>
#include <pybind11/cast.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

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

void py_griddedfield(py::module_& m) {
  m.add_object("_cleanupGriddedField", py::capsule([]() {
                 details::GriddedField::extract_slice = details::three_args;
                 details::GriddedField::refine_grid = details::five_args;
                 details::GriddedField::to_xarray = details::one_arg;
                 details::GriddedField::to_dict = details::one_arg;
                 details::GriddedField::from_xarray = details::two_args;
               }));

  py::class_<details::GriddedField>(m, "GriddedField::details")
      .def_readwrite_static("extract_slice",
                            &details::GriddedField::extract_slice)
      .def_readwrite_static("refine_grid", &details::GriddedField::refine_grid)
      .def_readwrite_static("to_xarray", &details::GriddedField::to_xarray)
      .def_readwrite_static("to_dict", &details::GriddedField::to_dict)
      .def_readwrite_static("from_xarray", &details::GriddedField::from_xarray);

  py::class_<GriddedField>(m, "GriddedField")
      .def_property_readonly("dim",
                             [](GriddedField& gf) { return gf.get_dim(); })
      .def_property("name", &GriddedField::get_name, &GriddedField::set_name)
      .def("get_grid",
           [](GriddedField& g, Index i) -> VectorOrArrayOfString {
             if (i < 0 or i > g.get_dim()) throw std::out_of_range("Bad axis");
             if (g.get_grid_type(i) == GRID_TYPE_NUMERIC)
               return g.get_numeric_grid(i);
             return g.get_string_grid(i);
           })
      .def("set_grid",
           [](GriddedField& g, Index i, const VectorOrArrayOfString& y) {
             if (i < 0 or i > g.get_dim()) throw std::out_of_range("Bad axis");
             std::visit([&](auto&& z) { g.set_grid(i, z); }, y);
           })
      .def("get_grid_name",
           [](GriddedField& g, Index i) {
             if (i < 0 or i > g.get_dim()) throw std::out_of_range("Bad axis");
             return g.get_grid_name(i);
           })
      .def("set_grid_name",
           [](GriddedField& g, Index i, const String& s) {
             if (i < 0 or i > g.get_dim()) throw std::out_of_range("Bad axis");
             g.set_grid_name(i, s);
           })
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
          })
      .def_property(
          "gridnames",
          [](GriddedField& g) {
            const Index n = g.get_dim();
            ArrayOfString o(n);
            for (Index i = 0; i < n; i++) o[i] = g.get_grid_name(i);
            return o;
          },
          [](GriddedField& g, const ArrayOfString& s) {
            const Index n = g.get_dim();
            ARTS_USER_ERROR_IF(std::size_t(n) not_eq s.size(), "Bad size");
            for (Index i = 0; i < n; i++) {
              g.set_grid_name(i, s[i]);
            }
          })
      .def("checksize_strict", &GriddedField::checksize_strict)
      .def("checksize", &GriddedField::checksize)
      .def(
          "extract_slice",
          [](py::object& g, py::object& s, py::object& i) {
            return details::GriddedField::extract_slice(g, s, i);
          },
          py::arg("slice"),
          py::arg("axis") = 0,
          py::doc(
              R"--(Return a new GriddedField containing a slice of the current one.

Parameters:
    s (slice): Slice.
    axis (int): Axis to slice along.

Returns:
    GriddedField containing sliced grids and data.
)--"))
      .def(
          "to_xarray",
          [](py::object& g) { return details::GriddedField::to_xarray(g); },
          py::doc(
              R"--(Convert GriddedField to xarray.DataArray object.

Convert a GriddedField object into a :func:`xarray.DataArray`
object.  The dataname is used as the DataArray name.

Returns:
    xarray.DataArray object corresponding to gridded field
)--"))
      .def(
          "to_dict",
          [](py::object& g) { return details::GriddedField::to_dict(g); },
          py::doc(
              R"--(Convert GriddedField to dictionary.

Converts a GriddedField object into a classic Python dictionary. The
gridname is used as dictionary key. If the grid is unnamed the key is
generated automatically ('grid1', 'grid2', ...). The data can be
accessed through the 'data' key.

Returns:
    Dictionary containing the grids and data.
)--"))
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
          py::arg("type") = py::str("linear"),
          py::doc(
              R"--(Interpolate GriddedField axis to a new grid.

This function replaces a grid of a GriddField and interpolates all
data to match the new coordinates. :func:`scipy.interpolate.interp1d`
is used for interpolation.

Parameters:
    new_grid (ndarray): The coordinates of the interpolated values.
    axis (int): Specifies the axis of data along which to interpolate.
        Interpolation defaults to the first axis of the GriddedField.
    type (str or function): Rescaling type for function if str or the
        actual rescaling function
    **kwargs:
        Keyword arguments passed to :func:`scipy.interpolate.interp1d`.

Returns: GriddedField
)--"));

  py::class_<GriddedField1, GriddedField>(m, "GriddedField1")
      .def(py::init([]() { return new GriddedField1{}; }))
      .def(py::init([](const String& s) { return new GriddedField1{s}; }))
      .def(py::init([](const std::array<VectorOrArrayOfString, 1>& grids,
                       const Vector& data,
                       const std::optional<std::array<String, 1>>& gridnames,
                       const String& name) {
             constexpr std::size_t n = 1;
             auto* out = new GriddedField1(name);
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
           py::arg_v("name", "", "None"))
      .PythonInterfaceCopyValue(GriddedField1)
      .PythonInterfaceWorkspaceVariableConversion(GriddedField1)
      .PythonInterfaceFileIO(GriddedField1)
      .PythonInterfaceBasicRepresentation(GriddedField1)
      .PythonInterfaceGriddedField(GriddedField1);

  py::class_<GriddedField2, GriddedField>(m, "GriddedField2")
      .def(py::init([]() { return new GriddedField2{}; }))
      .def(py::init([](const String& s) { return new GriddedField2{s}; }))
      .def(py::init([](const std::array<VectorOrArrayOfString, 2>& grids,
                       const Matrix& data,
                       const std::optional<std::array<String, 2>>& gridnames,
                       const String& name) {
             constexpr std::size_t n = 2;
             auto* out = new GriddedField2(name);
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
           py::arg_v("name", "", "None"))
      .PythonInterfaceCopyValue(GriddedField2)
      .PythonInterfaceWorkspaceVariableConversion(GriddedField2)
      .PythonInterfaceFileIO(GriddedField2)
      .PythonInterfaceBasicRepresentation(GriddedField2)
      .PythonInterfaceGriddedField(GriddedField2);

  py::class_<GriddedField3, GriddedField>(m, "GriddedField3")
      .def(py::init([]() { return new GriddedField3{}; }))
      .def(py::init([](const String& s) { return new GriddedField3{s}; }))
      .def(py::init([](const std::array<VectorOrArrayOfString, 3>& grids,
                       const Tensor3& data,
                       const std::optional<std::array<String, 3>>& gridnames,
                       const String& name) {
             constexpr std::size_t n = 3;
             auto* out = new GriddedField3(name);
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
           py::arg_v("name", "", "None"))
      .PythonInterfaceCopyValue(GriddedField3)
      .PythonInterfaceWorkspaceVariableConversion(GriddedField3)
      .PythonInterfaceFileIO(GriddedField3)
      .PythonInterfaceBasicRepresentation(GriddedField3)
      .PythonInterfaceGriddedField(GriddedField3);

  py::class_<GriddedField4, GriddedField>(m, "GriddedField4")
      .def(py::init([]() { return new GriddedField4{}; }))
      .def(py::init([](const String& s) { return new GriddedField4{s}; }))
      .def(py::init([](const std::array<VectorOrArrayOfString, 4>& grids,
                       const Tensor4& data,
                       const std::optional<std::array<String, 4>>& gridnames,
                       const String& name) {
             constexpr std::size_t n = 4;
             auto* out = new GriddedField4(name);
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
           py::arg_v("name", "", "None"))
      .PythonInterfaceCopyValue(GriddedField4)
      .PythonInterfaceWorkspaceVariableConversion(GriddedField4)
      .PythonInterfaceFileIO(GriddedField4)
      .PythonInterfaceBasicRepresentation(GriddedField4)
      .PythonInterfaceGriddedField(GriddedField4);

  py::class_<GriddedField5, GriddedField>(m, "GriddedField5")
      .def(py::init([]() { return new GriddedField5{}; }))
      .def(py::init([](const String& s) { return new GriddedField5{s}; }))
      .def(py::init([](const std::array<VectorOrArrayOfString, 5>& grids,
                       const Tensor5& data,
                       const std::optional<std::array<String, 5>>& gridnames,
                       const String& name) {
             constexpr std::size_t n = 5;
             auto* out = new GriddedField5(name);
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
           py::arg_v("name", "", "None"))
      .PythonInterfaceCopyValue(GriddedField5)
      .PythonInterfaceWorkspaceVariableConversion(GriddedField5)
      .PythonInterfaceFileIO(GriddedField5)
      .PythonInterfaceBasicRepresentation(GriddedField5)
      .PythonInterfaceGriddedField(GriddedField5);

  py::class_<GriddedField6, GriddedField>(m, "GriddedField6")
      .def(py::init([]() { return new GriddedField6{}; }))
      .def(py::init([](const String& s) { return new GriddedField6{s}; }))
      .def(py::init([](const std::array<VectorOrArrayOfString, 6>& grids,
                       const Tensor6& data,
                       const std::optional<std::array<String, 6>>& gridnames,
                       const String& name) {
             constexpr std::size_t n = 6;
             auto* out = new GriddedField6(name);
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
           py::arg_v("name", "", "None"))
      .PythonInterfaceCopyValue(GriddedField6)
      .PythonInterfaceWorkspaceVariableConversion(GriddedField6)
      .PythonInterfaceFileIO(GriddedField6)
      .PythonInterfaceBasicRepresentation(GriddedField6)
      .PythonInterfaceGriddedField(GriddedField6);

  PythonInterfaceWorkspaceArray(GriddedField1);
  PythonInterfaceWorkspaceArray(GriddedField2);
  PythonInterfaceWorkspaceArray(GriddedField3);
  PythonInterfaceWorkspaceArray(GriddedField4);

  PythonInterfaceWorkspaceArray(ArrayOfGriddedField1)
      .def(py::init([](const std::vector<std::vector<GriddedField1>>& x) {
        ArrayOfArrayOfGriddedField1 y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }));
  py::implicitly_convertible<std::vector<std::vector<GriddedField1>>,
                             ArrayOfArrayOfGriddedField1>();

  PythonInterfaceWorkspaceArray(ArrayOfGriddedField2)
      .def(py::init([](const std::vector<std::vector<GriddedField2>>& x) {
        ArrayOfArrayOfGriddedField2 y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }));
  py::implicitly_convertible<std::vector<std::vector<GriddedField2>>,
                             ArrayOfArrayOfGriddedField2>();

  PythonInterfaceWorkspaceArray(ArrayOfGriddedField3)
      .def(py::init([](const std::vector<std::vector<GriddedField3>>& x) {
        ArrayOfArrayOfGriddedField3 y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }));
  py::implicitly_convertible<std::vector<std::vector<GriddedField3>>,
                             ArrayOfArrayOfGriddedField3>();
}
}  // namespace Python