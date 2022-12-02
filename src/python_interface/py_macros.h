#ifndef py_macros_h
#define py_macros_h

#include <debug.h>
#include <global_data.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xml_io.h>

#include <functional>

#include "python_interface_value_type.h"

namespace Python {
namespace py = pybind11;
}  // namespace Python

constexpr Index negative_clamp(const Index i, const Index n) noexcept {
  return (i < 0) ? i + n : i;
}

#define PythonInterfaceFileIO(Type)                                        \
  def(                                                                     \
      "savexml",                                                           \
      [](const Type& x,                                                    \
         const char* const file,                                           \
         const char* const type,                                           \
         bool clobber) {                                                   \
        xml_write_to_file(                                                 \
            file, x, string2filetype(type), clobber ? 0 : 1, Verbosity()); \
      },                                                                   \
      py::arg("file").none(false),                                         \
      py::arg("type").none(false) = "ascii",                               \
      py::arg("clobber") = true,                                           \
      py::doc("Saves :class:`" #Type "` to file\n"                         \
              "\n"                                                         \
              "Parameters:\n"                                              \
              "    file (str): The path to which the file is written."     \
              " Note that several of the options might modify the"         \
              " name or write more files\n"                                \
              "    type (str): Type of file to save (ascii. zascii,"       \
              " or binary)\n"                                              \
              "    clobber (bool): Overwrite existing files or add new"    \
              " file with modified name?\n"                                \
              "\n"                                                         \
              "On Error:\n"                                                \
              "    Throws RuntimeError for any failure to save"))          \
      .def(                                                                \
          "readxml",                                                       \
          [](Type& x, const char* const file) {                            \
            xml_read_from_file(file, x, Verbosity());                      \
          },                                                               \
          py::arg("file").none(false),                                     \
          py::doc("Read :class:`" #Type "` from file\n"                    \
                  "\n"                                                     \
                  "Parameters:\n"                                          \
                  "    file (str): A file that can be read\n"              \
                  "\n"                                                     \
                  "On Error:\n"                                            \
                  "    Throws RuntimeError for any failure to read"))

#define PythonInterfaceIndexItemAccess(Type)                                \
  def("__len__", [](const Type& x) { return x.nelem(); })                   \
      .def(                                                                 \
          "__getitem__",                                                    \
          [](Type& x, Index i) -> decltype(x[i])& {                         \
            i = negative_clamp(i, x.nelem());                               \
            if (x.nelem() <= i or i < 0)                                    \
              throw std::out_of_range(var_string("Bad index access: ",      \
                                                 i,                         \
                                                 " in object of size [0, ", \
                                                 x.size(),                  \
                                                 ")"));                     \
            return as_ref(x[i]);                                            \
          },                                                                \
          py::return_value_policy::reference_internal)                      \
      .def(                                                                 \
          "__setitem__",                                                    \
          [](Type& x, Index i, decltype(x[i]) y) {                          \
            i = negative_clamp(i, x.nelem());                               \
            if (x.nelem() <= i or i < 0)                                    \
              throw std::out_of_range(var_string("Bad index access: ",      \
                                                 i,                         \
                                                 " in object of size [0, ", \
                                                 x.size(),                  \
                                                 ")"));                     \
            x[i] = std::move(y);                                            \
          },                                                                \
          py::return_value_policy::reference_internal)

#define PythonInterfaceBasicRepresentation(Type)       \
  def(                                                 \
      "__str__",                                       \
      [](const Type& x) { return var_string(x); },     \
      py::is_operator())                               \
      .def(                                            \
          "__repr__",                                  \
          [](const Type&) { return #Type; }, \
          py::is_operator())

#define PythonInterfaceCopyValue(Type)                                  \
  def(                                                                  \
      "__copy__", [](Type& t) -> Type { return t; }, py::is_operator()) \
      .def(                                                             \
          "__deepcopy__",                                               \
          [](Type& t, py::dict&) -> Type { return t; },                 \
          py::is_operator())

#define PythonInterfaceArrayEquality(Type)                              \
  def(                                                                  \
      "__eq__",                                                         \
      [](const Type& a, const Type& b) {                                \
        std::size_t n = a.size();                                       \
        if (n not_eq b.size()) return false;                            \
        for (size_t i = 0; i < n; i++) {                                \
          auto v1 = py::cast(a[i]);                                     \
          auto v2 = py::cast(b[i]);                                     \
          auto ev = v1.attr("__eq__")(v2);                              \
          if (py::hasattr(ev, "all")) ev = py::bool_(ev.attr("all")()); \
          if (not py::isinstance<py::bool_>(ev)) return false;          \
          if (not ev.cast<bool>()) return false;                        \
        }                                                               \
        return true;                                                    \
      },                                                                \
      py::is_operator())

/** Gives the basic Array<X> interface of Arts
 * 
 * The type should be the full name of the class
 *
 * Python operators provided includes:
 *
 * x = Array<X>()
 * x = Array<X>(INDEX)
 * x = Array<X>(INDEX, X())
 * x.append(X())
 * y = x.pop()
 * list(x)
 * x[INDEX]
 * x[INDEX] = X()
 * len(x)
 * for a in x: DO_SOMETHING(a)
 */
#define PythonInterfaceArrayDefault(BaseType)                                  \
  PythonInterfaceIndexItemAccess(Array<BaseType>)                              \
      .PythonInterfaceArrayEquality(Array<BaseType>)                           \
      .PythonInterfaceCopyValue(Array<BaseType>)                               \
      .def(py::init([]() { return new Array<BaseType>{}; }))                   \
      .def(py::init([](Index n, const BaseType& v) {                           \
             return new Array<BaseType>(n, v);                                 \
           }),                                                                 \
           py::arg("size"),                                                    \
           py::arg("cval"))                                                    \
      .def(py::init([](const std::vector<BaseType>& v) {                       \
             return new Array<BaseType>{v};                                    \
           }),                                                                 \
           py::arg("arr"))                                                     \
      .def(                                                                    \
          "append",                                                            \
          [](Array<BaseType>& x, BaseType y) {                                 \
            x.emplace_back(std::move(y));                                      \
          },                                                                   \
          py::doc("Appends a :class:`" #BaseType "` at the end of the Array")) \
      .def(                                                                    \
          "pop",                                                               \
          [](Array<BaseType>& x, Index i) -> BaseType {                        \
            if (x.size() < 1) throw std::out_of_range("pop from empty list");  \
            i = negative_clamp(i, x.nelem());                                  \
            if (x.nelem() <= i or i < 0)                                       \
              throw std::out_of_range(var_string("pop index out of range"));   \
            std::rotate(x.begin() + i, x.begin() + 1 + i, x.end());            \
            BaseType copy = std::move(x.back());                               \
            x.pop_back();                                                      \
            return copy;                                                       \
          },                                                                   \
          py::doc("Pops a :class:`" #BaseType "` from the Array"),             \
          py::arg("i") = -1)                                                   \
      .def("reverse",                                                          \
           [](Array<BaseType>& x) { std::reverse(x.begin(), x.end()); })       \
      .def("clear", [](Array<BaseType>& x) { x.clear(); })                     \
      .def("copy", [](Array<BaseType>& x) { return x; })                       \
      .def("extend",                                                           \
           [](Array<BaseType>& x, const Array<BaseType>& y) {                  \
             for (auto& z : y) x.push_back(z);                                 \
           })                                                                  \
      .def("insert",                                                           \
           [](Array<BaseType>& x, Index i, BaseType y) {                       \
             auto pos = (i < x.nelem() and i >= 0) ? x.begin() + i : x.end();  \
             x.insert(pos, std::move(y));                                      \
           })                                                                  \
      .def(py::pickle(                                                         \
          [](const Array<BaseType>& v) {                                       \
            auto n = v.size();                                                 \
            std::vector<BaseType> out(n);                                      \
            std::copy(v.begin(), v.end(), out.begin());                        \
            return py::make_tuple(std::move(out));                             \
          },                                                                   \
          [](const py::tuple& t) {                                             \
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")                \
            return new Array<BaseType>{t[0].cast<std::vector<BaseType>>()};    \
          }))

#define PythonInterfaceCommonMath(Type) \
  def(py::self + Type())                \
      .def(py::self - Type())           \
      .def(py::self * Type())           \
      .def(py::self / Type())           \
      .def(py::self += Type())          \
      .def(py::self -= Type())          \
      .def(py::self *= Type())          \
      .def(py::self /= Type())          \
      .def(py::self == Type())          \
      .def(py::self != Type())          \
      .def(py::self <= Type())          \
      .def(py::self < Type())           \
      .def(py::self >= Type())          \
      .def(py::self > Type())           \
      .def(Type() + py::self)           \
      .def(Type() - py::self)           \
      .def(Type() * py::self)           \
      .def(Type() / py::self)           \
      .def(Type() == py::self)          \
      .def(Type() != py::self)          \
      .def(Type() <= py::self)          \
      .def(Type() < py::self)           \
      .def(Type() >= py::self)          \
      .def(Type() > py::self)

#define PythonInterfaceCommonMathSelf \
  def(+py::self)                      \
      .def(-py::self)                 \
      .def(py::self + py::self)       \
      .def(py::self - py::self)       \
      .def(py::self * py::self)       \
      .def(py::self / py::self)       \
      .def(py::self += py::self)      \
      .def(py::self -= py::self)      \
      .def(py::self *= py::self)      \
      .def(py::self /= py::self)      \
      .def(py::self == py::self)      \
      .def(py::self != py::self)      \
      .def(py::self <= py::self)      \
      .def(py::self < py::self)       \
      .def(py::self >= py::self)      \
      .def(py::self > py::self)

#define PythonInterfaceWorkspaceVariableConversion(Type)                   \
  def(py::init([](const Type& x) { return new Type{x}; }), py::arg("val")) \
      .def(py::init([](const WorkspaceVariable* w) {                       \
             Type& v = *w;                                                 \
             return new Type(v);                                           \
           }),                                                             \
           py::arg("wsv").none(false).noconvert(),                         \
           py::doc("Automatic conversion from a workspace variable"))

//! Place at the end!
#define PythonInterfaceWorkspaceDocumentation(Type)                      \
  doc() = global_data::wsv_groups.at(global_data::WsvGroupMap.at(#Type)) \
              .desc.c_str()

//! Place at the end!  Add a few line-breaks to fit in extra documentation!
#define PythonInterfaceWorkspaceDocumentationExtra(Type, extra)              \
  doc() =                                                                    \
      (global_data::wsv_groups.at(global_data::WsvGroupMap.at(#Type)).desc + \
       extra)                                                                \
          .c_str()

/*! The workspace array interface

It takes the base type as input, e.g., ArrayOfIndex (Array<Index>) takes
"Index" as input.  This allows setting up the default conversion from
std::vector<BaseType> to be created.

The main limitation for this is that the type must be names ArrayOfBaseType
in Arts.  Without this limitation, the representation and automatic documentation
that is generated by this macro would show the C++ type name rather than the
desired python name.  "ArrayOfBaseType" is the class exposed to python
*/
#define PythonInterfaceWorkspaceArray(BaseType)                                    \
  auto auto_impl_name##BaseType =                                                  \
      py::class_<ArrayOf##BaseType>(m, "ArrayOf" #BaseType);                       \
  auto_impl_name##BaseType.PythonInterfaceWorkspaceDocumentationExtra(             \
      ArrayOf##BaseType,                                                           \
      "\n\nThis class is partly compatible with python list of :class:`" #BaseType \
      "`");                                                                        \
                                                                                   \
  py::implicitly_convertible<std::vector<BaseType>, ArrayOf##BaseType>();          \
                                                                                   \
  auto_impl_name##BaseType.PythonInterfaceFileIO(ArrayOf##BaseType)                \
      .def(                                                                        \
          "__repr__",                                                              \
          [](const py::object& arr) {                                              \
            return py::list{arr}.attr("__repr__")();                               \
          },                                                                       \
          py::is_operator())                                                       \
      .def(                                                                        \
          "__str__",                                                               \
          [](const py::object& arr) {                                              \
            return py::list{arr}.attr("__str__")();                                \
          },                                                                       \
          py::is_operator())                                                       \
      .PythonInterfaceArrayDefault(BaseType)                                       \
      .PythonInterfaceWorkspaceVariableConversion(ArrayOf##BaseType)

#define PythonInterfaceGriddedField(Type)                                       \
  def_readwrite("data", &Type::data)                                            \
      .def(py::pickle(                                                          \
          [](const Type& self) {                                                \
            Index n = self.get_dim();                                           \
            std::vector<std::variant<Vector, ArrayOfString>> outgrid;           \
            ArrayOfString outname;                                              \
                                                                                \
            for (Index i = 0; i < n; i++) {                                     \
              outname.emplace_back(self.get_grid_name(i));                      \
              if (self.get_grid_type(i) == GRID_TYPE_NUMERIC)                   \
                outgrid.emplace_back(self.get_numeric_grid(i));                 \
              else                                                              \
                outgrid.emplace_back(self.get_string_grid(i));                  \
            }                                                                   \
                                                                                \
            return py::make_tuple(                                              \
                self.get_name(), self.data, outname, outgrid);                  \
          },                                                                    \
          [](const py::tuple& t) {                                              \
            ARTS_USER_ERROR_IF(t.size() != 4, "Invalid state!")                 \
                                                                                \
            auto name = t[0].cast<String>();                                    \
            auto data = t[1].cast<decltype(Type::data)>();                      \
            auto grid_names = t[2].cast<ArrayOfString>();                       \
            auto grids =                                                        \
                t[3].cast<std::vector<std::variant<Vector, ArrayOfString>>>();  \
                                                                                \
            auto* out = new Type{name};                                         \
            out->data = data;                                                   \
            for (Index i = 0; i < out->get_dim(); i++) {                        \
              out->set_grid_name(i, grid_names[i]);                             \
              std::visit([&](auto&& v) { out->set_grid(i, v); }, grids[i]);     \
            }                                                                   \
                                                                                \
            return out;                                                         \
          }))                                                                   \
      .def(                                                                     \
          "__eq__",                                                             \
          [](Type& g1, Type& g2) {                                              \
            const auto n = g1.get_dim();                                        \
                                                                                \
            if (not g1.checksize()) return false;                               \
            if (not g2.checksize()) return false;                               \
                                                                                \
            for (Index i = 0; i < n; i++) {                                     \
              if (g1.get_grid_size(i) not_eq g2.get_grid_size(i)) {             \
                return false;                                                   \
              }                                                                 \
                                                                                \
              if (g1.get_grid_type(i) not_eq g2.get_grid_type(i)) {             \
                return false;                                                   \
              }                                                                 \
                                                                                \
              if (g1.get_grid_type(i) == GRID_TYPE_NUMERIC) {                   \
                auto o1 = py::cast(g1.get_numeric_grid(i));                     \
                auto o2 = py::cast(g2.get_numeric_grid(i));                     \
                if (not py::bool_((o1.attr("__eq__")(o2)).attr("all")())        \
                            .cast<bool>()) {                                    \
                  return false;                                                 \
                }                                                               \
              } else {                                                          \
                auto o1 = py::cast(g1.get_string_grid(i));                      \
                auto o2 = py::cast(g2.get_string_grid(i));                      \
                if (not o1.attr("__eq__")(o2).cast<bool>()) return false;       \
              }                                                                 \
                                                                                \
              if (g1.get_grid_name(i) not_eq g2.get_grid_name(i)) {             \
                return false;                                                   \
              }                                                                 \
            }                                                                   \
            if (g1.get_name() not_eq g2.get_name()) return false;               \
                                                                                \
            auto o1 = py::cast(g1);                                             \
            auto o2 = py::cast(g2);                                             \
            if (not py::bool_(                                                  \
                        (o1.attr("data").attr("__eq__")(o2.attr("data")))       \
                            .attr("all")())                                     \
                        .cast<bool>()) {                                        \
              return false;                                                     \
            }                                                                   \
            return true;                                                        \
          },                                                                    \
          py::is_operator())                                                    \
      .def_property_readonly(                                                   \
          "shape",                                                              \
          [](py::object& g) {                                                   \
            return g.attr("data").attr("value").attr("shape");                  \
          })                                                                    \
      .def(                                                                     \
          "__getitem__",                                                        \
          [](py::object& g, py::object& i) {                                    \
            return g.attr("data").attr("value").attr("__getitem__")(i);         \
          },                                                                    \
          py::is_operator())                                                    \
      .def(                                                                     \
          "__setitem__",                                                        \
          [](py::object& g, py::object& i, py::object& n) {                     \
            g.attr("data").attr("value").attr("__setitem__")(i, n);             \
          },                                                                    \
          py::is_operator())                                                    \
      .def(                                                                     \
          "get",                                                                \
          [](Type& g,                                                           \
             const String& key,                                                 \
             const py::object& def,                                             \
             bool keep_dims) -> py::object {                                    \
            ARTS_USER_ERROR_IF(                                                 \
                g.get_grid_type(0) == GRID_TYPE_NUMERIC,                        \
                "Method only works, if the first grid is an \"ArrayOfString\"") \
            auto& grid = g.get_string_grid(0);                                  \
            auto ptr = std::find(grid.begin(), grid.end(), key);                \
            if (ptr == grid.end()) return def;                                  \
                                                                                \
            auto o = py::cast(g);                                               \
            auto v = o.attr("data").attr("value").attr("__getitem__")(          \
                std::array<Index, 1>{std::distance(grid.begin(), ptr)});        \
                                                                                \
            if (keep_dims) return v;                                            \
            return v.attr("squeeze")();                                         \
          },                                                                    \
          py::arg("key"),                                                       \
          py::arg("default") = py::none(),                                      \
          py::arg("keep_dims") = true)                                          \
      .def(                                                                     \
          "set",                                                                \
          [](Type& g, const String& key, const py::object& data) {              \
            ARTS_USER_ERROR_IF(                                                 \
                g.get_grid_type(0) == GRID_TYPE_NUMERIC,                        \
                "Method only works, if the first grid is an \"ArrayOfString\"") \
            auto& grid = g.get_string_grid(0);                                  \
            auto ptr = std::find(grid.begin(), grid.end(), key);                \
            ARTS_USER_ERROR_IF(ptr == grid.end(), "Cannot find: ", key)         \
                                                                                \
            auto o = py::cast(g);                                               \
            auto v = o.attr("data").attr("value").attr("__setitem__")(          \
                std::array<Index, 1>{std::distance(grid.begin(), ptr)}, data);  \
          },                                                                    \
          py::arg("key"),                                                       \
          py::arg("data"))                                                      \
      .def(                                                                     \
          "scale",                                                              \
          [](py::object& g, const String& key, const py::object& fac) {         \
            g.attr("set")(key, (g.attr("get")(key)).attr("__mul__")(fac));      \
          },                                                                    \
          py::arg("key"),                                                       \
          py::arg("factor"))                                                    \
      .def_static(                                                              \
          "from_xarray",                                                        \
          [](py::object& v) {                                                   \
            auto g = py::type::of<Type>();                                      \
            return details::GriddedField::from_xarray(g, v);                    \
          },                                                                    \
          py::doc(                                                              \
              R"--(Create GriddedField from a xarray.DataArray object.\
\
The data and its dimensions are returned as a :class:`GriddedField` object.\
The DataArray name is used as name for the gridded field. If the attribute\
`data_name` is present, it is used as `dataname` on the :class:`GriddedField`.\
\
Parameters:\
    da (xarray.DataArray): xarray.DataArray containing the dimensions and data.\
\
Returns:\
    GriddedField object.\
)--"))

#define PythonInterfaceReadWriteData(Type, data) \
  def_readwrite(#data, &Type::data, py::return_value_policy::reference_internal)

#define PythonInterfaceBasicReferenceProperty(                               \
    Type, PropertyName, ReadFunction, WriteFunction)                         \
  def_property(                                                              \
      #PropertyName,                                                         \
      py::cpp_function(                                                      \
          [](Type& x)                                                        \
              -> std::remove_cv_t<                                           \
                  std::add_lvalue_reference_t<decltype(x.ReadFunction())>> { \
            return x.ReadFunction();                                         \
          },                                                                 \
          py::return_value_policy::reference_internal),                      \
      [](Type& x, std::remove_reference_t<decltype(x.ReadFunction())> y) {   \
        x.WriteFunction() = std::move(y);                                    \
      })

#define PythonInterfaceSelfAttribute(ATTR) \
  def_property_readonly(                   \
      #ATTR, [](py::object& x) { return x.attr("value").attr(#ATTR); })

#define PythonInterfaceSelfOperator(ATTR)                          \
  def(                                                             \
      #ATTR,                                                       \
      [](py::object& x) { return x.attr("value").attr(#ATTR)(); }, \
      py::is_operator())

#define PythonInterfaceValueOperator(ATTR)     \
  def(                                         \
      #ATTR,                                   \
      [](py::object& x, py::object& y) {       \
        return x.attr("value").attr(#ATTR)(y); \
      },                                       \
      py::is_operator())

#define PythonInterfaceTwoValueOperator(ATTR)           \
  def(                                                  \
      #ATTR,                                            \
      [](py::object& x, py::object& y, py::object& z) { \
        return x.attr("value").attr(#ATTR)(y, z);       \
      },                                                \
      py::is_operator())

#define PythonInterfaceInternalValueOperator(ATTR) \
  def(                                             \
      #ATTR,                                       \
      [](py::object& x, py::object& y) {           \
        x.attr("value").attr(#ATTR)(y);            \
        return x;                                  \
      },                                           \
      py::is_operator())

#define PythonInterfaceInternalTwoValueOperator(ATTR)   \
  def(                                                  \
      #ATTR,                                            \
      [](py::object& x, py::object& y, py::object& z) { \
        x.attr("value").attr(#ATTR)(y, z);              \
        return x;                                       \
      },                                                \
      py::is_operator())

#define PythonInterfaceNumpyValueProperties \
  PythonInterfaceSelfAttribute(ndim)        \
      .PythonInterfaceSelfAttribute(shape)  \
      .PythonInterfaceSelfAttribute(size)   \
      .PythonInterfaceSelfAttribute(dtype)

#define PythonInterfaceValueOperators                      \
  PythonInterfaceSelfOperator(__len__)                     \
      .PythonInterfaceSelfOperator(__iter__)               \
      .PythonInterfaceSelfOperator(__pos__)                \
      .PythonInterfaceSelfOperator(__neg__)                \
      .PythonInterfaceSelfOperator(__abs__)                \
      .PythonInterfaceSelfOperator(__invert__)             \
      .PythonInterfaceSelfOperator(__complex__)            \
      .PythonInterfaceSelfOperator(__int__)                \
      .PythonInterfaceSelfOperator(__float__)              \
                                                           \
      .PythonInterfaceValueOperator(__add__)               \
      .PythonInterfaceValueOperator(__sub__)               \
      .PythonInterfaceValueOperator(__mul__)               \
      .PythonInterfaceValueOperator(__matmul__)            \
      .PythonInterfaceValueOperator(__truediv__)           \
      .PythonInterfaceValueOperator(__floordiv__)          \
      .PythonInterfaceValueOperator(__mod__)               \
      .PythonInterfaceValueOperator(__divmod__)            \
      .PythonInterfaceValueOperator(__pow__)               \
      .PythonInterfaceTwoValueOperator(__pow__)            \
                                                           \
      .PythonInterfaceValueOperator(__radd__)              \
      .PythonInterfaceValueOperator(__rsub__)              \
      .PythonInterfaceValueOperator(__rmul__)              \
      .PythonInterfaceValueOperator(__rmatmul__)           \
      .PythonInterfaceValueOperator(__rtruediv__)          \
      .PythonInterfaceValueOperator(__rfloordiv__)         \
      .PythonInterfaceValueOperator(__rmod__)              \
      .PythonInterfaceValueOperator(__rdivmod__)           \
      .PythonInterfaceValueOperator(__rpow__)              \
      .PythonInterfaceTwoValueOperator(__rpow__)           \
                                                           \
      .PythonInterfaceInternalValueOperator(__iadd__)      \
      .PythonInterfaceInternalValueOperator(__isub__)      \
      .PythonInterfaceInternalValueOperator(__imul__)      \
      .PythonInterfaceInternalValueOperator(__imatmul__)   \
      .PythonInterfaceInternalValueOperator(__itruediv__)  \
      .PythonInterfaceInternalValueOperator(__ifloordiv__) \
      .PythonInterfaceInternalValueOperator(__imod__)      \
      .PythonInterfaceInternalValueOperator(__idivmod__)   \
      .PythonInterfaceInternalValueOperator(__ipow__)      \
      .PythonInterfaceInternalTwoValueOperator(__ipow__)   \
                                                           \
      .PythonInterfaceValueOperator(__lt__)                \
      .PythonInterfaceValueOperator(__le__)                \
      .PythonInterfaceValueOperator(__eq__)                \
      .PythonInterfaceValueOperator(__ne__)                \
      .PythonInterfaceValueOperator(__ge__)                \
      .PythonInterfaceValueOperator(__gt__)                \
                                                           \
      .PythonInterfaceValueOperator(__contains__)          \
                                                           \
      .PythonInterfaceValueOperator(__getitem__)           \
      .PythonInterfaceTwoValueOperator(__setitem__)        \
                                                           \
      .PythonInterfaceSelfOperator(flatten)

#endif
