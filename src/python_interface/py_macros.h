#ifndef py_macros_h
#define py_macros_h

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
      py::arg("file"),                                                     \
      py::arg("type") = "ascii",                                           \
      py::arg("clobber") = true)                                           \
      .def("loadxml", [](Type& x, const char* const file) {                \
        xml_read_from_file(file, x, Verbosity());                          \
      })

#define PythonInterfaceIndexItemAccess(Type)                                 \
  def("__len__", [](const Type& x) { return x.nelem(); })                    \
      .def("__getitem__",                                                    \
           [](const Type& x, Index i) {                                      \
             if (x.nelem() <= i or i < 0)                                    \
               throw std::out_of_range(var_string("Bad index access: ",      \
                                                  i,                         \
                                                  " in object of size [0, ", \
                                                  x.size(),                  \
                                                  ")"));                     \
             return x[i];                                                    \
           })                                                                \
      .def("__setitem__", [](Type& x, Index i, decltype(x[i]) y) {           \
             if (x.nelem() <= i or i < 0)                                    \
               throw std::out_of_range(var_string("Bad index access: ",      \
                                                  i,                         \
                                                  " in object of size [0, ", \
                                                  x.size(),                  \
                                                  ")"));                     \
        x[i] = std::move(y);                                                 \
      })

#define PythonInterfaceBasicRepresentation(Type) \
  def("__str__", [](const Type& x) {             \
    return var_string(x);                        \
  }).def("__repr__", [](const Type& x) { return var_string(x); })

#define PythonInterfaceBasicIteration(Type)                          \
  def(                                                               \
      "__iter__",                                                    \
      [](Type& x) { return py::make_iterator(x.begin(), x.end()); }, \
      py::keep_alive<0, 1>())

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
 * list(x)
 * x[INDEX]
 * len(x)
 * for a in x: DO_SOMETHING(a)
 */
#define PythonInterfaceArrayDefault(Type)          \
  PythonInterfaceIndexItemAccess(Type)             \
      .PythonInterfaceBasicIteration(Type)         \
      .def(py::init<>())                           \
      .def(py::init<Index>())                      \
      .def(py::init<Index, decltype(Type{}[0])>()) \
      .def("append",                               \
           [](Type& x, decltype(x[0]) y) { x.emplace_back(std::move(y)); })

/** Provides -=, +=, /=. and *= for all LHS Type */
#define PythonInterfaceInPlaceMathOperators(Type, Other)  \
  def(                                                    \
      "__imul__",                                         \
      [](Type& x, const Other& y) { return x *= y; },     \
      py::is_operator(),                                  \
      py::keep_alive<0, 1>())                             \
      .def(                                               \
          "__itruediv__",                                 \
          [](Type& x, const Other& y) { return x /= y; }, \
          py::is_operator(),                              \
          py::keep_alive<0, 1>())                         \
      .def(                                               \
          "__iadd__",                                     \
          [](Type& x, const Other& y) { return x += y; }, \
          py::is_operator(),                              \
          py::keep_alive<0, 1>())                         \
      .def(                                               \
          "__isub__",                                     \
          [](Type& x, const Other& y) { return x -= y; }, \
          py::is_operator(),                              \
          py::keep_alive<0, 1>())

/** Provides -, +, /. and * for Type (right and left)
 */
#define PythonInterfaceMathOperators(Type, Other)        \
  def(                                                   \
      "__mul__",                                         \
      [](Type& x, const Other& y) { return x * y; },     \
      py::is_operator())                                 \
      .def(                                              \
          "__truediv__",                                 \
          [](Type& x, const Other& y) { return x / y; }, \
          py::is_operator())                             \
      .def(                                              \
          "__add__",                                     \
          [](Type& x, const Other& y) { return x + y; }, \
          py::is_operator())                             \
      .def(                                              \
          "__sub__",                                     \
          [](Type& x, const Other& y) { return x - y; }, \
          py::is_operator())                             \
      .def(                                              \
          "__rmul__",                                    \
          [](Type& x, const Other& y) { return y * x; }, \
          py::is_operator())                             \
      .def(                                              \
          "__rtruediv__",                                \
          [](Type& x, const Other& y) { return y / x; }, \
          py::is_operator())                             \
      .def(                                              \
          "__radd__",                                    \
          [](Type& x, const Other& y) { return y + x; }, \
          py::is_operator())                             \
      .def(                                              \
          "__rsub__",                                    \
          [](Type& x, const Other& y) { return y - x; }, \
          py::is_operator())

#define PythonInterfaceWorkspaceArray(Type) \
  py::class_<Type>(m, #Type)                \
      .PythonInterfaceFileIO(Type)          \
      .PythonInterfaceArrayDefault(Type)    \
      .PythonInterfaceBasicRepresentation(Type)

#define PythonInterfaceGriddedField(Type)                                \
  def_readwrite("data", &Type::data)                                     \
      .def("get_dim", [](Type& gf) { return gf.get_dim(); })             \
      .def("get_name", [](Type& gf) { return gf.get_name(); })           \
      .def("set_name",                                                   \
           [](Type& gf, const String& str) { return gf.set_name(str); }) \
      .def("get_grid_name",                                              \
           [](Type& gf, Index i) {                                       \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.get_grid_name(i);                                 \
           })                                                            \
      .def("get_grid_size",                                              \
           [](Type& gf, Index i) {                                       \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.get_grid_size(i);                                 \
           })                                                            \
      .def("get_grid_type",                                              \
           [](Type& gf, Index i) {                                       \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.get_grid_type(i);                                 \
           })                                                            \
      .def("get_numeric_grid",                                           \
           [](Type& gf, Index i) {                                       \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.get_numeric_grid(i);                              \
           })                                                            \
      .def("get_string_grid",                                            \
           [](Type& gf, Index i) {                                       \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.get_string_grid(i);                               \
           })                                                            \
      .def("set_grid_name",                                              \
           [](Type& gf, Index i, const String& str) {                    \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.set_grid_name(i, str);                            \
           })                                                            \
      .def("set_grid",                                                   \
           [](Type& gf, Index i, const Vector& vec) {                    \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.set_grid(i, vec);                                 \
           })                                                            \
      .def("set_grid", [](Type& gf, Index i, const ArrayOfString& str) { \
        ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")           \
        return gf.set_grid(i, str);                                      \
      });

#endif