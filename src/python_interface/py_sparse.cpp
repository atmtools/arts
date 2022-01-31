#include <auto_md.h>
#include <pybind11/pybind11.h>
#include <xml_io.h>

#include "debug.h"
#include "matpackII.h"
#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_sparse(py::module_& m) {
  py::class_<Sparse>(m, "Sparse")
      .def(py::init<>())
      .PythonInterfaceFileIO(Sparse)
      .PythonInterfaceBasicRepresentation(Sparse);

  py::enum_<Block::MatrixType>(m, "BlockMatrixType")
      .value("dense", Block::MatrixType::dense)
      .value("sparse", Block::MatrixType::sparse);

  py::class_<Block>(m, "Block")
      .def("get_matrix_type", &Block::get_matrix_type)
      .def_property_readonly("dense",
                             [](Block& x) {
                               ARTS_USER_ERROR_IF(x.get_matrix_type() not_eq
                                                      Block::MatrixType::dense,
                                                  "Block is not dense")
                               return x.get_dense();
                             })
      .def_property_readonly("sparse", [](Block& x) {
        ARTS_USER_ERROR_IF(x.get_matrix_type() not_eq Block::MatrixType::sparse,
                           "Block is not sparse")
        return x.get_sparse();
      });

  py::class_<CovarianceMatrix>(m, "CovarianceMatrix")
      .def(py::init<>())
      .def_property(
          "blocks",
          [](CovarianceMatrix& x) { return x.get_blocks(); },
          [](CovarianceMatrix& x, std::vector<Block> y) {
            x.get_blocks() = std::move(y);
          })
      .PythonInterfaceFileIO(CovarianceMatrix)
      .PythonInterfaceBasicRepresentation(CovarianceMatrix);

  PythonInterfaceWorkspaceArray(Sparse);
}
}  // namespace Python