#include <py_auto_interface.h>
#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include <memory>
#include <utility>
#include <variant>

#include "covariance_matrix.h"
#include "matpack.h"
#include "py_macros.h"

namespace Python {
void py_sparse(py::module_& m) {
  py::class_<Sparse>(m, "Sparse")
      .def(py::init([]() { return new Sparse{}; }))
      .def(py::init([](Eigen::SparseMatrix<Numeric, Eigen::RowMajor> es) {
        Sparse s;
        s.matrix.swap(es);
        return s;
      }))
      .PythonInterfaceCopyValue(Sparse)
      .PythonInterfaceWorkspaceVariableConversion(Sparse)
      .PythonInterfaceFileIO(Sparse)
      .PythonInterfaceBasicRepresentation(Sparse)
      .def(
          "__getitem__",
          [](Sparse& x, std::tuple<Index, Index> ind) -> Numeric& {
            const auto [r, c] = ind;
            if (r < 0 or r >= x.nrows())
              throw std::out_of_range(var_string("row ", r));
            if (c < 0 or c >= x.ncols())
              throw std::out_of_range(var_string("col ", c));
            return x.rw(r, c);
          },
          py::return_value_policy::reference_internal)
      .def("__setitem__",
           [](Sparse& x, std::tuple<Index, Index> ind, Numeric y) {
             const auto [r, c] = ind;
             if (r < 0 or r >= x.nrows())
               throw std::out_of_range(var_string("row ", r));
             if (c < 0 or c >= x.ncols())
               throw std::out_of_range(var_string("col ", c));
             x.rw(r, c) = y;
           })
      .def_property(
          "value",
          [](Sparse& s) { return s.matrix; },
          [](Sparse& s, Eigen::SparseMatrix<Numeric, Eigen::RowMajor> ns) {
            s.matrix.swap(ns);
          })
      .def("toarray", [](Sparse& sp) { return Eigen::MatrixXd(sp.matrix); })
      .def(py::pickle(
          [](const Sparse& self) { return py::make_tuple(self.matrix); },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")

            auto* out = new Sparse{};
            out->matrix = t[0].cast<decltype(out->matrix)>();

            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(Sparse);
  py::implicitly_convertible<Eigen::SparseMatrix<Numeric, Eigen::RowMajor>,
                             Sparse>();

  py::enum_<Block::MatrixType>(m, "BlockMatrixType")
      .value("dense", Block::MatrixType::dense)
      .value("sparse", Block::MatrixType::sparse)
      .PythonInterfaceCopyValue(Block::MatrixType)
      .def(py::pickle(
          [](const Block::MatrixType& self) {
            return py::make_tuple(static_cast<Index>(self));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")

            return static_cast<Block::MatrixType>(t[0].cast<Index>());
          }));

  py::class_<Block>(m, "Block")
      .def(py::init([](Range row_range,
                       Range column_range,
                       IndexPair indices,
                       Matrix mat) {
        return new Block(row_range,
                         column_range,
                         indices,
                         std::make_shared<Matrix>(std::move(mat)));
      }))
      .def(py::init([](Range row_range,
                       Range column_range,
                       IndexPair indices,
                       Sparse mat) {
        return new Block(row_range,
                         column_range,
                         indices,
                         std::make_shared<Sparse>(std::move(mat)));
      }))
      .PythonInterfaceCopyValue(Block)
      .def_property("matrix",
                    py::cpp_function(
                        [](Block& x) -> std::variant<Matrix*, Sparse*> {
                          if (x.get_matrix_type() == Block::MatrixType::dense)
                            return &x.get_dense();
                          return &x.get_sparse();
                        },
                        py::return_value_policy::reference_internal),
                    [](Block& x, std::variant<Matrix*, Sparse*> y) {
                      if (std::holds_alternative<Matrix*>(y)) {
                        x.set_matrix(std::make_shared<Matrix>(
                            **std::get_if<Matrix*>(&y)));
                      } else {
                        x.set_matrix(std::make_shared<Sparse>(
                            **std::get_if<Sparse*>(&y)));
                      }
                    })
      .def(py::pickle(
          [](const Block& self) {
            if (self.get_matrix_type() == Block::MatrixType::sparse)
              return py::make_tuple(self.get_row_range(),
                                    self.get_column_range(),
                                    self.get_indices(),
                                    self.get_matrix_type(),
                                    self.get_sparse());
            return py::make_tuple(self.get_row_range(),
                                  self.get_column_range(),
                                  self.get_indices(),
                                  self.get_matrix_type(),
                                  self.get_dense());
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 5, "Invalid state!")

            auto row_range = t[0].cast<Range>();
            auto column_range = t[1].cast<Range>();
            auto indices = t[2].cast<IndexPair>();

            if (t[3].cast<Block::MatrixType>() == Block::MatrixType::sparse)
              return new Block{
                  row_range,
                  column_range,
                  indices,
                  std::make_shared<Sparse>(Sparse{t[4].cast<Sparse>()})};
            return new Block{
                row_range,
                column_range,
                indices,
                std::make_shared<Matrix>(Matrix{t[4].cast<Matrix>()})};
          }));

  py::class_<CovarianceMatrix>(m, "CovarianceMatrix")
      .def(py::init([]() { return new CovarianceMatrix{}; }))
      .PythonInterfaceCopyValue(CovarianceMatrix)
      .PythonInterfaceWorkspaceVariableConversion(CovarianceMatrix)
      .def_property(
          "blocks",
          [](CovarianceMatrix& x) { return x.get_blocks(); },
          [](CovarianceMatrix& x, std::vector<Block> y) {
            x.get_blocks() = std::move(y);
          })
      .PythonInterfaceFileIO(CovarianceMatrix)
      .PythonInterfaceBasicRepresentation(CovarianceMatrix)
      .def(py::pickle(
          [](CovarianceMatrix& self) {
            return py::make_tuple(self.get_blocks(), self.get_inverse_blocks());
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 2, "Invalid state!")

            auto b = t[0].cast<std::vector<Block>>();
            auto i = t[0].cast<std::vector<Block>>();

            auto* out = new CovarianceMatrix{};
            out->get_blocks() = std::move(b);
            out->get_inverse_blocks() = std::move(i);
            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(CovarianceMatrix);

  PythonInterfaceWorkspaceArray(Sparse);
}
}  // namespace Python
