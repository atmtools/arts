#include <nanobind/eigen/dense.h>
#include <nanobind/eigen/sparse.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/variant.h>

#include <memory>
#include <utility>
#include <variant>

#include "covariance_matrix.h"
#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "nanobind/nanobind.h"
#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_sparse(py::module_& m) try {
  py::class_<Sparse> sp(m, "Sparse");
  workspace_group_interface(sp);
  sp.def(
        "__init__",
        [](Sparse* s, Eigen::SparseMatrix<Numeric, Eigen::RowMajor> es) {
          new (s) Sparse{};
          s->matrix.swap(es);
        },
        "From :class:`scipy.sparse.csr_matrix`")
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
          py::rv_policy::reference_internal)
      .def("__setitem__",
           [](Sparse& x, std::tuple<Index, Index> ind, Numeric y) {
             const auto [r, c] = ind;
             if (r < 0 or r >= x.nrows())
               throw std::out_of_range(var_string("row ", r));
             if (c < 0 or c >= x.ncols())
               throw std::out_of_range(var_string("col ", c));
             x.rw(r, c) = y;
           })
      .def_prop_rw(
          "value",
          [](Sparse& s) { return s.matrix; },
          [](Sparse& s, Eigen::SparseMatrix<Numeric, Eigen::RowMajor> ns) {
            s.matrix.swap(ns);
          },
          "Treat as :class:`scipy.sparse.csr_matrix`")
      .def(
          "toarray",
          [](Sparse& sp) { return Eigen::MatrixXd(sp.matrix); },
          R"(Make a dense :class:`numpy.ndarray` from the sparse matrix

Returns
-------
arr : numpy.ndarray
    A dense array
)")
      .def("__getstate__",
           [](const Sparse& self) {
             return std::tuple<Eigen::SparseMatrix<Numeric, Eigen::RowMajor>>{
                 self.matrix};
           })
      .def("__setstate__",
           [](Sparse* self,
              const std::tuple<Eigen::SparseMatrix<Numeric, Eigen::RowMajor>>&
                  state) {
             new (self) Sparse{};
             self->matrix = std::get<0>(state);
           });
  py::implicitly_convertible<Eigen::SparseMatrix<Numeric, Eigen::RowMajor>,
                             Sparse>();

  auto a1 = py::bind_vector<ArrayOfSparse, py::rv_policy::reference_internal>(
      m, "ArrayOfSparse");
  workspace_group_interface(a1);
  vector_interface(a1);

  py::class_<Block>(m, "Block")
      .def(py::init<Range, Range, IndexPair, std::shared_ptr<Matrix>>(),
           "By value, dense")
      .def(py::init<Range, Range, IndexPair, std::shared_ptr<Sparse>>(),
           "By value, sparse")
      .PythonInterfaceCopyValue(Block)
      .def_prop_rw(
          "matrix",
          [](Block& x) -> std::variant<Matrix*, Sparse*> {
            if (x.is_dense()) return &x.get_dense();
            return &x.get_sparse();
          },
          [](Block& x, std::variant<Matrix*, Sparse*> y) {
            if (std::holds_alternative<Matrix*>(y)) {
              x.set_matrix(
                  std::make_shared<Matrix>(**std::get_if<Matrix*>(&y)));
            } else {
              x.set_matrix(
                  std::make_shared<Sparse>(**std::get_if<Sparse*>(&y)));
            }
          },
          ":class:`~pyarts.arts.Matrix` or :class:`~pyarts.arts.Sparse` The matrix held inside the instance")
      .def("__getstate__",
           [](const Block& self) {
             if (self.is_sparse())
               return std::
                   tuple<Range, Range, IndexPair, std::variant<Matrix, Sparse>>{
                       self.get_row_range(),
                       self.get_column_range(),
                       self.get_indices(),
                       self.get_sparse()};
             return std::
                 tuple<Range, Range, IndexPair, std::variant<Matrix, Sparse>>{
                     self.get_row_range(),
                     self.get_column_range(),
                     self.get_indices(),
                     self.get_dense()};
           })
      .def(
          "__setstate__",
          [](Block* self,
             const std::
                 tuple<Range, Range, IndexPair, std::variant<Matrix, Sparse>>&
                     state) {
            auto row_range    = std::get<0>(state);
            auto column_range = std::get<1>(state);
            auto indices      = std::get<2>(state);
            auto mat          = std::get<3>(state);

            if (std::holds_alternative<Sparse>(mat))
              new (self) Block(row_range,
                               column_range,
                               indices,
                               std::make_shared<Sparse>(std::get<Sparse>(mat)));
            else
              new (self) Block(row_range,
                               column_range,
                               indices,
                               std::make_shared<Matrix>(std::get<Matrix>(mat)));
          })
      .doc() = "A single block matrix";

  py::class_<BlockMatrix> bm(m, "BlockMatrix");
  bm.def(py::init_implicit<Matrix>());
  bm.def(py::init_implicit<Sparse>());
  bm.def(
      "__init__",
      [](BlockMatrix* s, Eigen::SparseMatrix<Numeric, Eigen::RowMajor> es) {
        new (s) BlockMatrix{};
        *s = py::cast<Sparse>(py::type<Sparse>()(es));
      },
      "From :class:`scipy.sparse.csr_matrix`");
  py::implicitly_convertible<Eigen::SparseMatrix<Numeric, Eigen::RowMajor>,
                             BlockMatrix>();
  bm.def(
      "__init__",
      [](BlockMatrix* v,
         const py::ndarray<py::numpy, Numeric, py::ndim<2>, py::c_contig>& a) {
        new (v) BlockMatrix{};
        *v = py::cast<Matrix>(py::type<Matrix>()(a));
      },
      "a"_a);
  py::implicitly_convertible<py::ndarray<py::numpy, Numeric, py::ndim<2>, py::c_contig>,
                             BlockMatrix>();
  bm.def_prop_rw(
      "matrix",
      [](BlockMatrix& bm) -> std::variant<Matrix, Sparse> {
        if (bm.not_null()) {
          if (bm.is_dense()) return bm.dense();
          return bm.sparse();
        }

        return Matrix{};
      },
      [](BlockMatrix& bm, const std::variant<Matrix, Sparse>& mat) {
        std::visit([&bm](auto& m) { bm = m; }, mat);
      });
  bm.def(
      "__array__",
      [](py::object& v, py::object dtype, py::object copy) {
        return v.attr("matrix").attr("__array__")(dtype, copy);
      },
      "dtype"_a.none() = py::none(),
      "copy"_a.none()  = py::none());
  bm.def_prop_rw(
      "value",
      [](py::object& x) { return x.attr("__array__")("copy"_a = false); },
      [](BlockMatrix& a, const std::variant<Matrix, Sparse>& b) {
        std::visit([&a](auto& c) { a = c; }, b);
      });
  common_ndarray(bm);
  workspace_group_interface(bm);

  py::class_<CovarianceMatrix> covm(m, "CovarianceMatrix");
  workspace_group_interface(covm);
  covm.def_prop_rw(
          "blocks",
          [](CovarianceMatrix& x) { return x.get_blocks(); },
          [](CovarianceMatrix& x, std::vector<Block> y) {
            x.get_blocks() = std::move(y);
          },
          ":class:`list` of :class:`~pyarts.arts.Block`")
      .def("__getstate__",
           [](CovarianceMatrix& self) {
             return std::tuple<std::vector<Block>, std::vector<Block>>(
                 self.get_blocks(), self.get_inverse_blocks());
           })
      .def("__setstate__",
           [](CovarianceMatrix* self,
              const std::tuple<std::vector<Block>, std::vector<Block>>& state) {
             new (self) CovarianceMatrix{};
             self->get_blocks()         = std::get<0>(state);
             self->get_inverse_blocks() = std::get<1>(state);
           });
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize sparse\n", e.what()));
}
}  // namespace Python
