#include <arts_omp.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <obsel.h>
#include <python_interface.h>

#include <algorithm>
#include <memory>
#include <stdexcept>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "rtepack.h"

namespace Python {
void py_sensor(py::module_& m) try {
  py::class_<SensorPosLos> splos(m, "SensorPosLos");
  generic_interface(splos);
  common_ndarray(splos);
  splos
      .def(
          "__array__",
          [](SensorPosLos& x, py::object dtype, py::object copy)
              -> std::variant<
                  py::ndarray<py::numpy, Numeric, py::shape<5>, py::c_contig>,
                  py::object> {
            std::array<size_t, 1> shape = {5};
            auto np                     = py::module_::import_("numpy");
            auto w =
                py::ndarray<py::numpy, Numeric, py::shape<5>, py::c_contig>(
                    &x, 1, shape.data(), py::cast(&x));

            if (not dtype.is_none()) {
              return np.attr("asarray")(w, "dtype"_a = dtype, "copy"_a = copy);
            }

            return w.cast((copy.is_none() or not py::bool_(copy))
                              ? py::rv_policy::automatic_reference
                              : py::rv_policy::copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none(),
          "Returns a :class:`~numpy.ndarray` of the object.")
      .def_prop_rw(
          "value",
          [](py::object& x) { return x.attr("__array__")(); },
          [](SensorPosLos& a, const SensorPosLos& b) { a = b; },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`")
      .def(py::init<Vector3, Vector2>(), "From pos and los")
      .def_rw("pos", &SensorPosLos::pos, "Position\n\n.. :class:`Vector3`")
      .def_rw(
          "los", &SensorPosLos::los, "Line of sight\n\n.. :class:`Vector2`");

  py::class_<SensorPosLosVector> vsplos(m, "SensorPosLosVector");
  generic_interface(vsplos);
  common_ndarray(vsplos);
  vsplos.def(
      "sort",
      [](SensorPosLosVector& v, const SensorKeyType x, bool reverse) {
        using enum SensorKeyType;

        switch (x) {
          case alt: stdr::sort(v, {}, &SensorPosLos::alt); break;
          case lat: stdr::sort(v, {}, &SensorPosLos::lat); break;
          case lon: stdr::sort(v, {}, &SensorPosLos::lon); break;
          case za:  stdr::sort(v, {}, &SensorPosLos::za); break;
          case aa:  stdr::sort(v, {}, &SensorPosLos::aa); break;
          case f:   break;
        }

        if (reverse) {
          if (f == x) throw std::runtime_error("Cannot reverse frequency");
          stdr::reverse(v);
        }
      },
      "key"_a,
      "reverse"_a = false,
      "Sorts the vector by the given key.");
  vsplos
      .def(
          "__array__",
          [](SensorPosLosVector& x,
             py::object dtype,
             py::object copy) -> std::variant<py::ndarray<py::numpy,
                                                          Numeric,
                                                          py::shape<-1, 5>,
                                                          py::c_contig>,
                                              py::object> {
            std::array<size_t, 2> shape = {static_cast<size_t>(x.size()), 5};
            auto np                     = py::module_::import_("numpy");
            auto w =
                py::ndarray<py::numpy, Numeric, py::shape<-1, 5>, py::c_contig>(
                    x.data_handle(), 2, shape.data(), py::cast(&x));

            if (not dtype.is_none()) {
              return np.attr("asarray")(w, "dtype"_a = dtype, "copy"_a = copy);
            }

            return w.cast((copy.is_none() or not py::bool_(copy))
                              ? py::rv_policy::automatic_reference
                              : py::rv_policy::copy);
          },
          "dtype"_a.none() = py::none(),
          "copy"_a.none()  = py::none(),
          "Returns a :class:`~numpy.ndarray` of the object.")
      .def_prop_rw(
          "value",
          [](py::object& x) { return x.attr("__array__")(); },
          [](SensorPosLosVector& x, Matrix& y) {
            if (y.ncols() != 5) {
              throw std::runtime_error("Bad shape");
            }
            x.resize(y.nrows());
            std::transform(y.begin(), y.end(), x.begin(), [](const auto& row) {
              return SensorPosLos{.pos = {row[0], row[1], row[2]},
                                  .los = {row[3], row[4]}};
            });
          },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`")
      .def("__getstate__",
           [](const py::object& self) {
             return py::make_tuple(self.attr("value"));
           })
      .def("__setstate__", [](py::object& self, const py::tuple& state) {
        self.attr("value") = state[0];
      });

  py::class_<sensor::SparseStokvecMatrix> ssm(m, "SparseStokvecMatrix");
  ssm.doc() = "A sparse matrix of Stokvec";
  generic_interface(ssm);
  ssm.def("__getitem__",
          [](const sensor::SparseStokvecMatrix& self, std::tuple<Size, Size> ij)
              -> Stokvec { return self[std::get<0>(ij), std::get<1>(ij)]; });
  ssm.def("__getitem__",
          [](const sensor::SparseStokvecMatrix& self,
             std::tuple<Size, Size, Size> ijk) -> Numeric {
            return self[std::get<0>(ijk), std::get<1>(ijk)][std::get<2>(ijk)];
          });
  ssm.def("__setitem__",
          [](sensor::SparseStokvecMatrix& self,
             std::tuple<Size, Size> ij,
             Stokvec v) { self[std::get<0>(ij), std::get<1>(ij)] = v; });
  ssm.def("__setitem__",
          [](sensor::SparseStokvecMatrix& self,
             std::tuple<Size, Size, Size> ijk,
             Numeric v) {
            self[std::get<0>(ijk), std::get<1>(ijk)][std::get<2>(ijk)] = v;
          });
  ssm.def_prop_ro(
      "dense",
      [](const sensor::SparseStokvecMatrix& self) {
        return StokvecMatrix{self};
      },
      "Get a dense copy.\n\n.. :class:`~pyarts3.arts.StokvecMatrix`");
  ssm.def_prop_ro(
      "value",
      [](const py::object& self) {
        return self.attr("dense").attr("__array__")("copy"_a = true);
      },
      "Allow using numpy on a copy of the class.\n\n.. :class:`~numpy.ndarray`");
  ssm.def(
      "__array__",
      [](const py::object& self,
         const py::object& dtype,
         const py::object& copy) {
        if (not copy.is_none() and not py::bool_(copy))
          throw std::invalid_argument("Must copy dense array");

        return self.attr("value").attr("__array__")("dtype"_a = dtype,
                                                    "copy"_a  = true);
      },
      "dtype"_a = py::none(),
      "copy"_a  = py::none(),
      "Allows :func:`~numpy.array` to be called with the object.");
  common_ndarray(ssm);
  ssm.def(
      "reduce",
      [](const sensor::SparseStokvecMatrix& self,
         Stokvec v,
         bool along_poslos,
         bool along_freq) -> std::variant<Matrix, Vector, Numeric> {
        if (along_poslos and along_freq) {
          Numeric sum = 0.0;
          for (const auto& elem : self.vector()) {
            sum += dot(elem.data, v);
          }
          return sum;
        }

        if (along_freq) {
          Vector out(self.nrows(), 0.0);
          for (const auto& elem : self.vector()) {
            out[elem.irow] += dot(elem.data, v);
          }
          return out;
        }

        if (along_poslos) {
          Vector out(self.ncols(), 0.0);
          for (const auto& elem : self.vector()) {
            out[elem.icol] += dot(elem.data, v);
          }
          return out;
        }

        Matrix out(self.nrows(), self.ncols(), 0.0);
        for (const auto& elem : self.vector()) {
          out[elem.irow, elem.icol] = dot(elem.data, v);
        }
        return out;
      },
      "pol"_a          = Stokvec{1., 0., 0., 0.},
      "along_poslos"_a = false,
      "along_freq"_a   = false,
      R"(Reduces the matrix by polarization.

The reduction happens either along the position/line-of-sight axis,
the frequency axis, or both.  The reduction is the equivalent of performing
``np.einsum("ijk,k->ij", weight_matrix, pol)``

where ``i`` is the position/line-of-sight index, ``j`` is the frequency index,
and ``k`` is the Stokes index if both ``along_poslos`` and ``along_freq``
are ``False``.  If both are ``True``, the result is a single number.
If only one of them is ``True``, the result is a vector the size
of the remaining axis.

Parameters
----------
pol : Stokvec
    The polarization to reduce with.  Default is unpolarized, i.e., ``[1, 0, 0, 0]``.
along_poslos : bool
    If ``True``, reduce along the position/line-of-sight axis.  Default is ``False``.
along_freq : bool
    If ``True``, reduce along the frequency axis.  Default is ``False``.

Returns
-------
Numeric, Vector, or Matrix
    The reduced value, depending on the parameters.
)");

  py::class_<SensorObsel> so(m, "SensorObsel");
  generic_interface(so);
  so.def(py::init<const AscendingGrid&,
                  const SensorPosLosVector&,
                  StokvecMatrix>())
      .def_prop_ro("f_grid",
                   &SensorObsel::f_grid,
                   "Frequency grid\n\n.. :class:`AscendingGrid`")
      .def_prop_ro(
          "weight_matrix",
          [](const SensorObsel& self) { return self.weight_matrix(); },
          "Weights matrix\n\n.. :class:`SparseStokvecMatrix`")
      .def_prop_ro(
          "poslos",
          &SensorObsel::poslos_grid,
          "Position and line of sight grid\n\n.. :class:`SensorPosLosVector`");

  auto a0 = py::bind_vector<Array<SensorPosLosVector>,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfSensorPosLosVector");
  a0.doc() = "Array of SensorPosLosVector";
  generic_interface(a0);
  vector_interface(a0);

  auto a1 =
      py::bind_vector<ArrayOfSensorObsel, py::rv_policy::reference_internal>(
          m, "ArrayOfSensorObsel");
  generic_interface(a1);
  vector_interface(a1);

  a1.def(
      "unique_freq_grids",
      [](const ArrayOfSensorObsel& x) {
        const SensorSimulations simuls = collect_simulations(x);
        std::vector<AscendingGrid> out;
        out.reserve(simuls.size());
        for (auto& k : simuls | std::views::keys) out.push_back(*k);
        return out;
      },
      "List of the unique frequency grids");

  a1.def(
      "unique_poslos_grids",
      [](const ArrayOfSensorObsel& x) {
        const SensorSimulations simuls = collect_simulations(x);
        std::unordered_set<std::shared_ptr<const SensorPosLosVector>> vecs;
        for (auto& k : simuls | std::views::values) {
          vecs.insert(k.begin(), k.end());
        }

        std::vector<SensorPosLosVector> out;
        out.reserve(vecs.size());
        for (auto& v : vecs) out.push_back(*v);
        return out;
      },
      "List of the unique poslos grids");

  a1.def(
      "collect_freq_grids",
      [](ArrayOfSensorObsel& x) {
        if (x.empty()) return;

        const auto tmp = [&x]() {
          std::set<std::shared_ptr<const AscendingGrid>> freq_set;
          for (const auto& obsel : x) freq_set.insert(obsel.f_grid_ptr());
          if (freq_set.size() == 1)
            return std::make_pair(true, *freq_set.begin());

          std::vector<Numeric> fs;
          fs.reserve(std::accumulate(
              freq_set.begin(),
              freq_set.end(),
              Size(0),
              [](Size a, const auto& b) { return a + b->size(); }));
          for (auto& freqs : freq_set) fs.append_range(*freqs);

          stdr::sort(fs);
          fs.erase(std::unique(fs.begin(), fs.end()), fs.end());
          return std::make_pair(
              false, std::make_shared<const AscendingGrid>(std::move(fs)));
        }();

        //! Clang bug workaround
        const auto& early    = tmp.first;
        const auto& freq_ptr = tmp.second;

        if (early) return;

#pragma omp parallel for if (not arts_omp_in_parallel())
        for (Size i = 0; i < x.size(); i++) {
          sensor::SparseStokvecMatrix w(x[i].poslos_grid().size(),
                                        freq_ptr->size());
          for (const auto& we : x[i].weight_matrix()) {
            const Index ifreq = stdr::distance(
                freq_ptr->begin(),
                stdr::lower_bound(*freq_ptr, x[i].f_grid()[we.icol]));
            w[we.irow, ifreq] = we.data;
          }

          x[i] = {freq_ptr, x[i].poslos_grid_ptr(), std::move(w)};
        }
      },
      R"(Collect all frequency grids in a single grid.

.. tip::
      This generally speeds up radiative transfer calculations.
      The internal weighting matrix is updated internally to
      give the same forward model result.  This matrix is stored
      sparsely, so memory usage should not be a big issue.
      Beware that accessing the weight matrix creates a dense matrix copy,
      which might be large if many frequency grids are combined.

.. note::
      If you need to consider sensor error corrections in frequency,
      such as frequency shifts, stretches or standing waves,
      this method might not have the desired effect.  All channels are considered
      as part of the same "physical" spectrometer, i.e., a standing wave in one channel
      will stretch to all the rest.
)");

  a1.def(
      "collect_poslos_grids",
      [](ArrayOfSensorObsel& x) {
        std::set<std::shared_ptr<const SensorPosLosVector>> poslos_set;
        for (const auto& obsel : x) poslos_set.insert(obsel.poslos_grid_ptr());

        const auto all_poslos_ptr = [&poslos_set]() {
          SensorPosLosVector all_poslos;
          for (const auto& pl : poslos_set) {
            for (auto& poslos : *pl) {
              if (stdr::find(all_poslos, poslos) == all_poslos.end()) {
                all_poslos.push_back(poslos);
              }
            }
          }
          return std::make_shared<const SensorPosLosVector>(
              std::move(all_poslos));
        }();

#pragma omp parallel for if (not arts_omp_in_parallel())
        for (Size i = 0; i < x.size(); i++) {
          sensor::SparseStokvecMatrix w(all_poslos_ptr->size(),
                                        x[i].f_grid().size());

          for (const auto& we : x[i].weight_matrix()) {
            const Index iposlos = stdr::distance(
                all_poslos_ptr->begin(),
                stdr::find(*all_poslos_ptr, x[i].poslos_grid()[we.irow]));
            w[iposlos, we.icol] = we.data;
          }

          x[i] = {x[i].f_grid_ptr(), all_poslos_ptr, std::move(w)};
        }
      },
      R"(Collect all poslos grids in a single grid.

.. tip::
      This generally speeds up radiative transfer calculations.
      The internal weighting matrix is updated internally to
      give the same forward model result.  This matrix is stored
      sparsely, so memory usage should not be a big issue.
      Beware that accessing the weight matrix creates a dense matrix copy,
      which might be large if many poslos grids are combined.

.. note::
      If you need to consider sensor error corrections in line-of-sight or position,
      this method might not have the desired effect.  All channels are considered
      as part of the same "physical" sensor, i.e., a correction in position or line-of-sight
      of one channel will affect all channels.
)");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize sensors\n{}", e.what()));
}
}  // namespace Python
