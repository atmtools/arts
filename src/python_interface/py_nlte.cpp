#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>
#include <workspace.h>

namespace Python {
namespace py = nanobind;

void py_nlte(py::module_ &m) try {
  py::class_<VibrationalEnergyLevels>(m, "VibrationalEnergyLevels")
      .def("__init__",
           [](VibrationalEnergyLevels *n,
              std::map<QuantumIdentifier, Numeric> &inv) {
             new (n) VibrationalEnergyLevels{};
             for (auto &x : inv) {
               n->operator[](x.first) = x.second;
             }
           })
      .def(
          "__getitem__",
          [](VibrationalEnergyLevels &x, const QuantumIdentifier &q) {
            if (x.data.find(q) == x.end())
              throw py::key_error(var_string(q).c_str());
            return x[q];
          },
          py::rv_policy::reference_internal)
      .def("__setitem__",
           [](VibrationalEnergyLevels &x,
              const QuantumIdentifier &q,
              Numeric y) { x[q] = y; })
      .def("__getstate__",
           [](const VibrationalEnergyLevels &t) {
             std::vector<QuantumIdentifier> qn;
             std::vector<Numeric> v;

             qn.reserve(t.size());
             v.reserve(t.size());

             for (auto &x : t) {
               qn.emplace_back(x.first);
               v.emplace_back(x.second);
             }

             return std::tuple<std::vector<QuantumIdentifier>,
                               std::vector<Numeric>>(qn, v);
           })
      .def("__setstate__",
           [](VibrationalEnergyLevels &t,
              const std::tuple<std::vector<QuantumIdentifier>,
                               std::vector<Numeric>> &s) {
             const auto qn = std::get<0>(s);
             const auto v  = std::get<1>(s);
             ARTS_USER_ERROR_IF(v.size() != qn.size(), "Invalid size!")

             for (std::size_t i = 0; i < v.size(); i++) {
               t[qn[i]] = v[i];
             }
           });
} catch (std::exception &e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize nlte\n", e.what()));
}
}  // namespace Python
