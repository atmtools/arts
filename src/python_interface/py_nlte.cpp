
#include <py_auto_wsg_init.h>
#include <python_interface.h>

#include <memory>
#include <vector>

#include "debug.h"
#include "py_macros.h"
#include "python_interface_value_type.h"
#include "quantum_numbers.h"

namespace Python {
void py_nlte(py::module_ &m) try {
  py_staticVibrationalEnergyLevels(m)
      .def(py::init([](std::map<QuantumIdentifier, Numeric> &in) {
        auto out = std::make_shared<VibrationalEnergyLevels>();
        for (auto &x : in) {
          out->operator[](x.first) = x.second;
        }
        return out;
      }))
      .def(
          "__getitem__",
          [](VibrationalEnergyLevels &x, const QuantumIdentifier &q) {
            if (x.data.find(q) == x.end()) throw py::key_error(var_string(q));
            return x[q];
          },
          py::return_value_policy::reference_internal)
      .def("__setitem__",
           [](VibrationalEnergyLevels &x,
              const QuantumIdentifier &q,
              Numeric y) { x[q] = y; })
      .def(py::pickle(
          [](const VibrationalEnergyLevels &t) {
            std::vector<QuantumIdentifier> qn;
            std::vector<Numeric> v;

            qn.reserve(t.size());
            v.reserve(t.size());

            for (auto &x : t) {
              qn.emplace_back(x.first);
              v.emplace_back(x.second);
            }

            return py::make_tuple(qn, v);
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 2, "Invalid state!")

            const auto qn = t[0].cast<std::vector<QuantumIdentifier>>();
            const auto v = t[1].cast<std::vector<Numeric>>();
            ARTS_USER_ERROR_IF(v.size() != qn.size(), "Invalid size!")

            auto out = std::make_shared<VibrationalEnergyLevels>();
            for (std::size_t i = 0; i < v.size(); i++) {
              out->operator[](qn[i]) = v[i];
            }
            return out;
          }));
} catch (std::exception &e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize nlte\n", e.what()));
}
}  // namespace Python
