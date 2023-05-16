#include <pybind11/pybind11.h>

#include <iostream>

void check_pyversion() {
  std::cout << "PyMajorVersion: " << PY_MAJOR_VERSION << "\n";
  std::cout << "PyMinorVersion: " << PY_MINOR_VERSION << "\n";
  std::cout << "PyMicroVersion: " << PY_MICRO_VERSION << "\n";
  std::cout << "Py_GetVersion()[4]: " << Py_GetVersion()[4] << "\n";
}

PYBIND11_MODULE(test_pyversion, m) {
  m.def("check_pyversion", &check_pyversion, "Output detected Python version");
}
