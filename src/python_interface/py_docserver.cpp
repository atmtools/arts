#include <docserver.h>
#include <pybind11/pybind11.h>

namespace Python {
namespace py = pybind11;

void py_docserver(py::module_& m) {
  auto docserver = m.def_submodule("docserver");
  docserver.doc() = R"--(Run the ARTS documentation server)--";
  docserver.def("run",
                &run_docserver,
                py::arg("port") = 9000,
                py::arg("baseurl") = "",
                py::arg("daemon") = false,
                R"--(Run the ARTS documentation server

The ARTS documentation server provides a web interface to the ARTS
documentation. It can be used to browse the documentation of variables,
methods and agendas.

By default, the server is started on port 9000 and is accessible via a
web browser at `http://localhost:9000`.

Parameters:
    port (int): Port to run the server on
    baseurl (str): Base URL to use for the server
    daemon (bool): Run the server as a daemon
)--");
}
}  // namespace Python

