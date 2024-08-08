#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
              << "<NUM VARIABLE FILES> <NUM METHOD FILES>\n";
    return EXIT_FAILURE;
  }

  const int num_variables = std::stoi(argv[1]);
  const int num_methods   = std::stoi(argv[2]);

  std::ofstream os("py_auto_interface.cpp");
  os << "#include <python_interface.h>\n\n";
  os << "namespace Python {\n";

  for (int i = 0; i < num_variables; i++) {
    os << "void py_auto_wsv_" << i << "(py::class_<Workspace>& ws);\n";
  }
  os << "void py_auto_wsv(py::class_<Workspace>& ws) {\n";
  for (int i = 0; i < num_variables; i++) {
    os << "  py_auto_wsv_" << i << "(ws);\n";
  }
  os << "}\n\n";

  for (int i = 0; i < num_methods; i++) {
    os << "void py_auto_wsm_" << i << "(py::class_<Workspace>& ws);\n";
  }
  os << "void py_auto_wsm(py::class_<Workspace>& ws) {\n";
  for (int i = 0; i < num_methods; i++) {
    os << "  py_auto_wsm_" << i << "(ws);\n";
  }
  os << "}\n}  // namespace Python\n";

  return EXIT_SUCCESS;
}
