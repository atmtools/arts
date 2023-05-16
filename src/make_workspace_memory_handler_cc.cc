#include <iostream>
#include <stdexcept>
#include "array.h"
#include "arts.h"
#include "file.h"
#include "global_data.h"
#include "matpack_data.h"

int main() {
  try {
    define_wsv_groups();
    using global_data::wsv_groups;

    ofstream ofs;
    open_output_file(ofs, "workspace_memory_handler.cc");

    ofs << "/*! \\file  workspace_memory_handler.cc\n"
        << " *\n"
        << " *  \\brief Defines global workspace_handler_objects and its \n"
        << " *   dispatch functions.\n\n"
        << " *   <b>DO NOT EDIT!</b>\n\n"
        << " *   \\date " << __DATE__ << ", " << __TIME__ << "\n"
        << " */\n\n";

    ofs << "#include \"workspace_memory_handler.h\"\n"
        << "#include <iostream>\n"
        << "#include \"matpack_sparse.h\"\n"
        << "#include \"m_general.h\"\n"
        << "#include \"supergeneric.h\"\n"
        << "#include \"artstime.h\"\n"
        << "#include \"ppath.h\"\n"
        << "#include \"gas_abs_lookup.h\"\n"
        << "#include \"linemixing_hitran.h\"\n"
        << "#include \"optproperties.h\"\n"
        << "#include \"gridded_fields.h\"\n"
        << "#include \"jacobian.h\"\n"
        << "#include \"agenda_class.h\"\n"
        << "#include \"mc_interp.h\"\n"
        << "#include \"mc_antenna.h\"\n"
        << "#include \"cia.h\"\n"
        << "#include \"propagationmatrix.h\"\n"
        << "#include \"transmissionmatrix.h\"\n"
        << "#include \"covariance_matrix.h\"\n"
        << "#include \"sun.h\"\n"
        << "#include \"telsem.h\"\n"
        << "#include \"tessem.h\"\n"
        << "#include \"xsec_fit.h\"\n"
        << "#include \"absorptionlines.h\"\n"
        << "#include \"linemixing.h\"\n"
        << "#include \"callback.h\"\n"
        << "\n";

    ofs << "// Allocation and deallocation routines for workspace groups\n";
    for (Index i = 0; i < wsv_groups.nelem(); ++i) {
      ofs << "std::shared_ptr<void> allocate_wsvg_" << wsv_groups[i] << "(){\n"
          << "  return std::make_shared<" << wsv_groups[i] << ">();\n}\n"
          << "std::shared_ptr<void> duplicate_wsvg_" << wsv_groups[i]
          << "(const std::shared_ptr<void>& vp) {"
          << "  return std::make_shared<" << wsv_groups[i] << ">(*static_cast<"
          << wsv_groups[i] << "*>(vp.get()));\n}\n\n";
    }

    ofs << "  /// Initialization dispatch functions.\n"
        << "void WorkspaceMemoryHandler::initialize() {\n"
        << "  allocation_ptrs_.resize(" << wsv_groups.size() << ");\n"
        << "  duplication_ptrs_.resize(" << wsv_groups.size() << ");\n\n";

    for (Index i = 0; i < wsv_groups.nelem(); ++i) {
      ofs << "  allocation_ptrs_[" << i << "] = allocate_wsvg_" << wsv_groups[i]
          << ";\n"
          << "  duplication_ptrs_[" << i << "] = duplicate_wsvg_" << wsv_groups[i]
          << ";\n";
    }
    ofs << "}\n";
  } catch (const std::runtime_error &x) {
    cout << "Something went wrong. Message text:\n";
    cout << x.what() << '\n';
    return 1;
  }

  return 0;
}
