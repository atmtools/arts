/* Copyright (C) 2020 Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA.
*/

#include <iostream>
#include <stdexcept>
#include "array.h"
#include "arts.h"
#include "file.h"
#include "global_data.h"
#include "matpackI.h"

int main() {
  try {
    define_wsv_group_names();
    using global_data::wsv_group_names;

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
        << "#include \"matpackII.h\"\n"
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
        << "#include \"telsem.h\"\n"
        << "#include \"tessem.h\"\n"
        << "#include \"hitran_xsec.h\"\n"
        << "#include \"absorptionlines.h\"\n"
        << "\n";

    ofs << "// Allocation and deallocation routines for workspace groups\n";
    for (Index i = 0; i < wsv_group_names.nelem(); ++i) {
      ofs << "void *allocate_wsvg_" << wsv_group_names[i] << "(){\n"
          << "  return (void *)new " << wsv_group_names[i] << ";\n}\n"
          << "void deallocate_wsvg_" << wsv_group_names[i]
          << "(void *vp)\n"
          << "  { delete (" << wsv_group_names[i] << " *)vp;\n}\n"
          << "void *duplicate_wsvg_" << wsv_group_names[i]
          << "(void *vp) {"
          << "  return (new " << wsv_group_names[i] << "(*("
          << wsv_group_names[i] << " *)vp));\n}\n\n";
    }

    ofs << "  /// Initialization dispatch functions.\n"
        << "void WorkspaceMemoryHandler::initialize() {\n"
        << "  allocation_ptrs_.resize(" << wsv_group_names.size() << ");\n"
        << "  deallocation_ptrs_.resize(" << wsv_group_names.size() << ");\n"
        << "  duplication_ptrs_.resize(" << wsv_group_names.size() << ");\n\n";

    for (Index i = 0; i < wsv_group_names.nelem(); ++i) {
      ofs << "  allocation_ptrs_[" << i << "] = allocate_wsvg_" << wsv_group_names[i]
          << ";\n"
          << "  deallocation_ptrs_[" << i << "] = deallocate_wsvg_" << wsv_group_names[i]
          << ";\n"
          << "  duplication_ptrs_[" << i << "] = duplicate_wsvg_" << wsv_group_names[i]
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
