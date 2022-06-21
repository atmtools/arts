/* Copyright (C) 2000-2012 Stefan Buehler <sbuehler@ltu.se>

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
   USA. */

#include "agenda_record.h"
#include "array.h"
#include "arts.h"
#include "file.h"
#include "global_data.h"
#include "methods.h"
#include "workspace_ng.h"

/* Adds commas and indentation to parameter lists. */
void align(ofstream& ofs, bool& is_first_parameter, const String& indent) {
  // Add comma and line break, if not first element:
  if (is_first_parameter)
    is_first_parameter = false;
  else {
    ofs << ",\n";
    // Make proper indentation:
    ofs << indent;
  }
}

int main() {
  try {
    // Make the global data visible:
    using global_data::md_data;
    using global_data::wsv_groups;
    const Array<WsvRecord>& wsv_data = Workspace::wsv_data;

    // Initialize the wsv group name array:
    define_wsv_groups();

    // Initialize wsv data.
    Workspace::define_wsv_data();

    // Initialize WsvMap.
    Workspace::define_wsv_map();

    // Initialize method data.
    define_md_data_raw();

    // Expand supergeneric methods:
    expand_md_data_raw_to_md_data();

    const Index n_md = md_data.nelem();

    // Write auto_md.cc:
    // -----------
    ofstream ofs;
    open_output_file(ofs, "auto_md.cc");

    ofs << "// This file was generated automatically by make_auto_md_cc.cc.\n";
    ofs << "// DO NOT EDIT !\n";
    ofs << "// Generated: " << __DATE__ << ", " << __TIME__ << "\n\n";

    ofs << "#include \"arts.h\"\n"
        << "#include \"auto_md.h\"\n"
        << "#include \"wsv_aux.h\"\n"
        << "#include \"mc_interp.h\"\n"
        << "#include \"m_append.h\"\n"
        << "#include \"m_delete.h\"\n"
        << "#include \"m_copy.h\"\n"
        << "#include \"m_conversion.h\"\n"
        << "#include \"m_extract.h\"\n"
        << "#include \"m_general.h\"\n"
        << "#include \"m_gridded_fields.h\"\n"
        << "#include \"m_ignore.h\"\n"
        << "#include \"m_nc.h\"\n"
        << "#include \"m_reduce.h\"\n"
        << "#include \"m_select.h\"\n"
        << "#include \"m_xml.h\"\n"
        << "#include \"m_basic_types.h\"\n"
        << "#include \"propagationmatrix.h\"\n"
        << "#include \"transmissionmatrix.h\"\n"
        << "#include \"agenda_record.h\"\n"
        << "#include \"workspace_ng.h\"\n"
        << "#include \"global_data.h\"\n"
        << "#include \"absorptionlines.h\"\n"
        << "\n";

    //ofs << "static Index agendacallcount = 0;\n";

    // Write all get-away functions:
    // -----------------------------
    for (Index i = 0; i < n_md; ++i) {
      const MdRecord& mdd = md_data[i];

      // This is needed to flag the first function parameter, which
      // needs no line break before being written:
      bool is_first_parameter = true;
      // The String indent is needed to achieve the correct
      // indentation of the functin parameters:
      String indent = String(mdd.Name().nelem() + 3, ' ');
      ;
      // Flag to pass the workspace to the WSM. Only true if the WSM has
      // an Agenda as input.
      bool pass_workspace = false;

      // There are four lists of parameters that we have to
      // write.
      ArrayOfIndex vo = mdd.Out();            // Output
      const ArrayOfIndex& vi = mdd.InOnly();  // Input
      ArrayOfIndex vgo = mdd.GOutType();      // Generic Output
      ArrayOfIndex vgi = mdd.GInType();       // Generic Input
      // vo and vi contain handles of workspace variables,
      // vgo and vgi handles of workspace variable groups.

      // There used to be a similar block here for the generic
      // input/output variables. However, this was a mistake. For
      // example, if a method has a vector as generic input and a
      // vector as generic output, this does not mean that it is
      // the same vector!

      {
        String ws, mr;

        // Use parameter name only if it is used inside the function
        // to avoid warnings
        ws = " ws";
        //              if (!mdd.AgendaMethod() && !mdd.PassWorkspace() && !vo.nelem () && !vi.nelem () && !vgo.nelem () && !vgi.nelem ())
        //              {
        //                ws = "";
        //              }

        // Find out if the WSM gets an agenda as input. If so, pass
        // the current workspace to this method
        for (Index j = 0; !pass_workspace && j < mdd.In().nelem(); j++) {
          if (is_agenda_group_id(wsv_data[mdd.In()[j]].Group())) {
            pass_workspace = true;
          }
        }

        // Find out if the WSM gets an agenda as input. If so, pass
        // the current workspace to this method
        for (Index j = 0; !pass_workspace && j < mdd.GInType().nelem(); j++) {
          if (is_agenda_group_id(mdd.GInType()[j])) {
            pass_workspace = true;
          }
        }

        // Use parameter name only if it is used inside the function
        // to avoid warnings
        if (vo.nelem() || vi.nelem() || vgo.nelem() || vgi.nelem() ||
            mdd.AgendaMethod()) {
          mr = " mr";
        }

        if (mdd.Supergeneric()) {
          ofs << "void " << mdd.Name() << "_sg_" << mdd.ActualGroups()
              << "_g(Workspace&" << ws << ", const MRecord&" << mr << ")\n"
              << "{\n";
        } else {
          ofs << "void " << mdd.Name() << "_g(Workspace&" << ws
              << ", const MRecord&" << mr << ")\n"
              << "{\n";
        }
      }

      // Erase all Output only variables to uncover if they are
      // misused as Input variables
#ifndef NDEBUG
      // Determine indexes of variables in vo that are only use as output
      ArrayOfIndex voutonly;  // Output only
      for (Index k = 0; k < vo.nelem(); ++k) {
        bool output_only = true;
        for (ArrayOfIndex::const_iterator j = mdd.In().begin();
             j != mdd.In().end();
             ++j)
          if (vo[k] == *j) {
            output_only = false;
            break;
          }

        if (output_only) voutonly.push_back(k);
      }

      // Set builtin-type variables to a certain value,
      // initialize objects with their default constructor
      for (Index j = 0; j < voutonly.nelem(); j++) {
        ostringstream docstr;
        docstr << "  // " << wsv_data[vo[voutonly[j]]].Name() << "\n";

        String gname = wsv_groups[wsv_data[vo[voutonly[j]]].Group()];
        ostringstream initstr;
        if (gname == "Numeric")
          initstr << " = NAN;";
        else if (gname == "Index")
          initstr << " = -1;";
        else
          initstr << " = " << gname << "();";

        // Only zero-out the variable if it is not passed as an input
        // variable for any other parameter
        ofs << "  if (mr.In().end() == find(mr.In().begin(), mr.In().end(),"
            << " mr.Out()[" << voutonly[j] << "]))\n";
        ofs << "    (*(static_cast<" << wsv_groups[wsv_data[vo[voutonly[j]]].Group()]
            << "*>(ws[mr.Out()[" << voutonly[j] << "]].get())))" << initstr.str();
        ofs << docstr.str();
      }
#endif /* NDEBUG */

      ofs << "  " << mdd.Name() << "(";

      if (pass_workspace || mdd.PassWorkspace() || mdd.AgendaMethod()) {
        ofs << "ws";
        is_first_parameter = false;
      }

      // Write the Output workspace variables:
      for (Index j = 0; j < vo.nelem(); ++j) {
        // Check by assert whether the group identifier is too
        // large to correspond to a group. This can easily
        // happen if somebody puts a variable identifier instead
        // of a group identifier in the argument of GOUTPUT:
        ARTS_ASSERT(wsv_data[vo[j]].Group() < wsv_groups.nelem());

        // Add comma and line break, if not first element:
        align(ofs, is_first_parameter, indent);

        ofs << "*(static_cast<" << wsv_groups[wsv_data[vo[j]].Group()]
            << "*>(ws[mr.Out()[" << j << "]].get()))";
      }

      // Write the Generic output workspace variables:
      for (Index j = 0; j < vgo.nelem(); ++j) {
        // Check by assert whether the group identifier is too
        // large to correspond to a group. This can easily
        // happen if somebody puts a variable identifier instead
        // of a group identifier in the argument of GOUTPUT:
        ARTS_ASSERT(vgo[j] < wsv_groups.nelem());

        // Add comma and line break, if not first element:
        align(ofs, is_first_parameter, indent);

        ofs << "*(static_cast<" << wsv_groups[vgo[j]] << "*>(ws[mr.Out()["
            << j + vo.nelem() << "]].get()))";
      }

      // Write the Generic output workspace variable names:
      if (mdd.PassWsvNames()) {
        for (Index j = 0; j < vgo.nelem(); ++j) {
          // Add comma and line break, if not first element:
          align(ofs, is_first_parameter, indent);

          ofs << "Workspace::wsv_data[mr.Out()[" << j + vo.nelem()
              << "]].Name()";
        }
      }

      // Write the Input workspace variables:
      for (Index j = 0; j < vi.nelem(); ++j) {
        // Add comma and line break, if not first element:
        align(ofs, is_first_parameter, indent);

        if (is_agenda_group_id(wsv_data[vi[j]].Group())) {
          ofs << "*(static_cast<" << wsv_groups[wsv_data[vi[j]].Group()]
              << "*>(ws[mr.In()[" << j << "]].get()))";
        } else {
          ofs << "*(static_cast<" << wsv_groups[wsv_data[vi[j]].Group()]
              << "*>(ws[mr.In()[" << j << "]].get()))";
        }
      }

      // Write the control parameters:
      {
        if (mdd.SetMethod()) {
          // Add comma and line break, if not first element:
          align(ofs, is_first_parameter, indent);

          ofs << "mr.SetValue()";
        } else {
          // Write the Generic input workspace variables:
          for (Index j = 0; j < vgi.nelem(); ++j) {
            // Check by assert whether the group identifier is too
            // large to correspond to a group. This can easily
            // happen if somebody puts a variable identifier instead
            // of a group identifier in the argument of GINPUT:
            ARTS_ASSERT(vgi[j] < wsv_groups.nelem());

            // Add comma and line break, if not first element:
            align(ofs, is_first_parameter, indent);

            ofs << "*(static_cast<" << wsv_groups[vgi[j]] << "*>(ws[mr.In()["
                << j + vi.nelem() << "]].get()))";
          }

          // Write the Generic input workspace variable names:
          if (mdd.PassWsvNames()) {
            for (Index j = 0; j < vgi.nelem(); ++j) {
              // Add comma and line break, if not first element:
              align(ofs, is_first_parameter, indent);

              ofs << "Workspace::wsv_data[mr.In()[" << j + vi.nelem()
                  << "]].Name()";
            }
          }
        }
      }

      // Write the agenda, if there is one.
      if (mdd.AgendaMethod()) {
        align(ofs, is_first_parameter, indent);
        ofs << "mr.Tasks()";
      }

      // Flag that is set to false if the WSM has verbosity as an input or
      // output already. Otherwise it's passed as the last parameter.
      bool pass_verbosity = true;

      // Find out if the WSM has the verbosity as input.
      for (Index j = 0; pass_verbosity && j < mdd.In().nelem(); j++) {
        if (wsv_data[mdd.In()[j]].Name() == "verbosity") {
          pass_verbosity = false;
        }
      }

      // Find out if the WSM has the verbosity as output.
      for (Index j = 0; pass_verbosity && j < mdd.Out().nelem(); j++) {
        if (wsv_data[mdd.Out()[j]].Name() == "verbosity") {
          pass_verbosity = false;
        }
      }

      if (pass_verbosity) {
        static Index verbosity_wsv_id = get_wsv_id("verbosity");
        static Index verbosity_group_id = get_wsv_group_id("Verbosity");
        align(ofs, is_first_parameter, indent);
        ofs << "*(static_cast<" << wsv_groups[verbosity_group_id] << "*>(ws["
            << verbosity_wsv_id << "].get()))";
      }

      ofs << ");\n";
      ofs << "}\n\n";
    }

    // Add getaways, the array that hold pointers to the getaway functions:
    {
      String indent = "     ";
      bool is_first_parameter = true;

      ofs << "// The array holding the pointers to the getaway functions.\n"
          << "void (*getaways[])(Workspace&, const MRecord&)\n"
          << "  = {";
      for (Index i = 0; i < n_md; ++i) {
        const MdRecord& mdd = md_data[i];

        // Add comma and line break, if not first element:
        align(ofs, is_first_parameter, indent);

        if (mdd.Supergeneric()) {
          ofs << mdd.Name() << "_sg_" << mdd.ActualGroups() << "_g";
        } else {
          ofs << mdd.Name() << "_g";
        }
      }
      // Add empty slot which can be used by the API for external callbacks.
      ofs << ", nullptr};\n\n";
    }

    // Agenda execute helper function

    ofs << "void auto_md_agenda_execute_helper(bool& agenda_failed, String& agenda_error_msg, Workspace& ws, const Agenda& input_agenda)\n";
    ofs << "{\n";
    ofs << "    const ArrayOfIndex& outputs_to_push = input_agenda.get_output2push();\n";
    ofs << "    const ArrayOfIndex& outputs_to_dup = input_agenda.get_output2dup();\n";
    ofs << "\n";
    ofs << "    for (auto&& i : outputs_to_push)\n";
    ofs << "    {\n";
    ofs << "        // Even if a variable is only used as WSM output inside this agenda,\n";
    ofs << "        // It is possible that it is used as input further down by another agenda,\n";
    ofs << "        // which we can't see here. Therefore initialized variables have to be\n";
    ofs << "        // duplicated.\n";
    ofs << "        if (ws.is_initialized(i))\n";
    ofs << "            ws.duplicate(i);\n";
    ofs << "        else\n";
    ofs << "            ws.emplace(i);\n";
    ofs << "    }\n";
    ofs << "\n";
    ofs << "    for (auto&& i : outputs_to_dup)\n";
    ofs << "        ws.duplicate(i);\n";
    ofs << "\n";
    ofs << "    agenda_failed = false;\n";
    ofs << "    try\n";
    ofs << "    {\n";
    ofs << "        input_agenda.execute(ws);\n";
    ofs << "    }\n";
    ofs << "    catch (const std::exception &e)\n";
    ofs << "    {\n";
    ofs << "        ostringstream os;\n";
    ofs << "        os << \"Run-time error in agenda: \"\n";
    ofs << "           << input_agenda.name() << \'\\n\' << e.what();\n";
    ofs << "        agenda_failed = true;\n";
    ofs << "        agenda_error_msg = os.str();\n";
    ofs << "    }\n";
    ofs << "\n";
    ofs << "    for (auto&& i : outputs_to_push)\n";
    ofs << "        ws.pop(i);\n";
    ofs << "\n";
    ofs << "    for (auto&& i : outputs_to_dup)\n";
    ofs << "        ws.pop(i);\n";
    ofs << "}\n\n";

    // Create implementation of the agenda wrappers

    // Initialize agenda data.
    Workspace::define_wsv_map();
    define_agenda_data();

    using global_data::agenda_data;
    for (Index i = 0; i < agenda_data.nelem(); i++) {
      const AgRecord& agr = agenda_data[i];
      const ArrayOfIndex& ago = agr.Out();
      const ArrayOfIndex& agi = agr.In();
      ostringstream ain_push_os, ain_pop_os;
      ostringstream aout_push_os, aout_pop_os;

      bool is_agenda_array = wsv_data[get_wsv_id(agr.Name())].Group() ==
                             get_wsv_group_id("ArrayOfAgenda");
      write_agenda_wrapper_header(ofs, agr, is_agenda_array);

      ofs << "\n";
      ofs << "{\n";

      if (ago.nelem() || agi.nelem()) {
        if (is_agenda_array) {
          ofs << "  if (agenda_array_index < 0 || agenda_array_index >= input_agenda_array.nelem())\n"
              << "  {\n"
              << "      std::ostringstream os;\n"
              << "      os << \"Agenda index \" << agenda_array_index\n"
              << "         << \" out of bounds. 0 <= index < \" << input_agenda_array.nelem();\n"
              << "      throw std::runtime_error(os.str());\n"
              << "  }\n\n"
              << "  const Agenda& input_agenda = input_agenda_array[agenda_array_index];\n\n";
        }
        ofs << "  using global_data::AgendaMap;\n"
            << "  using global_data::agenda_data;\n"
            << "\n"
            << "  if (!input_agenda.checked())\n"
            << "    throw std::runtime_error(\"" << agr.Name()
            << " is uninitialized. Use *AgendaSet* to add methods to it.\\n"
            << "Agenda variables defined in the controlfile cannot be executed.\\nThe agenda must "
            << "be copied to a workspace variable for execution.\");\n"
            << "\n"
            << "  const AgRecord& agr =\n"
            << "    agenda_data[AgendaMap.find (input_agenda.name ())->second];\n"
            << "\n";
      }
      if (ago.nelem()) {
        for (Index j = 0; j < ago.nelem(); j++) {
          // Mark agenda output-only variables as uninitialized
          auto it = agi.begin();
          while (it != agi.end() && *it != ago[j]) it++;
          if (it == agi.end()) {
            aout_push_os << "  auto borrowed_out_"<<j<<"{ws.borrow_uninitialized (aout[" << j << "], "
                         << wsv_data[ago[j]].Name() << ")};\n";
          } else {
            aout_push_os << "  auto borrowed_out_"<<j<<"{ws.borrow (aout[" << j << "], "
                         << wsv_data[ago[j]].Name() << ")};\n";
          }
        }
      }
      if (agi.nelem()) {
        for (Index j = 0; j < agi.nelem(); j++) {
          // Ignore Input parameters that are also output
          auto it = ago.begin();
          while (it != ago.end() && *it != agi[j]) it++;
          if (it == ago.end()) {
            ain_push_os << "  auto borrowed_in_"<<j<<"{ws.borrow (ain[" << j << "], "
                        << wsv_data[agi[j]].Name() << ")};\n";
          }
        }
      }

      if (aout_push_os.str().length()) {
        ofs << "  const ArrayOfIndex& aout = agr.Out();\n";
        ofs << aout_push_os.str() << "\n";
      }
      if (ain_push_os.str().length()) {
        ofs << "  const ArrayOfIndex& ain = agr.In();\n";
        ofs << ain_push_os.str() << "\n";
      }

      ofs << "  bool agenda_failed;\n";
      ofs << "  String agenda_error_msg;\n";
      ofs << "  auto_md_agenda_execute_helper(agenda_failed, agenda_error_msg, ws, input_agenda);\n\n";

      if (aout_pop_os.str().length()) {
        ofs << aout_pop_os.str() << "\n";
      }
      if (ain_pop_os.str().length()) {
        ofs << ain_pop_os.str() << "\n";
      }

      ofs << "  if (agenda_failed) throw runtime_error (agenda_error_msg);\n";
      ofs << "}\n\n";
    }

    // Create implementation of the GroupCreate WSMs
    //
    for (auto&& it : wsv_groups) {
      if (it != "Any") {
        ofs << "/* Workspace method: Doxygen documentation will be auto-generated */\n"
            << "void " << it << "Create(" << it << "& var, const Verbosity&)\n"
            << "{ ";

        // Treat atomic types separately.
        // For objects the default constructor is used.
        if (it == "Index")
          ofs << "var = 0;";
        else if (it == "Numeric")
          ofs << "var = 0.;";
        else
          ofs << "var = " << it << "();";

        ofs << " }\n\n";

        if (it not_eq "Agenda" and it not_eq "ArrayOfAgenda") {
          ofs << "/* Workspace method: Doxygen documentation will be auto-generated */\n"
              << "void " << it << "Set(" << it << "& out, const " << it
              << "& value, const Verbosity&) { out = value; }\n\n";
        }
      }
    }
  } catch (const std::runtime_error& x) {
    cout << "Something went wrong. Message text:\n";
    cout << x.what() << '\n';
    return 1;
  }

  return 0;
}
