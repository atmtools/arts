/* Copyright (C) 2000-2008 Stefan Buehler <sbuehler@ltu.se>

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

#include "arts.h"
#include "array.h"
#include "file.h"
#include "methods.h"
#include "workspace_ng.h"
#include "agenda_record.h"

/* Adds commas and indentation to parameter lists. */
void align(ofstream& ofs, bool& is_first_parameter, const String& indent)
{
  // Add comma and line break, if not first element:
  if (is_first_parameter)
    is_first_parameter = false;
  else
    {
      ofs << ",\n";
      // Make proper indentation:
      ofs << indent;
    }
}

int main()
{
  try
    {
      // Make the global data visible:
      extern Array<MdRecord> md_data;
      extern const ArrayOfString wsv_group_names;
      const Array<WsvRecord>& wsv_data = Workspace::wsv_data;

      // Initialize the wsv group name array:
      define_wsv_group_names();

      // Initialize wsv data.
      Workspace::define_wsv_data();
  
      // Initialize WsvMap.
      Workspace::define_wsv_map();

      // Initialize method data.
      define_md_data_raw();

      // Expand supergeneric methods:
      expand_md_data_raw_to_md_data();


      const Index n_md  = md_data.nelem();

      // Write auto_md.cc:
      // -----------
      ofstream ofs;
      open_output_file(ofs,"auto_md.cc");
  
      ofs << "// This file was generated automatically by make_auto_md_cc.cc.\n";
      ofs << "// DO NOT EDIT !\n";
      ofs << "// Generated: "
          << __DATE__ << ", "
          << __TIME__ << "\n\n";

      ofs << "#include \"arts.h\"\n"
          << "#include \"make_array.h\"\n"
          << "#include \"auto_md.h\"\n"
          << "#include \"wsv_aux.h\"\n"
          << "#include \"mc_interp.h\"\n"
          << "#include \"m_append.h\"\n"
          << "#include \"m_delete.h\"\n"
          << "#include \"m_copy.h\"\n"
          << "#include \"m_general.h\"\n"
          << "#include \"m_ignore.h\"\n"
          << "#include \"m_xml.h\"\n"
          << "#include \"m_basic_types.h\"\n"
          << "#include \"agenda_record.h\"\n"
          << "#include \"workspace_ng.h\"\n"
          << "\n";

      //ofs << "static Index agendacallcount = 0;\n";

      // Write all get-away functions:
      // -----------------------------
      for (Index i=0; i<n_md; ++i)
        {
          const MdRecord& mdd = md_data[i];

          // This is needed to flag the first function parameter, which 
          // needs no line break before being written:
          bool is_first_parameter = true;
          // The String indent is needed to achieve the correct
          // indentation of the functin parameters:
          String indent = String(mdd.Name().nelem()+3,' ');;
          
          // There are four lists of parameters that we have to
          // write. 
          ArrayOfIndex  vo=mdd.Output();   // Output 
          ArrayOfIndex  vi;                // Input
          ArrayOfIndex  vgo=mdd.GOutput(); // Generic Output 
          ArrayOfIndex  vgi=mdd.GInput();  // Generic Input
          ArrayOfIndex  kgi=mdd.Types();   // Keyword Input
          // vo and vi contain handles of workspace variables, 
          // vgo and vgi handles of workspace variable groups.

          mdd.input_only(vi);

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
              if (!mdd.PassWorkspace() && !vo.nelem () && !vi.nelem () && !vgo.nelem () && !vgi.nelem () && !kgi.nelem())
              {
                ws = "";
              }

              // Use parameter name only if it is used inside the function
              // to avoid warnings
              if ( vo.nelem () || vi.nelem () || vgo.nelem () || vgi.nelem ()
                   || mdd.Keywords().nelem() || mdd.AgendaMethod())
                {
                  mr = " mr";
                }

              if ( mdd.Supergeneric() )
                {
                  ofs << "void " << mdd.Name()
                    << "_sg_" << wsv_group_names[mdd.ActualGroup()]
                    << "_g(Workspace&" << ws
                    << ", const MRecord&" << mr << ")\n"
                    << "{\n";
                }
              else
                {
                  ofs << "void " << mdd.Name()
                    << "_g(Workspace&" << ws
                    << ", const MRecord&" << mr << ")\n"
                    << "{\n";
                }
            }

          // Create copy of input agendas with private workspace
          for (Index j=0; j<vi.nelem(); ++j)
            {
              if (wsv_data[vi[j]].Group() == get_wsv_group_id("Agenda"))
                {
#ifdef _OPENMP
                  ofs << "  Agenda AI" << j
#else
                  ofs << "  Agenda& AI" << j
#endif
                    << " = *(("
                    << wsv_group_names[wsv_data[vi[j]].Group()]
                    << " *)ws[mr.Input()[" << j
                    << "]]);\n";
                  ofs << "  AI" << j << ".set_workspace(&ws);\n";
                }
            }

          ofs << "  " << mdd.Name() << "(";

          if (mdd.PassWorkspace())
            {
              ofs << "ws";
              is_first_parameter = false;
            }

          // Write the Output workspace variables:
          for (Index j=0; j<vo.nelem(); ++j)
            {
              // Check by assert whether the group identifier is too
              // large to correspond to a group. This can easily
              // happen if somebody puts a variable identifier instead
              // of a group identifier in the argument of GOUTPUT:
              assert( wsv_data[vo[j]].Group() < wsv_group_names.nelem() );

              // Add comma and line break, if not first element:
              align(ofs,is_first_parameter,indent);

              ofs << "*(("
                  << wsv_group_names[wsv_data[vo[j]].Group()]
                  << " *)ws[mr.Output()[" << j
                  << "]])";
            }

          // Write the Generic output workspace variables:
          for (Index j=0; j<vgo.nelem(); ++j)
            {
              // Check by assert whether the group identifier is too
              // large to correspond to a group. This can easily
              // happen if somebody puts a variable identifier instead
              // of a group identifier in the argument of GOUTPUT:
              assert( vgo[j] < wsv_group_names.nelem() );

              // Add comma and line break, if not first element:
              align(ofs,is_first_parameter,indent);

              ofs << "*((" << wsv_group_names[vgo[j]]
                  << " *)ws[mr.Output()[" << j+vo.nelem()
                  << "]])";
            }

          // Write the Generic output workspace variable names:
          if (mdd.PassWsvNames())
            {
              for (Index j=0; j<vgo.nelem(); ++j)
                {
                  // Add comma and line break, if not first element:
                  align(ofs,is_first_parameter,indent);

                  ofs << "Workspace::wsv_data[mr.Output()["
                    << j+vo.nelem()
                    << "]].Name()";
                }
            }

          // Write the Input workspace variables:
          for (Index j=0; j<vi.nelem(); ++j)
            {
              // Add comma and line break, if not first element:
              align(ofs,is_first_parameter,indent);

              if (wsv_data[vi[j]].Group() == get_wsv_group_id("Agenda"))
                {
                  ofs << "AI" << j;
                }
              else
                {
                  ofs << "*(("
                    << wsv_group_names[wsv_data[vi[j]].Group()]
                    << " *)ws[mr.Input()[" << j
                    << "]])";
                }
            }

          // Write the Generic input workspace variables:
          for (Index j=0; j<vgi.nelem(); ++j)
            {
              // Check by assert whether the group identifier is too
              // large to correspond to a group. This can easily
              // happen if somebody puts a variable identifier instead
              // of a group identifier in the argument of GINPUT:
              assert( vgi[j] < wsv_group_names.nelem() );

              // Add comma and line break, if not first element:
              align(ofs,is_first_parameter,indent);

              ofs << "*((" << wsv_group_names[vgi[j]]
                  << " *)ws[mr.Input()[" << j+vi.nelem()
                  << "]])";
            }

          // Write the Generic input workspace variable names:
          if (mdd.PassWsvNames())
            {
              for (Index j=0; j<vgi.nelem(); ++j)
                {
                  // Add comma and line break, if not first element:
                  align(ofs,is_first_parameter,indent);

                  ofs << "Workspace::wsv_data[mr.Input()["
                    << j+vi.nelem()
                    << "]].Name()";
                }
            }

          // Write the control parameters:
          {
            if (mdd.SetMethod())
              {
                // Add comma and line break, if not first element:
                align(ofs,is_first_parameter,indent);

                ofs << "mr.SetValue()";
              }
            else
              {
                // Write the Keyword input workspace variables:
                for (Index j=0; j<kgi.nelem(); ++j)
                  {
                    // Check by assert whether the group identifier is too
                    // large to correspond to a group. This can easily
                    // happen if somebody puts a variable identifier instead
                    // of a group identifier in the argument of GINPUT:
                    assert( kgi[j] < wsv_group_names.nelem() );

                    // Add comma and line break, if not first element:
                    align(ofs,is_first_parameter,indent);

                    ofs << "*((" << wsv_group_names[kgi[j]]
                      << " *)ws[mr.Values()[" << j
                      << "]])";
                  }
              }
          }

          // Write the agenda, if there is one.
          if ( mdd.AgendaMethod() )
            {
              align(ofs,is_first_parameter,indent);
              ofs << "mr.Tasks()";
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
        for (Index i=0; i<n_md; ++i)
          {
            const MdRecord& mdd = md_data[i];
          
            // Add comma and line break, if not first element:
            align(ofs,is_first_parameter,indent);

            if ( mdd.Supergeneric() )
              {
                ofs << mdd.Name()
                    << "_sg_" << wsv_group_names[mdd.ActualGroup()]
                    << "_g";
              }
            else
              {
                ofs << mdd.Name() << "_g";
              }
          }
        ofs << "};\n\n";
      }

      // Create implementation of the agenda wrappers

      // Initialize agenda data.
      Workspace::define_wsv_map ();
      define_agenda_data ();

      extern const Array<AgRecord> agenda_data;
      for (Index i = 0; i < agenda_data.nelem (); i++)
        {
          const AgRecord& agr = agenda_data[i];
          const ArrayOfIndex& ago = agr.Output ();
          const ArrayOfIndex& agi = agr.Input ();
          ostringstream ain_push_os, ain_pop_os;
          ostringstream aout_push_os, aout_pop_os;

          write_agenda_wrapper_header (ofs, agr);

          ofs << "\n";
          ofs << "{\n";
          ofs << "  Agenda new_input_agenda (input_agenda);\n";
          ofs << "  if (safe_workspace)\n"
              << "    new_input_agenda.set_workspace (new Workspace (*new_input_agenda.workspace()));\n";
          ofs << "  Workspace& ws = *(new_input_agenda.workspace());\n";

          if (ago.nelem () || agi.nelem ())
            {
              ofs << "  extern map<String, Index> AgendaMap;\n"
                << "  extern const Array<AgRecord> agenda_data;\n"
                << "\n"
                << "  const AgRecord& agr =\n"
                << "    agenda_data[AgendaMap.find (new_input_agenda.name ())->second];\n"
                << "\n";
            }
          if (ago.nelem ())
            {
              for (Index j = 0; j < ago.nelem (); j++)
                {
                  // Mark agenda output-only variables as uninitialized
                  ArrayOfIndex::const_iterator it = agi.begin ();
                  while (it != agi.end () && *it != ago[j]) it++;
                  if (it == agi.end ())
                    {
                      aout_push_os << "  ws.push_uninitialized (aout[" << j << "], "
                        << "(void *)&" << wsv_data[ago[j]].Name () << ");\n";
                    }
                  else
                    {
                      aout_push_os << "  ws.push (aout[" << j << "], "
                        << "(void *)&" << wsv_data[ago[j]].Name () << ");\n";
                    }
                  aout_pop_os << "  ws.pop (aout[" << j << "]);\n";
                }
            }
          if (agi.nelem ())
            {
              for (Index j = 0; j < agi.nelem (); j++)
                {
                  // Ignore Input parameters that are also output
                  ArrayOfIndex::const_iterator it = ago.begin ();
                  while (it != ago.end () && *it != agi[j]) it++;
                  if (it == ago.end ())
                    {
                      ain_push_os << "  ws.push (ain[" << j << "], "
                        << "(void *)&" << wsv_data[agi[j]].Name () << ");\n";
                      ain_pop_os << "  ws.pop (ain[" << j << "]);\n";
                    }
                }
            }

          if (aout_push_os.str().length())
            {
              ofs << "  const ArrayOfIndex& aout = agr.Output ();\n";
              ofs << aout_push_os.str () << "\n";
            }
          if (ain_push_os.str().length())
            {
              ofs << "  const ArrayOfIndex& ain = agr.Input ();\n";
              ofs << ain_push_os.str () << "\n";
            }

          ofs << "  const ArrayOfIndex& outputs_to_push = new_input_agenda.get_output2push();\n"
              << "  const ArrayOfIndex& outputs_to_dup = new_input_agenda.get_output2dup();\n"
              << "\n"
              << "  for (ArrayOfIndex::const_iterator it = outputs_to_push.begin ();\n"
              << "       it != outputs_to_push.end (); it++)\n"
              << "  { ws.push (*it, NULL); }\n"
              << "\n"
              << "  for (ArrayOfIndex::const_iterator it = outputs_to_dup.begin ();\n"
              << "       it != outputs_to_dup.end (); it++)\n"
              << "  { ws.duplicate (*it); }\n"
              << "\n";

          ofs << "  String agenda_error_msg;\n"
              << "  bool agenda_failed = false;\n\n"
              << "  try {\n"
              << "    new_input_agenda.execute (silent);\n"
              << "  } catch (runtime_error e) {\n"
              << "    ostringstream os;\n"
              << "    os << \"Run-time error in agenda: \"\n"
              << "       << new_input_agenda.name() << \'\\n\' << e.what();\n"
              << "    agenda_failed = true;\n"
              << "    agenda_error_msg = os.str();\n"
              << "  }\n";

          ofs << "  for (ArrayOfIndex::const_iterator it = outputs_to_push.begin ();\n"
              << "       it != outputs_to_push.end (); it++)\n"
              << "    { ws.pop_free (*it); }\n"
              << "\n"
              << "  for (ArrayOfIndex::const_iterator it = outputs_to_dup.begin ();\n"
              << "       it != outputs_to_dup.end (); it++)\n"
              << "    { ws.pop_free (*it); }\n\n";

          if (aout_pop_os.str().length())
            {
              ofs << aout_pop_os.str () << "\n";
            }
          if (ain_pop_os.str().length())
            {
              ofs << ain_pop_os.str () << "\n";
            }

          ofs << "  if (safe_workspace)\n"
              << "    delete new_input_agenda.workspace();\n";

          ofs << "  if (agenda_failed) throw runtime_error (agenda_error_msg);\n\n";

          ofs << "}\n\n";
        }
    }
  catch (runtime_error x)
    {
      cout << "Something went wrong. Message text:\n";
      cout << x.what() << '\n';
      return 1;
    }

  return 0;
}
