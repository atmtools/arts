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

/*!
  \file   make_auto_md_h.cc
  \brief  This is a little C++ program that generates the file auto_md.h from the
          workspace methods data md_data. 

  The file auto_md.h declares the enum
  type MdHandle that is used to access the method data, so it has
  to be made sure that the two are allways consistent.

  A second file is produced: auto_md.cc.
  This contains the `get-away' functions that provided the interface
  between the engine and the workspace methods. The get-functions all 
  have the same arguments:

  \code
     void get_away_example_g(Workspace& ws,
     const MRecord& mr);
  \endcode

  Their names all have the extension _g

  Pointers to the get-away functions are stored in the array
  `getaway'. 

  Each get-away function simply contains a function call to the
  matching workspace method. The parameters are arranged similar to
  the follwing example:
   
  \code
  void SomeMethod(owsv1,iwsv1,iwsv2,iwsv3,c1,c2,c3,...)
  \endcode

  First come the output workspace variables, then the input workspace 
  variables, and then the control parameters. There can be an
  arbitrary number of parameters of each type, but the most usual
  case is to have only one output workspace variable. 

  The same variable may be both in the list of input and in the list
  of output workspace variables. This case makes good sense,
  actually, if you think for example of a method that adds an offset
  to the absorption coefficients. IN THAT CASE THE VARIABLE IS ADDED
  TO THE LIST ONLY ONCE, namely among the OUTPUT variables.

  For generic methods the names of the actual workspace variables are
  also passed on to the method function.

  \author Stefan Buehler
  \date   2000-07-29
*/

#include "arts.h"
#include "array.h"
#include "file.h"
#include "methods.h"
#include "workspace_ng.h"
#include "agenda_record.h"
#include "global_data.h"

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

//! Write method header documentation.
/*!
  \param[in,out] ofs The stream to write to.
  \param mdd Method lookup data.
*/
void write_method_header_documentation (ofstream& ofs, const MdRecord& mdd)
{
  const Array<WsvRecord>& wsv_data = Workspace::wsv_data;

  String fullname = mdd.Name();

  // This is needed to flag the first function parameter, which
  // needs no line break before being written:
  bool is_first_parameter = true;

  // The String indent is needed to achieve the correct
  // indentation of the functin parameters:
  String indent("  ");

  // Flag to pass the workspace to the WSM. Only true if the WSM has
  // an Agenda as input.
  bool pass_workspace = false;

  // There are four lists of parameters that we have to
  // write.
  ArrayOfIndex  vo=mdd.Out();   // Output
  const ArrayOfIndex &vi = mdd.InOnly(); // Input
  ArrayOfIndex  vgo=mdd.GOutType(); // Generic Output
  ArrayOfIndex  vgi=mdd.GInType();  // Generic Input
  // vo and vi contain handles of workspace variables,
  // vgo and vgi handles of workspace variable groups.

  // Find out if the WSM gets an agenda as input. If so, pass
  // the current workspace to this method
  for (Index j = 0; !pass_workspace && j < mdd.In().nelem(); j++)
    {
      if (is_agenda_group_id(wsv_data[mdd.In()[j]].Group()))
        {
          pass_workspace = true;
        }
    }

  // Find out if the WSM gets an agenda as input. If so, pass
  // the current workspace to this method
  for (Index j = 0; !pass_workspace && j < mdd.GInType().nelem(); j++)
    {
      if (is_agenda_group_id(mdd.GInType()[j]))
        {
          pass_workspace = true;
        }
    }

  // Start with the name of the one line description
  ofs << "//! WORKSPACE METHOD: " << fullname << ".\n";

  ofs << "/*!\n";

  String DoxyDescription = mdd.Description();
  size_t start_pos = 0;

  while (start_pos != string::npos)
    {
      start_pos = DoxyDescription.find ("\n ", start_pos);
      if (start_pos && start_pos != string::npos
          && DoxyDescription[start_pos]-1 != '\n')
        {
          DoxyDescription.insert (start_pos + 1, "<br>");
          start_pos += 5;
        }
    }

  // Some characters have to be masqueraded to appear properly in doxygen
  // documentation
  DoxyDescription.insert_substr ("<", "\\");
  DoxyDescription.insert_substr (">", "\\");

  ofs << DoxyDescription << "\n";

  // Write the authors:
  for (Index j=0; j<mdd.Authors().nelem(); ++j)
    {
      ofs << indent << "\\author " << mdd.Authors ()[j] << "\n";
    }

  ofs << "\n";

  if (pass_workspace || mdd.PassWorkspace() || mdd.AgendaMethod())
    {
      ofs << indent << "\\param[in,out] " << "ws Workspace\n";
    }

  // Write the Output workspace variables:
  for (Index j=0; j<vo.nelem(); ++j)
    {
      bool inout = (std::find(mdd.In().begin(), mdd.In().end(), vo[j])
                    != mdd.In().end());

      if (inout)
        {
          ofs << indent << "\\param[in,out] "
              << wsv_data[vo[j]].Name() << " WS Input/Output\n";
        }
      else
        {
          ofs << indent << "\\param[out] "
              << wsv_data[vo[j]].Name() << " WS Output\n";
        }
    }

  // Write the Generic output workspace variables:
  for (Index j=0; j<vgo.nelem(); ++j)
    {
      ofs << indent << "\\param[out]    ";

      if (mdd.GOut()[j].length())
        ofs << mdd.GOut()[j];
      else
        ofs << "genericoutput" << j+1;

      if (mdd.Supergeneric ()) ofs << " Supergeneric output\n";
      else ofs << " Generic output\n";
    }

  // Write the Generic output workspace variable names:
  if (mdd.PassWsvNames())
    {
      for (Index j=0; j<vgo.nelem(); ++j)
        {
          ofs << indent << "\\param[in]     ";
          if (mdd.GOut()[j].length())
            ofs << mdd.GOut()[j] << "_wsvname";
          else
            ofs << "genericoutput" << j+1 << "_wsvname";

          ofs << " Generic Output Name" << endl;
        }
    }

  // Write the Input workspace variables:
  for (Index j=0; j<vi.nelem(); ++j)
    {
      ofs << indent << "\\param[in]     "
        << wsv_data[vi[j]].Name() << " WS Input\n";
    }

  // Write the Generic input workspace variables:
  for (Index j=0; j<vgi.nelem(); ++j)
    {
      ofs << indent << "\\param[in]     ";
      if (mdd.GIn()[j] != "")
        ofs << mdd.GIn()[j];
      else
        ofs << "genericinput" << j+1;

      ofs << " Generic Input";

      if (mdd.GInDefault()[j] != NODEF)
        {
          ofs << " (Default: \"" << mdd.GInDefault()[j] << "\")";
        }
      ofs << endl;
    }

  // Write the Generic input workspace variable names:
  if (mdd.PassWsvNames())
    {
      for (Index j=0; j<vgi.nelem(); ++j)
        {
          ofs << indent << "\\param[in]     ";
          if (mdd.GIn()[j].length())
            ofs << mdd.GIn()[j] << "_wsvname";
          else
            ofs << "genericinput" << j+1 << "_wsvname";

          ofs << " Generic Input Name" << endl;
        }
    }


  // Write agenda, if there is one:
  if ( mdd.AgendaMethod() )
    {
      align(ofs,is_first_parameter,indent);
      ofs << indent << "\\param[in]     " << "input_agenda Agenda from controlfile\n";
    }

  ofs << "*/\n";
}

//! Write a method header.
/*!
  \param ofs The stream to write to.
  \param mdd Method lookup data.
*/
void write_method_header( ofstream& ofs,
                          const MdRecord& mdd )
{
  using global_data::wsv_group_names;
  const Array<WsvRecord>& wsv_data = Workspace::wsv_data;

//   // Work out the full name to use:
//   String fullname;
//   {
//     ostringstream os;
//     os << mdd.Name() << add_to_name;
//     fullname = os.str();
//   }

  String fullname = mdd.Name();

  // This is needed to flag the first function parameter, which
  // needs no line break before being written:
  bool is_first_parameter = true;

  // The String indent is needed to achieve the correct
  // indentation of the functin parameters:
  String indent(fullname.nelem()+6,' ');

  // Flag to pass the workspace to the WSM. Only true if the WSM has
  // an Agenda as input.
  bool pass_workspace = false;

  // There are four lists of parameters that we have to
  // write.
  ArrayOfIndex  vo=mdd.Out();   // Output
  const ArrayOfIndex &vi = mdd.InOnly(); // Input
  ArrayOfIndex  vgo=mdd.GOutType(); // Generic Output
  ArrayOfIndex  vgi=mdd.GInType();  // Generic Input
  // vo and vi contain handles of workspace variables,
  // vgo and vgi handles of workspace variable groups.

  // Find out if the WSM gets an agenda as input. If so, pass
  // the current workspace to this method
  for (Index j = 0; !pass_workspace && j < mdd.In().nelem(); j++)
    {
      if (is_agenda_group_id(wsv_data[mdd.In()[j]].Group()))
        {
          pass_workspace = true;
        }
    }

  // Find out if the WSM gets an agenda as input. If so, pass
  // the current workspace to this method
  for (Index j = 0; !pass_workspace && j < mdd.GInType().nelem(); j++)
    {
      if (is_agenda_group_id(mdd.GInType()[j]))
        {
          pass_workspace = true;
        }
    }

  // There used to be a similar block here for the generic
  // input/output variables. However, this was a mistake. For
  // example, if a method has a vector as generic input and a
  // vector as generic output, this does not mean that it is
  // the same vector!

  if (mdd.Supergeneric() && mdd.UsesTemplates())
    {
      ofs << "template <typename T>" << endl;
    }

  // Start with the name of the method:
  ofs << "void " << fullname << "(";

  if (pass_workspace || mdd.PassWorkspace() || mdd.AgendaMethod())
    {
      ofs << "// Workspace reference:\n";
      ofs << indent << "Workspace& ws";
      is_first_parameter = false;
    }

  // Write the Output workspace variables:
  {
    // Flag first parameter of this sort:
    bool is_first_of_these = true;

    for (Index j=0; j<vo.nelem(); ++j)
      {
        // Add comma and line break, if not first element:
        align(ofs,is_first_parameter,indent);

        // Add comment if this is the first of this sort
        if (is_first_of_these)
          {
            ofs << "// WS Output:\n";
            ofs << indent;
            is_first_of_these = false;
          }

        ofs << wsv_group_names[Workspace::wsv_data[vo[j]].Group()] << "& "
          << Workspace::wsv_data[vo[j]].Name();
      }
  }

  // Write the Generic output workspace variables:
  {
    // Flag first parameter of this sort:
    bool is_first_of_these = true;

    for (Index j=0; j<vgo.nelem(); ++j)
      {
        // Add comma and line break, if not first element:
        align(ofs,is_first_parameter,indent);

        // Add comment if this is the first of this sort
        if (is_first_of_these)
          {
            ofs << "// WS Generic Output:\n";
            ofs << indent;
            is_first_of_these = false;
          }

        if (wsv_group_names[mdd.GOutType()[j]] == "Any") ofs << "T& ";
        else ofs << wsv_group_names[mdd.GOutType()[j]] << "& ";

        if (mdd.GOut()[j].length()) ofs << mdd.GOut()[j];
        else ofs << "genericoutput" << j+1;
      }
  }

  // Write the Generic output workspace variable names:
  if (mdd.PassWsvNames())
    {
      // Flag first parameter of this sort:
      bool is_first_of_these = true;

      for (Index j=0; j<vgo.nelem(); ++j)
        {
          // Add comma and line break, if not first element:
          align(ofs,is_first_parameter,indent);

          // Add comment if this is the first of this sort
          if (is_first_of_these)
            {
              ofs << "// WS Generic Output Names:\n";
              ofs << indent;
              is_first_of_these = false;
            }

          ofs << "const String& ";
          if (mdd.GOut()[j].length())
            ofs << mdd.GOut()[j] << "_wsvname";
          else
            ofs << "genericoutput" << j+1 << "_wsvname";
        }
    }

  // Write the Input workspace variables:
  {
    // Flag first parameter of this sort.
    bool is_first_of_these = true;

    for (Index j=0; j<vi.nelem(); ++j)
      {
        // Add comma and line break, if not first element:
        align(ofs,is_first_parameter,indent);

        // Add type if this is the first of this sort.
        if (is_first_of_these)
          {
            ofs << "// WS Input:\n";
            ofs << indent;
            is_first_of_these = false;
          }

        ofs << "const "
          << wsv_group_names[Workspace::wsv_data[vi[j]].Group()] << "& "
          << Workspace::wsv_data[vi[j]].Name();
      }
  }

  // Write the Generic input workspace variables:
  {
    // Flag first parameter of this sort.
    bool is_first_of_these = true;

    for (Index j=0; j<vgi.nelem(); ++j)
      {
        // Add comma and line break, if not first element:
        align(ofs,is_first_parameter,indent);

        // Add type if this is the first of this sort.
        if (is_first_of_these)
          {
            ofs << "// WS Generic Input:\n";
            ofs << indent;
            is_first_of_these = false;
          }

        if (wsv_group_names[mdd.GInType()[j]] == "Any")
          {
            ofs << "const T& ";
            if (mdd.GIn()[j].length()) ofs << mdd.GIn()[j];
            else ofs << "genericinput" << j+1;
          }
        else
          {
            ofs << "const " << wsv_group_names[mdd.GInType()[j]] << "& ";
            if (mdd.GIn()[j].length()) ofs << mdd.GIn()[j];
            else ofs << "genericinput" << j+1;
          }
      }
  }

  // Write the Generic input workspace variable names:
  if (mdd.PassWsvNames())
    {
      // Flag first parameter of this sort:
      bool is_first_of_these = true;

      for (Index j=0; j<vgi.nelem(); ++j)
        {
          // Add comma and line break, if not first element:
          align(ofs,is_first_parameter,indent);

          // Add comment if this is the first of this sort
          if (is_first_of_these)
            {
              ofs << "// WS Generic Input Names:\n";
              ofs << indent;
              is_first_of_these = false;
            }

          ofs << "const String& ";
          if (mdd.GIn()[j].length())
            ofs << mdd.GIn()[j] << "_wsvname";
          else
            ofs << "genericinput" << j+1 << "_wsvname";
        }
    }

  // Write agenda, if there is one:
  if ( mdd.AgendaMethod() )
    {
      align(ofs,is_first_parameter,indent);
      ofs << "// Agenda from controlfile:\n";
      ofs << indent;
      ofs << "const Agenda& input_agenda";
    }

  // Flag that is set to false if the WSM has verbosity as an input or
  // output already. Otherwise it's passed as the last parameter.
  bool pass_verbosity = true;

  // Find out if the WSM has the verbosity as input.
  for (Index j = 0; pass_verbosity && j < mdd.In().nelem(); j++)
  {
    if (wsv_data[mdd.In()[j]].Name() == "verbosity")
    {
      pass_verbosity = false;
    }
  }

  // Find out if the WSM has the verbosity as output.
  for (Index j = 0; pass_verbosity && j < mdd.Out().nelem(); j++)
  {
    if (wsv_data[mdd.Out()[j]].Name() == "verbosity")
    {
      pass_verbosity = false;
    }
  }

  if (pass_verbosity)
  {
    align(ofs,is_first_parameter,indent);
    ofs << "// Verbosity object:\n";
    ofs << indent;
    ofs << "const Verbosity& verbosity";
  }

  ofs << ");\n\n";
}


bool md_sanity_checks (const Array<MdRecord>& md_data)
{
  ostringstream os;

  bool is_sane = true;
  for (Array<MdRecord>::const_iterator i = md_data.begin ();
       i < md_data.end (); ++i)
    {
      bool invalid_author = false;
      for (ArrayOfString::const_iterator j = i->Authors ().begin ();
           !invalid_author && j < i->Authors ().end (); ++j)
        {
          if (*j == "" || *j == "unknown")
            invalid_author = true;
        }

      if (invalid_author)
        {
          os << i->Name () << ": Missing or invalid author.\n";
          is_sane = false;
        }

      switch (check_newline (i->Description ())) {
        case 1:
          os << i->Name () << ": Empty description.\n";
          is_sane = false;
          break;
        case 2:
          os << i->Name () << ": Missing newline at the end of description.\n";
          is_sane = false;
          break;
        case 3:
          os << i->Name () << ": Extra newline at the end of description.\n";
          is_sane = false;
          break;
      }
    }

  if (!is_sane)
    {
      cerr << "Error(s) found in workspace method documentation (check methods.cc):\n"
        << os.str ();
    }

  return is_sane;
}


int main()
{
  try
    {
      // Make the global data visible:
      using global_data::md_data_raw;
      using global_data::md_data;

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

      if (!md_sanity_checks (md_data))
        return 1;

      const Index n_md  = md_data.nelem();

      // Write auto_md.h:
      // -----------
      ofstream ofs;
      open_output_file(ofs,"auto_md.h");

      ofs << "// This file was generated automatically by make_auto_md_h.cc.\n";
      ofs << "// DO NOT EDIT !\n";
      ofs << "// Generated: "
          << __DATE__ << ", "
          << __TIME__ << "\n\n";

      ofs << "#ifndef auto_md_h\n";
      ofs << "#define auto_md_h\n\n";

      ofs << "#include \"matpackI.h\"\n"
          << "#include \"matpackII.h\"\n"
          << "#include \"abs_species_tags.h\"\n"
          << "#include \"gas_abs_lookup.h\"\n"
          << "#include \"gridded_fields.h\"\n"
          << "#include \"optproperties.h\"\n"
          << "#include \"jacobian.h\"\n"
          << "#include \"mc_antenna.h\"\n"
          << "#include \"m_general.h\"\n"
          << "#include \"parser.h\"\n"
          << "#include \"workspace_ng.h\"\n"
          << "#include \"cia.h\"\n"
          << "#include \"covariance_matrix.h\"\n"
          << "#include \"propagationmatrix.h\"\n"
          << "#include \"transmissionmatrix.h\"\n"
          << "#include \"telsem.h\"\n"
          << "#include \"tessem.h\"\n"
          << "#include \"hitran_xsec.h\"\n"
          << "\n";

      ofs << "// This is only used for a consistency check. You can get the\n"
          << "// number of WSMs from md_data.nelem().\n"
          << "#define N_MD " << n_md << "\n\n";


      // Add all the method function declarations
      ofs << "// Method function declarations:\n\n";
      for (Index i=0; i<n_md; ++i)
        {
          const MdRecord& mdd = md_data[i];
          if ( !mdd.UsesTemplates() )
            {
              write_method_header_documentation( ofs, mdd );
              write_method_header( ofs, mdd );
            }
        }

      // Add all the method function declarations
      ofs << "// Supergeneric template function declarations:\n\n";
      for (Index i=0; i<md_data_raw.nelem (); ++i)
        {
          const MdRecord& mdd = md_data_raw[i];
          if ( mdd.Supergeneric() && mdd.UsesTemplates() )
            {
              write_method_header_documentation( ofs, mdd );
              write_method_header( ofs, mdd );
            }
        }

      // Add all the get-away function declarations:
      ofs << "// Get-away function declarations:\n\n";
      for (Index i=0; i<n_md; ++i)
        {
          const MdRecord& mdd = md_data[i];
          if ( mdd.Supergeneric() )
                {
                  ofs << "void " << mdd.Name()
                      << "_sg_" << mdd.ActualGroups()
                      << "_g(Workspace& ws, const MRecord& mr);\n";
                }
              else
                {
                  ofs << "void " << mdd.Name()
                      << "_g(Workspace& ws, const MRecord& mr);\n";
                }
        }

      ofs << "\n";

      // Create prototypes for the agenda wrappers

      // Initialize agenda data.
      Workspace::define_wsv_map ();
      define_agenda_data ();

      using global_data::agenda_data;
      const Array<WsvRecord>& wsv_data = Workspace::wsv_data;
      for (Index i = 0; i < agenda_data.nelem (); i++)
        {
          bool is_agenda_array =
                  wsv_data[get_wsv_id(agenda_data[i].Name())].Group() ==
                      get_wsv_group_id("ArrayOfAgenda");
            write_agenda_wrapper_header (ofs, agenda_data[i], is_agenda_array);

          ofs << ";\n\n";
        }

      ofs << "\n#endif  // auto_md_h\n";

      // Close auto_md.h.
      ofs.close();

    }
  catch (const std::runtime_error &x)
    {
      cout << "Something went wrong. Message text:\n";
      cout << x.what() << '\n';
      return 1;
    }

  return 0;
}
