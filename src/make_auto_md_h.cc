/* Copyright (C) 2000-2007 Stefan Buehler <sbuehler@ltu.se>

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
#include "token.h"
#include "array.h"
#include "file.h"
#include "auto_wsv.h"
#include "methods.h"
#include "wsv_aux.h"
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

//! Write method header documentation.
/*!
  \param[in,out] ofs The stream to write to.
  \param mdd Method lookup data.
*/
void write_method_header_documentation (ofstream& ofs, const MdRecord& mdd)
{
  extern const Array<WsvRecord> wsv_data;

  String fullname = mdd.Name();

  // This is needed to flag the first function parameter, which 
  // needs no line break before being written:
  bool is_first_parameter = true;

  // The String indent is needed to achieve the correct
  // indentation of the functin parameters:
  String indent("  ");

  // There are four lists of parameters that we have to
  // write. 
  ArrayOfIndex  vo=mdd.Output();   // Output 
  ArrayOfIndex  vi;                // Input
  ArrayOfIndex  vgo=mdd.GOutput(); // Generic Output 
  ArrayOfIndex  vgi=mdd.GInput();  // Generic Input
  // vo and vi contain handles of workspace variables, 
  // vgo and vgi handles of workspace variable groups.

  mdd.input_only(vi);

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

  if ( mdd.PassWorkspace() )
    {
      ofs << indent << "\\param[in,out] " << "ws Workspace\n";
    }

  // Write the Output workspace variables:
  for (Index j=0; j<vo.nelem(); ++j)
    {
      ofs << indent << "\\param[out] "
        << wsv_data[vo[j]].Name() << " WS Output\n";
    }

  // Write the Generic output workspace variables:
  for (Index j=0; j<vgo.nelem(); ++j)
    {
      if (mdd.Supergeneric ())
        {
          ofs << indent << "\\param[out] supergenericoutput" << j+1
            << " Supergeneric output\n";
        }
      else
        {
          ofs << indent << "\\param[out] genericoutput" << j+1
            << " Generic output\n";
        }
    }

  // Write the Input workspace variables:
  for (Index j=0; j<vi.nelem(); ++j)
    {
      ofs << indent << "\\param[in] "
        << wsv_data[vi[j]].Name() << " WS Input\n";
    }

  // Write the Generic input workspace variables:
  for (Index j=0; j<vgi.nelem(); ++j)
    {
      if (mdd.Supergeneric ())
        {
          ofs << indent << "\\param[in] supergenericinput" << j+1
            << " Supergeneric Input\n";
        }
      else
        {
          ofs << indent << "\\param[in] genericinput" << j+1
            << " Generic Input\n";
        }
    }

  // Write the control parameters:
  {
    // Number of keyword parameters.
    Index n_mr = mdd.Keywords().nelem();

    for (Index j=0; j!=n_mr; ++j)
      {
        ofs << indent << "\\param[in] "
          << mdd.Keywords()[j] << " Control parameter";

        if (mdd.Defaults()[j] != NODEF)
          {
            ofs << " (Default: " << mdd.Defaults()[j] << ")";
          }
        
        ofs << "\n";
      }
  }

  // Write agenda, if there is one:
  if ( mdd.AgendaMethod() )
    {
      align(ofs,is_first_parameter,indent);
      ofs << indent << "\\param[in] " << "input_agenda Agenda from controlfile\n";
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
  extern const ArrayOfString wsv_group_names;
  extern const Array<WsvRecord> wsv_data;

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

  // There are four lists of parameters that we have to
  // write. 
  ArrayOfIndex  vo=mdd.Output();   // Output 
  ArrayOfIndex  vi;                // Input
  ArrayOfIndex  vgo=mdd.GOutput(); // Generic Output 
  ArrayOfIndex  vgi=mdd.GInput();  // Generic Input
  // vo and vi contain handles of workspace variables, 
  // vgo and vgi handles of workspace variable groups.

  mdd.input_only(vi);

  // There used to be a similar block here for the generic
  // input/output variables. However, this was a mistake. For
  // example, if a method has a vector as generic input and a
  // vector as generic output, this does not mean that it is
  // the same vector!

  if (mdd.Supergeneric ())
    {
      ofs << "template <typename T>" << endl;
    }

  // Start with the name of the method:
  ofs << "void " << fullname << "(";

  if (mdd.PassWorkspace ())
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

        ofs << wsv_group_names[wsv_data[vo[j]].Group()] << "& "
          << wsv_data[vo[j]].Name();
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

        if (wsv_group_names[mdd.GOutput()[j]] == "Any")
          {
            ofs << "T& supergenericoutput" << j+1;
          }
        else
          {
            ofs << wsv_group_names[mdd.GOutput()[j]] << "& genericoutput"
              << j+1;
          }
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

          ofs << "const String& genericoutputname" << j+1;
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
          << wsv_group_names[wsv_data[vi[j]].Group()] << "& "
          << wsv_data[vi[j]].Name();
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
                
        if (wsv_group_names[mdd.GInput()[j]] == "Any")
          {
            ofs << "const T& supergenericinput" << j+1;
          }
        else
          {
            ofs << "const "
              << wsv_group_names[mdd.GInput()[j]]
              << "& genericinput" << j+1;
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

          ofs << "const String& genericinputname" << j+1;
        }
    }

  // Write the control parameters:
  {
    // Flag first parameter of this sort.
    bool is_first_of_these = true;

    // Number of keyword parameters.
    Index n_mr = mdd.Keywords().nelem();

    for (Index j=0; j!=n_mr; ++j)
      {
        // Add comma and line break, if not first element:
        align(ofs,is_first_parameter,indent);
                    
        // Add type if this is the first of this sort.
        if (is_first_of_these)
          {
            ofs << "// Control Parameters:\n";
            ofs << indent;                
            is_first_of_these = false;
          }

        extern String TokValTypeName[];
        ofs << "const " << TokValTypeName[mdd.Types()[j]] << "& "
            << mdd.Keywords()[j];
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
      extern Array<MdRecord> md_data_raw;
      extern Array<MdRecord> md_data;
      extern const ArrayOfString wsv_group_names;
      extern const Array<WsvRecord> wsv_data;

      // Initialize method data.
      define_md_data_raw();

      // Initialize the wsv group name array:
      define_wsv_group_names();

      // Expand supergeneric methods:
      expand_md_data_raw_to_md_data();

        // Initialize wsv data.
      define_wsv_data();

      if (!md_sanity_checks (md_data))
        return 1;

      const Index n_md  = md_data.nelem();
      const Index n_wsv = wsv_data.nelem();

      // For safety, check if n_wsv and N_WSV have the same value. If not, 
      // then the file wsv.h is not up to date.
      if (N_WSV != n_wsv)
        {
          cout << "The file wsv.h is not up to date!\n";
          cout << "(N_WSV = " << N_WSV << ", n_wsv = " << n_wsv << ")\n";
          cout << "Make wsv.h first. Check if Makefile is correct.\n";
          return 1;
        }

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
          << "#include \"auto_wsv.h\"\n"
          << "#include \"parser.h\"\n"
          << "#include \"workspace_ng.h\"\n"
          << "\n";

      ofs << "// This is only used for a consistency check. You can get the\n"
          << "// number of WSMs from md_data.nelem().\n"
          << "#define N_MD " << n_md << "\n\n";

     
      // We don't really need these handles, do we?

//       ofs << "enum MdHandle{\n";
//       for (Index i=0; i<n_md-1; ++i)
//      {
//        ofs << "  " << md_data[i].Name() << "_,\n";
//      }
//       ofs << "  " << md_data[n_md-1].Name() << "_\n";
//       ofs << "};\n\n";

      // Add all the method function declarations
      ofs << "// Method function declarations:\n\n";
      for (Index i=0; i<n_md; ++i)
        {
          const MdRecord& mdd = md_data[i];
          if ( !mdd.SuppressHeader() )
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
          if ( mdd.Supergeneric() )
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
                      << "_sg_" << wsv_group_names[mdd.ActualGroup()]
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
      define_wsv_map ();
      define_agenda_data ();

      extern const Array<AgRecord> agenda_data;
      for (Index i = 0; i < agenda_data.nelem (); i++)
        {
          write_agenda_wrapper_header (ofs, agenda_data[i], true);

          ofs << ";\n\n";
        }

      ofs << "\n#endif  // auto_md_h\n";

      // Close auto_md.h.
      ofs.close();

    }
  catch (exception x)
    {
      cout << "Something went wrong. Message text:\n";
      cout << x.what() << '\n';
      return 1;
    }

  return 0;
}
