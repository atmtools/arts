/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>

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
  \file   make_auto_wsv_h.cc
  \brief  Generates the file auto_wsv.h.

  This is a little C++ program that generates the header file from the
  workspace data in workspace.cc. The file auto_wsv.h declares the enum
  type that is used to access the workspace data.

  \author Stefan Buehler
  \date 1999-07-29 */


#include <iostream>
#include "arts.h"
#include "matpackI.h"
#include "array.h"
#include "file.h"
#include "absorption.h"
#include "auto_wsv_groups.h"
#include "wsv_aux.h"
#include "mystring.h"

int main()
{
  try
    {
      // Initialize wsv data and wsv group names.
      define_wsv_data();
      define_wsv_group_names();

      // Make the data visible.
      extern const Array<WsvRecord> wsv_data;
      extern const ArrayOfString wsv_group_names;

      const Index n_wsv = wsv_data.nelem();

      //      cout << "size = " << wsv_data.nelem() << '\n';

      ofstream ofs,ofs2;
      open_output_file(ofs,"auto_wsv.h");
      open_output_file(ofs2,"auto_wsv.txt");

      ofs << "/*! \\file  auto_wsv.h\n"
          << "    \\brief Declares the enum type that acts as a\n"
          << "    handle for workspace variables. Also declares the\n"
          << "    workspace itself.\n\n"

          << "    This file was generated automatically by make_auto_wsv_h.cc.\n"

          << "    <b>DO NOT EDIT!</b>\n\n"

          << "    \\date "
          << __DATE__ << ", "
          << __TIME__ << " */\n\n";

      ofs << "#ifndef auto_wsv_h\n";
      ofs << "#define auto_wsv_h\n\n";

      ofs << "#include \"absorption.h\"\n"
          << "#include \"agenda_class.h\"\n"
          << "#include \"ppath.h\"\n"
          << "#include \"matpackIII.h\"\n"
          << "#include \"gas_abs_lookup.h\"\n"
          << "\n";
      
      ofs << "//! Total number of workspace variables.\n"
          << "/*! \n"
          << "  This is mainly used for a consistency check. You can get the\n"
          << "  number of workspace variables from wsv_data.nelem(). \n"
          << "  \n"
          << "  However, this is used by member functions of Workspace. For\n"
          << "  example by the default constructor to initialize the occupied\n"
          << "  array. Also for range checking in is_occupied, set, and free\n"
          << "  member functions. \n"
          << "*/\n"
          << "#define N_WSV " << n_wsv << "\n\n";

      ofs << "enum WsvHandle{\n";
      for (Index i=0; i<n_wsv-1; ++i)
        {
          ofs << "  " << wsv_data[i].Name() << "_,\n";
        }
      ofs << "  " << wsv_data[n_wsv-1].Name() << "_\n";
      ofs << "};\n\n";

      // Now the workspace itself:

      ofs << "/** The declaration of the (great) workspace. */\n"
          << "class WorkSpace {\n"
          << "public:\n"
          << "  WorkSpace();\n"
          << "  bool is_occupied(Index i) const;\n"
          << "  void set(Index i);\n"
          << "  void free(Index i);\n"
          << "private:\n"
          << "  /** Keeps track which elements are occupied.\n"
          << "      For some strange reason, Array<bool> does not work, hence\n"
          << "      this is implemented internally as an Array of Index.*/\n"     
          << "  Array<Index> moccupied;\n"
          << "public:\n";
      
      for (Index i=0; i<n_wsv; ++i)
        {
          // First of all write the comment as a doxygen header.
          // For this we have to make some small replacements for
          // indendation and put everything starting from the second
          // sentence into a verbatim environment.  
          {
            // Local copy of the description String:
            String s = wsv_data[i].Description();

            // Add indentation:
            replace_all(s,"\n","\n    "); 

            // Look for the end of the first sentence. There we have
            // to include the beginning of the verbatim
            // environment. Not earlier, because the first sentence
            // has a special meaning.
            Index full_stop = s.find('.');

            // We have to check against the case that the point was
            // not found. In that case we set full_stop to the length
            // of the entire test. If it was found, we increase the
            // value by one, since we want to include the point in the
            // first part.
            if ( full_stop==s.npos )
              full_stop = s.nelem();
            else
              full_stop += 1;

            String first(s,0,full_stop);

            // Rest will contain the rest of the documentation. It could
            // be empty!
            String rest = "";

            if ( full_stop!=s.nelem() )
              rest = String(s,full_stop);

            // Remove leading whitespace and linebreaks in rest:
            while (
                   0 < rest.nelem() &&
                   ( ' ' == rest[0] || '\n' == rest[0] )
                   )
              {
                rest.erase(0,1);
              }

            ofs << "/** " << first;
            if ( 0==rest.nelem() )
              {
                ofs << " */\n";
              }
            else
              {
                ofs << '\n'
                    << "    \\verbatim\n"
                    << "    " << rest << '\n'
                    << "    \\endverbatim */\n";
              }
          }

          ofs << "  "
              << wsv_group_names[wsv_data[i].Group()]
              << " "
              << wsv_data[i].Name() << ";\n";

        }
      ofs << "};\n\n";

      ofs << "#endif  // auto_wsv_h\n";

      // Write text information file
      // If you make changes here for the "VARIABLE" and "DATA TYPE" rows,
      // you probably need to change get_artstype.m in AMI.
      for ( Index i=0; i<n_wsv; i++ )
      {
        ofs2 << "VARIABLE : " << wsv_data[i].Name() << "\n"
             << "DATA TYPE: " << wsv_group_names[wsv_data[i].Group()] <<"\n"
             << "DESCRIPTION:\n" 
             << wsv_data[i].Description() << "\n\n";
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

