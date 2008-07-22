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
#include "array.h"
#include "file.h"
#include "workspace_ng.h"
#include "mystring.h"

bool wsv_sanity_checks (const Array<WsvRecord>& wsv_data)
{
  ostringstream os;

  bool is_sane = true;
  for (Array<WsvRecord>::const_iterator i = wsv_data.begin ();
       i < wsv_data.end (); ++i)
    {
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
      cerr << "Error(s) found in workspace variable documentation (check methods.cc):\n"
        << os.str ();
    }

  return is_sane;
}


int main()
{
  try
    {

      // Initialize the wsv group name array:
      define_wsv_group_names();

      // Initialize wsv data and wsv group names.
      Workspace::define_wsv_data();

      // Make the data visible.
      const Array<WsvRecord>& wsv_data = Workspace::wsv_data;

      if (!wsv_sanity_checks (wsv_data))
        return 1;

      const Index n_wsv = wsv_data.nelem();

      //      cout << "size = " << wsv_data.nelem() << '\n';

      ofstream ofs,ofs2;
      open_output_file(ofs,"auto_wsv.h");

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

      ofs << "#include \"agenda_class.h\"\n"
          << "#include \"ppath.h\"\n"
          << "#include \"matpackII.h\"\n"
          << "#include \"matpackIII.h\"\n"
          << "#include \"gas_abs_lookup.h\"\n"
          << "#include \"optproperties.h\"\n"
          << "#include \"m_general.h\"\n"
          << "#include \"gridded_fields.h\"\n"
          << "#include \"jacobian.h\"\n"
          << "#include \"mc_interp.h\"\n"
          << "#include \"mc_antenna.h\"\n"
          << "#include \"supergeneric.h\"\n"
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

      ofs << "#endif  // auto_wsv_h\n";

    }
  catch (runtime_error x)
    {
      cout << "Something went wrong. Message text:\n";
      cout << x.what() << '\n';
      return 1;
    }

  return 0;
}

