/* Copyright (C) 2001 Stefan Buehler <sbuehler@uni-bremen.de>

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

#include <stdexcept>
#include <iostream>
#include "arts.h"
#include "matpackI.h"
#include "file.h"
#include "array.h"

int main()
{
  try
    {
      // Initialize the group names.
      define_wsv_group_names();

      // Make the names visible.
      extern const ArrayOfString wsv_group_names;

      const Index n_wsv_groups = wsv_group_names.nelem();

      ofstream ofs;
      open_output_file(ofs,"auto_wsv_groups.h");

      ofs << "/*! \\file  auto_wsv_groups.h\n"
          << "    \\brief Defines the enum type that acts as a\n"
          << "    handle for workspace variables groups.\n\n"

          << "    Also defined here is a special pointer class that can hold\n"
          << "    a pointer to any workspace variable.\n\n"

          << "    This file was generated automatically by make_auto_wsv_groups_h.cc.\n"

          << "    <b>DO NOT EDIT!</b>\n\n"

          << "    \\date "
          << __DATE__ << ", "
          << __TIME__ << " */\n\n";

      ofs << "#ifndef auto_wsv_groups_h\n"
          << "#define auto_wsv_groups_h\n\n";

      ofs << "#include <iostream>\n"
          << "#include \"supergeneric.h\"\n"
          << "#include \"ppath.h\"\n"
          << "#include \"absorption.h\"\n\n"
          << "#include \"gas_abs_lookup.h\"\n\n";
      
      ofs << "// Declare existence of class Agenda. We cannot include agenda.h here,\n"
          << "// because that would generate a dependency loop.\n"
          << "class Agenda;\n";

      ofs << "/*! This is only used for a consistency check. You can get the\n"
          << "    number of groups from wsv_group_names.nelem(). */\n"
          << "#define N_WSV_GROUPS " << n_wsv_groups << "\n\n";

      ofs << "/*! The enum type that identifies wsv groups.\n"
          << "    This is used to group workspace variables of the same type\n"
          << "    together, so that generic methods can operate on any of them. */\n";

      ofs << "enum WsvGroup{\n";
      // Now write the group handles one by one:
      for (Index i=0; i<n_wsv_groups; ++i)
        {
          ofs << (i?",\n":"") << "  " << wsv_group_names[i] << "_";
        }
      ofs << "\n};\n\n";

      
      // Now write the declaration of the WsvP class.

      ofs << "/*! Base class for the different Wsv pointers.\n"
          << "    This contains a virtual function for the\n"
          << "    conversion operator for each group.\n\n"
          << "    \\author Stefan Buehler */\n";
      
      ofs << "class WsvP {\n"
          << "public:\n";
      for (Index i=0; i<n_wsv_groups; ++i)
        {
          ofs << "  virtual operator "
              << wsv_group_names[i]
              << "*(){safety();return NULL;};\n";
        }

      ofs << "\nprivate:\n";

      ofs << "/*! Safety check. This is called by all the virtual conversion\n"
          << "    operators. It just stops the program with an error message. This\n"
          << "    should never happen, because conversion should only be attempted\n"
          << "    to the correct type, for which an overloaded conversion operator\n"
          << "    exists. */\n";

      ofs << "  void safety() {\n"
          << "    cerr << \"Internal error: Tried to convert a WsvP \"\n"
          << "         << \"pointer to the wrong type.\\n\";\n"
          << "    exit(1);\n"
          << "  };\n";

      ofs << "};\n\n";


      ofs << "#endif  // auto_wsv_groups_h\n";      
    }
  catch (runtime_error x)
    {
      cout << "Something went wrong. Message text:\n";
      cout << x.what() << '\n';
      return 1;
    }

  return 0;
}
