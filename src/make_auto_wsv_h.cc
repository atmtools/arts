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
          << "#include \"matpackII.h\"\n"
          << "#include \"matpackIII.h\"\n"
          << "#include \"gas_abs_lookup.h\"\n"
          << "#include \"optproperties.h\"\n"
          << "#include \"m_general.h\"\n"
          << "#include \"gridded_fields.h\"\n"
          << "#include \"jacobian.h\"\n"
          << "#include \"mc_interp.h\"\n"
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

      ////////////////////////////////////////////////////////////////////
      // WorkspaceMemoryHandler class
      //
      ofs << "class WorkspaceMemoryHandler {\n"
        <<   "private:\n"
        <<   "  // List of function pointers to allocation routines\n"
        <<   "  void *(*allocfp[" << wsv_data.nelem () << "])();\n"
        <<   "  // List of function pointers to deallocation routines\n"
        <<   "  void (*deallocfp[" << wsv_data.nelem () << "])(void *);\n\n"
        <<   "  // Allocation and deallocation routines for workspace groups\n";
      for (Index i = 0; i < wsv_group_names.nelem (); ++i)
        {
          ofs << "  static void *allocate_wsvg_" << wsv_group_names[i] << "()\n"
            <<   "    { return (void *)new " << wsv_group_names[i] << "; }\n\n"
            <<   "  static void deallocate_wsvg_" << wsv_group_names[i] << "(void *vp)\n"
            <<   "    { delete (" << wsv_group_names[i] << " *)vp; }\n\n";
        }

      ofs << "public:\n"
        <<   "  /** Default constructor. Initialize allocation and "
        <<   "deallocation\n"
        <<   "      function pointer lists.\n"
        <<   "  */\n"
        <<   "  WorkspaceMemoryHandler ()\n"
        <<   "    {\n";

      for (Index i = 0; i < wsv_data.nelem (); ++i)
        {
          ofs << "      allocfp[" << i << "] = allocate_wsvg_"
            <<            wsv_group_names[wsv_data[i].Group()] << ";\n"
            <<   "      deallocfp[" << i << "] = deallocate_wsvg_"
            <<            wsv_group_names[wsv_data[i].Group()] << ";\n";
        }

      ofs << "    }\n\n"
        <<   "  /** Getaway function to call the allocation function for the\n"
        <<   "      workspace variable with the given Index.\n"
        <<   "  */\n"
        <<   "  void *allocate (Index wsv)\n"
        <<   "    {\n"
        <<   "      return allocfp[wsv]();\n"
        <<   "    }\n\n"
        <<   "  /** Getaway function to call the deallocation function for the\n"
        <<   "      workspace variable with the given Index.\n"
        <<   "  */\n"
        <<   "  void deallocate (Index wsv, void *vp)\n"
        <<   "    {\n"
        <<   "      deallocfp[wsv](vp);\n"
        <<   "    }\n\n";

      ofs << "};\n\n";
      //
      ////////////////////////////////////////////////////////////////////

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

