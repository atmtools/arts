/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>

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

/* This is a little C++ program that generates the file wsv.h from the 
   workspace data. The file wsv.h declares the enum type that is used
   to access the workspace data, so it has to be made sure that the
   two are allways consistent. 

   History:
   SAB 29.07.99 Created.
*/

#include "arts.h"
#include "vecmat.h"
#include "file.h"
#include "workspace.h"

int main()
{
  try
    {
      // Initialize wsv data.
      define_wsv_data();

      // Make wsv data visible.
      extern const ARRAY<WsvRecord> wsv_data;

      const size_t n_wsv = wsv_data.size();

      //      cout << "size = " << wsv_data.size() << '\n';

      ofstream ofs;
      open_output_file(ofs,"wsv.h");

      ofs << "// This file was generated automatically by make_wsv_h.cc.\n";
      ofs << "// DO NOT EDIT !\n";
      ofs << "// Generated: "
	  << __DATE__ << ", "
	  << __TIME__ << "\n\n";

      ofs << "#ifndef wsv_h\n";
      ofs << "#define wsv_h\n\n";
      
      ofs << "// This is only used for a consistency check. You can get the\n"
	  << "// number of workspace variables from wsv_data.size().\n"
	  << "#define N_WSV " << n_wsv << "\n\n";

      ofs << "enum WsvHandle{\n";
      for (size_t i=0; i<n_wsv-1; ++i)
	{
	  ofs << "  " << wsv_data[i].Name() << "_,\n";
	}
      ofs << "  " << wsv_data[n_wsv-1].Name() << "_\n";
      ofs << "};\n\n";


      ofs << "#endif  // wsv_h\n";
    }
  catch (runtime_error x)
    {
      cout << "Something went wrong. Message text:\n";
      cout << x.what() << '\n';
      return 1;
    }

  return 0;
}
