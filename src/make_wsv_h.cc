/* This is a little C++ program that generates the file wsv.h from the 
   workspace data. The file wsv.h declares the enum type that is used
   to access the workspace data, so it has to be made sure that the
   two are allways consistent. 

   History:
   SAB 29.07.99 Created.
*/

#include "arts.h"
#include "vecmat.h"
#include "workspace.h"

int main()
{
  try
    {
      // Initialize the wsv data.
      define_wsv_data();

      // The lookup information for the workspace variables. 
      extern ARRAY<WsvRecord> wsv_data;

      const size_t n_wsv = wsv_data.size();

      //      cout << "size = " << wsv_data.size() << '\n';

      ofstream ofs("wsv.h");
      if (!ofs)
	{
	  cout << "Cannot open output file wsv.h.\n";
	  return 1;
	}

      ofs << "// This file was generated automatically by make_wsv_h.cc.\n";
      ofs << "// DO NOT EDIT !\n";
      ofs << "// Generated: "
	  << __DATE__ << ", "
	  << __TIME__ << "\n\n";

      ofs << "#ifndef wsv_h\n";
      ofs << "#define wsv_h\n\n";

      ofs << "#define N_WSV " << n_wsv << "\n\n";

      ofs << "enum WsvHandle{\n";
      for (size_t i=0; i<n_wsv-1; ++i)
	{
	  ofs << "  " << wsv_data[i].Name() << "_,\n";
	}
      ofs << "  " << wsv_data[n_wsv-1].Name() << "_\n";
      ofs << "};\n\n";

      ofs << "#endif  // wsv_h\n";
    }
  catch (exception x)
    {
      cout << "Something went wrong, I don't know what.\n";
      cout << "Message text:\n";
      cout << x.what();
      return 1;
    }

  return 0;
}
