#include "arts.h"
#include "vecmat.h"
#include "file.h"
#include "wsv_aux.h"

int main()
{
  try
    {
      // We need group names and WSV data:
      extern const ARRAY<string> wsv_group_names;
      extern const ARRAY<WsvRecord> wsv_data;

      // Initialize:
      define_wsv_group_names();
      define_wsv_data();

      const size_t n_wsv = wsv_data.size();

      ofstream ofs;
      open_output_file(ofs,"wsv_pointers.cc");

      ofs << "/*! \\file  wsv_pointers.cc\n"
	  << "    \\brief Defines the smart pointers that are used by\n"
	  << "            the engine to access workspace variables.\n\n"

	  << "    This file was generated automatically by make_wsv_pointers_cc.cc.\n"

	  << "    <b>DO NOT EDIT!</b>\n\n"

	  << "    \\date "
	  << __DATE__ << ", "
	  << __TIME__ << " */\n\n";

      ofs << "#include \"arts.h\"\n"
	  << "#include \"vecmat.h\"\n"
	  << "#include \"wsv_groups.h\"\n"
	  << "#include \"wsv_aux.h\"\n"
	  << "#include \"wsv.h\"\n\n";

      ofs << "ARRAY<WsvP*> wsv_pointers;\n\n";
      
      ofs << "void define_wsv_pointers()\n"
	  << "{\n"
	  << "  /* Make global data visible. */\n"
	  << "  extern WorkSpace workspace;\n\n";

      // Now write the pointers one by one:
      for (size_t i=0; i<n_wsv; ++i)
	{
	  const WsvRecord& wr = wsv_data[i];

	  ofs << "  {\n"
	      << "    static WsvPointer<"
	      << wsv_group_names[wr.Group()]
	      << "> p(&workspace."
	      << wr.Name()
	      << ");\n";

	  ofs << "    wsv_pointers.push_back(&p);\n"
	      << "  }\n\n";
	}

      ofs << "};\n";

    }
  catch (runtime_error x)
    {
      cout << "Something went wrong. Message text:\n";
      cout << x.what() << '\n';
      return 1;
    }

  return 0;
}
