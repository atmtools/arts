#include "arts.h"
#include "vecmat.h"
#include "file.h"

int main()
{
  try
    {
      // Initialize the group names.
      define_wsv_group_names();

      // Make the names visible.
      extern const ARRAY<string> wsv_group_names;

      const size_t n_wsv_groups = wsv_group_names.size();

      ofstream ofs;
      open_output_file(ofs,"wsv_group.h");

      ofs << "/*! \\file  wsv_group.h\n"
	  << "    \\brief Defines the enum type that acts as a\n"
	  << "    handle for workspace variables groups.\n\n"

	  << "    Also defined here is a special pointer class that can hold\n"
	  << "    a pointer to any workspace variable.\n\n"

	  << "    This file was generated automatically by make_wsv_groups_h.cc.\n"

	  << "    <b>DO NOT EDIT!</b>\n\n"

	  << "    \\date "
	  << __DATE__ << ", "
	  << __TIME__ << " */\n\n";

      ofs << "#ifndef wsv_group_h\n"
	  << "#define wsv_group_h\n\n";
      
      ofs << "/*! This is only used for a consistency check. You can get the\n"
	  << "    number of groups from wsv_group_names.size(). */\n"
	  << "#define N_WSV_GROUPS " << n_wsv_groups << "\n\n";

      ofs << "/*! The enum type that identifies wsv groups.\n"
	  << "    This is used to group workspace variables of the same type\n"
	  << "    together, so that generic methods can operate on any of them. */\n";

      ofs << "enum WsvGroup{\n";
      // Now write the group handles one by one:
      for (size_t i=0; i<n_wsv_groups; ++i)
	{
	  ofs << "  " << wsv_group_names[i] << "_,\n";
	}
      ofs << "};\n\n";

      
      // Now write the declaration of the WsvP class.

      ofs << "/*! Base class for the different Wsv pointers.\n"
	  << "    This contains a virtual function for the\n"
	  << "    conversion operator for each group.\n\n"
	  << "    \\author Stefan Buehler */\n";
      
      ofs << "class WsvP {\n"
	  << "public:\n";
      for (size_t i=0; i<n_wsv_groups; ++i)
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


      ofs << "#endif  // wsv_group_h\n";      
    }
  catch (runtime_error x)
    {
      cout << "Something went wrong. Message text:\n";
      cout << x.what() << '\n';
      return 1;
    }

  return 0;
}
