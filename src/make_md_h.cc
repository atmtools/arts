// FIXME: Should I also generate headers for the methods themselves?

/* This is a little C++ program that generates the file md.h from the
   workspace methods data md_data. The file md.h declares the enum
   type MdHandle that is used to access the method data, so it has
   to be made sure that the two are allways consistent.

   A second file is produced: md.cc.
   This contains the `get-away' functions that provided the interface
   between the engine and the workspace methods. The get-functions all 
   have the same arguments:

   void get_away_example_g(WorkSpace& ws,
			   const ARRAY<TokVal>& tv);

   Their names all have the extension _g

   Pointers to the get-away functions are stored in the array
   `getaway'. 

   Each get-away function simply contains a function call to the
   matching workspace method. The parameters are arranged similar to
   the follwing example:
   
   void SomeMethod(owsv1,iwsv1,iwsv2,iwsv3,c1,c2,c3,...)

   First come the output workspace variables, then the input workspace 
   variables, and then the control parameters. There can be an
   arbitrary number of parameters of each type, but the most usual
   case is to have only one output workspace variable. 

   The same variable may be both in the list of input and in the list
   of output workspace variables. This case makes good sense,
   actually, if you think for example of a method that adds an offset
   to the absorption coefficients. IN THAT CASE THE VARIABLE IS ADDED
   TO THE LIST ONLY ONCE, namely among the OUTPUT variables.

   History:
   SAB 29.07.99 Created.  */

#include "arts.h"
#include "token.h"
#include "vecmat.h"
#include "file.h"
#include "wsv.h"
#include "workspace.h"
#include "methods.h"

int main()
{
  try
    {
      // Initialize method data.
      define_md_data();

      // Make method data visible.
      extern ARRAY<MdRecord> md_data;

      // Initialize wsv data.
      define_wsv_data();
  
      // Make wsv data visible.
      extern ARRAY<WsvRecord> wsv_data;

      const size_t n_md  = md_data.size();
      const size_t n_wsv = wsv_data.size();

      // For safety, check if n_wsv and N_WSV have the same value. If not, 
      // then the file wsv.h is not up to date.
      if (N_WSV != n_wsv)
	{
	  cout << "The file wsv.h is not up to date!\n";
	  cout << "(N_WSV = " << N_WSV << ", n_wsv = " << n_wsv << ")\n";
	  cout << "Make wsv.h first. Check if Makefile is correct.\n";
	  return 1;
	}

      // Write md.h:
      // -----------
      ofstream ofs;
      open_output_file(ofs,"md.h");

      ofs << "// This file was generated automatically by make_md_h.cc.\n";
      ofs << "// DO NOT EDIT !\n";
      ofs << "// Generated: "
	  << __DATE__ << ", "
	  << __TIME__ << "\n\n";

      ofs << "#ifndef md_h\n";
      ofs << "#define md_h\n\n";

      ofs << "#include \"vecmat.h\"\n"
	  << "#include \"workspace.h\"\n"
	  << "#include \"token.h\"\n"
	  << "\n";

      ofs << "// This is only used for a consistency check. You can get the\n"
	  << "// number of workspace variables from wsv_data.size().\n"
	  << "#define N_MD " << n_md << "\n\n";

      ofs << "enum MdHandle{\n";
      for (size_t i=0; i<n_md-1; ++i)
	{
	  ofs << "  " << md_data[i].Name() << "_,\n";
	}
      ofs << "  " << md_data[n_md-1].Name() << "_\n";
      ofs << "};\n\n";

      // Add all the method function declarations
      ofs << "// Method function declarations:\n\n";
      for (size_t i=0; i<n_md; ++i)
	{
	  // The names of the types of workspace variables.
	  extern ARRAY<string> wsv_group_names;

	  // This is needed to flag the first function parameter, which 
	  // needs no line break before being written:
	  bool is_first_parameter = true;

	  // The string indent is needed to achieve the correct
	  // indentation of the functin parameters:
	  string indent(md_data[i].Name().size()+6,' ');

	  // There are two lists of parameters that we have to
	  // write: 
	  ARRAY<size_t>  vo=md_data[i].Output();   // Output 
	  ARRAY<size_t>  vi=md_data[i].Input();    // Input

	  // Check, if some workspace variables are in both the
	  // input and the output list, and erase those from the input 
	  // list:
	  for (ARRAY<size_t>::const_iterator j=vo.begin(); j!=vo.end(); ++j)
	    // It is important that the condition is k<vi.end(), not
	    // k!=vi.end, because if erase is called, vi.end() is
	    // decreased. Since k is increased at the same time, the
	    // case k=vi.end() can be missed!
	    for (ARRAY<size_t>::iterator k=vi.begin(); k<vi.end(); ++k)
	      if ( *j == *k )
		{
		  vi.erase(k);
		}

	  ofs << "void " << md_data[i].Name() << "(";

	  // Write the Output workspace variables:
	  {
	    // Flag first parameter of this sort:
	    bool is_first_of_these = true;

	    for (ARRAY<size_t>::const_iterator j=vo.begin(); j!=vo.end(); ++j)
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

		// Add comment if this is the first of this sort
		if (is_first_of_these)
		  {
		    ofs << "// WS Output:\n";
		    ofs << indent;
		    is_first_of_these = false;
		  }

		ofs << wsv_group_names[wsv_data[*j].Group()] << "&";
	      }
	  }

	  // Write the Input workspace variables:
	  {
	    // Flag first parameter of this sort.
	    bool is_first_of_these = true;

	    for (ARRAY<size_t>::const_iterator j=vi.begin(); j!=vi.end(); ++j)
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
		    
		// Add type if this is the first of this sort.
		if (is_first_of_these)
		  {
		    ofs << "// WS Input:\n";
		    ofs << indent;		  
		    is_first_of_these = false;
		  }

		ofs << "const ws." << wsv_group_names[wsv_data[*j].Group()] << "&";
	      }
	  }

	  // Write the control parameters:
	  {
	    // Flag first parameter of this sort.
	    	    bool is_first_of_these = true;

	    // Number of keyword parameters.
	    size_t n_tv = md_data[i].Keywords().size();

	    for (size_t j=0; j!=n_tv; ++j)
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
		    
		// Add type if this is the first of this sort.
		if (is_first_of_these)
		  {
		    ofs << "// Control Parameters:\n";
		    ofs << indent;		  
		    is_first_of_these = false;
		  }

		extern string TokValTypeName[];
		ofs << "const " << TokValTypeName[md_data[i].Types()[j]] << "&";
	      }
	  }

	  ofs << ");\n\n";
	}

      // Add all the get-away function declarations:
      ofs << "// Get-away function declarations:\n\n";
      for (size_t i=0; i<n_md; ++i)
	ofs << "void " << md_data[i].Name()
	    << "_g(WorkSpace& ws, const ARRAY<TokVal>& tv);\n";

      ofs << "\n";

      ofs << "\n#endif  // md_h\n";

      // Close md.h.
      ofs.close();

      // Write md.cc:
      // -----------
      open_output_file(ofs,"md.cc");
  
      ofs << "// This file was generated automatically by make_md_h.cc.\n";
      ofs << "// DO NOT EDIT !\n";
      ofs << "// Generated: "
	  << __DATE__ << ", "
	  << __TIME__ << "\n\n";

      ofs << "#include \"arts.h\"\n"
	  << "#include \"md.h\"\n"
	  << "\n";

      // Write all get-away functions:
      // -----------------------------
      for (size_t i=0; i<n_md; ++i)
	{
	  // This is needed to flag the first function parameter, which 
	  // needs no line break before being written:
	  bool is_first_parameter = true;
	  // The string indent is needed to achieve the correct
	  // indentation of the functin parameters:
	  string indent(md_data[i].Name().size()+3,' ');

	  
	  // There are two lists of parameters that we have to
	  // write: 
	  ARRAY<size_t>  vo=md_data[i].Output();   // Output 
	  ARRAY<size_t>  vi=md_data[i].Input();    // Input

	  // Check, if some workspace variables are in both the
	  // input and the output list, and erase those from the input 
	  // list:
	  for (ARRAY<size_t>::const_iterator j=vo.begin(); j!=vo.end(); ++j)
	    // It is important that the condition is k<vi.end(), not
	    // k!=vi.end, because if erase is called, vi.end() is
	    // decreased. Since k is increased at the same time, the
	    // case k=vi.end() can be missed!
	    for (ARRAY<size_t>::iterator k=vi.begin(); k<vi.end(); ++k)
	      if ( *j == *k )
		{
		  vi.erase(k);
		}

	  ofs << "void " << md_data[i].Name()
	      << "_g(WorkSpace& ws, const ARRAY<TokVal>& tv)\n";
	  ofs << "{\n";
	  ofs << "  " << md_data[i].Name() << "(";

	  // Write the Output workspace variables:
	  {
	    // Flag first parameter of this type:
	    // 	    bool is_first_of_these = true;

	    for (ARRAY<size_t>::const_iterator j=vo.begin(); j!=vo.end(); ++j)
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

		// Add type if this is the first of this type
		// 		if (is_first_of_these)
		// 		  {
		// 		    ofs << "// WS Output:\n";
		// 		    ofs << indent;
		// 		    is_first_of_these = false;
		// 		  }
		    
		ofs << "ws." << wsv_data[*j].Name();
	      }
	  }

	  // Write the Input workspace variables:
	  {
	    // Flag first parameter of this type:
	    // 	    bool is_first_of_these = true;

	    for (ARRAY<size_t>::const_iterator j=vi.begin(); j!=vi.end(); ++j)
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
		    
		// Add type if this is the first of this type
		// 		if (is_first_of_these)
		// 		  {
		// 		    ofs << "// WS Input:\n";
		// 		    ofs << indent;		  
		// 		    is_first_of_these = false;
		// 		  }

		ofs << "ws." << wsv_data[*j].Name();
	      }
	  }

	  // Write the control parameters:
	  {
	    // Flag first parameter of this type:
	    // 	    bool is_first_of_these = true;

	    // The tv parameters look all the same (tv[i]), so we just
	    // need to know the number of them: 
	    size_t n_tv = md_data[i].Keywords().size();

	    for (size_t j=0; j!=n_tv; ++j)
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
		    
		// Add type if this is the first of this type
		// 		if (is_first_of_these)
		// 		  {
		// 		    ofs << "// Control Parameters:\n";
		// 		    ofs << indent;		  
		// 		    is_first_of_these = false;
		// 		  }

		ofs << "tv[" << j << "]";
	      }
	  }

	  ofs << ");\n";
	  ofs << "}\n\n";
	}

      // Add getaway, the array that hold pointers to the getaway functions:
      {
	string indent = "     ";
	bool is_first_parameter = true;

	ofs << "/** The array holding the pointers to the getaway functions. */"
	    << "void (*getaways[])(WorkSpace&, const ARRAY<TokVal>&)\n"
	    << "  = {";
	for (size_t i=0; i<n_md; ++i)
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
	    ofs << md_data[i].Name() << "_g";
	  }
	ofs << "};\n\n";
      }

    }
  catch (exception x)
    {
      cout << "Something went wrong. Message text:\n";
      cout << x.what() << '\n';
      return 1;
    }

  return 0;
}
