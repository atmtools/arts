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

/*!
  \file   make_md_h.cc
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
     void get_away_example_g(WorkSpace& ws,
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
#include "vecmat.h"
#include "file.h"
#include "auto_wsv.h"
#include "methods.h"


/* Adds commas and indentation to parameter lists. */
void align(ofstream& ofs, bool& is_first_parameter, const string& indent)
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

int main()
{
  try
    {
      // Make the global data visible:
      extern ARRAY<MdRecord> md_data;
      extern const ARRAY<string> wsv_group_names;
      extern const ARRAY<WsvRecord> wsv_data;

      // Initialize method data.
      define_md_data();

      // Initialize the wsv group name array:
      define_wsv_group_names();

      // Initialize wsv data.
      define_wsv_data();
  

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

      // Write auto_md.h:
      // -----------
      ofstream ofs;
      open_output_file(ofs,"auto_md.h");

      ofs << "// This file was generated automatically by make_md_h.cc.\n";
      ofs << "// DO NOT EDIT !\n";
      ofs << "// Generated: "
	  << __DATE__ << ", "
	  << __TIME__ << "\n\n";

      ofs << "#ifndef auto_md_h\n";
      ofs << "#define auto_md_h\n\n";

      ofs << "#include \"vecmat.h\"\n"
	  << "#include \"auto_wsv.h\"\n"
	  << "#include \"parser.h\"\n"
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

	  // This is needed to flag the first function parameter, which 
	  // needs no line break before being written:
	  bool is_first_parameter = true;

	  // The string indent is needed to achieve the correct
	  // indentation of the functin parameters:
	  string indent(md_data[i].Name().size()+6,' ');

	  // There are four lists of parameters that we have to
	  // write. 
	  ARRAY<size_t>  vo=md_data[i].Output();   // Output 
	  ARRAY<size_t>  vi=md_data[i].Input();    // Input
	  ARRAY<size_t>  vgo=md_data[i].GOutput();   // Generic Output 
	  ARRAY<size_t>  vgi=md_data[i].GInput();    // Generic Input
	  // vo and vi contain handles of workspace variables, 
	  // vgo and vgi handles of workspace variable groups.

	  // Check, if some workspace variables are in both the
	  // input and the output list, and erase those from the input 
	  // list:
	  for (ARRAY<size_t>::const_iterator j=vo.begin(); j<vo.end(); ++j)
	    for (ARRAY<size_t>::iterator k=vi.begin(); k<vi.end(); ++k)
	      {
		if ( *j == *k )
		  {
		    //		  erase_vector_element(vi,k);
		    k = vi.erase(k) - 1;
		    // We need the -1 here, otherwise due to the
		    // following increment we would miss the element
		    // behind the erased one, which is now at the
		    // position of the erased one.
		  }
	      }

	  // There used to be a similar block here for the generic
	  // input/output variables. However, this was a mistake. For
	  // example, if a method has a vector as generic input and a
	  // vector as generic output, this does not mean that it is
	  // the same vector!


	  // Start with the name of the method:
	  ofs << "void " << md_data[i].Name() << "(";

	  // Write the Output workspace variables:
	  {
	    // Flag first parameter of this sort:
	    bool is_first_of_these = true;

	    for (size_t j=0; j<vo.size(); ++j)
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

		ofs << wsv_group_names[wsv_data[vo[j]].Group()] << "&";
	      }
	  }

	  // Write the Generic output workspace variables:
	  {
	    // Flag first parameter of this sort:
	    bool is_first_of_these = true;

	    for (size_t j=0; j<vgo.size(); ++j)
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

		  ofs << wsv_group_names[md_data[i].GOutput()[j]]   << "&";
	      }
	  }

	  // Write the Generic output workspace variable names:
	  {
	    // Flag first parameter of this sort:
	    bool is_first_of_these = true;

	    for (size_t j=0; j<vgo.size(); ++j)
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

		  ofs << "const string&";
	      }
	  }

	  // Write the Input workspace variables:
	  {
	    // Flag first parameter of this sort.
	    bool is_first_of_these = true;

	    for (size_t j=0; j<vi.size(); ++j)
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
		    << wsv_group_names[wsv_data[vi[j]].Group()] << "&";
	      }
	  }

	  // Write the Generic input workspace variables:
	  {
	    // Flag first parameter of this sort.
	    bool is_first_of_these = true;

	    for (size_t j=0; j<vgi.size(); ++j)
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
		
		ofs << "const "
		    << wsv_group_names[md_data[i].GInput()[j]]   << "&";
	      }
	  }

	  // Write the Generic input workspace variable names:
	  {
	    // Flag first parameter of this sort:
	    bool is_first_of_these = true;

	    for (size_t j=0; j<vgi.size(); ++j)
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

		  ofs << "const string&";
	      }
	  }

	  // Write the control parameters:
	  {
	    // Flag first parameter of this sort.
	    bool is_first_of_these = true;

	    // Number of keyword parameters.
	    size_t n_mr = md_data[i].Keywords().size();

	    for (size_t j=0; j!=n_mr; ++j)
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

		extern string TokValTypeName[];
		ofs << "const " << TokValTypeName[md_data[i].Types()[j]] << "& "
		    << md_data[i].Keywords()[j];
	      }
	  }

	  ofs << ");\n\n";
	}

      // Add all the get-away function declarations:
      ofs << "// Get-away function declarations:\n\n";
      for (size_t i=0; i<n_md; ++i)
	ofs << "void " << md_data[i].Name()
	    << "_g(WorkSpace& ws, const MRecord& mr);\n";

      ofs << "\n";

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
