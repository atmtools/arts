/*-----------------------------------------------------------------------
FILE:      methods_aux.cc

INCLUDES:  This file contains auxiliary material for the workspace
           methods, which used to be in methods.cc. The reason for the
           separation is that the stuff here hardly ever should be
           changed, whereas methods.cc has to be edited each time a
           new method is added. See methods.h for more
	   documentation. 

FUNCTIONS: void define_md_map()
           ostream& MdRecord::PrintTemplate(ostream& os,
	                                    bool show_description=true) const
	   ostream& operator<<(ostream& os, const MdRecord& mdr)

HISTORY:   10.06.2000 Created by Stefan Buehler
-----------------------------------------------------------------------*/

#include "arts.h"
#include "make_array.h"
#include "workspace.h"
#include "methods.h"


void define_md_map()
{
  extern const ARRAY<MdRecord> md_data;
  extern std::map<string, size_t> MdMap;

  for ( size_t i=0 ; i<md_data.size() ; ++i)
    {
      MdMap[md_data[i].Name()] = i;
    }
}


ostream& MdRecord::PrintTemplate(ostream& os,
				 bool show_description=true) const
{
  extern const  ARRAY<string> wsv_group_names;

  if (show_description)
    {
      // FIXME: Print description string!
    }
  
  os << Name();

  // Is this a generic method? -- Then we need round braces.
  if ( 0 != GOutput().size()+GInput().size() )
    {
      // First entry needs to comma before:
      bool first=true;

      os << '(';

      for (size_t i=0; i<GOutput().size(); ++i)
	{
	  if (first) first=false;
	  else os << ",\n";

	  os << wsv_group_names[GOutput()[i]];
	}

      for (size_t i=0; i<GInput().size(); ++i)
	{
	  if (first) first=false;
	  else os << ",\n";

	  os << wsv_group_names[GInput()[i]];
	}

      os << ')';
    }

  // Now the keywords:

  os << '{';

  // Determine the length of the longest keyword:
  size_t maxsize = 0;
  for (size_t i=0; i<Keywords().size(); ++i)
    if ( Keywords()[i].size() > maxsize )
      maxsize = Keywords()[i].size();


  for (size_t i=0; i<Keywords().size(); ++i)
    {
      os << "\t" << setw(maxsize)
	 << Keywords()[i] << " = \n";
    }

  os << '}';

  return os;
}

ostream& operator<<(ostream& os, const MdRecord& mdr)
{
  extern const ARRAY<WsvRecord> wsv_data;
  extern const ARRAY<string> wsv_group_names;
  extern const string TokValTypeName[];
  bool first;

  os << "Name = " << mdr.Name() << '\n'
     << "Description =\n" << mdr.Description() << "\n\n";

  // Output:
  first = true;
  os << "Output = ";
  for ( size_t i=0; i<mdr.Output().size(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << wsv_data[mdr.Output()[i]].Name();
    }
  os << '\n';

  // Input:
  first = true;
  os << "Input = ";
  for ( size_t i=0; i<mdr.Input().size(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << wsv_data[mdr.Input()[i]].Name();
    }
  os << '\n';
      
  // GOutput:
  first = true;
  os << "GOutput = ";
  for ( size_t i=0; i<mdr.GOutput().size(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << wsv_group_names[mdr.GOutput()[i]];
    }
  os << '\n';

  // GInput:
  first = true;
  os << "GInput = ";
  for ( size_t i=0; i<mdr.GInput().size(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << wsv_group_names[mdr.GInput()[i]];
    }
  os << '\n';

  // Keywords:
  first = true;
  os << "Keywords = ";
  for ( size_t i=0; i<mdr.Keywords().size(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << mdr.Keywords()[i];
    }
  os << '\n';

  // Types:
  first = true;
  os << "Types = ";
  for ( size_t i=0; i<mdr.Types().size(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << TokValTypeName[mdr.Types()[i]];
    }
  os << '\n';

  return os;
}

