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
  \file   methods_aux.cc
  \brief  Auxiliary material for the workspace
          methods, which used to be in methods.cc. 

  The reason for the separation is that the stuff here hardly ever
  should be changed, whereas methods.cc has to be edited each time a
  new method is added. See methods.h for more documentation.

  \author Stefan Buehler
  \date 2000-06-10 */

#include "arts.h"
#include "make_array.h"
#include "auto_wsv.h"
#include "methods.h"


void define_md_map()
{
  extern const Array<MdRecord> md_data;
  extern std::map<String, size_t> MdMap;

  for ( size_t i=0 ; i<md_data.size() ; ++i)
    {
      MdMap[md_data[i].Name()] = i;
    }
}


ostream& MdRecord::PrintTemplate(ostream& os,
				 bool show_description=true) const
{
  extern const  Array<String> wsv_group_names;

  if (show_description)
    {
      // FIXME: Print description String!
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
  extern const Array<WsvRecord> wsv_data;
  extern const Array<String> wsv_group_names;
  extern const String TokValTypeName[];
  bool first;

  os << "\n*--------------------------------------------------------------*\n"
     << "Workspace method = " << mdr.Name() << 
        "\n----------------------------------------------------------------\n"
     << "\n" << mdr.Description() << "\n" << 
        "\n----------------------------------------------------------------\n";

  //  os << "\n-----\nName = " << mdr.Name() << '\n\n'
  //     << "Description =\n" << mdr.Description() << "\n\n";

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
  os << "\n*--------------------------------------------------------------*\n";

  return os;
}

