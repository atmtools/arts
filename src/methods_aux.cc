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
  \file   methods_aux.cc
  \brief  Auxiliary material for the workspace
          methods, which used to be in methods.cc. 

  The reason for the separation is that the stuff here hardly ever
  should be changed, whereas methods.cc has to be edited each time a
  new method is added. See methods.h for more documentation.

  \author Stefan Buehler
  \date 2000-06-10 */

#include "arts.h"
#include <map>
#include "make_array.h"
#include "auto_wsv.h"
#include "methods.h"
#include "wsv_aux.h"

void define_md_map()
{
  extern const Array<MdRecord> md_data;
  extern std::map<String, Index> MdMap;

  for ( Index i=0 ; i<md_data.nelem() ; ++i)
    {
      MdMap[md_data[i].Name()] = i;
    }
}


ostream& MdRecord::PrintTemplate(ostream& os,
                                 bool show_description) const
{
  extern const  ArrayOfString wsv_group_names;

  if (show_description)
    {
      // FIXME: Print description String!
    }
  
  os << Name();

  // Is this a generic method? -- Then we need round braces.
  if ( 0 != GOutput().nelem()+GInput().nelem() )
    {
      // First entry needs to comma before:
      bool first=true;

      os << '(';

      for (Index i=0; i<GOutput().nelem(); ++i)
        {
          if (first) first=false;
          else os << ",\n";

          os << wsv_group_names[GOutput()[i]];
        }

      for (Index i=0; i<GInput().nelem(); ++i)
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
  Index maxsize = 0;
  for (Index i=0; i<Keywords().nelem(); ++i)
    if ( Keywords()[i].nelem() > maxsize )
      maxsize = Keywords()[i].nelem();


  for (Index i=0; i<Keywords().nelem(); ++i)
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
  extern const ArrayOfString wsv_group_names;
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
  for ( Index i=0; i<mdr.Output().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << wsv_data[mdr.Output()[i]].Name();
    }
  os << '\n';

  // Input:
  first = true;
  os << "Input = ";
  for ( Index i=0; i<mdr.Input().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << wsv_data[mdr.Input()[i]].Name();
    }
  os << '\n';
      
  // GOutput:
  first = true;
  os << "GOutput = ";
  for ( Index i=0; i<mdr.GOutput().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << wsv_group_names[mdr.GOutput()[i]];
    }
  os << '\n';

  // GInput:
  first = true;
  os << "GInput = ";
  for ( Index i=0; i<mdr.GInput().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << wsv_group_names[mdr.GInput()[i]];
    }
  os << '\n';

  // Keywords:
  first = true;
  os << "Keywords = ";
  for ( Index i=0; i<mdr.Keywords().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << mdr.Keywords()[i];
    }
  os << '\n';

  // Types:
  first = true;
  os << "Types = ";
  for ( Index i=0; i<mdr.Types().nelem(); ++i )
    {
      if (first) first=false;
      else os << ", ";

      os << TokValTypeName[mdr.Types()[i]];
    }
  os << "\n*--------------------------------------------------------------*\n";

  return os;
}

