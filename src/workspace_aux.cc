/*-----------------------------------------------------------------------
FILE:      workspace_aux.cc

INCLUDES:  This file contains auxiliary material for the
           workspace, which used to be in workspace.cc. The reason for
	   the separation is that the stuff here hardly ever should be
	   changed, whereas workspace.cc has to be edited each time
	   a new variable is added. See workspace.h for more
	   documentation. 

FUNCTIONS: void define_wsv_map()
           ostream& operator<<(ostream& os, const WsvRecord& wr)

HISTORY:   10.06.2000 Created by Stefan Buehler
-----------------------------------------------------------------------*/

#include "arts.h"
#include "vecmat.h"
#include "workspace.h"


void define_wsv_map()
{
  extern const ARRAY<WsvRecord> wsv_data;
  extern std::map<string, size_t> WsvMap;

  for ( size_t i=0 ; i<wsv_data.size() ; ++i)
    {
      WsvMap[wsv_data[i].Name()] = i;
    }
}


ostream& operator<<(ostream& os, const WsvRecord& wr)
{
  extern const ARRAY<string> wsv_group_names;

  os << "Name = "        << wr.Name() << '\n'
     << "Description = " << wr.Description() << '\n'
     <<	"Group = "       << wsv_group_names[wr.Group()];

  return os;
}

