/*-----------------------------------------------------------------------
FILE:      workspace_aux.cc

INCLUDES:  This file contains auxiliary material for the
           workspace, which used to be in workspace.cc. The reason for
	   the separation is that the stuff here hardly ever should be
	   changed, whereas workspace.cc has to be edited each time
	   a new variable is added. See workspace.h for more
	   documentation. 

GLOBALS:   workspace
           wsv_group_names
	   wsv_data
	   WsvMap

FUNCTIONS: void define_wsv_map()
           ostream& operator<<(ostream& os, const WsvRecord& wr)

HISTORY:   10.06.2000 Created by Stefan Buehler
-----------------------------------------------------------------------*/

#include "arts.h"
#include "vecmat.h"
#include "workspace.h"


/** The workspace itself. */
WorkSpace workspace;

/** The names associated with Wsv groups as strings. 
    
    \begin{verbatim}
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Must be consistent with the enum in workspace.h! Later on,
    these enums could also be generated automatically, but that would
    have to take place before the wsv_data is defined, since that
    needs these enums.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \end{verbatim}  */
ARRAY<string> wsv_group_names;

/** The lookup information for the workspace variables. */
ARRAY<WsvRecord> wsv_data;

/** The map assiciated with wsv_data. */
std::map<string, size_t> WsvMap;


void define_wsv_map()
{
  // These two declarations are unnecessary because the variables are
  // defined in this file. They lead to errors with g++-2.95.2.
  //  extern const ARRAY<WsvRecord> wsv_data;
  //  extern std::map<string, size_t> WsvMap;

  for ( size_t i=0 ; i<wsv_data.size() ; ++i)
    {
      WsvMap[wsv_data[i].Name()] = i;
    }
}


ostream& operator<<(ostream& os, const WsvRecord& wr)
{
  os << "Name = "        << wr.Name() << '\n'
     << "Description = " << wr.Description() << '\n'
     <<	"Group = "       << wsv_group_names[wr.Group()];

  return os;
}

