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
    \end{verbatim}

 */
ARRAY<string> wsv_group_names;

/** The lookup information for the workspace variables. */
ARRAY<WsvRecord> wsv_data;

/** Initializes the workspace lookup data. */
void define_wsv_data()
{

  //--------------------< Build the group names array >--------------------
  // Initialize to empty, just in case.
  wsv_group_names.clear();

  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     Must be consistent with the enum in workspace.h! Later on,
     these enums could also be generated automatically, but that would
     have to take place before the wsv_data is defined, since that
     needs these enums.
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  wsv_group_names.push_back("VECTOR");
  wsv_group_names.push_back("MATRIX");

  // As a primitive consistency check, compare the size of
  // wsv_group_names with N_WSV_GROUPS:  
  assert(N_WSV_GROUPS==wsv_group_names.size());

  //--------------------< Build the wsv data >--------------------
  // Initialize to empty, just in case.
  wsv_data.clear();

  {
    static WsvPointer<VECTOR> p(&workspace.p_abs);
    wsv_data.push_back
      (WsvRecord
       ("p_abs",
	"The pressure grid for the absorption coefficients [hPa].",
	VECTOR_,
	&p));
  }
  
  {
    static WsvPointer<VECTOR> p(&workspace.f_abs);
    wsv_data.push_back
      (WsvRecord
       ("f_abs",
	"The frequency grid for the absorption coefficients [GHz].",
	VECTOR_,
	&p));
  }
  
  {
    static WsvPointer<VECTOR> p(&workspace.t_abs);
    wsv_data.push_back
      (WsvRecord
       ("t_abs",
	"Temperature associated with the pressures in p_abs [K]",
	VECTOR_,
	&p));
  }
  
  {
    static WsvPointer<MATRIX> p(&workspace.abs);
    wsv_data.push_back
      (WsvRecord
       ("abs",
	"The matrix of absorption coefficients.",
	MATRIX_,
	&p));
  }

  //  cout << "size = " << wsv_data.size() << '\n';
}


ostream& operator<<(ostream& os, const WsvRecord& wr)
{
  os << "Name = "        << wr.Name() << '\n'
     << "Description = " << wr.Description() << '\n'
     <<	"Group = "       << wsv_group_names[wr.Group()];

  return os;
}

