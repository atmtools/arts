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
  wsv_group_names.push_back("string");
  wsv_group_names.push_back("int");
  wsv_group_names.push_back("Numeric");
  wsv_group_names.push_back("VECTOR");
  wsv_group_names.push_back("MATRIX");
  wsv_group_names.push_back("Los");

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
    static WsvPointer<VECTOR> p(&workspace.t_abs);
    wsv_data.push_back
      (WsvRecord
       ("t_abs",
	"Temperature associated with the pressures in p_abs [K]",
	VECTOR_,
	&p));
  }

  {
    static WsvPointer<VECTOR> p(&workspace.z_abs);
    wsv_data.push_back
      (WsvRecord
       ("z_abs",
	"Vertical altitudes associated with the pressures in p_abs [m]",
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
    static WsvPointer<MATRIX> p(&workspace.abs);
    wsv_data.push_back
      (WsvRecord
       ("abs",
	"The matrix of absorption coefficients.",
	MATRIX_,
	&p));
  }

  {
    static WsvPointer<VECTOR> p(&workspace.view1);
    wsv_data.push_back
      (WsvRecord
       ("view1",
	"Viewing angle 1, the angle between zenith and the LOS [deg].",
	VECTOR_,
	&p));
  }

  {
    static WsvPointer<Numeric> p(&workspace.z_plat);
    wsv_data.push_back
      (WsvRecord
       ("z_plat",
	"The vertical altitude, above the geiod, of the platform [m].",
	Numeric_,
	&p));
  }

  {
    static WsvPointer<Numeric> p(&workspace.l_step);
    wsv_data.push_back
      (WsvRecord
       ("l_step",
	"The length (along the LOS) between the points of LOS [m].",
	Numeric_,
	&p));
  }

  {
    static WsvPointer<int> p(&workspace.refr);
    wsv_data.push_back
      (WsvRecord
       ("refr",
	"Boolean to consider refraction (0=no refraction).",
	int_,
	&p));
  }

  {
    static WsvPointer<Numeric> p(&workspace.l_step_refr);
    wsv_data.push_back
      (WsvRecord
       ("l_step_refr",
	"The step length (along LOS) when determining the LOS with refraction [m].",
	Numeric_,
	&p));
  }

  {
    static WsvPointer<int> p(&workspace.cbgr);
    wsv_data.push_back
      (WsvRecord
       ("cbgr",
	"Boolean to consider cosmic background radiation (0=no cbgr).",
	int_,
	&p));
  }

  {
    static WsvPointer<Numeric> p(&workspace.z_ground);
    wsv_data.push_back
      (WsvRecord
       ("z_ground",
	"The vertical altitude above the geiod of the ground [m].",
	Numeric_,
	&p));
  }

  {
    static WsvPointer<Numeric> p(&workspace.t_ground);
    wsv_data.push_back
      (WsvRecord
       ("t_ground",
	"The physical temperature of the ground [K].",
	Numeric_,
	&p));
  }

  {
    static WsvPointer<VECTOR> p(&workspace.e_ground);
    wsv_data.push_back
      (WsvRecord
       ("e_ground",
	"The emission factor for the ground at the frequencies in f_abs [0-1].",
	VECTOR_,
	&p));
  }

  {
    static WsvPointer<Los> p(&workspace.los);
    wsv_data.push_back
      (WsvRecord
       ("los",
	"Structure to define the line of sight (LOS) for 1d cases.", 
	Los_,
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

