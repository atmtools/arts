/*-----------------------------------------------------------------------
FILE:      workspace.cc

INCLUDES:  This file contains only the function define_wsv_data, which
           sets the WSV group names and the lookup data for the WSVs.
	   You have to edit this function whenever you add a new
	   workspace variable. See workspace.h for more documentation.

GLOBALS:   None defined

FUNCTIONS: void define_wsv_data()

HISTORY:   10.06.2000 Created by Stefan Buehler
-----------------------------------------------------------------------*/

#include "arts.h"
#include "vecmat.h"
#include "workspace.h"


void define_wsv_data()
{
  /* The variables workspace, wsv_group_names, and wsv_data are defined
     in file workspace_aux.cc. */
  extern WorkSpace workspace;
  extern ARRAY<string> wsv_group_names;
  extern ARRAY<WsvRecord> wsv_data;

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
  wsv_group_names.push_back("ARRAYofMATRIX");
  wsv_group_names.push_back("ARRAYofVECTOR");
  wsv_group_names.push_back("Los");
  wsv_group_names.push_back("ARRAYofLineRecord");
  wsv_group_names.push_back("TagGroups");

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
    static WsvPointer<VECTOR> p(&workspace.refr_index);
    wsv_data.push_back
      (WsvRecord
       ("refr_index",
	"The refractive index associated with the pressures in p_abs [-].",
	VECTOR_,
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
	"The ground emission factor for the frequencies in f_abs [0-1].",
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

  {
    static WsvPointer<ARRAYofMATRIX> p(&workspace.source);
    wsv_data.push_back
      (WsvRecord
       ("source",
	"Mean source functions between the points of the LOS [W/(m3Hzsr)].",
	ARRAYofMATRIX_,
	&p));
  }

  {
    static WsvPointer<ARRAYofMATRIX> p(&workspace.trans);
    wsv_data.push_back
      (WsvRecord
       ("trans",
	"The transmissions between the points of the LOS [-].",
	ARRAYofMATRIX_,
	&p));
  }

  {
    static WsvPointer<VECTOR> p(&workspace.y_space);
    wsv_data.push_back
      (WsvRecord
       ("y_space",
	"Radiation entering the atmosphere at the start of the LOS.",
	VECTOR_,
	&p));
  }

  {
    static WsvPointer<VECTOR> p(&workspace.y);
    wsv_data.push_back
      (WsvRecord
       ("y",
	"The final spectrum, including sensor effects and data reduction.",
	VECTOR_,
	&p));
  }

  {
    static WsvPointer<ARRAYofMATRIX> p(&workspace.klos);
    wsv_data.push_back
      (WsvRecord
       ("klos",
	"Line of sight weighting functions.",
	ARRAYofMATRIX_,
	&p));
  }

  {
    static WsvPointer<ARRAYofLineRecord> p(&workspace.lines);
    wsv_data.push_back
      (WsvRecord
       ("lines",
	"A list of spectral line data.", 
	ARRAYofLineRecord_,
	&p));
  }

  {
    static WsvPointer<TagGroups> p(&workspace.tag_groups);
    wsv_data.push_back
      (WsvRecord
       ("tag_groups",
	"This is an array of arrays of OneTag tag definitions.\n"
	"It defines the available tag groups for the calculation\n"
	"of absorption coefficients and weighting functions.\n"
	"Contrary to the original Bredbeck definition, tags within a\n"
	"group must belong to the same species, because one VMR profile\n"
	"is associated with each tag group.", 
	TagGroups_,
	&p));
  }

  {
    static WsvPointer<int> p(&workspace.n_profiles);
    wsv_data.push_back
      (WsvRecord
       ("n_profiles",
	"The number of absorption profiles which should be calculated\n"
	"simultaneously. For 1-D calculations this is always 1.\n"
	"This should be useful for 2 and 3-D calculations. It affects\n"
	"the array dimension of pzT and t_and_all_vmrs, as well as the\n"
	"matrix dimension of raw_vmr_profiles", 
	int_,
	&p));
  }

  {
    static WsvPointer<ARRAYofMATRIX> p(&workspace.ptz);
    wsv_data.push_back
      (WsvRecord
       ("ptz",
	"Matrix has columns:\n"
	"1. Pressure in hPa\n"
	"2. Temperature in K\n"
	"3. Altitude in km\n"
	"\n"
	"The array dimension is determined by n_profiles (for 2 and 3-D)\n"
	"calculations.", 
	ARRAYofMATRIX_,
	&p));
  }

  {
    static WsvPointer<ARRAYofMATRIX> p(&workspace.raw_vmr_profiles);
    wsv_data.push_back
      (WsvRecord
       ("raw_vmr_profiles",
	"The individual VMR profiles. Each species VMR profile comes with a\n"
	"pressure profile. The different species can hence be on different grids.\n"
	"\n"
	"Matrix has columns:\n"
	"1. Pressure in hPa\n"
	"2. VMR profile (absolute number)\n"
	"3. more VMR profiles.\n"
	"\n"
	"The number of VMR profiles is given by n_profiles (1 in the 1-D case).\n"
	"The array dimension is determined by the number of tag groups.", 
	ARRAYofMATRIX_,
	&p));
  }

  {
    static WsvPointer<ARRAYofMATRIX> p(&workspace.t_and_all_vmrs);
    wsv_data.push_back
      (WsvRecord
       ("t_and_all_vmrs",
	"Matrix has columns:\n"
	"1. Pressure in hPa\n"
	"2. Temperature in K\n"
	"2. VMR profile (absolute number)\n"
	"3. more VMR profiles.\n"
	"\n"
	"The number of VMR profiles is determined by the number of tag groups.\n"
	"The array dimension is determined by n_profiles (for 2 and 3-D)\n"
	"calculations.", 
	ARRAYofMATRIX_,
	&p));
  }

  //  cout << "size = " << wsv_data.size() << '\n';
}

