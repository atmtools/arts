/*-----------------------------------------------------------------------
FILE:      globals_1.cc

INCLUDES:  This file contains all global variable definitions that do
	   NOT depend on the automatically generated header file
	   wsv.h. It is necessary to have these in a separate file for
	   compiler technical reasons. (With g++-2.95.2 it does not
	   work to declare them constant in the same file where they
	   are defined.)

FUNCTIONS: None

HISTORY:   10.06.2000 Created by Stefan Buehler
-----------------------------------------------------------------------*/

#include "arts.h"
#include "vecmat.h"
#include "workspace.h"
#include "absorption.h"

//                     -----------------
//--------------------< Workspace Stuff >--------------------
//                     -----------------

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


//                     ------------------
//--------------------< Absorption Stuff >--------------------
//                     ------------------

/** The lookup information for all the different species. */
ARRAY<SpeciesRecord> species_data;

/** The map associated with species_data. */
std::map<string, size_t> SpeciesMap;

