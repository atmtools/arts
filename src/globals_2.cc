/*-----------------------------------------------------------------------
FILE:      globals_2.cc

INCLUDES:  This file contains all global variable definitions that DO
	   depend on the automatically generated header file wsv.h. It
	   is necessary to have these in a separate file for compiler
	   technical reasons. (With g++-2.95.2 it does not work to
	   declare them constant in the same file where they are
	   defined.)

FUNCTIONS: None

HISTORY:   10.06.2000 Created by Stefan Buehler
-----------------------------------------------------------------------*/

#include "arts.h"
#include "vecmat.h"
//#include "workspace.h"
#include "methods.h"


//                     ---------------
//--------------------< Methods Stuff >--------------------
//                     ---------------

/** The lookup information for the workspace methods. */
ARRAY<MdRecord> md_data;

/** The map associated with md_data. */
std::map<string, size_t> MdMap;


