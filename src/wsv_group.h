/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>
                      Patrick Eriksson <patrick@rss.chalmers.se>

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

/*-----------------------------------------------------------------------
FILE:      wsv_group.h

INCLUDES:  Header stuff related to workspace variable groups. You have
           to edit this file if you are adding a new WSV group.

FUNCTIONS: None

CLASSES:   WsvGroup
           WsvP

HISTORY:   10.06.2000 Created by Stefan Buehler
-----------------------------------------------------------------------*/

#ifndef wsv_group_h
#define wsv_group_h


/** The enum type that identifies wsv groups.
    This is used to group workspace variables of the same type
    together, so that generic methods can operate on any of them. 

    \begin{verbatim}
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Has to be consistent with the group names in workspace.cc
    And with the WsvP classes!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \end{verbatim}

    @author Stefan Buehler */
enum WsvGroup{
  string_,
  int_,
  Numeric_,
  VECTOR_,
  MATRIX_,
  ARRAYofMATRIX_,
  ARRAYofVECTOR_,
  Los_,
  ARRAYofLineRecord_,
  ARRAYofARRAYofLineRecord_,
  TagGroups_
};

// For consistency check:
#define N_WSV_GROUPS 11




/** Base class for the different Wsv pointers. A virtual function for
    the conversion operator must be added here for each new group.  

    @author Stefan Buehler */ 
class WsvP {
public:
  virtual operator string*()        	    { safety(); return NULL; };
  virtual operator int*()           	    { safety(); return NULL; };
  virtual operator Numeric*()       	    { safety(); return NULL; };
  virtual operator VECTOR*()        	    { safety(); return NULL; };
  virtual operator MATRIX*()        	    { safety(); return NULL; };
  virtual operator ARRAYofMATRIX*() 	    { safety(); return NULL; };
  virtual operator ARRAYofVECTOR*() 	    { safety(); return NULL; };
  virtual operator Los*()           	    { safety(); return NULL; };
  virtual operator ARRAYofLineRecord*()        { safety(); return NULL; };
  virtual operator ARRAYofARRAYofLineRecord*() { safety(); return NULL; };
  virtual operator TagGroups*()        	    { safety(); return NULL; };
      
private:
  /** Safety check. This is called by all the virtual conversion
      operators. It just stops the program with an error message. This
      should never happen, because conversion should only be attempted
      to the correct type, for which an overloaded conversion operator
      exists. */
  void safety() {
    cerr << "Internal error: Tried to convert a WsvP "
	 << "pointer to the wrong type.\n";
    exit(1);
  };
};



#endif  // wsv_group_h
