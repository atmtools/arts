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
  \file   wsv_aux.h
  \brief  Auxiliary header stuff related to workspace variable
          groups. Normally you should not need to edit this file. 


  \author Stefan Buehler
  \date   2000-06-10
*/

#ifndef wsv_aux_h
#define wsv_aux_h

#include "wsv_groups.h"

/*! Template for Wsv Pointers. This defines for each pointer class the
    conversion operator back to the type that it is pointing
    to. 

    This makes it possible to store arbitrary pointers in an array of
    pointers to WsvP.

    /author Stefan Buehler */
template<class T>
class WsvPointer : public WsvP {
public:
  WsvPointer(T* x) : mx(x) { /* Nothing to do here. */ };
  operator T*() { return mx; }
private:
  T* mx;
};



/*! This class contains all static information for one workspace
    variable.

    The program make_wsv_h.cc uses these records to generate the file
    wsv.h, which contains both the declaration of the wsv handles and
    the declaration of the workspace itself.

    \author Stefan Buehler */
class WsvRecord {
public:
  /*! Initializing constructor.

    This is used by define_wsv_data() to set the information for each
    workspace variable. */
  WsvRecord(const char name[],
	    const char description[],
	    const size_t group)
    : mname(name),
      mdescription(description),
      mgroup(group)
  {
    //    Nothing to do here.
  }
  /*! Name of this workspace variable. */
  const string&  Name()        const { return mname;        }   
  /*! A text describing this workspace variable. */
  const string&  Description() const { return mdescription; }
  /*! The wsv group to which this variable belongs. */
  const size_t   Group()       const { return mgroup;       }
private:
  string mname;
  string mdescription;
  size_t mgroup;
};

/*! Output operator for WsvRecord.
  \author Stefan Buehler */
ostream& operator<<(ostream& os, const WsvRecord& wr);


/*! Define the lookup data for the workspace variables. The array
    wsv_data contains all that we need to know about each workspace
    variable. The array WsvGroupName contains the names of the work
    space variable groups. These two lookup tables are global
    variables. They can be made visible anywhere with an extern
    declaration.

    \author Stefan Buehler */
void define_wsv_data();

/*! Define WsvMap. WsvMap can be used to find workspace variable data
    by name.

    \author Stefan Buehler */ 
void define_wsv_map();


#endif   // wsv_aux_h
