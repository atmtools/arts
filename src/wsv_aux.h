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



/** This class contains all static information for one workspace
    variable. */
class WsvRecord {
public:
  WsvRecord(const char name[],
	    const char description[],
	    const size_t group,
	    WsvP* const pointer)
    : mname(name),
      mdescription(description),
      mgroup(group),
      mpointer(pointer)
  {
    // Assign mtotal to mid and increase mtotal by 1:
    //    mid = mtotal++;
  }
  //  const int      Total()       const { return mtotal;       }
  //  const int      Id()          const { return mid;          }
  const string&  Name()        const { return mname;        }   
  const string&  Description() const { return mdescription; }
  const size_t   Group()       const { return mgroup;       }
  WsvP* const    Pointer()     const { return mpointer;     }
private:
  // Total number of WsvRecords around:
  //  static int mtotal;
  // Id of this record:
  // (The Id is not really used, only for debugging. Eventually, the
  // Handle of this record should have the same value as Id.)
  //  int mid;
  /** Name of this workspace variable. */
  string mname;
  /** A text describing this workspace variable. */
  string mdescription;
  /** The wsv group to which this variable belongs. */
  size_t mgroup;
  /** Pointer to smart pointer to the variable itself. */
  WsvP* mpointer;
};

/** Output operator for WsvRecord.
    @author Stefan Buehler */
ostream& operator<<(ostream& os, const WsvRecord& wr);


/** Define the lookup data for the workspace variables. The array
    wsv_data contains all that we need to know about each workspace
    variable. The array WsvGroupName contains the names of the work
    space variable groups. These two lookup tables are global
    variables. They can be made visible anywhere with an extern
    declaration. */
void define_wsv_data();

/** Define WsvMap. WsvMap can be used to find workspace variable data
    by name. */ 
void define_wsv_map();


#endif   // wsv_aux_h
