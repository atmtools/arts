#ifndef workspace_h
#define workspace_h

#include <iostream>
#include "vecmat.h"
#include "absorption.h"


/** Define the enum type that identifies wsv groups.
    This is used to group workspace variables of the same type
    together, so that generic methods can operate on any of them. 

    \begin{verbatim}
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Has to be consistent with the group names in workspace.cc
    And with the WsvP classes!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \end{verbatim} */
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
  TagGroups_
};

// For consistency check:
#define N_WSV_GROUPS 10



// ======================================================================
// === Definition of workspace variables of not basic types
// ======================================================================


/** The line of sight (LOS). The LOS structure has the fields:
       ARRAYofVECTOR  p;
       VECTOR         l_step;
       ARRAY<int>     ground;
       ARRAY<int>     start;
       ARRAY<int>     stop;
    where 
       p        The pressures along LOS
       l_step   The geometrical length along LOS between the points.
       start    start index for the iteration
       stop     stop index for the iteration
       ground   O if no intersection with the ground. Else, GROUND
                gives the index for the ground.  

    The LOS is defined in equal long geometrical steps along the path.
    This step length (L_STEP) can vary between the viewing angles.

    Spectra are calculated in the following way (by RTE_ITERATE in m_los):
       1. Iteration from START down to 1 or GROUND
       2. If GROUND, including the effect of the ground reflection.
       3. Iteration from 1 or GROUND-1 to STOP

    The START and STOP variables make it possible to use a possible symmetry
    for 1D calculations. For example, for limb sounding from space, START
    and STOP are both set to the length of P. The GROUND variable is for
    1D calculations either 0 or 1.

    For cases without symmetry (upward looking and 2D), STOP is always 1
    and corresponds to the point closest to the sensor. Accordingly, START
    corresponds to the point of LOS furthest away from the sensor.

    The GROUND variable is used both as a flag to indicate ground 
    intersections of the LOS, and a variable to give the position of the
    ground. As mentioned, for 1D cases, the ground is always placed at 
    index 1. For 2D cases, GROUND gives the index for the ground point, 
    that is, the point of LOS with index GROUND corresponds to the ground 
    level.

    Written by Patrick Eriksson 07.06.00
                                                                    */
 struct Los {
  ARRAYofVECTOR  p;
  VECTOR         l_step;
  ARRAY<int>     ground;
  ARRAY<int>     start;
  ARRAY<int>     stop;
};



// ======================================================================
// === Definition of the (great) workspace
// ======================================================================

class WorkSpace {
public:
  VECTOR              p_abs;
  VECTOR              t_abs;
  VECTOR              z_abs;
  VECTOR              f_abs;
  MATRIX              abs;
  VECTOR              view1;
  Numeric             z_plat;
  Numeric             l_step;
  int                 refr;
  Numeric             l_step_refr;
  VECTOR              refr_index;
  Numeric             z_ground;
  Numeric             t_ground;
  VECTOR              e_ground;
  Los                 los; 
  ARRAYofMATRIX       source;
  ARRAYofMATRIX       trans;
  VECTOR              y_space;
  VECTOR              y;
  ARRAYofMATRIX       klos;
  ARRAYofLineRecord   lines;
  TagGroups           tag_groups;
  int                 n_profiles;
  ARRAYofMATRIX       ptz;
  ARRAYofMATRIX       raw_vmr_profiles;
  ARRAYofMATRIX	      t_and_all_vmrs;
};




/** Base class for the different Wsv pointers. A virtual function for
    the conversion operator must be added here. You can put variables
    of a new type in the group Other_ if you don't want to define a
    group for them. However, then you can not use generic methods on
    these variables! A conversion operator to other is on purpose not
    defined. */ 
class WsvP {
public:
  virtual operator string*()        	{ safety(); return NULL; };
  virtual operator int*()           	{ safety(); return NULL; };
  virtual operator Numeric*()       	{ safety(); return NULL; };
  virtual operator VECTOR*()        	{ safety(); return NULL; };
  virtual operator MATRIX*()        	{ safety(); return NULL; };
  virtual operator ARRAYofMATRIX*() 	{ safety(); return NULL; };
  virtual operator ARRAYofVECTOR*() 	{ safety(); return NULL; };
  virtual operator Los*()           	{ safety(); return NULL; };
  virtual operator ARRAYofLineRecord*() { safety(); return NULL; };
  virtual operator TagGroups*()        	{ safety(); return NULL; };
      
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

/** Template for Wsv Pointers. All you have to do when you want to add
    a new Wsv group `Smurf' is add a virtual method

    \begin{verbatim}
    virtual operator MATRIX*() = 0;
    \end{verbatim}

    to the common base class WsvP. */
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

#endif  // workshpace_h
