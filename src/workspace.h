#ifndef workspace_h
#define workspace_h

#include <iostream>
#include "vecmat.h"

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
  VECTOR_,
  MATRIX_,
  Numeric_
};

// For consistency check:
#define N_WSV_GROUPS 3


// Yes, this is the declaration of the great workspace itself.
class WorkSpace {
public:
  VECTOR p_abs;
  VECTOR f_abs;
  VECTOR t_abs;
  MATRIX abs;
  Numeric dummy1;
};




/** Base class for the different Wsv pointers. A virtual function for
    the conversion operator must be added here. You can put variables
    of a new type in the group Other_ if you don't want to define a
    group for them. However, then you can not use generic methods on
    these variables! A conversion operator to other is on purpose not
    defined. */ 
class WsvP {
public:
  virtual operator string*() { safety(); return NULL; };
  virtual operator VECTOR*() { safety(); return NULL; };
  virtual operator MATRIX*() { safety(); return NULL; };
  virtual operator Numeric*() { safety(); return NULL; };
      
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

#endif  // workshpace_h
