/*! \file  wsv_group.h
    \brief Defines the enum type that acts as a
    handle for workspace variables groups.

    Also defined here is a special pointer class that can hold
    a pointer to any workspace variable.

    This file was generated automatically by make_wsv_groups_h.cc.
    <b>DO NOT EDIT!</b>

    \date Aug  5 2000, 14:30:14 */

#ifndef wsv_group_h
#define wsv_group_h

/*! This is only used for a consistency check. You can get the
    number of groups from wsv_group_names.size(). */
#define N_WSV_GROUPS 11

/*! The enum type that identifies wsv groups.
    This is used to group workspace variables of the same type
    together, so that generic methods can operate on any of them. */
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
  TagGroups_,
};

/*! Base class for the different Wsv pointers.
    This contains a virtual function for the
    conversion operator for each group.

    \author Stefan Buehler */
class WsvP {
public:
  virtual operator string*(){safety();return NULL;};
  virtual operator int*(){safety();return NULL;};
  virtual operator Numeric*(){safety();return NULL;};
  virtual operator VECTOR*(){safety();return NULL;};
  virtual operator MATRIX*(){safety();return NULL;};
  virtual operator ARRAYofMATRIX*(){safety();return NULL;};
  virtual operator ARRAYofVECTOR*(){safety();return NULL;};
  virtual operator Los*(){safety();return NULL;};
  virtual operator ARRAYofLineRecord*(){safety();return NULL;};
  virtual operator ARRAYofARRAYofLineRecord*(){safety();return NULL;};
  virtual operator TagGroups*(){safety();return NULL;};

private:
/*! Safety check. This is called by all the virtual conversion
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
