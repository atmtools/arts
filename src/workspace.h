/*-----------------------------------------------------------------------
FILE:      workspace.h

INCLUDES:  The declaration of the workspace data type. You have to
           edit this file if you want to add a new workspace variable.

FUNCTIONS: None

CLASSES:   Workspace

HISTORY:   ??.??.???? Created by Stefan Buehler
           10.06.2000 Cleaned up by Stefan Buehler
-----------------------------------------------------------------------*/

#ifndef workspace_h
#define workspace_h

#include "arts.h"
#include "absorption.h"
#include "los.h"
#include "wsv_group.h"
#include "wsv_aux.h"


/** The declaration of the (great) workspace.
    @author Stefan Buehler */
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



#endif  // workshpace_h
