#ifndef ppath_struct_h
#define ppath_struct_h

#include "matpackI.h"
#include "mystring.h"
#include "interpolation.h"

/*===========================================================================
  === The Ppath structure
  ===========================================================================*/

/** The structure to describe a propagation path and releated quantities.
 * 
 *  The fields of the structure are described more in detail inside the ARTS
 *  user guide (AUG).
 */
struct Ppath {
  /** Atmospheric dimensionality */
  Index dim;
  /** Number of points describing the ppath */
  Index np;
  /** The propagation path constant (only used for 1D) */
  Numeric constant;
  /** Radiative background */
  String background;
  /** Start position */
  Vector start_pos;
  /** Start line-of-sight */
  Vector start_los;
  /** Length between sensor and atmospheric boundary */
  Numeric start_lstep;
  /** The distance between start pos and the last position in pos */
  Matrix pos;
  /** Line-of-sight at each ppath point */
  Matrix los;
  /** Radius of each ppath point */
  Vector r;
  /** The length between ppath points */
  Vector lstep;
  /** End position */
  Vector end_pos;
  /** End line-of-sight */
  Vector end_los;
  /** The distance between end pos and the first position in pos */
  Numeric end_lstep;
  /** The real part of the refractive index at each path position */
  Vector nreal;
  /** The group index of refraction */  
  Vector ngroup;
  /** Index position with respect to the pressure grid */
  ArrayOfGridPos gp_p;
  /** Index position with respect to the latitude grid */
  ArrayOfGridPos gp_lat;
  /** Index position with respect to the longitude grid */
  ArrayOfGridPos gp_lon;
};

/** An array of propagation paths. */
typedef Array<Ppath> ArrayOfPpath;

#endif
