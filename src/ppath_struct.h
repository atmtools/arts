#ifndef ppath_struct_h
#define ppath_struct_h

#include "matpack_data.h"
#include "mystring.h"
#include "interpolation.h"

enum PpathBackground {
  PPATH_BACKGROUND_UNDEFINED = 0,
  PPATH_BACKGROUND_SPACE = 1,
  PPATH_BACKGROUND_SURFACE = 2,
  PPATH_BACKGROUND_CLOUDBOX = 3,
  PPATH_BACKGROUND_TRANSMITTER = 4,
  PPATH_BACKGROUND_STOP_DISTANCE = 9  
};

/*===========================================================================
  === The Ppath structure
  ===========================================================================*/

/** The structure to describe a propagation path and related quantities. 
 *  
 *  The path is described in the viewing direction, that is starting
 *  from the observation side. This is the reversed direction compared
 *  to the photon direction. 
 * 
 *  Members start_pos and start_los give the position and LOS of the
 *  sensor (either real or imaginary). 
 *
 *  The path through the atmosphere is described by the members pos
 *  and los, that have np rows. If the sensor is totally outside of
 *  the atmosphere, np is 0. Otherwise the first row in pos and los
 *  describes the starting point of the path inside the atmosphere. If
 *  the sensor is inside the atmosphere, this equals start_pos/los.
 *  
 *  The members end_pos and end_los equal normally the last row of
 *  pos and los. The exception is if the path goes to a source or
 *  transmitter in space, when end_pos/los then give the position and
 *  LOS at the source/transmitter.
 */
struct Ppath {
  /** Atmospheric dimensionality */
  Index dim;                    // To be removed
  /** Number of points describing the ppath */
  Index np;
  /** The propagation path constant (only used for 1D) */
  Numeric constant;            // To be removed
  /** Radiative background */
  String background;
  enum PpathBackground backgroundZZZ;   // New version
  /** Start position */
  Vector start_pos;
  /** Start line-of-sight */
  Vector start_los;
  /** Length between sensor and atmospheric boundary */
  Numeric start_lstep;
  /** The positions representing the propagation path */
  Matrix pos;
  /** Line-of-sights at pos */
  Matrix los;
  /** Radius of each ppath point */
  Vector r;           // To be removed
  /** The length between pos (length np-1) */
  Vector lstep;
  /** End position */
  Vector end_pos;
  /** End line-of-sight */
  Vector end_los;
  /** The distance between end pos and the last position in pos */
  Numeric end_lstep;
  /** The real part of the refractive index at each path position */
  Vector nreal;
  /** The group index of refraction */  
  Vector ngroup;
  /** Index position with respect to the pressure grid */
  ArrayOfGridPos gp_p;          // To be removed
  /** Index position with respect to the latitude grid */
  ArrayOfGridPos gp_lat;        // To be removed
  /** Index position with respect to the longitude grid */
  ArrayOfGridPos gp_lon;        // To be removed

  friend std::ostream& operator<<(std::ostream& os, const Ppath& x);
};

/** An array of propagation paths. */
typedef Array<Ppath> ArrayOfPpath;

#endif
