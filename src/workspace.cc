/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>
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

/*!
  \file   workspace.cc
  \brief  Definition of function wsv_data.

  This file contains the function define_wsv_data, which
  sets the WSV group names and the lookup data for the WSVs.
  You have to edit this function whenever you add a new
  workspace variable. 

  \author Stefan Buehler
  \date   2000-06-10

  Removed the wsv pointers. They are now in a separate place.
  \author Stefan Buehler
  \date   2000-08-10 */


#include "arts.h"
#include "matpackI.h"
#include "matpackIII.h"
#include "array.h"
#include "auto_wsv_groups.h"
#include "wsv_aux.h"
#include "ppath.h"

/*! The lookup information for the workspace variables. */
Array<WsvRecord> wsv_data;

void define_wsv_data()
{

  //--------------------< Build the wsv data >--------------------
  // Initialize to empty, just in case.
  wsv_data.resize(0);

  /* Template for description strings:

     ----------------------------------------------------------------------
     Brief description of the variable (1 line).
     
     Detailed description of the variable. Don't be too short here,
     this is the main place where your documentation should be. I
     really recommend to edit this in a text buffer, so that you can
     do some re-formatting until it looks nice. Only at the end put it
     in quotes and add the line breaks.

     Use blank lines to separate paragraphs.  There really should be a
     detailed descriptions of all component of your variable, if it
     has a complicated type. Also some detailed discussion of the
     dimensions if necessary. Also some detailed discussion of the
     members if your variable is a structure.

     Usage:      Set by user (or "Function output.")

     Units:      E.g., kg/m

     Dimensions: [ first dimension, second dimension, ... ]

     Members:    Here you would list the members if your
                 variable is a structure.
     ----------------------------------------------------------------------

     The dashed lines are not part of the template, they just show you
     where it starts and ends. 

     Give the keywords above only if they apply, i.e., Members only
     for a structure, Units only for a physical variable. 
  */


  /////////////////////////////////////////////////////////////////////////////
  // Let's put in the variables in alphabetical order. This gives a clear rule
  // for where to place a new variable and this gives a nicer results when
  // the functions are listed by "arts -w all".
  // No distinction is made between uppercase and lowercase letters. The sign
  // "_" comes after all letters.
  // Patrick Eriksson 2002-05-08
  /////////////////////////////////////////////////////////////////////////////


  wsv_data.push_back
   (WsvRecord
    ("atmosphere_dim",
     "The atmospheric dimensionality (1-3).\n"
     "\n"
     "This variable defines the complexity of the atmospheric structure.\n"
     "The dimensionality is given by an integer between 1 and 3, where 1\n"
     "means 1D etc. This is the master variable for the atmospheric\n"
     "dimensionality, variables which size changes with the dimensionality\n"
     "are checked to match this variable. \n"
     "\n"
     "Functions adapt automatically to this variable. That is, it should\n"
     "not be needed to change any functions if the dimensionality is\n"
     "changed.\n"
     "\n"
     "The atmospheric dimensionalities (1D, 2D and 3D) are defined in the\n"
     "user guide (look for \"atmospheric dimensionality\" in the index).\n" 
     "\n"
     "Usage:      Set by the user.",
     Index_ ));

  wsv_data.push_back
   (WsvRecord
    ("azimuth_angle_1d",
     "The azimuthal angle for a 1D observation.\n"
     "\n"
     "This variable is onlyu used when calculating a curvature radius of the\n"
     "geoid for 1D cases. The given value shall be the angle between the\n"
     "observation and meridian plane where 0 degrees means an observation in\n"
     "the N-S direction.\n"
     "\n"
     "Values shall be in the range [-180,180].\n"
     "\n"
     "Usage:      Set by the user.\n"
     "\n"
     "Units:      degrees",
     Numeric_ ));

  wsv_data.push_back
   (WsvRecord
    ("a_los",
     "A line-of-sight (for test purposes).\n"
     "\n"
     "The purpose of this variable and *a_pos* are to enable calling of\n"
     "*ppathCalc* and maybe other functions from the workspace. This can be\n"
     "of interest both for testing of the functions and to display\n"
     "intermediate results as the propagation path. \n"
     "\n"
     "For 1D and 2D cases, *a_los* is a vector of length 1 holding the \n"
     "zenith angle. For 3D, the length of the vector is 2, where the\n"
     "additional element is the azimuthal angle. These angles are defined\n"
     "in the ARTS user guide (AUG). Look in the index for \"zenith angle\"\n"
     "and \"azimuthal angle\".\n"
     "\n"
     "Usage:      To call *ppathCalc* as an individual function. Other\n"
     "            purposes can be possible.\n"
     "\n"
     "Units:      [ degree, degree ]\n"
     "\n"
     "Dimensions: [ 1 for 1D and 2D, 2 for 3D ]",
     Vector_ ));

  wsv_data.push_back
   (WsvRecord
    ("a_pos",
     "A geographical position (for test purposes).\n"
     "\n"
     "The purpose of this variable and *a_los* are to enable calling of\n"
     "*ppathCalc* and maybe other functions from the workspace. This can be\n"
     "of interest both for testing of the functions and to display\n"
     "intermediate results as the propagation path. \n"
     "\n"
     "This variable is a vector with a length equalling the atmospheric\n"
     "dimensionality. The first element is the radius (from the coordinate\n"
     "system centre) of the position. Element 2 is the latitude and element\n"
     "3 is the longitude. Please note that the vertical position is given\n"
     "as the radius, not the altitude above the geoid.\n"
     "\n"
     "Usage:      To call *ppathCalc* as an individual function. Other\n"
     "            purposes can be possible.\n"
     "\n"
     "Units:      [ m, degree, degree ]\n"
     "\n"
     "Dimensions: [ atmosphere_dim ]",
     Vector_ ));

  wsv_data.push_back
   (WsvRecord
    ("blackbody_ground",
     "Flag to force the ground to be treated as a blackbody.\n"
     "\n"
     "This variable must have the value 0 or 1, where 1 means that the\n"
     "ground is treated as a blackbody (an emissivity of 1). This choice\n"
     "affects the starting point of propagation paths (including an inter-\n"
     "section with the ground) and the allowed choice for the scattering\n"
     "box.\n"
     "\n"
     "When *blackbody_ground* is set to 1, all values of *e_ground* must\n"
     "equal 1.\n"
     "\n"
     "See further the ARTS user guide (AUG). Use the index to find where\n"
     "this variable is discussed. The variable is listed as a subentry to\n"
     "\"workspace variables\".\n"
     "\n"
     "Usage:      Set by user.\n"
     "\n"
     "Units:      boolean",
     Index_ ));

  wsv_data.push_back
   (WsvRecord
    ("cloudbox_on",
     "Flag to activate the cloud box.\n"
     "\n"
     "Scattering calculations are confined to a part of the atmosphere\n"
     "denoted as the cloud box. The extension of the cloud box is given by\n"
     "*cloudbox_limits*. This variable tells functions if a cloud box is\n"
     "activated or not. \n"
     "\n"
     "See further the ARTS user guide (AUG). Use the index to find where\n"
     "this variable is discussed. The variable is listed as a subentry to\n"
     "\"workspace variables\".\n"
     "\n"
     "Usage:      Set by the user.\n"
     "\n"
     "Units:      boolean",
     Index_ ));

  wsv_data.push_back
   (WsvRecord
    ("cloudbox_limits",
     "The limits of the cloud box.\n"
     "\n"
     "This variable defines the extension of the cloud box. The cloud box is\n"
     "defined to be rectangular in the used coordinate system, with limits\n"
     "exactly at points of the involved grids. This means, for example, that\n"
     "the vertical limits of the cloud box are two pressure surfaces. For\n"
     "2D, the angular extension of the cloud box is between two points of\n"
     "the latitude grid (Figure~\ref{fig:fm_defs:2d}), and likewise for 3D\n"
     "but then also with a longitude extension between two grid points.  The\n"
     "latitude and longitude limits for the cloud box cannot be placed at\n"
     "the end points of the corresponding grid as it must be possible to\n"
     "calculate the incoming intensity field.\n"
     "\n"
     "In general the lower vertical limit must be the lowest pressure\n"
     "surface (that is, the first element of *cloudbox_limits is 0). This\n"
     "means that the practical lower limit of the cloud box is the ground.\n"
     "However, if the ground is treated to a blackbody (*blackbody_ground* =\n"
     "1) the lower limit can be above the ground. \n"
     " \n"
     "The variable *cloudbox_limits* is an array of index value with length\n"
     "twice *atmosphere_dim*. For each dimension there is a lower limit and\n"
     "an upper limit. The order of the dimensions is as usual pressure,\n"
     "latitude and longitude. The upper limit index must be greater then the\n"
     "lower limit index. For example, *cloudbox_limits* = [0 5 4 11] means\n"
     "that cloud box extends between pressure levels 0 and 5, and latitude\n"
     "points 4 and 11.\n"
     "\n"
     "If *cloudbox_on* = 0, the content of this variable is neglected, but\n"
     "it must be initiated to some dummy values.\n"
     "\n"
     "See further the ARTS user guide (AUG). Use the index to find where\n"
     "this variable is discussed. The variable is listed as a subentry to\n"
     "\"workspace variables\".\n"
     "\n"
     "Usage:      Set by the user, either directly or using a function\n"
     "            checking the extension of scattering particles.\n"
     "\n"
     "Units:      Index values.\n"
     "\n"
     "Dimensions: [ 2 * atmosphere_dim ]",
     ArrayOfIndex_ ));

  wsv_data.push_back
   (WsvRecord
    ("e_ground",
     "Ground emissivity.\n"
     "\n"
     "The emission of the ground is calculated as the product between\n"
     "*e_ground* and the Planck function for the temperature given by\n"
     "*t_ground*. The ground emissivity and temperature are interpolated to\n"
     "the point of interest before the product is calculated.\n"
     "\n"
     "The ground emissivity is a value in the range [0,1].\n"
     "\n"
     "The emissivity for a point between the grid crossings is obtained by\n"
     "linear (1D) or bi-linear (2D) interpolation of the *e_ground*.\n"
     "No interpolation in frequency is performed.\n"
     "\n"
     "If the ground is assumed to be a blackbody (*blackbody_ground* = 1),\n"
     "all values of *e_ground* must equal 1. For such cases, *e_ground*\n"
     "can be set to have the size 1x1x1.\n"
     "\n"
     "See further the ARTS user guide (AUG). Use the index to find where\n"
     "this variable is discussed. The variable is listed as a subentry to\n"
     "\"workspace variables\".\n"
     "\n"
     "Usage:      Set by user.\n"
     "\n"
     "Units:      - (0-1)\n"
     "\n"
     "Dimensions: [ f_grid, lat_grid, lon_grid ]",
     Tensor3_ ));

  wsv_data.push_back
   (WsvRecord
    ("lat_grid",
     "The latitude grid.\n"
     "\n"
     "The latitudes for which the atmospheric fields are defined. The\n"
     "atmosphere is undefined outside the range covered by the grid.\n"
     "The grid must be sorted in increasing order, with no repetitions.\n"
     "\n"
     "Geocentric latitudes shall be used.\n"
     "\n"
     "For 1D calculations this vector shall be set to be empty, but the\n"
     "number of latitudes shall be considered to be 1 when examining the\n"
     "size of variables.\n"
     "\n" 
     "In the case of 2D, the latitudes shall be interpreted as the angular\n"
     "distance inside the orbit plane from an arbitrary zero point. Any\n"
     "latitude values are accepted for 2D.\n"
     "\n"
     "For 3D, the valid latitude range is [-90,90].\n"
     "\n"
     "See further the ARTS user guide (AUG). Use the index to find where\n"
     "this variable is discussed. The variable is listed as a subentry to\n"
     "\"workspace variables\".\n"
     "\n"
     "Usage:      Set by the user.\n"
     "\n"
     "Units:      degrees",
     Vector_ ));

  wsv_data.push_back
   (WsvRecord
    ("latitude_1d",
     "The latitude for a 1D observation.\n"
     "\n"
     "This variable is used to assign a latitude valid for a 1D case. A\n"
     "special variable is needed for 1D as there exists no latitude grid for\n"
     "such cases. The variable can be used, for example, to set the geoid\n"
     "radius or select atmospheric profiles from a 2D/3D climatology. \n"
     "\n"
     "For limb sounding, the choosen latitude should normally be selected\n"
     "to be valid for the tangent points, rather than for the sensor.\n"
     "\n"
     "Values shall be inside the range [-90,90].\n"
     "\n"
     "Usage:      Set by the user.\n"
     "\n"
     "Units:      degrees",
     Numeric_ ));

  wsv_data.push_back
   (WsvRecord
    ("lon_grid",
     "The longitude grid.\n"
     "\n"
     "The longitudes for which the atmospheric fields are defined. The\n"
     "atmosphere is undefined outside the range covered by the grid.\n"
     "The grid must be sorted in increasing order, with no repetitions.\n"
     "\n"
     "For 1D and 2D calculations this vector shall be set to be empty, but\n"
     "the number of longitudes shall be considered to be 1 when examining\n"
     "the size of variables.\n"
     "\n"
     "Allowed values for longidudes is the range [-360,360].\n"
     "\n"
     "See further the ARTS user guide (AUG). Use the index to find where\n"
     "this variable is discussed. The variable is listed as a subentry to\n"
     "\"workspace variables\".\n"
     "\n"
     "Usage:      Set by the user.\n"
     "\n"
     "Units: degrees",
     Vector_ ));

  wsv_data.push_back
   (WsvRecord
    ("main_agenda",
     "The agenda corresponding to the entire controlfile. This is executed\n"
     "when ARTS is run.",
     Agenda_));

  wsv_data.push_back
   (WsvRecord
    ("p_grid",
     "The pressure grid.\n"
     "\n"
     "The pressure surfaces on which the atmospheric fields are defined.\n"
     "This variable must always be defined. The grid must be sorted in\n"
     "decreasing order, with no repetitions.\n"
     "\n"
     "No gap between the lowermost pressure level and the ground is \n"
     "allowed. The uppermost pressure level defines the practical upper\n"
     "limit of the atmosphere as vacuum is assumed above.\n"
     "\n"
     "See further the ARTS user guide (AUG). Use the index to find where\n"
     "this variable is discussed. The variable is listed as a subentry to\n"
     "\"workspace variables\".\n"
     "\n"
     "Usage:      Set by the user.\n"
     "\n"
     "Units:      Pa",
     Vector_ ));

  wsv_data.push_back
   (WsvRecord
    ("ppath",
     "The propagation path for one line-of-sight.\n"
     "\n"
     "This variable described the total (pencil beam) propagation path for\n"
     "a given combination of starting point and line-of-sight. The path is\n"
     "described by a data structure of type Ppath. This structure contains\n"
     "also additional fields to faciliate the calculation of spectra and\n"
     "interpolation of the atmospheric fields.\n"
     "\n"
     "The data struture is to extensive to be described here, but it is\n"
     "described carefully in the ARTS user guide (AUG). Use the index to\n"
     "find where the data structure, Ppath, for propagation paths is \n"
     "discussed. It is listed as a subentry to \"data structures\".\n"
     "\n"
     "Usage:      Output from functions of the type *ppathCalc*.\n"
     "\n"
     "Members:    To be written.",
     Ppath_ ));

  wsv_data.push_back
   (WsvRecord
    ("ppath_partial",
     "A part of a propagation path.\n"
     "\n"
     "The main intention of this variable is communication with the agenda\n"
     "*ppath_step*.\n"
     "\n"
     "Type \"arts -d ppath_step\" and \"arts -d ppathCalc\" for more\n"
     "information on this variable and the calculation of propagation\n"
     "paths. Or read the chapter on propgation paths in the ARTS user\n"
     "guide.\n"
     "\n"
     "Usage:      In/output to/from the agenda *ppath_step*.\n"
     "\n"
     "Members:    To be written.",
     Ppath_ ));

  wsv_data.push_back
   (WsvRecord
    ("ppath_step",
     "Calculation of a propagation path step.\n"
     "\n"
     "A propagation path step is defined as the path between some point to\n"
     "a crossing with either the pressure, latitude or longitude grid,\n"
     "and this agenda performs the calculations to determine such a\n"
     "partial propagation path. The starting point is normally a grid\n"
     "crossing point, but can also be an arbitrary point inside the\n"
     "atmosphere, such as the sensor position. Only points inside the model\n"
     "atmosphere are handled.\n"
     "\n"
     "The communication between this agenda and the calling function is\n"
     "handled by *ppath_partial*. That variable is used both as input and\n"
     "output to the *ppath_step* agenda. The agenda gets back\n"
     "*ppath_partial* as returned to the calling function and the last path\n"
     "point hold by the structure is accordingly the starting point for the\n"
     "new calculations. If a total propagation path shall be determined,\n"
     "this agenda is called repeatedly until the starting point of the\n"
     "propagation path is found and *ppath_partial* will hold all path\n"
     "steps that together make up *ppath*. The starting point is not\n"
     "included in the returned structure.\n"
     "\n"
     "The path is determined by starting at the end point and moving\n"
     "backwards to the starting point (the forward direction is defined as\n"
     "the way the photons are moving). The calculations are initiated by\n"
     "filling *ppath_partial* with the practical end point of the path.\n"
     "This is either the position of the sensor (true or hypothetical), or\n"
     "some point at the top of the atmosphere (determined by geometrical\n"
     "calculations starting at the sensor). This initialisation is not\n"
     "handled by *ppath_step*. All fields of *ppath_partial* are set by\n"
     "*ppath_step*. All the fields of *ppath_partial*, beside *l_step*, act\n"
     "as for *ppath*. In the *ppath_partial* initiated by the calling\n"
     "function, the *l_step* vector has length 0. For the later structures\n"
     "the vector has the length given by the field *np* (one value more than\n"
     "in *ppath*), where the additional value is the distance between the\n"
     "starting point and the first point in the structure.\n"
     "\n"
     "The zenith and azimuthal angles in *ppath_partial* are valid for\n"
     "the direction the path is calculated. The direction is later reversed\n"
     "when values are moved to *ppath*.\n"
     "\n"
     "The agenda performs only calculations to next crossing of a grid, all\n"
     "other tasks must be performed by the calling function, with one\n"
     "exception. If there is an intersection of a blackbody ground, the\n"
     "calculations stop at this point. This is flagged by *ppath_step* by\n"
     "setting the background field of *ppath_partial*. Beside this,\n"
     "the calling function must check if the starting point of the\n"
     "calculations is inside the scattering box or below the ground level,\n"
     "and check if the starting point of the path has been reached. The\n"
     "starting point of a full propagation path can be the top of the\n"
     "atmosphere, a blackbody ground (if *blackbody_ground* = 1) and the\n"
     "cloud box.\n"
     "\n"
     "A *ppath_step* agenda shall put in points along the propagation path\n"
     "at all crossings with the grids, tangent points and points of ground\n"
     "reflection. It is also allowed to make agendas that put in additional\n"
     "points to fulfil some criterion, such as a maximum distance along\n"
     "the path between the points. Accordingly, the number of points of\n"
     "each step can exceed one.\n"
     "\n"
     "For more information read the chapter on propagation paths in the\n"
     "ARTS user guide.\n"
     "\n"
     "Usage:      Called from *ppathCalc* and scattering functions.",
     Agenda_ ));

  wsv_data.push_back
   (WsvRecord
    ("r_geoid",
     "Geoid radius.\n"
     "\n"
     "Geometrical altitudes are defined as the vertical distance above the\n"
     "geoid, and the geoid is the reference surface used when giving, for\n"
     "example, *z_ground* and *z_field*. \n"
     "\n"
     "The geoid is defined by giving the radius from the coordinate centre\n"
     "to the geoid surface for each crossing of the latitude and longitude\n"
     "grids. The geoid should normally be selected to be an ellipsoid but\n"
     "any shape is allowed. For 1D calculations, the geoid is be definition\n"
     "a sphere.\n"
     "\n"
     "The radius for a point between the grid crossings is obtained by\n"
     "linear (2D) or bi-linear (3D) interpolation of the *r_geoid*. For 1D\n"
     "cases the geoid radius is constant.\n"
     "\n"
     "See further the ARTS user guide (AUG). Use the index to find where\n"
     "this variable is discussed. The variable is listed as a subentry to\n"
     "\"workspace variables\".\n"
     "\n"
     "Usage:      Set by using a function for a geodetic datum.\n"
     "\n"
     "Units:      m\n"
     "\n"
     "Dimensions: [ lat_grid, lon_grid ]",
     Matrix_ ));

  wsv_data.push_back
   (WsvRecord
    ("sensor_pos",
     "The sensor position for each measurement block.\n"
     "\n"
     "Text will be written. Right now I just want to make WSF for this\n"
     "variable, that also will be used for other purposes.\n"
     "\n"
     "Range for zenith angles is [0,180] beside for 2D where it\n"
     "is [-180,180].\n"
     "\n"
     "Dimensions: [ number of measurement blocks, atmospheric dimensionality ]",
     Matrix_ ));

  wsv_data.push_back
   (WsvRecord
    ("t_field",
     "The field of atmospheric temperatures.\n"
     "\n"
     "This variable gives the atmospheric temperature at each crossing of\n"
     "the pressure, latitude and longitude grids.\n"
     "\n"
     "The temperature for a point between the grid crossings is obtained by\n"
     "(multi-)linear interpolation of the *t_field*.\n"
     "\n"
     "See further the ARTS user guide (AUG). Use the index to find where\n"
     "this variable is discussed. The variable is listed as a subentry to\n"
     "\"workspace variables\".\n"
     "\n"
     "Usage:      Set by the user by interpolation of a climatology.\n"
     "\n"
     "Units:      K\n"
     "\n"
     "Dimensions: [ p_grid, lat_grid, lon_grid ]",
     Tensor3_ ));

  wsv_data.push_back
   (WsvRecord
    ("t_ground",
     "Effective emission temperature of the ground.\n"
     "\n"
     "The emission of the ground is calculated as the product between\n"
     "*e_ground* and the Planck function for the temperature given by\n"
     "*t_ground*. The ground emissivity and temperature are interpolated to\n"
     "the point of interest before the product is calculated.\n"
     "\n"
     "The temperatures can be set to any values >= 0 K. \n"
     "\n"
     "The temperature for a point between the grid crossings is obtained by\n"
     "linear (1D) or bi-linear (2D) interpolation of the *t_ground*.\n"
     "\n"
     "See further the ARTS user guide (AUG). Use the index to find where\n"
     "this variable is discussed. The variable is listed as a subentry to\n"
     "\"workspace variables\".\n"
     "\n"
     "Usage:      Set by user.\n"
     "\n"
     "Units:      K\n"
     "\n"
     "Dimensions: [ lat_grid, lon_grid ]",
     Matrix_ ));

  wsv_data.push_back
   (WsvRecord
    ("vector1",
     "An arbitrary vector.\n"
     "\n"
     "This variable shall be treated as a general variable of type Vector.\n"
     "It can be used, for example, when some intermediate data must be\n"
     "generated or to copy some data.\n"
     "\n"
     "Usage:      Set by user.",
     Vector_ ));

  wsv_data.push_back
   (WsvRecord
    ("vector2",
     "An arbitrary vector.\n"
     "\n"
     "This variable shall be treated as a general variable of type Vector.\n"
     "It can be used, for example, when some intermediate data must be\n"
     "generated or to copy some data.\n"
     "\n"
     "Usage:      Set by user.",
     Vector_ ));

  wsv_data.push_back
   (WsvRecord
    ("z_field",
     "The field of geometrical altitudes.\n"
     "\n"
     "This variable gives the geometrical altitude, above the geoid, of\n"
     "each crossing of the pressure, latitude and longitude grids. For 1D\n"
     "cases the altitudes give the geometrical position of the pressure\n"
     "surfaces.\n"
     "\n"
     "For each geographical position (lat,lon), the values must be sorted\n"
     "in increasing order, with no repetitions. Otherwise the altitudes can\n"
     "be set to arbitrary values. Hydrostatic equilibrium is not applied \n"
     "automatically. If hydrostatic equilibrium applies, *z_field* must be\n"
     "set by a function ensuring that this criterium is fulfilled.\n"
     "\n"
     "The radius (from the coordinate centre) for a point between the grid\n"
     "crossings is obtained by a (multi-)linear interpolation of the sum of\n"
     "*r_geoid* and *z_field*.\n" 
     "\n"
     "See further the ARTS user guide (AUG). Use the index to find where\n"
     "this variable is discussed. The variable is listed as a subentry to\n"
     "\"workspace variables\".\n"
     "\n"
     "Usage:      Set by the user by interpolation of a climatology, or set\n"
     "            by a function for applying hydrostatic equilibrium.\n"
     "\n"
     "Units:      m\n"
     "\n"
     "Dimensions: [ p_grid, lat_grid, lon_grid ]",
     Tensor3_ ));

  wsv_data.push_back
   (WsvRecord
    ("z_ground",
     "The ground altitude.\n"
     "\n"
     "This variable defines the shape of the ground, by giving the\n"
     "geometrical altitude above the geiod for each crossing of the \n"
     "latitude and longitude grids. Any shape of the ground is accepted.\n"
     "No gap between the ground and the lowermost pressure level is \n"
     "allowed.\n"
     "\n"
     "The radius (from the coordinate centre) for a point between the grid\n"
     "crossings is obtained by a linear (1D or bi-linear (2D) interpolation\n"
     "of the sum of *r_geoid* and *r_ground*. With other words, the radius\n"
     "for the ground is assumed to vary linear along the latitudes and \n"
     "longitudes in *lat_grid* and *lon_grid*.\n"
     "\n"
     "See further the ARTS user guide (AUG). Use the index to find where\n"
     "this variable is discussed. The variable is listed as a subentry to\n"
     "\"workspace variables\".\n"
     "\n"
     "Usage:      Set by user.\n"
     "\n"
     "Units:      m\n"
     "\n"
     "Dimensions: [ lat_grid, lon_grid ]",
     Matrix_ ));













  // Below this line are the definitions from ARTS-1. Shall be removed or
  // updated with better derscriptions.
  // =======================================================================
  // 


  //--------------------< Spectroscopy Stuff >--------------------
  //                     --------------------
  wsv_data.push_back
    (WsvRecord
     ("lines",
      "A list of spectral line data.", 
      ArrayOfLineRecord_));

  wsv_data.push_back
    (WsvRecord
     ("lines_per_tg",
      "A list of spectral line data for each tag.\n"
      "Dimensions: (tag_groups.nelem()) (# of lines for this tag)", 
      ArrayOfArrayOfLineRecord_));

  wsv_data.push_back
    (WsvRecord
     ("tgs",
      "This is an array of arrays of OneTag tag definitions.\n"
      "It defines the available tag groups for the calculation\n"
      "of absorption coefficients and weighting functions.\n"
      "Contrary to the original Bredbeck definition, tags within a\n"
      "group must belong to the same species, because one VMR profile\n"
      "is associated with each tag group.", 
      TagGroups_));

  wsv_data.push_back
    (WsvRecord
     ("wfs_tgs",
      "This is an array of arrays of tag group definitions.\n"
      "It defines the tag groups for the calculation of weighting\n"
      "functions. The selected tag groups must be a subgroup of the\n"
      "tag groups defined for the absorption coefficient calculation.", 
      TagGroups_));

  wsv_data.push_back
    (WsvRecord
     ("lineshape",
      "Lineshape specification: function, norm, cutoff. There is one entry for\n"
      "each abs_tag, not for each species. This means if you have several\n"
      "abs_tags for different isotopes or transitions of a species, you\n"
      "may use different lineshapes.",
      ArrayOfLineshapeSpec_));


  //--------------------< Continuum Stuff >--------------------
  //                     -----------------
  wsv_data.push_back
    (WsvRecord
     ("cont_description_names",
      "Continuum / full model absorption tag names. This variable should\n"
      "contain a list of tag names of continuum and full models, respectively.\n"
      "Associated with this WSV is the WSV\n"
      "`cont_description_models' which contains the specific model version of\n"
      "each continuum / full model absorption tag and the WSV\n"
      "`cont_description_parameters' which should contain the continuum / full model\n"
      "user defined parameters. The user defined parameters are only used when\n"
      "the specified model is 'user'. See also the online documentation in\n"
      "arts/doc/doxygen/html/continua_cc.html.\n"
      "\n"
      "The following full water vapor models are implemented:\n"
      "'H2O-MPM87': absorption model (line and continuum) according to \n"
      "   H. J. Liebe,\n" 
      "   A contribution to modeling atmospheric millimeter-wave properties,\n"
      "   Frequenz,  41, 1987, 31-36\n"
      "   and\n"
      "   H. J. Liebe and D. H. Layton,\n"
      "   Millimeter-wave properties of the atmosphere:\n"
      "   Laboratory studies and propagation modeling,\n"
      "   U.S. Dept. of Commerce, National Telecommunications and Information\n"
      "   Administration, Institute for Communication Sciences,\n"
      "   325 Broadway, Boulder, CO 80303-3328, report 87224.\n"
      "'H2O-MPM89': absorption model (line and continuum) according to \n"
      "   H. J. Liebe,\n Int. J. Infrared and Millimeter Waves, 10(6), 1989, 631\n"
      "'H2O-MPM93': absorption model (line and continuum) according to \n"
      "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
      "   Propagation modeling of moist air and suspended water/ice\n"
      "   particles at frequencies below 1000 GHz,\n"
      "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
      "   Propagation Panel,\n Palma de Mallorca, Spain, 1993, May 17-21 \n" 
      "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"
      "'H2O-CP98': absorption model (line and continuum) according to \n"
      "   S. L. Cruz-Pol et al.,\n Radio Science, 33(5), 1319, 1998"
      "   (ece.uprm.edu/~pol/Atmosphere.html)\n"
      "'H2O-PWR98': absorption model (line and continuum) according to \n"
      "   P. W. Rosenkranz,\n "
      "   Radio Science, 33(4),  919, 1998, Radio Science, 34(4), 1025, 1999\n"
      "   (ftp: mesa.mit.edu/phil/lbl_rt).\n"
      "\n"
      "The following full oxygen models are implemented:\n"
      "'O2-MPM93': absorption model (line and continuum) according to\n"
      "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
      "   Propagation modeling of moist air and suspended water/ice\n"
      "   particles at frequencies below 1000 GHz,\n"
      "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
      "   Propagation Panel,\n Palma de Mallorca, Spain, 1993, May 17-21\n"
      "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"
      "'O2-PWR93': absorption model (line and continuum) according to \n"
      "   P. W. Rosenkranz,\n Chapter 2, in M. A. Janssen, \n"
      "   Atmospheric Remote Sensing by Microwave Radiometry\n"
      "   John Wiley & Sons, Inc., 1993 (mesa.mit.edu/phil/lbl_rt)\n"
      "\n"
      "The following continuum parameterizations are implemented:\n"
      "H2O-H2O ('H2O-SelfContStandardType'):\n" 
      "   P. W. Rosenkranz, \n"
      "   Radio Science, Vol. 33, No 4, Pages 919-928, 1998 and \n"
      "   Radio Science, Vol. 34, No 4, Page 1025, 1999 (mesa.mit.edu/phil/lbl_rt)\n"
      "H2O-air ('H2O-ForeignContStandardType'): \n"
      "   P. W. Rosenkranz, \n"
      "   Radio Science, Vol. 33, No 4, Pages 919-928, 1998 and \n"
      "   Radio Science, Vol. 34, No 4, Page 1025, 1999 (mesa.mit.edu/phil/lbl_rt)\n"
      "H2O-air ('H2O-ContMPM93'): \n"
      "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
      "   Propagation modeling of moist air and suspended water/ice\n"
      "   particles at frequencies below 1000 GHz,\n"
      "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
      "   Propagation Panel,\n Palma de Mallorca, Spain, 1993, May 17-21\n"
      "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"      
      "O2-air ('O2-SelfContStandardType'):\n"
      "   P. W. Rosenkranz,\n"
      "   Chapter 2, in M. A. Janssen,\n"
      "   Atmospheric Remote Sensing by Microwave Radiometry,\n"
      "   John Wiley & Sons, Inc., 1993\n"
      "   (mesa.mit.edu/phil/lbl_rt)\n"
      "   and also described in \n"
      "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
      "   Propagation modeling of moist air and suspended water/ice\n"
      "   particles at frequencies below 1000 GHz,\n"
      "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
      "   Propagation Panel,\n Palma de Mallorca, Spain, 1993, May 17-21\n"
      "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"
      "N2-N2 ('N2-SelfContStandardType'):\n"
      "   The functional form of Rosenkranz but with more input parameters.\n"
      "   P. W. Rosenkranz,\n"
      "   Chapter 2, in M. A. Janssen,\n"
      "   Atmospheric Remote Sensing by Microwave Radiometry,\n"
      "   John Wiley & Sons, Inc., 1993 (mesa.mit.edu/phil/lbl_rt)\n"
      "N2-N2 ('N2-SelfContMPM93'):\n"
      "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
      "   Propagation modeling of moist air and suspended water/ice\n"
      "   particles at frequencies below 1000 GHz,\n"
      "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
      "   Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21 \n"
      "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"
      "CO2-CO2 ('CO2-SelfContPWR93'):\n"
      "   P. W. Rosenkranz,\n"
      "   Chapter 2, in M. A. Janssen,\n"
      "   Atmospheric Remote Sensing by Microwave Radiometry,\n"
      "   John Wiley & Sons, Inc., 1993 (mesa.mit.edu/phil/lbl_rt)\n"
      "CO2-N2 ('CO2-ForeignContPWR93'):\n"
      "   P. W. Rosenkranz,\n"
      "   Chapter 2, in M. A. Janssen,\n"
      "   Atmospheric Remote Sensing by Microwave Radiometry,\n"
      "   John Wiley & Sons, Inc., 1993 (mesa.mit.edu/phil/lbl_rt)\n"
      "\n"
      "The following cloud absorption models are implemented:\n"
      "Suspended water droplet ('liquidcloud-MPM93') \n"
      "   absorption parameterization from the MPM93 model:\n"
      "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
      "   Propagation modeling of moist air and suspended water/ice\n"
      "   particles at frequencies below 1000 GHz,\n"
      "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
      "   Propagation Panel,\n Palma de Mallorca, Spain, 1993, May 17-21\n"
      "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"
      "Ice water droplet absorption ('icecloud-MPM93') \n"
      "   parameterization from MPM93 model:\n"
      "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
      "   Propagation modeling of moist air and suspended water/ice\n"
      "   particles at frequencies below 1000 GHz,\n"
      "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
      "   Propagation Panel,\n Palma de Mallorca, Spain, 1993, May 17-21\n"
      "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"
      "\n",
      ArrayOfString_));
  
  wsv_data.push_back
    (WsvRecord
     ("cont_description_models",
      "Continuum / full model absorption model description parameter.\n"
      "See the WSV `cont_description_names' for a detailed description\n"
      "of the allowed continuum models. There should be one string here\n"
      "for each entry in `cont_description_names'.See also the online" 
      "documentation in arts/doc/doxygen/html/continua_cc.html.\n",
      ArrayOfString_));

  wsv_data.push_back
    (WsvRecord
     ("cont_description_parameters",
      "Continuum model parameters. See the WSV `cont_description_names'\n"
      "for a detailed description of the allowed continuum models. There\n"
      "should be one parameter vector here for each entry in\n"
      "`cont_description_names'. See also the online documentation in\n"
      "arts/doc/doxygen/html/continua_cc.html.\n",
      ArrayOfVector_));


  //--------------------< 1D Input Atmosphere Stuff >--------------------
  //                     ---------------------------
  wsv_data.push_back
    (WsvRecord
     ("raw_ptz",
      "Matrix has rows:\n"
      "1. Pressure in Pa\n"
      "2. Temperature in K\n"
      "3. Altitude in m", 
      Matrix_));

  wsv_data.push_back
    (WsvRecord
     ("raw_vmrs",
      "The individual VMR profiles. Each species VMR profile comes with a\n"
      "pressure profile. The different species can hence be on different\n"
      "grids.\n"
      "The matrix has rows:\n"
      "1. Pressure in Pa\n"
      "2. VMR profile (absolute number)\n"
      "The array dimension is determined by the number of tag groups.", 
      ArrayOfMatrix_));


  //--------------------< General Absorption Stuff >--------------------
  //                     --------------------------
  wsv_data.push_back
    (WsvRecord
     ("p_abs",
      "The pressure grid for the absorption coefficients [Pa]. This\n"
      "is the basic independent grid for the absorption calculation, both\n"
      "in the 1D and 2D case. Therefore it remains a vector, even in 2D.\n"
      "The \"raw\" atmospheric data shall be interpolated to p_abs before\n"
      "the absorption calculations starts.",
      Vector_));
  
  wsv_data.push_back
    (WsvRecord
     ("f_mono",
      "The monochromatic frequency grid [Hz]. ",
      Vector_));
    

  //--------------------< 1D Absorption Stuff >--------------------
  //                     ---------------------
  wsv_data.push_back
    (WsvRecord
     ("t_abs",
      "Temperature associated with the pressures in p_abs [K]",
      Vector_));

  wsv_data.push_back
    (WsvRecord
     ("z_abs",
      "Vertical altitudes associated with the pressures in p_abs [m]",
      Vector_));

  wsv_data.push_back
    (WsvRecord
     ("h2o_abs",
      "The total water profile associated with the pressures in p_abs [-]",
      Vector_));

  wsv_data.push_back
    (WsvRecord
     ("n2_abs",
      "The total nitrogen profile associated with the pressures in p_abs [-]",
      Vector_));

  wsv_data.push_back
    (WsvRecord
     ("vmrs",
      "The VMRs (unit: absolute number) on the p_abs grid.\n"
      "Dimensions: [tag_groups.nelem(), p_abs.nelem()]",
      Matrix_));

  wsv_data.push_back
    (WsvRecord
     ("abs",
      "The matrix of absorption coefficients (in units of [1/m]).\n"
      "Dimensions: [f_mono.nelem(), p_abs.nelem()]",
      Matrix_));

  wsv_data.push_back
    (WsvRecord
     ("abs_per_tg",
      "These are the absorption coefficients individually for each\n"
      "tag group. The Array contains one matrix for each tag group,\n"
      "the matrix format is the same as that of abs",
      ArrayOfMatrix_));

  wsv_data.push_back
    (WsvRecord
     ("xsec_per_tg",
      "These are the cross sections individually for each tag\n"
      "group. The Array contains one matrix for each tag group,\n"
      "the matrix format is the same as that of abs",
      ArrayOfMatrix_));


}
