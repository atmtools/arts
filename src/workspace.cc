/* Copyright (C) 2000, 2001, 2002 Stefan Buehler <sbuehler@uni-bremen.de>
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
*/


#include "arts.h"
#include "matpackI.h"
#include "matpackIII.h"
#include "array.h"
#include "auto_wsv_groups.h"
#include "wsv_aux.h"
#include "ppath.h"

// Some #defines to make the records better readable:
#define NAME(x)        x 
#define DESCRIPTION(x) x
#define GROUP(x)       x 


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

  Usage:      Set by user (or "Method output.")

  Units:      E.g., kg/m

  Dimensions: [ first dimension, second dimension, ... ]
  or
  Size:       [ .., nrows, ncols ]

  Members:    Here you would list the members if your
              variable is a structure.
  ----------------------------------------------------------------------

  The dashed lines are not part of the template, they just show you
  where it starts and ends. 

  Give the keywords above only if they apply, i.e., Members only
  for a structure, Units only for a physical variable. 
  Use either Dimensions or Size, depending on what is most appropiate for
  the variable.
  */

  
  /*----------------------------------------------------------------------
    Let's put in the variables in alphabetical order. This gives a clear
    rule for where to place a new variable and this gives a nicer
    results when the methods are listed by "arts -w all".  No
    distinction is made between uppercase and lowercase letters. The
    sign "_" comes after all letters.
    Patrick Eriksson 2002-05-08
  ----------------------------------------------------------------------*/

  wsv_data.push_back
    (WsvRecord
    ( NAME( "abs_vec" ),
      DESCRIPTION
      (
       "Total absorption vector.\n"
       "\n"
       "This variable contains the absorption coefficient vector which \n"
       "is used in the scattering RT calculation. It is \n"
       "the physical absorption which includes particle absorption \n"
       "for all chosen particle types as well as gaseous absorption for\n"
       "all chosen gaseous species.\n" 
       "The vector is calculated by the agenda *abs_vec_agenda* \n"
       "The dimensision of the variable adapts to *stokes_dim*.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      Output of the agenda *ext_mat_agenda*. \n"
       "\n"
       "Unit:        m^2 \n"
       "\n"
       "Dimensions: [stokes_dim]"
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_vec_agenda" ),
       DESCRIPTION
       (
	"See agendas.cc."
	),
       GROUP(  Agenda_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_vec_gas" ),
       DESCRIPTION
       (
	"Absorption Vector for gaseous species.\n"
	"\n"
	"This workspace variable represents the absorption vector of the\n"
	"gaseous species specified for the study.  It can either be\n"
	"calculated in the method *abs_vec_gasCalc* or can be read in from\n"
	"from a file.  The method *abs_vec_gasCalc* is not yet implemented.\n"
	"\n"
	"There will be for sure more documentation on this.\n"
	"\n"
	"Usage:      Input to the method abs_vecCalc\n"
	"\n"
	"Unit:        m^2\n"
	"\n"
	"Dimensions: [ stokes_dim ]\n"
	),
      GROUP( Vector_ )));
  

  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_vec_part" ),
       DESCRIPTION
       (
	"Physical absorption vector for particle.\n"
	"\n"
	"This workspace variable represents the actual physical absorption\n"
	"vector of the particles chosen for the study for given propagation \n"
	"directions.  This vector is calculated by the method *abs_vec_partCalc*\n"
	"\n"
	"ARTS user guide (AUG) gives the formula used for computing this \n"
	"variable. Use the index to find where this variable is discussed.\n"
	"The variable is listed as a subentry to \"workspace variables\".\n"
        ),
      GROUP( Vector_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME("abs_vec_spt"),
       DESCRIPTION
       (
	"Absorption Vector for a single particle type.\n"
	"\n"
	"This variable contains the elements of absorption vector of a \n"
	"single particle, given *ext_mat_spt* and *pha_mat_spt*. It is the\n"
	"input as well as the output of the method *abs_vec_sptCalc*. This \n"
	"variable is a matrix where the first dimension part_types indicate \n"
	"the particle type under consideration and the second dimension is \n"
        "the *stokes_dim*.\n"
        "*stokes_dim* can be 1,2,3 or 4 as as set by the user.\n"
	"\n"
	"ARTS user guide (AUG) gives the formulas used for computing all \n"
	"the elements of absorption vector.\n"
	"\n"
	"Usage:      Input and Output of the method abs_vec_sptCalc\n"
	"\n"
	"Unit:        m^2\n"
	"\n"
	"Dimensions: [part_types,stokes_dim]"
	),
       GROUP( Matrix_ ) ));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "amp_mat" ),
      DESCRIPTION
      (
       "Monochromatic amplitude matrix.\n"
       "\n"
       "The amplitude matrix is required for scattering calculations.\n"
       "It contains all optical properties of the scattering particles.\n"
       "It depends on the frequency, the particle type, the propagation \n"
       "direction and the scattered direction.\n"
       "*amp_mat* is generated inside the frequency loop, so the variable \n"
       "itself does not include the frequency information. \n"
       "The amplitude matrix is a 2x2 complex matrix. The workspace variable\n"
       "*amp_mat* stores the real and imaginary elements (i.e. 8 elements)\n"
       "separately. \n" 
       "\n"
       "Usage: Input to *ext_mat_agenda*, *abs_vec_agenda*, \n"
       "*sca_mat_agenda*\n"
       "Output of *get_amp*. \n"    
       "\n"
       "Unit:       m\n"
       "\n"
       "Dimensions: [ part_types, scat_za_grid, scat_aa_grid, \n"
       "              scat_za_grid, scat_aa_grid, amplitude matrix element]"
       ), 
      GROUP( Tensor6_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "amp_mat_raw" ),
      DESCRIPTION
      (
       "Amplitude matrix data.\n"
       "\n"
       "The amplitude matrix is required for scattering calculations.\n"
       "It contains all optical properties of the scattering particles.\n"
       "It depends on the frequency, the particle type, the propagation \n"
       "direction and the scattered direction.\n"
       "*amp_mat_raw* is an Array of Tensor6. Each Tensor6 corresponds \n"
       "to one particle type and contains the amplitude matrix for all \n"
       "frequencies defined in *f_grid* and all prpagation and \n"
       "scattering directions.\n" 
       "Furthermore it contains the grids where the data is stored on, i.e. \n"
       "*f_grid*, *scat_za_grid*, *scat_aa_grid*, *scat_za_grid* and \n"
       "*scat_aa_grid*. These are stored in the first element in each \n" 
       "dimension, therefore the size of the frequency dimension for \n"
       "example is the number of frequencies plus one. \n"
       "The amplitude matrix is a 2x2 complex matrix. The workspace variable\n"
       "*amp_mat_raw* stores the real and imaginary elements\n"
       "(i.e. 8 elements) separately. \n" 
       "\n"
       "Usage: not clear yet \n"     
       "\n"
       "Unit: m \n"
       "\n"
       "Size: Array[N_pt] \n "
       "      [N_f+1, N_za+1, N_aa+1, N_za+1, N_aa+1, 8] \n"
       "\n"
       "Dimensions: Array [part_types] \n"
       "            [f_grid, scat_za_grid, scat_aa_grid, scat_za_grid, \n"
       "             scat_aa_grid, amplitude matrix element]\n"
       "\n"
       ),
      GROUP(ArrayOfTensor6_ )));

  
  wsv_data.push_back
   (WsvRecord
    ( NAME( "antenna_dim" ),
      DESCRIPTION
      (
       "The dimensionality of the antenna pattern (1-2).\n"
       "\n"
       "Text will be written (PE).\n"
       "\n"
       "Usage:      Set by the user."
       ),
      GROUP( Index_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "atmosphere_dim" ),
      DESCRIPTION
      (
       "The atmospheric dimensionality (1-3).\n"
       "\n"
       "This variable defines the complexity of the atmospheric structure.\n"
       "The dimensionality is given by an integer between 1 and 3, where 1\n"
       "means 1D etc. This is the master variable for the atmospheric\n"
       "dimensionality, variables which size changes with the dimensionality\n"
       "are checked to match this variable. \n"
       "\n"
       "Methods adapt automatically to this variable. That is, it should\n"
       "not be needed to change any methods if the dimensionality is\n"
       "changed.\n"
       "\n"
       "The atmospheric dimensionalities (1D, 2D and 3D) are defined in the\n"
       "user guide (look for \"atmospheric dimensionality\" in the index).\n" 
       "\n"
       "Usage:      Set by the user."
       ),
      GROUP( Index_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "meridian_angle_1d" ),
      DESCRIPTION
      (
       "The meridian angle for a 1D observation.\n"
       "\n"
       "This variable is only used when calculating a curvature radius of the\n"
       "geoid for 1D cases. The given value shall be the angle between the\n"
       "observation and meridian plane where 0 degrees means an observation in\n"
       "the N-S direction.\n"
       "\n"
       "Values shall be in the range [-180,180].\n"
       "\n"
       "Usage:      Set by the user.\n"
       "\n"
       "Unit:       degrees"
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "a_los" ),
      DESCRIPTION
      (
       "A line-of-sight (for test purposes).\n"
       "\n"
       "The purpose of this variable and *a_pos* are to enable calling of\n"
       "*ppathCalc* and maybe other methods from the workspace. This can be\n"
       "of interest both for testing of the methods and to display\n"
       "intermediate results as the propagation path. \n"
       "\n"
       "For 1D and 2D cases, *a_los* is a vector of length 1 holding the \n"
       "zenith angle. For 3D, the length of the vector is 2, where the\n"
       "additional element is the azimuthal angle. These angles are defined\n"
       "in the ARTS user guide (AUG). Look in the index for \"zenith angle\"\n"
       "and \"azimuthal angle\".\n"
       "\n"
       "Usage: To call *ppathCalc* as an individual method. Other\n"
       "       purposes can be possible.\n"
       "\n"
       "Units: [ degree, degree ]\n"
       "\n"
       "Size:  [ 1 or 2 ]"
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "a_pos" ),
      DESCRIPTION
      (
       "A geographical position (for test purposes).\n"
       "\n"
       "The purpose of this variable and *a_los* are to enable calling of\n"
       "*ppathCalc* and maybe other methods from the workspace. This can be\n"
       "of interest both for testing of the methods and to display\n"
       "intermediate results as the propagation path. \n"
       "\n"
       "This variable is a vector with a length equalling the atmospheric\n"
       "dimensionality. The first element is the radius (from the coordinate\n"
       "system centre) of the position. Element 2 is the latitude and element\n"
       "3 is the longitude. Please note that the vertical position is given\n"
       "as the radius, not the altitude above the geoid.\n"
       "\n"
       "Usage: To call *ppathCalc* as an individual method. Other\n"
       "       purposes can be possible.\n"
       "\n"
       "Units: [ m, degree, degree ]\n"
       "\n"
       "Size:  [ atmosphere_dim ]"
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "blackbody_ground" ),
      DESCRIPTION
      (
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
       "Usage:      Set by user."
       ),
      GROUP( Index_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "cloudbox_on" ),
      DESCRIPTION
      (
       "Flag to activate the cloud box.\n"
       "\n"
       "Scattering calculations are confined to a part of the atmosphere\n"
       "denoted as the cloud box. The extension of the cloud box is given by\n"
       "*cloudbox_limits*. This variable tells methods if a cloud box is\n"
       "activated or not. \n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      Set by the user."
       ),
      GROUP( Index_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "cloudbox_limits" ),
      DESCRIPTION
      (
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
       "Usage: Set by the user, either directly or using a method\n"
       "       checking the extension of scattering particles.\n"
       "\n"
       "Unit:  Index values.\n"
       "\n"
       "Size:  [ 2 * atmosphere_dim ]"
       ),
      GROUP( ArrayOfIndex_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "convergence_flag" ),
      DESCRIPTION
      (
       "Flag for the convergence test.\n"
       "\n"
       "This variable is initialized with 0 inside the method \n"
       "*i_fieldIterate*.\n"
       "If after an iteration the convergence test is fulfilled, 1 is \n"
       "assigned which means that the iteration is completed. \n"
      ), 
      GROUP( Index_ ))); 

 wsv_data.push_back
   (WsvRecord
    ( NAME( "convergence_test_agenda" ),
      DESCRIPTION
      (
	"See agendas.cc."
       ),
      GROUP( Agenda_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "e_ground" ),
      DESCRIPTION
      (
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
       "Unit:       - (0-1)\n"
       "\n"
       "Dimensions: [ f_grid, lat_grid, lon_grid ]"
       ),
      GROUP( Tensor3_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "els" ),
      DESCRIPTION
      (
       "The elementary lineshape function.\n"
       "\n"
       "Holds the result of a method calculating the elementary lineshape\n"
       "function for a vector of frequencies. The difference to the lineshape\n"
       "functions ls is that the frequency grid associated with els is always\n"
       "relative to the line center.\n"
       "\n"
       "Usage: Agenda output, set by elementary lineshape\n"
       "       functions, e.g., elsLorentz.\n"
       "\n"
       "Unit: 1/Hz."
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "els_agenda" ),
       DESCRIPTION
       (
	"See agendas.cc."
	),
       GROUP(  Agenda_ )));
  
  wsv_data.push_back
    (WsvRecord
    ( NAME( "ext_mat" ),
      DESCRIPTION
      (
       "Total extinction matrix.\n"
       "\n"
       "This variable contains the extinction coefficient matrix which \n"
       "is used in the scattering RT calculation. It is \n"
       "the physical extinction matrix which includes particles extinction \n"
       "for all chosen particle types and gaseous extinction for all chosen \n"
       "gaseous species.\n" 
       "The matrix is calculated by the agenda *ext_mat_agenda* \n"
       "The dimensision of the variable adapts to *stokes_dim*.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      Output of the agenda *ext_mat_agenda*. \n"
       "\n"
       "Unit:        m^2 \n"
       "\n"
       "Dimensions: [stokes_dim, stokes_dim]"
       ),
      GROUP( Matrix_ )));
  
  wsv_data.push_back
    (WsvRecord
     ( NAME( "ext_mat_agenda" ),
       DESCRIPTION
       (
	"See agendas.cc."
	),
       GROUP(  Agenda_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "ext_mat_gas" ),
       DESCRIPTION
       (
	"Extinction matrix for gaseous species.\n"
	"\n"
	"This workspace variable represents the extinction matrix of the\n"
	"gaseous species specified for the study.  It can either be\n"
	"calculated in the method *ext_mat_gasCalc* or can be read in from\n"
	"a file.  The method *ext_mat_gasCalc* is not yet implemented. \n"
	"\n"
	"There will be for sure more documentation on this.\n"
	"\n"
	"Usage:      Input to the method ext_matCalc\n"
	"\n"
	"Unit:        m^2\n"
	"\n"
	"Dimensions: [ stokes_dim, stokes_dim ]\n"
	),
      GROUP( Matrix_ )));
	
  wsv_data.push_back
    (WsvRecord
     ( NAME( "ext_mat_part" ),
       DESCRIPTION
       (
	"Physical extinction matrix for particle.\n"
	"\n"
	"This workspace variable represents the actual physical extinction\n"
	"matrix of the particles chosen for the study for given propagation \n"
	"directions.  This is calculated by the method *ext_mat_partCalc*\n"
	"\n"
	"ARTS user guide (AUG) gives the formula used for computing this \n"
	"variable. Use the index to find where this variable is discussed.\n"
	"The variable is listed as a subentry to \"workspace variables\".\n"
	"\n"
	"Usage:      Output of the method ext_mat_partCalc\n"
	"\n"
	"Unit:        m^2\n"
	"\n"
	"Dimensions: [ stokes_dim, stokes_dim ]\n"
        ),
       GROUP( Matrix_ )));

  wsv_data.push_back
    (WsvRecord
    ( NAME( "ext_mat_spt" ),
      DESCRIPTION
      (
       "Extinction matrix for a single particle type.\n"
       "\n"
       "This variable contains the elements for extinction matrix of a  \n"
       "single particle for propagation direction given by *za_index* \n"
       "and *aa_index*. It is the input as well as the output of the \n"
       "method *ext_mat_sptCalc*.  The elements of extinction matrix \n"
       "are calculated from the elements of *amp_mat*. This variable \n"
       "comes under Tensor3 where the first dimension *part_types* \n"
       "indicate the particle type under consideration and the  \n"
       "second and third dimension indicates the *stokes_dim*.  \n"
       "*stokes_dim* can be 1,2,3 or 4 as set by the user. \n"
       "\n"
       "ARTS user guide (AUG) gives the formulae used for computing all \n"
       "the elements of the extinction matrix for a given particle  \n"
       "type. \n"
       "\n"
       "Usage:      Input and Output of the method ext_mat_sptCalc \n"
       "\n"
       "Unit:        m^2 \n"
       "\n"
       "Dimensions: [part_types,stokes_dim, stokes_dim]"
       ),
      GROUP( Tensor3_ )));
    
  wsv_data.push_back
    (WsvRecord
     ( NAME( "f_grid" ),
       DESCRIPTION
       (
	"The frequency grid for monochromatic pencil beam\n"
	"calculations.\n"
	"\n"
	"Text will be written (PE).\n" 
	"\n" 
	"Usage:      Set by the user.\n "
	"\n"   
	"Unit:        Hz"
	), 
       GROUP( Vector_ )));
 
  wsv_data.push_back
    (WsvRecord
     ( NAME( "i_field" ), 
       DESCRIPTION
      (
       "Radiation field.\n" 
       "\n"
       "This variable is used to store the intensity field inside the\n"
       "cloudbox which is found by an iterative solution.\n"
       "More decription will be written (CE).\n"
       "\n"
       "Usage: Input and output of *scat_mono_agenda*. \n"    
       "\n"
       "Unit: W / (m^2 Hz sr) for each Stokes component.\n"
       "\n"
       "Dimensions: [ p_grid, lat_grid, lon_grid, scat_za_grid, \n"
       "              scat_aa_grid, stokes_dim ]"
       ),
      GROUP( Tensor6_ )));
  
  wsv_data.push_back
    (WsvRecord
     ( NAME( "i_field_dim" ), 
       DESCRIPTION
      (
       "Dimension of the radiation field.\n" 
       "\n"
       "This variable is important if a 1D atmosphere is considered.\n"
       "1D atmosphere means that all profiles only depend on altitude or \n"
       "equivalently pressure. Nevertheless it can be useful to allow the \n"
       "radiation field to be 2D or 3D. Solar radiation or scattering can \n"
       "be sources of inhomogeneities in the radiation field.\n" 
       ),
      GROUP( Index_ )));
 
 wsv_data.push_back
   (WsvRecord
    ( NAME( "i_field_old" ),
      DESCRIPTION
      (
       "Intensity field inside the cloudbox.\n"
       "\n"
       "This variable is used to store the intensity field inside the\n"
       "cloudbox while performing the iteration.\n"
       "More decription will be written (CE).\n"
       "\n"
       "Usage: Input of *i_fieldUpdate1D*. \n"    
       "\n"
       "Unit: W / (m^2 Hz sr) for each Stokes component.\n"
       "\n"
       "Dimensions: [ p_grid, lat_grid, lon_grid, scat_za_grid, \n"
       "              scat_aa_grid, stokes_dim ]"
       ),
      GROUP( Tensor6_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "i_rte" ),
      DESCRIPTION
      (
       "The spectrum produced by *rte_agenda*.\n"
       "\n"
       "Text will be written (PE).\n"
       "\n"
       "Usage:      Output from *rte_agenda*.\n"
       "\n"
       "Unit:       W / (m^2 Hz sr) or optical thickness \n"
       "\n"
       "Dimensions: [ f_grid, stokes_dim ]"
       ),
      GROUP( Matrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "i_space" ),
      DESCRIPTION
      (
       "The monochromatic intensity entering at the top of the atmosphere.\n"
       "\n"
       "Text will be written (PE).\n"
       "\n"
       "This radiation is assumed to un-polarised, and a vector suffice.\n"
       "\n"
       "Usage:      Set by the user.\n"
       "\n"
       "Unit:       W / (m^2 Hz sr) \n"
       "\n"
       "Dimensions: [ f_grid ]"
       ),
      GROUP( Vector_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "l_step" ),
      DESCRIPTION
      (
       "The pathlegth through one grid cell/ layer.\n"
       "\n"
       "The pathlength is required in the methods for RT step calculations, \n"
       "which are *sto_vecGeneral* and *sto_vecScalar*.\n"
       "It can be calculated using the *ppath_step_agenda*.\n"
       "\n"
       "Usage:      Calculated in *i_fieldUpdate1D*."
       "\n"
       "Unit:       m \n"
       ),
      GROUP( Numeric_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "lat_grid" ),
      DESCRIPTION
      (
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
       "Unit:       degrees"
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "lat_1d" ),
      DESCRIPTION
      (
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
       "Unit:       degrees"
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "lon_grid" ),
      DESCRIPTION
      (
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
       "Unit:       degrees"
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ls" ),
      DESCRIPTION
      (
       "The lineshape function.\n"
       "\n"
       "Holds the result of a method calculating the lineshape function for a\n"
       "Vector of frequencies.\n"
       "\n"
       "Usage: Agenda output, set by lineshape functions, e.g., lsLorentz.\n"
       "\n"
       "Unit: 1/Hz."
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ls_cutoff" ),
      DESCRIPTION
      (
       "The lineshape cutoff frequency.\n"
       "\n"
       "The cutoff is meant to limit the frequency range where the lineshape\n"
       "is not zero. Method lsWithCutoffAdd uses this variable to set the\n"
       "lineshape to zero for frequencies for which\n"
       "\n"
       "abs(f-f0) > ls_cutoff\n"
       "\n"
       "Unit: Hz\n"
       "\n"
       "Usage: Set by the user. A typical value is 750e9 Hz."
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "els_f_grid" ),
      DESCRIPTION
      (
       "The frequency grid for the lineshape calculation.\n"
       "\n"
       "This is a local copy of the global f_grid. The copy is necessary,\n"
       "because in cases with cutoff we have to add the cutoff frequency to\n"
       "the frequency grid, so that we can subtract the lineshape value at the\n"
       "cutoff. \n"
       "\n"
       "Usage: Agenda input, set automatically by calling method.\n"
       "\n"
       "Unit: Hz"
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ls_f0" ),
      DESCRIPTION
      (
       "The line center frequency.\n"
       "\n"
       "Used as input by methods calculating the lineshape function.\n"
       "\n"
       "Usage: Agenda input, set automatically by calling method.\n"
       "\n"
       "Unit: Hz."
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ls_gamma" ),
      DESCRIPTION
      (
       "The pressure broadened line width.\n"
       "\n"
       "Used as input by methods calculating the lineshape function.\n"
       "\n"
       "Usage: Agenda input, set automatically by calling method.\n"
       "\n"
       "Unit: Hz."
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ls_sigma" ),
      DESCRIPTION
      (
       "The Doppler broadened line width.\n"
       "\n"
       "Used as input by methods calculating the lineshape function.\n"
       "\n"
       "Usage: Agenda input, set automatically by calling method.\n"
       "\n"
       "Unit: Hz."
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "main_agenda" ),
      DESCRIPTION
      (
	"See agendas.cc."
       ),
      GROUP( Agenda_)));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "mblock_aa_grid" ),
      DESCRIPTION
      (
       "The azimuthal angle grid for each measurement block.\n"
       "\n"
       "Text will be written (PE).\n"
       "\n"
       "Usage:      Set by the user.\n"
       "\n"
       "Unit:       degrees "
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "mblock_za_grid" ),
      DESCRIPTION
      (
       "The zenith angle grid for each measurement block.\n"
       "\n"
       "Text will be written (PE).\n"
       "\n"
       "Usage:      Set by the user.\n"
       "\n"
       "Unit:       degrees "
       ),
      GROUP( Vector_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "planck_function" ),
      DESCRIPTION
      (
       "The Planck function.\n"
       "\n"
       "The Planck function is used in the methods for RT step calculations, \n"
       "which are *sto_vecGeneral* and *sto_vecScalar*. \n"
       "\n"
       "Usage:      Calculated in *i_fieldUpdate1D*.\n"
       "\n"
       "Unit:       W\n "
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "p_grid" ),
      DESCRIPTION
      (
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
       "Unit:       Pa"
       ),
      GROUP( Vector_ )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "part_types" ),
      DESCRIPTION
      (
       "Particle types.\n"
       "\n"
       "A vector containing all particle types which shall be considered."
       "\n"
       "Usage:      Set by the user.\n"
       "\n"
       "Unit:       Pa"
       ),
      GROUP( Vector_ )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "pha_mat" ),
      DESCRIPTION
      (
       "Physical phase matrix for particle.\n"
       "\n"
       "This workspace variable represents the actual physical phase\n"
       "matrix of the particles chosen for the study for given propagation\n"
       "directions.  This is calculated by the method *pha_matCalc*\n"
       "\n"
       "ARTS user guide (AUG) gives the formula used for computing this\n"
       "variable. Use the index to find where this variable is discussed.\n"
       "The variable is listed as a subentry to \"workspace variables\".\n"
       "\n"
       "Usage:      Output of the method pha_matCalc\n"
       "\n"
       "Unit:        m^2\n"
       "\n"
       "Dimensions: [ scat_za_grid, scat_aa_grid, stokes_dim, stokes_dim ]\n"
       ),
      GROUP( Tensor4_ )));
   
   wsv_data.push_back
   (WsvRecord
    ( NAME( "pha_mat_spt" ),
      DESCRIPTION
      (
       "Phase matrix for a single particle type.\n"
       "\n"
       "This variable contains the elements of phase matrix for a single \n"
       "particle for given propagation direction *za_index* and *aa_index*. \n"
       "It is the input as well as the output of the method *pha_mat_sptCalc*.\n"
       "The elements of the phase matrix are calculated from the elements of \n"
       "*amp_mat*. This variable is a Tensor5 where the first dimension \n"
       "*part_types* indicate the particle type under consideration.The \n"
       "second and third dimension gives the incident zenith and azimuth \n"
       "angles respectively, and the fourth and the fifth dimension gives\n" 
       "the *stokes_dim*. *stokes_dim* can be 1,2,3 or 4 as set by the user.\n"
       "\n"
       "ARTS user guide (AUG) gives the formulas used for computing all \n"
       "elements of the phase matrix for a given particle type.\n"
       "\n"
       "Usage:      Input and Output of the method pha_mat_sptCalc\n"
       "\n"
       "Unit:        m^2\n"
       "\n"
       "Dimensions: [part_types,scat_za_grid,scat_aa_grid,stokes_dim, stokes_dim]"
       ),
      GROUP( Tensor5_ )));


   wsv_data.push_back
   (WsvRecord
    ( NAME( "pnd_field" ),
      DESCRIPTION
      (
       "The field representing particle number densities.\n"
       "\n"
       "This variable gives the particle number density of the chosen particle\n"
       "types as a function of p_grid, lat_grid, lon_grid. \n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "\n"
       "Usage:      Calculated internally.\n"
       "\n"
       "Unit:        m^-3\n"
       "\n"
       "Dimensions: [part_types, p_grid, lat_grid, lon_grid]\n"
        ),
      GROUP( Tensor4_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "pnd_field_raw" ),
      DESCRIPTION
      (
       "The particle number density field data.\n"
       "\n"
       "This variable contains the particle number densities for all \n"
       "chosen particle types. It includes the grids corresponding to the \n"
       "grids in the database. \n"
       "\n"
       "\n"
       "Usage:      Reading routine.\n"
       "\n"
       "Unit:        m^-3\n"
       "\n"
       "Size:  Array[N_pt] \n "
       "       [N_p, N_lat, N_lon] \n"
       "\n"
       "Dimensions: Array[part_types]\n" 
       "           [p_grid, lat_grid, lon_grid]\n"
        ),
      GROUP( ArrayOfTensor3_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath" ),
      DESCRIPTION
      (
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
       "Usage:      Output from the method *ppathCalc*.\n"
       "\n"
       "Members:    To be written."
       ),
      GROUP( Ppath_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath_step" ),
      DESCRIPTION
      (
       "A propagation path step.\n"
       "\n"
       "The main intention of this variable is communication with the agenda\n"
       "*ppath_step_agenda*.\n"
       "\n"
       "Type \"arts -d ppath_step\" and \"arts -d ppathCalc\" for more\n"
       "information on this variable and the calculation of propagation\n"
       "paths. Or read the chapter on propgation paths in the ARTS user\n"
       "guide.\n"
       "\n"
       "Usage:      In/output to/from *ppath_step_agenda*.\n"
       "\n"
       "Members:    To be written."
       ),
      GROUP( Ppath_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath_step_agenda" ),
      DESCRIPTION
      (
	"See agendas.cc."
       ),
      GROUP( Agenda_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_agenda" ),
      DESCRIPTION
      (
	"See agendas.cc."
       ),
      GROUP( Agenda_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "r_geoid" ),
      DESCRIPTION
      (
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
       "Usage:      Set by using a method for a geodetic datum.\n"
       "\n"
       "Unit:       m\n"
       "\n"
       "Dimensions: [ lat_grid, lon_grid ]"
       ),
      GROUP( Matrix_ )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_aa_grid" ),
      DESCRIPTION
      (
       "Azimuthal angle grid.\n"
       "\n"
       "The azimutal angle grid, on which the intensity field and the \n"
       "optical scattering properties are stored. The grid has to be defined\n"
       "if the cloudbox is activated by the flag *cloudbox_on*.\n"
       "The grid must be sorted in decreasing order, with no repetitions.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      Set by the user.\n"
       "\n"
       "Unit:       degrees "
       ),
      GROUP( Vector_ )));

 wsv_data.push_back
    (WsvRecord
    ( NAME( "sca_vec" ),
      DESCRIPTION
      (
       "Scattered field vector.\n"
       "\n"
       "This variable contains the scattered field vector which \n"
       "is used in the scattering RT calculation. Physically it is the \n"
       "amount of radiation which is scattered into the propagation \n"
       "direction. \n"
       "The vector is calculated by the agenda *sca_vec_agenda* \n"
       "The dimensision of the variable adapts to *stokes_dim*.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\"\n"
       "\n"
       "Usage:      Output of the agenda *sca_int_agenda*. \n"
       "\n"
       "Unit:       W / (m^2 Hz sr) \n"
       "\n"
       "Dimensions: [stokes_dim]"
       ),
      GROUP( Vector_ )));
 
   wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_aa_index" ),
      DESCRIPTION
      (
       "Azimuth angle index for scattering calculations.\n"
       "\n"
       "This variable is used in methods used for computing scattering\n"
       "properties of particles like ext_mat_sptCalc and pha_mat_sptCalc.\n"
       "This holds the information about the azimuth angles for which the \n"
       "scattering calculations are done.  The angles used for computing \n"
       "scattering properties of particle can be different from that used \n"
       "for radiative transfer calculation. \n"
       "\n"
       "Usage:    Input to the methods *ext_mat_sptCalc*, *pha_mat_sptCalc*\n"
       "\n"
       ),
     GROUP( Index_ ))); 

   wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_field" ),
      DESCRIPTION
      (
       "Scattered field field inside the cloudbox.\n"
       "\n"
       "This variable holds the value of the scattering integral.\n"
       "for all points inside the cloudbox. \n"
       "More decription will be written (CE).\n"
       "\n"
       "Usage: Input to *i_fieldUpdate1D*. \n"    
       "\n"
       "Unit: W / (m^2 Hz sr) for each Stokes component.\n"
       "\n"
       "Dimensions: [ p_grid, lat_grid, lon_grid, scat_za_grid, \n"
       "              scat_aa_grid, stokes_dim ]"
       ),
      GROUP( Tensor6_ )));   

 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_f_index" ),
      DESCRIPTION
      (
       "Frequency index for scattering calculations. \n"
       "\n"
       "The calculations inside the cloudbox are only done for one frequency\n"
       "at a time. Some methods used for scattering calculation require the \n"
       "frequency. *f_index* holds the information, for which frequency the \n"
       "scattering calcultations are performed.\n"
       "\n"
       "Usage:      Output of *scat_mono_agenda*.\n"
       ),
      GROUP( Index_ )));


 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_i_lat" ),
      DESCRIPTION
      (
       "Intensity field on cloudbox boundary (equal latitude surfaces).\n"
       "\n"
       "This variable gives the intensity field from all directions defined \n"
       "in *scat_aa_grid* and *scat_za_grid* on each grid point on the two \n"
       "equal \n"
       "latitude surfaces of the boundary of the cloudbox, which is defined \n"
       "by the workspace variable *cloudbox_limits*. It contains all four \n"
       "components of the Stokes vector.\n"
       "\n"
       "This variable is used as interface between the clear sky and the \n"
       "scattering calculations. \n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      In/output to/from *scat_iterateCalc* \n"
       "\n"
       "Unit:        W / (m^2 Hz sr) \n"
       "\n"
       "Dimensions: [ f_grid, p_grid, latitude surface, lon_grid, \n"
       "              scat_za_grid \n  scat_aa_grid, stokes_dim ]"
       ),
      GROUP( Tensor7_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_i_lon" ),
      DESCRIPTION
      (
       "Intensity field on cloudbox boundary (equal longitude surfaces).\n"
       "\n"
       "This variable gives the intensity field from all directions defined \n"
       "in *scat_aa_grid* and *scat_za_grid* on each grid point on the equal\n"
       "latitude surfaces of the boundary of the cloudbox, which is defined \n"
       "by the workspace variable *cloudbox_limits*. It contains all four \n"
       "components of the Stokes vector.\n"
       "\n"
       "This variable is used as interface between the clear sky and the \n"
       "scattering calculations. \n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      In/output to/from *scat_iterateCalc* \n"
       "\n"
       "Unit:        W / (m^2 Hz sr) \n"
       "\n"
       "Dimensions: [ f_grid, p_grid, lat_grid, latitude surface, \n"
       "              scat_za_grid, scat_aa_grid, stokes_dim]"
       ),
      GROUP( Tensor7_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_i_p" ),
      DESCRIPTION
      (
       "Intensity field on cloudbox boundary (equal pressure surfaces).\n"
       "\n"
       "This variable gives the intensity field from all directions defined \n"
       "in *scat_aa_grid* and *scat_za_grid* on each grid point on the equal\n"
       "latitude surfaces of the boundary of the cloudbox, which is defined \n"
       "by the workspace variable *cloudbox_limits*. It contains all four \n"
       "components of the Stokes vector.\n"
       "\n"
       "This variable is used as interface between the clear sky and the \n"
       "scattering calculations. \n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      In/output to/from *scat_iterateCalc* \n"
       "\n"
       "Unit:        W / (m^2 Hz sr) \n"
       "\n"
       "Dimensions: [ f_grid, pressure surfaces, lat_grid, lon_grid, \n" 
       "              scat_za_grid, scat_aa_grid, stokes_dim]"
       ),
      GROUP( Tensor7_ )));


 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_lat_index" ),
      DESCRIPTION
      (
       "Latitude index for scattering calculations.\n"
       "\n"
       "This variable is used in methods used for computing scattering\n"
       "properties of particles like ext_mat_partCalc and pha_matCalc.\n"
       "This holds the information about the position for which the \n"
       "scattering calculations are done. \n"
       "\n"
       "Usage:    Input to the methods *ext_mat_partCalc*, *pha_matCalc*\n"
       "\n"
       ),
     GROUP( Index_ ))); 

 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_lon_index" ),
      DESCRIPTION
      (
       "Longitude index for scattering calculations.\n"
       "\n"
       "This variable is used in methods used for computing scattering\n"
       "properties of particles like ext_mat_partCalc and pha_matCalc.\n"
       "This holds the information about the position for which the \n"
       "scattering calculations are done.  \n"
       "\n"
       "Usage:    Input to the methods *ext_mat_agenda*, *pha_mat_sptCalc*\n"
       "\n"
       ),
     GROUP( Index_ ))); 


 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_p_index" ),
      DESCRIPTION
      (
       "Pressure index for scattering calculations.\n"
       "\n"
       "This variable is used in methods used for computing scattering\n"
       "properties of particles like ext_mat_partCalc and pha_matCalc.\n"
       "This holds the information about the location for which the \n"
       "scattering calculations are done.\n"  
       "\n"
       "Usage:    Input to the methods *ext_mat_partCalc*, *pha_matCalc*\n"
       "\n"
       ),
     GROUP( Index_ ))); 
  wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_rte_agenda" ),
      DESCRIPTION
      (
	"See agendas.cc."
       ),
      GROUP( Agenda_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_za_grid" ),
      DESCRIPTION
      (
       "Zenith angle grid.\n"
       "\n"
       "The zenith angle grid, on which the intensity field and the \n"
       "optical scattering properties are stored. The grid has to be defined\n"
       "if the cloudbox is activated by the flag *cloudbox_on*.\n"
       "The grid must be sorted in decreasing order, with no repetitions.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      Set by the user.\n"
       "\n"
       "Unit:       degrees "
       ),
      GROUP( Vector_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_za_index" ),
      DESCRIPTION
      (
       "Zenith angle index for scattering calculations.\n"
       " \n"
       "This variable is used in methods used for computing scattering \n"
       "properties of particles like ext_mat_sptCalc and pha_mat_sptCalc.\n"
       "This holds the information about the zenith angles for which the \n"
       "scattering calculations are done.  The angles used for computing \n"
       "scattering properties of particle can be different from that used\n"
       "for radiative transfer calculation. \n"
       " \n"
       "Usage:    Input to the methods *ext_mat_sptCalc*, *pha_mat_sptCalc*\n"
       " \n"
       ),
      GROUP( Index_ )));


  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_los" ),
      DESCRIPTION
      (
       "The sensor line-of-sight for each measurement block.\n"
       "\n"
       "Text will be written (PE).\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  [ degrees, degrees ]\n"
       "\n"
       "Size:  [ number of measurement blocks, 1 or 2 ]"
       ),
      GROUP( Matrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_pos" ),
      DESCRIPTION
      (
       "The sensor position for each measurement block.\n"
       "\n"
       "Text will be written (PE).\n"
       "\n"
       "Range for zenith angles is [0,180] beside for 2D where it\n"
       "is [-180,180].\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  [ m, degrees, degrees ]\n"
       "\n"
       "Size:  [ number of measurement blocks, atmosphere_dim ]"
       ),
      GROUP( Matrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "stokes_dim" ),
      DESCRIPTION
      (
       "The dimensionality of the Stokes vector (1-4).\n"
       "\n"
       "Text will be written (Claudia?).\n"
       "\n"
       "Usage:      Set by the user."
       ),
      GROUP( Index_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "stokes_vec" ),
      DESCRIPTION
      (
       "The Stokes vector.\n"
       "\n"
       "The Stokes Vector is a 4 dimensional vector which stores the \n"
       "Stokes components during the RT calculation, for example in the \n"
       "methods *sto_vecGeneral* and *sto_vecScalar*. \n"
       "\n"
       "Usage:      Calculated internally."
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "t_field" ),
      DESCRIPTION
      (
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
       "Unit:       K\n"
       "\n"
       "Dimensions: [ p_grid, lat_grid, lon_grid ]"
       ),
      GROUP( Tensor3_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "t_ground" ),
      DESCRIPTION
      (
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
       "Unit:       K\n"
       "\n"
       "Dimensions: [ lat_grid, lon_grid ]"
       ),
      GROUP( Matrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "vector1" ),
      DESCRIPTION
      (
       "An arbitrary vector.\n"
       "\n"
       "This variable shall be treated as a general variable of type Vector.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage:      Set by user."
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "vector2" ),
      DESCRIPTION
      (
       "An arbitrary vector.\n"
       "\n"
       "This variable shall be treated as a general variable of type Vector.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage:      Set by user."
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "xml_output_type" ),
      DESCRIPTION
      (
       "Flag to determine whether XML output is binary or ascii\n"
       "\n"
       "This flag has to be set using the workspace method SetXMLOutputBinary\n"
       "or SetXMLOutputAscii. One of these methods MUST be called before\n"
       "writing the first output file."
       "\n"
       "Usage:      Set by user."
       ),
      GROUP( Index_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "y" ),
      DESCRIPTION
      (
       "The measurement vector.\n"
       "\n"
       "Text will be written (PE).\n"
       "\n"
       "Usage:      Output from *rteCalc* and methods to apply sensor\n"
       "            characteristics.\n"
       "\n"
       "Unit:       Undefined. Possibilities include: K, W/(m^2 Hz sr) and\n "
       "            optical thickness."
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "z_field" ),
      DESCRIPTION
      (
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
       "set by a method ensuring that this criterium is fulfilled.\n"
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
       "            by a method for applying hydrostatic equilibrium.\n"
       "\n"
       "Unit:       m\n"
       "\n"
       "Dimensions: [ p_grid, lat_grid, lon_grid ]"
       ),
      GROUP( Tensor3_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "z_ground" ),
      DESCRIPTION
      (
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
       "Unit:       m\n"
       "\n"
       "Dimensions: [ lat_grid, lon_grid ]"
       ),
      GROUP( Matrix_ )));





  // Below this line are the definitions from ARTS-1. Shall be removed or
  // updated with better descriptions.
  // =======================================================================
  // 


//   //--------------------< Spectroscopy Stuff >--------------------
//   //                     --------------------
//    wsv_data.push_back
//      (WsvRecord
//       ( NAME( "lines" ),
// 	DESCRIPTION
// 	(
// 	 "A list of spectral line data."
// 	 ), 
// 	GROUP( ArrayOfLineRecord_)));

//    wsv_data.push_back
//      (WsvRecord
//       ( NAME( "lines_per_tg" ),
// 	DESCRIPTION
// 	(
// 	 "A list of spectral line data for each tag.\n"
// 	 "Dimensions: (tag_groups.nelem()) (# of lines for this tag)"
// 	 ), 
// 	GROUP( ArrayOfArrayOfLineRecord_)));

//    wsv_data.push_back
//      (WsvRecord
//       ( NAME( "tgs" ),
// 	DESCRIPTION
// 	(
// 	 "This is an array of arrays of OneTag tag definitions.\n"
// 	 "It defines the available tag groups for the calculation\n"
// 	 "of absorption coefficients and weighting functions.\n"
// 	 "Contrary to the original Bredbeck definition, tags within a\n"
// 	 "group must belong to the same species, because one VMR profile\n"
// 	 "is associated with each tag group."
// 	 ), 
// 	GROUP( TagGroups_)));

//    wsv_data.push_back
//      (WsvRecord
//       ( NAME( "wfs_tgs" ),
// 	DESCRIPTION
// 	(
// 	 "This is an array of arrays of tag group definitions.\n"
// 	 "It defines the tag groups for the calculation of weighting\n"
// 	 "functions. The selected tag groups must be a subgroup of the\n"
// 	 "tag groups defined for the absorption coefficient calculation."
// 	 ), 
// 	GROUP( TagGroups_)));

//    wsv_data.push_back
//      (WsvRecord
//       ( NAME( "lineshape" ),
// 	DESCRIPTION
// 	(
// 	 "Lineshape specification: function, norm, cutoff. There is one entry for\n"
// 	 "each abs_tag, not for each species. This means if you have several\n"
// 	 "abs_tags for different isotopes or transitions of a species, you\n"
// 	 "may use different lineshapes."
// 	 ),
// 	GROUP( ArrayOfLineshapeSpec_)));


//    //--------------------< Continuum Stuff >--------------------
//    //                     -----------------
//    wsv_data.push_back
//      (WsvRecord
//       ( NAME( "cont_description_names" ),
// 	DESCRIPTION
// 	(
// 	 "Continuum / full model absorption tag names. This variable should\n"
// 	 "contain a list of tag names of continuum and full models, respectively.\n"
// 	 "Associated with this WSV is the WSV\n"
// 	 "`cont_description_models' which contains the specific model version of\n"
// 	 "each continuum / full model absorption tag and the WSV\n"
// 	 "`cont_description_parameters' which should contain the continuum / full model\n"
// 	 "user defined parameters. The user defined parameters are only used when\n"
// 	 "the specified model is 'user'. See also the online documentation in\n"
// 	 "arts/doc/doxygen/html/continua_cc.html.\n"
// 	 "\n"
// 	 "The following full water vapor models are implemented:\n"
// 	 "'H2O-MPM87': absorption model (line and continuum) according to \n"
// 	 "   H. J. Liebe,\n" 
// 	 "   A contribution to modeling atmospheric millimeter-wave properties,\n"
// 	 "   Frequenz,  41, 1987, 31-36\n"
// 	 "   and\n"
// 	 "   H. J. Liebe and D. H. Layton,\n"
// 	 "   Millimeter-wave properties of the atmosphere:\n"
// 	 "   Laboratory studies and propagation modeling,\n"
// 	 "   U.S. Dept. of Commerce, National Telecommunications and Information\n"
// 	 "   Administration, Institute for Communication Sciences,\n"
// 	 "   325 Broadway, Boulder, CO 80303-3328, report 87224.\n"
// 	 "'H2O-MPM89': absorption model (line and continuum) according to \n"
// 	 "   H. J. Liebe,\n Int. J. Infrared and Millimeter Waves, 10(6), 1989, 631\n"
// 	 "'H2O-MPM93': absorption model (line and continuum) according to \n"
// 	 "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
// 	 "   Propagation modeling of moist air and suspended water/ice\n"
// 	 "   particles at frequencies below 1000 GHz,\n"
// 	 "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
// 	 "   Propagation Panel,\n Palma de Mallorca, Spain, 1993, May 17-21 \n" 
// 	 "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"
// 	 "'H2O-CP98': absorption model (line and continuum) according to \n"
// 	 "   S. L. Cruz-Pol et al.,\n Radio Science, 33(5), 1319, 1998"
// 	 "   (ece.uprm.edu/~pol/Atmosphere.html)\n"
// 	 "'H2O-PWR98': absorption model (line and continuum) according to \n"
// 	 "   P. W. Rosenkranz,\n "
// 	 "   Radio Science, 33(4),  919, 1998, Radio Science, 34(4), 1025, 1999\n"
// 	 "   (ftp: mesa.mit.edu/phil/lbl_rt).\n"
// 	 "\n"
// 	 "The following full oxygen models are implemented:\n"
// 	 "'O2-MPM93': absorption model (line and continuum) according to\n"
// 	 "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
// 	 "   Propagation modeling of moist air and suspended water/ice\n"
// 	 "   particles at frequencies below 1000 GHz,\n"
// 	 "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
// 	 "   Propagation Panel,\n Palma de Mallorca, Spain, 1993, May 17-21\n"
// 	 "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"
// 	 "'O2-PWR93': absorption model (line and continuum) according to \n"
// 	 "   P. W. Rosenkranz,\n Chapter 2, in M. A. Janssen, \n"
// 	 "   Atmospheric Remote Sensing by Microwave Radiometry\n"
// 	 "   John Wiley & Sons, Inc., 1993 (mesa.mit.edu/phil/lbl_rt)\n"
// 	 "\n"
// 	 "The following continuum parameterizations are implemented:\n"
// 	 "H2O-H2O ('H2O-SelfContStandardType'):\n" 
// 	 "   P. W. Rosenkranz, \n"
// 	 "   Radio Science, Vol. 33, No 4, Pages 919-928, 1998 and \n"
// 	 "   Radio Science, Vol. 34, No 4, Page 1025, 1999 (mesa.mit.edu/phil/lbl_rt)\n"
// 	 "H2O-air ('H2O-ForeignContStandardType'): \n"
// 	 "   P. W. Rosenkranz, \n"
// 	 "   Radio Science, Vol. 33, No 4, Pages 919-928, 1998 and \n"
// 	 "   Radio Science, Vol. 34, No 4, Page 1025, 1999 (mesa.mit.edu/phil/lbl_rt)\n"
// 	 "H2O-air ('H2O-ContMPM93'): \n"
// 	 "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
// 	 "   Propagation modeling of moist air and suspended water/ice\n"
// 	 "   particles at frequencies below 1000 GHz,\n"
// 	 "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
// 	 "   Propagation Panel,\n Palma de Mallorca, Spain, 1993, May 17-21\n"
// 	 "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"      
// 	 "O2-air ('O2-SelfContStandardType'):\n"
// 	 "   P. W. Rosenkranz,\n"
// 	 "   Chapter 2, in M. A. Janssen,\n"
// 	 "   Atmospheric Remote Sensing by Microwave Radiometry,\n"
// 	 "   John Wiley & Sons, Inc., 1993\n"
// 	 "   (mesa.mit.edu/phil/lbl_rt)\n"
// 	 "   and also described in \n"
// 	 "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
// 	 "   Propagation modeling of moist air and suspended water/ice\n"
// 	 "   particles at frequencies below 1000 GHz,\n"
// 	 "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
// 	 "   Propagation Panel,\n Palma de Mallorca, Spain, 1993, May 17-21\n"
// 	 "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"
// 	 "N2-N2 ('N2-SelfContStandardType'):\n"
// 	 "   The functional form of Rosenkranz but with more input parameters.\n"
// 	 "   P. W. Rosenkranz,\n"
// 	 "   Chapter 2, in M. A. Janssen,\n"
// 	 "   Atmospheric Remote Sensing by Microwave Radiometry,\n"
// 	 "   John Wiley & Sons, Inc., 1993 (mesa.mit.edu/phil/lbl_rt)\n"
// 	 "N2-N2 ('N2-SelfContMPM93'):\n"
// 	 "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
// 	 "   Propagation modeling of moist air and suspended water/ice\n"
// 	 "   particles at frequencies below 1000 GHz,\n"
// 	 "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
// 	 "   Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21 \n"
// 	 "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"
// 	 "CO2-CO2 ('CO2-SelfContPWR93'):\n"
// 	 "   P. W. Rosenkranz,\n"
// 	 "   Chapter 2, in M. A. Janssen,\n"
// 	 "   Atmospheric Remote Sensing by Microwave Radiometry,\n"
// 	 "   John Wiley & Sons, Inc., 1993 (mesa.mit.edu/phil/lbl_rt)\n"
// 	 "CO2-N2 ('CO2-ForeignContPWR93'):\n"
// 	 "   P. W. Rosenkranz,\n"
// 	 "   Chapter 2, in M. A. Janssen,\n"
// 	 "   Atmospheric Remote Sensing by Microwave Radiometry,\n"
// 	 "   John Wiley & Sons, Inc., 1993 (mesa.mit.edu/phil/lbl_rt)\n"
// 	 "\n"
// 	 "The following cloud absorption models are implemented:\n"
// 	 "Suspended water droplet ('liquidcloud-MPM93') \n"
// 	 "   absorption parameterization from the MPM93 model:\n"
// 	 "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
// 	 "   Propagation modeling of moist air and suspended water/ice\n"
// 	 "   particles at frequencies below 1000 GHz,\n"
// 	 "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
// 	 "   Propagation Panel,\n Palma de Mallorca, Spain, 1993, May 17-21\n"
// 	 "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"
// 	 "Ice water droplet absorption ('icecloud-MPM93') \n"
// 	 "   parameterization from MPM93 model:\n"
// 	 "   H. J. Liebe and G. A. Hufford and M. G. Cotton,\n"
// 	 "   Propagation modeling of moist air and suspended water/ice\n"
// 	 "   particles at frequencies below 1000 GHz,\n"
// 	 "   AGARD 52nd Specialists Meeting of the Electromagnetic Wave\n"
// 	 "   Propagation Panel,\n Palma de Mallorca, Spain, 1993, May 17-21\n"
// 	 "   (ftp.its.bldrdoc.gov/pub/mpm93/)\n"
// 	 "\n"
// 	 ),
// 	GROUP( ArrayOfString_)));


//    wsv_data.push_back
//      (WsvRecord
//       ( NAME( "cont_description_models" ),
// 	DESCRIPTION
// 	(
// 	 "Continuum / full model absorption model description parameter.\n"
// 	 "See the WSV `cont_description_names' for a detailed description\n"
// 	 "of the allowed continuum models. There should be one string here\n"
// 	 "for each entry in `cont_description_names'.See also the online" 
// 	 "documentation in arts/doc/doxygen/html/continua_cc.html.\n"
// 	 ),
// 	GROUP( ArrayOfString_)));

//    wsv_data.push_back
//      (WsvRecord
//       ( NAME( "cont_description_parameters" ),
// 	DESCRIPTION
// 	(
// 	 "Continuum model parameters. See the WSV `cont_description_names'\n"
// 	 "for a detailed description of the allowed continuum models. There\n"
// 	 "should be one parameter vector here for each entry in\n"
// 	 "`cont_description_names'. See also the online documentation in\n"
// 	 "arts/doc/doxygen/html/continua_cc.html.\n"
// 	 ),
// 	GROUP( ArrayOfVector_)));


//    //--------------------< 1D Input Atmosphere Stuff >--------------------
//    //                     ---------------------------
//    wsv_data.push_back
//      (WsvRecord
//       ( NAME( "raw_ptz" ),
// 	DESCRIPTION
// 	(
// 	 "Matrix has rows:\n"
// 	 "1. Pressure in Pa\n"
// 	 "2. Temperature in K\n"
// 	 "3. Altitude in m"
// 	 ), 
// 	GROUP( Matrix_)));

//    wsv_data.push_back
//      (WsvRecord
//       ( NAME( "raw_vmrs" ),
// 	DESCRIPTION
// 	(
// 	 "The individual VMR profiles. Each species VMR profile comes with a\n"
// 	 "pressure profile. The different species can hence be on different\n"
// 	 "grids.\n"
// 	 "The matrix has rows:\n"
// 	 "1. Pressure in Pa\n"
// 	 "2. VMR profile (absolute number)\n"
// 	 "The array dimension is determined by the number of tag groups."
// 	 ), 
// 	GROUP( ArrayOfMatrix_)));

//    wsv_data.push_back
//      (WsvRecord
//       ( NAME( "vmrs" ),
// 	DESCRIPTION
// 	(
// 	 "The VMRs (unit: absolute number) on the p_abs grid.\n"
// 	 "Dimensions: [tag_groups.nelem(), p_abs.nelem()]"
// 	 ),
// 	GROUP( Matrix_)));


//    wsv_data.push_back
//      (WsvRecord
//        ( NAME( "abs" ),
// 	 DESCRIPTION
// 	 (
// 	  "\n"
// 	  "FIXME: This is the old abs. What should it be now?\n"
// 	  "\n"
// 	  "The matrix of absorption coefficients (in units of [1/m]).\n"
// 	  "Dimensions: [f_mono.nelem(), p_abs.nelem()]"
// 	  ),
// 	 GROUP( Matrix_)));

//    wsv_data.push_back
//      (WsvRecord
//       ( NAME( "abs_per_tg" ),
// 	DESCRIPTION
// 	(
// 	 "These are the absorption coefficients individually for each\n"
// 	 "tag group. The Array contains one matrix for each tag group,\n"
// 	 "the matrix format is the same as that of abs"
// 	 ),
// 	GROUP( ArrayOfMatrix_)));

//    wsv_data.push_back
//      (WsvRecord
//       ( NAME( "xsec_per_tg" ),
// 	DESCRIPTION
// 	(
// 	 "These are the cross sections individually for each tag\n"
// 	 "group. The Array contains one matrix for each tag group,\n"
// 	 "the matrix format is the same as that of abs"
// 	 ),
// 	GROUP( ArrayOfMatrix_)));


 
}
