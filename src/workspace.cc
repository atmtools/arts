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
#include "array.h"
#include "auto_wsv_groups.h"
#include "wsv_aux.h"

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


  //--------------------< Hydrostatic equilibrium >--------------------
  //                     -------------------------
  wsv_data.push_back
    (WsvRecord
     ("hse",
      "This vector holds the parameters for calculating hydrostatic \n"
      "equilibrium (HSE). The length of the vector is either 1 or 5, where\n"
      "the values are: \n "
      "  1: On/off flag. 0 = ignore HSE, 1 = consider HSE.\n " 
      "  2: The pressure of the reference point [Pa]. \n " 
      "  3: The altitude of the reference point [m]. \n " 
      "  4: Gravitational acceleration at the geoid surface [m/s2]. \n " 
      "  5: Number of iterations of the calculations.\n"
      "If the on/off flag is set to 1, the length of the vector must be 5,\n"
      "while if the flag is 0 a length of 1 is OK.\n"
      "See the function hseCalc for some more details.", 
      Vector_));


  //--------------------< RT Stuff >--------------------
  //                     ----------

  wsv_data.push_back
    (WsvRecord
     ("emission",
      "Boolean to include emssion in the calculation of spectra.\n"
      "If this variable is set to 0 (zero) pure transmission calculations \n"
      "will be simulated and, for example, yCalc will give optical \n"
      "thicknesses instead of radiation intensities.",
      Index_));

  wsv_data.push_back
    (WsvRecord
     ("za_pencil",
      "Pencil beam zenith angle grid [deg]. \n"
      "The observation direction is specified by the angle between zenith \n"
      "and the LOS.",
      Vector_));

  wsv_data.push_back
    (WsvRecord
     ("z_tan",
      "Tangent altitude for each spectrum [m].\n"
      "These tangent altitudes include the effect of refraction (if set). \n"
      "In the case of a ground intersection, a geometrical prolongation \n"
      "below the ground is applied to determine the tangent altitude. \n"
      "For upward observations where there are no tangent altitudes, \n" 
      "*z_tan* is set to 999 km. \n"
      "It should be noted that the LOS calculations take *za_pencil* as \n" 
      "input, not *z_tan*. However, *za_pencil* can be calculated from \n"
      "*z_tan* by the function *zaFromZtan*. ",
      Vector_));

  wsv_data.push_back
    (WsvRecord
     ("z_plat",
      "Vertical altitude, above the geiod, of the observation platform [m].",
      Numeric_));

  wsv_data.push_back
    (WsvRecord
     ("l_step",
      "The maximum length, along the LOS, between the points of LOS [m].\n"
      "The final step length will in most cases equal the selected length.\n"
      "There are two rare exceptions:\n"
      "  1. Downward observations from within the atmosphere, where the step\n"
      "     length is adjusted downwards to get an integer number of steps\n"
      "     between the sensor and the tangent or ground point.\n"
      "  2. Limb sounding and the distance from the tangent point to the\n"
      "     atmospheric limit (the highest absorption altitude) is smaller\n"
      "     the selected length. The length is then adjusted to this\n"
      "     distance",
      Numeric_));

  wsv_data.push_back
    (WsvRecord
     ("refr",
      "Boolean for inclusion of refraction (0=no refraction, 1=refraction).",
      Index_));

  wsv_data.push_back
    (WsvRecord
     ("refr_lfac",
      "This factor determines the step length used during the ray tracing \n"
      "performed when considering refraction. \n"
      "The step length applied is *l_step* divided by *refr_lfac*. \n" 
      "Accordingly, this factor gives how many ray tracing steps that are \n"
      "performed for each step of the LOS.",
      Index_));

  wsv_data.push_back
    (WsvRecord
     ("refr_model",
      "A string giving what parameterization to use for the calculation of \n"
      "refractive index. See *refrCalc* for existing choices.",
      String_));

  wsv_data.push_back
    (WsvRecord
     ("refr_index",
      "The refractive index at the pressure levels in p_abs [-].\n",
      Vector_));

  wsv_data.push_back
    (WsvRecord
     ("r_geoid",
      "The local curvature radius of the geoid along the LOS [m].",
      Numeric_));

  wsv_data.push_back
    (WsvRecord
     ("z_ground",
      "The vertical altitude above the geiod of the ground [m].",
      Numeric_));

  wsv_data.push_back
    (WsvRecord
     ("t_ground",
      "The physical temperature of the ground [K].",
      Numeric_));

  wsv_data.push_back
    (WsvRecord
     ("e_ground",
      "The ground emission factor for the frequencies in f_mono [0-1].",
      Vector_));

  wsv_data.push_back
    (WsvRecord
     ("los",
      "Structure to define the line of sight (LOS). See los.h for \n"
      "definition of the structure.", 
      Los_));

  wsv_data.push_back
    (WsvRecord
     ("source",
      "Mean source function between the points of the LOS.",
      ArrayOfMatrix_));

  wsv_data.push_back
    (WsvRecord
     ("trans",
      "The transmissions between the points of the LOS [-].",
      ArrayOfMatrix_));

  wsv_data.push_back
    (WsvRecord
     ("y_space",
      "Radiation entering the atmosphere at the top of the atmosphere, \n"
      "typically cosmic background radiation. This variable is most easily \n"
      "set by the function *y_spaceStd*.",
      Vector_));

  wsv_data.push_back
    (WsvRecord
     ("y",
      "The working set of spectra. \n"
      "The spectra from the different zenith angles are appended to form *y*.",
      Vector_));

  wsv_data.push_back
    (WsvRecord
     ("y0",
      "A reference spectrum. This variable can be used e.g. to save a copy\n"
      "of *y* or to compare the spectra before and after some operation(s).",
      Vector_));

  wsv_data.push_back
    (WsvRecord
     ("h",
      "The H matrix.\n"
      "\n"
      "Can be used to apply the sensor model to monochromatic pencil beam\n"
      "spectra and weighting functions. \n",
      Matrix_));


  //--------------------< WF Stuff >--------------------
  //                     ----------
  wsv_data.push_back
    (WsvRecord
     ("absloswfs",
      "Line of sight weighting functions. \n"
      "See AUG for definition of this quantity. ",
      ArrayOfMatrix_));

  wsv_data.push_back
    (WsvRecord
     ("k_grid",
      "Retrieval grid to be used in calculation of weighting functions (WFs)\n"
      "For example, *k_grid* is the pressure altitude grid for species WFs. \n"
      "Not all WFs need 'k_grid* as input.",
      Vector_));

  wsv_data.push_back
    (WsvRecord
     ("k",
      "The weighting functions (WFs) for a single retrieval/error group.",
      Matrix_));

  wsv_data.push_back
    (WsvRecord
     ("k_names",
      "Names of the retrieval identies in *k*.",
      ArrayOfString_));

  wsv_data.push_back
    (WsvRecord
     ("k_aux",
      "Auxiliary data for *k*. The number of rows of this matrix equals the\n"
      "length of the state vector for the retrieval group (the number of\n"
      "columns of k).\n"
      "The columns hold different quantities:\n"
      "  Col 1: retrieval grid (or correspondingly)\n"
      "  Col 2: a priori values",
      Matrix_));

  wsv_data.push_back
    (WsvRecord
     ("kx",
      "The state weighting function matrix.",
      Matrix_));

  wsv_data.push_back
    (WsvRecord
     ("kx_names",
      "Names of the retrieval identities in *kx*.",
      ArrayOfString_));

  wsv_data.push_back
    (WsvRecord
     ("kx_lengths",
      "The length of the state vector for each retrieval identity in *kx*.",
      ArrayOfIndex_));

  wsv_data.push_back
    (WsvRecord
     ("kx_aux",
      "Auxiliary data for *kx*. As *k_aux* but with the data of the \n"
      "different retrieval groups appended vertically.",
      Matrix_));

  wsv_data.push_back
    (WsvRecord
     ("kb",
      "The model parameters weighting function matrix.",
      Matrix_));

  wsv_data.push_back
    (WsvRecord
     ("kb_names",
      "Names of the model parameter identities in *kb*.",
      ArrayOfString_));

  wsv_data.push_back
    (WsvRecord
     ("kb_lengths",
      "The length of the model vector for each retrieval identity in *kb*.",
      ArrayOfIndex_));

  wsv_data.push_back
    (WsvRecord
     ("kb_aux",
      "Auxiliary data for *kb*. As *k_aux* but with the data of the \n"
      "different forward model groups appended vertically.",
      Matrix_));



  //-------------------< Batch calculation stuff >-----------------------
  //                     -----------------------

  wsv_data.push_back
    (WsvRecord
     ("batchname",
      "Default basename for batch data.",
      String_));

  wsv_data.push_back
    (WsvRecord
     ("ybatch",
      "A batch of spectra.\n"
      "The spectra are stored as columns in a matrix",
      Matrix_));

  wsv_data.push_back
    (WsvRecord
     ("radiosonde_data",
      "An array of Matrix holding data for many radiosonde launches. The\n"
      "dimension of the Array is the number of radiosonde launches. Each\n"
      "Matrix within the Array has dimension nx4, where n is the number of\n"
      "pressure levels. The columns of the Matrix are:\n"
      "\n"
      "pressure [Pa] temperature [K] altitude [m] VMR [1]",
      ArrayOfMatrix_));

  //-------------------< Methods as variables >-----------------------
  //                     --------------------

  wsv_data.push_back
    (WsvRecord
     ("main_agenda",
      "The agenda corresponding to the entire controlfile. This is executed\n"
      "when ARTS is run.",
      Agenda_));

}
