/* Copyright (C) 2000-2012
   Stefan Buehler <sbuehler@ltu.se>
   Patrick Eriksson <patrick.eriksson@chalmers.se>

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
#include "matpackII.h"
#include "matpackIII.h"
#include "matpackVI.h"
#include "array.h"
#include "wsv_aux.h"
#include "ppath.h"
#include "workspace_ng.h"

// Some #defines to make the records better readable:
#define NAME(x)        x 
#define DESCRIPTION(x) x
#define GROUP(x)       x 


void Workspace::define_wsv_data()
{
  
  //--------------------< Build the wsv data >--------------------
  // Initialize to empty, just in case.
  wsv_data.resize(0);

/* Templace record entry:

  wsv_data.push_back
    (WsvRecord
     ( NAME( "workspace_variable_name" ),
       DESCRIPTION
       (
        "Brief description of the variable (1 line).\n"
        "\n"
        "Detailed description of the variable. Don't be too short here,\n"
        "this is the main place where your documentation should be. I\n"
        "really recommend to edit this in a text buffer, so that you can\n"
        "do some re-formatting until it looks nice. Only at the end put it\n"
        "in quotes and add the line breaks.\n"
        "\n"
        "Use blank lines to separate paragraphs.  There really should be a\n"
        "detailed descriptions of all component of your variable, if it\n"
        "has a complicated type. Also some detailed discussion of the\n"
        "dimensions if necessary. Also some detailed discussion of the\n"
        "members if your variable is a structure.\n"
        "\n"
        "Usage:      Set by user (or "Method output.")\n"
        "\n"
        "Units:      E.g., kg/m\n"
        "\n"
        "Dimensions: [ first dimension, second dimension, ... ]\n"
        "or\n"
        "Size:       [ .., nrows, ncols ]\n"
        "\n"
        "Members:    Here you would list the members if your\n"
        "            variable is a structure.\n"
        "\n"
        "Dimensions: [x, y]\n"
        "\n"
        "Unit: Which unit this variable uses\n"
        "\n"
        "Give the keywords above only if they apply, i.e., Members only\n"
        "for a structure, Units only for a physical variable.\n"
        "Use either Dimensions or Size, depending on what is most appropiate\n"
        "for the variable.\n"
        ),
      GROUP( "VariableType" )));

*/


 
  /*----------------------------------------------------------------------
    Let's put in the variables in alphabetical order. This gives a clear
    rule for where to place a new variable and this gives a nicer
    results when the methods are listed by "arts -w all".  No
    distinction is made between uppercase and lowercase letters. The
    sign "_" comes after all letters.
    Patrick Eriksson 2002-05-08
  ----------------------------------------------------------------------*/

  // New name: abs_coef
  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_coef" ),
       DESCRIPTION
       (
        "The matrix of total absorption coefficients.\n"
        "\n"
	"This variable is not used explicitly in a standard calculation, where\n"
	"absorption comes from the lookup table *abs_lookup*. However, it is\n"
	"useful for testing the methods that actually calculate line-by-line\n"
	"absorption, which have this variable as output. These methods are\n"
	"called internally by the method *abs_lookupCreate*, which generates\n"
	"the lookup table.\n"
        "\n"
        "Dimensions: [f_grid, abs_p]\n"
        "\n"
        "Unit: 1/m\n"
        ),
      GROUP( "Matrix" )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "abs_coef_per_species" ),
      DESCRIPTION
      (
       "These are the absorption coefficients individually for each\n"
       "tag group. The Array contains one matrix for each tag group,\n"
       "the matrix format is the same as that of abs_coef\n"
      ),
      GROUP( "ArrayOfMatrix" )));

//   wsv_data.push_back
//    (WsvRecord
//     ( NAME( "abs_coef_agenda" ),
//       DESCRIPTION
//       (
//         "See agendas.cc.\n"
//        ),
//       GROUP( "Agenda" )));
  
  wsv_data.push_back
    (WsvRecord
     (NAME( "abs_cont_models" ),
      DESCRIPTION
      (
       "Continuum / full model absorption model description parameter.\n"
       "See the WSV `abs_cont_names' for a detailed description\n"
       "of the allowed continuum models. There should be one string here\n"
       "for each entry in `abs_cont_names'.See also the online\n" 
       "documentation in arts/doc/doxygen/html/continua_cc.html.\n"
      ),
      GROUP( "ArrayOfString" )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "abs_cont_names" ),
      DESCRIPTION
      (
       "Continuum / full model absorption tag names. This variable should\n"
       "contain a list of tag names of continuum and full models, respectively.\n"
       "Associated with this WSV is the WSV\n"
       "`abs_cont_models' which contains the specific model version of\n"
       "each continuum / full model absorption tag and the WSV\n"
       "`abs_cont_parameters' which should contain the continuum / full model\n"
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
       "\n"
       "The following rain extinction model is implemented:\n"
       "Rain extinction parameterization ('rain-MPM93') from the\n"
       "   MPM93 model, described in:\n" 
       "   H. J. Liebe,\n" 
       "   MPM - An Atmospheric Millimeter-Wave Propagation Model,\n"
       "   Int. J. Infrared and Millimeter Waves, vol. 10(6),\n"
       "   pp. 631-650, 1989;\n"
       "   and based on:\n" 
       "   Olsen, R.L., D.V. Rogers, and D. B. Hodge,\n"
       "   The aR^b relation in the calculation of rain attenuation,\n"
       "   IEEE Trans. Antennas Propagat., vol. AP-26, pp. 318-329, 1978.\n"
       "   IMPORTANT NOTE: rain-MPM93 parameterizes the EXTINCTION by rain,\n"
       "    not just the absorption. Therefore it is not suitable for \n"
       "    calculating thermal emission by rain!\n"
       "    Please use rain-MPM93 only for calculation of attenuation.\n"
      ),
      GROUP( "ArrayOfString" )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "abs_cont_parameters" ),
      DESCRIPTION
      (
       "Continuum model parameters. See the WSV *abs_cont_names*\n"
       "for a detailed description of the allowed continuum models. There\n"
       "should be one parameter vector here for each entry in\n"
       "*abs_cont_names*. See also the online documentation in\n"
       "arts/doc/doxygen/html/continua_cc.html.\n"
      ),
      GROUP( "ArrayOfVector" )));
    
    wsv_data.push_back
    (WsvRecord
     (NAME( "abs_h2o" ),
      DESCRIPTION
      (
       "The total water profile associated with the pressures in *abs_p* [-]\n"
      ),
      GROUP( "Vector" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_lines" ),
       DESCRIPTION
       (
        "A list of spectral line data.\n"
       ), 
       GROUP( "ArrayOfLineRecord" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_lineshape" ),
       DESCRIPTION
       (
        "Lineshape specification: function, norm, cutoff. There is one entry for\n"
        "each abs_tag, not for each species. This means if you have several\n"
        "abs_tags for different isotopes or transitions of a species, you\n"
        "may use different lineshapes.\n"
       ),
       GROUP( "ArrayOfLineshapeSpec" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_lines_per_species" ),
       DESCRIPTION
       (
        "A list of spectral line data for each tag.\n"
        "Dimensions: (tag_groups.nelem()) (# of lines for this tag)\n"
       ), 
       GROUP( "ArrayOfArrayOfLineRecord" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_lookup" ),
       DESCRIPTION
       (
        "An absorption lookup table.\n"
        "\n"
        "This holds an absorption lookup table, as well as all information that\n"
        "is necessary to use the table to extract absorption. Extraction\n"
        "routines are implemented as member functions. \n"
        "\n"
        "This has quite a complicated structure. See Doxygen documentation for\n"
        "class GasAbsLookup for details. See also the Arts User Guide, \n"
        "Chapter 3.1 (\"The gas absorption lookup table\").\n"
        ), 
       GROUP( "GasAbsLookup" )));

    wsv_data.push_back
      (WsvRecord
       ( NAME( "abs_nls" ),
         DESCRIPTION
         (
          "Nonlinear species for absorption lookup table generation.\n"
          "\n"
          "A list of absorption species that should be treated non-linearly.\n"
          "This means that the H2O VMR should be varied when calculating the\n"
          "lookup table for those species.\n"
          "\n"
          "A typical example is for this to containt the Rosenkranz full\n"
          "absorption model species for water vapor and oxygen \n"
          "([\"H2O-PWR98\", \"O2-PWR93\"]).\n"
          "\n"
          "It only makes sense to put a species here if is either a water vapor\n"
          "species, or some other species that uses *abs_h2o*, that is, for which\n"
          "the absorption coefficient depends directly on water vapor.\n"
          "\n"
          "See user guide and online documentation of *abs_pts* and *abs_lookupCreate*\n"
          "for more details and usage examples.\n"
          ), 
         GROUP( "ArrayOfArrayOfSpeciesTag" )));

    wsv_data.push_back
      (WsvRecord
       ( NAME( "abs_nls_pert" ),
         DESCRIPTION
         (
          "Fractional perturbations for the nonlinear species in the absorption\n"
          "lookup table.\n"
          "\n"
          "This is a vector of fractional perturbations that should contain 1\n"
          "(the unperturbed reference profile). A value of 0 may lead to error\n"
          "messages from some absorption routines, so a possible content for this\n"
          "variable is: [1e-24, 1, 2].\n"
          "(This is similar to *abs_t_pert*, but multiplicative, not additive.)\n"
          ), 
         GROUP( "Vector" )));

    wsv_data.push_back
      (WsvRecord
       ( NAME( "abs_nls_interp_order" ),
         DESCRIPTION
         (
          "The interpolation order to use when interpolating absorption between\n"
          "the H2O values given by *abs_nls_pert*. This is used by methods\n"
          "extracting absorption coefficients from the lookup table, and by\n"
          "methods setting up parameters for lookup table generation. Has a\n"
          "default value, which is set in general.arts.\n"
          "\n"
          "Note that the number of points used in the interpolation scheme is\n"
          "interpolation order + 1 (e.g., two for first order interpolation).\n"
          ), 
         GROUP( "Index" )));

    wsv_data.push_back
      (WsvRecord
       ( NAME( "abs_p_interp_order" ),
         DESCRIPTION
         (
          "The interpolation order to use when interpolating absorption\n"
          "between pressure levels. This is used by methods extracting\n"
          "absorption coefficients from the lookup table, and by methods\n"
          "setting up parameters for lookup table generation. Has a\n"
          "default value, which is set in general.arts.\n"
          "\n"
          "Note that the number of points used in the interpolation scheme is\n"
          "interpolation order + 1 (e.g., two for first order interpolation).\n"
          ), 
         GROUP( "Index" )));

    wsv_data.push_back
      (WsvRecord
       ( NAME( "abs_t_pert" ),
         DESCRIPTION
         (
          "Temperature perturbations for the absorption lookup table.\n"
          "\n"
          "This is a vector containing temperature perturbations (in Kelvin) that\n"
          "should be added to the reference temperature profile. (Similar to\n"
          "*abs_nls_pert*, but additive, not multiplicative.) Should normally\n"
          "contain 0, to include the reference profile itself. Example content:\n"
          "[-5, 0, 5].\n"
          ), 
         GROUP( "Vector" )));

    wsv_data.push_back
      (WsvRecord
       ( NAME( "abs_t_interp_order" ),
         DESCRIPTION
         (
          "The interpolation order to use when interpolating absorption between\n"
          "the temperature values given by *abs_t_pert*. This is used by methods\n"
          "extracting absorption coefficients from the lookup table, and by\n"
          "methods setting up parameters for lookup table generation. Has a\n"
          "default value, which is set in general.arts.\n"
          "\n"
          "Note that the number of points used in the interpolation scheme is\n"
          "interpolation order + 1 (e.g., two for first order interpolation).\n"
          ), 
         GROUP( "Index" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_lookup_is_adapted" ),
       DESCRIPTION
       (
        "Flag to indicate whether *abs_lookupAdapt* has already been\n"
        "called.\n"
        "\n"
        "Values: 0=false, 1=true.\n"
        ), 
       GROUP( "Index" )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "abs_n2" ),
      DESCRIPTION
      (
       "The total nitrogen profile associated with the pressures in *abs_p* [-]\n"
      ),
      GROUP( "Vector" )));

  wsv_data.push_back
    (WsvRecord
    ( NAME( "abs_p" ),
      DESCRIPTION
      (
       "List of pressures to be used for the calculation of absorption\n"
       "coefficients. \n"
       "\n"
       "This can be copied from the global *p_grid*, but could also be\n"
       "different. \n"
       "\n"
       "Any absorption method should check that the length of this vector\n"
       "is the same as that of *abs_t*\n"
       "\n"
       "Dimension: [number of pressures]\n"
       "\n"
       "Unit: Pa\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
    (WsvRecord
    ( NAME( "abs_scalar_gas" ),
      DESCRIPTION
      (
       "Scalar gas absorption coefficients.\n"
       "\n"
       "This contains the absorption coefficients for one point in the\n"
       "atmosphere (one set of pressure, temperature, and VMR values). There\n"
       "are two distinct cases:\n"
       "\n"
       "Case a:    For all frequencies and all species:\n"
       "Dimension: [ f_grid, abs_species ]\n"
       "\n"
       "Case b:    For a single frequency for all species:\n"
       "Dimension: [ 1,      abs_species ]\n"
       "\n"
       "Unit: 1/m\n"
       ),
      GROUP( "Matrix" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "abs_scalar_gas_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));
  
  wsv_data.push_back
   (WsvRecord
    ( NAME( "abs_field" ),
      DESCRIPTION
      (
       "Scalar gas absorption field.\n"
       "\n"
       "Contains the scalar gas absorption for all species as a function of\n"
       "*f_grid*, *p_grid*, *lat_grid*, and *lon_grid*. \n"
       "\n"
       "This is mainly for testing and plotting gas absorption. For RT\n"
       "calculations, gas absorption is calculated or extracted locally,\n"
       "therefore there is no need to store a global field. But this variable\n"
       "is handy for easy plotting of absorption vs. pressure, for example.\n"
       "\n"
       "Unit:       1/m\n"
       "\n"
       "Dimensions: [species, f_grid, p_grid, lat_grid, lon_grid]\n"
        ),
      GROUP( "Tensor5" ))); 

  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_species" ),
       DESCRIPTION
       (
        "Tag groups for scalar gas absorption.\n"
        "\n"
        "This is an array of arrays of SpeciesTag tag definitions. It defines the\n"
        "available tag groups for the calculation of scalar gas absorption\n"
        "coefficients.  See online documentation of method *SpeciesSet* for\n"
        "more detailed information how tag groups work and some examples.\n"
        ), 
       GROUP( "ArrayOfArrayOfSpeciesTag" )));

  wsv_data.push_back
    (WsvRecord
    ( NAME( "abs_t" ),
      DESCRIPTION
      (
       "List of temperatures to be used for the calculation of absorption\n"
       "coefficients.\n"
       "\n"
       "In contrast to the global *t_field*, this is just a vector. Any\n"
       "absorption method should check that the length of this vector is the\n"
       "same as that of *abs_p*\n"
       "\n"
       "Dimension: [number of pressures]\n"
       "\n"
       "Unit: K\n"
       ),
      GROUP( "Vector" )));

 wsv_data.push_back
    (WsvRecord
    ( NAME( "abs_vec" ),
      DESCRIPTION
      (
       "Total absorption vector.\n"
       "\n"
       "This variable contains the absorption coefficient vector which \n"
       "is used in the RTE calculation. It is \n"
       "the physical absorption which includes particle absorption \n"
       "for all chosen particle types as well as gaseous absorption for\n"
       "all chosen gaseous species.\n" 
       "The vector is calculated by the agendas *opt_prop_gas_agenda* \n"
       "and, if scattering calculations are performed, \n"
       "*opt_prop_part_agenda* \n"
       "The dimensision of the variable adapts to *stokes_dim*.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      Output of the agendas *opt_prop_gas_agenda* \n"
       "                             and *opt_prop_part_agenda* \n"
       "\n"
       "Unit:        [Hz, m^2]\n"
       "\n"
       "Dimensions: [f_grid, stokes_dim]\n"
        ),
       GROUP( "Matrix" )));

 wsv_data.push_back
    (WsvRecord
     ( NAME("abs_vec_spt"),
       DESCRIPTION
       (
        "Absorption vector for a single particle type.\n"
        "\n"
        "This variable contains the elements of absorption vector of a \n"
        "single particle, given  It is calculated in the agenda \n"
        "*spt_calc_agenda*.\n"
        "\n"
        "ARTS user guide (AUG) gives the formulas used for computing all \n"
        "the elements of absorption vector.\n"
        "\n"
        "Usage:      Input and Output of the method abs_vec_sptCalc\n"
        "\n"
        "Unit:        m^2\n"
        "\n"
        "Dimensions: [N_particletypes,stokes_dim]\n"
        ),
       GROUP( "Matrix" ) ));

  wsv_data.push_back
    (WsvRecord
     (NAME( "abs_vmrs" ),
      DESCRIPTION
      (
       "The VMRs (unit: absolute number) on the abs_p grid.\n"
       "Dimensions: [tag_groups.nelem(), abs_p.nelem()]\n"
      ),
      GROUP( "Matrix" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_xsec_per_species" ),
       DESCRIPTION
       (
        "Absorption cross sections.\n"
        "\n"
        "This variable contains absorption cross section xsec individually for\n"
        "each tag group. The Array contains one matrix for each tag group, the\n"
        "matrix format is the same as that of abs_coef.\n"
        "\n"
        "Dimensions: [abs_species](f_grid, abs_p)\n"
        "\n"
        "Unit:       m^2 (alpha = xsec * n * VMR),\n"
        "            where n is total density.\n"
        ),
       GROUP( "ArrayOfMatrix" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "antenna_dim" ),
      DESCRIPTION
      (
       "The dimensionality of the antenna pattern (1-2).\n"
       "\n"
       "A dimensionality of 1 means that only the respons variation in the\n"
       "zenith direction is considered. The provided respons shall then be the\n"
       "integrated in the azimuth direction. For 2D, the respons of the\n"
       "antenna has both a zenith and azimuth variation.\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  Integer value [1-2].\n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "antenna_los" ),
      DESCRIPTION
      (
       "The relative line-of-sight of each antenna pattern.\n"
       "\n"
       "This variable describes the line-of-sight of the individual antennae\n"
       "relative to *sensor_los*. If only one antenna is present the matrix\n"
       "should contain a row of zero(s). The number of columns corresponds to\n"
       "the *antenna_dim*, with the first column containing zenith angles\n"
       "and the second azimuth angles. If each measurement block corresponds\n"
       "to a single antenna pattern, the normal choice is to set the angle(s)\n"
       "of this variable to zero.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  [ degrees, degrees ]\n"
       "\n"
       "Size:  [ number of antennae, 1 or 2 ]\n"
       ),
      GROUP( "Matrix" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "antenna_response" ),
      DESCRIPTION
      (
       "The antenna pattern/response.\n"
       "\n"
       "This WSV describes the antenna response as a function of polarisation\n"
       "(pol), frequencue (f), zenith angle (za) and azimuth angle (aa).\n"
       "\n"
       "Polarisation dimension: If this dimension has size 1, the data are\n"
       "applied for all polarisations of concern. The data are otherwise used\n"
       "in sequential order. This signifies that, in general, the first\n"
       "polarisation \"layer\" corresponds to the first stokes dimension\n"
       "etc. An exception is if a polarisation rotation has been applied. In\n"
       "any case, it is up to the user to ensure that polarisations are\n"
       "consistently defined.\n"
       "\n"
       "Frequency dimension: If this dimension has size 1, the data are\n"
       "applied for all frequencies of concern. The given frequency must be\n"
       "inside the frequency range of concern. A linear interpolation is\n"
       "otherwise applied.\n"
       "\n"
       "Zenith angle dimension: This dimension must always have a size >= 2\n"
       "The response outside covered grid range is treated as zero. If\n"
       "*antenna_dim* equals 1, the data should correspond to the response\n"
       "integrated in the azimuthal direction.\n"
       "\n"
       "Azimuth angle dimension: If *antenna_dim* equals 1, this dimension\n"
       "must have size 1. A size >= 2 is otherwise required. The response\n"
       "outside covered grid range is treated as zero.\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Dimensions: \n"
       "   GriddedField4:\n"
       "      ArrayOfString field_names[N_pol]\n"
       "      Vector f_grid[N_f]\n"
       "      Vector za_grid[N_za]\n"
       "      Vector aa_grid[N_aa]\n"
       "      Tensor4 data[N_pol][N_f][N_za][N_aa]\n"
       ),
      GROUP( "GriddedField4" )));

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
       "changed. However, not all methods are working for higher dimensions.\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit: Integer value.\n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "atm_fields_compact" ),
      DESCRIPTION
      (
       "A compact combination of all atmospheric fields for a clear-sky\n"
       "calculation on a common set of grids.\n"
       "\n"
       "This concerns temperature, altitude, and gas VMRs.\n"
       "\n"
       "The data are stored in a *GriddedField4*.\n"
       "\n"
       "The order of the fields must be:\n"
       "T[K] z[m] VMR_1[1] ... VMR_n[1]\n"
       "(the field names for the gases do not have to start with VMR_)\n"
       "\n"
       "Usage: Used inside batch calculations, to hold successive atmospheric\n"
       "       states from an *ArrayOfGriddedField4*.\n"
       "\n"
       "Possible future extensions: Add a similar variable\n"
       "particle_fields_compact for hydrometeors?\n"
       "\n"
       "Dimensions: \n"
       "   GriddedField4:\n"
       "      ArrayOfString field_names[N_fields]\n"
       "      Vector p_grid[N_p]\n"
       "      Vector lat_grid[N_lat]\n"
       "      Vector lon_grid[N_lon]\n"
       "      Tensor4 data[N_fields][N_p][N_lat][N_lon]\n"
       ),
      GROUP( "GriddedField4" )));
   
    wsv_data.push_back
   (WsvRecord
    ( NAME( "atm_fields_compact_all" ),
      DESCRIPTION
      (
       "A compact combination of all atmospheric fields for a clear-sky\n"
       "and cloud particle scattering calculation on a common set of grids.\n"
       "\n"
       "This concerns temperature, altitude, scattering particles and gas VMRs.\n"
       "\n"
       "The data are stored in a *GriddedField4*.\n"
       "\n"
       "The order of the fields must be:\n"
       "T[K] z[m] LWC[kg/m3] IWC[kg/m3] Rain[kg/m2/s] Snow[kg/m2/s] VMR_1[1] ... VMR_n[1]\n"
       "(the field names for the gases do not have to start with VMR_)\n"
       "\n"
       "Usage: Used inside batch calculations, to hold successive atmospheric\n"
       "       states from an *ArrayOfGriddedField4*.\n"
       "\n"
       "Dimensions: \n"
       "   GriddedField4:\n"
       "      ArrayOfString field_names[N_fields]\n"
       "      Vector p_grid[N_p]\n"
       "      Vector lat_grid[N_lat]\n"
       "      Vector lon_grid[N_lon]\n"
       "      Tensor4 data[N_fields][N_p][N_lat][N_lon]\n"
       ),
      GROUP( "GriddedField4" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "backend_channel_response" ),
      DESCRIPTION
      (
       "The response of each backend channel.\n"
       "\n"
       "The response is given as an *ArrayOfGriddedField1*. The grid consists of\n"
       "relative frequencies. These relative frequencies are added to \n"
       "*f_backend* to obtain the absolute frequency for each response value.\n"
       "The actual data are the response at each frequency grid point.\n"
       "\n"
       "There are here two options. If the array has length 1, the same\n"
       "response is applied for all channels. Accordingly, this assumes that\n"
       "all channels have the same response function. The second option is to\n"
       "specify the response for each channel seperately. This signifies that\n"
       "the *backend_channel_response* array has either 1 or n elements, where\n"
       "n is the length of *f_backend*\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Size:  Array[N_ch]\n"
       "       GriddedField1 \n "
       "       [N_f] \n"
       "       [N_f] \n"
       ),
      GROUP( "ArrayOfGriddedField1" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "backend_channel_response_multi" ),
       DESCRIPTION
       (
        "As *backend_channel_response* but describes an instrument with\n"
        "muliple mixer/reciever chains.\n"
        "\n"
        "See *f_backend_multi* for when to use this variable and size\n"
        "constraints.\n"
        "\n"
        "Usage: Set by the user.\n "
        ),
      GROUP( "ArrayOfArrayOfGriddedField1" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "basics_checked" ),
      DESCRIPTION
      (
       "OK-flag for basic variables and grids.\n"
       "\n"
       "For example, this variable flags that the (clear-sky) atmosphere\n"
       "and the surface are defined in a formally correct way, e.g. that\n"
       "the sizes of atmospheric fields match the atmospheric grids.\n"
       "\n"
       "Shall be set by *basics_checkedCalc*. See that WSM for treated WSVs.\n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "batch_atm_fields_compact" ),
      DESCRIPTION
      (
       "An array of compact atmospheric states.\n"
       "\n"
       "This is used to hold a set of *atm_fields_compact* for batch\n"
       "calculations. \n"
       ),
      GROUP( "ArrayOfGriddedField4" )));
   
  wsv_data.push_back
   (WsvRecord
    ( NAME( "batch_atm_fields_compact_all" ),
      DESCRIPTION
      (
       "An array of compact atmospheric states, including scattering particles.\n"
       "\n"
       "This is used to hold a set of *atm_fields_compact_all* for batch\n"
       "calculations. \n"
       ),
      GROUP( "ArrayOfGriddedField4" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "batch_cloudbox_limits" ),
      DESCRIPTION
      (
       "An array of cloudbox_limits.\n"
       "\n"
       "This is used to hold a set of *cloudbox-limits* for batch\n"
       "calculations. \n"
       ),
      GROUP( "ArrayOfArrayOfIndex" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "batch_pnd_fields" ),
      DESCRIPTION
      (
       "An array of compact pnd states.\n"
       "\n"
       "This is used to hold a set of 1D *pnd_field* for batch\n"
       "calculations. \n"
       ),
      GROUP( "ArrayOfTensor4" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "blackbody_radiation" ),
      DESCRIPTION
      (
       "The blackbody radiation for a given temperature.\n"
       "\n"
       "This variable can be seen as the source term for thermal emission.\n"
       "Normally, this variable should match the Planck function. The version\n"
       "of the Planck function taking frequency as input is considered as\n"
       "default for ARTS. The unit for radiance is then W / [m2 Hz sr].\n"
       "For frequencies where the Rayleigh-Jeans approximation is valid\n"
       "(but not recommended) option is to set this variable to be equal\n"
       "to the physical temperature, resulting in K as unit.\n"
       "\n"
       "Inside some methods, such as DOIT, the calculation of this source\n"
       "term can be hard-coded.\n"
       "\n"
       "Usage:      Set by *blackbody_radiation_agenda*.\n"
       "\n"
       "Unit:       See above.\n"
       "\n"
       "Dimensions: [ f_grid ]\n"
       ),
      GROUP( "Vector" )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "blackbody_radiation_agenda" ),
       DESCRIPTION
       (
        "See agendas.cc.\n"
        ),
       GROUP( "Agenda" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "cloudbox_checked" ),
      DESCRIPTION
      (
       "OK-flag for variables associated with the cloudbox.\n"
       "\n"
       "This variable flags that cloudbox variables are defined in a formally\n"
       "and practically correct way. For example, that there is sufficient\n"
       "space between the cloudbox and edges of the model atmosphere (for\n"
       "2D and 3D). Pure clear-sky variables are covered by *basics_checked*.\n"
       "\n"
       "Relevant checks are performed by *cloudbox_checkedCalc.\n"
       ),
      GROUP( "Index" )));

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
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  Boolean.\n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "cloudbox_limits" ),
      DESCRIPTION
      (
       "The limits of the cloud box.\n"
       "\n"
       "This variable defines the extension of the cloud box. The cloud box \n"
       "is defined to be rectangular in the used coordinate system, with \n"
       "limits exactly at points of the involved grids. This means, for \n"
       "example, that the vertical limits of the cloud box are two pressure \n"
       "levels. For 2D, the angular extension of the cloud box is between \n"
       "two points of the latitude grid, and likewise for 3D but then also \n"
       "with a longitude extension between two grid points. The latitude and\n"
       "longitude limits for the cloud box cannot be placed at the end \n"
       "points of the corresponding grid as it must be possible to calculate\n"
       "the incoming intensity field.\n"
       "\n"
       "The variable *cloudbox_limits* is an array of index value with\n"
       "length twice *atmosphere_dim*. For each dimension there is a lower \n"
       "limit and an upper limit. The order of the dimensions is as usual \n"
       "pressure, latitude and longitude. The upper limit index must be \n"
       "greater then the lower limit index. For example, \n"
       "*cloudbox_limits* = [0 5 4 11 4 11] means that cloud box extends\n"
       "between pressure levels 0 and 5, and latitude and longitude points 4\n"
       "and 11.\n"
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
       "Size:  [ 2 * atmosphere_dim ]\n"
       ),
      GROUP( "ArrayOfIndex" )));

 wsv_data.push_back
    (WsvRecord
     ( NAME("complex_n"),
       DESCRIPTION
       (
        "Complex refractive index (n).\n"
        "\n"
        "This matrix describes the dielectric properties of a medium. \n"
        "A typical usage of this variable is to describe the properties of\n"
        "the surface (if it is assumed to be flat).\n"
        "\n"
        "This is a two-column matrix. The first column holds the real part\n"
        "of n, and the second column the imaginary part. The number of rows\n"
        "should match *f_grid*. The case of one row is also accepted,\n"
        "interpreted as that the refractive index is constant with \n"
        "frequency.\n"
        "\n"
        "Unit:       -\n"
        "\n"
        "Dimensions: [f_grid or 1, 2]\n"
        ),
       GROUP( "Matrix" ) ));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "diy_dx" ),
      DESCRIPTION
      (
       "Derivative of *iy* with respect to retrieval quantities.\n"
       "\n"
       "The variable gives the derivative if *iy* with respect to some\n"
       "variables (but not all jacobian variables). Handled are only variables\n"
       "affecting monochromatic pencil beam radiances where an (semi-)\n"
       "analytical expression can be applied (and that this calculation way\n"
       "has been selected when the jacobian has been set-up).\n"
       "\n"
       "Usage:      Output of *iy_clearsky_agenda*.\n"
       "\n"
       "Dimensions: \n"
       "     [n_quantities][ n_retrieval_points, f_grid, stokes_dim ]\n"
       ),
      GROUP( "ArrayOfTensor3" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "dispersion_do" ),
      DESCRIPTION
      (
       "Flag to consider dispersion.\n"
       "\n"
       "Inside ARTS, the general assumption is that all frequency components\n"
       "share a common propagation path, i.e. no dispersion. If a WSM has\n"
       "this variable as input, it is possible to override this assumption.\n"
       "In that case a propagation path is calculated for each frequency\n"
       "component. A correct treatment of dispersion requires further that\n"
       "*refr_index_agenda* contains methods returning a frequency dependent\n"
       "refractive index.\n"
       "\n"
       "Usage:   Set by the user..\n"
       ),
      GROUP( "Index" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_conv_flag" ),
      DESCRIPTION
      (
       "Flag for the convergence test.\n"
       "\n"
       "This variable is initialized with 0 inside the method \n"
       "*doit_i_fieldIterate*.\n"
       "If after an iteration the convergence test is fulfilled, 1 is \n"
       "assigned which means that the iteration is completed. \n"
       "\n"
       "Usage: Method output. \n"
      ), 
      GROUP( "Index" ))); 

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_conv_test_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_i_field" ), 
      DESCRIPTION
      (
       "Radiation field.\n" 
       "\n"
       "This variable is used to store the monochromatic radiation field \n"
       "inside the cloudbox which is found by an iterative solution (DOIT).\n"
       "Refer to AUG for further information.\n"
       "\n"
       "Usage: Method output. \n"    
       "\n"
       "Unit: W / (m^2 Hz sr) for each Stokes component.\n"
       "\n"
       "Size: [(cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
       "       (cloudbox_limits[3] - cloudbox_limits[2]) +1, \n"
       "       (cloudbox_limits[5] - cloudbox_limits[4]) +1, \n"
       "        N_za, N_aa, N_i ]\n"
       ),
       GROUP( "Tensor6" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_i_field1D_spectrum" ),
      DESCRIPTION
      (
       "Radiation field for the whole frequency spectrum. \n"
       "\n"
       "This variable holds the radiation field. In contrast to \n"
       "*doit_i_field* this variable has an additional freqeuncy \n"
       "dimension. This variable is only used for 1D DOIT \n"
       "calculations.\n"
       "\n"
       "Usage: Output of *DoitCloudboxFieldPut*\n"
       "\n"
       "Unit: W / (m^2 Hz sr)\n"
       "\n"
        "Size: [N_f \n"
       "       (cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
       "        N_za, N_aa, N_i ]\n"
       ),
      GROUP( "Tensor4" )));
 
 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_i_field_old" ),
      DESCRIPTION
      (
       "Intensity field inside the cloudbox.\n"
       "\n"
       "This variable is used to store the intensity field inside the\n"
       "cloudbox while performing the iteration. One has to store the\n"
       "intensity field of the previous iteration to be able to do the \n"
       "convergence test after each iteration.\n"
       "Refer to AUG for more information.\n"
       "\n"
       "Usage: Method output. \n"    
       "\n"
       "Unit: W / (m^2 Hz sr) for each Stokes component.\n"
       "\n"
       "Size: [(cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
       "       (cloudbox_limits[3] - cloudbox_limits[2]) +1, \n"
       "       (cloudbox_limits[5] - cloudbox_limits[4]) +1, \n"
       "        N_za, N_aa, N_i ]\n"
       ),
      GROUP( "Tensor6" )));
 
 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_is_initialized" ),
      DESCRIPTION
      (
       "Flag to determine if *DoitInit* was called.\n"
       "\n"
       "This flag is checked by *ScatteringDoit* to make sure that\n"
       "*DoitInit* was called before.\n"
       ),
      GROUP( "Index" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_iteration_counter" ),
      DESCRIPTION
      (
       "Counter for number of iterations.\n"
       "\n"
       "This variable holds the number of iterations \n"
       "while solving the VRTE using the DOIT method. \n"
       ),
      GROUP( "Index" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_mono_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));
 
 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_rte_agenda" ),
      DESCRIPTION
      (
       "See agendas.cc.\n"
       ),
      GROUP( "Agenda" ))); 

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_scat_field_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" ))); 

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_scat_field" ),
      DESCRIPTION
      (
       "Scattered field field inside the cloudbox.\n"
       "\n"
       "This variable holds the value of the scattering integral for all\n"
       "points inside the cloudbox. For more information refer to AUG.\n"
       "\n"
       "Usage: Input to *doit_i_fieldUpdate...*. \n"    
       "\n"
       "Unit: W / (m^2 Hz sr) for each Stokes component.\n"
       "\n"
       "Size: [(cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
       "       (cloudbox_limits[3] - cloudbox_limits[2]) +1, \n"
       "       (cloudbox_limits[5] - cloudbox_limits[4]) +1, \n"
       "        N_za, N_aa, N_i ]\n"
       ),
      GROUP( "Tensor6" )));   

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_za_grid_opt" ),
      DESCRIPTION
       (
        "Optimized zenith angle grid.\n"
        "\n"
        "Output of the method *doit_za_grid_optCalc*.\n"
        "\n"
        "Usage:   Output of *doit_za_grid_optCalc*   \n"
        "\n"
        "Unit:    degrees \n"
        ),
      GROUP( "Vector" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_za_grid_size" ),
      DESCRIPTION
      (
       "Number of equidistant grid points of the zenith angle grid, \n"
       "defined from 0 to 180 deg, for the scattering integral calculation. \n"
       "\n"
       "Usage: Output of *DoitAngularGridsSet*.\n"
       ),
      GROUP( "Index" )));
 
 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_za_interp" ),
      DESCRIPTION
      (
       "Flag for interplation method in zenith angle dimension.\n"
       "\n"
       "0 - linear interpolation \n"
       "1 - cubic interpolation \n"
       "\n"
       "Usage: Set by user in *doit_za_interpSet*. \n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "edensity_field" ),
      DESCRIPTION
      (
       "The field of free electrons.\n"
       "\n"
       "This variable gives the density of free electrons at each crossing of\n"
       "the pressure, latitude and longitude grids. The density for a point\n"
       "between the grid crossings is obtained by multi-linear interpolation\n"
       "of the field.\n"
       "\n"
       "Can be set to be empty, which is interpreted as zero wind speed\n"
       "everywhere.\n"
       "\n"
       "Unit:       1/m3\n"
       "\n"
       "Dimensions: [ p_grid, lat_grid, lon_grid ]  or [ 0 0 0 ].\n"
       ),
      GROUP( "Tensor3" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "ext_mat" ),
       DESCRIPTION
      (
       "Total extinction matrix.\n"
       "\n"
       "This variable contains the extinction coefficient matrix which\n"
       "is used in the RT calculation in the cloudbox . It is the physical\n"
       "extinction matrix which includes particles extinction for all chosen\n"
       "particle types and gaseous extinction for all chosen gaseous species.\n" 
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      Output of the agendas *opt_prop_gas_agenda* \n"
       "                             and *opt_prop_part_agenda* \n" 
       "\n"
       "Unit:       [Hz, m^2, m^2] "
       "\n"
       "Dimensions: [f_grid, stokes_dim, stokes_dim]\n"
       ),
       GROUP( "Tensor3" )));

  wsv_data.push_back
     (WsvRecord
    ( NAME( "ext_mat_spt" ),
      DESCRIPTION
      (
       "Extinction matrix for a single particle type.\n"
       "\n"
       "This variable contains the elements for extinction matrix of a  \n"
       "single particle for a given propagation direction. It is calculated\n"
       "input as well as the output of the agenda *spt_calc_agenda*.  \n"
       "\n"
       "Usage:      Output of *spt_calc_agenda* \n"
       "\n"
       "Unit:        m^2 \n"
       "\n"
       "Dimensions: [N_particletypes, stokes_dim, stokes_dim]\n"
       ),
      GROUP( "Tensor3" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "file_index" ),
       DESCRIPTION
       (
        "Index number for files.\n"
        "\n"
        "See *WriteXMLIndexed* for further information.\n"
        "\n"
        "Usage:   Input to *WriteXMLIndexed* and *ReadXMLIndexed*. \n"
        ),
        GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "forloop_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));
  
  wsv_data.push_back
    (WsvRecord
     ( NAME( "forloop_index" ),
       DESCRIPTION
       (
        "The index for for-loops.\n"
        "\n"
        "This is the index that is used by method *ForLoop* to loop over\n"
        "*forloop_agenda*. \n"
        ),
        GROUP( "Index" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "f_backend" ),
       DESCRIPTION
       (
        "The frequency position of each backend (spectrometer) channel.\n"
        "\n"
        "Usage: Set by the user.\n "
        "\n"
        "Unit:  Hz\n"
        ),
        GROUP( "Vector" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "f_backend_multi" ),
       DESCRIPTION
       (
        "As *f_backend* but describes an instrument with muliple\n"
        "mixer/reciever chains.\n"
        "\n"
        "This variable is needed when e.g. the reciever has several mixers\n"
        "or the the reciever measures several polarisation and the channels\n"
        "differ in position or response function. \n"
        "\n"
        "The array has one element for each \"reciever chain\". The array\n"
        "length must match *backend_channel_response_multi*, and possibly\n"
        "also *lo_multi*.\n"
        "\n"
        "Usage: Set by the user.\n "
        "\n"
        "Unit:  Hz\n"
        ),
        GROUP( "ArrayOfVector" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "f_grid" ),
       DESCRIPTION
       (
        "The frequency grid for monochromatic pencil beam calculations.\n"
        "\n"
        "Usage: Set by the user.\n "
        "\n"
        "Unit:  Hz\n"
        ),
        GROUP( "Vector" )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "f_index" ),
      DESCRIPTION
      (
       "Frequency index. \n"
       "\n"
       "Not all methods handle all monochromatic frequencies (of *f_grid*) in\n"
       "parellel and this variable is used for communication between methods,\n"
       "holding the index of the frequency treated presently.\n"
       "\n"
       "In some contexts, a negative f_index means all frequencies.\n"
       "\n"
       "Usage: Method output.\n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "fos_angles" ),
      DESCRIPTION
      (
       "Angles and weights for calculation of the scattering integral.\n"
       "\n"
       "The approach for calculating the scattering integral is described\n"
       "in *iyFOS*. This variable gives the angles and weights to apply.\n"
       "\n"
       "This is a matrix with three columns. The first column holds the\n"
       "zenith angles. These are normally in the range [0,180], but values\n"
       "down to -180 are accepted for 2D calculations. Azimuth angles are\n"
       "found in the second column, that shall be in the range [-180,180].\n"
       "The third columns holds the wieght to apply for each direction. The\n"
       "sum of these weights shall theoretically be 4*pi, the solid angle\n"
       "of the integration sphere. However, for increased flexibility this\n"
       "sum is allowed to deviate with about 50% from that value.\n" ),
      GROUP( "Matrix" )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "fos_i" ),
      DESCRIPTION
      (
       "Scattering order handled presently.\n"
       "\n"
       "This variable must be set to 0 by the user. This variable is later\n"
       "used internally by the FOS methods to handle recursive calls, to\n" 
       "know when *fos_n* has been reached.\n" ),
      GROUP( "Index" )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "fos_n" ),
      DESCRIPTION
      (
       "The scattering order to apply.\n"
       "\n"
       "This variable determines the maximum scattering order that will be\n"
       "considered in the FOS scheme. For example, 1 corresponds to that only\n"
       "single scattering is considered. The value 0 is accepted and results\n"
       "in \"pure clearsky\" calculations.\n" ),
      GROUP( "Index" )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "fos_y" ),
      DESCRIPTION
      (
       "Incoming radiation, for estimating the scattering integral.\n"
       "\n"
       "Output of *fos_y_agenda*. See further that agenda.\n"
       "\n"
       "Dimensions: [fos_angles, f_grid, stokes_dim]\n" ),
      GROUP( "Tensor3" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "fos_y_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));
  
//  wsv_data.push_back
//    (WsvRecord
//     ( NAME( "geomag_los_calc_agenda" ),
//       DESCRIPTION
//       (
//         "See agendas.cc.\n"
//        ),
//       GROUP( "Agenda" )));

//   wsv_data.push_back
//     (WsvRecord
//     ( NAME( "geomag_los" ),
//       DESCRIPTION
//       (
//        "Magnetic field along the line of sight\n"
//        "\n"
//        "more text by Nikolay \n"
//        "\n"
//        "Unit: ..."
//        "\n"
//        "Dimensions: [Magnetic field B, angle between B and los] \n"
//        ),
//       GROUP( "Matrix" )));

   
  wsv_data.push_back
    (WsvRecord
     (NAME( "g0" ),
      DESCRIPTION
      (
       "Gravity at zero altitude.\n"
       "\n"
       "This variable is \"little g\" at the reference ellipsiod. That is,\n"
       "for Earth this is a value around 9.81 m/s2\n"
       ),
      GROUP( "Numeric" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "g0_agenda" ),
       DESCRIPTION
       (
        "See agendas.cc.\n"
        ),
       GROUP( "Agenda" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "iy" ),
      DESCRIPTION
      (
       "Monochromatic pencil beam radiance spectrum.\n"
       "\n"
       "This variable holds a single spectrum, with values corresponding\n"
       "to infinite frequency and spatial resolution (compare to *y*).\n"
       "\n"
       "The variable is used to represent spectra at all positions of the\n"
       "propagation path and can e.g. temporarily hold radiation entering\n"
       "the atmpophere from space. The unit depends on if emission is \n"
       "considered or not (no conversion to e.g. brightness temperature shall\n"
       "be applied).\n"
       "\n"
       "Usage:      Used by radiative transfer methods.\n"
       "\n"
       "Unit:       W / (m^2 Hz sr) or transmission.\n"
       "\n"
       "Dimensions: [ f_grid, stokes_dim ]\n"
       ),
      GROUP( "Matrix" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "iyb" ),
      DESCRIPTION
      (
       "Monochromatic pencil beam data for one measurement block.\n"
       "\n"
       "The data for all *iy* of a measurement block appended to a vector,\n"
       "following the sorting order used for *y*.\n"
       "\n"
       "Usage:      Used internally.\n"
       "\n"
       "Unit:       W / (m^2 Hz sr) or transmission.\n"
       "\n"
       "Dimensions: [ naa*nza*nf*stokes_dim ] where naa is length of\n"
       "            mblock_aa_grid, za length of mblock_za_grid and nf is\n"
       "            length of f_grid.\n"
       ),
      GROUP( "Vector" )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "iy_agenda_call1" ),
       DESCRIPTION
       (
        "Flag to handle recursive calls of *iy_clearsky_agenda*\n"
        "\n"
        "The agenda *iy_clearsky_agenda* can be used recursively and this flag\n"
        "is used to tell the methods inside the agenda which is the primary\n"
        " call. This is handled automatically for methods using\n"
        "*iy_clearsky_agenda*, such as *yCalc*, but the user must set this\n"
        "variable to 1 if the agenda is called directly inside the control\n"
        "file (which should be a rare case).\n"
        ),
       GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "iy_aux" ),
      DESCRIPTION
      (
       "Data auxiliary to *iy*.\n"
       "\n"
       "This variable makes it possible to provide auxiliary information for\n"
       "each value in *iy*. The standard usage is to provide transmission,\n"
       "then to complement radiance data. See the relevant WSMs for details.\n"
       "\n"
       "The data are not used for any calculations. The data are just passed\n"
       "on to create *y_aux*. A possibility here is to set *iy_aux* to e.g.\n"
       "-1000 inside *iy_cloudbox_agenda*, to obtain a way to catch\n"
       "interceptions with the cloudbox.\n"
       "\n"
       "Usage:      Used by radiative transfer methods.\n"
       "\n"
       "Dimensions: [ f_grid, stokes_dim ]\n"
       ),
      GROUP( "Matrix" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "iy_aux_variable" ),
      DESCRIPTION
      (
       "Selection of variable for *iy_aux*.\n"
       "\n"
       "Some radiative transfer methods can return several variables. This\n"
       "string is used for selecting what shall be the auxiliary output, ie.\n"
       "what quantity to put in *iy_aux*.\n"
       "\n"
       "This WSV is not used by all WSMs returning *iy_aux*. And allowed\n"
       "options vary between the WSM that consider the WSV.\n"
       "\n"
       "Usage:   Set by the user.\n"
       ),
      GROUP( "String" )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "iy_clearsky_agenda" ),
       DESCRIPTION
       (
        "See agendas.cc.\n"
        ),
       GROUP( "Agenda" )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "iy_clearsky_basic_agenda" ),
       DESCRIPTION
       (
        "See agendas.cc.\n"
        ),
       GROUP( "Agenda" )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "iy_cloudbox_agenda" ),
       DESCRIPTION
       (
        "See agendas.cc.\n"
        ),
       GROUP( "Agenda" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "iy_error" ),
      DESCRIPTION
      (
       "Estimation of calculation errors in *iy*.\n"
       "\n"
       "As *y_error*, but treats *iy* and can be left empty if\n"
       "*iy_error_type* is 0.\n"
       "\n"
       "Usage:      Used by radiative transfer methods.\n"
       "\n"
       "Unit:       W / (m^2 Hz sr) or transmission.\n"
       "\n"
       "Dimensions: [ f_grid, stokes_dim ]\n"
       ),
      GROUP( "Matrix" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "iy_error_type" ),
      DESCRIPTION
      (
       "Characteristics of error values in *iy_error*.\n"
       "\n"
       "These options are defined:\n"
       "   0: The error is zero. *iy_error* can then be left undefined or\n"
       "      empty.\n"
       "   1: The error values are totally uncorrelated.\n"
       "   2: The error values are totally correlated.\n"
       "The distinction between case 1 and 2 is important, as the weighting\n"
       "with sensor response data must be performed differently for the two\n"
       "cases.\n"
       ),
      GROUP( "Index" )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "iy_space_agenda" ),
       DESCRIPTION
       (
        "See agendas.cc.\n"
        ),
       GROUP( "Agenda" )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "iy_surface_agenda" ),
       DESCRIPTION
       (
        "See agendas.cc.\n"
        ),
       GROUP( "Agenda" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "iy_transmission" ),
      DESCRIPTION
      (
       "Transmission to be included in *iy*.\n"
       "\n"
       "The calculation of *iy* can be performed over several propation path\n"
       "branches, and there can be recursive calls of *iy_clearsky_agenda*.\n"
       "This variable gives the transmission from the end point of the present\n"
       "branch and the sensor for such recursive cases.\n"
       "\n"
       "This variable is used purely internally. The exact usage can vary\n"
       "between different RT integration schemes.\n"
       "\n"
       "Usage:      Internally inside iy_clearsky_agenda.\n"
       "\n"
       "Unit:       1\n"
       "\n"
       "Dimensions: [ f_grid, stokes_dim, stokes_dim ]\n"
       ),
      GROUP( "Tensor3" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "iy_variable" ),
      DESCRIPTION
      (
       "Selection of variable for *iy*.\n"
       "\n"
       "Some radiative transfer methods can return several variables. This\n"
       "string is used for selecting what shall be the main output, ie. what\n"
       "quantity to put in *iy*.\n"
       "\n"
       "This WSV is not used by all WSMs returning *iy*. And allowed\n"
       "options vary between the WSM that consider the WSV.\n"
       "\n"
       "Usage:   Set by the user.\n"
       ),
      GROUP( "String" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "jacobian" ),
      DESCRIPTION
      (
       "The Jacobian matrix.\n"
       "\n"
       "The matrix holding the Jacobians of the retrieval quantities. Each\n"
       "quantity, and its subdivision into atmospheric grids, are stored as\n"
       "columns in the matrix. The matrix has to be initialised before the\n"
       "quantities can be defined. Initialisation WSM is *jacobianInit*.\n"
       "Retrieval quantities are then added with *jacobianAdd...* methods.\n"
       "See the online help. Pure numerical calculation is described by\n"
       "*jacobian_calc_agenda* and are performed by *jacobianCalc*.\n"
       "\n"
       "Units:   See the different retrieval quantities.\n"
       "\n"
       "Dimension: [ y, number of retrieval quantities and grids ]\n"
      ),
      GROUP( "Matrix" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "jacobian_do" ),
      DESCRIPTION
      (
       "Flag to activate jacobian calculations.\n"
       "\n"
       "If this variable is set to 0, no jacobian calculations will be\n"
       "even if such calculations have been set-up (through the jacobianAddXxx\n"
       "methods).\n"
      ),
      GROUP( "Index" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "jacobian_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "jacobian_indices" ),
      DESCRIPTION
      (
       "First and last column index in *jacobian* for each retrieval quantity."
       "\n"
       "This variable tells which part of *jacobian* that corresponds to \n"
       "each jacobian quantity.\n"
       "\n"
       "Usage:      Set by *jacobianClose*.\n"
      ),
      GROUP( "ArrayOfArrayOfIndex" )));

//  Keep this (PE 081113)
//
//  wsv_data.push_back
//    (WsvRecord
//     ( NAME( "jacobian_particle_update_agenda" ),
//       DESCRIPTION
//       (
//         "See agendas.cc.\n"
//        ),
//       GROUP( "Agenda" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "jacobian_quantities" ),
      DESCRIPTION
      (
       "The retrieval quantities in the Jacobian matrix.\n"
       "\n"
       "An array of retrieval quantities for which the jacobians are\n"
       "calculated.\n"
       "\n"
       "Usage: Quantities are added by the jacobianAdd WSMs.\n"
      ),
      GROUP( "ArrayOfRetrievalQuantity" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "jacobian_y_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "lat" ),
      DESCRIPTION
      (
       "A latitude.\n"
       "\n"
       "Unit:  degrees\n"
       ),
      GROUP( "Numeric" )));

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
       "Geocentric latitudes are used.\n"
       "\n"
       "For 1D calculations this vector shall be set to be empty.\n"
       "\n"
       "For 2D cases the latitudes shall be interpreted as the angular\n"
       "distance inside the orbit plane from the equator (values\n"
       "outside +-90 deg are allowed).\n"
       "\n"
       "For 3D, the valid latitude range is [-90,90].\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  degrees\n"
       ),
      GROUP( "Vector" )));
 
 wsv_data.push_back
   (WsvRecord
    ( NAME( "lat_true" ),
      DESCRIPTION
      (
       "Latitudinal geolocation for 1D and 2D data.\n"
       "\n"
       "The variables *lat_grid* and *lon_grid* contain true positions only\n"
       "for 3D. For 1D and 2D, the geographical position is given by\n"
       "*lat_true* and *lon_true*. Can be left empty when not used.\n"
       "Otherwise:\n"
       "\n"
       "   1D: *lat_true* shall have length 1\n"
       "\n"
       "   2D: Both *lat_true* and *lon_true* shall have a length matching\n"
       "       *lat_grid*. That is, *lat_true* and *lon_true* shall not be\n"
       "       seen as grids, they are vectors giving the actual lat or lon\n"
       "       for each point corresponding to *lat_grid*.\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  degrees\n"
       ),
      GROUP( "Vector" )));
 
  wsv_data.push_back
   (WsvRecord
    ( NAME( "lo" ),
      DESCRIPTION
      (
       "The local oscillator frequency.\n"
       "\n"
       "A local oscillator frequency is used in a heterodyne system when\n"
       "the mixer folds the spectra from from radio frequencies (RF) to\n"
       "intermediate frequencies (IF).\n"
       "\n"
       "Unit:  Hz\n"
       "\n"
       "Usage: Set by the user.\n"
       ),
      GROUP( "Numeric" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "lo_multi" ),
      DESCRIPTION
      (
       "Local oscillator frequencies.\n"
       "\n"
       "As *lo* but describes an instrument with multiple mixers. A vector\n"
       "element for each LO. The size of this variable and\n"
       "*sideband_response_multi* shall match, and probably also\n"
       "*sideband_mode_multi*.\n"
       "\n"
       "Unit:  Hz\n"
       "\n"
       "Usage: Set by the user.\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "lon" ),
      DESCRIPTION
      (
       "A longitude.\n"
       "\n"
       "Unit:  degrees\n"
       ),
      GROUP( "Numeric" )));

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
       "For 1D and 2D, this WSV shall be set to be empty.\n"
       "\n"
       "Allowed values for longitudes is the range [-360,360]. The difference\n"
       "between last and first value can not exceed 360 degrees. A difference\n"
       "of exactly 360 deg. means that the complete globe is covered and no\n"
       "propagation paths will reach a longitude edge.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  degrees\n"
       ),
      GROUP( "Vector" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "lon_true" ),
      DESCRIPTION
      (
       "Longitudinal geolocation for 1D and 2D data.\n"
       "\n"
       "The variables *lat_grid* and *lon_grid* contain true positions only\n"
       "for 3D. For 1D and 2D, the geographical position is given by\n"
       "*lat_true* and *lon_true*. Can be left empty when not used.\n"
       "Otherwise:\n"
       "\n"
       "   1D: *lon_true* shall have length 1\n"
       "\n"
       "   2D: Both *lat_true* and *lon_true* shall have a length matching\n"
       "       *lat_grid*. That is, *lat_true* and *lon_true* shall not be\n"
       "       seen as grids, they are vectors giving the actual lat or lon\n"
       "       for each point corresponding to *lat_grid*.\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  degrees\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "main_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));
   
   wsv_data.push_back
    (WsvRecord
     (NAME( "massdensity_field" ),
      DESCRIPTION
      (
       "The field of atmospheric scattering particle types in unit of massdensity\n"
       "(like IWC/LWC/Rain/Aerosol etc.).\n"
       "\n"
       "If *massdensity_field* is obtained from Chevallier91L data, the scat. particles are\n"
       "hydrometeors of type: CLW[kg/m^3] CIW[kg/m^3] Rain[kg/(m2*s)] Snow[kg/(m2*s)]\n"
       "\n"
       "NOTE: naming discussion: *massdensity_field* should be renamed in the future\n"
       "to a more general term (not limited to massdensities).\n"
       "\n"
       "Possible future extension: In the future a *massdensity_field_raw* might be needed,\n"
       "which contains the not-interpolated mass concentrations of scattering particles.\n"
       "This is not needed in the moment, since *massdensity_field* is only used with \n"
       "batch profile data (where all atmospheric variables are on the same grid).\n"
       "\n"
       "Usage:\tmassdensity data is used to calculate pnd_fields\n"
       "\n"
       "Unit:\tdepending on what is read, usually [kg/m3]\n"
       "\n"
       "Dimension:\t[ type, p_grid, lat_grid, lon_grid ]\n"
      
       ),
      GROUP( "Tensor4" ))); 
    
/*    wsv_data.push_back
    (WsvRecord
     (NAME( "massdensity_threshold" ),
      DESCRIPTION
      (
       "Threshold value for minimum massdensity concentration in *massdensity_field*.\n"
       "\n"
       "Check WSM *Massdensity_cleanup* for more information!\n"
       "\n"
       "Default value:\t1e-15\n"      
       ),
      GROUP( "Numeric" )));   
*/     

  wsv_data.push_back
   (WsvRecord
    ( NAME( "mblock_aa_grid" ),
      DESCRIPTION
      (
       "The azimuthal angle grid for each measurement block.\n"
       "\n"
       "This variable should normally contain the azimuth grid of the\n"
       "antenna pattern. The grid is given as an angular off-set with\n"
       "respect to the angles in *sensor_los*.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  degrees\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "mblock_index" ),
      DESCRIPTION
      (
       "Measurement block index. \n"
       "\n"
       "Used to tell agendas the index of present measurement block.\n"
       "\n"
       "Usage: Used internally.\n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "mblock_za_grid" ),
      DESCRIPTION
      (
       "The zenith angle grid for each measurement block.\n"
       "\n"
       "This variable should normally contain the zenith grid of the\n"
       "antenna pattern. The grid is given as an angular off-set with\n"
       "respect to the angles in *sensor_los*.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  degrees\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_antenna" ),
       DESCRIPTION
       (
        "Antenna pattern description for dedicated MC calculaions.\n"
        "\n"
        "Usage: Input to MCGeneral. Set by *mc_antennaSetGaussian* and similar\n"
        "       methods.\n"
        ), 
       GROUP( "MCAntenna" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_cloud_opt_path" ),
       DESCRIPTION
       (
        "The cloud optical path integrated over the field of view.\n"
        "\n"
        "Usage: Output from mc_IWP_cloud_opt_pathCalc \n"
        ), 
       GROUP( "Numeric" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_cloud_opt_path_error" ),
       DESCRIPTION
       (
        "Standad error in the cloud optical path integrated over the field\n"
        "of view.\n"
        "\n"
        "Usage: Output from mc_IWP_cloud_opt_pathCalc \n"
        ), 
       GROUP( "Numeric" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_error" ),
       DESCRIPTION
       (
        "Error in simulated *y* when using a Monte Carlo approach.\n"
        "\n"
        "Usage: Output from Monte Carlo functions. \n"
        "\n"
        "Units: Depends on *y_unit*.\n"
        "\n"
        "Size:  [ stokes_dim ]\n"
        ), 
       GROUP( "Vector" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_iteration_count" ),
       DESCRIPTION
       (
        "Counts the number of iterations (or photons) used in the MC\n "
        "scattering algorithm.\n"
        "\n"
        "Usage: Set by MCGeneral and other MC methods.\n"
       ),
       GROUP( "Index" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_IWP" ),
       DESCRIPTION
       (
        "The ice water path integrated over the field of view\n"
        "\n"
        "Usage: Output from mc_IWP_cloud_opt_pathCalc \n"
        ), 
       GROUP( "Numeric" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_IWP_error" ),
       DESCRIPTION
       (
        "The standard error of ice water path integrated over the field of view\n"
        "\n"
        "Usage: Output from mc_IWP_cloud_opt_pathCalc \n"
        ), 
       GROUP( "Numeric" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_points" ),
       DESCRIPTION
       (
        "Counts the number of MC endpoints in each grid cell.\n"
        "\n"
        "Usage: Set by MCGeneral and other MC methods.\n"
        ),
       GROUP( "Tensor3" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_seed" ),
       DESCRIPTION
       (
        "The integer seed for the random number generator used by\n"
        "Monte Carlo methods.\n"
        "\n"
        "Usage: Set by MCSetSeed.\n"
        ),
       GROUP( "Index" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_std_err" ),
       DESCRIPTION
       (
        "Target precision (1 std. dev.) for Monte Carlo calculations.\n"
        "\n"
        "Usage: Set by the user.\n"
        ),
       GROUP( "Numeric" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_max_time" ),
       DESCRIPTION
       (
        "The maximum time allowed for Monte Carlo calculations.\n"
        "\n"
        "Usage: Set by the user.\n"
        "\n"
        "Unit: s\n"
        ),
       GROUP( "Index" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_max_iter" ),
       DESCRIPTION
       (
        "The maximum number of iterations allowed for Monte Carlo\n"
        "calculations.\n"
        "\n"
        "Usage: Set by the user.\n"
        ),
       GROUP( "Index" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_z_field_is_1D" ),
       DESCRIPTION
       (
        "Flag for MCGeneral and associated methods.\n"
        "\n"
        "If fields outside the cloud box are 1D, this flag can be set to 1\n"
        "and the calculations will be more rapid.\n"
        "\n"
        "Usage: Set by the user.\n"
        ),
       GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "met_amsu_data" ),
      DESCRIPTION
      (
       "The AMSU data set.\n"
       "\n"
       "This is intended as input for the method ybatchMetProfiles. It holds the\n"
       "latitude, longitude, satellite zenith angle and amsu-b corrected and \n"
       "uncorrected brightness temperatures.  It also has information about \n"
       "the particular pixel corresponds to a land or sea point.  This will be \n"
       "read in the method ybatchMetProfiles and the profiles corresponding to \n"
       "each latitude and longitude will be read in.\n"
       "\n"
       "See documentation of WSM *ybatchMetProfiles* for more information.\n"
       ),
      GROUP( "Matrix" )));
 
 wsv_data.push_back
   (WsvRecord
    ( NAME( "nelem" ),
      DESCRIPTION
      (
        "This variable is used by the VectorSet workspace method.\n"
       ),
      GROUP( "Index" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "ncols" ),
      DESCRIPTION
      (
        "This variable is used by the MatrixSet, Tensor3Set, etc. \n"
        "workspace methods.\n"
       ),
      GROUP( "Index" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "nrows" ),
      DESCRIPTION
      (
        "See *ncols*.\n"
       ),
      GROUP( "Index" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "npages" ),
      DESCRIPTION
      (
        "See *ncols*.\n"
       ),
      GROUP( "Index" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "nbooks" ),
      DESCRIPTION
      (
        "See *ncols*.\n"
       ),
      GROUP( "Index" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "nshelves" ),
      DESCRIPTION
      (
        "See *ncols*.\n"
       ),
      GROUP( "Index" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "nvitrines" ),
      DESCRIPTION
      (
        "See *ncols*.\n"
       ),
      GROUP( "Index" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "nlibraries" ),
      DESCRIPTION
      (
        "See *ncols*.\n"
       ),
      GROUP( "Index" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "met_profile_calc_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "molarmass_dry_air" ),
      DESCRIPTION
      (
       "The average molar mass of dry air.\n"
       "\n"
       "This could also be referred to as the average molecular weight for\n"
       "dry air. The definition of \"dry air\" can differ between planets and\n"
       "methods using the WSV. For Earth, this should be a value around\n"
       "28.97.\n"
       ),
      GROUP( "Numeric" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "opt_prop_gas_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "opt_prop_part_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "output_file_format" ),
      DESCRIPTION
      (
       "Output file format. \n"
       "\n"
       "This variable sets the format for output files. It could be set to\n"
       "\"ascii\" for plain xml files, \"zascii\" for zipped xml files, or\n"
       "\"binary\".\n"
       "\n"
       "To change the value of this variable use the workspace methods\n"
       "*output_file_formatSetAscii*, *output_file_formatSetZippedAscii*, and\n"
       "*output_file_formatSetBinary*\n"
       ),
      GROUP( "String" )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "particle_masses" ),
      DESCRIPTION
      (
       "The mass of each particle type stored in a vector \n"
       "\n"
       "Usage: Set by the user.\n"
       ),
      GROUP( "Vector" )));
    
   wsv_data.push_back
   (WsvRecord
    ( NAME( "part_species" ),
      DESCRIPTION
      (
       "Array of Strings containing user defined input for particle number\n"
       "density calculations performed by *pnd_fieldSetup*.\n"
       "\n"
       "Each String has to have the following structure with elements separated\n"
       "by dashes:\n"
       "\n"
       "- profile type [*String*]\n"
       "\t the type of particle bulk profile (mass content, precip rate, or\n"
       "\t similar) to act on. Free form, but needs to match (name, order)\n"
       "\t *atm_fields_compact* field names.\n"
       "- particle size distribution [*String*]:\n"
       "\t the size distribution function/parameterisation to apply. For currently\n"
       "\t possible PSDs see *pnd_fieldSetup*.\n"
       "- particle type [*String*]:\n"
       "\t the type (material/phase) of the individual particles to select from\n"
       "\t *scat_data_raw* (Ice, Water, or similar). Free form, but will select\n"
       "\t particles with matching *scat_data_meta*.type.\n"
       "- sizemin and sizemax [*Numeric*]:\n"
       "\t the minimum and maximum size (volume equivalent sphere radius in um) of\n"
       "\t the individual particles to consider. Wildcards ('*') select all\n"
       "\t particles on lower and upper size ends, respectively. When left out,\n"
       "\t whole size range is considered.\n"
       "\n"
       "Example: [''IWC-MH97-Ice-2-1000'', ''LWC-HM98_STCO-Water-0.1-10'', ...]\n"
       ),
      GROUP( "ArrayOfString" )));

    wsv_data.push_back
   (WsvRecord
    ( NAME( "pha_mat" ),
      DESCRIPTION
      (
       "Ensemble averaged phase matrix.\n"
       "\n"
       "This workspace variable represents the actual physical phase\n"
       "matrix (averaged over all particle types) for given propagation\n"
       "directions. It is calculated in the method *pha_matCalc*\n"
       "\n"
       "ARTS user guide (AUG) gives the formula used for computing this\n"
       "variable. Use the index to find where this variable is discussed.\n"
       "The variable is listed as a subentry to \"workspace variables\".\n"
       "\n"
       "Usage:      Output of the method *pha_matCalc*\n"
       "\n"
       "Unit:        m^2\n"
       "\n"
       "Dimensions: [ scat_za_grid, scat_aa_grid, stokes_dim, stokes_dim ]\n"
       ),
      GROUP( "Tensor4" )));
   
   wsv_data.push_back
   (WsvRecord
    ( NAME( "pha_mat_spt" ),
      DESCRIPTION
      (
       "Phase matrix for a single particle type.\n"
       "\n"
       "This variable contains the elements of phase matrix for a single \n"
       "particle for given propagation direction. \n"
       "It is the calculated in the agenda *pha_mat_spt_agenda*. \n"
       "The elements of the phase matrix are calculated from   \n"
       "the single scattering data. "
       "ARTS user guide (AUG) gives the formulas used for computing all \n"
       "elements of the phase matrix for a given particle type.\n"
// Commented out by Gerrit 2011-04-12, this is incorrect and I cannot
// correct it.
// 
//       "\n"
//       "Usage:      Input and Output of the method pha_mat_sptCalc\n"
       "\n"
       "Unit:        m^2\n"
       "\n"
       "Dimensions: [N_particletypes, *scat_za_grid*, *scat_aa_grid*, *stokes_dim*, *stokes_dim*]\n"
       ),
      GROUP( "Tensor5" )));

    wsv_data.push_back
   (WsvRecord
    ( NAME( "pha_mat_spt_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" ))); 

   wsv_data.push_back
   (WsvRecord
    ( NAME( "pha_mat_sptDOITOpt" ),
      DESCRIPTION
      (
       "Interpolated phase matrix.\n"
       "\n"
       "This variable contains the data of the phase matrix in the \n"
       "scattering frame interpolated on the actual frequency (the variable\n"
       "is used inside *doit_mono_agenda*) and also interpolated on all \n"
       "possible scattering angles following from all combinations of \n"
       "*scat_za_grid* and *scat_aa_grid*. \n"
       "\n"
       "Usage:      Input of the method *pha_mat_sptFromDataDOITOpt*\n"
       "\n"
       "Unit:        m^2\n"
       "\n"
       "Dimensions: \n"
       "[particle types]\n"
       "[T, scat_za_grid,scat_aa_grid, scat_za_grid, scat_aa_grid,\n"
       "stokes_dim, stokes_dim]\n"
       ),
      GROUP( "ArrayOfTensor7" )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "pnd_field" ),
      DESCRIPTION
      (
       "Particle number density field.\n"
       "\n"
       "This veriable corresponds to the particle number density fields \n"
       "for all particle types being read in the WSMs *ParticleTypeAdd* \n"
       "or *ParticleTypeAddAll*. \n"
       "Another alternative is to create *pnd_field* with *pnd_fieldSetup*.\n"
       "\n"
       "Usage:      Calculated internally.\n"
       "\n"
       "Unit:        m^-3\n"
       "\n"
       "Size: [N_particletypes, \n"
       "       (cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
       "       (cloudbox_limits[3] - cloudbox_limits[2]) +1, \n"
       "       (cloudbox_limits[5] - cloudbox_limits[4]) +1 ] \n"
        ),
      GROUP( "Tensor4" )));

//    wsv_data.push_back
//    (WsvRecord
//     ( NAME( "pnd_field_perturb" ),
//       DESCRIPTION
//       (
//        "The field representing particle number density perturbations.\n"
//        "\n"
//        "This variable gives the perturbation of particle number density\n"
//        "of the chosen particle types as a function of p_grid, lat_grid,\n"
//        "lon_grid for each retrieval quantity. The variable has to be\n"
//        "prepared outside ARTS and it has to be setup prior to calling\n"
//        "*jacobianAddParticle*. Since it is added to *pnd_field* during\n"
//        "the calculation of the Jacobian, the perturbations are absolute\n"
//        "and as such should have the same unit as *pnd_field*\n"
//        "\n"
//        "See further the ARTS user guide (AUG). Use the index to find where\n"
//        "this variable is discussed. The variable is listed as a subentry to\n"
//        "\"workspace variables\".\n"
//        "\n"
//        "Usage:      Set by the user.\n"
//        "\n"
//        "Unit:       m^-3\n"
//        "\n"
//        "Dimensions: [N_retrieval_quantities, as *pnd_field* ]\n"
//         ),
//       GROUP( "Tensor5" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "pnd_field_raw" ),
      DESCRIPTION
      (
       "The particle number density field data.\n"
       "\n"
       "This variable contains the particle number density data for all \n"
       "chosen particle types. It includes the grids corresponding to the \n"
       "grids in the database. \n"
       "*pnd_field_raw* is an Array of GriddedField3. It includes a\n"
       "GriddedField3 for each particle type which contains the data and \n"
       "also the grids.\n"
       "\n"
       "Usage: Used in the methods *ParticleTypeAdd* and \n"
       "       *ParticleTypeAddAll*\n"
       "\n"
       "Unit:  m^-3\n"
       "\n"
       "Size:  Array[N_pt]\n"
       "       GriddedField3 \n "
       "       [N_p] \n"
       "       [N_lat] \n"
       "       [N_lon] \n"
       "       [N_p, N_lat, N_lon] \n"
       ),
      GROUP( "ArrayOfGriddedField3" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath" ),
      DESCRIPTION
      (
       "The propagation path for one line-of-sight.\n"
       "\n"
       "This variable describes the total (pencil beam) propagation path for\n"
       "a given combination of starting point and line-of-sight. The path is\n"
       "described by a data structure of type Ppath. This structure contains\n"
       "also additional fields to faciliate the calculation of spectra and\n"
       "interpolation of the atmospheric fields.\n"
       "\n"
       "The data struture is too extensive to be described here, but it is\n"
       "described carefully in the ARTS user guide (AUG). Use the index to\n"
       "find where the data structure, Ppath, for propagation paths is \n"
       "discussed. It is listed as a subentry to \"data structures\".\n"
       "\n"
       "Usage: Output from *ppath_agenda*.\n"
       ),
      GROUP( "Ppath" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath_inside_cloudbox_do" ),
      DESCRIPTION
      (
       "Flag to perform ray tracing inside the cloudbox.\n"
       "\n"
       "Standard propagation path calculations stop at the bounday of the\n"
       "clouxbox, or stop directly if started inside the cloudbox. This WSV\n"
       "allows scatting methods to obtain propagation paths inside the\n"
       "cloudbox. Hence, this variable is for internal usage primarily.\n"
       "\n"
       "Usage: For communication between modules of arts.\n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath_lmax" ),
      DESCRIPTION
      (
       "Maximum length between points describing propagation paths.\n"
       "\n"
       "See *ppath_stepGeometric* for a description of this variable.\n"
       "\n"
       "Usage: Ppath methods such as *ppath_stepGeometric*.\n"
       ),
      GROUP( "Numeric" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath_lraytrace" ),
      DESCRIPTION
      (
       "Maximum length of ray tracing steps when determining propagation\n"
       "paths.\n"
       "\n"
       "See *ppath_stepRefractionEuler* for a description of this variable.\n"
       "\n"
       "Usage: Refraction ppath methods such as *ppath_stepRefractionEuler*.\n"
       ),
      GROUP( "Numeric" )));

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
       "See *ppath_step_agenda* for more information on this variable and\n"
       "the calculation of propagation paths. Or read the chapter on\n"
       "propagation paths in the ARTS user guide.\n"
       "\n"
       "Usage:   In/output to/from *ppath_step_agenda*.\n"
       "\n"
       "Members: See AUG.\n"
       ),
      GROUP( "Ppath" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath_step_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "p_grid" ),
      DESCRIPTION
      (
       "The pressure grid.\n"
       "\n"
       "The pressure levels on which the atmospheric fields are defined.\n"
       "This variable must always be defined. The grid must be sorted in\n"
       "decreasing order, with no repetitions.\n"
       "\n"
       "No gap between the lowermost pressure level and the surface is \n"
       "allowed. The uppermost pressure level defines the practical upper\n"
       "limit of the atmosphere as vacuum is assumed above.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  Pa\n"
       ),
      GROUP( "Vector" )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "p_hse" ),
      DESCRIPTION
      (
       "Reference pressure calculation of hydrostatic equilibrium.\n"
       "\n"
       "The altitude specified by this pressure is used as the reference\n"
       "when calculating hydrostatic equilibrium. That is, the geometrical\n"
       "altitude at this pressure is not changed.\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  Pa\n"
       ),
      GROUP( "Numeric" )));

  wsv_data.push_back
    (WsvRecord
    ( NAME( "refr_index" ),
      DESCRIPTION
      (
       "Real part of the refractive index of air.\n"
       "\n"
       "The variable contains the refractive index summed over all relevant\n"
       "constituents, at one position in the atmosphere. This refractive\n"
       "is related to the phase velocity. See also *refr_index_group*.\n"
       "\n"
       "Unit: 1\n"
       ),
      GROUP( "Numeric" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "refr_index_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));

  wsv_data.push_back
    (WsvRecord
    ( NAME( "refr_index_group" ),
      DESCRIPTION
      (
       "Group index of refractivity.\n"
       "\n"
       "This variable is defined as the ratio between group velocity and the\n"
       "speed of ligh in vacuum. That is, it is defined as the \"standard\"\n"
       "refractive index, but refers to the group velocity instead of the\n"
       "phase velocity. See also *refr_index*.\n"
       "\n"
       "Unit: 1\n"
       ),
      GROUP( "Numeric" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "refellipsoid" ),
      DESCRIPTION
      (
       "Reference ellipsoid.\n"
       "\n"
       "This vector specifies the shape of the reference ellipsoid. The\n"
       "vector must have length 2, where the two elements are:\n"
       "  1: Equatorial radius.\n"
       "  2: The eccentricity.\n"
       "The eccentricity is sqrt(1-b*b/a*a) where a and b are equatorial and\n"
       "polar radius, respectively. If the eccentricity is set to 0, an\n"
       "average radius should be used instead of the equatorial one.\n"
       "\n"
       "The eccentricity must be 0 for 1D calculations, as a spherical Earth\n"
       "is implied by setting *atmosphere_dim* to 1. For 2D, the selected\n"
       "ellipsoid parameters should be selected according to cross-section\n"
       "between the real ellipsoid and the 2D plane considered. That is\n"
       "the applied ellipsoid shall have een converted to match the internal\n"
       "treatment of 2D cases. For 3D, models can be used, such as WGS84.\n"
       "\n"
       "Usage:  Set by the user.\n"
       "\n"
       "Size:   [ 2 ]\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
  (WsvRecord
   ( NAME( "rte_doppler" ),
    DESCRIPTION
    (
     "A doppler shift for radiative transfer calculations.\n"
     "\n"
     "This scalar variable can hold the local doppler shift. It is intended\n"
     "mainly for communication with various methods and agendas, such as\n"
     "methods and agendas calculating absorption coefficients.\n"
     "\n"
     "Positive values mean that the spectrum is shifted towards higher\n"
     "frequencies\n"
     "\n"
     "Usage: Communication variable.\n"
     "\n"
     "Units: [ Hz ]\n"
     ),
    GROUP( "Numeric" )));
  
 wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_do_vmr_jacs" ),
      DESCRIPTION
      (
       "Index of gas species for which *rte_agenda* shall calculate "
       "VMR jacobians (with respect to changes along the propagation path).\n"
       "\n"
       "These indexes refer to the position in *abs_species*.\n"
       "\n"
       "Usage:   Set internally, by *RteCalc*.\n"
      ),
      GROUP( "ArrayOfIndex" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_do_t_jacs" ),
      DESCRIPTION
      (
       "Flag to *rte_agenda* to calculate jacobians for temperature.\n"
       "\n"
       "Usage:   Set internally, by *RteCalc*.\n"
      ),
      GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_edensity" ),
      DESCRIPTION
      (
       "An electron density for radiative transfer calculations.\n"
       "\n"
       "This scalar variable can hold the local electron density. It is\n"
       "intended mainly for communication with various methods and agendas,\n"
       "such as methods and agendas calculating refractive index.\n"
       "\n"
       "Usage: Communication variable.\n"
       "\n"
       "Units: [ 1/m3 ]\n"
       ),
      GROUP( "Numeric" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_los" ),
      DESCRIPTION
      (
       "A line-of-sight for radiative transfer calculations.\n"
       "\n"
       "The main purpose of this WSV and *rte_pos* is communication with \n"
       "different agendas involved in the RTE calculations.\n"
       "\n"
       "For 1D and 2D cases, *rte_los* is a vector of length 1 holding the \n"
       "zenith angle. For 3D, the length of the vector is 2, where the\n"
       "additional element is the azimuthal angle. These angles are defined\n"
       "in the ARTS user guide (AUG). Look in the index for \"zenith angle\"\n"
       "and \"azimuthal angle\".\n"
       "\n"
       "Usage: See above.\n"
       "\n"
       "Units: [ degree, degree ]\n"
       "\n"
       "Size:  [ 1 or 2 ]\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_pos" ),
      DESCRIPTION
      (
       "A geographical position for radiative transfer calculations.\n"
       "\n"
       "The main purpose of this WSV and *rte_los* is communication with \n"
       "different agendas involved in the RTE calculations.\n"
        "\n"
       "This variable is a vector with a length equalling the atmospheric\n"
       "dimensionality. The first element is the geomtrical altitude.\n"
       "Element 2 is the latitude and element 3 is the longitude.\n"
       "\n"
       "Usage: See above. \n"
       "\n"
       "Units: [ m, degree, degree ]\n"
       "\n"
       "Size:  [ atmosphere_dim ]\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_pos2" ),
      DESCRIPTION
      (
       "A second geographical position for radiative transfer calculations.\n"
       "\n"
       "Used when the propagation path is defined by two positions, instead\n"
       "of a position and a line-of-sight. That is, this variable basically\n"
       "replaces *rte_los* for the cases of consideration.\n"
       "\n"
       "As *rte_pos* with the exception that a \"latitude\" must also be\n"
       "specified for 1D. This is the angular distance to *rtte_pos*, where\n"
       "this distance is defined as the 2D-\"latitude\".\n"
       "\n"
       "Usage: See above. \n"
       "\n"
       "Units: [ m, degree, degree ]\n"
       "\n"
       "Size:  [ atmosphere_dim ]\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_pressure" ),
      DESCRIPTION
      (
       "A pressure for radiative transfer calculations.\n"
       "\n"
       "This scalar variable can hold the local pressure. It is intended\n"
       "mainly for communication with various methods and agendas, such as\n"
       "methods and agendas calculating absorption coefficients.\n"
       "\n"
       "Usage: Communication variable.\n"
       "\n"
       "Units: [ Pa ]\n"
       ),
      GROUP( "Numeric" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_temperature" ),
      DESCRIPTION
      (
       "A temperature for radiative transfer calculations.\n"
       "\n"
       "This scalar variable can hold the local temperature. It is intended\n"
       "mainly for communication with various methods and agendas, such as\n"
       "methods and agendas calculating absorption coefficients.\n"
       "\n"
       "Usage: Communication variable.\n"
       "\n"
       "Units: [ K ]\n"
       ),
      GROUP( "Numeric" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_vmr_list" ),
      DESCRIPTION
      (
       "A list of VMR values for radiative transfer calculations.\n"
       "\n"
       "This vector variable holds the local VMR value for all used species\n"
       "(as given by *abs_species*). It is intended mainly for communication\n"
       "with various methods and agendas, such as methods and agendas \n"
       "calculating absorption coefficients.\n"
       "\n"
       "Usage: Communication variable.\n"
       "\n"
       "Units: [ Absolute value ]\n"
       "\n"
       "Size:  Should match abs_species.nelem()\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "scat_aa_grid" ),
       DESCRIPTION
       (
        "Azimuthal angle grid.\n"
        "\n"
        "The azimutal angle grid, on which the intensity field is stored. \n"
        "This grid is used for RT calculations inside the cloudbox, \n"
        "therefore one has to define it if the cloudbox is activated by \n"
        "the flag *cloudbox_on*.\n"
        "The grid must be sorted in increasing order, with no repetitions.\n"
        "\n"
        "Usage:      Set by the user.\n"
        "\n"
        "Unit:       degrees \n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_aa_index" ),
      DESCRIPTION
      (
       "Azimuth angle index for scattering calculations.\n"
       "\n"
       "This variable is used in methods used for computing scattering\n"
       "properties. \n"
       "It holds the information about the azimuth angles for which the \n"
       "scattering calculations are done.  The angles used for computing \n"
       "scattering properties of particles can be different from that used \n"
       "for radiative transfer calculation. \n"
       "\n"
       "Usage:    Method output.\n"
       ),
     GROUP( "Index" ))); 

   wsv_data.push_back
     (WsvRecord
      ( NAME( "scat_data_meta" ),
        DESCRIPTION
        (
	 "Structure for the scattering meta data.\n"
	 "\n"
	 "This variable holds the scattering meta data for a single particle.\n"
	 "This data is needed for particle size distribution calculations.\n"
	 "\n"
	 "Currently some of the meta data entries are not used by ARTS,\n" 
	 "but were included for future extensions of pnd calculations in ARTS.\n"
	 "This applies to: \"shape\", \"A_proje\", \"asratio\"\n"
	 "\n"	 
	 "Usage: Set by the user.\n"
	 "\n"
	 "Dimensions/Units:\n"
	 "\tString[description]\t[particle description]\n"
	 "\tString[type]\t\t[''Ice'', ''Water''...]\n"
	 "\tString[shape]\t\t[''Droxtal'', ''Aggregate'', ''Spherical'', ...]\n"
	 "\tNumeric[density]\t[kg/m3]\n"
	 "\tNumeric[d_max]\t\t[m]\n"
	 "\tNumeric[V]\t\t[m3]\n"
	 "\tNumeric[A_proje]\t[m2]\n"
	 "\tNumeric[asratio]\t[]\n"
         ),
        GROUP( "ScatteringMetaData" ))); 

   wsv_data.push_back
     (WsvRecord
      ( NAME( "scat_data_meta_array" ),
        DESCRIPTION
        (
	 "An Array of ScatteringMetaData.\n"
	 "\n"
	 "The elements of the array hold the scattering meta data for a single particle.\n"
	 "This data is needed for particle size distribution calculations, with *pnd_fieldSetup*.\n"
	 "\n"
	 "Currently some of the meta data entries are not used by ARTS, but were included for \n"
	 "future extensions of pnd calculations in ARTS.\n"
	 "This applies to: \"shape\", \"A_proje\", \"asratio\"\n"
	 "\n"
	 "Note: This array must contain as many elements as *scat_data_raw*\n"
	 "\n"
	 "Usage: Set by the user.\n"
	 "\n"
	 "Dimensions/Units: Array[particle types]\n"
	 "\tString[description]\t[particle description]\n"
	 "\tString[type]\t\t[''Ice'', ''Water''...]\n"
	 "\tString[shape]\t\t[''Droxtal'', ''Aggregate'', ''Spherical'', ...]\n"
	 "\tNumeric[density]\t[kg/m3]\n"
	 "\tNumeric[d_max]\t\t[m]\n"
	 "\tNumeric[V]\t\t[m3]\n"
	 "\tNumeric[A_proje]\t[m2]\n"
	 "\tNumeric[asratio]\t[]\n"
         ),
        GROUP( "ArrayOfScatteringMetaData" ))); 


   wsv_data.push_back
     (WsvRecord
      ( NAME( "scat_data_mono" ),
        DESCRIPTION
        (
         "Monochromatic single scattering data.\n"
         "\n"
         "This variable holds the single scattering properties for all \n"
         "hydrometeor species. It is calculated from *scat_data_raw* by \n"
         "*scat_data_monoCalc*, which interpolates *scat_data_raw* for the \n"
         "required frequency.\n"
         ),
        GROUP( "ArrayOfSingleScatteringData" ))); 

   wsv_data.push_back
     (WsvRecord
      ( NAME( "scat_data_nelem" ),
        DESCRIPTION
        (
	  "This Array has the same size as *part_species*. \n"
	  "\n"
	  "Each element holds the number of scattering data files, associated\n"
	  "with the criteria given by *part_species*.\n"
	  "\n"
          "Usage: WSM *ScatteringParticlesSelect* creates *scat_data_nelem*.\n"
              ),
        GROUP( "ArrayOfIndex" ))); 
     
   wsv_data.push_back
     (WsvRecord
      ( NAME( "scat_data_raw" ),
        DESCRIPTION
        (
         "Raw data of single scattering data.\n"
         "\n"
         "This variable holds the single scattering properties for all \n"
         "hydrometeor species included in a calculation by using the \n"
         "methods *ParticleTypeAdd* or *ParticleTypeAddAll*. \n" 
         "For more information refer to ArtsWiki.\n"
	 "\n"
	 "Option: See also *ScatteringParticleTypeAndMetaRead* for reading \n"
	 "*SingleScatteringData* and *ScatteringMetaData*.\n"
	 "\n"
         "The unit of the single scattering properties is m^2.\n"
         "\n"
         "This may be used in combination with *scat_data_meta_array*\n"
         "\n"
         "Usage: Method ouput.\n"
         "\n"
         "Dimensions: Array[particle types] \n"
         "  SingleScatteringData \n"
         "  Enum[particle type attribute]\n"
         "  String[description] \n"
         "  Vector[f_grid]\n"
         "  Vector[T_grid]\n"
         "  Vector[za_grid]\n"
         "  Vector[aa_grid]\n"
         "  Tensor7[pha_mat_data]\n"
         "      [f_grid, T_grid, za_grid, aa_grid, za_grid, aa_grid, matrix_element]\n"
         "                       ^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^\n"
         "                       scattered         incoming\n"
         "  Tensor5[ext_mat_data]\n"
         "      [f_grid, T_grid, za_grid, aa_grid, matrix_element]\n"
         "  Tensor5[abs_vec_data]\n"
         "      [f_grid, T_grid, za_grid, aa_grid, matrix_element]\n"
         ),
        GROUP( "ArrayOfSingleScatteringData" ))); 
   
 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_i_lat" ),
      DESCRIPTION
      (
       "Intensity field on cloudbox boundary (constant latitude slice).\n"
       "\n"
       "This variable gives the intensity field from all directions defined \n"
       "in *scat_aa_grid* and *scat_za_grid* on each grid point on the two \n"
       "equal \n"
       "latitude levels of the cloudbox boundary. It contains all four \n"
       "components of the Stokes vector.\n"
       "\n"
       "This variable is used as interface between the clear sky and the \n"
       "scattering calculations. \n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      In/Output from/to *ScatteringDoit* and *ScatteringDisort*\n"
       "\n"
       "Unit:        W / (m^2 Hz sr) \n"
       "\n"
       "Dimensions: [ f_grid, p_grid, 2, lon_grid, \n"
       "              scat_za_grid \n  scat_aa_grid, stokes_dim ]\n"
       ),
      GROUP( "Tensor7" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_i_lon" ),
      DESCRIPTION
      (
       "Intensity field on cloudbox boundary (equal longitude slice).\n"
       "\n"
       "This variable gives the intensity field from all directions defined \n"
       "in *scat_aa_grid* and *scat_za_grid* on each grid point on the equal\n"
       "longitude level of the boundary of the cloudbox, which is defined \n"
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
       "Usage:      Output from *ScatteringDoit* and *ScatteringDisort* \n"
       "\n"
       "Unit:        W / (m^2 Hz sr) \n"
       "\n"
       "Dimensions: [ f_grid, p_grid, lat_grid, 2, \n"
       "              scat_za_grid, scat_aa_grid, stokes_dim]\n"
       ),
      GROUP( "Tensor7" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_i_p" ),
      DESCRIPTION
      (
       "Intensity field on cloudbox boundary (equal pressure slice).\n"
       "\n"
       "This variable gives the intensity field from all directions defined \n"
       "in *scat_aa_grid* and *scat_za_grid* on each grid point on the equal\n"
       "pressure levels of the cloudbox boundary. It contains all four \n"
       "components of the Stokes vector.\n"
       "\n"
       "This variable is used as interface between the clear sky and the \n"
       "scattering calculations. \n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      In/Output from *ScatteringDoit* and *ScatteringDisort* \n"
       "\n"
       "Unit:        W / (m^2 Hz sr) \n"
       "\n"
       "Dimensions: [ f_grid, 2, lat_grid, lon_grid, \n" 
       "              scat_za_grid, scat_aa_grid, stokes_dim]\n"
       ),
      GROUP( "Tensor7" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_lat_index" ),
      DESCRIPTION
      (
       "Latitude index for scattering calculations.\n"
       "\n"
       "This variable is used in methods used for computing scattering\n"
       "properties of particles like *opt_prop_sptFromData* and *pha_matCalc*.\n"
       "It holds the information about the position for which the \n"
       "scattering calculations are done. \n"
       "\n"
       "Usage:    Input to the methods *spt_calc_agenda*,\n"
       "                               *pha_mat_spt_agenda*\n"
       ),
     GROUP( "Index" ))); 

 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_lon_index" ),
      DESCRIPTION
      (
       "Longitude index for scattering calculations.\n"
       "\n"
       "This variable is used in methods used for computing scattering\n"
       "properties of particles like *opt_prop_sptFromData* and *pha_matCalc*.\n"
       "It holds the information about the position for which the \n"
       "scattering calculations are done.  \n"
       "\n"
       "Usage:    Input to the methods *spt_calc_agenda*,\n"
       "                               *pha_mat_spt_agenda*\n"
       ),
     GROUP( "Index" ))); 

  wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_p_index" ),
      DESCRIPTION
      (
       "Pressure index for scattering calculations.\n"
       "\n"
       "This variable is used in methods used for computing scattering\n"
       "properties of particles like *opt_prop_sptFromData* and *pha_matCalc*.\n"
       "It holds the information about the location for which the \n"
       "scattering calculations are done.\n"  
       "\n"
       "Usage:    Input to the methods *spt_calc_agenda*,\n"
       "                               *pha_mat_spt_agenda*\n"
       ),
     GROUP( "Index" ))); 
  
  wsv_data.push_back
    (WsvRecord
     ( NAME( "scat_za_grid" ),
       DESCRIPTION
       (
        "Zenith angle grid.\n"
        "\n"
        "The zenith angle grid, on which the intensity field is stored. \n"
        "This grid is used for RT calculations inside the cloudbox, therefore\n"
        "the grid has to be defined\n"
        "if the cloudbox is activated by the flag *cloudbox_on*.\n"
        "The grid must be sorted in increasing order, with no repetitions.\n"
        "\n"
        "Usage:      Set by the user.\n"
        "\n"
        "Unit:       degrees \n"
        ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_za_index" ),
      DESCRIPTION
      (
       "Zenith angle index for scattering calculations.\n"
       " \n"
       "This variable is used internally in WSMs for computing scattering \n"
       "properties. \n"
       "\n"
       "Usage:    Input to the agendas *spt_calc_agenda*, \n "
       "                               *pha_mat_spt_agenda*.\n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_description_amsu" ),
      DESCRIPTION
      (
	   "Sensor description for simple AMSU setup.\n"
	   "\n"
	   "This is a compact description of an AMSU-type sensor. The matrix\n"
	   "contains one row for each instrument channel. Each row contains three\n"
	   "elements: LO position [Hz], offset of the channel center from the LO\n"
	   "[Hz], and channel width [Hz].\n"
	   "\n"
	   "Usage: Set by the user.\n"
	   "\n"
	   "Unit: All entries in Hz.\n"
	   "\n"
	   "Size: [number of channels, 3]\n"
       ),
      GROUP( "Matrix" )));

	wsv_data.push_back
	(WsvRecord
	 ( NAME( "sensor_los" ),
      DESCRIPTION
      (
       "The sensor line-of-sight (LOS) for each measurement block.\n"
       "\n"
       "Line-of-sights are specified by giving the zenith and azimuth angles.\n"
       "Column 1 holds the zenith angle. This angle is simply the angle \n"
       "between the zenith and LOS directions. For 1D and 3D the valid\n"
       "range is [0 180], while for 2D angles down to -180 degrees are\n" 
       "allowed. Negative angles signifies for 2D observations towards\n"
       "lower latitudes, while positive angles means observations towards\n"
       "higher latitudes. Nadir corresponds throughout to 180 degrees.\n"
       "\n"
       "The azimuth angle is given with respect to the meridian plane. That\n"
       "is, the plane going through the north and south poles. The valid \n"
       "range is [-180,180] where angles are counted clockwise; 0 means\n"
       "that the viewing or propagation direction is north-wise and +90 means\n"
       "that the direction of concern goes eastward.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  [ degrees, degrees ]\n"
       "\n"
       "Size:  [ number of measurement blocks, 1 or 2 ]\n"
       ),
      GROUP( "Matrix" )));
	
  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_norm" ),
      DESCRIPTION
      (
       "Flag if sensor response should be normalised or not (0 or 1).\n"
       "\n"
       "If the flag is set to 1 each sensor response is normalised (where\n"
       "applicable). If set to 0 the sensor responses are left as provided.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a sub-entry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_pol" ),
      DESCRIPTION
      (
       "Sensor polarisations.\n"
       "\n"
       "The default for output is to give data for the selected Stokes\n"
       "elements (1:stokes_dim). This variable defines the polarisations\n"
       "that are actually measured (= or just what shall be outputted).\n"
       "This variable is used as input for WSM that handles the extraction\n"
       "of polarisation components. This in contrast to *sensor_response_pol*\n"
       "and that are used for internal bookkeeping\n"
       "\n"
       "The polarisation states/components are coded as\n"
       "   0 = Undefined.\n"
       "   1 = I, total intensity.\n"
       "   2 = Q, second Stokes component, Iv - Ih.\n"
       "   3 = U, third Stokes component, I+45 - I-45.\n"
       "   4 = V, forth Stokes component, Irc - Ilc\n"
       "   5 = Iv, intensity of vertically polarised component.\n"
       "   6 = Ih, intensity of horizontally polarised component.\n"
       "   7 = I+45, intensity of +45 deg linearly polarised component.\n"
       "   8 = I-45, intensity of -45 deg linearly polarised component.\n"
       "   9 = Irhc, intensity of right-hand circularly polarised component.\n"
       "  10 = Ilhc, intensity of left-hand circularly polarised component.\n"
       "\n"
       "See the documentation for definition of the Stokes vector and the\n"
       "different components.\n"
       "\n"
       "If the sensor measures the vertical and horisontal componenets, this\n"
       "variable shall accordingly be set to [5,6].\n"
       "\n"
       "The conversion to Planck-BT of components 2-4 requires that component\n"
       "1 is kept, and is put as first element.\n"
       "\n"
       "The shift from the Stokes vector can be made at any stage when of the\n"
       "sensor response set-up. The responses used must of course be adopted\n"
       "correspondingly. Or reversed, if the antenna response is defined for\n"
       "Iv and Ih it could be useful to shift polarisation as first sensor\n"
       "operation.\n"
       "\n"
       "Usage: Set by the user.\n"
       ),
      GROUP( "ArrayOfIndex" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_pos" ),
      DESCRIPTION
      (
       "The sensor position for each measurement block.\n"
       "\n"
       "The sensor positions are specified as a matrix, where the number of\n"
       "columns shall be equal to *atmosphere_dim*. Column 1 shall contain\n"
       "the altitude of observation posotion, column 2 the latitude and the \n"
       "last column the longitude. The number of rows corresponds to the\n"
       "number of measurement blocks.\n" 
       "\n"
       "Valid range for latitudes in 3D is [-90,90], while for 2D any value\n"
       "is accepted. Accepted range for longitudes are [-360,360].\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  [ m, degrees, degrees ]\n"
       "\n"
       "Size:  [ number of measurement blocks, atmosphere_dim ]\n"
       ),
      GROUP( "Matrix" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response" ),
      DESCRIPTION
      (
        "The matrix modelling the total sensor response.\n"
        "\n"
        "This matrix describes the sensor respons for one measurement block\n"
        "The response is assumed to be identical for each such block.\n"
        "\n"
        "The matrix is the product of all the individual sensor response\n"
        "matrices. Therefore its dimensions are depending on the total sensor\n"
        "configuration. The *sensor_response* has to initialised by the \n"
        "*sensor_responseInit* method.\n"
        "\n"
        "Usage: Output/input to the *sensor_response...* methods.\n"
        "\n"
        "Units: -\n"
        "\n"
        "Dimension: See the individual *sensor_response...* method. \n"
       ),
      GROUP( "Sparse" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_aa" ),
      DESCRIPTION
      (
       "The relative azimuth angles associated with the output of\n"
       "*sensor_response*.\n"
       "\n"
       "Definition of angle matches *mblock_aa_grid*. Works otherwise as\n"
       "*sensor_response_f*.\n"
       "\n"
       "The variable shall not be set manually, it will be set together with\n"
       "*sensor_response* by sensor response WSMs.\n"
       "\n"
       "Usage: Set by sensor response methods.\n"
       "\n"
       "Unit:  [ degrees ]\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_aa_array" ),
      DESCRIPTION
      (
       "As *sensor_response_aa* refers to *sensor_response_array*.\n"
       ),
      GROUP( "ArrayOfVector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_aa_grid" ),
      DESCRIPTION
      (
       "The azimuth angle grid associated with *sensor_response*.\n"
       "\n"
       "A variable for communication between sensor response WSMs. Matches\n"
       "initially *mblock_aa_grid*, but is later adjusted according to the\n"
       "sensor specifications. Only defined when a common grid exists. Values\n"
       "are here not repeated as in *sensor_response_aa*\n"
       "\n"
       "The zenith and azimuth dimensions are joined into a single dimension\n"
       "after the antenna. The variables *sensor_response_za_grid* and \n"
       "*sensor_response_aa_grid* have then the same length after the antenna\n"
       "(if antenna_dim = 2), holding data taken from the columns of \n"
       "*antenna_los*.\n"
       "\n"
       "Usage: Set by sensor response methods.\n"
       "\n"
       "Unit:  [ degrees ]\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_array" ),
      DESCRIPTION
      (
        "A set of *sensor_response* matrices.\n"
        "\n"
        "The main application of this variable should be to describe the\n"
        "sensor characterstics for different measurement blocks (then\n"
        "together with *sensor_response_index*.)\n"
       ),
      GROUP( "ArrayOfSparse" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_f" ),
      DESCRIPTION
      (
       "The frequencies associated with the output of *sensor_response*.\n"
       "\n"
       "This vector gives the frequency for each element of the measurement\n"
       "vector produced inside one measurement block. The frequencies of\n"
       "the total measurement vector, *y*, are obtained by repeating these\n"
       "frequencies n times, where n is the number of measurement blocks\n"
       "(e.g. the number of rows in *sensor_pos*).\n"
       "\n"
       "The variable shall not be set manually, it will be set together with\n"
       "*sensor_response* by sensor response WSMs.\n"
       "\n"
       "Usage: Set by sensor response methods.\n"
       "\n"
       "Unit:  [ Hz ]\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_f_array" ),
      DESCRIPTION
      (
       "As *sensor_response_f* refers to *sensor_response_array*.\n"
       ),
      GROUP( "ArrayOfVector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_f_grid" ),
      DESCRIPTION
      (
       "The frequency grid associated with *sensor_response*.\n"
       "\n"
       "A variable for communication between sensor response WSMs. Matches\n"
       "initially *f_grid*, but is later adjusted according to the sensor\n"
       "specifications. Only defined when a common grid exists. Values are\n"
       "here not repeated as in *sensor_response_f*\n"
       "\n"
       "Usage: Set by sensor response methods.\n"
       "\n"
       "Unit:  [ Hz ]\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_index" ),
      DESCRIPTION
      (
       "Identification of *sensor_response* for each measurement block.\n"
       "\n"
       "This variable shall hold an index number for each measurement block.\n"
       "Typically this index identifies which element of\n"
       "*sensor_response_array* that shall be applied.\n"
       "\n"
       "Zero based indexing is applied.\n"
       ),
      GROUP( "ArrayOfIndex" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_pol" ),
      DESCRIPTION
      (
       "The polarisation states associated with the output of\n"
       "*sensor_response*.\n"
       "\n"
       "Works basically as *sensor_response_f*.\n"
       "\n"
       "See *sensor_pol* for coding of polarisation states.\n"
       "\n"
       "The variable shall not be set manually, it will be set together with\n"
       "*sensor_response* by sensor response WSMs.\n"
       "\n"
       "Usage: Set by sensor response methods.\n"
       "\n"
       "Unit:  [ - ]\n"
       ),
      GROUP( "ArrayOfIndex" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_pol_array" ),
      DESCRIPTION
      (
       "As *sensor_response_pol* refers to *sensor_response_array*.\n"
       ),
      GROUP( "ArrayOfArrayOfIndex" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_pol_grid" ),
      DESCRIPTION
      (
       "The \"polarisation grid\" associated with *sensor_response*.\n"
       "\n"
       "A variable for communication between sensor response WSMs. It is\n"
       "initially 1:stokes_dim, but can later adjusted according to the \n"
       "sensor specifications. Only defined when a common grid exists. \n"
       "\n"
       "See *sensor_pol* for coding of polarisation states.\n"
       "\n"
       "Usage: Set by sensor response methods.\n"
       "\n"
       "Unit:  [ - ]\n"
       ),
      GROUP( "ArrayOfIndex" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_za" ),
      DESCRIPTION
      (
       "The relative zenith angles associated with the output of\n"
       "*sensor_response*.\n"
       "\n"
       "Definition of angle matches *mblock_za_grid*. Works otherwise as\n"
       "*sensor_response_f*.\n"
       "\n"
       "The variable shall not be set manually, it will be set together with\n"
       "*sensor_response* by sensor response WSMs.\n"
       "\n"
       "Usage: Set by sensor response methods.\n"
       "\n"
       "Unit:  [ degrees ]\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_za_array" ),
      DESCRIPTION
      (
       "As *sensor_response_za* refers to *sensor_response_array*.\n"
       ),
      GROUP( "ArrayOfVector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_za_grid" ),
      DESCRIPTION
      (
       "The zenith angle grid associated with *sensor_response*.\n"
       "\n"
       "A variable for communication between sensor response WSMs. Matches\n"
       "initially *mblock_za_grid*, but is later adjusted according to the\n"
       "sensor specifications. Only defined when a common grid exists. Values\n"
       "are here not repeated as in *sensor_response_za*\n"
       "\n"
       "The zenith and azimuth dimensions are joined into a single dimension\n"
       "after the antenna. The variables *sensor_response_za_grid* and \n"
       "*sensor_response_aa_grid* have then the same length after the antenna\n"
       "(if antenna_dim = 2), holding data taken from the columns of \n"
       "*antenna_los*.\n"
       "\n"
       "Usage: Set by sensor response methods.\n"
       "\n"
       "Unit:  [ degrees ]\n"
       ),
      GROUP( "Vector" )));

  /* Sensor rotation not yet updated 
  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_rot" ),
      DESCRIPTION
      (
       "The rotation of the sensor for each antenna line-of-sight.\n"
       "\n"
       "The rotation is the angle between the atmospheric and sensor frames\n"
       "for polarisation. The angle increases with clockwise rotation of the\n"
       "sensor when looking along the line-of-sight of the sensor. \n"
       "\n"
       "If the purpose of the simulations is to extract the polarisation\n"
       "of the radiation coming from the atmosphere (no sensor), the angles\n"
       "shall be set to 0.\n"
       "\n"
       "The size of the vector shall either be equal to the number of rows in\n"
       "*antenna_los* or be one. In the latter case the constant rotation\n"
       "will be applied for all antennae line-of-sight.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  [ degrees ]\n"
       "\n"
       "Size:  [ number of antennae or one ]\n"
       ),
      GROUP( "Vector" )));
  */

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_time" ),
      DESCRIPTION
      (
       "The time for each measurement block.\n"
       "\n"
       "This WSV is used when a time must be assigned to the measurements.\n"
       "No specific time format has (yet) been specified.\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  [ arbitrary ]\n"
       "\n"
       "Size:  [ number of measurement blocks ]\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sideband_mode" ),
      DESCRIPTION
      (
        "Description of target sideband.\n"
        "\n"
        "A text string describing which of the two sidebands (of a heterodyne\n"
        "instrument) that can be seen as \"main\" band. Possible choices are:\n"
        " \"lower\" : Low frequency sideband shall be considered as target.\n"
        " \"upper\" : High frequency sideband shall be considered as target.\n"
        "\n"
        "Usage: Set by the user.\n"
       ),
      GROUP( "String" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sideband_mode_multi" ),
      DESCRIPTION
      (
        "Description of target sideband for a multiple LO reciever.\n"
        "\n"
        "As *sideband_mode* but handles an instrument with several LO chains.\n"
        "See further *lo_multi* and *sideband_response_multi*. This length of\n"
        "this array must match the size of those WSVs.\n"
        "\n"
        "Usage: Set by the user.\n"
       ),
      GROUP( "ArrayOfString" )));


  wsv_data.push_back
   (WsvRecord
    ( NAME( "sideband_response" ),
      DESCRIPTION
      (
       "Description of (mixer) sideband response.\n"
       "\n"
       "This variable describes the response of each sideband of a heterodyne\n"
       "receiver. The response is given as a GriddedField1, with frequency as the\n"
       "grid. The actual data describe the sideband filter function at each\n"
       "frequency grid point. An interpolation is applied to obtain the\n"
       "response for other frequencies.\n"
       "\n"
       "The frequency grid should be given in terms of IF, with end points\n"
       "symmetrically placed around zero. That is, the grid must contain\n"
       "both negative and positive values. The sideband response (after \n"
       "summation with *lo*) is not allowed to extend outside the range\n"
       "for which spectral data exist (normally determined by *f_grid*).\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Dimensions: \n"
       "   GriddedField1:\n"
       "      Vector f_grid[N_f]\n"
       "      Vector data[N_f]\n"
       ),
      GROUP( "GriddedField1" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sideband_response_multi" ),
      DESCRIPTION
      (
       "Description of multiple (mixer) sideband responses.\n"
       "\n"
       "As *sideband_response* but describes an instrument with multiple\n"
       "mixers. An array element for each LO. The size of this variable and\n"
       "*lo_multi* shall match.\n"
       "\n"
       "Unit: Hz\n"
       "\n"
       "Usage: Set by the user.\n"
       ),
      GROUP( "ArrayOfGriddedField1" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "single_scattering_data" ),
      DESCRIPTION
      (
       "Structure for the  single scattering data.\n"
       "\n"
       "See futher the ArtsWiki documentation were the SingleScatteringData\n"
       "structure is disussed.\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Dimensions:  SingleScatteringData \n"
       "  Enum[particle type attribute]\n"
       "  String[description] \n"
       "  Vector[f_grid]\n"
       "  Vector[T_grid]\n"
       "  Vector[za_grid]\n"
       "  Vector[aa_grid]\n"
       "  Tensor7[pha_mat_data]\n"
       "      [f_grid, T_grid, za_grid, aa_grid, za_grid, aa_grid, matrix_element]\n"
       "                       ^^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^\n"
       "                       scattered         incoming\n"
       "  Tensor5[ext_mat_data]\n"
       "      [f_grid, T_grid, za_grid, aa_grid, matrix_element]\n"
       "  Tensor5[abs_vec_data]\n"
       "      [f_grid, T_grid, za_grid, aa_grid, matrix_element]\n"
       ),
      GROUP( "SingleScatteringData" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "stokes_dim" ),
      DESCRIPTION
      (
       "The dimensionality of the Stokes vector (1-4).\n"
       "\n"
       "Usage:      Set by the user.\n"
       ),
      GROUP( "Index" )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "spt_calc_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));

   wsv_data.push_back
     (WsvRecord
      ( NAME( "surface_emission" ),
        DESCRIPTION
        ( "The emission from the surface.\n"
          "\n"
          "See specific methods generating *surface_emission* and the user\n"
          "guide for more information.\n"
          "\n"
          "Dimensions: [ f_grid, stokes_dim ]\n"
         ), 
        GROUP( "Matrix" )));

   wsv_data.push_back
     (WsvRecord
      ( NAME( "surface_emissivity_DISORT" ),
        DESCRIPTION
        ( "The surface emissivity specified on lat_grid and lon_grid.\n"
          "\n"
          "Remnant from a first solution for surface emissivity fields.\n"
          "Should be replaced with more flexible solution allowing emissivity\n"
          "to vary with incidence angle.\n"
          "\n"
          "Dimensions: [ lat_grid, lon_grid ]\n"
         ), 
        GROUP( "Matrix" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "surface_los" ),
       DESCRIPTION
       (
        "Directions for which to calculate downwelling radiation when\n"
        "considering a surface reflection.\n"
        "\n"
        "See further the user guide.\n"
        "\n"
        "Units: degrees\n"
        "\n"
        "Size:  [ any number, 1 or 2 ]\n"
        ), 
       GROUP( "Matrix" )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "surface_rtprop_agenda" ),
       DESCRIPTION
       (
        "See agendas.cc.\n"
        ),
       GROUP( "Agenda" )));
  
  wsv_data.push_back
    (WsvRecord
     ( NAME( "surface_rmatrix" ),
       DESCRIPTION
       (
        "The reflection coefficients for the directions given by\n"
        "*surface_los* to the direction of interest.\n"
        "\n"
        "The rows and columns of this tensor holds the reflection\n"
        "coefficient matrix for one frequency and one LOS. The reflection\n"
        "coefficients shall take into accound the angular weighting if the\n"
        "downwelling radiation.\n"
        "\n"
        "See specific methods generating *surface_rmatrix* and the user guide\n"
        "for more information.\n"
        "\n"
        "Usage:      Input to methods for *surface_rtprop_agenda*."
        "\n"
        "Units:      -\n"
        "\n"
        "Dimensions: [ surface_los, f_grid, stokes_dim, stokes_dim ]\n"
        ), 
       GROUP( "Tensor4" )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "surface_skin_t" ),
      DESCRIPTION
      (
       "Surface skin temperature.\n"
       "\n"
       "This temperature shall be selected considering the radiative\n"
       "properties of the surface, and can differ from the \"bulk\"\n"
       "temperature.\n"
       "\n"
       "Usage:   Input to methods for *surface_rtprop_agenda*.\n"
       ),
      GROUP( "Numeric" )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "surface_reflectivity" ),
      DESCRIPTION
      (
       "Surface reflectivity, for a given position and angle.\n"
       "\n"
       "This variable describes the surface reflectivity at one position\n"
       "and one incidence angle. It works as *surface_scalar_reflectivity*\n"
       "but is also defined for vector radiative transfer.\n"
       "\n"
       "The first dimension of the variable shall either match *f_grid* or\n"
       "be 1. The later case is interpreted as the reflectivity is the same\n"
       "for all frequencies.\n"
       "\n"
       "Usage:   Input to some surface properties methods.\n"
       "\n"
       "Dimensions: [ f_grid or 1, stokes_dim, stokes_dim]\n"
       ),
      GROUP( "Tensor3" )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "surface_scalar_reflectivity" ),
      DESCRIPTION
      (
       "Surface reflectivity, assuming it can be described as a scalar value.\n"
       "\n"
       "This variable describes the surface reflectivity at one position\n"
       "and one incidence angle.\n"
       "\n"
       "The variable can only be used for scalar radiative transfer, that\n"
       "is, *stokes_dim* equals 1. Use *surface_reflectivity* for vector\n"
       "radiative transfer.\n"
       "\n"
       "The length of the vector shall either match *f_grid* or be 1. The \n"
       "later case is interpreted as the reflectivity is the same for all\n"
       "frequencies (ie. matches a constant vector).\n"
       "\n"
       "Usage:   Input to some surface properties methods.\n"
       "\n"
       "Dimensions: [ f_grid or 1]\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "test_agenda" ),
     DESCRIPTION
     (
      "See agendas.cc.\n"
     ),
     GROUP( "Agenda" )));
  
  wsv_data.push_back
   (WsvRecord
    ( NAME( "timer" ),
      DESCRIPTION
      (
       "Stores the starting time for time measurements.\n"
       ),
      GROUP( "Timer" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "transmitter_pos" ),
      DESCRIPTION
      (
       "Transmitter positions.\n"
       "\n"
       "Used for radio link calculations and gives then the position of the\n"
       "transmitting device. The corresponding positions of the reciever are\n"
       "given by *sensor_pos. The number of rows in *transitter_pos* and\n"
       "*sensor_pos* must be equal.\n"
       "\n" 
       "This WSV is also defined as *sensor_pos* regarding the content of the\n"
       "columns, accepted range for latitudes etc. With one exception, this\n"
       "WSV is demanded to have two columns also for 1D. The additional\n"
       "second value is the angular distance between the transmitter and the\n"
       "reciver. This angle is defined as \"latitude\" for 2D, with the\n"
       "sensor fixed at the angle of 0 degree. This definition fits\n"
       "*rte_pos2* and single transmitter positions should notmally be put\n"
       "into *rte_pos2*.\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  [ m, degrees, degrees ]\n"
       ),
      GROUP( "Matrix" )));

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
       "The temperature for a point between the grid crossings is obtained \n"
       "by (multi-)linear interpolation of the *t_field*.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      Output of *AtmFieldsCalc*.\n"
       "\n"
       "Unit:       K\n"
       "\n"
       "Dimensions: [ p_grid, lat_grid, lon_grid ]\n"
       ),
      GROUP( "Tensor3" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "t_field_raw" ),
      DESCRIPTION
      (
       "Raw data for atmospheric temperatures.\n"
       "\n"
       "This variable gives the atmospheric temperature as stored in the \n"
       "database for the atmospheric scenarios.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user by choosing a climatology.\n"
       "\n"
       "Unit:  K\n"
       "\n"
       "Size   GriddedField3 \n "
       "       [N_p] \n"
       "       [N_lat] \n"
       "       [N_lon] \n"
       "       [N_p, N_lat, N_lon] \n"
       ),
      GROUP( "GriddedField3" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "t_surface" ),
      DESCRIPTION
      (
       "The surface temperature.\n"
       "\n"
       "This variable holds the temperature of the surface at each latitude\n"
       "and longitude grid crossing. The normal case should be that this \n"
       "temperature field is interpolated to obtain *surface_skin_t*.\n"
       "Accordingly, for 1D cases it could be a better idea to specify\n"
       "*surface_skin_t* directly.\n"
       "\n"
       "These temperature shall be selected considering the radiative\n"
       "properties of the surface, and can differ from the \"bulk\"\n"
       "temperatures.\n"
       "\n"
       "Usage:      Set by user.\n"
       "\n"
       "Unit:       K\n"
       "\n"
       "Dimensions: [ lat_grid, lon_grid ]\n"
       ),
      GROUP( "Matrix" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "use_mean_scat_data" ),
      DESCRIPTION
      (
       "Flag to use same scattering properties for all frequencies.\n"
       "\n"
       "This flag is not considered y all scattering methods, but (at least)\n"
       "by *iyFOS* and *iyBeerLambertStandardCloudbox*.\n"
       "\n"
       "If set to 1, the scattering properties are extracted for a single\n"
       "frequency, and these properties are applied for all frequncies when\n"
       "performing the actual radiative transfer calculations. This can save\n"
       "considerable time. The option can be when the width of the band is\n"
       "small comapred to the mean frequency. The properties are extracted\n"
       "for the mean of min and max of *f_grid*.\n"
       "\n"
       "If set to 0, standard calculations are made (scattering properties\n"
       "extracted for each frequency).\n"
       "\n"
       "Usage:      Set by user.\n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "verbosity" ),
      DESCRIPTION
      (
       "ARTS verbosity.\n"
       "\n"
       "!!! UNDER CONSTRUCTION !!! Currently unused\n"
       "\n"
       "Usage:      Set by user.\n"
       ),
      GROUP( "Verbosity" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "vmr_field" ),
      DESCRIPTION
      (
       "VMR field.\n"
       "\n"
       "This variable gives the volume mixing ratio of the chosen gaseous \n"
       "species as a function of p_grid, lat_grid, lon_grid. \n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "\n"
       "Usage:      Calculated internally.\n"
       "\n"
       "Unit:        absolute numbers \n"
       "\n"
       "Dimensions: [species, p_grid, lat_grid, lon_grid]\n"
        ),
      GROUP( "Tensor4" ))); 

  wsv_data.push_back
   (WsvRecord
    ( NAME( "vmr_field_raw" ),
      DESCRIPTION
      (
       "VMR data for the chosen gaseous species.\n"
       "\n"
       "This variable contains the volume mixing ratios (VMR) for all \n"
       "chosen gaseous species. It includes the grids corresponding to the \n"
       "grids in the database. \n"
       "*vmr_field_raw* is an Array of Array of Tensor3. It contains one \n"
       "gridded field for each species which contains the data and \n"
       "also the grids.\n"
       "For the calculation the data is \n"
       "interpolated on *p_grid*, *lat_grid* and *lon_grid*\n"  
       "\n"
       "Usage: Output of *AtmRawRead*\n"
       "       Input to *AtmFieldsCalc*.\n"
       "\n"
       "Unit:  absolute number\n"
       "\n"
       "Size:  Array[N_pt]\n"
       "       GriddedField3 \n "
       "       [N_p] \n"
       "       [N_lat] \n"
       "       [N_lon] \n"
       "       [N_p, N_lat, N_lon] \n"
       ),
      GROUP( "ArrayOfGriddedField3" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "wmrf_channels" ),
      DESCRIPTION
      (
       "Channel selection for WMRF fast calculation.\n"
       "\n"
       "This variable can be used to select one or several instrument channels\n"
       "from the list of all possible channels. Zero-based indexing is used, so\n"
       "Channel 0 is the first instrument channel!\n"
       ),
      GROUP( "ArrayOfIndex" ))); 

  wsv_data.push_back
   (WsvRecord
    ( NAME( "wind_u_field" ),
      DESCRIPTION
      (
       "Meridional (second horizontal wind) component field for 3D.\n"
       "\n"
       "The east-west wind component. Air moving towards higher\n"
       "longitudes is a positive wind.\n"
       "\n"       
       "Can be set to be empty, which is interpreted as zero wind speed\n"
       "everywhere.\n"
       "\n"
       "Unit:       m/s\n"
       "\n"
       "Dimensions: [ p_grid, lat_grid, lon_grid ]  or [ 0 0 0 ].\n"
       ),
      GROUP( "Tensor3" ))); 

  wsv_data.push_back
   (WsvRecord
    ( NAME( "wind_v_field" ),
      DESCRIPTION
      (
       "Horizontal 1D, 2D or 3D wind components.\n"
       "\n"
       "This wind component is defined differently for the atmospheric\n"
       "dimensionalities.\n"
       "\n"       
       " 1D: Complete vertical wind. The wind is treated to always be aligned\n"
       "azimuthally with the observation direction. That is, the vertical\n"
       "winds converge/start at geographical position of the sensor. A wind\n"
       "going away from the sensor is treated as positive.\n"
       "\n"       
       " 2D: Complete vertical wind. The wind is treated to be aligned with\n"
       "the 2D cross-section the atmosphere represents. Air moving towards\n"
       "higher latitudes is considered as a positive wind.\n"
       "\n"       
       " 3D: The north-south wind component. Air moving towards higher\n"
       "latitudes is a positive wind.\n"
       "\n"       
       "Can be set to be empty, which is interpreted as zero wind speed\n"
       "everywhere.\n"
       "\n"
       "Unit:       m/s\n"
       "\n"
       "Dimensions: [ p_grid, lat_grid, lon_grid ] or [ 0 0 0 ]\n"
       ),
      GROUP( "Tensor3" ))); 

  wsv_data.push_back
   (WsvRecord
    ( NAME( "wind_w_field" ),
      DESCRIPTION
      (
       "Vertical wind component field.\n"
       "\n"
       "Upward moving air corresponds to a positive wind speed.\n"
       "\n"       
       "Can be set to be empty, which is interpreted as zero wind speed\n"
       "everywhere.\n"
       "\n"
       "Unit:       m/s\n"
       "\n"
       "Dimensions: [ p_grid, lat_grid, lon_grid ] or [ 0 0 0 ]\n"
       ),
      GROUP( "Tensor3" ))); 

  wsv_data.push_back
   (WsvRecord
    ( NAME( "wmrf_weights" ),
      DESCRIPTION
      (
       "The weights for a WMRF fast calculation.\n"
       "\n"
       "Weights are stored in a sparse matrix. This can be used as a\n"
       "sensor_response matrix.\n"
       "\n"
       "The dimension of the matrix is (nchan, nfreq), where nchan\n"
       "is the number of instrument channels and nfreq is the number\n"
       "of monochromatic frequencies.\n"
       ),
      GROUP( "Sparse" ))); 

  wsv_data.push_back
   (WsvRecord
    ( NAME( "xml_output_type" ),
      DESCRIPTION
      (
       "Flag to determine whether XML output shall be binary or ascii.\n"
       "\n"
       "This flag has to be set using the workspace method\n"
       "*output_file_formatSetAscii* or *output_file_formatSetBinary*.\n"
       "One of these methods MUST be called before writing the first\n"
       "output file.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "y" ),
      DESCRIPTION
      (
       "The measurement vector.\n"
       "\n"
       "This vector holds radiances averaged in frequency and spatially,\n"
       "and can contain many spectra appended. \n"
       "\n"
       "Usage: Output from radiative transfer calculations considering\n"
       "       sensor response.\n"
       "\n"
       "Unit:  Undefined. Possibilities include: K, W/(m^2 Hz sr) and\n "
       "       optical thickness.\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "y_aux" ),
      DESCRIPTION
      (
       "Data auxilary to *y*\n"
       "\n"
       "A general variable to provide data auxilary to *y*. The content of\n"
       "this variable depends on method selected for *iy_clearsky_agenda*\n"
       "and *iy_cloudbox_agenda*. For methods dealing with radiance data,\n"
       "this variable is used to report the transmission of the atmosphere.\n"
       "This transmission can be calculated in a simplified manner, or be\n"
       "set to special values to indicate e.g. an interception with the\n"
       "cloudbox. See the relevant WSMs for details, and also *iy_aux*. \n"
       "\n"
       "If created through *yCalc*, the weighting with *sensor_response*\n"
       "is included. A value of 999 indicates that no data have been\n"
       "provided by the agenda methods.\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "y_error" ),
      DESCRIPTION
      (
       "Estimate of calculation errors in *y*.\n"
       "\n"
       "This variable is used for providing an error estimate. The estimate\n"
       "covers only the actual calculation approach, and e.g. uncertainties\n"
       "in spectroscopic data or representation errors due to coarse grids\n"
       "are not considered. This means that the error should normally be\n"
       "zero for clear-sky cases where calculations are performed for the\n"
       "complete propagation path.\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "y_f" ),
      DESCRIPTION
      (
       "The frequencies associated with *y*.\n"
       "\n"
       "A value is returned for each element of *y*. Depending on the sensor\n"
       "set-up and number of measurement blocks, this can be a copy of\n"
       "*sensor_response_f*, sveral copies of this vector appended, or some\n"
       "other frequenices.\n"
       "\n"
       "Usage: Output from radiative transfer calculations considering\n"
       "       sensor response.\n"
       "\n"
       "Unit:  [ Hz ]\n"
       ),
      GROUP( "Vector" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "y_los" ),
      DESCRIPTION
      (
       "The line-of-sights associated with *y*.\n"
       "\n"
       "Definition of angles matches *sensor_los* (such as first column holds\n"
       "zenith angles), but gives actual observed LOS. That is, the values of\n"
       "both *sensor_los* and *antenna_los* are considered. Data are provided\n"
       "for each element of *y*, following y_f, and the number of rows equals\n"
       "the length of *y*.\n"
       "\n"
       "Usage: Output from radiative transfer calculations considering\n"
       "       sensor response.\n"
       "\n"
       "Unit:  [ degrees, degrees ]\n"
        ),
      GROUP( "Matrix" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "y_pol" ),
      DESCRIPTION
      (
       "The polarisation states associated with *y*.\n"
       "\n"
       "Data are provided for each element of *y*, following y_f, and the\n"
       "length of this variable and *y* is equal.\n"
       "\n"
       "See *sensor_pol* for coding of polarisation components.\n"
       "\n"
       "Usage: Output from radiative transfer calculations considering\n"
       "       sensor response.\n"
       "\n"
       "Unit:  [ - ]\n"
       ),
      GROUP( "ArrayOfIndex" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "y_pos" ),
      DESCRIPTION
      (
       "The sensor positions associated with *y*.\n"
       "\n"
       "Definition of positions matches *sensor_pos* (such as first column\n"
       "holds the altitude). Data are provided for each element of *y*,\n"
       "following y_f, and the number of rows equals the length of *y*.\n"
       "\n"
       "Usage: Output from radiative transfer calculations considering\n"
       "       sensor response.\n"
       "\n"
       "Unit:  [ degrees, degrees ]\n"
        ),
      GROUP( "Matrix" )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "y_unit" ),
       DESCRIPTION
       (
        "Unit for spectral values returned by radiative transfer methods.\n"
        "\n"
        "The basic unit is determined by the definition of background\n"
        "radiation and atmospheric and surface source terms. The standard\n"
        "choices give radiances in unit of [W/m2/Hz/sr]. This variable\n"
        "allows conversion to other units.\n"
        "\n"
        "Please note that there is no way for ARTS to check the actual unit\n"
        "of a spectral value. In principle, the unit can differ between the\n"
        "elements. The user must makes sure that any unit conversion is\n"
        "applied correctly, and in accordance with the calibration of the\n"
        "instrument of concern. Note especially that the brightness\n"
        "temperature of pure Stokes elements is obtained with a common power\n"
        "reference. In short, the ratio is Q/I is not changed by the\n"
        "conversion (beside a non-linear effect for PlanckBT). This in\n"
        "contrast to Iv and other single \"measurement\" single polarisation\n"
        "componenets, where a Tb of e.g. Iv=200K corresponds to half of the\n"
        "power for I=200K (total power). See further the documentation.\n"
        "See *y_pol* for defined polarisation componenets.\n"
        "\n"
        "Possible choices are:\n"
        " \"1\"             : No conversion.\n"
        " \"RJBT\"          : Conversion to Rayleigh-Jean brightness temperature.\n"
        " \"PlanckBT\"      : Conversion to Planck brightness temperature.\n"
        " \"W/(m^2 m sr)\"  : Conversion to [W/(m^2 m sr)] (radiance per wavelength unit).\n"
        " \"W/(m^2 m-1 sr)\": Conversion to [W/(m^2 m-1 sr)] (radiance per wavenumber unit).\n"
        "\n"
        "Usage: Set by the user.\n"
        ),
       GROUP( "String" )));
  
  wsv_data.push_back
   (WsvRecord
    ( NAME( "yb" ),
      DESCRIPTION
      (
       "The measurement vector for a single measurement block.\n"
       "\n"
       "Exactly as *y*, but holds data only for a single measurement block.\n"
       "\n"
       "Usage: Used internally.\n"
       ),
      GROUP( "Vector" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "ybatch" ),
      DESCRIPTION
      (
       "Batch of spectra.\n"
       "\n"
       "Each column of *ybatch* corresponds to a spectrum vector *y*. \n"
       "See further *ybatchCalc*.\n"
       "\n"
       "Usage: Most commonly produced by *ybatch*.\n"
       "\n"
       "Unit:  Undefined. Possibilities include: K, W/(m^2 Hz sr) and\n "
       "       optical thickness.\n"
       "\n"
       "Dimensions: (length(y), number of batch cases)\n"
       ),
      GROUP( "Matrix" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "ybatch_calc_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( "Agenda" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "ybatch_index" ),
      DESCRIPTION
      (
       "Index of batch case.\n"
       "\n"
       "See further *ybatchCalc*.\n"
       "\n"
       "Usage: Set by *ybatchCalc*, for communication with\n"
       "       *ybatch_calc_agenda*.\n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
  (WsvRecord
   ( NAME( "ybatch_jacobians" ),
    DESCRIPTION
    (
     "All the Jacobians associated with ybatch.\n"
     "\n"
     "The batch index here is the array dimension.\n"
     "\n"
     "Usage: Most commonly produced by *ybatch*.\n"
     "\n"
     "Unit:  Depends on unit of y and on Jacobian type.\n"
     "\n"
     "Dimensions: [number of batch cases] \n"
     "             (length(y), \n"
     "             number of retrieval quantities and grids)\n" 
     ),
    GROUP( "ArrayOfMatrix" )));
  
  wsv_data.push_back
   (WsvRecord
    ( NAME( "ybatch_n" ),
      DESCRIPTION
      (
       "Number of batch cases for *ybatchCalc*.\n"
       "\n"
       "See further *ybatchCalc*.\n"
       "\n"
       "Usage: Input to *ybatchCalc*.\n"
       ),
      GROUP( "Index" )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "ybatch_start" ),
      DESCRIPTION
      (
       "Start index for *ybatchCalc*.\n"
       "\n"
       "This is set to a default of zero in *general.arts*.\n"
       "\n"
       "See further *ybatchCalc*.\n"
       "\n"
       "Usage: Input to *ybatchCalc*.\n"
       ),
      GROUP( "Index" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "z_field" ),
      DESCRIPTION
      (
       "The field of geometrical altitudes.\n"
       "\n"
       "This variable gives the geometrical altitude, above the ellipsoid, of\n"
       "each crossing of the pressure, latitude and longitude grids. For 1D\n"
       "cases the altitudes give the geometrical position of the pressure\n"
       "levels.\n"
       "\n"
       "For each geographical position (lat,lon), the values must be sorted\n"
       "in increasing order, with no repetitions. Otherwise the altitudes\n"
       "can be set to arbitrary values. Hydrostatic equilibrium is not\n"
       "applied automatically. If hydrostatic equilibrium applies, *z_field*\n"
       "must be set by a method ensuring that this criterium is fulfilled.\n"
       "\n"
       "The radius (from the coordinate centre) for a point between the grid\n"
       "crossings is obtained by a (multi-)linear interpolation of the sum\n"
       "of the ellipsoid radius and *z_field*.\n" 
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      Output of *AtmFieldsCalc*\n"
       "\n"
       "Unit:       m\n"
       "\n"
       "Dimensions: [ p_grid, lat_grid, lon_grid ]\n"
       ),
      GROUP( "Tensor3" )));

//    wsv_data.push_back
//    (WsvRecord
//     ( NAME( "gfield3" ),
//       DESCRIPTION
//       (
//        "Variable for testing the new gridded fields implementation.\n"
//        ),
//       GROUP( "GriddedField3" )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "z_field_raw" ),
      DESCRIPTION
      (
       "Raw data for geometrical altitudes.\n"
       "\n"
       "This variable gives the geometrical altitudes as stored in the \n"
       "database for atmospheric scenarios.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user by choosing a climatology.\n"
       "\n"
       "Unit:  K\n"
       "\n"
       "Size   GriddedField3 \n "
       "       [N_p] \n"
       "       [N_lat] \n"
       "       [N_lon] \n"
       "       [N_p, N_lat, N_lon] \n"
       ),
      GROUP( "GriddedField3" )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "z_hse_accuracy" ),
      DESCRIPTION
      (
       "Minimum accuracy for calculation of hydrostatic equilibrium.\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  m\n"
       ),
      GROUP( "Numeric" )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "z_surface" ),
      DESCRIPTION
      (
       "The surface altitude.\n"
       "\n"
       "This variable defines the shape of the surface, by giving the\n"
       "geometrical altitude above the geiod for each crossing of the \n"
       "latitude and longitude grids. Any shape of the surface is accepted.\n"
       "No gap between the surface and the lowermost pressure level is \n"
       "allowed.\n"
       "\n"
       "The radius (from the coordinate centre) for a point between the grid\n"
       "crossings is obtained by a linear (1D) or bi-linear (2D) \n"
       "interpolation of the sum of the ellipsoid radius and *z_surface*.\n"
       "That is, the radius for the surface is assumed to vary linear along \n"
       "the latitudes and longitudes in *lat_grid* and *lon_grid*.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      Set by user.\n"
       "\n"
       "Unit:       m\n"
       "\n"
       "Dimensions: [ lat_grid, lon_grid ]\n"
       ),
      GROUP( "Matrix" )));
}


//! Get index of WSV
/** 
 Returns the index the Workspace of the given WSV.
 
 \param[in]  name   WSV name
 \returns           Index in Workspace
 
 \author Oliver Lemke
 */
Index get_wsv_id(const String& name)
{
  map<String, Index>::const_iterator it = Workspace::WsvMap.find (name);
  if (it == Workspace::WsvMap.end())
    return -1;
  else
    return it->second;
}


//! Get index of WSV
/** 
 Returns the index the Workspace of the given WSV.
 
 Convenience function which can be called from within the debugger because it
 takes a plain char pointer instead of a String object as input.
 
 \param[in]  name   WSV name
 \returns           Index in Workspace
 
 \author Oliver Lemke
 */
Index get_wsv_id(const char *name)
{
  return get_wsv_id(String(name));
}
