/* Copyright (C) 2000-2007
   Stefan Buehler <sbuehler@ltu.se>
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
#include "matpackII.h"
#include "matpackIII.h"
#include "matpackVI.h"
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

  // New name: abs_coef
  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_coef" ),
       DESCRIPTION
       (
        "The matrix of total absorption coefficients.\n"
        "\n"
        "FIXME: Is this used much?\n"
        "\n"
        "Dimensions: [f_grid, abs_p]\n"
        "\n"
        "Unit: 1/m\n"
        ),
      GROUP( Matrix_ )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "abs_coef_per_species" ),
      DESCRIPTION
      (
       "These are the absorption coefficients individually for each\n"
       "tag group. The Array contains one matrix for each tag group,\n"
       "the matrix format is the same as that of abs_coef\n"
      ),
      GROUP( ArrayOfMatrix_ )));

//   wsv_data.push_back
//    (WsvRecord
//     ( NAME( "abs_coef_agenda" ),
//       DESCRIPTION
//       (
//         "See agendas.cc.\n"
//        ),
//       GROUP( Agenda_)));
  
  wsv_data.push_back
    (WsvRecord
     (NAME( "abs_cont_models" ),
      DESCRIPTION
      (
       "Continuum / full model absorption model description parameter.\n"
       "See the WSV `abs_cont_names' for a detailed description\n"
       "of the allowed continuum models. There should be one string here\n"
       "for each entry in `abs_cont_names'.See also the online" 
       "documentation in arts/doc/doxygen/html/continua_cc.html.\n"
      ),
      GROUP( ArrayOfString_ )));

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
      GROUP( ArrayOfString_ )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "abs_cont_parameters" ),
      DESCRIPTION
      (
       "Continuum model parameters. See the WSV `abs_cont_names'\n"
       "for a detailed description of the allowed continuum models. There\n"
       "should be one parameter vector here for each entry in\n"
       "`abs_cont_names'. See also the online documentation in\n"
       "arts/doc/doxygen/html/continua_cc.html.\n"
      ),
      GROUP( ArrayOfVector_ )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "abs_h2o" ),
      DESCRIPTION
      (
       "The total water profile associated with the pressures in abs_p [-]\n"
      ),
      GROUP( Vector_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_lines" ),
       DESCRIPTION
       (
        "A list of spectral line data.\n"
       ), 
       GROUP( ArrayOfLineRecord_ )));

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
       GROUP( ArrayOfLineshapeSpec_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_lines_per_species" ),
       DESCRIPTION
       (
        "A list of spectral line data for each tag.\n"
        "Dimensions: (tag_groups.nelem()) (# of lines for this tag)\n"
       ), 
       GROUP( ArrayOfArrayOfLineRecord_ )));

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
        "class GasAbsLookup for details. FIXME: Add here a reference to AUG,\n"
        "once the chapter on the lookup table has been written.\n"
        ), 
       GROUP( GasAbsLookup_ )));

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
          "species, or some other species that uses *h2o_abs*, that is, for which\n"
          "the absorption coefficient depends directly on water vapor.\n"
          "\n"
          "See user guide and online documentation of *abs_pts* and *abs_lookupCreate*\n"
          "for more details and usage examples.\n"
          ), 
         GROUP( ArrayOfArrayOfSpeciesTag_ )));

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
         GROUP( Vector_ )));

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
         GROUP( Vector_ )));

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
       GROUP( Index_ )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "abs_n2" ),
      DESCRIPTION
      (
       "The total nitrogen profile associated with the pressures in abs_p [-]\n"
      ),
      GROUP( Vector_ )));

  wsv_data.push_back
    (WsvRecord
    ( NAME( "abs_p" ),
      DESCRIPTION
      (
       "List of pressures to be used for the calculation of absorption\n"
       "coefficients. \n"
       "\n"
       "This can be copied from the global p_grid, but could also be\n"
       "different. \n"
       "\n"
       "Any absorption method should check that the length of this vector\n"
       "is the same as that of abs_t\n"
       "\n"
       "Dimension: [number of pressures]\n"
       "\n"
       "Unit: Pa\n"
       ),
      GROUP( Vector_ )));

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
      GROUP( Matrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "abs_scalar_gas_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_)));
  
  wsv_data.push_back
   (WsvRecord
    ( NAME( "abs_field" ),
      DESCRIPTION
      (
       "Scalar gas absorption field.\n"
       "\n"
       "Contains the scalar gas absorption for all species as a function of\n"
       "f_grid, p_grid, lat_grid, and lon_grid. \n"
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
      GROUP( Tensor5_ ))); 

  wsv_data.push_back
    (WsvRecord
     ( NAME( "abs_species" ),
       DESCRIPTION
       (
        "Tag groups for scalar gas absorption.\n"
        "\n"
        "This is an array of arrays of SpeciesTag tag definitions. It defines the\n"
        "available tag groups for the calculation of scalar gas absorption\n"
        "coefficients.  See online documentation of method *abs_speciesSet* for\n"
        "more detailed information how tag groups work and some examples.\n"
        ), 
       GROUP( ArrayOfArrayOfSpeciesTag_ )));

  wsv_data.push_back
    (WsvRecord
    ( NAME( "abs_t" ),
      DESCRIPTION
      (
       "List of temperatures to be used for the calculation of absorption\n"
       "coefficients.\n"
       "\n"
       "In contrast to the global t_field, this is just a vector. Any\n"
       "absorption method should check that the length of this vector is the\n"
       "same as that of abs_p\n"
       "\n"
       "Dimension: [number of pressures]\n"
       "\n"
       "Unit: K\n"
       ),
      GROUP( Vector_ )));

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
       GROUP( Matrix_ )));

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
        "Dimensions: [part_types,stokes_dim]\n"
        ),
       GROUP( Matrix_ ) ));

  wsv_data.push_back
    (WsvRecord
     (NAME( "abs_vmrs" ),
      DESCRIPTION
      (
       "The VMRs (unit: absolute number) on the abs_p grid.\n"
       "Dimensions: [tag_groups.nelem(), abs_p.nelem()]\n"
      ),
      GROUP( Matrix_ )));

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
       GROUP( ArrayOfMatrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "antenna_diagram" ),
      DESCRIPTION
      (
       "The antenna diagram.\n"
       "\n"
       "The diagram is described by an ArrayOfArrayOfMatrix. The highest\n"
       "level corresponds to different viewing angles of an multiple antenna\n"
       "array or multiple beam antenna. The next level, i.e. the elements\n"
       "of the different ArrayOfMatrix corresponds to different polarisations.\n"
       "The individual matrices then describes the antenna diagrams, with the\n"
       "first column describing a relative angle grid and the antenna gain\n"
       "given for different frequencies in the consecutive columns.\n"
       "\n"
       "For each level in the antenna diagram there is a choice to provide\n"
       "a single ArrayOfMatrix/Matrix/column that will be used for all\n"
       "existing viewing angles/polarisations/frequencies, or to provide a\n"
       "full description for each viewing angle/polarisation/frequency.\n"
       "\n"
       "Usage:      Set by the user.\n"
       ),
      GROUP( ArrayOfArrayOfMatrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "antenna_dim" ),
      DESCRIPTION
      (
       "The dimensionality of the antenna pattern (1-2).\n"
       "\n"
       "A dimensionality of 1 means that only the respons variation in the\n"
       "zenith direction is considered. The respons is then integrated in \n"
       "the azimuth direction. For 2D, the respons of the antenna has both a\n"
       "zenith and azimuth variation.\n"
       "\n"
       "Usage:      Set by the user.\n"
       ),
      GROUP( Index_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "antenna_los" ),
      DESCRIPTION
      (
       "The line-of-sight of individual antennae.\n"
       "\n"
       "This variable describes the line-of-sight of the individual antennae\n"
       "relative to *sensor_los*. If only one antenna is present the matrix\n"
       "should contain a row of zero(s). The number of columns corresponds to\n"
       "the *antenna_dim*, with the first column containing zenith angles\n"
       "and the second azimuth angles.\n"
       "\n"
       "The number of rows also describes the number of antennae to consider\n"
       "and shall be consistent with the number of elements in *sensor_rot*.\n"
       "A special case is when *sensor_rot* only has one element, see\n"
       "*sensor_rot* for more details.\n"
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
      GROUP( Matrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "arrayofmatrix_1" ),
      DESCRIPTION
      (
       "An arbitrary array of matrices.\n"
       "\n"
       "This variable shall be treated as a general variable of type\n"
       "ArrayOfMatrix. It can be used, for example, when some intermediate\n"
       "data must be generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( ArrayOfMatrix_ )));


  wsv_data.push_back
    (WsvRecord
     ( NAME( "arrayofstring_1" ),
       DESCRIPTION
       (
        "An arbitrary array of strings.\n"
        "\n"
        "This variable shall be treated as a general variable of type\n"
        "ArrayOfString. \n"
        "\n"
        "Usage: Set by user.\n"
        ),
       GROUP( ArrayOfString_ )));

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
       "changed. However, not all methods are working for higher dimesnions.\n"
       "\n"
       "The atmospheric dimensionalities (1D, 2D and 3D) are defined in the\n"
       "user guide (look for \"atmospheric dimensionality\" in the index).\n"
       "\n"
       "Usage:      Set by the user.\n"
       ),
      GROUP( Index_ )));

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
       "The data is stored in a GriddedField4.\n"
       "\n"
       "The order of the fields must be:\n"
       "T[K] z[m] VMR_1[1] ... VMR[2]\n"
       "\n"
       "It is up to the user to decide which VMR belongs to which species in\n"
       "an ARTS calculation. The method *AtmFieldsFromCompact* will\n"
       "simply split up the data found here, and create a *vmr_field* from the\n"
       "various VMR profiles.\n"
       "\n"
       "Usage: Used inside batch calculations, to hold successive atmospheric\n"
       "       states from an ArrayOfGriddedField4.\n"
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
      GROUP( GriddedField4_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "ybatch_calc_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_ )));

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
       "Usage:      Set by the user.\n"
       ),
      GROUP( Index_ )));

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
       "surfaces. For 2D, the angular extension of the cloud box is between \n"
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
       "*cloudbox_limits* = [0 5 4 11] means that cloud box extends between \n"
       "pressure levels 0 and 5, and latitude points 4 and 11.\n"
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
      GROUP( ArrayOfIndex_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "diy_dt" ),
      DESCRIPTION
      (
       "Derivative of *iy* with respect to temperature along the ppaths.\n"
       "\n"
       "This variable holds the derivative of monochromatic pencil beam\n"
       "radiances with respect to the temperature at each point along each \n"
       "propagation path part.\n"
       "\n"
       "The number of books is here always 1. This extra dimension is\n"
       "included in order to have the same data type for all diy-variables.\n"
       "\n"
       "Usage:      Set by *iy_calc* and *rte_agenda*.\n"
       "\n"
       "Dimensions: [ppath_array][ 1 ppath.np, f_grid, stokes_dim ]\n"
       ),
      GROUP( ArrayOfTensor4_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "diy_dvmr" ),
      DESCRIPTION
      (
       "Derivative of *iy* with respect to VMR along the propagation path.\n"
       "\n"
       "This variable holds the derivative of monochromatic pencil beam\n"
       "radiances with respect to the VMR of each species at each point\n"
       "along each propagation path part.\n"
       "\n"
       "Usage:      Set by *iy_calc* and *rte_agenda*.\n"
       "\n"
       "Dimensions: \n"
       "     [ppath_array][ rte_do_vmr_species, ppath.np, f_grid, stokes_dim ]\n"
       ),
      GROUP( ArrayOfTensor4_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_conv_flag" ),
      DESCRIPTION
      (
       "Flag for the convergence test.\n"
       "\n"
       "This variable is initialized with 0 inside the method \n"
       "*i_fieldIterate*.\n"
       "If after an iteration the convergence test is fulfilled, 1 is \n"
       "assigned which means that the iteration is completed. \n"
       "\n"
       "Usage: Method output. \n"
      ), 
      GROUP( Index_ ))); 

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_conv_test_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_ )));

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
       GROUP( Tensor6_ )));

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
       "Usage: Output of *CloudboxFieldPut*\n"
       "\n"
       "Unit: W / (m^2 Hz sr)\n"
       "\n"
        "Size: [N_f \n"
       "       (cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
       "        N_za, N_aa, N_i ]\n"
       ),
      GROUP( Tensor4_ )));
 
 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_i_field_old" ),
      DESCRIPTION
      (
       "Intensity field inside the cloudbox.\n"
       "\n"
       "This variable is used to store the intensity field inside the\n"
       "cloudbox while performing the iteration. One has to store the\n"
       "intensity field of the previos iteration to be able to do the \n"
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
      GROUP( Tensor6_ )));
 
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
      GROUP( Index_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_mono_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_ )));
 
 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_rte_agenda" ),
      DESCRIPTION
      (
       "See agendas.cc.\n"
       ),
      GROUP( Agenda_ ))); 

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_scat_field_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_ ))); 

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_scat_field" ),
      DESCRIPTION
      (
       "Scattered field field inside the cloudbox.\n"
       "\n"
       "This variable holds the value of the scattering integral.\n"
       "for all points inside the cloudbox. \n"
       "For more information refer to AUG.\n"
       "\n"
       "Usage: Input to *doit_i_fieldUpdateXXX*. \n"    
       "\n"
       "Unit: W / (m^2 Hz sr) for each Stokes component.\n"
       "\n"
       "Size: [(cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
       "       (cloudbox_limits[3] - cloudbox_limits[2]) +1, \n"
       "       (cloudbox_limits[5] - cloudbox_limits[4]) +1, \n"
       "        N_za, N_aa, N_i ]\n"
       ),
      GROUP( Tensor6_ )));   

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_za_grid_opt" ),
      DESCRIPTION
       (
        "Optimized zenith angle grid.\n"
        "\n"
        "Output of the method *DoitGridOptimization*. For limb simulations \n"
        "with scattering it is very \n"
        "important to use optimized grids for both, accuracy and speed of \n"
        "the calculations.\n"
        "\n"
        "Usage:   Output of *DoitGridOptimization*   \n"
        "\n"
        "Unit:    degrees \n"
        ),
      GROUP( Vector_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_za_grid_size" ),
      DESCRIPTION
      (
       "This vector contains the discretisation of the zenith angle grid \n"
       "for the scattering integral calculation. \n"
       "\n"
       "The zenith angle grid is defined from 0 to 180.\n"
       "*doit_za_grid_size* is the number of points of the zenith angle grid.\n"
       "\n"
       "Usage: Output of *DoitAngularGridsSet*.\n"
       ),
      GROUP( Index_ )));
 
 wsv_data.push_back
   (WsvRecord
    ( NAME( "doit_za_interp" ),
      DESCRIPTION
      (
       "Flag for interplation method in zenith angle dimension.\n"
       "\n"
       "0 - linear interpolation \n"
       "1 - cubic interpolation \n"
       "Default is linear interpolation. \n"
       "\n"
       "Usage: Set by user in *doit_za_interpSet*. \n"
       ),
      GROUP( Index_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "emission" ),
      DESCRIPTION
      (
       "Thermal emission source term.\n"
       "\n"
       "This variable holds the emission at one position along\n"
       "the propagation path. This source term is used for calculations\n"
       "inside *rte_agenda*. Inside scattering methods, such as DOIT,\n"
       "the calculation of the source term can be hard coded.\n"
       "\n"
       "Usage:      Set by *emission_agenda.\n"
       "\n"
       "Unit:       W / (m^2 Hz sr) or optical thickness \n"
       "\n"
       "Dimensions: [ f_grid ]\n"
       ),
      GROUP( Vector_ )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "emission_agenda" ),
       DESCRIPTION
       (
        "See agendas.cc.\n"
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
       "is used in the RT calculation in the cloudbox . It is \n"
       "the physical extinction matrix which includes particles extinction \n"
       "for all chosen particle types and gaseous extinction for all chosen \n"
       "gaseous species.\n" 
       "The matrix is calculated by the agendas *opt_prop_gas_agenda* \n"
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
       "Unit:       [Hz, m^2, m^2] "
       "\n"
       "Dimensions: [f_grid, stokes_dim, stokes_dim]\n"
       ),
       GROUP( Tensor3_ )));

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
       "Dimensions: [part_types, stokes_dim, stokes_dim]\n"
       ),
      GROUP( Tensor3_ )));

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
        GROUP( Index_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "f_backend" ),
       DESCRIPTION
       (
        "The frequency grid of the backend channels.\n"
        "\n"
        "Usage:      Input to *sensor_responseBackend*.\n "
        "\n"
        "Unit:       Hz\n"
        ),
        GROUP( Vector_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "f_grid" ),
       DESCRIPTION
       (
        "The frequency grid for monochromatic pencil beam calculations.\n"
        "\n"
        "Usage:      Set by the user.\n "
        "\n"
        "Unit:       Hz\n"
        ),
        GROUP( Vector_ )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "f_index" ),
      DESCRIPTION
      (
       "Frequency index. \n"
       "\n"
       "The calculations inside the cloudbox are only done for one frequency\n"
       "at a time. Some methods used for scattering calculation require the \n"
       "frequency. *f_index* holds the information, for which frequency the \n"
       "scattering calcultations are performed.\n"
       "\n"
       "*f_index* is an input to *opt_prop_gas_agenda* and \n"
       "*opt_prop_part_agenda*. In the clearsky case *opt_prop_gas_agenda* \n"
       "has to calculate optical properties for all frequencies defined\n"
       "in *f_grid*, in this case *f_index* has to be set to -1.\n "
       "\n"
       "Usage:      Input and output to *scat_mono_agenda*,\n"
       "                     *opt_prop_gas_agenda\n"
       "                     *opt_prop_part_agenda*.\n"
       ),
      GROUP( Index_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "f_mixer" ),
       DESCRIPTION
       (
        "The frequency grid of output from the mixer.\n"
        "\n"
        "Usage:      Input and output in *sensor_responseMixer*.\n "
        "\n"
        "Unit:       Hz\n"
        ),
        GROUP( Vector_ )));
  
 wsv_data.push_back
   (WsvRecord
    ( NAME( "geomag_los_calc_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_ )));


  wsv_data.push_back
    (WsvRecord
    ( NAME( "geomag_los" ),
      DESCRIPTION
      (
       "Magnetic field along the line of sight\n"
       "\n"
       "more text by Nikolay \n"
       "\n"
       "Unit: ..."
       "\n"
       "Dimensions: [Magnetic field B, angle between B and los] \n"
       ),
      GROUP( Matrix_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "index_1" ),
      DESCRIPTION
      (
       "An arbitrary index.\n"
       "\n"
       "This variable shall be treated as a general variable of type Index.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( Index_ )));

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
       "considered or not.\n"
       "\n"
       "Usage:      Used by radiative transfer methods.\n"
       "\n"
       "Unit:       W / (m^2 Hz sr) or transmission.\n"
       "\n"
       "Dimensions: [ f_grid, stokes_dim ]\n"
       ),
      GROUP( Matrix_ )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "iy_cloudbox_agenda" ),
       DESCRIPTION
       (
        "See agendas.cc.\n"
        ),
       GROUP(  Agenda_ )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "iy_space_agenda" ),
       DESCRIPTION
       (
        "See agendas.cc.\n"
        ),
       GROUP(  Agenda_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "i_space" ),
      DESCRIPTION
      (
       "Monochromatic intensities at the top of the atmosphere.\n"
       "\n"
       "The matrix holds the intensity entering the atmosphere from above\n"
       "along a propagation path. This should normally correspond to cosmic\n"
       "background radiation, but could also be radiation from the sun or \n"
       "a transmitting satellite.\n"
       "\n"
       "Usage:      Input to *iy_space_agenda*.\n"
       "\n"
       "Unit:       W / (m^2 Hz sr) \n"
       "\n"
       "Dimensions: [ f_grid, stokes_dim ]\n"
       ),
      GROUP( Matrix_ )));

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
        "See the online help. The calculation is controlled by an agenda,\n"
        "and is performed by *jacobianCalc*.\n"
        "\n"
        "Units:   See the different retrieval quantities.\n"
        "\n"
        "Dimension:   [ y, number of retrieval quantities and grids ]\n"
      ),
      GROUP( Matrix_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "jacobian_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_ )));

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
      GROUP( ArrayOfArrayOfIndex_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "jacobian_particle_update_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "jacobian_lat_grid" ),
      DESCRIPTION
      (
       "The latitude grid for the Jacobian matrix.\n"
       "\n"
       "The latitudes for which the Jacobian is determined. The jacobian\n"
       "is undefined outside the range covered by the grid. The grid must\n"
       "be sorted in increasing order, with no repetitions.\n"
       "\n"
       "Geocentric latitudes shall be used.\n"
       "\n"
       "The latitude grid can differ between retrieval quantities, and it\n"
       "is the grid defined when the retrieval quantity is added that is\n"
       "used.\n"
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
       "Usage:      Set by the user.\n"
       "\n"
       "Unit:       degrees\n"
      ),
      GROUP( Vector_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "jacobian_lon_grid" ),
      DESCRIPTION
      (
       "The longitude grid for the Jacobian matrix.\n"
       "\n"
       "The longitudes for which the Jacobian is determined. The jacobian\n"
       "is undefined outside the range covered by the grid. The grid must\n"
       "be sorted in increasing order, with no repetitions.\n"
       "\n"
       "Geocentric latitudes shall be used.\n"
       "\n"
       "The latitude grid can differ between retrieval quantities, and it\n"
       "is the grid defined when the retrieval quantity is added that is\n"
       "used.\n"
       "\n"
       "For 1D and 2D calculations this vector shall be set to be empty, but\n"
       "the number of latitudes shall be considered to be 1 when examining\n"
       "the size of variables.\n"
       "\n"
       "Allowed values for 3D, is [-360,360].\n"
       "\n"
       "Usage:      Set by the user.\n"
       "\n"
       "Unit:       degrees\n"
      ),
      GROUP( Vector_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "jacobian_p_grid" ),
      DESCRIPTION
      (
       "The pressure grid for the Jacobian matrix.\n"
       "\n"
       "The pressure surfaces on which the Jacobian matrix is determined.\n"
       "This variable has to be defined when calculating the the Jacobian.\n"
       "It must be sorted in decreasing order, with no repetitions.\n"
       "\n"
       "Usage:       Set by the user.\n"
       "\n"
       "Unit:        Pa\n"
      ),
      GROUP( Vector_ )));

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
       "Usage: Quantities are added by specific WSM;\n"
       "         jacobianAddGas,\n"
       "         jacobianAddTemp,\n"
       "         ...\n"
      ),
      GROUP( ArrayOfRetrievalQuantity_ )));

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
       "Unit:       degrees\n"
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
       "special variable is needed for 1D as there exists no latitude grid \n"
       "for such cases. The variable can be used, for example, to set the \n"
       "geoid radius or select atmospheric profiles from a 2D/3D \n"
       "climatology. \n"
       "\n"
       "For limb sounding, the choosen latitude should normally be selected\n"
       "to be valid for the tangent points, rather than for the sensor.\n"
       "\n"
       "Values shall be inside the range [-90,90].\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  degrees\n"
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "lo" ),
      DESCRIPTION
      (
       "The local oscillator frequency.\n"
       "\n"
       "The local oscillator frequency is used in a heterodyne system when\n"
       "the mixer folds the spectra from from radio frequencies (RF) to\n"
       "intermediate frequencies (IF).\n"
       "\n"
       "Unit: Hz\n"
       "\n"
       "Usage: Set by the user.\n"
       ),
      GROUP( Vector_ )));

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
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  degrees\n"
       ),
      GROUP( Vector_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "l_step" ),
      DESCRIPTION
      (
       "The pathlegth through one grid cell/layer.\n"
       "\n"
       "The pathlength is required in the methods for RT step calculations,\n"
       "which are *sto_vecGeneral* and *sto_vecScalar*.\n"
       "It can be calculated using the *ppath_step_agenda*.\n"
       "\n"
       "Usage:      Used in *doit_i_fieldUpdateXXX.\n"
       "\n"
       "Unit:       m \n"
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "matrix_1" ),
      DESCRIPTION
      (
       "An arbitrary matrix.\n"
       "\n"
       "This variable shall be treated as a general variable of type Matrix.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( Matrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "matrix_2" ),
      DESCRIPTION
      (
       "An arbitrary matrix.\n"
       "\n"
       "This variable shall be treated as a general variable of type Matrix.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( Matrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "main_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_)));

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
      GROUP( Vector_ )));

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
      GROUP( Vector_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_antenna" ),
       DESCRIPTION
       (
        "MCAntenna object used by MCGeneral to sample the field of view."
        "Possible antenna types include Pencil Beam, Gaussian (2D), and"
        "Antenna pattern lookup (yet to be implemented).\n"
        "\n"
        "Usage: Input to MCGeneral \n"
        ), 
       GROUP( MCAntenna_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_cloud_opt_path" ),
       DESCRIPTION
       (
        "The cloud optical path integrated over the field of view\n"
        "\n"
        "Usage: Output from mc_IWP_cloud_opt_pathCalc \n"
        ), 
       GROUP( Numeric_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_cloud_opt_path_error" ),
       DESCRIPTION
       (
        "standrad error in the cloud optical path integrated over the field of view\n"
        
        "\n"
        "Usage: Output from mc_IWP_cloud_opt_pathCalc \n"
        ), 
       GROUP( Numeric_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_error" ),
       DESCRIPTION
       (
        "Error in simulated *y* from cloudy-sky arts simulations using ScatteringMonteCarlo.\n"
        "\n"
        "Usage: Output from ScatteringMonteCarlo.. \n"
        "\n"
        "Units: W / (m^2 Hz sr)\n"
        "\n"
        "Size:  [ stokes_dim ]\n"
        ), 
       GROUP( Vector_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_incoming" ),
       DESCRIPTION
       (
        "Incoming radiance lookup table for Monte Carlo calculations\n"
        "\n"
        "SLIData2 object with x1=r (in m), x2=LOS zentiah angle, and\n"
        "y=incoming radiance in W / (m^2 Hz sr)]\n"
        "Usage: Input for ScatteringMonteCarlo.. \n"
        "\n"
        "Units: [m,degrees,W / (m^2 Hz sr)]\n"
        ), 
       GROUP( SLIData2_ )));


  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_iteration_count" ),
       DESCRIPTION
       (
        "Counts the number of iterations (or photons) used in the MC\n "
        "scattering algorithm.\n"
        "\n"
        "Usage: Set by ScatteringMonteCarlo.\n"
        ),
       GROUP( Index_)));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_IWP" ),
       DESCRIPTION
       (
        "The ice water path integrated over the field of view\n"
        "\n"
        "Usage: Output from mc_IWP_cloud_opt_pathCalc \n"
        ), 
       GROUP( Numeric_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_IWP_error" ),
       DESCRIPTION
       (
        "The standard error of ice water path integrated over the field of view\n"
        "\n"
        "Usage: Output from mc_IWP_cloud_opt_pathCalc \n"
        ), 
       GROUP( Numeric_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_points" ),
       DESCRIPTION
       (
        "Counts the number of MC endpoints in each grid cell\n"
        "\n"
        "Usage: Set by MCGeneral.\n"
        ),
       GROUP( Tensor3_)));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_seed" ),
       DESCRIPTION
       (
        "The integer seed for the random number generator used by\n"
        "ScatteringMonteCarlo and MCGeneral.\n"
        "\n"
        "Usage: Set by MCSetSeed.\n"
        ),
       GROUP( Index_)));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "mc_unit" ),
       DESCRIPTION
       (
        "Determines the unit used for the radiance returned by MCGeneral\n"
        "Possible values are \"RJBT\", which returns the Rayleigh Jeans\n" 
        "Brightness temperature in Kelvin, and \"radiance\", which returns\n"
        "the radiance in (Watts per meter squared per steradian)\n"
        "\n"
        "Usage: Set by the user.\n"
        ),
       GROUP( String_)));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "meridian_angle_1d" ),
      DESCRIPTION
      (
       "The meridian angle for a 1D observation.\n"
       "\n"
       "This variable is only used when calculating a curvature radius of \n"
       "the geoid for 1D cases. The given value shall be the angle between \n"
       "the observation and meridian plane where 0 degrees means an \n" 
       "observation in the N-S direction.\n"
       "\n"
       "Values shall be in the range [-180,180].\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit: degrees\n"
       ),
      GROUP( Numeric_ )));

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
      GROUP( Matrix_ )));
 
 wsv_data.push_back
   (WsvRecord
    ( NAME( "nelem" ),
      DESCRIPTION
      (
        "This variable is used by the VectorSet workspace method.\n"
       ),
      GROUP( Index_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "ncols" ),
      DESCRIPTION
      (
        "This variable is used by the MatrixSet, Tensor3Set, etc. \n"
        "workspace methods.\n"
       ),
      GROUP( Index_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "nrows" ),
      DESCRIPTION
      (
        "See *ncols*.\n"
       ),
      GROUP( Index_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "npages" ),
      DESCRIPTION
      (
        "See *ncols*.\n"
       ),
      GROUP( Index_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "nbooks" ),
      DESCRIPTION
      (
        "See *ncols*.\n"
       ),
      GROUP( Index_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "nshelves" ),
      DESCRIPTION
      (
        "See *ncols*.\n"
       ),
      GROUP( Index_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "nvitrines" ),
      DESCRIPTION
      (
        "See *ncols*.\n"
       ),
      GROUP( Index_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "nlibraries" ),
      DESCRIPTION
      (
        "See *ncols*.\n"
       ),
      GROUP( Index_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "numeric_1" ),
      DESCRIPTION
      (
      "An arbitrary numeric scalar.\n"
      "\n"
      "This variable shall be treated as a general variable of type Numeric.\n"
      "It can be used, for example, when some intermediate data must be\n"
      "generated or to copy some data.\n"
      "\n"
      "Usage: Set by user.\n"
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "numeric_2" ),
      DESCRIPTION
      (
      "An arbitrary numeric scalar.\n"
      "\n"
      "This variable shall be treated as a general variable of type Numeric.\n"
      "It can be used, for example, when some intermediate data must be\n"
      "generated or to copy some data.\n"
      "\n"
      "Usage: Set by user.\n"
       ),
      GROUP( Numeric_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "met_profile_calc_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_)));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "opt_prop_gas_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_)));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "opt_prop_part_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_)));

  wsv_data.push_back
    (WsvRecord
     (NAME( "output_file_format" ),
      DESCRIPTION
      (
       "Output file format. \n"
       "\n"
       "This variable sets the format for output files. It could be set to\n"
       "\"ascii\" or \"binary\".\n"
       "\n"
       "To change the value of this variable use the workspace methods\n"
       "*output_file_formatSetAscii* and *output_file_formatSetBinary*\n"
       ),
      GROUP( String_ )));

  wsv_data.push_back
    (WsvRecord
     (NAME( "particle_masses" ),
      DESCRIPTION
      (
       "Output file format.The mass of each particle type stored in a vector \n"
       "\n"
       "Usage: Set by the user\n"
       ),
      GROUP( Vector_ )));

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
      GROUP( Tensor4_ )));
   
   wsv_data.push_back
   (WsvRecord
    ( NAME( "pha_mat_spt" ),
      DESCRIPTION
      (
       "Phase matrix for a single particle type.\n"
       "\n"
       "This variable contains the elements of phase matrix for a single \n"
       "particle for given propagation direction. \n"
       "It is the calculated in the agenda *pha_mat_spt_agenda* .\n"
       "The elements of the phase matrix are calculated from   \n"
       "single_scattering_data*."
       "ARTS user guide (AUG) gives the formulas used for computing all \n"
       "elements of the phase matrix for a given particle type.\n"
       "\n"
       "Usage:      Input and Output of the method pha_mat_sptCalc\n"
       "\n"
       "Unit:        m^2\n"
       "\n"
       "Dimensions: [part_types, scat_za_grid, scat_aa_grid, stokes_dim, stokes_dim]\n"
       ),
      GROUP( Tensor5_ )));

    wsv_data.push_back
   (WsvRecord
    ( NAME( "pha_mat_spt_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_ ))); 

   wsv_data.push_back
   (WsvRecord
    ( NAME( "pha_mat_sptDOITOpt" ),
      DESCRIPTION
      (
       "Interpolated phase matrix.\n"
       "\n"
       "This variable contains the data of the phase matrix in the \n"
       "scattering frame interpolated on the actual frequency (the variable\n"
       "is used inside *scat_mono_agenda*) and also interpolated on all \n"
       "possible scattering angles following from all combinations of \n"
       "*scat_za_grid* and *scat_aa_grid*. \n"
       "\n"
       "Usage:      Input of the method *pha_mat_sptFromDataDOITOpt\n"
       "\n"
       "Unit:        m^2\n"
       "\n"
       "Dimensions: \n"
       "[particle types]\n"
       "[T, scat_za_grid,scat_aa_grid, scat_za_grid, scat_aa_grid,\n"
       "stokes_dim, stokes_dim]\n"
       ),
      GROUP( ArrayOfTensor7_ )));

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
       "\n"
       "Usage:      Calculated internally.\n"
       "\n"
       "Unit:        m^-3\n"
       "\n"
       "Size: [N_part_types, \n"
       "       (cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
       "       (cloudbox_limits[3] - cloudbox_limits[2]) +1, \n"
       "       (cloudbox_limits[5] - cloudbox_limits[4]) +1 ] \n"
        ),
      GROUP( Tensor4_ )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "pnd_field_perturb" ),
      DESCRIPTION
      (
       "The field representing particle number density perturbations.\n"
       "\n"
       "This variable gives the perturbation of particle number density\n"
       "of the chosen particle types as a function of p_grid, lat_grid,\n"
       "lon_grid for each retrieval quantity. The variable has to be\n"
       "prepared outside ARTS and it has to be setup prior to calling\n"
       "*jacobianAddParticle*. Since it is added to *pnd_field* during\n"
       "the calculation of the Jacobian, the perturbations are absolute\n"
       "and as such should have the same unit as *pnd_field*\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      Set by the user.\n"
       "\n"
       "Unit:       m^-3\n"
       "\n"
       "Dimensions: [N_retrieval_quantities, as *pnd_field* ]\n"
        ),
      GROUP( Tensor5_ )));

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
       "Usage:      Used in the methods *ParticleTypeAdd* and \n"
       "                  *ParticleTypeAddAll*\n"
       "\n"
       "Unit:        m^-3\n"
       "\n"
       "Size:  Array[N_pt]\n"
       "       GriddedField3 \n "
       "       [N_p] \n"
       "       [N_lat] \n"
       "       [N_lon] \n"
       "       [N_p, N_lat, N_lon] \n"
       ),
      GROUP( ArrayOfGriddedField3_ )));

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
       "The data struture is too extensive to be described here, but it is\n"
       "described carefully in the ARTS user guide (AUG). Use the index to\n"
       "find where the data structure, Ppath, for propagation paths is \n"
       "discussed. It is listed as a subentry to \"data structures\".\n"
       "\n"
       "Usage: Output from the method *ppathCalc*.\n"
       ),
      GROUP( Ppath_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath_array" ),
      DESCRIPTION
      (
       "The complete set of propagation paths.\n"
       "\n"
       "In the case of scattering (either inside the atmosphere or by the \n"
       "surface) the propagation path is divided into parts, with one part\n"
       "between each scattering event. This variable describes this complete\n"
       "set of propagation paths, in contrast to *ppath* that describes only\n"
       "one part.\n"
       "\n"
       "This variable is not always filled. It used as part of analytical \n"
       "jacobian calculations for gases and temperature. This variable can\n"
       "also be used for making plots. To force this variable to be filled,\n"
       "activate *ppath_array_do*.\n"
       "\n"
       "See the user guide for further details.\n"
       "\n"
       "Usage: See above.\n"
       ),
      GROUP( ArrayOfPpath_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath_array_do" ),
      DESCRIPTION
      (
       "Flag to fill *ppath_array*.\n"
       "\n"
       "Include FlagOn(ppath_array){} to fill *ppath_array* even if this is\n"
       "not needed for internal purposes.\n"
       "\n"
       "Note that this variable is set to 0 by *jacobianOff/Init*.\n"
       "\n"
       "Usage: Set by *RteCalc*.\n"
       ),
      GROUP( Index_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath_array_index" ),
      DESCRIPTION
      (
       "Index in *ppath_array* of present propagation path part.\n"
       "\n"
       "This variable shall point to correct position in *ppath_array* when\n"
       "*iy_calc* is called. A negative value means that there is no\n"
       "previous path part. The variable is then modified by *iy_calc* as\n"
       "soon as a new propagation path is calculated.\n"
       "\n"
       "Usage: Communication with *iy_calc*.\n"
       ),
      GROUP( Index_ )));

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
       "Usage:   In/output to/from *ppath_step_agenda*.\n"
       "\n"
       "Members: See AUG.\n"
       ),
      GROUP( Ppath_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath_step_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "ppath_transmissions" ),
      DESCRIPTION
      (
        "?\n"
       ),
      GROUP( Tensor4_ )));

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
      GROUP( Vector_ )));

  wsv_data.push_back
    (WsvRecord
    ( NAME( "refr_index" ),
      DESCRIPTION
      (
       "The total refractive index at some point in the atmosphere.\n"
       "\n"
       "This variable contains the refractive index summed over all relevant\n"
       "constituents, at one position in the atmosphere. The standard set of\n"
       "functions assumes that all frequency components propagate along the\n"
       "same path. That is, dispersion is neglected and this variable has\n"
       "frequency dimension.\n"
       "\n"
       "Unit: 1\n"
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "refr_index_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_ )));

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
      GROUP( ArrayOfIndex_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_do_t_jacs" ),
      DESCRIPTION
      (
       "Flag to *rte_agenda* to calculate jacobians for temperature.\n"
       "\n"
       "Usage:   Set internally, by *RteCalc*.\n"
      ),
      GROUP( Index_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_gp_p" ),
      DESCRIPTION
      (
       "The pressure grid position of *rte_pos*.\n"
       "\n"
       "This variable is used to give the grid position for an end point\n"
       "of the propagation path to some workspace method part of the\n"
       "radiative transfer calculations.\n"
       "\n"
       "Usage:   Set internally.\n"
       ),
      GROUP( GridPos_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_gp_lat" ),
      DESCRIPTION
      (
       "The latitude grid position of *rte_pos*.\n"
       "\n"
       "This variable is used to give the grid position for an end point\n"
       "of the propagation path to some workspace method part of the\n"
       "radiative transfer calculations.\n"
       "\n"
       "Usage:   Set internally.\n"
       ),
      GROUP( GridPos_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_gp_lon" ),
      DESCRIPTION
      (
       "The longitude grid position of *rte_pos*.\n"
       "\n"
       "This variable is used to give the grid position for an end point\n"
       "of the propagation path to some workspace method part of the\n"
       "radiative transfer calculations.\n"
       "\n"
       "Usage:   Set internally.\n"
       ),
      GROUP( GridPos_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_los" ),
      DESCRIPTION
      (
       "A line-of-sight for radiative transfer calculations.\n"
       "\n"
       "The main purpose of this WSV and *rte_pos* is communication with \n"
       "different agendas involved in the RTE calculations. These variables \n"
       "can also be used to enable calling of *ppathCalc* (and maybe other \n"
       "methods) from the workspace. \n"
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
      GROUP( Vector_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_planck_value" ),
      DESCRIPTION
      (
       "A Planck function value for radiative transfer calculations.\n"
       "\n"
       "The Planck function is used in the methods for RT step calculations,\n"
       "which are *sto_vecGeneral* and *sto_vecScalar*. \n"
       "\n"
       "Usage:      Calculated in *i_fieldUpdate1D*.\n"
       "\n"
       "Unit:       W / (m^2 Hz sr)\n "
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_pos" ),
      DESCRIPTION
      (
       "A geographical position for radiative transfer calculations.\n"
       "\n"
       "The main purpose of this WSV and *rte_los* is communication with \n"
       "different agendas involved in the RTE calculations. These variables \n"
       "can also be used to enable calling of *ppathCalc* (and maybe other \n"
       "methods) from the workspace. \n"
        "\n"
       "This variable is a vector with a length equalling the atmospheric\n"
       "dimensionality. The first element is the radius (from the coordinate\n"
       "system centre) of the position. Element 2 is the latitude and \n"
       "element 3 is the longitude. Please note that the vertical position \n"
       "is given as the radius, not the altitude above the geoid.\n"
       "\n"
       "Usage: See above. \n"
       "\n"
       "Units: [ m, degree, degree ]\n"
       "\n"
       "Size:  [ atmosphere_dim ]\n"
       ),
      GROUP( Vector_ )));

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
      GROUP( Numeric_ )));

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
      GROUP( Numeric_ )));

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
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "r_geoid" ),
      DESCRIPTION
      (
       "Geoid radius.\n"
       "\n"
       "Geometrical altitudes are defined as the vertical distance above the\n"
       "geoid, and the geoid is the reference surface used when giving, for\n"
       "example, *z_surface* and *z_field*. \n"
       "\n"
       "The geoid is defined by giving the radius from the coordinate centre\n"
       "to the geoid surface for each crossing of the latitude and longitude\n"
       "grids. The geoid should normally be selected to be an ellipsoid but\n"
       "any shape is allowed. For 1D calculations, the geoid is by\n"
       "definition a sphere.\n"
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
       "Dimensions: [ lat_grid, lon_grid ]\n"
       ),
      GROUP( Matrix_ )));

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
      GROUP( Vector_ )));

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
     GROUP( Index_ ))); 

   wsv_data.push_back
     (WsvRecord
      ( NAME( "scat_data_mono" ),
        DESCRIPTION
        (
         "Monochromatic single scattering data.\n"
         "\n"
         "This variable holds the single scattering properties for all \n"
         "hydrometeor species. It is calculated from scat_data_raw by \n"
         "*scat_data_monoCalc*, which interpolates scat_data_raw for the \n"
         "required frequency.\n"
         ),
        GROUP( ArrayOfSingleScatteringData_ ))); 

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
         "The unit of the single scattering properties is m^2.\n"
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
         "  Tensor6[pha_mat_data]\n"
         "      [f_grid, za_grid, aa_grid, za_grid, aa_grid, matrix_element]\n"
         "  Tensor4[ext_mat_data]\n"
         "      [f_grid, za_grid, aa_grid, matrix_element]\n"
         "  Tensor4[abs_vec_data]\n"
         "      [f_grid, za_grid, aa_grid, matrix_element]\n"
         ),
        GROUP( ArrayOfSingleScatteringData_ ))); 
   
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
       "latitude surfaces of the cloudbox boundary. It contains all four \n"
       "components of the Stokes vector.\n"
       "\n"
       "This variable is used as interface between the clear sky and the \n"
       "scattering calculations. \n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      In/Output from/to *ScatteringMain* \n"
       "\n"
       "Unit:        W / (m^2 Hz sr) \n"
       "\n"
       "Dimensions: [ f_grid, p_grid, latitude surface, lon_grid, \n"
       "              scat_za_grid \n  scat_aa_grid, stokes_dim ]\n"
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
       "Usage:      Output from *ScatteringMain* \n"
       "\n"
       "Unit:        W / (m^2 Hz sr) \n"
       "\n"
       "Dimensions: [ f_grid, p_grid, lat_grid, latitude surface, \n"
       "              scat_za_grid, scat_aa_grid, stokes_dim]\n"
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
       "latitude surfaces of the cloudbox boundary. It contains all four \n"
       "components of the Stokes vector.\n"
       "\n"
       "This variable is used as interface between the clear sky and the \n"
       "scattering calculations. \n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage:      In/Output from *ScatteringMain* \n"
       "\n"
       "Unit:        W / (m^2 Hz sr) \n"
       "\n"
       "Dimensions: [ f_grid, pressure surfaces, lat_grid, lon_grid, \n" 
       "              scat_za_grid, scat_aa_grid, stokes_dim]\n"
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
       "properties of particles like *ext_mat_partCalc* and *pha_matCalc*.\n"
       "It holds the information about the position for which the \n"
       "scattering calculations are done. \n"
       "\n"
       "Usage:    Input to the methods *spt_calc_agenda*,\n"
       "                               *pha_mat_spt_agenda*\n"
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
       "properties of particles like *ext_mat_partCalc* and *pha_matCalc*.\n"
       "It holds the information about the position for which the \n"
       "scattering calculations are done.  \n"
       "\n"
       "Usage:    Input to the methods *spt_calc_agenda*,\n"
       "                               *pha_mat_spt_agenda*\n"
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
       "properties of particles like *ext_mat_partCalc* and *pha_matCalc*.\n"
       "It holds the information about the location for which the \n"
       "scattering calculations are done.\n"  
       "\n"
       "Usage:    Input to the methods *spt_calc_agenda*,\n"
       "                               *pha_mat_spt_agenda*\n"
       ),
     GROUP( Index_ ))); 
  
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
      GROUP( Vector_ )));

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
      GROUP( Index_ )));

 
 
  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_los" ),
      DESCRIPTION
      (
       "The sensor line-of-sight for each measurement block.\n"
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
      GROUP( Matrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_norm" ),
      DESCRIPTION
      (
       "Flag if sensor response should be normalised (0 or 1).\n"
       "\n"
       "If the flag is set to 1 each sensor response block will be\n"
       "normalised. For a flag equal 0 the sensor response is left as is.\n"
       "The polarisation response is the exception, it is not affected by\n"
       "this flag.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       ),
      GROUP( Index_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_pol" ),
      DESCRIPTION
      (
       "The sensor polarisation response.\n"
       "\n"
       "The number of rows of this matrix equals the number of polarisation\n"
       "channels of the sensor. For example, if horisontal and vertical \n"
       "polarisation are measured simultaneously, the matrix has two rows. \n"
       "This matrix is multiplicated with the Stokes vector, converted to \n"
       "the sensor frame, before the sensor response is applied. \n"
       "\n"
       "If the purpose of the simulations is to extract the polarisation\n"
       "of the radiation coming from the atmosphere (no sensor), the matrix\n"
       "shall be set to the identity matrix with a size matching\n"
       "*stokes_dim*. \n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  [ - (0-1) ]\n"
       "\n"
       "Size:  [ number of polarisation values, stokes_dim ]\n"
       ),
      GROUP( Matrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_pos" ),
      DESCRIPTION
      (
       "The sensor position for each measurement block.\n"
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
      GROUP( Matrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response" ),
      DESCRIPTION
      (
        "The response block matrix modelling the total sensor response.\n"
        "\n"
        "The matrix is the product of all the individual sensor response\n"
        "matrices. Therefore its dimension are depending on the sensor\n"
        "configuration and where in the calculations we are.\n"
        "The *sensor_response* has to initialised by the \n"
        "*sensor_responseInit* method.\n"
        "\n"
        "Usage:   Output/input to the *sensor_response...* methods.\n"
        "\n"
        "Units:   -\n"
        "\n"
        "Dimension:     See the individual *sensor_response...* method \n"
        "               documentation.\n"
       ),
      GROUP( Sparse_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_f" ),
      DESCRIPTION
      (
       "The frequencies associated with *sensor_response*.\n"
       "\n"
       "This vector gives the output frequencies for the sensor parts\n"
       "considered in *sensor_response*.\n"
       "\n"
       "The variable shall not be set manually, it will be set together with\n"
       "*sensor_response* by sensor response WSMs.\n"
       "\n"
       "Usage: Set by sensor response methods.\n"
       "\n"
       "Unit:  [ Hz ]\n"
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_za" ),
      DESCRIPTION
      (
       "The relative zenith angles associated with *sensor_response*.\n"
       "\n"
       "This vector summed with the zenith angles in *sensor_los* gives the\n"
       "observation directions observed. For example, if the measurement \n"
       "blocks consist of one spectrum, this vector has length 1, and most \n"
       "likely with the value 0.\n"
       "\n"
       "The variable shall not be set manually, it will be set together with\n"
       "*sensor_response* by sensor response WSMs.\n"
       "\n"
       "Usage: Set by sensor response methods.\n"
       "\n"
       "Unit:  [ degrees ]\n"
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_aa" ),
      DESCRIPTION
      (
       "The azimuth angles associated with *sensor_response*.\n"
       "\n"
       "The relative azimuth angles associated with *sensor_response*.\n"
       "The variable shall be empty if *antenna_dim* = 1.\n"
       "\n"
       "This vector summed with the azimuth angles in *sensor_los* gives the\n"
       "observation directions observed. For example, if the measurement \n"
       "blocks consist of one spectrum, this vector has length 1, and most \n"
       "likely with the value 0.\n"
       "\n"
       "The variable shall not be set manually, it will be set together with\n"
       "*sensor_response* by sensor response WSMs.\n"
       "\n"
       "Usage: Set by sensor response methods.\n"
       "\n"
       "Unit:  [ degrees ]\n"
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_response_pol" ),
      DESCRIPTION
      (
       "The number of polarisations measured by the sensor.\n"
       "\n"
       "This variable is set to be equal to *stokes_dim* when the sensor\n"
       "is initialised. If a polarisation response is calculated the\n"
       "variable is updated to match the output polarisations from\n"
       "the sensor response.\n"
       "\n"
       "The variable shall not be set manually, it will be set together with\n"
       "*sensor_response* by sensor response WSMs.\n"
       "\n"
       "Usage: Set by sensor response methods.\n"
       ),
      GROUP( Index_ )));

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
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_time" ),
      DESCRIPTION
      (
       "The time for each measurement block.\n"
       "\n"
       "This WSV is used when calculating pointing offset Jacobians.\n"
       "It can either be used to store actual times, in any desired unit,\n"
       "for real measurements or it can store a relative time series\n"
       "for the measurement blocks.\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  [ arbitrary ]\n"
       "\n"
       "Size:  [ number of measurement blocks ]\n"
       ),
      GROUP( Vector_ )));

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
       "      [f_grid, T_grid, za_grid, aa_grid, za_grid, aa_grid,"
       "matrix_element]\n"
       "  Tensor5[ext_mat_data]\n"
       "      [f_grid, T_grid, za_grid, aa_grid, matrix_element]\n"
       "  Tensor5[abs_vec_data]\n"
       "      [f_grid, T_grid, za_grid, aa_grid, matrix_element]\n"
       ),
      GROUP( SingleScatteringData_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sparse_1" ),
      DESCRIPTION
      (
       "An arbitrary sparse matrix.\n"
       "\n"
       "This variable shall be treated as a general variable of type Sparse.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( Sparse_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "species_index" ),
      DESCRIPTION
      (
       "This ArrayOfIndex yields the tag positions of key species like\n"
       "N2 (=0), O2 (=1), H2O (=2), O3 (=3), CO2 (=4).\n"
       "For example species_index[2] gives the first H2Otag position in the\n"
       "controle file specified list of tags for which calculations should\n"
       "be performed.\n"
       ),
      GROUP( ArrayOfIndex_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "stokes_dim" ),
      DESCRIPTION
      (
       "The dimensionality of the Stokes vector (1-4).\n"
       "\n"
       "Usage:      Set by the user.\n"
       ),
      GROUP( Index_ )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "spt_calc_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc.\n"
       ),
      GROUP( Agenda_ )));

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
       "Usage:      Calculated internally.\n"
       ),
      GROUP( Vector_ )));
 
   wsv_data.push_back
     (WsvRecord
      ( NAME( "surface_emission" ),
        DESCRIPTION
        ( "The emission from the surface at a specified position.\n"
          "\n"
      "The position is normally specified by *rte_pos* or the combination of\n"
      "*rte_gp_p*, *rte_gp_lat* and *rte_gp_lon*.\n"
          "\n"
          "See further *surfaceCalc* and the user guide.\n"
          "\n"
      "Usage:      Input to methods for *iy_surface_agenda*."
          "\n"
          "Unit:       W / (m^2 Hz sr)\n"
          "\n"
          "Dimensions: [ f_grid, stokes_dim ]\n"
         ), 
        GROUP( Matrix_ )));

   wsv_data.push_back
     (WsvRecord
      ( NAME( "surface_emissivity" ),
        DESCRIPTION
        ( "The surface emissivity at position of interest.\n"
          "\n"
          "Usage: Input to surfaceSingleEmissivity.\n"
          "\n"
          "Unit: a value between 0 and 1\n"
         ), 
        GROUP( Numeric_ )));

   wsv_data.push_back
     (WsvRecord
      ( NAME( "surface_emissivity_field" ),
        DESCRIPTION
        ( "The surface emissivity specified on lat_grid and lon_grid.\n"
          "\n"
          "Dimensions: [ lat_grid, lon_grid ]\n"
         ), 
        GROUP( Matrix_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "surface_los" ),
       DESCRIPTION
       (
        "Directions for which to calculate downwelling radiation when \n"
        "considering a surface reflection.\n"
        "\n"
        "See further the user guide.\n"
        "\n"
        "Units: degrees\n"
        "\n"
        "Size:  [ any number, 1 or 2 ]\n"
        ), 
       GROUP( Matrix_ )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "surface_prop_agenda" ),
       DESCRIPTION
       (
        "See agendas.cc.\n"
        ),
       GROUP(  Agenda_ )));
  
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
        "downwelling radiation. For example, if the surface has isotropic\n"
        "scattering, without absorbtion of incoming radiation, and the\n"
        "downwelling radiation is calculated at ten angles (i.e. the length\n"
        "of *surface_los* is ten), the surface reflection coefficients \n"
        "are 0.1.\n"
        "\n"
        "See further *surfaceCalc* and the user guide.\n"
        "\n"
        "Usage:      Input to methods for *iy_surface_agenda*."
        "\n"
        "Units:      -\n"
        "\n"
        "Dimensions: [ surface_los, f_grid, stokes_dim, stokes_dim ]\n"
        ), 
       GROUP( Tensor4_ )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "surface_skin_t" ),
      DESCRIPTION
      (
       "Surface skin temperature.\n"
       "\n"
       "Usage:   Input to methods for *iy_surface_agenda*.\n"
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "tensor3_1" ),
      DESCRIPTION
      (
       "An arbitrary Tensor3.\n"
       "\n"
       "This variable is a general variable of type Tensor3.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( Tensor3_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "tensor3_2" ),
      DESCRIPTION
      (
       "An arbitrary Tensor3.\n"
       "\n"
       "This variable is a general variable of type Tensor3.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( Tensor3_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "tensor4_1" ),
      DESCRIPTION
      (
       "An arbitrary Tensor4.\n"
       "\n"
       "This variable is a general variable of type Tensor4.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( Tensor4_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "tensor4_2" ),
      DESCRIPTION
      (
       "An arbitrary Tensor4.\n"
       "\n"
       "This variable is a general variable of type Tensor4.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( Tensor4_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "tensor5_1" ),
      DESCRIPTION
      (
       "An arbitrary Tensor5.\n"
       "\n"
       "This variable is a general variable of type Tensor5.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( Tensor5_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "tensor5_2" ),
      DESCRIPTION
      (
       "An arbitrary Tensor4.\n"
       "\n"
       "This variable is a general variable of type Tensor5.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( Tensor5_ )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "tensor6_1" ), 
       DESCRIPTION
      (
       "An arbitrary tensor6. \n"
       "\n"
       "This variable is a general variable of type Tensor6.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
       GROUP( Tensor6_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "timer" ),
      DESCRIPTION
      (
       "Stores the starting time for time measurements.\n"
       ),
      GROUP( Timer_ )));

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
      GROUP( Tensor3_ )));

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
       "Usage:      Set by the user by choosing a climatology.\n"
       "\n"
       "Unit:       K\n"
       "\n"
       "Size   GriddedField3 \n "
       "       [N_p] \n"
       "       [N_lat] \n"
       "       [N_lon] \n"
       "       [N_p, N_lat, N_lon] \n"
       ),
      GROUP( GriddedField3_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "vector_1" ),
      DESCRIPTION
      (
       "An arbitrary vector.\n"
       "\n"
       "This variable shall be treated as a general variable of type Vector.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "vector_2" ),
      DESCRIPTION
      (
       "An arbitrary vector.\n"
       "\n"
       "This variable shall be treated as a general variable of type Vector.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "vector_3" ),
      DESCRIPTION
      (
       "An arbitrary vector.\n"
       "\n"
       "This variable shall be treated as a general variable of type Vector.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user.\n"
       ),
      GROUP( Vector_ )));

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
      GROUP( Tensor4_ ))); 

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
       "Usage:      Output of *AtmRawRead*\n"
       "            Input to *AtmFieldsCalc*.\n"
       "\n"
       "Unit:        absolute number\n"
       "\n"
       "Size:  Array[N_pt]\n"
       "       GriddedField3 \n "
       "       [N_p] \n"
       "       [N_lat] \n"
       "       [N_lon] \n"
       "       [N_p, N_lat, N_lon] \n"
       ),
      GROUP( ArrayOfGriddedField3_ )));


  wsv_data.push_back
   (WsvRecord
    ( NAME( "xml_output_type" ),
      DESCRIPTION
      (
       "Flag to determine whether XML output is binary or ascii\n"
       "\n"
       "This flag has to be set using the workspace method *output_file_formatSetAscii*\n"
       "or *output_file_formatSetBinary*. One of these methods MUST be called before\n"
       "writing the first output file."
       "\n"
       "Usage:      Set by user.\n"
       ),
      GROUP( Index_ )));

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
      GROUP( Vector_ )));

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
       ),
      GROUP( Matrix_ )));

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
       "       *bach_update_agenda*.\n"
       ),
      GROUP( Index_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "ybatch_n" ),
      DESCRIPTION
      (
       "Number of batch cases defined.\n"
       "\n"
       "See further *ybatchCalc*.\n"
       "\n"
       "Usage: Output from *batch_pre_agenda*.\n"
       ),
      GROUP( Index_ )));

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
       "in increasing order, with no repetitions. Otherwise the altitudes\n"
       "can be set to arbitrary values. Hydrostatic equilibrium is not\n"
       "applied automatically. If hydrostatic equilibrium applies, *z_field*\n"
       "must be set by a method ensuring that this criterium is fulfilled.\n"
       "\n"
       "The radius (from the coordinate centre) for a point between the grid\n"
       "crossings is obtained by a (multi-)linear interpolation of the sum\n"
       "of *r_geoid* and *z_field*.\n" 
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
      GROUP( Tensor3_ )));

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
       "Usage:      Set by the user by choosing a climatology.\n"
       "\n"
       "Unit:       K\n"
       "\n"
       "Size   GriddedField3 \n "
       "       [N_p] \n"
       "       [N_lat] \n"
       "       [N_lon] \n"
       "       [N_p, N_lat, N_lon] \n"
       ),
      GROUP( GriddedField3_ )));

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
       "interpolation of the sum of *r_geoid* and *z_surface*. With other \n"
       "words, the radius for the surface is assumed to vary linear along \n"
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
      GROUP( Matrix_ )));
 

  //--------------------------------------------------------------------------------
  // Zeeman WSVs, commented out after back-porting absorption from arts-1-0.

//  wsv_data.push_back
//     (WsvRecord
//     ( NAME( "abs_vec_zee" ),
//       DESCRIPTION
//       (
//        "Zeeman absorption vector.\n"
//        "\n"
//        "This variable contains the total absorption coefficient vector \n"
//        "used in the RTE calculation for user specified O2 (currently) lines, \n"
//        "affected by the Zeeman effect. The physical absorption includes  \n"
//        "O2 absorption from all 3 polarization components of the absorbed radiation. \n"
//        "\n"
//        "The vector is calculated by the agenda *zeeman_prop_agenda* \n"
//        "\n"
//        "The dimensision of the variable adapts to *stokes_dim*.\n"
//        "\n"
//        "See further the ARTS user guide (AUG). Use the index to find where\n"
//        "this variable is discussed. The variable is listed as a subentry to\n"
//        "\"workspace variables\".\n"
//        "\n"
//        "Usage:      Output of the agenda *zeeman_prop_agenda* \n"
//        "\n"
//        "Unit:        [Hz, m^2]\n"
//        "\n"
//        "Dimensions: [f_grid, stokes_dim]\n"
//         ),
//        GROUP( Matrix_ )));

// wsv_data.push_back
//     (WsvRecord
//      ( NAME( "ext_mat_zee" ),
//        DESCRIPTION
//       (
//        "Zeeman extinction matrix.\n"
//        "\n"
//        "This variable contains the total extinction matrix used \n"
//        "in the RTE calculation for user specified O2 (currently) lines, \n"
//        "affected by the Zeeman effect. It is the physical extinction matrix  \n"
//        "of all 3 polarization components of the radiation. \n"
//        "\n"
//        "Usage:      Output of the agendas *zeeman_prop_agenda* \n"
//        "\n"
//        "Unit:       [Hz, m^2, m^2] "
//        "\n"
//        "Dimensions: [f_grid, stokes_dim, stokes_dim]\n"
//        ),
//        GROUP( Tensor3_ )));

//  wsv_data.push_back
//    (WsvRecord
//     ( NAME( "zeeman_prop_agenda" ),
//       DESCRIPTION
//       (
//         "See agendas.cc.\n"
//        ),
//       GROUP( Agenda_ )));

//   wsv_data.push_back
//     (WsvRecord
//      (NAME( "zeeman_o2_onoff" ),
//       DESCRIPTION
//       (
//        "Make the Zeeman specific settings for O2 Zeeman spectral line\n"
//        "splitting for the microwave range (1-1000 GHz).\n"
//        "If zeeman_o2_onoff=1 the Zeeman effect is considered,\n"
//        "and if zeeman_o2_onoff=0 the Zeeman effect is omitted.\n"
//        ),
//       GROUP( Index_ )));
 
//   wsv_data.push_back
//     (WsvRecord
//      (NAME( "zeeman_o2_pressure_limit" ),
//       DESCRIPTION
//       (
//        "Make the Zeeman specific settings for O2 Zeeman spectral line\n"
//        "splitting for the microwave range (1-1000 GHz).\n"
//        "This variable sets the upper pressure limit [Pa] at which the\n"
//        " Zeeman splitting is taken into account\n"
//        ),
//       GROUP( Numeric_ )));

//   wsv_data.push_back
//     (WsvRecord
//      (NAME( "zeeman_o2_line" ),
//       DESCRIPTION
//       (
//        "Calculate this line with Zeeman splitting, lines are\n"
//        "identified by their upper rotational angular momentum quantum\n"
//        "number N, postive values of zeeman_o2_line are N+ transitions,\n"
//        "negative are N-.\n"
//        ),
//       GROUP( Index_ )));






}
