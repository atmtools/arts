/* Copyright (C) 2000, 2001, 2002
   Stefan Buehler <sbuehler@uni-bremen.de>
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
  \file   methods.cc
  \brief  Definition of method description data.

  This file contains only the definition of the function
  define_md_data, which sets the WSV lookup data. You have to change
  this function each time you add a new method. See methods.h for more
  documentation.

  \author Stefan Buehler
  \date 2000-06-10 */

#include "arts.h"
#include "make_array.h"
#include "auto_wsv.h"
#include "methods.h"
#include "auto_wsv_groups.h"

// Some #defines and typedefs to make the records better readable:
#define NAME(x) x
#define DESCRIPTION(x) x
#define OUTPUT   MakeArray<Index>
#define INPUT    MakeArray<Index>
#define GOUTPUT  MakeArray<Index>
#define GINPUT   MakeArray<Index>
#define KEYWORDS MakeArray<String>
#define TYPES    MakeArray<TokValType>
#define AGENDAMETHOD(x) x
#define SUPPRESSHEADER(x) x


/* Here's a template record entry:  (PE 2001-09-18)

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "FunctionName" ),
        DESCRIPTION
        (
         "A summary of the function in one sentence.\n"
         "\n"
         "A detailed description of the function. Please, try to be as \n"
         "clear and detailed as possible, this will help both you and \n"
         "others in the long run. \n"
         "\n"
         "Paragraphs are seperated with blank lines.\n"
         "\n"
         "The names of workspace variables and other methods are marked by\n"
         "stars, for example *z_plat*.\n"
         "\n"
         "Generic input and output, and keywords shall be described \n"
         "as exemplified below. If there is no variables of a group, \n"
         "(e.g. generic input) remove that part totally. Note that the \n"
         "on-line help just gives the type of generic input/output and the \n"
         "keyword names, and additional information is for sure needed.\n"
         "\n"
         "Leave space and brake lines when listing input and output \n"
         "variabales to make the code easier to read. See example below. \n"
         "\n"
         "Generic input: \n"
         "   Vector : Vector giving some very important input. Don't \n"
         "            be too short. Use the type of indention used here. \n"
         "\n"
         "Generic output: \n"
         "   Vector : Return vector for the zenith angles. The normal \n"
         "            options are ZA_PENCIL and ZA_SENSOR. \n"
         "\n"
         "Keywords:\n"
         "   delta_t   : Time increment between observations.\n"
         "   z_tan_lim : Vector with start and stop tangent altitudes."
        ),
        OUTPUT(),
        INPUT( z_plat_, p_abs_, z_abs_, l_step_, refr_, refr_lfac_,
               refr_index_, r_geoid_, z_ground_ ),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "delta_t", "z_tan_lim" ),
        TYPES(    Numeric_t, Vector_t    )));
  */

  /* Here's an empty record entry:  (PE 2001-09-18)

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "" ),
        DESCRIPTION
        (
         "\n"
         "\n"
         "Generic input: \n"
         "   \n"
         "\n"
         "Generic output: \n"
         "   \n"
         "\n"
         "Keywords:\n"
         "   "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));
  */




void define_md_data_raw()
{
  // The variable md_data is defined in file methods_aux.cc.
  extern Array<MdRecord> md_data_raw;


  // Initialize to zero, just in case:
  md_data_raw.resize(0);


  /////////////////////////////////////////////////////////////////////////////
  // Let's put in the functions in alphabetical order. This gives a clear rule
  // for where to place a new function and this gives a nicer results when
  // the functions are listed by "arts -m all".
  // No distinction is made between uppercase and lowercase letters. The sign
  // "_" comes after all letters.
  // Patrick Eriksson 2002-05-08
  /////////////////////////////////////////////////////////////////////////////


 //  md_data_raw.push_back
//     ( MdRecord
//       ( NAME("abs_vec_gasExample"),
//      DESCRIPTION
//      (
//       "This is only an example method created to perform test \n"
//       "calculations.\n"
//       ),
//      OUTPUT(abs_vec_),
//      INPUT(p_grid_,  atmosphere_dim_, stokes_dim_, scat_p_index_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("absMPM02_H2O"),
        DESCRIPTION
        (
         "Calculate water vapor absorption\n"
         "\n"
         "Corrected version of MPM93 by TKS, iup, 2002 The H2O lines at\n"
         "547.676440 GHz and 552.020960 GHz are isotopic lines: 547 GHz is from\n"
         "the isotope 1-8-1 (HITRAN code 181, JPL code 20003) with an isotopic\n"
         "ratio of 0.00199983 and 552 GHz is from the isotope 1-7-1 (HITRAN code\n"
         "171, JPL code 19003) with an isotopic ratio of 0.00037200.\n"
         "\n"
         "The original source code of MPM93 has these isotopic ratios not\n"
         "included in the line strength parameter b1, which is an error. In the\n"
         "arts implementation the line strength parameter b1 of these two lines\n"
         "is multiplied with the appropriate isotopic ratio.\n"
         "\n"
         "Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton, Propagation\n"
         "modeling of moist air and suspended water/ice particles at frequencies\n"
         "below 1000 GHz, AGARD 52nd Specialists Meeting of the Electromagnetic\n"
         "Wave Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21.\n"
         "\n"
         "Output:\n"
         "   abs    : absorption coefficients [1/m], \n"
         "            dimension: [ f_grid, abs_p (=abs_t) ]\n"
         "\n"
         "Input:\n"
         "   f_grid : Frequency grid [Hz].\n"
         "   abs_p  : List of pressures [Pa].\n"
         "   abs_t  : List of temperatures [K]. Must have same length as abs_p!\n"
         "   abs_vmr   : List of volume mixing ratios [absolute number].\n"
         "               Must have same length as abs_p!\n"
         "   abs_model : String specifying the model to use. Allowed options:\n"
         "               \"MPM02\"          - Calculate lines and continuum.\n"
         "               \"MPM02Lines\"     - Calculate only lines.\n"
         "               \"MPM02Continuum\" - Calculate only continuum.\n"
         "               \"user\"           - Use parameters given by abs_user_parameters,\n"
         "                                    instead of the predefined settings.\n"
         "   abs_user_parameters : Only used if abs_model==\"user\". In that case,\n"
         "                         abs_user_parameters must have 3 elements:\n"
         "                         1. Continuum scaling factor\n"
         "                         2. Line strength scaling factor\n"
         "                         3. Line broadening scaling factor\n"
         "                         Setting all scaling factors to 1 gives the\n"
         "                         same behavior as abs_model==\"MPM02\".\n"
         "                         Must be empty if one of the predefined models is used."
        ),
        OUTPUT( abs_ ),
        INPUT(  f_grid_, abs_p_, abs_t_, abs_vmr_,
                abs_model_, abs_user_parameters_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));


  md_data_raw.push_back     
    ( MdRecord
      ( NAME("absCKDMT_H2O_H2O"),
        DESCRIPTION
        (
         "Calculate water vapor self continuum absorption\n"
         "\n"
         "CKD version MT 1.00 self continuum absorption coefficient.\n"
         "The original  code is taken from the FORTRAN77 code of\n"
	 "Atmospheric and Environmental Research Inc. (AER),\n"
         "    Radiation and Climate Group\n"
         "    131 Hartwell Avenue\n"
         "    Lexington, MA 02421, USA\n"
         "    http://www.rtweb.aer.com/continuum_frame.html\n"
         "\n"
         "Output:\n"
         "   abs    : absorption coefficients [1/m], \n"
         "            dimension: [ f_grid, abs_p (=abs_t) ]\n"
         "\n"
         "Input:\n"
         "   f_grid              : Frequency grid [Hz].\n"
         "   abs_p               : List of pressures [Pa].\n"
         "   abs_t               : List of temperatures [K].\n"
         "                         (Must have same length as abs_p!)\n"
         "   abs_vmr             : List of H2O volume mixing ratios [absolute number].\n"
         "                         Must have same length as abs_p!\n"
         "   abs_model           : String specifying the model to use.\n"
         "                         Allowed options are:\n"
         "                         \"CKDMT100\": calculate self continuum absorption\n"
         "                                       according to CKD_MT version 1.00\n"
         "                         \"user\"    : use user defined parameter given by\n"
         "                                       abs_user_parameters, instead of \n"
         "                                       the predefined settings.\n"
         "   abs_user_parameters : Only used if abs_model==\"user\". In that case,\n"
         "                         abs_user_parameters must have 1 element:\n"
         "                         1. Continuum scaling factor\n"
         "                         abs_user_parameters must be empty if one of the\n"
         "                         pre-defined models is used."
        ),
        OUTPUT( abs_ ),
        INPUT(  f_grid_, abs_p_, abs_t_, abs_vmr_,
                abs_model_, abs_user_parameters_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("absCKDMT_H2O_AIR"),
        DESCRIPTION
        (
         "Calculate water vapor foreign continuum absorption\n"
         "\n"
         "CKD version MT 1.00 foreign continuum absorption coefficient.\n"
         "The original  code is taken from the FORTRAN77 code of\n"
	 "Atmospheric and Environmental Research Inc. (AER),\n"
         "    Radiation and Climate Group\n"
         "    131 Hartwell Avenue\n"
         "    Lexington, MA 02421, USA\n"
         "    http://www.rtweb.aer.com/continuum_frame.html\n"
         "\n"
         "Output:\n"
         "   abs    : absorption coefficients [1/m], \n"
         "            dimension: [ f_grid, abs_p (=abs_t) ]\n"
         "\n"
         "Input:\n"
         "   f_grid              : Frequency grid [Hz].\n"
         "   abs_p               : List of pressures [Pa].\n"
         "   abs_t               : List of temperatures [K].\n"
         "                         (Must have same length as abs_p!)\n"
         "   abs_vmr             : List of H2O volume mixing ratios [absolute number].\n"
         "                         Must have same length as abs_p!\n"
         "   abs_model           : String specifying the model to use.\n"
         "                         Allowed options are:\n"
         "                         \"CKDMT100\": calculate self continuum absorption\n"
         "                                       according to CKD_MT version 1.00\n"
         "                         \"user\"    : use user defined parameter given by\n"
         "                                       abs_user_parameters, instead of \n"
         "                                       the predefined settings.\n"
         "   abs_user_parameters : Only used if abs_model==\"user\". In that case,\n"
         "                         abs_user_parameters must have 1 element:\n"
         "                         1. Continuum scaling factor\n"
         "                         abs_user_parameters must be empty if one of the\n"
         "                         pre-defined models is used."
        ),
        OUTPUT( abs_ ),
        INPUT(  f_grid_, abs_p_, abs_t_, abs_vmr_,
                abs_model_, abs_user_parameters_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));


  md_data_raw.push_back     
    ( MdRecord
      ( NAME("absO2Model"),
        DESCRIPTION
        (
         "Calculate oxygen absorption in the 1-1000GHz range from  the absorption"
	 " model based on P.W.Rosenkranz and H. J. Liebe (MPM).\n"
         "Output:\n"
         "   abs    : absorption coefficients [1/m], \n"
         "            dimension: [ f_grid, abs_p (=abs_t) ]\n"
         "\n"
         "Input:\n"
         "   geomag_los          : magnetic filed strength and angle with respect to the LOS\n"
         "   f_grid              : Frequency grid [Hz].\n"
         "   abs_p               : List of pressures [Pa].\n"
         "   abs_t               : List of temperatures [K].\n"
         "                         (Must have same length as abs_p!)\n"
         "   abs_vmr             : List of H2O volume mixing ratios [absolute number].\n"
         "                         Must have same length as abs_p!\n"
         "   abs_model           : String specifying the model to use.\n"
         "                         Allowed options are:\n"
         "   abs_user_parameters : Only used if abs_model==\"user\" or \"*Scaling\". \n"
         "                         In that case, abs_user_parameters must have \n"
         "                         4 or 1 element(s) respectively.\n"
	 "                         abs_model==\"user\": \n"
         "                         1. O2 Continuum scaling factor\n"
         "                         2. O2 line strength scaling factor\n"
         "                         3. O2 line pressure broadening scaling factor\n"
         "                         4. O2 line mixing scaling factor\n"
	 "                         abs_model==\"*Scaling\": \n"
         "                         1. O2 absorption overall scaling factor\n"
         "                         Note:\n"
         "                         abs_user_parameters must be empty if one of the\n"
         "                         pre-defined models is used."
        ),
        OUTPUT( abs_ ),
        INPUT(  geomag_los_, f_grid_, abs_p_, abs_t_, abs_vmr_,
                abs_model_, abs_user_parameters_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));
  


  md_data_raw.push_back     
    ( MdRecord
      ( NAME("abs_scalar_gasExtractFromLookup"),
        DESCRIPTION
        (
         "Extract scalar gas absorption coefficients from lookup table.\n"
         "\n"
         "This extracts the absorption coefficient for all species in the\n"
         "current calculation from the lookup table. Extraction is for one\n"
         "specific atmospheric condition, i.e., a set of pressure, temperature,\n"
         "and VMR values.\n"
         "\n"
         "Extraction can be either for a single frequency (f_index>=0), or for\n"
         "all frequencies (f_index<0). The dimension of the output\n"
         "abs_scalar_gas is adjusted accordingly."
        ),
        OUTPUT( abs_scalar_gas_ ),
        INPUT(  gas_abs_lookup_, gas_abs_lookup_is_adapted_,
                f_index_, 
                rte_pressure_, rte_temperature_, rte_vmr_list_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("abs_scalar_gas_fieldCalc"),
        DESCRIPTION
        (
         "Calculate scalar gas absorption for all points in the atmosphere.\n"
         "\n"
         "This is useful in two different contexts:\n"
         "\n"
         "1. For testing and plotting gas absorption. (For RT calculations, gas\n"
         "absorption is calculated or extracted locally, therefore there is no\n"
         "need to calculate a global field. But this method is handy for easy\n"
         "plotting of absorption vs. pressure, for example.)\n"
         "\n"
         "2. Inside the scattering region, monochromatic absorption is\n"
         "pre-calculated for the entire atmospheric field.\n"
         "\n"
         "Because of the different contexts, the method can calculate absorption\n"
         "either for all frequencies in the frequency grid (f_index<0), or just\n"
         "for the frequency indicated by f_index (f_index>=0).\n"
         "\n"
         "The calculation itself is performed by the\n"
         "*scalar_gas_absorption_agenda*, which needs the input variables\n"
         "*rte_pressure*, *rte_temperature*, and *rte_vmr_list*, and returns the\n"
         "output variable *abs_scalar_gas*."
        ),
        OUTPUT( abs_scalar_gas_field_,
                abs_scalar_gas_,
                rte_pressure_, rte_temperature_, rte_vmr_list_),
        INPUT(  scalar_gas_absorption_agenda_,
                f_index_,
                f_grid_,
                atmosphere_dim_,
                p_grid_, lat_grid_, lon_grid_,
                t_field_, vmr_field_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_vecAddGas"),
        DESCRIPTION
        (
         "Add gas absorption to first element of absorption vector.\n"
         "\n"
         "The task of this method is to sum up the gas absorption of the\n"
         "different gas species and add the result to the first element of the\n"
         "absorption vector."
         ),
        OUTPUT(abs_vec_),
        INPUT(abs_vec_, abs_scalar_gas_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_vecAddGasZeeman"),
        DESCRIPTION
        (
         "Add zeeman absorption to the elements of absorption vector.\n"
         "\n"
	 ),
        OUTPUT(abs_vec_),
        INPUT(abs_vec_, abs_vec_zee_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_vecAddPart"),
        DESCRIPTION
        (
         "The particle absorption is added to *abs_vec* \n"
         "\n"
         "This function sums up the absorption vectors for all particle \n"
         "types weighted with particle number density.\n"
         "The resluling absorption vector is added to the workspace \n"
         "variable *abs_vec* \n"
         "Output and input of this method is *abs_vec* (stokes_dim).\n"
         "The inputs are the absorption vector for the single particle type \n"
         "*abs_vec_spt* (part_types, stokes_dim) and the local particle\n"
         " number densities for all particle types namely the \n"
         "*pnd_field* (part_types, p_grid, lat_grid, lon_grid, ) for given \n"
         "*p_grid*, *lat_grid*, and *lon_grid*. The particle types required \n"
         "are specified in the control file.  \n"
         ),
        OUTPUT(abs_vec_),
        INPUT(abs_vec_, abs_vec_spt_, pnd_field_, atmosphere_dim_,
              scat_p_index_,  scat_lat_index_, scat_lon_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_vecInit"),
        DESCRIPTION
        (
         "Initialize absorption vector.\n"
         "\n"
         "This method is necessary, because all other absorption methods just\n"
         "add to the existing absorption vector.\n"
         "\n"
         "So, here we have to make it the right size and fill it with 0.\n"
         "\n"
         "Note, that the vector is not really a vector, because it has a\n"
         "leading frequency dimension."
         ),
        OUTPUT(abs_vec_),
        INPUT(f_grid_, stokes_dim_, f_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("abs_vec_sptCalc"),
        DESCRIPTION
        (
         "This method calculates the absorption vector for a single particle\n"
         "type.\n"
         "\n"
         "All elements of absorption vector for a particle can be expressed\n"
         "in terms of extinction matrix and phase matrix. So it is necessary\n"
         "that the methods ext_mat_sptCalc and pha_mat_sptCalc are done \n"
         "before calling this method. \n"
         "\n"
         "The output of the method *abs_vec_sptCalc* is *abs_vec_spt* \n"
         "(Matrix, size: [Npt,stokes_dim]). The input to the method \n"
         "*abs_vec_sptCalc are abs_vec_spt, *pha_mat_spt*(Tensor 5,*\n"
         "size: [Npt,Nza,Naa,stokes_dim,stokes_dim]), *ext_mat_spt* \n"
         "(Tensor 3,size = [Npt,stokes_dim,stokes_dim]), *scat_za_grid*\n"
         "(Vector,size = [Nza]) and *scat_aa_grid* (Vector,size = [Naa]).\n"
         "This method calls the function amp2abs which does the actual\n"
         "physics, that of computing the elements of absorption vector \n"
         "from the elements of extinction matrix and phase matrix.\n"
         ),
        OUTPUT(abs_vec_spt_ ),
        INPUT(abs_vec_spt_, ext_mat_spt_, pha_mat_spt_, scat_za_grid_,
              scat_aa_grid_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AgendaSet"),
        DESCRIPTION
        ( 
         "Set up an agenda.\n"
         "\n"
         "An agenda is used to store a list of methods that are meant to be\n"
         "executed sequentially.\n"
         "\n"
         "This method takes the methods given in the body (in the curly braces)\n"
         "and puts them in the agenda given by the output argument (in the round\n"
         "braces).\n"
         "\n"
         "It also uses the agenda lookup data (defined in file agendas.cc) to\n"
         "check, whether the given methods use the right input WSVs and produce\n"
         "the right output WSVs.\n"
         " \n"
         "Generic output:\n"
         "   Agenda : The new agenda.\n"
         "\n"
         "Keywords:\n"
         "   No keywords, but other methods can appear in the method body."      
        ),
        OUTPUT(  ),
        INPUT(  ),
        GOUTPUT( Agenda_ ),
        GINPUT(),
        KEYWORDS(),
        TYPES(),
        AGENDAMETHOD( true )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("amp_matCalc"),
        DESCRIPTION
        (
         "This method converts the raw amplitude matrix data namely \n"
         "amp_mat_raw to workspace variable amp_mat which can be directly\n"
         "used for calculating scattering properties of particles.  \n"
         "\n"
         "The data type of amp_mat_raw is an ArrayOfArrayOfTensor6 which \n"
         "contains one gridded field for each particle type the user wants.\n"
         "One data file contains amp_mat_raw for one particle type.  One \n"
         "data file has amp_mat_raw for one particle for a set of \n"
         "frequencies, incoming and outgoing zenith angles and azimuth \n"
         "angles. The frequencies, angles, and the amplitude matrix are \n"
         "each a Tensor 6. The size of amp_mat_raw is amp_mat_raw[Npt][7],\n"
         "where Npt is the number of particle types. amp_mat_raw[Npt] [0]\n"
         "gives the frequency tensor [N_f, 1, 1, 1, 1, 1] where N_f gives\n"
         "the number of frequencies considered in that particular database \n"
         "file. Similarly, amp_mat_raw[ Npt ][ 1 ] gives the outgoing zenith\n"
         "angle tensor [1, Nza, 1, 1, 1, 1], amp_mat_raw[ Npt ][ 2 ] gives \n"
         "the outgoing azimuth angle tensor [1, 1, Naa, 1, 1, 1], \n"
         "amp_mat_raw[ Npt ][ 3 ] gives the incoming zentih angle tensor\n"
         "[1, 1, 1, Nza, 1, 1], amp_mat_raw[ Npt ][ 4 ] gives the incoming\n"
         "azimuth angle tensor [1, 1, 1, 1, Naa, 1], amp_mat_raw[ Npt ][ 5 ]\n"
         "is a dummy tensor6 and amp_mat_raw[ Npt ][ 6 ] gives amplitude\n"
         "matrix which is also a tensor6 of size \n"
         "[N_f, N_za, N_aa, N_za, N_aa, 8]. Here, Nza is the number of \n"
         "zenith angles, Naa is the number of azimuth angles and 8 denotes \n"
         "the amplitude matrix elements.  \n"
         "\n"
         "In this method, we have to interpolate the raw data calculated on \n"
         "specific angular and frequency grids onto a grid which is set by \n"
         "the user. Since we decide that frequency should be the outermost \n"
         "loop for our radiative transfer calculation, the frequency grid \n"
         "contains just one value specified by the index scat_f_index. The\n"
         "angles for which the calculations are to be done are specified by\n"
         "scat_za_grid and scat_aa_grid. The interpolation has to be done \n"
         "for the frequency grid, zenith angle grid and azimuth angle grid. \n"
         "Since this interpolation is from a gridded field to a new field, \n"
         "we have to perform a green interpolation. For more insight into \n"
         "the interpolation schemes refer to Chapter 8-Interpolation of AUG.\n"
         "\n"
         "The output of this method is amp_mat has to be a Tensor6 with the \n"
         "first dimension being that of the particle type, then the angles \n"
         "and finally the amplitude matrix element 8. The size of amp_mat is\n"
         "(Npt, Nza, Naa, Nza, Naa, 8).  Note that the dimension frequency \n"
         "is taken out.\n"
         ),
        OUTPUT(amp_mat_),
        INPUT(amp_mat_raw_, f_index_, f_grid_, scat_za_grid_, 
              scat_aa_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 

  md_data_raw.push_back
    ( MdRecord
      ( NAME("antenna_diagramAppendArray"),
        DESCRIPTION
        (
         "Appends a ArrayOfMatrix to *antenna_diagram*.\n"
         "\n"
         "This method can be used both to initialise and expand\n"
         "the viewing angles of *antenna_diagram*. At least one viewing\n"
         "angle must be given in *antenna_diagram*, and the array that\n"
         "is appended must have at least one element but not more than\n"
         "the number of polarisation given by *sensor_pol*."
        ),
        OUTPUT( antenna_diagram_ ),
        INPUT( sensor_pol_ ),
        GOUTPUT(),
        GINPUT( ArrayOfMatrix_ ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AntennaSet1D"),
        DESCRIPTION
        (
         "To set the antenna dimension to be 1D.\n"
         "\n"
         "Sets *antenna_dim* to 1 and sets *mblock_aa_grid* to be empty."
        ),
        OUTPUT( antenna_dim_, mblock_aa_grid_ ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AntennaSet2D"),
        DESCRIPTION
        (
         "To set the antenna dimension to be 2D.\n"
         "\n"
         "Sets *antenna_dim* to 2.\n"
         "\n"
         "It is only allowed to set *antenna_dim* to 2 when *atmosphere_dim*\n"
         "equals 3."
        ),
        OUTPUT( antenna_dim_ ),
        INPUT( atmosphere_dim_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixSet"),
        DESCRIPTION
        (
         "Sets a element in an ArrayOfMatrix to a Matrix.\n"
         "\n"
         "The keyword can be used to chose which element will be set, If a\n"
         "negative number is given, the matrix will be appended to the array.\n"
         "Note that zero-based indexing is used.\n"
         "\n"
         "Generic output:\n"
         "  ArrayOfMatrix : The array to be expanded.\n"
         "\n"
         "Generic input:\n"
         "         Matrix : The matrix to append.\n"
         "Keywords:\n"
         "        element : The index to be set."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfMatrix_ ),
        GINPUT( Matrix_ ),
        KEYWORDS( "element" ),
        TYPES( Index_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfStringSet"),
        DESCRIPTION
        (
         "Sets a String array according the given text.\n"
         "The format is text = [\"String1\",\"String2\",...]"
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfString_ ),
        GINPUT(),
        KEYWORDS( "text"         ),
        TYPES(    Array_String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmosphereSet1D"),
        DESCRIPTION
        (
         "To set the atmosheric dimension to be 1D.\n"
         "\n"
         "Sets *atmosphere_dim* to 1 and gives some variables dummy values.\n"
         "\n"
         "The latitude and longitude grids are set to be empty."
        ),
        OUTPUT( atmosphere_dim_, lat_grid_, lon_grid_ ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmosphereSet2D"),
        DESCRIPTION
        (
         "To set the atmosheric dimension to be 2D.\n"
         "\n"
         "Sets *atmosphere_dim* to 2 and gives some variables dummy values.\n"
         "\n"
         "The longitude grid is set to be empty. The variables *lat_1d*\n"
         "and *meridian_angle_1d* are given values that cause an error\n"
         "message if used."
        ),
        OUTPUT( atmosphere_dim_, lon_grid_, lat_1d_, meridian_angle_1d_),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmosphereSet3D"),
        DESCRIPTION
        (
         "To set the atmosheric dimension to be 3D.\n"
         "\n"
         "Sets *atmosphere_dim* to 3 and gives some variables dummy values.\n"
         "\n"
         "The variables *lat_1d* and *meridian_angle_1d* are given\n"
         "values that cause an error message if used."
        ),
        OUTPUT( atmosphere_dim_, lat_1d_, meridian_angle_1d_ ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmFieldsCalc"),
        DESCRIPTION
        (
         "Interpolated the atmospheric fields.\n"
         "\n"
         "An atmospheric scenario includes the following data for each \n"
         "position (pressure, latitude, longitude) in the atmosphere: \n"
         "           1. temperature field \n"
         "           2. the corresponding altitude field \n"
         "           3. vmr fields for the gaseous species \n"
         "This methods interpolates the fields from the raw data\n"
         "(*t_field_raw*, *z_field_raw*) which can be stored on \n"
         "arbitrary grids on the grids for the calculation\n"
         "(*p_grid*, *lat_grid*, *lon_grid*). "
        ),
        OUTPUT(t_field_, z_field_, vmr_field_),
        INPUT(p_grid_, lat_grid_, lon_grid_, t_field_raw_, z_field_raw_, 
              vmr_field_raw_, atmosphere_dim_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmFieldsCalcExpand1D"),
        DESCRIPTION
        (
         "Interpolate 1D raw atmospheric fields to create  2D or 3D \n"
         "homogenous atmospheric fields.\n"
         "\n"
         "The method works as *AtmFieldsCalc* but accepts only raw 1D\n"
         "atmsopheres. The raw atmsophere is interpolated to *p_grid* and \n"
         "the obtained values are applied for all latitudes, and also \n"
         "longitudes for 3D, to create a homogenous atmsophere. \n"
         "\n"
         "Note that the method only deals with the atmospheric fields, and\n"
         "to create a 2D or 3D version of a 1D case, a demand is also that\n"
         "the geoid radius is set to be constant for all latitudes/longitudes."
        ),
        OUTPUT( t_field_, z_field_, vmr_field_ ),
        INPUT( p_grid_, lat_grid_, lon_grid_, t_field_raw_, z_field_raw_, 
               vmr_field_raw_, atmosphere_dim_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

   md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmFieldsFromAscii2Xml"),
        DESCRIPTION
        (
         "Converts atmospheric data profiles (ascii-format) to xml format.\n"
         "\n"
         "This function can be used to convert data into the right format.\n"
         "The input files are matrices, holding the pressure grid and the \n"
         "data and output is a 3D gridded field (ArrayOfTensor3). We need \n"
         "this format, because generally we store the data not only on the \n"
         "pressure grid, but also on latitude and longitude grid. \n"
         "\n"
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("path", "basename"),
        TYPES(String_t, String_t)));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("AtmRawRead"),
        DESCRIPTION
        (
         "Reads atmospheric data from a scenario.\n"
         "\n"
         "An atmospheric scenario includes the following data for each \n"
         "position (pressure, latitude, longitude) in the atmosphere: \n"
         "           1. temperature field \n"
         "           2. the corresponding altitude field \n"
         "           3. vmr fields for the gaseous species \n"
         "The data is stored in different files. This methods reads all \n"
         "files and creates the variables *t_field_raw*, *z_field_raw* \n"
         "\n"
         "Different atmospheric scenarios are available in arts data:\n"
         "For example tropical and midlatitude-summer. 3D scenarios are \n"
         "not available yet.\n"
         "\n"
         "Files in the scenarios look like this: tropical.H2O.xml \n"
         "\n"
         "The basename must include the path, i.e., the files can be \n"
         "anywhere, but they must be all in the same directory.\n"
         "The profile is chosen by the species name. If you have more than \n"
         "one tag group for the same species, the same profile will be \n"
         "used.\n"
         "\n"
         "Keywords: \n"
         "basename :The name and path of a particular atmospheric scenario.\n"
         "For example:\n"
         "/smiles_local/arts-data/atmosphere/fascod/tropical \n"
         "\n"
        ),
        OUTPUT(t_field_raw_, z_field_raw_, vmr_field_raw_),
        INPUT(gas_species_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("basename"),
        TYPES(String_t)));


  md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_losSet" ),
        DESCRIPTION
        (
         "Sets *rte_los* to the given angles.\n"
         "\n"
         "The keyword argument *za* is put in as first element of *rte_los*\n"
         "and *aa* as the second element. However, when *atmosphere_dim* is\n"
         "set to 1D or 2D, the length of *rte_los* is set to 1 and only the\n"
         "given zenith angle is considered.\n"
         "\n"
         "Keywords: \n"
         "   za : Zenith angle of sensor line-of-sight.\n"
         "   aa : Azimuth angle of sensor line-of-sight."
        ),
        OUTPUT( rte_los_ ),
        INPUT( atmosphere_dim_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "za",      "aa"      ),
        TYPES(    Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_posAddGeoidWGS84" ),
        DESCRIPTION
        (
         "Adds a geoid radius according to WGS-84 to a geometric altitude.\n"
         "\n"
         "This function assumes that the first element of *rte_pos* is set\n"
         "to the geometric altitude for the position of the sensor. \n"
         "The variable *rte_pos* shall contain the radius instead of the\n"
         "altitude and that can be achieved by this function. The function\n"
         "adds a geoid radius to the given altitude. The geoid radius is\n"
         "taken from the WGS-84 reference ellipsiod.\n"
         "\n"
         "For 1D, the geoid radius is set to the radius of curvature of the\n"
         "WGS-84 ellipsiod for the position and observation direction \n"
         "described with *lat_1d* and *meridian_angle_1d*.\n"
         "For 2D and 3D, the geoid radius is set to the radius of the WGS-84\n"
         "ellipsiod for the latitude value in *rte_pos*."
        ),
        OUTPUT( rte_pos_ ),
        INPUT( rte_pos_, atmosphere_dim_, lat_1d_, meridian_angle_1d_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_posAddRgeoid" ),
        DESCRIPTION
        (
         "Adds a geoid radius by interpolating *r_geoid*.\n"
         "\n"
         "This function assumes that the first element of *rte_pos* is set\n"
         "to the geometric altitude for the position of the sensor. \n"
         "The variable *rte_pos* shall contain the radius instead of the\n"
         "altitude and that can be achieved by this function. The function\n"
         "adds a geoid radius to the given altitude. The geoid radius is\n"
         "obtained by interpolation of *r_geoid*. There is an error if the\n"
         "given position is outside the latitude and longitude grids."
        ),
        OUTPUT( rte_pos_ ),
        INPUT( rte_pos_, atmosphere_dim_, lat_grid_, lon_grid_, r_geoid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_posSet" ),
        DESCRIPTION
        (
         "Sets *rte_pos* to the given co-ordinates.\n"
         "\n"
         "The keyword argument *r_or_z* is put in as first element of\n"
         "*rte_pos*, *lat* as the second element and *lon* as third element.\n"
         "However, the length of *rte_pos* is set to *atmosphere_dim* and\n"
         "keyword arguments for dimensions not used are ignored.\n"
         "\n"
         "The first keyword argument can either be a radius, or an altitude\n"
         "above the geoid. In the latter case, a function such as\n"
         "*rte_posAddGeoidWGS84* should also be called to obtain a radius as\n"
         "first element of *rte_pos*.\n"
         "\n"
         "Keywords: \n"
         "   r_or_z : Radius or geometrical altitude of sensor position.\n"
         "   lat : Latitude of sensor position.\n"
         "   lon : Longitude of sensor position."
        ),
        OUTPUT( rte_pos_ ),
        INPUT( atmosphere_dim_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "r_or_z",  "lat",     "lon"     ),
        TYPES(    Numeric_t, Numeric_t, Numeric_t )));

md_data_raw.push_back
    ( MdRecord
      ( NAME( "Cloudbox_ppathCalc" ),
        DESCRIPTION
        (
         "Main function for calculation of propagation paths within the \n"
           "cloud box. This function will be used in *ScatteringMonteCarlo*\n"
         "\n"
         "\n"
         "The definition of a propgation path cannot be accomodated here.\n"
         "For more information read the chapter on propagation paths in the\n"
         "ARTS user guide and read the  on-line information for\n"
         "*ppath_step_agenda* (type \"arts -d ppath_step_agenda\")."
        ),
        OUTPUT( ppath_, ppath_step_ ),
        INPUT( ppath_step_agenda_, atmosphere_dim_, p_grid_, lat_grid_, 
               lon_grid_, z_field_, r_geoid_, z_ground_, 
               cloudbox_limits_, rte_pos_, rte_los_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));


  md_data_raw.push_back
    ( MdRecord
      ( NAME( "CloudboxGetIncoming" ),
        DESCRIPTION
        (
         "This method gives the intensity field at the boundary of the\n"
         "cloudbox.\n"
         "\n"
         "The method uses *PpathCalc* and *rte_agenda*.  The input to\n"
         "this method is the position of the cloudbox given by the \n"
         "variable *cloudbox_limits*. Then for each propagation direction\n"
         "it calls the function *PpathCalc* and executes the agenda \n"
         "*rte_agenda*.  This gives *i_rte* which holds the Stokes vector\n"
         "for all frequencies btained by the monochromatic pencil beam \n"
         "calculations performed by *rte_agenda*. Then this is copied to\n"
         "the interface variable. \n"
         "\n"
         ),
        OUTPUT(scat_i_p_, scat_i_lat_, scat_i_lon_, ppath_, ppath_step_,
               i_rte_, i_space_, ground_emission_, ground_los_,
               ground_refl_coeffs_,
               rte_los_, rte_pos_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_),
        INPUT( cloudbox_limits_, atmosphere_dim_, stokes_dim_, scat_za_grid_,
                scat_aa_grid_, f_grid_, ppath_step_agenda_,  rte_agenda_,
                i_space_agenda_, ground_refl_agenda_, p_grid_, lat_grid_,
                lon_grid_, z_field_, t_field_, r_geoid_, z_ground_),

        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "CloudboxGetIncoming1DAtm" ),
        DESCRIPTION
        (
         "This method gives the intensity field on the boundary of the\n"
         "cloudbox for a spherically symmetric clearsky atmosphere.\n"
         "\n"
         "The method uses *PpathCalc* and *rte_agenda*.  The input to\n"
         "this method is the position of the cloudbox given by the \n"
         "variable *cloudbox_limits*. Then for each propagation direction\n"
         "it calls the function *PpathCalc* and executes the agenda \n"
         "*rte_agenda*.  This gives *i_rte* which holds the Stokes vector\n"
         "for all frequencies btained by the monochromatic pencil beam \n"
         "calculations performed by *rte_agenda*. Then this is copied to\n"
         "the interface variable. \n"
         "\n"
         ),
        OUTPUT(scat_i_p_, scat_i_lat_, scat_i_lon_, ppath_, ppath_step_,
               i_rte_, i_space_, ground_emission_, ground_los_,
               ground_refl_coeffs_,
               rte_los_, rte_pos_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_),
        INPUT( cloudbox_limits_, atmosphere_dim_, stokes_dim_, scat_za_grid_,
                scat_aa_grid_, f_grid_, ppath_step_agenda_,  rte_agenda_,
                i_space_agenda_, ground_refl_agenda_, p_grid_, lat_grid_,
                lon_grid_, z_field_, t_field_, r_geoid_, z_ground_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "CloudboxGetOutgoing" ),
        DESCRIPTION
        (
         "Scattered radiance on the cloudbox boundary.\n"
         "\n"
         "This method returns the radiances for a given direction and \n"
         "position on the boundary of the cloudbox. The variable *y_scat* \n"
         "is a matrix with the dimensions [f_grid, stokes_dim].\n"
         "\n"
         ),
        OUTPUT(),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, rte_gp_p_, rte_gp_lat_,
               rte_gp_lon_, rte_los_,  cloudbox_on_, cloudbox_limits_,
               atmosphere_dim_, stokes_dim_,
               scat_za_grid_, scat_aa_grid_, f_grid_),
        GOUTPUT(Matrix_),
        GINPUT(),
        KEYWORDS(),
        TYPES()));


  md_data_raw.push_back
    ( MdRecord
      ( NAME( "CloudboxSetIncomingForTauCalc1D" ),
        DESCRIPTION
        (
         "This method sets the intensity field at the pressure boundaries of\n"
         "the cloudbox.\n"
         "\n"
         "The method sets the upper intensity field to zero and the lower\n"
         "is calculated as 1/-cos(za), where angles za are taken from the\n"
         "scat_za_grid vector. This gives intensity 1 from the downlooking\n"
         "direction and increasing values for decreasing zenith angle.\n"
         "\n"
         "The incoming clearsky radiation is set to zero at the vertical\n"
         "limits of the cloudbox.\n"
         "\n"
         "This is implemented only for the case of 1D atmosphere at present\n"
         "\n"
         "Keywords: \n"
         "   za_low : Lower zenith angle limit for incoming intensity.\n"
         ),
        OUTPUT(scat_i_p_, scat_i_lat_, scat_i_lon_ ),
        INPUT( cloudbox_limits_, atmosphere_dim_, stokes_dim_, scat_za_grid_,
               f_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "za_low" ),
        TYPES( Numeric_t )));


  md_data_raw.push_back
    ( MdRecord
      ( NAME("cloudboxOff"),
        DESCRIPTION
        (
         "Deactivates the cloud box. \n"
         "\n"
         "The function sets *cloudbox_on* to 0, and *cloudbox_limits* to be\n"
         "an empty vector. The variables *scat_i_p*, *scat_i_lat*,  \n"
         "*scat_i_lon*, *scat_za_grid* and *scat_aa_grid* are also set to be\n"
         "empty."
        ),
        OUTPUT( cloudbox_on_, cloudbox_limits_, scat_i_p_, scat_i_lat_,
                scat_i_lon_, scat_za_grid_, scat_aa_grid_ ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "cloudboxSetManually" ),
        DESCRIPTION
        (
         "Sets the cloud box to encompass the given positions.\n"
         "\n"
         "The function sets *cloudbox_on* to 1 and sets *cloudbox_limits*\n"
         "following the given pressure, latitude and longitude positions.\n"
         "The index limits in *cloudbox_limits* are selected to give the\n" 
         "smallest possible cloud box that encompass the given points. \n"
         "\n"
         "The points must be given in the same order as used in\n"
         "*cloudbox_limits*. That means that the first keyword argument \n"
         "shall be a higher pressure than argument two, while the latitude\n"
         "and longitude points are given in increasing order. Positions\n"
         "given for dimensions not used by the selected atmospheric\n"
         "dimensionality are ignored.\n"
         "\n"
         "The given pressure points can be outside the range of *p_grid*.\n"
         "The pressure limit is then set to the end point of *p_grid*.\n"
         "The given latitude and longitude points must be inside the range\n"
         "of the corresponding grid. In addition, the latitude and longitude\n"
         "points cannot be inside the outermost grid ranges as the latitude\n"
         "and longitude limits in *cloudbox_limits* are not allowed to be\n"
         "grid end points.\n"
         "\n"
         "Keywords: \n"
         "   p1   : Upper pressure point.\n"
         "   p2   : Lower pressure point.\n"
         "   lat1 : Lower latitude point.\n"
         "   lat2 : Upper latitude point.\n"
         "   lon1 : Lower longitude point.\n"
         "   lon2 : Upper longitude point."
        ),
        OUTPUT( cloudbox_on_, cloudbox_limits_ ),
        INPUT( atmosphere_dim_, p_grid_, lat_grid_, lon_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( "p1", "p2", "lat1", "lat2", "lon1", "lon2" ),
        TYPES( Numeric_t, Numeric_t, Numeric_t, Numeric_t, Numeric_t, 
               Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ConvertIFToRF"),
        DESCRIPTION
        (
         "Convert *sensor_response_f* from IF to RF, the function also\n"
         "unfolds the measurement spectra *y*.\n"
         "\n"
         "This function should be used when the sensor configuration contains\n"
         "a mixer and the spectra should be given in brightness temperature.\n"
         "The reason is that the mixer converts the RF to IF, and to be able\n"
         "to perform the Rayleigh-Jeans conversion from radiance to\n"
         "brightness temperature, the radiance needs to be given in RF.\n"
         "\n"
         "Note that the number of elements in both *sensor_response_f* and\n"
         "*y* will potentially increase, since the IF is mapped to both the\n"
         "lower and upper sidebands.\n"
         "\n"
         "Keyword: \n"
         "   output : Which sideband(s) to output, \"lower\", \"upper\" or\n"
         "            \"double\""
         ),
         OUTPUT( sensor_response_f_, y_ ),
         INPUT( sensor_pol_, sensor_response_za_, sensor_response_aa_, lo_,
                atmosphere_dim_ ),
         GOUTPUT( ),
         GINPUT( ),
         KEYWORDS( "output" ),
         TYPES( String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Copy"),
        DESCRIPTION
        (
         "Copy a workspace variable.\n"
         "\n"
         "This is a supergeneric method. It can copy any workspace variable\n"
         "to another workspace variable of the same group. (E.g., a Matrix to\n"
         "another Matrix.)\n"
         "\n"
         "As allways, output comes first in the argument list!\n"
         "\n"
         "Usage example:\n"
         "\n"
         "Copy(f_grid,p_grid){}\n"
         "\n"
         "Will copy the content of *p_grid* to *f_grid*. The size of *f_grid*\n"
         "is adjusted automatically (the normal behaviour for workspace\n"
         "methods).\n"
         "\n"
         "Supergeneric output:\n"
         "   Any_ : The output variable.\n"
         "\n"
         "Supergeneric input:\n"
         "   Any_ : The input variable."
         ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Any_ ),
        GINPUT(  Any_ ),
        KEYWORDS(),
        TYPES(),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("convergence_flagAbs"),
        DESCRIPTION
        (
         "Convergence test (maximum absolute difference).\n"
         "\n"
         "The function calculates the absolute differences for two successive\n"
         "iteration fields. It picks out the maximum values for each Stokes \n"
         "component separately. The convergence test is fullfilled under the\n"
         "following conditions: \n"
         "|I(m+1) - I(m)| < epsilon_1     Intensity.\n"
         "|Q(m+1) - Q(m)| < epsilon_2     The other Stokes components.\n" 
         "|U(m+1) - U(m)| < epsilon_3    \n"
         "|V(m+1) - V(m)| < epsilon_4    \n" 
         "\n"
         "The limits for convergence have to be set in the controlfile by \n"
         "setting the vector *epsilon* to appropriate values.\n"
         "\n"
         "The conditions have to be valid for all positions in the cloudbox \n"
         "and for all directions.\n"  
         "Then *convergence_flag* is set to 1.\n"
         "\n"
         "Unit of *epsilon* is that of radiance.\n"
        ),
        OUTPUT(convergence_flag_),
        INPUT(i_field_, i_field_old_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS("epsilon"),
        TYPES(Vector_t)));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("convergence_flagAbs_BT"),
        DESCRIPTION
        (
         "Convergence test (maximum absolute difference in Rayleigh-Jeans\n"
         "Brightness temperature units)\n"
         "\n"
         "This method work exactly similar to convergence_flagAbs except that\n"
         "now we se the convergence criteria keyword *epsilon* to Brightness\n"
         "temperature units. /n"
         "\n"
         "Note that we use Rayleigh Jeans Brightness temperature for epsilon.\n"
         "This is because epsilon is a difference of intensity and Planck BT\n"
         "is not linear for small radiance values.  For higher stokes components\n"
         "also Planck BT cannot be used because of the same reason.\n"
         "\n"
         "The function calculates the absolute differences for two successive\n"
         "iteration fields. It picks out the maximum values for each Stokes \n"
         "component separately. The convergence test is fullfilled under the\n"
         "following conditions: \n"
         "|I(m+1) - I(m)| < epsilon_1     Intensity.\n"
         "|Q(m+1) - Q(m)| < epsilon_2     The other Stokes components.\n" 
         "|U(m+1) - U(m)| < epsilon_3    \n"
         "|V(m+1) - V(m)| < epsilon_4    \n" 
         "\n"
         "The limits for convergence have to be set in the controlfile by \n"
         "setting the vector *epsilon* to appropriate values.\n"
         "\n"
         "The conditions have to be valid for all positions in the cloudbox \n"
         "and for all directions.\n"  
         "Then *convergence_flag* is set to 1.\n"
         "\n"
         "Unit of *epsilon* is that of brightness temperature(RJ).\n"
        ),
        OUTPUT(convergence_flag_),
        INPUT(i_field_, i_field_old_, f_grid_, f_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS("epsilon"),
        TYPES(Vector_t)));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Exit"),
        DESCRIPTION
        (
         "Stops the execution and exits ARTS.\n"
         "\n"
         "This method is handy if you want to debug one of your control\n"
         "files. You can insert it anywhere in the control file. When\n"
         "it is reached, it will terminate the program."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_matAddGas"),
        DESCRIPTION
        (
         "Add gas absorption to all diagonal elements of extinction matrix.\n"
         " \n"
         "The task of this method is to sum up the gas absorption of the\n"
         "different gas species and add the result to the extinction matrix."
         ),
        OUTPUT(ext_mat_),
        INPUT(ext_mat_, abs_scalar_gas_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_matAddGasZeeman"),
        DESCRIPTION
        (
         "Add Zeeman extinction  to the elements of extinction matrix.\n"
         " \n"
	 ),
        OUTPUT(ext_mat_),
        INPUT(ext_mat_, ext_mat_zee_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_matAddPart"),
        DESCRIPTION
        (
         "The particle extinction is added to *ext_mat* \n"
         "\n"
         "This function sums up the extinction matrices for all particle \n"
         "types weighted with particle number density.\n"
         "The resluling extinction matrix is added to the workspace \n"
         "variable *ext_mat* \n"
         "The output of this method is *ext_mat* (stokes_dim, stokes_dim).\n"
         "The inputs are the extinction matrix for the single particle type \n"
         "*ext_mat_spt* (part_types, stokes_dim, stokes_dim) and the local \n"
         "particle number densities for all particle types namely the \n"
         "*pnd_field* (part_types, p_grid, lat_grid, lon_grid ) for given \n"
         "*p_grid*, *lat_grid*, and *lon_grid*. The particle types required \n"
         "are specified in the control file.  \n"
         ),
        OUTPUT( ext_mat_  ),
        INPUT( ext_mat_, ext_mat_spt_, pnd_field_, atmosphere_dim_, 
               scat_p_index_, scat_lat_index_, scat_lon_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_matInit"),
        DESCRIPTION
        (
         "Initialize extinction matrix.\n"
         "\n"
         "This method is necessary, because all other extinction methods just\n"
         "add to the existing extinction matrix. \n"
         "\n"
         "So, here we have to make it the right size and fill it with 0.\n"
         "\n"
         "Note, that the matrix is not really a matrix, because it has a\n"
         "leading frequency dimension."
         ),
        OUTPUT(ext_mat_ ),
        INPUT(f_grid_, stokes_dim_, f_index_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("ext_mat_partScat"),
//      DESCRIPTION
//      (
//       "This function sums up the convergence extinction matrices \n"
//          "for all particle types weighted with particle number density.\n"
//       "\n"
//          "This method is only used for the convergence test, where particle\n"
//          "absorption is set to be 0. \n"
//       "The output of this method is *ext_mat_part* (stokes_dim, stokes_dim).\n"
//       "The inputs are the convergence extinction matrix for the single\n"
//          "particle type *ext_mat_conv_spt* \n"
//          "(part_types, stokes_dim, stokes_dim) and the local \n"
//       "particle number densities for all particle types namely the \n"
//       "*pnd_field* (part_types, p_grid, lat_grid, lon_grid ) for given \n"
//       "*p_grid*, *lat_grid*, and *lon_grid*. The particle types required \n"
//       "are specified in the control file.  \n"
//       ),
//      OUTPUT( ext_mat_  ),
//         INPUT( ext_mat_spt_, pnd_field_, atmosphere_dim_, scat_p_index_, 
//             scat_lat_index_, scat_lon_index_),
//      GOUTPUT( ),
//      GINPUT( ),
//      KEYWORDS( ),
//      TYPES( )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_mat_sptScat"),
        DESCRIPTION
        ( 
         "This method calculates the convergence extinction matrix.\n"
         "\n"
         "This function calculates extinction due to scattering (without \n"
         "absorption). It is used only for testing, if the iterative method \n"
         "converges towards the right solution.\n"
         "\n"
         ),
        OUTPUT( ext_mat_spt_  ),
        INPUT( pha_mat_spt_, scat_za_grid_, scat_aa_grid_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));
 
 md_data_raw.push_back
    ( MdRecord
      ( NAME("ext_mat_sptCalc"),
        DESCRIPTION
        ( 
         "This method calculates the extinction matrix for a single particle\n"
         "type.\n"
         "\n"
         "Extinction matrix describes the total attenuation of the incident \n"
         "radiation resulting from the combined effect of scattering and\n"
         "absorption by the particle.  It is a 4X4 matrix and all the\n"
         "elements of extinction matrix for a particle can be expressed in\n"
         "terms of the elements of the forward scattering amplitude matrix.\n"
         "\n"
         "The output of this method is \n"
         "*ext_mat_spt*(Tensor 3, size = [Npt,stokes_dim,stokes_dim])and the\n"
         "inputs are *ext_mat_spt*,*amp_mat*(Tensor 6, Size=[Npt,Nza,Naa,Nza,Naa,8]), \n"
         "*scat_za_index*,*scat_aa_index*,*f_index* and *scat_f_grid*. \n"
         "\n"
         "The variables *scat_za_index* and *scat_aa_index picks the right \n"
         "element of the Tensor *amp_mat*. *f_grid* and *f_index* picks \n"
         "the right frequeny for calculation. Frequeny is needed because the\n"
         "computation of extinction matrix from amplitude matrix involves \n"
         "multiplication by wavelength.  Then this method calls the \n"
         "function amp2ext which does the actual physics, that of computing\n"
         "the elements of extinction matrix from the elements of amplitude \n"
         "matrix.\n"
         ),
        OUTPUT( ext_mat_spt_  ),
        INPUT( ext_mat_spt_,amp_mat_, scat_za_index_, scat_aa_index_,
               f_index_, f_grid_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));


  md_data_raw.push_back     
    ( MdRecord
      ( NAME("FlagOff"),
        DESCRIPTION
        (
         "Sets an index variable that acts as an on/off flag to 0."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Index_ ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("FlagOn"),
        DESCRIPTION
        (
         "Sets an index variable that acts as an on/off flag to 1."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Index_ ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("gas_abs_lookupAdapt"),
        DESCRIPTION
        (
         "Adapts a gas absorption lookup table to the current calculation.\n"
         "\n"
         "The lookup table can contain more species and more frequencies than\n"
         "are needed for the current calculation. This method cuts down the\n"
         "table in memory, so that it contains just what is needed. Also, the\n"
         "species in the table are brought in the same order as the species in\n"
         "the current calculation.\n"
         "\n"
         "Of course, the method also performs quite a lot of checks on the\n"
         "table. If something is not ok, a runtime error is thrown."
        ),
        OUTPUT( gas_abs_lookup_, gas_abs_lookup_is_adapted_ ),
        INPUT(  gas_abs_lookup_, gas_species_, f_grid_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("gas_abs_lookupInit"),
        DESCRIPTION
        (
         "Creates an empty gas absorption lookup table.\n"
         "\n"
         "This is mainly there to help developers. For example, you can write\n"
         "the empty table to an XML file, to see the file format."
        ),
        OUTPUT( gas_abs_lookup_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("gas_speciesSet"),
        DESCRIPTION
        (
         "Set up the list of gas species tag groups.\n"
         "\n"
         "The workspace variable *gas_species* contains several tag groups. Each\n"
         "tag group contain one or more tags. This method converts descriptions\n"
         "of tag groups given in the keyword to the internal representation of\n"
         "*gas_species*. A tag group selects spectral features which belong to\n"
         "the same species.\n"
         "\n"
         "A tag is defined in terms of the name of the species, isotope, and a\n"
         "range of frequencies. Species are named after the standard chemical\n"
         "names, e.g., \"O3\".  Isotopes are given by the last digit of the atomic\n"
         "weight, i.g., \"O3-668\" for the asymmetric ozone molecule including an\n"
         "oxygen 18 atom. Groups of transitions are specified by giving a lower\n"
         "and upper limit of a frequency range, e.g., \"O3-666-500e9-501e9\".\n"
         "\n"
         "The symbol \"*\" acts as a wild card. Furthermore, frequency range or\n"
         "frequency range and isotope may be omitted.\n"
         "\n"
         "Finally, instead of the isotope the special letter \"nl\" may be given,\n"
         "e.g., \"H2O-nl\". This means that no lines of this species should be\n"
         "included in the general line-by-line calculation. This feature is\n"
         "useful if you want to define a tag group just for a continuum, or for\n"
         "a complete absorption model.\n"
         "\n"
         "Keywords:\n"
         "   species : Specify one String for each tag group that you want to\n"
         "             create. Inside the String, separate the tags by commas\n"
         "             (plus optional blanks).\n"
         "\n"
         "Example:\n"
         "\n"
         "   species = [ \"O3-666-500e9-501e9, O3-686\",\n"
         "               \"O3\",\n"
         "               \"H2O-nl\" ]\n"
         "\n"
         "   The first tag group selects all O3-666 lines between 500 and\n"
         "   501 GHz plus all O3-686 lines.  \n"
         "\n"
         "   The second tag group selects all remaining O3 transitions.\n"
         "\n"
         "   The third tag group selects H2O, but will not put any lines in the\n"
         "   line list for this species. Presumably, we are using a complete\n"
         "   absorption model like MPM89 for H2O in this case."
         ),
        OUTPUT( gas_species_ ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "species" ),
        TYPES(    Array_String_t   )));
 
  md_data_raw.push_back
    ( MdRecord
      ( NAME("GaussianResponse"),
        DESCRIPTION
        (
         "Creates a matrix with a relative grid and a gaussian response.\n"
         "\n"
         "The grid is a relative grid so therefore it is centred around\n"
         "zero, the TotWidth keyword describes the difference between the\n"
         "maximum grid point and the minimum. The grid range is then divided\n"
         "into grid points equally spaced with max distance equal or less\n"
         "than MaxSpacing. The results are the stored in columns in the matrix.\n"
         "Grid points in the first and values in the second column.\n"
         "\n"
         "Generic output: \n"
         "   Matrix     : The matrix containing the grid and response values.\n"
         "\n"
         "Keywords:\n"
         "   fwhm       : The Full Width at Half Mean value for the response.\n"
         "   tot_width   : The total width of the relative grid.\n"
         "   max_spacing : The maximum step between grid points."
         ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT(),
        KEYWORDS( "fwhm", "tot_width", "max_spacing" ),
        TYPES( Numeric_t, Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Ignore"),
        DESCRIPTION
        (
         "Ignore a workspace variable.\n"
         "\n"
         "This method is handy for use in agendas in order to suppress warnings\n"
         "about unused input workspace variables. What it does is: Nothing!\n"
         "In other words, it just ignores the variable it is called on.\n"
         "\n"
         "This is a supergeneric method. It can ignore any workspace variable\n"
         "you want.\n"
         "\n"
         "Usage example:\n"
         "\n"
         "AgendaSet(els_agenda){\n"
         "  Ignore(ls_sigma){}\n"
         "  elsLorentz{}\n"
         "}\n"
         "\n"
         "Without Ignore you would get an error message, because els_agenda is\n"
         "supposed to use the Doppler width *ls_sigma*, but the Lorentz lineshape\n"
         "*elsLorentz* does not need it.\n"
         "\n"
         "Supergeneric input:\n"
         "   Any_ : The input variable."
         ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT(  Any_ ),
        KEYWORDS(),
        TYPES(),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldIterate" ),
        DESCRIPTION
        (
         "Iterative solution of the RTE.\n"
         "\n"
         "A solution for the RTE with scattering is found using an iterative\n"
         "scheme:\n"
         "\n"
         "1. Calculate scattering integral.\n"
         "   To calculate the scattering integral the total scattering \n"
         "   matrix is required. This is obtained using the functions \n"
         "   *pha_mat_sptCalc* and  *pha_matCalc*. The method *.....* \n"
         "   performs the integration.\n"
         "2. Calculate RT with fixed scattered field.\n"
         "   The radiative transfer equation with fixed scattering integral\n"
         "   can be solved analytically if the coefficients are assumed to \n"
         "   be constant.\n"
         "   According to *atmosphere_dim* either *i_fieldUpdate1D* or \n"
         "   *i_fieldUpdate2D* are called to perform the calculation. Inside\n"
         "   these methods the agenda *scat_rte_agenda* is executed. \n"
         "3. Convergence test.\n"
         "   Here the *convergence_test_agenda* is executed.\n"
         "\n"
         "Note: The atmospheric dimensionality *atmosphere_dim* can be \n"
         "      either 1 or 3. To these dimensions the method adapts \n"
         "      automatically. \n"
         "      If *atmosphere_dim* equals 2, it returns an error message,\n"
         "      as 2D scattering calculations can not be performed.\n"
         ),
        OUTPUT(i_field_, i_field_old_, convergence_flag_),
        INPUT( scat_field_agenda_, scat_rte_agenda_, 
               convergence_test_agenda_, f_grid_,
               f_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldSetClearsky" ),
        DESCRIPTION
        (
         "This method uses a linear 3D interpolation scheme to obtain the\n"
         "radiation field on all grid points inside the cloud box from the\n"
         "clear sky field on the cloud bod boundary. This radiation field\n"
         "is taken as the first guess radiation field for the iterative \n"
         "solution method of the RTE.\n"
         "\n"
         "The inputs to this method are *scat_i_p*[ N_f, 2, N_lat, N_lon,\n"
         "N_za, N_aa, stokes_dim], *scat_i_lat*[ N_f, N_p, 2, N_lon, \n"
         "N_za, N_aa, stokes_dim] and *scat_i_lon*[ N_f, N_p, N_lat, 2,\n"
         "N_za, N_aa, stokes_dim].  The method has to pick the \n"
         "monochromatic radiation field out of these variables.  The \n"
         "output of the method is the initial field stored in the \n"
         "workspace variable *i_field*[ N_p, N_lat, N_lon, N_za, N_aa,\n"
         "stokes_dim].\n"
         ),
        OUTPUT(i_field_),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, f_grid_, 
               f_index_, p_grid_, lat_grid_, lon_grid_, 
               cloudbox_limits_, atmosphere_dim_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldSetConst" ),
        DESCRIPTION
        (
         "This method can be used to set the initial field inside the \n"
         "cloudbox to a constant value. \n"
         "The value can be given by the user .\n"
         "\n"
         "The inputs to this method are *scat_i_p*[ N_f, 2, N_lat, N_lon,\n"
         "N_za, N_aa, stokes_dim], *scat_i_lat*[ N_f, N_p, 2, N_lon, \n"
         "N_za, N_aa, stokes_dim] and *scat_i_lon*[ N_f, N_p, N_lat, 2,\n"
         "N_za, N_aa, stokes_dim].  The method has to pick the \n"
         "monochromatic radiation field out of these variables.  The \n"
         "output of the method is the initial field stored in the \n"
         "workspace variable *i_field*[ N_p, N_lat, N_lon, N_za, N_aa,\n"
         "stokes_dim].\n"
         ),
        OUTPUT(i_field_),
        INPUT( scat_i_p_, scat_i_lat_, scat_i_lon_, f_grid_, 
               f_index_, p_grid_, lat_grid_, lon_grid_, 
               cloudbox_limits_, atmosphere_dim_, stokes_dim_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("value"),
        TYPES(Vector_t)));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldUpdate1D" ),
        DESCRIPTION
        (
         "Updates the i_field during the iteration. It performs the RT \n"
         "calculation using a fixed value for the scattering integral stored\n"
         "in *scat_field*.\n"
         "\n" 
        ),
        OUTPUT(i_field_, rte_pressure_, rte_temperature_,
               rte_vmr_list_, scat_za_index_, ext_mat_, abs_vec_,
               scat_p_index_, ppath_step_),
        INPUT(i_field_old_, scat_field_, cloudbox_limits_, 
              scalar_gas_absorption_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, 
              opt_prop_part_agenda_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, z_field_, r_geoid_, t_field_,
              f_grid_, f_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldUpdateSeq1D" ),
        DESCRIPTION
        (
         "Updates the i_field during the iteration sequentially.\n"
         "It performs the RT \n"
         "calculation using a fixed value for the scattering integral stored \n"
         "in *scat_field*.\n"
         "\n" 
         ),
        OUTPUT(i_field_, rte_pressure_, rte_temperature_,
               rte_vmr_list_, scat_za_index_, ext_mat_, abs_vec_,
               scat_p_index_, ppath_step_, ground_los_, ground_emission_,
               ground_refl_coeffs_, rte_los_, rte_pos_, rte_gp_p_),
        INPUT(scat_field_, cloudbox_limits_, 
              scalar_gas_absorption_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, 
              opt_prop_part_agenda_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, z_field_, r_geoid_, t_field_,
              f_grid_, f_index_, ground_refl_agenda_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));


  md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldUpdate3D" ),
        DESCRIPTION
        (
         "Updates the i_field during the iteration. It performs the RT \n"
         "calculation using a fixed value for the scattering integral stored \n"
         "in *scat_field*.\n"
         "\n " 
        ),
        OUTPUT(i_field_, rte_pressure_, rte_temperature_,
               rte_vmr_list_, scat_za_index_, scat_aa_index_, ext_mat_, abs_vec_,
               scat_p_index_, scat_lat_index_, scat_lon_index_,  ppath_step_),
        INPUT(i_field_old_, scat_field_, cloudbox_limits_, 
              scalar_gas_absorption_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, scat_aa_grid_,
              opt_prop_part_agenda_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, lat_grid_, lon_grid_, z_field_,
              r_geoid_, t_field_,
              f_grid_, f_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

   md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldUpdateSeq3D" ),
        DESCRIPTION
        (
         "Updates the i_field during the iteration sequentially.\n"
         "It performs the RT \n"
         "calculation using a fixed value for the scattering integral stored \n"
         "in *scat_field*.\n"
         "\n " 
        ),
        OUTPUT(i_field_, rte_pressure_, rte_temperature_,
               rte_vmr_list_, scat_za_index_, scat_aa_index_, ext_mat_, abs_vec_,
               scat_p_index_, scat_lat_index_, scat_lon_index_,  ppath_step_),
        INPUT(scat_field_, cloudbox_limits_, 
              scalar_gas_absorption_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, scat_aa_grid_,
              opt_prop_part_agenda_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, lat_grid_, lon_grid_, z_field_,
              r_geoid_, t_field_,
              f_grid_, f_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "i_fieldUpdateSeq1D_PlaneParallel" ),
        DESCRIPTION
        (
         "Updates the i_field during the iteration sequentially.\n"
         "It performs the RT \n"
         "calculation using a fixed value for the scattering integral stored \n"
         "in *scat_field*.\n"
         "\n" 
         ),
        OUTPUT(i_field_, rte_pressure_, rte_temperature_,
               rte_vmr_list_, scat_za_index_, ext_mat_, abs_vec_,
               scat_p_index_, ppath_step_, ground_los_, ground_emission_,
               ground_refl_coeffs_, rte_los_, rte_pos_, rte_gp_p_),
        INPUT(scat_field_, cloudbox_limits_, 
              scalar_gas_absorption_agenda_,
              vmr_field_, spt_calc_agenda_, scat_za_grid_, 
              opt_prop_part_agenda_, opt_prop_gas_agenda_,
              ppath_step_agenda_, p_grid_, z_field_, r_geoid_, t_field_,
              f_grid_, f_index_, ground_refl_agenda_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("IndexSet"),
        DESCRIPTION
        (
         "Sets an index workspace variable to the given value. \n"
         "\n"
         "Generic output: \n"
         "   Index : The index variable to be set. \n"
         "\n"
         "Keywords:\n"
         "   value : A positive integer." 
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Index_ ),
        GINPUT( ),
        KEYWORDS( "value" ),
        TYPES(     Index_t   )));
  
  md_data_raw.push_back     
    ( MdRecord
      ( NAME("iteration_counterIncrease"),
        DESCRIPTION
        (
         "Increase iteration counter. \n"
         "\n"
         "This function can be used for writing the separate iteration \n"
         "fields into differtent files using *Tensor6WriteIteration*.\n"
         ),
        OUTPUT(iteration_counter_),
        INPUT(iteration_counter_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));
  
   md_data_raw.push_back
    ( MdRecord
      ( NAME("elsDoppler"),
        DESCRIPTION
        (
         "The Doppler lineshape.\n"
         "\n"
         "This computes the Doppler lineshape as:\n"
         "\n"
         "els[i] = 1/(sqrt(PI) * ls_sigma) /\n"
         "         exp(-(els_f_grid[i]/ls_sigma)^2)\n"
         "\n"
         "Note that the frequency grid els_f_grid must hold\n"
         "offset frequencies from line center. Hence, the\n"
         "line center frequency is not needed as input.\n"
         "\n"
         "Output:\n"
         "   els        : The lineshape function [1/Hz]\n"
         "\n"
         "Input:\n"
         "   ls_sigma   : Line width [Hz].\n"
         "   els_f_grid : Frequency grid [Hz]."
        ),
        OUTPUT( els_ ),
        INPUT(  ls_sigma_, els_f_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 

  md_data_raw.push_back
    ( MdRecord
      ( NAME("elsLorentz"),
        DESCRIPTION
        (
         "The Lorentz lineshape.\n"
         "\n"
         "This computes the simple Lorentz lineshape as:\n"
         "\n"
         "els[i] = 1/PI * ls_gamma /\n"
         "         ( (els_f_grid[i])^2 + ls_gamma^2 )\n"
         "\n"
         "Note that the frequency grid els_f_grid must hold\n"
         "offset frequencies from line center. Hence, the\n"
         "line center frequency is not needed as input.\n"
         "\n"
         "Output:\n"
         "   els        : The lineshape function [1/Hz]\n"
         "\n"
         "Input:\n"
         "   ls_gamma   : Line width [Hz].\n"
         "   els_f_grid : Frequency grid [Hz]."
        ),
        OUTPUT( els_ ),
        INPUT(  ls_gamma_, els_f_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("elsVoigt_Drayson"),
        DESCRIPTION
        (
         "The Voigt Drayson linshape.\n"
         "\n"
         "This computes the Voigt profile\n"
         "through the Drayson approximation.\n"
         "\n"
         "Note that the frequency grid els_f_grid must hold\n"
         "offset frequencies from line center. Hence, the\n"
         "line center frequency is not needed as input.\n"
         "\n"
         "Output:\n"
         "   els        : The lineshape function [1/Hz]\n"
         "\n"
         "Input:\n"
         "   ls_sigma   : Lorentz width [Hz].\n"
         "   ls_gamma   : Doppler width [Hz].\n"
         "   els_f_grid : Frequency grid [Hz]."
        ),
        OUTPUT( els_ ),
        INPUT(  ls_sigma_, ls_gamma_, els_f_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("elsVoigt_Kuntz6"),
        DESCRIPTION
        (
         "The Voigt-Kuntz linshape.\n"
         "\n"
         "This computes the Voigt profile\n"
         "through the Kuntz algorithm with\n"
         "a relative accuracy better than 2*10-6.\n"
         "\n"
         "Note that the frequency grid els_f_grid must hold\n"
         "offset frequencies from line center. Hence, the\n"
         "line center frequency is not needed as input.\n"
         "\n"
         "Output:\n"
         "   els        : The lineshape function [1/Hz]\n"
         "\n"
         "Input:\n"
         "   ls_sigma   : Lorentz width [Hz].\n"
         "   ls_gamma   : Doppler width [Hz].\n"
         "   els_f_grid : Frequency grid [Hz]."
        ),
        OUTPUT( els_ ),
        INPUT(  ls_sigma_, ls_gamma_, els_f_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("grid_stepsizeCheck"),
        DESCRIPTION
        (
         "Check grid_stepsize."
         "\n"
         "This method calculates the value of the constant gridstepsize of\n"
         "*scat_za_grid* and *scat_aa_grid* and stores it\n"
         "in the Vector *grid_stepsize*\n"
         "\n"
         "If one of the grid stepsizes isn't constant, -1 will be stored\n"
         "\n"
         "grid_stepsize[0] <-> scat_za_grid\n"
         "grid_stepsize[1] <-> scat_aa_grid\n"
         ),
        OUTPUT( grid_stepsize_ ),
        INPUT( scat_za_grid_, scat_aa_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("grid_sizeSet"),
        DESCRIPTION
        (
         "Set *grid_size* for scattering integral."
         "\n"
         "In this method the grid sizes (number of points of zenith grid \n"
         "and azimuth  angle grid) for the scattering integral calculation \n"
         "can be specified. For the calculation of the scattering integral\n"
         "equidistant grids are appropriate as the peak of the phase \n"
         "function can be anywhere, depending on incoming and scattered \n"
         "directions." 
         "\n"
         "Keywords: \n"
         "  N_za_grid : Number of points in zenith angle grid. \n"
         "  N_aa_grid : Number of points in azimuth angle grid. \n"
         ),
        OUTPUT( za_grid_size_, scat_aa_grid_),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("N_za_grid", "N_aa_grid"),
        TYPES(Index_t, Index_t)));


  md_data_raw.push_back     
    ( MdRecord
      ( NAME("groundNoScatteringSingleEmissivity"),
        DESCRIPTION
        (
         "Treats the ground to not cause any scattering, and to have a\n"
         "reflection coefficient of 1-e. \n"
         "\n"
         "The size of *ground_emission* is set to [ nf, stokes_dim ] where \n"
         "nf is the length of *f_grid*. Columns 2-4 are set to zero.\n"
         "The temperature of the ground is obtained by interpolating \n"
         "*t_field* to the position of the ground reflection. The obtained \n"
         "temperature and *f_grid* are then used as input to the Planck\n"
         "function. The emission from the ground is then calculated as eB,\n"
         "where B is the Planck function.\n"
         "\n"
         "It is here assumed that the downwelling radiation to consider\n"
         "comes from a single direction and the returned *ground_los*\n"
         "contains only one LOS. The slope of the ground is considered\n"
         "when calculating the LOS for the downwelling radiation. The\n"
         "reflection matrices in *ground_refl_coeffs* are all set to be\n"
         "diagonal matrices, where all diagonal elements are 1-e.\n"
         "\n"
         "Keywords: \n"
         "   e : Ground emissivity. Must be a value in the range [0,1].\n"
         "       All frequencies are assumed to have the same e."
        ),
        OUTPUT( ground_emission_, ground_los_, ground_refl_coeffs_ ),
        INPUT( f_grid_, stokes_dim_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_, 
               rte_los_, atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, 
               r_geoid_,z_ground_, t_field_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS(    "e"    ),
        TYPES(    Numeric_t )));
 
  md_data_raw.push_back     
    ( MdRecord
      ( NAME("groundTreatAsBlackbody"),
        DESCRIPTION
        (
         "Sets the ground variables (see below) to model a blackbdoy ground.\n"
         "\n"
         "The function creates the variables *ground_emission*, *ground_los*\n"
         "and *ground_refl_coeffs*. When the ground is treated to act as a\n"
         "blackbody, no downwelling radiation needs to be calculated and\n"
         "*ground_los* and *ground_refl_coeffs* are set to be empty.\n"
         "\n"
         "The size of *ground_emission* is set to [ nf, stokes_dim ] where \n"
         "nf is the length of *f_grid*. Columns 2-4 are set to zero.\n"
         "\n"
         "The temperature of the ground is obtained by interpolating \n"
         "*t_field* to the position of the ground reflection. The obtained \n"
         "temperature and *f_grid* are then used as input to the Planck\n"
         "function and the calculated blackbody radiation is put into the\n"
         "first column of *ground_emission*.\n"
         "\n"
         "Note that this function does not use *rte_los*, *r_geoid* and\n"
         "*z_ground* as input, and if used inside *ground_refl_agenda*,\n"
         "ignore commands for those variables must be added to the agenda."
        ),
        OUTPUT( ground_emission_, ground_los_, ground_refl_coeffs_ ),
        INPUT( f_grid_, stokes_dim_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_,
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, t_field_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("lsWithCutoffAdd"),
        DESCRIPTION
        (
         "Adds lineshape with cutoff to ls.\n"
         "\n"
         "The output variable ls has to exist before and has to have the correct\n"
         "size. Use this method only if you want to have a cutoff, otherwise use\n"
         "the simpler method lsAdd. What this method does is the following:\n"
         "\n"
         "0. Check if any part of f_grid is inside the cutoff. If not, return\n"
         "   immediately.\n"
         "\n"
         "1. Create els_f_grid from f_grid. Only frequencies inside the cutoff\n"
         "   are used, but the cutoff frequency itself is also added. Note that\n"
         "   els_f_grid is relative to ls_f0.\n"
         "\n"
         "2. Use els_agenda to compute a lineshape. (This\n"
         "   should be something very simple, e.g., Lorentz or Voigt.)\n"
         "\n"
         "3. Find the generated lineshape in WSV els.\n"
         "\n"
         "4. Subtract the value at the cutoff.\n"
         "\n"
         "5. Add to ls.\n"
         "\n"
         "Output:\n"
         "   ls         : The lineshape added to previous content.\n"
         "   els        : The output of els_agenda.\n"
         "                (Comunication variable with called agenda.)\n"
         "   els_f_grid : The frequency grid for els_agenda.\n"
         "                (Comunication variable with called agenda.)\n"
         "\n"
         "Input:\n"
         "   ls         : The lineshape to add to.\n"
         "   els_agenda : Used to comute lineshape.\n"
         "   ls_cutoff  : Cutoff frequency (must be > 0).\n"
         "   ls_f0      : Line center frequency (may be < 0).\n"
         "   ls_gamma   : The pressure broadening parameter (must be > 0).\n"
         "   ls_sigma   : The Doppler broadening parameter (must be > 0).\n"
         "   f_grid     : Global frequency grid."
        ),
        OUTPUT( ls_,
                els_,
                els_f_grid_ ),
        INPUT(  ls_,
                els_agenda_,
                ls_cutoff_,
                ls_f0_, ls_gamma_, ls_sigma_, f_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Main"),
        DESCRIPTION
        ( 
         "Run the agenda that is specified inside the curly braces. ARTS\n"
         "controlfiles must define this method. It is executed automatically\n"
         "when ARTS is run on the controlfile." 
        ),
        OUTPUT( ),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES(),
        AGENDAMETHOD( true )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("MatrixCBR"),
        DESCRIPTION
        (
         "Sets a matrix to hold cosmic background radiation (CBR).\n"
         "\n"
         "The CBR is assumed to be un-polarized and Stokes components 2-4\n"
         "are zero. Number of Stokes components, that equals the number \n"
         "of columns in the created matrix, is determined by *stokes_dim. \n"
         "The number of rows in the created matrix equals the length of the \n"
         "given frequency vector. \n"
         "\n"
         "The cosmic radiation is modelled as blackbody radiation for the \n"
         "temperature given by the global constant COSMIC_BG_TEMP, set in \n"
         "the file constants.cc. The frequencies are taken from the generic \n"
         "input vector. A standard way to use the function should be (as \n"
         "part of *i_space_agenda*): \n"
         "   MatrixCBR(i_space,f_grid){} \n"
         "\n"
         "Generic output: \n"
         "   Matrix : Matrix with cosmic background radiation. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A set of frequencies. "
        ),
        OUTPUT(),
        INPUT( stokes_dim_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Matrix1ColFromVector"),
        DESCRIPTION
        (
         "Forms a matrix containing 1 column from a vector.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Matrix2ColFromVectors"),
        DESCRIPTION
        (
         "Forms a matrix containing 2 columns from two vectors.\n"
         "\n"
         "The vectors are put as columns in the matrix in the same order\n"
         "as they are given.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied into the first column. \n"
         "   Vector : The vector to be copied into the second column."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Vector_ ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Matrix3ColFromVectors"),
        DESCRIPTION
        (
         "Forms a matrix containing 3 columns from three vectors.\n"
         "\n"
         "The vectors are put as columns in the matrix in the same order\n"
         "as they are given.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied into the first column. \n"
         "   Vector : The vector to be copied into the second column. \n"
         "   Vector : The vector to be copied into the third column."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Matrix1RowFromVector"),
        DESCRIPTION
        (
         "Forms a matrix containing 1 row from a vector.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_ ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Matrix2RowFromVectors"),
        DESCRIPTION
        (
         "Forms a matrix containing 2 rows from two vectors.\n"
         "\n"
         "The vectors are put as rows in the matrix in the same order\n"
         "as they are given.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied into the first row. \n"
         "   Vector : The vector to be copied into the second row."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Vector_ ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Matrix3RowFromVectors"),
        DESCRIPTION
        (
         "Forms a matrix containing 3 rows from three vectors.\n"
         "\n"
         "The vectors are put as rows in the matrix in the same order\n"
         "as they are given.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied into the first row. \n"
         "   Vector : The vector to be copied into the second row. \n"
         "   Vector : The vector to be copied into the third row."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_, Vector_, Vector_ ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixPlanck"),
        DESCRIPTION
        (
         "Sets a matrix to hold blackbody radiation.\n"
         "\n"
         "This function works as MatrixCBR but the temperature for which \n"
         "(unpolarised) balckbody radiation shall be calculated is selected \n"
         "as a keyword argument.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : Matrix with cosmic background radiation. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A set of frequencies. \n"
         "\n"
         "Keyword: \n"
         "   t : Temperature for the balckbody radiation. "
        ),
        OUTPUT(),
        INPUT( stokes_dim_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "t" ),
        TYPES(    Numeric_t )));


  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixScale"),
        DESCRIPTION
        (
         "Scales all elements of a matrix with the same value. \n"
         "\n"
         "The result can either be stored in the same or another matrix. \n"
         "\n"
         "Generic output: \n"
         "   Matrix : Return matrix. \n"
         "\n"
         "Generic input: \n"
         "   Matrix : Original matrix. \n"
         "\n"
         "Keywords: \n"
         "   value : The value to be multiplicated with the matrix."  
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_ ),
        KEYWORDS( "value" ),
        TYPES(    Numeric_t   )));

   md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixSet"),
        DESCRIPTION
        (
         "Creates a workspace matrix with the specified size and sets \n"
         "all values of the matrix to the specified value. \n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Keywords:\n"
         "   nrows : The number of rows of the matrix to create. \n"
         "   ncols : The number of columns of the matrix to create. \n"
         "   value : The value of the matrix elements. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT(),
        KEYWORDS( "nrows", "ncols", "value"   ),
        TYPES(    Index_t, Index_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixSetTakingSizeFromMatrix"),
        DESCRIPTION
        (
         "Creates a workspace vector with the same size as another matrix,\n"
         "and sets all values of the new matrix to the specified value. \n"
         "\n"
         "Generic output: \n"
         "   Matrix : The matrix to be created. \n"
         "\n"
         "Generic input: \n"
         "   Matrix : The matrix specifying the size. \n"
         "\n"
         "Keywords:\n"
         "   value  : The value of the matrix elements. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_ ),
        KEYWORDS( "value"   ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixToTbByPlanck" ),
        DESCRIPTION
        (
         "Converts a matrix of radiances to brightness temperatures by \n"
         "inverting the Planck function.\n"
         "\n"
         "This function works as *MatrixToTbByRJ*. However, this function \n"
         "is not recommended in connection with inversions, but can be used \n"
         "to display calculated spectra in a temperature scale.\n"
         "\n"
         "Generic output: \n"
         "   Matrix : A matrix with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Matrix : A matrix with radiance values."
        ),
        OUTPUT(),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_, sensor_pol_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_ ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "MatrixToTbByRJ" ),
        DESCRIPTION
        (
         "Converts a matrix of radiances to brightness temperatures by \n"
         "the Rayleigh-Jeans approximation of the Planck function.\n"
         "\n"
         "This function works as *VectorToTbByRJ*, but operates on a matrix.\n"
         "Each column of the matrix is treated to contain a spectral vector,\n"
         "with frequencies repeated as assumed in *VectorToTbByRJ*. \n"
         "\n"
         "Generic output: \n"
         "   Matrix : A matrix with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Matrix : A matrix with radiance values."
        ),
        OUTPUT(),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_, sensor_pol_ ),
        GOUTPUT( Matrix_ ),
        GINPUT( Matrix_ ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "montecarlo_p_from_belowCscaAdapt" ),
        DESCRIPTION
        (
         "Reduces montecarlo_p_from_belowCscaAdapt to one frequency" 
        ),
        OUTPUT(montecarlo_p_from_belowCsca_),
        INPUT(f_grid_,f_index_,scat_data_raw_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("NumericSet"),
        DESCRIPTION
        (
         "Sets a numeric workspace variable to the given value. \n"
         "\n"
         "Generic output: \n"
         "   Numeric : The numeric variable to be set. \n"
         "\n"
         "Keywords:\n"
         "   value : The value." 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Numeric_ ),
        GINPUT(),
        KEYWORDS( "value"   ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("output_file_formatSetAscii"),
        DESCRIPTION
        (
         "Sets the output file format to ASCII."
        ),
        OUTPUT( output_file_format_),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("output_file_formatSetBinary"),
        DESCRIPTION
        (
         "Sets the output file format to binary."
        ),
        OUTPUT( output_file_format_),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("opt_prop_sptFromData"),
        DESCRIPTION
        (
         "Calculates opticle properties for the single particle types.\n"
         "\n"
         "In this function extinction matrix and absorption vector are \n"
         "calculated in the laboratory frame. These properties are required\n"
         "for the RT calculation, inside the the i_fieldUpdateXXX \n"
         "functions.\n" 
         "\n"
         "The interpolation of the data on the actual frequency is the \n"
         "first step in this function. \n"
         "\n"
         "Then the transformation from the database coordinate system to to \n"
         "laboratory coordinate system is done.\n"
         "\n"
         "Output of the function are *ext_mat_spt*, and *abs_vec_spt* which\n"
         "hold the optical properties for a specified propagation direction\n"
         "for each particle type. \n"
         "\n"
        ),
        OUTPUT( ext_mat_spt_, abs_vec_spt_ ),
        INPUT(  ext_mat_spt_, abs_vec_spt_, scat_data_raw_,
               scat_za_grid_, 
               scat_aa_grid_, scat_za_index_, scat_aa_index_, 
               f_index_, f_grid_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "ParticleTypeAdd" ),
        DESCRIPTION
        (
         "This method reads single scattering data and particle number\n"
         "density fields from a data base. \n"
         "\n"
         "The method allows the user to chose hydro-meteor species and \n"
         "particle \n"
         "number density fields. The methods reads from the chosen files \n"
         "and appends the variables *scat_data_raw* and *pnd_field_raw*. \n"
         "There is one database for particle number density fields ( ....),\n"
         "which includes the following particle types:\n"
         "\n"
         "Another database (....) containes the scattering properties for \n"
         "those particle types. \n"
         ),
        OUTPUT(scat_data_raw_, pnd_field_raw_),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("filename_scat_data", "filename_pnd_field"),
        TYPES(String_t, String_t)));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "ParticleTypeAddAmpl" ),
        DESCRIPTION
        (
         "This method reads the amplitute matrix and the particle number\n"
         "density field from a data base. \n"
         "\n"
         "The method allows the user to chose particle types and particle \n"
         "number density fields. The methods reads from the chosen files \n"
         "and appends the variables *amp_mat_raw* and *pnd_field_raw*. \n"
         "There is one database for particle number density fields ( ....),\n"
         "which includes the following particle types:\n"
         "\n"
         "Another database (....) containes the amplitude matrices for \n"
         "those particle types from which all optical properties can be \n"
         "derived.\n"
         ),
        OUTPUT(amp_mat_raw_, pnd_field_raw_),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("filename_amp_mat", "filename_pnd_field"),
        TYPES(String_t, String_t)));

 md_data_raw.push_back
   ( MdRecord
     ( NAME( "ParticleTypeInit" ),
       DESCRIPTION
       (
        "This method initializes variables containing data about the \n"
        "optical properties of particles (*scat_data_raw*) and about the \n"
        "particle number distribution (*pnd_field_raw*)\n"
        "\n"
        "*ParticleTypeInit* has to be executed before executing \n"
        "*ParticleTypeAdd*.\n"
        ),
       OUTPUT(scat_data_raw_, pnd_field_raw_),
       INPUT(),
       GOUTPUT(),
       GINPUT(),
       KEYWORDS(), 
       TYPES())); 


   md_data_raw.push_back
    ( MdRecord
      ( NAME( "ParticleTypeInitAmpl" ),
        DESCRIPTION
        (
         "This method initializes variables containing data about the \n"
         "optical properties of particles (*amp_mat_raw*) and about the \n"
         "particle number distribution (*pnd_field_raw*)\n"
         "\n"
         "*ParticleTypeInit* has to be executed before executing \n"
         "*ParticleTypeAdd*.\n"
         ),
        OUTPUT(amp_mat_raw_, pnd_field_raw_),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(), 
        TYPES())); 

   md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_matCalc" ),
        DESCRIPTION
        (
         "This function sums up the phase matrices for all particle\n"
         "types weighted with particle number density.\n"
         "\n"
         "The output of this method is *pha_mat* (Nza, Naa, stokes_dim,\n"
         "stokes_dim). The inputs are the phase matrix for the single particle\n"
         "type *pha_mat_spt* (part_types, Nza, Naa, stokes_dim, stokes_dim)\n"
         "and the local particle  number densities for all particle types namely \n"
         "the *pnd_field* (part_types, p_grid, lat_grid, lon_grid ) for given\n"
         "*p_grid*, *lat_grid*, and *lon_grid*. The particle types required \n"
         "are specified in the control file.\n"
         ),
        OUTPUT(pha_mat_),
        INPUT(pha_mat_spt_, pnd_field_, atmosphere_dim_, scat_p_index_,
              scat_lat_index_, scat_lon_index_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_mat_sptCalc" ),
        DESCRIPTION
        (
         "This method calculates the phase matrix for a single particle type.\n"
         "\n"
         "Phase matrix explains the tranformation of Stokes parameters of \n"
         "incident plane wave into those of the scattered spherical wave \n"
         "due to scattering of radiation by the particle. All the elements \n"
         "of phase matrix can be expressed in terms of the elements of the\n"
         "amplitude matrix.\n"
         "\n"
         "The output of the method pha_mat_sptCalc is pha_mat_spt(Tensor 5,\n"
         "size: [Npt,Nza,Naa,stokes_dim,stokes_dim]). The input to the method\n"
         "pha_mat_sptCalc are *pha_mat_spt*,*amp_mat*(Tensor 6,\n"
         "Size=[Npt,Nza,Naa,Nza,Naa,8]), *scat_za_index* and scat_aa_index.\n"
         "\n"
         "The variables *scat_za_index* and *scat_aa_index picks the right \n"
         "element of the Tensor *amp_mat*.Then this method calls the \n"
         "function amp2pha which does the actual physics, that of computing\n"
         "the elements of phase matrix from the elements of amplitude \n"
         "matrix."
         ),
        OUTPUT(pha_mat_spt_),
        INPUT(pha_mat_spt_, amp_mat_, scat_za_index_, scat_aa_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 
  
  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_mat_sptFromData" ),
        DESCRIPTION
        (
         "Calculation of the phase matrix for the single particle types.\n"
         "\n"
         "In this function the phase matrix calculated \n"
         "in the laboratory frame. This function is used in the calculation\n"
         "of the scattering integral (*scat_fieldCalc*).\n"
         "\n"
         "The interpolation of the data on the actual frequency is the first\n"
         "step in this function. \n"
         "\n"
         "Then the transformation from the database coordinate system to to \n"
         "laboratory coordinate system is done.\n"
         "\n"
         ),
        OUTPUT(pha_mat_spt_),
        INPUT(pha_mat_spt_, scat_data_raw_, scat_za_grid_, scat_aa_grid_, 
              scat_za_index_, scat_aa_index_, f_index_, f_grid_, scat_theta_,
              scat_theta_gps_, scat_theta_itws_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "pha_mat_sptFromDataDOITOpt" ),
        DESCRIPTION
        (
         "Calculation of the phase matrix for the single particle types.\n"
         "\n"
         "In this function the phase matrix calculated \n"
         "in the laboratory frame. This function can be used in the agenda\n"
         "*pha_mat_spt_agenda* which is part of the calculation of the scattering\n"
         "integral (*scat_fieldCalc*).\n"
         "\n"
         "The function uses *pha_mat_sptDOITOpt*, which is the phase matrix \n"
         "interpolated\n"
         "on the right frequency for all scattering angles.\n"
         "\n"
         "The transformation from the database coordinate system to to \n"
         "laboratory coordinate system is done in this function.\n"
         "\n"
         "NOTE: This is a special function for 1 particle type, including \n"
         "randomly oriented particles. \n"
         "\n"
         "In this function the phase matrix calculated \n"
         "in the laboratory frame. This function is used in the calculation\n"
         "of the scattering integral (*scat_fieldCalc*).\n"
         "\n"
         ),
	OUTPUT(pha_mat_spt_),
        INPUT(pha_mat_spt_, pha_mat_sptDOITOpt_, scat_theta_, za_grid_size_,
              scat_aa_grid_, 
              scat_za_index_, scat_aa_index_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "pnd_fieldCalc" ),
        DESCRIPTION
        ("Interpolate the particle number density fields.\n"
         "\n"
         "This methods interpolates the particle number density field\n"
         "from the raw data *pnd_field_raw* to pnd_field* onto the grids for the\n"
         "calculation namely, *p_grid*, *lat_grid*, *lon_grid*.  The method is\n"
         "similar to *AtmFieldsCalc* where temperature, altitude and vmr are\n"
         "interpolated onto *p_grid*, *lat_grid*, *lon_grid*. \n"
         "\n"
         "The method takes in as input *pnd_field_raw* which is an\n"
         "ArrayOfArrayOfTensor3 which contains one gridded field for each\n"
         "particle type. See the online documentation of *pnd_field_raw* for\n"
         "more information about this variable.  The output *pnd_field* is a\n"
         "Tensor4 with dimension pnd_field(part_types, p_grid, lat_grid,\n"
         "lon_grid). \n"
         ),
        OUTPUT(pnd_field_),
        INPUT(p_grid_, lat_grid_, lon_grid_, pnd_field_raw_, atmosphere_dim_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES())); 

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ppathCalc" ),
        DESCRIPTION
        (
         "Main function for calculation of propagation paths.\n"
         "\n"
         "There exists only one function to calculate total propagation\n"
         "paths and this is that function. The function is normally not\n"
         "visible in the control file, it is called from inside *RteCalc*.\n"
         "A reason to call this function directly would be to plot a\n"
         "propgation path.\n"
         "\n"
         "The definition of a propgation path cannot be accomodated here.\n"
         "For more information read the chapter on propagation paths in the\n"
         "ARTS user guide and read the  on-line information for\n"
         "*ppath_step_agenda* (type \"arts -d ppath_step_agenda\")."
        ),
        OUTPUT( ppath_, ppath_step_ ),
        INPUT( ppath_step_agenda_, atmosphere_dim_, p_grid_, lat_grid_, 
               lon_grid_, z_field_, r_geoid_, z_ground_, 
               cloudbox_on_, cloudbox_limits_, rte_pos_, rte_los_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));


  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ppath_stepGeometric" ),
        DESCRIPTION
        (
         "Calculates a geometrical propagation path step.\n"
         "\n"
         "This function determines a propagation path step by pure\n"
         "geometrical calculations. That is, refraction is neglected. Path\n"
         "points are always included for crossings with the grids, tangent\n"
         "points and points of ground intersections. The keyword *lmax* \n"
         "gives the option to include additional points to ensure that the\n"
         "distance along the path between the points does not exceed the \n"
         "selected maximum length. No additional points are included if\n"
         "*lmax* is set to <= 0.\n"
         "\n"
         "As functions of this kind should very seldom be called directly,\n"
         "and that the functions can be called a high number of times, these\n"
         "functions do not perform any checks of the input that give\n" 
         "detailed error messages, but asserts are performed (if turned on).\n"
         "\n"
         "For further information, type see the on-line information for\n"
         "*ppath_step_agenda* (type \"arts -d ppath_step_agenda\").\n"
         "\n"
         "Keywords: \n"
         "   lmax      : Maximum allowed length between path points."
        ),
        OUTPUT( ppath_step_ ),
        INPUT( ppath_step_, atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, 
               z_field_, r_geoid_, z_ground_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "lmax" ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ppath_stepRefractionEuler" ),
        DESCRIPTION
        (
         "Calculates a propagation path step, considering refraction by a\n"
         "straightforward Euler approach.\n"
         "\n"
         "Refraction is taken into account by probably the simplest approach\n"
         "possible. The path is treated to consist of piece-wise geometric \n"
         "steps. A geometric path step is calculated from each point by \n"
         "using the local line-of-sight. Except for 1D zenith angles, the\n"
         "path quantities are propagated by solving the differential \n"
         "equations by the Euler method. Snell's law for a case with \n"
         "spherical symmetry is used for 1D to update the zenith angles. \n"
         "\n"
         "See further the on-line information for *ppath_stepGeometric*\n"
         "(type \"arts -d ppath_stepGeometric\") and the user guide for more\n"
         "details on the algorithms used.\n"
         "\n"
         "The maximum length of each ray tracing step is given by the \n"
         "keyword argument *lraytrace*. The length will never exceed the \n" 
         "given maximum value, but can be smaller. The ray tracing steps are\n"
         "only used to determine the path, points to describe the path for \n" 
         "*RteCalc* are included as for *ppath_stepGeometric*, this\n"
         "including the functionality for the keyword *lmax*.\n"
         "\n"
         "Keywords: \n"
         "   lraytrace : Maximum length of ray tracing steps.\n"
         "   lmax      : Maximum allowed length between path points."
        ),
        OUTPUT( ppath_step_, rte_pressure_, rte_temperature_, rte_vmr_list_, 
                refr_index_ ),
        INPUT( refr_index_agenda_, ppath_step_, atmosphere_dim_, p_grid_, 
               lat_grid_, lon_grid_, z_field_, t_field_, vmr_field_, r_geoid_,
               z_ground_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "lraytrace", "lmax"    ),
        TYPES(    Numeric_t,   Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ReadXML"),
        DESCRIPTION
        (
         "Reads a workspace variable from an XML file.\n"
         "\n"
         "This is a supergeneric method. It can read variables of any group.\n"
         "\n"
         "If the filename is omitted, the variable is read\n"
         "from <basename>.<variable_name>.xml.\n"
         "\n"
         "Usage example:\n"
         "\n"
         "ReadXML(f_grid){\"frequencies.xml\"}\n"
         "Will read the frequency grid *f_grid* from the specified file.\n"
         "\n"
         "Supergeneric output:\n"
         "   Any_     : The variable to read.\n"
         "\n"
         "Keywords:\n"
         "   filename : Name of the input file."
         ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Any_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("refr_indexFieldAndGradients"),
        DESCRIPTION
        (
         "Calculates the field and gradients of the refractive index.\n"
         "\n"
         "This function calculates the refractive index and its gradients\n"
         "for a rectangular grid. \n"
         "\n"
         "Calculations are performed for all combinations of the given \n"
         "vectors, where the first vector shall contain pressure values, the\n"
         "second latitude values, and the last longitude values. For \n"
         "dimensions not used, the corresponding position vector is ignored.\n"
         "\n"
         "The calculated values for a Tensor4, with size:\n"
         "   [atmosphere_dim+1, np, nlat, nlon] \n"
         "where np is the number of pressures given etc. The book of the\n"
         "tensor with the following index holds:\n"
         "   0: the refractive index \n"
         "   1: radial gradient of the refractive index \n"
         "   2: latitude gradient of the refractive index \n"
         "   3: longitude gradient of the refractive index \n"
         "\n"
         "To calculate these quantities for the atmsopheric mesh, execute:\n"
         "   RefrIndexFieldAndGradients(tensor4_1,p_grid,lat_grid,lon_grid)"
        ),
        OUTPUT( refr_index_, rte_pressure_, rte_temperature_, rte_vmr_list_ ),
        INPUT( refr_index_agenda_, atmosphere_dim_, p_grid_, lat_grid_, 
               lon_grid_, r_geoid_, z_field_, t_field_, vmr_field_ ),
        GOUTPUT( Tensor4_ ),
        GINPUT( Vector_, Vector_, Vector_  ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("refr_indexIR"),
        DESCRIPTION
        (
         "Calculates the microwave refractive index due to gases in the\n"
         "Earth's atmosphere. \n"
         "\n"
         "Only refractivity of dry air is considered. All other gases has\n"
                 "a negligible contribution.  \n"
         "\n"
         "The formula used is contributed by Michael Höpfner,\n"
                 "Forschungszentrum Karlsruhe."
        ),
        OUTPUT( refr_index_ ),
        INPUT( rte_pressure_, rte_temperature_, rte_vmr_list_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("refr_indexThayer"),
        DESCRIPTION
        (
         "Calculates the microwave refractive index due to gases in the\n"
         "Earth's atmosphere. \n"
         "\n"
         "The refractivity of dry air and water vapour is summed. All\n"
         "other gases has a negligible contribution.  \n"
         "\n"
         "The parameterisation of Thayer (Radio Science, 9, 803-807, 1974)\n"
         "is used. See also Eq. 3 and 5 of Solheim et al. (JGR, 104, \n"
         "pp. 9664). "
        ),
        OUTPUT( refr_index_ ),
        INPUT( rte_pressure_, rte_temperature_, rte_vmr_list_, gas_species_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("refr_indexUnit"),
        DESCRIPTION
        (
         "Sets the refractive index to 1.\n"
         "\n"
         "If this method is used, the obtained path should be identical to\n"
         "the geomtrical path.\n"
         "\n"
         "As this function does not need any input, you have to include call\n"
         "of *Ignore* for all variables expected to be used by\n"
         "*refr_index_agenda*."
        ),
        OUTPUT( refr_index_ ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "RteCalc" ),
        DESCRIPTION
        (
         "Main function for calculation of spectra.\n"
         "\n"
         "Spectra are calculated in a fixed manner, by this function.\n"
         "The function calculates monochromatic spectra for all pencil beam\n"
         "directions and applies the sensor response WSVs on obtained\n"
         "values.\n"
         "\n"
         "The first step is to calculate the propagation path through the\n"
         "atmosphere for the considered viewing direction. The next step is\n"
         "to determine the spectrum at the starting point of the propagation\n"
         "path. The start point of the propagation path can be found at the\n"
         "top of the atmosphere, the ground, or at the boundary or inside\n"
         "the cloud box. To determine the start spectrum can involve a\n"
         "recursive call of *RteCalc (for example to calculate the radiation\n"
         "refelected by the ground). After this, the vector radiative\n"
         "transfer equation is solved to the end point of the propagation\n"
         "path. Finally, the polarisation and intensity response of the \n"
         "sensor are applied.\n"
         "\n"
         "See further the user guide."
        ),
        OUTPUT( y_, ppath_, ppath_step_, i_rte_,
                rte_pos_, rte_los_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_, 
                i_space_, ground_emission_, ground_los_, ground_refl_coeffs_ ),
        INPUT( ppath_step_agenda_, rte_agenda_, i_space_agenda_,
               ground_refl_agenda_,
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, z_field_, 
               t_field_, r_geoid_, z_ground_, cloudbox_on_, cloudbox_limits_, 
               scat_i_p_, scat_i_lat_, scat_i_lon_, 
               scat_za_grid_, scat_aa_grid_, sensor_response_,
               sensor_pos_, sensor_los_, sensor_pol_, sensor_rot_, 
               f_grid_, stokes_dim_,
               antenna_dim_, mblock_za_grid_, mblock_aa_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "RteEmissionStd" ),
        DESCRIPTION
        (
         "Standard RTE function with emission.\n"
         "\n"
         "This function does a clearsky radiative transfer calculation for\n"
         "a given propagation path. Designed to be part of *rte_agenda*.\n"
         "\n"
         "Gaseous emission and absorption is calculated for each propagation\n"
         "path point using the agenda *gas_absorption_agenda*. \n"
         "Absorption vector and extinction matrix are created using \n" 
         "*opt_prop_part_agenda*.\n"
         "Absorption and extinction variables are averaged between two\n" 
         "successive propagation path points. The same applies to the\n"
         "Planck function values.\n"
         "\n"
         "See futher the user guide."
        ),
        OUTPUT( i_rte_, abs_vec_, ext_mat_, rte_pressure_, rte_temperature_,
                rte_vmr_list_, f_index_ ),
        INPUT( i_rte_, ppath_, f_grid_, stokes_dim_, 
               atmosphere_dim_, p_grid_, lat_grid_, lon_grid_, t_field_,
               vmr_field_, scalar_gas_absorption_agenda_, 
               opt_prop_gas_agenda_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "r_geoidSpherical" ),
        DESCRIPTION
        (
         "Sets the geoid to be a perfect sphere.\n"
         "\n"
         "The radius of the sphere is selected by the keyword argument *r*.\n"
         "If the keyword is set to be negative, the radius is set to the\n"
         "global internal variable *EARTH_RADIUS*, defined in constants.cc.\n"
         "\n"
         "Keywords:\n"
         "   r : Radius of geoid sphere. See further above."
        ),
        OUTPUT( r_geoid_ ),
        INPUT( atmosphere_dim_, lat_grid_, lon_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS( "r" ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "r_geoidWGS84" ),
        DESCRIPTION
        (
         "Sets the geoid radius to match the WGS-84 reference ellipsiod.\n"
         "\n"
         "For 1D, the geoid radius is set to the radius of curvature of the\n"
         "WGS-84 ellipsiod for the position and observation direction \n"
         "described with *lat_1d* and *meridian_angle_1d*.\n"
         "For 2D and 3D, *r_geoid* is set to the radius of the WGS-84\n"
         "ellipsiod for the crossing points of the latitude and longitude\n"
         "grids.\n"
         "\n"
         "Please note that the latitude grid must contain true latitudes\n"
         "if the function shall give correct result, and not just arbitrary\n"
         "orbit angles which is allowed for 2D cases."
        ),
        OUTPUT( r_geoid_ ),
        INPUT( atmosphere_dim_, lat_grid_, lon_grid_, lat_1d_,
               meridian_angle_1d_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_posShift" ),
        DESCRIPTION
        (
         "shifts rte_pos and rte_los, and rte_gp_XXX to the end of ppath."
        ),
        OUTPUT( rte_pos_, rte_los_, rte_gp_p_, rte_gp_lat_, rte_gp_lon_ ),
        INPUT( ppath_, atmosphere_dim_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "rte_pos_and_losFromTangentPressure" ),
        DESCRIPTION
        (
         "If you are doing limb calculations it can be useful to specify\n"
         "viewing direction and sensor position by the tangent pressure.\n"
         "This function takes tan_p as a keyword argument and sets rte_los\n"
         "and rte_pos to the apropriate position on the edge of the modelled\n"
         "atmosphere\n\n"
         "This function is a work in progress. Only 1D is currently supported"
        ),
        OUTPUT( rte_pos_, rte_los_, ppath_, ppath_step_ ),
        INPUT( atmosphere_dim_, p_grid_, z_field_, lat_grid_, lon_grid_,
               ppath_step_agenda_, r_geoid_, z_ground_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("tan_p"),
        TYPES(Numeric_t)));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "scat_data_rawCheck" ),
        DESCRIPTION
        (
         "Method for checking the consistency of the optical properties\n"
         "in the database. \n"
         "\n"
         "This function can be used to check datafiles containing data for \n"
         "randomly oriented scattering media.\n"
         "It is checked whether the data is consistent. The integral over \n"
         "the phase matrix should result the scattering cross section \n"
         "<C_sca>.\n"
         "\n"
         "The check is if:\n"
         "<C_ext> - <C_sca> = <C_abs>\n"
         "\n"
         "The result is printed on the screen.\n"
         "\n"
         ),
        OUTPUT( ),
        INPUT( scat_data_raw_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "scat_fieldCalc" ),
        DESCRIPTION
        (
         "This method calculates the scattering integral.\n"
         "\n"
         "By scattering integral we mean the field generated by integrating\n"
         "the product of intensity field and phase matrix over all incident \n"
         "angles. This term results solely due to the scattering properties\n"
         "of ice particles in cirrus clouds.  \n"
         "\n"
         "The output of this method is the scattered field *scat_field*\n"
         "which is used in the radiative transfer equation to give a new\n"
         "radiation field *i_field*. The dimensions of *scat_field* and \n"
         "*i_field* are the same. This resultant field is again given as \n"
         "input to this method, which calculates a new *scat_field*.  The \n"
         "iteration continues till the field converges.  Another important\n"
         "requirement for this method is the phase matrix.  For this we \n"
         "give as input to this method *pha_mat_spt* and *pnd_field*. From\n"
         "these two workspace variables we calculate *pha_mat* with\n"
         "the method *pha_matCalc*.  \n"
         ),
        OUTPUT( scat_field_, pha_mat_, pha_mat_spt_, scat_za_index_,
                scat_aa_index_),
        INPUT( pha_mat_spt_agenda_, i_field_, pnd_field_, scat_za_grid_, 
               scat_aa_grid_, p_grid_, lat_grid_, lon_grid_, 
               atmosphere_dim_, cloudbox_limits_, za_grid_size_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

md_data_raw.push_back
    ( MdRecord
      ( NAME( "scat_fieldCalcLimb" ),
        DESCRIPTION
        (
         "This method calculates the scattering integral (Limb).\n"
         "\n"
         "By scattering integral we mean the field generated by integrating\n"
         "the product of intensity field and phase matrix over all incident \n"
         "angles. This term results solely due to the scattering properties\n"
         "of ice particles in cirrus clouds.  \n"
         "\n"
         "The output of this method is the scattered field *scat_field*\n"
         "which is used in the radiative transfer equation to give a new\n"
         "radiation field *i_field*. If different zenith angle grids are\n"
         "required for the calculation of the scattering integral and the  \n"
         "radiative transfer (this is the case for Limb observations) this \n"
         "method has to be used instead of *scat_fieldCalc*.\n"
         "\n"
         "The dimensions of *scat_field* and \n"
         "*i_field* are the same. This resultant field is again given as \n"
         "input to this method, which calculates a new *scat_field*.  The \n"
         "iteration continues till the field converges.  Another important\n"
         "requirement for this method is the phase matrix.  For this we \n"
         "give as input to this method *pha_mat_spt* and *pnd_field*. From\n"
         "these two workspace variables we calculate *pha_mat* with\n"
         "the method *pha_matCalc*.  \n"
         ),
        OUTPUT( scat_field_, pha_mat_, pha_mat_spt_, scat_za_index_,
                scat_aa_index_),
        INPUT( pha_mat_spt_agenda_, i_field_, pnd_field_, scat_za_grid_,
               scat_aa_grid_, p_grid_, lat_grid_, lon_grid_, 
               atmosphere_dim_, cloudbox_limits_, za_grid_size_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));


  md_data_raw.push_back
    ( MdRecord
      ( NAME( "scat_fieldCalcFromAmpMat" ),
        DESCRIPTION
        (
         "This method calculates the scattering integral.\n"
         "\n"
         "By scattering integral we mean the field generated by integrating\n"
         "the product of intensity field and phase matrix over all incident \n"
         "angles. This term results solely due to the scattering properties\n"
         "of ice particles in cirrus clouds.  \n"
         "\n"
         "The output of this method is the scattered field *scat_field*\n"
         "which is used in the radiative transfer equation to give a new\n"
         "radiation field *i_field*. The dimensions of *scat_field* and \n"
         "*i_field* are the same. This resultant field is again given as \n"
         "input to this method, which calculates a new *scat_field*.  The \n"
         "iteration continues till the field converges.  Another important\n"
         "requirement for this method is the phase matrix.  For this we \n"
         "give as input to this method *pha_mat_spt* and *pnd_field*. From\n"
         "these two workspace variables we calculate *pha_mat* with\n"
         "the method *pha_matCalc*.  \n"
         ),
        OUTPUT( scat_field_, pha_mat_, pha_mat_spt_ ),
        INPUT( amp_mat_, i_field_, pnd_field_, scat_za_grid_, 
               scat_aa_grid_, p_grid_, lat_grid_, lon_grid_, 
               atmosphere_dim_, cloudbox_limits_,
               grid_stepsize_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "scat_iPut" ),
        DESCRIPTION
        (
         "Method for the communication between cloudbox and clearsky.\n"
         "\n"
         "This method puts the scattered radiation field into the interface\n"
         "variables between the cloudbox and the clearsky, which are \n"
         "*scat_i_p*, *scat_i_lat* and *scat_i_lon*. As i_field is only\n"
         "stored for one frequency given by *f_index* this method has\n" 
         "to be\n"
         "executed after each scattering calculation to store the scattered\n"
         "field on the boundary of the cloudbox.\n"
         "\n"
         "The best way to calculate spectra including the influence of\n" 
         "scattering is to set up the *scat_mono_agenda* where this method \n"
         "can be included.\n"
         "\n"
         ),
        OUTPUT( scat_i_p_, scat_i_lat_, scat_i_lon_ ),
        INPUT( i_field_, f_grid_, f_index_,   p_grid_, lat_grid_, 
               lon_grid_, scat_za_grid_, scat_aa_grid_, stokes_dim_,
               atmosphere_dim_, cloudbox_limits_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

md_data_raw.push_back
    ( MdRecord
      ( NAME( "scat_iPutMonteCarlo" ),
        DESCRIPTION
        (
         "Method for the communication between cloudbox and clearsky.\n"
         "\n"
         "This is the equivalent of scat_iPut for use with ScatteringMonteCarlo.\n"
         "To fit in with pre-existing code the Stokes vector I is simply copied \n"
         "several times to make the sizes of scat_i_p, scat_i_lat, and scat_i_lon\n"
         "the same as they would be of successive order of scattering calculations.\n"
         "This means that after *scat_iPutMonteCarlo* the radiative transfer \n"
         "calculation can be completed by simply calling *RteCalc*\n"
         "\n"
         ),
        OUTPUT( scat_i_p_, scat_i_lat_, scat_i_lon_ ),
        INPUT( i_rte_, stokes_dim_, f_grid_, cloudbox_limits_, scat_za_grid_,
               scat_aa_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
   ( MdRecord
      ( NAME( "ScatteringInit" ),
        DESCRIPTION
        (
         "Initialize variables for a scattering calculation. \n"
         "\n"
         "Variables needed in the scattering calculations are initialzed\n"
         "here. This method has to be executed before using \n"
         "*ScatteringMain*.\n"
         "\n"
         ),
        OUTPUT(scat_p_index_, scat_lat_index_, scat_lon_index_, 
               scat_za_index_, scat_aa_index_, pha_mat_,
               pha_mat_spt_, ext_mat_spt_, abs_vec_spt_, scat_field_,
               i_field_),
        INPUT(stokes_dim_, atmosphere_dim_, scat_za_grid_, scat_aa_grid_,
              za_grid_size_, 
              cloudbox_limits_, scat_data_raw_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
   ( MdRecord
      ( NAME( "ScatteringInitAmpMat" ),
        DESCRIPTION
        (
         "Initialize variables for a scattering calculation. \n"
         "\n"
         "This function has to be used, when the scattering data is stored \n"
         "in amplitude matrix format. \n"
         "Variables needed in the scattering calculations are initialzed\n"
         "here. This method has to be executed before using \n"
         "*ScatteringMain*.\n"
         "\n"
         ),
        OUTPUT(scat_p_index_, scat_lat_index_, scat_lon_index_, 
               scat_za_index_, scat_aa_index_, iteration_counter_, pha_mat_,
               pha_mat_spt_, ext_mat_spt_, abs_vec_spt_, scat_field_,
               i_field_),
        INPUT(stokes_dim_, atmosphere_dim_, scat_za_grid_, scat_aa_grid_,
              cloudbox_limits_, amp_mat_raw_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "ScatteringMain" ),
        DESCRIPTION
        (
         "This method executes *scat_mono_agenda* for each frequency defined\n"
         "in *f_grid* \n"
         "\n"
         "If the number of frequencies is only two, it is assumed that the \n"
         "user is only interested in a monochromatic scattering calculation\n"
         "and executes the agenda only for the first frequency.\n"
         "\n"
         ),
        OUTPUT(f_index_),
        INPUT(f_grid_, scat_mono_agenda_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));
 
 md_data_raw.push_back
    ( MdRecord
      ( NAME( "ScatteringMonteCarlo" ),
        DESCRIPTION
        (
         "This method performs a single pencil beam monochromatic scattering\n"
         "calculation using a Monte Carlo algorithm \n"
         "\n"
         "\n"
         "The main output variables *i_rte* and *i_montecarlo_error* represent the \n"
         "Stokes vector, and the estimated error in the Stokes vector respectively.\n"
         "The keyword parameter `maxiter\' describes the number of `photons\'\n"
         "used in the simulation (more photons means smaller *Ierror*).\n"
         "Non-zero values of keyword parameters record_ppathcloud and record_ppath\n"
         "enable the saving of internal and external ppath data for diagnostic purposes.\n"
         "  record_ppathcloud and record_ppath should be set to 0 for large values of\n"
         " max_iter.\n Negative values of rng_seed seed the random number generator \n "
         "according to system time, positive rng_seed values are taken literally.\n"
          "if keyword parameter silent is non-zero iterative output showing the photon\n"
         "number and scattering order are suppressed" ),
        OUTPUT(ppath_, ppath_step_, i_montecarlo_error_, rte_pos_, rte_los_,
               rte_gp_p_, rte_gp_lat_, rte_gp_lon_, i_space_, ground_emission_,
               ground_los_, ground_refl_coeffs_, i_rte_, 
               scat_za_grid_,scat_aa_grid_, rte_pressure_, rte_temperature_, 
               rte_vmr_list_, ext_mat_, abs_vec_, f_index_, scat_za_index_, 
               scat_aa_index_, ext_mat_spt_, abs_vec_spt_),
        INPUT(ppath_step_agenda_, atmosphere_dim_, p_grid_, lat_grid_,
              lon_grid_, z_field_, r_geoid_, z_ground_, cloudbox_limits_,
              stokes_dim_, rte_agenda_, i_space_agenda_, ground_refl_agenda_,
              t_field_, scat_za_grid_,
              scat_aa_grid_, f_grid_, opt_prop_gas_agenda_,
              spt_calc_agenda_,scalar_gas_absorption_agenda_, vmr_field_,
              scat_data_raw_, pnd_field_, scat_theta_, scat_theta_gps_,
              scat_theta_itws_,montecarlo_p_from_belowCsca_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("maxiter","rng_seed","record_ppathcloud","record_ppath","silent", 
                 "record_histdata", "histdata_filename", "los_sampling_method",
		 "strat_sampling"),
        TYPES( Index_t, Index_t, Index_t, Index_t, Index_t, Index_t, String_t, 
	       Index_t, Index_t)));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensorOff" ),
        DESCRIPTION
        (
         "Sets sensor WSVs to obtain monochromatic pencil beam values.\n"
         "\n"
         "The variables are set as follows:\n"
         "   sensor_response : As returned by *sensor_responseInit*.\n"
         "   sensor_pol      : Identity matrix, with size following\n"
         "                     *stokes_dim*.\n"
         "   sensor_rot      : Length matching *sensor_pos/los*. All values\n"
         "                     are set 0.\n"
         "   antenna_dim     : 1.\n"
         "   mblock_za_grid  : Length 1, value 0."
         "   mblock_aa_grid  : Empty."
        ),
        OUTPUT( sensor_response_, sensor_response_f_, sensor_response_za_,
                sensor_response_aa_, sensor_pol_, sensor_rot_,
                antenna_dim_, mblock_za_grid_, mblock_aa_grid_ ),
        INPUT( atmosphere_dim_, stokes_dim_, sensor_pos_, sensor_los_, 
               f_grid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_posAddGeoidWGS84" ),
        DESCRIPTION
        (
         "Adds a geoid radius according to WGS-84 to a geometric altitude.\n"
         "\n"
         "This function assumes that the first element of *sensor_pos* is\n"
         "set to the geometric altitude for the positions of the sensor. \n"
         "The variable *sensor_pos* shall contain the radius instead of the\n"
         "altitude and that can be achieved by this function. The function\n"
         "adds a geoid radius to the given altitude. The geoid radius is\n"
         "taken from the WGS-84 reference ellipsiod.\n"
         "\n"
         "For 1D, the geoid radius is set to the radius of curvature of the\n"
         "WGS-84 ellipsiod for the position and observation direction \n"
         "described with *lat_1d* and *meridian_angle_1d*.\n"
         "For 2D and 3D, the geoid radius is set to the radius of the WGS-84\n"
         "ellipsiod for the latitude values in *sensor_pos*."
        ),
        OUTPUT( sensor_pos_ ),
        INPUT( sensor_pos_, atmosphere_dim_, lat_1d_, meridian_angle_1d_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "sensor_posAddRgeoid" ),
        DESCRIPTION
        (
         "Adds a geoid radius by interpolating *r_geoid*.\n"
         "\n"
         "This function assumes that the first element of *rte_pos* is set\n"
         "to the geometric altitude for the position of the sensor. \n"
         "The variable *rte_pos* shall contain the radius instead of the\n"
         "altitude and that can be achieved by this function. The function\n"
         "adds a geoid radius to the given altitude. The geoid radius is\n"
         "obtained by interpolation of *r_geoid*. There is an error if the\n"
         "given position is outside the latitude and longitude grids."
        ),
        OUTPUT( sensor_pos_ ),
        INPUT( sensor_pos_, atmosphere_dim_, lat_grid_, lon_grid_, r_geoid_ ),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responseAntenna1D"),
        DESCRIPTION
        (
         "Returns the response block matrix after it has been modified by\n"
         "a 1D antenna response.\n"
         "\n"
         "The antenna diagram patterns are given as the generic input\n"
         "ArrayOfArrayOfMatrix. The structure of this variable is that\n"
         "the Matrix describes the antenna diagram values by a relative\n"
         "zenith angle grid, the ArrayOfMatrix then contains antenna\n"
         "diagrams for each polarisation given by the rows of *sensor_pol*\n"
         "and at the top level, the ArrayOfArrayOfMatrix contains antenna\n"
         "diagrams for each viewing angle of the antennas/beams described\n"
         "by the generic input Vector.\n"
         "\n"
         "The individual antenna diagrams, described by the matrices,\n"
         "contain at least two columns where the first column describes a\n"
         "relative grid of angles and the following column(s) describe\n"
         "the antenna diagram.\n"
         "\n"
         "For each level in the antenna diagram there exist two cases,\n"
         "either only one element/column of antenna gain values are given,\n"
         "this element/column will then be used for all directions/-\n"
         "polarisations/frequencies. Or else each direction/polarisation/-\n"
         "frequency is given its individual element/column.\n"
         "\n"
         "Generic Input: \n"
         "ArrayOfArrayOfMatrix : The antenna diagram(s).\n"
         "              Vector : The antenna/beam zenith angle grid."
        ),
        OUTPUT( sensor_response_, sensor_response_za_ ),
        INPUT( f_grid_, mblock_za_grid_, antenna_dim_, sensor_pol_ ),
        GOUTPUT( ),
        GINPUT( ArrayOfArrayOfMatrix_, Vector_ ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responseBackend"),
        DESCRIPTION
        (
         "Returns the response block matrix after it has been modified by\n"
         "a spectrometer backend response.\n"
         "\n"
         "The channel response is given as the generic input array of\n"
         "matrices, where each element in the array represent different\n"
         "polarisations given by *sensor_pol*. The individual matrices\n"
         "describe the channel responses as function of frequency, where the\n"
         "first column describes a relative grid of frequencies and the rest\n"
         "of the columns describe the backend response.\n"
         "\n"
         "For each level, the response can be described in two ways. Either\n"
         "one single array element/matrix column is given and then used for\n"
         "each polarisation/frequency. Or a complete set of array\n"
         "elements/matrix columns covering all polarisations/frequencies are\n"
         "given and in each case a individual response will be used.\n"
         "Note that for both cases there must allways be a column in the\n"
         "matrices, the first, of a relative frequency grid.\n"
         "\n"
         "Generic Input: \n"
         "  ArrayOfMatrix : The backend channel response."
        ),
        OUTPUT( sensor_response_, sensor_response_f_ ),
        INPUT( f_backend_, sensor_pol_, sensor_response_za_ ),
        GOUTPUT( ),
        GINPUT( ArrayOfMatrix_ ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responseInit"),
        DESCRIPTION
        (
         "Initialises the response block matrix to an identity matrix.\n"
         "\n"
         "The initialised matrix is a quadratic matrix with sidelength equal\n"
         "to the product of the length of *f_grid*, *mblock_za_grid*,\n"
         "*mblock_aa_grid* and the columns of *sensor_pol*."
        ),
        OUTPUT( sensor_response_, sensor_response_f_, sensor_response_za_,
                sensor_response_aa_  ),
        INPUT( f_grid_, mblock_za_grid_, mblock_aa_grid_, antenna_dim_,
               sensor_pol_, atmosphere_dim_, stokes_dim_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("sensor_responseMixer"),
        DESCRIPTION
        (
         "Returns the response block matrix after it has been modified by\n"
         "the mixer and sideband filter. The returned matrix converts RF to\n"
         "IF.\n"
         "\n"
         "The generic input matrix is a two-column matrix where the first\n"
         "column should be equal to *f_grid* and the second column desrcibes\n"
         "the sideband filter function.\n"
         "\n"
         "The local oscillator frequency is set by the keyword *lo*.\n"
         "\n"
         "Generic Input: \n"
         "       Matrix : The sideband filter response matrix."
        ),
        OUTPUT( sensor_response_, sensor_response_f_, f_mixer_ ),
        INPUT( sensor_pol_, sensor_response_za_, lo_ ),
        GOUTPUT( ),
        GINPUT( Matrix_ ),
        KEYWORDS( ),
        TYPES( )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("ScatteringDataPrepareDOIT"),
        DESCRIPTION
        (
         "Prepare single scattering data for a DOIT scattering calculation.\n"
         "\n"
         "This function has to be used for scattering calculations using the\n"
         "DOIT method. It prepares the data in such a way that the code can be\n"
         "speed-optimized.\n"
         "\n"
         "For different hydrometeor species different preparations are required.\n" 
         "So far, we only the case of randomly oriented particles is implemented.\n" 
         "\n"
         "For randomly oriented hydrometeor species the scattering angles in the\n"
         "scattering frame are precalculated for all possible incoming and\n" 
         "scattered directions and stored the the WSV *scat_theta*. In the program\n" 
         "all phase matrix emements have to be interpolated on the scattering angle.\n"
         "This has to be done for each iteration and for all frequencies.\n"
         "\n" 
         ),
        OUTPUT( scat_theta_, scat_theta_gps_,  scat_theta_itws_),
        INPUT( scat_za_grid_, scat_aa_grid_, scat_data_raw_ ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ScatteringDataPrepareDOITOpt"),
        DESCRIPTION
        (
          "Prepare single scattering data for a DOIT scattering calculation.\n"
         "\n"
         "This function has can be used for scatttering calcualtions using the\n" 
         " DOIT method. It is a method optimized for randomly oriented \n"
         "scattering media and it can only be used for such cases. \n"
         "\n"
         "It has to be used in *scat_mono_agenda*. The phase matrix data is \n"
         "interpolated on the actual frequency and on all scattering angles \n"
         "following from all possible combinations in *scat_za_grid* and \n"
         "*scat_aa_grid*. \n"
         "\n" 
          ),
        OUTPUT( scat_theta_, pha_mat_sptDOITOpt_),
        INPUT( za_grid_size_, scat_aa_grid_, scat_data_raw_, f_grid_, f_index_,
               atmosphere_dim_),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));


md_data_raw.push_back
    ( MdRecord
      ( NAME("ScatteringDataPrepareOFF"),
        DESCRIPTION
        (
        "No preparation of single scattering data.\n"
        "\n"
        "The parameters below are set to be empty. If this function is used\n" 
        "scattering angles, grid positions and interpolation weights are\n"
        "calculated inside the WSM *pha_mat_sptFromData*.\n"
        "The output variables are set to be empty. \n"
        "\n" 
         ),
        OUTPUT( scat_theta_, scat_theta_gps_,  scat_theta_itws_),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("StringSet"),
        DESCRIPTION
        (
         "Sets a String to the given text String."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( String_ ),
        GINPUT(),
        KEYWORDS( "text"   ),
        TYPES(    String_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor3FillWithVector"),
        DESCRIPTION
        (
         "Forms a tensor of order 3 by repeating a vector.\n"
         "\n"
         "The direction of the vector inside the tensor is selected by\n"
         "setting the size determined by the vector length to 0. For \n"
         "example, if the keyword *ncols* is set to 0, the vector will be\n"
         "put in as rows on every page. The remaining sizes are taken from \n"
         "the keyword arguments. \n"
         "\n"
         "One, but only one, keyword argument must be 0.\n"
         "\n"
         "Generic output: \n"
         "   Tensor3 : The tensor to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector to be copied. \n"
         "Keyword: \n"
         "   npages : Number of pages in the tensor.\n"
         "   nrows  : Number of rows in the tensor.\n"
         "   ncols  : Number of columns in the tensor. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor3_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "npages", "nrows", "ncols"   ),
        TYPES(    Index_t,  Index_t, Index_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor3Scale"),
        DESCRIPTION
        (
         "Scales a workspace tensor3 with the specified value. \n"
         "\n"
         "The result can either be stored in the input tensor3 or\n"
         "in a new tensor3.\n"
         "\n"
         "Generic output: \n"
         "   Tensor3 : The scaled tensor3. \n"
         "\n"
         "Generic input: \n"
         "   Tensor3 : The tensor3 to be scaled.\n"
         "\n"
         "Keywords:\n"
         "   value  : The scale factor. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor3_ ),
        GINPUT( Tensor3_ ),
        KEYWORDS( "value"   ),
        TYPES( Numeric_t )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor3Set"),
        DESCRIPTION
        (
         "Creates a workspace tensor3 with the specified size and sets \n"
         "all values of the tensor3 to the specified value. \n"
         "\n"
         "Generic output: \n"
         "   Tensor3 : The tensor3 to be created. \n"
         "\n"
         "Keywords:\n"
         "   npages : The number of pages of the tensor3 to create. \n"
         "   nrows  : The number of rows of the tensor3 to create. \n"
         "   ncols  : The number of columns of the tensor3 to create. \n"
         "   value  : The value of the tensor3 elements. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor3_ ),
        GINPUT(),
        KEYWORDS( "npages", "nrows", "ncols", "value"   ),
        TYPES( Index_t, Index_t, Index_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor4Scale"),
        DESCRIPTION
        (
         "Scales a workspace tensor4 with the specified value. \n"
         "\n"
         "The result can either be stored in the input tensor4 or\n"
         "in a new tensor4.\n"
         "\n"
         "Generic output: \n"
         "   Tensor4 : The scaled tensor4. \n"
         "\n"
         "Generic input: \n"
         "   Tensor4 : The tensor4 to be scaled.\n"
         "\n"
         "Keywords:\n"
         "   value  : The scale factor. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor4_ ),
        GINPUT( Tensor4_ ),
        KEYWORDS( "value"   ),
        TYPES( Numeric_t )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor4Set"),
        DESCRIPTION
        (
         "Creates a workspace tensor4 with the specified size and sets \n"
         "all values of the tensor4 to the specified value. \n"
         "\n"
         "Generic output: \n"
         "   Tensor4 : The tensor4 to be created. \n"
         "\n"
         "Keywords:\n"
         "   nbooks : The number of books of the tensor4 to create. \n"
         "   npages : The number of pages of the tensor4 to create. \n"
         "   nrows  : The number of rows of the tensor4 to create. \n"
         "   ncols  : The number of columns of the tensor4 to create. \n"
         "   value  : The value of the tensor4 elements. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor4_ ),
        GINPUT(),
        KEYWORDS( "nbooks", "npages", "nrows", "ncols", "value"   ),
        TYPES( Index_t, Index_t, Index_t, Index_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor5Scale"),
        DESCRIPTION
        (
         "Scales a workspace tensor5 with the specified value. \n"
         "\n"
         "The result can either be stored in the input tensor5 or\n"
         "in a new tensor5.\n"
         "\n"
         "Generic output: \n"
         "   Tensor5 : The scaled tensor5. \n"
         "\n"
         "Generic input: \n"
         "   Tensor5 : The tensor5 to be scaled.\n"
         "\n"
         "Keywords:\n"
         "   value  : The scale factor. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor5_ ),
        GINPUT( Tensor5_ ),
        KEYWORDS( "value"   ),
        TYPES( Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor5Set"),
        DESCRIPTION
        (
         "Creates a workspace tensor5 with the specified size and sets \n"
         "all values of the tensor5 to the specified value. \n"
         "\n"
         "Generic output: \n"
         "   Tensor5 : The tensor5 to be created. \n"
         "\n"
         "Keywords:\n"
         "   nshelfs : The number of shelfs of the tensor5 to create. \n"
         "   nbooks  : The number of books of the tensor5 to create. \n"
         "   npages  : The number of pages of the tensor5 to create. \n"
         "   nrows   : The number of rows of the tensor5 to create. \n"
         "   ncols   : The number of columns of the tensor5 to create. \n"
         "   value   : The value of the tensor5 elements. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor5_ ),
        GINPUT(),
        KEYWORDS( "nshelfs", "nbooks", "npages", "nrows", "ncols", "value" ),
        TYPES( Index_t, Index_t, Index_t, Index_t, Index_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor6Scale"),
        DESCRIPTION
        (
         "Scales a workspace tensor6 with the specified value. \n"
         "\n"
         "The result can either be stored in the input tensor6 or\n"
         "in a new tensor6.\n"
         "\n"
         "Generic output: \n"
         "   Tensor6 : The scaled tensor6. \n"
         "\n"
         "Generic input: \n"
         "   Tensor6 : The tensor6 to be scaled.\n"
         "\n"
         "Keywords:\n"
         "   value  : The scale factor. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor6_ ),
        GINPUT( Tensor6_ ),
        KEYWORDS( "value"   ),
        TYPES( Numeric_t )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor6Set"),
        DESCRIPTION
        (
         "Creates a workspace tensor6 with the specified size and sets \n"
         "all values of the tensor6 to the specified value. \n"
         "\n"
         "Generic output: \n"
         "   Tensor6 : The tensor6 to be created. \n"
         "\n"
         "Keywords:\n"
         "   nvitrines : The number of vitrines of the tensor6 to create. \n"
         "   nshelfs   : The number of shelfs of the tensor6 to create. \n"
         "   nbooks    : The number of books of the tensor6 to create. \n"
         "   npages    : The number of pages of the tensor6 to create. \n"
         "   nrows     : The number of rows of the tensor6 to create. \n"
         "   ncols     : The number of columns of the tensor6 to create. \n"
         "   value     : The value of the tensor6 elements. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor6_ ),
        GINPUT(),
        KEYWORDS( "nvitrines", "nshelfs", "nbooks", "npages", "nrows",
                  "ncols", "value" ),
        TYPES( Index_t, Index_t, Index_t, Index_t, Index_t, Index_t,
               Numeric_t )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "Tensor6ToTbByPlanck" ),
        DESCRIPTION
        (
         "Converts a Tensor6 of radiances to brightness temperatures by \n"
         "inverting the Planck function. \n"
         "\n"
         "Generic output: \n"
         "   Tensor6 : A Tensor6 with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Tenosr6 : A Tensor6 with radiance values. \n"
         ),
        OUTPUT(),
        INPUT(f_index_, f_grid_),
        GOUTPUT( Tensor6_ ),
        GINPUT( Tensor6_ ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor7Scale"),
        DESCRIPTION
        (
         "Scales a workspace tensor7 with the specified value. \n"
         "\n"
         "The result can either be stored in the input tensor7 or\n"
         "in a new tensor7.\n"
         "\n"
         "Generic output: \n"
         "   Tensor7 : The scaled tensor7. \n"
         "\n"
         "Generic input: \n"
         "   Tensor7 : The tensor7 to be scaled.\n"
         "\n"
         "Keywords:\n"
         "   value  : The scale factor. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor7_ ),
        GINPUT( Tensor7_ ),
        KEYWORDS( "value"   ),
        TYPES( Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor7Set"),
        DESCRIPTION
        (
         "Creates a workspace tensor7 with the specified size and sets \n"
         "all values of the tensor7 to the specified value. \n"
         "\n"
         "Generic output: \n"
         "   Tensor7 : The tensor7 to be created. \n"
         "\n"
         "Keywords:\n"
         "   nlibraries : The number of libraries of the tensor7 to create. \n"
         "   nvitrines  : The number of vitrines of the tensor7 to create. \n"
         "   nshelfs    : The number of shelfs of the tensor7 to create. \n"
         "   nbooks     : The number of books of the tensor7 to create. \n"
         "   npages     : The number of pages of the tensor7 to create. \n"
         "   nrows      : The number of rows of the tensor7 to create. \n"
         "   ncols      : The number of columns of the tensor7 to create. \n"
         "   value      : The value of the tensor7 elements. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Tensor7_ ),
        GINPUT(),
        KEYWORDS( "nlibraries", "nvitrines", "nshelfs", "nbooks", "npages",
                  "nrows", "ncols", "value" ),
        TYPES( Index_t, Index_t, Index_t, Index_t, Index_t, Index_t,
               Index_t, Numeric_t )));

 md_data_raw.push_back
    ( MdRecord
      ( NAME("Tensor6WriteIteration"),
        DESCRIPTION
        (
       "Write iterated fields.\n"
       "\n"
       "This function writes intermediate resultes, the iterations of \n"
       "fields to xml files. It can be used to check the solution method \n"
       "for the RTE with scattering integral, which is an iterative \n"
       "numerical method. It is useful to look how the radiation field \n"
       "*i_field* and the scattered field *scat_field* behave. \n"
       "\n"
       "The user can give an array containing the iterations which shall \n"
       "be written to files as a keyword to the method. E.g. if \n"
       "'iterations = [3, 6, 9]', the 3rd, 6th and 9th iterations are \n"
       "stored in the files 'iteration_field_3.xml', \n"
       "'iteration_field_6.xml' ...\n"
       "\n"
       "If you want to save all the iterations the array has to contain \n"
       "just one element set to 0: 'iterations = [0]'.\n"
       "\n"
       "Note: The workspace variable iteration_counter has to be set as 0 \n"
       "in the control file before using this method.\n"
       "\n"
        ),
        OUTPUT( ),
        INPUT(iteration_counter_),
        GOUTPUT( ),
        GINPUT(Tensor6_),
        KEYWORDS("iterations"),
        TYPES(Array_Index_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("Test"),
        DESCRIPTION
        (
         "A dummy method that can be used for test purposes.\n"
         "\n"
         "This method can be used by ARTS developers to quickly test stuff.\n"
         "The implementation is in file m_io.cc. This just saves you the \n"
         "trouble of adding a dummy method everytime you want to try \n"
         "something out quickly."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( ),
        KEYWORDS( ),
        TYPES( )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("timerStart"),
        DESCRIPTION
        (
         "Initializes the CPU timer."
         "\n"
         "Use *timerStop* to output the consumed cpu time\n"
         "since *timerStart*.\n"
         "\n"
         "Usage example:\n"
         "\n"
         "timerStart()"
         "ReadXML(f_grid){\"frequencies.xml\"}\n"
         "timerStop()\n"
         "Prints the CPU time spent for reading the XML file"
         ),
        OUTPUT(timer_),
        INPUT(),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));


  md_data_raw.push_back
    ( MdRecord
      ( NAME("timerStop"),
        DESCRIPTION
        (
         "Stops the CPU timer."
         "\n"
         "Use *timerStop* to output the consumed cpu time\n"
         "since *timerStart*. See *timerStart* for example"
         ),
        OUTPUT(),
        INPUT(timer_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS(),
        TYPES()));


  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorAddScalar"),
        DESCRIPTION
        (
         "Adds a scalar to all elements of a vector. \n"
         "\n"
         "The result can either be stored in the same or another vector. \n"
         "\n"
         "Generic output: \n"
         "   Vector : Return vector. \n"
         "\n"
         "Generic input: \n"
         "   Vector : Original vector. \n"
         "\n"
         "Keywords:\n"
         "   value : The value to be added to the vector."
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "value" ),
        TYPES( Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorLinSpace"),
        DESCRIPTION
        (
         "Creates a vector with linear spacing.\n"
         "\n"
         "The first element equals always the start value, and the spacing\n"
         "equlas always the step value, but note that the last value can  \n"
         "deviate from the stop value. The keyword step can be both positive\n"
         "and negative. \n"
         "   The vector is [start, start+step, start+2*step, ...]\n "
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Keywords:\n"
         "   start : The start value. \n"
         "    stop : The maximum value of the end value. \n"
         "    step : The spacing of the vector. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "start",   "stop",    "step"    ),
        TYPES(    Numeric_t, Numeric_t, Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorNLinSpace"),
        DESCRIPTION
        (
         "Creates a vector with defined length, equally spaced between the \n"
         "given end values. \n"
         "\n"
         "The length must be larger than 1. \n"
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Keywords:\n"
         "   start : The start value. \n"
         "    stop : The end value. \n"  
         "       n : Number of elements of the vector. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(Vector_),
        GINPUT(),
        KEYWORDS( "start",   "stop",    "n"   ),
        TYPES(    Numeric_t, Numeric_t, Index_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorNLogSpace"),
        DESCRIPTION
        (
         "Creates a vector with defined length, equally logarithmically \n"
         "spaced between the given end values. \n"
         "\n"
         "The length must be larger than 1. \n"
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Keywords:\n"
         "   start : The start value. \n"
         "    stop : The end value. \n"  
         "       n : Number of elements of the vector. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(Vector_),
        GINPUT(),
        KEYWORDS( "start",   "stop",    "n"   ),
        TYPES(    Numeric_t, Numeric_t, Index_t )));

  md_data_raw.push_back     
    ( MdRecord
      ( NAME("Print"),
        DESCRIPTION
        (
         "Prints a variable on the screen."
         "\n"
         "Keywords:\n"
         "   level : Output level to use. \n"
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( ),
        GINPUT( Any_ ),
        KEYWORDS( "level" ),
        TYPES( Index_t ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));


  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorScale"),
        DESCRIPTION
        (
         "Scales all elements of a vector with the same value. \n"
         "\n"
         "The result can either be stored in the same or another vector. \n"
         "\n"
         "Generic output: \n"
         "   Vector : Return vector. \n"
         "\n"
         "Generic input: \n"
         "   Vector : Original vector. \n"
         "\n"
         "Keywords:\n"
         "   value : The value to be multiplicated with the vector."  
        ),
        OUTPUT( ),
        INPUT( ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "value" ),
        TYPES( Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorSet"),
        DESCRIPTION
        (
         "Creates a workspace vector with the specified length and sets \n"
         "all values of the vector to the specified value. \n"
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Keywords:\n"
         "   length : The length of the new vector. \n"
         "   value  : The value of the vector elements. " 
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "length", "value"   ),
        TYPES(    Index_t,    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorSetExplicitly"),
        DESCRIPTION
        (
         "Create a vector from the given list of numbers.\n"
         "\n"
         "Generic output:\n"
         "   Vector : The vector to be created.\n"
         "\n"
         "Keywords:\n"
         "   values  : The vector elements.\n"
         "\n"
         "Usage:\n"
         "   VectorSetExplicitly(p_grid){[1000 100 10]}\n"
         "   Will create a p_grid vector with these three elements."
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "values"   ),
        TYPES(    Vector_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorSetTakingLengthFromVector"),
        DESCRIPTION
        (
         "Creates a workspace vector with the same length as another vector,\n"
         "and sets all values of the new vector to the specified value. \n"
         "\n"
         "A possible usage of the function is: \n"
         "  VectorSetLengthFromVector(vector1,f_grid){value=0.75} \n"
         "where *vector1* then can be used to set *e_ground*. \n"
         "\n"
         "Generic output: \n"
         "   Vector : The vector to be created. \n"
         "\n"
         "Generic input: \n"
         "   Vector : The vector specifying the length. \n"
         "\n"
         "Keywords:\n"
         "   value  : The value of the vector elements. "
        ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "value"   ),
        TYPES(    Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorToTbByPlanck" ),
        DESCRIPTION
        (
         "Converts a vector of radiances to brightness temperatures by \n"
         "inverting the Planck function.\n"
         "\n"
         "This function works as *VectorToTbByRJ*. However, this function \n"
         "is not recommended in connection with inversions, but can be used \n"
         "to display calculated spectra in a temperature scale.\n"
         "\n"
         "Generic output: \n"
         "   Vector : A vector with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A vector with radiance values. "
        ),
        OUTPUT(),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_,
               sensor_response_za_, sensor_response_aa_, sensor_pol_ ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS(),
        TYPES()));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorToTbByRJ" ),
        DESCRIPTION
        (
         "Converts a vector of radiances to brightness temperatures by \n"
         "the Rayleigh-Jeans approximation of the Planck function.\n"
         "\n"
         "This function performs a linear transformation of spectral \n"
         "radiances to an approximative temperature scale. The advantage \n"
         "of this linear transformation is that the obtained values can be \n"
         "used for retrievals if the weighting functions are handled \n"
         "likewise (by *MatrixToTbByRJ*). This is not the case if the \n"
         "radiances are converted to temparatures by the Planck function \n"
         "directly. \n"
         "\n"
         "The conversion assumes that the elements of the input vector are\n"
         "stored in standard order and correspond to *sensor_response_f*,\n"
         "*sensor_response_za*, *sensor_response_aa* and *sensor_pol*.\n"
         "The standard option shall accordingly be to perform this \n"
         "conversion directly after *RteCalc*.\n"
         "\n"
         "If *y* shall be converted from radiances to brightness \n"
         "temperatures: \n"
         "   VectorToTbByRJ(y,y){} \n"
         "\n"
         "Generic output: \n"
         "   Vector : A vector with brightness temperature values. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A vector with radiance values."
        ),
        OUTPUT(),
        INPUT( sensor_pos_, sensor_los_, sensor_response_f_, 
               sensor_response_za_, sensor_response_aa_, sensor_pol_ ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS(),
        TYPES()));

 md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorZtanToZaRefr" ),
        DESCRIPTION
        (
         "Converts a set of true tangent altitudes to zenith angles.\n"
         "\n"
         "The tangent altitudes are given to the function as a vector, which\n"
         "are converted to a generic vector of zenith angles. The position of\n"
         "the sensor is given by the WSV *sensor_pos*. The function works only\n"
         "for a spherical geoid. The zenith angles are always set to be\n"
                 "positive.\n"
         "The tangent altitudes are given as the altitude above the geoid.\n"
         "\n"
         "Generic output: \n"
         "   Vector : A vector with zenith angles. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A vector with true tangent altitudes\n"
        ),
        OUTPUT( refr_index_, rte_pressure_, rte_temperature_, rte_vmr_list_ ),
        INPUT( refr_index_agenda_, sensor_pos_, p_grid_, t_field_, z_field_,
                           vmr_field_, r_geoid_, atmosphere_dim_ ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS(),
        TYPES()));

md_data_raw.push_back
    ( MdRecord
      ( NAME( "VectorZtanToZa" ),
        DESCRIPTION
        (
         "Converts a set of geometrical tangent altitudes to zenith angles.\n"
         "\n"
         "The tangent altitudes are given to the function as a vector, which\n"
         "are converted to a generic vector of zenith angles. The position of\n"
         "the sensor is given by the WSV *sensor_pos*. The function works only\n"
         "for a spherical geoid, where the geoid radius is given as a keyword\n"
         "argument (*r_geoid*). The zenith angles are always set to be positive.\n"
         "The tangent altitudes are given as the altitude above the geoid.\n"
         "\n"
         "Generic output: \n"
         "   Vector : A vector with zenith angles. \n"
         "\n"
         "Generic input: \n"
         "   Vector : A vector with geometric tangent altitudes\n"
         "\n"
         "Keywords:\n"
         "   r_geoid : The geoid radius for the given tangent altitudes."
        ),
        OUTPUT(),
        INPUT( sensor_pos_ ),
        GOUTPUT( Vector_ ),
        GINPUT( Vector_ ),
        KEYWORDS( "r_geoid" ),
        TYPES( Numeric_t )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("WriteXML"),
        DESCRIPTION
        (
         "Writes a workspace variable to an XML file.\n"
         "\n"
         "This is a supergeneric method. It can write variables of any group.\n"
         "\n"
         "If the filename is omitted, the variable is written\n"
         "to <basename>.<variable_name>.xml.\n"
         "\n"
         "Usage example:\n"
         "\n"
         "WriteXML(f_grid){\"\"}\n"
         "Will write the frequency grid *f_grid* to the default file.\n"
         "\n"
         "Supergeneric input:\n"
         "   Any_     : The variable to write.\n"
         "\n"
         "Keywords:\n"
         "   filename : Name of the output file."
         ),
        OUTPUT(),
        INPUT( output_file_format_ ),
        GOUTPUT( ),
        GINPUT(  Any_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   ),
        AGENDAMETHOD(   false ),
        SUPPRESSHEADER( true  )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ybatchMetProfiles" ),
        DESCRIPTION
        (
         "This method is used for simulating ARTS for metoffice model fields"
         "\n"
         "This method reads in *met_amsu_data* which contains the\n"
         "lat-lon of the metoffice profile files as a Matrix. It then \n"
         "loops over the number of profiles and corresponding to each \n"
         "longitude create the appropriate profile basename. Then, \n"
         "Corresponding to each basename we have temperature field, altitude\n"
         "field, humidity field and particle number density field.  The\n"
         "temperature field and altitude field are stored in the same dimensions\n"
         "as *t_field_raw* and *z_field_raw*.  The oxygen and nitrogen VMRs are\n"
         "set to constant values of 0.209 and 0.782, respectively and are used\n"
         "along with humidity field to generate *vmr_field_raw*.  \n"
         "\n"
         "The three fields *t_field_raw*, *z_field_raw*, and *vmr_field_raw* are\n"
         "given as input to *met_profile_calc_agenda* which is called in this\n"
         "method.  See documentation of WSM *met_profile_calc_agenda* for more\n"
         "information on this agenda.  \n"
         "\n"
         "The method also converts satellite zenith angle to appropriate \n"
         "*sensor_los*.  It also sets the *p_grid* and *cloudbox_limits* \n"
         "from the profiles inside the function\n"
         ),
        OUTPUT( ybatch_, y_, t_field_raw_, z_field_raw_, vmr_field_raw_, 
		pnd_field_raw_,	 p_grid_, sensor_los_,cloudbox_on_, 
		cloudbox_limits_, z_ground_),
        INPUT(gas_species_, met_profile_calc_agenda_, f_grid_, met_amsu_data_,
	      sensor_pos_, r_geoid_, lat_grid_, lon_grid_, atmosphere_dim_,
	      scat_data_raw_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("nelem_p_grid", "met_profile_path", "met_profile_pnd_path"),
        TYPES(Index_t, String_t, String_t)));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ybatchMetProfilesClear" ),
        DESCRIPTION
        (
         "This method is used for simulating ARTS for metoffice model fields\n"
         "for clear sky conditions.\n"
         "\n"
         "This method reads in *met_amsu_data* which contains the\n"
         "lat-lon of the metoffice profile files as a Matrix. It then \n"
         "loops over the number of profiles and corresponding to each \n"
         "longitude create the appropriate profile basename. Then, \n"
         "Corresponding to each basename we have temperature field, altitude\n"
         "field, humidity field and particle number density field.  The\n"
         "temperature field and altitude field are stored in the same dimensions\n"
         "as *t_field_raw* and *z_field_raw*.  The oxygen and nitrogen VMRs are\n"
         "set to constant values of 0.209 and 0.782, respectively and are used\n"
         "along with humidity field to generate *vmr_field_raw*.  \n"
         "\n"
         "The three fields *t_field_raw*, *z_field_raw*, and *vmr_field_raw* are\n"
         "given as input to *met_profile_calc_agenda* which is called in this\n"
         "method.  See documentation of WSM *met_profile_calc_agenda* for more\n"
         "information on this agenda.  \n"
         "\n"
         "The method also converts satellite zenith angle to appropriate \n"
         "*sensor_los*.  It also sets the *p_grid* and *cloudbox_limits* \n"
         "from the profiles inside the function\n"
         ),
        OUTPUT( ybatch_, t_field_raw_, z_field_raw_, vmr_field_raw_, 
                 y_, p_grid_, sensor_los_, z_ground_),
        INPUT(gas_species_, met_profile_calc_agenda_, 
              f_grid_, met_amsu_data_, sensor_pos_, r_geoid_),
        GOUTPUT(),
        GINPUT(),
        KEYWORDS("nelem_p_grid","met_profile_path" ),
        TYPES(Index_t, String_t)));

  md_data_raw.push_back
    ( MdRecord
      ( NAME( "ZaSatOccultation" ),
        DESCRIPTION
        (
         "Calculates zenith angles for satellite occultations.\n"
         "\n"
         "The zenith angles are calculated with an interval of *t_sample\n"
         "with the recieving satellite at height *z_recieve* above the geoid\n"
         "and the transmitting satellite at height *z_send*.\n"
         "The zenith angles are restricted by the two tangent altitudes\n"
         "*z_scan_low* and *z_scan_high*."
        ),
        OUTPUT( ppath_step_),
        INPUT( ppath_step_agenda_, atmosphere_dim_, p_grid_, lat_grid_,
               lon_grid_, z_field_, r_geoid_, z_ground_ ),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "z_recieve", "z_send", "t_sample", 
                  "z_scan_low", "z_scan_high" ),
        TYPES( Numeric_t, Numeric_t, Numeric_t,
               Numeric_t, Numeric_t )));








  //
  // Below this line you find methods not touched for ARTS 2. 
  // Please revise the documentation etc. before a methods is moved up.
  // Place functions in alphabetical order. 
  // Finally, all methods below the line will be removed.
  //--------------------------------------------------------------------

//======================================================================
//=== IO methods
//======================================================================

//=== Index ============================================================

  md_data_raw.push_back
    ( MdRecord
      ( NAME("IndexWriteAscii"),
        DESCRIPTION(
                    "Writes an index value to an ASCII file.\n"
                    "\n"
                    "The index value of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the index is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Index : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Index_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("IndexReadAscii"),
        DESCRIPTION(
                    "Reads a index value from an ASCII file.\n"
                    "\n"
                    "The index value is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the index is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Index : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Index_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("IndexWriteBinary"),
        DESCRIPTION(
                    "Writes an index to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Index : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename: Name of the output file.\n"
                    ),
   
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Index_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("IndexReadBinary"),
        DESCRIPTION(
                    "Reads an index from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Index : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Index_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT

//=== NUMERIC ==========================================================

  md_data_raw.push_back
    ( MdRecord
      ( NAME("NumericWriteAscii"),
        DESCRIPTION(
                    "Writes a numeric value to an ASCII file.\n"
                    "\n"
                    "The numeric value of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the numeric is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Numeric : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Numeric_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("NumericReadAscii"),
        DESCRIPTION(
                    "Reads a numeric value from an ASCII file.\n"
                    "\n"
                    "The numeric value is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the numeric is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Numeric : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Numeric_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("NumericWriteBinary"),
        DESCRIPTION(
                    "Writes a numeric value to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Numeric : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Numeric_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("NumericReadBinary"),
        DESCRIPTION(
                    "Reads a numeric from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Numeric : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Numeric_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT


//=== Vector ==========================================================



  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorWriteAscii"),
        DESCRIPTION(
                    "Writes a vector to an ASCII file.\n"
                    "\n"
                    "The vector of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the vector is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Vector : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Vector_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorReadAscii"),
        DESCRIPTION(
                    "Reads a vector from an ASCII file.\n"
                    "\n"
                    "The vector is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the vector is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Vector : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorWriteBinary"),
        DESCRIPTION(
                    "Writes a vector to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Vector : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Vector_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("VectorReadBinary"),
        DESCRIPTION(
                    "Reads a vector from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n" 
                    "Generic output: \n"
                    "   Vector : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Vector_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT




//=== Matrix ==========================================================


  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixWriteAscii"),
        DESCRIPTION(
                    "Writes a matrix to an ASCII file.\n"
                    "\n"
                    "The matrix of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the matrix is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Matrix : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Matrix_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixReadAscii"),
        DESCRIPTION(
                    "Reads a matrix from an ASCII file.\n"
                    "\n"
                    "The matrix is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the matrix is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Matrix : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixWriteBinary"),
        DESCRIPTION(
                    "Writes a matrix to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   Matrix : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( Matrix_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("MatrixReadBinary"),
        DESCRIPTION(
                    "Reads a matrix from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   Matrix : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( Matrix_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT





//=== ArrayOfIndex =====================================================

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexWriteAscii"),
        DESCRIPTION(
                    "Writes a index array to an ASCII file.\n"
                    "\n"
                    "The index array of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the index array is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfIndex : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfIndex_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexReadAscii"),
        DESCRIPTION(
                    "Reads a index array from an ASCII file.\n"
                    "\n"
                    "The index array is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the index array is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfIndex : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"                    
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfIndex_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexWriteBinary"),
        DESCRIPTION(
                    "Writes an index array to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfIndex : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfIndex_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfIndexReadBinary"),
        DESCRIPTION(
                    "Reads an index array from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfIndex : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfIndex_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT


//=== ArrayOfVector ====================================================

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfVectorWriteAscii"),
        DESCRIPTION(
                    "Writes an array of vectors to an ASCII file.\n"
                    "\n"
                    "The array of vectors of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the array of vectors is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfVector : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfVector_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfVectorReadAscii"),
        DESCRIPTION(
                    "Reads an array of vectors from an ASCII file.\n"
                    "\n"
                    "The array of vectors is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the array of vectors is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfVector : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfVector_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfVectorWriteBinary"),
        DESCRIPTION(
                    "Writes a vector array to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfVector : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfVector_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfVectorReadBinary"),
        DESCRIPTION(
                    "Reads a vector array from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output:  \n"
                    "   ArrayOfVector : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfVector_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT



//=== ArrayOfMatrix ====================================================

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixWriteAscii"),
        DESCRIPTION(
                    "Writes an array of matrices to an ASCII file.\n"
                    "\n"
                    "The array of matrices of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the array of matrices is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "The format is as follows:\n"
                    "\n"
                    "# <comments>\n"
                    "<n_array_elements>\n"
                    "<n_rows> <n_columns>\n"
                    "<elements>\n"
                    "<n_rows> <n_columns>\n"
                    "<elements>\n"
                    "...\n"
                    "\n"
                    "Example:\n"
                    "# Generated by arts-0.0.16, Apr 29 2000, 17:38:44\n"
                    "2\n"
                    "3 4\n"
                    "xx xx xx xx\n"
                    "xx xx xx xx\n"
                    "xx xx xx xx\n"
                    "2 2\n"
                    "yy yy\n"
                    "yy yy"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfMatrix : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfMatrix_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixReadAscii"),
        DESCRIPTION(
                    "Reads an array of matrices from an ASCII file.\n"
                    "\n"
                    "The array of matrices is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the array of matrices is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfMatrixWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfMatrix : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfMatrix_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixWriteBinary"),
        DESCRIPTION(
                    "Writes a matrix array to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfMatrix : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfMatrix_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfMatrixReadBinary"),
        DESCRIPTION(
                    "Reads a matrix array from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfMatrix : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfMatrix_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT

//=== STRING ============================================================


  md_data_raw.push_back
    ( MdRecord
      ( NAME("StringWriteAscii"),
        DESCRIPTION(
                    "Writes a string to an ASCII file.\n"
                    "\n"
                    "The string of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the string is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfStringWriteAscii* for file format.\n"
                    "\n"
                    "Generic input: \n"
                    "   String : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( String_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("StringReadAscii"),
        DESCRIPTION(
                    "Reads a string from an ASCII file.\n"
                    "\n"
                    "The string is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the string is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfStringWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   String : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( String_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("StringWriteBinary"),
        DESCRIPTION(
                    "Writes a String to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   String : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( String_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("StringReadBinary"),
        DESCRIPTION(
                    "Reads a String from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   String : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( String_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT


//=== ArrayOfSTRING =========================================================


  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfStringWriteAscii"),
        DESCRIPTION(
                    "Writes an array of strings to an ASCII file.\n"
                    "\n"
                    "The array of strings of the given workspace variable\n"
                    "is written to the file with the specified name.\n"
                    "If the filename is omitted, the string is written\n"
                    "to <basename>.<variable_name>.aa.\n"
                    "\n"
                    "The format is as follows:\n"
                    "\n"
                    "# <comments>\n"
                    "<n_Strings>\n"
                    "<String 1>\n"
                    "<String 2>\n"
                    "...\n"
                    "\n"
                    "Example:\n"
                    "# Generated by arts-0.0.16, Apr 29 2000, 17:38:44\n"
                    "2\n"
                    "A String\n"
                    "Another String\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfString : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfString_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfStringReadAscii"),
        DESCRIPTION(
                    "Reads an array of strings from an ASCII file.\n"
                    "\n"
                    "The array of strings is read from the file with the\n"
                    "specified name and stored in the given workspace\n"
                    "variable.\n"
                    "If the filename is omitted, the string is read\n"
                    "from <basename>.<variable_name>.aa.\n"
                    "\n"
                    "See *ArrayOfStringWriteAscii* for file format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfString : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords: \n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfString_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

#ifdef HDF_SUPPORT
  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfStringWriteBinary"),
        DESCRIPTION(
                    "Writes a String array to a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic input: \n"
                    "   ArrayOfString : Name of the workspace variable to write.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the output file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT(),
        GINPUT( ArrayOfString_ ),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));

  md_data_raw.push_back
    ( MdRecord
      ( NAME("ArrayOfStringReadBinary"),
        DESCRIPTION(
                    "Reads a String array from a binary file.\n"
                    "\n"
                    "The filename can be specified or an empty String.\n"
                    "If empty, it is set to <basename>.<variable_name>.ab.\n"
                    "File is in HDF4 format.\n"
                    "\n"
                    "Generic output: \n"
                    "   ArrayOfString : Name of the workspace variable to read.\n"
                    "\n"
                    "Keywords:\n"
                    "   filename : Name of the input file.\n"
                    ),
        OUTPUT(),
        INPUT(),
        GOUTPUT( ArrayOfString_ ),
        GINPUT(),
        KEYWORDS( "filename" ),
        TYPES(    String_t   )));
#endif // HDF_SUPPORT



//======================================================================
//=== Absorption methods
//======================================================================

//=== Spectroscopic methods ============================================

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgSetEmpty"),
//      DESCRIPTION
//      (
//       "Sets lines_per_tg to empty line lists.\n"
//       "\n"
//       "You can use this method to set lines per tag if you do not reall want\n"
//       "to compute line spectra. Formally, absCalc will still require\n"
//       "lines_per_tg to be set.\n"
//       ),
//      OUTPUT(   lines_per_tg_      ),
//      INPUT(    tgs_        ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(  ),
//      TYPES(    )));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgReadFromCatalogues"),
//      DESCRIPTION(
//                  "This method can read lines from different line \n"
//                  "catalogues.\n"
//                  "\n"
//                  "For each tag group, you can specify which catalogue\n"
//                  "to use. Because the method creates lines_per_tg directly,\n"
//                  "it replaces for example thefollowing two method calls:\n"
//                  "  - linesReadFromHitran\n"
//                  "  - lines_per_tgCreateFromLines\n"
//                  "   This method needs as input WSVs the list of tag \n"
//                  "groups. Keyword parameters must specify the names of\n"
//                  "the catalogue files to use and the matching formats.\n"
//                  "Names can be anything, formats can currently be \n"
//                  "HITRAN96, MYTRAN2, JPL, or ARTS. Furthermore, keyword\n"
//                  "parameters have to specify minimum and maximum \n"
//                  "frequency for each tag group. To safe typing, if there\n"
//                  "are less elements in the keyword parameters than there\n"
//                  "are tag groups, the last parameters are applied to all\n"
//                  "following tag groups.\n"
//                  "\n"
//                  "Example usage:\n"
//                  "\n"
//                  "lines_per_tgReadFromCatalogues{\n"
//                  "  filenames = [ \"../data/cat1.dat\", \"../data/cat2.dat\" ]\n"
//                  "  formats   = [ \"MYTRAN2\",          \"HITRAN96\"         ]\n"
//                  "  fmin      = [ 0,                  0                  ]\n"
//                  "  fmax      = [ 2000e9,             100e9              ]\n"
//                  "}\n"
//                  "   In this example, lines for the first tag group will\n"
//                  "be taken from cat1, lines for all other tag groups \n"
//                  "will be taken from cat2.\n"
//                  "   This methods allows you for example to use a \n"
//                  "special line file just for water vapor lines. This\n"
//                  "could be the  improved water vapor line file \n"
//                  "generated by Thomas Kuhn.\n"
//                  "   Catalogues are only read once, even if several tag\n"
//                  "groups have the same catalogue. However, in that case\n"
//                  "the frequency ranges MUST be the same. (If you want \n"
//                  "to do fine-tuning of the frequency ranges, you can do \n"
//                  "this inside the tag definitions, e.g., \"H2O-*-0-2000e9\".)\n"
//                  "   This function uses the various reading routines\n"
//                  "(linesReadFromHitran, etc.), as well as\n"
//                  "lines_per_tgCreateFromLines.\n"
//                  "\n"
//                  "Keywords: \n"
//                  "   filenames = Name (and path) of the catalogue files.\n"
//                  "   formats   = allowed formats are HITRAN96,MYTRAN2,JPL,ARTS \n"
//                  "   fmin      = Minimum frequency for lines to read in Hz.\n"
//                  "   fmax      = Maximum frequency for lines to read in Hz.\n"),
//      OUTPUT(   lines_per_tg_      ),
//      INPUT(    tgs_        ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filenames",    "formats",      "fmin",   "fmax" ),
//      TYPES(    Array_String_t, Array_String_t, Vector_t, Vector_t)));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesReadFromHitran"),
//      DESCRIPTION(
//                  "Read all the lines from a HITRAN catalogue file in the \n"
//                  "given frequency range. Otherwise a runtime error will be\n"
//                  "thrown\n"
//                  "\n"
//                  "Please note that all lines must correspond\n"
//                  "to the legal species / isotope combinations\n"
//                  "\n"
//                  "Keywords: \n"
//                  "   filename = Name (and path) of the catalogue file.\n"
//                  "   fmin     = Minimum frequency for lines to read in Hz.\n"
//                  "   fmax     = Maximum frequency for lines to read in Hz."),
//      OUTPUT(   lines_   ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filename",  "fmin",    "fmax"),
//      TYPES(    String_t,    Numeric_t, Numeric_t)));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesReadFromMytran2"),
//      DESCRIPTION(
//                  "Read all the lines from a MYTRAN2 catalogue file in the \n"
//                  "given frequency range. Otherwise a runtime error will be\n"
//                  "thrown\n"
//                  "\n"
//                  "Please note that all lines must correspond\n"
//                  "to the legal species / isotope combinations\n"
//                  "\n"
//                  "Keywords: \n"
//                  "   filename = Name (and path) of the catalogue file.\n"
//                  "   fmin     = Minimum frequency for lines to read in Hz.\n"
//                  "   fmax     = Maximum frequency for lines to read in Hz."),
//      OUTPUT(   lines_   ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filename",  "fmin",    "fmax"),
//      TYPES(    String_t,    Numeric_t, Numeric_t)));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesReadFromJpl"),
//      DESCRIPTION(
//                  "Read all the lines from a JPL catalogue file in the \n"
//                  "given frequency range. Otherwise a runtime error will be\n"
//                  "thrown\n"
//                  "\n"
//                  "Please note that all lines must correspond\n"
//                  "to the legal species / isotope combinations.\n"
//                  "\n"
//                  "Keywords: \n"
//                  "   filename = Name (and path) of the catalogue file.\n"
//                  "   fmin     = Minimum frequency for lines to read in Hz.\n"
//                  "   fmax     = Maximum frequency for lines to read in Hz."),
//      OUTPUT(   lines_   ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filename",  "fmin",    "fmax"),
//      TYPES(    String_t,    Numeric_t, Numeric_t)));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesReadFromArts"),
//      DESCRIPTION(
//                  "Read all the lines from an Arts catalogue file in the \n"
//                  "given frequency range. Otherwise a runtime error will be\n"
//                  "thrown \n"
//                  "\n"
//                  "Please note that all lines must correspond\n"
//                  "to the legal species / isotope combinations\n"
//                  "\n"
//                  "Keywords: \n"
//                  "   filename = Name (and path) of the catalogue file.\n"
//                  "   fmin     = Minimum frequency for lines to read in Hz.\n"
//                  "   fmax     = Maximum frequency for lines to read in Hz."),
//      OUTPUT(   lines_   ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filename",  "fmin",    "fmax"),
//      TYPES(    String_t,    Numeric_t, Numeric_t)));
  
//   // FIXME: Remove this one.
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesElowToJoule"),
//      DESCRIPTION(
//                  "Just a little helper to convert the lower state energy from cm^-1\n"
//                  "(ARTSCAT-2) to Joule (ARTSCAT-3). This should be removed soon\n"),
//      OUTPUT(   lines_   ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( ),
//      TYPES(    )));
      
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgCreateFromLines"),
//      DESCRIPTION(
//                  "Split lines up into the different tag groups.\n"
//                  "\n"
//                  "The tag groups are tested in the order in which they are\n" 
//                  "specified in the controlfile. The lines are assigned to \n"
//                  "the tag groups in the order as the groups  are specified.\n"
//                  "That means if you do [\"O3-666\",\"O3\"],the last group O3 \n"
//                  "gets assigned all the O3 lines that do not fit in the first group."),
//      OUTPUT(   lines_per_tg_      ),
//      INPUT(    lines_, tgs_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgAddMirrorLines"),
//      DESCRIPTION(
//                  "Adds mirror lines at negative frequencies to the *lines_per_tg*.\n"
//                  "\n"
//                  "For each line at frequency +f in *lines_per_tg* a corresponding\n"
//                  "entry at frequency -f is added to *lines_per_tg*.The mirror \n"
//                  "lines are appended to the line lists after the original lines."),
//      OUTPUT(   lines_per_tg_      ),
//      INPUT(    lines_per_tg_      ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgCompact"),
//      DESCRIPTION(
//                  "Removes all lines outside the defined lineshape cutoff frequency\n"
//                  "from the *lines_per_tg*. This can save computation time.\n"
//                  "It should be particularly useful to call this method after\n"
//                  "*lines_per_tgAddMirrorLines*."),
//      OUTPUT(   lines_per_tg_      ),
//      INPUT(    lines_per_tg_, lineshape_, f_mono_  ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("linesWriteAscii"),
//      DESCRIPTION(
//                     "Writes the workspace variable *lines* to an ASCII file.\n"
//                     "\n"
//                  "The content of the workspace variable 'lines`\n"
//                  "The content of the workspace variable *lines*\n"
//                  "is written in ARTS line format to the file with\n"
//                     "the specified name. If the filename is omitted, the\n"
//                     "lines are written to <basename>.lines.aa.\n"
//                     "\n"
//                     "Keywords: \n"
//                     "   filename : Name of the output file.\n"
//                     ), 
//      OUTPUT(),
//      INPUT( lines_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filename" ),
//      TYPES(    String_t   )));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lines_per_tgWriteAscii"),
//      DESCRIPTION(
//                     "Writes the workspace variable *lines_per_tg* to an ASCII file.\n"
//                     "\n"
//                     "The content of the workspace variable *lines_per_tg*\n"
//                     "is written in ARTS line format to the file with\n"
//                     "the specified name. If the filename is omitted, the\n"
//                     "lines are written to <basename>.lines_per_tg.aa.\n"
//                     "\n"
//                     "The array dimension is handled in a similar way as by the\n"
//                     "array of vector and matrix output functions:\n"
//                     "First an integer stating the number of tag groups.\n"
//                     "Then an integer specifying the number of lines for the\n"
//                     "first group. Then the other groups in similar fashion."
//                     "\n"
//                     "Keywords: \n"
//                     "   filename : Name of the output file.\n"
//                     ),
//      OUTPUT(),
//      INPUT( lines_per_tg_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "filename" ),
//      TYPES(    String_t   )));
  
//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("tgsDefineAllInScenario"),
//      DESCRIPTION
//      (
//       "Define one tag group for each species known to ARTS and included in an\n"
//       "atmospheric scenario.\n"
//       "\n"
//       "You can use this as an alternative to tgsDefine if you want to make an\n"
//       "absorption calculation that is as complete as possible. The method\n"
//       "goes through all defined species and tries to open the VMR file. If\n"
//       "this works the tag is included, otherwise it is skipped.\n"
//       "\n"
//       "Keywords:\n"
//       "   basename : The name and path of a particular atmospheric scenario.\n"
//       "              For example: /pool/lookup2/arts-data/atmosphere/fascod/tropical"
//            ),
//      OUTPUT( tgs_ ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "basename" ),
//      TYPES(    String_t   )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lineshapeDefine"),
//      DESCRIPTION(
//           "Sets the lineshape for all calculated lines.\n\n"
//           "\n"
//           "   A general lineshape profile is specified, according to a given  \n"
//           "approximation. Alongside a normalization factor is to be set - a  \n"
//           "multiplicative forefactor through which the profile can be \n"
//           "modified. This factor is just the 0th or 1st, or 2nd power of the \n"
//           "ratio between the frequency of calculation f and the center frequency\n"
//           "for a specific line f0. A cutoff frequency must also be specified in\n"
//           "order to restrict the calculation within a desired frequency region or\n"
//           "not, when there's no such region.\n"
//           "   The general lineshape profile is given by the keyword shape,\n"
//           "while the normalization factor and the cutoff frequency by\n"
//           "normalizationfactor and cutoff respectively.\n"
//           "\n"
//           "   The available values for these keywords are given below.\n"
//           "shape - \"no_shape\" : no specified shape\n"
//           "        \"Doppler\" : Doppler lineshape\n"
//           "        \"Lorentz\" : Lorentz lineshape\n"
//           "        \"Voigt_Kuntz3\" : Kuntz approximation to the Voigt profile,\n"
//           "                         accuracy > 2x10^(-3)\n"
//           "        \"Voigt_Kuntz4\" : Kuntz approximation to the Voigt profile,\n"
//           "                         accuracy > 2x10^(-4)\n"
//           "        \"Voigt_Kuntz6\" : Kuntz approximation to the Voigt profile,\n"
//           "                         accuracy > 2x10^(-6)\n"   
//           "        \"Voigt_Drayson\" : Drayson approximation to the Voigt profile \n"
//           "        \"Rosenkranz_Voigt_Drayson\" : Rosenkrantz oxygen absortion with overlap correction\n" 
//           "                                     on the basis of Drayson routine\n"                                    
//           "        \"Rosenkranz_Voigt_Kuntz6\" : Rosenkrantz oxygen absortion with overlap correction\n"
//           "                                    on the basis of Kuntz routine, accuracy > 2x10^(-6)\n"
//           "normalizationfactor - \"no_norm\": 1\n"
//           "                      \"linear\": f/f0\n" 
//           "                      \"quadratic\": (f/f0)^2.\n"
//           "cutoff - \" -1\" : no cutoff\n"
//           "         \"Number\": positive cutoff frequency in Hz.\n"
//           "\n"
//           "Example usage:\n"
//        "shape=[\"Lorentz\"]\n"
//           "normalizationfactor=[\"linear\"]\n"
//           "cutoff= [650e9]"
//           "\n"
//           "Keywords:\n"
//           "   shape               : The general profile according to an approximation.\n"
//           "   normalizationfactor : The multiplicative forefactor for the general profile.\n"
//           "   cutoff              : The frequency at which a cutoff can be made.\n"),
//      OUTPUT( lineshape_ ),
//      INPUT( tgs_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(  "shape",    "normalizationfactor",  "cutoff" ),
//      TYPES(     String_t,        String_t,         Numeric_t )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("lineshape_per_tgDefine"),
//      DESCRIPTION(
//           "Sets the lineshape per tag group for all calculated lines.\n\n"
//           "\n" 
//           "   A general lineshape profile is specified, according to a given  \n"
//           "approximation for each tag group. Alongside a normalization factor\n" 
//           "is to be set also for each tag group - a multiplicative forefactor through\n"
//           "which the profile can be modified. This factor is just the 0th or 1st,\n"
//           "or 2nd power of the ratio between the frequency of calculation f and\n"
//           "the center frequency for a specific line f0. A cutoff frequency must also be\n"
//           "specified for each of the tags in  order to restrict the calculation within\n" 
//           "a desired region or not, when there's no such region.\n"
//           "   The general lineshape profile is given by the keyword shape,\n"
//           "while the normalization factor and the cutoff frequency by\n"
//           "normalizationfactor and cutoff respectively.\n"
//           "\n"
//           "   The available values for these keywords are given below.\n"
//           "shape - \"no_shape\" : no specified shape\n"
//           "        \"Doppler\" : Doppler lineshape\n"
//           "        \"Lorentz\" : Lorentz lineshape\n"
//           "        \"Voigt_Kuntz3\" : Kuntz approximation to the Voigt profile,\n"
//           "                        accuracy > 2x10^(-3)\n"
//           "        \"Voigt_Kuntz4\" : Kuntz approximation to the Voigt profile,\n"
//           "                         accuracy > 2x10^(-4)\n"
//           "        \"Voigt_Kuntz6\" : Kuntz approximation to the Voigt profile,\n"
//           "                         accuracy > 2x10^(-6)\n"   
//           "        \"Voigt_Drayson\" : Drayson approximation to the Voigt profile \n"
//           "        \"Rosenkranz_Voigt_Drayson\" : Rosenkrantz oxygen absortion with overlap correction\n" 
//           "                                     on the basis of Drayson routine\n"                                    
//           "        \"Rosenkranz_Voigt_Kuntz6\" : Rosenkrantz oxygen absortion with overlap correction\n"
//           "                                    on the basis of Kuntz routine, accuracy > 2x10^(-6)\n"
//           "normalizationfactor - \"no_norm\": 1\n"
//           "                      \"linear\": f/f0\n" 
//           "                      \"quadratic\": (f/f0)^2.\n"
//           "cutoff - \" -1\" : no cutoff\n"
//           "           \"Number\": positive cutoff frequency in Hz.\n"
//           "\n"
//           "Example usage:\n"
//        "shape = [\"Lorentz\",\"Voigt_Kuntz6\"] \n"
//        "normalizationfactor= [\"linear\", \"quadratic\"] \n"
//        "cutoff = [ 650e9, -1 ]"
//           "\n"
//           "Keywords:\n"
//           "   shape               : The general profile according to an approximation.\n"
//           "   normalizationfactor : The multiplicative forefactor for the general profile.\n"
//           "   cutoff              : The frequency at which a cutoff can be made.\n"),
//      OUTPUT( lineshape_ ),
//      INPUT( tgs_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(  "shape",           "normalizationfactor",    "cutoff" ),
//      TYPES(   Array_String_t,         Array_String_t,        Vector_t )));


// //=== Continuum methods ============================================

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("cont_descriptionInit"),
//      DESCRIPTION
//      (
//       "Initializes the two workspace variables for the continuum description,\n"
//       "*cont_description_names* and *cont_description_parameters*.\n"
//       " \n"
//       "This method does not really do anything, except setting the two\n"
//       "variables to empty Arrays. It is just necessary because the method\n"
//       "*cont_descriptionAppend* wants to append to the variables.\n"
//       "   Formally, the continuum description workspace variables are required\n"
//       "by the absorption calculation methods (e.g., *absCalc*). Therefore you\n"
//       "always have to call at least *cont_descriptionInit*, even if you do\n"
//       "not want to use any continua."
//       ),
//      OUTPUT( cont_description_names_, 
//                 cont_description_models_,
//                 cont_description_parameters_ ),
//      INPUT(),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("cont_descriptionAppend"),
//      DESCRIPTION
//      (
//       "Appends the description of a continuum model or a complete absorption\n"
//       "model to *cont_description_names* and *cont_description_parameters*.\n"
//       "\n"
//       "See online documentation for *cont_description_names* for a list of\n"
//       "allowed models and for information what parameters they require. See\n"
//       "file cont.arts in the doc/examples directory for usage examples and\n"
//       "default parameters for the various models. \n"
//       "\n"
//       "Keywords:\n"
//       "   name       : The name of a continuum model. Must match one of the models\n"
//       "                implemented in ARTS. \n"
//          "   option     : give here the option of this continuum/full model.\n"
//       "   parameters : A Vector containing the required number of parameters\n"
//       "                for the model given. The meaning of the parameters and\n"
//       "                how many parameters are required depends on the model.\n"
//       ),
//      OUTPUT( cont_description_names_, 
//                 cont_description_models_,
//                 cont_description_parameters_ ),
//      INPUT(  cont_description_names_, 
//                 cont_description_models_,
//                 cont_description_parameters_),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "tagname",  "model",   "userparameters" ),
//      TYPES(    String_t,   String_t,   Vector_t         )));


//=== Input Atmosphere methods ===========================================

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("raw_vmrsReadFromFiles"),
//         DESCRIPTION(
//           "Reads the individual VMR profile for each TAGS from file.\n"
//           "\n"
//           "Using this function one can read VMRs of specific TAGS from\n"
//           "explicitly specified files and the remaing from a scenario.\n"
//           "The filenames and the base name of atmospheric scenario\n"
//           "should be specified as keywords. One file name must\n"
//           "be specified for each tag group(each element of *tgs*).\n"
//           "The name may include a path.\n"
//        "\n"
//        "Keywords:\n"
//        "   seltags   : Must be a sub group of tags which should be read from files.\n"
//        "   filenames : Names of the files containing VMR profiles of seltags.\n"
//        "   basename  : The name of a particular atmospheric scenario.\n"
//        "               See *raw_vmrsReadFromScenario* for details. Remaining\n"
//        "               VMRs will be read from the scenario.\n"
//        "\n"
//           ),
//         OUTPUT(   raw_vmrs_         ),
//         INPUT(    tgs_                 ),
//         GOUTPUT(                       ),
//         GINPUT(                        ),
//         KEYWORDS( "seltags",       "filenames",    "basename"),
//         TYPES(    Array_String_t,  Array_String_t, String_t)));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("raw_vmrsReadFromScenario"),
//      DESCRIPTION(
//        "Reads the individual VMR profile for each tag group from a standard\n"
//        "atmospheric scenario.\n" 
//        "\n"
//           "Five different atmospheric scenarios are available in arts data:\n"
//           "tropical, midlatitude-summer, midlatitude-winter, subartic-summer\n"
//           "and subartic-winter.\n"
//        "\n"
//        "   Files in the scenarios look like this: tropical.H2O.aa\n"
//        "\n"
//        "   The basename must include the path, i.e., the files can be anywhere,\n"
//        "but they must be all in the same directory.\n"
//        "   The profile is chosen by the species name. If you have more than one\n"
//        "tag group for the same species, the same profile will be used.\n"
//        "\n"
//        "Keywords:\n"
//        "   basename :The name and path of a particular atmospheric scenario.\n"
//        "   For example:\n"
//        "   /pool/lookup2/arts-data/atmosphere/fascod/tropical\n"
//        "\n"
//        ),
//      OUTPUT(   raw_vmrs_    ),
//      INPUT(    tgs_                 ),
//      GOUTPUT(                       ),
//      GINPUT(                        ),
//      KEYWORDS( "basename"           ),
//      TYPES(    String_t             )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("AtmFromRaw"),
//      DESCRIPTION(
//        "Interpolates temperature, altitude, and VMRs to the pressure grid\n"
//        "given by p_abs.\n" 
//        "\n"
//           "The altitude is not used by the absorption routines,\n"
//        "but later on by the RT routines.\n"
//        "\n"
//        "Interpolations used: \n"
//        "\n"
//        "Temperature      : Linear interpolation in ln(p)\n"
//        "Altitude         : Linear interpolation in ln(p)\n"
//        "VMRs             : Linear interpolation in ln(p)\n"
//        "Cloud Parameters : Linear interpolation in ln(p)\n"
//        "\n"
//        ),
//      OUTPUT(   t_abs_    , z_abs_   , vmrs_           ),
//      INPUT(    tgs_, p_abs_    , raw_ptz_ , raw_vmrs_ ),
//      GOUTPUT(                                         ),         
//      GINPUT(                                          ),
//      KEYWORDS(                                        ),
//      TYPES(                                           )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("WaterVaporSaturationInClouds"),
//      DESCRIPTION(
//        "Calculates the water vapor saturation volume mixing ratio (VMR) in the\n"
//        "vertical range where liquid or ice clouds are in the atmosphere.\n"
//        "At the pressure/altitude grid points where the liquid water content (LWC)\n"
//        "or ice water content (IWC) of the clouds (tags 'liquidcloud' and 'icecloud')\n"
//           "is larger than zero the H2O-VMR is set to liquid water/ice saturation VMR.\n"
//           "The saturation pressure is calculated according to Goff-Gratch equations.\n"
//        ),
//      OUTPUT(   vmrs_ , p_abs_                         ),
//      INPUT(    vmrs_ , p_abs_ , t_abs_ , tgs_         ),
//      GOUTPUT(                                         ),         
//      GINPUT(                                          ),
//      KEYWORDS(                                        ),
//      TYPES(                                           )));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("vmrsScale"),
//      DESCRIPTION(
//           "Scales the vmr input of the tgs given in scaltgs by the\n"
//        "factors given in scalfac.\n"
//        "\n"
//        "Keywords:\n"
//        "   scaltgs : subgroup of tags which has to be scaled.\n"
//        "   scalfac : the factor with which vmr to be scaled.\n"
//        "\n"
//        ),
//      OUTPUT( vmrs_ ),
//      INPUT(  tgs_, vmrs_  ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS( "scaltgs", "scalfac"),
//      TYPES( Array_String_t, Vector_t)));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("h2o_absSet"),
//      DESCRIPTION(
//           "Sets h2o_abs to the profile of the first tag group containing\n"
//        "water.\n" 
//        "\n"
//           "This is necessary, because for example *absCalc* requires h2o_abs\n"
//        "to contain the water vapour profile(the reason for this is the\n"
//           "calculation of oxygen line brodening requires water vapour profile).\n"
//        "Then this function can be used to copy the profile of the first tag\n"
//           "group of water.\n"
//        "\n"
//        ),
//      OUTPUT( h2o_abs_ ),
//      INPUT(  tgs_, vmrs_  ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("n2_absSet"),
//      DESCRIPTION(
//           "Sets n2_abs to the profile of the first tag group containing\n"
//        "molecular nitrogen. See *h2o_absSet* for more details.\n"
//        "\n"
//        ),
//      OUTPUT(     n2_abs_ ),
//      INPUT(  tgs_, vmrs_  ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));




// //=== Absorption methods ===============================================

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME( "absCalc" ),
//      DESCRIPTION(
//         "Calculate absorption coefficients. \n"
//         "\n"
//         "This function calculates both, the total absorption (*abs*)\n"
//         "and the absorption per tag group (*abs_per_tg*).\n"
//             ) ,
//      OUTPUT(abs_  , abs_per_tg_ ),
//      INPUT(tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_, 
//               lines_per_tg_, lineshape_,
//            cont_description_names_, cont_description_models_, 
//               cont_description_parameters_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("absCalcFromXsec"),
//      DESCRIPTION(
//                  "Calculate absorption coefficients from cross sections.\n"
//                  "\n"
//                  "This calculates both the total absorption and the\n"
//                  "absorption per tag group. \n"
//                  "This method calls three other  methods:\n"
//                  "1. *xsec_per_tgInit* - initialize *xsec_per_tg* \n"
//                  "2. *xsec_per_tgAddLine* - calculate cross sections per \n"
//                  "                   tag group for line spectra.\n"
//                  "3. *xsec_per_tgAddConts* - calculate cross sections per \n"
//                  "                   tag group for continua.\n"
//                  "Then it calculates the absorption coefficient by multiplying\n"
//                  "the cross section by VMR.\n"
//                     "This is done once for each tag group (output: *abs_per_tg*)\n"
//                  "and for the sum of all tag group to get the total absorption\n"
//                  "coefficient (output: *abs*)\n"
//                  ),
//      OUTPUT(     abs_  , abs_per_tg_ ),
//      INPUT(      xsec_per_tg_, vmrs_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME( "xsec_per_tgInit" ),
//      DESCRIPTION(
//         "Initialize *xsec_per_tg*.\n"
//         "\n"
//         "The initialization is\n"
//         "necessary, because methods *xsec_per_tgAddLines*\n"
//         "and *xsec_per_tgAddConts* just add to *xsec_per_tg*.\n"
//         "The size is determined from *tgs*.\n"
//         ),
//      OUTPUT( xsec_per_tg_ ),
//      INPUT(tgs_, f_mono_, p_abs_),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("xsec_per_tgAddLines"),
//      DESCRIPTION(
//                  "Calculate cross sections per tag group for line spectra.\n"
//                 ),
//      OUTPUT(     xsec_per_tg_                             ),
//      INPUT(      tgs_, f_mono_, p_abs_, t_abs_, h2o_abs_, vmrs_, 
//                  lines_per_tg_, lineshape_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("xsec_per_tgAddConts"),
//      DESCRIPTION(
//                  "Calculate cross sections per tag group for continua.\n"
//                      ),
//      OUTPUT(     xsec_per_tg_                             ),
//      INPUT(      tgs_, f_mono_, p_abs_, t_abs_, n2_abs_, h2o_abs_, vmrs_,
//                  cont_description_names_, cont_description_parameters_,
//                     cont_description_models_),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));


// //=== Methods operating on absorption ========================================

//   md_data_raw.push_back
//     ( MdRecord
//       ( NAME("abs_per_tgReduce"),
//      DESCRIPTION(
//                  "Reduces absorption coefficients. Only absorption\n"
//                  "coefficients for which weighting functions are\n"
//                  "calculated are kept in memory.\n"
//                  ),
//      OUTPUT(     abs_per_tg_ ),
//      INPUT(      abs_per_tg_, tgs_, wfs_tgs_ ),
//      GOUTPUT(),
//      GINPUT(),
//      KEYWORDS(),
//      TYPES()));



}

