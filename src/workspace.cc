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
       "Dimension: [ f_grid, gas_species ]\n"
       "\n"
       "Case b:    For a single frequency for all species:\n"
       "Dimension: [ 1,      gas_species ]\n"
       "\n"
       "Unit: 1/m"
       ),
      GROUP( Matrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "abs_scalar_gas_field" ),
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
       "Dimensions: [species, f_grid, p_grid, lat_grid, lon_grid]"
        ),
      GROUP( Tensor5_ ))); 

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
       "Dimensions: [f_grid, stokes_dim]"
        ),
       GROUP( Matrix_ )));

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
       "Output of *amp_matCalc*. \n"    
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
       "\n"
       "*amp_mat_raw* is an Array of Array of Tensor6. It contains one \n"
       "gridded field for each particle type which contains the data and \n"
       "also the grids on which the data is stored.\n"
       "For the calculation the data is \n"
       "interpolated on *f_grid*, *scat_za_grid*, *scat_aa_grid*,\n"
       "*scat_za_grid* and *scat_aa_grid*. The interpolated data is stored \n"
       "in *amp_mat*. \n"
       "The amplitude matrix is a 2x2 complex matrix. The workspace variable\n"
       "*amp_mat_raw* stores the real and imaginary elements\n"
       "(i.e. 8 elements) separately. \n" 
       "\n"
       "Usage: Created by *ParticleTypesAdd\n"     
       "\n"
       "Unit: m \n"
       "\n"
       "Size: Array[N_pt] \n" 
       "      Array[7] \n "
       "      [N_f, 1, 1, 1, 1, 1] \n"
       "      [1, N_za, 1, 1, 1, 1] \n"
       "      [1, 1, N_aa, 1, 1, 1] \n"
       "      [1, 1, 1, N_za, 1, 1] \n"
       "      [1, 1, 1, 1, N_aa, 1] \n"
       "      [1, 1, 1, 1, 1, 8] \n"
       "      [N_f, N_za, N_aa, N_za, N_aa, 8] \n"
       "\n"
       ),
      GROUP(ArrayOfArrayOfTensor6_ )));

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
       "Usage:      Set by the user."
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
       "Usage:      Set by the user."
       ),
      GROUP( Index_ )));

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
       "Usage: Set by user."
       ),
      GROUP( ArrayOfMatrix_ )));

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
       "Usage:      Set by the user."
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
     ( NAME( "ext_mat" ),
       DESCRIPTION
      (
       "Total extinction matrix.\n"
       "\n"
       "This variable contains the extinction coefficient matrix which \n"
       "is used in the RTE calculation. It is \n"
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
       "Dimensions: [f_grid, stokes_dim, stokes_dim]"
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
     ( NAME( "f_backend" ),
       DESCRIPTION
       (
        "The frequency grid of the backend channels.\n"
        "\n"
        "Usage:      Input to *sensor_responseBackend*.\n "
        "\n"
        "Unit:       Hz"
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
        "Unit:        Hz"
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
       "\n"
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
        "Unit:       Hz"
        ),
        GROUP( Vector_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "gas_abs_lookup" ),
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
        "once the chapter on the lookup table has been written."
        ), 
       GROUP( GasAbsLookup_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "gas_abs_lookup_is_adapted" ),
       DESCRIPTION
       (
        "Flag to indicate whether *gas_abs_lookupAdapt* has already been\n"
        "called.\n"
        "\n"
        "Values: 0=false, 1=true."
        ), 
       GROUP( Index_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "gas_species" ),
       DESCRIPTION
       (
        "Tag groups for scalar gas absorption.\n"
        "\n"
        "This is an array of arrays of SpeciesTag tag definitions. It defines the\n"
        "available tag groups for the calculation of scalar gas absorption\n"
        "coefficients.  See online documentation of method *gas_speciesSet* for\n"
        "more detailed information how tag groups work and some examples."
        ), 
       GROUP( ArrayOfArrayOfSpeciesTag_ )));

  wsv_data.push_back
    (WsvRecord
    ( NAME( "grid_stepsize" ),
      DESCRIPTION
      (
       "Vector of Numeric with two components:\n"
       "\n"
       "grid_stepsize[0] is the stepsize of the Workspacevariable scat_za_grid,\n"
       "grid_stepsize[1] is the stepsize of the Workspacevariable scat_aa_grid\n"
       "\n"
       "If the value is equal to -1, there is no constant gridspace"
       "\n"
       "Unit: none"
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "ground_emission" ),
       DESCRIPTION
       (
        "The emission from the ground at a specified position.\n"
        "\n"
        "See further *ground_refl_agenda* and the user guide.\n"
        "\n"
        "Usage:      Output from *ground_refl_agenda*.. \n"
        "\n"
        "Unit:       W / (m^2 Hz sr)\n"
        "\n"
        "Dimensions: [ f_grid, stokes_dim ]"
        ), 
       GROUP( Matrix_ )));
 
 wsv_data.push_back
   (WsvRecord
    ( NAME( "ground_refl_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc."
       ),
      GROUP( Agenda_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "ground_refl_coeffs" ),
       DESCRIPTION
       (
        "The reflection coefficients from the directions given by\n"
        "*ground_los* to the direction of interest.\n"
        "\n"
        "The rows and columns of this tensor holds the reflection\n"
        "coefficient matrix for one frequency and one LOS.\n"
        "\n"
        "See further *ground_refl_agenda* and the user guide.\n"
        "\n"
        "Usage:      Output from *ground_refl_agenda*. \n"
        "\n"
        "Units:      -\n"
        "\n"
        "Dimensions: [ ground_los, f_grid, stokes_dim, stokes_dim ]"
        ), 
       GROUP( Tensor4_ )));
 
  wsv_data.push_back
    (WsvRecord
     ( NAME( "ground_los" ),
       DESCRIPTION
       (
        "Directions for which to calculate downwelling radiation when \n"
        "considerin g a ground reflection.\n"
        "\n"
        "See further *ground_refl_agenda* and the user guide.\n"
        "\n"
        "Usage: Output from *ground_refl_agenda*. \n"
        "\n"
        "Units: degrees\n"
        "\n"
        "Size:  [ any number, 1 or 2 ]"
        ), 
       GROUP( Matrix_ )));

  wsv_data.push_back
    (WsvRecord
     ( NAME( "i_montecarlo_error" ),
       DESCRIPTION
       (
        "Error in *I* from ScatteringMonteCarlo.\n"
        "\n"
        "\n"
        "Usage: Output from ScatteringMonteCArlo.. \n"
        "\n"
        "Units: W / (m^2 Hz sr)\n"
        "\n"
        "Size:  [ stokes_dim ]"
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
       "Size: [(cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
       "       (cloudbox_limits[3] - cloudbox_limits[2]) +1, \n"
       "       (cloudbox_limits[5] - cloudbox_limits[4]) +1, \n"
       "        N_za, N_aa, N_i ]"
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
       "Size: [(cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
       "       (cloudbox_limits[3] - cloudbox_limits[2]) +1, \n"
       "       (cloudbox_limits[5] - cloudbox_limits[4]) +1, \n"
       "        N_za, N_aa, N_i ]"
       ),
      GROUP( Tensor6_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "i_rte" ),
      DESCRIPTION
      (
       "One spectrum calculated by *rte_agenda*.\n"
       "\n"
       "This matrix holds the Stokes vector for all frequencies obtained\n"
       "by the monochromatic pencil beam calculations performed by\n"
       "*rte_agenda*. The unit depends on if emission is considered or not."
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
       "Monochromatic intensities at the top of the atmosphere.\n"
       "\n"
       "The matrix holds the intensity entering the atmosphere from above\n"
       "along a propagation path. This should normally correspond to cosmic\n"
       "background radiation, but could also be radiation from the sun or \n"
       "a transmitting satellite.\n"
       "\n"
       "Usage:      Output from *i_space_agenda*.\n"
       "\n"
       "Unit:       W / (m^2 Hz sr) \n"
       "\n"
       "Dimensions: [ f_grid, stokes_dim ]"
       ),
      GROUP( Matrix_ )));

 wsv_data.push_back
    (WsvRecord
     ( NAME( "i_space_agenda" ),
       DESCRIPTION
       (
        "See agendas.cc."
        ),
       GROUP(  Agenda_ )));
  
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
       "Usage:      Calculated in *i_fieldUpdate1D*."
       "\n"
       "Unit:       m \n"
       ),
      GROUP( Numeric_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "iteration_counter" ),
      DESCRIPTION
      (
       "Counter for iterations.\n"
       "\n"
       "This variable holds the number of iterations which have been \n"
       "while solving the RTE with scattering. \n"
       "It is used in the method *Tensor6WriteIteration* and has to be set \n"
       "0 in the control file if this method is used.\n"
       "\n"
  ),
      GROUP( Index_ )));


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
       "Unit:  degrees"
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
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  degrees"
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
    ( NAME( "matrix_1" ),
      DESCRIPTION
      (
       "An arbitrary matrix.\n"
       "\n"
       "This variable shall be treated as a general variable of type Matrix.\n"
       "It can be used, for example, when some intermediate data must be\n"
       "generated or to copy some data.\n"
       "\n"
       "Usage: Set by user."
       ),
      GROUP( Matrix_ )));

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
    ( NAME( "met_profile_calc_agenda" ),
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
       "Unit:  degrees "
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
       "Unit:  degrees "
       ),
      GROUP( Vector_ )));

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
       "Unit:  degrees"
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
       "See documentation of WSM *ybatchMetProfiles* for more information."
       ),
      GROUP( Matrix_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "met_profile_basenames" ),
      DESCRIPTION
      (
       "A list of met office profile basenames.\n"
       "\n"
       "This is intended as input for the method ybatchMetProfiles. It holds a\n"
       "list of profiles basenames. For each basename, there should exist\n"
       "files with different extensions for the different profile parameters,\n"
       "such as humidity, temperature, etc.\n"
       "\n"
       "See documentation of WSM *ybatchMetProfiles* for more information."
       ),
      GROUP( ArrayOfString_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "met_profile_path" ),
      DESCRIPTION
      (
       "Path of the metoffice data.\n"
       "\n"
       "This is intended as input for the method ybatchMetProfiles. It holds \n"
       "the path of the list of profiles. "
       "\n"
       "See documentation of WSM *ybatchMetProfiles* for more information."
       ),
      GROUP( String_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "opt_prop_gas_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc."
       ),
      GROUP( Agenda_)));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "opt_prop_part_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc."
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
       "\n"
       ),
      GROUP( String_ )));

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
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  Pa"
       ),
      GROUP( Vector_ )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "part_types" ),
      DESCRIPTION
      (
       "Particle types.\n"
       "\n"
       "An ArrayOfString containing the filenames for  all particle types\n"
       "which shall be considered. Theses files must be part of the Single \n"
       "Scattering Database. \n"
       "\n"
       "Usage:      Set by the user.\n"
       "\n"
       ),
      GROUP( ArrayOfString_ )));

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
    ( NAME( "pha_mat_spt_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc."
       ),
      GROUP( Agenda_ ))); 


   wsv_data.push_back
   (WsvRecord
    ( NAME( "pha_mat_sptDOITOpt" ),
      DESCRIPTION
      (
       "Interpolated phase matrix for a single particle type.\n"
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
       "Dimensions: [scat_za_grid,scat_aa_grid, scat_za_grid, scat_aa_grid, 6]"
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
       "*pnd_field_raw* is an Array of Array of Tensor3. It contains one \n"
       "gridded field for each particle type which contains the data and \n"
       "also the grids.\n"
       "For the calculation the data is \n"
       "interpolated on *p_grid*, *lat_grid* and *lon_grid*\n"  
       "\n"
       "Usage:      Used in the method *ParticleTypesAdd*.\n"
       "\n"
       "Unit:        m^-3\n"
       "\n"
       "Size:  Array[N_pt]\n"
       "       Array[4] \n "
       "       [N_p, 1, 1] \n"
       "       [1, N_lat, 1] \n"
       "       [1, 1, N_lon] \n"
       "       [N_p, N_lat, N_lon] \n"
       "\n"
       ),
      GROUP( ArrayOfArrayOfTensor3_ )));

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
       "Usage: Output from the method *ppathCalc*."
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
       "Usage:   In/output to/from *ppath_step_agenda*.\n"
       "\n"
       "Members: See AUG."
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
       "Unit: 1"
       ),
      GROUP( Numeric_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "refr_index_agenda" ),
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
    ( NAME( "rte_gp_p" ),
      DESCRIPTION
      (
       "A pressure grid position for radiative transfer calculations.\n"
       "\n"
       "This variable is used to give the grid position for an end point\n"
       "of the propagation path to some workspace method part of the.\n"
       "radiative transfer calculations.\n"
       "\n"
       "Usage:      Set by the calling function, or by the user when calling"
       "            the method directly in the control file."
       ),
      GROUP( GridPos_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_gp_lat" ),
      DESCRIPTION
      (
       "A latitude grid position for radiative transfer calculations.\n"
       "\n"
       "This variable is used to give the grid position for an end point\n"
       "of the propagation path to some workspace method part of the.\n"
       "radiative transfer calculations.\n"
       "\n"
       "Usage:      Set by the calling function, or by the user when calling"
       "            the method directly in the control file."
       ),
      GROUP( GridPos_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "rte_gp_lon" ),
      DESCRIPTION
      (
       "A longitude grid position for radiative transfer calculations.\n"
       "\n"
       "This variable is used to give the grid position for an end point\n"
       "of the propagation path to some workspace method part of the.\n"
       "radiative transfer calculations.\n"
       "\n"
       "Usage:      Set by the calling function, or by the user when calling"
       "            the method directly in the control file."
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
       "Size:  [ 1 or 2 ]"
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
       "Size:  [ atmosphere_dim ]"
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
       "Units: [ Pa ]"
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
       "Units: [ K ]"
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
       "(as given by *gas_species*). It is intended mainly for communication\n"
       "with various methods and agendas, such as methods and agendas \n"
       "calculating absorption coefficients.\n"
       "\n"
       "Usage: Communication variable.\n"
       "\n"
       "Units: [ Absolute value ]\n"
       "\n"
       "Size:  Should match gas_species.nelem()"
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
       "example, *z_ground* and *z_field*. \n"
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
       "Dimensions: [ lat_grid, lon_grid ]"
       ),
      GROUP( Matrix_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "scalar_gas_absorption_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc."
       ),
      GROUP( Agenda_)));
  
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
      ( NAME( "scat_data_raw" ),
        DESCRIPTION
        (
         "Raw data of single scattering data.\n"
         "\n"
         "This variable holds the single scattering properties for all \n"
         "hydrometeor species included in a calculation by using the \n"
         "method *ParticleTypeAdd*.\n"
         "\n"
         ),
        GROUP( ArrayOfSingleScatteringData_ ))); 

  wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_field_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc."
       ),
      GROUP( Agenda_ ))); 


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
       "Size: [(cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
       "       (cloudbox_limits[3] - cloudbox_limits[2]) +1, \n"
       "       (cloudbox_limits[5] - cloudbox_limits[4]) +1, \n"
       "        N_za, N_aa, N_i ]"
       ),
      GROUP( Tensor6_ )));   

 
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
       "Usage:      Output from CloudboxGetIncoming \n"
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
    ( NAME( "scat_mono_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc."
       ),
      GROUP( Agenda_ )));

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
    ( NAME( "scat_theta" ),
      DESCRIPTION
      (
       "Scattering angles.\n"
       "\n"
       "The scattering angles (angle between incident and scattered direction) \n"
       "for all angle combinations of *scat_za_grid* and *scat_aa_grid. \n"
       "The scattering angles are used for the case of randomly oriented \n"
       "particles. \n" 
       "\n"
       "Usage: Output of *single_scattering_dataPrepare* and used in the DOIT\n"
       "       method. \n"
       "\n"
       ),
      GROUP( Tensor4_ ))); 

 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_theta_gps" ),
      DESCRIPTION
      (
       "Grid positions of scattering angles.\n"
       "\n"
       "The scattering angles (angle between incident and scattered direction) \n"
       "for all angle combinations of *scat_za_grid* and *scat_aa_grid. \n"
       "The scattering angles are used for the case of randomly oriented \n"
       "particles. \n" 
       "\n"
       "Usage: Output of *single_scattering_dataPrepare* and used in the DOIT\n"
       "       method. \n"
       "\n"
       ),
      GROUP( ArrayOfArrayOfArrayOfArrayOfGridPos_ ))); 
 
 wsv_data.push_back
   (WsvRecord
    ( NAME( "scat_theta_itws" ),
      DESCRIPTION
      (
       "Interpolation weights of scattering angles.\n"
       "\n"
       "The scattering angles (angle between incident and scattered direction) \n"
       "for all angle combinations of *scat_za_grid* and *scat_aa_grid. \n"
       "The scattering angles are used for the case of randomly oriented \n"
       "particles \n" 
       "\n"
       "Usage: Output of *single_scattering_dataPrepare* and used in the DOIT\n"
       "       method. \n"
       "\n"
       ),
      GROUP( Tensor5_ ))); 


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
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
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
       "Size:  [ number of polarisation values, stokes_dim ]"
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
       "Size:  [ number of measurement blocks, atmosphere_dim ]"
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
	"Usage:		Output/input to the *sensor_response...* methods.\n"
	"\n"
	"Units:		-\n"
	"\n"
	"Dimension:     See the individual *sensor_response...* method \n"
        "documentation."
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
       "Unit:  [ Hz ]"
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
       "Unit:  [ degrees ]"
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
       "Unit:  [ degrees ]"
       ),
      GROUP( Vector_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "sensor_rot" ),
      DESCRIPTION
      (
       "The rotation of the sensor for each measurement block.\n"
       "\n"
       "The rotation is the angle between the atmospheric and sensor frames\n"
       "for polarisation. The angle increases with clockwise rotation of the\n"
       "sensor when looking along the line-of-sight of the sensor. \n"
       "\n"
       "If the purpose of the simulations is to extract the polarisation\n"
       "of the radiation coming from the atmosphere (no sensor), the angles\n"
       "shall be set to 0.\n"
       "\n"       
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "Unit:  [ degrees ]\n"
       "\n"
       "Size:  [ number of measurement blocks ]"
       ),
      GROUP( Vector_ )));

wsv_data.push_back
   (WsvRecord
    ( NAME( "single_scattering_data" ),
      DESCRIPTION
      (
       "Structure for the  single scattering data.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user.\n"
       "\n"
       "More will be written (Claudia)\n"
       "\n"
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
       "Usage: Set by user."
       ),
      GROUP( Sparse_ )));

  wsv_data.push_back
   (WsvRecord
    ( NAME( "stokes_dim" ),
      DESCRIPTION
      (
       "The dimensionality of the Stokes vector (1-4).\n"
       "\n"
       "Usage:      Set by the user."
       ),
      GROUP( Index_ )));

   wsv_data.push_back
   (WsvRecord
    ( NAME( "spt_calc_agenda" ),
      DESCRIPTION
      (
        "See agendas.cc."
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
       "Usage:      Calculated internally."
       ),
      GROUP( Vector_ )));

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
       "Usage: Set by user."
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
       "Usage: Set by user."
       ),
      GROUP( Tensor4_ )));

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
       "Usage: Set by user."
       ),
       GROUP( Tensor6_ )));

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
       "Dimensions: [ p_grid, lat_grid, lon_grid ]"
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
       "Size   Array[4] \n "
       "       [N_p, 1, 1] \n"
       "       [1, N_lat, 1] \n"
       "       [1, 1, N_lon] \n"
       "       [N_p, N_lat, N_lon] \n"
       "\n"
       ),
      GROUP( ArrayOfTensor3_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "timer" ),
      DESCRIPTION
      (
       "Stores the starting time for time measurements.\n"
       "\n"
       ),
      GROUP( Timer_ )));

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
       "Usage: Set by user."
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
       "Usage: Set by user."
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
       "Usage: Set by user."
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
       "       Array[4] \n "
       "       [N_p, 1, 1] \n"
       "       [1, N_lat, 1] \n"
       "       [1, 1, N_lon] \n"
       "       [N_p, N_lat, N_lon] \n"
       "\n"
       ),
      GROUP( ArrayOfArrayOfTensor3_ )));


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
       "Usage: Output from radiative transfer calculations considering\n"
       "       sensor response.\n"
       "\n"
       "Unit:  Undefined. Possibilities include: K, W/(m^2 Hz sr) and\n "
       "       optical thickness."
       ),
      GROUP( Vector_ )));

 wsv_data.push_back
   (WsvRecord
    ( NAME( "ybatch" ),
      DESCRIPTION
      (
       "Spectra for a batch of metoffice profiles.\n"
       "The spectra is in the columns of the matrix and the rows\n"
       "corresponds to the number of profiles.\n"
       "\n"
       "Usage: Output from batch methods.\n"
       "\n"
       "Unit:  Undefined. Possibilities include: K, W/(m^2 Hz sr) and\n "
       "       optical thickness."
       ),
      GROUP( Matrix_ )));

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
       "Dimensions: [ p_grid, lat_grid, lon_grid ]"
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
       "Size   Array[4] \n "
       "       [N_p, 1, 1] \n"
       "       [1, N_lat, 1] \n"
       "       [1, 1, N_lon] \n"
       "       [N_p, N_lat, N_lon] \n"
       "\n"
       ),
      GROUP( ArrayOfTensor3_ )));

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
       "crossings is obtained by a linear (1D) or bi-linear (2D) \n"
       "interpolation of the sum of *r_geoid* and *r_ground*. With other \n"
       "words, the radius for the ground is assumed to vary linear along the\n"
       "latitudes and longitudes in *lat_grid* and *lon_grid*.\n"
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
 
}
