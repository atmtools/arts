/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>
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
#include "vecmat.h"
#include "wsv_aux.h"

/*! The lookup information for the workspace variables. */
ARRAY<WsvRecord> wsv_data;

void define_wsv_data()
{

  //--------------------< Build the wsv data >--------------------
  // Initialize to empty, just in case.
  wsv_data.clear();


  //--------------------< Spectroscopy Stuff >--------------------
  //                     --------------------
  wsv_data.push_back
    (WsvRecord
     ("lines",
      "A list of spectral line data.", 
      ARRAYofLineRecord_));

  wsv_data.push_back
    (WsvRecord
     ("lines_per_tg",
      "A list of spectral line data for each tag.\n"
      "Dimensions: (tag_groups.dim()) (# of lines for this tag)", 
      ARRAYofARRAYofLineRecord_));

  wsv_data.push_back
    (WsvRecord
     ("tag_groups",
      "This is an array of arrays of OneTag tag definitions.\n"
      "It defines the available tag groups for the calculation\n"
      "of absorption coefficients and weighting functions.\n"
      "Contrary to the original Bredbeck definition, tags within a\n"
      "group must belong to the same species, because one VMR profile\n"
      "is associated with each tag group.", 
      TagGroups_));

  wsv_data.push_back
    (WsvRecord
     ("lineshape",
      "Indices of lineshape functions. There is one entry for\n"
      "each abs_tag, not for each species. This means if you have several\n"
      "abs_tags for different isotopes or transitions of a species, you\n"
      "may use different lineshapes.",
      ARRAYofsizet_));

  wsv_data.push_back
    (WsvRecord
     ("lineshape_norm",
      "Indices of normalizations to the lineshape functions. There is one\n"
      "entry for each abs_tag, not for each species. This means if you have\n"
      "several abs_tags for different isotopes or transitions of a species, you\n"
      "may use different lineshapes and normalizations.",
      ARRAYofsizet_));


  //--------------------< 1D Input Atmosphere Stuff >--------------------
  //                     ---------------------------
  wsv_data.push_back
    (WsvRecord
     ("raw_ptz_1d",
      "Matrix has rows:\n"
      "1. Pressure in Pa\n"
      "2. Temperature in K\n"
      "3. Altitude in m", 
      MATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("raw_vmrs_1d",
      "The individual VMR profiles. Each species VMR profile comes with a\n"
      "pressure profile. The different species can hence be on different\n"
      "grids.\n"
      "The matrix has rows:\n"
      "1. Pressure in Pa\n"
      "2. VMR profile (absolute number)\n"
      "The array dimension is determined by the number of tag groups.", 
      ARRAYofMATRIX_));


  //--------------------< General Absorption Stuff >--------------------
  //                     --------------------------
  wsv_data.push_back
    (WsvRecord
     ("p_abs",
      "The pressure grid for the absorption coefficients [Pa]. This\n"
      "is the basic independent grid, both in the 1D and 2D\n"
      "case. Therefore it remains a vector, even in 2D.",
      VECTOR_));
  
  wsv_data.push_back
    (WsvRecord
     ("f_mono",
      "The monochromatic frequency grid [Hz]. This grid is used when\n"
      "calculating absorption and pencil b",
      VECTOR_));
    

  //--------------------< 2D Absorption Stuff >--------------------
  //                     ---------------------
  wsv_data.push_back
    (WsvRecord
     ("t_abs_2d",
      "2D temperatures associated with the pressures in p_abs [K].\n"
      "Array coordinate is the profile index, i.e., the horizontal\n"
      "dimension. This dimension must be consistent with z_abs_2d,\n"
      "vmr_2d, and abs_2d.", 
      ARRAYofVECTOR_));

  wsv_data.push_back
    (WsvRecord
     ("z_abs_2d",
      "2D vertical altitudes associated with the pressures in p_abs [m]."
      "Array coordinate is the profile index, i.e., the horizontal\n"
      "dimension. This dimension must be consistent with t_abs_2d,\n"
      "vmr_2d, and abs_2d.", 
      ARRAYofVECTOR_));

  wsv_data.push_back
    (WsvRecord
     ("vmrs_2d",
      "2D VMRs associated with the pressures in p_abs [absolute number]."
      "Array coordinate is the profile index, i.e., the horizontal\n"
      "dimension. This dimension must be consistent with t_abs_2d\n"
      "z_abs_2d, and abs_2d."
      "\n"
      "The matrix dimensions are [tag_groups.dim(),p_abs.dim()].", 
      ARRAYofMATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("abs_2d",
      "The array of absorption coefficient matrices."
      "Array coordinate is the profile index, i.e., the horizontal\n"
      "dimension. This dimension must be consistent with t_abs_2d\n"
      "and z_abs_2d."
      "\n"
      "The matrix dimensions are [f_mono.dim(),p_abs.dim()].", 
      ARRAYofMATRIX_));


  //--------------------< 1D Absorption Stuff >--------------------
  //                     ---------------------
  wsv_data.push_back
    (WsvRecord
     ("t_abs",
      "Temperature associated with the pressures in p_abs [K]",
      VECTOR_));

  wsv_data.push_back
    (WsvRecord
     ("z_abs",
      "Vertical altitudes associated with the pressures in p_abs [m]",
      VECTOR_));

  wsv_data.push_back
    (WsvRecord
     ("vmrs",
      "The VMRs (unit: absolute number) on the p_abs grid.\n"
      "Dimensions: [tag_groups.dim(), p_abs.dim()]",
      ARRAYofVECTOR_));

  wsv_data.push_back
    (WsvRecord
     ("abs",
      "The matrix of absorption coefficients (in units of [1/m]).\n"
      "Dimensions: [f_mono.dim(), p_abs.dim()]",
      MATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("abs_per_tg",
      "These are the absorption coefficients individually for each\n"
      "tag group. The ARRAY contains one matrix for each tag group,\n"
      "the matrix format is the same as that of abs",
      ARRAYofMATRIX_));


  //--------------------< RT Stuff >--------------------
  //                     ----------
  wsv_data.push_back
    (WsvRecord
     ("za_pencil",
      "Pencil beam zenith angle, the angle between zenith and the LOS [deg].\n"
      "This grid is applied when calculating pencil beam spectra.",
      VECTOR_));

  wsv_data.push_back
    (WsvRecord
     ("z_plat",
      "The vertical altitude, above the geiod, of the platform [m].",
      Numeric_));

  wsv_data.push_back
    (WsvRecord
     ("l_step",
      "The maximum length, along the LOS, between the points of LOS [m].\n"
      "The final step length will in most cases equal the selected length.\n"
      "There are two rare exceptions:\n"
      "  1. Downward observations from within the atmsophere, where the step\n"
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
      "Boolean to consider refraction (0=no refraction).",
      int_));

  wsv_data.push_back
    (WsvRecord
     ("l_step_refr",
      "The step length (along LOS) when performing the calculations to\n"
      "determining the LOS with refraction [m].\n"
      "Note that the final step length between the LOS points is l_step.\n"
      "The step length here is only applied during the calculations.",
      Numeric_));

  wsv_data.push_back
    (WsvRecord
     ("refr_index",
      "The refractive index associated with the pressures in p_abs [-].\n",
      VECTOR_));

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
      VECTOR_));

  wsv_data.push_back
    (WsvRecord
     ("los",
      "Structure to define the line of sight (LOS). See los.h.", 
      LOS_));

  wsv_data.push_back
    (WsvRecord
     ("source",
      "Mean source functions between the points of the LOS.",
      ARRAYofMATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("trans",
      "The transmissions between the points of the LOS [-].",
      ARRAYofMATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("y_space",
      "Radiation entering the atmosphere at the start of the LOS,\n"
      "typically cosmic background radiation.",
      VECTOR_));

  wsv_data.push_back
    (WsvRecord
     ("y",
      "The working spectrum.",
      VECTOR_));

  wsv_data.push_back
    (WsvRecord
     ("y0",
      "A reference spectrum. This variable can be used e.g. to save a copy\n"
      "of y or to compare the spectra before and after some operation(s).",
      VECTOR_));


  //--------------------< WF Stuff >--------------------
  //                     ----------
  wsv_data.push_back
    (WsvRecord
     ("absloswfs",
      "Line of sight weighting functions.",
      ARRAYofMATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("k_grid",
      "Grid for the retrieval identity for which weighting functions (WFS)\n"
      "shall be calculated (when applicable). For example, pressure altitude\n"
      "grid for species\n"
      "WFs.",
      VECTOR_));

  wsv_data.push_back
    (WsvRecord
     ("k",
      "The weighting functions (WFs) for a single retrieval identity.",
      MATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("k_names",
      "Name(s) on the retrieval identity associated with k.",
      ARRAYofstring_));

  wsv_data.push_back
    (WsvRecord
     ("k_aux",
      "Auxiliary data for k. The number of rows of this matrix equals the\n"
      "length of the state vector for the retrieval identity (the number of\n"
      "columns of k).\n"
      "The columns hold different quantities:\n"
      "  Col 1: retrieval grid (or correspondingly)\n"
      "  Col 2: a priori values\n"
      "  Col 3: volume mixing ratios",
      MATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("kx",
      "The state weighting function matrix.",
      MATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("kx_names",
      "Names on the retrieval identities associated with kx.",
      ARRAYofstring_));

  wsv_data.push_back
    (WsvRecord
     ("kx_index",
      "This is a two-column matrix holding first and last index of the state\n"
      "vector for each retrieval identity. That is, each row corresponds to\n"
      "a retrieval identity as [i_first,i_last].",
      MATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("kx_aux",
      "Auxiliary data for kx. As k_aux but with the data of the different\n"
      "retrieval identies appended vertically.",
      MATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("kb",
      "The model parameters weighting function matrix.",
      MATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("kb_names",
      "Names on the model parameter identities associated with kb.",
      ARRAYofstring_));

  wsv_data.push_back
    (WsvRecord
     ("kb_index",
      "This is a two-column matrix holding first and last index of the model\n"
      "parameter vector (b) for each  identity. That is, each row \n"
      "corresponds to a forward model identity as [i_first,i_last].",
      MATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("kb_aux",
      "Auxiliary data for kb. As k_aux but with the data of the different\n"
      "forward model identies appended vertically.",
      MATRIX_));


  //--------------------< H matrices >--------------------
  //                     ------------
  wsv_data.push_back
    (WsvRecord
     ("h1",
      "A transfer matrix for sensor effects and data reduction.\n"
      "This is typically the total transfer matrix and includedes effects\n"
      "of both the sensor and data reduction.\n"
      "This matrix and h1 can also be used to split the sensor calculations\n"
      "in two parts.",
      Hmatrix_));

  wsv_data.push_back
    (WsvRecord
     ("h2",
      "A second transfer matrix for sensor effects and data reduction.\n"
      "This matrix includes typically only effects of the data reduction.\n"
      "See further h1.",
      Hmatrix_));


  //  cout << "size = " << wsv_data.size() << '\n';
}
