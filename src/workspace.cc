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

  <b>History:</b>
  <dl>
  <dt> 2000-08-10 Stefan Buehler</dt>
  <dd> Removed the wsv pointers. They are now in a separate place.</dd>
  </dl>
*/

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
      "pressure profile. The different species can hence be on different grids.\n"
      "\n"
      "Matrix has rows:\n"
      "1. Pressure in Pa\n"
      "2. VMR profile (absolute number)\n"
      "\n"
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
      "This grid is applied when calculating pecil beam spectra",
      VECTOR_));

  wsv_data.push_back
    (WsvRecord
     ("z_plat",
      "The vertical altitude, above the geiod, of the platform [m].",
      Numeric_));

  wsv_data.push_back
    (WsvRecord
     ("l_step",
      "The length (along the LOS) between the points of LOS [m].",
      Numeric_));

  wsv_data.push_back
    (WsvRecord
     ("refr",
      "Boolean to consider refraction (0=no refraction).",
      int_));

  wsv_data.push_back
    (WsvRecord
     ("l_step_refr",
      "The step length (along LOS) when determining the LOS with refraction [m].",
      Numeric_));

  wsv_data.push_back
    (WsvRecord
     ("refr_index",
      "The refractive index associated with the pressures in p_abs [-].",
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
      "Structure to define the line of sight (LOS) for 1d cases.", 
      Los_));

  wsv_data.push_back
    (WsvRecord
     ("source",
      "Mean source functions between the points of the LOS [W/(m3Hzsr)].",
      ARRAYofMATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("trans",
      "The transmissions between the points of the LOS [-].",
      ARRAYofMATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("y_space",
      "Radiation entering the atmosphere at the start of the LOS.",
      VECTOR_));

  wsv_data.push_back
    (WsvRecord
     ("y",
      "The final spectrum, including sensor effects and data reduction.",
      VECTOR_));


  //--------------------< WF Stuff >--------------------
  //                     ----------
  wsv_data.push_back
    (WsvRecord
     ("klos",
      "Line of sight weighting functions.",
      ARRAYofMATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("k_grid",
      "Grid for the weighting function matrix to be calculated.",
      VECTOR_));

  wsv_data.push_back
    (WsvRecord
     ("k",
      "A weighting function matrix.",
      MATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("h",
      "The instrument matrix. Patrick, please put more information here.",
      SPARSEMATRIX_));

  wsv_data.push_back
    (WsvRecord
     ("lineshapes",
      "Indices of lineshape functions. There is one entry for\n"
      "each abs_tag, not for each species. This means if you have several\n"
      "abs_tags for different isotopes or transitions of a species, you\n"
      "may use different lineshapes.",
      ARRAYofsizet_));


  //  cout << "size = " << wsv_data.size() << '\n';
}
