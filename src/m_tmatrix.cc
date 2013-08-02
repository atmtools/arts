/* Copyright (C) 2013 Oliver Lemke
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*!
  \file   m_tmatrix.cc
  \author Oliver Lemke
  \date   2013-06-25

  \brief  T-Matrix related workspace methods.
*/

#include <stdexcept>
#include "messages.h"
#include "tmatrix.h"
#include "check_input.h"


void TMatrixTest(const Verbosity& verbosity)
{
    tmatrix_tmd_test(verbosity);
    tmatrix_ampld_test(verbosity);
    calc_ssp_random_test(verbosity);
    calc_ssp_fixed_test(verbosity);
}
//-----------------------------

void scat_data_meta_arrayAdd(// WS Output:
                 ArrayOfScatteringMetaData& scat_data_meta_array,
                 // WS Generic input:
                 const String& description,
                 const String& material,
                 const String& shape,
                 const String& p_type,
                 const Numeric& density,
                 const Numeric& aspect_ratio,
                 const Vector& r_grid,
                 const Vector& f_gridScat,
                 const Vector& T_gridScat,
                 const GriddedField3& ref_index,
                 const Verbosity&) 

{ chk_if_equal("f_gridScat", "f_grid from ref_index", f_gridScat, ref_index.get_numeric_grid(0));
  chk_if_equal("T_gridScat", "T_grid from ref_index", T_gridScat, ref_index.get_numeric_grid(1));
    
    for(Index k=0; k<r_grid.nelem(); k++)
        
        
    {

        //Works for spheres only atm
        
        extern const Numeric PI;
        
       
        
        
        ScatteringMetaData smd;
        smd.description = description;
        smd.type = material;
        smd.shape = shape;
        smd.p_type = p_type;
        smd.density = density;
        smd.d_max =2*r_grid[k];
        smd.V = 4/(float)3*PI*r_grid[k]*r_grid[k]*r_grid[k];
        smd.A_projec = 2*r_grid[k]*PI;
        smd.asratio = aspect_ratio;
        smd.f_grid = f_gridScat;
        smd.T_grid = T_gridScat;
        smd.ref_index = ref_index.data;
        
        
        
        scat_data_meta_array.push_back(smd);
    }
    
    
    
}

//-----------------------------------


void single_scattering_dataCalcTMatrixTest(// WS Output:
                                           ArrayOfSingleScatteringData& scat_data_raw,
                                           // WS Generic input:
                                           const String& p_type,
                                           const Vector& f_grid,
                                           const Vector& T_grid,
                                           const Vector& za_grid,
                                           const Vector& aa_grid,
                                           const Matrix& ref_index_real,
                                           const Matrix& ref_index_imag,
                                           const Vector& equiv_radius,
                                           const Index& np,
                                           const String& phase,
                                           const Numeric& aspect_ratio,
                                           const Numeric& precision,
                                           const Verbosity&)
{
    SingleScatteringData sdd;
    sdd.f_grid = f_grid;
    sdd.T_grid = T_grid;
    sdd.za_grid = za_grid;
    sdd.aa_grid = aa_grid;

    if (p_type == "MACROS_ISO")
        sdd.ptype = PARTICLE_TYPE_MACROS_ISO;
    else if (p_type == "HORIZ_AL")
        sdd.ptype = PARTICLE_TYPE_HORIZ_AL;
    else
    {
        ostringstream os;
        os << "Unknown particle type: " << p_type << "\n"
        << "Must be MACROS_ISO or HORIZ_AL";
        throw std::runtime_error(os.str());
    }
    
    for(Index ii=0;ii<equiv_radius.nelem();ii++)
    {
        calcSingleScatteringDataProperties(sdd,
                                       ref_index_real,
                                       ref_index_imag,
                                       equiv_radius[ii],
                                       np,
                                       phase,
                                       aspect_ratio,
                                       precision);
        scat_data_raw.push_back(sdd);
    }
}
