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
#include <cmath>
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
void scat_data_meta_arrayInit(// WS Output:
                 ArrayOfScatteringMetaData& scat_data_meta_array,
                 //WS Input
                 const Verbosity&)
{
    scat_data_meta_array.resize(0);
}
                
                 
void scat_data_meta_arrayAddTmatrix(// WS Output:
                             ArrayOfScatteringMetaData& scat_data_meta_array,
                             // WS Generic input:
                             const String& description,
                             const String& material,
                             const String& shape,
                             const String& particle_type,
                             const Numeric& density,
                             const Numeric& aspect_ratio,
                             const Vector& diameter_grid,
                             const Vector& scat_f_grid,
                             const Vector& scat_T_grid,
                             const GriddedField3& complex_refr_index,
                             const Verbosity&) 

{
  chk_if_equal("f_gridScat", "f_grid from complex_refr_index", scat_f_grid, complex_refr_index.get_numeric_grid(0));
  chk_if_equal("T_gridScat", "T_grid from complex_refr_index", scat_T_grid, complex_refr_index.get_numeric_grid(1));

  for(Index k=0; k<diameter_grid.nelem(); k++)
    {
      extern const Numeric PI;

      Numeric diameter_max;

      if (shape == "spheroidal")

        if (aspect_ratio<1)
          diameter_max=diameter_grid[k]*pow(aspect_ratio, -2./3.);

        else if (aspect_ratio>1)
          diameter_max=diameter_grid[k]*pow(aspect_ratio, 1./3.);

        else
          {
            ostringstream os;
            os << "Incorrect aspect ratio: " << aspect_ratio << "\n"
              << "Can not be equal to one";
            throw std::runtime_error(os.str());
          }

      else if (shape == "cylindrical")

          diameter_max=pow(pow(16./3., 2./3.)*pow(diameter_grid[k]/.2, 2.)*
                    (pow(aspect_ratio, -4./3.)+pow(aspect_ratio, 2./3.)), 1./2.);
 

      else
        {
          ostringstream os;
          os << "Unknown particle shape: " << shape << "\n"
            << "Must be spheroidal or cylindrical";
          throw std::runtime_error(os.str());
        }

      ScatteringMetaData smd;
      if (description=="")
        {   
          ostringstream os;
          os << shape<< " "<< material << " particle of type " << particle_type<<
            ", with volume equivalent diameter "
            <<diameter_grid[k]<<" meters.";
          smd.description=os.str();
        }
      else 
        smd.description = description;
 
      smd.material = material;
      smd.shape = shape;
      smd.particle_type = particle_type;
      smd.density = density;
      smd.diameter_max =diameter_max;
      smd.volume = 4./3.*PI*pow(diameter_grid[k]/2., 3);
      smd.area_projected = 0;
      smd.aspect_ratio = aspect_ratio;
      smd.scat_f_grid = scat_f_grid;
      smd.scat_T_grid = scat_T_grid;
      smd.complex_refr_index = complex_refr_index.data;

      scat_data_meta_array.push_back(smd);
    }
}

//-----------------------------------


void scat_data_rawFromMeta(// WS Output:
                              ArrayOfSingleScatteringData& scat_data_raw,
                              //WS Input
                              const ArrayOfScatteringMetaData& scat_data_meta_array,
                              // WS Generic input:
                              const Vector& za_grid,
                              const Vector& aa_grid,
                              const Numeric& precision,
                              const Verbosity&)
{
  for(Index ii=0;ii<scat_data_meta_array.nelem();ii++)
    {

      extern const Numeric PI;  
      const String& particle_type = scat_data_meta_array[ii].particle_type;
      Index  np;

      SingleScatteringData sdd;
      sdd.f_grid = scat_data_meta_array[ii].scat_f_grid;
      sdd.T_grid = scat_data_meta_array[ii].scat_T_grid;
      sdd.za_grid = za_grid;
      sdd.aa_grid = aa_grid;


      if (particle_type == "macros_iso")
        sdd.particle_type = PARTICLE_TYPE_MACROS_ISO;
      else if (particle_type == "horiz_al")
        sdd.particle_type = PARTICLE_TYPE_HORIZ_AL;
      else
        {
          ostringstream os;
          os << "Unknown particle type: " << particle_type << "\n"
            << "Must be macros_iso or horiz_al";
          throw std::runtime_error(os.str());
        }

      if (scat_data_meta_array[ii].shape == "spheroidal" )
        np=-1;

      else if (scat_data_meta_array[ii].shape == "cylindrical")
        np=-2;
      else
        {
          ostringstream os;
          os << "Unknown particle shape: " << scat_data_meta_array[ii].shape << "\n"
            << "Must be spheroidal or cylindrical";
          throw std::runtime_error(os.str());
        }

      calcSingleScatteringDataProperties(sdd,
                                         scat_data_meta_array[ii].complex_refr_index(joker,joker,0),
                                         scat_data_meta_array[ii].complex_refr_index(joker,joker,1),
                                         pow(scat_data_meta_array[ii].volume*1e18*3./(4.*PI),1./3.),
                                         np,
                                         scat_data_meta_array[ii].aspect_ratio,
                                         precision);

      scat_data_raw.push_back(sdd);
    }
}


void single_scattering_dataFromTmatrix( 
         SingleScatteringData& single_scattering_data,
   const Vector&        f_grid,
   const Vector&        T_grid,
   const Vector&        za_grid,
   const Vector&        aa_grid,
   const GriddedField3& complex_refr_index,
   const Numeric&       equiv_radius,
   const String&        orientation,
   const String&        shape,
   const Numeric&       aspect_ratio,
   const Numeric&       precision,
   const Verbosity&)
{
    single_scattering_data.f_grid  = f_grid;
    single_scattering_data.T_grid  = T_grid;
    single_scattering_data.za_grid = za_grid;
    single_scattering_data.aa_grid = aa_grid;

    if( orientation == "macros_iso" )
      { single_scattering_data.particle_type = PARTICLE_TYPE_MACROS_ISO; }
    else if( orientation == "horiz_al" )
      { single_scattering_data.particle_type = PARTICLE_TYPE_HORIZ_AL; }
    else
      {
        ostringstream os;
        os << "Unknown particle orientation: " << orientation << "\n"
        << "Must be \"macros_iso\" or \"horiz_al\".";
        throw std::runtime_error(os.str());
      }

    Index ishape = 999;
    if( shape == "spheroidal" )
      { ishape = -1; }
    else if( shape == "cylindrical" )
      { ishape = -2; }
    else
      {
        ostringstream os;
        os << "Unknown particle shape: " << shape << "\n"
        << "Must be \"spheroidal\" or \"cylindrical\".";
        throw std::runtime_error(os.str());
      }
    
    // Extract real and imaginary part of n
    // So far we just set to first value
    const Index nf = f_grid.nelem();
    const Index nt = T_grid.nelem();
    //
    Matrix n_real(nt,nf), n_imag(nt,nf);
    //
    n_real = complex_refr_index.data(0,0,0);
    n_imag = complex_refr_index.data(0,0,1);

    calcSingleScatteringDataProperties( single_scattering_data,
                                        n_real, n_imag, equiv_radius,
                                        ishape, aspect_ratio,
                                        precision );
}
