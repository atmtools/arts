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

#include <cfloat>
#include <cmath>
#include <stdexcept>
#include "check_input.h"
#include "logic.h"
#include "make_array.h"
#include "math_funcs.h"
#include "messages.h"
#include "refraction.h"
#include "special_interp.h"
#include "tmatrix.h"

extern const Numeric PI;  



/* Workspace method: Doxygen documentation will be auto-generated */
void diameter_maxFromDiameter_volume_equ(
          Numeric&   diameter_max,
          Numeric&   aspect_area_max,
    const String&    shape,
    const Numeric&   diameter_volume_equ,
    const Numeric&   aspect_ratio,
    const Verbosity& )
{
  const Numeric volume = ( PI*pow(diameter_volume_equ,3) ) / 6.;

  if( shape == "spheroidal" )
    {
      if( aspect_ratio > 1 ) // oblate spheriod
        {   
          //rotational axis
          const Numeric a = pow( (3.*volume*aspect_ratio)/(4*PI), 1./3. );
          diameter_max   = 2.*a;
          aspect_area_max = PI * a * a;
        }
      else if( aspect_ratio < 1) // prolate spheroid
        {  
          //non-rotational axis (perpendicular to a)
          const Numeric b = pow( (3.*volume)/(4.*PI*pow(aspect_ratio,2)), 1./3.);
          diameter_max   = 2.*b;
          aspect_area_max = PI * b * b;
        }
      else
        {
          throw runtime_error( "For spheriodal particles, the aspect ratio "
                               "is not allowed to be exactly 1 (due to "
                               "potential numerical problems)." );
        }
    }
  
  else if (shape == "cylindrical")
    {
      //aspect_ratio=D/L
      const Numeric D = pow( (volume*4*aspect_ratio)/PI, 1./3. );
      const Numeric L = D / aspect_ratio;
      diameter_max   = pow( pow(D,2) + pow(L,2), 1./2. );
      aspect_area_max = max( PI*D*D/4., D*L );
    }
  
  else
    {
      ostringstream os;
      os << "Unknown particle shape: " << shape << "\n"
         << "Must be spheroidal or cylindrical";
      throw runtime_error( os.str() );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void diameter_volume_equFromDiameter_max(
          Numeric&   diameter_volume_equ,
          Numeric&   volume,
    const String&    shape,
    const Numeric&   diameter_max,
    const Numeric&   aspect_ratio,
    const Verbosity&)
{
  if( shape == "spheroidal" )
    {
      if( aspect_ratio > 1 ) // oblate spheriod
        {   
          const Numeric a = diameter_max/2.;
          volume = (pow(a,3.)*4.*PI)/(3.*aspect_ratio);
        }
      else if( aspect_ratio < 1 ) // prolate spheroid
        {  
          const Numeric b = diameter_max/2.;
          volume = (pow(b,3.)*4.*PI*pow(aspect_ratio,2.))/3.;
        }
      else
        {
          throw runtime_error( "For spheriodal particles, the aspect ratio "
                               "is not allowed to be exactly 1 (due to "
                               "numerical problems)." );
        }
    }
  
  else if (shape == "cylindrical")
    {
      volume = pow( diameter_max/pow((pow(aspect_ratio,2.)+1.),1./2.),3)*
               pow(aspect_ratio,2.)*PI/4.;
    }
  
  else
    {
      ostringstream os;
      os << "Unknown particle shape: " << shape << "\n"
         << "Must be spheroidal or cylindrical";
      throw runtime_error( os.str() );
    }

  diameter_volume_equ = pow( (6.*volume) / PI, 1./3. );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void scat_data_singleTmatrix(
          SingleScatteringData&   scat_data_single,
          ScatteringMetaData&     scat_meta_single,
    const GriddedField3&          complex_refr_index,
    const String&                 shape,
    const Numeric&                diameter_volume_equ,
    const Numeric&                aspect_ratio,
    const Numeric&                mass,
    const String&                 ptype,
    const Vector&                 data_f_grid,
    const Vector&                 data_t_grid,
    const Vector&                 data_za_grid,
    const Vector&                 data_aa_grid,
    const Numeric&                precision,
    const String&                 cri_source,
    const Index&                  ndgs,
    const Index&                  robust,
    const Index&                  quiet,
    const Verbosity&              verbosity )
{
  // Get internal coding for ptype
  scat_data_single.ptype = PTypeFromString( ptype );

  // Set description
  {   
    ostringstream os;
    os << "T-matrix calculation for a " << shape << " particle, with " 
       << "diameter_volume_equ = " << 1e6*diameter_volume_equ << "um and "
       << "aspect ratio = " << aspect_ratio << ".";
    scat_data_single.description = os.str();
  }


  // Add grids to scat_data_single
  //
  scat_data_single.f_grid  = data_f_grid;
  scat_data_single.T_grid  = data_t_grid;

  if( scat_data_single.ptype == PTYPE_MACROS_ISO )
    {
      // tmatrix random orient requires equidistant angular grid. checking here
      // that given data_za_grid fulfills this requirement
      if( !(is_same_within_epsilon(data_za_grid[0],0.,2*DBL_EPSILON) &&
            is_same_within_epsilon(last(data_za_grid),180.,2*DBL_EPSILON)) )
        {
          ostringstream os;
          os << "Zenith angle (=scattering angle) grid needs to include\n"
             << "0 deg and 180 deg as first and last grid points, respectively.\n"
             << "At least one of them does not fit.";
          throw runtime_error( os.str() );
        }
      Index nza = data_za_grid.nelem();
      Numeric dza = 180./((Numeric)nza-1.);
      for( Index iza=1; iza<nza; iza++ )
        {
          if( !(is_same_within_epsilon(data_za_grid[iza],(Numeric)iza*dza,2*DBL_EPSILON)) )
            {
              ostringstream os;
              os << "Input zenith angle grid *data_za_grid* is required to be\n"
                 << "equidistant for randomly oriented particles, but it is not.";
              throw runtime_error( os.str() );
            }
        }
    }
  scat_data_single.za_grid = data_za_grid;

  if( scat_data_single.ptype == PTYPE_MACROS_ISO )
    {
      // in case of random orientation, azimuth grid should be empty. We just
      // set that here, ignoring whatever is in data_aa_grid (but giving a
      // warning that we do so).
      Vector empty_grid(0);
      scat_data_single.aa_grid = empty_grid;
    }
  else
    {
      // For horizontally-aligned particles, the azimuth angle grid must cover
      // 0-180 degrees.
      if (scat_data_single.ptype == PTYPE_HORIZ_AL && data_aa_grid.nelem() == 0.)
        {
          ostringstream os;
          os << "For ptype = \"horizontally_aligned\""
          << " the azimuth angle grid can not be empty.";
          throw runtime_error( os.str() );
        }
      if (scat_data_single.ptype == PTYPE_HORIZ_AL && data_aa_grid[0] != 0.)
        {
          ostringstream os;
          os << "For ptype = \"horizontally_aligned\""
          << " the first value of the aa grid must be 0.";
          throw runtime_error( os.str() );
        }

      if (scat_data_single.ptype == PTYPE_HORIZ_AL && last(data_aa_grid) != 180.)
        {
          ostringstream os;
          os << "For ptype = \"horizontally_aligned\""
          << " the last value of the aa grid must be 180.";
          throw runtime_error( os.str() );
        }

      scat_data_single.aa_grid = data_aa_grid;
    }


  // Index coding for shape
  Index np;
  if( shape == "spheroidal" )
    { 
      np=-1; 
      if( aspect_ratio == 1 )
        throw runtime_error( "For spheroidal particles, the aspect ratio "
                             "is not allowed to be exactly 1 (due to "
                             "numerical problems)." );
    }
  else if( shape == "cylindrical" )
    { np=-2; }
  else
    {
      ostringstream os;
      os << "Unknown particle shape: " << shape << "\n"
         << "Must be \"spheroidal\" or \"cylindrical\".";
      throw runtime_error( os.str() );
    }

  // Interpolate refractive index to relevant grids
  //
  const Index nf = data_f_grid.nelem();
  const Index nt = data_t_grid.nelem();
  //
  Tensor3 ncomp( nf, nt, 2 );
  complex_n_interp( ncomp(joker,joker,0), ncomp(joker,joker,1),
                    complex_refr_index, "complex_refr_index", 
                    data_f_grid, data_t_grid );

  // Run T-matrix and we are ready (T-matrix takes size as equiv radius(!) )
  calcSingleScatteringDataProperties( scat_data_single,
                                      ncomp(joker,joker,0), 
                                      ncomp(joker,joker,1),
                                      0.5*diameter_volume_equ, 
                                      np, aspect_ratio, precision, ndgs,
                                      robust, quiet
                                    );


  // Meta data
  scat_meta_single.description  = 
                  "Meta data for associated file with single scattering data.";
  scat_meta_single.source       = 
                        "ARTS interface to T-matrix code by Mishchenko et al.";
  scat_meta_single.refr_index   = cri_source;
  //
  Numeric diameter_max, area_max;
  diameter_maxFromDiameter_volume_equ( diameter_max, area_max, shape, 
                                       diameter_volume_equ, aspect_ratio,
                                       verbosity );
  //
  scat_meta_single.mass                            = mass;
  scat_meta_single.diameter_max                    = diameter_max;
  scat_meta_single.diameter_volume_equ             = diameter_volume_equ;
  scat_meta_single.diameter_area_equ_aerodynamical = area_max;
}



void TMatrixTest(const Verbosity& verbosity)
{
    tmatrix_tmd_test(verbosity);
    tmatrix_ampld_test(verbosity);
    calc_ssp_random_test(verbosity);
    calc_ssp_fixed_test(verbosity);
}



//-----------------------------------

                 
/* commenting out. in ARTS 2.3 we change the interface-TMatrix approach:
   Instead of setting the meta data and doing TMatrix based on that, now we
   first do the TMatrix (per scattering element; use scat_data_singleTmatrix),
   then or along with this the meta data is dervied and set.
   Keeping old code, but commented out, for now in case we want to re-use
   something...

void scat_metaAddTmatrix(// WS Output:
                         ArrayOfArrayOfScatteringMetaData& scat_meta,
                         // WS Input:
                         const GriddedField3& complex_refr_index,
                         // WS Generic input:
                         const String& description,
                         const String& material,
                         const String& shape,
                         const String& ptype,
                         const Numeric& density,
                         const Vector& aspect_ratio_grid,
                         const Vector& diameter_max_grid,
                         const Vector& scat_f_grid,
                         const Vector& scat_T_grid,
                         const Verbosity&)

{
    // FIXME OLE: Where to add to scat_meta?
  Index last_species = scat_meta.nelem()-1;
  if (last_species < 0)
  {
      scat_meta.resize(1);
      last_species = 0;
  }

  for(Index k=0; k<diameter_max_grid.nelem(); k++)
    {
     for ( Index i=0; i < aspect_ratio_grid.nelem(); i++ )
      {
        Numeric volume;
        
        if (shape == "spheroidal")
        {
            if (aspect_ratio_grid[i]>1) // If oblate spheroid
            {   
                Numeric a; //rotational axis (perpendicular to a)
                
                a=diameter_max_grid[k]/2.;
                volume=pow(a,3.)*4.*PI/(3.*aspect_ratio_grid[i]);
            }
        
            else if (aspect_ratio_grid[i]<1) // If prolate spheroid
            {  
                Numeric b; //non-rotational axis (perpendicular to a)
                
                b=diameter_max_grid[k]/2.;
                volume=pow(b,3.)*4.*PI*pow(aspect_ratio_grid[i],2.)/3.;
            }
        
            else
            {
                ostringstream os;
                os << "Incorrect aspect ratio: " << aspect_ratio_grid[i] << "\n"
                << "Can not be equal to one";
                throw std::runtime_error(os.str());
            }
        }
        else if (shape == "cylindrical")
        {
            // Formulas to determine cylindrical volume from diamter_max:
            
            // D/L=aspect_ratio
            // diameter_max=sqrt(D²+L²)
            // volume=D³*pi/(4*aspect_ratio)
            
            volume=pow(diameter_max_grid[k]/pow((pow(aspect_ratio_grid[i], 2.)+1.), 1./2.), 3)*pow(aspect_ratio_grid[i], 2.)*PI/4.;
         }

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
            os << shape<< " "<< material << " particle of type " << ptype<<
            ", with volume equivalent diameter "
            <<diameter_max_grid[k]<<" meters.";
            smd.description=os.str();
        }
        else 
            smd.description = description;
        
        smd.material = material;
        smd.shape = shape;
        smd.ptype = PTypeFromString(ptype);
        smd.ssd_method = PARTICLE_SSDMETHOD_TMATRIX;
        smd.density = density;
        smd.diameter_max =diameter_max_grid[k];
        smd.volume = volume;
        smd.area_projected = 0;
        smd.aspect_ratio = aspect_ratio_grid[i];
        smd.scat_f_grid = scat_f_grid;
        smd.scat_T_grid = scat_T_grid;
        smd.complex_refr_index = complex_refr_index;

        scat_meta[last_species].push_back(smd);
      }
    }
}
*/



/* commenting out. in ARTS 2.3 we change the interface-TMatrix approach:
   Instead of setting the meta data and doing TMatrix based on that, now we
   first do the TMatrix (per scattering element; use scat_data_singleTmatrix),
   then or along with this the meta data is dervied and set.
   Keeping old code, but commented out, for now in case we want to re-use
   something...

void scat_dataFromMeta(// WS Output:
                       ArrayOfArrayOfSingleScatteringData& scat_data,
                       //WS Input
                       const ArrayOfArrayOfScatteringMetaData& scat_meta,
                       // WS Generic input:
                       const Vector& za_grid,
                       const Vector& aa_grid,
                       const Numeric& precision,
                       const Verbosity&)
{
    for(Index i_ss=0;i_ss<scat_meta.nelem();i_ss++)
    {
        ArrayOfSingleScatteringData assd;

        for(Index i_se=0;i_se<scat_meta.nelem();i_se++)
        {

            //Interpolation

            Tensor3 complex_refr_index;

            const Index nf = scat_meta[i_ss][i_se].scat_f_grid.nelem();
            const Index nt = scat_meta[i_ss][i_se].scat_T_grid.nelem();

            complex_refr_index.resize(nf,nt,2);

            complex_n_interp(complex_refr_index(joker, joker,0), complex_refr_index(joker, joker,1),
                             scat_meta[i_ss][i_se].complex_refr_index, "complex_refr_index",
                             scat_meta[i_ss][i_se].scat_f_grid, scat_meta[i_ss][i_se].scat_T_grid);


            Index  np;

            SingleScatteringData ssd;
            ssd.f_grid = scat_meta[i_ss][i_se].scat_f_grid;
            ssd.T_grid = scat_meta[i_ss][i_se].scat_T_grid;
            ssd.za_grid = za_grid;
            ssd.aa_grid = aa_grid;
            ssd.ptype = scat_meta[i_ss][i_se].ptype;

            if (scat_meta[i_ss][i_se].shape == "spheroidal" )
                np=-1;

            else if (scat_meta[i_ss][i_se].shape == "cylindrical")
                np=-2;
            else
            {
                ostringstream os;
                os << "Unknown particle shape: " << scat_meta[i_ss][i_se].shape << "\n"
                << "Must be spheroidal or cylindrical";
                throw std::runtime_error(os.str());
            }

            calcSingleScatteringDataProperties(ssd,
                                               complex_refr_index(joker,joker,0),
                                               complex_refr_index(joker,joker,1),
                                               pow(scat_meta[i_ss][i_se].volume*1e18*3./(4.*PI),1./3.),
                                               np,
                                               scat_meta[i_ss][i_se].aspect_ratio,
                                               precision);
            
            assd.push_back(ssd);
        }

        scat_data.push_back(assd);
    }
}

*/



/* commenting out. in ARTS 2.3 we change the interface-TMatrix approach:
   Instead of setting the meta data and doing TMatrix based on that, now we
   first do the TMatrix (per scattering element; use scat_data_singleTmatrix),
   then or along with this the meta data is dervied and set.
   Keeping old code, but commented out, for now in case we want to re-use
   something...

void scat_metaAddTmatrixOldVersion(// WS Output:
                             ArrayOfScatteringMetaData& scat_meta,
                             // WS Input:
                             const GriddedField3& complex_refr_index,
                             // WS Generic input:
                             const String& description,
                             const String& material,
                             const String& shape,
                             const String& ptype,
                             const Numeric& density,
                             const Numeric& aspect_ratio,
                             const Vector& diameter_grid,
                             const Vector& scat_f_grid,
                             const Vector& scat_T_grid,
                             const Verbosity&) 

{
  for(Index k=0; k<diameter_grid.nelem(); k++)
    {

      Numeric volume;
      
      volume=4./3.*PI*pow(diameter_grid[k]/2.,3);
      
      
      Numeric diameter_max;

      if (shape == "spheroidal")
        {
       if (aspect_ratio>1) // If oblate spheroid
            {   
            Numeric a; //rotational axis
            
            a=pow(3*volume*aspect_ratio/(4*PI), 1./3.);
            diameter_max=2*a;
            }
        
       else if (aspect_ratio<1) // If prolate spheroid
            {  
            Numeric b; //non-rotational axis (perpendicular to a)
            
            b=pow(3*volume/(4*PI*pow(aspect_ratio, 2)),1./3.);
            diameter_max=2*b;
            }
        
        else
            {
            ostringstream os;
            os << "Incorrect aspect ratio: " << aspect_ratio << "\n"
              << "Can not be equal to one";
            throw std::runtime_error(os.str());
            }
        }
      else if (shape == "cylindrical")
        {
            Numeric D;
            Numeric L;
            
            //aspect_ratio=D/L
            
            D=pow(volume*4./PI*aspect_ratio, 1./3.);
            L=D/aspect_ratio;
            diameter_max=pow((pow(D,2)+pow(L,2)), 1./2.);
         }

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
          os << shape<< " "<< material << " particle of type " << ptype<<
            ", with volume equivalent diameter "
            <<diameter_grid[k]<<" meters.";
          smd.description=os.str();
        }
      else 
        smd.description = description;
 
      smd.material = material;
      smd.shape = shape;
      smd.ptype = PTypeFromString(ptype);
      smd.ssd_method = PARTICLE_SSDMETHOD_TMATRIX;
      smd.density = density;
      smd.diameter_max =diameter_max;
      smd.volume = volume;
      smd.area_projected = 0;
      smd.aspect_ratio = aspect_ratio;
      smd.scat_f_grid = scat_f_grid;
      smd.scat_T_grid = scat_T_grid;
      smd.complex_refr_index = complex_refr_index;

      scat_meta.push_back(smd);
    }
}

*/
