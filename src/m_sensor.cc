/* Copyright (C) 2003-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler(at)ltu.se>

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



/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   m_sensor.cc
  \author Mattias Ekström/Patrick Eriksson <patrick.eriksson@chalmers.se>
  \date   2003-02-13

  \brief  Workspace functions related to sensor modelling variables.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include <string>
#include "arts.h"
#include "check_input.h"
#include "interpolation_poly.h"
#include "math_funcs.h"
#include "messages.h"
#include "ppath.h"
#include "special_interp.h"
#include "xml_io.h"
#include "sensor.h"
#include "sorting.h"
#include "auto_md.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Index GFIELD1_F_GRID;
extern const Index GFIELD4_FIELD_NAMES;
extern const Index GFIELD4_F_GRID;
extern const Index GFIELD4_ZA_GRID;
extern const Index GFIELD4_AA_GRID;
extern const Numeric SPEED_OF_LIGHT;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void AntennaConstantGaussian1D(Index&           antenna_dim,
                               Vector&          mblock_za_grid,
                               Vector&          mblock_aa_grid,
                               GriddedField4&   r,
                               Matrix&          antenna_los,
                               const Index&     n_za_grid,
                               const Numeric&   fwhm,
                               const Numeric&   xwidth_si,
                               const Numeric&   dx_si,
                               const Verbosity& verbosity)
{
  if( dx_si > xwidth_si )
    throw runtime_error( "It is demanded that dx_si <= xwidth_si." );

  AntennaSet1D( antenna_dim, mblock_aa_grid, verbosity );
  antenna_responseGaussian( r, fwhm, xwidth_si, dx_si, verbosity );

  antenna_los.resize(1,1);
  antenna_los(0,0) = 0.0;

  // za grid for response
  ConstVectorView r_za_grid =  r.get_numeric_grid(GFIELD4_ZA_GRID);
  const Index nr = r_za_grid.nelem();

  // Cumulative integral of response (factor /2 skipped, but does not matter)
  Vector cumtrapz(nr);
  cumtrapz[0] = 0;
  for( Index i=1; i<nr; i++ )  
    { cumtrapz[i] = cumtrapz[i-1] + r.data(0,0,i-1,0) + r.data(0,0,i,0); }

  // Equally spaced vector between end points of cumulative sum
  Vector csp;
  nlinspace( csp, cumtrapz[0], cumtrapz[nr-1], n_za_grid );

  // Get mblock_za_grid by interpolation
  mblock_za_grid.resize(n_za_grid);   
  ArrayOfGridPos gp(n_za_grid);
  gridpos( gp, cumtrapz, csp );
  Matrix itw(n_za_grid,2);
  interpweights( itw, gp );
  interp( mblock_za_grid, itw, r_za_grid, gp );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void AntennaMultiBeamsToPencilBeams(Matrix&          sensor_pos,
                                    Matrix&          sensor_los,
                                    Matrix&          antenna_los,
                                    Index&           antenna_dim,
                                    Vector&          mblock_za_grid,
                                    Vector&          mblock_aa_grid,
                                    const Index&     atmosphere_dim,
                                    const Verbosity& verbosity)
{
  // Sizes
  const Index nmblock = sensor_pos.nrows();
  const Index nant    = antenna_los.nrows();

  // Input checks
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_if_in_range( "antenna_dim", antenna_dim, 1, 2 );
  //
  if( sensor_pos.ncols() != atmosphere_dim )
    throw runtime_error( "The number of columns of sensor_pos must be "
                              "equal to the atmospheric dimensionality." );
  if( atmosphere_dim <= 2  &&  sensor_los.ncols() != 1 )
    throw runtime_error( 
                      "For 1D and 2D, sensor_los shall have one column." );
  if( atmosphere_dim == 3  &&  sensor_los.ncols() != 2 )
    throw runtime_error( "For 3D, sensor_los shall have two columns." );
  if( sensor_los.nrows() != nmblock )
    {
      ostringstream os;
      os << "The number of rows of sensor_pos and sensor_los must be "
         << "identical, but sensor_pos has " << nmblock << " rows,\n"
         << "while sensor_los has " << sensor_los.nrows() << " rows.";
      throw runtime_error( os.str() );
    }
  if( antenna_dim == 2  &&  atmosphere_dim < 3 )
    {
      throw runtime_error(
                          "If *antenna_dim* is 2, *atmosphere_dim* must be 3." );
    }
  if( antenna_dim != antenna_los.ncols() ) 
    {
      throw runtime_error(
             "The number of columns of *antenna_los* must be *antenna_dim*.\n" );
    }

  // Claculate new sensor_pos and sensor_los
  const Matrix pos_copy = sensor_pos;
  const Matrix los_copy = sensor_los;
  //
  sensor_pos.resize( nmblock*nant, pos_copy.ncols() );
  sensor_los.resize( nmblock*nant, los_copy.ncols() );
  //
  for( Index ib=0; ib<nmblock; ib++ )
    {
      for( Index ia=0; ia<nant; ia++ )
        {
          const Index i = ib*nant + ia;

          sensor_pos(i,joker) = pos_copy(ib,joker);
          sensor_los(i,joker) = los_copy(ib,joker);
          
          sensor_los(i,0) += antenna_los(ia,0);
          if( antenna_dim == 2 )
            sensor_los(i,1) += antenna_los(ia,1);
        }
    }

  // Set other variables
  AntennaOff( antenna_dim, mblock_za_grid, mblock_aa_grid, verbosity );
  //
  antenna_los.resize( 1, 1 );
  antenna_los = 0;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void AntennaOff(// WS Output:
                Index&  antenna_dim,
                Vector& mblock_za_grid,
                Vector& mblock_aa_grid,
                const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  out2 << "  Sets the antenna dimensionality to 1.\n";
  antenna_dim = 1;

  out2 << "  Sets *mblock_za_grid* to have length 1 with value 0.\n";
  mblock_za_grid.resize(1);
  mblock_za_grid[0] = 0;

  out2 << "  Sets *mblock_aa_grid* to be an empty vector.\n";
  mblock_aa_grid.resize(0);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void AntennaSet1D(// WS Output:
                  Index&   antenna_dim,
                  Vector&  mblock_aa_grid,
                  const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  out2 << "  Sets the antenna dimensionality to 1.\n";
  out3 << "    antenna_dim = 1\n";
  out3 << "    mblock_aa_grid is set to be an empty vector\n";
  antenna_dim = 1;
  mblock_aa_grid.resize(0);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void AntennaSet2D(// WS Output:
                  Index&   antenna_dim,
                  // WS Input:
                  const Index&   atmosphere_dim,
                  const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  if( atmosphere_dim != 3 )
    throw runtime_error("Antenna dimensionality 2 is only allowed when the "
                                          "atmospheric dimensionality is 3." );
  out2 << "  Sets the antenna dimensionality to 1.\n";
  out3 << "    antenna_dim = 2\n";
  antenna_dim = 2;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void antenna_responseGaussian(GriddedField4&   r,
                              const Numeric&   fwhm,
                              const Numeric&   xwidth_si,
                              const Numeric&   dx_si,
                              const Verbosity&)
{
  if( dx_si > xwidth_si )
    throw runtime_error( "It is demanded that dx_si <= xwidth_si." );

  Vector x, y;
  gaussian_response_autogrid( x, y, 0, fwhm, xwidth_si, dx_si );

  r.set_name( "Antenna response" );

  r.set_grid_name( 0, "Polarisation" );
  r.set_grid( 0, MakeArray<String>( "NaN" ) ); 

  r.set_grid_name( 1, "Frequency" );
  r.set_grid( 1, Vector(1,-999) );

  r.set_grid_name( 2, "Zenith angle" );
  r.set_grid( 2, x );

  r.set_grid_name( 3, "Azimuth angle" );
  r.set_grid( 3, Vector(1,0) );

  const Index n = y.nelem();
  r.data.resize( 1, 1, n, 1 );
  r.data(0,0,joker,0) = y;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void antenna_responseVaryingGaussian(
         GriddedField4&   r,
   const Numeric&         leff,
   const Numeric&         xwidth_si,
   const Numeric&         dx_si,
   const Index&           nf,
   const Numeric&         fstart,
   const Numeric&         fstop,
   const Verbosity&       verbosity )
{
  r.set_name( "Antenna response" );

  r.set_grid_name( 0, "Polarisation" );
  r.set_grid( 0, MakeArray<String>( "NaN" ) ); 

  r.set_grid_name( 3, "Azimuth angle" );
  r.set_grid( 3, Vector(1,0) );

  Vector f_grid;
  VectorNLogSpace( f_grid, nf, fstart, fstop, verbosity );
  r.set_grid_name( 1, "Frequency" );
  r.set_grid( 1, f_grid );

  // Calculate response for highest frequency, with xwidth_si scaled from
  // fstart
  Vector x, y;
  Numeric fwhm = RAD2DEG * SPEED_OF_LIGHT / ( leff * fstop );
  gaussian_response_autogrid( x, y, 0, fwhm, (fstop/fstart)*xwidth_si, dx_si );

  r.set_grid_name( 2, "Zenith angle" );
  r.set_grid( 2, x );

  const Index n = y.nelem();
  r.data.resize( 1, nf, n, 1 );
  //
  r.data(0,nf-1,joker,0) = y;
  //
  for( Index i=0; i<nf-1; i++ )
    {
      fwhm = RAD2DEG * SPEED_OF_LIGHT / ( leff * f_grid[i] );
      gaussian_response( y, x, 0, fwhm );
      r.data(0,i,joker,0) = y;
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void backend_channel_responseFlat(ArrayOfGriddedField1& r,
                                  const Numeric&        resolution,
                                  const Verbosity&)
{
  r.resize( 1 );
  r[0].set_name( "Backend channel response function" );

  Vector x(2);

  r[0].set_grid_name( 0, "Frequency" );
  x[1] = resolution / 2.0;
  x[0] = -x[1];
  r[0].set_grid( 0, x );

  r[0].data.resize( 2 );
  r[0].data[0] = 1/ resolution;
  r[0].data[1] = r[0].data[0];
}



/* Workspace method: Doxygen documentation will be auto-generated */
void backend_channel_responseGaussian(ArrayOfGriddedField1&   r,
                                      const Numeric&   fwhm,
                                      const Numeric&   xwidth_si,
                                      const Numeric&   dx_si,
                                      const Verbosity&)
{
  r.resize( 1 );
  Vector x, y;

  gaussian_response_autogrid( x, y, 0, fwhm, xwidth_si, dx_si );

  r[0].set_name( "Backend channel response function" );

  r[0].set_grid_name( 0, "Frequency" );
  r[0].set_grid( 0, x );

  const Index n = y.nelem();
  r[0].data.resize( n );
  for( Index i=0; i<n; i++ )
    r[0].data[i] = y[i];
}



/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridFromSensorAMSU(// WS Output:
                          Vector& f_grid,
                          // WS Input:
                          const Vector& lo,
                          const ArrayOfVector& f_backend,
                          const ArrayOfArrayOfGriddedField1& backend_channel_response,
                          // Control Parameters:
                          const Numeric& spacing,
                          const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  // Find out how many channels we have in total:
  // f_backend is an array of vectors, containing the band frequencies for each Mixer.
  Index n_chan = 0;
  for (Index i=0; i<f_backend.nelem(); ++i)
    for (Index s=0; s<f_backend[i].nelem(); ++s)
      ++n_chan;

  // Checks on input quantities:

  // There must be at least one channel:
  if (n_chan < 1)
    {
      ostringstream os;
      os << "There must be at least one channel.\n"
         << "(The vector *lo* must have at least one element.)";
      throw runtime_error(os.str());
    }

  // Is number of LOs consistent in all input variables?
  if ( (f_backend.nelem()                != lo.nelem()) ||
       (backend_channel_response.nelem() != lo.nelem()) )
    {
      ostringstream os;
      os << "Variables *lo_multi*, *f_backend_multi* and *backend_channel_response_multi*\n"
         << "must have same number of elements (number of LOs).";
      throw runtime_error(os.str());
    }

  // Is number of bands consistent for each LO?
  for (Index i=0; i<f_backend.nelem(); ++i)
    if (f_backend[i].nelem() != backend_channel_response[i].nelem())
    {
      ostringstream os;
      os << "Variables *f_backend_multi* and *backend_channel_response_multi*\n"
         << "must have same number of bands for each LO.";
      throw runtime_error(os.str());
    }

  // Start the actual work.
  
  // We construct the necessary input for function find_effective_channel_boundaries, 
  // which will identify channel boundaries, taking care of overlaping channels.

  // A "flat" vector of nominal band frequencies (two for each AMSU channel):
  Vector f_backend_flat(2*n_chan);

  // A "flat" list of channel response functions (two for each AMSU channel)
  ArrayOfGriddedField1 backend_channel_response_flat(2*n_chan);
  
  // Counts position inside the flat arrays during construction:
  Index j=0;
  
  for (Index i=0; i<f_backend.nelem(); ++i)
    for (Index s=0; s<f_backend[i].nelem(); ++s)
    {      
      const GriddedField1& this_grid = backend_channel_response[i][s];
      const Numeric this_f_backend = f_backend[i][s];

      // Signal sideband:
      f_backend_flat[j] = this_f_backend;
      backend_channel_response_flat[j] = this_grid;
      j++;
      
      // Image sideband:
      Numeric offset  = this_f_backend - lo[i];
      Numeric f_image = lo[i] - offset;
      f_backend_flat[j] = f_image;
      backend_channel_response_flat[j] = this_grid;      
      j++;
    }  
  
  // We build up a total list of absolute frequency ranges for
  // both signal and image sidebands:
  Vector fmin(2*n_chan), fmax(2*n_chan);

  // We have to add some additional margin at the band edges, 
  // otherwise the instrument functions are not happy. Define 
  // this in terms of the grid spacing:
  Numeric delta = 1*spacing;
  
  // Call subfunction to do the actual work of merging overlapping 
  // channels and identifying channel boundaries:
  find_effective_channel_boundaries(fmin,
                                    fmax,
                                    f_backend_flat,
                                    backend_channel_response_flat,
                                    delta, verbosity);
  
  // Create f_grid_array. This is an array of Numeric, so that we
  // can use the STL push_back function.
  ArrayOfNumeric f_grid_array;

  for (Index i=0; i<fmin.nelem(); ++i)
    {
      // Band width:
      const Numeric bw = fmax[i] - fmin[i];

      // How many grid intervals do I need?
      const Numeric npf = ceil(bw/spacing);

      // How many grid points to store? - Number of grid intervals
      // plus 1.
      const Index   npi = (Index) npf + 1;

      // Create the grid for this band:
      Vector grid;
      nlinspace(grid,fmin[i],fmax[i],npi);
      

      out3 << "  Band range " << i << ": " << grid << "\n";

      // Append to f_grid_array:
      f_grid_array.reserve(f_grid_array.nelem()+npi);
      for (Index s=0; s<npi; ++s)
        f_grid_array.push_back(grid[s]);
    }

  // Copy result to output vector:
  f_grid = f_grid_array;
  
  out2 << "  Total number of frequencies in f_grid: " << f_grid.nelem() << "\n";
  //  cout << "f_grid: " << f_grid << "\n";
}

void f_gridFromSensorAMSUgeneric(// WS Output:
    Vector& f_grid,
    // WS Input:
    const ArrayOfVector& f_backend_multi,  // Center frequency for each channel 
    const ArrayOfArrayOfGriddedField1& backend_channel_response_multi,
    // Control Parameters:
    const Numeric& spacing,
    const Vector& verbosityVect,
    const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  Index numFpoints = 0;
  // Calculate the total number of frequency points
  for(Index idx = 0;idx<backend_channel_response_multi.nelem();++idx)
  {
    for(Index idy = 0;idy<backend_channel_response_multi[idx].nelem();++idy)
    {
      numFpoints += backend_channel_response_multi[idx][idy].get_grid_size(0);
    }
  }

  Index numChan = backend_channel_response_multi.nelem(); // number of channels

  // Checks on input quantities:

  // There must be at least one channel:
  if (numChan < 1)
  {
    ostringstream os;
    os << "There must be at least one channel.\n"
      << "(The vector *lo* must have at least one element.)";
    throw runtime_error(os.str());
  }


  // Start the actual work.
  // We construct the necessary input for function find_effective_channel_boundaries, 
  // which will identify channel boundaries, taking care of overlaping channels.


  // A "flat" vector of nominal band frequencies (one for each AMSU channel):
  Vector f_backend_flat(numChan);
  // A "flat" list of channel response functions (one for each AMSU channel)
  ArrayOfGriddedField1 backend_channel_response_nonflat(numChan);

  // Save the correct verbosity value to each passband 
  Vector FminVerbosityVect(numFpoints);
  Vector FmaxVerbosityVect(numFpoints);
  Vector VerbosityValVect(numFpoints);
  Index VerbVectIdx =0; 


  for(Index i =0;i<f_backend_multi.nelem();++i)
  {
    for(Index ii =0;ii<f_backend_multi[i].nelem();++ii)
    {
      const GriddedField1& this_grid = backend_channel_response_multi[i][ii];
      const Numeric this_f_backend = f_backend_multi[i][ii];
      // Signal passband :
      f_backend_flat[i] = this_f_backend;
      backend_channel_response_nonflat[i] = this_grid;
      for(Index idy=1; idy < backend_channel_response_multi[i][ii].get_grid_size(0);++idy)
      {
        if((backend_channel_response_multi[i][ii].data[idy-1]==0) && (backend_channel_response_multi[i][ii].data[idy]>0))
        {
          FminVerbosityVect[VerbVectIdx] = f_backend_multi[i][ii]+backend_channel_response_multi[i][ii].get_numeric_grid(0)[idy];
          VerbosityValVect[VerbVectIdx] = verbosityVect[i];
        }
        if((backend_channel_response_multi[i][ii].data[idy-1]>0) && (backend_channel_response_multi[i][ii].data[idy]==0))
        {
          FmaxVerbosityVect[VerbVectIdx] = f_backend_multi[i][ii]+backend_channel_response_multi[i][ii].get_numeric_grid(0)[idy];
          VerbVectIdx ++;
        }
      }
    }
  }
  // We build up a total list of absolute frequency ranges for all passbands 
  // Code reused from the function "f_gridFromSensorAMSU"
  Vector fmin,fmax;  // - these variables will be resized, therefore len(1) is enough for now.,

  // We have to add some additional margin at the band edges, 
  // otherwise the instrument functions are not happy. Define 
  // this in terms of the grid spacing:
  const Numeric delta = 10;

  // Call subfunction to do the actual work of merging overlapping 
  // channels and identifying channel boundaries:
  find_effective_channel_boundaries(fmin,
      fmax,
      f_backend_flat,
      backend_channel_response_nonflat,
      delta, verbosity);

  // Create f_grid_array. This is an array of Numeric, so that we
  // can use the STL push_back function.
  ArrayOfNumeric f_grid_array;

  for (Index i=0; i<fmin.nelem(); ++i)
  {
    // Bandwidth:
    const Numeric bw = fmax[i] - fmin[i];
    Numeric npf  = ceil(bw/spacing); // Set a default value

    // How many grid intervals do I need?
    Index verbIdx = 0;
    if(verbosityVect.nelem()>0)
    { 
      // find the grid needed for the particular part of passband
      for(Index ii =0;ii<VerbVectIdx;++ii)
      {
        if((FminVerbosityVect[ii]>=fmin[i]) &&(FmaxVerbosityVect[ii]<=fmax[i]))        
        {
          if(verbIdx ==0)
          {
            verbIdx = ii;
          } 
          else 
          {
            if(VerbosityValVect[ii]<VerbosityValVect[verbIdx])
            {
              verbIdx = ii;
            }
          }
        }
      }
      if(spacing > VerbosityValVect[verbIdx])
      {
        npf = ceil(bw/VerbosityValVect[verbIdx]); // is the default value to coarse?
      }
      else 
      {
        npf  = ceil(bw/spacing); // Default value
      }
    }


    // How many grid points to store? - Number of grid intervals
    // plus 1.
    const Index   npi = (Index) npf + 1;

    // What is the actual grid spacing inside the band?
    const Numeric gs = bw/npf;

    // Create the grid for this band:
    Vector grid(fmin[i], npi, gs);

    out3 << "  Band range " << i << ": " << grid << "\n";

    // Append to f_grid_array:
    f_grid_array.reserve(f_grid_array.nelem()+npi);
    for (Index s=0; s<npi; ++s)
      f_grid_array.push_back(grid[s]);
  }

  // Copy result to output vector:
  f_grid = f_grid_array;

  out2 << "  Total number of frequencies in f_grid: " << f_grid.nelem() << "\n";

}





/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridFromSensorHIRS(// WS Output:
                          Vector& f_grid,
                          // WS Input:
                          const Vector& f_backend,
                          const ArrayOfGriddedField1& backend_channel_response,
                          // Control Parameters:
                          const Numeric& spacing,
                          const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  // Check input
  if (spacing <= 0) {
    ostringstream os;
    os << "Expected positive spacing. Found spacing to be: "
       << spacing << "\n";
    throw runtime_error(os.str());
  }
  // Call subfunction to get channel boundaries. Also does input
  // consistency checking for us.
  Vector fmin, fmax;

  // We have to add some additional margin at the band edges, 
  // otherwise the instrument functions are not happy. Define 
  // this in terms of the grid spacing:
  Numeric delta = 1*spacing;
  
  find_effective_channel_boundaries(fmin,
                                    fmax,
                                    f_backend,
                                    backend_channel_response,
                                    delta, verbosity);

  // Ok, now we just have to create a frequency grid for each of the
  // fmin/fmax ranges.

  // Create f_grid_array. This is an array of Numeric, so that we
  // can use the STL push_back function.
  ArrayOfNumeric f_grid_array;

  for (Index i=0; i<fmin.nelem(); ++i)
    {
      // Band width:
      const Numeric bw = fmax[i] - fmin[i];

      // How many grid intervals do I need?
      const Numeric npf = ceil(bw/spacing);

      // How many grid points to store? - Number of grid intervals
      // plus 1.
      const Index   npi = (Index) npf + 1;

      // Create the grid for this band:
      Vector grid;
      nlinspace(grid,fmin[i],fmax[i],npi);
     
      out3 << "  Band range " << i << ": " << grid << "\n";

      // Append to f_grid_array:
      f_grid_array.reserve(f_grid_array.nelem()+npi);
      for (Index s=0; s<npi; ++s)
        f_grid_array.push_back(grid[s]);
    }

  // Copy result to output vector:
  f_grid = f_grid_array;

  out2 << "  Total number of frequencies in f_grid: " << f_grid.nelem() << "\n";

}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseAntenna(// WS Output:
                            Sparse&       sensor_response,
                            Vector&       sensor_response_f,
                            ArrayOfIndex& sensor_response_pol,
                            Vector&       sensor_response_za,
                            Vector&       sensor_response_aa,
                            Vector&       sensor_response_za_grid,
                            Vector&       sensor_response_aa_grid,
                            // WS Input:
                            const Vector&         sensor_response_f_grid,
                            const ArrayOfIndex&   sensor_response_pol_grid,
                            const Index&          atmosphere_dim,
                            const Index&          antenna_dim,
                            const Matrix&         antenna_los,
                            const GriddedField4&  antenna_response,
                            const Index&   sensor_norm,
                            const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Basic checks
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_if_in_range( "antenna_dim",    antenna_dim,    1, 2 );
  chk_if_bool(     "sensor_norm",    sensor_norm          );

  // Some sizes
  const Index nf   = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nza  = sensor_response_za_grid.nelem();
  const Index naa  = max( Index(1), sensor_response_aa_grid.nelem() );
  const Index nin  = nf * npol * nza * naa;

  // Initialise a output stream for runtime errors and a flag for errors
  ostringstream os;
  bool          error_found = false;


  // Check that sensor_response variables are consistent in size
  if( sensor_response_f.nelem() != nin )
  {
    os << "Inconsistency in size between *sensor_response_f* and the sensor\n"
       << "grid variables (sensor_response_f_grid etc.).\n";
    error_found = true;
  }
  if( sensor_response.nrows() != nin )
  {
    os << "The sensor block response matrix *sensor_response* does not have\n"
       << "right size compared to the sensor grid variables\n"
       << "(sensor_response_f_grid etc.).\n";
    error_found = true;
  }


  // Checks related to antenna dimension
  if( antenna_dim == 2  &&  atmosphere_dim < 3 )
  {
    os << "If *antenna_dim* is 2, *atmosphere_dim* must be 3.\n";
    error_found = true;
  }
  if( antenna_dim == 1  &&  sensor_response_aa_grid.nelem() )
  {
    os << "If *antenna_dim* is 1, *sensor_response_aa_grid* (and\n"
       << "*mblock_aa_grid*) must be empty.";
    error_found = true;
  }


  // Check of antenna_los
  if( antenna_dim != antenna_los.ncols() ) 
  {
    os << "The number of columns of *antenna_los* must be *antenna_dim*.\n";
    error_found = true;
  }
  // We allow angles in antenna_los to be unsorted


  // Checks of antenna_response polarisation dimension
  //
  const Index lpolgrid = 
                  antenna_response.get_string_grid(GFIELD4_FIELD_NAMES).nelem();
  //
  if( lpolgrid != 1  &&  lpolgrid != npol ) 
  {
    os << "The number of polarisation in *antenna_response* must be 1 or be\n"
       << "equal to the number of polarisations used (determined by\n"
       << "*stokes_dim* or *sensor_pol*).\n";
    error_found = true;
  }


  // Checks of antenna_response frequency dimension
  //
  ConstVectorView aresponse_f_grid = 
                              antenna_response.get_numeric_grid(GFIELD4_F_GRID);
  //
  chk_if_increasing( "f_grid of antenna_response", aresponse_f_grid );
  //
  Numeric f_dlow  = 0.0;
  Numeric f_dhigh = 0.0;
  //
  f_dlow  = min(sensor_response_f_grid) - aresponse_f_grid[0];
  f_dhigh = last(aresponse_f_grid) - max(sensor_response_f_grid);
  //
  if( aresponse_f_grid.nelem() > 1 )
  {
    if( f_dlow < 0 ) 
    {
      os << "The frequency grid of *antenna_response is too narrow. It must\n"
         << "cover all considered frequencies (*f_grid*), if the length\n"
         << "is > 1. The grid needs to be expanded with "<<-f_dlow<<" Hz in\n"
         << "the lower end.\n";
      error_found = true;
    }
    if( f_dhigh < 0 ) 
    {
      os << "The frequency grid of *antenna_response is too narrow. It must\n"
         << "cover all considered frequencies (*f_grid*), if the length\n"
         << "is > 1. The grid needs to be expanded with "<<-f_dhigh<<" Hz in\n"
         << "the upper end.\n";
      error_found = true;
    }
  }


  // Checks of antenna_response za dimension
  //
  ConstVectorView aresponse_za_grid = 
                             antenna_response.get_numeric_grid(GFIELD4_ZA_GRID);
  //
  chk_if_increasing( "za_grid of *antenna_response*", aresponse_za_grid );
  //
  if( aresponse_za_grid.nelem() < 2 )
  {
    os << "The zenith angle grid of *antenna_response* must have >= 2 values.\n";
    error_found = true;
    
  }
  //
  // Check if the relative grid added to the antenna_los za angles
  // outside sensor_response_za_grid.
  //
  Numeric za_dlow  = 0.0;
  Numeric za_dhigh = 0.0;
  //
  za_dlow = min(antenna_los(joker,0)) + aresponse_za_grid[0] -
                                                   min(sensor_response_za_grid);
  za_dhigh = max(sensor_response_za_grid) - ( max(antenna_los(joker,0)) +
                                                      last(aresponse_za_grid) );
  //
  if( za_dlow < 0 ) 
  {
    os << "The WSV *sensor_response_za_grid* is too narrow. It should be\n"
       << "expanded with "<<-za_dlow<<" deg in the lower end. This change\n"
       << "should be probably applied to *mblock_za_grid*.\n";
    error_found = true;
  }
  if( za_dhigh < 0 ) 
  {
    os << "The WSV *sensor_response_za_grid* is too narrow. It should be\n"
       << "expanded with "<<-za_dhigh<<" deg in the higher end. This change\n"
       << "should be probably applied to *mblock_za_grid*.\n";
    error_found = true;
  }


  // Checks of antenna_response aa dimension
  //
  ConstVectorView aresponse_aa_grid = 
                             antenna_response.get_numeric_grid(GFIELD4_AA_GRID);
  //
  if( antenna_dim == 1 )
  {
    if( aresponse_aa_grid.nelem() != 1 )
    {
      os << "The azimuthal dimension of *antenna_response* must be 1 if\n"
         << "*antenna_dim* equals 1.\n";
      error_found = true;    
    }
  }
  else
  {
    chk_if_increasing( "aa_grid of antenna_response", aresponse_aa_grid );
    //
    if( aresponse_za_grid.nelem() < 2 )
      {
        os << "The zenith angle grid of *antenna_response* must have >= 2\n"
           << "values.\n";
        error_found = true;
      }
    // Check if the relative grid added to the antena_los aa angles
    // outside sensor_response_aa_grid.
    //
    Numeric aa_dlow  = 0.0;
    Numeric aa_dhigh = 0.0;
    //
    aa_dlow = min(antenna_los(joker,1)) + aresponse_aa_grid[0] -
                                                   min(sensor_response_aa_grid);
    aa_dhigh = max(sensor_response_aa_grid) - ( max(antenna_los(joker,1)) +
                                                      last(aresponse_aa_grid) );
    //
    if( aa_dlow < 0 ) 
    {
      os << "The WSV *sensor_response_aa_grid* is too narrow. It should be\n"
         << "expanded with "<<-aa_dlow<<" deg in the lower end. This change\n"
         << "should be probably applied to *mblock_aa_grid*.\n";
      error_found = true;
    }
    if( aa_dhigh < 0 ) 
    {
      os << "The WSV *sensor_response_aa_grid* is too narrow. It should be\n"
         << "expanded with "<<-aa_dhigh<<" deg in the higher end. This change\n"
         << "should be probably applied to *mblock_aa_grid*.\n";
      error_found = true;
    }
  }


  // If errors where found throw runtime_error with the collected error
  // message.
  if (error_found)
    throw runtime_error(os.str());

  // And finally check if grids and data size match
  antenna_response.checksize_strict();


  // Call the core function 
  //
  Sparse hantenna;
  //
  if( antenna_dim == 1 )
    antenna1d_matrix( hantenna, antenna_dim, antenna_los, antenna_response, 
                      sensor_response_za_grid, sensor_response_f_grid, 
                      npol, sensor_norm );
  else
    antenna2d_simplified( hantenna, antenna_dim, antenna_los, antenna_response,
                          sensor_response_za_grid, sensor_response_aa_grid, 
                          sensor_response_f_grid, npol, sensor_norm ); 

  // Here we need a temporary sparse that is copy of the sensor_response
  // sparse matrix. We need it since the multiplication function can not
  // take the same object as both input and output.
  Sparse htmp = sensor_response;
  sensor_response.resize( hantenna.nrows(), htmp.ncols());
  mult( sensor_response, hantenna, htmp );

  // Some extra output.
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  // Update sensor_response_za_grid
  sensor_response_za_grid = antenna_los(joker,0);

  // Update sensor_response_aa_grid
  if( antenna_dim == 2 )
    sensor_response_aa_grid = antenna_los(joker,1);

  // Set aux variables
  sensor_aux_vectors( sensor_response_f,       sensor_response_pol, 
                      sensor_response_za,      sensor_response_aa, 
                      sensor_response_f_grid,  sensor_response_pol_grid, 
                      sensor_response_za_grid, sensor_response_aa_grid, 0 );
}





/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseBackend(// WS Output:
                            Sparse&   sensor_response,
                            Vector&   sensor_response_f,
                            ArrayOfIndex&   sensor_response_pol,
                            Vector&   sensor_response_za,
                            Vector&   sensor_response_aa,
                            Vector&   sensor_response_f_grid,
                            // WS Input:
                            const ArrayOfIndex&   sensor_response_pol_grid,
                            const Vector&   sensor_response_za_grid,
                            const Vector&   sensor_response_aa_grid,
                            const Vector&   f_backend,
                            const ArrayOfGriddedField1&   backend_channel_response,
                            const Index&   sensor_norm,
                            const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Some sizes
  const Index nf   = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nza  = sensor_response_za_grid.nelem();
  const Index naa  = sensor_response_aa_grid.nelem();
  const Index nin  = nf * npol * nza;
  // Note that there is no distinction between za and aa grids after the antenna

  // Initialise an output stream for runtime errors and a flag for errors
  ostringstream os;
  bool          error_found = false;

  // Check that sensor_response variables are consistent in size
  if( sensor_response_f.nelem() != nin )
  {
    os << "Inconsistency in size between *sensor_response_f* and the sensor\n"
       << "grid variables (sensor_response_f_grid etc.).\n";
    error_found = true;
  }
  if( naa  &&  naa != nza )
  {
    os << "Incorrect size of *sensor_response_aa_grid*.\n";
    error_found = true;
  }
  if( sensor_response.nrows() != nin )
  {
    os << "The sensor block response matrix *sensor_response* does not have\n"
       << "right size compared to the sensor grid variables\n"
       << "(sensor_response_f_grid etc.).\n";
    error_found = true;
  }

  // We allow f_backend to be unsorted, but must be inside sensor_response_f_grid
  if( min(f_backend) < min(sensor_response_f_grid) )
    {
      os << "At least one value in *f_backend* (" << min(f_backend) 
         << ") below range\ncovered by *sensor_response_f_grid* ("
         << min(sensor_response_f_grid) << ").\n";
      error_found = true;
    }
  if( max(f_backend) > max(sensor_response_f_grid) )
    {
      os << "At least one value in *f_backend* (" << max(f_backend) 
         << ") above range\ncovered by *sensor_response_f_grid* ("
         << max(sensor_response_f_grid) << ").\n";
      error_found = true;
    }

  // Check number of columns in backend_channel_response
  //
  const Index nrp = backend_channel_response.nelem();
  //
  if( nrp != 1  &&  nrp != f_backend.nelem() ) 
    {
      os << "The WSV *backend_channel_response* must have 1 or n elements,\n"
         << "where n is the length of *f_backend*.\n"; 
      error_found = true;
    }

  // If errors where found throw runtime_error with the collected error
  // message (before error message gets too long).
  if( error_found )
    throw runtime_error(os.str());

  Numeric f_dlow  = 0.0;
  Numeric f_dhigh = 0.0;

  for( Index i=0; i<nrp; i++ )
    {
      ConstVectorView bchr_f_grid =     
                   backend_channel_response[i].get_numeric_grid(GFIELD1_F_GRID);

      if( bchr_f_grid.nelem() != backend_channel_response[i].data.nelem() )
        {
          os << "Mismatch in size of grid and data in element " << i
             << "\nof *sideband_response*.\n"; 
          error_found = true;
        }

      if( !is_increasing( bchr_f_grid ) ) 
        {
          os << "The frequency grid of element " << i
             << " in *backend_channel_response*\nis not strictly increasing.\n"; 
          error_found = true;
        }

      // Check if the relative grid added to the channel frequencies expands
      // outside sensor_response_f_grid.
      //
      Numeric f1 = f_backend[i] + bchr_f_grid[0] - min(sensor_response_f_grid);
      Numeric f2 = (max(sensor_response_f_grid) - 
                    f_backend[i]) - last(bchr_f_grid);
      //
      f_dlow  = min( f_dlow, f1 );
      f_dhigh = min( f_dhigh, f2 );
    }

  if( f_dlow < 0 ) 
  {
    os << "The WSV *sensor_response_f_grid* is too narrow. It should be\n"
       << "expanded with "<<-f_dlow<<" Hz in the lower end. This change\n"
       << "should be applied to either *f_grid* or the sensor part in\n"
       << "front of *sensor_responseBackend*.\n";
    error_found = true;
  }
  if( f_dhigh < 0 ) 
  {
    os << "The WSV *sensor_response_f_grid* is too narrow. It should be\n"
       << "expanded with "<<-f_dhigh<<" Hz in the higher end. This change\n"
       << "should be applied to either *f_grid* or the sensor part in\n"
       << "front of *sensor_responseBackend*.\n";
    error_found = true;
  }

  // If errors where found throw runtime_error with the collected error
  // message.
  if (error_found)
    throw runtime_error(os.str());
  

  // Call the core function 
  //
  Sparse hbackend;
  //
  spectrometer_matrix( hbackend, f_backend, backend_channel_response,
                       sensor_response_f_grid, npol, nza, sensor_norm );

  // Here we need a temporary sparse that is copy of the sensor_response
  // sparse matrix. We need it since the multiplication function can not
  // take the same object as both input and output.
  Sparse htmp = sensor_response;
  sensor_response.resize( hbackend.nrows(), htmp.ncols());
  mult( sensor_response, hbackend, htmp );

  // Some extra output.
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  // Update sensor_response_f_grid
  sensor_response_f_grid = f_backend;

  // Set aux variables
  sensor_aux_vectors( sensor_response_f,       sensor_response_pol, 
                      sensor_response_za,      sensor_response_aa, 
                      sensor_response_f_grid,  sensor_response_pol_grid, 
                      sensor_response_za_grid, sensor_response_aa_grid, 0 );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseBackendFrequencySwitching(Sparse&         sensor_response,
                                              Vector&         sensor_response_f,
                                              ArrayOfIndex&   sensor_response_pol,
                                              Vector&         sensor_response_za,
                                              Vector&         sensor_response_aa,
                                              Vector&         sensor_response_f_grid,
                                              const ArrayOfIndex&   sensor_response_pol_grid,
                                              const Vector&   sensor_response_za_grid,
                                              const Vector&   sensor_response_aa_grid,
                                              const Vector&   f_backend,
                                              const ArrayOfGriddedField1&   backend_channel_response,
                                              const Index&    sensor_norm,
                                              const Numeric&  df1,
                                              const Numeric&  df2,
                                              const Verbosity& verbosity)
{
  // All needed checks are done in sensor_responseBackend

  Sparse H1=sensor_response, H2=sensor_response;

  // Some needed vectors
  Vector f_backend_shifted;
  Vector fdummy=sensor_response_f, fdummy_grid=sensor_response_f_grid;

  // Cycle 1
  f_backend_shifted  = f_backend;
  f_backend_shifted += df1;
  //
  sensor_responseBackend(H1, fdummy, sensor_response_pol, 
                         sensor_response_za, sensor_response_aa, 
                         fdummy_grid, sensor_response_pol_grid, 
                         sensor_response_za_grid, sensor_response_aa_grid, 
                         f_backend_shifted, backend_channel_response, 
                         sensor_norm, verbosity);
  // Cycle 2
  f_backend_shifted  = f_backend;
  f_backend_shifted += df2;
  //
  sensor_responseBackend(H2, sensor_response_f, sensor_response_pol, 
                         sensor_response_za, sensor_response_aa, 
                         sensor_response_f_grid, sensor_response_pol_grid, 
                         sensor_response_za_grid, sensor_response_aa_grid, 
                         f_backend_shifted, backend_channel_response, 
                         sensor_norm, verbosity);

  // Total response
  sub( sensor_response, H2, H1 );

  // sensor_response_f_grid shall be f_backend
  sensor_response_f_grid = f_backend;

  // Set aux variables
  sensor_aux_vectors( sensor_response_f,       sensor_response_pol, 
                      sensor_response_za,      sensor_response_aa, 
                      sensor_response_f_grid,  sensor_response_pol_grid, 
                      sensor_response_za_grid, sensor_response_aa_grid, 0 );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseBeamSwitching(// WS Output:
                                  Sparse&         sensor_response,
                                  Vector&         sensor_response_f,
                                  ArrayOfIndex&   sensor_response_pol,
                                  Vector&         sensor_response_za,
                                  Vector&         sensor_response_aa,
                                  Vector&         sensor_response_za_grid,
                                  Vector&         sensor_response_aa_grid,
                                  // WS Input:
                                  const Vector&   sensor_response_f_grid,
                                  const ArrayOfIndex&   sensor_response_pol_grid,
                                  const Numeric&   w1,
                                  const Numeric&   w2,
                                  const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  if( sensor_response_za_grid.nelem() != 2 )
    throw runtime_error( 
       "This method requires that the number of observation directions is 2." );

  if( sensor_response_pol_grid.nelem() != 1 )
    throw runtime_error( 
               "This method handles (so far) only single polarisation cases." );

  const Index n = sensor_response_f_grid.nelem();

  // Form H matrix representing beam switching
  Sparse Hbswitch( n, 2*n );
  Vector hrow( 2*n, 0.0 );
  //
  for( Index i=0; i<n; i++ )
    {
      hrow[i]   = w1;
      hrow[i+n] = w2;
      //
      Hbswitch.insert_row( i, hrow );
      //
      hrow = 0;
    }

  // Here we need a temporary sparse that is copy of the sensor_response
  // sparse matrix. We need it since the multiplication function can not
  // take the same object as both input and output.
  Sparse Htmp = sensor_response;
  sensor_response.resize( Hbswitch.nrows(), Htmp.ncols() );
  mult( sensor_response, Hbswitch, Htmp );

  // Some extra output.
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  // Update sensor_response_za_grid
  const Numeric za = sensor_response_za_grid[1];
  sensor_response_za_grid.resize(1);
  sensor_response_za_grid[0] = za;

  // Update sensor_response_aa_grid
  if( sensor_response_aa_grid.nelem() > 0 )
    {
      const Numeric aa = sensor_response_aa_grid[1];
      sensor_response_aa_grid.resize(1);
      sensor_response_aa_grid[0] = aa;
    }

  // Set aux variables
  sensor_aux_vectors( sensor_response_f,       sensor_response_pol, 
                      sensor_response_za,      sensor_response_aa, 
                      sensor_response_f_grid,  sensor_response_pol_grid, 
                      sensor_response_za_grid, sensor_response_aa_grid, 0 );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseFrequencySwitching(// WS Output:
                                       Sparse&         sensor_response,
                                       Vector&         sensor_response_f,
                                       ArrayOfIndex&   sensor_response_pol,
                                       Vector&         sensor_response_za,
                                       Vector&         sensor_response_aa,
                                       Vector&         sensor_response_f_grid,
                                       // WS Input:
                                       const ArrayOfIndex&   sensor_response_pol_grid,
                                       const Vector&   sensor_response_za_grid,
                                       const Vector&   sensor_response_aa_grid,
                                       const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  if( sensor_response_za_grid.nelem() != 1 )
    throw runtime_error( 
       "This method requires that the number of observation directions is 1." );

  if( sensor_response_pol_grid.nelem() != 1 )
    throw runtime_error( 
               "This method handles (so far) only single polarisation cases." );

  const Index n  = sensor_response_f_grid.nelem();
  const Index n2 = n/2;

  if( sensor_response.nrows() != n )
    throw runtime_error( "Assumptions of method are not fulfilled, "
                         "considering number of rows in *sensor_response* "
                         "and length of *sensor_response_f_grid*." );

  if( !is_multiple(n,2) )
    throw runtime_error( "There is an odd number of total frequencies, "
                         "which is not consistent with the assumptions of "
                         "the method." );


  // Form H matrix representing frequency switching
  Sparse Hbswitch( n2, n );
  Vector hrow( n, 0.0 );
  //
  for( Index i=0; i<n2; i++ )
    {
      hrow[i]    = -1;
      hrow[i+n2] = 1;
      //
      Hbswitch.insert_row( i, hrow );
      //
      hrow = 0;
    }

  // Here we need a temporary sparse that is copy of the sensor_response
  // sparse matrix. We need it since the multiplication function can not
  // take the same object as both input and output.
  Sparse Htmp = sensor_response;
  sensor_response.resize( Hbswitch.nrows(), Htmp.ncols() );
  mult( sensor_response, Hbswitch, Htmp );

  // Some extra output.
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  // Update sensor_response_f_grid
  const Vector f = sensor_response_f_grid;
  sensor_response_f_grid.resize(n2);
  sensor_response_f_grid = f[Range(n2,n2)];

  // Set aux variables
  sensor_aux_vectors( sensor_response_f,       sensor_response_pol, 
                      sensor_response_za,      sensor_response_aa, 
                      sensor_response_f_grid,  sensor_response_pol_grid, 
                      sensor_response_za_grid, sensor_response_aa_grid, 0 );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseIF2RF(// WS Output:
                          Vector&         sensor_response_f,
                          Vector&         sensor_response_f_grid,
                          // WS Input:
                          const Numeric&  lo,
                          const String&   sideband_mode,
                          const Verbosity&)
{
  // Check that frequencies are not too high. This might be a floating limit.
  // For this we use the variable f_lim, given in Hz.
  Numeric f_lim = 30e9;
  if( max(sensor_response_f_grid) > f_lim )
    throw runtime_error( "The frequencies seem to already be given in RF." );


  // Lower band
  if( sideband_mode == "lower" ) 
    {
      sensor_response_f      *= -1;
      sensor_response_f_grid *= -1;
      sensor_response_f      += lo;
      sensor_response_f_grid += lo;
    }

  // Upper band
  else if( sideband_mode=="upper" ) 
    {
      sensor_response_f      += lo;
      sensor_response_f_grid += lo;
    }

  // Unknown option
  else
    {
      throw runtime_error(
      "Only allowed options for *sideband _mode* are \"lower\" and \"upper\"." );
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseFillFgrid(// WS Output:
                              Sparse&         sensor_response,
                              Vector&         sensor_response_f,
                              ArrayOfIndex&   sensor_response_pol,
                              Vector&         sensor_response_za,
                              Vector&         sensor_response_aa,
                              Vector&         sensor_response_f_grid,
                              // WS Input:
                              const ArrayOfIndex&   sensor_response_pol_grid,
                              const Vector&   sensor_response_za_grid,
                              const Vector&   sensor_response_aa_grid,
                              const Index&    polyorder,
                              const Index&    nfill,
                              const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Some sizes
  const Index nf   = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nza  = sensor_response_za_grid.nelem();
  const Index naa  = max( Index(1), sensor_response_aa_grid.nelem() );
  const Index nin  = nf * npol * nza * naa;

  // Initialise a output stream for runtime errors and a flag for errors
  ostringstream os;
  bool          error_found = false;

  // Check that sensor_response variables are consistent in size
  if( sensor_response_f.nelem() != nin )
  {
    os << "Inconsistency in size between *sensor_response_f* and the sensor\n"
       << "grid variables (sensor_response_f_grid etc.).\n";
    error_found = true;
  }
  if( sensor_response.nrows() != nin )
  {
    os << "The sensor block response matrix *sensor_response* does not have\n"
       << "right size compared to the sensor grid variables\n"
       << "(sensor_response_f_grid etc.).\n";
    error_found = true;
  }

  // Check polyorder and nfill
  if( polyorder < 2  ||  polyorder > 7 )
  {
    os << "Accepted range for *polyorder* is [3,7].\n";
    error_found = true;
  }
  if( nfill < 1 )
  {
    os << "The argument *nfill* must be > 1.\n";
    error_found = true;
  }

  // If errors where found throw runtime_error with the collected error
  // message.
  if (error_found)
    throw runtime_error(os.str());


  // New frequency grid
  //
  const Index n1   = nfill+1;
  const Index n2   = nfill+2;
  const Index nnew = (nf-1)*n1 + 1;
  //
  Vector fnew( nnew );
  //
  for( Index i=0; i<nf-1; i++ )
    {
      Vector fp(n2);
      nlinspace( fp, sensor_response_f_grid[i], sensor_response_f_grid[i+1], n2 );
      fnew[Range(i*n1,n2)] = fp;
    }
  
  // Find interpolation weights
  //
  ArrayOfGridPosPoly   gp( nnew );
  Matrix               itw( nnew, polyorder+1 );
  //
  gridpos_poly( gp, sensor_response_f_grid, fnew, polyorder );
  interpweights( itw, gp );

  // Set up H for this part
  //
  Sparse hpoly( nnew * npol * nza * naa, nin );
  Vector hrow( nin, 0.0 );
  Index  row = 0;
  //
  for( Index iza=0; iza<nza; iza++ )
    {
      for( Index iaa=0; iaa<naa; iaa++ )
        { 
          for( Index iv=0; iv<nnew; iv++ )
            { 
              for( Index ip=0; ip<npol; ip++ )
                {  
                  const Index col0 = (iza*naa+iaa)*nf*npol;
                  for( Index i=0; i<gp[iv].idx.nelem(); i++ )
                    { 
                      const Numeric w = gp[iv].w[i];
                      if( abs(w) > 1e-5 )
                        { 
                          hrow[col0+gp[iv].idx[i]*npol+ip] = w; 
                        }
                    }
                  hpoly.insert_row( row, hrow );
                  for( Index i=0; i<gp[iv].idx.nelem(); i++ )
                    { hrow[col0+gp[iv].idx[i]*npol+ip] = 0; }
                  row += 1;
                }
            }
        }
    }

  // Here we need a temporary sparse that is copy of the sensor_response
  // sparse matrix. We need it since the multiplication function can not
  // take the same object as both input and output.
  Sparse htmp = sensor_response;
  sensor_response.resize( hpoly.nrows(), htmp.ncols());
  mult( sensor_response, hpoly, htmp );

  // Some extra output.
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  // Update sensor_response_za_grid
  sensor_response_f_grid = fnew;

  // Set aux variables
  sensor_aux_vectors( sensor_response_f,       sensor_response_pol, 
                      sensor_response_za,      sensor_response_aa, 
                      sensor_response_f_grid,  sensor_response_pol_grid, 
                      sensor_response_za_grid, sensor_response_aa_grid, 1 );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseInit(// WS Output:
                         Sparse&         sensor_response,
                         Vector&         sensor_response_f,
                         ArrayOfIndex&   sensor_response_pol,
                         Vector&         sensor_response_za,
                         Vector&         sensor_response_aa,
                         Vector&         sensor_response_f_grid,
                         ArrayOfIndex&   sensor_response_pol_grid,
                         Vector&         sensor_response_za_grid,
                         Vector&         sensor_response_aa_grid,
                         // WS Input:
                         const Vector&   f_grid,
                         const Vector&   mblock_za_grid,
                         const Vector&   mblock_aa_grid,
                         const Index&    antenna_dim,
                         const Index&    atmosphere_dim,
                         const Index&    stokes_dim,
                         const Index&    sensor_norm,
                         const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  // Check input

  // Basic variables
  chk_if_in_range( "stokes_dim",  stokes_dim,  1, 4 );
  chk_if_in_range( "antenna_dim", antenna_dim, 1, 2 );
  chk_if_bool(     "sensor_norm", sensor_norm       );

  // f_grid (could in fact be decreasing, but an increasing grid is
  // demanded in other parts).
  chk_if_increasing( "f_grid", f_grid );
  
  // mblock_za_grid
  if( mblock_za_grid.nelem() == 0 )
    throw runtime_error( "The measurement block zenith angle grid is empty." );
  if( !is_increasing(mblock_za_grid)  &&  !is_decreasing(mblock_za_grid) )  
    throw runtime_error( 
        "The WSV *mblock_za_grid* must be strictly increasing or decreasing." );

  // mblock_aa_grid
  if( antenna_dim == 1 )
    {
      if( mblock_aa_grid.nelem() != 0 )
        throw runtime_error( 
              "For antenna_dim = 1, the azimuthal angle grid must be empty." );
    }
  else
    {
      if( atmosphere_dim < 3 )
        throw runtime_error( "2D antennas (antenna_dim=2) can only be "
                                                 "used with 3D atmospheres." );
      if( mblock_aa_grid.nelem() == 0 )
        {
          ostringstream os;
          os << "The measurement block azimuthal angle grid is empty despite"
             << "a 2D antenna pattern is flagged (*antenna_dim*).";
          throw runtime_error( os.str() );
        }
      if( !is_increasing(mblock_aa_grid)  &&  !is_decreasing(mblock_aa_grid) )  
        throw runtime_error( 
        "The WSV *mblock_aa_grid* must be strictly increasing or decreasing." );
    }


  // Set grid variables
  sensor_response_f_grid   = f_grid;
  sensor_response_za_grid  = mblock_za_grid;
  sensor_response_aa_grid  = mblock_aa_grid;
  //
  sensor_response_pol_grid.resize(stokes_dim);
  //
  for( Index is=0; is<stokes_dim; is++ )
    {
      sensor_response_pol_grid[is] = is + 1;
    }


  // Set aux variables
  sensor_aux_vectors( sensor_response_f,       sensor_response_pol, 
                      sensor_response_za,      sensor_response_aa, 
                      sensor_response_f_grid,  sensor_response_pol_grid, 
                      sensor_response_za_grid, sensor_response_aa_grid, 1 );

  //Set response matrix to identity matrix
  //
  const Index   n = sensor_response_f.nelem();
  //
  out2 << "  Initialising *sensor_reponse* as a identity matrix.\n";
  out3 << "  Size of *sensor_response*: " << n << "x" << n << "\n";
  //
  sensor_response.make_I( n, n );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensorOff(// WS Output:
               Sparse&         sensor_response,
               Vector&         sensor_response_f,
               ArrayOfIndex&   sensor_response_pol,
               Vector&         sensor_response_za,
               Vector&         sensor_response_aa,
               Vector&         sensor_response_f_grid,
               ArrayOfIndex&   sensor_response_pol_grid,
               Vector&         sensor_response_za_grid,
               Vector&         sensor_response_aa_grid,
               Index&          antenna_dim,
               Vector&         mblock_za_grid,
               Vector&         mblock_aa_grid,
               const Index&    stokes_dim,
               const Vector&   f_grid,
               const Verbosity& verbosity)
{
  // Checks are done in sensor_responseInit.

  AntennaOff( antenna_dim, mblock_za_grid, mblock_aa_grid, verbosity );

  // Dummy variables (The method is independent of atmosphere_dim.
  // atmosphere_dim used below just for some checks when antenna is active, not
  // relevant here. ).
  const Index sensor_norm = 1, atmosphere_dim = 1;

  sensor_responseInit( sensor_response, sensor_response_f, 
                  sensor_response_pol, sensor_response_za, sensor_response_aa, 
                  sensor_response_f_grid, sensor_response_pol_grid, 
                  sensor_response_za_grid, sensor_response_aa_grid, f_grid, 
                  mblock_za_grid, mblock_aa_grid, antenna_dim, atmosphere_dim, 
                  stokes_dim, sensor_norm, verbosity );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseMixer(// WS Output:
                          Sparse&         sensor_response,
                          Vector&         sensor_response_f,
                          ArrayOfIndex&   sensor_response_pol,
                          Vector&         sensor_response_za,
                          Vector&         sensor_response_aa,
                          Vector&         sensor_response_f_grid,
                          // WS Input:
                          const ArrayOfIndex&   sensor_response_pol_grid,
                          const Vector&   sensor_response_za_grid,
                          const Vector&   sensor_response_aa_grid,
                          const Numeric&  lo,
                          const GriddedField1&   sideband_response,
                          const Index&    sensor_norm,
                          const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Some sizes
  const Index nf   = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nza  = sensor_response_za_grid.nelem();
  const Index naa  = sensor_response_aa_grid.nelem();
  const Index nin  = nf * npol * nza;
  // Note that there is no distinction between za and aa grids after the antenna

  // Frequency grid of for sideband response specification
  ConstVectorView sbresponse_f_grid = 
                             sideband_response.get_numeric_grid(GFIELD1_F_GRID);

  // Initialise a output stream for runtime errors and a flag for errors
  ostringstream os;
  bool          error_found = false;

  // Check that sensor_response variables are consistent in size
  if( sensor_response_f.nelem() != nin )
  {
    os << "Inconsistency in size between *sensor_response_f* and the sensor\n"
       << "grid variables (sensor_response_f_grid etc.).\n";
    error_found = true;
  }
  if( naa  &&  naa != nza )
  {
    os << "Incorrect size of *sensor_response_aa_grid*.\n";
    error_found = true;
  }
  if( sensor_response.nrows() != nin )
  {
    os << "The sensor block response matrix *sensor_response* does not have\n"
       << "right size compared to the sensor grid variables\n"
       << "(sensor_response_f_grid etc.).\n";
    error_found = true;
  }

  // Check that the lo frequency is within the sensor_response_f_grid
  if( lo <= sensor_response_f_grid[0]  ||  lo >= last(sensor_response_f_grid) )
  {
    os << "The given local oscillator frequency is outside the sensor\n"
       << "frequency grid. It must be within the *sensor_response_f_grid*.\n";
    error_found = true;
  }

  // Checks of sideband_response, partly in combination with lo
  if( sbresponse_f_grid.nelem() != sideband_response.data.nelem() )
    {
      os << "Mismatch in size of grid and data in *sideband_response*.\n"; 
      error_found = true;
    }
  if( sbresponse_f_grid.nelem() < 2 )
    {
      os << "At least two data points must be specified in "
         << "*sideband_response*.\n"; 
      error_found = true;
    }
  if( !is_increasing( sbresponse_f_grid ) ) 
    {
      os << "The frequency grid of *sideband_response* must be strictly\n"
         << "increasing.\n"; 
      error_found = true;
    }
  if( fabs(last(sbresponse_f_grid)+sbresponse_f_grid[0]) > 0 )
    {
      os << "The end points of the *sideband_response* frequency grid must be\n"
         << "symmetrically placed around 0. That is, the grid shall cover a\n"
         << "a range that can be written as [-df,df]. \n";
      error_found = true;      
    }

  // Check that response function does not extend outside sensor_response_f_grid
  Numeric df_high = lo + last(sbresponse_f_grid) - last(sensor_response_f_grid);
  Numeric df_low  = sensor_response_f_grid[0] - lo - sbresponse_f_grid[0];
  if( df_high > 0  &&  df_low > 0 )
  {
    os << "The *sensor_response_f* grid must be extended by at least\n"
       << df_low << " Hz in the lower end and " << df_high << " Hz in the\n"
       << "upper end to cover frequency range set by *sideband_response*\n"
       << "and *lo*. Or can the frequency grid of *sideband_response* be\n"
       << "decreased?";
    error_found = true;
  }
  else if( df_high > 0 )
  {
    os << "The *sensor_response_f* grid must be extended by at " << df_high 
       << " Hz\nin the upper end to cover frequency range set by\n"
       << "*sideband_response* and *lo*. Or can the frequency grid of\n"
       << "*sideband_response* be decreased?";
    error_found = true;
  }
  else if( df_low > 0 )
  {
    os << "The *sensor_response_f* grid must be extended by at " << df_low
       << " Hz\nin the lower end to cover frequency range set by\n"
       << "*sideband_response* and *lo*. Or can the frequency grid of\n"
       << "*sideband_response* be decreased?";
    error_found = true;
  }

  // If errors where found throw runtime_error with the collected error
  // message.
  if (error_found)
    throw runtime_error(os.str());


  //Call the core function
  //
  Sparse hmixer;
  Vector f_mixer;
  //
  mixer_matrix( hmixer, f_mixer, lo, sideband_response, 
                sensor_response_f_grid, npol, nza, sensor_norm );

  // Here we need a temporary sparse that is copy of the sensor_response
  // sparse matrix. We need it since the multiplication function can not
  // take the same object as both input and output.
  Sparse htmp = sensor_response;
  sensor_response.resize( hmixer.nrows(), htmp.ncols() );
  mult( sensor_response, hmixer, htmp );

  // Some extra output.
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  // Update sensor_response_f_grid
  sensor_response_f_grid = f_mixer;

  // Set aux variables
  sensor_aux_vectors( sensor_response_f,       sensor_response_pol, 
                      sensor_response_za,      sensor_response_aa, 
                      sensor_response_f_grid,  sensor_response_pol_grid, 
                      sensor_response_za_grid, sensor_response_aa_grid, 0 );
}



void sensor_responseMultiMixerBackend(// WS Output:
                                      Sparse&         sensor_response,
                                      Vector&         sensor_response_f,
                                      ArrayOfIndex&   sensor_response_pol,
                                      Vector&         sensor_response_za,
                                      Vector&         sensor_response_aa,
                                      Vector&         sensor_response_f_grid,
                                      // WS Input:
                                      const ArrayOfIndex&   sensor_response_pol_grid,
                                      const Vector&   sensor_response_za_grid,
                                      const Vector&   sensor_response_aa_grid,
                                      const Vector&   lo_multi,
                                      const ArrayOfGriddedField1&   sideband_response_multi,
                                      const ArrayOfString&   sideband_mode_multi,
                                      const ArrayOfVector&   f_backend_multi,
                                      const ArrayOfArrayOfGriddedField1&   backend_channel_response_multi,
                                      const Index&   sensor_norm,
                                      const Verbosity& verbosity)
{
  // Some sizes
  const Index nf   = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nza  = sensor_response_za_grid.nelem();
  const Index naa  = sensor_response_aa_grid.nelem();
  const Index nin  = nf * npol * nza;
  // Note that there is no distinction between za and aa grids after the antenna
  const Index nlo  = lo_multi.nelem();

  // Initialise a output stream for runtime errors and a flag for errors
  ostringstream os;
  bool          error_found = false;

  // Check that sensor_response variables are consistent in size
  if( sensor_response_f.nelem() != nin )
  {
    os << "Inconsistency in size between *sensor_response_f* and the sensor\n"
       << "grid variables (sensor_response_f_grid etc.).\n";
    error_found = true;
  }
  if( naa  &&  naa != nza )
  {
    os << "Incorrect size of *sensor_response_aa_grid*.\n";
    error_found = true;
  }
  if( sensor_response.nrows() != nin )
  {
    os << "The sensor block response matrix *sensor_response* does not have\n"
       << "right size compared to the sensor grid variables\n"
       << "(sensor_response_f_grid etc.).\n";
    error_found = true;
  }

  // Check that response data are consistent with respect to number of
  // mixer/reciever chains.
  if( sideband_response_multi.nelem() != nlo )
  {
    os << "Inconsistency in length between *lo_mixer* and "
       << "*sideband_response_multi*.\n";
    error_found = true;
  }
  if( sideband_mode_multi.nelem() != nlo )
  {
    os << "Inconsistency in length between *lo_mixer* and "
       << "*sideband_mode_multi*.\n";
    error_found = true;
  }
  if( f_backend_multi.nelem() != nlo )
  {
    os << "Inconsistency in length between *lo_mixer* and "
       << "*f_backend_multi*.\n";
    error_found = true;
  }
  if( backend_channel_response_multi.nelem() != nlo )
  {
    os << "Inconsistency in length between *lo_mixer* and "
       << "*backend_channel_response_multi*.\n";
    error_found = true;
  }

  // If errors where found throw runtime_error with the collected error
  // message. Data for each mixer and reciever chain are checked below.
  if (error_found)
    throw runtime_error(os.str());


  // Variables for data to be appended
  Array<Sparse> sr;
  ArrayOfVector srfgrid;
  ArrayOfIndex  cumsumf(nlo+1,0);

  for( Index ilo=0; ilo<nlo; ilo++ )
    {
      // Copies of variables that will be changed, but must be
      // restored for next loop
      Sparse       sr1      = sensor_response;
      Vector       srf1     = sensor_response_f;
      ArrayOfIndex srpol1   = sensor_response_pol;
      Vector       srza1    = sensor_response_za;
      Vector       sraa1    = sensor_response_aa;
      Vector       srfgrid1 = sensor_response_f_grid;

      // Call single reciever methods. Try/catch for improved error message.
      try
        {
          sensor_responseMixer(sr1, srf1, srpol1, srza1, sraa1, srfgrid1,
                               sensor_response_pol_grid,
                               sensor_response_za_grid, 
                               sensor_response_aa_grid,
                               lo_multi[ilo], 
                               sideband_response_multi[ilo], 
                               sensor_norm,
                               verbosity);

          sensor_responseIF2RF(srf1, srfgrid1,
                               lo_multi[ilo], 
                               sideband_mode_multi[ilo],
                               verbosity);

          sensor_responseBackend(sr1, srf1, srpol1, srza1, sraa1, srfgrid1,
                                 sensor_response_pol_grid,
                                 sensor_response_za_grid, 
                                 sensor_response_aa_grid,
                                 f_backend_multi[ilo],
                                 backend_channel_response_multi[ilo],
                                 sensor_norm, verbosity);
        } 
      catch( runtime_error e ) 
        {
          ostringstream os2;
          os2 << "Error when dealing with receiver/mixer chain (1-based index) " 
              << ilo+1 << ":\n" << e.what();
          throw runtime_error(os2.str());
        }

      // Store in temporary arrays
      sr.push_back( sr1 );
      srfgrid.push_back( srfgrid1 );
      //
      cumsumf[ilo+1] = cumsumf[ilo] + srfgrid1.nelem();
    }

  // Append data to create sensor_response_f_grid
  //
  const Index  nfnew = cumsumf[nlo];
  sensor_response_f_grid.resize( nfnew );
  //
  for( Index ilo=0; ilo<nlo; ilo++ )
    {
      for( Index i=0; i<srfgrid[ilo].nelem(); i++ )
        {
          sensor_response_f_grid[cumsumf[ilo]+i] = srfgrid[ilo][i];
        }
    }

  // Append data to create total sensor_response
  //
  const Index  ncols   = sr[0].ncols();
  const Index  npolnew = sensor_response_pol_grid.nelem();
  const Index  nfpolnew = nfnew * npolnew;
  //
  sensor_response.resize( nza*nfpolnew, ncols );
  //
  Vector dummy( ncols, 0.0 );
  //
  for( Index ilo=0; ilo<nlo; ilo++ )
    {
      const Index nfpolthis = (cumsumf[ilo+1]-cumsumf[ilo]) * npolnew;

      assert( sr[ilo].nrows() == nza*nfpolthis );
      assert( sr[ilo].ncols() == ncols );

      for( Index iz=0; iz<nza; iz++ )
        {
          for( Index i=0; i<nfpolthis; i++ )
            {
              // "Poor mans" transfer of a row from one sparse to another 
              for( Index ic=0; ic<ncols; ic++ )
                { dummy[ic] = sr[ilo](iz*nfpolthis+i,ic); }

              sensor_response.insert_row( iz*nfpolnew+cumsumf[ilo]*npolnew+i, 
                                                                       dummy );
            }
        }
    }  

  // Set aux variables
  sensor_aux_vectors( sensor_response_f,       sensor_response_pol, 
                      sensor_response_za,      sensor_response_aa, 
                      sensor_response_f_grid,  sensor_response_pol_grid, 
                      sensor_response_za_grid, sensor_response_aa_grid, 0 );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responsePolarisation(// WS Output:
                                 Sparse&         sensor_response,
                                 Vector&         sensor_response_f,
                                 ArrayOfIndex&   sensor_response_pol,
                                 Vector&         sensor_response_za,
                                 Vector&         sensor_response_aa,
                                 ArrayOfIndex&   sensor_response_pol_grid,
                                 // WS Input:
                                 const Vector&   sensor_response_f_grid,
                                 const Vector&   sensor_response_za_grid,
                                 const Vector&   sensor_response_aa_grid,
                                 const Index&    stokes_dim,
                                 const String&   y_unit,
                                 const ArrayOfIndex&   sensor_pol,
                                 const Verbosity&)
{
  // Vectors for extracting polarisation components
  //
  Numeric w = 0.5;
  if( y_unit == "PlanckBT"  ||  y_unit == "RJBT"  )
    { w = 1.0; }
  //
  ArrayOfVector pv;
  stokes2pol( pv, w );

  // Some sizes
  const Index nnew = sensor_pol.nelem();
  const Index nf   = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nza  = sensor_response_za_grid.nelem();

  // Initialise an output stream for runtime errors and a flag for errors
  ostringstream os;
  bool          error_found = false;

  // For this method, the za and aa can be both "independent" and "dependent".
  // The size of sensor_response resolves this
  //
  bool za_aa_independent = true;
  //
  const Index naa  = max( Index(1), sensor_response_aa_grid.nelem() );
        Index nfz  = nf * nza;
        Index nin  = nfz *npol;
  //
  if( sensor_response.nrows() == nin )
    { za_aa_independent = false; }
  else if( sensor_response.nrows() == nin*naa )
    { nfz *= naa;  nin *= naa; }
  else
    {
      os << "The sensor block response matrix *sensor_response* does not have\n"
         << "right size compared to the sensor grid variables\n"
         << "(sensor_response_f_grid etc.).\n";
      error_found = true;
    }

  // Check that sensor_response variables are consistent in size
  if( sensor_response_f.nelem() != nin )
    {
      os << "Inconsistency in size between *sensor_response_f* and the sensor\n"
         << "grid variables (sensor_response_f_grid etc.).\n";
      error_found = true;
    }
  if( naa  &&  naa != nza )
    {
      os << "Incorrect size of *sensor_response_aa_grid*.\n";
      error_found = true;
    }
  if( npol != stokes_dim )
    {
      os << "Number of input polarisation does not match *stokes_dim*.\n";
      error_found = true;
    }
  if( nnew == 0 )
    {
      os << "The WSV *sensor_pol* can not be empty.\n";
      error_found = true;
    }
  // If errors where found throw runtime_error with the collected error
  // message (before it gets too long)
  if( error_found )
    throw runtime_error(os.str());

  // Check polarisation data more in detail
  for( Index i=0; i<npol && !error_found; i++ )
    {
      if( sensor_response_pol_grid[i] != i+1 )
        {
          os << "The input polarisations must be I, Q, U and V (up to "
             << "stokes_dim). It seems that input data are for other "
             << "polarisation components.";
          error_found = true;
        }      
    }
  for( Index i=0; i<nnew && !error_found; i++ )
    {
      if( sensor_pol[i] < 1  || sensor_pol[i] > 10 )
        {
          os << 
             "The elements of *sensor_pol* must be inside the range [1,10].\n";
          error_found = true;
        }
    }
  // If errors where found throw runtime_error with the collected error
  // message (before it gets too long)
  if( error_found )
    throw runtime_error(os.str());

  for( Index i=0; i<nnew && !error_found; i++ )
    {
      if( pv[sensor_pol[i]-1].nelem() > stokes_dim )
        {
          os << "You have selected an output polarisation that is not covered "
             << "by present value of *stokes_dim* (the later has to be "
             << "increased).";
          error_found = true;
        }
    }  
  // If errors where found throw runtime_error with the collected error
  // message 
  if( error_found )
    throw runtime_error(os.str());

  // Form H matrix representing polarisation response
  //
  Sparse Hpol( nfz*nnew, nin );
  Vector hrow( nin, 0.0 );
  Index row = 0;
  //
  for( Index i=0; i<nfz; i++ )
    {
      Index col = i*npol;
      for( Index in=0; in<nnew; in++ )
        {
          Index p = sensor_pol[in] - 1;
          //
          for( Index iv=0; iv<pv[p].nelem(); iv++ )
            { hrow[col+iv] = pv[p][iv]; }
          //
          Hpol.insert_row( row, hrow );
          //
          hrow = 0;
          row += 1;
        }
    }

  // Here we need a temporary sparse that is copy of the sensor_response
  // sparse matrix. We need it since the multiplication function can not
  // take the same object as both input and output.
  Sparse Htmp = sensor_response;
  sensor_response.resize( Hpol.nrows(), Htmp.ncols());
  mult( sensor_response, Hpol, Htmp );

  // Update sensor_response_pol_grid
  sensor_response_pol_grid = sensor_pol;

  // Set aux variables
  sensor_aux_vectors( sensor_response_f,       sensor_response_pol, 
                      sensor_response_za,      sensor_response_aa, 
                      sensor_response_f_grid,  sensor_response_pol_grid, 
                      sensor_response_za_grid, sensor_response_aa_grid, 
                      za_aa_independent );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseStokesRotation(
         Sparse&         sensor_response,
   const Vector&         sensor_response_f_grid,
   const ArrayOfIndex&   sensor_response_pol_grid,
   const Vector&         sensor_response_za_grid,
   const Vector&         sensor_response_aa_grid,
   const Index&          stokes_dim,
   const Matrix&         stokes_rotation,
   const Verbosity& )
{
  // Basic checks
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );

  // Some sizes
  const Index nf   = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nza  = sensor_response_za_grid.nelem();
  const Index naa  = max( Index(1), sensor_response_aa_grid.nelem() );
  const Index nin  = nf * npol * nza * naa;


  //---------------------------------------------------------------------------
  // Initialise a output stream for runtime errors and a flag for errors
  ostringstream os;
  bool          error_found = false;

  // Check that sensor_response variables are consistent in size
  if( sensor_response.nrows() != nin )
  {
    os << "The sensor block response matrix *sensor_response* does not have\n"
       << "right size compared to the sensor grid variables\n"
       << "(sensor_response_f_grid etc.).\n";
    error_found = true;
  }
  
  // Check special stuff for this method
  if( stokes_dim < 3 )
  {
    os << "To perform a rotation of the Stokes coordinate system,\n"
       << "*stokes_dim* must be >= 3.\n";
    error_found = true;
  }
  if( stokes_rotation.nrows() != nza  ||  stokes_rotation.ncols() != naa )  
  {
    os << "Incorrect number of angles in *stokes_rotation*. The size of this\n"
       << "matrix must match *sensor_response_za_grid* and "
       << "*sensor_response_aa_grid*.\n";
    error_found = true;
  }

  // If errors where found throw runtime_error with the collected error
  // message.
  if (error_found)
    throw runtime_error(os.str());
  //---------------------------------------------------------------------------


  // Set up the H matrix for applying rotation
  //
  Sparse hrot( sensor_response.nrows(), sensor_response.ncols() );
  Vector row( hrot.ncols(), 0 );
  Index  irow = -1;
  //
  for( Index iaa=0; iaa<naa; iaa++ )
    {
      for( Index iza=0; iza<nza; iza++ )
        {
          const Numeric a = cos( 2 * DEG2RAD * stokes_rotation(iza,iaa) );
          const Numeric b = sin( 2 * DEG2RAD * stokes_rotation(iza,iaa) );

          for( Index ifr=0; ifr<nf; ifr++ ) 
            {
              for( Index ip=0; ip<npol; ip++ )
                {
                  irow += 1;
                  if( ip==0  ||  ip==3 )  // Row 1 and 4 of Mueller matrix
                    { row[irow] = 1; }
                  else if( ip==1 )        // Row 2
                    { row[irow] = a; row[irow+1] = b; }
                  else                    // Row 3
                    { row[irow] = a; row[irow-1] = -b; }
                  hrot.insert_row( irow, row );
                  // Make sure that set values are not affecting next row
                  // (note that column irow+1 always will be overwritten)
                  row[irow] = 0;
                  if( ip == 2 )
                    { row[irow-1] = 0; }
                }
            }
        }
    }

  // Here we need a temporary sparse that is copy of the sensor_response
  // sparse matrix. We need it since the multiplication function can not
  // take the same object as both input and output.
  Sparse htmp = sensor_response;
  sensor_response.resize( htmp.nrows(), htmp.ncols());  //Just in case!
  mult( sensor_response, hrot, htmp );  
}



void sensor_responseGenericAMSU(// WS Output:
    Vector& f_grid,
    Index& antenna_dim,
    Vector& mblock_za_grid,
    Vector& mblock_aa_grid,
    Sparse& sensor_response,
    Vector& sensor_response_f,
    ArrayOfIndex& sensor_response_pol,
    Vector& sensor_response_za,
    Vector& sensor_response_aa,
    Vector& sensor_response_f_grid,
    ArrayOfIndex& sensor_response_pol_grid,
    Vector& sensor_response_za_grid,
    Vector& sensor_response_aa_grid,
    Index& sensor_norm,
    // WS Input:
    const Index& atmosphere_dim,
    const Index& stokes_dim,
    const Matrix& sensor_description_amsu,
    // WS Generic Input:
    const Numeric& spacing,
    const Verbosity& verbosity)
    {
      // Number of instrument channels:
      const Index n = sensor_description_amsu.nrows(); 
      const Index m = sensor_description_amsu.ncols(); 

      // FIXME  - Oscar Isoz-090413
      // add error checks 
      //

      // The meaning of the columns in sensor_description_amsu is:
      // LO frequency, channel center offset from LO, channel width, second offset from LO (to be extended if needed);
      // (All in Hz.)
      //
      if(5>sensor_description_amsu.ncols() )
      {
        ostringstream os;
        os << "Input variable sensor_description_amsu must have atleast five columns, but it has "
          << sensor_description_amsu.ncols() << ".";
        throw runtime_error( os.str() );
      }

      ConstVectorView lo_multi = sensor_description_amsu(Range(joker),0);
      ConstMatrixView offset   = sensor_description_amsu(Range(joker),Range(1,m-2)); // Remember to ignore column 2..
      ConstVectorView verbosityVectIn = sensor_description_amsu(Range(joker),m-1);
      ConstVectorView width    = sensor_description_amsu(Range(joker),m-2);



      //is there any undefined verbosity values in the vector?
      //Set the verbosity to one third of the bandwidth to make sure that the passband flanks does't overlap
      const Numeric minRatioVerbosityVsFdiff = 10;  // To be used when to passbands are closer then one verbosity value  
      Index totNumPb =0;
      Vector verbosityVect(n);

      for(Index idx = 0;idx<n;++idx)
      {
        if((verbosityVectIn[idx] ==0) || (verbosityVectIn[idx]> width[idx]))
        {
          verbosityVect[idx] = ((Numeric)width[idx])/3;
        }
        else 
        {
          verbosityVect[idx] = verbosityVectIn[idx];
        }	
      }

      // Create a vector to store the number of passbands (PB) for each channel 
      Vector numPBpseudo(n); // store values used for calc
      Vector numPB(n); // Store the true values 
      // Find the number of IFs for each channel and calculate the number of passbands 
      for (Index i=0;i<n;++i)
      {
        numPB[i] = 0 ; // make sure that it is zero
        for(Index j=0;j<(m-2);++j)
        {
          if(j!=2)
          {
            if (offset(i,j)>0)
            {
              numPB[i]++;
            }
          }
        }  
        numPB[i] = 1 << (int)numPB[i] ;// number of passbands= 2^numLO
        if(numPB[i] ==1)
        {
          numPBpseudo[i] = 2;
        } 
        else
        {
          numPBpseudo[i] = numPB[i];
        }

        totNumPb += (int)numPB[i];

        if(numPB[i]>4)
        {
          ostringstream os;
          os << "This function does currently not support more than 4 passbands per channel"
            << numPB[i]<< ".";
          throw runtime_error( os.str() );
        }
      }  

      // Find the center frequencies for all sub-channels
      // Create one center frequency for each passband
      ArrayOfArrayOfGriddedField1 backend_channel_response_multi(n);
      ArrayOfVector f_backend_multi(n); // changed !!!
      for (Index i=0; i<n; ++i)
      {
        // Channel frequencies
        Vector &f = f_backend_multi[i];
        f.resize(1);
        f[0] = lo_multi[i]+0.0*width[i];//(offset(i,0)+offset(i,1));

        //channel response 
        const Index numVal = 4;
        backend_channel_response_multi[i].resize(1);
        GriddedField1& b_resp = backend_channel_response_multi[i][0];
        b_resp.set_name("Backend channel response function");
        b_resp.resize(numVal*(Index)numPBpseudo[i]);
        Vector f_range(numVal*(Index)numPBpseudo[i]);
        Numeric pbOffset = 0;
        b_resp.set_grid_name(0,"Frequency");

        Numeric slope =0; // 1900; 
        // To avoid overlapping passbands in the AMSU-A sensor, reduce the passbands of each channel by a few Hz
        for(Index pbOffsetIdx = 0;pbOffsetIdx<numPBpseudo[i];++pbOffsetIdx)
        { // Filter response 
          slope = width[i]/100;
          f_range[pbOffsetIdx*numVal+0] = -0.5*width[i]-0*slope;
          f_range[pbOffsetIdx*numVal+1] = -0.5*width[i]+1*slope;
          f_range[pbOffsetIdx*numVal+2] = +0.5*width[i]-1*slope;
          f_range[pbOffsetIdx*numVal+3] = +0.5*width[i]+0*slope;

          b_resp.data[pbOffsetIdx*numVal+0] = 0.0/numPB[i];;
          b_resp.data[pbOffsetIdx*numVal+1] = 1.0/numPB[i];
          b_resp.data[pbOffsetIdx*numVal+2] = 1.0/numPB[i];
          b_resp.data[pbOffsetIdx*numVal+3] = 0.0/numPB[i];

          if(numPB[i] ==1)
          {
            if(pbOffsetIdx==0)
            {
              pbOffset = -0.0*width[i];
              b_resp.data[pbOffsetIdx*numVal+0] = 0;
              b_resp.data[pbOffsetIdx*numVal+1] = 1.0/1;

              b_resp.data[pbOffsetIdx*numVal+2] = 1.0/1;
              b_resp.data[pbOffsetIdx*numVal+3] = 1.0/1;
              f_range[pbOffsetIdx*numVal+0] = -0.5*width[i]-2*slope;
              f_range[pbOffsetIdx*numVal+1] = -0.5*width[i]-1*slope;
              f_range[pbOffsetIdx*numVal+2] = -0.5*width[i]+1*slope;
              f_range[pbOffsetIdx*numVal+3] = -0.5*width[i]+2*slope;
            }
            if(pbOffsetIdx ==1)
            {
              pbOffset = 0.0*width[i]; //just a dummy band 
              b_resp.data[pbOffsetIdx*numVal+0] = 1.0/1;
              b_resp.data[pbOffsetIdx*numVal+1] = 1.0/1;
              b_resp.data[pbOffsetIdx*numVal+2] = 1.0/1;
              b_resp.data[pbOffsetIdx*numVal+3] = 0;
              f_range[pbOffsetIdx*numVal+0] = +0.5*width[i]-3*slope;
              f_range[pbOffsetIdx*numVal+1] = +0.5*width[i]-2*slope;
              f_range[pbOffsetIdx*numVal+2] = +0.5*width[i]-1*slope;
              f_range[pbOffsetIdx*numVal+3] = +0.5*width[i]+0*slope-10; 
              // without the extra '-10' it will fail due to too narrow backend sensor response.
            }
          }
          else if (numPB[i]==2) 
          { // move the passband in frequency to the correct frequency 
            if(pbOffsetIdx==0)
            {
              pbOffset = -offset(i,0); 
            }
            if(pbOffsetIdx==1)
            {
              pbOffset = offset(i,0); 
            }
          } 
          if(numPB[i]==4)
          {
            if(pbOffsetIdx==0)
            {
              pbOffset = -offset(i,0)-offset(i,1); 
            }
            if(pbOffsetIdx==1)
            {
              pbOffset = -offset(i,0)+offset(i,1); 
            }
            if(pbOffsetIdx==2)
            {
              pbOffset = offset(i,0)-offset(i,1); 
            }
            if(pbOffsetIdx==3)
            {
              pbOffset = offset(i,0)+offset(i,1); 
            }
          }
          for(Index iii=0;iii<numVal;++iii)
          {
            f_range[pbOffsetIdx*numVal+iii] += 1*pbOffset;
          }
        }
        // Are any passbands overlapping?
        for(Index ii = 2;ii<(f_range.nelem()-2);++ii)
        {
          if(((b_resp.data[ii-1]==1) && (b_resp.data[ii]==0) &&(b_resp.data[ii+1]==0) && (b_resp.data[ii+2]==1)))
          {
            if((f_range[ii]>=f_range[ii+1]))  // Overlapping passbands
            { 
              if((f_range[ii+2]-f_range[ii-1])>verbosityVectIn[i])  // difference in frequency between passbands 
              {
                f_range[ii+1] = f_range[ii+2]-verbosityVectIn[i]/2;
                f_range[ii] = f_range[ii-1]+verbosityVectIn[i]/2;
              } 
              else
              {
                f_range[ii-1] = (f_range[ii]+f_range[ii+2])/2-2*verbosityVectIn[i]/minRatioVerbosityVsFdiff;
                f_range[ii+1] = (f_range[ii]+f_range[ii+2])/2+verbosityVectIn[i]/minRatioVerbosityVsFdiff;
                f_range[ii] = f_range[ii-1]+verbosityVectIn[i]/minRatioVerbosityVsFdiff ;
                f_range[ii+2] = f_range[ii+1]+verbosityVectIn[i]/minRatioVerbosityVsFdiff;
              }
            }
          }
        }
        b_resp.set_grid(0,f_range); 
      }

      // construct sideband response 
      ArrayOfGriddedField1 sideband_response_multi(n);
      for (Index i=0;i<n;++i)
      {
        GriddedField1& r = sideband_response_multi[i];
        r.set_name("Sideband response function");
        r.resize((Index)numPBpseudo[i]);
        Vector f((Index)numPBpseudo[i]);
        if(numPB[i]==1)
        {
          r.data[0]=0.5;
          f[0] = -.0*width[i];
          r.data[1]=0.5;
          f[1] = .0*width[i];
        }
        else if(numPB[i]==2)
        {
          r.data[0]=1/numPB[i];
          r.data[1]=1/numPB[i];
          f[0] = -1*offset(i,0)-0.5*width[i];
          f[1] = +1*offset(i,0)+0.5*width[i];
        }
        else if(numPB[i]==4)
        {
          r.data[0]=1/numPB[i];
          r.data[1]=1/numPB[i];
          r.data[2]=1/numPB[i];
          r.data[3]=1/numPB[i];
          f[0] = -offset(i,0)-offset(i,1)-0.5*width[i];;
          f[1] = -offset(i,0)+offset(i,1)-0.5*width[i];;
          f[2] = +offset(i,0)-offset(i,1)+0.5*width[i];;
          f[3] = +offset(i,0)+offset(i,1)+0.5*width[i];;

        }
        r.set_grid_name(0, "Frequency");
        r.set_grid(0,f);
      }

      sensor_norm =1;
      f_gridFromSensorAMSUgeneric(
          // out
          f_grid,
          // in
          f_backend_multi, 
          backend_channel_response_multi,
          spacing,
          verbosityVect,
          verbosity);

      // do some final work...
      AntennaOff(         // out 
          antenna_dim,
          mblock_za_grid, 
          mblock_aa_grid, 
          // in 
          verbosity);

      sensor_responseInit(
          // out 
          sensor_response, 
          sensor_response_f, sensor_response_pol, 
          sensor_response_za, sensor_response_aa, 
          sensor_response_f_grid, 
          sensor_response_pol_grid, 
          sensor_response_za_grid, 
          sensor_response_aa_grid,
          // in
          f_grid, 
          mblock_za_grid, 
          mblock_aa_grid, 
          antenna_dim,
          atmosphere_dim, 
          stokes_dim,
          sensor_norm,
          verbosity);  

      Index numLO = lo_multi.nelem();
      // Variables for data to be appended
      // Based on code from m_sensor->sensor_responseMultiMixerBackend()
      Array<Sparse> sr;
      ArrayOfVector srfgrid;
      ArrayOfIndex  cumsumf(numLO+1,0);
      const Index nza  = sensor_response_za_grid.nelem();

      // Do this for all channels .... 
      for( Index idxLO = 0; idxLO < numLO ; idxLO++)
      {
        Sparse       sr1      = sensor_response;
        Vector       srf1     = sensor_response_f;
        ArrayOfIndex srpol1   = sensor_response_pol;
        Vector       srza1    = sensor_response_za;
        Vector       sraa1    = sensor_response_aa;
        Vector       srfgrid1 = sensor_response_f_grid;


        sensor_responseBackend( // out
            sr1, 
            srf1,
            srpol1,
            srza1, 
            sraa1, 
            srfgrid1, 

            //in
            sensor_response_pol_grid,
            sensor_response_za_grid,
            sensor_response_aa_grid, 
            f_backend_multi[idxLO],
            backend_channel_response_multi[idxLO],
            sensor_norm,
            verbosity); 

        // Store in temporary arrays
        sr.push_back( sr1 );
        srfgrid.push_back( srfgrid1 );
        //
        cumsumf[idxLO+1] = cumsumf[idxLO] + srfgrid1.nelem();

      }

      // Append data to create sensor_response_f_grid
      //
      const Index  nfnew = cumsumf[numLO];
      sensor_response_f_grid.resize( nfnew );
      //
      for( Index ilo=0; ilo<numLO; ilo++ )
      {
        for( Index i=0; i<srfgrid[ilo].nelem(); i++ )
        {
          sensor_response_f_grid[cumsumf[ilo]+i] = srfgrid[ilo][i];

        }
      }


      const Index  ncols   = sr[0].ncols();
      const Index  npolnew = sensor_response_pol_grid.nelem();
      const Index  nfpolnew = nfnew * npolnew;
      //
      sensor_response.resize( nza*nfpolnew, ncols );
      //
      Vector dummy( ncols, 0.0 );
      //
      for( Index ilo=0; ilo<numLO; ilo++ )
      {
        const Index nfpolthis = (cumsumf[ilo+1]-cumsumf[ilo]) * npolnew;

        assert( sr[ilo].nrows() == nza*nfpolthis );
        assert( sr[ilo].ncols() == ncols );

        for( Index iz=0; iz<nza; iz++ )
        {
          for( Index i=0; i<nfpolthis; i++ )
          {
            // "Poor mans" transfer of a row from one sparse to another 
            for( Index ic=0; ic<ncols; ic++ )
            { dummy[ic] = sr[ilo](iz*nfpolthis+i,ic); }

            sensor_response.insert_row( iz*nfpolnew+cumsumf[ilo]*npolnew+i, 
                dummy );
          }
        }
      }  

      sensor_aux_vectors( // out
          sensor_response_f, 
          sensor_response_pol,
          sensor_response_za,
          sensor_response_aa, 
          sensor_response_f_grid,
          //in
          sensor_response_pol_grid,
          sensor_response_za_grid,
          sensor_response_aa_grid, 
          0);
  }


void sensor_responseSimpleAMSU(// WS Output:
                               Vector& f_grid,
                               Index& antenna_dim,
                               Vector& mblock_za_grid,
                               Vector& mblock_aa_grid,
                               Sparse& sensor_response,
                               Vector& sensor_response_f,
                               ArrayOfIndex& sensor_response_pol,
                               Vector& sensor_response_za,
                               Vector& sensor_response_aa,
                               Vector& sensor_response_f_grid,
                               ArrayOfIndex& sensor_response_pol_grid,
                               Vector& sensor_response_za_grid,
                               Vector& sensor_response_aa_grid,
                               Index& sensor_norm,
                               // WS Input:
                               const Index& atmosphere_dim,
                               const Index& stokes_dim,
                               const Matrix& sensor_description_amsu,
                               // WS Generic Input:
                               const Numeric& spacing,
                               const Verbosity& verbosity)
{
  // Check that sensor_description_amsu has the right dimension:
  if ( 3 != sensor_description_amsu.ncols() )
  {
    ostringstream os;
    os << "Input variable sensor_description_amsu must have three columns, but it has "
       << sensor_description_amsu.ncols() << ".";
    throw runtime_error( os.str() );
  }
  
  // Number of instrument channels:
  const Index n = sensor_description_amsu.nrows(); 

  // The meaning of the columns in sensor_description_amsu is:
  // LO frequency, channel center offset from LO, channel width.
  // (All in Hz.)
  ConstVectorView lo_multi = sensor_description_amsu(Range(joker),0);
  ConstVectorView offset   = sensor_description_amsu(Range(joker),1);
  ConstVectorView width    = sensor_description_amsu(Range(joker),2);

  // Channel frequencies:
  ArrayOfVector f_backend_multi(n);
  for (Index i=0; i<n; ++i) {
    Vector& f = f_backend_multi[i];
    f.resize(1);
    f[0] = lo_multi[i] + offset[i];
  }  
  
  // Construct channel response
  ArrayOfArrayOfGriddedField1 backend_channel_response_multi(n);
  for (Index i=0; i<n; ++i) {
    backend_channel_response_multi[i].resize(1);
    GriddedField1& r = backend_channel_response_multi[i][0];
    r.set_name("Backend channel response function");
    r.resize(2);
    
    // Frequency range:
    Vector f(2);
    f[0] = - 0.5 * width[i];
    f[1] = + 0.5 * width[i];
    r.set_grid_name(0, "Frequency");
    r.set_grid(0,f);
    
    // Response:
    r.data[0] = 1;
    r.data[1] = 1;
  }
  
  // Construct sideband response:
  ArrayOfGriddedField1 sideband_response_multi(n);
  for (Index i=0; i<n; ++i) {
    GriddedField1& r = sideband_response_multi[i];
    r.set_name("Sideband response function");
    r.resize(2);
    
    // Frequency range:
    Vector f(2);
    f[0] = - (offset[i] + 0.5*width[i]);
    f[1] = + (offset[i] + 0.5*width[i]);
    r.set_grid_name(0, "Frequency");
    r.set_grid(0,f);
    
    // Response:
    r.data[0] = 0.5;
    r.data[1] = 0.5;
  }

  // Set sideband mode:
  ArrayOfString sideband_mode_multi(n,"upper");
  
  // We want to automatically normalize the sensor response data, so set sensor_norm to 1:
  sensor_norm = 1;
   
  // Now the rest is just to use some workspace methods:
  // ---------------------------------------------------

  f_gridFromSensorAMSU(f_grid, lo_multi, 
                       f_backend_multi, backend_channel_response_multi, 
                       spacing, verbosity);
  
  AntennaOff(antenna_dim, mblock_za_grid, mblock_aa_grid, verbosity);
  
  sensor_responseInit(sensor_response, 
                      sensor_response_f, sensor_response_pol, 
                      sensor_response_za, sensor_response_aa, 
                      sensor_response_f_grid, 
                      sensor_response_pol_grid, 
                      sensor_response_za_grid, 
                      sensor_response_aa_grid, f_grid, mblock_za_grid, 
                      mblock_aa_grid, antenna_dim, atmosphere_dim, 
                      stokes_dim, sensor_norm,
                      verbosity);  
  
  sensor_responseMultiMixerBackend(sensor_response, sensor_response_f, 
                                   sensor_response_pol, 
                                   sensor_response_za, 
                                   sensor_response_aa, 
                                   sensor_response_f_grid, 
                                   sensor_response_pol_grid, 
                                   sensor_response_za_grid, 
                                   sensor_response_aa_grid, lo_multi, 
                                   sideband_response_multi, 
                                   sideband_mode_multi, 
                                   f_backend_multi, 
                                   backend_channel_response_multi, 
                                   sensor_norm,
                                   verbosity);
  
}

// Declare select functions needed by WMRFSelectChannels:

void Select(// WS Generic Output:
            Vector& needles,
            // WS Generic Input:
            const Vector& haystack,
            const ArrayOfIndex& needleind,
            const Verbosity& verbosity);

template< class T >
void Select(// WS Generic Output:
            Array<T>& needles,
            // WS Generic Input:
            const Array<T>& haystack,
            const ArrayOfIndex& needleind,
            const Verbosity& verbosity);

void Select(// WS Generic Output:
            Sparse& needles,
            // WS Generic Input:
            const Sparse& haystack,
            const ArrayOfIndex& needleind,
            const Verbosity& verbosity);

/* Workspace method: Doxygen documentation will be auto-generated */
void WMRFSelectChannels(// WS Output:
                        Vector& f_grid,
                        Sparse& wmrf_weights,
                        Vector& f_backend,
                        // WS Input:
                        const ArrayOfIndex& wmrf_channels,
                        const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  // For error messages:
  ostringstream os;

  // Some checks of input parameters:

  // wmrf_weights must have same number of rows as f_backend, and same
  // number of columns as f_grid.
  if ( (wmrf_weights.nrows() != f_backend.nelem()) ||
       (wmrf_weights.ncols() != f_grid.nelem()) )
    {
      os << "The WSV *wmrf_weights* must have same number of rows as\n"
         << "*f_backend*, and same number of columns as *f_grid*.\n"
         << "wmrf_weights.nrows() = " << wmrf_weights.nrows() << "\n"
         << "f_backend.nelem()    = " << f_backend.nelem()    << "\n"
         << "wmrf_weights.ncols() = " << wmrf_weights.ncols() << "\n"
         << "f_grid.nelem()       = " << f_grid.nelem();
      throw runtime_error(os.str());
    }

  // wmrf_channels must be strictly increasing (no repetitions).
  chk_if_increasing("wmrf_channels", wmrf_channels);

  // All selected channels must be within the original set of
  // channels.
  if ( min(wmrf_channels)<0 )
    {
      os << "Min(wmrf_channels) must be >= 0, but it is "
         << min(wmrf_channels) << ".";
    }
  if ( max(wmrf_channels)>=f_backend.nelem() )
    {
      os << "Max(wmrf_channels) must be less than the total number of channels.\n"
         << "(We use zero-based indexing!)\n"
         << "The actual value you have is "
         << max(wmrf_channels) << ".";
    }

  if (wmrf_channels.nelem()==f_backend.nelem())
    {
      // No channels to be removed, I can return the original grid.
      out2 << "  Retaining all channels.\n";
    }
  else
    {
      out2 << "  Reducing number of channels from "
           << f_backend.nelem() << " to " << wmrf_channels.nelem() << ".\n";
    }


  // Now the real work starts:

  //  1. Remove unwanted channels from f_backend:  
  Select(f_backend, f_backend, wmrf_channels, verbosity);

  // 2. Remove unwanted channels from wmrf_weights. (We also have to
  // do something about the frequency dimension of wmrf_weights, but
  // we'll do that later.)
  Select(wmrf_weights, wmrf_weights, wmrf_channels, verbosity);

  // 3. Identify, which frequencies are still needed, and which are
  // now obsolete. We store the still needed frequencies in an
  // ArrayOfIndex. 

  // Create f_grid_array. We do not store the frequencies themselves,
  // but the indices of the frequencies to use.
  ArrayOfIndex selection;
  // Make sure that selection does not have to be reallocated along
  // the way. (This is purely to improve performance a bit.)
  selection.reserve(f_grid.nelem());

  // Go through f_grid, and check for each frequency whether it is in
  // the set of WMRF frequencies for any of the channels. 
  assert( wmrf_weights.nrows() == f_backend.nelem() ); 
  assert( wmrf_weights.ncols() == f_grid.nelem() );
  for (Index fi=0; fi<wmrf_weights.ncols(); ++fi)
    {
      Index i;
      for (i=0; i<wmrf_weights.nrows(); ++i)
        {
          if ( wmrf_weights(i,fi) != 0 )
            {
              selection.push_back(fi);
              break;
            }
        }
      if (i==wmrf_weights.nrows())
        {
          out3 << "  The frequency with index " << fi
               << " is not used by any channel.\n";
        }
    }

  if (selection.nelem()==f_grid.nelem())
    {
      // No frequencies were removed, I can return the original grid.
      out2 << "  No unnecessary frequencies, leaving f_grid untouched.\n";
    }
  else if (selection.nelem() == 0)
    {
      throw runtime_error("No frequencies found for any selected channels.\n");
    }
  else
    {
      out2 << "  Reducing number of frequency grid points from "
           << f_grid.nelem() << " to " << selection.nelem() << ".\n";
    }

  // 4. Select the right frequencies in f_grid:
  Select(f_grid, f_grid, selection, verbosity);

  // 5. Select the right frequencies in wmrf_weights. This is a bit
  // tricky, since Select works on the row dimension. So we have to
  // take the transpose.
  Sparse wt(wmrf_weights.ncols(), wmrf_weights.nrows());
  transpose(wt, wmrf_weights);
  Select(wt, wt, selection, verbosity);
  wmrf_weights.resize(wt.ncols(), wt.nrows());
  transpose(wmrf_weights, wt);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseWMRF(// WS Output:
                         Sparse& sensor_response,
                         Vector& sensor_response_f,
                         ArrayOfIndex& sensor_response_pol,
                         Vector& sensor_response_za,
                         Vector& sensor_response_aa,
                         Vector& sensor_response_f_grid,
                         // WS Input:
                         const ArrayOfIndex& sensor_response_pol_grid,
                         const Vector& sensor_response_za_grid,
                         const Vector& sensor_response_aa_grid,
                         const Sparse& wmrf_weights,
                         const Vector& f_backend,
                         const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Some sizes
  const Index nf   = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nza  = sensor_response_za_grid.nelem();
  const Index naa  = sensor_response_aa_grid.nelem();
  const Index nin  = nf * npol * nza;
  // Note that there is no distinction between za and aa grids after the antenna

  // Initialise output stream for runtime errors and a flag for errors
  ostringstream os;
  bool          error_found = false;

  // Check that sensor_response variables are consistent in size
  if( sensor_response_f.nelem() != nin )
  {
    os << "Inconsistency in size between *sensor_response_f* and the sensor\n"
       << "grid variables (sensor_response_f_grid etc.).\n";
    error_found = true;
  }
  if( naa  &&  naa != nza )
  {
    os << "Incorrect size of *sensor_response_aa_grid*.\n";
    error_found = true;
  }
  if( sensor_response.nrows() != nin )
  {
    os << "The sensor block response matrix *sensor_response* does not have\n"
       << "right size compared to the sensor grid variables\n"
       << "(sensor_response_f_grid etc.).\n";
    error_found = true;
  }

  if (nin == 0)
    {
      os << "One of f_grid, pol_grid, za_grid are empty. Sizes are: ("
         << nf << ", " << npol << ", " << nza << ")" << "\n";
      error_found = true;
    }

  // Check number of rows in WMRF weight matrix
  //
  const Index nrw = wmrf_weights.nrows();
  //
  if( nrw != f_backend.nelem() ) 
    {
      os << "The WSV *wmrf_weights* must have as many rows\n"
         << "as *f_backend* has elements.\n"
         << "wmrf_weights.nrows() = " << nrw << "\n"
         << "f_backend.nelem()    = " << f_backend.nelem() << "\n";
      error_found = true;
    }

  // Check number of columns in WMRF weight matrix
  //
  const Index ncw = wmrf_weights.ncols();
  //
  if( ncw != sensor_response_f_grid.nelem() ) 
    {
      os << "The WSV *wmrf_weights* must have as many columns\n"
         << "as *sensor_response_f_grid* has elements.\n"
         << "wmrf_weights.ncols()           = " << ncw << "\n"
         << "sensor_response_f_grid.nelem() = " << sensor_response_f_grid.nelem() << "\n";
      error_found = true;
    }

  // If errors where found throw runtime_error with the collected error
  // message (before error message gets too long).
  if( error_found )
    throw runtime_error(os.str());


  // Ok, now the actual work.

  // Here we need a temporary sparse that is copy of the sensor_response
  // sparse matrix. We need it since the multiplication function can not
  // take the same object as both input and output.
  Sparse htmp = sensor_response;
  sensor_response.resize( wmrf_weights.nrows(), htmp.ncols());
  mult( sensor_response, wmrf_weights, htmp );

  // Some extra output.
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  // Update sensor_response_f_grid
  sensor_response_f_grid = f_backend;

  // Set aux variables
  sensor_aux_vectors( sensor_response_f,       sensor_response_pol, 
                      sensor_response_za,      sensor_response_aa, 
                      sensor_response_f_grid,  sensor_response_pol_grid, 
                      sensor_response_za_grid, sensor_response_aa_grid, 0 );
  
}



/* Workspace method: Doxygen documentation will be auto-generated */
void ySimpleSpectrometer(
         Vector&         y,
         Vector&         y_f,
   const Matrix&         iy,
   const Index&          stokes_dim,
   const Vector&         f_grid,
   const Numeric&        df,
   const Verbosity&      verbosity )
{
  // Some dummy values
  const Index sensor_norm = 1, atmosphere_dim = 1;

  // Init sensor reponse
  //
  Index        antenna_dim;
  Vector       mblock_za_grid, mblock_aa_grid, sensor_response_f;
  Vector       sensor_response_za, sensor_response_aa, sensor_response_f_grid;
  Vector       sensor_response_za_grid, sensor_response_aa_grid;
  ArrayOfIndex sensor_response_pol, sensor_response_pol_grid;
  Sparse       sensor_response;
  //
  AntennaOff( antenna_dim, mblock_za_grid, mblock_aa_grid, verbosity );


  sensor_responseInit( sensor_response, sensor_response_f, 
                  sensor_response_pol, sensor_response_za, sensor_response_aa, 
                  sensor_response_f_grid, sensor_response_pol_grid, 
                  sensor_response_za_grid, sensor_response_aa_grid, f_grid, 
                  mblock_za_grid, mblock_aa_grid, antenna_dim, atmosphere_dim, 
                  stokes_dim, sensor_norm, verbosity );

  // Center position of "channels"
  Vector f_backend;
  linspace( f_backend, f_grid[0]+df/2, last(f_grid), df );

  // Create channel response
  ArrayOfGriddedField1 r;
  backend_channel_responseFlat( r, df, verbosity );

  // New sensor response
  sensor_responseBackend( sensor_response, sensor_response_f, 
                          sensor_response_pol, sensor_response_za,
                          sensor_response_aa, sensor_response_f_grid,
                          sensor_response_pol_grid, sensor_response_za_grid,
                          sensor_response_aa_grid, f_backend, r, sensor_norm,
                          verbosity );

  // Some sizes
  const Index nf = f_grid.nelem();
  const Index n = sensor_response.nrows();

  // Convert iy to a vector
  //
  Vector iyb( nf*stokes_dim );
  //
  for( Index is=0; is<stokes_dim; is++ )
    { iyb[Range(is,nf,stokes_dim)] = iy(joker,is); }

  // y and y_f
  //
  y_f = sensor_response_f;
  y.resize( n );
  mult( y, sensor_response, iyb );
}

