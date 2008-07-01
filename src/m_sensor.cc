/* Copyright (C) 2003-2008
   Mattias Ekström <ekstrom@rss.chalmers.se>
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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
  \author Mattias Ekström <ekstrom@rss.chalmers.se>
  \date   2003-02-13

  \brief  Workspace functions related to sensor modelling variables.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <string>
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "messages.h"
#include "ppath.h"
#include "special_interp.h"
#include "xml_io.h"
#include "sensor.h"
#include "make_vector.h"
#include "sorting.h"

extern const Numeric PI;
extern const Index GFIELD4_FIELD_NAMES;
extern const Index GFIELD4_F_GRID;
extern const Index GFIELD4_ZA_GRID;
extern const Index GFIELD4_AA_GRID;


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void AntennaOff(
        // WS Output:
        Index&          antenna_dim,
        Vector&         mblock_za_grid,
        Vector&         mblock_aa_grid )
{
  out2 << "  Sets the antenna dimensionality to 1.\n";
  antenna_dim = 1;

  out2 << "  Sets *mblock_za_grid* to have length 1 with value 0.\n";
  mblock_za_grid.resize(1);
  mblock_za_grid[0] = 0;

  out2 << "  Sets *mblock_aa_grid* to be an empty vector.\n";
  mblock_aa_grid.resize(0);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void AntennaSet1D(
        // WS Output:
        Index&    antenna_dim,
        Vector&   mblock_aa_grid )
{
  out2 << "  Sets the antenna dimensionality to 1.\n";
  out3 << "    antenna_dim = 1\n";
  out3 << "    mblock_aa_grid is set to be an empty vector\n";
  antenna_dim = 1;
  mblock_aa_grid.resize(0);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void AntennaSet2D(
        // WS Output:
              Index&   antenna_dim,
        // WS Input:
        const Index&   atmosphere_dim )
{
  if( atmosphere_dim != 3 )
    throw runtime_error("Antenna dimensionality 2 is only allowed when the "
                                          "atmospheric dimensionality is 3." );
  out2 << "  Sets the antenna dimensionality to 1.\n";
  out3 << "    antenna_dim = 2\n";
  antenna_dim = 2;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensorOff_NEW(
        // WS Output:
              Sparse&         sensor_response,
              Vector&         sensor_response_f,
              ArrayOfIndex&   sensor_response_pol,
              Vector&         sensor_response_za,
              Vector&         sensor_response_aa,
              Index&          antenna_dim,
              Vector&         mblock_za_grid,
              Vector&         mblock_aa_grid,
        const Index&          atmosphere_dim,
        const Index&          stokes_dim,
        const Vector&         f_grid )
{
  // Checks are done in sensor_responseInit.

  AntennaOff( antenna_dim, mblock_za_grid, mblock_aa_grid );

  // Dummy variables
  Index         sensor_norm = 1;
  Vector        sensor_response_f_grid;
  ArrayOfIndex  sensor_response_pol_grid;
  Vector        sensor_response_za_grid;
  Vector        sensor_response_aa_grid;

  sensor_responseInit_NEW( sensor_response, sensor_response_f, 
    sensor_response_pol, sensor_response_za, sensor_response_aa, 
    sensor_response_f_grid, sensor_response_pol_grid, 
    sensor_response_za_grid, sensor_response_aa_grid, f_grid, mblock_za_grid, 
    mblock_aa_grid, antenna_dim, atmosphere_dim, stokes_dim, sensor_norm );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseAntenna_NEW(
        // WS Output:
              Sparse&           sensor_response,
              Vector&           sensor_response_f,
              ArrayOfIndex&     sensor_response_pol,
              Vector&           sensor_response_za,
              Vector&           sensor_response_aa,
              Vector&           sensor_response_za_grid,
              Vector&           sensor_response_aa_grid,
        // WS Input:
        const Vector&           sensor_response_f_grid,
        const ArrayOfIndex&     sensor_response_pol_grid,
        const Index&            atmosphere_dim,
        const Index&            antenna_dim,
        const Matrix&           antenna_los,
        const GField4&          aresponse,
        const Index&            sensor_norm )
{
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
  const Index lpolgrid = aresponse.get_string_grid(GFIELD4_FIELD_NAMES).nelem();
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
  ConstVectorView aresponse_f_grid = aresponse.get_numeric_grid(GFIELD4_F_GRID);
  //
  chk_if_increasing( "f_grid of antenna_response", aresponse_f_grid );
  //
  Numeric f_dlow  = 0.0;
  Numeric f_dhigh = 0.0;
  //
  f_dlow  = min(sensor_response_f_grid) - aresponse_f_grid[0];
  f_dhigh = last(aresponse_f_grid) - max(sensor_response_f_grid);
  //
  if( f_dlow < 0 ) 
  {
    os << "The frequency grid of *antenna_response is too narrow. It must\n"
       << "cover all considered frequencies (*f_grid*). The grid needs to be\n"
       << "expanded with "<<-f_dlow<<" Hz in the lower end.\n";
    error_found = true;
  }
  if( f_dhigh < 0 ) 
  {
    os << "The frequency grid of *antenna_response* is too narrow. It must\n"
       << "cover all considered frequencies (*f_grid*). The grid needs to be\n"
       << "expanded with "<<-f_dhigh<<" Hz in the higher end.\n";
    error_found = true;
  }
  

  // Checks of antenna_response za dimension
  //
  ConstVectorView aresponse_za_grid = 
                                    aresponse.get_numeric_grid(GFIELD4_ZA_GRID);
  //
  chk_if_increasing( "za_grid of antenna_response", aresponse_za_grid );
  //
  // Check if the relative grid added to the antena_los za angles
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
  if( f_dhigh < 0 ) 
  {
    os << "The WSV *sensor_response_za_grid* is too narrow. It should be\n"
       << "expanded with "<<-za_dhigh<<" deg in the higher end. This change\n"
       << "should be probably applied to *mblock_za_grid*.\n";
    error_found = true;
  }


  // Checks of antenna_response aa dimension
  //
  ConstVectorView aresponse_aa_grid = 
                                    aresponse.get_numeric_grid(GFIELD4_AA_GRID);
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
    if( f_dhigh < 0 ) 
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



  // Call the core function 
  //
  Sparse hantenna;
  //
  // antenna_matrix_NEW( hantenna

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
                      sensor_response_za_grid, sensor_response_aa_grid );
}





/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseBackend_NEW(
        // WS Output:
              Sparse&       sensor_response,
              Vector&       sensor_response_f,
              ArrayOfIndex& sensor_response_pol,
              Vector&       sensor_response_za,
              Vector&       sensor_response_aa,
              Vector&       sensor_response_f_grid,
        // WS Input:
        const ArrayOfIndex& sensor_response_pol_grid,
        const Vector&       sensor_response_za_grid,
        const Vector&       sensor_response_aa_grid,
        const Vector&       f_backend,
        const Matrix&       backend_channel_response,
        const Index&        sensor_norm )
{
  // Some sizes
  const Index nf   = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nza  = sensor_response_za_grid.nelem();
  const Index naa  = sensor_response_aa_grid.nelem();
  const Index nin  = nf * npol * nza;
  // Note that there is no distinction between za and aa grids after the antenna

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

  // Check number of columns in backend_channel_response and that frequencies
  // are strictly increasing
  if( backend_channel_response.ncols() != 2  &&
      backend_channel_response.ncols() != f_backend.nelem()+1 ) 
    {
      os << "The WSV *backend_channel_response* must have 2 or n+1 columns,\n"
         << "where n is the length of *f_backend*.\n"; 
      error_found = true;
    }
  if( !is_increasing( backend_channel_response(joker,0) ) ) 
    {
      os << "The frequency grid of *backend_channel_response* must be strictly\n"
         << "increasing.\n"; 
      error_found = true;
    }
  // We allow f_backend to be unsorted

  // Check if the relative grid added to the channel frequencies expands
  // outside sensor_response_f_grid.
  //
  Numeric f_dlow  = 0.0;
  Numeric f_dhigh = 0.0;
  //
  f_dlow  = min(f_backend) + backend_channel_response(0,0) -
                                                    min(sensor_response_f_grid);
  f_dhigh = max(sensor_response_f_grid) - ( max(f_backend) +
               backend_channel_response(backend_channel_response.nrows()-1,0) );
  //
  if( f_dlow < 0 ) 
  {
    os << "The WSV *sensor_response_f_grid* is too narrow. It should be\n"
       << "expanded with "<<-f_dlow<<" Hz in the lower end. This change\n"
       << "should be applied to either *f_grid* or the sensor part in\n"
       << "front of *sensor_responseBackend*\n";
    error_found = true;
  }
  if( f_dhigh < 0 ) 
  {
    os << "The WSV *sensor_response_f_grid* is too narrow. It should be\n"
       << "expanded with "<<-f_dhigh<<" Hz in the higher end. This change\n"
       << "should be applied to either *f_grid* or the sensor part in\n"
       << "front of *sensor_responseBackend*\n";
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
  spectrometer_matrix_NEW( hbackend, backend_channel_response,
                           f_backend, sensor_response_f_grid, npol, nza,
                           sensor_norm );

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
                      sensor_response_za_grid, sensor_response_aa_grid );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseInit_NEW(
        // WS Output:
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
        const Vector&         f_grid,
        const Vector&         mblock_za_grid,
        const Vector&         mblock_aa_grid,
        const Index&          antenna_dim,
        const Index&          atmosphere_dim,
        const Index&          stokes_dim,
        const Index&          sensor_norm )
{
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
      if( !is_increasing(mblock_za_grid)  &&  !is_decreasing(mblock_za_grid) )  
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
                      sensor_response_za_grid, sensor_response_aa_grid );

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
void sensor_responseMixer_NEW(
        // WS Output:
              Sparse&         sensor_response,
              Vector&         sensor_response_f,
              ArrayOfIndex&   sensor_response_pol,
              Vector&         sensor_response_za,
              Vector&         sensor_response_aa,
              Vector&         sensor_response_f_grid,
        // WS Input:
        const ArrayOfIndex&   sensor_response_pol_grid,
        const Vector&         sensor_response_za_grid,
        const Vector&         sensor_response_aa_grid,
        const Numeric&        lo,
        const Matrix&         sideband_response,
        const Index&          sensor_norm )
{
  // Some sizes
  const Index nf   = sensor_response_f_grid.nelem();
  const Index npol = sensor_response_pol_grid.nelem();
  const Index nza  = sensor_response_za_grid.nelem();
  const Index naa  = sensor_response_aa_grid.nelem();
  const Index nin  = nf * npol * nza;
  // Note that there is no distinction between za and aa grids after the antenna
  const Index nrp  = sideband_response.nrows();


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

  // Check that the sideband filter matrix has OK size, that frequencies
  // are strictly increasing and that grid end poinst are symmetric around 0
  if( sideband_response.ncols()!= 2  ||  nrp < 2  )
  {
    os << "The sideband filter response matrix has not been\n"
       << "correctly initialised. A two column matrix with at least two rows\n"
       << "is expected.\n";
    error_found = true;
  }
  if( !is_increasing( sideband_response(joker,0) ) ) 
    {
      os << "The frequency grid of *sideband_response* must be strictly\n"
         << "increasing.\n"; 
      error_found = true;
    }
  if( fabs(sideband_response(nrp-1,0)+sideband_response(0,0)) > 1e3 )
    {
      os << "The end points of the *sideband_response* frequency grid must be\n"
         << "symmetrically placed around 0. That is, the grid shall cover a\n"
         << "a range that can be written as [-df,df]. \n";
      error_found = true;      
    }

  // Check that response function does not extend outside sensor_response_f_grid
  Numeric df_high = lo + sideband_response(sideband_response.nrows()-1,0) - 
                                                   last(sensor_response_f_grid);
  Numeric df_low  = sensor_response_f_grid[0] - lo - sideband_response(0,0);
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
  mixer_matrix_NEW( hmixer, f_mixer, sensor_response_f_grid, lo,
                    sideband_response, npol, nza, sensor_norm );

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
                      sensor_response_za_grid, sensor_response_aa_grid );
}












//--- Old stuff --------------------------------------------------------------

/* Workspace method: Doxygen documentation will be auto-generated */
void antenna_diagramAppendArray(// WS Output:
                                ArrayOfArrayOfMatrix&   antenna_diagram,
                                // WS Input:
                                const Matrix&           sensor_pol,
                                // WS Generic Input:
                                const ArrayOfMatrix&    a)
{
  // Check the size of the array of matrices
  if (a.nelem()==0) {
    ostringstream os;
    os << "Input array must at least contain one element.";
    throw runtime_error(os.str());
  }
  if (a.nelem()>sensor_pol.nrows()) {
    ostringstream os;
    os << "Input array can not contain more elements than\n"
       << "the number of polarisations given by *sensor_pol*";
    throw runtime_error(os.str());
  }

  // Append the array to antenna_diagram
  antenna_diagram.push_back(a);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void ConvertIFToRF(
                   // WS Output:
                   Vector&          sensor_response_f,
                   Vector&          y,
                   // WS Input:
                   const Matrix&    sensor_pol,
                   const Vector&    sensor_response_za,
                   const Vector&    sensor_response_aa,
                   const Vector&    lo,
                   const Index&     atmosphere_dim,
                   const Matrix&    sensor_pos,
                   const String&    sideband_mode )
{
  // This function only supports single mixer spectra, so far. So first
  // check that *lo* only has one element.
  if( lo.nelem()!=1 )
    throw runtime_error(
      "The *ConvertIFToRF* function only supports single mixer setups.");

  Numeric l = lo[0];

  // Check that frequencies are not too high. This might be a floating limit.
  // For this we use the variable f_lim, given in Hz.
  Numeric f_lim = 25e9;
  if( min(sensor_response_f) > f_lim )
    throw runtime_error("The frequencies seems to already be given in RF.");

  // Check that *y* has the right size
  Index n_mb = sensor_pos.nrows();
  Index n_az = sensor_response_za.nelem();
  Index n_pol = sensor_pol.nrows();
  Index n_freq = sensor_response_f.nelem();
  if( atmosphere_dim>2 )
    n_az *= sensor_response_aa.nelem();
  if( y.nelem()!=n_freq*n_pol*n_az*n_mb )
    throw runtime_error(
      "The measurement vector *y* does not have the correct size.");

  // Check what sideband_mode option is wanted
  if( sideband_mode=="lower" ) {
    // Only lower sideband will be output, the sizes of the output vectors
    // will not increase, only the values are changed.
    // First we reverse *sensor_response_f* and subtract it from *lo*
    Vector f_out =
      sensor_response_f[Range(sensor_response_f.nelem()-1, joker, -1)];
    f_out *= -1;
    f_out += l;

    // The measurement vector *y* then also has to be reversed for each
    // measurement block, zenith and azimuth angle, and polarisation to 
    // match the frequency vector.
    Vector y_out( y.nelem() );
    for( Index mb=0; mb<n_mb; mb++ ) {
      for( Index az=0; az<n_az; az++ ) {
        Index mb_az = (mb*n_az+az)*n_pol*n_freq;
        for( Index p=0; p<n_pol; p++ ) {
          y_out[Range(mb_az+p,n_freq,n_pol)] =
            y[Range(mb_az+(n_freq-1)*n_pol+p,n_freq,-n_pol)];
        }
      }
    }

    // Copy temporary vectors to output vectors
    sensor_response_f = f_out;
    y = y_out;

  } else if( sideband_mode=="upper") {
    // Only upper sideband will be output, as above the sizes will not change
    // This time we only need to add *lo* to *sensor_response_f*
    // The measurement vector *y* is untouched.
    sensor_response_f += l;

  } else if( sideband_mode=="double") {
    // Both upper and lower sideband will be output. The size of both the
    // frequency and measurement vectors will be doubled.
    // First copy the content of *sensor_response_f* twice to f_out, first
    // reversed and then in normal order
    Vector f_out( n_freq*2 );
    f_out[Range(0,n_freq)] = sensor_response_f[Range(n_freq-1,joker,-1)];
    f_out[Range(n_freq,n_freq)] = sensor_response_f;

    // Subtract f_out from *lo* for the first half, and then add *lo* to
    // f_out for the second half
    f_out[Range(0,n_freq)] *= -1;
    f_out += l;

    // The measurement vector *y* is also unfolded in a similar manner.
    Vector y_out( y.nelem()*2 );
    for( Index mb=0; mb<n_mb; mb++ ) {
      for( Index az=0; az<n_az; az++ ) {
        // Primary band, just moved down the vector
        Index mb_az = (mb*n_az+az)*n_pol*n_freq;
        y_out[Range(2*mb_az+n_az*n_freq*n_pol,n_az*n_freq*n_pol)] =
          y[Range(mb_az,n_az*n_freq*n_pol)];
        // Image band, here we need to resort the vector        
        for( Index p=0; p<n_pol; p++ ) {
          Index mb2_az = (2*mb*n_az+az)*n_pol*n_freq;
          y_out[Range(mb2_az+p,n_freq,n_pol)] =
            y[Range(mb_az+(n_freq-1)*n_pol+p,n_freq,-n_pol)];
        }
      }
    }
    
    // Copy temporary vectors to output vectors
    sensor_response_f = f_out;
    y = y_out;

  } else {
    // If keyword is wrong, throw runtime error
    throw runtime_error(
      "The keyword \"output\" is either wrong or misspelled.");
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridFromSensor(// WS Output:
                      Vector& f_grid,
                      // WS Input:
                      const Vector& lo,
                      const Vector& f_backend,
                      const Matrix& backend_channel_response,
                      // Control Parameters:
                      const Numeric& spacing)
{
  const Index n_chan = lo.nelem();

  // Checks on input quantities:

  // There must be at least one channel:
  if (n_chan < 1)
    {
      ostringstream os;
      os << "There must be at least one channel.\n"
         << "(The vector *lo* must have at least one element.)";
      throw runtime_error(os.str());
    }

  // Does length of f_backend match lo?
  if (f_backend.nelem() != n_chan)
    {
      ostringstream os;
      os << "The vectors *lo* and *f_backend* must have same number of elements.";
      throw runtime_error(os.str());
    }

  // Does backend_channel_response have the right dimension?
  if (backend_channel_response.ncols() != n_chan+1)
    {
      ostringstream os;
      os << "The number of columns in matrix *backend_channel_response*\n"
         << "Must be 1 plus the number of channels. (First column is frequency.)";
      throw runtime_error(os.str());
    }

  // Check that the frequency grid in backend_channel_response is strictly
  // increasing:
  if (!is_increasing(backend_channel_response(joker,0)))
    {
      ostringstream os;
      os << "The frequency grid in backend_channel_response must be "
         << "strictly increasing.";
      throw runtime_error(os.str());
    }

  // Construct image bands:
  Vector f_image(n_chan);
  for (Index i=0; i<n_chan; ++i)
    {
      Numeric offset = f_backend[i] - lo[i];
      f_image[i] = lo[i] - offset;
    }
  out3 << "  Image band nominal frequencies: " << f_image << "\n";

  // Find out the non-zero frequency range for each band:

  Vector f_min(n_chan), f_max(n_chan);
  const Index nf = backend_channel_response.nrows(); // Number of
                                                     // frequency grid points in 
                                                     // backend_channel_response

  for (Index i=0; i<n_chan; ++i)
    {
      // Make sure that not all response values are zero: 
      if (max(backend_channel_response(joker,1+i)) == 0)
        {
          ostringstream os;
          os << "The response for one of the channels seems to be all zero!";
          throw runtime_error(os.str());
        }
          
      // To go through backend_channel_response frequency grid:
      Index imin=0;                
      while (backend_channel_response(imin,1+i)==0) ++imin;

      Index imax=nf-1;                
      while (backend_channel_response(imax,1+i)==0) --imax;
      
      if (imax == imin)
        {
          ostringstream os;
          os << "For one channel there is only a single respons value in"
             << "backend_channel_response, but we need at least two.";
          throw runtime_error(os.str());
        }

      // Actually, f_grid must cover not only the non-zero response
      // frequencies of the channel response matrix, but one point
      // more on both sides, because the response is assumed to vary linearly
      // between the grid points.
      if (imin>0)    --imin;
      if (imax<nf-1) ++imax;

      f_min[i] = backend_channel_response(imin,0);
      f_max[i] = backend_channel_response(imax,0);      
    }  

  // Now we build up a total list of absolute frequency ranges for
  // both signal and image sidebands:
  Vector fabs_min(2*n_chan), fabs_max(2*n_chan);
  Index ifabs=0;
  for (Index i=0; i<n_chan; ++i)
    {
      // Signal sideband:
      fabs_min[ifabs] = f_backend[i] + f_min[i];
      fabs_max[ifabs] = f_backend[i] + f_max[i];
      ++ifabs;

      // Image sideband:
      fabs_min[ifabs] = f_image[i] + f_min[i];
      fabs_max[ifabs] = f_image[i] + f_max[i];
      ++ifabs;
    }

//   cout << "fabs_min: " << fabs_min << "\n";
//   cout << "fabs_max: " << fabs_max << "\n";

  // Check for overlap:
  for (Index i=1; i<fabs_min.nelem(); ++i)
    {
      for (Index s=0; s<i; ++s)
        {
          // We check if either fabs_min[i] or fabs_max[i] are inside
          // the interval of any fabs with a smaller index.
          if (((fabs_min[i]>=fabs_min[s]) && (fabs_min[i]<=fabs_max[s])) ||
              ((fabs_max[i]>=fabs_min[s]) && (fabs_max[i]<=fabs_max[s])) )
            {
              ostringstream os;
              os << "Your instrument bands overlap. This case it not (yet) handled.";
              throw runtime_error(os.str());
            }
        }
    }  

  // Create f_grid_unsorted. This is an array of Numeric, so that we
  // can use the STL push_back function.
  ArrayOfNumeric f_grid_unsorted;
  for (Index i=0; i<fabs_min.nelem(); ++i)
    {
      // Band width:
      const Numeric bw = fabs_max[i] - fabs_min[i];

      // How many grid intervals do I need?
      const Numeric npf = ceil(bw/spacing);

      // How many grid points to store? - Number of grid intervals
      // plus 1.
      const Index   npi = (Index) npf + 1;

      // What is the actual grid spacing inside the band?
      const Numeric gs = bw/npf;

      // Create the grid for this band:
      Vector grid(fabs_min[i], npi, gs);

      out3 << "  Band range " << i << ": " << grid << "\n";

      // Append to f_grid_unsorted:
      f_grid_unsorted.reserve(f_grid_unsorted.nelem()+npi);
      for (Index s=0; s<npi; ++s)
        f_grid_unsorted.push_back(grid[s]);
    }

  // Sort the entire f_grid by increasing frequency:
  f_grid.resize(f_grid_unsorted.nelem());
  ArrayOfIndex si;
  get_sorted_indexes(si, f_grid_unsorted);
  for (Index i=0; i<f_grid_unsorted.nelem(); ++i)
    f_grid[i] = f_grid_unsorted[si[i]];

  //  cout << "Sorted indices: " << si << "\n";

  //   cout << "Created f_grid:\n"
  //        << "  " << f_grid << "\n";

  // That's it, we're done!
}



/* Workspace method: Doxygen documentation will be auto-generated */
void GaussianResponse(// WS Generic Output:
                      Matrix&           r_matrix,
                      // Control Parameters:
                      const Numeric&    FWHM,
                      const Numeric&    TotWidth,
                      const Numeric&    MaxSpacing)
{
  //Calculate new size of matrix
  Index nrows = Index (ceil(TotWidth / MaxSpacing)+1);
  r_matrix.resize(nrows,2);

  //Set up grid column, using temporary vector since nlinspace resizes the
  //input vector.
  Vector tmp(nrows);
  nlinspace(tmp, -TotWidth/2, TotWidth/2, nrows);
  r_matrix(joker,0) = tmp;

  //Calculate standard deviation from Full Width at Half Mean
  Numeric sigma = FWHM / (2.*sqrt(2.*log(2.)));

  //Calculate the normalised gaussian response
  for( Index i=0; i<nrows; i++) {
    r_matrix(i,1) = 1/(sigma*sqrt(2*PI)) *
                    exp(-pow(r_matrix(i,0),(Numeric)2.0)/(2*pow(sigma,(Numeric)2.0)));
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensorOff(
        // WS Output:
              Sparse&   sensor_response,
              Vector&   sensor_response_f,
              Vector&   sensor_response_za,
              Vector&   sensor_response_aa,
              Index&    sensor_response_pol,
              Index&    antenna_dim,
              Vector&   mblock_za_grid,
              Vector&   mblock_aa_grid,
        const Index&    atmosphere_dim,
        const Index&    stokes_dim,
        const Matrix&   sensor_pos,
        const Matrix&   sensor_los,
        const Vector&   f_grid )
{
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );

  if( sensor_los.nrows() != sensor_pos.nrows() )
    {
      ostringstream os;
      os << "The number of rows of sensor_pos and sensor_los must be "
         << "identical,but sensor_pos has " << sensor_pos.nrows()
         << " rows, while sensor_los has " << sensor_los.nrows() << " rows.";
      throw runtime_error( os.str() );
    }

  out2 << "  Sets the antenna dimensionality to 1.\n";
  antenna_dim = 1;

  out2 << "  Sets *mblock_za_grid* to have length 1 with value 0.\n";
  mblock_za_grid.resize(1);
  mblock_za_grid[0] = 0;

  out2 << "  Sets *mblock_aa_grid* to be an empty vector.\n";
  mblock_aa_grid.resize(0);

  // Dummy value
  Index sensor_norm = 1;

  sensor_responseInit( sensor_response, sensor_response_f, sensor_response_za,
    sensor_response_aa, sensor_response_pol, f_grid, mblock_za_grid,
    mblock_aa_grid, antenna_dim, atmosphere_dim, stokes_dim, sensor_norm );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseAntenna1D(
       // WS Output:
       Sparse&                      sensor_response,
       Vector&                      sensor_response_za,
       // WS Input:
       const Vector&                sensor_response_f,
       const Index&                 sensor_response_pol,
       const Vector&                mblock_za_grid,
       const Index&                 antenna_dim,
       const ArrayOfArrayOfMatrix&  diag,
       const Index&                 sensor_norm,
       const Matrix&                antenna_los)
{
  // Check that the antenna has the right dimension, this implies that the
  // mblock_aa_grid is empty (as set by AntennaSet1D).
  if( antenna_dim!=1 ) throw runtime_error( "Antenna dimension must be 1." );
  if( antenna_los.ncols()>1 ) throw runtime_error(
    "Too many columns in *antenna_los*, only zenith angles is considered." );

  // Initialise ostringstream and error flag to collect error messages.
  // This way, if several errors occur, they can all be displayed at the
  // same time. The za_dlow and za_dhigh will be used to store maximum
  // differences between the relative grid in the antenna diagrams and
  // mblock_za_grid,
  ostringstream os;
  bool error_found = false;
  Numeric za_dlow = 0.0;
  Numeric za_dhigh = 0.0;

  // Check that sensor_response has the right size for multiplication,
  // i.e. at least has been initialised by sensor_responseInit.
  Index n = sensor_response_za.nelem()*sensor_response_f.nelem()*
            sensor_response_pol;
  if (sensor_response.nrows() != n ) {
    os << "The sensor block response matrix *sensor_response* does not have the\n"
       << "right size. Check that at least *sensor_responseInit* has been run\n"
       << "prior to this method.\n";
    error_found = true;
  }

  // Check that the number of elements in diag equals 1 or the number of
  // elements in antenna_za. That is if the same values are used for all
  // directions (or there exist only one viewing angle) or if each angle
  // has its individual values.
  if (diag.nelem()==0 || antenna_los.nrows()==0) {
    ostringstream os2;
    os2 << "The antenna response array *antenna_diagram* and the viewing\n"
       << "angle matrix *antenna_los* must contain at least one element.\n";
    error_found = true;
  } else if (diag.nelem()==1 && antenna_los.nrows()>1) {
    // The same antenna diagram is used for all viewing directions
  } else if (diag.nelem()==antenna_los.nrows()) {
    // Each viewing direction uses individual values
  } else {
    ostringstream os2;
    os2 << "The antenna response array *antenna_diagram* does not have the"
       << " right\n size. It should either has one element or as many elements"
       << " as\n the number of viewing angles given by *antenna_los*.\n";
    error_found = true;
  }

  // Check each ArrayOfMatrix in diag, it should contain either one Matrix
  // or one per polarisation
  for (Index i=0; i<diag.nelem(); i++) {
    if (diag[i].nelem()!=1 && diag[i].nelem()!=sensor_response_pol) {
      ostringstream os2;
      os2 << "The number of Matrix in element " << i << " in *antenna_diagram*"
         << "\nmust be equal one or the number of polarisations.\n";
      error_found = true;
    }
    // Check each Matrix in diag[i], it should contain either on column
    // or one per frequency. Also check the difference between the antenna
    // diagram zenith angle grid and mblock_za_grid. This is to make sure
    // that the antenna diagram is covered.
    for (Index j=0; j<diag[i].nelem(); j++) {
      if ((diag[i])[j].ncols()!=2 &&
          (diag[i])[j].ncols()!=sensor_response_f.nelem()+1) {
        ostringstream os2;
        os2 << "The number of columns in Matrix " << j << " in array element "
           << i << " in  *antenna_diagram*\nmust equal two or the number of "
           << "frequencies plus one.\n";
        error_found = true;
      }

      // Get the difference between the relative zenith angle grid
      // (modified by the antenna viewing directions) and mblock_za_grid,
      // store the value if it is lower than previous differences.
      za_dlow = min(za_dlow,
        ((diag[i])[j](0,0)+antenna_los(i,0))-mblock_za_grid[0]);
      za_dhigh = min(za_dhigh,
        last(mblock_za_grid)-(diag[i])[j]((diag[i])[j].nrows()-1,0)
        -antenna_los(antenna_los.nrows()-1,0));
    }
  }

  // Check if za_dlow and za_dhigh are negative. If so the relative zenith
  // angle grid expands outside the mblock_za_grid.
  if (za_dlow<0) {
    os << "The *mblock_za_grid* is too narrow, it has to be expanded in the\n"
       << "lower end by " << -za_dlow << " degree(s).\n";
    error_found = true;
  }
  if (za_dhigh<0) {
    os << "The *mblock_za_grid* is too narrow, it has to be expanded in the\n"
       << "upper end by " << -za_dhigh << " degree(s).\n";
    error_found = true;
  }

  // If errors where found throw runtime_error with the collected error
  // message.
  if (error_found)
    throw runtime_error(os.str());

  // Tell the user what is happening
  out2 << "  Calculating the antenna response.\n";

  // Create the response matrix for the antenna, this matrix will later be
  // multiplied with the original sensor_response matrix.
  Index nout = sensor_response_f.nelem()*sensor_response_pol*antenna_los.nrows();
  Sparse antenna_response( nout, n );
  antenna_matrix( antenna_response, sensor_response_za, diag,
                  sensor_response_f, antenna_los(joker,0), 
                  sensor_response_pol, sensor_norm );

  // It's forbidden to have same matrix as input twice to mult and we
  // want to multiply antenna_response with sensor_response and store the
  // result in sensor_response. So we need to create a temporary copy
  // of sensor_response matrix.
  Sparse sensor_response_tmp(sensor_response);
  sensor_response.resize( nout, sensor_response_tmp.ncols());
  mult( sensor_response, antenna_response, sensor_response_tmp);

  // Some extra information to the user
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  // Update some descriptive variables
  sensor_response_za = antenna_los(joker,0);
  out3 << "  *sensor_response_za* set to *sensor_los* with *antenna_los* "
       << "added to it.\n";
}






/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseBackend(// WS Output:
                            Sparse&               sensor_response,
                            Vector&               sensor_response_f,
                            // WS Input:
                            const Vector&         f_backend,
                            const Index&          sensor_response_pol,
                            const Vector&         sensor_response_za,
                            const Index&          sensor_norm,
                            // WS Generic Input:
                            const ArrayOfMatrix&  ch_response)
{
  // Initialise a output stream for runtime errors, a flag for errors
  // and counters for difference between sensor_response_f and the
  // relative frequency grid.
  ostringstream os;
  bool error_found = false;
  Numeric f_dlow = 0.0;
  Numeric f_dhigh = 0.0;
  Index n_za_pol = sensor_response_za.nelem()*sensor_response_pol;

  // Check that sensor_response has the right size for the multiplication with
  // the spectrometer response.
  if( sensor_response.nrows() != sensor_response_f.nelem()*n_za_pol) {
    os << "The sensor block response matrix *sensor_response* does not have\n"
       << "right size. Either it has not been initialised or some sensor in\n"
       << "front of the backend has not been considered.\n";
    error_found = true;
  }

  // Check that the number of elements in ch_response is equal one or
  // the number of polarisations
  if (ch_response.nelem()!=1 && ch_response.nelem()==sensor_response_pol) {
    os << "The ArrayOfMatrix can only contain 1 or "
       << sensor_response_pol <<" elements.\n";
    error_found = true;
  }

  // Check number of columns in each element of ch_response and also calculate
  // the difference between channel frequencies added with the relative
  // frequency grid and the frequency grid of the sensor_response_f.
  for (Index i=0;i<ch_response.nelem();i++) {
    if (ch_response[i].ncols()!=2 &&
        ch_response[i].ncols()==f_backend.nelem()+1) {
      os << "Matrix number "<<i+1<<" in "<<ch_response<<" must have 2 or "
         << f_backend.nelem()+1<<" columns.";
      error_found = true;
    }
    f_dlow =
      min(f_dlow,(ch_response[i](0,0)+f_backend[0])-sensor_response_f[0]);
    f_dhigh =
      min(f_dhigh,last(sensor_response_f)-
      (ch_response[i](ch_response[i].nrows()-1,0)+last(f_backend)));
  }

  // Check if the relative grid added to the channel frequencies expands
  // outside the sensor_response_f grid.
  if (f_dlow<0) {
    os << "The *sensor_response_f* grid is too narrow. It should be\n"
       << "expanded with "<<-f_dlow<<" Hz in the lower end. This change\n"
       << "should be applied to either *f_grid* or the sensor part in\n"
       << "front of *sensor_responseBackend*\n";
    error_found = true;
  }
  if (f_dhigh<0) {
    os << "The *sensor_response_f* grid is too narrow. It should be\n"
       << "expanded with "<<-f_dhigh<<" Hz in the higher end. This change\n"
       << "should be applied to either *f_grid* or the sensor part in\n"
       << "front of *sensor_responseBackend*\n";
    error_found = true;
  }

  // If errors where found throw runtime_error with the collected error
  // message.
  if (error_found)
    throw runtime_error(os.str());

  // Call the function that calculates the sensor transfer matrix.
  Sparse backend_response(f_backend.nelem()*n_za_pol,
    sensor_response_f.nelem()*n_za_pol);
  spectrometer_matrix(backend_response,ch_response,
    f_backend,sensor_response_f,sensor_response_za.nelem(),
    sensor_response_pol, sensor_norm);

  // Here we need a temporary sparse that is copy of the sensor_response
  // sparse matrix. We need it since the multiplication function can not
  // take the same object as both input and output.
  Sparse sensor_response_tmp = sensor_response;
  sensor_response.resize(f_backend.nelem()*n_za_pol,
    sensor_response_tmp.ncols());
  mult(sensor_response,backend_response,sensor_response_tmp);

  // Some extra output.
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  // Update the sensor_response_f variable
  sensor_response_f = f_backend;
  out3 << "  *sensor_response_f* set to *f_backend*\n";
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseIF2RF(
        // WS Output:
              Vector&         sensor_response_f,
              Vector&         sensor_response_f_grid,
        // WS Input:
        const Numeric&        lo,
        const String&         sideband_mode )
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
void sensor_responseInit(// WS Output:
                               Sparse&      sensor_response,
                               Vector&      sensor_response_f,
                               Vector&      sensor_response_za,
                               Vector&      sensor_response_aa,
                               Index&       sensor_response_pol,
                         // WS Input:
                         const Vector&      f_grid,
                         const Vector&      mblock_za_grid,
                         const Vector&      mblock_aa_grid,
                         const Index&       antenna_dim,
                         const Index&       atmosphere_dim,
                         const Index&       stokes_dim,
                         const Index&       sensor_norm )
{
  // Check input
  chk_if_in_range( "antenna_dim", antenna_dim, 1, 2 );
  if( mblock_za_grid.nelem() == 0 )
    throw runtime_error(
                         "The measurement block zenith angle grid is empty." );
  chk_if_increasing( "mblock_za_grid", mblock_za_grid );
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
        throw runtime_error(
                      "The measurement block azimuthal angle grid is empty." );
      chk_if_increasing( "mblock_aa_grid", mblock_aa_grid );
    }

  if (sensor_norm!=1 && sensor_norm!=0)
    throw runtime_error(
      "The normalisation flag, *sensor_norm*, has to be either 0 or 1." );

  // Set description variables
  sensor_response_f   = f_grid;
  sensor_response_za  = mblock_za_grid;
  sensor_response_aa  = mblock_aa_grid;
  sensor_response_pol = stokes_dim;
  Index n = sensor_response_f.nelem()*sensor_response_za.nelem()
                                                          *sensor_response_pol;
  if ( antenna_dim == 2)
    n *= sensor_response_aa.nelem();

  out2 << "  Initialising *sensor_reponse* as a identity matrix.\n";
  out3 << "  Size of *sensor_response*: " << n << "x" << n << "\n";

  //Set matrix to identity matrix
  sensor_response.make_I(n,n);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseMixer(// WS Output:
                          Sparse&           sensor_response,
                          Vector&           sensor_response_f,
                          Vector&           f_mixer,
                          // WS Input:
                          const Index&      sensor_response_pol,
                          const Vector&     sensor_response_za,
                          const Vector&     lo,
                          const Index&      sensor_norm,
                          // WS Generic Input:
                          const Matrix&     filter)
{
  // Initialise a output stream for runtime errors, a flag for errors
  // and counters for difference between sensor_response_f and the
  // relative frequency grid.
  ostringstream os;
  bool error_found = false;
  Index n_za_pol = sensor_response_za.nelem()*sensor_response_pol;

  // Check that sensor_response has the right size, i.e. has been initialised
  if( sensor_response.nrows() != sensor_response_f.nelem()*n_za_pol)
  {
    os << "The sensor block response matrix *sensor_response* does not have\n"
       << "right size. Either it has not been initialised or some part in\n"
       << "front of the mixer has not been considered.\n";
    error_found = true;
  }

  // Check that the sideband filter matrix has been initialised...
  if( filter.ncols()!=2 )
  {
    os << "The sideband filter response matrix has not been\n"
      << "correctly initialised. A two column matrix is expected.\n";
    error_found = true;
  }

  // and that it covers the whole sensor_response_f range.
  Numeric df_high = filter(filter.nrows()-1,0)-last(sensor_response_f);
  Numeric df_low = sensor_response_f[0]-filter(0,0);
  if( df_high<0 && df_low<0 )
  {
    os << "The frequency grid of the sideband filter matrix\n"
       << "must be extended by at least " << -df_low << " Hz in the "
       << "lower\n end and " << -df_high << " Hz in the upper end to cover"
       << "the *sensor_response_f* grid.\n";
    error_found = true;
  }
  else if( df_high<0 )
  {
    os << "The frequency grid of the sideband filter matrix\n"
       << "must be extended by at least " << -df_high << " Hz in the "
       << "upper\n end to cover the *sensor_response_f* grid.\n";
    error_found = true;
  }
  else if( df_low<0 )
  {
   os << "The frequency grid of the sideband filter matrix\n"
       << "must be extended by at least " << -df_low << " Hz in the "
       << "lower\n end to cover the *sensor_response_f* grid.\n";
    error_found = true;
  }

  // Check that there is only one lo frequency and that it is within the
  // sensor_response_f grid
  if( lo.nelem()!=1 )
  {
    os << "Only one local oscillator frequency should be given when using\n"
       << "the *sensor_responseMixer* method. This method only simulates\n"
       << "one mixer/sideband filter unit.\n";
    error_found = true;
  }
  if( lo[0]<sensor_response_f[0] || lo[0]>last(sensor_response_f) )
  {
    os << "The given local oscillator frequency is outside the sensor\n"
       << "frequency grid. It must be within the *sensor_response_f* grid.\n";
    error_found = true;
  }

  // If errors where found throw runtime_error with the collected error
  // message.
  if (error_found)
    throw runtime_error(os.str());

  //Call to calculating function
  Sparse mixer_response;
  mixer_matrix( mixer_response, f_mixer, sensor_response_f, lo[0],
    filter, sensor_response_pol, sensor_response_za.nelem(), sensor_norm );

  // Here we need a temporary sparse that is copy of the sensor_response
  // sparse matrix. We need it since the multiplication function can not
  // take the same object as both input and output.
  Sparse sensor_response_tmp = sensor_response;
  sensor_response.resize( f_mixer.nelem()*n_za_pol,
    sensor_response_tmp.ncols());
  mult( sensor_response, mixer_response, sensor_response_tmp);

  // Some extra output.
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  // Update the sensor_response_f variable
  sensor_response_f = f_mixer;
  out3 << "  *sensor_response_f* set to *f_mixer*\n";
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseMultiMixerBackend(
     // WS Output:
     Sparse&                sensor_response,
     Vector&                sensor_response_f,
     Vector&                f_mixer,
     Index&                 sensor_response_pol,
     // WS Input:
     const Vector&          sensor_response_za,
     const Vector&          sensor_response_aa,
     const Vector&          lo,
     const Index&           sensor_norm,
     const Vector&          f_backend,
     const Matrix&          sensor_pol,
     const Matrix&          ch_resp,
     // WS Generic Input:
     const Matrix&          sb_filter)
{
  // Check if a *lo* is given for each polarisation (row in *sensor_pol*).
  Vector lo_tmp(sensor_pol.nrows());
  if( lo.nelem()==1 )
  {
    // Here we expand *lo* to fit call to multi_mixer_matrix
    lo_tmp = lo[0];
  }
  else if( lo.nelem()==sensor_pol.nrows() )
  {
    lo_tmp = lo;
  }
  else
  {
    ostringstream os;
    os << "The number of elements in *lo* has to be either one or equal the\n"
       << "number of polarisations given by *sensor_pol*.";
    throw runtime_error(os.str());
  }

  // Check size of *sensor_response*
  Index n_pza = sensor_response_za.nelem()*sensor_response_pol;
  Index n_aa = 1;
  if( sensor_response_aa.nelem()!=0 )
    n_aa = sensor_response_aa.nelem();
  if( sensor_response.nrows()!=n_pza*n_aa*sensor_response_f.nelem() )
  {
    ostringstream os;
    os << "The sensor block response matrix *sensor_response* does not have\n"
       << "the right size. Check that at least *sensor_responseInit* has been\n"
       << "run prior to this method.\n";
    throw runtime_error(os.str());
  }

  // Check that *f_backend* combined with the channel response frequency
  // grid is covered by *sensor_response_f*
  if( max(f_backend)+max(ch_resp(joker,0)) > max(sensor_response_f) ||
      min(f_backend)-min(ch_resp(joker,0)) < min(sensor_response_f) )
  {
    ostringstream os;
    os << "The combination of the backend channel frequencies and the "
      << "frequency grid are outside the current sensor response\n"
      << "frequency grid\n. "
      << "No weighting can be performed.";
    throw runtime_error(os.str());
  }

  // Check that the sideband filter covers *sensor_response_f*
  if( min(sb_filter(joker,0)) > min(sensor_response_f) ||
      max(sb_filter(joker,0)) < max(sensor_response_f) )
  {
    ostringstream os;
    os << "The sideband filter has to cover the current sensor response "
      << "frequency grid.\n"
      << "The frequencies in the sideband filter response matrix has\n"
      << "to be expanded to cover the frequencies from the previous\n"
      << "sensor parts.";
    throw runtime_error(os.str());
  }

  //Call to calculating function
  Sparse mixer_response(n_pza*n_aa, sensor_response_f.nelem()*n_pza*n_aa);
  multi_mixer_matrix( mixer_response, sensor_response_f, f_backend, lo_tmp,
                      sb_filter, ch_resp, sensor_response_za.nelem(), n_aa, 
                      sensor_response_pol, sensor_norm );

  // Here we need a temporary sparse that is copy of the sensor_response
  // sparse matrix. We need it since the multiplication function can not
  // take the same object as both input and output.
  Sparse sensor_response_tmp = sensor_response;
  sensor_response.resize( n_pza*n_aa,
    sensor_response_tmp.ncols());
  mult( sensor_response, mixer_response, sensor_response_tmp);

  // Some extra output.
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  // Update the sensor_response_f variable
  sensor_response_f = f_backend;
  f_mixer = f_backend;
  out3 << "  *sensor_response_f* set to *f_backend*\n";
  // Since this method implies a special use of the sensor polaristion response
  // the sensor_response_pol value is set to one.
  sensor_response_pol = 1;
  out3 << "  *sensor_response_pol* set to one due to use of " 
       << "multi-mixer configuration.\n";
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responsePolarisation(// WS Output:
                                 Sparse&          sensor_response,
                                 Index&           sensor_response_pol,
                                 // WS Input:
                                 const Matrix&    sensor_pol,
                                 const Vector&    sensor_response_za,
                                 const Vector&    sensor_response_aa,
                                 const Vector&    sensor_response_f,
                                 const Index&     stokes_dim )
{
  Index n_aa;
  if( sensor_response_aa.nelem()==0 )
    n_aa = 1;
  else
    n_aa = sensor_response_aa.nelem();
  Index n_f_a = sensor_response_f.nelem()*sensor_response_za.nelem()*n_aa;

  // Check that the initial sensor_response has the right size.
  if( sensor_response.nrows() != sensor_response_pol*n_f_a ) {
    ostringstream os;
    os << "The sensor block response matrix *sensor_response* does not have\n"
       << "the right size. Check that at least *sensor_responseInit* has been\n"
       << "run prior to this method.\n";
    throw runtime_error(os.str());
  }

  // Check that sensor_pol and stokes_dim are consistent??
  if ( sensor_pol.ncols()!=stokes_dim ) {
    ostringstream os;
    os << "The number of columns in *sensor_pol* does not match *stokes_dim*.";
    throw runtime_error(os.str());
  }

  /*
  // Check that *sensor_pol* is not a identity matrix. If so this method is
  // not just unnecessary but also gives wrong output.
  if( sensor_pol.nrows()==sensor_pol.ncols() ) {
    bool is_I = true;
    for( Index it=0; it<stokes_dim; it++ ) {
      for( Index jt=0; jt<stokes_dim; jt++ ) {
        if( it==jt && sensor_pol(it,jt)!=1 )
          is_I = false;
        else if( it!=jt && sensor_pol(it,jt)!=0 )
          is_I = false;
      }
    }
    if( is_I ) {
      ostringstream os;
      os << "The matrix *sensor_pol* is an identity matrix and this method is\n"
         << "therfor unnecessary.";
      throw runtime_error(os.str());
    }
  }

  // Check each row of *sensor_pol* so that the first element is 1 and the
  // sum of the squares of the others also equal 1.
  bool input_error = false;
  for( Index it=0; it<sensor_pol.nrows(); it++ ) {
    if( sensor_pol(it,1)!=1 )
      input_error = true;
    Numeric row_sum = 0.0;
    for( Index jt=1; jt<sensor_pol.ncols(); jt++ )
      row_sum += pow(sensor_pol(it,jt),2.0);
    if( row_sum!=1.0 )
      input_error = true;
  }
  if( input_error ) {
    ostringstream os;
    os << "The elements in *sensor_pol* are not correct. The first element\n"
       << "has to be 1 and the sum of the squares of the following should\n"
       << "also be 1.";
    throw runtime_error(os.str());
  }
  */

  // Output to the user.
  out2 << "  Calculating the polarisation response using *sensor_pol*.\n";

  // Call to calculating function
  Sparse pol_response( sensor_pol.nrows()*n_f_a, stokes_dim*n_f_a );
  polarisation_matrix( pol_response, sensor_pol,
    sensor_response_f.nelem(), sensor_response_za.nelem(), stokes_dim );

  // Multiply with sensor_response
  Sparse tmp = sensor_response;
  sensor_response.resize( sensor_pol.nrows()*n_f_a, tmp.ncols());
  mult( sensor_response, pol_response, tmp);

  // Some extra output.
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  // Update sensor_response variable
  sensor_response_pol = sensor_pol.nrows();
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseRotation(// WS Output:
                             Sparse&        sensor_response,
                             // WS Input:
                             const Vector&  sensor_rot,
                             const Matrix&  antenna_los,
                             const Index&   antenna_dim,
                             const Index&   stokes_dim,
                             const Vector&  sensor_response_f,
                             const Vector&  sensor_response_za)
{
  // Check that at least 3 stokes components are simulated. This since no 2
  // and 3 are weighted together.
  if( stokes_dim<3 ) {
    ostringstream os;
    os << "For a rotating sensor *stokes_dim* has to be at least 3.";
    throw runtime_error(os.str());
  }

  // Check that the antenna dimension and the columns of antenna_los is ok.
  if( antenna_dim!=antenna_los.ncols() ) {
    ostringstream os;
    os << "The antenna line-of-sight is not defined in consistency with the\n"
       << "antenna dimension. The number of columns in *antenna_los* should be\n"
       << "equal to *antenna_dim*";
    throw runtime_error(os.str());
  }

  // Check that the incoming sensor response matrix has the right size. Here
  // we use *stokes_dim* instead of *sensor_response_pol* since this function
  // should be used on the 'raw' stokes components. Also check that a antenna
  // response has been run prior to this method.
  // NOTE: Here we test that sensor_response_za has the same length as
  // antenna_los number of rows. Because for 1D, 2D and 3D cases the zenith
  // angles are allways represented (so far).
  Index n = stokes_dim*antenna_los.nrows()*sensor_response_f.nelem();
  if( sensor_response.nrows()!=n ||
      sensor_response_za.nelem()!=antenna_los.nrows() ) {
    ostringstream os;
    os << "A sensor_response antenna function has to be run prior to\n"
       << "*sensor_responseRotation*.";
    throw runtime_error(os.str());
  }

  // Check the size of *sensor_rot* vs. the number rows of *antenna_los*.
  if( sensor_rot.nelem()!=1 && sensor_rot.nelem()!=antenna_los.nrows() ) {
    ostringstream os;
    os << "The size of *sensor_rot* and number of rows in *antenna_los* has\n"
       << "to be equal.";
    throw runtime_error(os.str());
  }

  // Output to the user.
  out2 << "  Calculating the rotation response with rotation from "
       << "*sensor_rot*.\n";

  // Setup L-matrix, iterate through rotation and insert in sensor_response
  Sparse rot_resp( n, n);
  rotation_matrix(rot_resp, sensor_rot, sensor_response_f.nelem(),
    stokes_dim );

  // Multiply with sensor_response
  Sparse tmp = sensor_response;
  sensor_response.resize( n, tmp.ncols());
  mult( sensor_response, rot_resp, tmp);

  // Some extra output.
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

}
