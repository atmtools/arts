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
extern const Index GFIELD1_F_GRID;
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
    Index&   antenna_dim,
   Vector&   mblock_za_grid,
   Vector&   mblock_aa_grid )
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
    Index&   antenna_dim,
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
void sensorOff(
   // WS Output:
         Sparse&   sensor_response,
         Vector&   sensor_response_f,
   ArrayOfIndex&   sensor_response_pol,
         Vector&   sensor_response_za,
         Vector&   sensor_response_aa,
          Index&   antenna_dim,
         Vector&   mblock_za_grid,
         Vector&   mblock_aa_grid,
    const Index&   atmosphere_dim,
    const Index&   stokes_dim,
   const Vector&   f_grid )
{
  // Checks are done in sensor_responseInit.

  AntennaOff( antenna_dim, mblock_za_grid, mblock_aa_grid );

  // Dummy variables
  Index         sensor_norm = 1;
  Vector        sensor_response_f_grid;
  ArrayOfIndex  sensor_response_pol_grid;
  Vector        sensor_response_za_grid;
  Vector        sensor_response_aa_grid;

  sensor_responseInit( sensor_response, sensor_response_f, 
                  sensor_response_pol, sensor_response_za, sensor_response_aa, 
                  sensor_response_f_grid, sensor_response_pol_grid, 
                  sensor_response_za_grid, sensor_response_aa_grid, f_grid, 
                  mblock_za_grid, mblock_aa_grid, antenna_dim, atmosphere_dim, 
                  stokes_dim, sensor_norm );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_responseAntenna(
   // WS Output:
               Sparse&   sensor_response,
               Vector&   sensor_response_f,
         ArrayOfIndex&   sensor_response_pol,
               Vector&   sensor_response_za,
               Vector&   sensor_response_aa,
               Vector&   sensor_response_za_grid,
               Vector&   sensor_response_aa_grid,
   // WS Input:
         const Vector&   sensor_response_f_grid,
   const ArrayOfIndex&   sensor_response_pol_grid,
          const Index&   atmosphere_dim,
          const Index&   antenna_dim,
         const Matrix&   antenna_los,
        const GField4&   antenna_response,
          const Index&   sensor_norm )
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
  if( aresponse_f_grid.nelem() == 1 )
  {
    if( f_dlow > 0  ||  f_dhigh > 0 )
    {
      os << "The frequency grid of *antenna_response has a single value. In \n"
         << "this case, the grid frequency point must be inside the range\n"
         << "of considered frequencies (*f_grid*).\n";
      error_found = true;
    } 
  }
  else
  {
    //
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
  if( antenna_dim == 1 )
    antenna1d_matrix( hantenna, antenna_dim, antenna_los, antenna_response, 
                      sensor_response_za_grid, sensor_response_f_grid, 
                      npol, sensor_norm );
  else
    throw runtime_error( "2D antenna cases are not yet handled." );

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
void sensor_responseBackend(
     // WS Output:
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
   const ArrayOfGField1&   backend_channel_response,
            const Index&   sensor_norm )
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

      if( bchr_f_grid.nelem() != backend_channel_response[i].nelem() )
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
      Numeric f2 = max(sensor_response_f_grid) - 
                                          ( f_backend[i] + last(bchr_f_grid) );
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
                      sensor_response_za_grid, sensor_response_aa_grid );
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
void sensor_responseInit(
   // WS Output:
         Sparse&   sensor_response,
         Vector&   sensor_response_f,
   ArrayOfIndex&   sensor_response_pol,
         Vector&   sensor_response_za,
         Vector&   sensor_response_aa,
         Vector&   sensor_response_f_grid,
   ArrayOfIndex&   sensor_response_pol_grid,
         Vector&   sensor_response_za_grid,
         Vector&   sensor_response_aa_grid,
   // WS Input:
   const Vector&   f_grid,
   const Vector&   mblock_za_grid,
   const Vector&   mblock_aa_grid,
    const Index&   antenna_dim,
    const Index&   atmosphere_dim,
    const Index&   stokes_dim,
    const Index&   sensor_norm )
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
void sensor_responseMixer(
    // WS Output:
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
         const Numeric&   lo,
         const GField1&   sideband_response,
           const Index&   sensor_norm )
{
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
  if( sbresponse_f_grid.nelem() != sideband_response.nelem() )
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
  if( fabs(last(sbresponse_f_grid)+sbresponse_f_grid[0]) > 1e3 )
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
    cout << sensor_response_f_grid[0] << "\n";
    cout << lo << "\n";
    cout << sbresponse_f_grid[0] << "\n";

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
                      sensor_response_za_grid, sensor_response_aa_grid );
}



void sensor_responseMultiMixerBackend(
   // WS Output:
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
                  const Vector&   lo_multi,
          const ArrayOfGField1&   sideband_response_multi,
           const ArrayOfString&   sideband_mode_multi,
           const ArrayOfVector&   f_backend_multi,
   const ArrayOfArrayOfGField1&   backend_channel_response_multi,
                   const Index&   sensor_norm )
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
  Index         ntot = 0, nftot = 0; 

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
          sensor_responseMixer( sr1, srf1, srpol1, srza1, sraa1, srfgrid1,
                                sensor_response_pol_grid,
                                sensor_response_za_grid, 
                                sensor_response_aa_grid,
                                lo_multi[ilo], 
                                sideband_response_multi[ilo], 
                                sensor_norm );

          sensor_responseIF2RF( srf1, srfgrid1,
                                lo_multi[ilo], 
                                sideband_mode_multi[ilo] );

          sensor_responseBackend( sr1, srf1, srpol1, srza1, sraa1, srfgrid1,
                                  sensor_response_pol_grid,
                                  sensor_response_za_grid, 
                                  sensor_response_aa_grid,
                                  f_backend_multi[ilo],
                                  backend_channel_response_multi[ilo],
                                  sensor_norm );
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
      ntot        += sr1.nrows();
      nftot       += srfgrid1.nelem();
    }

  // Append data to create total sensor_response and sensor_response_f_grid
  //
  const Index ncols = sr[0].ncols();
  Index row0 = 0, if0 = 0;
  Vector dummy( ncols, 0.0 );
  //
  sensor_response.resize( ntot, ncols );
  sensor_response_f_grid.resize( nftot );
  //
  for( Index ilo=0; ilo<nlo; ilo++ )
    {
      for( Index row=0; row<sr[ilo].nrows(); row++ )
        {
          // "Poor mans" transfer of a row from one sparse to another 
          for( Index ic=0; ic<ncols; ic++ )
            { dummy[ic] = sr[ilo](row,ic); }

          sensor_response.insert_row( row0+row, dummy );
        }
      row0 += sr[ilo].nrows();
      
      for( Index i=0; i<srfgrid[ilo].nelem(); i++ )
        {
          sensor_response_f_grid[if0+i] = srfgrid[ilo][i];
        }
      if0 += srfgrid[ilo].nelem();
    }  

  // Set aux variables
  sensor_aux_vectors( sensor_response_f,       sensor_response_pol, 
                      sensor_response_za,      sensor_response_aa, 
                      sensor_response_f_grid,  sensor_response_pol_grid, 
                      sensor_response_za_grid, sensor_response_aa_grid );
}




//--- Stuff to be updated ----------------------------------------------------



/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridFromSensor(// WS Output:
                      Vector& f_grid,
                      // WS Input:
                      const Vector& lo,
                      const ArrayOfVector& f_backend,
                      const ArrayOfArrayOfGField1& backend_channel_response,
                      // Control Parameters:
                      const Numeric& spacing)
{
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

  // We build up a total list of absolute frequency ranges for
  // both signal and image sidebands:
  Vector fabs_min(2*n_chan), fabs_max(2*n_chan);
  Index ifabs=0;
  for (Index i=0; i<f_backend.nelem(); ++i)
    for (Index s=0; s<f_backend[i].nelem(); ++s)
      {
        ConstVectorView this_grid = backend_channel_response[i][s].get_numeric_grid(0);
        const Numeric this_f_backend = f_backend[i][s];

        // We need to add a bit of extra margin at both sides,
        // otherwise there is a numerical problem in the sensor
        // WSMs. We take 1% of the requested frequncy grid spacing. 
        const Numeric delta = 0.001 * spacing;

        // Signal sideband:
        fabs_min[ifabs] = this_f_backend + this_grid[0] - delta;
        fabs_max[ifabs] = this_f_backend + this_grid[this_grid.nelem()-1] + delta;
        ++ifabs;
        
        // Image sideband:
        Numeric offset  = this_f_backend - lo[i];
        Numeric f_image = lo[i] - offset;

        fabs_min[ifabs] = f_image + this_grid[0] - delta;                  
        fabs_max[ifabs] = f_image + this_grid[this_grid.nelem()-1] + delta;
        ++ifabs;
      }

  //  cout << "fabs_min: " << fabs_min << "\n";
  //  cout << "fabs_max: " << fabs_max << "\n";

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



/* Workspace method: Doxygen documentation will be auto-generated 
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

  //// Check that *sensor_pol* is not a identity matrix. If so this method is
  //// not just unnecessary but also gives wrong output.
  //  if( sensor_pol.nrows()==sensor_pol.ncols() ) {
  //  bool is_I = true;
  //  for( Index it=0; it<stokes_dim; it++ ) {
  //    for( Index jt=0; jt<stokes_dim; jt++ ) {
  //      if( it==jt && sensor_pol(it,jt)!=1 )
  //        is_I = false;
  //      else if( it!=jt && sensor_pol(it,jt)!=0 )
  //        is_I = false;
  //    }
  //  }
  //  if( is_I ) {
  //    ostringstream os;
  //    os << "The matrix *sensor_pol* is an identity matrix and this method is\n"
  //       << "therfor unnecessary.";
  //    throw runtime_error(os.str());
  //  }
  //}

  //// Check each row of *sensor_pol* so that the first element is 1 and the
  //// sum of the squares of the others also equal 1.
  //bool input_error = false;
  //for( Index it=0; it<sensor_pol.nrows(); it++ ) {
  //  if( sensor_pol(it,1)!=1 )
  //    input_error = true;
  //  Numeric row_sum = 0.0;
  //  for( Index jt=1; jt<sensor_pol.ncols(); jt++ )
  //    row_sum += pow(sensor_pol(it,jt),2.0);
  //  if( row_sum!=1.0 )
  //    input_error = true;
  //}
  //if( input_error ) {
  //  ostringstream os;
  //  os << "The elements in *sensor_pol* are not correct. The first element\n"
  //     << "has to be 1 and the sum of the squares of the following should\n"
  //     << "also be 1.";
  //  throw runtime_error(os.str());
  //}


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
*/



/* Workspace method: Doxygen documentation will be auto-generated 
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
*/
