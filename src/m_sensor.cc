/* Copyright (C) 2003 Mattias Ekström <ekstrom@rss.chalmers.se>

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

extern const Numeric PI;



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! AntennaSet1D
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-28
*/
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



//! AntennaSet2D
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-28
*/
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



//! GaussianResponse
/*!
   See the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2003-05-20
*/
void GaussianResponse(// WS Generic Output:
                      Matrix&           r_matrix,
                      // WS Generic Output Names:
                      const String&     r_matrix_name,
                      // Control Parameters:
                      const Numeric&    FWHM,
                      const Numeric&    TotWidth,
                      const Numeric&    MaxSpacing)
{
  //Calculate new size of matrix
  Index nrows = Index (ceil(TotWidth / MaxSpacing)+1);
  r_matrix.resize(nrows,2);

  out2 << "  Setting up a sensor response matrix in *"
       << r_matrix_name << "*\n  with gaussian distribution.\n";

  //Set up grid column
  Vector tmp(nrows);
  nlinspace(tmp, -TotWidth/2, TotWidth/2, nrows);
  r_matrix(joker,0) = tmp;

  //Calculate standard deviation from Full Width at Half Mean
  Numeric sigma = FWHM / (2.*sqrt(2.*log(2.)));

  //Calculate the normalised gaussian response
  for( Index i=0; i<nrows; i++) {
    r_matrix(i,1) = 1/(sigma*sqrt(2*PI)) *
                    exp(-pow(r_matrix(i,0),2.0)/(2*pow(sigma,2.0)));
  }
}



//! sensorOff
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2003-07-12
*/
void sensorOff(
        // WS Output:
              Sparse&   sensor_response,
              Vector&   sensor_response_f,
              Vector&   sensor_response_za,
              Vector&   sensor_response_aa,
              Matrix&   sensor_pol,
              Vector&   sensor_rot,
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
         << "identical,\Nbut sensor_pos has " << sensor_pos.nrows() 
         << " rows, while sensor_los has " << sensor_los.nrows() << " rows.";
      throw runtime_error( os.str() );
    }

  out2 << "  Sets *sensor_pol* to be the identity matrix.\n";
  sensor_pol.resize( stokes_dim, stokes_dim );
  sensor_pol = 0;
  for( Index i=0; i<stokes_dim; i++ )
    { sensor_pol(i,i) = 1; }

  out2 << "  Sets *sensor_rot* to a vector of zeros.\n";
  sensor_rot.resize( sensor_pos.nrows() );
  sensor_rot = 0;

  out2 << "  Sets the antenna dimensionality to 1.\n";
  antenna_dim = 1;

  out2 << "  Sets *mblock_za_grid* to have length 1 with value 0.\n";
  mblock_za_grid.resize(1);
  mblock_za_grid[0] = 0;

  out2 << "  Sets *mblock_aa_grid* to be an empty vector.\n";
  mblock_aa_grid.resize(0);

  sensor_responseInit( sensor_response, sensor_response_f, sensor_response_za,
                   sensor_response_aa, f_grid, mblock_za_grid, mblock_aa_grid,
                         antenna_dim, sensor_pol, atmosphere_dim, stokes_dim );
}



void sensor_responseAntenna1D(// WS Output:
                              Sparse&           sensor_response,
                              // WS Input:
                              const Vector&     f_grid,
                              const Vector&     mblock_za_grid,
                              const Index&      antenna_dim,
                              // WS Generic Input:
                              const Matrix&     srm,
                              // WS Generic Input Names:
                              const String&     srm_name )
{
  //Check that the antenna has the right dimension
  assert(antenna_dim==1);

  //Check that sensor_response has the right size, i.e. has been initialised
  Index n = f_grid.nelem() * mblock_za_grid.nelem();
  if( sensor_response.nrows() != n ) {
    ostringstream os;
    os << "The sensor block response matrix *sensor_response* does not have the\n"
       << "right size. The number of rows has to equal the product of the number\n"
       << "of elements of the frequency grid and the measurement block zenith angle\n"
       << "grid. Check that at least *sensor_responseInit* has been run prior to\n"
       << "this method.";
    throw runtime_error( os.str() );
  }


  //Check that the antenna response matrix has been initialised
  if( srm.ncols()!=2 ) {
    ostringstream os;
    os << "The antenna response matrix *" << srm_name << "* has not"
       << " been\n correctly initialised. A two column matrix is expected.\n";
    throw runtime_error( os.str() );
  }

  out2 << "  Calculating the antenna response using values and grid from *"
       << srm_name << "*.\n";

  Sparse antenna_response( f_grid.nelem(), n);
  antenna_transfer_matrix( antenna_response, mblock_za_grid, srm, f_grid);

  //cout << "sensor_response (pre):\n" << sensor_response << "\n";
  //xml_write_to_file ("sensor_response.xml", sensor_response);
  //Sparse sr_t( sensor_response.ncols(), sensor_response.nrows() );
  //transpose( sr_t, sensor_response );

  Sparse sensor_response_tmp(sensor_response);
  sensor_response.resize( f_grid.nelem(), n);
  mult( sensor_response, antenna_response, sensor_response_tmp);

  //cout << "sensor_response_tmp:\n" << sensor_response_tmp << "\n";

  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  //cout << "antenna_response:\n" << antenna_response << "\n";
  //xml_write_to_file ("antenna_response.xml", antenna_response);
  //xml_write_to_file ("sensor_response_tmp.xml", sensor_response_tmp);
}

void sensor_responseBackend(// WS Output:
                            Sparse&         sensor_response,
                            // WS Input:
                            const Vector&   f_backend,
                            const Vector&   f_mixer,
                            // WS Generic Input:
                            const Matrix&   srm,
                            // WS Generic Input Names:
                            const String&   srm_name)
{
  //Check that sensor_response has the right size, i.e. has been initialised
  if( sensor_response.nrows() != f_mixer.nelem()) {
    ostringstream os;
    os << "The sensor block response matrix *sensor_response* has not been\n"
       << "initialised or some sensor in front of the backend has not been\n"
       << "considered.";
    throw runtime_error( os.str());
  }

  //Check that the backend response matrix has been initialised
  if( srm.ncols()!=2 ) {
    ostringstream os;
    os << "The backend response response matrix *" << srm_name << "* has not"
       << " been\n correctly initialised. A two column matrix is expected,\n"
       << "and can be created by *GaussianResponse*.\n";
    throw runtime_error( os.str() );
  }

  out2 << "  Calculating the backend response using values and grid from *"
       << srm_name << "*.\n";

  Sparse backend_response( f_backend.nelem(), f_mixer.nelem());
  spectrometer_transfer_matrix( backend_response, srm, f_backend, f_mixer);

  Sparse sensor_response_tmp = sensor_response;
  sensor_response.resize( f_backend.nelem(), f_mixer.nelem());
  mult( sensor_response, backend_response, sensor_response_tmp);
  
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";
}

void sensor_responseInit(// WS Output:
                               Sparse&      sensor_response,
                               Vector&      sensor_response_f,
                               Vector&      sensor_response_za,
                               Vector&      sensor_response_aa,
                         // WS Input:
                         const Vector&      f_grid,
                         const Vector&      mblock_za_grid,
                         const Vector&      mblock_aa_grid,
                         const Index&       antenna_dim,
                         const Matrix&      sensor_pol,
                         const Index&       atmosphere_dim,
                         const Index&       stokes_dim )
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
  if( sensor_pol.ncols() != stokes_dim )
    throw runtime_error( 
         "The number of columns in *sensor_pol* must be equal *stokes_dim*." );


  Index n = f_grid.nelem() * mblock_za_grid.nelem() * sensor_pol.nrows();
  if ( antenna_dim == 2)
    n *= mblock_aa_grid.nelem();

  out2 << "  Initialising *sensor_reponse* as a identity matrix.\n";
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";

  //Set matrix to identity matrix
  sensor_response.resize(n,n);
  for( Index i=0; i<n; i++) {
    sensor_response.rw(i,i) = 1.0;
  }

  // Set description variables
  sensor_response_f.resize( f_grid.nelem() );
  sensor_response_f = f_grid;
  sensor_response_za.resize( mblock_za_grid.nelem() );
  sensor_response_za = mblock_za_grid;
  sensor_response_aa.resize( mblock_aa_grid.nelem() );
  sensor_response_aa = mblock_aa_grid;
}


void sensor_responseMixer(// WS Output:
                          Sparse&           sensor_response,
                          Vector&           f_mixer,
                          // WS Input:
                          const Vector&     f_grid,
                          // WS Generic Input:
                          const Matrix&     sfrm,
                          // WS Generic Input Names:
                          const String&     sfrm_name,
                          // Control Parameters:
                          const Numeric&    lo,
                          const String&     primary_band)
{
  //Check that sensor_response has the right size, i.e. has been initialised
  if( sensor_response.nrows() != f_grid.nelem()) {
    ostringstream os;
    os << "The sensor block response matrix *sensor_response* has not been\n"
       << "initialised or some sensor in front of the mixer has not been\n"
       << "considered.";
    throw runtime_error( os.str() );
  }

  //Check that the sideband filter matrix has been initialised
  if( sfrm.ncols()!=2 ) {
    ostringstream os;
    os << "The sideband filter response matrix *" << sfrm_name << "* has not"
       << " been\n correctly initialised. A two column matrix is expected.";
    throw runtime_error( os.str() );
  }

  //Determine if primary band is upper or not.
  bool is_upper;
  if( primary_band=="upper" ) {
    is_upper = true;
  } else if ( primary_band=="lower" ) {
    is_upper = false;
  } else {
    ostringstream os;
    os << "The primary band has to be specified by either \"upper\" or \"lower\".\n";
    throw runtime_error( os.str());
  }

  out2 << "  Calculating the mixer and sideband filter response using values\n"
       << "  and grid from *" << sfrm_name << "*.\n";

  //Call to calculating function
  Sparse mixer_response(1,1);
  mixer_transfer_matrix( mixer_response, f_mixer, f_grid, is_upper, lo, sfrm );

  Sparse sensor_response_tmp = sensor_response;
  sensor_response.resize( f_mixer.nelem(), f_grid.nelem());
  mult( sensor_response, mixer_response, sensor_response_tmp);

  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";
}


/*===========================================================================
  === Obsolete functions?
  ===========================================================================*/

//! SensorIntegrationVector
/*!
   See the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2003-02-13
*/
/*
void SensorIntegrationVector(
        // WS Generic Output:
              Vector&   h_out,
        // WS Generic Output Names:
        const String&   h_out_name,
        // WS Generic Input:
        const Vector&   f_values,
        const Vector&   f_grid,
        const Vector&   g_grid,
        // WS Generic Input Names:
        const String&   f_values_name,
        const String&   f_grid_name,
        const String&   g_grid_name )
{
  //Call sensor_integration_vector in sensor.cc
  h_out.resize(g_grid.nelem());
  sensor_integration_vector( h_out, f_values, f_grid, g_grid);

}
*/

//! AntennaTransferMatrix
/*!
   See the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2003-03-06
*/
/*
void AntennaTransferMatrix(// WS Generic Output:
                                 Matrix&    Hb,
                           // WS Generic Output Names:
                           const String&    Hb_name,
                           // WS Generic Input:
                           const Vector&    mblock_za_grid,
                           const Matrix&    a,
                           const Vector&    a_grid,
                           const Vector&    f_grid,
                           // WS Generic Input Names:
                           const String&    mblock_za_grid_name,
                           const String&    a_name,
                           const String&    a_grid_name,
                           const String&    f_grid_name)
{
  //Call antenna_transfer_matrix in sensor.cc
  Hb.resize( f_grid.nelem(), mblock_za_grid.nelem() * f_grid.nelem() );
  antenna_transfer_matrix( Hb, mblock_za_grid, a, a_grid, f_grid );

}
*/

//! AntennaTest
/*!
   See the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date  2003-03-10
*/
/*
void AntennaTest(// WS Generic Output:
                 Vector&          a,
                 // WS Generic Output Names:
                 const String&    a_name)
{
  //Set variables
  Vector a_grid(181);
  nlinspace(a_grid,-90,90,181);
  Numeric theta = 3;

  //Set size of a
  a.resize( a_grid.nelem() );

  //Call function
  antenna_diagram_gaussian(a, a_grid, theta);

  //Set up new vector and scale antenna diagram
  Vector a_new = scale_antenna_diagram(a, 1, 100);

  a = a_new;
}
*/
