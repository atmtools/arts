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

//! antenna_diagramAppendArray
/*!
   See the online help (arts -d FUNCTION_NAME)
   
   \author Mattias Ekström
   \date   2003-08-25
*/
void antenna_diagramAppendArray(// WS Output:
                                ArrayOfArrayOfMatrix&   antenna_diagram,
                                // WS Input:
                                const Matrix&           sensor_pol,
                                // WS Generic Input:
                                const ArrayOfMatrix&    a,
                                // WS Generic Input Names:
                                const String&           a_name)
{
  // Check the size of the array of matrices
  if (a.nelem()==0) {
    ostringstream os;
    os << "Input "<<a_name<<" must at least contain one element.";
    throw runtime_error(os.str());
  }
  if (a.nelem()>sensor_pol.nrows()) {
    ostringstream os;
    os << "Input "<<a_name<<" can not contain more elements than\n"
       << "the number of polarisations given by *sensor_pol*";
    throw runtime_error(os.str());
  }

  // Append the array to antenna_diagram
  antenna_diagram.push_back(a);

  // Output info to user
  out2 << "  Appending "<<a_name<<" to *antenna_diagram*\n";
}


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
         << "identical,but sensor_pos has " << sensor_pos.nrows()
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


//! sensor_responseAntenna1D
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2003-08-14
*/
void sensor_responseAntenna1D(
       // WS Output:
       Sparse&                      sensor_response,
       Vector&                      sensor_response_za,
       // WS Input:
       const Vector&                f_grid,
       const Vector&                mblock_za_grid,
       const Index&                 antenna_dim,
       const Matrix&                sensor_pol,
       // WS Generic Input:
       const ArrayOfArrayOfMatrix&  diag,
       const Vector&                antenna_za,
       // WS Generic Input Names:
       const String&                diag_name,
       const String&                antenna_za_name )
{
  // Check that the antenna has the right dimension, this implies that the
  // mblock_aa_grid is empty (as set by AntennaSet1D).
  assert(antenna_dim==1);

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
  Index n = f_grid.nelem() * mblock_za_grid.nelem() * sensor_pol.nrows();
  if( sensor_response.nrows() != n ) {
    os << "The sensor block response matrix *sensor_response* does not have the\n"
       << "right size. Check that at least *sensor_responseInit* has been run\n"
       << "prior to this method.\n";
    error_found = true;
  }

  // Check that the number of elements in diag equals 1 or the number of
  // elements in antenna_za. That is if the same values are used for all
  // directions (or there exist only on viewing angle) or if each angle
  // has its individual values.
  if (diag.nelem()==0 || antenna_za.nelem()==0) {
    ostringstream os;
    os << "The antenna response array " << diag_name << " and the viewing\n"
       << "angle vector " << antenna_za_name << " must contain at least\n"
       << "one element.\n";
    error_found = true;
  } else if (diag.nelem()==1 && antenna_za.nelem()>1) {
    //FIXME: Give output about the same antenna diagram is used for all
    // viewing directions
  } else if (diag.nelem()==antenna_za.nelem()) {
    //FIXME: Give output that each viewing direction uses individual values
  } else {
    ostringstream os;
    os << "The antenna response array " << diag_name << " does not have the"
       << " right\n size. It should either has one element or as many elements"
       << " as\n the number of viewing angles given by " << antenna_za_name
       << ".\n";
    error_found = true;
  }

  // Check each ArrayOfMatrix in diag, it should contain either one Matrix
  // or one per polarisation
  for (Index i=0; i<diag.nelem(); i++) {
    if (diag[i].nelem()!=1 && diag[i].nelem()!=sensor_pol.nrows()) {
      ostringstream os;
      os << "The number of Matrix in element " << i << " in "
         << diag_name << "\nmust be equal one or the number of "
         << "polarisations.\n";
      error_found = true;
    }
    // Check each Matrix in diag[i], it should contain either on column
    // or one per frequency. Also get the difference between the relative
    // zenith angle grid and mblock_za_grid, store the value if it is
    // lower than previous differences.
    for (Index j=0; j<diag[i].nelem(); j++) {
      if ((diag[i])[j].ncols()!=2 &&
          (diag[i])[j].ncols()!=f_grid.nelem()+1) {
        ostringstream os;
        os << "The number of columns in Matrix " << j << " in array element "
           << i << " in  " << diag_name << "\nmust equal two or the number of "
           << "frequencies plus one.\n";
        error_found = true;
      }
      // Also get the difference between the relative zenith angle grid
      // (modified by the antenna viewing directions) and mblock_za_grid,
      // store the value if it is lower than previous differences.
      za_dlow = min(za_dlow,
        ((diag[i])[j](0,0)+antenna_za[i])-mblock_za_grid[0]);
      za_dhigh = min(za_dhigh,
        last(mblock_za_grid)-(diag[i])[j]((diag[i])[j].nrows()-1,0)
        -antenna_za[i]);
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
  out2 << "  Calculating the antenna response using values and grids from *"
       << diag_name << "*.\n";

  // Create the response matrix for the antenna, this matrix will later be
  // multiplied with the original sensor_response matrix.
  Index nout = f_grid.nelem()*sensor_pol.nrows()*antenna_za.nelem();
  Sparse antenna_response(nout,n);
  antenna_transfer_matrix(antenna_response,mblock_za_grid,diag,f_grid,
    antenna_za,sensor_pol.nrows());

  // It's forbidden to have same matrix as input twice to mult and we
  // want to multiply antenna_response with sensor_response and store the
  // result in sensor_response. So we need to create a temporary copy
  // of sensor_response matrix.
  Sparse sensor_response_tmp(sensor_response);
  sensor_response.resize(nout, n);
  mult( sensor_response, antenna_response, sensor_response_tmp);

  // Some extra information to the user
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";
       
  // Update some descriptive variables
  sensor_response_za = antenna_za;
  out3 << "  *sensor_response_za* set to "<<antenna_za_name<<".\n";
}

//! sensor_responseBackend
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2003-08-15
*/
void sensor_responseBackend(// WS Output:
                            Sparse&               sensor_response,
                            Vector&               sensor_response_f,
                            // WS Input:
                            const Vector&         f_backend,
                            const Matrix&         sensor_pol,
                            const Vector&         sensor_response_za,
                            // WS Generic Input:
                            const ArrayOfMatrix&  ch_response,
                            // WS Generic Input Names:
                            const String&         ch_response_name)
{
  // Initialise a output stream for runtime errors, a flag for errors
  // and counters for difference between sensor_response_f and the
  // relative frequency grid.
  ostringstream os;
  bool error_found = false;
  Numeric f_dlow = 0.0;
  Numeric f_dhigh = 0.0;
  Index n_za_pol = sensor_response_za.nelem()*sensor_pol.nrows();

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
  if (ch_response.nelem()!=1 && ch_response.nelem()==sensor_pol.nrows()) {
    os << "The ArrayOfMatrix "<<ch_response_name<<" can only contain 1 or "
       << sensor_pol.nrows()<<" elements.\n";
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

  // Give some output to the user.
  out2 << "  Calculating the backend response using values and grid from *"
       << ch_response_name << "*.\n";

  // Call the function that calculates the sensor transfer matrix.
  Sparse backend_response(f_backend.nelem()*n_za_pol,
    sensor_response_f.nelem()*n_za_pol);
  spectrometer_transfer_matrix(backend_response,ch_response,
    f_backend,sensor_response_f,sensor_response_za.nelem(),
    sensor_pol.nrows());

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

