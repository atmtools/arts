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
void sensor_responseAntenna1D(// WS Output:
                              Sparse&           sensor_response,
                              // WS Input:
                              const Vector&     f_grid,
                              const Vector&     mblock_za_grid,
                              const Index&      antenna_dim,
                              const Matrix&     sensor_pol,
                              // WS Generic Input:
                              const Matrix&     srm,
                              // WS Generic Input Names:
                              const String&     srm_name,
                              // Control Parameters:
                              const String&     diagram_type,
                              const Numeric&    f_ref )
{
  // Check that the antenna has the right dimension, this implies that the
  // mblock_aa_grid is empty (as set by AntennaSet1D).
  assert(antenna_dim==1);

  // Check that sensor_response has the right size for multiplication,
  // i.e. at least has been initialised by sensor_responseInit.
  Index n = f_grid.nelem() * mblock_za_grid.nelem() * sensor_pol.nrows();
  if( sensor_response.nrows() != n ) {
    ostringstream os;
    os << "The sensor block response matrix *sensor_response* does not have the\n"
       << "right size. Check that at least *sensor_responseInit* has been run\n"
       << "prior to this method.";
    throw runtime_error(os.str());
  }

  // Check that the 'type' keyword is correct and that the antenna response
  // matrix has the right size in combination with the keyword.
  if (diagram_type=="single" || diagram_type=="scale") {
    if (srm.ncols()!=2) {
      ostringstream os;
      os << "When diagram_type is \"single\" or \"scale\", the antenna diagram\n"
         << "matrix, " << srm_name << ", must have two columns.\n";
      throw runtime_error(os.str());
    }
  } else if (diagram_type=="full") {
    if (srm.ncols()!=f_grid.nelem()+1) {
      ostringstream os;
      os << "When diagram_type is \"full\" the antenna diagram matrix must be\n"
         << "defined for every frequency, i.e. " << srm_name << " must have\n"
         << f_grid.nelem()+1 << " columns.";
      throw runtime_error(os.str());
    }
  } else {
    ostringstream os;
    os << "The keyword diagram_type = " << diagram_type << " is not a valid option.\n"
       << "Only \"single\", \"scale\" and \"full\" are valid.";
    throw runtime_error(os.str());
  }

  // FIXME: If a antenna zenith angle offset from the center angle should be
  // applied, this should be done here before we check if the antenna grid
  // extends outside the mblock_za_grid.

  // Check that the antenna diagram angles does not extend beyond
  // mblock_za_grid min and max. First calculate the difference between
  // the lower and upper bounds.
  // NOTE: this assumes that mblock_za_grid and srm(joker,0) are increasing
  Numeric diff_min = srm(0,0)-mblock_za_grid[0];
  Numeric diff_max = last(mblock_za_grid)-srm(srm.nrows()-1,0);
  if (diff_min<0) {
    ostringstream os;
    os << "The antenna diagram is too wide, it has to be shortened in the\n"
       << "lower end by " << -diff_min << " degrees.";
    throw runtime_error( os.str() );
  } else if (diff_max<0) {
    ostringstream os;
    os << "The antenna diagram is too wide, it has to be shortened in the\n"
       << "upper end by " << -diff_max << " degrees.";
    throw runtime_error( os.str() );
  }

  // Now we need to set up the antenna diagram matrix that will be sent to
  // antenna_transfer_matrix. In the cases of 'single' or 'full' we only
  // need to copy the matrix, but for 'scale' we will create a full matrix
  // by scaling the diagram using the reference frequency.
  Matrix srm_tmp;
  if (diagram_type=="scale") {
    // Resize the matrix to the size of 'full', copy the grid column and
    // then loop through the frequencies and scale the gain.
    srm_tmp.resize(srm.nrows(), f_grid.nelem()+1);
    srm_tmp(joker,0) = srm(joker,0);
    for (Index i=0; i<f_grid.nelem(); i++) {
      scale_antenna_diagram(srm_tmp(joker,i+1),srm,f_ref,f_grid[i]);
    }
  } else {
    // In this case diagram_type is either 'single' or 'full', just copy the
    // matrix.
    srm_tmp = srm;
  }

  // Tell the user what is happening
  out2 << "  Calculating the antenna response using values and grid from *"
       << srm_name << "*.\n";

  // Create the response matrix for the antenna, this matrix will later be
  // multiplied with the original sensor_response matrix.
  Sparse antenna_response(f_grid.nelem(),n);
  antenna_transfer_matrix(antenna_response,mblock_za_grid,srm_tmp,f_grid);

  // It's forbidden to have same matrix as input twice to mult and we
  // want to multiply antenna_response with sensor_response and store the
  // result in sensor_response. So we need to create a temporary copy
  // of sensor_response matrix.
  Sparse sensor_response_tmp(sensor_response);
  sensor_response.resize( f_grid.nelem(), n);
  mult( sensor_response, antenna_response, sensor_response_tmp);

  // Some extra information to the user
  out3 << "  Size of *sensor_response*: " << sensor_response.nrows()
       << "x" << sensor_response.ncols() << "\n";
}

//! sensor_responseBackend
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2003-08-15
*/
void sensor_responseBackend(// WS Output:
                            Sparse&         sensor_response,
                            // WS Input:
                            const Vector&   f_backend,
                            const Vector&   f_mixer,
                            // WS Generic Input:
                            const Matrix&   srm,
                            // WS Generic Input Names:
                            const String&   srm_name,
                            // Control Parameters:
                            const String&   response_type )
{
  // Check that sensor_response has the right size for the multiplication with
  // the spectrometer response.
  if( sensor_response.nrows() != f_mixer.nelem()) {
    ostringstream os;
    os << "The sensor block response matrix *sensor_response* does not have\n"
       << "right size. Either it has not been initialised or some sensor in\n"
       << "front of the backend has not been considered.";
    throw runtime_error( os.str());
  }

  // Check that the backend channel response matrix has the right size in 
  // combination with the keyword 'response_type'.
  if (response_type=="single") {
    if (srm.ncols()!=2) {
      ostringstream os;
      os << "The backend channel response *" << srm_name << "* must have "
         << "two columns.";
      throw runtime_error( os.str() );
    }
  } else if (response_type=="full") {
    if (srm.ncols()!=f_backend.nelem()+1) {
      ostringstream os;
      os << "The backend channel response *" << srm_name << "* must have "
         << f_backend.nelem()+1 << " columns,\n the first containing a "
         << "relative frequency grid and the rest containing\n response "
         << "elements for each backend channel frequency.";
      throw runtime_error(os.str());
    }
  } else {
    ostringstream os;
    os << "The keyword response_type must be either \"single\" or \"full\".";
    throw runtime_error(os.str());
  }

  // Check that the channel frequencies together with the relative frequency
  // grid does not expand outside the frequency grid of the sensor_response
  // matrix (FIXME: Fuzzy explanation, change to sensor_response_f?)
  Numeric diff_low = (f_backend[0]+srm(0,0))-f_mixer[0];
  Numeric diff_high = f_mixer[0]-(last(f_backend)+last(srm(0,joker)));
  if (diff_low<0) {
    ostringstream os;
    os << "There is a " << -diff_low << " Hz overlap in the lower end between the\n"
       << "channel response " << srm_name << " and the frequency grid of\n"
       << "*sensor_reponse*. The frequency grid must be expanded by this amount to\n"
       << "meet the channel response.";
    throw runtime_error(os.str());
  } else if (diff_high<0) {
    ostringstream os;
    os << "There is a " << -diff_high << " Hz overlap in the higher end between the\n"
       << "channel response " << srm_name << " and the frequency grid of\n"
       << "*sensor_reponse*. The frequency grid must be expanded by this amount to\n"
       << "meet the channel response.";
    throw runtime_error(os.str());
  }

  // Give some output to the user.
  out2 << "  Calculating the backend response using values and grid from *"
       << srm_name << "*.\n";

  // Call the function that calculates the sensor transfer matrix.
  Sparse backend_response(f_backend.nelem(),f_mixer.nelem());
  spectrometer_transfer_matrix(backend_response,srm,f_backend,f_mixer);

  // Here we need a temporary sparse that is copy of the sensor_response
  // sparse matrix. We need it since the multiplication function can not
  // take the same object as both input and output.
  Sparse sensor_response_tmp = sensor_response;
  sensor_response.resize(f_backend.nelem(),f_mixer.nelem());
  mult(sensor_response,backend_response,sensor_response_tmp);

  // Some extra output.
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

