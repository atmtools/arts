/* Copyright (C) 2003-2007
   Mattias Ekström <ekstrom@rss.chalmers.se>

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


/* Workspace method: Doxygen documentation will be auto-generated */
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
                   // Control Parameters:
                   const String&    output )
{
  // This function only supports single mixer spectra, so far. So first
  // check that *lo* only has one element.
  if( lo.nelem()!=1 )
    throw runtime_error(
      "The *ConvertIFToRF* function only supports single mixer setups.");
  Numeric l = lo[0];

  // Check that frequencies are not too high. This might be a floating limit.
  // For this we use the variable f_lim, given in Hz.
  Numeric f_lim = 20e9;
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

  // Check what output option is wanted
  if( output=="lower" ) {
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

  } else if( output=="upper") {
    // Only upper sideband will be output, as above the sizes will not change
    // This time we only need to add *lo* to *sensor_response_f*
    // The measurement vector *y* is untouched.
    sensor_response_f += l;

  } else if( output=="double") {
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
    os << "The ArrayOfMatrix "<<ch_response_name<<" can only contain 1 or "
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

  // Give some output to the user.
  out2 << "  Calculating the backend response using values and grid from "
       << "*" << ch_response_name << "*.\n";

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
                          const Matrix&     filter,
                          // WS Generic Input Names:
                          const String&     filter_name)
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
    os << "The sideband filter response matrix *" << filter_name << "* has not"
       << " been\n correctly initialised. A two column matrix is expected.\n";
    error_found = true;
  }

  // and that it covers the whole sensor_response_f range.
  Numeric df_high = filter(filter.nrows()-1,0)-last(sensor_response_f);
  Numeric df_low = sensor_response_f[0]-filter(0,0);
  if( df_high<0 && df_low<0 )
  {
    os << "The frequency grid of the sideband filter matrix *" << filter_name
       << "*\n must be extended by at least " << -df_low << " Hz in the "
       << "lower\n end and " << -df_high << " Hz in the upper end to cover"
       << "the *sensor_response_f* grid.\n";
    error_found = true;
  }
  else if( df_high<0 )
  {
    os << "The frequency grid of the sideband filter matrix *" << filter_name
       << "*\n must be extended by at least " << -df_high << " Hz in the "
       << "upper\n end to cover the *sensor_response_f* grid.\n";
    error_found = true;
  }
  else if( df_low<0 )
  {
   os << "The frequency grid of the sideband filter matrix *" << filter_name
       << "*\n must be extended by at least " << -df_low << " Hz in the "
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

  // Give some output to the user.
  out2 << "  Calculating the mixer and sideband filter response using values\n"
       << "  and grid from *" << filter_name << "*.\n";

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
     // WS Input:
     const Index&           sensor_response_pol,
     const Vector&          sensor_response_za,
     const Vector&          sensor_response_aa,
     const Vector&          lo,
     const Index&           sensor_norm,
     const Vector&          f_backend,
     const Matrix&          sensor_pol,
     // WS Generic Input:
     const Matrix&          sb_filter,
     const Matrix&          ch_resp,
     // WS Generic Input Names:
     const String&          sb_filter_name,
     const String&          ch_resp_name)
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
          "frequency grid of\n*" << ch_resp_name <<
          "* are outside the current sensor response frequency grid\n. "
          "No weighting can be performed.";
    throw runtime_error(os.str());
  }

  // Check that the sideband filter covers *sensor_response_f*
  if( min(sb_filter(joker,0)) > min(sensor_response_f) ||
      max(sb_filter(joker,0)) < max(sensor_response_f) )
  {
    ostringstream os;
    os << "The sideband filter has to cover the current sensor response "
          "frequency grid.\nThe frequencies in *" << sb_filter_name <<
          "* has to be expanded to cover the frequencies\n"
          "from the previous sensor parts.";
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

