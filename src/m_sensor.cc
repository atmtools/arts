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

//! GaussianResponse
/*!
   See the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2003-05-20
*/
void GaussianResponse(// WS Generic Output:
                      Matrix&			r_matrix,
                      // WS Generic Output Names:
                      const String& 	r_matrix_name,
                      // Control Parameters:
                      const Numeric& 	FWHM,
                      const Numeric& 	TotWidth,
                      const Numeric& 	MaxSpacing)
{
  //Calculate new size of matrix
  Index nrows = Index (ceil(TotWidth / MaxSpacing)+1);
  r_matrix.resize(nrows,2);

  out2 << "Setting up a sensor response matrix in *"
       << r_matrix_name << "*\n with gaussian distribution.\n";

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

void sensor_responseAntenna1D(// WS Output:
                   			  Sparse&			sensor_response,
                   			  // WS Input:
                   			  const Vector&		f_grid,
                   			  const Vector&		mblock_za_grid,
							  const Index&		antenna_dim,
         					  // WS Generic Input:
                   			  const Matrix&		srm,
                   			  // WS Generic Input Names:
                   			  const String&		srm_name )
{
  //Check that the antenna has the right dimension
  assert(antenna_dim==1);

  //Check that sensor_response has the right size, i.e. has been initialised
  Index n = f_grid.nelem() * mblock_za_grid.nelem();
  if( sensor_response.ncols() != n ) {
    ostringstream os;
    os << "The sensor block response matrix *sensor_response* has not been\n"
	   << "initialised or some sensor in front of the antenna has not been\n"
       << "considered.";
    throw runtime_error( os.str() );
  }

  out2 << "   Calculating the antenna response using values and grid from *"
       << srm_name << "*.\n";

  //FIXME: Temporary solution, change when Sparse has been added to groups
  Sparse antenna_response( f_grid.nelem(), n);
  antenna_transfer_matrix( antenna_response, mblock_za_grid, srm, f_grid);

  Sparse sensor_response_tmp = sensor_response;
  sensor_response.resize( f_grid.nelem(), n);
  mult( sensor_response, antenna_response, sensor_response_tmp);
}

void sensor_responseBackend(// WS Output:
                 			Sparse&			sensor_response,
                 			// WS Input:
                 			const Vector&	f_grid,
                 			const Vector&	f_grid_backend,
                 			// WS Generic Input:
                 			const Matrix&	srm,
                 			// WS Generic Input Names:
                 			const String&	srm_name)
{
  //Check that sensor_response has the right size, i.e. has been initialised
  if( sensor_response.nrows() != f_grid_backend.nelem() ||
      sensor_response.ncols() != f_grid.nelem()) {
    ostringstream os;
    os << "The sensor block response matrix *sensor_response* has not been\n"
	   << "initialised or some sensor in front of the backend has not been\n"
	   << "considered.";
    throw runtime_error( os.str());
  }

  out2 << "   Calculating the backend response using values and grid from *"
       << srm_name << "*.\n";

  //FIXME: Temporary solution, change when Sparse has been added to groups
  Sparse backend_response( f_grid_backend.nelem(), f_grid.nelem());
  spectrometer_transfer_matrix( backend_response, srm, f_grid_backend, f_grid);

  Sparse sensor_response_tmp = sensor_response;
  sensor_response.resize( f_grid_backend.nelem(), f_grid.nelem());
  mult( sensor_response, backend_response, sensor_response_tmp);
}

void sensor_responseInit(// WS Output:
              			 Sparse&			sensor_response,
              			 // WS Input:
              			 const Vector&		f_grid,
              			 const Vector&		mblock_za_grid,
              			 const Vector&		mblock_aa_grid,
              			 const Index&		antenna_dim )
{
  //FIXME: Check the antenna dimension??
  assert( antenna_dim == 1);

  //Resize matrix and store values
  //FIXME: Keep this as close to Sparse procedure as possible for future switch
  Index n = f_grid.nelem() * mblock_za_grid.nelem() * mblock_aa_grid.nelem();
  sensor_response.resize(n,n);

  for( Index i=0; i<n; i++) {
  	sensor_response(i,i) = 1;
  }
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
