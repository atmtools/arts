/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                      Stefan Buehler   <sbuehler@uni-bremen.de>
                            
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
  \file   m_physics.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-08-20 

  \brief  Workspace methods of physical character.

  This file includes workspace methods for operations that have some
  connection to basic physics. Example of methods are:  <br>
  1. Setting WSV to hold blackbody radiation. <br>
  2. Conversion to brightness temperature.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts.h"
#include "auto_md.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "physics_funcs.h"

extern const Numeric COSMIC_BG_TEMP;



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/



//! MatrixCBR
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-28
*/
void MatrixCBR(
        // WS Output:
              Matrix&   m,
	// WS Generic Output Names:
        const String&   m_name,
        // WS Input:
        const Index&    stokes_dim,
        // WS Generic Input:
        const Vector&   f,
        // WS Generic Input Names:
        const String&   f_name )
{
  const Index n = f.nelem();

  if( n == 0 )
    throw runtime_error( "The given frequency vector is empty." );

  out2 << "  Setting *" << m_name << "* to hold cosmic background "
       << "radiation.\n";

  m.resize(n,stokes_dim);
  m = 0;

  for( Index i=0; i<n; i++ )
    { m(i,0) = planck( f[i], COSMIC_BG_TEMP ); }
}



//! MatrixPlanck
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-28
*/
void MatrixPlanck(
        // WS Output:
              Matrix&   m,
	// WS Generic Output Names:
        const String&   m_name,
        // WS Input:
        const Index&    stokes_dim,
        // WS Generic Input:
        const Vector&   f,
        // WS Generic Input Names:
        const String&   f_name,
        // Control Parameters:
        const Numeric&  t )
{
  const Index n = f.nelem();

  if( n == 0 )
    throw runtime_error( "The given frequency vector is empty." );

  out2 << "  Setting *" << m_name << "* to hold blackbody radiation for a\n"
       << "temperature of " << t << " K.\n";

  m.resize(n,stokes_dim);
  m = 0;

  for( Index i=0; i<n; i++ )
    { m(i,0) = planck( f[i], t ); }
}



//! MatrixToTbByPlanck
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-08-11
*/
void MatrixToTbByPlanck(
        // WS Generic Output:
              Matrix&   y_out,
        // WS Generic Output Names:
        const String&   y_out_name,
        // WS Generic Input:
        const Matrix&   y_in,
        const Vector&   f_grid,
        // WS Generic Input Names:
        const String&   y_in_name,
        const String&   f_grid_name )
{
  // Some lengths
  const Index nin   = y_in.nrows();
  const Index ncols = y_in.ncols();
  const Index nf    = f_grid.nelem();

  // Check sizes
  if( nf == 0 )
    {
      ostringstream os;
      os << "The vector *" << f_grid << "* is empty.";
      throw runtime_error( os.str() );
    }
  if( nin == 0  ||  ncols == 0 )
    {
      ostringstream os;
      os << "The matrix *" << y_in << "* is empty.";
      throw runtime_error( os.str() );
    }
  if( !is_multiple( nin, nf ) )
    {
      ostringstream os;
      os << "The length of *" << f_grid_name << "* is not an integer multiple "
	 << "of the\nnumber of rows of *" << y_in_name << "*.";
      throw runtime_error( os.str() );
    }

  out2 << "   " << y_out_name << " = inv_of_planck(" << y_in_name << "," 
       << f_grid_name << ")\n" ;

  // Note that y_in and y_out can be the same matrix
  if ( &y_out != &y_in )
    { y_out.resize(nin,ncols); }

  // Nummber of repitions of the frequency values
  const Index nrep = integer_div(nin,nf);
        Index ii, irep, icol;

  for( Index iv=0; iv<nf; iv++ )
    {
      for( irep=0; irep<nrep; irep++ )
	{  
	  ii = irep*nf + iv;

	  for( icol=0; icol<ncols; icol++ )
	    { y_out(ii,icol) = invplanck( y_in(ii,icol), f_grid[iv] ); }
	}
    }
}



//! MatrixToTbByRJ
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-08-11
*/
void MatrixToTbByRJ(
        // WS Generic Output:
              Matrix&   y_out,
        // WS Generic Output Names:
        const String&   y_out_name,
        // WS Generic Input:
        const Matrix&   y_in,
        const Vector&   f_grid,
        // WS Generic Input Names:
        const String&   y_in_name,
        const String&   f_grid_name )
{
  // Some lengths
  const Index nin   = y_in.nrows();
  const Index ncols = y_in.ncols();
  const Index nf    = f_grid.nelem();

  // Check sizes
  if( nf == 0 )
    {
      ostringstream os;
      os << "The vector *" << f_grid << "* is empty.";
      throw runtime_error( os.str() );
    }
  if( nin == 0  ||  ncols == 0 )
    {
      ostringstream os;
      os << "The matrix *" << y_in << "* is empty.";
      throw runtime_error( os.str() );
    }
  if( !is_multiple( nin, nf ) )
    {
      ostringstream os;
      os << "The length of *" << f_grid_name << "* is not an integer multiple "
	 << "of the\nnumber of rows of *" << y_in_name << "*.";
      throw runtime_error( os.str() );
    }

  out2 << "   " << y_out_name << " = inv_of_rj(" << y_in_name << "," 
       << f_grid_name << ")\n" ;

  // Note that y_in and y_out can be the same matrix
  if ( &y_out != &y_in )
    { y_out.resize(nin,ncols); }

  // Nummber of repitions of the frequency values
  const Index nrep = integer_div(nin,nf);
        Index irep, ii, icol;

  // To be as fast as possible, there is a special versions for nrep=1 
  // and ncols=1 
  
  if( nrep == 1  &&  ncols == 1 )
    {
      for( Index iv=0; iv<nf; iv++ )
	{
	  for( irep=0; irep<nrep; irep++ )
	    {  
	      ii = irep*nf + iv;
	      
	      for( icol=0; icol<ncols; icol++ )
		{ y_out(ii,icol) = invrayjean( y_in(ii,icol), f_grid[iv] ); }
	    }
	}
    }

  else
    {
      // Here we try to save time by determining the scaling from radiances
      // to brightness temperature for each frequency, and applying this
      // scaling on each repition for that frequency. Note that relationship
      // between radiance and Tb is linear for a given frequency.

      Numeric scfac;
      Index   irep, ii, icol;

      for( Index iv=0; iv<nf; iv++ )
	{
	  scfac = invrayjean( 1, f_grid[iv] );

	  for( irep=0; irep<nrep; irep++ )
	    {  
	      ii = irep*nf + iv;

	      for( icol=0; icol<ncols; icol++ )
	        { y_out(ii,icol) = scfac * y_in(ii,icol); }
	    }
	}
    }
}



//! Converts i_field(Tensor6) in radiance units to brightness temperature unit
/*! 
 This is used to convert intensity in radiance units to brightness temperature 
units for a Tensor6. It uses the function invplanck from physics_funcs.cc.  
The frequency grid is specific input since inside the cloudbox we have only
one frequency and we need Tensor6 conversions probably only inside the cloudbox.  
 
\param y_out        Output : Tensor6 in Brightness temperatures
\param scat_f_index Input  : frequency index
\param f_grid       Input  : frequency grid
\param y_in         Input  : intensity in radiance units
*/
void Tensor6ToTbByPlanck( // WS Generic Output:
			 Tensor6&   y_out,
			 // WS Generic Output Names:
			 const String&   y_out_name,
			 // WS Specific Input: 
			 const Index& scat_f_index,
			 const Vector&   f_grid,
			 // WS Generic Input:
			 const Tensor6&   y_in,
			 // WS Generic Input Names:
			 const String&   y_in_name)
		
{
  // Some lengths
  const Index nv    = y_in.nvitrines();
  const Index ns    = y_in.nshelves();
  const Index nb    = y_in.nbooks();
  const Index np    = y_in.npages();
  const Index nr    = y_in.nrows();
  const Index nc    = y_in.ncols();
  
  Numeric f = f_grid[scat_f_index];
   
  // Note that y_in and y_out can be the same matrix
  if ( &y_out != &y_in )
    { 
      y_out.resize(nv, ns, nb, np, nr, nc);
    }
  
   
  for( Index iv = 0; iv < nv; ++ iv )
    {
      for( Index is = 0; is < ns; ++ is )
	{
	  for( Index ib = 0; ib < nb; ++ ib )
	    {
	      for( Index ip = 0; ip < np; ++ ip )
		{
		  for( Index ir = 0; ir < nr; ++ ir )
		    {
		      for( Index ic = 0; ic < nc; ++ ic)
			{
			  y_out(iv, is, ib, ip, ir, ic) = 
			   invplanck( y_in (iv, is, ib, ip, ir, ic ), f);  
			}
		    }
		}
	    }
	}
    }
}



//! VectorToTbByPlanck
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-08-11
*/
void VectorToTbByPlanck(
        // WS Generic Output:
              Vector&   y_out,
        // WS Generic Output Names:
        const String&   y_out_name,
        // WS Generic Input:
        const Vector&   y_in,
        const Vector&   f_grid,
        // WS Generic Input Names:
        const String&   y_in_name,
        const String&   f_grid_name )
{
  // Some lengths
  const Index nin = y_in.nelem();
  const Index nf  = f_grid.nelem();

  // Check sizes
  if( nf == 0 )
    {
      ostringstream os;
      os << "The vector *" << f_grid << "* is empty.";
      throw runtime_error( os.str() );
    }
  if( nin == 0 )
    {
      ostringstream os;
      os << "The vector *" << y_in << "* is empty.";
      throw runtime_error( os.str() );
    }
  if( !is_multiple( nin, nf ) )
    {
      ostringstream os;
      os << "The length of *" << f_grid_name << "* is not an integer multiple "
	 << "of the\nnumber of rows of *" << y_in_name << "*.";
      throw runtime_error( os.str() );
    }

  out2 << "   " << y_out_name << " = inv_of_planck(" << y_in_name << "," 
       << f_grid_name << ")\n" ;

  // Note that y_in and y_out can be the same vector
  if ( &y_out != &y_in )
    { y_out.resize(nin); }

  // Nummber of repitions of the frequency values
  const Index nrep = integer_div(nin,nf);
        Index irep, ii;

  for( Index iv=0; iv<nf; iv++ )
    {
      for( irep=0; irep<nrep; irep++ )
	{  
	  ii = irep*nf + iv;
	  y_out[ii] = invplanck( y_in[ii], f_grid[iv] );
	}
    }
}



//! VectorToTbByRJ
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-08-09
*/
void VectorToTbByRJ(
        // WS Generic Output:
              Vector&   y_out,
        // WS Generic Output Names:
        const String&   y_out_name,
        // WS Generic Input:
        const Vector&   y_in,
        const Vector&   f_grid,
        // WS Generic Input Names:
        const String&   y_in_name,
        const String&   f_grid_name )
{
  // Some lengths
  const Index nin = y_in.nelem();
  const Index nf  = f_grid.nelem();

  // Check sizes
  if( nf == 0 )
    {
      ostringstream os;
      os << "The vector *" << f_grid << "* is empty.";
      throw runtime_error( os.str() );
    }
  if( nin == 0 )
    {
      ostringstream os;
      os << "The vector *" << y_in << "* is empty.";
      throw runtime_error( os.str() );
    }
  if( !is_multiple( nin, nf ) )
    {
      ostringstream os;
      os << "The length of *" << f_grid_name << "* is not an integer multiple "
	 << "of the\nnumber of rows of *" << y_in_name << "*.";
      throw runtime_error( os.str() );
    }

  out2 << "   " << y_out_name << " = inv_of_rj(" << y_in_name << "," 
       << f_grid_name << ")\n" ;

  // Note that y_in and y_out can be the same vector
  if ( &y_out != &y_in )
    { y_out.resize(nin); }

  // Nummber of repitions of the frequency values
  const Index nrep = integer_div(nin,nf);

  // To be as fast as possible, there are special versions for nrep=1 
  // and nrep>1.
  
  if( nrep == 1 )
    {
      // Here we just loop the frequencies and call invrayjean
      for( Index iv=0; iv<nf; iv++ )
	{ y_out[iv] = invrayjean( y_in[iv], f_grid[iv] ); }
    }

  else
    {
      // Here we try to save time by determining the scaling from radiances
      // to brightness temperature for each frequency, and applying this
      // scaling on each repition for that frequency. Note that relationship
      // between radiance and Tb is linear for a given frequency.

      Numeric scfac;
      Index   irep, ii;

      for( Index iv=0; iv<nf; iv++ )
	{
	  scfac = invrayjean( 1, f_grid[iv] );

	  for( irep=0; irep<nrep; irep++ )
	    {  
	      ii = irep*nf + iv;
	      y_out[ii] = scfac * y_in[ii];
	    }
	}
    }
}



