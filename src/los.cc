/*
Outline of interpolation functions:

There are two types of interpolations that we want to do:
1. Interpolation of some field at a set of points, for example, along the PPATH.
2. To interpolate an input field to the calculation grids used, for example,
   to determine the temperatures at the pressure and latitude grid selected
   from a climatology. The interpolation is here done for all grid crossings.

To distinguish clearly these different cases, the interpolation of the kind in
case 2 will be denoted as re-sampling, and function names for case 1 will start
with "interp", while for case the names will start with "resample".

Both interp and resample functions will take index positions as input
(instead of physical positions). For example, the index position 6.5 means that
a position is exactly between the points with index 6 and 7, index position 5 
is exactly at the point with index 5 etc. [! A more stringent description
should be added. !]

Only linear interpolation will be implemented.

No resizing of vectors, matrices and tensors will be made inside the
interpolation functions. The functions assert that the arguments have 
consistent sizes.

Details for interp-functions (case 1):
The interpolation positions are given as a number of vectors, one for each
dimension that can exist (for example, pressure/altitude, latitude and
longitude). The length of these vectors must be equal, with the exception
that an empty vector (length zero) means that that dimension is not specified.
   For cases there the interpolation is performed along all dimensions, general
functions can be made, and these are called interp_1D, interp_2D, interp_3D.
The output of these functions are throughout a vector, where the length equals
the number of positions for which interpolation is performed. For example,
the function interp_3D can be used to interpolate the temperature field
(stored as a Tensor3) to the PPATH. 
  On the other hand, there are cases when the interpolation is not performed 
along all dimensions. A typical example is interpolation of the absorption
tensor to get the absorption at the points along the PPATH, where no 
interpolation is made in the frequency dimension. For such interpolation cases
special functions will be made (e.g. interp_abs2ppath). To make general 
functions would be too messy and inefficient.
  Both for the general and special cases, the interpolation functions adapt
automatically to the present dimensionality of the simulations. For example, 
if the simulations are 2D, indicated by that the vector with longitude 
positions is empty and the number of pages of the temperature tensor is 1
(consistency between these two criteria is checked), the longitude dimension is
ignored during the interpolation.

Details for resample-functions (case 2):
Resampling corresponds to an interpolation over all involved dimensions. The
grids are given as individual vectors. An empty vector indicates (as above)
that the corresponding dimension is not specified. Interpolation is performed
for all possible combinations between the vectors, that is, all grid crossings
(as expected). The functions adapt automatically to the dimensionality of the
simulations (as above).
*/


/*
Names and symbols:

pressure        : p
altitude        : z
latitude        : alpha
longitude       : beta
zenith angle    : psi
azimuthal angle : omega

*/

#include <math.h>
#include "arts.h"
#include "matpackI.h"
#include "matpackIII.h"
#include "messages.h"          
#include "auto_wsv.h"          

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric EARTH_RADIUS;


struct Ppath {
  Index     dim;
  Index     np;
  Index     i_start;
  Index     i_stop;
  Matrix    pos;
  Matrix    ip_pos;
  Vector    z;
  Vector    l_step;
  Matrix    los;
  String    background;
  Index     ground;
  Index     i_ground;
  Vector    tan_pos;
};



Index is_bool( 
        const Index&    x )
{
  return ( x==0 || x==1 );
}

Index is_sorted( 
        ConstVectorView&    x )
{
  if( x.nelem() > 1 )
    {
      for( Index i=1; i<x.nelem(); i++ )
	{
	  if( x[i] < x[i-1] )
	    return 0;
	}
    }
  return 1;
}

Index is_increasing( 
        ConstVectorView&    x )
{
  if( x.nelem() > 1 )
    {
      for( Index i=1; i<x.nelem(); i++ )
	{
	  if( x[i] <= x[i-1] )
	    return 0;
	}
    }
  return 1;
}


void chk_if_bool( 
        const String&   x_name,
        const Index&    x )
{
  if ( !is_bool(x) )
    {
      ostringstream os;
      os << "The variable *" << x_name <<  "* must be a boolean (0 or 1).\n" 
	 << "The present value of *"<< x_name <<  "* is " << x << ".";
      throw runtime_error( os.str() );
    }
}

void chk_if_in_range( 
	const String&   x_name,
        const Index&    x, 
        const Index&    x_low, 
        const Index&    x_high )
{
  if ( (x<x_low) || (x>x_high) )
    {
      ostringstream os;
      os << "The variable *" << x_name <<  "* must fulfill:\n"
	 << "   " << x_low << " <= " << x_name << " <= " << x_high << "\n" 
	 << "The present value of *"<< x_name <<  "* is " << x << ".";
      throw runtime_error( os.str() );
    }
}



void chk_if_over_0( 
	const String&    x_name,
        const Numeric&   x ) 
{
  if ( x <= 0 )
    {
      ostringstream os;
      os << "The variable *" << x_name <<  "* must exceed 0.\n"
	 << "The present value of *"<< x_name <<  "* is " << x << ".";
      throw runtime_error( os.str() );
    }
}

void chk_if_in_range( 
	const String&    x_name,
        const Numeric&   x, 
        const Numeric&   x_low, 
        const Numeric&   x_high )
{
  if ( (x<x_low) || (x>x_high) )
    {
      ostringstream os;
      os << "The variable *" << x_name <<  "* must fulfill:\n"
	 << "   " << x_low << " <= " << x_name << " <= " << x_high << "\n" 
	 << "The present value of *"<< x_name <<  "* is " << x << ".";
      throw runtime_error( os.str() );
    }
}



void chk_vector_length( 
	const String&      x_name,
        ConstVectorView&   x,
        const Index&       l ) 
{
  if ( x.nelem() != l )
    {
      ostringstream os;
      os << "The vector *" << x_name <<  "* must have the length " << l 
         << ".\n" 
         << "The present length of *"<< x_name <<  "* is " << x.nelem() << ".";
      throw runtime_error( os.str() );
    }
}

void chk_vector_length( 
	const String&      x1_name,
	const String&      x2_name,
        ConstVectorView&   x1, 
        ConstVectorView&   x2 ) 
{
  if ( x1.nelem() != x2.nelem() )
    {
      ostringstream os;
      os << "The vectors *" << x1_name <<  "* and *" << x1_name 
         <<  "* must have the same length.\n"
         << "The length of *"<< x1_name <<  "* is " << x1.nelem() << ".\n"
         << "The length of *"<< x2_name <<  "* is " << x2.nelem() << ".";
      throw runtime_error( os.str() );
    }
}

void chk_if_increasing( 
	const String&      x_name,
        ConstVectorView&   x ) 
{
  if ( !is_increasing(x) )
    {
      ostringstream os;
      os << "The vector *" << x_name <<  "* must have strictly\nincreasing "
         << "values, but this is not the case.";
      throw runtime_error( os.str() );
    }
}



void chk_atm_grids( 
	const Index&      dim,
	ConstVectorView   p_grid,
	ConstVectorView   alpha_grid,
	ConstVectorView   beta_grid )
{
  if( p_grid.nelem() < 2 )
    throw runtime_error( "The length of *p_grid* must be >= 2." );
  chk_if_increasing( "p_grid", p_grid );

  if( dim == 1 )
    {
      if( alpha_grid.nelem() != 0 )
	throw runtime_error(
                          "For dim=1, the length of *alpha_grid* must be 0." );
    }
  else
    {
      if( alpha_grid.nelem() < 2 )
	throw runtime_error(
                         "For dim>1, the length of *alpha_grid* must be >=2.");
      chk_if_increasing( "alpha_grid", alpha_grid );
    }

  if( dim < 3 )
    { 
      if( beta_grid.nelem() != 0 )
	throw runtime_error(
                           "For dim<3, the length of *beta_grid* must be 0." );
    }
  else
    {
      if( beta_grid.nelem() < 2 )
	throw runtime_error(
                          "For dim=3, the length of *beta_grid* must be >=2.");
      chk_if_increasing( "beta_grid", beta_grid );
    }
}

void chk_atm_field( 
	const String&     x_name,
        const Tensor3&    x, 
	const Index&      dim,
	ConstVectorView   p_grid,
	ConstVectorView   alpha_grid,
	ConstVectorView   beta_grid )
{
  Index ncols=p_grid.nelem(), nrows=1, npages=1;
  if( dim > 1 )
    nrows = alpha_grid.nelem();
  if( dim > 2 )
    npages = beta_grid.nelem();
  if( x.ncols()!=ncols || x.nrows()!=nrows || x.npages()!=npages ) 
    {
      ostringstream os;
      os << "The atmospheric field *" << x_name <<  "* has wrong size.\n"
         << "Expected size is " << npages << " x " << nrows << " x " 
         << ncols << ",while actual size is " << x.npages() << " x " 
         << x.nrows() << " x " << x.ncols() << ".";
      throw runtime_error( os.str() );
    }
}

void chk_atm_surface( 
	const String&     x_name,
        const Matrix&     x, 
	const Index&      dim,
	ConstVectorView   alpha_grid,
	ConstVectorView   beta_grid )
{
  Index ncols=1, nrows=1;
  if( dim > 1 )
    ncols = alpha_grid.nelem();
  if( dim > 2 )
    nrows = beta_grid.nelem();
  if( x.ncols()!=ncols || x.nrows()!=nrows ) 
    {
      ostringstream os;
      os << "The atmospheric surface *" << x_name <<  "* has wrong size.\n"
         << "Expected size is " << nrows << " x " << ncols << ","
         << "while actual size is " << x.nrows() << " x " << x.ncols() << ".";
      throw runtime_error( os.str() );
    }
}


/** Returns the first value of a vector.

    \return       The first value of x.
    \param    x   A vector.

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric first2( ConstVectorView x )
{
  assert( x.nelem() > 0 );
  return x[0]; 
}



/** Returns the last value of a vector.

    \return      The last value of x.
    \param   x   A vector.

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric last2( ConstVectorView x )
{
  assert( x.nelem() > 0 );
  return x[x.nelem()-1]; 
}



void truncate_vector( 
	      Vector&   x, 
	const Index&    n )        
{
  assert( n >= 0 );
  if( n > 0 ) 
    {
      Vector   dummy(n);	  
      dummy = x[ Range(0,n) ];
      x.resize(n); 
      x = dummy;
    }
  else
    x.resize(0);
}

void truncate_vector( 
	      Vector&   x, 
	const Index&    i1,
	const Index&    i2 )        
{
  assert( i1 >= 0 );
  assert( i2 >= 0 );
  if( i2 >= i1 ) 
    {
      Index    n = i2-i1+1;
      Vector   dummy(n);	  
      dummy = x[ Range(i1,n) ];
      x.resize(n); 
      x = dummy;
    }
  else
    x.resize(0);
}



void assert_size( 
	ConstVectorView   x,
	const Index&      l ) 
{
  assert( x.nelem() == l );
}



void assert_size( 
	ConstMatrixView   x,
	const Index&      nrows,
	const Index&      ncols ) 
{
  assert( x.ncols() == ncols );
  assert( x.nrows() == nrows );
}



void assert_size( 
	const Tensor3&    x,
	const Index&      npages,
	const Index&      nrows,
	const Index&      ncols ) 
{
  assert( x.ncols() == ncols );
  assert( x.nrows() == nrows );
  assert( x.npages() == npages );
}



void assert_maxdim_of_tensor(
	const Tensor3&   x,
	const Index&     dim )
{
  assert( dim >= 1 );
  assert( dim <= 3 );

  if( dim == 1 )
    {
      assert( x.nrows() == 1 );
      assert( x.npages() == 1 );
    }
  else if( dim == 2 )
    assert( x.npages() == 1 );
}



void get_interp_weights (
	      ArrayOfIndex&   ii,
	      Vector&         w,
	ConstVectorView       y,
	ConstVectorView       ip_x )
{
  // Sizes
  const Index   n_y   = y.nelem();
  const Index   n_out = ip_x.nelem();

  // Asserts
  assert( ii.nelem() == n_out );
  assert( w.nelem() == n_out );

  for( Index ix=0; ix<n_out; ix++ )
    {
      ii[ix] = Index( floor( ip_x[ix] ) );
      assert( ii[ix] >= 0 );
      w[ix]  = ip_x[ix] - Numeric( ii[ix] );
      
      // If w is very small (< 1e-6), treat it to be 0.
      // If w=0, we can be at the end point.
      if( w[ix] < 1e-6 ) 
	{
	  assert( ii[ix] < n_y );
	  w[ix] = 0;
	}
      else
	  assert( ii[ix] < n_y-1 );
    }
}	     



Index get_position_dim (
	ConstVectorView   ip_p,
	ConstVectorView   ip_alpha,
	ConstVectorView   ip_beta )
{
  if( ip_alpha.nelem() == 0 )
    {
      assert( ip_beta.nelem() == 0 );
      return 1;
    }
  else if( ip_beta.nelem() == 0 )
    {
      assert( ip_p.nelem() == ip_alpha.nelem() );
      return 2;
    }
  else
    {
      assert( ip_p.nelem() == ip_alpha.nelem() );
      assert( ip_p.nelem() == ip_beta.nelem() );
      return 3;
    }
}



void interp_1D (
	      Vector&     y,
	ConstVectorView   yi,
	ConstVectorView   ip_x )
{
  // Sizes
  const Index   n_in = yi.nelem();
  const Index   n_out = ip_x.nelem();

  // Asserts
  assert( y.nelem() == n_out );

  // Interpolatation weights
  ArrayOfIndex   ii(n_out); 
  Vector         w(n_out);
  get_interp_weights( ii, w, yi, ip_x );

  for( Index ix=0; ix<n_out; ix++ )
    {
      if( w[ix] == 0 )   // No interpolation, just copy data
	y[ix] = yi[ii[ix]];
      else
	y[ix] = (1-w[ix])*yi[ii[ix]] + w[ix]*yi[ii[ix]+1];
    }
}	     



void interp_abs2ppath (
	      Matrix&     abs_out,
	const Tensor3&    abs_in,
	ConstVectorView   ip_p,
	ConstVectorView   ip_alpha,
	ConstVectorView   ip_beta )
{
  // Get dimension and check consistency of index vectors
  Index   dim = get_position_dim( ip_p, ip_alpha, ip_beta );

  // Check that the input absorption tensor match the found dimension
  assert_maxdim_of_tensor( abs_in, dim+1 );

  // Check that return matrix has the correct size
  const Index   n  = ip_p.nelem();
  const Index   nv = abs_in.ncols();
  assert_size( abs_out, n, nv );

  // Get interpolation weights for pressure/altitude dimension
  ArrayOfIndex   ii1(n); 
  Vector         w1(n);
  get_interp_weights( ii1, w1, abs_in(0,Range(joker),0), ip_p );

  // Interpolate
  Index iv;
  Numeric wv;
  if( dim == 1 )
    {
      for( Index ix=0; ix<n; ix++ )
	{
	  wv = w1[ix];
	  if( wv == 0 )   // No interpolation, just copy data
	    {
	      for( iv=0; iv<nv; iv++ )
		abs_out(ix,iv) = abs_in(0,ii1[ix],iv);
	    }
	  else
	    {
	      for( iv=0; iv<nv; iv++ )
	        abs_out(ix,iv) = (1-wv) * abs_in(0,ii1[ix],iv) +
	                             wv * abs_in(0,ii1[ix]+1,iv);
	    }
	}
    }
  else
    throw runtime_error("Absorption interpolation only implemented for 1D." );
}



void empty_ppath( 
	      Ppath&      ppath,
	const Index&      dim )
{
  ppath.dim        = dim;
  ppath.np         = 0;
  ppath.i_start    = 0;
  ppath.i_stop     = 0;
  ppath.background = "Cosmic background radiation";
  ppath.ground     = 0;
  ppath.i_ground   = 0;
  ppath.pos.resize(0,dim);
  ppath.ip_pos.resize(0,dim);
  ppath.z.resize(0);
  ppath.l_step.resize(0);
  ppath.los.resize(0,0);
  ppath.tan_pos.resize(0);
}



void ppath_1d_geom (
              Index&      nz,
	      Vector&     z,
	      Vector&     ip_z,
	      Vector&     alpha,
	      Vector&     l_step,
	      Vector&     psi,
 	      Index&      ground,
              Vector&     tan_pos,
	ConstVectorView   z_field,
	const Numeric&    r_geoid,
	const Numeric&    z_ground,
	const Numeric&    ppath_lmax,
	const Numeric&    r_s,
	const Numeric&    psi_s )
{
  // Asserts
  assert( z_ground >= z_field[0] );
  assert( z_ground < last2(z_field) );
  assert( z_field.nelem() > 1 );
  assert( ppath_lmax > 0 );
  assert( r_s >= r_geoid + z_ground );
  assert( psi_s >= 0 );
  assert( psi_s <= 180 );

  // Guess a value for index return arguments.
  nz      = 0;
  ground  = 0;  

  // Get highest altitude and length of z_field
  const Index     n_zabs = z_field.nelem();
  const Numeric   z_max  = z_field[n_zabs-1];

  // Determine the lowest point of the PPATH, z1, and the zenith angle 
  // at this point, psi1. The tangent altitude, ground flag and some other 
  // variables are set on the same time.
  // The latitude distance from the sensor and z1 is denoted as alpha0.
  // The case with downward observation from inside the atmosphere needs 
  // special treatment, and it is handled by the flag do_down.
  //
  Numeric z1, psi1, alpha0;
  Index   do_down = 0;
  //
  if( psi_s <= 90 )   // Upward observation (no tangent point)
    {
      z1     = r_s - r_geoid;
      psi1   = psi_s;
      alpha0 = 0;
      tan_pos.resize(0);
    }
  else              // Downward observation (limb sounding)
    {
      tan_pos.resize(2);
      tan_pos[0] = r_s * sin(DEG2RAD*psi_s) - r_geoid; 
      tan_pos[1] = psi_s - 90.0;
      if( tan_pos[0] >= z_ground )   // No intersection with the ground
	{
	  z1     = tan_pos[0];
	  psi1   = -90;
          alpha0 = psi_s - 90.0;
	}
      else                       // Intersection with ground
	{
	  ground = 1;
	  z1     = z_ground;
	  psi1   = -RAD2DEG*asin( (r_geoid+tan_pos[0]) / (r_geoid+z_ground) );
          alpha0 = psi_s - psi1 - 180.0;
	}
      if( r_s < r_geoid + z_max )
	do_down = 1;
    }

  // The return vectors are set to be empty if z1 >= z_max
  if( z1 >= z_max )
    {
      z.resize(0);
      ip_z.resize(0);
      l_step.resize(0);
      alpha.resize(0);
    }

  // Do anything only if z1 is inside the atmosphere
  else
    {
      // Create vectors for special points.
      // Special points are z1 and r_s (if do_down).
      // This vector shall start with z1 and end with a dummy value > z_max.
      const Index    n_special = 1 + do_down;  // The dummy value is ignored
            Vector   z_special(n_special+1); 
      //
      z_special[0] = z1;
      if( do_down )
	z_special[1] = r_s - r_geoid;
      z_special[do_down+1] = z_max*2;

      // Determine index of first z_field above z1, i_above.
      Index   i_above;
      for( i_above=0; z_field[i_above]<=z1; i_above++ ) {}

      // Create a vector containing z_special and z_field levels above z1 (zs).
      // The altitudes shall be sorted and there shall be no duplicates.
      // It is assumed that the first shall be taken from z_special.
      Index    n_zs = n_zabs - i_above + n_special; 
      Vector   zs(n_zs);
      Index    i1, i2=1;
      //
      zs[0] = z_special[0];
      n_zs  = 1;            // n_zs counts now the number of values moved to zs
      for( i1=i_above; i1<n_zabs; i1++  )
	{
	  // Check if values from z_special shall be copied
	  for( ; z_special[i2] <= z_field[i1]; i2++ )
	    {
	      if( zs[n_zs] != z_special[i2] )
		{
		  zs[n_zs] = z_special[i2];
		  n_zs++;
		}
	    }
	  // Copy next z_field
	  if( zs[n_zs] != z_field[i1] )
	    {
	      zs[n_zs] = z_field[i1];
	      n_zs++;
	    }
	}

      // Calculate the length along the PPATH from z1 (ls) and the number of
      // PPATH steps needed to reach next altitude in zs (ns).
      // The geometric calculations can go wrong if float is used, and to
      // avoid this, the variables a-e are hard-coded to be double.
      Vector         ls(n_zs);
      ArrayOfIndex   ns(n_zs-1);               
      double         a, b, c, d, e;
      Index          n_sum = 0;       // n_sum is sum(ns)
      //
      // Handle first point seperately
      zs[0] = zs[0]; 
      ls[0] = 0; 
      //
      // Loop zs
      a  = cos( DEG2RAD * psi1 );
      b  = sin( DEG2RAD * psi1 );
      c  = ( r_geoid + z1 ) * b;
      c *= c;
      d  = ( r_geoid + z1 ) * a;
      for( i1=1; i1<n_zs; i1++ )
	{
	  e        = r_geoid + zs[i1];
	  ls[i1]   = sqrt( e*e - c ) - d;
          ns[i1-1] = Index( ceil( (ls[i1]-ls[i1-1])/ppath_lmax ) );
          n_sum   += ns[i1-1];
	}
      n_sum += 1;   // To account for the point at the last z_field level

      // The length of the z and alpha vectors is n_sum
      nz = n_sum;
      
      // Create the return vectors.
      z.resize(nz);
      ip_z.resize(nz);
      l_step.resize(nz-1);
      alpha.resize(nz);
      psi.resize(nz);
      //
      Numeric   dl, l;              // Different lengts along the PPATH
      Numeric   zv;                 // Next altitude
      c = r_geoid + z1;  
      d = c * c;
      n_sum = 0;                 // Re-use this variable
      //
      for( i1=0; i1<n_zs-1; i1++ )
	{
          dl = (ls[i1+1] - ls[i1]) / ns[i1];
	  for( i2=0; i2<ns[i1]; i2++ )
	    {
	      l = ls[i1] + i2*dl;
	      if( i2 == 0 )
		zv = zs[i1];
	      else
		zv = sqrt( d + l*l + 2 * c * l * a ) - r_geoid;
	      z[n_sum+i2]      = zv;
	      l_step[n_sum+i2] = dl;
	      ip_z[n_sum+i2]   = i_above - 1 + ( zv - z_field[i_above-1]) /
                                       (z_field[i_above] - z_field[i_above-1]);
	      alpha[n_sum+i2]  = alpha0 + RAD2DEG*asin( l*b / (r_geoid+zv) );
	      psi[n_sum+i2]    = psi_s - alpha[n_sum+i2];
	    }

          // Increase n_sum with the points done
          n_sum += ns[i1];

          // Check if i_above shall be increased
	  if( zs[i1] >= z_field[i_above] )
	    i_above++;
	}
      
      // Put in uppermost z_field level that is not covered above
      z[n_sum]     = z_field[n_zabs-1];
      ip_z[n_sum]  = n_zabs-1;
      alpha[n_sum] = alpha0 + RAD2DEG*asin( ls[n_zs-1]*b / (r_geoid+z[n_sum]));
      psi[n_sum]   = psi_s - alpha[n_sum];
    }
}



void ppath_1d (
	      Ppath&          ppath,
	ConstVectorView       p_grid,
	ConstVectorView       z_field,
	const Numeric&        r_geoid,
	const Numeric&        z_ground,
        const Index&          refr_on,
        const Index&          blackbody_ground,
	const Index&          cloudbox_on,
	const ArrayOfIndex&   cloudbox_limits,
	const Numeric&        ppath_lmax,
	const Numeric&        r_s,
        const Numeric&        psi_s )
{
  // Asserts (remaining variables are asserted in the sub-functions).
  assert( p_grid.nelem() == z_field.nelem() );
  assert( is_bool( refr_on ) );
  assert( is_bool( blackbody_ground ) );
  assert( is_bool( cloudbox_on ) );
  if( cloudbox_on )
    {
      assert( cloudbox_limits.nelem() == 1 );
      assert( cloudbox_limits[0] >= 1 );
      assert( cloudbox_limits[0] < p_grid.nelem() );
    }

  // Check that the sensor not is below the ground.
  if( r_s <= r_geoid + z_ground )
    {
      ostringstream os;
      os << "The sensor position cannot be below, or equal to, the ground "
         << "altitude.\nThe ground altitude is here " << z_ground/1e3 
         << " km and the sensor altitude is " << (r_s-r_geoid)/1e3 << " km."; 
      throw runtime_error( os.str() );
    }

  // Start by setting PPATH to be empty 
  empty_ppath( ppath, 1 );

  // For simplicity, a full calculation of the path is always done. For cases
  // with intersection with a blackbody ground or cloud box the path is
  // truncated later.

  // Do stuff that differ between with and without refraction.
  // Stuff going to ppath.pos/ip_pos/los is stored in temporary vectors.
  Vector ip_p, alpha, psi;
  if( refr_on )
    throw runtime_error(
   	         "1D PPATH calculations with refraction not yet implemented" );
  else
    ppath_1d_geom ( ppath.np, ppath.z, ip_p, alpha, ppath.l_step, psi, 
    		    ppath.ground, ppath.tan_pos, z_field, r_geoid,
    			      z_ground, ppath_lmax, r_s, fabs(psi_s) );

  // Set i_start and i_stop assuming that blackbody_ground and cloudbox are 0.
  if( fabs(psi_s) <= 90 )
    ppath.i_stop = 0;
  else
    ppath.i_stop = ppath.np - 1;
  ppath.i_start = ppath.np - 1;

  // Downward observation from inside the atmosphere is a special case.
  // We must here set i_stop to the index corresponding to the sensor.
  if( fabs(psi_s)>90 && r_s<r_geoid+last2(z_field) )
    for( ppath.i_stop=0; ppath.z[ppath.i_stop]!=(r_s-r_geoid); 
                                                           ppath.i_stop++ ) {}

  // Intersection with cloud box? 
  if( cloudbox_on )
    {
      // Start point inside cloud box?
      if( ip_p[ppath.i_stop]<cloudbox_limits[0] || 
           ( ip_p[ppath.i_stop]==cloudbox_limits[0] && fabs(psi_s)>90) )
	{
    	  truncate_vector( ppath.z, ppath.i_stop, ppath.i_stop );
    	  truncate_vector( ip_p, ppath.i_stop, ppath.i_stop );
    	  truncate_vector( alpha, ppath.i_stop, ppath.i_stop );
    	  truncate_vector( psi, ppath.i_stop, ppath.i_stop );
    	  truncate_vector( ppath.l_step, 0 );
	  ppath.background = "Inside cloud box";
	  ppath.i_start    = 0;
	  ppath.i_stop     = 0;
	  ppath.ground     = 0;
    	  ppath.np         = 1;
	}

      // Intersection with cloud box?
      else if( ip_p[0]<cloudbox_limits[0] )
	{
          // Find index of first point above cloud box limit
          Index i1, np;
          for( i1=0; ip_p[i1]<cloudbox_limits[0]; i1++ ) {} 

	  ppath.background = "Surface of cloud box";
    	  np               = ppath.np;
	  ppath.np         = np - i1;
	  ppath.i_stop     = ppath.np - 1;
	  ppath.i_start    = 0;
	  ppath.ground     = 0;
    	  truncate_vector( ppath.z, i1, np-1 );
    	  truncate_vector( ip_p, i1, np-1 );
    	  truncate_vector( alpha, i1, np-1 );
    	  truncate_vector( psi, i1, np-1 );
    	  truncate_vector( ppath.l_step, i1, np-2 );
	}
    }

  // Ignore part of the PPATH before ground reflection if blackbody ground.
  // The start index is then always 0. As the PPATH has then the ground as
  // background, there is no ground intersection along the PPATH (ground=0).
  // If i_stop deviates from np-1, the vectors shall be truncated 
  // (corresponds to downward observation from within the atmosphere and 
  // blackbody ground).
  // Note the "else if" here with respect to cloud box.
  else if( ppath.ground && blackbody_ground )
    {
      ppath.background = "Blackbody ground";
      ppath.i_start    = 0;
      ppath.ground     = 0;

      // Truncate vectors
      if( ppath.i_stop < ppath.np-1 )
    	{
    	  const Index    n = ppath.i_stop+1;
    	  ppath.np = n;
    	  truncate_vector( ppath.z, n );
    	  truncate_vector( ip_p, n );
    	  truncate_vector( alpha, n );
    	  truncate_vector( psi, n );
    	  truncate_vector( ppath.l_step, n-1 );
    	}
    }


  // Convert altitudes to pressures
  Vector p(ppath.np);
  interp_1D ( p, p_grid, ip_p ); 

  // Create ppath.pos, ppath.ip_pos and ppath.los
  ppath.pos.resize( ppath.np, 2 );
  ppath.pos(Range(joker),0) = p;
  ppath.pos(Range(joker),1) = alpha;
  ppath.ip_pos.resize( ppath.np, 1 );
  ppath.pos(Range(joker),0) = ip_p; 
  ppath.los.resize( ppath.np, 1 );
  ppath.los(Range(joker),0) = psi; 
}



void ppathCalc(
	      Ppath&          ppath,
	const Index&          dim,
	ConstVectorView       p_grid,
	ConstVectorView       alpha_grid,
	ConstVectorView       beta_grid,
	const Tensor3&        z_field,
 	ConstMatrixView       r_geoid,
 	ConstMatrixView       z_ground,
	const Tensor3&        e_ground,
        const Index&          refr_on,
	const Index&          cloudbox_on,
	const ArrayOfIndex&   cloudbox_limits,
	const Numeric&        ppath_lmax,
	ConstVectorView       sensor_pos,
	ConstVectorView       sensor_los )
{
  // Check input:
  // That atm. grids are strictly increasing is checked inside chk_atm_grids.
  // Consistency between the grids and sensor position checked in sub-funs.
  chk_if_in_range( "dim", dim, 1, 3 );
  chk_atm_grids( dim, p_grid, alpha_grid, beta_grid );
  chk_atm_field( "z_field", z_field, dim, p_grid, alpha_grid, beta_grid );
  for( Index page=0; page<z_field.npages(); page++ )
    {
      for( Index row=0; row<z_field.nrows(); row++ )
	{
	  // I could not get the compliler to accept a solution without dummy!!
	  Vector dummy(z_field.ncols());
	  dummy = z_field(page,row,Range(joker));
	  ostringstream os;
	  os << "z_field for latitude nr " << row << " and longitude nr " 
             << page;
	  chk_if_increasing( os.str(), dummy ); 
	}
    }
  chk_atm_surface( "r_geoid", r_geoid, dim, alpha_grid, beta_grid );
  chk_atm_surface( "z_ground", z_ground, dim, alpha_grid, beta_grid );
  for( Index row=0; row<z_ground.nrows(); row++ )
    {
      for( Index col=0; col<z_ground.ncols(); col++ )
	{
	  if( z_ground(row,col) < z_field(row,col,0) )
	    {
	      ostringstream os;
	      os << "The ground altitude (*z_ground*) cannot be below the "
		 << "lowest altitude in *z_field*.";
		if( dim > 1 )
		  os << "\nThis was found to be the case for:\nlatitude " 
		     << alpha_grid[col];
		if( dim > 2 )
		  os << "\nlongitude " << beta_grid[row];
	      throw runtime_error( os.str() );
	    }
	}
    }
  chk_atm_surface( "e_ground (for one frequency)", 
           e_ground(Range(joker),Range(joker),0), dim, alpha_grid, beta_grid );
  chk_if_bool(  "refr_on", refr_on );
  chk_if_bool(  "cloudbox_on", cloudbox_on );
  if( cloudbox_on )
    {
      if( cloudbox_limits.nelem() != dim*2-1 )
	{
	  ostringstream os;
	  os << "The array *cloudbox_limits* has incorrect length.\n"
	     << "For dim = " << dim << "the length shall be " << dim*2-1
	     << "but it is " << cloudbox_limits.nelem() << ".";
	  throw runtime_error( os.str() );
	}
      if( cloudbox_limits[0]<1 || cloudbox_limits[0]>=p_grid.nelem() )
	{
	  ostringstream os;
	  os << "Incorrect value for cloud box pressure limit found.\n"
	     << "With present length of *p_grid*, OK values are 1 - " 
	     << p_grid.nelem() << ",\nbut the limit is set to " 
	     << cloudbox_limits[0] << ".";
	  throw runtime_error( os.str() );
	}
      if( dim >= 2 )
	{
	  if( cloudbox_limits[2]<=cloudbox_limits[1] || cloudbox_limits[1]<1 ||
                                cloudbox_limits[2]>=alpha_grid.nelem()-1 )
	    {
	      ostringstream os;
	      os << "Incorrect value(s) for cloud box latitude limit(s) found."
		 << "\nValues are either out of range or upper limit is not "
		 << "greater than lower limit.\nWith present length of "
                 << "*alpha_grid*, OK values are 1 - " << alpha_grid.nelem()-2
                 << ".\nThe latitude index limits are set to " 
		 << cloudbox_limits[1] << " - " << cloudbox_limits[2] << ".";
	      throw runtime_error( os.str() );
	    }
	}
      if( dim >= 3 )
	{
	  if( cloudbox_limits[4]<=cloudbox_limits[3] || cloudbox_limits[3]<1 ||
                                cloudbox_limits[4]>=beta_grid.nelem()-1 )
	    {
	      ostringstream os;
	      os << "Incorrect value(s) for cloud box longitude limit(s) found"
		 << ".\nValues are either out of range or upper limit is not "
		 << "greater than lower limit.\nWith present length of "
                 << "*beta_grid*, OK values are 1 - " << beta_grid.nelem()-2
                 << ".\nThe longitude limits are set to " 
		 << cloudbox_limits[3] << " - " << cloudbox_limits[4] << ".";
	      throw runtime_error( os.str() );
	    }
	}
    }
  chk_if_over_0( "ppath_lmax", ppath_lmax );
  chk_vector_length( "sensor_pos", sensor_pos, dim );
  if( dim<3 )
    chk_vector_length( "sensor_los", sensor_los, 1 );
  else
    {
      chk_vector_length( "sensor_los", sensor_los, 2 );
      chk_if_in_range( "sensor_los azimuth angle", sensor_los[1], -180, 180 );
    }
  chk_if_in_range( "sensor_los zenith angle", sensor_los[0], -180, 180 );


  // Is the ground is perfect blackbody everywhere?
  Index blackbody_ground = 1;
  for( Index i=0; blackbody_ground && i<e_ground.npages(); i++ )
    {
      if( min( e_ground(i,Range(joker),Range(joker)) ) < 0.99999 )
	blackbody_ground = 0;
    }


  // The symmetry for 1D limb sounding makes it not suitable to use the same
  // algorithm for all dimensions, and 1D is treated seperately.
  if( dim == 1 )
    ppath_1d( ppath, p_grid, z_field(0,0,Range(joker)), r_geoid(0,0),
                 z_ground(0,0), refr_on, blackbody_ground, cloudbox_on, 
                   cloudbox_limits, ppath_lmax, sensor_pos[0], sensor_los[0] );
  else
    throw runtime_error(
                          "2D and 3D PPATH calculations not yet implemented" );
}



void test_new_los()
{
  Ppath ppath;
 
  // Dimension
  Index dim = 1;

  // Grids
  const Index n_p = 11;
  Vector p_grid( 0, n_p, 1 );
  Vector alpha_grid(0), beta_grid(0);

  // z_field
  Tensor3 z_field(1,1,n_p);
  for( Index i=0; i<n_p; i++ )
    z_field(0,0,i) = i*1e3;

  // Geoid and ground
  Tensor3 e_ground(1,1,5);
  Matrix r_geoid(1,1), z_ground(1,1);
  r_geoid(0,0)  = EARTH_RADIUS;
  z_ground(0,0) = 0;
  e_ground      = 0.8;

  // Booleans
  Index refr_on = 0;

  // Cloud box
  Index cloudbox_on;
  ArrayOfIndex   cloudbox_limits;
  cloudbox_on = 0;
  cloudbox_limits.resize(1);
  cloudbox_limits[0] = 4;

  // Path step length
  Numeric ppath_lmax;
  ppath_lmax = 50e3;

  // Sensor  
  Vector sensor_pos(1), sensor_los(1);
  sensor_pos[0] = EARTH_RADIUS + 20e3;
  sensor_los[0] = 98;

  ppathCalc( ppath, dim, p_grid, alpha_grid, beta_grid, z_field, r_geoid, 
 	         z_ground, e_ground, refr_on, cloudbox_on, cloudbox_limits,
                                          ppath_lmax, sensor_pos, sensor_los );

  cout << "z = \n" << ppath.z << "\n";
  //cout << "p = \n" << ppath.p << "\n";
  //cout << "ip_p = \n" << ppath.ip_p << "\n";
  cout << "alpha = \n" << ppath.pos(Range(joker),1) << "\n";
  cout << "psi = \n" << ppath.los(Range(joker),0) << "\n";
  cout << "l_step = \n" << ppath.l_step << "\n";
  cout << "nz      = " << ppath.np << "\n";
  cout << "i_start = " << ppath.i_start << "\n";
  cout << "i_stop  = " << ppath.i_stop << "\n";
  cout << "bground = " << ppath.background << "\n";
  cout << "ground  = " << ppath.ground << "\n";
  if( ppath.ground )
    cout << "i_ground= " << ppath.i_ground << "\n";
  if( ppath.tan_pos.nelem() > 0 )
  {
    cout << "z_tan   = " << ppath.tan_pos[0]/1e3 << " km\n";
    cout << "a_tan   = " << ppath.tan_pos[1] << " degs\n";
  }  

  /*
  Tensor3 abs_in(1,z_field.nelem(),2,0.0);
  for ( Index j=0; j<abs_in.nrows(); ++j )
    for ( Index k=0; k<abs_in.ncols()-1; ++k )
      abs_in(0,j,k) = j;
  //  cout << "abs_in = " << abs_in << "\n";

  Matrix abs_out(ppath.p.nelem(),2);
  
  interp_abs2ppath(abs_out,abs_in,ppath.ip_p,ppath.ip_alpha,ppath.ip_beta);
  //  cout << "abs_out = \n" << abs_out << "\n";  
  */
}
