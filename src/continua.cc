/**
   \file   continua.cc

   This contains the physical functions for the various continua. 

   \author Stefan Buehler
   \date   2001-01-17
*/

#include "vecmat.h"
#include "absorption.h"

/**
   Calculates the self broadened continuum according to the simple
   empirical function as formulated by Rosenkranz98:

   xsec += C * (300/t_abs)^x * f_mono^2 * (300/t_abs)^3 * p_abs^2 * vmr

   See the list of parameters below for the meaning of these
   variables. The equation is slightly modified from Rosenkranz's
   original, because we want to return only the cross section, not the
   absorption coefficient. (Only vmr, not vmr^2 at the end of the
   equation.) 

   This xsec is understood to parameterize the difference between
   observed absorption cross section and the one calculated from the
   pure line spectrum. The continuum is added to the previous content
   of xsec!

   \retval xsec  Absorption cross section, defined such that the
                 absorption coefficient alpha is:<br>
                 alpaha [1/m] = xsec * VMR.<br>
		 The functions adds to xsec, rather than replacing the
		 previous content. 

   \param  C       Continuum coefficient.
   \param  x       Temperature exponenet.
   \param  f_mono  Frequency grid.
   \param  p_abs   Pressure grid.
   \param  t_abs   Temperatures associated with p_abs.
   \param  vmr     Volume mixing ratio of the calculated species (H2O).

   References:
   
   P.W. Rosenkranz, `Water vapor microwave continuum absorption: A
   comparison of measurements and models', Radio Science, Vol. 33, No
   4, Pages 919-928, July-August 1998.


   \author Stefan Buehler
   \date   2001-01-17
*/
void simple_empirical_self_continuum( MATRIX&           xsec,
				      Numeric		C,
				      Numeric		x,
				      const VECTOR&  	f_mono,
				      const VECTOR&  	p_abs,
				      const VECTOR&  	t_abs,
				      const VECTOR&     vmr	 )
{
  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop pressure/temperature:
  for ( size_t i=0; i<n_p; ++i )
    {
      // Dummy scalar holds everything except the quadratic frequency dependence:
      Numeric dummy =
	C * pow( 300./t_abs[i], x+3. ) * p_abs[i] * p_abs[i] * vmr[i];

      // Loop frequency:
      for ( size_t s=0; s<n_f; ++s )
	{
	  xsec[s][i] += dummy * f_mono[s] * f_mono[s];
	  //	  cout << "xsec[" << s << "][" << i << "]: " << xsec[s][i] << "\n";
	}
    }
}


/** Calculates continuum absorption for one continuum tag. Note, that
    only one tag can be taken at a time. That means for water vapor
    you will have to call this function two times, once with the
    self-continuum tag and once with the foreign continuum tag.

    Calculated is the absorption cross section, that means you have to
    multiply this with the VMR in order to get the absorption
    coefficient:

    alpaha [1/m] = xsec * VMR

    \retval xsec    Cross section of one continuum tag.
    \param  tag     The tag for which to do the calculation.
    \param  f_mono  Frequency grid.
    \param  p_abs   Pressure grid.
    \param  t_abs   Temperatures associated with p_abs.
    \param  h2o_abs Total volume mixing ratio of water vapor. This
                    will be needed only for the oxygen continuum
    \param  vmr     Volume mixing ratio of the calculated species.

    \date   2001-01-16
    \author Stefan Buehler */
void xsec_continuum_tag( MATRIX&                    xsec,
			 const OneTag&              tag,
			 const VECTOR&  	    f_mono,
			 const VECTOR&  	    p_abs,
			 const VECTOR&  	    t_abs,
			 const VECTOR&  	    h2o_abs,
			 const VECTOR&              vmr )
{
  extern const ARRAY<SpeciesRecord> species_data;

  // Get only the continuum name. The full tag name is something like:
  // H2O-HITRAN96Self-*-*. We want only the `H2O-HITRAN96Self' part:
  const string name =
    species_data[tag.Species()].Name() + "-"
    + species_data[tag.Species()].Isotope()[tag.Isotope()].Name();
  
  // A simple switch statement does not work here, because the
  // switching condition is not a simple value. So we need to use a
  // chain of if-else statements.
  if ( "H2O-HITRAN96Self"==name )
    {
      // FIXME: Thomas, the values for C and x are just my guess from the
      // Rosenkranz paper. You have to put in the right ones. Note
      // that xsec is in [1/m], f_mono is in [Hz], p_abs is in [Pa], and t_abs is in
      // [K]. Therefor, the exponent in C is very different from
      // Rosenkranz's paper. I probably have not got it right, I did
      // this very quickly. You have to find out the right exponent yourself.
      simple_empirical_self_continuum( xsec,
				       7.8e-27,
				       4.5,
				       f_mono,
				       p_abs,
				       t_abs,
				       vmr );
    }
  else if ( "H2O-HITRAN96Foreign"==name )
    {
      // FIXME: Thomas, you have to fill in this one! You should implement a
      // function simple_empirical_foreign_continuum, similar to my
      // function for the self continuum.
    }
  else
    {
      ostringstream os;
      os << "Continuum tag `" << name << "' not implemented.";
      throw runtime_error(os.str());
    }
}

