/**
   \file   continua.h

   Contains declarations of continuum functions.

   \author Stefan Buehler
   \date   2001-01-17
*/

#ifndef contiua_h
#define contiua_h


void xsec_continuum_tag( MATRIX&                    xsec,
			 const OneTag&              tag,
			 const VECTOR&  	    f_mono,
			 const VECTOR&  	    p_abs,
			 const VECTOR&  	    t_abs,
			 const VECTOR&  	    h2o_abs,
			 const VECTOR&              vmr );


#endif // contiua_h
