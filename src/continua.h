/**
   \file   continua.h

   Contains declarations of continuum functions.

   \author Stefan Buehler
   \date   2001-01-17
*/

#ifndef contiua_h
#define contiua_h

const Numeric VMRCalcLimit = 1.000e-25;


void xsec_continuum_tag( MATRIX&                    xsec,
			 const string&              name,
			 const VECTOR&              parameters,
			 const VECTOR&  	    f_mono,
			 const VECTOR&  	    p_abs,
			 const VECTOR&  	    t_abs,
			 const VECTOR&  	    n2_abs,
			 const VECTOR&  	    h2o_abs,
			 const VECTOR&              vmr );


void check_continuum_model(const string& name);


#endif // contiua_h
