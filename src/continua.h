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


//////////////////////////////////////////////////////////////////////////// 
// water vapor absorption continua / models
//////////////////////////////////////////////////////////////////////////// 

void Rosenkranz_H2O_self_continuum( MATRIX&           xsec,
				    Numeric	      C,
				    Numeric	      x,
				    const VECTOR&     f_mono,
				    const VECTOR&     p_abs,
				    const VECTOR&     t_abs,
				    const VECTOR&     vmr	 );

void Rosenkranz_H2O_foreign_continuum( MATRIX&           xsec,
				       Numeric	         C,
				       Numeric	         x,
				       const VECTOR&     f_mono,
				       const VECTOR&     p_abs,
				       const VECTOR&     t_abs,
				       const VECTOR&     vmr	 );

void MPM93_H2O_continuum( MATRIX&           xsec,
			  const VECTOR&     f_mono,
			  const VECTOR&     p_abs,
			  const VECTOR&     t_abs,
			  const VECTOR&     vmr	 );

//////////////////////////////////////////////////////////////////////////// 
// water droplet and ice particle absorption
//////////////////////////////////////////////////////////////////////////// 

void MPM93WaterDropletAbs( MATRIX&           xsec,
			   const VECTOR&   f_mono,  // frequency vector
			   const VECTOR&    p_abs,  // pressure vector
			   const VECTOR&    t_abs,  // temperature vector
			   const VECTOR&      vmr); // suspended water droplet density vector

void MPM93IceCrystalAbs( MATRIX&           xsec,
			 const VECTOR&   f_mono,    // frequency vector
			 const VECTOR&    p_abs,    // pressure vector
			 const VECTOR&    t_abs,    // temperature vector
			 const VECTOR&      vmr	 ); // suspended ice particle density vector, 
                                                    // valid range: 0-10.0e-3 kg/m3

//////////////////////////////////////////////////////////////////////////// 
// oxygen absorption continua / models
//////////////////////////////////////////////////////////////////////////// 

void MPM93O2AbsModel( MATRIX&           xsec,
		      const VECTOR&     f_mono,
		      const VECTOR&     p_abs,
		      const VECTOR&     t_abs,
		      const VECTOR&     h2o_abs,
		      const VECTOR&     vmr );

void MPM93_O2_continuum( MATRIX&           xsec,
			 const VECTOR&     f_mono,
			 const VECTOR&     p_abs,
			 const VECTOR&     t_abs,
			 const VECTOR&     h2o_abs,
			 const VECTOR&     vmr	 );

void Rosenkranz_O2_continuum( MATRIX&           xsec,
			      const VECTOR&  	f_mono,
			      const VECTOR&  	p_abs,
			      const VECTOR&  	t_abs,
			      const VECTOR&     h2o_abs,
			      const VECTOR&     vmr	 );

//////////////////////////////////////////////////////////////////////////// 
// nitrogen absorption continua / models
//////////////////////////////////////////////////////////////////////////// 

void MPM93_N2_continuum( MATRIX&           xsec,
			 const VECTOR&     f_mono,
			 const VECTOR&     p_abs,
			 const VECTOR&     t_abs,
			 const VECTOR&     h2o_abs,
			 const VECTOR&     vmr	 );

void Rosenkranz_N2_self_continuum( MATRIX&           xsec,
				   const VECTOR&     f_mono,
				   const VECTOR&     p_abs,
				   const VECTOR&     t_abs,
				   const VECTOR&     vmr    );

//////////////////////////////////////////////////////////////////////////// 
// carbon dioxide absorption continua / models
//////////////////////////////////////////////////////////////////////////// 

void Rosenkranz_CO2_self_continuum( MATRIX&           xsec,
				    const VECTOR&     f_mono,
				    const VECTOR&     p_abs,
				    const VECTOR&     t_abs,
				    const VECTOR&     vmr	 );

void Rosenkranz_CO2_foreign_continuum( MATRIX&           xsec,
				       const VECTOR&     f_mono,
				       const VECTOR&     p_abs,
				       const VECTOR&     t_abs,
				       const VECTOR&     n2_abs,
				       const VECTOR&     vmr	 );

//////////////////////////////////////////////////////////////////////////// 
// help functions
//////////////////////////////////////////////////////////////////////////// 

Numeric MPMLineShapeFunction( Numeric gamma, 
			      Numeric fl, 
			      Numeric f);

Numeric MPMLineShapeO2Function( Numeric gamma, 
				Numeric fl, 
				Numeric f,
                                Numeric delta);

Numeric WVSatPressureLiquidWater(Numeric t);

Numeric WVSatPressureIce(Numeric t);

#endif // contiua_h
