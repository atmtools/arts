/**
   \file   continua.h

   Contains declarations of continuum functions.

   \author Stefan Buehler
   \date   2001-01-17
*/

#ifndef contiua_h
#define contiua_h

const Numeric VMRCalcLimit = 1.000e-25;


void xsec_continuum_tag( Matrix&                    xsec,
			 const string&              name,
			 const Vector&              parameters,
			 const Vector&  	    f_mono,
			 const Vector&  	    p_abs,
			 const Vector&  	    t_abs,
			 const Vector&  	    n2_abs,
			 const Vector&  	    h2o_abs,
			 const Vector&              vmr );


void check_continuum_model(const string& name);


//////////////////////////////////////////////////////////////////////////// 
// water vapor absorption continua / models
//////////////////////////////////////////////////////////////////////////// 

void Rosenkranz_H2O_self_continuum( Matrix&           xsec,
				    Numeric	      C,
				    Numeric	      x,
				    const Vector&     f_mono,
				    const Vector&     p_abs,
				    const Vector&     t_abs,
				    const Vector&     vmr	 );

void Rosenkranz_H2O_foreign_continuum( Matrix&           xsec,
				       Numeric	         C,
				       Numeric	         x,
				       const Vector&     f_mono,
				       const Vector&     p_abs,
				       const Vector&     t_abs,
				       const Vector&     vmr	 );

void MPM93_H2O_continuum( Matrix&           xsec,
			  const Vector&     f_mono,
			  const Vector&     p_abs,
			  const Vector&     t_abs,
			  const Vector&     vmr	 );

//////////////////////////////////////////////////////////////////////////// 
// water droplet and ice particle absorption
//////////////////////////////////////////////////////////////////////////// 

void MPM93WaterDropletAbs( Matrix&           xsec,
			   const Vector&   f_mono,  // frequency vector
			   const Vector&    p_abs,  // pressure vector
			   const Vector&    t_abs,  // temperature vector
			   const Vector&      vmr); // suspended water droplet density vector

void MPM93IceCrystalAbs( Matrix&           xsec,
			 const Vector&   f_mono,    // frequency vector
			 const Vector&    p_abs,    // pressure vector
			 const Vector&    t_abs,    // temperature vector
			 const Vector&      vmr	 ); // suspended ice particle density vector, 
                                                    // valid range: 0-10.0e-3 kg/m3

//////////////////////////////////////////////////////////////////////////// 
// oxygen absorption continua / models
//////////////////////////////////////////////////////////////////////////// 

void MPM93O2AbsModel( Matrix&           xsec,
		      const Vector&     f_mono,
		      const Vector&     p_abs,
		      const Vector&     t_abs,
		      const Vector&     h2o_abs,
		      const Vector&     vmr );

void MPM93_O2_continuum( Matrix&           xsec,
			 const Vector&     f_mono,
			 const Vector&     p_abs,
			 const Vector&     t_abs,
			 const Vector&     h2o_abs,
			 const Vector&     vmr	 );

void Rosenkranz_O2_continuum( Matrix&           xsec,
			      const Vector&  	f_mono,
			      const Vector&  	p_abs,
			      const Vector&  	t_abs,
			      const Vector&     h2o_abs,
			      const Vector&     vmr	 );

//////////////////////////////////////////////////////////////////////////// 
// nitrogen absorption continua / models
//////////////////////////////////////////////////////////////////////////// 

void MPM93_N2_continuum( Matrix&           xsec,
			 const Vector&     f_mono,
			 const Vector&     p_abs,
			 const Vector&     t_abs,
			 const Vector&     h2o_abs,
			 const Vector&     vmr	 );

void Rosenkranz_N2_self_continuum( Matrix&           xsec,
				   const Vector&     f_mono,
				   const Vector&     p_abs,
				   const Vector&     t_abs,
				   const Vector&     vmr    );

void General_N2_self_continuum(    Matrix&           xsec,
                                   Numeric           C,
                                   Numeric           xf,
                                   Numeric           xt,
                                   Numeric           xp,
				   const Vector&     f_mono,
				   const Vector&     p_abs,
				   const Vector&     t_abs,
				   const Vector&     vmr    );

//////////////////////////////////////////////////////////////////////////// 
// carbon dioxide absorption continua / models
//////////////////////////////////////////////////////////////////////////// 

void Rosenkranz_CO2_self_continuum( Matrix&           xsec,
				    const Vector&     f_mono,
				    const Vector&     p_abs,
				    const Vector&     t_abs,
				    const Vector&     vmr	 );

void Rosenkranz_CO2_foreign_continuum( Matrix&           xsec,
				       const Vector&     f_mono,
				       const Vector&     p_abs,
				       const Vector&     t_abs,
				       const Vector&     n2_abs,
				       const Vector&     vmr	 );

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
