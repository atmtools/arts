/* Copyright (C) 2001 Thomas Kuhn    <tkuhn@uni-bremen.de>
                      Stefan Buehler <sbuehler@uni-bremen.de>

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

/**
   \file   continua.h

   This header file contains all the declarations of the implemented 
   continua and full absorption (lines+continuum) models.

   \author Thomas Kuhn
   \date   2001-11-05
*/

#ifndef continua_h
#define continua_h

#include "matpackI.h"


//////////////////////////////////////////////////////////////////////////// 
// entry function to all continua and full model functions
//////////////////////////////////////////////////////////////////////////// 

void xsec_continuum_tag( MatrixView         xsec,       // calculated x-section
			 const String&      name,       // model name
			 ConstVectorView    parameters, // model 
			 const String&      model,      // model option
			 ConstVectorView    f_mono,     // frequency vector
			 ConstVectorView    p_abs,      // pressure vector
			 ConstVectorView    t_abs,      // temperature vector 
			 ConstVectorView    n2_abs,     // N2 vmr profile
			 ConstVectorView    h2o_abs,    // H2O vmr profile
			 ConstVectorView    vmr );      // species vmr profile

//////////////////////////////////////////////////////////////////////////// 
// check of consistency of all full and continua absorption models
//////////////////////////////////////////////////////////////////////////// 

void check_continuum_model(const String& name);


//////////////////////////////////////////////////////////////////////////// 
// water vapor line+continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void MPM87H2OAbsModel( MatrixView        xsec,       // calculated x-section
                       const Numeric	 CC,         // continuum scale factor 
		       const Numeric	 CL,         // line strength scale factor
		       const Numeric	 CW,         // line broadening scale factor
                       const String&     model,      // model option
		       ConstVectorView   f_mono,     // frequency vector
		       ConstVectorView   p_abs,      // pressure vector
		       ConstVectorView   t_abs,      // temperature vector
		       ConstVectorView   vmr );      // H2O vmr profile

void MPM89H2OAbsModel( MatrixView        xsec,       // calculated x-section
                       const Numeric	 CCin,       // continuum scale factor 
		       const Numeric	 CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
		       const String&     model,      // model option
		       ConstVectorView   f_mono,     // frequency vector
		       ConstVectorView   p_abs,      // pressure vector
		       ConstVectorView   t_abs,      // temperature vector
		       ConstVectorView   vmr );      // H2O vmr profile

void MPM93H2OAbsModel( MatrixView        xsec,
                       const Numeric     CCin,       // continuum scale factor 
		       const Numeric	 CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
		       const String&     model,      // model option
		       ConstVectorView   f_mono,     // frequency vector
		       ConstVectorView   p_abs,      // pressure vector
		       ConstVectorView   t_abs,      // temperature vector
		       ConstVectorView   vmr );      // H2O vmr profile

void PWR98H2OAbsModel( MatrixView        xsec,       // calculated x-section
		       const Numeric	 CCin,       // continuum scale factor 
		       const Numeric     CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
		       const String&     model,      // model option
		       ConstVectorView   f_mono,     // frequency vector
		       ConstVectorView   p_abs,      // pressure vector
		       ConstVectorView   t_abs,      // temperature vector
		       ConstVectorView   vmr );      // H2O vmr profile

void CP98H2OAbsModel( MatrixView        xsec,        // calculated x-section
                      const Numeric     CCin,        // continuum scale factor 
		      const Numeric     CLin,        // line strength scale factor
		      const Numeric     CWin,        // line broadening scale factor
		      const String&     model,       // model option
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   vmr );       // H2O vmr profile

//////////////////////////////////////////////////////////////////////////// 
// water vapor continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void Pardo_ATM_H2O_ForeignContinuum( MatrixView          xsec,   // calculated x-section
				     const Numeric       Cin,    // model parameter
				     const String&       model,  // model option
				     ConstVectorView     f_mono, // frequency vector
				     ConstVectorView     p_abs,  // pressure vector
				     ConstVectorView     t_abs,  // temperature vector 
				     ConstVectorView     vmr);   // H2O vmr profile

void Standard_H2O_self_continuum( MatrixView        xsec,        // calculated x-section
				  const Numeric     C,           // model parameter
				  const Numeric     x,           // model parameter
				  const String&     model,       // model option
				  ConstVectorView   f_mono,      // frequency vector
				  ConstVectorView   p_abs,       // pressure vector
				  ConstVectorView   t_abs,       // temperature vector 
				  ConstVectorView   vmr);        // H2O vmr profile

void Standard_H2O_foreign_continuum( MatrixView        xsec,     // calculated x-section
				     const Numeric	 C,      // model parameter
				     const Numeric	 x,      // model parameter
				     const String&     model,    // model option
				     ConstVectorView   f_mono,   // frequency vector
				     ConstVectorView   p_abs,    // pressure vector
				     ConstVectorView   t_abs,    // temperature vector 
				     ConstVectorView   vmr);     // H2O vmr profile

void MaTipping_H2O_foreign_continuum( MatrixView        xsec,     // calculated x-section
				     const Numeric	 C,      // model parameter
				     const Numeric	 x,      // model parameter
				     const String&     model,    // model option
				     ConstVectorView   f_mono,   // frequency vector
				     ConstVectorView   p_abs,    // pressure vector
				     ConstVectorView   t_abs,    // temperature vector 
				     ConstVectorView   vmr);     // H2O vmr profile

void MPM93_H2O_continuum( MatrixView        xsec,                // calculated x-section
			  const Numeric	    fcenter,             // model parameter
			  const Numeric	    b1,                  // model parameter
			  const Numeric	    b2,                  // model parameter
			  const Numeric	    b3,                  // model parameter
			  const Numeric	    b4,                  // model parameter
			  const Numeric	    b5,                  // model parameter
			  const Numeric	    b6,                  // model parameter
			  const String&     model,               // model option
			  ConstVectorView   f_mono,              // frequency vector
			  ConstVectorView   p_abs,               // pressure vector
			  ConstVectorView   t_abs,               // temperature vector
			  ConstVectorView   vmr	 );              // H2O vmr profile

void CKD24_H20( MatrixView          xsec,      // calculated x-section
		int                 isf,       // flag if self or foreign cont.
		const Numeric       Cin,       // model scaling factor
		const String&       model,     // model option
		ConstVectorView     f_mono,    // frequency vector
		ConstVectorView     p_abs,     // pressure vector
		ConstVectorView     t_abs,     // temperature vector
		ConstVectorView     vmr,       // H2O vmr profile
                ConstVectorView     n2_abs );  // N2 vmr profile

//////////////////////////////////////////////////////////////////////////// 
// oxygen line+continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void MPM85O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   h2o_abs,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile


void MPM87O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   h2o_abs,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile


void MPM89O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   h2o_abs,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile


void MPM92O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   h2o_abs,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile


void MPM93O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   h2o_abs,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile

void PWR93O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      const String&     version,     // model version 1993 or 1988
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   h2o_abs,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile

//////////////////////////////////////////////////////////////////////////// 
// oxygen continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void MPM93_O2_continuum( MatrixView        xsec,             // calculated x-section
			 const Numeric     S0in,             // model parameter
			 const Numeric     G0in,             // model parameter
			 const Numeric     XSOin,            // model parameter
			 const Numeric     XG0in,            // model parameter
			 const String&     model,            // model option
			 ConstVectorView   f_mono,           // frequency vector
			 ConstVectorView   p_abs,            // pressure vector
			 ConstVectorView   t_abs,            // temperature vector
			 ConstVectorView   h2o_abs,          // H2O vmr profile
			 ConstVectorView   vmr	 );          // O2 vmr profile

void Rosenkranz_O2_continuum( MatrixView        xsec,        // calculated x-section
			      const Numeric     S0in,        // model parameter
			      const Numeric     G0in,        // model parameter
			      const Numeric     XSOin,       // model parameter
			      const Numeric     XG0in,       // model parameter
			      const String&     model,       // model option
			      ConstVectorView  	f_mono,      // frequency vector
			      ConstVectorView  	p_abs,       // pressure vector
			      ConstVectorView  	t_abs,       // temperature vector
			      ConstVectorView   h2o_abs,     // H2O vmr profile
			      ConstVectorView   vmr);        // O2 vmr profile

//////////////////////////////////////////////////////////////////////////// 
// nitrogen continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void BF86_CIA_N2( MatrixView          xsec,      // calculated x-section
		  const Numeric       Cin,       // model parameter
		  const String&       model,     // model option 
		  ConstVectorView     f_mono,    // frequency vector
		  ConstVectorView     p_abs,     // pressure vector
		  ConstVectorView     t_abs,     // temperature vector
		  ConstVectorView     vmr   );   // N2 vmr profile 

void MPM93_N2_continuum( MatrixView        xsec,                  // calculated x-section
			 const Numeric     Cin,                   // model parameter
			 const Numeric     Gin,                   // model parameter
			 const Numeric     xTin,                  // model parameter
			 const Numeric     xfin,                  // model parameter
			 const String&     model,                 // model option
			 ConstVectorView   f_mono,                // frequency vector
			 ConstVectorView   p_abs,                 // pressure vector
			 ConstVectorView   t_abs,                 // temperature vector
			 ConstVectorView   h2o_abs,               // H2O vmr profile
			 ConstVectorView   vmr	 );               // N2 vmr profile

void Rosenkranz_N2_self_continuum( MatrixView        xsec,        // calculated x-section
				   const Numeric     Cin,         // model parameter
				   const Numeric     xin,         // model parameter
				   const String&     model,       // model option
				   ConstVectorView   f_mono,      // frequency vector
				   ConstVectorView   p_abs,       // pressure vector
				   ConstVectorView   t_abs,       // temperature vector
				   ConstVectorView   vmr );       // N2 vmr profile

void Standard_N2_self_continuum(   MatrixView        xsec,        // calculated x-section
                                   const Numeric     Cin,         // model parameter
                                   const Numeric     xfin,        // model parameter
                                   const Numeric     xtin,        // model parameter
                                   const Numeric     xpin,        // model parameter
				   const String&     model,       // model option
				   ConstVectorView   f_mono,      // frequency vector
				   ConstVectorView   p_abs,       // pressure vector
				   ConstVectorView   t_abs,       // temperature vector
				   ConstVectorView   vmr );       // N2 vmr profile

void Pardo_ATM_N2_dry_continuum( MatrixView        xsec,          // calculated x-section
				  const Numeric     Cin,          // model parameter
				  const String&     model,        // model option
				  ConstVectorView   f_mono,       // frequency vector
				  ConstVectorView   p_abs,        // pressure vector
				  ConstVectorView   t_abs,        // temperature vector
				  ConstVectorView   vmr,	  // N2 vmr profile
				  ConstVectorView   h2o_abs);     // H2O vmr profile

//////////////////////////////////////////////////////////////////////////// 
// carbon dioxide continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void Rosenkranz_CO2_self_continuum( MatrixView        xsec,       // calculated x-section
				    const Numeric     C,          // model parameter
				    const Numeric     x,          // model parameter
				    const String&     model,      // model option
				    ConstVectorView   f_mono,     // frequency vector
				    ConstVectorView   p_abs,      // pressure vector
				    ConstVectorView   t_abs,      // temperature vector
				    ConstVectorView   vmr );      // CO2 vmr profile

void Rosenkranz_CO2_foreign_continuum( MatrixView        xsec,    // calculated x-section
				       const Numeric     C,       // model parameter
				       const Numeric     x,       // model parameter
				       const String&     model,   // model option
				       ConstVectorView   f_mono,  // frequency vector
				       ConstVectorView   p_abs,   // pressure vector
				       ConstVectorView   t_abs,   // temperature vector
				       ConstVectorView   n2_abs,  // N2 vmr profile
				       ConstVectorView   vmr );   // CO2 vmr profile

//////////////////////////////////////////////////////////////////////////// 
// water droplet and ice particle absorption (clouds)
//////////////////////////////////////////////////////////////////////////// 

void MPM93WaterDropletAbs( MatrixView        xsec,     // calculated x-section
			   const Numeric     CC,       // model parameter
			   const Numeric     CG,       // model parameter
			   const Numeric     CE,       // model parameter
			   const String&     model,    // model option
			   ConstVectorView   f_mono,   // frequency vector
			   ConstVectorView   p_abs,    // pressure vector
			   ConstVectorView   t_abs,    // temperature vector
			   ConstVectorView   vmr);     // suspended water droplet density vector

void MPM93IceCrystalAbs( MatrixView        xsec,       // calculated x-section
			 const Numeric     CC,         // model parameter
			 const Numeric     CA,         // model parameter
			 const Numeric     CB,         // model parameter
			 const String&     model,      // model option
			 ConstVectorView   f_mono,     // frequency vector
			 ConstVectorView   p_abs,      // pressure vector
			 ConstVectorView   t_abs,      // temperature vector
			 ConstVectorView   vmr	 );    // suspended ice particle density vector, 

void MPM93RainExt( MatrixView        xsec,       // calculated x-section
		   const Numeric     CE,         // model parameter
		   const Numeric     CA,         // model parameter
		   const Numeric     CB,         // model parameter
		   const String&     model,      // model option
		   ConstVectorView   f_mono,     // frequency vector
		   ConstVectorView   p_abs,      // pressure vector
		   ConstVectorView   t_abs,      // temperature vector
		   ConstVectorView   vmr   );    // rain rate vector, 

//////////////////////////////////////////////////////////////////////////// 
// help functions
//////////////////////////////////////////////////////////////////////////// 

Numeric MPMLineShapeFunction( const Numeric gamma,     // line width
			      const Numeric fl,        // line center frequency
			      const Numeric f);        // frequency

Numeric MPMLineShapeO2Function( const Numeric gamma,   // line width
				const Numeric fl,      // line center frequency
				const Numeric f,       // frequency
                                const Numeric delta);  // line coupling

Numeric WVSatPressureLiquidWater(const Numeric t);     // temperature

Numeric WVSatPressureIce(const Numeric t);             // temperature

#endif // continua_h
