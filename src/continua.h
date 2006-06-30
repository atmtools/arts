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

   This is the file from arts-1-0, back-ported to arts-1-1.

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
			 ConstVectorView    f_grid,     // frequency vector
			 ConstVectorView    abs_p,      // pressure vector
			 ConstVectorView    abs_t,      // temperature vector 
			 ConstVectorView    abs_n2,     // N2 vmr profile
			 ConstVectorView    abs_h2o,    // H2O vmr profile
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
		       ConstVectorView   f_grid,     // frequency vector
		       ConstVectorView   abs_p,      // pressure vector
		       ConstVectorView   abs_t,      // temperature vector
		       ConstVectorView   vmr );      // H2O vmr profile

void MPM89H2OAbsModel( MatrixView        xsec,       // calculated x-section
                       const Numeric	 CCin,       // continuum scale factor 
		       const Numeric	 CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
		       const String&     model,      // model option
		       ConstVectorView   f_grid,     // frequency vector
		       ConstVectorView   abs_p,      // pressure vector
		       ConstVectorView   abs_t,      // temperature vector
		       ConstVectorView   vmr );      // H2O vmr profile

void MPM93H2OAbsModel( MatrixView        xsec,
                       const Numeric     CCin,       // continuum scale factor 
		       const Numeric	 CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
		       const String&     model,      // model option
		       ConstVectorView   f_grid,     // frequency vector
		       ConstVectorView   abs_p,      // pressure vector
		       ConstVectorView   abs_t,      // temperature vector
		       ConstVectorView   vmr );      // H2O vmr profile

void PWR98H2OAbsModel( MatrixView        xsec,       // calculated x-section
		       const Numeric	 CCin,       // continuum scale factor 
		       const Numeric     CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
		       const String&     model,      // model option
		       ConstVectorView   f_grid,     // frequency vector
		       ConstVectorView   abs_p,      // pressure vector
		       ConstVectorView   abs_t,      // temperature vector
		       ConstVectorView   vmr );      // H2O vmr profile

void CP98H2OAbsModel( MatrixView        xsec,        // calculated x-section
                      const Numeric     CCin,        // continuum scale factor 
		      const Numeric     CLin,        // line strength scale factor
		      const Numeric     CWin,        // line broadening scale factor
		      const String&     model,       // model option
		      ConstVectorView   f_grid,      // frequency vector
		      ConstVectorView   abs_p,       // pressure vector
		      ConstVectorView   abs_t,       // temperature vector
		      ConstVectorView   vmr );       // H2O vmr profile

//////////////////////////////////////////////////////////////////////////// 
// water vapor continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void Pardo_ATM_H2O_ForeignContinuum( MatrixView          xsec,   // calculated x-section
				     const Numeric       Cin,    // model parameter
				     const String&       model,  // model option
				     ConstVectorView     f_grid, // frequency vector
				     ConstVectorView     abs_p,  // pressure vector
				     ConstVectorView     abs_t,  // temperature vector 
				     ConstVectorView     vmr);   // H2O vmr profile

void Standard_H2O_self_continuum( MatrixView        xsec,        // calculated x-section
				  const Numeric     C,           // model parameter
				  const Numeric     x,           // model parameter
				  const String&     model,       // model option
				  ConstVectorView   f_grid,      // frequency vector
				  ConstVectorView   abs_p,       // pressure vector
				  ConstVectorView   abs_t,       // temperature vector 
				  ConstVectorView   vmr);        // H2O vmr profile

void Standard_H2O_foreign_continuum( MatrixView        xsec,     // calculated x-section
				     const Numeric	 C,      // model parameter
				     const Numeric	 x,      // model parameter
				     const String&     model,    // model option
				     ConstVectorView   f_grid,   // frequency vector
				     ConstVectorView   abs_p,    // pressure vector
				     ConstVectorView   abs_t,    // temperature vector 
				     ConstVectorView   vmr);     // H2O vmr profile

void MaTipping_H2O_foreign_continuum( MatrixView        xsec,     // calculated x-section
				     const Numeric	 C,      // model parameter
				     const Numeric	 x,      // model parameter
				     const String&     model,    // model option
				     ConstVectorView   f_grid,   // frequency vector
				     ConstVectorView   abs_p,    // pressure vector
				     ConstVectorView   abs_t,    // temperature vector 
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
			  ConstVectorView   f_grid,              // frequency vector
			  ConstVectorView   abs_p,               // pressure vector
			  ConstVectorView   abs_t,               // temperature vector
			  ConstVectorView   vmr	 );              // H2O vmr profile


void CKD_222_self_h2o( MatrixView          xsec,
		       const Numeric       Cin,
		       const String&       model,
		       ConstVectorView     f_grid,
		       ConstVectorView     abs_p,
		       ConstVectorView     abs_t,
		       ConstVectorView     vmr,
		       ConstVectorView     abs_n2 );

void CKD_222_foreign_h2o( MatrixView          xsec,
			  const Numeric       Cin,
			  const String&       model,
			  ConstVectorView     f_grid,
			  ConstVectorView     abs_p,
			  ConstVectorView     abs_t,
			  ConstVectorView     vmr,
			  ConstVectorView     abs_n2 );

void CKD_242_self_h2o( MatrixView          xsec,
		       const Numeric       Cin,
		       const String&       model,
		       ConstVectorView     f_grid,
		       ConstVectorView     abs_p,
		       ConstVectorView     abs_t,
		       ConstVectorView     vmr,
		       ConstVectorView     abs_n2 );

void CKD24_H20( MatrixView          xsec,      // calculated x-section
		int                 isf,       // flag if self or foreign cont.
		const Numeric       Cin,       // model scaling factor
		const String&       model,     // model option
		ConstVectorView     f_grid,    // frequency vector
		ConstVectorView     abs_p,     // pressure vector
		ConstVectorView     abs_t,     // temperature vector
		ConstVectorView     vmr,       // H2O vmr profile
                ConstVectorView     abs_n2 );  // N2 vmr profile

void CKD_242_foreign_h2o( MatrixView          xsec,
			  const Numeric       Cin,
			  const String&       model,
			  ConstVectorView     f_grid,
			  ConstVectorView     abs_p,
			  ConstVectorView     abs_t,
			  ConstVectorView     vmr,
			  ConstVectorView     abs_n2 );

void CKD_mt_100_self_h2o( MatrixView          xsec,
			  const Numeric       Cin,
			  const String&       model,
			  ConstVectorView     f_grid,
			  ConstVectorView     abs_p,
			  ConstVectorView     abs_t,
			  ConstVectorView     vmr,
			  ConstVectorView     abs_n2 );

void CKD_mt_100_foreign_h2o( MatrixView          xsec,
			     const Numeric       Cin,
			     const String&       model,
			     ConstVectorView     f_grid,
			     ConstVectorView     abs_p,
			     ConstVectorView     abs_t,
			     ConstVectorView     vmr,
			     ConstVectorView     abs_n2 );

//////////////////////////////////////////////////////////////////////////// 
// oxygen line+continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void MPM85O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_grid,      // frequency vector
		      ConstVectorView   abs_p,       // pressure vector
		      ConstVectorView   abs_t,       // temperature vector
		      ConstVectorView   abs_h2o,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile


void MPM87O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_grid,      // frequency vector
		      ConstVectorView   abs_p,       // pressure vector
		      ConstVectorView   abs_t,       // temperature vector
		      ConstVectorView   abs_h2o,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile


void MPM89O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_grid,      // frequency vector
		      ConstVectorView   abs_p,       // pressure vector
		      ConstVectorView   abs_t,       // temperature vector
		      ConstVectorView   abs_h2o,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile


void MPM92O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_grid,      // frequency vector
		      ConstVectorView   abs_p,       // pressure vector
		      ConstVectorView   abs_t,       // temperature vector
		      ConstVectorView   abs_h2o,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile


void MPM93O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_grid,      // frequency vector
		      ConstVectorView   abs_p,       // pressure vector
		      ConstVectorView   abs_t,       // temperature vector
		      ConstVectorView   abs_h2o,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile

void PWR93O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      const String&     version,     // model version 1993 or 1988
		      ConstVectorView   f_grid,      // frequency vector
		      ConstVectorView   abs_p,       // pressure vector
		      ConstVectorView   abs_t,       // temperature vector
		      ConstVectorView   abs_h2o,     // H2O vmr profile
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
			 ConstVectorView   f_grid,           // frequency vector
			 ConstVectorView   abs_p,            // pressure vector
			 ConstVectorView   abs_t,            // temperature vector
			 ConstVectorView   abs_h2o,          // H2O vmr profile
			 ConstVectorView   vmr	 );          // O2 vmr profile

void Rosenkranz_O2_continuum( MatrixView        xsec,        // calculated x-section
			      const Numeric     S0in,        // model parameter
			      const Numeric     G0in,        // model parameter
			      const Numeric     XSOin,       // model parameter
			      const Numeric     XG0in,       // model parameter
			      const String&     model,       // model option
			      ConstVectorView  	f_grid,      // frequency vector
			      ConstVectorView  	abs_p,       // pressure vector
			      ConstVectorView  	abs_t,       // temperature vector
			      ConstVectorView   abs_h2o,     // H2O vmr profile
			      ConstVectorView   vmr);        // O2 vmr profile

void CKD_mt_CIAfun_o2( MatrixView         xsec,        // calculated x-section
		      const Numeric       Cin,         // scaling factor
		      const String&       model,       // model option
		      ConstVectorView     f_grid,      // frequency vector
		      ConstVectorView     abs_p,       // pressure vector
		      ConstVectorView     abs_t,       // temperature vector
		      ConstVectorView     vmr );       // O2 vmr profile

void CKD_mt_v0v0_o2( MatrixView          xsec,        // calculated x-section
		     const Numeric       Cin,         // scaling factor
		     const String&       model,       // model option
		     ConstVectorView     f_grid,      // frequency vector
		     ConstVectorView     abs_p,       // pressure vector
		     ConstVectorView     abs_t,       // temperature vector
		     ConstVectorView     vmr,         // O2 vmr profile
		     ConstVectorView     abs_n2 );    // N2 vmr profile

void CKD_mt_v1v0_o2( MatrixView          xsec,        // calculated x-section
		     const Numeric       Cin,         // scaling factor
		     const String&       model,       // model option
		     ConstVectorView     f_grid,      // frequency vector
		     ConstVectorView     abs_p,       // pressure vector
		     ConstVectorView     abs_t,       // temperature vector
		     ConstVectorView     vmr );       // O2 vmr profile

//////////////////////////////////////////////////////////////////////////// 
// nitrogen continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void CKD_mt_CIArot_n2( MatrixView         xsec,        // calculated x-section
		       const Numeric      Cin,         // scaling factor
		       const String&      model,       // model option
		       ConstVectorView    f_grid,      // frequency vector
		       ConstVectorView    abs_p,       // pressure vector
		       ConstVectorView    abs_t,       // temperature vector
		       ConstVectorView    vmr );       // N2 vmr profile

void CKD_mt_CIAfun_n2( MatrixView         xsec,        // calculated x-section
		      const Numeric       Cin,         // scaling factor
		      const String&       model,       // model option
		      ConstVectorView     f_grid,      // frequency vector
		      ConstVectorView     abs_p,       // pressure vector
		      ConstVectorView     abs_t,       // temperature vector
		      ConstVectorView     vmr );       // N2 vmr profile

void BF86_CIA_N2( MatrixView              xsec,        // calculated x-section
		  const Numeric           Cin,         // model parameter
		  const String&           model,       // model option 
		  ConstVectorView         f_grid,      // frequency vector
		  ConstVectorView         abs_p,       // pressure vector
		  ConstVectorView         abs_t,       // temperature vector
		  ConstVectorView         vmr   );     // N2 vmr profile 

void MPM93_N2_continuum( MatrixView       xsec,        // calculated x-section
			 const Numeric    Cin,         // model parameter
			 const Numeric    Gin,         // model parameter
			 const Numeric    xTin,        // model parameter
			 const Numeric    xfin,        // model parameter
			 const String&    model,       // model option
			 ConstVectorView  f_grid,      // frequency vector
			 ConstVectorView  abs_p,       // pressure vector
			 ConstVectorView  abs_t,       // temperature vector
			 ConstVectorView  abs_h2o,     // H2O vmr profile
			 ConstVectorView  vmr	 );    // N2 vmr profile

void Rosenkranz_N2_self_continuum( MatrixView        xsec,        // calculated x-section
				   const Numeric     Cin,         // model parameter
				   const Numeric     xin,         // model parameter
				   const String&     model,       // model option
				   ConstVectorView   f_grid,      // frequency vector
				   ConstVectorView   abs_p,       // pressure vector
				   ConstVectorView   abs_t,       // temperature vector
				   ConstVectorView   vmr );       // N2 vmr profile

void Standard_N2_self_continuum(   MatrixView        xsec,        // calculated x-section
                                   const Numeric     Cin,         // model parameter
                                   const Numeric     xfin,        // model parameter
                                   const Numeric     xtin,        // model parameter
                                   const Numeric     xpin,        // model parameter
				   const String&     model,       // model option
				   ConstVectorView   f_grid,      // frequency vector
				   ConstVectorView   abs_p,       // pressure vector
				   ConstVectorView   abs_t,       // temperature vector
				   ConstVectorView   vmr );       // N2 vmr profile

void Pardo_ATM_N2_dry_continuum( MatrixView         xsec,         // calculated x-section
				  const Numeric     Cin,          // model parameter
				  const String&     model,        // model option
				  ConstVectorView   f_grid,       // frequency vector
				  ConstVectorView   abs_p,        // pressure vector
				  ConstVectorView   abs_t,        // temperature vector
				  ConstVectorView   vmr,	  // N2 vmr profile
				  ConstVectorView   abs_h2o);     // H2O vmr profile

//////////////////////////////////////////////////////////////////////////// 
// carbon dioxide continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void CKD_241_co2( MatrixView          xsec,        // calculated x-section
		  const Numeric       Cin,         // scaling factor
		  const String&       model,       // model option
		  ConstVectorView     f_grid,      // frequency vector
		  ConstVectorView     abs_p,       // pressure vector
		  ConstVectorView     abs_t,       // temperature vector
		  ConstVectorView     vmr );	   // CO2 vmr profile

void CKD_mt_co2( MatrixView          xsec,         // calculated x-section
		 const Numeric       Cin,          // scaling factor
		 const String&       model,        // model option
		 ConstVectorView     f_grid,       // frequency vector
		 ConstVectorView     abs_p,        // pressure vector
		 ConstVectorView     abs_t,        // temperature vector
		 ConstVectorView     vmr );        // CO2 vmr profile

void Rosenkranz_CO2_self_continuum( MatrixView        xsec,       // calculated x-section
				    const Numeric     C,          // model parameter
				    const Numeric     x,          // model parameter
				    const String&     model,      // model option
				    ConstVectorView   f_grid,     // frequency vector
				    ConstVectorView   abs_p,      // pressure vector
				    ConstVectorView   abs_t,      // temperature vector
				    ConstVectorView   vmr );      // CO2 vmr profile

void Rosenkranz_CO2_foreign_continuum( MatrixView        xsec,    // calculated x-section
				       const Numeric     C,       // model parameter
				       const Numeric     x,       // model parameter
				       const String&     model,   // model option
				       ConstVectorView   f_grid,  // frequency vector
				       ConstVectorView   abs_p,   // pressure vector
				       ConstVectorView   abs_t,   // temperature vector
				       ConstVectorView   abs_n2,  // N2 vmr profile
				       ConstVectorView   vmr );   // CO2 vmr profile

//////////////////////////////////////////////////////////////////////////// 
// water droplet and ice particle absorption (clouds)
//////////////////////////////////////////////////////////////////////////// 

void MPM93WaterDropletAbs( MatrixView        xsec,     // calculated x-section
			   const Numeric     CC,       // model parameter
			   const Numeric     CG,       // model parameter
			   const Numeric     CE,       // model parameter
			   const String&     model,    // model option
			   ConstVectorView   f_grid,   // frequency vector
			   ConstVectorView   abs_p,    // pressure vector
			   ConstVectorView   abs_t,    // temperature vector
			   ConstVectorView   vmr);     // suspended water droplet density vector

void MPM93IceCrystalAbs( MatrixView        xsec,       // calculated x-section
			 const Numeric     CC,         // model parameter
			 const Numeric     CA,         // model parameter
			 const Numeric     CB,         // model parameter
			 const String&     model,      // model option
			 ConstVectorView   f_grid,     // frequency vector
			 ConstVectorView   abs_p,      // pressure vector
			 ConstVectorView   abs_t,      // temperature vector
			 ConstVectorView   vmr	 );    // suspended ice particle density vector, 

void MPM93RainExt( MatrixView        xsec,       // calculated x-section
		   const Numeric     CE,         // model parameter
		   const Numeric     CA,         // model parameter
		   const Numeric     CB,         // model parameter
		   const String&     model,      // model option
		   ConstVectorView   f_grid,     // frequency vector
		   ConstVectorView   abs_p,      // pressure vector
		   ConstVectorView   abs_t,      // temperature vector
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




//////////////////////////////////////////////////////////////////////////// 
// arrays of the CKD H2O, CO2, N2, O2 absorption models
//////////////////////////////////////////////////////////////////////////// 

// additional array fields due to different numbering schemes of F77 and C/C++
const int addF77fields = 1;


/*
      ---------------------------------------------------------------------

c     29 August 2002
c
c     This  version, ckd_2.4.2, is composed of the updates to all continua 
c     including water vapor (ckd_2.2).  This is to retain consistency between the 
c     continuum and the HITRAN line parameters.
c

      further information can be found under
      http://www.rtweb.aer.com/continuum_frame.html

      Transformation from original F77 code to C/C++ by 
      T. Kuhn, iup Bremen, October 2003 
 
      ---------------------------------------------------------------------
 */
// H2O self continuum parameters at T=296 K for the CKD versions 0 to 2.4.2
// date of last update: 06/28/82
// units of (CM**3/MOL)*1.E-20
const Numeric SL296_ckd_0_v1  =   -20.0;
const Numeric SL296_ckd_0_v2  = 20000.0;
const Numeric SL296_ckd_0_dv  =    10.0;
const int     SL296_ckd_0_npt =  2003;
const double  SL296_ckd_0[SL296_ckd_0_npt+addF77fields] = {
  0.0000e-00, 1.1109e-01, 1.0573e-01, 
  1.0162e-01, 1.0573e-01, 1.1109e-01, 1.2574e-01, 1.3499e-01,  // f02740
  1.4327e-01, 1.5065e-01, 1.5164e-01, 1.5022e-01, 1.3677e-01,  // f02750
  1.3115e-01, 1.2253e-01, 1.1271e-01, 1.0070e-01, 8.7495e-02,  // f02760
  8.0118e-02, 6.9940e-02, 6.2034e-02, 5.6051e-02, 4.7663e-02,  // f02770
  4.2450e-02, 3.6690e-02, 3.3441e-02, 3.0711e-02, 2.5205e-02,  // f02780
  2.2113e-02, 1.8880e-02, 1.6653e-02, 1.4626e-02, 1.2065e-02,  // f02790
  1.0709e-02, 9.1783e-03, 7.7274e-03, 6.7302e-03, 5.6164e-03,  // f02800
  4.9089e-03, 4.1497e-03, 3.5823e-03, 3.1124e-03, 2.6414e-03,  // f02810
  2.3167e-03, 2.0156e-03, 1.7829e-03, 1.5666e-03, 1.3928e-03,  // f02820
  1.2338e-03, 1.0932e-03, 9.7939e-04, 8.8241e-04, 7.9173e-04,  // f02830
  7.1296e-04, 6.4179e-04, 5.8031e-04, 5.2647e-04, 4.7762e-04,  // f02850
  4.3349e-04, 3.9355e-04, 3.5887e-04, 3.2723e-04, 2.9919e-04,  // f02860
  2.7363e-04, 2.5013e-04, 2.2876e-04, 2.0924e-04, 1.9193e-04,  // f02870
  1.7618e-04, 1.6188e-04, 1.4891e-04, 1.3717e-04, 1.2647e-04,  // f02880
  1.1671e-04, 1.0786e-04, 9.9785e-05, 9.2350e-05, 8.5539e-05,  // f02890
  7.9377e-05, 7.3781e-05, 6.8677e-05, 6.3993e-05, 5.9705e-05,  // f02900
  5.5788e-05, 5.2196e-05, 4.8899e-05, 4.5865e-05, 4.3079e-05,  // f02910
  4.0526e-05, 3.8182e-05, 3.6025e-05, 3.4038e-05, 3.2203e-05,  // f02920
  3.0511e-05, 2.8949e-05, 2.7505e-05, 2.6170e-05, 2.4933e-05,  // f02930
  2.3786e-05, 2.2722e-05, 2.1736e-05, 2.0819e-05, 1.9968e-05,  // f02940
  1.9178e-05, 1.8442e-05, 1.7760e-05, 1.7127e-05, 1.6541e-05,  // f02960
  1.5997e-05, 1.5495e-05, 1.5034e-05, 1.4614e-05, 1.4230e-05,  // f02970
  1.3883e-05, 1.3578e-05, 1.3304e-05, 1.3069e-05, 1.2876e-05,  // f02980
  1.2732e-05, 1.2626e-05, 1.2556e-05, 1.2544e-05, 1.2604e-05,  // f02990
  1.2719e-05, 1.2883e-05, 1.3164e-05, 1.3581e-05, 1.4187e-05,  // f03000
  1.4866e-05, 1.5669e-05, 1.6717e-05, 1.8148e-05, 2.0268e-05,  // f03010
  2.2456e-05, 2.5582e-05, 2.9183e-05, 3.3612e-05, 3.9996e-05,  // f03020
  4.6829e-05, 5.5055e-05, 6.5897e-05, 7.5360e-05, 8.7213e-05,  // f03030
  1.0046e-04, 1.1496e-04, 1.2943e-04, 1.5049e-04, 1.6973e-04,  // f03040
  1.8711e-04, 2.0286e-04, 2.2823e-04, 2.6780e-04, 2.8766e-04,  // f03050
  3.1164e-04, 3.3640e-04, 3.6884e-04, 3.9159e-04, 3.8712e-04,  // f03070
  3.7433e-04, 3.4503e-04, 3.1003e-04, 2.8027e-04, 2.5253e-04,  // f03080
  2.3408e-04, 2.2836e-04, 2.4442e-04, 2.7521e-04, 2.9048e-04,  // f03090
  3.0489e-04, 3.2646e-04, 3.3880e-04, 3.3492e-04, 3.0987e-04,  // f03100
  2.9482e-04, 2.8711e-04, 2.6068e-04, 2.2683e-04, 1.9996e-04,  // f03110
  1.7788e-04, 1.6101e-04, 1.3911e-04, 1.2013e-04, 1.0544e-04,  // f03120
  9.4224e-05, 8.1256e-05, 7.3667e-05, 6.2233e-05, 5.5906e-05,  // f03130
  5.1619e-05, 4.5140e-05, 4.0273e-05, 3.3268e-05, 3.0258e-05,  // f03140
  2.6440e-05, 2.3103e-05, 2.0749e-05, 1.8258e-05, 1.6459e-05,  // f03150
  1.4097e-05, 1.2052e-05, 1.0759e-05, 9.1400e-06, 8.1432e-06,  // f03160
  7.1460e-06, 6.4006e-06, 5.6995e-06, 4.9372e-06, 4.4455e-06,  // f03180
  3.9033e-06, 3.4740e-06, 3.1269e-06, 2.8059e-06, 2.5558e-06,  // f03190
  2.2919e-06, 2.0846e-06, 1.8983e-06, 1.7329e-06, 1.5929e-06,  // f03200
  1.4631e-06, 1.3513e-06, 1.2461e-06, 1.1519e-06, 1.0682e-06,  // f03210
  9.9256e-07, 9.2505e-07, 8.6367e-07, 8.0857e-07, 7.5674e-07,  // f03220
  7.0934e-07, 6.6580e-07, 6.2580e-07, 5.8853e-07, 5.5333e-07,  // f03230
  5.2143e-07, 4.9169e-07, 4.6431e-07, 4.3898e-07, 4.1564e-07,  // f03240
  3.9405e-07, 3.7403e-07, 3.5544e-07, 3.3819e-07, 3.2212e-07,  // f03250
  3.0714e-07, 2.9313e-07, 2.8003e-07, 2.6777e-07, 2.5628e-07,  // f03260
  2.4551e-07, 2.3540e-07, 2.2591e-07, 2.1701e-07, 2.0866e-07,  // f03270
  2.0082e-07, 1.9349e-07, 1.8665e-07, 1.8027e-07, 1.7439e-07,  // f03290
  1.6894e-07, 1.6400e-07, 1.5953e-07, 1.5557e-07, 1.5195e-07,  // f03300
  1.4888e-07, 1.4603e-07, 1.4337e-07, 1.4093e-07, 1.3828e-07,  // f03310
  1.3569e-07, 1.3270e-07, 1.2984e-07, 1.2714e-07, 1.2541e-07,  // f03320
  1.2399e-07, 1.2102e-07, 1.1878e-07, 1.1728e-07, 1.1644e-07,  // f03330
  1.1491e-07, 1.1305e-07, 1.1235e-07, 1.1228e-07, 1.1224e-07,  // f03340
  1.1191e-07, 1.1151e-07, 1.1098e-07, 1.1068e-07, 1.1109e-07,  // f03350
  1.1213e-07, 1.1431e-07, 1.1826e-07, 1.2322e-07, 1.3025e-07,  // f03360
  1.4066e-07, 1.5657e-07, 1.7214e-07, 1.9449e-07, 2.2662e-07,  // f03370
  2.6953e-07, 3.1723e-07, 3.7028e-07, 4.4482e-07, 5.3852e-07,  // f03380
  6.2639e-07, 7.2175e-07, 7.7626e-07, 8.7248e-07, 9.6759e-07,  // f03400
  1.0102e-06, 1.0620e-06, 1.1201e-06, 1.2107e-06, 1.2998e-06,  // f03410
  1.3130e-06, 1.2856e-06, 1.2350e-06, 1.1489e-06, 1.0819e-06,  // f03420
  1.0120e-06, 9.4795e-07, 9.2858e-07, 9.8060e-07, 1.0999e-06,  // f03430
  1.1967e-06, 1.2672e-06, 1.3418e-06, 1.3864e-06, 1.4330e-06,  // f03440
  1.4592e-06, 1.4598e-06, 1.4774e-06, 1.4726e-06, 1.4820e-06,  // f03450
  1.5050e-06, 1.4984e-06, 1.5181e-06, 1.5888e-06, 1.6850e-06,  // f03460
  1.7690e-06, 1.9277e-06, 2.1107e-06, 2.3068e-06, 2.5347e-06,  // f03470
  2.8069e-06, 3.1345e-06, 3.5822e-06, 3.9051e-06, 4.3422e-06,  // f03480
  4.8704e-06, 5.5351e-06, 6.3454e-06, 7.2690e-06, 8.2974e-06,  // f03490
  9.7609e-06, 1.1237e-05, 1.3187e-05, 1.5548e-05, 1.8784e-05,  // f03510
  2.1694e-05, 2.5487e-05, 3.0092e-05, 3.5385e-05, 4.2764e-05,  // f03520
  4.9313e-05, 5.5800e-05, 6.2968e-05, 7.1060e-05, 7.7699e-05,  // f03530
  8.7216e-05, 8.9335e-05, 9.2151e-05, 9.2779e-05, 9.4643e-05,  // f03540
  9.7978e-05, 1.0008e-04, 1.0702e-04, 1.1026e-04, 1.0828e-04,  // f03550
  1.0550e-04, 1.0432e-04, 1.0428e-04, 9.8980e-05, 9.4992e-05,  // f03560
  9.5159e-05, 1.0058e-04, 1.0738e-04, 1.1550e-04, 1.1229e-04,  // f03570
  1.0596e-04, 1.0062e-04, 9.1742e-05, 8.4492e-05, 6.8099e-05,  // f03580
  5.6295e-05, 4.6502e-05, 3.8071e-05, 3.0721e-05, 2.3297e-05,  // f03590
  1.8688e-05, 1.4830e-05, 1.2049e-05, 9.6754e-06, 7.9192e-06,  // f03600
  6.6673e-06, 5.6468e-06, 4.8904e-06, 4.2289e-06, 3.6880e-06,  // f03620
  3.2396e-06, 2.8525e-06, 2.5363e-06, 2.2431e-06, 1.9949e-06,  // f03630
  1.7931e-06, 1.6164e-06, 1.4431e-06, 1.2997e-06, 1.1559e-06,  // f03640
  1.0404e-06, 9.4300e-07, 8.4597e-07, 7.6133e-07, 6.8623e-07,  // f03650
  6.2137e-07, 5.6345e-07, 5.1076e-07, 4.6246e-07, 4.1906e-07,  // f03660
  3.8063e-07, 3.4610e-07, 3.1554e-07, 2.8795e-07, 2.6252e-07,  // f03670
  2.3967e-07, 2.1901e-07, 2.0052e-07, 1.8384e-07, 1.6847e-07,  // f03680
  1.5459e-07, 1.4204e-07, 1.3068e-07, 1.2036e-07, 1.1095e-07,  // f03690
  1.0237e-07, 9.4592e-08, 8.7530e-08, 8.1121e-08, 7.5282e-08,  // f03700
  6.9985e-08, 6.5189e-08, 6.0874e-08, 5.6989e-08, 5.3530e-08,  // f03710
  5.0418e-08, 4.7745e-08, 4.5367e-08, 4.3253e-08, 4.1309e-08,  // f03730
  3.9695e-08, 3.8094e-08, 3.6482e-08, 3.4897e-08, 3.3500e-08,  // f03740
  3.2302e-08, 3.0854e-08, 2.9698e-08, 2.8567e-08, 2.7600e-08,  // f03750
  2.6746e-08, 2.5982e-08, 2.5510e-08, 2.5121e-08, 2.4922e-08,  // f03760
  2.4909e-08, 2.5013e-08, 2.5216e-08, 2.5589e-08, 2.6049e-08,  // f03770
  2.6451e-08, 2.6978e-08, 2.7687e-08, 2.8600e-08, 2.9643e-08,  // f03780
  3.0701e-08, 3.2058e-08, 3.3695e-08, 3.5558e-08, 3.7634e-08,  // f03790
  3.9875e-08, 4.2458e-08, 4.5480e-08, 4.8858e-08, 5.2599e-08,  // f03800
  5.7030e-08, 6.2067e-08, 6.7911e-08, 7.4579e-08, 8.1902e-08,  // f03810
  8.9978e-08, 9.9870e-08, 1.1102e-07, 1.2343e-07, 1.3732e-07,  // f03820
  1.5394e-07, 1.7318e-07, 1.9383e-07, 2.1819e-07, 2.4666e-07,  // f03840
  2.8109e-07, 3.2236e-07, 3.7760e-07, 4.4417e-07, 5.2422e-07,  // f03850
  6.1941e-07, 7.4897e-07, 9.2041e-07, 1.1574e-06, 1.4126e-06,  // f03860
  1.7197e-06, 2.1399e-06, 2.6266e-06, 3.3424e-06, 3.8418e-06,  // f03870
  4.5140e-06, 5.0653e-06, 5.8485e-06, 6.5856e-06, 6.8937e-06,  // f03880
  6.9121e-06, 6.9005e-06, 6.9861e-06, 6.8200e-06, 6.6089e-06,  // f03890
  6.5809e-06, 7.3496e-06, 8.0311e-06, 8.3186e-06, 8.4260e-06,  // f03900
  9.0644e-06, 9.4965e-06, 9.4909e-06, 9.0160e-06, 9.1494e-06,  // f03910
  9.3629e-06, 9.5944e-06, 9.5459e-06, 8.9919e-06, 8.6040e-06,  // f03920
  7.8613e-06, 7.1567e-06, 6.2677e-06, 5.1899e-06, 4.4188e-06,  // f03930
  3.7167e-06, 3.0636e-06, 2.5573e-06, 2.0317e-06, 1.6371e-06,  // f03950
  1.3257e-06, 1.0928e-06, 8.9986e-07, 7.4653e-07, 6.1111e-07,  // f03960
  5.1395e-07, 4.3500e-07, 3.7584e-07, 3.2633e-07, 2.8413e-07,  // f03970
  2.4723e-07, 2.1709e-07, 1.9294e-07, 1.7258e-07, 1.5492e-07,  // f03980
  1.3820e-07, 1.2389e-07, 1.1189e-07, 1.0046e-07, 9.0832e-08,  // f03990
  8.2764e-08, 7.4191e-08, 6.7085e-08, 6.0708e-08, 5.4963e-08,  // f04000
  4.9851e-08, 4.5044e-08, 4.0916e-08, 3.7220e-08, 3.3678e-08,  // f04010
  3.0663e-08, 2.7979e-08, 2.5495e-08, 2.3286e-08, 2.1233e-08,  // f04020
  1.9409e-08, 1.7770e-08, 1.6260e-08, 1.4885e-08, 1.3674e-08,  // f04030
  1.2543e-08, 1.1551e-08, 1.0655e-08, 9.8585e-09, 9.1398e-09,  // f04040
  8.4806e-09, 7.8899e-09, 7.3547e-09, 6.8670e-09, 6.4131e-09,  // f04060
  5.9930e-09, 5.6096e-09, 5.2592e-09, 4.9352e-09, 4.6354e-09,  // f04070
  4.3722e-09, 4.1250e-09, 3.9081e-09, 3.7118e-09, 3.5372e-09,  // f04080
  3.3862e-09, 3.2499e-09, 3.1324e-09, 3.0313e-09, 2.9438e-09,  // f04090
  2.8686e-09, 2.8050e-09, 2.7545e-09, 2.7149e-09, 2.6907e-09,  // f04100
  2.6724e-09, 2.6649e-09, 2.6642e-09, 2.6725e-09, 2.6871e-09,  // f04110
  2.7056e-09, 2.7357e-09, 2.7781e-09, 2.8358e-09, 2.9067e-09,  // f04120
  2.9952e-09, 3.1020e-09, 3.2253e-09, 3.3647e-09, 3.5232e-09,  // f04130
  3.7037e-09, 3.9076e-09, 4.1385e-09, 4.3927e-09, 4.6861e-09,  // f04140
  5.0238e-09, 5.4027e-09, 5.8303e-09, 6.3208e-09, 6.8878e-09,  // f04150
  7.5419e-09, 8.3130e-09, 9.1952e-09, 1.0228e-08, 1.1386e-08,  // f04170
  1.2792e-08, 1.4521e-08, 1.6437e-08, 1.8674e-08, 2.1160e-08,  // f04180
  2.4506e-08, 2.8113e-08, 3.2636e-08, 3.7355e-08, 4.2234e-08,  // f04190
  4.9282e-08, 5.7358e-08, 6.6743e-08, 7.8821e-08, 9.4264e-08,  // f04200
  1.1542e-07, 1.3684e-07, 1.6337e-07, 2.0056e-07, 2.3252e-07,  // f04210
  2.6127e-07, 2.9211e-07, 3.3804e-07, 3.7397e-07, 3.8205e-07,  // f04220
  3.8810e-07, 3.9499e-07, 3.9508e-07, 3.7652e-07, 3.5859e-07,  // f04230
  3.6198e-07, 3.7871e-07, 4.0925e-07, 4.2717e-07, 4.8241e-07,  // f04240
  5.2008e-07, 5.6530e-07, 5.9531e-07, 6.1994e-07, 6.5080e-07,  // f04250
  6.6355e-07, 6.9193e-07, 6.9930e-07, 7.3058e-07, 7.4678e-07,  // f04260
  7.9193e-07, 8.3627e-07, 9.1267e-07, 1.0021e-06, 1.1218e-06,  // f04280
  1.2899e-06, 1.4447e-06, 1.7268e-06, 2.0025e-06, 2.3139e-06,  // f04290
  2.5599e-06, 2.8920e-06, 3.3059e-06, 3.5425e-06, 3.9522e-06,  // f04300
  4.0551e-06, 4.2818e-06, 4.2892e-06, 4.4210e-06, 4.5614e-06,  // f04310
  4.6739e-06, 4.9482e-06, 5.1118e-06, 5.0986e-06, 4.9417e-06,  // f04320
  4.9022e-06, 4.8449e-06, 4.8694e-06, 4.8111e-06, 4.9378e-06,  // f04330
  5.3231e-06, 5.7362e-06, 6.2350e-06, 6.0951e-06, 5.7281e-06,  // f04340
  5.4585e-06, 4.9032e-06, 4.3009e-06, 3.4776e-06, 2.8108e-06,  // f04350
  2.2993e-06, 1.7999e-06, 1.3870e-06, 1.0750e-06, 8.5191e-07,  // f04360
  6.7951e-07, 5.5336e-07, 4.6439e-07, 4.0243e-07, 3.5368e-07,  // f04370
  3.1427e-07, 2.7775e-07, 2.4486e-07, 2.1788e-07, 1.9249e-07,  // f04390
  1.7162e-07, 1.5115e-07, 1.3478e-07, 1.2236e-07, 1.1139e-07,  // f04400
  1.0092e-07, 9.0795e-08, 8.2214e-08, 7.4691e-08, 6.7486e-08,  // f04410
  6.0414e-08, 5.4584e-08, 4.8754e-08, 4.3501e-08, 3.8767e-08,  // f04420
  3.4363e-08, 3.0703e-08, 2.7562e-08, 2.4831e-08, 2.2241e-08,  // f04430
  1.9939e-08, 1.8049e-08, 1.6368e-08, 1.4863e-08, 1.3460e-08,  // f04440
  1.2212e-08, 1.1155e-08, 1.0185e-08, 9.3417e-09, 8.5671e-09,  // f04450
  7.8292e-09, 7.1749e-09, 6.5856e-09, 6.0588e-09, 5.5835e-09,  // f04460
  5.1350e-09, 4.7395e-09, 4.3771e-09, 4.0476e-09, 3.7560e-09,  // f04470
  3.4861e-09, 3.2427e-09, 3.0240e-09, 2.8278e-09, 2.6531e-09,  // f04480
  2.4937e-09, 2.3511e-09, 2.2245e-09, 2.1133e-09, 2.0159e-09,  // f04500
  1.9330e-09, 1.8669e-09, 1.8152e-09, 1.7852e-09, 1.7752e-09,  // f04510
  1.7823e-09, 1.8194e-09, 1.8866e-09, 1.9759e-09, 2.0736e-09,  // f04520
  2.2083e-09, 2.3587e-09, 2.4984e-09, 2.6333e-09, 2.8160e-09,  // f04530
  3.0759e-09, 3.3720e-09, 3.6457e-09, 4.0668e-09, 4.4541e-09,  // f04540
  4.7976e-09, 5.0908e-09, 5.4811e-09, 6.1394e-09, 6.3669e-09,  // f04550
  6.5714e-09, 6.8384e-09, 7.1918e-09, 7.3741e-09, 7.2079e-09,  // f04560
  7.2172e-09, 7.2572e-09, 7.3912e-09, 7.6188e-09, 8.3291e-09,  // f04570
  8.7885e-09, 9.2412e-09, 1.0021e-08, 1.0752e-08, 1.1546e-08,  // f04580
  1.1607e-08, 1.1949e-08, 1.2346e-08, 1.2516e-08, 1.2826e-08,  // f04590
  1.3053e-08, 1.3556e-08, 1.4221e-08, 1.5201e-08, 1.6661e-08,  // f04610
  1.8385e-08, 2.0585e-08, 2.3674e-08, 2.7928e-08, 3.3901e-08,  // f04620
  4.1017e-08, 4.9595e-08, 6.0432e-08, 7.6304e-08, 9.0764e-08,  // f04630
  1.0798e-07, 1.2442e-07, 1.4404e-07, 1.6331e-07, 1.8339e-07,  // f04640
  2.0445e-07, 2.2288e-07, 2.3083e-07, 2.3196e-07, 2.3919e-07,  // f04650
  2.3339e-07, 2.3502e-07, 2.3444e-07, 2.6395e-07, 2.9928e-07,  // f04660
  3.0025e-07, 3.0496e-07, 3.1777e-07, 3.4198e-07, 3.4739e-07,  // f04670
  3.2696e-07, 3.4100e-07, 3.5405e-07, 3.7774e-07, 3.8285e-07,  // f04680
  3.6797e-07, 3.5800e-07, 3.2283e-07, 2.9361e-07, 2.4881e-07,  // f04690
  2.0599e-07, 1.7121e-07, 1.3641e-07, 1.1111e-07, 8.9413e-08,  // f04700
  7.3455e-08, 6.2078e-08, 5.2538e-08, 4.5325e-08, 3.9005e-08,  // f04720
  3.4772e-08, 3.1203e-08, 2.8132e-08, 2.5250e-08, 2.2371e-08,  // f04730
  2.0131e-08, 1.7992e-08, 1.6076e-08, 1.4222e-08, 1.2490e-08,  // f04740
  1.1401e-08, 1.0249e-08, 9.2279e-09, 8.5654e-09, 7.6227e-09,  // f04750
  6.9648e-09, 6.2466e-09, 5.7252e-09, 5.3800e-09, 4.6960e-09,  // f04760
  4.2194e-09, 3.7746e-09, 3.3813e-09, 3.0656e-09, 2.6885e-09,  // f04770
  2.4311e-09, 2.1572e-09, 1.8892e-09, 1.7038e-09, 1.4914e-09,  // f04780
  1.3277e-09, 1.1694e-09, 1.0391e-09, 9.2779e-10, 8.3123e-10,  // f04790
  7.4968e-10, 6.8385e-10, 6.2915e-10, 5.7784e-10, 5.2838e-10,  // f04800
  4.8382e-10, 4.4543e-10, 4.1155e-10, 3.7158e-10, 3.3731e-10,  // f04810
  3.0969e-10, 2.8535e-10, 2.6416e-10, 2.4583e-10, 2.2878e-10,  // f04830
  2.1379e-10, 2.0073e-10, 1.8907e-10, 1.7866e-10, 1.6936e-10,  // f04840
  1.6119e-10, 1.5424e-10, 1.4847e-10, 1.4401e-10, 1.4068e-10,  // f04850
  1.3937e-10, 1.3943e-10, 1.4281e-10, 1.4766e-10, 1.5701e-10,  // f04860
  1.7079e-10, 1.8691e-10, 2.0081e-10, 2.1740e-10, 2.4847e-10,  // f04870
  2.6463e-10, 2.7087e-10, 2.7313e-10, 2.8352e-10, 2.9511e-10,  // f04880
  2.8058e-10, 2.7227e-10, 2.7356e-10, 2.8012e-10, 2.8034e-10,  // f04890
  2.9031e-10, 3.1030e-10, 3.3745e-10, 3.8152e-10, 4.0622e-10,  // f04900
  4.2673e-10, 4.3879e-10, 4.5488e-10, 4.7179e-10, 4.6140e-10,  // f04910
  4.6339e-10, 4.6716e-10, 4.7024e-10, 4.7931e-10, 4.8503e-10,  // f04920
  4.9589e-10, 4.9499e-10, 5.0363e-10, 5.3184e-10, 5.6451e-10,  // f04940
  6.0932e-10, 6.6469e-10, 7.4076e-10, 8.3605e-10, 9.4898e-10,  // f04950
  1.0935e-09, 1.2593e-09, 1.4913e-09, 1.8099e-09, 2.1842e-09,  // f04960
  2.7284e-09, 3.2159e-09, 3.7426e-09, 4.5226e-09, 5.3512e-09,  // f04970
  6.1787e-09, 6.8237e-09, 7.9421e-09, 9.0002e-09, 9.6841e-09,  // f04980
  9.9558e-09, 1.0232e-08, 1.0591e-08, 1.0657e-08, 1.0441e-08,  // f04990
  1.0719e-08, 1.1526e-08, 1.2962e-08, 1.4336e-08, 1.6150e-08,  // f05000
  1.8417e-08, 2.0725e-08, 2.3426e-08, 2.5619e-08, 2.7828e-08,  // f05010
  3.0563e-08, 3.3438e-08, 3.6317e-08, 4.0400e-08, 4.4556e-08,  // f05020
  5.0397e-08, 5.3315e-08, 5.9185e-08, 6.5311e-08, 6.9188e-08,  // f05030
  7.7728e-08, 7.9789e-08, 8.6598e-08, 8.7768e-08, 9.1773e-08,  // f05050
  9.7533e-08, 1.0007e-07, 1.0650e-07, 1.0992e-07, 1.0864e-07,  // f05060
  1.0494e-07, 1.0303e-07, 1.0031e-07, 1.0436e-07, 1.0537e-07,  // f05070
  1.1184e-07, 1.2364e-07, 1.3651e-07, 1.4881e-07, 1.4723e-07,  // f05080
  1.4118e-07, 1.3371e-07, 1.1902e-07, 1.0007e-07, 7.9628e-08,  // f05090
  6.4362e-08, 5.0243e-08, 3.8133e-08, 2.9400e-08, 2.3443e-08,  // f05100
  1.9319e-08, 1.6196e-08, 1.4221e-08, 1.2817e-08, 1.1863e-08,  // f05110
  1.1383e-08, 1.1221e-08, 1.1574e-08, 1.1661e-08, 1.2157e-08,  // f05120
  1.2883e-08, 1.3295e-08, 1.4243e-08, 1.4240e-08, 1.4614e-08,  // f05130
  1.4529e-08, 1.4685e-08, 1.4974e-08, 1.4790e-08, 1.4890e-08,  // f05140
  1.4704e-08, 1.4142e-08, 1.3374e-08, 1.2746e-08, 1.2172e-08,  // f05160
  1.2336e-08, 1.2546e-08, 1.3065e-08, 1.4090e-08, 1.5215e-08,  // f05170
  1.6540e-08, 1.6144e-08, 1.5282e-08, 1.4358e-08, 1.2849e-08,  // f05180
  1.0998e-08, 8.6956e-09, 7.0881e-09, 5.5767e-09, 4.2792e-09,  // f05190
  3.2233e-09, 2.5020e-09, 1.9985e-09, 1.5834e-09, 1.3015e-09,  // f05200
  1.0948e-09, 9.4141e-10, 8.1465e-10, 7.1517e-10, 6.2906e-10,  // f05210
  5.5756e-10, 4.9805e-10, 4.3961e-10, 3.9181e-10, 3.5227e-10,  // f05220
  3.1670e-10, 2.8667e-10, 2.5745e-10, 2.3212e-10, 2.0948e-10,  // f05230
  1.8970e-10, 1.7239e-10, 1.5659e-10, 1.4301e-10, 1.3104e-10,  // f05240
  1.2031e-10, 1.1095e-10, 1.0262e-10, 9.5130e-11, 8.8595e-11,  // f05250
  8.2842e-11, 7.7727e-11, 7.3199e-11, 6.9286e-11, 6.5994e-11,  // f05270
  6.3316e-11, 6.1244e-11, 5.9669e-11, 5.8843e-11, 5.8832e-11,  // f05280
  5.9547e-11, 6.1635e-11, 6.4926e-11, 7.0745e-11, 7.8802e-11,  // f05290
  8.6724e-11, 1.0052e-10, 1.1575e-10, 1.3626e-10, 1.5126e-10,  // f05300
  1.6751e-10, 1.9239e-10, 2.1748e-10, 2.2654e-10, 2.2902e-10,  // f05310
  2.3240e-10, 2.4081e-10, 2.3930e-10, 2.2378e-10, 2.2476e-10,  // f05320
  2.2791e-10, 2.4047e-10, 2.5305e-10, 2.8073e-10, 3.1741e-10,  // f05330
  3.6592e-10, 4.1495e-10, 4.6565e-10, 5.0990e-10, 5.5607e-10,  // f05340
  6.1928e-10, 6.6779e-10, 7.3350e-10, 8.1434e-10, 8.9635e-10,  // f05350
  9.9678e-10, 1.1256e-09, 1.2999e-09, 1.4888e-09, 1.7642e-09,  // f05360
  1.9606e-09, 2.2066e-09, 2.4601e-09, 2.7218e-09, 3.0375e-09,  // f05380
  3.1591e-09, 3.2852e-09, 3.2464e-09, 3.3046e-09, 3.2710e-09,  // f05390
  3.2601e-09, 3.3398e-09, 3.7446e-09, 4.0795e-09, 4.0284e-09,  // f05400
  4.0584e-09, 4.1677e-09, 4.5358e-09, 4.4097e-09, 4.2744e-09,  // f05410
  4.5449e-09, 4.8147e-09, 5.2656e-09, 5.2476e-09, 5.0275e-09,  // f05420
  4.7968e-09, 4.3654e-09, 3.9530e-09, 3.2447e-09, 2.6489e-09,  // f05430
  2.1795e-09, 1.7880e-09, 1.4309e-09, 1.1256e-09, 9.1903e-10,  // f05440
  7.6533e-10, 6.3989e-10, 5.5496e-10, 4.9581e-10, 4.5722e-10,  // f05450
  4.3898e-10, 4.3505e-10, 4.3671e-10, 4.5329e-10, 4.6827e-10,  // f05460
  4.9394e-10, 5.1122e-10, 5.1649e-10, 5.0965e-10, 4.9551e-10,  // f05470
  4.8928e-10, 4.7947e-10, 4.7989e-10, 4.9071e-10, 4.8867e-10,  // f05490
  4.7260e-10, 4.5756e-10, 4.5400e-10, 4.5993e-10, 4.4042e-10,  // f05500
  4.3309e-10, 4.4182e-10, 4.6735e-10, 5.0378e-10, 5.2204e-10,  // f05510
  5.0166e-10, 4.6799e-10, 4.3119e-10, 3.8803e-10, 3.3291e-10,  // f05520
  2.6289e-10, 2.1029e-10, 1.7011e-10, 1.3345e-10, 1.0224e-10,  // f05530
  7.8207e-11, 6.2451e-11, 5.0481e-11, 4.1507e-11, 3.5419e-11,  // f05540
  3.0582e-11, 2.6900e-11, 2.3778e-11, 2.1343e-11, 1.9182e-11,  // f05550
  1.7162e-11, 1.5391e-11, 1.3877e-11, 1.2619e-11, 1.1450e-11,  // f05560
  1.0461e-11, 9.6578e-12, 8.9579e-12, 8.3463e-12, 7.8127e-12,  // f05570
  7.3322e-12, 6.9414e-12, 6.6037e-12, 6.3285e-12, 6.1095e-12,  // f05580
  5.9387e-12, 5.8118e-12, 5.7260e-12, 5.6794e-12, 5.6711e-12,  // f05600
  5.7003e-12, 5.7670e-12, 5.8717e-12, 6.0151e-12, 6.1984e-12,  // f05610
  6.4232e-12, 6.6918e-12, 7.0065e-12, 7.3705e-12, 7.7873e-12,  // f05620
  8.2612e-12, 8.7972e-12, 9.4009e-12, 1.0079e-11, 1.0840e-11,  // f05630
  1.1692e-11, 1.2648e-11, 1.3723e-11, 1.4935e-11, 1.6313e-11,  // f05640
  1.7905e-11, 1.9740e-11, 2.1898e-11, 2.4419e-11, 2.7426e-11,  // f05650
  3.0869e-11, 3.4235e-11, 3.7841e-11, 4.1929e-11, 4.6776e-11,  // f05660
  5.2123e-11, 5.8497e-11, 6.5294e-11, 7.4038e-11, 8.4793e-11,  // f05670
  9.6453e-11, 1.1223e-10, 1.2786e-10, 1.4882e-10, 1.7799e-10,  // f05680
  2.0766e-10, 2.4523e-10, 2.8591e-10, 3.3386e-10, 4.0531e-10,  // f05690
  4.7663e-10, 5.4858e-10, 6.3377e-10, 7.1688e-10, 8.4184e-10,  // f05710
  9.5144e-10, 1.0481e-09, 1.1356e-09, 1.2339e-09, 1.3396e-09,  // f05720
  1.4375e-09, 1.5831e-09, 1.7323e-09, 1.9671e-09, 2.2976e-09,  // f05730
  2.6679e-09, 3.0777e-09, 3.4321e-09, 3.8192e-09, 4.2711e-09,  // f05740
  4.4903e-09, 4.8931e-09, 5.2253e-09, 5.4040e-09, 5.6387e-09,  // f05750
  5.6704e-09, 6.0345e-09, 6.1079e-09, 6.2576e-09, 6.4039e-09,  // f05760
  6.3776e-09, 6.1878e-09, 5.8616e-09, 5.7036e-09, 5.5840e-09,  // f05770
  5.6905e-09, 5.8931e-09, 6.2478e-09, 6.8291e-09, 7.4528e-09,  // f05780
  7.6078e-09, 7.3898e-09, 6.7573e-09, 5.9827e-09, 5.0927e-09,  // f05790
  4.0099e-09, 3.1933e-09, 2.4296e-09, 1.8485e-09, 1.4595e-09,  // f05800
  1.2017e-09, 1.0164e-09, 8.7433e-10, 7.7108e-10, 7.0049e-10,  // f05820
  6.5291e-10, 6.1477e-10, 5.9254e-10, 5.8150e-10, 5.7591e-10,  // f05830
  5.8490e-10, 5.8587e-10, 5.9636e-10, 6.2408e-10, 6.5479e-10,  // f05840
  7.0480e-10, 7.2313e-10, 7.5524e-10, 8.0863e-10, 8.3386e-10,  // f05850
  9.2342e-10, 9.6754e-10, 1.0293e-09, 1.0895e-09, 1.1330e-09,  // f05860
  1.2210e-09, 1.2413e-09, 1.2613e-09, 1.2671e-09, 1.2225e-09,  // f05870
  1.1609e-09, 1.0991e-09, 1.0600e-09, 1.0570e-09, 1.0818e-09,  // f05880
  1.1421e-09, 1.2270e-09, 1.3370e-09, 1.4742e-09, 1.4946e-09,  // f05890
  1.4322e-09, 1.3210e-09, 1.1749e-09, 1.0051e-09, 7.8387e-10,  // f05900
  6.1844e-10, 4.6288e-10, 3.4164e-10, 2.5412e-10, 1.9857e-10,  // f05910
  1.5876e-10, 1.2966e-10, 1.0920e-10, 9.4811e-11, 8.3733e-11,  // f05930
  7.3906e-11, 6.7259e-11, 6.1146e-11, 5.7119e-11, 5.3546e-11,  // f05940
  4.8625e-11, 4.4749e-11, 4.1089e-11, 3.7825e-11, 3.4465e-11,  // f05950
  3.1018e-11, 2.8109e-11, 2.5610e-11, 2.2859e-11, 2.0490e-11,  // f05960
  1.8133e-11, 1.5835e-11, 1.3949e-11, 1.2295e-11, 1.0799e-11,  // f05970
  9.6544e-12, 8.7597e-12, 7.9990e-12, 7.3973e-12, 6.9035e-12,  // f05980
  6.4935e-12, 6.1195e-12, 5.8235e-12, 5.5928e-12, 5.4191e-12,  // f05990
  5.2993e-12, 5.2338e-12, 5.2272e-12, 5.2923e-12, 5.4252e-12,  // f06000
  5.6523e-12, 5.9433e-12, 6.3197e-12, 6.9016e-12, 7.5016e-12,  // f06010
  8.2885e-12, 9.4050e-12, 1.0605e-11, 1.2257e-11, 1.3622e-11,  // f06020
  1.5353e-11, 1.7543e-11, 1.9809e-11, 2.2197e-11, 2.4065e-11,  // f06040
  2.6777e-11, 2.9751e-11, 3.2543e-11, 3.5536e-11, 3.9942e-11,  // f06050
  4.6283e-11, 5.4556e-11, 6.5490e-11, 7.6803e-11, 9.0053e-11,  // f06060
  1.0852e-10, 1.2946e-10, 1.4916e-10, 1.7748e-10, 2.0073e-10,  // f06070
  2.2485e-10, 2.5114e-10, 2.7715e-10, 3.1319e-10, 3.3305e-10,  // f06080
  3.5059e-10, 3.5746e-10, 3.6311e-10, 3.7344e-10, 3.6574e-10,  // f06090
  3.7539e-10, 3.9434e-10, 4.3510e-10, 4.3340e-10, 4.2588e-10,  // f06100
  4.3977e-10, 4.6062e-10, 4.7687e-10, 4.6457e-10, 4.8578e-10,  // f06110
  5.2344e-10, 5.6752e-10, 5.8702e-10, 5.6603e-10, 5.3784e-10,  // f06120
  4.9181e-10, 4.3272e-10, 3.5681e-10, 2.8814e-10, 2.3320e-10,  // f06130
  1.8631e-10, 1.4587e-10, 1.1782e-10, 9.8132e-11, 8.2528e-11,  // f06150
  6.9174e-11, 6.1056e-11, 5.3459e-11, 4.7116e-11, 4.1878e-11,  // f06160
  3.8125e-11, 3.6347e-11, 3.5071e-11, 3.3897e-11, 3.3541e-11,  // f06170
  3.3563e-11, 3.5469e-11, 3.8111e-11, 3.8675e-11, 4.1333e-11,  // f06180
  4.3475e-11, 4.6476e-11, 4.9761e-11, 5.1380e-11, 5.4135e-11,  // f06190
  5.3802e-11, 5.5158e-11, 5.6864e-11, 5.9311e-11, 6.3827e-11,  // f06200
  6.7893e-11, 6.8230e-11, 6.6694e-11, 6.6018e-11, 6.4863e-11,  // f06210
  6.5893e-11, 6.3813e-11, 6.4741e-11, 6.8630e-11, 7.0255e-11,  // f06220
  7.0667e-11, 6.8810e-11, 6.4104e-11, 5.8136e-11, 4.7242e-11,  // f06230
  3.7625e-11, 3.1742e-11, 2.5581e-11, 1.8824e-11, 1.3303e-11,  // f06240
  9.6919e-12, 7.5353e-12, 6.0986e-12, 5.0742e-12, 4.3094e-12,  // f06260
  3.7190e-12, 3.2520e-12, 2.8756e-12, 2.5680e-12, 2.3139e-12,  // f06270
  2.1025e-12, 1.9257e-12, 1.7777e-12, 1.6539e-12, 1.5508e-12,  // f06280
  1.4657e-12, 1.3966e-12, 1.3417e-12, 1.2998e-12, 1.2700e-12,  // f06290
  1.2514e-12, 1.2437e-12, 1.2463e-12, 1.2592e-12, 1.2823e-12,  // f06300
  1.3157e-12, 1.3596e-12, 1.4144e-12, 1.4806e-12, 1.5588e-12,  // f06310
  1.6497e-12, 1.7544e-12, 1.8738e-12, 2.0094e-12, 2.1626e-12,  // f06320
  2.3354e-12, 2.5297e-12, 2.7483e-12, 2.9941e-12, 3.2708e-12,  // f06330
  3.5833e-12, 3.9374e-12, 4.3415e-12, 4.8079e-12, 5.3602e-12,  // f06340
  5.9816e-12, 6.7436e-12, 7.6368e-12, 8.6812e-12, 9.8747e-12,  // f06350
  1.1350e-11, 1.3181e-11, 1.5406e-11, 1.7868e-11, 2.0651e-11,  // f06370
  2.4504e-11, 2.9184e-11, 3.4159e-11, 3.9979e-11, 4.8704e-11,  // f06380
  5.7856e-11, 6.7576e-11, 7.9103e-11, 9.4370e-11, 1.1224e-10,  // f06390
  1.3112e-10, 1.5674e-10, 1.8206e-10, 2.0576e-10, 2.3187e-10,  // f06400
  2.7005e-10, 3.0055e-10, 3.3423e-10, 3.6956e-10, 3.8737e-10,  // f06410
  4.2630e-10, 4.5154e-10, 4.8383e-10, 5.3582e-10, 5.8109e-10,  // f06420
  6.3741e-10, 6.3874e-10, 6.3870e-10, 6.5818e-10, 6.5056e-10,  // f06430
  6.5291e-10, 6.3159e-10, 6.3984e-10, 6.4549e-10, 6.5444e-10,  // f06440
  6.7035e-10, 6.7665e-10, 6.9124e-10, 6.8451e-10, 6.9255e-10,  // f06450
  6.9923e-10, 7.0396e-10, 6.7715e-10, 6.0371e-10, 5.3774e-10,  // f06460
  4.6043e-10, 3.7635e-10, 2.9484e-10, 2.2968e-10, 1.8185e-10,  // f06480
  1.4191e-10, 1.1471e-10, 9.4790e-11, 7.9613e-11, 6.7989e-11,  // f06490
  5.9391e-11, 5.2810e-11, 4.7136e-11, 4.2618e-11, 3.8313e-11,  // f06500
  3.4686e-11, 3.1669e-11, 2.9110e-11, 2.6871e-11, 2.5074e-11,  // f06510
  2.4368e-11, 2.3925e-11, 2.4067e-11, 2.4336e-11, 2.4704e-11,  // f06520
  2.5823e-11, 2.7177e-11, 2.9227e-11, 3.1593e-11, 3.5730e-11,  // f06530
  4.0221e-11, 4.3994e-11, 4.8448e-11, 5.3191e-11, 5.8552e-11,  // f06540
  6.3458e-11, 6.6335e-11, 7.2457e-11, 7.9091e-11, 8.2234e-11,  // f06550
  8.7668e-11, 8.7951e-11, 9.2952e-11, 9.6157e-11, 9.5926e-11,  // f06560
  1.0120e-10, 1.0115e-10, 9.9577e-11, 9.6633e-11, 9.2891e-11,  // f06570
  9.3315e-11, 9.5584e-11, 1.0064e-10, 1.0509e-10, 1.1455e-10,  // f06590
  1.2443e-10, 1.2963e-10, 1.2632e-10, 1.1308e-10, 1.0186e-10,  // f06600
  8.5880e-11, 6.7863e-11, 5.1521e-11, 3.7780e-11, 2.8842e-11,  // f06610
  2.2052e-11, 1.7402e-11, 1.4406e-11, 1.1934e-11, 1.0223e-11,  // f06620
  8.9544e-12, 7.9088e-12, 7.0675e-12, 6.2222e-12, 5.6051e-12,  // f06630
  5.0502e-12, 4.5578e-12, 4.2636e-12, 3.9461e-12, 3.7599e-12,  // f06640
  3.5215e-12, 3.2467e-12, 3.0018e-12, 2.6558e-12, 2.3928e-12,  // f06650
  2.0707e-12, 1.7575e-12, 1.5114e-12, 1.2941e-12, 1.1004e-12,  // f06660
  9.5175e-13, 8.2894e-13, 7.3253e-13, 6.5551e-13, 5.9098e-13,  // f06670
  5.3548e-13, 4.8697e-13, 4.4413e-13, 4.0600e-13, 3.7188e-13,  // f06680
  3.4121e-13, 3.1356e-13, 2.8856e-13, 2.6590e-13, 2.4533e-13,  // f06700
  2.2663e-13, 2.0960e-13, 1.9407e-13, 1.7990e-13, 1.6695e-13,  // f06710
  1.5512e-13, 1.4429e-13, 1.3437e-13, 1.2527e-13, 1.1693e-13,  // f06720
  1.0927e-13, 1.0224e-13, 9.5767e-14, 8.9816e-14, 8.4335e-14,  // f06730
  7.9285e-14, 7.4626e-14, 7.0325e-14, 6.6352e-14, 6.2676e-14,  // f06740
  5.9274e-14, 5.6121e-14, 5.3195e-14, 5.0479e-14, 4.7953e-14,  // f06750
  4.5602e-14, 4.3411e-14, 4.1367e-14, 3.9456e-14, 3.7670e-14,  // f06760
  3.5996e-14, 3.4427e-14, 3.2952e-14, 3.1566e-14, 3.0261e-14,  // f06770
  2.9030e-14, 2.7868e-14, 2.6770e-14, 2.5730e-14, 2.4745e-14,  // f06780
  2.3809e-14, 2.2921e-14, 2.2076e-14, 2.1271e-14, 2.0504e-14,  // f06790
  1.9772e-14, 1.9073e-14, 1.8404e-14, 1.7764e-14, 1.7151e-14,  // f06810
  1.6564e-14, 1.6000e-14, 1.5459e-14, 1.4939e-14, 1.4439e-14,  // f06820
  1.3958e-14, 1.3495e-14, 1.3049e-14, 1.2620e-14, 1.2206e-14,  // f06830
  1.1807e-14, 1.1422e-14, 1.1050e-14, 1.0691e-14, 1.0345e-14,  // f06840
  1.0010e-14, 9.6870e-15, 9.3747e-15, 9.0727e-15, 8.7808e-15,  // f06850
  8.4986e-15, 8.2257e-15, 7.9617e-15, 7.7064e-15, 7.4594e-15,  // f06860
  7.2204e-15, 6.9891e-15, 6.7653e-15, 6.5488e-15, 6.3392e-15,  // f06870
  6.1363e-15, 5.9399e-15, 5.7499e-15, 5.5659e-15, 5.3878e-15,  // f06880
  5.2153e-15, 5.0484e-15, 4.8868e-15, 4.7303e-15, 4.5788e-15,  // f06890
  4.4322e-15, 4.2902e-15, 4.1527e-15, 4.0196e-15, 3.8907e-15,  // f06900
  3.7659e-15, 3.6451e-15, 3.5281e-15, 3.4149e-15, 3.3052e-15,  // f06920
  3.1991e-15, 3.0963e-15, 2.9967e-15, 2.9004e-15, 2.8071e-15,  // f06930
  2.7167e-15, 2.6293e-15, 2.5446e-15, 2.4626e-15, 2.3833e-15,  // f06940
  2.3064e-15, 2.2320e-15, 2.1600e-15, 2.0903e-15, 2.0228e-15,  // f06950
  1.9574e-15, 1.8942e-15, 1.8329e-15, 1.7736e-15, 1.7163e-15,  // f06960
  1.6607e-15, 1.6069e-15, 1.5548e-15, 1.5044e-15, 1.4557e-15,  // f06970
  1.4084e-15, 1.3627e-15, 1.3185e-15, 1.2757e-15, 1.2342e-15,  // f06980
  1.1941e-15, 1.1552e-15, 1.1177e-15, 1.0813e-15, 1.0461e-15,  // f06990
  1.0120e-15, 9.7900e-16, 9.4707e-16, 9.1618e-16, 8.8628e-16,  // f07000
  8.5734e-16, 8.2933e-16, 8.0223e-16, 7.7600e-16, 7.5062e-16,  // f07010
  7.2606e-16, 7.0229e-16, 6.7929e-16, 6.5703e-16, 6.3550e-16,  // f07030
  6.1466e-16, 5.9449e-16, 5.7498e-16, 5.5610e-16, 5.3783e-16,  // f07040
  5.2015e-16, 5.0305e-16, 4.8650e-16, 4.7049e-16, 4.5500e-16,  // f07050
  4.4002e-16, 4.2552e-16, 4.1149e-16, 3.9792e-16, 3.8479e-16,  // f07060
  3.7209e-16, 3.5981e-16, 3.4792e-16, 3.3642e-16, 3.2530e-16,  // f07070
  3.1454e-16, 3.0413e-16, 2.9406e-16, 2.8432e-16, 2.7490e-16,  // f07080
  2.6579e-16, 2.5697e-16, 2.4845e-16, 2.4020e-16, 2.3223e-16,  // f07090
  2.2451e-16, 2.1705e-16, 2.0984e-16, 2.0286e-16, 1.9611e-16,  // f07100
  1.8958e-16, 1.8327e-16, 1.7716e-16, 1.7126e-16, 1.6555e-16,  // f07110
  1.6003e-16, 1.5469e-16, 1.4952e-16, 1.4453e-16, 1.3970e-16,  // f07120
  1.3503e-16};



// H2O foreign continuum parameters at T=296 K for the CKD versions 0 to 2.4.2
// date of last update: 06/28/82
// units of (CM**3/MOL)*1.E-20
const Numeric FH2O_ckd_0_v1  =   -20.0;
const Numeric FH2O_ckd_0_v2  = 20000.0;
const Numeric FH2O_ckd_0_dv  =    10.0;
const int     FH2O_ckd_0_npt =  2003;
const double  FH2O_ckd_0[FH2O_ckd_0_npt+addF77fields] = {
  0.0000e-00, 1.2859E-02, 1.1715E-02, 
  1.1038E-02, 1.1715E-02, 1.2859E-02, 1.5326E-02, 1.6999E-02,  // F12660
  1.8321E-02, 1.9402E-02, 1.9570E-02, 1.9432E-02, 1.7572E-02,  // F12670
  1.6760E-02, 1.5480E-02, 1.3984E-02, 1.2266E-02, 1.0467E-02,  // F12680
  9.4526E-03, 8.0485E-03, 6.9484E-03, 6.1416E-03, 5.0941E-03,  // F12690
  4.4836E-03, 3.8133E-03, 3.4608E-03, 3.1487E-03, 2.4555E-03,  // F12700
  2.0977E-03, 1.7266E-03, 1.4920E-03, 1.2709E-03, 9.8081E-04,  // F12710
  8.5063E-04, 6.8822E-04, 5.3809E-04, 4.4679E-04, 3.3774E-04,  // F12720
  2.7979E-04, 2.1047E-04, 1.6511E-04, 1.2993E-04, 9.3033E-05,  // F12730
  7.4360E-05, 5.6428E-05, 4.5442E-05, 3.4575E-05, 2.7903E-05,  // F12740
  2.1374E-05, 1.6075E-05, 1.3022E-05, 1.0962E-05, 8.5959E-06,  // F12750
  6.9125E-06, 5.3808E-06, 4.3586E-06, 3.6394E-06, 2.9552E-06,  // F12770
  2.3547E-06, 1.8463E-06, 1.6036E-06, 1.3483E-06, 1.1968E-06,  // F12780
  1.0333E-06, 8.4484E-07, 6.7195E-07, 5.0947E-07, 4.2343E-07,  // F12790
  3.4453E-07, 2.7830E-07, 2.3063E-07, 1.9951E-07, 1.7087E-07,  // F12800
  1.4393E-07, 1.2575E-07, 1.0750E-07, 8.2325E-08, 5.7524E-08,  // F12810
  4.4482E-08, 3.8106E-08, 3.4315E-08, 2.9422E-08, 2.5069E-08,  // F12820
  2.2402E-08, 1.9349E-08, 1.6152E-08, 1.2208E-08, 8.9660E-09,  // F12830
  7.1322E-09, 6.1028E-09, 5.2938E-09, 4.5350E-09, 3.4977E-09,  // F12840
  2.9511E-09, 2.4734E-09, 2.0508E-09, 1.8507E-09, 1.6373E-09,  // F12850
  1.5171E-09, 1.3071E-09, 1.2462E-09, 1.2148E-09, 1.2590E-09,  // F12860
  1.3153E-09, 1.3301E-09, 1.4483E-09, 1.6944E-09, 2.0559E-09,  // F12880
  2.2954E-09, 2.6221E-09, 3.2606E-09, 4.2392E-09, 5.2171E-09,  // F12890
  6.2553E-09, 8.2548E-09, 9.5842E-09, 1.1280E-08, 1.3628E-08,  // F12900
  1.7635E-08, 2.1576E-08, 2.4835E-08, 3.0014E-08, 3.8485E-08,  // F12910
  4.7440E-08, 5.5202E-08, 7.0897E-08, 9.6578E-08, 1.3976E-07,  // F12920
  1.8391E-07, 2.3207E-07, 2.9960E-07, 4.0408E-07, 5.9260E-07,  // F12930
  7.8487E-07, 1.0947E-06, 1.4676E-06, 1.9325E-06, 2.6587E-06,  // F12940
  3.4534E-06, 4.4376E-06, 5.8061E-06, 7.0141E-06, 8.4937E-06,  // F12950
  1.0186E-05, 1.2034E-05, 1.3837E-05, 1.6595E-05, 1.9259E-05,  // F12960
  2.1620E-05, 2.3681E-05, 2.7064E-05, 3.2510E-05, 3.5460E-05,  // F12970
  3.9109E-05, 4.2891E-05, 4.7757E-05, 5.0981E-05, 5.0527E-05,  // F12990
  4.8618E-05, 4.4001E-05, 3.7982E-05, 3.2667E-05, 2.7794E-05,  // F13000
  2.4910E-05, 2.4375E-05, 2.7316E-05, 3.2579E-05, 3.5499E-05,  // F13010
  3.8010E-05, 4.1353E-05, 4.3323E-05, 4.3004E-05, 3.9790E-05,  // F13020
  3.7718E-05, 3.6360E-05, 3.2386E-05, 2.7409E-05, 2.3626E-05,  // F13030
  2.0631E-05, 1.8371E-05, 1.5445E-05, 1.2989E-05, 1.1098E-05,  // F13040
  9.6552E-06, 8.0649E-06, 7.2365E-06, 5.9137E-06, 5.2759E-06,  // F13050
  4.8860E-06, 4.1321E-06, 3.5918E-06, 2.7640E-06, 2.4892E-06,  // F13060
  2.1018E-06, 1.7848E-06, 1.5855E-06, 1.3569E-06, 1.1986E-06,  // F13070
  9.4693E-07, 7.4097E-07, 6.3443E-07, 4.8131E-07, 4.0942E-07,  // F13080
  3.3316E-07, 2.8488E-07, 2.3461E-07, 1.7397E-07, 1.4684E-07,  // F13100
  1.0953E-07, 8.5396E-08, 6.9261E-08, 5.4001E-08, 4.5430E-08,  // F13110
  3.2791E-08, 2.5995E-08, 2.0225E-08, 1.5710E-08, 1.3027E-08,  // F13120
  1.0229E-08, 8.5277E-09, 6.5249E-09, 5.0117E-09, 3.9906E-09,  // F13130
  3.2332E-09, 2.7847E-09, 2.4570E-09, 2.3359E-09, 2.0599E-09,  // F13140
  1.8436E-09, 1.6559E-09, 1.4910E-09, 1.2794E-09, 9.8229E-10,  // F13150
  8.0054E-10, 6.0769E-10, 4.5646E-10, 3.3111E-10, 2.4428E-10,  // F13160
  1.8007E-10, 1.3291E-10, 9.7974E-11, 7.8271E-11, 6.3833E-11,  // F13170
  5.4425E-11, 4.6471E-11, 4.0209E-11, 3.5227E-11, 3.1212E-11,  // F13180
  2.8840E-11, 2.7762E-11, 2.7935E-11, 3.2012E-11, 3.9525E-11,  // F13190
  5.0303E-11, 6.8027E-11, 9.3954E-11, 1.2986E-10, 1.8478E-10,  // F13210
  2.5331E-10, 3.4827E-10, 4.6968E-10, 6.2380E-10, 7.9106E-10,  // F13220
  1.0026E-09, 1.2102E-09, 1.4146E-09, 1.6154E-09, 1.7510E-09,  // F13230
  1.8575E-09, 1.8742E-09, 1.8700E-09, 1.8582E-09, 1.9657E-09,  // F13240
  2.1204E-09, 2.0381E-09, 2.0122E-09, 2.0436E-09, 2.1213E-09,  // F13250
  2.0742E-09, 1.9870E-09, 2.0465E-09, 2.1556E-09, 2.2222E-09,  // F13260
  2.1977E-09, 2.1047E-09, 1.9334E-09, 1.7357E-09, 1.5754E-09,  // F13270
  1.4398E-09, 1.4018E-09, 1.5459E-09, 1.7576E-09, 2.1645E-09,  // F13280
  2.9480E-09, 4.4439E-09, 5.8341E-09, 8.0757E-09, 1.1658E-08,  // F13290
  1.6793E-08, 2.2694E-08, 2.9468E-08, 3.9278E-08, 5.2145E-08,  // F13300
  6.4378E-08, 7.7947E-08, 8.5321E-08, 9.7848E-08, 1.0999E-07,  // F13320
  1.1489E-07, 1.2082E-07, 1.2822E-07, 1.4053E-07, 1.5238E-07,  // F13330
  1.5454E-07, 1.5018E-07, 1.4048E-07, 1.2359E-07, 1.0858E-07,  // F13340
  9.3486E-08, 8.1638E-08, 7.7690E-08, 8.4625E-08, 1.0114E-07,  // F13350
  1.1430E-07, 1.2263E-07, 1.3084E-07, 1.3380E-07, 1.3573E-07,  // F13360
  1.3441E-07, 1.2962E-07, 1.2638E-07, 1.1934E-07, 1.1371E-07,  // F13370
  1.0871E-07, 9.8843E-08, 9.1877E-08, 9.1050E-08, 9.3213E-08,  // F13380
  9.2929E-08, 1.0155E-07, 1.1263E-07, 1.2370E-07, 1.3636E-07,  // F13390
  1.5400E-07, 1.7656E-07, 2.1329E-07, 2.3045E-07, 2.5811E-07,  // F13400
  2.9261E-07, 3.4259E-07, 4.0770E-07, 4.8771E-07, 5.8081E-07,  // F13410
  7.2895E-07, 8.7482E-07, 1.0795E-06, 1.3384E-06, 1.7208E-06,  // F13430
  2.0677E-06, 2.5294E-06, 3.1123E-06, 3.7900E-06, 4.7752E-06,  // F13440
  5.6891E-06, 6.6261E-06, 7.6246E-06, 8.7730E-06, 9.6672E-06,  // F13450
  1.0980E-05, 1.1287E-05, 1.1670E-05, 1.1635E-05, 1.1768E-05,  // F13460
  1.2039E-05, 1.2253E-05, 1.3294E-05, 1.4005E-05, 1.3854E-05,  // F13470
  1.3420E-05, 1.3003E-05, 1.2645E-05, 1.1715E-05, 1.1258E-05,  // F13480
  1.1516E-05, 1.2494E-05, 1.3655E-05, 1.4931E-05, 1.4649E-05,  // F13490
  1.3857E-05, 1.3120E-05, 1.1791E-05, 1.0637E-05, 8.2760E-06,  // F13500
  6.5821E-06, 5.1959E-06, 4.0158E-06, 3.0131E-06, 2.0462E-06,  // F13510
  1.4853E-06, 1.0365E-06, 7.3938E-07, 4.9752E-07, 3.4148E-07,  // F13520
  2.4992E-07, 1.8363E-07, 1.4591E-07, 1.1380E-07, 9.0588E-08,  // F13540
  7.3697E-08, 6.0252E-08, 5.1868E-08, 4.2660E-08, 3.6163E-08,  // F13550
  3.2512E-08, 2.9258E-08, 2.4238E-08, 2.1209E-08, 1.6362E-08,  // F13560
  1.3871E-08, 1.2355E-08, 9.6940E-09, 7.7735E-09, 6.2278E-09,  // F13570
  5.2282E-09, 4.3799E-09, 3.5545E-09, 2.7527E-09, 2.0950E-09,  // F13580
  1.6344E-09, 1.2689E-09, 1.0403E-09, 8.4880E-10, 6.3461E-10,  // F13590
  4.7657E-10, 3.5220E-10, 2.7879E-10, 2.3021E-10, 1.6167E-10,  // F13600
  1.1732E-10, 8.9206E-11, 7.0596E-11, 5.8310E-11, 4.4084E-11,  // F13610
  3.1534E-11, 2.5068E-11, 2.2088E-11, 2.2579E-11, 2.2637E-11,  // F13620
  2.5705E-11, 3.2415E-11, 4.6116E-11, 6.5346E-11, 9.4842E-11,  // F13630
  1.2809E-10, 1.8211E-10, 2.4052E-10, 3.0270E-10, 3.5531E-10,  // F13650
  4.2402E-10, 4.6730E-10, 4.7942E-10, 4.6813E-10, 4.5997E-10,  // F13660
  4.5788E-10, 4.0311E-10, 3.7367E-10, 3.3149E-10, 2.9281E-10,  // F13670
  2.5231E-10, 2.1152E-10, 1.9799E-10, 1.8636E-10, 1.9085E-10,  // F13680
  2.0786E-10, 2.2464E-10, 2.3785E-10, 2.5684E-10, 2.7499E-10,  // F13690
  2.6962E-10, 2.6378E-10, 2.6297E-10, 2.6903E-10, 2.7035E-10,  // F13700
  2.5394E-10, 2.5655E-10, 2.7184E-10, 2.9013E-10, 3.0585E-10,  // F13710
  3.0791E-10, 3.1667E-10, 3.4343E-10, 3.7365E-10, 4.0269E-10,  // F13720
  4.7260E-10, 5.6584E-10, 6.9791E-10, 8.6569E-10, 1.0393E-09,  // F13730
  1.2067E-09, 1.5047E-09, 1.8583E-09, 2.2357E-09, 2.6498E-09,  // F13740
  3.2483E-09, 3.9927E-09, 4.6618E-09, 5.5555E-09, 6.6609E-09,  // F13760
  8.2139E-09, 1.0285E-08, 1.3919E-08, 1.8786E-08, 2.5150E-08,  // F13770
  3.3130E-08, 4.5442E-08, 6.3370E-08, 9.0628E-08, 1.2118E-07,  // F13780
  1.5927E-07, 2.1358E-07, 2.7825E-07, 3.7671E-07, 4.4894E-07,  // F13790
  5.4442E-07, 6.2240E-07, 7.3004E-07, 8.3384E-07, 8.7933E-07,  // F13800
  8.8080E-07, 8.6939E-07, 8.6541E-07, 8.2055E-07, 7.7278E-07,  // F13810
  7.5989E-07, 8.6909E-07, 9.7945E-07, 1.0394E-06, 1.0646E-06,  // F13820
  1.1509E-06, 1.2017E-06, 1.1915E-06, 1.1259E-06, 1.1549E-06,  // F13830
  1.1938E-06, 1.2356E-06, 1.2404E-06, 1.1716E-06, 1.1149E-06,  // F13840
  1.0073E-06, 8.9845E-07, 7.6639E-07, 6.1517E-07, 5.0887E-07,  // F13850
  4.1269E-07, 3.2474E-07, 2.5698E-07, 1.8893E-07, 1.4009E-07,  // F13870
  1.0340E-07, 7.7724E-08, 5.7302E-08, 4.2178E-08, 2.9603E-08,  // F13880
  2.1945E-08, 1.6301E-08, 1.2806E-08, 1.0048E-08, 7.8970E-09,  // F13890
  6.1133E-09, 4.9054E-09, 4.1985E-09, 3.6944E-09, 3.2586E-09,  // F13900
  2.7362E-09, 2.3647E-09, 2.1249E-09, 1.8172E-09, 1.6224E-09,  // F13910
  1.5158E-09, 1.2361E-09, 1.0682E-09, 9.2312E-10, 7.9220E-10,  // F13920
  6.8174E-10, 5.6147E-10, 4.8268E-10, 4.1534E-10, 3.3106E-10,  // F13930
  2.8275E-10, 2.4584E-10, 2.0742E-10, 1.7840E-10, 1.4664E-10,  // F13940
  1.2390E-10, 1.0497E-10, 8.5038E-11, 6.7008E-11, 5.6355E-11,  // F13950
  4.3323E-11, 3.6914E-11, 3.2262E-11, 3.0749E-11, 3.0318E-11,  // F13960
  2.9447E-11, 2.9918E-11, 3.0668E-11, 3.1315E-11, 3.0329E-11,  // F13980
  2.8259E-11, 2.6065E-11, 2.3578E-11, 2.0469E-11, 1.6908E-11,  // F13990
  1.4912E-11, 1.1867E-11, 9.9730E-12, 8.1014E-12, 6.7528E-12,  // F14000
  6.3133E-12, 5.8599E-12, 6.0145E-12, 6.5105E-12, 7.0537E-12,  // F14010
  7.4973E-12, 7.8519E-12, 8.5039E-12, 9.1995E-12, 1.0694E-11,  // F14020
  1.1659E-11, 1.2685E-11, 1.3087E-11, 1.3222E-11, 1.2634E-11,  // F14030
  1.1077E-11, 9.6259E-12, 8.3202E-12, 7.4857E-12, 6.8069E-12,  // F14040
  6.7496E-12, 7.3116E-12, 8.0171E-12, 8.6394E-12, 9.2659E-12,  // F14050
  1.0048E-11, 1.0941E-11, 1.2226E-11, 1.3058E-11, 1.5193E-11,  // F14060
  1.8923E-11, 2.3334E-11, 2.8787E-11, 3.6693E-11, 4.8295E-11,  // F14070
  6.4260E-11, 8.8269E-11, 1.1865E-10, 1.5961E-10, 2.0605E-10,  // F14090
  2.7349E-10, 3.7193E-10, 4.8216E-10, 6.1966E-10, 7.7150E-10,  // F14100
  1.0195E-09, 1.2859E-09, 1.6535E-09, 2.0316E-09, 2.3913E-09,  // F14110
  3.0114E-09, 3.7495E-09, 4.6504E-09, 5.9145E-09, 7.6840E-09,  // F14120
  1.0304E-08, 1.3010E-08, 1.6441E-08, 2.1475E-08, 2.5892E-08,  // F14130
  2.9788E-08, 3.3820E-08, 4.0007E-08, 4.4888E-08, 4.5765E-08,  // F14140
  4.6131E-08, 4.6239E-08, 4.4849E-08, 4.0729E-08, 3.6856E-08,  // F14150
  3.6164E-08, 3.7606E-08, 4.1457E-08, 4.3750E-08, 5.1150E-08,  // F14160
  5.6054E-08, 6.1586E-08, 6.4521E-08, 6.6494E-08, 6.9024E-08,  // F14170
  6.8893E-08, 7.0901E-08, 6.9760E-08, 7.1485E-08, 7.0740E-08,  // F14180
  7.3764E-08, 7.6618E-08, 8.4182E-08, 9.3838E-08, 1.0761E-07,  // F14200
  1.2851E-07, 1.4748E-07, 1.8407E-07, 2.2109E-07, 2.6392E-07,  // F14210
  2.9887E-07, 3.4493E-07, 4.0336E-07, 4.3551E-07, 4.9231E-07,  // F14220
  5.0728E-07, 5.3781E-07, 5.3285E-07, 5.4496E-07, 5.5707E-07,  // F14230
  5.6944E-07, 6.1123E-07, 6.4317E-07, 6.4581E-07, 6.1999E-07,  // F14240
  6.0191E-07, 5.7762E-07, 5.7241E-07, 5.7013E-07, 6.0160E-07,  // F14250
  6.6905E-07, 7.4095E-07, 8.2121E-07, 8.0947E-07, 7.6145E-07,  // F14260
  7.2193E-07, 6.3722E-07, 5.4316E-07, 4.2186E-07, 3.2528E-07,  // F14270
  2.5207E-07, 1.8213E-07, 1.2658E-07, 8.6746E-08, 6.0216E-08,  // F14280
  4.1122E-08, 2.8899E-08, 2.1740E-08, 1.7990E-08, 1.5593E-08,  // F14290
  1.3970E-08, 1.2238E-08, 1.0539E-08, 9.2386E-09, 7.8481E-09,  // F14310
  6.8704E-09, 5.7615E-09, 5.0434E-09, 4.6886E-09, 4.3770E-09,  // F14320
  3.9768E-09, 3.5202E-09, 3.1854E-09, 2.9009E-09, 2.5763E-09,  // F14330
  2.2135E-09, 1.9455E-09, 1.6248E-09, 1.3368E-09, 1.0842E-09,  // F14340
  8.4254E-10, 6.7414E-10, 5.4667E-10, 4.5005E-10, 3.4932E-10,  // F14350
  2.6745E-10, 2.2053E-10, 1.8162E-10, 1.4935E-10, 1.1618E-10,  // F14360
  9.1888E-11, 8.0672E-11, 6.8746E-11, 6.2668E-11, 5.5715E-11,  // F14370
  4.5074E-11, 3.7669E-11, 3.2082E-11, 2.8085E-11, 2.4838E-11,  // F14380
  1.9791E-11, 1.6964E-11, 1.3887E-11, 1.1179E-11, 9.7499E-12,  // F14390
  7.8255E-12, 6.3698E-12, 5.3265E-12, 4.6588E-12, 4.4498E-12,  // F14400
  3.9984E-12, 3.7513E-12, 3.7176E-12, 3.9148E-12, 4.2702E-12,  // F14420
  5.0090E-12, 6.5801E-12, 8.7787E-12, 1.2718E-11, 1.8375E-11,  // F14430
  2.5304E-11, 3.5403E-11, 4.8842E-11, 6.4840E-11, 8.0911E-11,  // F14440
  1.0136E-10, 1.2311E-10, 1.4203E-10, 1.5869E-10, 1.8093E-10,  // F14450
  2.1370E-10, 2.5228E-10, 2.8816E-10, 3.4556E-10, 3.9860E-10,  // F14460
  4.4350E-10, 4.7760E-10, 5.2357E-10, 6.0827E-10, 6.3635E-10,  // F14470
  6.5886E-10, 6.8753E-10, 7.2349E-10, 7.2789E-10, 6.8232E-10,  // F14480
  6.6081E-10, 6.4232E-10, 6.3485E-10, 6.4311E-10, 7.2235E-10,  // F14490
  7.7263E-10, 8.1668E-10, 9.0324E-10, 9.7643E-10, 1.0535E-09,  // F14500
  1.0195E-09, 1.0194E-09, 1.0156E-09, 9.6792E-10, 9.2725E-10,  // F14510
  8.7347E-10, 8.4484E-10, 8.2647E-10, 8.4363E-10, 9.1261E-10,  // F14530
  1.0051E-09, 1.1511E-09, 1.4037E-09, 1.8066E-09, 2.4483E-09,  // F14540
  3.2739E-09, 4.3194E-09, 5.6902E-09, 7.7924E-09, 9.7376E-09,  // F14550
  1.2055E-08, 1.4303E-08, 1.6956E-08, 1.9542E-08, 2.2233E-08,  // F14560
  2.5186E-08, 2.7777E-08, 2.8943E-08, 2.8873E-08, 2.9417E-08,  // F14570
  2.7954E-08, 2.7524E-08, 2.7040E-08, 3.1254E-08, 3.6843E-08,  // F14580
  3.7797E-08, 3.8713E-08, 4.0135E-08, 4.2824E-08, 4.3004E-08,  // F14590
  4.0279E-08, 4.2781E-08, 4.5220E-08, 4.8948E-08, 5.0172E-08,  // F14600
  4.8499E-08, 4.7182E-08, 4.2204E-08, 3.7701E-08, 3.0972E-08,  // F14610
  2.4654E-08, 1.9543E-08, 1.4609E-08, 1.1171E-08, 8.3367E-09,  // F14620
  6.3791E-09, 5.0790E-09, 4.0655E-09, 3.3658E-09, 2.7882E-09,  // F14640
  2.4749E-09, 2.2287E-09, 2.0217E-09, 1.8191E-09, 1.5897E-09,  // F14650
  1.4191E-09, 1.2448E-09, 1.0884E-09, 9.3585E-10, 7.9429E-10,  // F14660
  7.3214E-10, 6.5008E-10, 5.7549E-10, 5.4300E-10, 4.7251E-10,  // F14670
  4.3451E-10, 3.8446E-10, 3.5589E-10, 3.4432E-10, 2.8209E-10,  // F14680
  2.4620E-10, 2.1278E-10, 1.8406E-10, 1.6314E-10, 1.3261E-10,  // F14690
  1.1696E-10, 9.6865E-11, 7.6814E-11, 6.6411E-11, 5.0903E-11,  // F14700
  4.0827E-11, 3.0476E-11, 2.3230E-11, 1.7707E-11, 1.3548E-11,  // F14710
  1.0719E-11, 9.3026E-12, 8.7967E-12, 8.3136E-12, 7.3918E-12,  // F14720
  6.5293E-12, 5.9243E-12, 5.3595E-12, 3.5266E-12, 2.2571E-12,  // F14730
  1.6150E-12, 1.1413E-12, 8.4998E-13, 7.0803E-13, 5.1747E-13,  // F14750
  4.0694E-13, 3.6528E-13, 3.3670E-13, 3.1341E-13, 2.9390E-13,  // F14760
  2.8680E-13, 3.1283E-13, 3.7294E-13, 5.0194E-13, 6.7919E-13,  // F14770
  1.0455E-12, 1.5230E-12, 2.3932E-12, 3.4231E-12, 5.0515E-12,  // F14780
  7.3193E-12, 9.9406E-12, 1.2193E-11, 1.4742E-11, 1.9269E-11,  // F14790
  2.1816E-11, 2.2750E-11, 2.2902E-11, 2.3888E-11, 2.4902E-11,  // F14800
  2.2160E-11, 2.0381E-11, 1.9903E-11, 2.0086E-11, 1.9304E-11,  // F14810
  2.0023E-11, 2.2244E-11, 2.5450E-11, 3.1228E-11, 3.4560E-11,  // F14820
  3.6923E-11, 3.7486E-11, 3.8124E-11, 3.8317E-11, 3.4737E-11,  // F14830
  3.3037E-11, 3.1724E-11, 2.9840E-11, 2.8301E-11, 2.5857E-11,  // F14840
  2.3708E-11, 1.9452E-11, 1.6232E-11, 1.5174E-11, 1.4206E-11,  // F14860
  1.4408E-11, 1.5483E-11, 1.8642E-11, 2.3664E-11, 3.0181E-11,  // F14870
  4.0160E-11, 5.2287E-11, 7.2754E-11, 1.0511E-10, 1.4531E-10,  // F14880
  2.0998E-10, 2.6883E-10, 3.3082E-10, 4.2638E-10, 5.3132E-10,  // F14890
  6.3617E-10, 7.1413E-10, 8.5953E-10, 9.9715E-10, 1.0796E-09,  // F14900
  1.0978E-09, 1.1052E-09, 1.1095E-09, 1.0641E-09, 9.7881E-10,  // F14910
  9.6590E-10, 1.0332E-09, 1.1974E-09, 1.3612E-09, 1.5829E-09,  // F14920
  1.8655E-09, 2.1465E-09, 2.4779E-09, 2.7370E-09, 2.9915E-09,  // F14930
  3.3037E-09, 3.6347E-09, 3.9587E-09, 4.4701E-09, 5.0122E-09,  // F14940
  5.8044E-09, 6.1916E-09, 6.9613E-09, 7.7863E-09, 8.2820E-09,  // F14950
  9.4359E-09, 9.7387E-09, 1.0656E-08, 1.0746E-08, 1.1210E-08,  // F14970
  1.1905E-08, 1.2194E-08, 1.3145E-08, 1.3738E-08, 1.3634E-08,  // F14980
  1.3011E-08, 1.2511E-08, 1.1805E-08, 1.2159E-08, 1.2390E-08,  // F14990
  1.3625E-08, 1.5678E-08, 1.7886E-08, 1.9933E-08, 1.9865E-08,  // F15000
  1.9000E-08, 1.7812E-08, 1.5521E-08, 1.2593E-08, 9.5635E-09,  // F15010
  7.2987E-09, 5.2489E-09, 3.5673E-09, 2.4206E-09, 1.6977E-09,  // F15020
  1.2456E-09, 9.3744E-10, 7.8379E-10, 6.9960E-10, 6.6451E-10,  // F15030
  6.8521E-10, 7.4234E-10, 8.6658E-10, 9.4972E-10, 1.0791E-09,  // F15040
  1.2359E-09, 1.3363E-09, 1.5025E-09, 1.5368E-09, 1.6152E-09,  // F15050
  1.6184E-09, 1.6557E-09, 1.7035E-09, 1.6916E-09, 1.7237E-09,  // F15060
  1.7175E-09, 1.6475E-09, 1.5335E-09, 1.4272E-09, 1.3282E-09,  // F15080
  1.3459E-09, 1.4028E-09, 1.5192E-09, 1.7068E-09, 1.9085E-09,  // F15090
  2.1318E-09, 2.1020E-09, 1.9942E-09, 1.8654E-09, 1.6391E-09,  // F15100
  1.3552E-09, 1.0186E-09, 7.8540E-10, 5.7022E-10, 3.9247E-10,  // F15110
  2.5441E-10, 1.6699E-10, 1.1132E-10, 6.8989E-11, 4.5255E-11,  // F15120
  3.1106E-11, 2.3161E-11, 1.7618E-11, 1.4380E-11, 1.1601E-11,  // F15130
  9.7148E-12, 8.4519E-12, 6.5392E-12, 5.4113E-12, 4.7624E-12,  // F15140
  4.0617E-12, 3.6173E-12, 2.8608E-12, 2.2724E-12, 1.7436E-12,  // F15150
  1.3424E-12, 1.0358E-12, 7.3064E-13, 5.4500E-13, 4.0551E-13,  // F15160
  2.8642E-13, 2.1831E-13, 1.6860E-13, 1.2086E-13, 1.0150E-13,  // F15170
  9.3550E-14, 8.4105E-14, 7.3051E-14, 6.9796E-14, 7.9949E-14,  // F15190
  1.0742E-13, 1.5639E-13, 2.1308E-13, 3.1226E-13, 4.6853E-13,  // F15200
  6.6917E-13, 1.0088E-12, 1.4824E-12, 2.2763E-12, 3.3917E-12,  // F15210
  4.4585E-12, 6.3187E-12, 8.4189E-12, 1.1302E-11, 1.3431E-11,  // F15220
  1.5679E-11, 1.9044E-11, 2.2463E-11, 2.3605E-11, 2.3619E-11,  // F15230
  2.3505E-11, 2.3805E-11, 2.2549E-11, 1.9304E-11, 1.8382E-11,  // F15240
  1.7795E-11, 1.8439E-11, 1.9146E-11, 2.1966E-11, 2.6109E-11,  // F15250
  3.1883E-11, 3.7872E-11, 4.3966E-11, 4.8789E-11, 5.3264E-11,  // F15260
  5.9705E-11, 6.3744E-11, 7.0163E-11, 7.9114E-11, 8.8287E-11,  // F15270
  9.9726E-11, 1.1498E-10, 1.3700E-10, 1.6145E-10, 1.9913E-10,  // F15280
  2.2778E-10, 2.6216E-10, 2.9770E-10, 3.3405E-10, 3.7821E-10,  // F15300
  3.9552E-10, 4.1322E-10, 4.0293E-10, 4.0259E-10, 3.8853E-10,  // F15310
  3.7842E-10, 3.8551E-10, 4.4618E-10, 5.0527E-10, 5.0695E-10,  // F15320
  5.1216E-10, 5.1930E-10, 5.5794E-10, 5.3320E-10, 5.2008E-10,  // F15330
  5.6888E-10, 6.1883E-10, 6.9006E-10, 6.9505E-10, 6.6768E-10,  // F15340
  6.3290E-10, 5.6753E-10, 5.0327E-10, 3.9830E-10, 3.1147E-10,  // F15350
  2.4416E-10, 1.8860E-10, 1.3908E-10, 9.9156E-11, 7.3779E-11,  // F15360
  5.6048E-11, 4.2457E-11, 3.4505E-11, 2.9881E-11, 2.7865E-11,  // F15370
  2.8471E-11, 3.1065E-11, 3.4204E-11, 3.9140E-11, 4.3606E-11,  // F15380
  4.9075E-11, 5.3069E-11, 5.5236E-11, 5.5309E-11, 5.3832E-11,  // F15390
  5.3183E-11, 5.1783E-11, 5.2042E-11, 5.4422E-11, 5.5656E-11,  // F15410
  5.4409E-11, 5.2659E-11, 5.1696E-11, 5.1726E-11, 4.9003E-11,  // F15420
  4.9050E-11, 5.1700E-11, 5.6818E-11, 6.3129E-11, 6.6542E-11,  // F15430
  6.4367E-11, 5.9908E-11, 5.4470E-11, 4.7903E-11, 3.9669E-11,  // F15440
  2.9651E-11, 2.2286E-11, 1.6742E-11, 1.1827E-11, 7.7739E-12,  // F15450
  4.8805E-12, 3.1747E-12, 2.0057E-12, 1.2550E-12, 8.7434E-13,  // F15460
  6.2755E-13, 4.9752E-13, 4.0047E-13, 3.5602E-13, 3.0930E-13,  // F15470
  2.4903E-13, 1.9316E-13, 1.4995E-13, 1.2059E-13, 8.7242E-14,  // F15480
  6.4511E-14, 5.3300E-14, 4.3741E-14, 3.4916E-14, 2.6560E-14,  // F15490
  1.6923E-14, 1.1816E-14, 6.7071E-15, 3.6474E-15, 2.0686E-15,  // F15500
  1.1925E-15, 6.8948E-16, 3.9661E-16, 2.2576E-16, 1.2669E-16,  // F15520
  6.9908E-17, 3.7896E-17, 2.0280E-17, 1.1016E-17, 6.7816E-18,  // F15530
  6.0958E-18, 8.9913E-18, 1.7201E-17, 3.4964E-17, 7.0722E-17,  // F15540
  1.4020E-16, 2.7167E-16, 5.1478E-16, 9.5500E-16, 1.7376E-15,  // F15550
  3.1074E-15, 5.4789E-15, 9.5640E-15, 1.6635E-14, 2.9145E-14,  // F15560
  5.2179E-14, 8.8554E-14, 1.4764E-13, 2.3331E-13, 3.5996E-13,  // F15570
  5.2132E-13, 6.3519E-13, 7.3174E-13, 8.3752E-13, 9.8916E-13,  // F15580
  1.1515E-12, 1.4034E-12, 1.6594E-12, 2.1021E-12, 2.7416E-12,  // F15590
  3.4135E-12, 4.5517E-12, 5.5832E-12, 7.2303E-12, 9.9484E-12,  // F15600
  1.2724E-11, 1.6478E-11, 2.0588E-11, 2.5543E-11, 3.3625E-11,  // F15610
  4.1788E-11, 5.0081E-11, 6.0144E-11, 6.9599E-11, 8.4408E-11,  // F15630
  9.7143E-11, 1.0805E-10, 1.1713E-10, 1.2711E-10, 1.3727E-10,  // F15640
  1.4539E-10, 1.6049E-10, 1.7680E-10, 2.0557E-10, 2.4967E-10,  // F15650
  3.0096E-10, 3.5816E-10, 4.0851E-10, 4.6111E-10, 5.2197E-10,  // F15660
  5.5043E-10, 6.0324E-10, 6.4983E-10, 6.7498E-10, 7.0545E-10,  // F15670
  7.0680E-10, 7.5218E-10, 7.5723E-10, 7.7840E-10, 8.0081E-10,  // F15680
  8.0223E-10, 7.7271E-10, 7.1676E-10, 6.7819E-10, 6.4753E-10,  // F15690
  6.5844E-10, 7.0163E-10, 7.7503E-10, 8.8152E-10, 9.9022E-10,  // F15700
  1.0229E-09, 9.9296E-10, 8.9911E-10, 7.7813E-10, 6.3785E-10,  // F15710
  4.7491E-10, 3.5280E-10, 2.4349E-10, 1.6502E-10, 1.1622E-10,  // F15720
  8.6715E-11, 6.7360E-11, 5.3910E-11, 4.5554E-11, 4.1300E-11,  // F15740
  3.9728E-11, 3.9000E-11, 3.9803E-11, 4.1514E-11, 4.3374E-11,  // F15750
  4.6831E-11, 4.8921E-11, 5.1995E-11, 5.7242E-11, 6.2759E-11,  // F15760
  7.0801E-11, 7.4555E-11, 7.9754E-11, 8.7616E-11, 9.1171E-11,  // F15770
  1.0349E-10, 1.1047E-10, 1.2024E-10, 1.2990E-10, 1.3725E-10,  // F15780
  1.5005E-10, 1.5268E-10, 1.5535E-10, 1.5623E-10, 1.5009E-10,  // F15790
  1.4034E-10, 1.3002E-10, 1.2225E-10, 1.1989E-10, 1.2411E-10,  // F15800
  1.3612E-10, 1.5225E-10, 1.7202E-10, 1.9471E-10, 1.9931E-10,  // F15810
  1.9079E-10, 1.7478E-10, 1.5259E-10, 1.2625E-10, 9.3332E-11,  // F15820
  6.8796E-11, 4.6466E-11, 2.9723E-11, 1.8508E-11, 1.2106E-11,  // F15830
  8.0142E-12, 5.4066E-12, 3.9329E-12, 3.1665E-12, 2.7420E-12,  // F15850
  2.3996E-12, 2.3804E-12, 2.3242E-12, 2.4476E-12, 2.5331E-12,  // F15860
  2.3595E-12, 2.2575E-12, 2.1298E-12, 2.0088E-12, 1.8263E-12,  // F15870
  1.6114E-12, 1.4422E-12, 1.2946E-12, 1.0837E-12, 9.1282E-13,  // F15880
  7.2359E-13, 5.3307E-13, 3.8837E-13, 2.6678E-13, 1.6769E-13,  // F15890
  1.0826E-13, 7.2364E-14, 4.5201E-14, 3.0808E-14, 2.2377E-14,  // F15900
  1.7040E-14, 9.2181E-15, 5.2934E-15, 3.5774E-15, 3.1431E-15,  // F15910
  3.7647E-15, 5.6428E-15, 9.5139E-15, 1.7322E-14, 2.8829E-14,  // F15920
  4.7708E-14, 6.9789E-14, 9.7267E-14, 1.4662E-13, 1.9429E-13,  // F15930
  2.5998E-13, 3.6636E-13, 4.7960E-13, 6.5129E-13, 7.7638E-13,  // F15940
  9.3774E-13, 1.1467E-12, 1.3547E-12, 1.5686E-12, 1.6893E-12,  // F15960
  1.9069E-12, 2.1352E-12, 2.3071E-12, 2.4759E-12, 2.8247E-12,  // F15970
  3.4365E-12, 4.3181E-12, 5.6107E-12, 7.0017E-12, 8.6408E-12,  // F15980
  1.0974E-11, 1.3742E-11, 1.6337E-11, 2.0157E-11, 2.3441E-11,  // F15990
  2.6733E-11, 3.0247E-11, 3.3737E-11, 3.8618E-11, 4.1343E-11,  // F16000
  4.3870E-11, 4.4685E-11, 4.4881E-11, 4.5526E-11, 4.3628E-11,  // F16010
  4.4268E-11, 4.6865E-11, 5.3426E-11, 5.4020E-11, 5.3218E-11,  // F16020
  5.4587E-11, 5.6360E-11, 5.7740E-11, 5.6426E-11, 6.0399E-11,  // F16030
  6.6981E-11, 7.4319E-11, 7.7977E-11, 7.5539E-11, 7.1610E-11,  // F16040
  6.4606E-11, 5.5498E-11, 4.3944E-11, 3.3769E-11, 2.5771E-11,  // F16050
  1.9162E-11, 1.3698E-11, 1.0173E-11, 7.8925E-12, 6.1938E-12,  // F16070
  4.7962E-12, 4.0811E-12, 3.3912E-12, 2.8625E-12, 2.4504E-12,  // F16080
  2.2188E-12, 2.2139E-12, 2.2499E-12, 2.2766E-12, 2.3985E-12,  // F16090
  2.5459E-12, 2.9295E-12, 3.4196E-12, 3.6155E-12, 4.0733E-12,  // F16100
  4.4610E-12, 4.9372E-12, 5.4372E-12, 5.7304E-12, 6.1640E-12,  // F16110
  6.1278E-12, 6.2940E-12, 6.4947E-12, 6.8174E-12, 7.5190E-12,  // F16120
  8.2608E-12, 8.4971E-12, 8.3484E-12, 8.1888E-12, 7.8552E-12,  // F16130
  7.8468E-12, 7.5943E-12, 7.9096E-12, 8.6869E-12, 9.1303E-12,  // F16140
  9.2547E-12, 8.9322E-12, 8.2177E-12, 7.3408E-12, 5.7956E-12,  // F16150
  4.4470E-12, 3.5881E-12, 2.6748E-12, 1.7074E-12, 9.6700E-13,  // F16160
  5.2645E-13, 2.9943E-13, 1.7316E-13, 1.0039E-13, 5.7859E-14,  // F16180
  3.2968E-14, 1.8499E-14, 1.0192E-14, 5.5015E-15, 2.9040E-15,  // F16190
  1.4968E-15, 7.5244E-16, 3.6852E-16, 1.7568E-16, 8.1464E-17,  // F16200
  3.6717E-17, 1.6076E-17, 6.8341E-18, 2.8195E-18, 1.1286E-18,  // F16210
   .0000E+00,  .0000E+00,  .0000E+00,  .0000E+00,  .0000E+00,  // F16220
   .0000E+00,  .0000E+00,  .0000E+00,  .0000E+00, 1.4070E-18,  // F16230
  3.0405E-18, 6.4059E-18, 1.3169E-17, 2.6443E-17, 5.1917E-17,  // F16240
  9.9785E-17, 1.8802E-16, 3.4788E-16, 6.3328E-16, 1.1370E-15,  // F16250
  2.0198E-15, 3.5665E-15, 6.3053E-15, 1.1309E-14, 2.1206E-14,  // F16260
  3.2858E-14, 5.5165E-14, 8.6231E-14, 1.2776E-13, 1.7780E-13,  // F16270
  2.5266E-13, 3.6254E-13, 5.1398E-13, 6.8289E-13, 8.7481E-13,  // F16290
  1.1914E-12, 1.6086E-12, 2.0469E-12, 2.5761E-12, 3.4964E-12,  // F16300
  4.4980E-12, 5.5356E-12, 6.7963E-12, 8.5720E-12, 1.0700E-11,  // F16310
  1.2983E-11, 1.6270E-11, 1.9609E-11, 2.2668E-11, 2.5963E-11,  // F16320
  3.0918E-11, 3.4930E-11, 3.9330E-11, 4.4208E-11, 4.6431E-11,  // F16330
  5.1141E-11, 5.4108E-11, 5.8077E-11, 6.5050E-11, 7.2126E-11,  // F16340
  8.1064E-11, 8.1973E-11, 8.1694E-11, 8.3081E-11, 8.0240E-11,  // F16350
  7.9225E-11, 7.6256E-11, 7.8468E-11, 8.0041E-11, 8.1585E-11,  // F16360
  8.3485E-11, 8.3774E-11, 8.5870E-11, 8.6104E-11, 8.8516E-11,  // F16370
  9.0814E-11, 9.2522E-11, 8.8913E-11, 7.8381E-11, 6.8568E-11,  // F16380
  5.6797E-11, 4.4163E-11, 3.2369E-11, 2.3259E-11, 1.6835E-11,  // F16400
  1.1733E-11, 8.5273E-12, 6.3805E-12, 4.8983E-12, 3.8831E-12,  // F16410
  3.2610E-12, 2.8577E-12, 2.5210E-12, 2.2913E-12, 2.0341E-12,  // F16420
  1.8167E-12, 1.6395E-12, 1.4890E-12, 1.3516E-12, 1.2542E-12,  // F16430
  1.2910E-12, 1.3471E-12, 1.4689E-12, 1.5889E-12, 1.6989E-12,  // F16440
  1.8843E-12, 2.0902E-12, 2.3874E-12, 2.7294E-12, 3.3353E-12,  // F16450
  4.0186E-12, 4.5868E-12, 5.2212E-12, 5.8856E-12, 6.5991E-12,  // F16460
  7.2505E-12, 7.6637E-12, 8.5113E-12, 9.4832E-12, 9.9678E-12,  // F16470
  1.0723E-11, 1.0749E-11, 1.1380E-11, 1.1774E-11, 1.1743E-11,  // F16480
  1.2493E-11, 1.2559E-11, 1.2332E-11, 1.1782E-11, 1.1086E-11,  // F16490
  1.0945E-11, 1.1178E-11, 1.2083E-11, 1.3037E-11, 1.4730E-11,  // F16510
  1.6450E-11, 1.7403E-11, 1.7004E-11, 1.5117E-11, 1.3339E-11,  // F16520
  1.0844E-11, 8.0915E-12, 5.6615E-12, 3.7196E-12, 2.5194E-12,  // F16530
  1.6569E-12, 1.1201E-12, 8.2335E-13, 6.0270E-13, 4.8205E-13,  // F16540
  4.1313E-13, 3.6243E-13, 3.2575E-13, 2.7730E-13, 2.5292E-13,  // F16550
  2.3062E-13, 2.1126E-13, 2.1556E-13, 2.1213E-13, 2.2103E-13,  // F16560
  2.1927E-13, 2.0794E-13, 1.9533E-13, 1.6592E-13, 1.4521E-13,  // F16570
  1.1393E-13, 8.3772E-14, 6.2077E-14, 4.3337E-14, 2.7165E-14,  // F16580
  1.6821E-14, 9.5407E-15, 5.3093E-15, 3.0320E-15, 1.7429E-15,  // F16590
  9.9828E-16, 5.6622E-16, 3.1672E-16, 1.7419E-16, 9.3985E-17,  // F16600
  4.9656E-17, 2.5652E-17, 1.2942E-17, 6.3695E-18, 3.0554E-18,  // F16620
  1.4273E-18, -0.       , -0.       , -0.       , -0.       ,  // F16630
  -0.       , 0.        , 0.        , 0.        , 0.        ,  // F16640
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16650
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16660
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16670
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16680
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16690
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16700
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16710
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16730
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16740
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16750
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16760
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16770
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16780
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16790
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16800
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16810
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16820
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16840
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16850
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16860
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16870
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16880
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16890
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16900
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16910
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16920
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16930
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16950
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16960
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16970
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16980
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F16990
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F17000
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F17010
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F17020
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F17030
  0.        , 0.        , 0.        , 0.        , 0.        ,  // F17040
  0.};






// H2O self continuum parameters at T=260 K for the CKD versions 0 to 2.4.2
// date of last update: 06/28/82
// units of (CM**3/MOL)*1.E-20
const Numeric SL260_ckd_0_v1  =   -20.0;
const Numeric SL260_ckd_0_v2  = 20000.0;
const Numeric SL260_ckd_0_dv  =    10.0;
const int     SL260_ckd_0_npt =  2003;
const double  SL260_ckd_0[SL260_ckd_0_npt+addF77fields] = {
  0.0000e-00, 1.7750e-01, 1.7045e-01, 
  1.6457e-01, 1.7045e-01, 1.7750e-01, 2.0036e-01, 2.1347e-01,  // f07700
  2.2454e-01, 2.3428e-01, 2.3399e-01, 2.3022e-01, 2.0724e-01,  // f07710
  1.9712e-01, 1.8317e-01, 1.6724e-01, 1.4780e-01, 1.2757e-01,  // f07720
  1.1626e-01, 1.0098e-01, 8.9033e-02, 7.9770e-02, 6.7416e-02,  // f07730
  5.9588e-02, 5.1117e-02, 4.6218e-02, 4.2179e-02, 3.4372e-02,  // f07740
  2.9863e-02, 2.5252e-02, 2.2075e-02, 1.9209e-02, 1.5816e-02,  // f07750
  1.3932e-02, 1.1943e-02, 1.0079e-02, 8.7667e-03, 7.4094e-03,  // f07760
  6.4967e-03, 5.5711e-03, 4.8444e-03, 4.2552e-03, 3.6953e-03,  // f07770
  3.2824e-03, 2.9124e-03, 2.6102e-03, 2.3370e-03, 2.1100e-03,  // f07780
  1.9008e-03, 1.7145e-03, 1.5573e-03, 1.4206e-03, 1.2931e-03,  // f07790
  1.1803e-03, 1.0774e-03, 9.8616e-04, 9.0496e-04, 8.3071e-04,  // f07810
  7.6319e-04, 7.0149e-04, 6.4637e-04, 5.9566e-04, 5.4987e-04,  // f07820
  5.0768e-04, 4.6880e-04, 4.3317e-04, 4.0037e-04, 3.7064e-04,  // f07830
  3.4325e-04, 3.1809e-04, 2.9501e-04, 2.7382e-04, 2.5430e-04,  // f07840
  2.3630e-04, 2.1977e-04, 2.0452e-04, 1.9042e-04, 1.7740e-04,  // f07850
  1.6544e-04, 1.5442e-04, 1.4425e-04, 1.3486e-04, 1.2618e-04,  // f07860
  1.1817e-04, 1.1076e-04, 1.0391e-04, 9.7563e-05, 9.1696e-05,  // f07870
  8.6272e-05, 8.1253e-05, 7.6607e-05, 7.2302e-05, 6.8311e-05,  // f07880
  6.4613e-05, 6.1183e-05, 5.8001e-05, 5.5048e-05, 5.2307e-05,  // f07890
  4.9761e-05, 4.7395e-05, 4.5197e-05, 4.3155e-05, 4.1256e-05,  // f07900
  3.9491e-05, 3.7849e-05, 3.6324e-05, 3.4908e-05, 3.3594e-05,  // f07920
  3.2374e-05, 3.1244e-05, 3.0201e-05, 2.9240e-05, 2.8356e-05,  // f07930
  2.7547e-05, 2.6814e-05, 2.6147e-05, 2.5551e-05, 2.5029e-05,  // f07940
  2.4582e-05, 2.4203e-05, 2.3891e-05, 2.3663e-05, 2.3531e-05,  // f07950
  2.3483e-05, 2.3516e-05, 2.3694e-05, 2.4032e-05, 2.4579e-05,  // f07960
  2.5234e-05, 2.6032e-05, 2.7119e-05, 2.8631e-05, 3.0848e-05,  // f07970
  3.3262e-05, 3.6635e-05, 4.0732e-05, 4.5923e-05, 5.3373e-05,  // f07980
  6.1875e-05, 7.2031e-05, 8.5980e-05, 9.8642e-05, 1.1469e-04,  // f07990
  1.3327e-04, 1.5390e-04, 1.7513e-04, 2.0665e-04, 2.3609e-04,  // f08000
  2.6220e-04, 2.8677e-04, 3.2590e-04, 3.8624e-04, 4.1570e-04,  // f08010
  4.5207e-04, 4.9336e-04, 5.4500e-04, 5.8258e-04, 5.8086e-04,  // f08030
  5.6977e-04, 5.3085e-04, 4.8020e-04, 4.3915e-04, 4.0343e-04,  // f08040
  3.7853e-04, 3.7025e-04, 3.9637e-04, 4.4675e-04, 4.7072e-04,  // f08050
  4.9022e-04, 5.2076e-04, 5.3676e-04, 5.2755e-04, 4.8244e-04,  // f08060
  4.5473e-04, 4.3952e-04, 3.9614e-04, 3.4086e-04, 2.9733e-04,  // f08070
  2.6367e-04, 2.3767e-04, 2.0427e-04, 1.7595e-04, 1.5493e-04,  // f08080
  1.3851e-04, 1.1874e-04, 1.0735e-04, 9.0490e-05, 8.1149e-05,  // f08090
  7.4788e-05, 6.5438e-05, 5.8248e-05, 4.8076e-05, 4.3488e-05,  // f08100
  3.7856e-05, 3.3034e-05, 2.9592e-05, 2.6088e-05, 2.3497e-05,  // f08110
  2.0279e-05, 1.7526e-05, 1.5714e-05, 1.3553e-05, 1.2145e-05,  // f08120
  1.0802e-05, 9.7681e-06, 8.8196e-06, 7.8291e-06, 7.1335e-06,  // f08140
  6.4234e-06, 5.8391e-06, 5.3532e-06, 4.9079e-06, 4.5378e-06,  // f08150
  4.1716e-06, 3.8649e-06, 3.5893e-06, 3.3406e-06, 3.1199e-06,  // f08160
  2.9172e-06, 2.7348e-06, 2.5644e-06, 2.4086e-06, 2.2664e-06,  // f08170
  2.1359e-06, 2.0159e-06, 1.9051e-06, 1.8031e-06, 1.7074e-06,  // f08180
  1.6185e-06, 1.5356e-06, 1.4584e-06, 1.3861e-06, 1.3179e-06,  // f08190
  1.2545e-06, 1.1951e-06, 1.1395e-06, 1.0873e-06, 1.0384e-06,  // f08200
  9.9250e-07, 9.4935e-07, 9.0873e-07, 8.7050e-07, 8.3446e-07,  // f08210
  8.0046e-07, 7.6834e-07, 7.3800e-07, 7.0931e-07, 6.8217e-07,  // f08220
  6.5648e-07, 6.3214e-07, 6.0909e-07, 5.8725e-07, 5.6655e-07,  // f08230
  5.4693e-07, 5.2835e-07, 5.1077e-07, 4.9416e-07, 4.7853e-07,  // f08250
  4.6381e-07, 4.5007e-07, 4.3728e-07, 4.2550e-07, 4.1450e-07,  // f08260
  4.0459e-07, 3.9532e-07, 3.8662e-07, 3.7855e-07, 3.7041e-07,  // f08270
  3.6254e-07, 3.5420e-07, 3.4617e-07, 3.3838e-07, 3.3212e-07,  // f08280
  3.2655e-07, 3.1865e-07, 3.1203e-07, 3.0670e-07, 3.0252e-07,  // f08290
  2.9749e-07, 2.9184e-07, 2.8795e-07, 2.8501e-07, 2.8202e-07,  // f08300
  2.7856e-07, 2.7509e-07, 2.7152e-07, 2.6844e-07, 2.6642e-07,  // f08310
  2.6548e-07, 2.6617e-07, 2.6916e-07, 2.7372e-07, 2.8094e-07,  // f08320
  2.9236e-07, 3.1035e-07, 3.2854e-07, 3.5481e-07, 3.9377e-07,  // f08330
  4.4692e-07, 5.0761e-07, 5.7715e-07, 6.7725e-07, 8.0668e-07,  // f08340
  9.3716e-07, 1.0797e-06, 1.1689e-06, 1.3217e-06, 1.4814e-06,  // f08360
  1.5627e-06, 1.6519e-06, 1.7601e-06, 1.9060e-06, 2.0474e-06,  // f08370
  2.0716e-06, 2.0433e-06, 1.9752e-06, 1.8466e-06, 1.7526e-06,  // f08380
  1.6657e-06, 1.5870e-06, 1.5633e-06, 1.6520e-06, 1.8471e-06,  // f08390
  1.9953e-06, 2.0975e-06, 2.2016e-06, 2.2542e-06, 2.3081e-06,  // f08400
  2.3209e-06, 2.2998e-06, 2.3056e-06, 2.2757e-06, 2.2685e-06,  // f08410
  2.2779e-06, 2.2348e-06, 2.2445e-06, 2.3174e-06, 2.4284e-06,  // f08420
  2.5290e-06, 2.7340e-06, 2.9720e-06, 3.2332e-06, 3.5392e-06,  // f08430
  3.9013e-06, 4.3334e-06, 4.9088e-06, 5.3428e-06, 5.9142e-06,  // f08440
  6.6106e-06, 7.4709e-06, 8.5019e-06, 9.6835e-06, 1.0984e-05,  // f08450
  1.2831e-05, 1.4664e-05, 1.7080e-05, 2.0103e-05, 2.4148e-05,  // f08470
  2.7948e-05, 3.2855e-05, 3.9046e-05, 4.6429e-05, 5.6633e-05,  // f08480
  6.6305e-05, 7.6048e-05, 8.7398e-05, 1.0034e-04, 1.1169e-04,  // f08490
  1.2813e-04, 1.3354e-04, 1.3952e-04, 1.4204e-04, 1.4615e-04,  // f08500
  1.5144e-04, 1.5475e-04, 1.6561e-04, 1.7135e-04, 1.6831e-04,  // f08510
  1.6429e-04, 1.6353e-04, 1.6543e-04, 1.5944e-04, 1.5404e-04,  // f08520
  1.5458e-04, 1.6287e-04, 1.7277e-04, 1.8387e-04, 1.7622e-04,  // f08530
  1.6360e-04, 1.5273e-04, 1.3667e-04, 1.2364e-04, 9.7576e-05,  // f08540
  7.9140e-05, 6.4241e-05, 5.1826e-05, 4.1415e-05, 3.1347e-05,  // f08550
  2.5125e-05, 2.0027e-05, 1.6362e-05, 1.3364e-05, 1.1117e-05,  // f08560
  9.4992e-06, 8.1581e-06, 7.1512e-06, 6.2692e-06, 5.5285e-06,  // f08580
  4.9000e-06, 4.3447e-06, 3.8906e-06, 3.4679e-06, 3.1089e-06,  // f08590
  2.8115e-06, 2.5496e-06, 2.2982e-06, 2.0861e-06, 1.8763e-06,  // f08600
  1.7035e-06, 1.5548e-06, 1.4107e-06, 1.2839e-06, 1.1706e-06,  // f08610
  1.0709e-06, 9.8099e-07, 8.9901e-07, 8.2394e-07, 7.5567e-07,  // f08620
  6.9434e-07, 6.3867e-07, 5.8845e-07, 5.4263e-07, 5.0033e-07,  // f08630
  4.6181e-07, 4.2652e-07, 3.9437e-07, 3.6497e-07, 3.3781e-07,  // f08640
  3.1292e-07, 2.9011e-07, 2.6915e-07, 2.4989e-07, 2.3215e-07,  // f08650
  2.1582e-07, 2.0081e-07, 1.8700e-07, 1.7432e-07, 1.6264e-07,  // f08660
  1.5191e-07, 1.4207e-07, 1.3306e-07, 1.2484e-07, 1.1737e-07,  // f08670
  1.1056e-07, 1.0451e-07, 9.9060e-08, 9.4135e-08, 8.9608e-08,  // f08690
  8.5697e-08, 8.1945e-08, 7.8308e-08, 7.4808e-08, 7.1686e-08,  // f08700
  6.8923e-08, 6.5869e-08, 6.3308e-08, 6.0840e-08, 5.8676e-08,  // f08710
  5.6744e-08, 5.5016e-08, 5.3813e-08, 5.2792e-08, 5.2097e-08,  // f08720
  5.1737e-08, 5.1603e-08, 5.1656e-08, 5.1989e-08, 5.2467e-08,  // f08730
  5.2918e-08, 5.3589e-08, 5.4560e-08, 5.5869e-08, 5.7403e-08,  // f08740
  5.8968e-08, 6.0973e-08, 6.3432e-08, 6.6245e-08, 6.9353e-08,  // f08750
  7.2686e-08, 7.6541e-08, 8.0991e-08, 8.5950e-08, 9.1429e-08,  // f08760
  9.7851e-08, 1.0516e-07, 1.1349e-07, 1.2295e-07, 1.3335e-07,  // f08770
  1.4488e-07, 1.5864e-07, 1.7412e-07, 1.9140e-07, 2.1078e-07,  // f08780
  2.3369e-07, 2.5996e-07, 2.8848e-07, 3.2169e-07, 3.5991e-07,  // f08800
  4.0566e-07, 4.5969e-07, 5.3094e-07, 6.1458e-07, 7.1155e-07,  // f08810
  8.3045e-07, 9.9021e-07, 1.2042e-06, 1.4914e-06, 1.8145e-06,  // f08820
  2.2210e-06, 2.7831e-06, 3.4533e-06, 4.4446e-06, 5.1989e-06,  // f08830
  6.2289e-06, 7.1167e-06, 8.3949e-06, 9.6417e-06, 1.0313e-05,  // f08840
  1.0485e-05, 1.0641e-05, 1.0898e-05, 1.0763e-05, 1.0506e-05,  // f08850
  1.0497e-05, 1.1696e-05, 1.2654e-05, 1.3029e-05, 1.3175e-05,  // f08860
  1.4264e-05, 1.4985e-05, 1.4999e-05, 1.4317e-05, 1.4616e-05,  // f08870
  1.4963e-05, 1.5208e-05, 1.4942e-05, 1.3879e-05, 1.3087e-05,  // f08880
  1.1727e-05, 1.0515e-05, 9.0073e-06, 7.3133e-06, 6.1181e-06,  // f08890
  5.0623e-06, 4.1105e-06, 3.3915e-06, 2.6711e-06, 2.1464e-06,  // f08910
  1.7335e-06, 1.4302e-06, 1.1847e-06, 9.9434e-07, 8.2689e-07,  // f08920
  7.0589e-07, 6.0750e-07, 5.3176e-07, 4.6936e-07, 4.1541e-07,  // f08930
  3.6625e-07, 3.2509e-07, 2.9156e-07, 2.6308e-07, 2.3819e-07,  // f08940
  2.1421e-07, 1.9366e-07, 1.7626e-07, 1.5982e-07, 1.4567e-07,  // f08950
  1.3354e-07, 1.2097e-07, 1.1029e-07, 1.0063e-07, 9.2003e-08,  // f08960
  8.4245e-08, 7.7004e-08, 7.0636e-08, 6.4923e-08, 5.9503e-08,  // f08970
  5.4742e-08, 5.0450e-08, 4.6470e-08, 4.2881e-08, 3.9550e-08,  // f08980
  3.6541e-08, 3.3803e-08, 3.1279e-08, 2.8955e-08, 2.6858e-08,  // f08990
  2.4905e-08, 2.3146e-08, 2.1539e-08, 2.0079e-08, 1.8746e-08,  // f09000
  1.7517e-08, 1.6396e-08, 1.5369e-08, 1.4426e-08, 1.3543e-08,  // f09020
  1.2724e-08, 1.1965e-08, 1.1267e-08, 1.0617e-08, 1.0010e-08,  // f09030
  9.4662e-09, 8.9553e-09, 8.4988e-09, 8.0807e-09, 7.7043e-09,  // f09040
  7.3721e-09, 7.0707e-09, 6.8047e-09, 6.5702e-09, 6.3634e-09,  // f09050
  6.1817e-09, 6.0239e-09, 5.8922e-09, 5.7824e-09, 5.7019e-09,  // f09060
  5.6368e-09, 5.5940e-09, 5.5669e-09, 5.5583e-09, 5.5653e-09,  // f09070
  5.5837e-09, 5.6243e-09, 5.6883e-09, 5.7800e-09, 5.8964e-09,  // f09080
  6.0429e-09, 6.2211e-09, 6.4282e-09, 6.6634e-09, 6.9306e-09,  // f09090
  7.2336e-09, 7.5739e-09, 7.9562e-09, 8.3779e-09, 8.8575e-09,  // f09100
  9.3992e-09, 1.0004e-08, 1.0684e-08, 1.1450e-08, 1.2320e-08,  // f09110
  1.3311e-08, 1.4455e-08, 1.5758e-08, 1.7254e-08, 1.8927e-08,  // f09130
  2.0930e-08, 2.3348e-08, 2.6074e-08, 2.9221e-08, 3.2770e-08,  // f09140
  3.7485e-08, 4.2569e-08, 4.8981e-08, 5.5606e-08, 6.2393e-08,  // f09150
  7.1901e-08, 8.2921e-08, 9.5513e-08, 1.1111e-07, 1.3143e-07,  // f09160
  1.5971e-07, 1.8927e-07, 2.2643e-07, 2.7860e-07, 3.2591e-07,  // f09170
  3.7024e-07, 4.2059e-07, 4.9432e-07, 5.5543e-07, 5.7498e-07,  // f09180
  5.9210e-07, 6.1005e-07, 6.1577e-07, 5.9193e-07, 5.6602e-07,  // f09190
  5.7403e-07, 6.0050e-07, 6.4723e-07, 6.7073e-07, 7.5415e-07,  // f09200
  8.0982e-07, 8.7658e-07, 9.1430e-07, 9.4459e-07, 9.8347e-07,  // f09210
  9.8768e-07, 1.0153e-06, 1.0066e-06, 1.0353e-06, 1.0353e-06,  // f09220
  1.0722e-06, 1.1138e-06, 1.1923e-06, 1.2947e-06, 1.4431e-06,  // f09240
  1.6537e-06, 1.8662e-06, 2.2473e-06, 2.6464e-06, 3.1041e-06,  // f09250
  3.4858e-06, 4.0167e-06, 4.6675e-06, 5.0983e-06, 5.7997e-06,  // f09260
  6.0503e-06, 6.4687e-06, 6.5396e-06, 6.7986e-06, 7.0244e-06,  // f09270
  7.2305e-06, 7.6732e-06, 7.9783e-06, 7.9846e-06, 7.7617e-06,  // f09280
  7.7657e-06, 7.7411e-06, 7.8816e-06, 7.8136e-06, 8.0051e-06,  // f09290
  8.5799e-06, 9.1659e-06, 9.8646e-06, 9.4920e-06, 8.7670e-06,  // f09300
  8.2034e-06, 7.2297e-06, 6.2324e-06, 4.9315e-06, 3.9128e-06,  // f09310
  3.1517e-06, 2.4469e-06, 1.8815e-06, 1.4627e-06, 1.1698e-06,  // f09320
  9.4686e-07, 7.8486e-07, 6.6970e-07, 5.8811e-07, 5.2198e-07,  // f09330
  4.6809e-07, 4.1671e-07, 3.7006e-07, 3.3066e-07, 2.9387e-07,  // f09350
  2.6415e-07, 2.3409e-07, 2.0991e-07, 1.9132e-07, 1.7519e-07,  // f09360
  1.5939e-07, 1.4368e-07, 1.3050e-07, 1.1883e-07, 1.0772e-07,  // f09370
  9.6884e-08, 8.7888e-08, 7.8956e-08, 7.1024e-08, 6.3824e-08,  // f09380
  5.7256e-08, 5.1769e-08, 4.7037e-08, 4.2901e-08, 3.8970e-08,  // f09390
  3.5467e-08, 3.2502e-08, 2.9827e-08, 2.7389e-08, 2.5111e-08,  // f09400
  2.3056e-08, 2.1267e-08, 1.9610e-08, 1.8133e-08, 1.6775e-08,  // f09410
  1.5491e-08, 1.4329e-08, 1.3265e-08, 1.2300e-08, 1.1420e-08,  // f09420
  1.0593e-08, 9.8475e-09, 9.1585e-09, 8.5256e-09, 7.9525e-09,  // f09430
  7.4226e-09, 6.9379e-09, 6.4950e-09, 6.0911e-09, 5.7242e-09,  // f09440
  5.3877e-09, 5.0821e-09, 4.8051e-09, 4.5554e-09, 4.3315e-09,  // f09460
  4.1336e-09, 3.9632e-09, 3.8185e-09, 3.7080e-09, 3.6296e-09,  // f09470
  3.5804e-09, 3.5776e-09, 3.6253e-09, 3.7115e-09, 3.8151e-09,  // f09480
  3.9804e-09, 4.1742e-09, 4.3581e-09, 4.5306e-09, 4.7736e-09,  // f09490
  5.1297e-09, 5.5291e-09, 5.9125e-09, 6.4956e-09, 7.0362e-09,  // f09500
  7.5318e-09, 7.9947e-09, 8.6438e-09, 9.7227e-09, 1.0130e-08,  // f09510
  1.0549e-08, 1.1064e-08, 1.1702e-08, 1.2043e-08, 1.1781e-08,  // f09520
  1.1838e-08, 1.1917e-08, 1.2131e-08, 1.2476e-08, 1.3611e-08,  // f09530
  1.4360e-08, 1.5057e-08, 1.6247e-08, 1.7284e-08, 1.8420e-08,  // f09540
  1.8352e-08, 1.8722e-08, 1.9112e-08, 1.9092e-08, 1.9311e-08,  // f09550
  1.9411e-08, 1.9884e-08, 2.0508e-08, 2.1510e-08, 2.3143e-08,  // f09570
  2.5050e-08, 2.7596e-08, 3.1231e-08, 3.6260e-08, 4.3410e-08,  // f09580
  5.2240e-08, 6.3236e-08, 7.7522e-08, 9.8688e-08, 1.1859e-07,  // f09590
  1.4341e-07, 1.6798e-07, 1.9825e-07, 2.2898e-07, 2.6257e-07,  // f09600
  2.9884e-07, 3.3247e-07, 3.4936e-07, 3.5583e-07, 3.7150e-07,  // f09610
  3.6580e-07, 3.7124e-07, 3.7030e-07, 4.1536e-07, 4.6656e-07,  // f09620
  4.6677e-07, 4.7507e-07, 4.9653e-07, 5.3795e-07, 5.4957e-07,  // f09630
  5.2238e-07, 5.4690e-07, 5.6569e-07, 5.9844e-07, 5.9835e-07,  // f09640
  5.6522e-07, 5.4123e-07, 4.7904e-07, 4.2851e-07, 3.5603e-07,  // f09650
  2.8932e-07, 2.3655e-07, 1.8592e-07, 1.4943e-07, 1.1971e-07,  // f09660
  9.8482e-08, 8.3675e-08, 7.1270e-08, 6.2496e-08, 5.4999e-08,  // f09680
  4.9821e-08, 4.5387e-08, 4.1340e-08, 3.7453e-08, 3.3298e-08,  // f09690
  3.0120e-08, 2.7032e-08, 2.4236e-08, 2.1500e-08, 1.8988e-08,  // f09700
  1.7414e-08, 1.5706e-08, 1.4192e-08, 1.3204e-08, 1.1759e-08,  // f09710
  1.0737e-08, 9.6309e-09, 8.8179e-09, 8.2619e-09, 7.2264e-09,  // f09720
  6.4856e-09, 5.8037e-09, 5.2093e-09, 4.7205e-09, 4.1749e-09,  // f09730
  3.7852e-09, 3.3915e-09, 3.0089e-09, 2.7335e-09, 2.4398e-09,  // f09740
  2.2031e-09, 1.9786e-09, 1.7890e-09, 1.6266e-09, 1.4830e-09,  // f09750
  1.3576e-09, 1.2518e-09, 1.1587e-09, 1.0726e-09, 9.9106e-10,  // f09760
  9.1673e-10, 8.5084e-10, 7.9147e-10, 7.2882e-10, 6.7342e-10,  // f09770
  6.2593e-10, 5.8294e-10, 5.4435e-10, 5.0997e-10, 4.7806e-10,  // f09790
  4.4931e-10, 4.2357e-10, 4.0023e-10, 3.7909e-10, 3.5999e-10,  // f09800
  3.4285e-10, 3.2776e-10, 3.1468e-10, 3.0377e-10, 2.9479e-10,  // f09810
  2.8877e-10, 2.8512e-10, 2.8617e-10, 2.8976e-10, 3.0001e-10,  // f09820
  3.1718e-10, 3.3898e-10, 3.5857e-10, 3.8358e-10, 4.3131e-10,  // f09830
  4.5741e-10, 4.6948e-10, 4.7594e-10, 4.9529e-10, 5.1563e-10,  // f09840
  4.9475e-10, 4.8369e-10, 4.8829e-10, 5.0047e-10, 5.0203e-10,  // f09850
  5.1954e-10, 5.5352e-10, 5.9928e-10, 6.7148e-10, 7.1121e-10,  // f09860
  7.4317e-10, 7.6039e-10, 7.8313e-10, 8.0684e-10, 7.8553e-10,  // f09870
  7.8312e-10, 7.8537e-10, 7.8872e-10, 8.0185e-10, 8.1004e-10,  // f09880
  8.2608e-10, 8.2525e-10, 8.3857e-10, 8.7920e-10, 9.2451e-10,  // f09900
  9.8661e-10, 1.0629e-09, 1.1659e-09, 1.2922e-09, 1.4387e-09,  // f09910
  1.6254e-09, 1.8425e-09, 2.1428e-09, 2.5477e-09, 3.0379e-09,  // f09920
  3.7570e-09, 4.4354e-09, 5.1802e-09, 6.2769e-09, 7.4894e-09,  // f09930
  8.7474e-09, 9.8037e-09, 1.1582e-08, 1.3293e-08, 1.4471e-08,  // f09940
  1.5025e-08, 1.5580e-08, 1.6228e-08, 1.6413e-08, 1.6020e-08,  // f09950
  1.6393e-08, 1.7545e-08, 1.9590e-08, 2.1449e-08, 2.3856e-08,  // f09960
  2.7050e-08, 3.0214e-08, 3.3733e-08, 3.6487e-08, 3.9353e-08,  // f09970
  4.2660e-08, 4.6385e-08, 4.9955e-08, 5.5313e-08, 6.0923e-08,  // f09980
  6.8948e-08, 7.3649e-08, 8.2602e-08, 9.2212e-08, 9.9080e-08,  // f09990
  1.1319e-07, 1.1790e-07, 1.2941e-07, 1.3199e-07, 1.3914e-07,  // f10010
  1.4843e-07, 1.5300e-07, 1.6419e-07, 1.7095e-07, 1.6988e-07,  // f10020
  1.6494e-07, 1.6327e-07, 1.6067e-07, 1.6909e-07, 1.7118e-07,  // f10030
  1.8106e-07, 1.9857e-07, 2.1696e-07, 2.3385e-07, 2.2776e-07,  // f10040
  2.1402e-07, 1.9882e-07, 1.7362e-07, 1.4308e-07, 1.1158e-07,  // f10050
  8.8781e-08, 6.8689e-08, 5.2062e-08, 4.0427e-08, 3.2669e-08,  // f10060
  2.7354e-08, 2.3200e-08, 2.0580e-08, 1.8676e-08, 1.7329e-08,  // f10070
  1.6621e-08, 1.6433e-08, 1.6953e-08, 1.7134e-08, 1.7948e-08,  // f10080
  1.9107e-08, 1.9875e-08, 2.1416e-08, 2.1556e-08, 2.2265e-08,  // f10090
  2.2171e-08, 2.2534e-08, 2.3029e-08, 2.2828e-08, 2.3143e-08,  // f10100
  2.2965e-08, 2.2223e-08, 2.1108e-08, 2.0265e-08, 1.9516e-08,  // f10120
  1.9941e-08, 2.0312e-08, 2.1080e-08, 2.2611e-08, 2.4210e-08,  // f10130
  2.6069e-08, 2.5097e-08, 2.3318e-08, 2.1543e-08, 1.8942e-08,  // f10140
  1.5960e-08, 1.2386e-08, 9.9340e-09, 7.7502e-09, 5.9462e-09,  // f10150
  4.5113e-09, 3.5523e-09, 2.8844e-09, 2.3394e-09, 1.9584e-09,  // f10160
  1.6749e-09, 1.4624e-09, 1.2809e-09, 1.1359e-09, 1.0087e-09,  // f10170
  9.0166e-10, 8.1079e-10, 7.2219e-10, 6.4922e-10, 5.8803e-10,  // f10180
  5.3290e-10, 4.8590e-10, 4.4111e-10, 4.0184e-10, 3.6644e-10,  // f10190
  3.3529e-10, 3.0789e-10, 2.8286e-10, 2.6089e-10, 2.4125e-10,  // f10200
  2.2355e-10, 2.0783e-10, 1.9370e-10, 1.8088e-10, 1.6948e-10,  // f10210
  1.5929e-10, 1.5013e-10, 1.4193e-10, 1.3470e-10, 1.2841e-10,  // f10230
  1.2307e-10, 1.1865e-10, 1.1502e-10, 1.1243e-10, 1.1099e-10,  // f10240
  1.1066e-10, 1.1216e-10, 1.1529e-10, 1.2171e-10, 1.3128e-10,  // f10250
  1.4153e-10, 1.5962e-10, 1.8048e-10, 2.0936e-10, 2.3165e-10,  // f10260
  2.5746e-10, 2.9600e-10, 3.3707e-10, 3.5267e-10, 3.5953e-10,  // f10270
  3.6822e-10, 3.8363e-10, 3.8286e-10, 3.5883e-10, 3.6154e-10,  // f10280
  3.6653e-10, 3.8507e-10, 4.0250e-10, 4.4435e-10, 4.9889e-10,  // f10290
  5.6932e-10, 6.3599e-10, 7.0281e-10, 7.5777e-10, 8.1279e-10,  // f10300
  8.8910e-10, 9.3400e-10, 1.0076e-09, 1.0945e-09, 1.1898e-09,  // f10310
  1.3108e-09, 1.4725e-09, 1.7028e-09, 1.9619e-09, 2.3527e-09,  // f10320
  2.6488e-09, 3.0327e-09, 3.4396e-09, 3.8797e-09, 4.4115e-09,  // f10340
  4.6853e-09, 4.9553e-09, 4.9551e-09, 5.1062e-09, 5.0996e-09,  // f10350
  5.1119e-09, 5.2283e-09, 5.8297e-09, 6.3439e-09, 6.2675e-09,  // f10360
  6.3296e-09, 6.5173e-09, 7.1685e-09, 7.0528e-09, 6.8856e-09,  // f10370
  7.3182e-09, 7.6990e-09, 8.3461e-09, 8.1946e-09, 7.7153e-09,  // f10380
  7.2411e-09, 6.4511e-09, 5.7336e-09, 4.6105e-09, 3.6962e-09,  // f10390
  2.9944e-09, 2.4317e-09, 1.9399e-09, 1.5331e-09, 1.2633e-09,  // f10400
  1.0613e-09, 9.0136e-10, 7.9313e-10, 7.1543e-10, 6.6485e-10,  // f10410
  6.4225e-10, 6.3980e-10, 6.4598e-10, 6.7428e-10, 7.0270e-10,  // f10420
  7.4694e-10, 7.7946e-10, 7.9395e-10, 7.8716e-10, 7.6933e-10,  // f10430
  7.6220e-10, 7.4825e-10, 7.4805e-10, 7.6511e-10, 7.6492e-10,  // f10450
  7.4103e-10, 7.1979e-10, 7.1686e-10, 7.3403e-10, 7.1142e-10,  // f10460
  7.0212e-10, 7.1548e-10, 7.5253e-10, 8.0444e-10, 8.2378e-10,  // f10470
  7.8004e-10, 7.1712e-10, 6.4978e-10, 5.7573e-10, 4.8675e-10,  // f10480
  3.7945e-10, 3.0118e-10, 2.4241e-10, 1.9100e-10, 1.4816e-10,  // f10490
  1.1567e-10, 9.4183e-11, 7.7660e-11, 6.5270e-11, 5.6616e-11,  // f10500
  4.9576e-11, 4.4137e-11, 3.9459e-11, 3.5759e-11, 3.2478e-11,  // f10510
  2.9419e-11, 2.6703e-11, 2.4365e-11, 2.2412e-11, 2.0606e-11,  // f10520
  1.9067e-11, 1.7800e-11, 1.6695e-11, 1.5729e-11, 1.4887e-11,  // f10530
  1.4135e-11, 1.3519e-11, 1.2992e-11, 1.2563e-11, 1.2223e-11,  // f10540
  1.1962e-11, 1.1775e-11, 1.1657e-11, 1.1605e-11, 1.1619e-11,  // f10560
  1.1697e-11, 1.1839e-11, 1.2046e-11, 1.2319e-11, 1.2659e-11,  // f10570
  1.3070e-11, 1.3553e-11, 1.4113e-11, 1.4754e-11, 1.5480e-11,  // f10580
  1.6298e-11, 1.7214e-11, 1.8236e-11, 1.9372e-11, 2.0635e-11,  // f10590
  2.2036e-11, 2.3590e-11, 2.5317e-11, 2.7242e-11, 2.9400e-11,  // f10600
  3.1849e-11, 3.4654e-11, 3.7923e-11, 4.1695e-11, 4.6055e-11,  // f10610
  5.0940e-11, 5.5624e-11, 6.0667e-11, 6.6261e-11, 7.2692e-11,  // f10620
  7.9711e-11, 8.7976e-11, 9.6884e-11, 1.0775e-10, 1.2093e-10,  // f10630
  1.3531e-10, 1.5404e-10, 1.7315e-10, 1.9862e-10, 2.3341e-10,  // f10640
  2.7014e-10, 3.1716e-10, 3.6957e-10, 4.3233e-10, 5.2566e-10,  // f10650
  6.2251e-10, 7.2149e-10, 8.3958e-10, 9.5931e-10, 1.1388e-09,  // f10670
  1.2973e-09, 1.4442e-09, 1.5638e-09, 1.6974e-09, 1.8489e-09,  // f10680
  1.9830e-09, 2.1720e-09, 2.3662e-09, 2.6987e-09, 3.1697e-09,  // f10690
  3.6907e-09, 4.2625e-09, 4.7946e-09, 5.3848e-09, 6.0897e-09,  // f10700
  6.4730e-09, 7.1483e-09, 7.7432e-09, 8.0851e-09, 8.5013e-09,  // f10710
  8.5909e-09, 9.1890e-09, 9.3124e-09, 9.5936e-09, 9.8787e-09,  // f10720
  9.9036e-09, 9.6712e-09, 9.2036e-09, 9.0466e-09, 8.9380e-09,  // f10730
  9.1815e-09, 9.5092e-09, 1.0027e-08, 1.0876e-08, 1.1744e-08,  // f10740
  1.1853e-08, 1.1296e-08, 1.0134e-08, 8.8245e-09, 7.3930e-09,  // f10750
  5.7150e-09, 4.4884e-09, 3.4027e-09, 2.6054e-09, 2.0790e-09,  // f10760
  1.7267e-09, 1.4724e-09, 1.2722e-09, 1.1234e-09, 1.0186e-09,  // f10780
  9.4680e-10, 8.8854e-10, 8.5127e-10, 8.3157e-10, 8.2226e-10,  // f10790
  8.3395e-10, 8.3294e-10, 8.4725e-10, 8.8814e-10, 9.3697e-10,  // f10800
  1.0112e-09, 1.0412e-09, 1.0948e-09, 1.1810e-09, 1.2267e-09,  // f10810
  1.3690e-09, 1.4512e-09, 1.5568e-09, 1.6552e-09, 1.7321e-09,  // f10820
  1.8797e-09, 1.9210e-09, 1.9686e-09, 1.9917e-09, 1.9357e-09,  // f10830
  1.8486e-09, 1.7575e-09, 1.7113e-09, 1.7163e-09, 1.7623e-09,  // f10840
  1.8536e-09, 1.9765e-09, 2.1334e-09, 2.3237e-09, 2.3259e-09,  // f10850
  2.1833e-09, 1.9785e-09, 1.7308e-09, 1.4596e-09, 1.1198e-09,  // f10860
  8.7375e-10, 6.5381e-10, 4.8677e-10, 3.6756e-10, 2.9155e-10,  // f10870
  2.3735e-10, 1.9590e-10, 1.6638e-10, 1.4549e-10, 1.2947e-10,  // f10890
  1.1511e-10, 1.0548e-10, 9.6511e-11, 9.0469e-11, 8.5170e-11,  // f10900
  7.7804e-11, 7.1971e-11, 6.6213e-11, 6.1063e-11, 5.5881e-11,  // f10910
  5.0508e-11, 4.5932e-11, 4.1997e-11, 3.7672e-11, 3.3972e-11,  // f10920
  3.0318e-11, 2.6769e-11, 2.3874e-11, 2.1336e-11, 1.9073e-11,  // f10930
  1.7313e-11, 1.5904e-11, 1.4684e-11, 1.3698e-11, 1.2873e-11,  // f10940
  1.2175e-11, 1.1542e-11, 1.1024e-11, 1.0602e-11, 1.0267e-11,  // f10950
  1.0012e-11, 9.8379e-12, 9.7482e-12, 9.7564e-12, 9.8613e-12,  // f10960
  1.0092e-11, 1.0418e-11, 1.0868e-11, 1.1585e-11, 1.2351e-11,  // f10970
  1.3372e-11, 1.4841e-11, 1.6457e-11, 1.8681e-11, 2.0550e-11,  // f10980
  2.2912e-11, 2.5958e-11, 2.9137e-11, 3.2368e-11, 3.4848e-11,  // f11000
  3.8462e-11, 4.2190e-11, 4.5629e-11, 4.9022e-11, 5.4232e-11,  // f11010
  6.1900e-11, 7.1953e-11, 8.5368e-11, 9.9699e-11, 1.1734e-10,  // f11020
  1.4185e-10, 1.7017e-10, 1.9813e-10, 2.3859e-10, 2.7304e-10,  // f11030
  3.0971e-10, 3.5129e-10, 3.9405e-10, 4.5194e-10, 4.8932e-10,  // f11040
  5.2436e-10, 5.4098e-10, 5.5542e-10, 5.7794e-10, 5.6992e-10,  // f11050
  5.8790e-10, 6.1526e-10, 6.8034e-10, 6.7956e-10, 6.6864e-10,  // f11060
  6.9329e-10, 7.2971e-10, 7.6546e-10, 7.5078e-10, 7.8406e-10,  // f11070
  8.3896e-10, 9.0111e-10, 9.1994e-10, 8.7189e-10, 8.1426e-10,  // f11080
  7.3097e-10, 6.3357e-10, 5.1371e-10, 4.0936e-10, 3.2918e-10,  // f11090
  2.6255e-10, 2.0724e-10, 1.6879e-10, 1.4165e-10, 1.1989e-10,  // f11110
  1.0125e-10, 8.9629e-11, 7.8458e-11, 6.8826e-11, 6.0935e-11,  // f11120
  5.5208e-11, 5.2262e-11, 5.0260e-11, 4.8457e-11, 4.7888e-11,  // f11130
  4.8032e-11, 5.0838e-11, 5.4668e-11, 5.5790e-11, 6.0056e-11,  // f11140
  6.3811e-11, 6.8848e-11, 7.4590e-11, 7.8249e-11, 8.3371e-11,  // f11150
  8.3641e-11, 8.6591e-11, 8.9599e-11, 9.3487e-11, 1.0066e-10,  // f11160
  1.0765e-10, 1.0851e-10, 1.0619e-10, 1.0557e-10, 1.0460e-10,  // f11170
  1.0796e-10, 1.0523e-10, 1.0674e-10, 1.1261e-10, 1.1431e-10,  // f11180
  1.1408e-10, 1.0901e-10, 9.9105e-11, 8.8077e-11, 6.9928e-11,  // f11190
  5.4595e-11, 4.5401e-11, 3.6313e-11, 2.6986e-11, 1.9463e-11,  // f11200
  1.4577e-11, 1.1583e-11, 9.5492e-12, 8.0770e-12, 6.9642e-12,  // f11220
  6.0966e-12, 5.4046e-12, 4.8431e-12, 4.3815e-12, 3.9987e-12,  // f11230
  3.6790e-12, 3.4113e-12, 3.1868e-12, 2.9992e-12, 2.8434e-12,  // f11240
  2.7153e-12, 2.6120e-12, 2.5311e-12, 2.4705e-12, 2.4290e-12,  // f11250
  2.4053e-12, 2.3988e-12, 2.4087e-12, 2.4349e-12, 2.4771e-12,  // f11260
  2.5355e-12, 2.6103e-12, 2.7019e-12, 2.8110e-12, 2.9383e-12,  // f11270
  3.0848e-12, 3.2518e-12, 3.4405e-12, 3.6527e-12, 3.8902e-12,  // f11280
  4.1555e-12, 4.4510e-12, 4.7801e-12, 5.1462e-12, 5.5539e-12,  // f11290
  6.0086e-12, 6.5171e-12, 7.0884e-12, 7.7357e-12, 8.4831e-12,  // f11300
  9.3096e-12, 1.0282e-11, 1.1407e-11, 1.2690e-11, 1.4148e-11,  // f11310
  1.5888e-11, 1.7992e-11, 2.0523e-11, 2.3342e-11, 2.6578e-11,  // f11330
  3.0909e-11, 3.6228e-11, 4.2053e-11, 4.9059e-11, 5.9273e-11,  // f11340
  7.0166e-11, 8.2298e-11, 9.7071e-11, 1.1673e-10, 1.4010e-10,  // f11350
  1.6621e-10, 2.0127e-10, 2.3586e-10, 2.7050e-10, 3.0950e-10,  // f11360
  3.6584e-10, 4.1278e-10, 4.6591e-10, 5.2220e-10, 5.5246e-10,  // f11370
  6.1500e-10, 6.5878e-10, 7.1167e-10, 7.9372e-10, 8.6975e-10,  // f11380
  9.6459e-10, 9.7368e-10, 9.8142e-10, 1.0202e-09, 1.0200e-09,  // f11390
  1.0356e-09, 1.0092e-09, 1.0269e-09, 1.0366e-09, 1.0490e-09,  // f11400
  1.0717e-09, 1.0792e-09, 1.1016e-09, 1.0849e-09, 1.0929e-09,  // f11410
  1.0971e-09, 1.0969e-09, 1.0460e-09, 9.2026e-10, 8.1113e-10,  // f11420
  6.8635e-10, 5.5369e-10, 4.2908e-10, 3.3384e-10, 2.6480e-10,  // f11440
  2.0810e-10, 1.6915e-10, 1.4051e-10, 1.1867e-10, 1.0158e-10,  // f11450
  8.8990e-11, 7.9175e-11, 7.0440e-11, 6.3453e-11, 5.7009e-11,  // f11460
  5.1662e-11, 4.7219e-11, 4.3454e-11, 4.0229e-11, 3.7689e-11,  // f11470
  3.6567e-11, 3.5865e-11, 3.5955e-11, 3.5928e-11, 3.6298e-11,  // f11480
  3.7629e-11, 3.9300e-11, 4.1829e-11, 4.4806e-11, 5.0534e-11,  // f11490
  5.6672e-11, 6.2138e-11, 6.8678e-11, 7.6111e-11, 8.4591e-11,  // f11500
  9.2634e-11, 9.8085e-11, 1.0830e-10, 1.1949e-10, 1.2511e-10,  // f11510
  1.3394e-10, 1.3505e-10, 1.4342e-10, 1.4874e-10, 1.4920e-10,  // f11520
  1.5872e-10, 1.5972e-10, 1.5821e-10, 1.5425e-10, 1.4937e-10,  // f11530
  1.5089e-10, 1.5521e-10, 1.6325e-10, 1.6924e-10, 1.8265e-10,  // f11550
  1.9612e-10, 2.0176e-10, 1.9359e-10, 1.7085e-10, 1.5197e-10,  // f11560
  1.2646e-10, 9.8552e-11, 7.4530e-11, 5.5052e-11, 4.2315e-11,  // f11570
  3.2736e-11, 2.6171e-11, 2.1909e-11, 1.8286e-11, 1.5752e-11,  // f11580
  1.3859e-11, 1.2288e-11, 1.1002e-11, 9.7534e-12, 8.8412e-12,  // f11590
  8.0169e-12, 7.2855e-12, 6.8734e-12, 6.4121e-12, 6.1471e-12,  // f11600
  5.7780e-12, 5.3478e-12, 4.9652e-12, 4.4043e-12, 3.9862e-12,  // f11610
  3.4684e-12, 2.9681e-12, 2.5791e-12, 2.2339e-12, 1.9247e-12,  // f11620
  1.6849e-12, 1.4863e-12, 1.3291e-12, 1.2021e-12, 1.0947e-12,  // f11630
  1.0015e-12, 9.1935e-13, 8.4612e-13, 7.8036e-13, 7.2100e-13,  // f11640
  6.6718e-13, 6.1821e-13, 5.7353e-13, 5.3269e-13, 4.9526e-13,  // f11660
  4.6093e-13, 4.2937e-13, 4.0034e-13, 3.7361e-13, 3.4895e-13,  // f11670
  3.2621e-13, 3.0520e-13, 2.8578e-13, 2.6782e-13, 2.5120e-13,  // f11680
  2.3581e-13, 2.2154e-13, 2.0832e-13, 1.9605e-13, 1.8466e-13,  // f11690
  1.7408e-13, 1.6425e-13, 1.5511e-13, 1.4661e-13, 1.3869e-13,  // f11700
  1.3131e-13, 1.2444e-13, 1.1803e-13, 1.1205e-13, 1.0646e-13,  // f11710
  1.0124e-13, 9.6358e-14, 9.1789e-14, 8.7509e-14, 8.3498e-14,  // f11720
  7.9735e-14, 7.6202e-14, 7.2882e-14, 6.9760e-14, 6.6822e-14,  // f11730
  6.4053e-14, 6.1442e-14, 5.8978e-14, 5.6650e-14, 5.4448e-14,  // f11740
  5.2364e-14, 5.0389e-14, 4.8516e-14, 4.6738e-14, 4.5048e-14,  // f11750
  4.3441e-14, 4.1911e-14, 4.0453e-14, 3.9063e-14, 3.7735e-14,  // f11770
  3.6467e-14, 3.5254e-14, 3.4093e-14, 3.2980e-14, 3.1914e-14,  // f11780
  3.0891e-14, 2.9909e-14, 2.8965e-14, 2.8058e-14, 2.7185e-14,  // f11790
  2.6344e-14, 2.5535e-14, 2.4755e-14, 2.4002e-14, 2.3276e-14,  // f11800
  2.2576e-14, 2.1899e-14, 2.1245e-14, 2.0613e-14, 2.0002e-14,  // f11810
  1.9411e-14, 1.8839e-14, 1.8285e-14, 1.7749e-14, 1.7230e-14,  // f11820
  1.6727e-14, 1.6240e-14, 1.5768e-14, 1.5310e-14, 1.4867e-14,  // f11830
  1.4436e-14, 1.4019e-14, 1.3614e-14, 1.3221e-14, 1.2840e-14,  // f11840
  1.2471e-14, 1.2112e-14, 1.1764e-14, 1.1425e-14, 1.1097e-14,  // f11850
  1.0779e-14, 1.0469e-14, 1.0169e-14, 9.8775e-15, 9.5943e-15,  // f11860
  9.3193e-15, 9.0522e-15, 8.7928e-15, 8.5409e-15, 8.2962e-15,  // f11880
  8.0586e-15, 7.8278e-15, 7.6036e-15, 7.3858e-15, 7.1742e-15,  // f11890
  6.9687e-15, 6.7691e-15, 6.5752e-15, 6.3868e-15, 6.2038e-15,  // f11900
  6.0260e-15, 5.8533e-15, 5.6856e-15, 5.5226e-15, 5.3642e-15,  // f11910
  5.2104e-15, 5.0610e-15, 4.9158e-15, 4.7748e-15, 4.6378e-15,  // f11920
  4.5047e-15, 4.3753e-15, 4.2497e-15, 4.1277e-15, 4.0091e-15,  // f11930
  3.8939e-15, 3.7820e-15, 3.6733e-15, 3.5677e-15, 3.4651e-15,  // f11940
  3.3655e-15, 3.2686e-15, 3.1746e-15, 3.0832e-15, 2.9944e-15,  // f11950
  2.9082e-15, 2.8244e-15, 2.7431e-15, 2.6640e-15, 2.5872e-15,  // f11960
  2.5126e-15, 2.4401e-15, 2.3697e-15, 2.3014e-15, 2.2349e-15,  // f11970
  2.1704e-15, 2.1077e-15, 2.0468e-15, 1.9877e-15, 1.9302e-15,  // f11990
  1.8744e-15, 1.8202e-15, 1.7675e-15, 1.7164e-15, 1.6667e-15,  // f12000
  1.6184e-15, 1.5716e-15, 1.5260e-15, 1.4818e-15, 1.4389e-15,  // f12010
  1.3971e-15, 1.3566e-15, 1.3172e-15, 1.2790e-15, 1.2419e-15,  // f12020
  1.2058e-15, 1.1708e-15, 1.1368e-15, 1.1037e-15, 1.0716e-15,  // f12030
  1.0405e-15, 1.0102e-15, 9.8079e-16, 9.5224e-16, 9.2451e-16,  // f12040
  8.9758e-16, 8.7142e-16, 8.4602e-16, 8.2136e-16, 7.9740e-16,  // f12050
  7.7414e-16, 7.5154e-16, 7.2961e-16, 7.0830e-16, 6.8761e-16,  // f12060
  6.6752e-16, 6.4801e-16, 6.2906e-16, 6.1066e-16, 5.9280e-16,  // f12070
  5.7545e-16, 5.5860e-16, 5.4224e-16, 5.2636e-16, 5.1094e-16,  // f12080
  4.9596e-16};







/*
      --------------------------- CKD_MT_1.00 -------------------------

      11 February 2003
 
      Release Notes for MT_CKD_1.00
 
      Prepared by S. A. Clough,
                  AER Inc.,
                  131 Harwell Avenue
                  Lexington,  MA  02421
                  clough@aer.com

      This is the initial release of the MT_CKD water vapor continuum and
      represents the  first recomputation of the entire self and foreign
      broadened continuum since the original model was developed in the
      1980s.  This version of the continuum is implemented in the line-by-line
      model LBLRTM v7.01 and will be utilized in all related AER Radiative
      Transfer models.

      further information can be found under
      http://www.rtweb.aer.com/continuum_frame.html

      Transformation from original F77 code to C/C++ by 
      T. Kuhn, iup Bremen, August 2003 

      --------------------------- CKD_MT_1.00 -------------------------
 */


// H2O self continuum parameters at T=296 K
// date of last update: 11/18/02
// units of (CM**3/MOL)*1.E-20
const Numeric SL296_ckd_mt_100_v1  =   -20.0;
const Numeric SL296_ckd_mt_100_v2  = 20000.0;
const Numeric SL296_ckd_mt_100_dv  =    10.0;
const int     SL296_ckd_mt_100_npt =  2003;
const double  SL296_ckd_mt_100[SL296_ckd_mt_100_npt+addF77fields] = {
        0.000e0,    1.720e-01,  1.695e-01,
        1.700e-01,  1.695e-01,  1.720e-01,  1.680e-01,  1.687e-01,
        1.624e-01,  1.606e-01,  1.508e-01,  1.447e-01,  1.344e-01,
        1.214e-01,  1.133e-01,  1.009e-01,  9.217e-02,  8.297e-02,
        6.989e-02,  6.513e-02,  5.469e-02,  5.056e-02,  4.417e-02,
        3.779e-02,  3.484e-02,  2.994e-02,  2.720e-02,  2.325e-02,
        2.063e-02,  1.818e-02,  1.592e-02,  1.405e-02,  1.251e-02,
        1.080e-02,  9.647e-03,  8.424e-03,  7.519e-03,  6.555e-03,
        5.880e-03,  5.136e-03,  4.511e-03,  3.989e-03,  3.509e-03,
        3.114e-03,  2.740e-03,  2.446e-03,  2.144e-03,  1.895e-03,
        1.676e-03,  1.486e-03,  1.312e-03,  1.164e-03,  1.031e-03,
        9.129e-04,  8.106e-04,  7.213e-04,  6.400e-04,  5.687e-04,
        5.063e-04,  4.511e-04,  4.029e-04,  3.596e-04,  3.220e-04,
        2.889e-04,  2.597e-04,  2.337e-04,  2.108e-04,  1.907e-04,
        1.728e-04,  1.570e-04,  1.430e-04,  1.305e-04,  1.195e-04,
        1.097e-04,  1.009e-04,  9.307e-05,  8.604e-05,  7.971e-05,
        7.407e-05,  6.896e-05,  6.433e-05,  6.013e-05,  5.631e-05,
        5.283e-05,  4.963e-05,  4.669e-05,  4.398e-05,  4.148e-05,
        3.917e-05,  3.702e-05,  3.502e-05,  3.316e-05,  3.142e-05,
        2.978e-05,  2.825e-05,  2.681e-05,  2.546e-05,  2.419e-05,
        2.299e-05,  2.186e-05,  2.079e-05,  1.979e-05,  1.884e-05,
        1.795e-05,  1.711e-05,  1.633e-05,  1.559e-05,  1.490e-05,
        1.426e-05,  1.367e-05,  1.312e-05,  1.263e-05,  1.218e-05,
        1.178e-05,  1.143e-05,  1.112e-05,  1.088e-05,  1.070e-05,
        1.057e-05,  1.050e-05,  1.051e-05,  1.059e-05,  1.076e-05,
        1.100e-05,  1.133e-05,  1.180e-05,  1.237e-05,  1.308e-05,
        1.393e-05,  1.483e-05,  1.614e-05,  1.758e-05,  1.930e-05,
        2.123e-05,  2.346e-05,  2.647e-05,  2.930e-05,  3.279e-05,
        3.745e-05,  4.152e-05,  4.813e-05,  5.477e-05,  6.203e-05,
        7.331e-05,  8.056e-05,  9.882e-05,  1.050e-04,  1.210e-04,
        1.341e-04,  1.572e-04,  1.698e-04,  1.968e-04,  2.175e-04,
        2.431e-04,  2.735e-04,  2.867e-04,  3.190e-04,  3.371e-04,
        3.554e-04,  3.726e-04,  3.837e-04,  3.878e-04,  3.864e-04,
        3.858e-04,  3.841e-04,  3.852e-04,  3.815e-04,  3.762e-04,
        3.618e-04,  3.579e-04,  3.450e-04,  3.202e-04,  3.018e-04,
        2.785e-04,  2.602e-04,  2.416e-04,  2.097e-04,  1.939e-04,
        1.689e-04,  1.498e-04,  1.308e-04,  1.170e-04,  1.011e-04,
        9.237e-05,  7.909e-05,  7.006e-05,  6.112e-05,  5.401e-05,
        4.914e-05,  4.266e-05,  3.963e-05,  3.316e-05,  3.037e-05,
        2.598e-05,  2.294e-05,  2.066e-05,  1.813e-05,  1.583e-05,
        1.423e-05,  1.247e-05,  1.116e-05,  9.760e-06,  8.596e-06,
        7.720e-06,  6.825e-06,  6.108e-06,  5.366e-06,  4.733e-06,
        4.229e-06,  3.731e-06,  3.346e-06,  2.972e-06,  2.628e-06,
        2.356e-06,  2.102e-06,  1.878e-06,  1.678e-06,  1.507e-06,
        1.348e-06,  1.210e-06,  1.089e-06,  9.806e-07,  8.857e-07,
        8.004e-07,  7.261e-07,  6.599e-07,  6.005e-07,  5.479e-07,
        5.011e-07,  4.595e-07,  4.219e-07,  3.885e-07,  3.583e-07,
        3.314e-07,  3.071e-07,  2.852e-07,  2.654e-07,  2.474e-07,
        2.311e-07,  2.162e-07,  2.026e-07,  1.902e-07,  1.788e-07,
        1.683e-07,  1.587e-07,  1.497e-07,  1.415e-07,  1.338e-07,
        1.266e-07,  1.200e-07,  1.138e-07,  1.080e-07,  1.027e-07,
        9.764e-08,  9.296e-08,  8.862e-08,  8.458e-08,  8.087e-08,
        7.744e-08,  7.429e-08,  7.145e-08,  6.893e-08,  6.664e-08,
        6.468e-08,  6.322e-08,  6.162e-08,  6.070e-08,  5.992e-08,
        5.913e-08,  5.841e-08,  5.796e-08,  5.757e-08,  5.746e-08,
        5.731e-08,  5.679e-08,  5.577e-08,  5.671e-08,  5.656e-08,
        5.594e-08,  5.593e-08,  5.602e-08,  5.620e-08,  5.693e-08,
        5.725e-08,  5.858e-08,  6.037e-08,  6.249e-08,  6.535e-08,
        6.899e-08,  7.356e-08,  7.918e-08,  8.618e-08,  9.385e-08,
        1.039e-07,  1.158e-07,  1.290e-07,  1.437e-07,  1.650e-07,
        1.871e-07,  2.121e-07,  2.427e-07,  2.773e-07,  3.247e-07,
        3.677e-07,  4.037e-07,  4.776e-07,  5.101e-07,  6.214e-07,
        6.936e-07,  7.581e-07,  8.486e-07,  9.355e-07,  9.942e-07,
        1.063e-06,  1.123e-06,  1.191e-06,  1.215e-06,  1.247e-06,
        1.260e-06,  1.271e-06,  1.284e-06,  1.317e-06,  1.323e-06,
        1.349e-06,  1.353e-06,  1.362e-06,  1.344e-06,  1.329e-06,
        1.336e-06,  1.327e-06,  1.325e-06,  1.359e-06,  1.374e-06,
        1.415e-06,  1.462e-06,  1.526e-06,  1.619e-06,  1.735e-06,
        1.863e-06,  2.034e-06,  2.265e-06,  2.482e-06,  2.756e-06,
        3.103e-06,  3.466e-06,  3.832e-06,  4.378e-06,  4.913e-06,
        5.651e-06,  6.311e-06,  7.169e-06,  8.057e-06,  9.253e-06,
        1.047e-05,  1.212e-05,  1.360e-05,  1.569e-05,  1.776e-05,
        2.020e-05,  2.281e-05,  2.683e-05,  2.994e-05,  3.488e-05,
        3.896e-05,  4.499e-05,  5.175e-05,  6.035e-05,  6.340e-05,
        7.281e-05,  7.923e-05,  8.348e-05,  9.631e-05,  1.044e-04,
        1.102e-04,  1.176e-04,  1.244e-04,  1.283e-04,  1.326e-04,
        1.400e-04,  1.395e-04,  1.387e-04,  1.363e-04,  1.314e-04,
        1.241e-04,  1.228e-04,  1.148e-04,  1.086e-04,  1.018e-04,
        8.890e-05,  8.316e-05,  7.292e-05,  6.452e-05,  5.625e-05,
        5.045e-05,  4.380e-05,  3.762e-05,  3.290e-05,  2.836e-05,
        2.485e-05,  2.168e-05,  1.895e-05,  1.659e-05,  1.453e-05,
        1.282e-05,  1.132e-05,  1.001e-05,  8.836e-06,  7.804e-06,
        6.922e-06,  6.116e-06,  5.429e-06,  4.824e-06,  4.278e-06,
        3.788e-06,  3.371e-06,  2.985e-06,  2.649e-06,  2.357e-06,
        2.090e-06,  1.858e-06,  1.647e-06,  1.462e-06,  1.299e-06,
        1.155e-06,  1.028e-06,  9.142e-07,  8.132e-07,  7.246e-07,
        6.451e-07,  5.764e-07,  5.151e-07,  4.603e-07,  4.121e-07,
        3.694e-07,  3.318e-07,  2.985e-07,  2.690e-07,  2.428e-07,
        2.197e-07,  1.992e-07,  1.810e-07,  1.649e-07,  1.506e-07,
        1.378e-07,  1.265e-07,  1.163e-07,  1.073e-07,  9.918e-08,
        9.191e-08,  8.538e-08,  7.949e-08,  7.419e-08,  6.940e-08,
        6.508e-08,  6.114e-08,  5.761e-08,  5.437e-08,  5.146e-08,
        4.890e-08,  4.636e-08,  4.406e-08,  4.201e-08,  4.015e-08,
        3.840e-08,  3.661e-08,  3.510e-08,  3.377e-08,  3.242e-08,
        3.130e-08,  3.015e-08,  2.918e-08,  2.830e-08,  2.758e-08,
        2.707e-08,  2.656e-08,  2.619e-08,  2.609e-08,  2.615e-08,
        2.630e-08,  2.675e-08,  2.745e-08,  2.842e-08,  2.966e-08,
        3.125e-08,  3.318e-08,  3.565e-08,  3.850e-08,  4.191e-08,
        4.590e-08,  5.059e-08,  5.607e-08,  6.239e-08,  6.958e-08,
        7.796e-08,  8.773e-08,  9.880e-08,  1.114e-07,  1.258e-07,
        1.422e-07,  1.610e-07,  1.822e-07,  2.060e-07,  2.337e-07,
        2.645e-07,  2.996e-07,  3.393e-07,  3.843e-07,  4.363e-07,
        4.935e-07,  5.607e-07,  6.363e-07,  7.242e-07,  8.230e-07,
        9.411e-07,  1.071e-06,  1.232e-06,  1.402e-06,  1.600e-06,
        1.820e-06,  2.128e-06,  2.386e-06,  2.781e-06,  3.242e-06,
        3.653e-06,  4.323e-06,  4.747e-06,  5.321e-06,  5.919e-06,
        6.681e-06,  7.101e-06,  7.983e-06,  8.342e-06,  8.741e-06,
        9.431e-06,  9.952e-06,  1.026e-05,  1.055e-05,  1.095e-05,
        1.095e-05,  1.087e-05,  1.056e-05,  1.026e-05,  9.715e-06,
        9.252e-06,  8.452e-06,  7.958e-06,  7.268e-06,  6.295e-06,
        6.003e-06,  5.000e-06,  4.591e-06,  3.983e-06,  3.479e-06,
        3.058e-06,  2.667e-06,  2.293e-06,  1.995e-06,  1.747e-06,
        1.517e-06,  1.335e-06,  1.165e-06,  1.028e-06,  9.007e-07,
        7.956e-07,  7.015e-07,  6.192e-07,  5.491e-07,  4.859e-07,
        4.297e-07,  3.799e-07,  3.380e-07,  3.002e-07,  2.659e-07,
        2.366e-07,  2.103e-07,  1.861e-07,  1.655e-07,  1.469e-07,
        1.309e-07,  1.162e-07,  1.032e-07,  9.198e-08,  8.181e-08,
        7.294e-08,  6.516e-08,  5.787e-08,  5.163e-08,  4.612e-08,
        4.119e-08,  3.695e-08,  3.308e-08,  2.976e-08,  2.670e-08,
        2.407e-08,  2.171e-08,  1.965e-08,  1.780e-08,  1.617e-08,
        1.470e-08,  1.341e-08,  1.227e-08,  1.125e-08,  1.033e-08,
        9.524e-09,  8.797e-09,  8.162e-09,  7.565e-09,  7.040e-09,
        6.560e-09,  6.129e-09,  5.733e-09,  5.376e-09,  5.043e-09,
        4.750e-09,  4.466e-09,  4.211e-09,  3.977e-09,  3.759e-09,
        3.558e-09,  3.373e-09,  3.201e-09,  3.043e-09,  2.895e-09,
        2.760e-09,  2.635e-09,  2.518e-09,  2.411e-09,  2.314e-09,
        2.230e-09,  2.151e-09,  2.087e-09,  2.035e-09,  1.988e-09,
        1.946e-09,  1.927e-09,  1.916e-09,  1.916e-09,  1.933e-09,
        1.966e-09,  2.018e-09,  2.090e-09,  2.182e-09,  2.299e-09,
        2.442e-09,  2.623e-09,  2.832e-09,  3.079e-09,  3.368e-09,
        3.714e-09,  4.104e-09,  4.567e-09,  5.091e-09,  5.701e-09,
        6.398e-09,  7.194e-09,  8.127e-09,  9.141e-09,  1.035e-08,
        1.177e-08,  1.338e-08,  1.508e-08,  1.711e-08,  1.955e-08,
        2.216e-08,  2.534e-08,  2.871e-08,  3.291e-08,  3.711e-08,
        4.285e-08,  4.868e-08,  5.509e-08,  6.276e-08,  7.262e-08,
        8.252e-08,  9.400e-08,  1.064e-07,  1.247e-07,  1.411e-07,
        1.626e-07,  1.827e-07,  2.044e-07,  2.284e-07,  2.452e-07,
        2.854e-07,  3.026e-07,  3.278e-07,  3.474e-07,  3.693e-07,
        3.930e-07,  4.104e-07,  4.220e-07,  4.439e-07,  4.545e-07,
        4.778e-07,  4.812e-07,  5.018e-07,  4.899e-07,  5.075e-07,
        5.073e-07,  5.171e-07,  5.131e-07,  5.250e-07,  5.617e-07,
        5.846e-07,  6.239e-07,  6.696e-07,  7.398e-07,  8.073e-07,
        9.150e-07,  1.009e-06,  1.116e-06,  1.264e-06,  1.439e-06,
        1.644e-06,  1.856e-06,  2.147e-06,  2.317e-06,  2.713e-06,
        2.882e-06,  2.990e-06,  3.489e-06,  3.581e-06,  4.033e-06,
        4.260e-06,  4.543e-06,  4.840e-06,  4.826e-06,  5.013e-06,
        5.252e-06,  5.277e-06,  5.306e-06,  5.236e-06,  5.123e-06,
        5.171e-06,  4.843e-06,  4.615e-06,  4.385e-06,  3.970e-06,
        3.693e-06,  3.231e-06,  2.915e-06,  2.495e-06,  2.144e-06,
        1.910e-06,  1.639e-06,  1.417e-06,  1.226e-06,  1.065e-06,
        9.290e-07,  8.142e-07,  7.161e-07,  6.318e-07,  5.581e-07,
        4.943e-07,  4.376e-07,  3.884e-07,  3.449e-07,  3.060e-07,
        2.712e-07,  2.412e-07,  2.139e-07,  1.903e-07,  1.689e-07,
        1.499e-07,  1.331e-07,  1.183e-07,  1.050e-07,  9.362e-08,
        8.306e-08,  7.403e-08,  6.578e-08,  5.853e-08,  5.216e-08,
        4.632e-08,  4.127e-08,  3.678e-08,  3.279e-08,  2.923e-08,
        2.612e-08,  2.339e-08,  2.094e-08,  1.877e-08,  1.686e-08,
        1.516e-08,  1.366e-08,  1.234e-08,  1.114e-08,  1.012e-08,
        9.182e-09,  8.362e-09,  7.634e-09,  6.981e-09,  6.406e-09,
        5.888e-09,  5.428e-09,  5.021e-09,  4.650e-09,  4.326e-09,
        4.033e-09,  3.770e-09,  3.536e-09,  3.327e-09,  3.141e-09,
        2.974e-09,  2.825e-09,  2.697e-09,  2.584e-09,  2.488e-09,
        2.406e-09,  2.340e-09,  2.292e-09,  2.259e-09,  2.244e-09,
        2.243e-09,  2.272e-09,  2.310e-09,  2.378e-09,  2.454e-09,
        2.618e-09,  2.672e-09,  2.831e-09,  3.050e-09,  3.225e-09,
        3.425e-09,  3.677e-09,  3.968e-09,  4.221e-09,  4.639e-09,
        4.960e-09,  5.359e-09,  5.649e-09,  6.230e-09,  6.716e-09,
        7.218e-09,  7.746e-09,  7.988e-09,  8.627e-09,  8.999e-09,
        9.442e-09,  9.820e-09,  1.015e-08,  1.060e-08,  1.079e-08,
        1.109e-08,  1.137e-08,  1.186e-08,  1.180e-08,  1.187e-08,
        1.194e-08,  1.192e-08,  1.224e-08,  1.245e-08,  1.246e-08,
        1.318e-08,  1.377e-08,  1.471e-08,  1.582e-08,  1.713e-08,
        1.853e-08,  2.063e-08,  2.270e-08,  2.567e-08,  2.891e-08,
        3.264e-08,  3.744e-08,  4.286e-08,  4.915e-08,  5.623e-08,
        6.336e-08,  7.293e-08,  8.309e-08,  9.319e-08,  1.091e-07,
        1.243e-07,  1.348e-07,  1.449e-07,  1.620e-07,  1.846e-07,
        1.937e-07,  2.040e-07,  2.179e-07,  2.298e-07,  2.433e-07,
        2.439e-07,  2.464e-07,  2.611e-07,  2.617e-07,  2.582e-07,
        2.453e-07,  2.401e-07,  2.349e-07,  2.203e-07,  2.066e-07,
        1.939e-07,  1.780e-07,  1.558e-07,  1.391e-07,  1.203e-07,
        1.048e-07,  9.464e-08,  8.306e-08,  7.239e-08,  6.317e-08,
        5.520e-08,  4.847e-08,  4.282e-08,  3.796e-08,  3.377e-08,
        2.996e-08,  2.678e-08,  2.400e-08,  2.134e-08,  1.904e-08,
        1.705e-08,  1.523e-08,  1.350e-08,  1.204e-08,  1.070e-08,
        9.408e-09,  8.476e-09,  7.470e-09,  6.679e-09,  5.929e-09,
        5.267e-09,  4.711e-09,  4.172e-09,  3.761e-09,  3.288e-09,
        2.929e-09,  2.609e-09,  2.315e-09,  2.042e-09,  1.844e-09,
        1.640e-09,  1.470e-09,  1.310e-09,  1.176e-09,  1.049e-09,
        9.377e-10,  8.462e-10,  7.616e-10,  6.854e-10,  6.191e-10,
        5.596e-10,  5.078e-10,  4.611e-10,  4.197e-10,  3.830e-10,
        3.505e-10,  3.215e-10,  2.956e-10,  2.726e-10,  2.521e-10,
        2.338e-10,  2.173e-10,  2.026e-10,  1.895e-10,  1.777e-10,
        1.672e-10,  1.579e-10,  1.496e-10,  1.423e-10,  1.358e-10,
        1.302e-10,  1.254e-10,  1.216e-10,  1.187e-10,  1.163e-10,
        1.147e-10,  1.145e-10,  1.150e-10,  1.170e-10,  1.192e-10,
        1.250e-10,  1.298e-10,  1.345e-10,  1.405e-10,  1.538e-10,
        1.648e-10,  1.721e-10,  1.872e-10,  1.968e-10,  2.089e-10,
        2.172e-10,  2.317e-10,  2.389e-10,  2.503e-10,  2.585e-10,
        2.686e-10,  2.800e-10,  2.895e-10,  3.019e-10,  3.037e-10,
        3.076e-10,  3.146e-10,  3.198e-10,  3.332e-10,  3.397e-10,
        3.540e-10,  3.667e-10,  3.895e-10,  4.071e-10,  4.565e-10,
        4.983e-10,  5.439e-10,  5.968e-10,  6.676e-10,  7.456e-10,
        8.405e-10,  9.478e-10,  1.064e-09,  1.218e-09,  1.386e-09,
        1.581e-09,  1.787e-09,  2.032e-09,  2.347e-09,  2.677e-09,
        3.008e-09,  3.544e-09,  4.056e-09,  4.687e-09,  5.331e-09,
        6.227e-09,  6.854e-09,  8.139e-09,  8.945e-09,  9.865e-09,
        1.125e-08,  1.178e-08,  1.364e-08,  1.436e-08,  1.540e-08,
        1.672e-08,  1.793e-08,  1.906e-08,  2.036e-08,  2.144e-08,
        2.292e-08,  2.371e-08,  2.493e-08,  2.606e-08,  2.706e-08,
        2.866e-08,  3.036e-08,  3.136e-08,  3.405e-08,  3.665e-08,
        3.837e-08,  4.229e-08,  4.748e-08,  5.320e-08,  5.763e-08,
        6.677e-08,  7.216e-08,  7.716e-08,  8.958e-08,  9.419e-08,
        1.036e-07,  1.108e-07,  1.189e-07,  1.246e-07,  1.348e-07,
        1.310e-07,  1.361e-07,  1.364e-07,  1.363e-07,  1.343e-07,
        1.293e-07,  1.254e-07,  1.235e-07,  1.158e-07,  1.107e-07,
        9.961e-08,  9.011e-08,  7.910e-08,  6.916e-08,  6.338e-08,
        5.564e-08,  4.827e-08,  4.198e-08,  3.695e-08,  3.276e-08,
        2.929e-08,  2.633e-08,  2.391e-08,  2.192e-08,  2.021e-08,
        1.890e-08,  1.772e-08,  1.667e-08,  1.603e-08,  1.547e-08,
        1.537e-08,  1.492e-08,  1.515e-08,  1.479e-08,  1.450e-08,
        1.513e-08,  1.495e-08,  1.529e-08,  1.565e-08,  1.564e-08,
        1.553e-08,  1.569e-08,  1.584e-08,  1.570e-08,  1.538e-08,
        1.513e-08,  1.472e-08,  1.425e-08,  1.349e-08,  1.328e-08,
        1.249e-08,  1.170e-08,  1.077e-08,  9.514e-09,  8.614e-09,
        7.460e-09,  6.621e-09,  5.775e-09,  5.006e-09,  4.308e-09,
        3.747e-09,  3.240e-09,  2.840e-09,  2.481e-09,  2.184e-09,
        1.923e-09,  1.710e-09,  1.504e-09,  1.334e-09,  1.187e-09,
        1.053e-09,  9.367e-10,  8.306e-10,  7.419e-10,  6.630e-10,
        5.918e-10,  5.277e-10,  4.717e-10,  4.222e-10,  3.783e-10,
        3.390e-10,  3.036e-10,  2.729e-10,  2.455e-10,  2.211e-10,
        1.995e-10,  1.804e-10,  1.635e-10,  1.485e-10,  1.355e-10,
        1.240e-10,  1.139e-10,  1.051e-10,  9.757e-11,  9.114e-11,
        8.577e-11,  8.139e-11,  7.792e-11,  7.520e-11,  7.390e-11,
        7.311e-11,  7.277e-11,  7.482e-11,  7.698e-11,  8.162e-11,
        8.517e-11,  8.968e-11,  9.905e-11,  1.075e-10,  1.187e-10,
        1.291e-10,  1.426e-10,  1.573e-10,  1.734e-10,  1.905e-10,
        2.097e-10,  2.280e-10,  2.473e-10,  2.718e-10,  2.922e-10,
        3.128e-10,  3.361e-10,  3.641e-10,  3.910e-10,  4.196e-10,
        4.501e-10,  4.932e-10,  5.258e-10,  5.755e-10,  6.253e-10,
        6.664e-10,  7.344e-10,  7.985e-10,  8.877e-10,  1.005e-09,
        1.118e-09,  1.251e-09,  1.428e-09,  1.610e-09,  1.888e-09,
        2.077e-09,  2.331e-09,  2.751e-09,  3.061e-09,  3.522e-09,
        3.805e-09,  4.181e-09,  4.575e-09,  5.167e-09,  5.634e-09,
        6.007e-09,  6.501e-09,  6.829e-09,  7.211e-09,  7.262e-09,
        7.696e-09,  7.832e-09,  7.799e-09,  7.651e-09,  7.304e-09,
        7.150e-09,  6.977e-09,  6.603e-09,  6.209e-09,  5.690e-09,
        5.432e-09,  4.764e-09,  4.189e-09,  3.640e-09,  3.203e-09,
        2.848e-09,  2.510e-09,  2.194e-09,  1.946e-09,  1.750e-09,
        1.567e-09,  1.426e-09,  1.302e-09,  1.197e-09,  1.109e-09,
        1.035e-09,  9.719e-10,  9.207e-10,  8.957e-10,  8.578e-10,
        8.262e-10,  8.117e-10,  7.987e-10,  7.875e-10,  7.741e-10,
        7.762e-10,  7.537e-10,  7.424e-10,  7.474e-10,  7.294e-10,
        7.216e-10,  7.233e-10,  7.075e-10,  6.892e-10,  6.618e-10,
        6.314e-10,  6.208e-10,  5.689e-10,  5.550e-10,  4.984e-10,
        4.600e-10,  4.078e-10,  3.879e-10,  3.459e-10,  2.982e-10,
        2.626e-10,  2.329e-10,  1.988e-10,  1.735e-10,  1.487e-10,
        1.297e-10,  1.133e-10,  9.943e-11,  8.736e-11,  7.726e-11,
        6.836e-11,  6.053e-11,  5.384e-11,  4.789e-11,  4.267e-11,
        3.804e-11,  3.398e-11,  3.034e-11,  2.710e-11,  2.425e-11,
        2.173e-11,  1.950e-11,  1.752e-11,  1.574e-11,  1.418e-11,
        1.278e-11,  1.154e-11,  1.044e-11,  9.463e-12,  8.602e-12,
        7.841e-12,  7.171e-12,  6.584e-12,  6.073e-12,  5.631e-12,
        5.254e-12,  4.937e-12,  4.679e-12,  4.476e-12,  4.328e-12,
        4.233e-12,  4.194e-12,  4.211e-12,  4.286e-12,  4.424e-12,
        4.628e-12,  4.906e-12,  5.262e-12,  5.708e-12,  6.254e-12,
        6.914e-12,  7.714e-12,  8.677e-12,  9.747e-12,  1.101e-11,
        1.256e-11,  1.409e-11,  1.597e-11,  1.807e-11,  2.034e-11,
        2.316e-11,  2.622e-11,  2.962e-11,  3.369e-11,  3.819e-11,
        4.329e-11,  4.932e-11,  5.589e-11,  6.364e-11,  7.284e-11,
        8.236e-11,  9.447e-11,  1.078e-10,  1.229e-10,  1.417e-10,
        1.614e-10,  1.843e-10,  2.107e-10,  2.406e-10,  2.728e-10,
        3.195e-10,  3.595e-10,  4.153e-10,  4.736e-10,  5.410e-10,
        6.088e-10,  6.769e-10,  7.691e-10,  8.545e-10,  9.621e-10,
        1.047e-09,  1.161e-09,  1.296e-09,  1.424e-09,  1.576e-09,
        1.739e-09,  1.893e-09,  2.080e-09,  2.336e-09,  2.604e-09,
        2.760e-09,  3.001e-09,  3.365e-09,  3.550e-09,  3.895e-09,
        4.183e-09,  4.614e-09,  4.846e-09,  5.068e-09,  5.427e-09,
        5.541e-09,  5.864e-09,  5.997e-09,  5.997e-09,  6.061e-09,
        5.944e-09,  5.855e-09,  5.661e-09,  5.523e-09,  5.374e-09,
        4.940e-09,  4.688e-09,  4.170e-09,  3.913e-09,  3.423e-09,
        2.997e-09,  2.598e-09,  2.253e-09,  1.946e-09,  1.710e-09,
        1.507e-09,  1.336e-09,  1.190e-09,  1.068e-09,  9.623e-10,
        8.772e-10,  8.007e-10,  7.420e-10,  6.884e-10,  6.483e-10,
        6.162e-10,  5.922e-10,  5.688e-10,  5.654e-10,  5.637e-10,
        5.701e-10,  5.781e-10,  5.874e-10,  6.268e-10,  6.357e-10,
        6.525e-10,  7.137e-10,  7.441e-10,  8.024e-10,  8.485e-10,
        9.143e-10,  9.536e-10,  9.717e-10,  1.018e-09,  1.042e-09,
        1.054e-09,  1.092e-09,  1.079e-09,  1.064e-09,  1.043e-09,
        1.020e-09,  9.687e-10,  9.273e-10,  9.208e-10,  9.068e-10,
        7.687e-10,  7.385e-10,  6.595e-10,  5.870e-10,  5.144e-10,
        4.417e-10,  3.804e-10,  3.301e-10,  2.866e-10,  2.509e-10,
        2.202e-10,  1.947e-10,  1.719e-10,  1.525e-10,  1.361e-10,
        1.210e-10,  1.084e-10,  9.800e-11,  8.801e-11,  7.954e-11,
        7.124e-11,  6.335e-11,  5.760e-11,  5.132e-11,  4.601e-11,
        4.096e-11,  3.657e-11,  3.250e-11,  2.909e-11,  2.587e-11,
        2.297e-11,  2.050e-11,  1.828e-11,  1.632e-11,  1.462e-11,
        1.314e-11,  1.185e-11,  1.073e-11,  9.760e-12,  8.922e-12,
        8.206e-12,  7.602e-12,  7.100e-12,  6.694e-12,  6.378e-12,
        6.149e-12,  6.004e-12,  5.941e-12,  5.962e-12,  6.069e-12,
        6.265e-12,  6.551e-12,  6.935e-12,  7.457e-12,  8.074e-12,
        8.811e-12,  9.852e-12,  1.086e-11,  1.207e-11,  1.361e-11,
        1.553e-11,  1.737e-11,  1.930e-11,  2.175e-11,  2.410e-11,
        2.706e-11,  3.023e-11,  3.313e-11,  3.657e-11,  4.118e-11,
        4.569e-11,  5.025e-11,  5.660e-11,  6.231e-11,  6.881e-11,
        7.996e-11,  8.526e-11,  9.694e-11,  1.106e-10,  1.222e-10,
        1.355e-10,  1.525e-10,  1.775e-10,  1.924e-10,  2.181e-10,
        2.379e-10,  2.662e-10,  2.907e-10,  3.154e-10,  3.366e-10,
        3.579e-10,  3.858e-10,  4.046e-10,  4.196e-10,  4.166e-10,
        4.457e-10,  4.466e-10,  4.404e-10,  4.337e-10,  4.150e-10,
        4.083e-10,  3.910e-10,  3.723e-10,  3.514e-10,  3.303e-10,
        2.847e-10,  2.546e-10,  2.230e-10,  1.994e-10,  1.733e-10,
        1.488e-10,  1.297e-10,  1.144e-10,  1.004e-10,  8.741e-11,
        7.928e-11,  7.034e-11,  6.323e-11,  5.754e-11,  5.250e-11,
        4.850e-11,  4.502e-11,  4.286e-11,  4.028e-11,  3.899e-11,
        3.824e-11,  3.761e-11,  3.804e-11,  3.839e-11,  3.845e-11,
        4.244e-11,  4.382e-11,  4.582e-11,  4.847e-11,  5.209e-11,
        5.384e-11,  5.887e-11,  6.371e-11,  6.737e-11,  7.168e-11,
        7.415e-11,  7.827e-11,  8.037e-11,  8.120e-11,  8.071e-11,
        8.008e-11,  7.851e-11,  7.544e-11,  7.377e-11,  7.173e-11,
        6.801e-11,  6.267e-11,  5.727e-11,  5.288e-11,  4.853e-11,
        4.082e-11,  3.645e-11,  3.136e-11,  2.672e-11,  2.304e-11,
        1.986e-11,  1.725e-11,  1.503e-11,  1.315e-11,  1.153e-11,
        1.014e-11,  8.942e-12,  7.901e-12,  6.993e-12,  6.199e-12,
        5.502e-12,  4.890e-12,  4.351e-12,  3.878e-12,  3.461e-12,
        3.094e-12,  2.771e-12,  2.488e-12,  2.241e-12,  2.025e-12,
        1.838e-12,  1.677e-12,  1.541e-12,  1.427e-12,  1.335e-12,
        1.262e-12,  1.209e-12,  1.176e-12,  1.161e-12,  1.165e-12,
        1.189e-12,  1.234e-12,  1.300e-12,  1.389e-12,  1.503e-12,
        1.644e-12,  1.814e-12,  2.017e-12,  2.255e-12,  2.534e-12,
        2.858e-12,  3.231e-12,  3.661e-12,  4.153e-12,  4.717e-12,
        5.360e-12,  6.094e-12,  6.930e-12,  7.882e-12,  8.966e-12,
        1.020e-11,  1.162e-11,  1.324e-11,  1.510e-11,  1.720e-11,
        1.965e-11,  2.237e-11,  2.560e-11,  2.927e-11,  3.371e-11,
        3.842e-11,  4.429e-11,  5.139e-11,  5.798e-11,  6.697e-11,
        7.626e-11,  8.647e-11,  1.022e-10,  1.136e-10,  1.300e-10,
        1.481e-10,  1.672e-10,  1.871e-10,  2.126e-10,  2.357e-10,
        2.583e-10,  2.997e-10,  3.289e-10,  3.702e-10,  4.012e-10,
        4.319e-10,  4.527e-10,  5.001e-10,  5.448e-10,  5.611e-10,
        5.760e-10,  5.965e-10,  6.079e-10,  6.207e-10,  6.276e-10,
        6.222e-10,  6.137e-10,  6.000e-10,  5.814e-10,  5.393e-10,
        5.350e-10,  4.947e-10,  4.629e-10,  4.117e-10,  3.712e-10,
        3.372e-10,  2.923e-10,  2.550e-10,  2.232e-10,  1.929e-10,
        1.679e-10,  1.460e-10,  1.289e-10,  1.130e-10,  9.953e-11,
        8.763e-11,  7.760e-11,  6.900e-11,  6.160e-11,  5.525e-11,
        4.958e-11,  4.489e-11,  4.072e-11,  3.728e-11,  3.438e-11,
        3.205e-11,  3.006e-11,  2.848e-11,  2.766e-11,  2.688e-11,
        2.664e-11,  2.670e-11,  2.696e-11,  2.786e-11,  2.861e-11,
        3.009e-11,  3.178e-11,  3.389e-11,  3.587e-11,  3.819e-11,
        4.054e-11,  4.417e-11,  4.703e-11,  5.137e-11,  5.460e-11,
        6.055e-11,  6.333e-11,  6.773e-11,  7.219e-11,  7.717e-11,
        8.131e-11,  8.491e-11,  8.574e-11,  9.010e-11,  9.017e-11,
        8.999e-11,  8.959e-11,  8.838e-11,  8.579e-11,  8.162e-11,
        8.098e-11,  7.472e-11,  7.108e-11,  6.559e-11,  5.994e-11,
        5.172e-11,  4.424e-11,  3.951e-11,  3.340e-11,  2.902e-11,
        2.541e-11,  2.215e-11,  1.945e-11,  1.716e-11,  1.503e-11,
        1.339e-11,  1.185e-11,  1.050e-11,  9.336e-12,  8.307e-12,
        7.312e-12,  6.550e-12,  5.836e-12,  5.178e-12,  4.600e-12,
        4.086e-12,  3.639e-12,  3.247e-12,  2.904e-12,  2.604e-12,
        2.341e-12,  2.112e-12,  1.914e-12,  1.744e-12,  1.598e-12,
        1.476e-12,  1.374e-12,  1.293e-12,  1.230e-12,  1.185e-12,
        1.158e-12,  1.147e-12,  1.154e-12,  1.177e-12,  1.219e-12,
        1.280e-12,  1.360e-12,  1.463e-12,  1.591e-12,  1.750e-12,
        1.940e-12,  2.156e-12,  2.430e-12,  2.748e-12,  3.052e-12,
        3.533e-12,  3.967e-12,  4.471e-12,  5.041e-12,  5.860e-12,
        6.664e-12,  7.522e-12,  8.342e-12,  9.412e-12,  1.072e-11,
        1.213e-11,  1.343e-11,  1.496e-11,  1.664e-11,  1.822e-11,
        2.029e-11,  2.233e-11,  2.457e-11,  2.709e-11,  2.928e-11,
        3.115e-11,  3.356e-11,  3.592e-11,  3.818e-11,  3.936e-11,
        4.061e-11,  4.149e-11,  4.299e-11,  4.223e-11,  4.251e-11,
        4.287e-11,  4.177e-11,  4.094e-11,  3.942e-11,  3.772e-11,
        3.614e-11,  3.394e-11,  3.222e-11,  2.791e-11,  2.665e-11,
        2.309e-11,  2.032e-11,  1.740e-11,  1.535e-11,  1.323e-11,
        1.151e-11,  9.803e-12,  8.650e-12,  7.540e-12,  6.619e-12,
        5.832e-12,  5.113e-12,  4.503e-12,  3.975e-12,  3.520e-12,
        3.112e-12,  2.797e-12,  2.500e-12,  2.240e-12,  2.013e-12,
        1.819e-12,  1.653e-12,  1.513e-12,  1.395e-12,  1.299e-12,
        1.225e-12,  1.168e-12,  1.124e-12,  1.148e-12,  1.107e-12,
        1.128e-12,  1.169e-12,  1.233e-12,  1.307e-12,  1.359e-12,
        1.543e-12,  1.686e-12,  1.794e-12,  2.028e-12,  2.210e-12,
        2.441e-12,  2.653e-12,  2.828e-12,  3.093e-12,  3.280e-12,
        3.551e-12,  3.677e-12,  3.803e-12,  3.844e-12,  4.068e-12,
        4.093e-12,  4.002e-12,  3.904e-12,  3.624e-12,  3.633e-12,
        3.622e-12,  3.443e-12,  3.184e-12,  2.934e-12,  2.476e-12,
        2.212e-12,  1.867e-12,  1.594e-12,  1.370e-12,  1.192e-12,
        1.045e-12,  9.211e-13,  8.170e-13,  7.290e-13,  6.550e-13,
        5.929e-13,  5.415e-13,  4.995e-13,  4.661e-13,  4.406e-13,
        4.225e-13,  4.116e-13,  4.075e-13,  4.102e-13,  4.198e-13,
        4.365e-13,  4.606e-13,  4.925e-13,  5.326e-13,  5.818e-13,
        6.407e-13,  7.104e-13,  7.920e-13,  8.868e-13,  9.964e-13,
        1.123e-12,  1.268e-12,  1.434e-12,  1.626e-12,  1.848e-12,
        2.107e-12,  2.422e-12,  2.772e-12,  3.145e-12,  3.704e-12,
        4.270e-12,  4.721e-12,  5.361e-12,  6.083e-12,  7.095e-12,
        7.968e-12,  9.228e-12,  1.048e-11,  1.187e-11,  1.336e-11,
        1.577e-11,  1.772e-11,  2.017e-11,  2.250e-11,  2.630e-11,
        2.911e-11,  3.356e-11,  3.820e-11,  4.173e-11,  4.811e-11,
        5.254e-11,  5.839e-11,  6.187e-11,  6.805e-11,  7.118e-11,
        7.369e-11,  7.664e-11,  7.794e-11,  7.947e-11,  8.036e-11,
        7.954e-11,  7.849e-11,  7.518e-11,  7.462e-11,  6.926e-11,
        6.531e-11,  6.197e-11,  5.421e-11,  4.777e-11,  4.111e-11,
        3.679e-11,  3.166e-11,  2.786e-11,  2.436e-11,  2.144e-11,
        1.859e-11,  1.628e-11,  1.414e-11,  1.237e-11,  1.093e-11,
        9.558e-12 };

// H2O self continuum parameters at T=260 K
// date of last update: 11/18/02
// units of (CM**3/MOL)*1.E-20
const Numeric SL260_ckd_mt_100_v1  =   -20.0;
const Numeric SL260_ckd_mt_100_v2  = 20000.0;
const Numeric SL260_ckd_mt_100_dv  =    10.0;
const int     SL260_ckd_mt_100_npt =  2003;
const double  SL260_ckd_mt_100[SL260_ckd_mt_100_npt+addF77fields] = {
        0.000e0,    2.749e-01,  2.732e-01,
        2.752e-01,  2.732e-01,  2.749e-01,  2.676e-01,  2.667e-01,
        2.545e-01,  2.497e-01,  2.327e-01,  2.218e-01,  2.036e-01,
        1.825e-01,  1.694e-01,  1.497e-01,  1.353e-01,  1.210e-01,
        1.014e-01,  9.405e-02,  7.848e-02,  7.195e-02,  6.246e-02,
        5.306e-02,  4.853e-02,  4.138e-02,  3.735e-02,  3.171e-02,
        2.785e-02,  2.431e-02,  2.111e-02,  1.845e-02,  1.640e-02,
        1.405e-02,  1.255e-02,  1.098e-02,  9.797e-03,  8.646e-03,
        7.779e-03,  6.898e-03,  6.099e-03,  5.453e-03,  4.909e-03,
        4.413e-03,  3.959e-03,  3.581e-03,  3.199e-03,  2.871e-03,
        2.583e-03,  2.330e-03,  2.086e-03,  1.874e-03,  1.684e-03,
        1.512e-03,  1.361e-03,  1.225e-03,  1.100e-03,  9.890e-04,
        8.916e-04,  8.039e-04,  7.256e-04,  6.545e-04,  5.918e-04,
        5.359e-04,  4.867e-04,  4.426e-04,  4.033e-04,  3.682e-04,
        3.366e-04,  3.085e-04,  2.833e-04,  2.605e-04,  2.403e-04,
        2.221e-04,  2.055e-04,  1.908e-04,  1.774e-04,  1.653e-04,
        1.544e-04,  1.443e-04,  1.351e-04,  1.267e-04,  1.190e-04,
        1.119e-04,  1.053e-04,  9.922e-05,  9.355e-05,  8.831e-05,
        8.339e-05,  7.878e-05,  7.449e-05,  7.043e-05,  6.664e-05,
        6.307e-05,  5.969e-05,  5.654e-05,  5.357e-05,  5.075e-05,
        4.810e-05,  4.560e-05,  4.322e-05,  4.102e-05,  3.892e-05,
        3.696e-05,  3.511e-05,  3.339e-05,  3.177e-05,  3.026e-05,
        2.886e-05,  2.756e-05,  2.636e-05,  2.527e-05,  2.427e-05,
        2.337e-05,  2.257e-05,  2.185e-05,  2.127e-05,  2.080e-05,
        2.041e-05,  2.013e-05,  2.000e-05,  1.997e-05,  2.009e-05,
        2.031e-05,  2.068e-05,  2.124e-05,  2.189e-05,  2.267e-05,
        2.364e-05,  2.463e-05,  2.618e-05,  2.774e-05,  2.937e-05,
        3.144e-05,  3.359e-05,  3.695e-05,  4.002e-05,  4.374e-05,
        4.947e-05,  5.431e-05,  6.281e-05,  7.169e-05,  8.157e-05,
        9.728e-05,  1.079e-04,  1.337e-04,  1.442e-04,  1.683e-04,
        1.879e-04,  2.223e-04,  2.425e-04,  2.838e-04,  3.143e-04,
        3.527e-04,  4.012e-04,  4.237e-04,  4.747e-04,  5.057e-04,
        5.409e-04,  5.734e-04,  5.944e-04,  6.077e-04,  6.175e-04,
        6.238e-04,  6.226e-04,  6.248e-04,  6.192e-04,  6.098e-04,
        5.818e-04,  5.709e-04,  5.465e-04,  5.043e-04,  4.699e-04,
        4.294e-04,  3.984e-04,  3.672e-04,  3.152e-04,  2.883e-04,
        2.503e-04,  2.211e-04,  1.920e-04,  1.714e-04,  1.485e-04,
        1.358e-04,  1.156e-04,  1.021e-04,  8.887e-05,  7.842e-05,
        7.120e-05,  6.186e-05,  5.730e-05,  4.792e-05,  4.364e-05,
        3.720e-05,  3.280e-05,  2.946e-05,  2.591e-05,  2.261e-05,
        2.048e-05,  1.813e-05,  1.630e-05,  1.447e-05,  1.282e-05,
        1.167e-05,  1.041e-05,  9.449e-06,  8.510e-06,  7.596e-06,
        6.961e-06,  6.272e-06,  5.728e-06,  5.198e-06,  4.667e-06,
        4.288e-06,  3.897e-06,  3.551e-06,  3.235e-06,  2.952e-06,
        2.688e-06,  2.449e-06,  2.241e-06,  2.050e-06,  1.879e-06,
        1.722e-06,  1.582e-06,  1.456e-06,  1.339e-06,  1.236e-06,
        1.144e-06,  1.060e-06,  9.830e-07,  9.149e-07,  8.535e-07,
        7.973e-07,  7.466e-07,  6.999e-07,  6.574e-07,  6.180e-07,
        5.821e-07,  5.487e-07,  5.180e-07,  4.896e-07,  4.631e-07,
        4.386e-07,  4.160e-07,  3.945e-07,  3.748e-07,  3.562e-07,
        3.385e-07,  3.222e-07,  3.068e-07,  2.922e-07,  2.788e-07,
        2.659e-07,  2.539e-07,  2.425e-07,  2.318e-07,  2.219e-07,
        2.127e-07,  2.039e-07,  1.958e-07,  1.885e-07,  1.818e-07,
        1.758e-07,  1.711e-07,  1.662e-07,  1.630e-07,  1.605e-07,
        1.580e-07,  1.559e-07,  1.545e-07,  1.532e-07,  1.522e-07,
        1.510e-07,  1.495e-07,  1.465e-07,  1.483e-07,  1.469e-07,
        1.448e-07,  1.444e-07,  1.436e-07,  1.426e-07,  1.431e-07,
        1.425e-07,  1.445e-07,  1.477e-07,  1.515e-07,  1.567e-07,
        1.634e-07,  1.712e-07,  1.802e-07,  1.914e-07,  2.024e-07,
        2.159e-07,  2.295e-07,  2.461e-07,  2.621e-07,  2.868e-07,
        3.102e-07,  3.394e-07,  3.784e-07,  4.223e-07,  4.864e-07,
        5.501e-07,  6.039e-07,  7.193e-07,  7.728e-07,  9.514e-07,
        1.073e-06,  1.180e-06,  1.333e-06,  1.472e-06,  1.566e-06,
        1.677e-06,  1.784e-06,  1.904e-06,  1.953e-06,  2.020e-06,
        2.074e-06,  2.128e-06,  2.162e-06,  2.219e-06,  2.221e-06,
        2.249e-06,  2.239e-06,  2.235e-06,  2.185e-06,  2.141e-06,
        2.124e-06,  2.090e-06,  2.068e-06,  2.100e-06,  2.104e-06,
        2.142e-06,  2.181e-06,  2.257e-06,  2.362e-06,  2.500e-06,
        2.664e-06,  2.884e-06,  3.189e-06,  3.480e-06,  3.847e-06,
        4.313e-06,  4.790e-06,  5.250e-06,  5.989e-06,  6.692e-06,
        7.668e-06,  8.520e-06,  9.606e-06,  1.073e-05,  1.225e-05,
        1.377e-05,  1.582e-05,  1.761e-05,  2.029e-05,  2.284e-05,
        2.602e-05,  2.940e-05,  3.483e-05,  3.928e-05,  4.618e-05,
        5.240e-05,  6.132e-05,  7.183e-05,  8.521e-05,  9.111e-05,
        1.070e-04,  1.184e-04,  1.264e-04,  1.475e-04,  1.612e-04,
        1.704e-04,  1.818e-04,  1.924e-04,  1.994e-04,  2.061e-04,
        2.180e-04,  2.187e-04,  2.200e-04,  2.196e-04,  2.131e-04,
        2.015e-04,  1.988e-04,  1.847e-04,  1.729e-04,  1.597e-04,
        1.373e-04,  1.262e-04,  1.087e-04,  9.439e-05,  8.061e-05,
        7.093e-05,  6.049e-05,  5.120e-05,  4.435e-05,  3.817e-05,
        3.340e-05,  2.927e-05,  2.573e-05,  2.291e-05,  2.040e-05,
        1.827e-05,  1.636e-05,  1.463e-05,  1.309e-05,  1.170e-05,
        1.047e-05,  9.315e-06,  8.328e-06,  7.458e-06,  6.665e-06,
        5.940e-06,  5.316e-06,  4.752e-06,  4.252e-06,  3.825e-06,
        3.421e-06,  3.064e-06,  2.746e-06,  2.465e-06,  2.216e-06,
        1.990e-06,  1.790e-06,  1.609e-06,  1.449e-06,  1.306e-06,
        1.177e-06,  1.063e-06,  9.607e-07,  8.672e-07,  7.855e-07,
        7.118e-07,  6.460e-07,  5.871e-07,  5.340e-07,  4.868e-07,
        4.447e-07,  4.068e-07,  3.729e-07,  3.423e-07,  3.151e-07,
        2.905e-07,  2.686e-07,  2.484e-07,  2.306e-07,  2.142e-07,
        1.995e-07,  1.860e-07,  1.738e-07,  1.626e-07,  1.522e-07,
        1.427e-07,  1.338e-07,  1.258e-07,  1.183e-07,  1.116e-07,
        1.056e-07,  9.972e-08,  9.460e-08,  9.007e-08,  8.592e-08,
        8.195e-08,  7.816e-08,  7.483e-08,  7.193e-08,  6.892e-08,
        6.642e-08,  6.386e-08,  6.154e-08,  5.949e-08,  5.764e-08,
        5.622e-08,  5.479e-08,  5.364e-08,  5.301e-08,  5.267e-08,
        5.263e-08,  5.313e-08,  5.410e-08,  5.550e-08,  5.745e-08,
        6.003e-08,  6.311e-08,  6.713e-08,  7.173e-08,  7.724e-08,
        8.368e-08,  9.121e-08,  9.986e-08,  1.097e-07,  1.209e-07,
        1.338e-07,  1.486e-07,  1.651e-07,  1.837e-07,  2.048e-07,
        2.289e-07,  2.557e-07,  2.857e-07,  3.195e-07,  3.587e-07,
        4.015e-07,  4.497e-07,  5.049e-07,  5.665e-07,  6.366e-07,
        7.121e-07,  7.996e-07,  8.946e-07,  1.002e-06,  1.117e-06,
        1.262e-06,  1.416e-06,  1.611e-06,  1.807e-06,  2.056e-06,
        2.351e-06,  2.769e-06,  3.138e-06,  3.699e-06,  4.386e-06,
        5.041e-06,  6.074e-06,  6.812e-06,  7.790e-06,  8.855e-06,
        1.014e-05,  1.095e-05,  1.245e-05,  1.316e-05,  1.390e-05,
        1.504e-05,  1.583e-05,  1.617e-05,  1.652e-05,  1.713e-05,
        1.724e-05,  1.715e-05,  1.668e-05,  1.629e-05,  1.552e-05,
        1.478e-05,  1.340e-05,  1.245e-05,  1.121e-05,  9.575e-06,
        8.956e-06,  7.345e-06,  6.597e-06,  5.612e-06,  4.818e-06,
        4.165e-06,  3.579e-06,  3.041e-06,  2.623e-06,  2.290e-06,
        1.984e-06,  1.748e-06,  1.534e-06,  1.369e-06,  1.219e-06,
        1.092e-06,  9.800e-07,  8.762e-07,  7.896e-07,  7.104e-07,
        6.364e-07,  5.691e-07,  5.107e-07,  4.575e-07,  4.090e-07,
        3.667e-07,  3.287e-07,  2.931e-07,  2.633e-07,  2.356e-07,
        2.111e-07,  1.895e-07,  1.697e-07,  1.525e-07,  1.369e-07,
        1.233e-07,  1.114e-07,  9.988e-08,  9.004e-08,  8.149e-08,
        7.352e-08,  6.662e-08,  6.030e-08,  5.479e-08,  4.974e-08,
        4.532e-08,  4.129e-08,  3.781e-08,  3.462e-08,  3.176e-08,
        2.919e-08,  2.687e-08,  2.481e-08,  2.292e-08,  2.119e-08,
        1.967e-08,  1.828e-08,  1.706e-08,  1.589e-08,  1.487e-08,
        1.393e-08,  1.307e-08,  1.228e-08,  1.156e-08,  1.089e-08,
        1.028e-08,  9.696e-09,  9.159e-09,  8.658e-09,  8.187e-09,
        7.746e-09,  7.340e-09,  6.953e-09,  6.594e-09,  6.259e-09,
        5.948e-09,  5.660e-09,  5.386e-09,  5.135e-09,  4.903e-09,
        4.703e-09,  4.515e-09,  4.362e-09,  4.233e-09,  4.117e-09,
        4.017e-09,  3.962e-09,  3.924e-09,  3.905e-09,  3.922e-09,
        3.967e-09,  4.046e-09,  4.165e-09,  4.320e-09,  4.522e-09,
        4.769e-09,  5.083e-09,  5.443e-09,  5.872e-09,  6.366e-09,
        6.949e-09,  7.601e-09,  8.371e-09,  9.220e-09,  1.020e-08,
        1.129e-08,  1.251e-08,  1.393e-08,  1.542e-08,  1.720e-08,
        1.926e-08,  2.152e-08,  2.392e-08,  2.678e-08,  3.028e-08,
        3.390e-08,  3.836e-08,  4.309e-08,  4.900e-08,  5.481e-08,
        6.252e-08,  7.039e-08,  7.883e-08,  8.849e-08,  1.012e-07,
        1.142e-07,  1.300e-07,  1.475e-07,  1.732e-07,  1.978e-07,
        2.304e-07,  2.631e-07,  2.988e-07,  3.392e-07,  3.690e-07,
        4.355e-07,  4.672e-07,  5.110e-07,  5.461e-07,  5.828e-07,
        6.233e-07,  6.509e-07,  6.672e-07,  6.969e-07,  7.104e-07,
        7.439e-07,  7.463e-07,  7.708e-07,  7.466e-07,  7.668e-07,
        7.549e-07,  7.586e-07,  7.384e-07,  7.439e-07,  7.785e-07,
        7.915e-07,  8.310e-07,  8.745e-07,  9.558e-07,  1.038e-06,
        1.173e-06,  1.304e-06,  1.452e-06,  1.671e-06,  1.931e-06,
        2.239e-06,  2.578e-06,  3.032e-06,  3.334e-06,  3.980e-06,
        4.300e-06,  4.518e-06,  5.321e-06,  5.508e-06,  6.211e-06,
        6.590e-06,  7.046e-06,  7.555e-06,  7.558e-06,  7.875e-06,
        8.319e-06,  8.433e-06,  8.590e-06,  8.503e-06,  8.304e-06,
        8.336e-06,  7.739e-06,  7.301e-06,  6.827e-06,  6.078e-06,
        5.551e-06,  4.762e-06,  4.224e-06,  3.538e-06,  2.984e-06,
        2.619e-06,  2.227e-06,  1.923e-06,  1.669e-06,  1.462e-06,
        1.294e-06,  1.155e-06,  1.033e-06,  9.231e-07,  8.238e-07,
        7.360e-07,  6.564e-07,  5.869e-07,  5.236e-07,  4.673e-07,
        4.174e-07,  3.736e-07,  3.330e-07,  2.976e-07,  2.657e-07,
        2.367e-07,  2.106e-07,  1.877e-07,  1.671e-07,  1.494e-07,
        1.332e-07,  1.192e-07,  1.065e-07,  9.558e-08,  8.586e-08,
        7.717e-08,  6.958e-08,  6.278e-08,  5.666e-08,  5.121e-08,
        4.647e-08,  4.213e-08,  3.815e-08,  3.459e-08,  3.146e-08,
        2.862e-08,  2.604e-08,  2.375e-08,  2.162e-08,  1.981e-08,
        1.817e-08,  1.670e-08,  1.537e-08,  1.417e-08,  1.310e-08,
        1.215e-08,  1.128e-08,  1.050e-08,  9.793e-09,  9.158e-09,
        8.586e-09,  8.068e-09,  7.595e-09,  7.166e-09,  6.778e-09,
        6.427e-09,  6.108e-09,  5.826e-09,  5.571e-09,  5.347e-09,
        5.144e-09,  4.968e-09,  4.822e-09,  4.692e-09,  4.589e-09,
        4.506e-09,  4.467e-09,  4.440e-09,  4.466e-09,  4.515e-09,
        4.718e-09,  4.729e-09,  4.937e-09,  5.249e-09,  5.466e-09,
        5.713e-09,  6.030e-09,  6.436e-09,  6.741e-09,  7.330e-09,
        7.787e-09,  8.414e-09,  8.908e-09,  9.868e-09,  1.069e-08,
        1.158e-08,  1.253e-08,  1.300e-08,  1.409e-08,  1.470e-08,
        1.548e-08,  1.612e-08,  1.666e-08,  1.736e-08,  1.763e-08,
        1.812e-08,  1.852e-08,  1.923e-08,  1.897e-08,  1.893e-08,
        1.888e-08,  1.868e-08,  1.895e-08,  1.899e-08,  1.876e-08,
        1.960e-08,  2.020e-08,  2.121e-08,  2.239e-08,  2.379e-08,
        2.526e-08,  2.766e-08,  2.994e-08,  3.332e-08,  3.703e-08,
        4.158e-08,  4.774e-08,  5.499e-08,  6.355e-08,  7.349e-08,
        8.414e-08,  9.846e-08,  1.143e-07,  1.307e-07,  1.562e-07,
        1.817e-07,  2.011e-07,  2.192e-07,  2.485e-07,  2.867e-07,
        3.035e-07,  3.223e-07,  3.443e-07,  3.617e-07,  3.793e-07,
        3.793e-07,  3.839e-07,  4.081e-07,  4.117e-07,  4.085e-07,
        3.920e-07,  3.851e-07,  3.754e-07,  3.490e-07,  3.229e-07,
        2.978e-07,  2.691e-07,  2.312e-07,  2.029e-07,  1.721e-07,
        1.472e-07,  1.308e-07,  1.132e-07,  9.736e-08,  8.458e-08,
        7.402e-08,  6.534e-08,  5.811e-08,  5.235e-08,  4.762e-08,
        4.293e-08,  3.896e-08,  3.526e-08,  3.165e-08,  2.833e-08,
        2.551e-08,  2.288e-08,  2.036e-08,  1.820e-08,  1.626e-08,
        1.438e-08,  1.299e-08,  1.149e-08,  1.030e-08,  9.148e-09,
        8.122e-09,  7.264e-09,  6.425e-09,  5.777e-09,  5.060e-09,
        4.502e-09,  4.013e-09,  3.567e-09,  3.145e-09,  2.864e-09,
        2.553e-09,  2.311e-09,  2.087e-09,  1.886e-09,  1.716e-09,
        1.556e-09,  1.432e-09,  1.311e-09,  1.202e-09,  1.104e-09,
        1.013e-09,  9.293e-10,  8.493e-10,  7.790e-10,  7.185e-10,
        6.642e-10,  6.141e-10,  5.684e-10,  5.346e-10,  5.032e-10,
        4.725e-10,  4.439e-10,  4.176e-10,  3.930e-10,  3.714e-10,
        3.515e-10,  3.332e-10,  3.167e-10,  3.020e-10,  2.887e-10,
        2.769e-10,  2.665e-10,  2.578e-10,  2.503e-10,  2.436e-10,
        2.377e-10,  2.342e-10,  2.305e-10,  2.296e-10,  2.278e-10,
        2.321e-10,  2.355e-10,  2.402e-10,  2.478e-10,  2.670e-10,
        2.848e-10,  2.982e-10,  3.263e-10,  3.438e-10,  3.649e-10,
        3.829e-10,  4.115e-10,  4.264e-10,  4.473e-10,  4.630e-10,
        4.808e-10,  4.995e-10,  5.142e-10,  5.313e-10,  5.318e-10,
        5.358e-10,  5.452e-10,  5.507e-10,  5.698e-10,  5.782e-10,
        5.983e-10,  6.164e-10,  6.532e-10,  6.811e-10,  7.624e-10,
        8.302e-10,  9.067e-10,  9.937e-10,  1.104e-09,  1.221e-09,
        1.361e-09,  1.516e-09,  1.675e-09,  1.883e-09,  2.101e-09,
        2.349e-09,  2.614e-09,  2.920e-09,  3.305e-09,  3.724e-09,
        4.142e-09,  4.887e-09,  5.614e-09,  6.506e-09,  7.463e-09,
        8.817e-09,  9.849e-09,  1.187e-08,  1.321e-08,  1.474e-08,
        1.698e-08,  1.794e-08,  2.090e-08,  2.211e-08,  2.362e-08,
        2.556e-08,  2.729e-08,  2.880e-08,  3.046e-08,  3.167e-08,
        3.367e-08,  3.457e-08,  3.590e-08,  3.711e-08,  3.826e-08,
        4.001e-08,  4.211e-08,  4.315e-08,  4.661e-08,  5.010e-08,
        5.249e-08,  5.840e-08,  6.628e-08,  7.512e-08,  8.253e-08,
        9.722e-08,  1.067e-07,  1.153e-07,  1.347e-07,  1.428e-07,
        1.577e-07,  1.694e-07,  1.833e-07,  1.938e-07,  2.108e-07,
        2.059e-07,  2.157e-07,  2.185e-07,  2.208e-07,  2.182e-07,
        2.093e-07,  2.014e-07,  1.962e-07,  1.819e-07,  1.713e-07,
        1.510e-07,  1.340e-07,  1.154e-07,  9.890e-08,  8.880e-08,
        7.673e-08,  6.599e-08,  5.730e-08,  5.081e-08,  4.567e-08,
        4.147e-08,  3.773e-08,  3.460e-08,  3.194e-08,  2.953e-08,
        2.759e-08,  2.594e-08,  2.442e-08,  2.355e-08,  2.283e-08,
        2.279e-08,  2.231e-08,  2.279e-08,  2.239e-08,  2.210e-08,
        2.309e-08,  2.293e-08,  2.352e-08,  2.415e-08,  2.430e-08,
        2.426e-08,  2.465e-08,  2.500e-08,  2.496e-08,  2.465e-08,
        2.445e-08,  2.383e-08,  2.299e-08,  2.165e-08,  2.113e-08,
        1.968e-08,  1.819e-08,  1.644e-08,  1.427e-08,  1.270e-08,
        1.082e-08,  9.428e-09,  8.091e-09,  6.958e-09,  5.988e-09,
        5.246e-09,  4.601e-09,  4.098e-09,  3.664e-09,  3.287e-09,
        2.942e-09,  2.656e-09,  2.364e-09,  2.118e-09,  1.903e-09,
        1.703e-09,  1.525e-09,  1.365e-09,  1.229e-09,  1.107e-09,
        9.960e-10,  8.945e-10,  8.080e-10,  7.308e-10,  6.616e-10,
        5.994e-10,  5.422e-10,  4.929e-10,  4.478e-10,  4.070e-10,
        3.707e-10,  3.379e-10,  3.087e-10,  2.823e-10,  2.592e-10,
        2.385e-10,  2.201e-10,  2.038e-10,  1.897e-10,  1.774e-10,
        1.667e-10,  1.577e-10,  1.502e-10,  1.437e-10,  1.394e-10,
        1.358e-10,  1.324e-10,  1.329e-10,  1.324e-10,  1.360e-10,
        1.390e-10,  1.424e-10,  1.544e-10,  1.651e-10,  1.817e-10,
        1.984e-10,  2.195e-10,  2.438e-10,  2.700e-10,  2.991e-10,
        3.322e-10,  3.632e-10,  3.957e-10,  4.360e-10,  4.701e-10,
        5.030e-10,  5.381e-10,  5.793e-10,  6.190e-10,  6.596e-10,
        7.004e-10,  7.561e-10,  7.934e-10,  8.552e-10,  9.142e-10,
        9.570e-10,  1.027e-09,  1.097e-09,  1.193e-09,  1.334e-09,
        1.470e-09,  1.636e-09,  1.871e-09,  2.122e-09,  2.519e-09,
        2.806e-09,  3.203e-09,  3.846e-09,  4.362e-09,  5.114e-09,
        5.643e-09,  6.305e-09,  6.981e-09,  7.983e-09,  8.783e-09,
        9.419e-09,  1.017e-08,  1.063e-08,  1.121e-08,  1.130e-08,
        1.201e-08,  1.225e-08,  1.232e-08,  1.223e-08,  1.177e-08,
        1.151e-08,  1.116e-08,  1.047e-08,  9.698e-09,  8.734e-09,
        8.202e-09,  7.041e-09,  6.074e-09,  5.172e-09,  4.468e-09,
        3.913e-09,  3.414e-09,  2.975e-09,  2.650e-09,  2.406e-09,
        2.173e-09,  2.009e-09,  1.861e-09,  1.727e-09,  1.612e-09,
        1.514e-09,  1.430e-09,  1.362e-09,  1.333e-09,  1.288e-09,
        1.249e-09,  1.238e-09,  1.228e-09,  1.217e-09,  1.202e-09,
        1.209e-09,  1.177e-09,  1.157e-09,  1.165e-09,  1.142e-09,
        1.131e-09,  1.138e-09,  1.117e-09,  1.100e-09,  1.069e-09,
        1.023e-09,  1.005e-09,  9.159e-10,  8.863e-10,  7.865e-10,
        7.153e-10,  6.247e-10,  5.846e-10,  5.133e-10,  4.360e-10,
        3.789e-10,  3.335e-10,  2.833e-10,  2.483e-10,  2.155e-10,
        1.918e-10,  1.709e-10,  1.529e-10,  1.374e-10,  1.235e-10,
        1.108e-10,  9.933e-11,  8.932e-11,  8.022e-11,  7.224e-11,
        6.520e-11,  5.896e-11,  5.328e-11,  4.813e-11,  4.365e-11,
        3.961e-11,  3.594e-11,  3.266e-11,  2.967e-11,  2.701e-11,
        2.464e-11,  2.248e-11,  2.054e-11,  1.878e-11,  1.721e-11,
        1.579e-11,  1.453e-11,  1.341e-11,  1.241e-11,  1.154e-11,
        1.078e-11,  1.014e-11,  9.601e-12,  9.167e-12,  8.838e-12,
        8.614e-12,  8.493e-12,  8.481e-12,  8.581e-12,  8.795e-12,
        9.131e-12,  9.601e-12,  1.021e-11,  1.097e-11,  1.191e-11,
        1.303e-11,  1.439e-11,  1.601e-11,  1.778e-11,  1.984e-11,
        2.234e-11,  2.474e-11,  2.766e-11,  3.085e-11,  3.415e-11,
        3.821e-11,  4.261e-11,  4.748e-11,  5.323e-11,  5.935e-11,
        6.619e-11,  7.418e-11,  8.294e-11,  9.260e-11,  1.039e-10,
        1.156e-10,  1.297e-10,  1.460e-10,  1.641e-10,  1.858e-10,
        2.100e-10,  2.383e-10,  2.724e-10,  3.116e-10,  3.538e-10,
        4.173e-10,  4.727e-10,  5.503e-10,  6.337e-10,  7.320e-10,
        8.298e-10,  9.328e-10,  1.059e-09,  1.176e-09,  1.328e-09,
        1.445e-09,  1.593e-09,  1.770e-09,  1.954e-09,  2.175e-09,
        2.405e-09,  2.622e-09,  2.906e-09,  3.294e-09,  3.713e-09,
        3.980e-09,  4.384e-09,  4.987e-09,  5.311e-09,  5.874e-09,
        6.337e-09,  7.027e-09,  7.390e-09,  7.769e-09,  8.374e-09,
        8.605e-09,  9.165e-09,  9.415e-09,  9.511e-09,  9.704e-09,
        9.588e-09,  9.450e-09,  9.086e-09,  8.798e-09,  8.469e-09,
        7.697e-09,  7.168e-09,  6.255e-09,  5.772e-09,  4.970e-09,
        4.271e-09,  3.653e-09,  3.154e-09,  2.742e-09,  2.435e-09,
        2.166e-09,  1.936e-09,  1.731e-09,  1.556e-09,  1.399e-09,
        1.272e-09,  1.157e-09,  1.066e-09,  9.844e-10,  9.258e-10,
        8.787e-10,  8.421e-10,  8.083e-10,  8.046e-10,  8.067e-10,
        8.181e-10,  8.325e-10,  8.517e-10,  9.151e-10,  9.351e-10,
        9.677e-10,  1.071e-09,  1.126e-09,  1.219e-09,  1.297e-09,
        1.408e-09,  1.476e-09,  1.517e-09,  1.600e-09,  1.649e-09,
        1.678e-09,  1.746e-09,  1.742e-09,  1.728e-09,  1.699e-09,
        1.655e-09,  1.561e-09,  1.480e-09,  1.451e-09,  1.411e-09,
        1.171e-09,  1.106e-09,  9.714e-10,  8.523e-10,  7.346e-10,
        6.241e-10,  5.371e-10,  4.704e-10,  4.144e-10,  3.683e-10,
        3.292e-10,  2.942e-10,  2.620e-10,  2.341e-10,  2.104e-10,
        1.884e-10,  1.700e-10,  1.546e-10,  1.394e-10,  1.265e-10,
        1.140e-10,  1.019e-10,  9.279e-11,  8.283e-11,  7.458e-11,
        6.668e-11,  5.976e-11,  5.330e-11,  4.794e-11,  4.289e-11,
        3.841e-11,  3.467e-11,  3.130e-11,  2.832e-11,  2.582e-11,
        2.356e-11,  2.152e-11,  1.970e-11,  1.808e-11,  1.664e-11,
        1.539e-11,  1.434e-11,  1.344e-11,  1.269e-11,  1.209e-11,
        1.162e-11,  1.129e-11,  1.108e-11,  1.099e-11,  1.103e-11,
        1.119e-11,  1.148e-11,  1.193e-11,  1.252e-11,  1.329e-11,
        1.421e-11,  1.555e-11,  1.685e-11,  1.839e-11,  2.054e-11,
        2.317e-11,  2.571e-11,  2.839e-11,  3.171e-11,  3.490e-11,
        3.886e-11,  4.287e-11,  4.645e-11,  5.047e-11,  5.592e-11,
        6.109e-11,  6.628e-11,  7.381e-11,  8.088e-11,  8.966e-11,
        1.045e-10,  1.120e-10,  1.287e-10,  1.486e-10,  1.662e-10,
        1.866e-10,  2.133e-10,  2.524e-10,  2.776e-10,  3.204e-10,
        3.559e-10,  4.028e-10,  4.448e-10,  4.882e-10,  5.244e-10,
        5.605e-10,  6.018e-10,  6.328e-10,  6.579e-10,  6.541e-10,
        7.024e-10,  7.074e-10,  7.068e-10,  7.009e-10,  6.698e-10,
        6.545e-10,  6.209e-10,  5.834e-10,  5.412e-10,  5.001e-10,
        4.231e-10,  3.727e-10,  3.211e-10,  2.833e-10,  2.447e-10,
        2.097e-10,  1.843e-10,  1.639e-10,  1.449e-10,  1.270e-10,
        1.161e-10,  1.033e-10,  9.282e-11,  8.407e-11,  7.639e-11,
        7.023e-11,  6.474e-11,  6.142e-11,  5.760e-11,  5.568e-11,
        5.472e-11,  5.390e-11,  5.455e-11,  5.540e-11,  5.587e-11,
        6.230e-11,  6.490e-11,  6.868e-11,  7.382e-11,  8.022e-11,
        8.372e-11,  9.243e-11,  1.004e-10,  1.062e-10,  1.130e-10,
        1.176e-10,  1.244e-10,  1.279e-10,  1.298e-10,  1.302e-10,
        1.312e-10,  1.295e-10,  1.244e-10,  1.211e-10,  1.167e-10,
        1.098e-10,  9.927e-11,  8.854e-11,  8.011e-11,  7.182e-11,
        5.923e-11,  5.212e-11,  4.453e-11,  3.832e-11,  3.371e-11,
        2.987e-11,  2.651e-11,  2.354e-11,  2.093e-11,  1.863e-11,
        1.662e-11,  1.486e-11,  1.331e-11,  1.193e-11,  1.071e-11,
        9.628e-12,  8.660e-12,  7.801e-12,  7.031e-12,  6.347e-12,
        5.733e-12,  5.182e-12,  4.695e-12,  4.260e-12,  3.874e-12,
        3.533e-12,  3.235e-12,  2.979e-12,  2.760e-12,  2.579e-12,
        2.432e-12,  2.321e-12,  2.246e-12,  2.205e-12,  2.196e-12,
        2.223e-12,  2.288e-12,  2.387e-12,  2.525e-12,  2.704e-12,
        2.925e-12,  3.191e-12,  3.508e-12,  3.876e-12,  4.303e-12,
        4.793e-12,  5.347e-12,  5.978e-12,  6.682e-12,  7.467e-12,
        8.340e-12,  9.293e-12,  1.035e-11,  1.152e-11,  1.285e-11,
        1.428e-11,  1.586e-11,  1.764e-11,  1.972e-11,  2.214e-11,
        2.478e-11,  2.776e-11,  3.151e-11,  3.591e-11,  4.103e-11,
        4.660e-11,  5.395e-11,  6.306e-11,  7.172e-11,  8.358e-11,
        9.670e-11,  1.110e-10,  1.325e-10,  1.494e-10,  1.736e-10,
        2.007e-10,  2.296e-10,  2.608e-10,  3.004e-10,  3.361e-10,
        3.727e-10,  4.373e-10,  4.838e-10,  5.483e-10,  6.006e-10,
        6.535e-10,  6.899e-10,  7.687e-10,  8.444e-10,  8.798e-10,
        9.135e-10,  9.532e-10,  9.757e-10,  9.968e-10,  1.006e-09,
        9.949e-10,  9.789e-10,  9.564e-10,  9.215e-10,  8.510e-10,
        8.394e-10,  7.707e-10,  7.152e-10,  6.274e-10,  5.598e-10,
        5.028e-10,  4.300e-10,  3.710e-10,  3.245e-10,  2.809e-10,
        2.461e-10,  2.154e-10,  1.910e-10,  1.685e-10,  1.487e-10,
        1.313e-10,  1.163e-10,  1.031e-10,  9.172e-11,  8.221e-11,
        7.382e-11,  6.693e-11,  6.079e-11,  5.581e-11,  5.167e-11,
        4.811e-11,  4.506e-11,  4.255e-11,  4.083e-11,  3.949e-11,
        3.881e-11,  3.861e-11,  3.858e-11,  3.951e-11,  4.045e-11,
        4.240e-11,  4.487e-11,  4.806e-11,  5.133e-11,  5.518e-11,
        5.919e-11,  6.533e-11,  7.031e-11,  7.762e-11,  8.305e-11,
        9.252e-11,  9.727e-11,  1.045e-10,  1.117e-10,  1.200e-10,
        1.275e-10,  1.341e-10,  1.362e-10,  1.438e-10,  1.450e-10,
        1.455e-10,  1.455e-10,  1.434e-10,  1.381e-10,  1.301e-10,
        1.276e-10,  1.163e-10,  1.089e-10,  9.911e-11,  8.943e-11,
        7.618e-11,  6.424e-11,  5.717e-11,  4.866e-11,  4.257e-11,
        3.773e-11,  3.331e-11,  2.958e-11,  2.629e-11,  2.316e-11,
        2.073e-11,  1.841e-11,  1.635e-11,  1.464e-11,  1.310e-11,
        1.160e-11,  1.047e-11,  9.408e-12,  8.414e-12,  7.521e-12,
        6.705e-12,  5.993e-12,  5.371e-12,  4.815e-12,  4.338e-12,
        3.921e-12,  3.567e-12,  3.265e-12,  3.010e-12,  2.795e-12,
        2.613e-12,  2.464e-12,  2.346e-12,  2.256e-12,  2.195e-12,
        2.165e-12,  2.166e-12,  2.198e-12,  2.262e-12,  2.364e-12,
        2.502e-12,  2.682e-12,  2.908e-12,  3.187e-12,  3.533e-12,
        3.946e-12,  4.418e-12,  5.013e-12,  5.708e-12,  6.379e-12,
        7.430e-12,  8.390e-12,  9.510e-12,  1.078e-11,  1.259e-11,
        1.438e-11,  1.630e-11,  1.814e-11,  2.055e-11,  2.348e-11,
        2.664e-11,  2.956e-11,  3.300e-11,  3.677e-11,  4.032e-11,
        4.494e-11,  4.951e-11,  5.452e-11,  6.014e-11,  6.500e-11,
        6.915e-11,  7.450e-11,  7.971e-11,  8.468e-11,  8.726e-11,
        8.995e-11,  9.182e-11,  9.509e-11,  9.333e-11,  9.386e-11,
        9.457e-11,  9.210e-11,  9.019e-11,  8.680e-11,  8.298e-11,
        7.947e-11,  7.460e-11,  7.082e-11,  6.132e-11,  5.855e-11,
        5.073e-11,  4.464e-11,  3.825e-11,  3.375e-11,  2.911e-11,
        2.535e-11,  2.160e-11,  1.907e-11,  1.665e-11,  1.463e-11,
        1.291e-11,  1.133e-11,  9.997e-12,  8.836e-12,  7.839e-12,
        6.943e-12,  6.254e-12,  5.600e-12,  5.029e-12,  4.529e-12,
        4.102e-12,  3.737e-12,  3.428e-12,  3.169e-12,  2.959e-12,
        2.798e-12,  2.675e-12,  2.582e-12,  2.644e-12,  2.557e-12,
        2.614e-12,  2.717e-12,  2.874e-12,  3.056e-12,  3.187e-12,
        3.631e-12,  3.979e-12,  4.248e-12,  4.817e-12,  5.266e-12,
        5.836e-12,  6.365e-12,  6.807e-12,  7.470e-12,  7.951e-12,
        8.636e-12,  8.972e-12,  9.314e-12,  9.445e-12,  1.003e-11,
        1.013e-11,  9.937e-12,  9.729e-12,  9.064e-12,  9.119e-12,
        9.124e-12,  8.704e-12,  8.078e-12,  7.470e-12,  6.329e-12,
        5.674e-12,  4.808e-12,  4.119e-12,  3.554e-12,  3.103e-12,
        2.731e-12,  2.415e-12,  2.150e-12,  1.926e-12,  1.737e-12,
        1.578e-12,  1.447e-12,  1.340e-12,  1.255e-12,  1.191e-12,
        1.146e-12,  1.121e-12,  1.114e-12,  1.126e-12,  1.156e-12,
        1.207e-12,  1.278e-12,  1.372e-12,  1.490e-12,  1.633e-12,
        1.805e-12,  2.010e-12,  2.249e-12,  2.528e-12,  2.852e-12,
        3.228e-12,  3.658e-12,  4.153e-12,  4.728e-12,  5.394e-12,
        6.176e-12,  7.126e-12,  8.188e-12,  9.328e-12,  1.103e-11,
        1.276e-11,  1.417e-11,  1.615e-11,  1.840e-11,  2.155e-11,
        2.429e-11,  2.826e-11,  3.222e-11,  3.664e-11,  4.140e-11,
        4.906e-11,  5.536e-11,  6.327e-11,  7.088e-11,  8.316e-11,
        9.242e-11,  1.070e-10,  1.223e-10,  1.341e-10,  1.553e-10,
        1.703e-10,  1.900e-10,  2.022e-10,  2.233e-10,  2.345e-10,
        2.438e-10,  2.546e-10,  2.599e-10,  2.661e-10,  2.703e-10,
        2.686e-10,  2.662e-10,  2.560e-10,  2.552e-10,  2.378e-10,
        2.252e-10,  2.146e-10,  1.885e-10,  1.668e-10,  1.441e-10,
        1.295e-10,  1.119e-10,  9.893e-11,  8.687e-11,  7.678e-11,
        6.685e-11,  5.879e-11,  5.127e-11,  4.505e-11,  3.997e-11,
        3.511e-11};


// H2O foreign continuum parameters at all temperatures
// date of last update: 11/18/02
// units of (CM**3/MOL)*1.E-20
const Numeric FH2O_ckd_mt_100_v1  =   -20.0;
const Numeric FH2O_ckd_mt_100_v2  = 20000.0;
const Numeric FH2O_ckd_mt_100_dv  =    10.0;
const int     FH2O_ckd_mt_100_npt =  2003;
const double  FH2O_ckd_mt_100[FH2O_ckd_mt_100_npt+addF77fields] = {
        0.000e0,    1.205e-02,  1.126e-02,
        1.095e-02,  1.126e-02,  1.205e-02,  1.322e-02,  1.430e-02,
        1.506e-02,  1.548e-02,  1.534e-02,  1.486e-02,  1.373e-02,
        1.262e-02,  1.134e-02,  1.001e-02,  8.702e-03,  7.475e-03,
        6.481e-03,  5.480e-03,  4.600e-03,  3.833e-03,  3.110e-03,
        2.543e-03,  2.049e-03,  1.680e-03,  1.374e-03,  1.046e-03,
        8.193e-04,  6.267e-04,  4.968e-04,  3.924e-04,  2.983e-04,
        2.477e-04,  1.997e-04,  1.596e-04,  1.331e-04,  1.061e-04,
        8.942e-05,  7.168e-05,  5.887e-05,  4.848e-05,  3.817e-05,
        3.170e-05,  2.579e-05,  2.162e-05,  1.768e-05,  1.490e-05,
        1.231e-05,  1.013e-05,  8.555e-06,  7.328e-06,  6.148e-06,
        5.207e-06,  4.387e-06,  3.741e-06,  3.220e-06,  2.753e-06,
        2.346e-06,  1.985e-06,  1.716e-06,  1.475e-06,  1.286e-06,
        1.122e-06,  9.661e-07,  8.284e-07,  7.057e-07,  6.119e-07,
        5.290e-07,  4.571e-07,  3.948e-07,  3.432e-07,  2.983e-07,
        2.589e-07,  2.265e-07,  1.976e-07,  1.704e-07,  1.456e-07,
        1.260e-07,  1.101e-07,  9.648e-08,  8.415e-08,  7.340e-08,
        6.441e-08,  5.643e-08,  4.940e-08,  4.276e-08,  3.703e-08,
        3.227e-08,  2.825e-08,  2.478e-08,  2.174e-08,  1.898e-08,
        1.664e-08,  1.458e-08,  1.278e-08,  1.126e-08,  9.891e-09,
        8.709e-09,  7.652e-09,  6.759e-09,  5.975e-09,  5.310e-09,
        4.728e-09,  4.214e-09,  3.792e-09,  3.463e-09,  3.226e-09,
        2.992e-09,  2.813e-09,  2.749e-09,  2.809e-09,  2.913e-09,
        3.037e-09,  3.413e-09,  3.738e-09,  4.189e-09,  4.808e-09,
        5.978e-09,  7.088e-09,  8.071e-09,  9.610e-09,  1.210e-08,
        1.500e-08,  1.764e-08,  2.221e-08,  2.898e-08,  3.948e-08,
        5.068e-08,  6.227e-08,  7.898e-08,  1.033e-07,  1.437e-07,
        1.889e-07,  2.589e-07,  3.590e-07,  4.971e-07,  7.156e-07,
        9.983e-07,  1.381e-06,  1.929e-06,  2.591e-06,  3.453e-06,
        4.570e-06,  5.930e-06,  7.552e-06,  9.556e-06,  1.183e-05,
        1.425e-05,  1.681e-05,  1.978e-05,  2.335e-05,  2.668e-05,
        3.022e-05,  3.371e-05,  3.715e-05,  3.967e-05,  4.060e-05,
        4.010e-05,  3.809e-05,  3.491e-05,  3.155e-05,  2.848e-05,
        2.678e-05,  2.660e-05,  2.811e-05,  3.071e-05,  3.294e-05,
        3.459e-05,  3.569e-05,  3.560e-05,  3.434e-05,  3.186e-05,
        2.916e-05,  2.622e-05,  2.275e-05,  1.918e-05,  1.620e-05,
        1.373e-05,  1.182e-05,  1.006e-05,  8.556e-06,  7.260e-06,
        6.107e-06,  5.034e-06,  4.211e-06,  3.426e-06,  2.865e-06,
        2.446e-06,  1.998e-06,  1.628e-06,  1.242e-06,  1.005e-06,
        7.853e-07,  6.210e-07,  5.071e-07,  4.156e-07,  3.548e-07,
        2.825e-07,  2.261e-07,  1.916e-07,  1.510e-07,  1.279e-07,
        1.059e-07,  9.140e-08,  7.707e-08,  6.170e-08,  5.311e-08,
        4.263e-08,  3.518e-08,  2.961e-08,  2.457e-08,  2.119e-08,
        1.712e-08,  1.439e-08,  1.201e-08,  1.003e-08,  8.564e-09,
        7.199e-09,  6.184e-09,  5.206e-09,  4.376e-09,  3.708e-09,
        3.157e-09,  2.725e-09,  2.361e-09,  2.074e-09,  1.797e-09,
        1.562e-09,  1.364e-09,  1.196e-09,  1.042e-09,  8.862e-10,
        7.648e-10,  6.544e-10,  5.609e-10,  4.791e-10,  4.108e-10,
        3.531e-10,  3.038e-10,  2.618e-10,  2.268e-10,  1.969e-10,
        1.715e-10,  1.496e-10,  1.308e-10,  1.147e-10,  1.008e-10,
        8.894e-11,  7.885e-11,  7.031e-11,  6.355e-11,  5.854e-11,
        5.534e-11,  5.466e-11,  5.725e-11,  6.447e-11,  7.943e-11,
        1.038e-10,  1.437e-10,  2.040e-10,  2.901e-10,  4.051e-10,
        5.556e-10,  7.314e-10,  9.291e-10,  1.134e-09,  1.321e-09,
        1.482e-09,  1.596e-09,  1.669e-09,  1.715e-09,  1.762e-09,
        1.817e-09,  1.828e-09,  1.848e-09,  1.873e-09,  1.902e-09,
        1.894e-09,  1.864e-09,  1.841e-09,  1.797e-09,  1.704e-09,
        1.559e-09,  1.382e-09,  1.187e-09,  1.001e-09,  8.468e-10,
        7.265e-10,  6.521e-10,  6.381e-10,  6.660e-10,  7.637e-10,
        9.705e-10,  1.368e-09,  1.856e-09,  2.656e-09,  3.954e-09,
        5.960e-09,  8.720e-09,  1.247e-08,  1.781e-08,  2.491e-08,
        3.311e-08,  4.272e-08,  5.205e-08,  6.268e-08,  7.337e-08,
        8.277e-08,  9.185e-08,  1.004e-07,  1.091e-07,  1.159e-07,
        1.188e-07,  1.175e-07,  1.124e-07,  1.033e-07,  9.381e-08,
        8.501e-08,  7.956e-08,  7.894e-08,  8.331e-08,  9.102e-08,
        9.836e-08,  1.035e-07,  1.064e-07,  1.060e-07,  1.032e-07,
        9.808e-08,  9.139e-08,  8.442e-08,  7.641e-08,  6.881e-08,
        6.161e-08,  5.404e-08,  4.804e-08,  4.446e-08,  4.328e-08,
        4.259e-08,  4.421e-08,  4.673e-08,  4.985e-08,  5.335e-08,
        5.796e-08,  6.542e-08,  7.714e-08,  8.827e-08,  1.040e-07,
        1.238e-07,  1.499e-07,  1.829e-07,  2.222e-07,  2.689e-07,
        3.303e-07,  3.981e-07,  4.840e-07,  5.910e-07,  7.363e-07,
        9.087e-07,  1.139e-06,  1.455e-06,  1.866e-06,  2.440e-06,
        3.115e-06,  3.941e-06,  4.891e-06,  5.992e-06,  7.111e-06,
        8.296e-06,  9.210e-06,  9.987e-06,  1.044e-05,  1.073e-05,
        1.092e-05,  1.106e-05,  1.138e-05,  1.171e-05,  1.186e-05,
        1.186e-05,  1.179e-05,  1.166e-05,  1.151e-05,  1.160e-05,
        1.197e-05,  1.241e-05,  1.268e-05,  1.260e-05,  1.184e-05,
        1.063e-05,  9.204e-06,  7.584e-06,  6.053e-06,  4.482e-06,
        3.252e-06,  2.337e-06,  1.662e-06,  1.180e-06,  8.150e-07,
        5.950e-07,  4.354e-07,  3.302e-07,  2.494e-07,  1.930e-07,
        1.545e-07,  1.250e-07,  1.039e-07,  8.602e-08,  7.127e-08,
        5.897e-08,  4.838e-08,  4.018e-08,  3.280e-08,  2.720e-08,
        2.307e-08,  1.972e-08,  1.654e-08,  1.421e-08,  1.174e-08,
        1.004e-08,  8.739e-09,  7.358e-09,  6.242e-09,  5.303e-09,
        4.567e-09,  3.940e-09,  3.375e-09,  2.864e-09,  2.422e-09,
        2.057e-09,  1.750e-09,  1.505e-09,  1.294e-09,  1.101e-09,
        9.401e-10,  8.018e-10,  6.903e-10,  5.965e-10,  5.087e-10,
        4.364e-10,  3.759e-10,  3.247e-10,  2.809e-10,  2.438e-10,
        2.123e-10,  1.853e-10,  1.622e-10,  1.426e-10,  1.260e-10,
        1.125e-10,  1.022e-10,  9.582e-11,  9.388e-11,  9.801e-11,
        1.080e-10,  1.276e-10,  1.551e-10,  1.903e-10,  2.291e-10,
        2.724e-10,  3.117e-10,  3.400e-10,  3.562e-10,  3.625e-10,
        3.619e-10,  3.429e-10,  3.221e-10,  2.943e-10,  2.645e-10,
        2.338e-10,  2.062e-10,  1.901e-10,  1.814e-10,  1.827e-10,
        1.906e-10,  1.984e-10,  2.040e-10,  2.068e-10,  2.075e-10,
        2.018e-10,  1.959e-10,  1.897e-10,  1.852e-10,  1.791e-10,
        1.696e-10,  1.634e-10,  1.598e-10,  1.561e-10,  1.518e-10,
        1.443e-10,  1.377e-10,  1.346e-10,  1.342e-10,  1.375e-10,
        1.525e-10,  1.767e-10,  2.108e-10,  2.524e-10,  2.981e-10,
        3.477e-10,  4.262e-10,  5.326e-10,  6.646e-10,  8.321e-10,
        1.069e-09,  1.386e-09,  1.743e-09,  2.216e-09,  2.808e-09,
        3.585e-09,  4.552e-09,  5.907e-09,  7.611e-09,  9.774e-09,
        1.255e-08,  1.666e-08,  2.279e-08,  3.221e-08,  4.531e-08,
        6.400e-08,  9.187e-08,  1.295e-07,  1.825e-07,  2.431e-07,
        3.181e-07,  4.009e-07,  4.941e-07,  5.880e-07,  6.623e-07,
        7.155e-07,  7.451e-07,  7.594e-07,  7.541e-07,  7.467e-07,
        7.527e-07,  7.935e-07,  8.461e-07,  8.954e-07,  9.364e-07,
        9.843e-07,  1.024e-06,  1.050e-06,  1.059e-06,  1.074e-06,
        1.072e-06,  1.043e-06,  9.789e-07,  8.803e-07,  7.662e-07,
        6.378e-07,  5.133e-07,  3.958e-07,  2.914e-07,  2.144e-07,
        1.570e-07,  1.140e-07,  8.470e-08,  6.200e-08,  4.657e-08,
        3.559e-08,  2.813e-08,  2.222e-08,  1.769e-08,  1.391e-08,
        1.125e-08,  9.186e-09,  7.704e-09,  6.447e-09,  5.381e-09,
        4.442e-09,  3.669e-09,  3.057e-09,  2.564e-09,  2.153e-09,
        1.784e-09,  1.499e-09,  1.281e-09,  1.082e-09,  9.304e-10,
        8.169e-10,  6.856e-10,  5.866e-10,  5.043e-10,  4.336e-10,
        3.731e-10,  3.175e-10,  2.745e-10,  2.374e-10,  2.007e-10,
        1.737e-10,  1.508e-10,  1.302e-10,  1.130e-10,  9.672e-11,
        8.375e-11,  7.265e-11,  6.244e-11,  5.343e-11,  4.654e-11,
        3.975e-11,  3.488e-11,  3.097e-11,  2.834e-11,  2.649e-11,
        2.519e-11,  2.462e-11,  2.443e-11,  2.440e-11,  2.398e-11,
        2.306e-11,  2.183e-11,  2.021e-11,  1.821e-11,  1.599e-11,
        1.403e-11,  1.196e-11,  1.023e-11,  8.728e-12,  7.606e-12,
        6.941e-12,  6.545e-12,  6.484e-12,  6.600e-12,  6.718e-12,
        6.785e-12,  6.746e-12,  6.724e-12,  6.764e-12,  6.995e-12,
        7.144e-12,  7.320e-12,  7.330e-12,  7.208e-12,  6.789e-12,
        6.090e-12,  5.337e-12,  4.620e-12,  4.037e-12,  3.574e-12,
        3.311e-12,  3.346e-12,  3.566e-12,  3.836e-12,  4.076e-12,
        4.351e-12,  4.691e-12,  5.114e-12,  5.427e-12,  6.167e-12,
        7.436e-12,  8.842e-12,  1.038e-11,  1.249e-11,  1.540e-11,
        1.915e-11,  2.480e-11,  3.256e-11,  4.339e-11,  5.611e-11,
        7.519e-11,  1.037e-10,  1.409e-10,  1.883e-10,  2.503e-10,
        3.380e-10,  4.468e-10,  5.801e-10,  7.335e-10,  8.980e-10,
        1.110e-09,  1.363e-09,  1.677e-09,  2.104e-09,  2.681e-09,
        3.531e-09,  4.621e-09,  6.106e-09,  8.154e-09,  1.046e-08,
        1.312e-08,  1.607e-08,  1.948e-08,  2.266e-08,  2.495e-08,
        2.655e-08,  2.739e-08,  2.739e-08,  2.662e-08,  2.589e-08,
        2.590e-08,  2.664e-08,  2.833e-08,  3.023e-08,  3.305e-08,
        3.558e-08,  3.793e-08,  3.961e-08,  4.056e-08,  4.102e-08,
        4.025e-08,  3.917e-08,  3.706e-08,  3.493e-08,  3.249e-08,
        3.096e-08,  3.011e-08,  3.111e-08,  3.395e-08,  3.958e-08,
        4.875e-08,  6.066e-08,  7.915e-08,  1.011e-07,  1.300e-07,
        1.622e-07,  2.003e-07,  2.448e-07,  2.863e-07,  3.317e-07,
        3.655e-07,  3.960e-07,  4.098e-07,  4.168e-07,  4.198e-07,
        4.207e-07,  4.289e-07,  4.384e-07,  4.471e-07,  4.524e-07,
        4.574e-07,  4.633e-07,  4.785e-07,  5.028e-07,  5.371e-07,
        5.727e-07,  5.955e-07,  5.998e-07,  5.669e-07,  5.082e-07,
        4.397e-07,  3.596e-07,  2.814e-07,  2.074e-07,  1.486e-07,
        1.057e-07,  7.250e-08,  4.946e-08,  3.430e-08,  2.447e-08,
        1.793e-08,  1.375e-08,  1.096e-08,  9.091e-09,  7.709e-09,
        6.631e-09,  5.714e-09,  4.886e-09,  4.205e-09,  3.575e-09,
        3.070e-09,  2.631e-09,  2.284e-09,  2.002e-09,  1.745e-09,
        1.509e-09,  1.284e-09,  1.084e-09,  9.163e-10,  7.663e-10,
        6.346e-10,  5.283e-10,  4.354e-10,  3.590e-10,  2.982e-10,
        2.455e-10,  2.033e-10,  1.696e-10,  1.432e-10,  1.211e-10,
        1.020e-10,  8.702e-11,  7.380e-11,  6.293e-11,  5.343e-11,
        4.532e-11,  3.907e-11,  3.365e-11,  2.945e-11,  2.558e-11,
        2.192e-11,  1.895e-11,  1.636e-11,  1.420e-11,  1.228e-11,
        1.063e-11,  9.348e-12,  8.200e-12,  7.231e-12,  6.430e-12,
        5.702e-12,  5.052e-12,  4.469e-12,  4.000e-12,  3.679e-12,
        3.387e-12,  3.197e-12,  3.158e-12,  3.327e-12,  3.675e-12,
        4.292e-12,  5.437e-12,  7.197e-12,  1.008e-11,  1.437e-11,
        2.035e-11,  2.905e-11,  4.062e-11,  5.528e-11,  7.177e-11,
        9.064e-11,  1.109e-10,  1.297e-10,  1.473e-10,  1.652e-10,
        1.851e-10,  2.079e-10,  2.313e-10,  2.619e-10,  2.958e-10,
        3.352e-10,  3.796e-10,  4.295e-10,  4.923e-10,  5.490e-10,
        5.998e-10,  6.388e-10,  6.645e-10,  6.712e-10,  6.549e-10,
        6.380e-10,  6.255e-10,  6.253e-10,  6.459e-10,  6.977e-10,
        7.590e-10,  8.242e-10,  8.920e-10,  9.403e-10,  9.701e-10,
        9.483e-10,  9.135e-10,  8.617e-10,  7.921e-10,  7.168e-10,
        6.382e-10,  5.677e-10,  5.045e-10,  4.572e-10,  4.312e-10,
        4.145e-10,  4.192e-10,  4.541e-10,  5.368e-10,  6.771e-10,
        8.962e-10,  1.210e-09,  1.659e-09,  2.330e-09,  3.249e-09,
        4.495e-09,  5.923e-09,  7.642e-09,  9.607e-09,  1.178e-08,
        1.399e-08,  1.584e-08,  1.730e-08,  1.816e-08,  1.870e-08,
        1.868e-08,  1.870e-08,  1.884e-08,  1.990e-08,  2.150e-08,
        2.258e-08,  2.364e-08,  2.473e-08,  2.602e-08,  2.689e-08,
        2.731e-08,  2.816e-08,  2.859e-08,  2.839e-08,  2.703e-08,
        2.451e-08,  2.149e-08,  1.787e-08,  1.449e-08,  1.111e-08,
        8.282e-09,  6.121e-09,  4.494e-09,  3.367e-09,  2.487e-09,
        1.885e-09,  1.503e-09,  1.249e-09,  1.074e-09,  9.427e-10,
        8.439e-10,  7.563e-10,  6.772e-10,  6.002e-10,  5.254e-10,
        4.588e-10,  3.977e-10,  3.449e-10,  3.003e-10,  2.624e-10,
        2.335e-10,  2.040e-10,  1.771e-10,  1.534e-10,  1.296e-10,
        1.097e-10,  9.173e-11,  7.730e-11,  6.547e-11,  5.191e-11,
        4.198e-11,  3.361e-11,  2.732e-11,  2.244e-11,  1.791e-11,
        1.509e-11,  1.243e-11,  1.035e-11,  8.969e-12,  7.394e-12,
        6.323e-12,  5.282e-12,  4.543e-12,  3.752e-12,  3.140e-12,
        2.600e-12,  2.194e-12,  1.825e-12,  1.511e-12,  1.245e-12,
        1.024e-12,  8.539e-13,  7.227e-13,  6.102e-13,  5.189e-13,
        4.430e-13,  3.774e-13,  3.236e-13,  2.800e-13,  2.444e-13,
        2.156e-13,  1.932e-13,  1.775e-13,  1.695e-13,  1.672e-13,
        1.704e-13,  1.825e-13,  2.087e-13,  2.614e-13,  3.377e-13,
        4.817e-13,  6.989e-13,  1.062e-12,  1.562e-12,  2.288e-12,
        3.295e-12,  4.550e-12,  5.965e-12,  7.546e-12,  9.395e-12,
        1.103e-11,  1.228e-11,  1.318e-11,  1.380e-11,  1.421e-11,
        1.390e-11,  1.358e-11,  1.336e-11,  1.342e-11,  1.356e-11,
        1.424e-11,  1.552e-11,  1.730e-11,  1.951e-11,  2.128e-11,
        2.249e-11,  2.277e-11,  2.226e-11,  2.111e-11,  1.922e-11,
        1.775e-11,  1.661e-11,  1.547e-11,  1.446e-11,  1.323e-11,
        1.210e-11,  1.054e-11,  9.283e-12,  8.671e-12,  8.670e-12,
        9.429e-12,  1.062e-11,  1.255e-11,  1.506e-11,  1.818e-11,
        2.260e-11,  2.831e-11,  3.723e-11,  5.092e-11,  6.968e-11,
        9.826e-11,  1.349e-10,  1.870e-10,  2.580e-10,  3.430e-10,
        4.424e-10,  5.521e-10,  6.812e-10,  8.064e-10,  9.109e-10,
        9.839e-10,  1.028e-09,  1.044e-09,  1.029e-09,  1.005e-09,
        1.002e-09,  1.038e-09,  1.122e-09,  1.233e-09,  1.372e-09,
        1.524e-09,  1.665e-09,  1.804e-09,  1.908e-09,  2.015e-09,
        2.117e-09,  2.219e-09,  2.336e-09,  2.531e-09,  2.805e-09,
        3.189e-09,  3.617e-09,  4.208e-09,  4.911e-09,  5.619e-09,
        6.469e-09,  7.188e-09,  7.957e-09,  8.503e-09,  9.028e-09,
        9.571e-09,  9.990e-09,  1.055e-08,  1.102e-08,  1.132e-08,
        1.141e-08,  1.145e-08,  1.145e-08,  1.176e-08,  1.224e-08,
        1.304e-08,  1.388e-08,  1.445e-08,  1.453e-08,  1.368e-08,
        1.220e-08,  1.042e-08,  8.404e-09,  6.403e-09,  4.643e-09,
        3.325e-09,  2.335e-09,  1.638e-09,  1.190e-09,  9.161e-10,
        7.412e-10,  6.226e-10,  5.516e-10,  5.068e-10,  4.831e-10,
        4.856e-10,  5.162e-10,  5.785e-10,  6.539e-10,  7.485e-10,
        8.565e-10,  9.534e-10,  1.052e-09,  1.115e-09,  1.173e-09,
        1.203e-09,  1.224e-09,  1.243e-09,  1.248e-09,  1.261e-09,
        1.265e-09,  1.250e-09,  1.217e-09,  1.176e-09,  1.145e-09,
        1.153e-09,  1.199e-09,  1.278e-09,  1.366e-09,  1.426e-09,
        1.444e-09,  1.365e-09,  1.224e-09,  1.051e-09,  8.539e-10,
        6.564e-10,  4.751e-10,  3.404e-10,  2.377e-10,  1.631e-10,
        1.114e-10,  7.870e-11,  5.793e-11,  4.284e-11,  3.300e-11,
        2.620e-11,  2.152e-11,  1.777e-11,  1.496e-11,  1.242e-11,
        1.037e-11,  8.725e-12,  7.004e-12,  5.718e-12,  4.769e-12,
        3.952e-12,  3.336e-12,  2.712e-12,  2.213e-12,  1.803e-12,
        1.492e-12,  1.236e-12,  1.006e-12,  8.384e-13,  7.063e-13,
        5.879e-13,  4.930e-13,  4.171e-13,  3.569e-13,  3.083e-13,
        2.688e-13,  2.333e-13,  2.035e-13,  1.820e-13,  1.682e-13,
        1.635e-13,  1.674e-13,  1.769e-13,  2.022e-13,  2.485e-13,
        3.127e-13,  4.250e-13,  5.928e-13,  8.514e-13,  1.236e-12,
        1.701e-12,  2.392e-12,  3.231e-12,  4.350e-12,  5.559e-12,
        6.915e-12,  8.519e-12,  1.013e-11,  1.146e-11,  1.240e-11,
        1.305e-11,  1.333e-11,  1.318e-11,  1.263e-11,  1.238e-11,
        1.244e-11,  1.305e-11,  1.432e-11,  1.623e-11,  1.846e-11,
        2.090e-11,  2.328e-11,  2.526e-11,  2.637e-11,  2.702e-11,
        2.794e-11,  2.889e-11,  2.989e-11,  3.231e-11,  3.680e-11,
        4.375e-11,  5.504e-11,  7.159e-11,  9.502e-11,  1.279e-10,
        1.645e-10,  2.098e-10,  2.618e-10,  3.189e-10,  3.790e-10,
        4.303e-10,  4.753e-10,  5.027e-10,  5.221e-10,  5.293e-10,
        5.346e-10,  5.467e-10,  5.796e-10,  6.200e-10,  6.454e-10,
        6.705e-10,  6.925e-10,  7.233e-10,  7.350e-10,  7.538e-10,
        7.861e-10,  8.077e-10,  8.132e-10,  7.749e-10,  7.036e-10,
        6.143e-10,  5.093e-10,  4.089e-10,  3.092e-10,  2.299e-10,
        1.705e-10,  1.277e-10,  9.723e-11,  7.533e-11,  6.126e-11,
        5.154e-11,  4.428e-11,  3.913e-11,  3.521e-11,  3.297e-11,
        3.275e-11,  3.460e-11,  3.798e-11,  4.251e-11,  4.745e-11,
        5.232e-11,  5.606e-11,  5.820e-11,  5.880e-11,  5.790e-11,
        5.661e-11,  5.491e-11,  5.366e-11,  5.341e-11,  5.353e-11,
        5.336e-11,  5.293e-11,  5.248e-11,  5.235e-11,  5.208e-11,
        5.322e-11,  5.521e-11,  5.725e-11,  5.827e-11,  5.685e-11,
        5.245e-11,  4.612e-11,  3.884e-11,  3.129e-11,  2.404e-11,
        1.732e-11,  1.223e-11,  8.574e-12,  5.888e-12,  3.986e-12,
        2.732e-12,  1.948e-12,  1.414e-12,  1.061e-12,  8.298e-13,
        6.612e-13,  5.413e-13,  4.472e-13,  3.772e-13,  3.181e-13,
        2.645e-13,  2.171e-13,  1.778e-13,  1.464e-13,  1.183e-13,
        9.637e-14,  7.991e-14,  6.668e-14,  5.570e-14,  4.663e-14,
        3.848e-14,  3.233e-14,  2.706e-14,  2.284e-14,  1.944e-14,
        1.664e-14,  1.430e-14,  1.233e-14,  1.066e-14,  9.234e-15,
        8.023e-15,  6.993e-15,  6.119e-15,  5.384e-15,  4.774e-15,
        4.283e-15,  3.916e-15,  3.695e-15,  3.682e-15,  4.004e-15,
        4.912e-15,  6.853e-15,  1.056e-14,  1.712e-14,  2.804e-14,
        4.516e-14,  7.113e-14,  1.084e-13,  1.426e-13,  1.734e-13,
        1.978e-13,  2.194e-13,  2.388e-13,  2.489e-13,  2.626e-13,
        2.865e-13,  3.105e-13,  3.387e-13,  3.652e-13,  3.984e-13,
        4.398e-13,  4.906e-13,  5.550e-13,  6.517e-13,  7.813e-13,
        9.272e-13,  1.164e-12,  1.434e-12,  1.849e-12,  2.524e-12,
        3.328e-12,  4.523e-12,  6.108e-12,  8.207e-12,  1.122e-11,
        1.477e-11,  1.900e-11,  2.412e-11,  2.984e-11,  3.680e-11,
        4.353e-11,  4.963e-11,  5.478e-11,  5.903e-11,  6.233e-11,
        6.483e-11,  6.904e-11,  7.569e-11,  8.719e-11,  1.048e-10,
        1.278e-10,  1.557e-10,  1.869e-10,  2.218e-10,  2.610e-10,
        2.975e-10,  3.371e-10,  3.746e-10,  4.065e-10,  4.336e-10,
        4.503e-10,  4.701e-10,  4.800e-10,  4.917e-10,  5.038e-10,
        5.128e-10,  5.143e-10,  5.071e-10,  5.019e-10,  5.025e-10,
        5.183e-10,  5.496e-10,  5.877e-10,  6.235e-10,  6.420e-10,
        6.234e-10,  5.698e-10,  4.916e-10,  4.022e-10,  3.126e-10,
        2.282e-10,  1.639e-10,  1.142e-10,  7.919e-11,  5.690e-11,
        4.313e-11,  3.413e-11,  2.807e-11,  2.410e-11,  2.166e-11,
        2.024e-11,  1.946e-11,  1.929e-11,  1.963e-11,  2.035e-11,
        2.162e-11,  2.305e-11,  2.493e-11,  2.748e-11,  3.048e-11,
        3.413e-11,  3.754e-11,  4.155e-11,  4.635e-11,  5.110e-11,
        5.734e-11,  6.338e-11,  6.990e-11,  7.611e-11,  8.125e-11,
        8.654e-11,  8.951e-11,  9.182e-11,  9.310e-11,  9.273e-11,
        9.094e-11,  8.849e-11,  8.662e-11,  8.670e-11,  8.972e-11,
        9.566e-11,  1.025e-10,  1.083e-10,  1.111e-10,  1.074e-10,
        9.771e-11,  8.468e-11,  6.958e-11,  5.470e-11,  4.040e-11,
        2.940e-11,  2.075e-11,  1.442e-11,  1.010e-11,  7.281e-12,
        5.409e-12,  4.138e-12,  3.304e-12,  2.784e-12,  2.473e-12,
        2.273e-12,  2.186e-12,  2.118e-12,  2.066e-12,  1.958e-12,
        1.818e-12,  1.675e-12,  1.509e-12,  1.349e-12,  1.171e-12,
        9.838e-13,  8.213e-13,  6.765e-13,  5.378e-13,  4.161e-13,
        3.119e-13,  2.279e-13,  1.637e-13,  1.152e-13,  8.112e-14,
        5.919e-14,  4.470e-14,  3.492e-14,  2.811e-14,  2.319e-14,
        1.948e-14,  1.660e-14,  1.432e-14,  1.251e-14,  1.109e-14,
        1.006e-14,  9.450e-15,  9.384e-15,  1.012e-14,  1.216e-14,
        1.636e-14,  2.305e-14,  3.488e-14,  5.572e-14,  8.479e-14,
        1.265e-13,  1.905e-13,  2.730e-13,  3.809e-13,  4.955e-13,
        6.303e-13,  7.861e-13,  9.427e-13,  1.097e-12,  1.212e-12,
        1.328e-12,  1.415e-12,  1.463e-12,  1.495e-12,  1.571e-12,
        1.731e-12,  1.981e-12,  2.387e-12,  2.930e-12,  3.642e-12,
        4.584e-12,  5.822e-12,  7.278e-12,  9.193e-12,  1.135e-11,
        1.382e-11,  1.662e-11,  1.958e-11,  2.286e-11,  2.559e-11,
        2.805e-11,  2.988e-11,  3.106e-11,  3.182e-11,  3.200e-11,
        3.258e-11,  3.362e-11,  3.558e-11,  3.688e-11,  3.800e-11,
        3.929e-11,  4.062e-11,  4.186e-11,  4.293e-11,  4.480e-11,
        4.643e-11,  4.704e-11,  4.571e-11,  4.206e-11,  3.715e-11,
        3.131e-11,  2.541e-11,  1.978e-11,  1.508e-11,  1.146e-11,
        8.700e-12,  6.603e-12,  5.162e-12,  4.157e-12,  3.408e-12,
        2.829e-12,  2.405e-12,  2.071e-12,  1.826e-12,  1.648e-12,
        1.542e-12,  1.489e-12,  1.485e-12,  1.493e-12,  1.545e-12,
        1.637e-12,  1.814e-12,  2.061e-12,  2.312e-12,  2.651e-12,
        3.030e-12,  3.460e-12,  3.901e-12,  4.306e-12,  4.721e-12,
        5.008e-12,  5.281e-12,  5.541e-12,  5.791e-12,  6.115e-12,
        6.442e-12,  6.680e-12,  6.791e-12,  6.831e-12,  6.839e-12,
        6.946e-12,  7.128e-12,  7.537e-12,  8.036e-12,  8.392e-12,
        8.526e-12,  8.110e-12,  7.325e-12,  6.329e-12,  5.183e-12,
        4.081e-12,  2.985e-12,  2.141e-12,  1.492e-12,  1.015e-12,
        6.684e-13,  4.414e-13,  2.987e-13,  2.038e-13,  1.391e-13,
        9.860e-14,  7.240e-14,  5.493e-14,  4.288e-14,  3.427e-14,
        2.787e-14,  2.296e-14,  1.909e-14,  1.598e-14,  1.344e-14,
        1.135e-14,  9.616e-15,  8.169e-15,  6.957e-15,  5.938e-15,
        5.080e-15,  4.353e-15,  3.738e-15,  3.217e-15,  2.773e-15,
        2.397e-15,  2.077e-15,  1.805e-15,  1.575e-15,  1.382e-15,
        1.221e-15,  1.090e-15,  9.855e-16,  9.068e-16,  8.537e-16,
        8.270e-16,  8.290e-16,  8.634e-16,  9.359e-16,  1.055e-15,
        1.233e-15,  1.486e-15,  1.839e-15,  2.326e-15,  2.998e-15,
        3.934e-15,  5.256e-15,  7.164e-15,  9.984e-15,  1.427e-14,
        2.099e-14,  3.196e-14,  5.121e-14,  7.908e-14,  1.131e-13,
        1.602e-13,  2.239e-13,  3.075e-13,  4.134e-13,  5.749e-13,
        7.886e-13,  1.071e-12,  1.464e-12,  2.032e-12,  2.800e-12,
        3.732e-12,  4.996e-12,  6.483e-12,  8.143e-12,  1.006e-11,
        1.238e-11,  1.484e-11,  1.744e-11,  2.020e-11,  2.274e-11,
        2.562e-11,  2.848e-11,  3.191e-11,  3.617e-11,  4.081e-11,
        4.577e-11,  4.937e-11,  5.204e-11,  5.401e-11,  5.462e-11,
        5.507e-11,  5.510e-11,  5.605e-11,  5.686e-11,  5.739e-11,
        5.766e-11,  5.740e-11,  5.754e-11,  5.761e-11,  5.777e-11,
        5.712e-11,  5.510e-11,  5.088e-11,  4.438e-11,  3.728e-11,
        2.994e-11,  2.305e-11,  1.715e-11,  1.256e-11,  9.208e-12,
        6.745e-12,  5.014e-12,  3.785e-12,  2.900e-12,  2.239e-12,
        1.757e-12,  1.414e-12,  1.142e-12,  9.482e-13,  8.010e-13,
        6.961e-13,  6.253e-13,  5.735e-13,  5.433e-13,  5.352e-13,
        5.493e-13,  5.706e-13,  6.068e-13,  6.531e-13,  7.109e-13,
        7.767e-13,  8.590e-13,  9.792e-13,  1.142e-12,  1.371e-12,
        1.650e-12,  1.957e-12,  2.302e-12,  2.705e-12,  3.145e-12,
        3.608e-12,  4.071e-12,  4.602e-12,  5.133e-12,  5.572e-12,
        5.987e-12,  6.248e-12,  6.533e-12,  6.757e-12,  6.935e-12,
        7.224e-12,  7.422e-12,  7.538e-12,  7.547e-12,  7.495e-12,
        7.543e-12,  7.725e-12,  8.139e-12,  8.627e-12,  9.146e-12,
        9.443e-12,  9.318e-12,  8.649e-12,  7.512e-12,  6.261e-12,
        4.915e-12,  3.647e-12,  2.597e-12,  1.785e-12,  1.242e-12,
        8.660e-13,  6.207e-13,  4.610e-13,  3.444e-13,  2.634e-13,
        2.100e-13,  1.725e-13,  1.455e-13,  1.237e-13,  1.085e-13,
        9.513e-14,  7.978e-14,  6.603e-14,  5.288e-14,  4.084e-14,
        2.952e-14,  2.157e-14,  1.593e-14,  1.199e-14,  9.267e-15,
        7.365e-15,  6.004e-15,  4.995e-15,  4.218e-15,  3.601e-15,
        3.101e-15,  2.692e-15,  2.360e-15,  2.094e-15,  1.891e-15,
        1.755e-15,  1.699e-15,  1.755e-15,  1.987e-15,  2.506e-15,
        3.506e-15,  5.289e-15,  8.311e-15,  1.325e-14,  2.129e-14,
        3.237e-14,  4.595e-14,  6.441e-14,  8.433e-14,  1.074e-13,
        1.383e-13,  1.762e-13,  2.281e-13,  2.831e-13,  3.523e-13,
        4.380e-13,  5.304e-13,  6.290e-13,  7.142e-13,  8.032e-13,
        8.934e-13,  9.888e-13,  1.109e-12,  1.261e-12,  1.462e-12,
        1.740e-12,  2.099e-12,  2.535e-12,  3.008e-12,  3.462e-12,
        3.856e-12,  4.098e-12,  4.239e-12,  4.234e-12,  4.132e-12,
        3.986e-12,  3.866e-12,  3.829e-12,  3.742e-12,  3.705e-12,
        3.694e-12,  3.765e-12,  3.849e-12,  3.929e-12,  4.056e-12,
        4.092e-12,  4.047e-12,  3.792e-12,  3.407e-12,  2.953e-12,
        2.429e-12,  1.931e-12,  1.460e-12,  1.099e-12,  8.199e-13,
        6.077e-13,  4.449e-13,  3.359e-13,  2.524e-13,  1.881e-13,
        1.391e-13,  1.020e-13,  7.544e-14,  5.555e-14,  4.220e-14,
        3.321e-14,  2.686e-14,  2.212e-14,  1.780e-14,  1.369e-14,
        1.094e-14,  9.130e-15,  8.101e-15,  7.828e-15,  8.393e-15,
        1.012e-14,  1.259e-14,  1.538e-14,  1.961e-14,  2.619e-14,
        3.679e-14,  5.049e-14,  6.917e-14,  8.880e-14,  1.115e-13,
        1.373e-13,  1.619e-13,  1.878e-13,  2.111e-13,  2.330e-13,
        2.503e-13,  2.613e-13,  2.743e-13,  2.826e-13,  2.976e-13,
        3.162e-13,  3.360e-13,  3.491e-13,  3.541e-13,  3.595e-13,
        3.608e-13,  3.709e-13,  3.869e-13,  4.120e-13,  4.366e-13,
        4.504e-13,  4.379e-13,  3.955e-13,  3.385e-13,  2.741e-13,
        2.089e-13,  1.427e-13,  9.294e-14,  5.775e-14,  3.565e-14,
        2.210e-14,  1.398e-14,  9.194e-15,  6.363e-15,  4.644e-15,
        3.550e-15,  2.808e-15,  2.274e-15,  1.871e-15,  1.557e-15,
        1.308e-15,  1.108e-15,  9.488e-16,  8.222e-16,  7.238e-16,
        6.506e-16,  6.008e-16,  5.742e-16,  5.724e-16,  5.991e-16,
        6.625e-16,  7.775e-16,  9.734e-16,  1.306e-15,  1.880e-15,
        2.879e-15,  4.616e-15,  7.579e-15,  1.248e-14,  2.030e-14,
        3.244e-14,  5.171e-14,  7.394e-14,  9.676e-14,  1.199e-13,
        1.467e-13,  1.737e-13,  2.020e-13,  2.425e-13,  3.016e-13,
        3.700e-13,  4.617e-13,  5.949e-13,  7.473e-13,  9.378e-13,
        1.191e-12,  1.481e-12,  1.813e-12,  2.232e-12,  2.722e-12,
        3.254e-12,  3.845e-12,  4.458e-12,  5.048e-12,  5.511e-12,
        5.898e-12,  6.204e-12,  6.293e-12,  6.386e-12,  6.467e-12,
        6.507e-12,  6.466e-12,  6.443e-12,  6.598e-12,  6.873e-12,
        7.300e-12,  7.816e-12,  8.368e-12,  8.643e-12,  8.466e-12,
        7.871e-12,  6.853e-12,  5.714e-12,  4.482e-12,  3.392e-12,
        2.613e-12,  2.008e-12,  1.562e-12,  1.228e-12,  9.888e-13,
        7.646e-13,  5.769e-13,  4.368e-13,  3.324e-13,  2.508e-13,
        1.916e-13};


// CO2 continuum Ridgeway 1982, implementation of CKD_MT_1.00
// UNITS OF (CM**3/MOL)*1.E-20
const Numeric FCO2_ckd_mt_100_v1  =   -20.0;
const Numeric FCO2_ckd_mt_100_v2  = 10000.0;
const Numeric FCO2_ckd_mt_100_dv  =    10.0;
const int     FCO2_ckd_mt_100_npt =  1003;
const double  FCO2_ckd_mt_100[FCO2_ckd_mt_100_npt+addF77fields] = {
           0.000e0,    1.1110E-11, 1.0188E-11, 
           9.3516E-12, 1.0188E-11, 1.1110E-11, 1.2127E-11, 1.3251E-11,    // F17590
           1.4495E-11, 1.5872E-11, 1.7400E-11, 1.9097E-11, 2.0985E-11,    // F17600
           2.3087E-11, 2.5431E-11, 2.8051E-11, 3.0982E-11, 3.4268E-11,    // F17610
           3.7956E-11, 4.2105E-11, 4.6779E-11, 5.2056E-11, 5.8025E-11,    // F17620
           6.4791E-11, 7.2477E-11, 8.1226E-11, 9.1209E-11, 1.0263E-10,    // F17630
           1.1572E-10, 1.3078E-10, 1.4814E-10, 1.6821E-10, 1.9148E-10,    // F17640
           2.1857E-10, 2.5019E-10, 2.8723E-10, 3.3080E-10, 3.8223E-10,    // F17650
           4.4321E-10, 5.1583E-10, 6.0274E-10, 7.0725E-10, 8.3363E-10,    // F17660
           9.8735E-10, 1.1755E-09, 1.4074E-09, 1.6953E-09, 2.0557E-09,    // F17670
           2.5107E-09, 3.0909E-09, 3.8391E-09, 4.8165E-09, 6.1117E-09,    // F17680
           7.8550E-09, 1.0241E-08, 1.3593E-08, 1.8344E-08, 2.5408E-08,    // F17700
           3.6386E-08, 5.4251E-08, 8.4262E-08, 1.3273E-07, 2.1867E-07,    // F17710
           3.5007E-07, 6.0011E-07, 1.0797E-06, 1.8254E-06, 3.1621E-06,    // F17720
           4.0293E-06, 4.3683E-06, 4.4552E-06, 4.2684E-06, 3.9341E-06,    // F17730
           2.5972E-06, 1.5617E-06, 8.9063E-07, 5.0360E-07, 3.0616E-07,    // F17740
           1.9066E-07, 1.1904E-07, 7.6078E-08, 4.9304E-08, 3.3335E-08,    // F17750
           2.3494E-08, 1.7114E-08, 1.2742E-08, 9.6068E-09, 7.3706E-09,    // F17760
           5.7386E-09, 4.5302E-09, 3.6223E-09, 2.9309E-09, 2.4001E-09,    // F17770
           1.9927E-09, 1.6877E-09, 1.4602E-09, 1.2764E-09, 1.1317E-09,    // F17780
           1.0273E-09, 9.1943E-10, 8.0353E-10, 6.8746E-10, 5.9354E-10,    // F17790
           5.1722E-10, 4.4975E-10, 4.2350E-10, 4.2282E-10, 4.2610E-10,    // F17810
           4.5465E-10, 4.6166E-10, 4.3149E-10, 3.7615E-10, 3.1576E-10,    // F17820
           2.6490E-10, 1.9143E-10, 1.2885E-10, 9.4954E-11, 7.6499E-11,    // F17830
           6.4581E-11, 5.5923E-11, 4.9200E-11, 4.3813E-11, 3.9533E-11,    // F17840
           3.6338E-11, 3.4320E-11, 3.3329E-11, 3.2400E-11, 3.1700E-11,    // F17850
           3.1267E-11, 2.9940E-11, 2.7628E-11, 2.4496E-11, 2.1764E-11,    // F17860
           1.9306E-11, 1.7352E-11, 1.7292E-11, 1.8733E-11, 2.0224E-11,    // F17870
           2.2396E-11, 2.4225E-11, 2.4890E-11, 2.3513E-11, 2.0824E-11,    // F17880
           1.8642E-11, 1.5676E-11, 1.2882E-11, 1.1054E-11, 1.0074E-11,    // F17890
           9.6324E-12, 9.4910E-12, 9.5134E-12, 9.6427E-12, 9.8552E-12,    // F17900
           1.0140E-11, 1.0494E-11, 1.0915E-11, 1.1405E-11, 1.1965E-11,    // F17920
           1.2601E-11, 1.3316E-11, 1.4116E-11, 1.5006E-11, 1.5997E-11,    // F17930
           1.7092E-11, 1.8305E-11, 1.9641E-11, 2.1121E-11, 2.2744E-11,    // F17940
           2.4503E-11, 2.6419E-11, 2.8221E-11, 3.0609E-11, 3.3260E-11,    // F17950
           3.6247E-11, 3.9581E-11, 4.3279E-11, 4.7376E-11, 5.1932E-11,    // F17960
           5.7001E-11, 6.2654E-11, 6.8973E-11, 7.6058E-11, 8.4037E-11,    // F17970
           9.3081E-11, 1.0344E-10, 1.1547E-10, 1.2970E-10, 1.4659E-10,    // F17980
           1.6724E-10, 1.9481E-10, 2.3520E-10, 2.9424E-10, 3.6319E-10,    // F17990
           4.2279E-10, 4.8494E-10, 5.2296E-10, 5.6111E-10, 5.8935E-10,    // F18000
           6.0807E-10, 6.4204E-10, 6.8457E-10, 7.6709E-10, 8.7664E-10,    // F18010
           1.0183E-09, 1.2116E-09, 1.4874E-09, 1.8596E-09, 2.2742E-09,    // F18030
           2.7577E-09, 3.1932E-09, 3.6381E-09, 4.1207E-09, 4.6458E-09,    // F18040
           5.3065E-09, 6.0741E-09, 7.1942E-09, 8.7103E-09, 1.0713E-08,    // F18050
           1.3344E-08, 1.6831E-08, 2.1524E-08, 2.7967E-08, 3.7047E-08,    // F18060
           5.0312E-08, 7.0566E-08, 1.0275E-07, 1.5419E-07, 2.3309E-07,    // F18070
           3.4843E-07, 5.3194E-07, 8.7207E-07, 1.5075E-06, 2.7077E-06,    // F18080
           4.7125E-06, 7.1734E-06, 9.2381E-06, 1.1507E-05, 1.3737E-05,    // F18090
           1.4004E-05, 1.2679E-05, 1.0478E-05, 8.5684E-06, 6.1472E-06,    // F18100
           3.2424E-06, 1.5291E-06, 8.0390E-07, 4.6767E-07, 2.9170E-07,    // F18110
           1.9148E-07, 1.3076E-07, 9.2156E-08, 6.6652E-08, 4.9265E-08,    // F18120
           3.7094E-08, 2.8380E-08, 2.2019E-08, 1.7297E-08, 1.3738E-08,    // F18140
           1.1019E-08, 8.9178E-09, 7.2762E-09, 5.9810E-09, 4.9500E-09,    // F18150
           4.1226E-09, 3.4534E-09, 2.9082E-09, 2.4611E-09, 2.0922E-09,    // F18160
           1.7864E-09, 1.5313E-09, 1.3176E-09, 1.1379E-09, 9.8612E-10,    // F18170
           8.5741E-10, 7.4782E-10, 6.5416E-10, 5.7384E-10, 5.0471E-10,    // F18180
           4.4503E-10, 3.9334E-10, 3.4841E-10, 3.0927E-10, 2.7510E-10,    // F18190
           2.4519E-10, 2.1893E-10, 1.9587E-10, 1.7555E-10, 1.5762E-10,    // F18200
           1.4178E-10, 1.2772E-10, 1.1524E-10, 1.0414E-10, 9.4248E-11,    // F18210
           8.5421E-11, 7.7530E-11, 7.0466E-11, 6.4134E-11, 5.8450E-11,    // F18220
           5.3342E-11, 4.8746E-11, 4.4607E-11, 4.0874E-11, 3.7507E-11,    // F18230
           3.4466E-11, 3.1719E-11, 2.9237E-11, 2.6993E-11, 2.4968E-11,    // F18250
           2.3139E-11, 2.1494E-11, 2.0022E-11, 1.8709E-11, 1.7541E-11,    // F18260
           1.6533E-11, 1.5690E-11, 1.5027E-11, 1.4560E-11, 1.4169E-11,    // F18270
           1.3796E-11, 1.3553E-11, 1.3526E-11, 1.3567E-11, 1.3399E-11,    // F18280
           1.3149E-11, 1.3049E-11, 1.3078E-11, 1.3093E-11, 1.3168E-11,    // F18290
           1.3572E-11, 1.4383E-11, 1.5698E-11, 1.7658E-11, 2.0197E-11,    // F18300
           2.2845E-11, 2.5944E-11, 3.0250E-11, 3.5900E-11, 4.1482E-11,    // F18310
           4.6602E-11, 5.2453E-11, 5.9754E-11, 6.9308E-11, 8.0696E-11,    // F18320
           9.5737E-11, 1.1733E-10, 1.4793E-10, 1.9119E-10, 2.5355E-10,    // F18330
           3.4588E-10, 4.8343E-10, 6.9378E-10, 1.0212E-09, 1.4858E-09,    // F18340
           2.0906E-09, 3.0576E-09, 4.6318E-09, 7.1585E-09, 1.1259E-08,    // F18360
           1.7954E-08, 2.9760E-08, 4.6693E-08, 6.2035E-08, 7.4399E-08,    // F18370
           9.1705E-08, 9.9448E-08, 9.5181E-08, 8.3050E-08, 7.1756E-08,    // F18380
           6.6261E-08, 6.0357E-08, 6.6988E-08, 8.3419E-08, 9.8834E-08,    // F18390
           1.2385E-07, 1.3962E-07, 1.3651E-07, 1.1963E-07, 9.7731E-08,    // F18400
           8.0083E-08, 5.1660E-08, 2.5778E-08, 1.2600E-08, 6.8779E-09,    // F18410
           4.1161E-09, 2.6276E-09, 1.7595E-09, 1.2225E-09, 8.7493E-10,    // F18420
           6.4179E-10, 4.7987E-10, 3.6491E-10, 2.8191E-10, 2.2084E-10,    // F18430
           1.7507E-10, 1.4025E-10, 1.1344E-10, 9.2580E-11, 7.6170E-11,    // F18440
           6.3142E-11, 5.2694E-11, 4.4260E-11, 3.7421E-11, 3.1847E-11,    // F18450
           2.7263E-11, 2.3352E-11, 2.0081E-11, 1.7332E-11, 1.5000E-11,    // F18470
           1.2978E-11, 1.1204E-11, 9.7513E-12, 8.5300E-12, 7.4888E-12,    // F18480
           6.5947E-12, 5.8231E-12, 5.1548E-12, 4.5739E-12, 4.0675E-12,    // F18490
           3.6250E-12, 3.2371E-12, 2.8963E-12, 2.5964E-12, 2.3316E-12,    // F18500
           2.0975E-12, 1.8902E-12, 1.7061E-12, 1.5425E-12, 1.3967E-12,    // F18510
           1.2665E-12, 1.1503E-12, 1.0463E-12, 9.5319E-13, 8.6963E-13,    // F18520
           7.9461E-13, 7.2718E-13, 6.6654E-13, 6.1201E-13, 5.6296E-13,    // F18530
           5.1894E-13, 4.7969E-13, 4.4494E-13, 4.1320E-13, 3.8529E-13,    // F18540
           3.6202E-13, 3.4320E-13, 3.2546E-13, 3.0741E-13, 2.9156E-13,    // F18550
           2.7819E-13, 2.6576E-13, 2.5327E-13, 2.4319E-13, 2.3770E-13,    // F18560
           2.3645E-13, 2.3967E-13, 2.4960E-13, 2.6858E-13, 2.9679E-13,    // F18580
           3.3247E-13, 3.8487E-13, 4.7576E-13, 6.1833E-13, 8.0740E-13,    // F18590
           1.0267E-12, 1.2291E-12, 1.4710E-12, 1.7211E-12, 1.8251E-12,    // F18600
           1.8982E-12, 1.9768E-12, 2.1877E-12, 2.5008E-12, 3.0545E-12,    // F18610
           4.1513E-12, 5.7469E-12, 7.7913E-12, 1.0873E-11, 1.5538E-11,    // F18620
           2.2838E-11, 3.4153E-11, 4.9751E-11, 7.0591E-11, 1.0794E-10,    // F18630
           1.7287E-10, 2.6554E-10, 3.5250E-10, 4.1952E-10, 5.1979E-10,    // F18640
           5.7649E-10, 5.6168E-10, 5.0014E-10, 4.3670E-10, 4.0057E-10,    // F18650
           3.5169E-10, 3.7578E-10, 5.5054E-10, 8.8962E-10, 1.2940E-09,    // F18660
           1.6293E-09, 2.0553E-09, 2.3945E-09, 2.3926E-09, 2.1385E-09,    // F18670
           1.7637E-09, 1.4623E-09, 1.0150E-09, 5.5612E-10, 3.5162E-10,    // F18690
           3.4009E-10, 4.1744E-10, 5.0009E-10, 6.0748E-10, 7.3258E-10,    // F18700
           7.6553E-10, 7.2066E-10, 6.1317E-10, 5.1585E-10, 3.9136E-10,    // F18710
           2.2991E-10, 1.2590E-10, 6.9549E-11, 3.8699E-11, 2.2976E-11,    // F18720
           1.4702E-11, 9.9989E-12, 7.1233E-12, 5.2612E-12, 4.0298E-12,    // F18730
           3.2395E-12, 2.7932E-12, 2.6331E-12, 2.7835E-12, 3.3167E-12,    // F18740
           3.3581E-12, 3.3404E-12, 3.1243E-12, 2.8459E-12, 2.4092E-12,    // F18750
           1.5349E-12, 9.7039E-13, 5.8611E-13, 3.9686E-13, 2.9332E-13,    // F18760
           2.2795E-13, 1.8432E-13, 1.5287E-13, 1.2898E-13, 1.1019E-13,    // F18770
           9.5041E-14, 8.2617E-14, 7.2310E-14, 6.3711E-14, 5.6561E-14,    // F18780
           5.0763E-14, 4.6525E-14, 4.4418E-14, 4.4681E-14, 4.7199E-14,    // F18800
           5.0389E-14, 5.3620E-14, 6.0817E-14, 6.0192E-14, 5.5878E-14,    // F18810
           4.9874E-14, 4.3955E-14, 3.9854E-14, 3.1697E-14, 3.1135E-14,    // F18820
           3.4683E-14, 3.8789E-14, 4.6932E-14, 5.0213E-14, 4.7156E-14,    // F18830
           4.2130E-14, 3.5554E-14, 3.0465E-14, 1.9216E-14, 1.1378E-14,    // F18840
           8.2878E-15, 6.8260E-15, 6.0960E-15, 5.8135E-15, 5.9618E-15,    // F18850
           6.8295E-15, 9.2943E-15, 1.2572E-14, 1.4837E-14, 1.8595E-14,    // F18860
           2.1533E-14, 2.2008E-14, 2.1305E-14, 1.9743E-14, 2.0413E-14,    // F18870
           2.1131E-14, 2.5346E-14, 3.3709E-14, 4.3995E-14, 5.8911E-14,    // F18880
           7.8451E-14, 1.0537E-13, 1.4559E-13, 2.0405E-13, 2.6734E-13,    // F18890
           3.5029E-13, 4.9788E-13, 7.3207E-13, 1.0979E-12, 1.4960E-12,    // F18910
           1.7906E-12, 2.2171E-12, 2.5369E-12, 2.5873E-12, 2.3871E-12,    // F18920
           2.0730E-12, 1.9095E-12, 1.6227E-12, 1.3981E-12, 1.5228E-12,    // F18930
           2.0956E-12, 3.2493E-12, 5.2740E-12, 8.6666E-12, 1.2672E-11,    // F18940
           1.5725E-11, 1.9496E-11, 2.2858E-11, 2.2939E-11, 2.0597E-11,    // F18950
           1.7021E-11, 1.4456E-11, 1.0794E-11, 7.1327E-12, 6.5438E-12,    // F18960
           8.8057E-12, 1.2311E-11, 1.5284E-11, 1.9273E-11, 2.2796E-11,    // F18970
           2.3156E-11, 2.0914E-11, 1.7298E-11, 1.4424E-11, 1.0127E-11,    // F18980
           5.2952E-12, 2.5759E-12, 1.4304E-12, 9.4758E-13, 7.9895E-13,    // F18990
           9.1124E-13, 1.2297E-12, 1.5898E-12, 1.9056E-12, 2.3905E-12,    // F19000
           2.6695E-12, 2.6297E-12, 2.3467E-12, 2.0058E-12, 1.6773E-12,    // F19020
           1.1327E-12, 6.7331E-13, 4.0954E-13, 2.5152E-13, 1.4491E-13,    // F19030
           9.0916E-14, 6.6510E-14, 5.9022E-14, 6.4403E-14, 8.3126E-14,    // F19040
           1.2409E-13, 1.5153E-13, 1.6909E-13, 1.7938E-13, 1.9169E-13,    // F19050
           2.1173E-13, 2.1941E-13, 2.6360E-13, 3.5956E-13, 4.8369E-13,    // F19060
           5.9657E-13, 7.4062E-13, 8.9452E-13, 8.7899E-13, 8.2012E-13,    // F19070
           7.4109E-13, 6.9845E-13, 6.3130E-13, 5.6538E-13, 6.9516E-13,    // F19080
           9.9486E-13, 1.5226E-12, 2.4155E-12, 3.9119E-12, 6.3541E-12,    // F19090
           1.0075E-11, 1.5903E-11, 2.5091E-11, 3.6282E-11, 4.6076E-11,    // F19100
           5.6240E-11, 7.1126E-11, 7.0230E-11, 6.3642E-11, 5.3722E-11,    // F19110
           4.4651E-11, 3.4409E-11, 1.5287E-11, 7.2479E-12, 3.9218E-12,    // F19130
           2.3172E-12, 1.4585E-12, 9.6297E-13, 6.6017E-13, 4.6655E-13,    // F19140
           3.3814E-13, 2.5034E-13, 1.8874E-13, 1.4457E-13, 1.1228E-13,    // F19150
           8.8284E-14, 7.0188E-14, 5.6365E-14, 4.5685E-14, 3.7357E-14,    // F19160
           3.0817E-14, 2.5674E-14, 2.1679E-14, 1.8780E-14, 1.7243E-14,    // F19170
           1.6273E-14, 1.5201E-14, 1.5091E-14, 1.4725E-14, 1.3668E-14,    // F19180
           1.1940E-14, 1.0097E-14, 8.8905E-15, 7.1475E-15, 5.8080E-15,    // F19190
           5.5216E-15, 5.9338E-15, 7.1932E-15, 9.9780E-15, 1.6167E-14,    // F19200
           2.9100E-14, 5.2355E-14, 8.4889E-14, 1.1311E-13, 1.4192E-13,    // F19210
           1.7648E-13, 1.8657E-13, 1.7498E-13, 1.4877E-13, 1.2578E-13,    // F19220
           1.0051E-13, 6.7213E-14, 5.4750E-14, 7.0454E-14, 1.1351E-13,    // F19240
           1.8015E-13, 2.4825E-13, 3.0875E-13, 3.9200E-13, 4.2550E-13,    // F19250
           4.0067E-13, 3.4438E-13, 2.8204E-13, 2.2432E-13, 1.3172E-13,    // F19260
           6.2820E-14, 3.6474E-14, 2.9409E-14, 3.4164E-14, 4.8300E-14,    // F19270
           6.4140E-14, 7.7284E-14, 9.7973E-14, 1.0969E-13, 1.0580E-13,    // F19280
           9.2070E-14, 7.5008E-14, 6.1722E-14, 3.8874E-14, 1.9007E-14,    // F19290
           9.6765E-15, 5.5169E-15, 3.5254E-15, 2.5012E-15, 2.0013E-15,    // F19300
           1.8810E-15, 2.2143E-15, 3.5332E-15, 5.7552E-15, 7.3359E-15,    // F19310
           8.3292E-15, 9.9174E-15, 1.0930E-14, 1.1185E-14, 1.0884E-14,    // F19320
           1.0577E-14, 1.1048E-14, 1.1611E-14, 1.1128E-14, 1.0729E-14,    // F19330
           1.0248E-14, 1.0630E-14, 1.1793E-14, 1.3977E-14, 1.9857E-14,    // F19350
           2.9182E-14, 4.2229E-14, 6.2710E-14, 9.0717E-14, 1.2561E-13,    // F19360
           1.6951E-13, 2.2520E-13, 3.2470E-13, 4.5178E-13, 6.3104E-13,    // F19370
           8.7521E-13, 1.1073E-12, 1.3534E-12, 1.6954E-12, 1.7005E-12,    // F19380
           1.5993E-12, 1.4416E-12, 1.3280E-12, 1.2760E-12, 1.1076E-12,    // F19390
           1.2850E-12, 1.6208E-12, 1.9527E-12, 2.4941E-12, 2.5077E-12,    // F19400
           2.3156E-12, 2.0069E-12, 1.6301E-12, 1.2885E-12, 5.9863E-13,    // F19410
           2.8012E-13, 1.5065E-13, 8.8802E-14, 5.5888E-14, 3.6951E-14,    // F19420
           2.5393E-14, 1.8001E-14, 1.3093E-14, 9.7308E-15, 7.3665E-15,    // F19430
           5.6662E-15, 4.4194E-15, 3.4897E-15, 2.7857E-15, 2.2457E-15,    // F19440
           1.8264E-15, 1.4973E-15, 1.2365E-15, 1.0280E-15, 8.5996E-16,    // F19460
           7.2345E-16, 6.1182E-16, 5.1994E-16, 4.4388E-16, 3.8055E-16,    // F19470
           3.2756E-16, 2.8300E-16, 2.4537E-16, 2.1347E-16, 1.8630E-16,    // F19480
           1.6307E-16, 1.4314E-16, 1.2599E-16, 1.1117E-16, 9.8344E-17,    // F19490
           8.7197E-17, 7.7487E-17, 6.9004E-17, 6.1577E-17, 5.5060E-17,    // F19500
           4.9325E-17, 4.4271E-17, 3.9810E-17, 3.5861E-17, 3.2361E-17,    // F19510
           2.9252E-17, 2.6487E-17, 2.4023E-17, 2.1826E-17, 1.9862E-17,    // F19520
           1.8107E-17, 1.6536E-17, 1.5129E-17, 1.3869E-17, 1.2739E-17,    // F19530
           1.1726E-17, 1.0820E-17, 1.0009E-17, 9.2846E-18, 8.6398E-18,    // F19540
           8.0682E-18, 7.5641E-18, 7.1229E-18, 6.7411E-18, 6.4161E-18,    // F19550
           6.1455E-18, 5.9290E-18, 5.7662E-18, 5.6574E-18, 5.6049E-18,    // F19570
           5.6112E-18, 5.6811E-18, 5.8200E-18, 6.0364E-18, 6.3405E-18,    // F19580
           6.7450E-18, 7.2674E-18, 7.9298E-18, 8.7613E-18, 9.8010E-18,    // F19590
           1.1086E-17, 1.2686E-17, 1.4679E-17, 1.7177E-17, 2.0335E-17,    // F19600
           2.4384E-17, 2.9538E-17, 3.6416E-17, 4.5520E-17, 5.7788E-17,    // F19610
           7.4676E-17, 9.8513E-17, 1.3323E-16, 1.8570E-16, 2.6897E-16,    // F19620
           4.0958E-16, 6.6785E-16, 1.2064E-15, 2.4023E-15, 4.3240E-15,    // F19630
           6.6353E-15, 8.6393E-15, 1.1433E-14, 1.3946E-14, 1.3611E-14,    // F19640
           1.2557E-14, 1.0934E-14, 1.0039E-14, 8.5099E-15, 7.9557E-15,    // F19650
           1.1346E-14, 1.8512E-14, 2.9285E-14, 4.1585E-14, 5.2809E-14,    // F19660
           7.0377E-14, 7.8094E-14, 7.3735E-14, 6.5845E-14, 5.5023E-14,    // F19680
           4.6866E-14, 2.7430E-14, 1.5975E-14, 1.4522E-14, 1.7075E-14,    // F19690
           2.0408E-14, 2.5119E-14, 3.1194E-14, 3.0280E-14, 2.7676E-14,    // F19700
           2.3344E-14, 1.9466E-14, 1.4140E-14, 6.2087E-15, 3.0307E-15,    // F19710
           1.6815E-15, 1.0169E-15, 6.5448E-16, 4.4162E-16, 3.0928E-16,    // F19720
           2.2320E-16, 1.6511E-16, 1.2471E-16, 9.5881E-17, 7.4850E-17,    // F19730
           5.9216E-17, 4.7400E-17, 3.8338E-17, 3.1298E-17, 2.5765E-17,    // F19740
           2.1371E-17, 1.7848E-17, 1.5000E-17, 1.2679E-17, 1.0774E-17,    // F19750
           9.2002E-18, 7.8922E-18, 6.7987E-18, 5.8800E-18, 5.1042E-18,    // F19760
           4.4461E-18, 3.8855E-18, 3.4060E-18, 2.9944E-18, 2.6397E-18,    // F19770
           2.3331E-18};

// CKD_MT 1.00 implementation of N2-N2 model of 
// Borysow, A, and L. Frommhold, 
//  "Collision-induced rototranslational absorption spectra of N2-N2
//  pairs for temperatures from 50 to 300 K", The
//  Astrophysical Journal, 311, 1043-1057, 1986.
// absorption coefficient in units of [cm^-1 Amagat^-2] 
// these data are for T=296K
const Numeric N2N2_CT296_ckd_mt_100_v1  = -10.0;
const Numeric N2N2_CT296_ckd_mt_100_v2  = 350.0;
const Numeric N2N2_CT296_ckd_mt_100_dv  =   5.0;
const int     N2N2_CT296_ckd_mt_100_npt =  73;
const double  N2N2_CT296_ckd_mt_100[N2N2_CT296_ckd_mt_100_npt+addF77fields] = {
           0.0000e0, 
           0.4303E-06, 0.4850E-06, 0.4979E-06, 0.4850E-06, 0.4303E-06,
           0.3715E-06, 0.3292E-06, 0.3086E-06, 0.2920E-06, 0.2813E-06,
           0.2804E-06, 0.2738E-06, 0.2726E-06, 0.2724E-06, 0.2635E-06,
           0.2621E-06, 0.2547E-06, 0.2428E-06, 0.2371E-06, 0.2228E-06,
           0.2100E-06, 0.1991E-06, 0.1822E-06, 0.1697E-06, 0.1555E-06,
           0.1398E-06, 0.1281E-06, 0.1138E-06, 0.1012E-06, 0.9078E-07,
           0.7879E-07, 0.6944E-07, 0.6084E-07, 0.5207E-07, 0.4540E-07,
           0.3897E-07, 0.3313E-07, 0.2852E-07, 0.2413E-07, 0.2045E-07,
           0.1737E-07, 0.1458E-07, 0.1231E-07, 0.1031E-07, 0.8586E-08,
           0.7162E-08, 0.5963E-08, 0.4999E-08, 0.4226E-08, 0.3607E-08,
           0.3090E-08, 0.2669E-08, 0.2325E-08, 0.2024E-08, 0.1783E-08,
           0.1574E-08, 0.1387E-08, 0.1236E-08, 0.1098E-08, 0.9777E-09,
           0.8765E-09, 0.7833E-09, 0.7022E-09, 0.6317E-09, 0.5650E-09,
           0.5100E-09, 0.4572E-09, 0.4115E-09, 0.3721E-09, 0.3339E-09,
           0.3005E-09, 0.2715E-09, 0.2428E-09};


// CKD_MT 1.00 implementation of N2-N2 model of 
// Borysow, A, and L. Frommhold, 
//  "Collision-induced rototranslational absorption spectra of N2-N2
//  pairs for temperatures from 50 to 300 K", The
//  Astrophysical Journal, 311, 1043-1057, 1986.
// absorption coefficient in units of [cm^-1 Amagat^-2] 
// these data are for T=220K
const Numeric N2N2_CT220_ckd_mt_100_v1  = -10.0;
const Numeric N2N2_CT220_ckd_mt_100_v2  = 350.0;
const Numeric N2N2_CT220_ckd_mt_100_dv  =   5.0;
const int     N2N2_CT220_ckd_mt_100_npt =  73;
const double  N2N2_CT220_ckd_mt_100[N2N2_CT220_ckd_mt_100_npt+addF77fields] = {
           0.0000e0,
           0.4946E-06, 0.5756E-06, 0.5964E-06, 0.5756E-06, 0.4946E-06,
           0.4145E-06, 0.3641E-06, 0.3482E-06, 0.3340E-06, 0.3252E-06,
           0.3299E-06, 0.3206E-06, 0.3184E-06, 0.3167E-06, 0.2994E-06,
           0.2943E-06, 0.2794E-06, 0.2582E-06, 0.2468E-06, 0.2237E-06,
           0.2038E-06, 0.1873E-06, 0.1641E-06, 0.1474E-06, 0.1297E-06,
           0.1114E-06, 0.9813E-07, 0.8309E-07, 0.7059E-07, 0.6068E-07,
           0.5008E-07, 0.4221E-07, 0.3537E-07, 0.2885E-07, 0.2407E-07,
           0.1977E-07, 0.1605E-07, 0.1313E-07, 0.1057E-07, 0.8482E-08,
           0.6844E-08, 0.5595E-08, 0.4616E-08, 0.3854E-08, 0.3257E-08,
           0.2757E-08, 0.2372E-08, 0.2039E-08, 0.1767E-08, 0.1548E-08,
           0.1346E-08, 0.1181E-08, 0.1043E-08, 0.9110E-09, 0.8103E-09,
           0.7189E-09, 0.6314E-09, 0.5635E-09, 0.4976E-09, 0.4401E-09,
           0.3926E-09, 0.3477E-09, 0.3085E-09, 0.2745E-09, 0.2416E-09,
           0.2155E-09, 0.1895E-09, 0.1678E-09, 0.1493E-09, 0.1310E-09,
           0.1154E-09, 0.1019E-09, 0.8855E-10};


// CKD_MT 1.00 implementation of N2-N2 model of 
// Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and J._M. Hartmann,
// Infrared collision-induced absorption by N2 near 4.3 microns for
// atmospheric applications: measurements and emprirical modeling, 
// Appl. Optics, 35, 5911-5917, (1996).
const Numeric N2N2_N2F_ckd_mt_100_v1  = 2085.000;
const Numeric N2N2_N2F_ckd_mt_100_v2  = 2670.000;
const Numeric N2N2_N2F_ckd_mt_100_dv  =    5.000;
const int     N2N2_N2F_ckd_mt_100_npt =  118;
const double  N2N2_N2F_ckd_mt_100[N2N2_N2F_ckd_mt_100_npt+addF77fields] = {
            0.000E+00, 
            0.000E+00,  2.000E-10,  5.200E-09,  1.020E-08,  1.520E-08,  
            2.020E-08,  2.520E-08,  3.020E-08,  4.450E-08,  5.220E-08,  
            6.460E-08,  7.750E-08,  9.030E-08,  1.060E-07,  1.210E-07,  
            1.370E-07,  1.570E-07,  1.750E-07,  2.010E-07,  2.300E-07,  
            2.590E-07,  2.950E-07,  3.260E-07,  3.660E-07,  4.050E-07,  
            4.470E-07,  4.920E-07,  5.340E-07,  5.840E-07,  6.240E-07,  
            6.670E-07,  7.140E-07,  7.260E-07,  7.540E-07,  7.840E-07,  
            8.090E-07,  8.420E-07,  8.620E-07,  8.870E-07,  9.110E-07,  
            9.360E-07,  9.760E-07,  1.030E-06,  1.110E-06,  1.230E-06,  
            1.390E-06,  1.610E-06,  1.760E-06,  1.940E-06,  1.970E-06,  
            1.870E-06,  1.750E-06,  1.560E-06,  1.420E-06,  1.350E-06,  
            1.320E-06,  1.290E-06,  1.290E-06,  1.290E-06,  1.300E-06,  
            1.320E-06,  1.330E-06,  1.340E-06,  1.350E-06,  1.330E-06,  
            1.310E-06,  1.290E-06,  1.240E-06,  1.200E-06,  1.160E-06,  
            1.100E-06,  1.040E-06,  9.960E-07,  9.380E-07,  8.630E-07,  
            7.980E-07,  7.260E-07,  6.550E-07,  5.940E-07,  5.350E-07,  
            4.740E-07,  4.240E-07,  3.770E-07,  3.330E-07,  2.960E-07,  
            2.630E-07,  2.340E-07,  2.080E-07,  1.850E-07,  1.670E-07,  
            1.470E-07,  1.320E-07,  1.200E-07,  1.090E-07,  9.850E-08,  
            9.080E-08,  8.180E-08,  7.560E-08,  6.850E-08,  6.140E-08,  
            5.830E-08,  5.770E-08,  5.000E-08,  4.320E-08,  3.140E-08,  
            2.890E-08,  2.640E-08,  2.390E-08,  2.140E-08,  1.890E-08,  
            1.640E-08,  1.390E-08,  1.140E-08,  8.900E-09,  6.400E-09,  
            3.900E-09,  1.400E-09,  0.000E+00};

//     temperature coefficients:
const double  N2N2_N2Ft_ckd_mt_100[N2N2_N2F_ckd_mt_100_npt+addF77fields] = {
            0.000E+00, 
            1.040E+03,  1.010E+03,  9.800E+02,  9.500E+02,  9.200E+02,  
            8.900E+02,  8.600E+02,  8.300E+02,  8.020E+02,  7.610E+02,  
            7.220E+02,  6.790E+02,  6.460E+02,  6.090E+02,  5.620E+02,  
            5.110E+02,  4.720E+02,  4.360E+02,  4.060E+02,  3.770E+02,  
            3.550E+02,  3.380E+02,  3.190E+02,  2.990E+02,  2.780E+02,  
            2.550E+02,  2.330E+02,  2.080E+02,  1.840E+02,  1.490E+02,  
            1.070E+02,  6.600E+01,  2.500E+01, -1.300E+01, -4.900E+01,  
           -8.200E+01, -1.040E+02, -1.190E+02, -1.300E+02, -1.390E+02,  
           -1.440E+02, -1.460E+02, -1.460E+02, -1.470E+02, -1.480E+02,  
           -1.500E+02, -1.530E+02, -1.600E+02, -1.690E+02, -1.810E+02,  
           -1.890E+02, -1.950E+02, -2.000E+02, -2.050E+02, -2.090E+02,  
           -2.110E+02, -2.100E+02, -2.100E+02, -2.090E+02, -2.050E+02,  
           -1.990E+02, -1.900E+02, -1.800E+02, -1.680E+02, -1.570E+02,  
           -1.430E+02, -1.260E+02, -1.080E+02, -8.900E+01, -6.300E+01,  
           -3.200E+01,  1.000E+00,  3.500E+01,  6.500E+01,  9.500E+01,  
            1.210E+02,  1.410E+02,  1.520E+02,  1.610E+02,  1.640E+02,  
            1.640E+02,  1.610E+02,  1.550E+02,  1.480E+02,  1.430E+02,  
            1.370E+02,  1.330E+02,  1.310E+02,  1.330E+02,  1.390E+02,  
            1.500E+02,  1.650E+02,  1.870E+02,  2.130E+02,  2.480E+02,  
            2.840E+02,  3.210E+02,  3.720E+02,  4.490E+02,  5.140E+02,  
            5.690E+02,  6.090E+02,  6.420E+02,  6.730E+02,  7.000E+02,  
            7.300E+02,  7.600E+02,  7.900E+02,  8.200E+02,  8.500E+02,  
            8.800E+02,  9.100E+02,  9.400E+02,  9.700E+02,  1.000E+03,  
            1.030E+03,  1.060E+03,  1.090E+03};


// CKD_MT 1.00 implementation of oxygen collision induced fundamental model of 
// F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, 
// J.-M. Hartmann, Ch. Boulet, 
// "Infrared collision-induced absorption by O2 near 6.4 microns for
// atmospheric applications: measurements and emprirical modeling", 
// Appl. Optics, 35, 5911-5917, (1996).
const Numeric O2O2_O2F_ckd_mt_100_v1  = 1340.000;
const Numeric O2O2_O2F_ckd_mt_100_v2  = 1850.000;
const Numeric O2O2_O2F_ckd_mt_100_dv  =    5.000;
const int     O2O2_O2F_ckd_mt_100_npt =  103;
const double  O2O2_O2Fo_ckd_mt_100[O2O2_O2F_ckd_mt_100_npt+addF77fields] = {
            0.000E+00, 
            0.000E+00,  9.744E-09,  2.256E-08,  3.538E-08,  4.820E-08,  
            6.100E-08,  7.400E-08,  8.400E-08,  9.600E-08,  1.200E-07,  
            1.620E-07,  2.080E-07,  2.460E-07,  2.850E-07,  3.140E-07,  
            3.800E-07,  4.440E-07,  5.000E-07,  5.710E-07,  6.730E-07,  
            7.680E-07,  8.530E-07,  9.660E-07,  1.100E-06,  1.210E-06,  
            1.330E-06,  1.470E-06,  1.590E-06,  1.690E-06,  1.800E-06,  
            1.920E-06,  2.040E-06,  2.150E-06,  2.260E-06,  2.370E-06,  
            2.510E-06,  2.670E-06,  2.850E-06,  3.070E-06,  3.420E-06,  
            3.830E-06,  4.200E-06,  4.450E-06,  4.600E-06,  4.530E-06,  
            4.280E-06,  3.960E-06,  3.680E-06,  3.480E-06,  3.350E-06,  
            3.290E-06,  3.250E-06,  3.230E-06,  3.230E-06,  3.210E-06,  
            3.190E-06,  3.110E-06,  3.030E-06,  2.910E-06,  2.800E-06,  
            2.650E-06,  2.510E-06,  2.320E-06,  2.130E-06,  1.930E-06,  
            1.760E-06,  1.590E-06,  1.420E-06,  1.250E-06,  1.110E-06,  
            9.900E-07,  8.880E-07,  7.910E-07,  6.780E-07,  5.870E-07,  
            5.240E-07,  4.640E-07,  4.030E-07,  3.570E-07,  3.200E-07,  
            2.900E-07,  2.670E-07,  2.420E-07,  2.150E-07,  1.820E-07,  
            1.600E-07,  1.460E-07,  1.280E-07,  1.030E-07,  8.700E-08,  
            8.100E-08,  7.100E-08,  6.400E-08,  5.807E-08,  5.139E-08,  
            4.496E-08,  3.854E-08,  3.212E-08,  2.569E-08,  1.927E-08,  
            1.285E-08,  6.423E-09,  0.000E+00};
                          
const double  O2O2_O2Ft_ckd_mt_100[O2O2_O2F_ckd_mt_100_npt+addF77fields] = {
            0.000E+00, 
            4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  
            4.670E+02,  4.000E+02,  3.150E+02,  3.790E+02,  3.680E+02,  
            4.750E+02,  5.210E+02,  5.310E+02,  5.120E+02,  4.420E+02,  
            4.440E+02,  4.300E+02,  3.810E+02,  3.350E+02,  3.240E+02,  
            2.960E+02,  2.480E+02,  2.150E+02,  1.930E+02,  1.580E+02,  
            1.270E+02,  1.010E+02,  7.100E+01,  3.100E+01, -6.000E+00,  
           -2.600E+01, -4.700E+01, -6.300E+01, -7.900E+01, -8.800E+01,  
           -8.800E+01, -8.700E+01, -9.000E+01, -9.800E+01, -9.900E+01,  
           -1.090E+02, -1.340E+02, -1.600E+02, -1.670E+02, -1.640E+02,  
           -1.580E+02, -1.530E+02, -1.510E+02, -1.560E+02, -1.660E+02,  
           -1.680E+02, -1.730E+02, -1.700E+02, -1.610E+02, -1.450E+02,  
           -1.260E+02, -1.080E+02, -8.400E+01, -5.900E+01, -2.900E+01,  
            4.000E+00,  4.100E+01,  7.300E+01,  9.700E+01,  1.230E+02,  
            1.590E+02,  1.980E+02,  2.200E+02,  2.420E+02,  2.560E+02,  
            2.810E+02,  3.110E+02,  3.340E+02,  3.190E+02,  3.130E+02,  
            3.210E+02,  3.230E+02,  3.100E+02,  3.150E+02,  3.200E+02,  
            3.350E+02,  3.610E+02,  3.780E+02,  3.730E+02,  3.380E+02,  
            3.190E+02,  3.460E+02,  3.220E+02,  2.910E+02,  2.900E+02,  
            3.500E+02,  3.710E+02,  5.040E+02,  4.000E+02,  4.000E+02,  
            4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  
            4.000E+02,  4.000E+02,  4.000E+02};



// CKD_MT 1.00 implementation of oxygen v0<-v0 band model of 
//   Mate et al. over the spectral region 7550-8486 cm-1: 
//   B. Mate, C. Lugez, G.T. Fraser, W.J. Lafferty,
//   "Absolute Intensities for the O2 1.27 micron
//   continuum absorption",  
//   J. Geophys. Res., 104, 30,585-30,590, 1999. 
//
// The units of these continua coefficients are  1 / (amagat_O2*amagat_air)
//
// Also, refer to the paper "Observed  Atmospheric
// Collision Induced Absorption in Near Infrared Oxygen Bands",
// Mlawer, Clough, Brown, Stephen, Landry, Goldman, & Murcray,
// Journal of Geophysical Research (1998).
//
// (comment: v0<--v0 band around 7536 to 8500 cm-1 (1.27mu))
const Numeric O2_00_ckd_mt_100_v1  = 7536.000e0;
const Numeric O2_00_ckd_mt_100_v2  = 8500.000e0;
const Numeric O2_00_ckd_mt_100_dv  =    2.000e0;
const int     O2_00_ckd_mt_100_npt =  483;
const double  O2_00_ckd_mt_100[O2_00_ckd_mt_100_npt+addF77fields] = {
            0.000E+00, 
            0.000E+00,  4.355E-11,  8.709E-11,  1.742E-10,  3.484E-10,  
            6.968E-10,  1.394E-09,  2.787E-09,  3.561E-09,  3.314E-09,  
            3.368E-09,  3.435E-09,  2.855E-09,  3.244E-09,  3.447E-09,  
            3.891E-09,  4.355E-09,  3.709E-09,  4.265E-09,  4.772E-09,  
            4.541E-09,  4.557E-09,  4.915E-09,  4.688E-09,  5.282E-09,  
            5.755E-09,  5.096E-09,  5.027E-09,  4.860E-09,  4.724E-09,  
            5.048E-09,  5.248E-09,  5.473E-09,  4.852E-09,  5.362E-09,  
            6.157E-09,  6.150E-09,  6.347E-09,  6.388E-09,  6.213E-09,  
            6.521E-09,  8.470E-09,  8.236E-09,  8.269E-09,  8.776E-09,  
            9.122E-09,  9.189E-09,  9.778E-09,  8.433E-09,  9.964E-09,  
            9.827E-09,  1.064E-08,  1.063E-08,  1.031E-08,  1.098E-08,  
            1.156E-08,  1.295E-08,  1.326E-08,  1.467E-08,  1.427E-08,  
            1.452E-08,  1.456E-08,  1.554E-08,  1.605E-08,  1.659E-08,  
            1.754E-08,  1.757E-08,  1.876E-08,  1.903E-08,  1.876E-08,  
            1.869E-08,  2.036E-08,  2.203E-08,  2.221E-08,  2.284E-08,  
            2.288E-08,  2.394E-08,  2.509E-08,  2.663E-08,  2.720E-08,  
            2.839E-08,  2.923E-08,  2.893E-08,  2.949E-08,  2.962E-08,  
            3.057E-08,  3.056E-08,  3.364E-08,  3.563E-08,  3.743E-08,  
            3.813E-08,  3.946E-08,  4.082E-08,  4.201E-08,  4.297E-08,  
            4.528E-08,  4.587E-08,  4.704E-08,  4.962E-08,  5.115E-08,  
            5.341E-08,  5.365E-08,  5.557E-08,  5.891E-08,  6.084E-08,  
            6.270E-08,  6.448E-08,  6.622E-08,  6.939E-08,  7.233E-08,  
            7.498E-08,  7.749E-08,  8.027E-08,  8.387E-08,  8.605E-08,  
            8.888E-08,  9.277E-08,  9.523E-08,  9.880E-08,  1.037E-07,  
            1.076E-07,  1.114E-07,  1.151E-07,  1.203E-07,  1.246E-07,  
            1.285E-07,  1.345E-07,  1.408E-07,  1.465E-07,  1.519E-07,  
            1.578E-07,  1.628E-07,  1.685E-07,  1.760E-07,  1.847E-07,  
            1.929E-07,  2.002E-07,  2.070E-07,  2.177E-07,  2.262E-07,  
            2.365E-07,  2.482E-07,  2.587E-07,  2.655E-07,  2.789E-07,  
            2.925E-07,  3.023E-07,  3.153E-07,  3.296E-07,  3.409E-07,  
            3.532E-07,  3.680E-07,  3.859E-07,  3.951E-07,  4.074E-07,  
            4.210E-07,  4.381E-07,  4.588E-07,  4.792E-07,  4.958E-07,  
            5.104E-07,  5.271E-07,  5.501E-07,  5.674E-07,  5.913E-07,  
            6.243E-07,  6.471E-07,  6.622E-07,  6.831E-07,  6.987E-07,  
            7.159E-07,  7.412E-07,  7.698E-07,  7.599E-07,  7.600E-07,  
            7.918E-07,  8.026E-07,  8.051E-07,  8.049E-07,  7.914E-07,  
            7.968E-07,  7.945E-07,  7.861E-07,  7.864E-07,  7.741E-07,  
            7.675E-07,  7.592E-07,  7.400E-07,  7.362E-07,  7.285E-07,  
            7.173E-07,  6.966E-07,  6.744E-07,  6.597E-07,  6.413E-07,  
            6.265E-07,  6.110E-07,  5.929E-07,  5.717E-07,  5.592E-07,  
            5.411E-07,  5.235E-07,  5.061E-07,  4.845E-07,  4.732E-07,  
            4.593E-07,  4.467E-07,  4.328E-07,  4.161E-07,  4.035E-07,  
            3.922E-07,  3.820E-07,  3.707E-07,  3.585E-07,  3.475E-07,  
            3.407E-07,  3.317E-07,  3.226E-07,  3.134E-07,  3.016E-07,  
            2.969E-07,  2.894E-07,  2.814E-07,  2.749E-07,  2.657E-07,  
            2.610E-07,  2.536E-07,  2.467E-07,  2.394E-07,  2.337E-07,  
            2.302E-07,  2.241E-07,  2.191E-07,  2.140E-07,  2.093E-07,  
            2.052E-07,  1.998E-07,  1.963E-07,  1.920E-07,  1.862E-07,  
            1.834E-07,  1.795E-07,  1.745E-07,  1.723E-07,  1.686E-07,  
            1.658E-07,  1.629E-07,  1.595E-07,  1.558E-07,  1.523E-07,  
            1.498E-07,  1.466E-07,  1.452E-07,  1.431E-07,  1.408E-07,  
            1.381E-07,  1.362E-07,  1.320E-07,  1.298E-07,  1.262E-07,  
            1.247E-07,  1.234E-07,  1.221E-07,  1.197E-07,  1.176E-07,  
            1.142E-07,  1.121E-07,  1.099E-07,  1.081E-07,  1.073E-07,  
            1.061E-07,  1.041E-07,  1.019E-07,  9.969E-08,  9.727E-08,  
            9.642E-08,  9.487E-08,  9.318E-08,  9.116E-08,  9.046E-08,  
            8.827E-08,  8.689E-08,  8.433E-08,  8.324E-08,  8.204E-08,  
            8.036E-08,  7.951E-08,  7.804E-08,  7.524E-08,  7.392E-08,  
            7.227E-08,  7.176E-08,  6.975E-08,  6.914E-08,  6.859E-08,  
            6.664E-08,  6.506E-08,  6.368E-08,  6.262E-08,  6.026E-08,  
            6.002E-08,  5.866E-08,  5.867E-08,  5.641E-08,  5.589E-08,  
            5.499E-08,  5.309E-08,  5.188E-08,  5.139E-08,  4.991E-08,  
            4.951E-08,  4.833E-08,  4.640E-08,  4.524E-08,  4.479E-08,  
            4.304E-08,  4.228E-08,  4.251E-08,  4.130E-08,  3.984E-08,  
            3.894E-08,  3.815E-08,  3.732E-08,  3.664E-08,  3.512E-08,  
            3.463E-08,  3.503E-08,  3.218E-08,  3.253E-08,  3.107E-08,  
            2.964E-08,  2.920E-08,  2.888E-08,  2.981E-08,  2.830E-08,  
            2.750E-08,  2.580E-08,  2.528E-08,  2.444E-08,  2.378E-08,  
            2.413E-08,  2.234E-08,  2.316E-08,  2.199E-08,  2.088E-08,  
            1.998E-08,  1.920E-08,  1.942E-08,  1.859E-08,  1.954E-08,  
            1.955E-08,  1.749E-08,  1.720E-08,  1.702E-08,  1.521E-08,  
            1.589E-08,  1.469E-08,  1.471E-08,  1.543E-08,  1.433E-08,  
            1.298E-08,  1.274E-08,  1.226E-08,  1.204E-08,  1.201E-08,  
            1.298E-08,  1.220E-08,  1.220E-08,  1.096E-08,  1.080E-08,  
            9.868E-09,  9.701E-09,  1.130E-08,  9.874E-09,  9.754E-09,  
            9.651E-09,  9.725E-09,  8.413E-09,  7.705E-09,  7.846E-09,  
            8.037E-09,  9.163E-09,  8.098E-09,  8.160E-09,  7.511E-09,  
            7.011E-09,  6.281E-09,  6.502E-09,  7.323E-09,  7.569E-09,  
            5.941E-09,  5.867E-09,  5.676E-09,  4.840E-09,  5.063E-09,  
            5.207E-09,  4.917E-09,  5.033E-09,  5.356E-09,  3.795E-09,  
            4.983E-09,  4.600E-09,  3.635E-09,  3.099E-09,  2.502E-09,  
            3.823E-09,  3.464E-09,  4.332E-09,  3.612E-09,  3.682E-09,  
            3.709E-09,  3.043E-09,  3.593E-09,  3.995E-09,  4.460E-09,  
            3.583E-09,  3.290E-09,  3.132E-09,  2.812E-09,  3.109E-09,  
            3.874E-09,  3.802E-09,  4.024E-09,  3.901E-09,  2.370E-09,  
            1.821E-09,  2.519E-09,  4.701E-09,  3.855E-09,  4.685E-09,  
            5.170E-09,  4.387E-09,  4.148E-09,  4.043E-09,  3.545E-09,  
            3.392E-09,  3.609E-09,  4.635E-09,  3.467E-09,  2.558E-09,  
            3.389E-09,  2.672E-09,  2.468E-09,  1.989E-09,  2.816E-09,  
            4.023E-09,  2.664E-09,  2.219E-09,  3.169E-09,  1.654E-09,  
            3.189E-09,  2.535E-09,  2.618E-09,  3.265E-09,  2.138E-09,  
            1.822E-09,  2.920E-09,  2.002E-09,  1.300E-09,  3.764E-09,  
            3.212E-09,  3.222E-09,  2.961E-09,  2.108E-09,  1.708E-09,  
            2.636E-09,  2.937E-09,  2.939E-09,  2.732E-09,  2.218E-09,  
            1.046E-09,  6.419E-10,  1.842E-09,  1.112E-09,  1.265E-09,  
            4.087E-09,  2.044E-09,  1.022E-09,  5.109E-10,  2.554E-10,  
            1.277E-10,  6.386E-11,  0.000E+00};

// CKD_MT 1.00 implementation of oxygen v0<-v0 band model of 
//   Mate et al. over the spectral region 7550-8486 cm-1: 
//   B. Mate, C. Lugez, G.T. Fraser, W.J. Lafferty,
//   "Absolute Intensities for the O2 1.27 micron
//   continuum absorption",  
//   J. Geophys. Res., 104, 30,585-30,590, 1999. 
//
// The units of these continua coefficients are  1 / (amagat_O2*amagat_air)
//
// Also, refer to the paper "Observed  Atmospheric
// Collision Induced Absorption in Near Infrared Oxygen Bands",
// Mlawer, Clough, Brown, Stephen, Landry, Goldman, & Murcray,
// Journal of Geophysical Research (1998).
//
// (comment: v1<--v0 band around 9100 to 11000 cm-1 (1mu))
const Numeric O2_10_ckd_mt_100_v1  =  9100.000e0;
const Numeric O2_10_ckd_mt_100_v2  = 11000.000e0;
const Numeric O2_10_ckd_mt_100_dv  =     2.000e0;

#endif // continua_h
