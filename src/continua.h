/* Copyright (C) 2001-2012
   Thomas Kuhn    <tkuhn@uni-bremen.de>
   Stefan Buehler <sbuehler@ltu.se>

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
#include "mystring.h"
#include "messages.h"


//////////////////////////////////////////////////////////////////////////// 
// entry function to all continua and full model functions
//////////////////////////////////////////////////////////////////////////// 

void xsec_continuum_tag (MatrixView         xsec,       // calculated x-section
                         const String&      name,       // model name
                         ConstVectorView    parameters, // model 
                         const String&      model,      // model option
                         ConstVectorView    f_grid,     // frequency vector
                         ConstVectorView    abs_p,      // pressure vector
                         ConstVectorView    abs_t,      // temperature vector 
                         ConstVectorView    abs_n2,     // N2 vmr profile
                         ConstVectorView    abs_h2o,    // H2O vmr profile
                         ConstVectorView    abs_o2,    // H2O vmr profile
                         ConstVectorView    vmr,        // species vmr profile
                         const Verbosity& verbosity);

//////////////////////////////////////////////////////////////////////////// 
// check of consistency of all full and continua absorption models
//////////////////////////////////////////////////////////////////////////// 

void check_continuum_model(const String& name);


//////////////////////////////////////////////////////////////////////////// 
// water vapor line+continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void MPM87H2OAbsModel (MatrixView        xsec,       // calculated x-section
                       const Numeric     CC,         // continuum scale factor 
                       const Numeric     CL,         // line strength scale factor
                       const Numeric     CW,         // line broadening scale factor
                       const String&     model,      // model option
                       ConstVectorView   f_grid,     // frequency vector
                       ConstVectorView   abs_p,      // pressure vector
                       ConstVectorView   abs_t,      // temperature vector
                       ConstVectorView   vmr,        // H2O vmr profile
                       const Verbosity& verbosity);

void MPM89H2OAbsModel (MatrixView        xsec,       // calculated x-section
                       const Numeric     CCin,       // continuum scale factor 
                       const Numeric     CLin,       // line strength scale factor
                       const Numeric     CWin,       // line broadening scale factor
                       const String&     model,      // model option
                       ConstVectorView   f_grid,     // frequency vector
                       ConstVectorView   abs_p,      // pressure vector
                       ConstVectorView   abs_t,      // temperature vector
                       ConstVectorView   vmr,        // H2O vmr profile
                       const Verbosity& verbosity);

void MPM93H2OAbsModel (MatrixView        xsec,
                       const Numeric     CCin,       // continuum scale factor 
                       const Numeric     CLin,       // line strength scale factor
                       const Numeric     CWin,       // line broadening scale factor
                       const String&     model,      // model option
                       ConstVectorView   f_grid,     // frequency vector
                       ConstVectorView   abs_p,      // pressure vector
                       ConstVectorView   abs_t,      // temperature vector
                       ConstVectorView   vmr,
                       const Verbosity& verbosity);      // H2O vmr profile

void PWR98H2OAbsModel (MatrixView        xsec,       // calculated x-section
                       const Numeric     CCin,       // continuum scale factor 
                       const Numeric     CLin,       // line strength scale factor
                       const Numeric     CWin,       // line broadening scale factor
                       const String&     model,      // model option
                       ConstVectorView   f_grid,     // frequency vector
                       ConstVectorView   abs_p,      // pressure vector
                       ConstVectorView   abs_t,      // temperature vector
                       ConstVectorView   vmr,        // H2O vmr profile
                       const Verbosity& verbosity);

void CP98H2OAbsModel (MatrixView        xsec,        // calculated x-section
                      const Numeric     CCin,        // continuum scale factor 
                      const Numeric     CLin,        // line strength scale factor
                      const Numeric     CWin,        // line broadening scale factor
                      const String&     model,       // model option
                      ConstVectorView   f_grid,      // frequency vector
                      ConstVectorView   abs_p,       // pressure vector
                      ConstVectorView   abs_t,       // temperature vector
                      ConstVectorView   vmr,         // H2O vmr profile
                      const Verbosity& verbosity);

//////////////////////////////////////////////////////////////////////////// 
// water vapor continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void Pardo_ATM_H2O_ForeignContinuum (MatrixView          xsec,   // calculated x-section
                                     const Numeric       Cin,    // model parameter
                                     const String&       model,  // model option
                                     ConstVectorView     f_grid, // frequency vector
                                     ConstVectorView     abs_p,  // pressure vector
                                     ConstVectorView     abs_t,  // temperature vector 
                                     ConstVectorView     vmr,    // H2O vmr profile
                                     const Verbosity& verbosity);

void Standard_H2O_self_continuum (MatrixView        xsec,        // calculated x-section
                                  const Numeric     C,           // model parameter
                                  const Numeric     x,           // model parameter
                                  const String&     model,       // model option
                                  ConstVectorView   f_grid,      // frequency vector
                                  ConstVectorView   abs_p,       // pressure vector
                                  ConstVectorView   abs_t,       // temperature vector 
                                  ConstVectorView   vmr,         // H2O vmr profile
                                  const Verbosity& verbosity);

void Standard_H2O_foreign_continuum (MatrixView        xsec,     // calculated x-section
                                     const Numeric     C,        // model parameter
                                     const Numeric     x,        // model parameter
                                     const String&     model,    // model option
                                     ConstVectorView   f_grid,   // frequency vector
                                     ConstVectorView   abs_p,    // pressure vector
                                     ConstVectorView   abs_t,    // temperature vector 
                                     ConstVectorView   vmr,      // H2O vmr profile
                                     const Verbosity& verbosity);

void MaTipping_H2O_foreign_continuum (MatrixView        xsec,     // calculated x-section
                                      const Numeric     C,        // model parameter
                                      const Numeric     x,        // model parameter
                                      const String&     model,    // model option
                                      ConstVectorView   f_grid,   // frequency vector
                                      ConstVectorView   abs_p,    // pressure vector
                                      ConstVectorView   abs_t,    // temperature vector 
                                      ConstVectorView   vmr,      // H2O vmr profile
                                      const Verbosity& verbosity);

void MPM93_H2O_continuum (MatrixView        xsec,                // calculated x-section
                          const Numeric     fcenter,             // model parameter
                          const Numeric     b1,                  // model parameter
                          const Numeric     b2,                  // model parameter
                          const Numeric     b3,                  // model parameter
                          const Numeric     b4,                  // model parameter
                          const Numeric     b5,                  // model parameter
                          const Numeric     b6,                  // model parameter
                          const String&     model,               // model option
                          ConstVectorView   f_grid,              // frequency vector
                          ConstVectorView   abs_p,               // pressure vector
                          ConstVectorView   abs_t,               // temperature vector
                          ConstVectorView   vmr,                 // H2O vmr profile
                          const Verbosity& verbosity);


void CKD_222_self_h2o (MatrixView          xsec,
                       const Numeric       Cin,
                       const String&       model,
                       ConstVectorView     f_grid,
                       ConstVectorView     abs_p,
                       ConstVectorView     abs_t,
                       ConstVectorView     vmr,
                       const Verbosity& verbosity);

void CKD_222_foreign_h2o (MatrixView          xsec,
                          const Numeric       Cin,
                          const String&       model,
                          ConstVectorView     f_grid,
                          ConstVectorView     abs_p,
                          ConstVectorView     abs_t,
                          ConstVectorView     vmr,
                          const Verbosity& verbosity);

void CKD_242_self_h2o (MatrixView          xsec,
                       const Numeric       Cin,
                       const String&       model,
                       ConstVectorView     f_grid,
                       ConstVectorView     abs_p,
                       ConstVectorView     abs_t,
                       ConstVectorView     vmr,
                       const Verbosity& verbosity);

void CKD24_H20 (MatrixView          xsec,      // calculated x-section
                int                 isf,       // flag if self or foreign cont.
                const Numeric       Cin,       // model scaling factor
                const String&       model,     // model option
                ConstVectorView     f_grid,    // frequency vector
                ConstVectorView     abs_p,     // pressure vector
                ConstVectorView     abs_t,     // temperature vector
                ConstVectorView     vmr,       // H2O vmr profile
                ConstVectorView     abs_n2,    // N2 vmr profile
                const Verbosity& verbosity);

void CKD_242_foreign_h2o (MatrixView          xsec,
                          const Numeric       Cin,
                          const String&       model,
                          ConstVectorView     f_grid,
                          ConstVectorView     abs_p,
                          ConstVectorView     abs_t,
                          ConstVectorView     vmr,
                          const Verbosity& verbosity);

void CKD_mt_100_self_h2o (MatrixView          xsec,
                          const Numeric       Cin,
                          const String&       model,
                          ConstVectorView     f_grid,
                          ConstVectorView     abs_p,
                          ConstVectorView     abs_t,
                          ConstVectorView     vmr,
                          const Verbosity& verbosity);

void CKD_mt_100_foreign_h2o (MatrixView          xsec,
                             const Numeric       Cin,
                             const String&       model,
                             ConstVectorView     f_grid,
                             ConstVectorView     abs_p,
                             ConstVectorView     abs_t,
                             ConstVectorView     vmr,
                             const Verbosity& verbosity);

void CKD_mt_250_self_h2o (MatrixView          xsec,
                          const Numeric       Cin,
                          const String&       model,
                          ConstVectorView     f_grid,
                          ConstVectorView     abs_p,
                          ConstVectorView     abs_t,
                          ConstVectorView     vmr,
                          const Verbosity& verbosity);


void CKD_mt_250_foreign_h2o (MatrixView          xsec,
                             const Numeric       Cin,
                             const String&       model,
                             ConstVectorView     f_grid,
                             ConstVectorView     abs_p,
                             ConstVectorView     abs_t,
                             ConstVectorView     vmr,
                             const Verbosity& verbosity);

//////////////////////////////////////////////////////////////////////////// 
// oxygen line+continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void MPM85O2AbsModel (MatrixView        xsec,        // calculated x-section
                      const Numeric     CC,          // model parameter
                      const Numeric     CL,          // model parameter
                      const Numeric     CW,          // model parameter
                      const Numeric     CO,          // model parameter
                      const String&     model,       // model option
                      ConstVectorView   f_grid,      // frequency vector
                      ConstVectorView   abs_p,       // pressure vector
                      ConstVectorView   abs_t,       // temperature vector
                      ConstVectorView   abs_h2o,     // H2O vmr profile
                      ConstVectorView   vmr,         // O2 vmr profile
                      const Verbosity& verbosity);


void MPM87O2AbsModel (MatrixView        xsec,        // calculated x-section
                      const Numeric     CC,          // model parameter
                      const Numeric     CL,          // model parameter
                      const Numeric     CW,          // model parameter
                      const Numeric     CO,          // model parameter
                      const String&     model,       // model option
                      ConstVectorView   f_grid,      // frequency vector
                      ConstVectorView   abs_p,       // pressure vector
                      ConstVectorView   abs_t,       // temperature vector
                      ConstVectorView   abs_h2o,     // H2O vmr profile
                      ConstVectorView   vmr,         // O2 vmr profile
                      const Verbosity& verbosity);


void MPM89O2AbsModel (MatrixView        xsec,        // calculated x-section
                      const Numeric     CC,          // model parameter
                      const Numeric     CL,          // model parameter
                      const Numeric     CW,          // model parameter
                      const Numeric     CO,          // model parameter
                      const String&     model,       // model option
                      ConstVectorView   f_grid,      // frequency vector
                      ConstVectorView   abs_p,       // pressure vector
                      ConstVectorView   abs_t,       // temperature vector
                      ConstVectorView   abs_h2o,     // H2O vmr profile
                      ConstVectorView   vmr,         // O2 vmr profile
                      const Verbosity& verbosity);


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


void MPM93O2AbsModel (MatrixView        xsec,        // calculated x-section
                      const Numeric     CC,          // model parameter
                      const Numeric     CL,          // model parameter
                      const Numeric     CW,          // model parameter
                      const Numeric     CO,          // model parameter
                      const String&     model,       // model option
                      ConstVectorView   f_grid,      // frequency vector
                      ConstVectorView   abs_p,       // pressure vector
                      ConstVectorView   abs_t,       // temperature vector
                      ConstVectorView   abs_h2o,     // H2O vmr profile
                      ConstVectorView   vmr,         // O2 vmr profile
                      const Verbosity& verbosity);

void TRE05O2AbsModel (MatrixView        xsec,        // calculated x-section
                      const Numeric     CC,          // model parameter
                      const Numeric     CL,          // model parameter
                      const Numeric     CW,          // model parameter
                      const Numeric     CO,          // model parameter
                      const String&     model,       // model option
                      ConstVectorView   f_grid,      // frequency vector
                      ConstVectorView   abs_p,       // pressure vector
                      ConstVectorView   abs_t,       // temperature vector
                      ConstVectorView   abs_h2o,     // H2O vmr profile
                      ConstVectorView   vmr,         // O2 vmr profile
                      const Verbosity& verbosity);

void PWR93O2AbsModel (MatrixView        xsec,        // calculated x-section
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
                      ConstVectorView   vmr,         // O2 vmr profile
                      const Verbosity& verbosity);

//////////////////////////////////////////////////////////////////////////// 
// oxygen continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void MPM93_O2_continuum (MatrixView        xsec,             // calculated x-section
                         const Numeric     S0in,             // model parameter
                         const Numeric     G0in,             // model parameter
                         const Numeric     XSOin,            // model parameter
                         const Numeric     XG0in,            // model parameter
                         const String&     model,            // model option
                         ConstVectorView   f_grid,           // frequency vector
                         ConstVectorView   abs_p,            // pressure vector
                         ConstVectorView   abs_t,            // temperature vector
                         ConstVectorView   abs_h2o,          // H2O vmr profile
                         ConstVectorView   vmr,              // O2 vmr profile
                         const Verbosity& verbosity);

void Rosenkranz_O2_continuum (MatrixView        xsec,        // calculated x-section
                              const Numeric     S0in,        // model parameter
                              const Numeric     G0in,        // model parameter
                              const Numeric     XSOin,       // model parameter
                              const Numeric     XG0in,       // model parameter
                              const String&     model,       // model option
                              ConstVectorView   f_grid,      // frequency vector
                              ConstVectorView   abs_p,       // pressure vector
                              ConstVectorView   abs_t,       // temperature vector
                              ConstVectorView   abs_h2o,     // H2O vmr profile
                              ConstVectorView   vmr,         // O2 vmr profile
                              const Verbosity& verbosity);

void CKD_mt_CIAfun_o2 (MatrixView          xsec,        // calculated x-section
                       const Numeric       Cin,         // scaling factor
                       const String&       model,       // model option
                       ConstVectorView     f_grid,      // frequency vector
                       ConstVectorView     abs_p,       // pressure vector
                       ConstVectorView     abs_t,       // temperature vector
                       ConstVectorView     vmr,         // O2 vmr profile
                       const Verbosity& verbosity);

void CKD_mt_v0v0_o2 (MatrixView          xsec,        // calculated x-section
                     const Numeric       Cin,         // scaling factor
                     const String&       model,       // model option
                     ConstVectorView     f_grid,      // frequency vector
                     ConstVectorView     abs_p,       // pressure vector
                     ConstVectorView     abs_t,       // temperature vector
                     ConstVectorView     vmr,         // O2 vmr profile
                     ConstVectorView     abs_n2,      // N2 vmr profile
                     const Verbosity& verbosity);

void CKD_mt_v1v0_o2 (MatrixView          xsec,        // calculated x-section
                     const Numeric       Cin,         // scaling factor
                     const String&       model,       // model option
                     ConstVectorView     f_grid,      // frequency vector
                     ConstVectorView     abs_p,       // pressure vector
                     ConstVectorView     abs_t,       // temperature vector
                     ConstVectorView     vmr,         // O2 vmr profile
                     const Verbosity& verbosity);

void CKD_mt_250_o2_vis (MatrixView          xsec,        // calculated x-section
                     const Numeric       Cin,         // scaling factor
                     const String&       model,       // model option
                     ConstVectorView     f_grid,      // frequency vector
                     ConstVectorView     abs_p,       // pressure vector
                     ConstVectorView     abs_t,       // temperature vector
                     ConstVectorView     vmr,         // O2 vmr profile
                     const Verbosity& verbosity);


//////////////////////////////////////////////////////////////////////////// 
// nitrogen continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void CKD_mt_CIArot_n2 (MatrixView         xsec,        // calculated x-section
                       const Numeric      Cin,         // scaling factor
                       const String&      model,       // model option
                       ConstVectorView    f_grid,      // frequency vector
                       ConstVectorView    abs_p,       // pressure vector
                       ConstVectorView    abs_t,       // temperature vector
                       ConstVectorView    vmr,         // N2 vmr profile
                       const Verbosity& verbosity);

void CKD_mt_CIAfun_n2 (MatrixView          xsec,        // calculated x-section
                       const Numeric       Cin,         // scaling factor
                       const String&       model,       // model option
                       ConstVectorView     f_grid,      // frequency vector
                       ConstVectorView     abs_p,       // pressure vector
                       ConstVectorView     abs_t,       // temperature vector
                       ConstVectorView     vmr,         // N2 vmr profile
                       const Verbosity& verbosity);

void CKD_mt_250_CIArot_n2 (MatrixView         xsec,        // calculated x-section
                       const Numeric      Cin,         // scaling factor
                       const String&      model,       // model option
                       ConstVectorView    f_grid,      // frequency vector
                       ConstVectorView    abs_p,       // pressure vector
                       ConstVectorView    abs_t,       // temperature vector
                       ConstVectorView    vmr,         // N2 vmr profile
                       ConstVectorView    abs_h2o,     // H2O vmr profile
                       ConstVectorView    abs_o2,      // O2 vmr profile
                       const Verbosity& verbosity);


void CKD_mt_250_CIAfun_n2 (MatrixView          xsec,        // calculated x-section
                       const Numeric       Cin,         // scaling factor
                       const String&       model,       // model option
                       ConstVectorView     f_grid,      // frequency vector
                       ConstVectorView     abs_p,       // pressure vector
                       ConstVectorView     abs_t,       // temperature vector
                       ConstVectorView     vmr,         // N2 vmr profile
                       ConstVectorView     abs_h2o,     // H2O vmr profile
                       ConstVectorView     abs_o2,      // O2 vmr profile
                       const Verbosity& verbosity);


void BF86_CIA_N2 (MatrixView              xsec,        // calculated x-section
                  const Numeric           Cin,         // model parameter
                  const String&           model,       // model option 
                  ConstVectorView         f_grid,      // frequency vector
                  ConstVectorView         abs_p,       // pressure vector
                  ConstVectorView         abs_t,       // temperature vector
                  ConstVectorView         vmr,         // N2 vmr profile
                  const Verbosity& verbosity);

void MPM93_N2_continuum (MatrixView       xsec,        // calculated x-section
                         const Numeric    Cin,         // model parameter
                         const Numeric    Gin,         // model parameter
                         const Numeric    xTin,        // model parameter
                         const Numeric    xfin,        // model parameter
                         const String&    model,       // model option
                         ConstVectorView  f_grid,      // frequency vector
                         ConstVectorView  abs_p,       // pressure vector
                         ConstVectorView  abs_t,       // temperature vector
                         ConstVectorView  abs_h2o,     // H2O vmr profile
                         ConstVectorView  vmr,         // N2 vmr profile
                         const Verbosity& verbosity);

void Rosenkranz_N2_self_continuum (MatrixView        xsec,        // calculated x-section
                                   const Numeric     Cin,         // model parameter
                                   const Numeric     xin,         // model parameter
                                   const String&     model,       // model option
                                   ConstVectorView   f_grid,      // frequency vector
                                   ConstVectorView   abs_p,       // pressure vector
                                   ConstVectorView   abs_t,       // temperature vector
                                   ConstVectorView   vmr,         // N2 vmr profile
                                   const Verbosity& verbosity);

void Standard_N2_self_continuum (MatrixView        xsec,        // calculated x-section
                                 const Numeric     Cin,         // model parameter
                                 const Numeric     xfin,        // model parameter
                                 const Numeric     xtin,        // model parameter
                                 const Numeric     xpin,        // model parameter
                                 const String&     model,       // model option
                                 ConstVectorView   f_grid,      // frequency vector
                                 ConstVectorView   abs_p,       // pressure vector
                                 ConstVectorView   abs_t,       // temperature vector
                                 ConstVectorView   vmr,         // N2 vmr profile
                                 const Verbosity& verbosity);

void Pardo_ATM_N2_dry_continuum (MatrixView        xsec,         // calculated x-section
                                 const Numeric     Cin,          // model parameter
                                 const String&     model,        // model option
                                 ConstVectorView   f_grid,       // frequency vector
                                 ConstVectorView   abs_p,        // pressure vector
                                 ConstVectorView   abs_t,        // temperature vector
                                 ConstVectorView   vmr,          // N2 vmr profile
                                 ConstVectorView   abs_h2o,      // H2O vmr profile
                                 const Verbosity& verbosity);

//////////////////////////////////////////////////////////////////////////// 
// carbon dioxide continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void CKD_241_co2 (MatrixView          xsec,        // calculated x-section
                  const Numeric       Cin,         // scaling factor
                  const String&       model,       // model option
                  ConstVectorView     f_grid,      // frequency vector
                  ConstVectorView     abs_p,       // pressure vector
                  ConstVectorView     abs_t,       // temperature vector
                  ConstVectorView     vmr,         // CO2 vmr profile
                  const Verbosity& verbosity);

void CKD_mt_co2 (MatrixView          xsec,         // calculated x-section
                 const Numeric       Cin,          // scaling factor
                 const String&       model,        // model option
                 ConstVectorView     f_grid,       // frequency vector
                 ConstVectorView     abs_p,        // pressure vector
                 ConstVectorView     abs_t,        // temperature vector
                 ConstVectorView     vmr,          // CO2 vmr profile
                 const Verbosity& verbosity);

void CKD_250_mt_co2 (MatrixView          xsec,         // calculated x-section
                 const Numeric       Cin,          // scaling factor
                 const String&       model,        // model option
                 ConstVectorView     f_grid,       // frequency vector
                 ConstVectorView     abs_p,        // pressure vector
                 ConstVectorView     abs_t,        // temperature vector
                 ConstVectorView     vmr,          // CO2 vmr profile
                 const Verbosity& verbosity);

void Rosenkranz_CO2_self_continuum (MatrixView        xsec,       // calculated x-section
                                    const Numeric     C,          // model parameter
                                    const Numeric     x,          // model parameter
                                    const String&     model,      // model option
                                    ConstVectorView   f_grid,     // frequency vector
                                    ConstVectorView   abs_p,      // pressure vector
                                    ConstVectorView   abs_t,      // temperature vector
                                    ConstVectorView   vmr,        // CO2 vmr profile
                                    const Verbosity& verbosity);

void Rosenkranz_CO2_foreign_continuum (MatrixView        xsec,    // calculated x-section
                                       const Numeric     C,       // model parameter
                                       const Numeric     x,       // model parameter
                                       const String&     model,   // model option
                                       ConstVectorView   f_grid,  // frequency vector
                                       ConstVectorView   abs_p,   // pressure vector
                                       ConstVectorView   abs_t,   // temperature vector
                                       ConstVectorView   abs_n2,  // N2 vmr profile
                                       ConstVectorView   vmr,     // CO2 vmr profile
                                       const Verbosity& verbosity);

//////////////////////////////////////////////////////////////////////////// 
// water droplet and ice particle absorption (clouds)
//////////////////////////////////////////////////////////////////////////// 

void MPM93WaterDropletAbs (MatrixView        xsec,     // calculated x-section
                           const Numeric     CC,       // model parameter
                           const Numeric     CG,       // model parameter
                           const Numeric     CE,       // model parameter
                           const String&     model,    // model option
                           ConstVectorView   f_grid,   // frequency vector
                           ConstVectorView   abs_p,    // pressure vector
                           ConstVectorView   abs_t,    // temperature vector
                           ConstVectorView   vmr,      // suspended water droplet density vector
                           const Verbosity& verbosity);

void MPM93IceCrystalAbs (MatrixView        xsec,       // calculated x-section
                         const Numeric     CC,         // model parameter
                         const Numeric     CA,         // model parameter
                         const Numeric     CB,         // model parameter
                         const String&     model,      // model option
                         ConstVectorView   f_grid,     // frequency vector
                         ConstVectorView   abs_p,      // pressure vector
                         ConstVectorView   abs_t,      // temperature vector
                         ConstVectorView   vmr,        // suspended ice particle density vector,
                         const Verbosity& verbosity);

void MPM93RainExt (MatrixView        xsec,       // calculated x-section
                   const Numeric     CE,         // model parameter
                   const Numeric     CA,         // model parameter
                   const Numeric     CB,         // model parameter
                   const String&     model,      // model option
                   ConstVectorView   f_grid,     // frequency vector
                   ConstVectorView   abs_p,      // pressure vector
                   ConstVectorView   abs_t,      // temperature vector
                   ConstVectorView   vmr,        // rain rate vector, 
                   const Verbosity& verbosity);

void ELL07WaterDropletAbs (MatrixView       xsec,     // calculatd x-section
                           const String&    model,    // model option
                           ConstVectorView  f_grid,   // frequency vector
                           ConstVectorView  abs_p,    // pressure vector
                           ConstVectorView  abs_t,    // temperature vector
                           ConstVectorView  vmr,      // suspended water droplet density vector
                           const Verbosity& verbosity);

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
