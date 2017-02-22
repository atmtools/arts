/* Copyright (C) 2000-2012
   Stefan Buehler <sbuehler@ltu.se>
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Oliver Lemke <olemke@core-dump.info>

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
  \file  arts.h

  The global header file for ARTS. This file is included directly or
  indirectly by each and every ARTS source file. It must therefor not
  contain stuff that should not always be present.

  Note that you do not have to include this file explicitly in many
  cases, since it is included directly or indirectly by most ARTS
  header files.

  \author Stefan Buehler
  \date 16.05.1999 */

/** \mainpage

    <center><img src="arts-splash.png" alt="ARTS"></center>

    <h2>What is ARTS?</h2>

    ARTS-2: 3D version for polarized radiative transfer calculations
    including particle scattering. It includes almost all features of the first
    version and also a number of additional functions.

    One of the main new features is the implementation of particle
    scattering as for many applications, scattering of microwave radiation
    by ice particles in the atmosphere emerges as an important issue. The
    new model should be able to simulate realistic cloud cases for microwave
    measurements in limb sounding geometry. Modeling radiative transfer
    through clouds is a complicated topic for various reasons.

    The cloud coverage is vertically and horizontally strongly inhomogeneous
    which implies that a 3D model is unavoidable for simulating realistic
    cases. Especially for simulating limb measurements, a 3D geometry is
    required as the observed region in the atmosphere has a horizontally
    large extent. Clouds consist of a variety of hydrometeors. There
    are liquid water clouds but also cirrus clouds which consist of ice
    particles of different sizes and shapes. Particle scattering leads to
    polarization effects, therefore modeling only the first component of
    the Stokes vector, the scalar intensity, is not sufficient. At least
    the first two components are required, in some cases, depending on the
    formation of the cloud, even all four components.

    The VRTE (Vector Radiative Transfer Equation) is an inhomogeneous vector
    differential equation for the Stokes vector. This equation can be
    solved numerically using an iterative method. So far gaseous absorption
    is pre-calculated using the first version of ARTS and stored in a
    lookup table. The single scattering properties of the particles, i.e.
    extinction, absorption and scattering, are calculated using the T-matrix
    method and stored in a data base.

    The Zeeman effect is also currently being implemented.

    As the program is modular the user can adjust the control file according
    to his/her requirements. The atmospheric dimensionality can be chosen to
    be 1D, 2D or 3D. If clearsky calculations are performed without Zeeman
    effect it does not make sense to calculate all 4 Stokes components. And
    for special symmetries in the scattering region, which can be switched
    on or off, the 3rd and 4th component of the Stokes vector are negligible
    small. Thus the user also can also decide how many Stokes components
    shall be simulated.

    <h2>Documentation</h2>

    You can use the HTML version to browse the source text. Just point and
    click, and eventually you will see the real implementation of functions
    and classes.

    If you are looking for a more comprehensive text, check out the
    Arts User Guide that is also distributed along with the
    program. Section `Documentation' in Chapter `The art of developing
    ARTS' there also tells you how you should add documentation
    headers to your code if you are an ARTS developer.
 */

#ifndef arts_h
#define arts_h

#include <cstddef>
#include <cstdlib>
#include "matpack.h"
#include "mystring.h"
#include "debug.h"

//----------< First of all, include the configuration header >----------
#include "config.h"

// C Assert macro:
#include <cassert>

#ifdef HAVE_NAMESPACES
  // We need those to support ansi-compliant compilers (gcc-3x)
  namespace std {}
  using namespace std;
#endif

//---------------< Global variable declarations >---------------
// See global_data.h

//---------------< Global function declarations: >---------------
// Documentations are with function definitions.
// FIXME: OLE: These should be moved to a separate header file.
class ArtsOut;

void define_wsv_group_names();  
Index get_wsv_id(const String& name);
Index get_wsv_id(const char *name);
bool is_valid_keyword_group(const Index name);
void define_species_data();
void define_species_map();
void define_lineshape_data();
void define_lineshape_norm_data();

void arts_exit (int status = EXIT_FAILURE);
void arts_exit_with_error_message(const String& m, ArtsOut &os);

//
// Physical constants are now in constants.cc
//

//---------------< Global macro definitions: >---------------

#endif // arts_h

