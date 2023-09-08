/*!
  \file   m_rt4.cc
  \author Jana Mendrok <jana.mendrok@gmail.com>
  \date   2016-05-24
  
  \brief  Workspace functions related to application of scattering solver RT4.
  
  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include <workspace.h>
#include "disort.h"
#include "m_xml.h"
#include "rt4.h"
#include "species_tags.h"


/* Workspace method: Doxygen documentation will be auto-generated */
#ifdef ENABLE_RT4
void RT4Test(Tensor4& out_rad,
             const String& datapath) {
  rt4_test(out_rad, datapath);
}
#else
void RT4Test(Tensor4&, const String&) {
  ARTS_USER_ERROR ("This version of ARTS was compiled without RT4 support.");
}
#endif
