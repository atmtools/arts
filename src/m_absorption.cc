/*!
  \file   m_absorption.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Wed Nov 20 18:04:20 2002
  
  \brief  Methods related to absorption, lookup table, etc.
*/

#include "arts.h"
#include "messages.h"
#include "gas_abs_lookup.h"

//! Creates an empty gas absorption lookup table.
/*! 
  This is mainly there to help developers. For example, you can write
  the empty table to an XML file, to see the file format.

  \param GasAbsLookup Absorption lookup table.
*/
void gas_abs_lookupInit(GasAbsLookup& x)
{
  // Nothing to do here.
  // That means, we rely on the default constructor.

  out2 << "  Created an empty gas absorption lookup table.\n";
}
