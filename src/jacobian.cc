/* Copyright (C) 2004 Mattias Ekström  <ekstrom@rss.chalmers.se>

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

/*!
  \file   jacobian.cc
  \author Mattias Ekström <ekstrom@rss.chalmers.se>
  \date   2004-09-14

  \brief  Routines for setting up the jacobian.
*/

#include "arts.h"
#include "jacobian.h"

ostream& operator << (ostream& os, const RetrievalQuantity& ot)
{
  return os << ot.MainTag() << " " << ot.Subtag();
}
