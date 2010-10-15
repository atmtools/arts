/* Copyright (C) 2010 Oliver Lemke <olemke@core-dump.info>

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
  \file   docserver.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2010-09-21

  \brief  Declarations for the arts documentation server.
*/

#ifndef docserver_h
#define docserver_h

#include "arts.h"

int docserver_start(Index port, bool daemon);

#endif /* docserver_h */

