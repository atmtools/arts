/* Copyright (C) 2004 Oliver Lemke
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */

#include <cstdlib>
#include <iostream>
#include "bifstream.h"
#include "exceptions.h"

int main (int argc, char *argv[])
{
  if (argc != 2)
    {
      cerr << "Usage: " << argv[0] << " PPFILE" << endl;
      exit (EXIT_FAILURE);
    }


  bifstream is (argv[1]);
  if (is.good())
    {
      long i;
      double d;

      is.setFlag (binio::BigEndian, true);
      cout << "Reading ppheader from " << argv[1] << endl;
      for (int j = 0; j < 65 && is.good (); j++)
        {
          i = is.readInt (4);
          cout << j << " " << i << " " << is.error () << endl;
        }
      cout << endl << endl;
      for (int j = 0; j < 300 && is.good (); j++)
        {
          d = is.readFloat (binio::Single);
          cout << j << " " << d << " " << is.error () << endl;
        }
    }
  else
    {
      cerr << "Error reading from " << argv[1] << endl;
    }

  return (EXIT_SUCCESS);
}

