/* Copyright (C) 2013 Oliver Lemke
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*!
  \file   tmatrix.h
  \author Oliver Lemke
  \date   2013-06-25

  \brief  Declarations for the T-Matrix interface.
*/

#ifndef tmatrix_h
#define tmatrix_h

/** T-Matrix validation test.

 Executes the standard test included with the double precision T-Matrix code
 for randomly oriented nonspherical particles.
 Should give the same as running the tmatrix_tmd executable in
 3rdparty/tmatrix/.

 \author Oliver Lemke
 */
void tmatrix_tmd_test();

#endif //  tmatrix_h

