
/* Copyright 2015 Haakan T. Johansson */

/*  This file is part of FASTWIGXJ.
 *
 *  FASTWIGXJ is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  FASTWIGXJ is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with FASTWIGXJ.  If not, see
 *  <http://www.gnu.org/licenses/>.
 */

void sort_deduplicate()
{
  if (!_nbuf)
    return;

  std::sort (_buf, _buf+_nbuf, compare<uint64_t>);
  
  // Remove duplicates
  
  size_t j = 0;
  
  uint64_t last = _buf[0];
  
  for (size_t i = 1; i < _nbuf; i++)
    {
      if (_buf[i] != last)
	j++;
      last = _buf[j] = _buf[i];
    }
  j++;
  /*
    fprintf (stderr,"%zd %zd %016x %016x\n",
    j,_nbuf,_buf[j-1],_buf[j-2]);
  */
  _nbuf = j;
}

