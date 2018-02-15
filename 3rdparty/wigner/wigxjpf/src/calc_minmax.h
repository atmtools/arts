
/* Copyright 2015 Haakan T. Johansson */

/*  This file is part of WIGXJPF.
 *
 *  WIGXJPF is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  WIGXJPF is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with WIGXJPF.  If not, see
 *  <http://www.gnu.org/licenses/>.
 */

#ifndef __WIGXJPF_CALC_MINMAX_H__
#define __WIGXJPF_CALC_MINMAX_H__

#define CHOOSE_MIN(a,b) ((a) < (b) ? (a) : (b))
#define CHOOSE_MAX(a,b) ((a) > (b) ? (a) : (b))

static inline int min3(int a, int b, int c)
{
  int ab = CHOOSE_MIN(a,b);
  return CHOOSE_MIN(ab,c);
}

static inline int max3(int a, int b, int c)
{
  int ab = CHOOSE_MAX(a,b);
  return CHOOSE_MAX(ab,c);
}

static inline int max4(int a, int b, int c, int d)
{
  int ab = CHOOSE_MAX(a,b);
  int cd = CHOOSE_MAX(c,d);
  return CHOOSE_MAX(ab,cd);
}

#endif/*__WIGXJPF_CALC_MINMAX_H__*/
