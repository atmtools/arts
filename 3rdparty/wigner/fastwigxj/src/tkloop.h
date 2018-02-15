
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

do {
  int tk_d04 = abs(two_j[0] - two_j[4]);
  int tk_s04 =     two_j[0] + two_j[4];
  int tk_d27 = abs(two_j[2] - two_j[7]);
  int tk_s27 =     two_j[2] + two_j[7];
  int tkmin0 = tk_d04 > tk_d27 ? tk_d04 : tk_d27;
  int tkmax0 = tk_s04 < tk_s27 ? tk_s04 : tk_s27;
  int tk_d56 = abs(two_j[5] - two_j[6]);
  int tk_s56  =    two_j[5] + two_j[6];
  tkmin0 = tkmin0 > tk_d56 ? tkmin0 : tk_d56;
  tkmax0 = tkmax0 < tk_s56 ? tkmax0 : tk_s56;
  int tkdiff0 = tkmax0 - tkmin0;
  {
    tkmin  = tkmin0;
    tkmax  = tkmax0;
    tkdiff = tkdiff0;
    order  = 0;
    if (tkdiff == 0) goto order0;
  }
  int tk_d05 = abs(two_j[0] - two_j[5]);
  int tk_s05 =     two_j[0] + two_j[5];
  int tk_d18 = abs(two_j[1] - two_j[8]);
  int tk_s18 =     two_j[1] + two_j[8];
  int tkmin1 = tk_d05 > tk_d18 ? tk_d05 : tk_d18;
  int tkmax1 = tk_s05 < tk_s18 ? tk_s05 : tk_s18;
  int tk_d46 = abs(two_j[4] - two_j[6]);
  int tk_s46  =    two_j[4] + two_j[6];
  tkmin1 = tkmin1 > tk_d46 ? tkmin1 : tk_d46;
  tkmax1 = tkmax1 < tk_s46 ? tkmax1 : tk_s46;
  int tkdiff1 = tkmax1 - tkmin1;
  if (tkdiff1 < tkdiff)   {
    tkmin  = tkmin1;
    tkmax  = tkmax1;
    tkdiff = tkdiff1;
    order  = 1;
    if (tkdiff == 0) goto order1;
  }
  int tk_d07 = abs(two_j[0] - two_j[7]);
  int tk_s07 =     two_j[0] + two_j[7];
  int tk_d24 = abs(two_j[2] - two_j[4]);
  int tk_s24 =     two_j[2] + two_j[4];
  int tkmin2 = tk_d07 > tk_d24 ? tk_d07 : tk_d24;
  int tkmax2 = tk_s07 < tk_s24 ? tk_s07 : tk_s24;
  int tk_d38 = abs(two_j[3] - two_j[8]);
  int tk_s38  =    two_j[3] + two_j[8];
  tkmin2 = tkmin2 > tk_d38 ? tkmin2 : tk_d38;
  tkmax2 = tkmax2 < tk_s38 ? tkmax2 : tk_s38;
  int tkdiff2 = tkmax2 - tkmin2;
  if (tkdiff2 < tkdiff)   {
    tkmin  = tkmin2;
    tkmax  = tkmax2;
    tkdiff = tkdiff2;
    order  = 2;
    if (tkdiff == 0) goto order2;
  }
  int tk_d08 = abs(two_j[0] - two_j[8]);
  int tk_s08 =     two_j[0] + two_j[8];
  int tk_d15 = abs(two_j[1] - two_j[5]);
  int tk_s15 =     two_j[1] + two_j[5];
  int tkmin3 = tk_d08 > tk_d15 ? tk_d08 : tk_d15;
  int tkmax3 = tk_s08 < tk_s15 ? tk_s08 : tk_s15;
  int tk_d37 = abs(two_j[3] - two_j[7]);
  int tk_s37  =    two_j[3] + two_j[7];
  tkmin3 = tkmin3 > tk_d37 ? tkmin3 : tk_d37;
  tkmax3 = tkmax3 < tk_s37 ? tkmax3 : tk_s37;
  int tkdiff3 = tkmax3 - tkmin3;
  if (tkdiff3 < tkdiff)   {
    tkmin  = tkmin3;
    tkmax  = tkmax3;
    tkdiff = tkdiff3;
    order  = 3;
    if (tkdiff == 0) goto order3;
  }
  int tk_d13 = abs(two_j[1] - two_j[3]);
  int tk_s13 =     two_j[1] + two_j[3];
  int tk_d26 = abs(two_j[2] - two_j[6]);
  int tk_s26 =     two_j[2] + two_j[6];
  int tkmin4 = tk_d13 > tk_d26 ? tk_d13 : tk_d26;
  int tkmax4 = tk_s13 < tk_s26 ? tk_s13 : tk_s26;
  int tk_d57 = abs(two_j[5] - two_j[7]);
  int tk_s57  =    two_j[5] + two_j[7];
  tkmin4 = tkmin4 > tk_d57 ? tkmin4 : tk_d57;
  tkmax4 = tkmax4 < tk_s57 ? tkmax4 : tk_s57;
  int tkdiff4 = tkmax4 - tkmin4;
  if (tkdiff4 < tkdiff)   {
    tkmin  = tkmin4;
    tkmax  = tkmax4;
    tkdiff = tkdiff4;
    order  = 4;
    if (tkdiff == 0) goto order4;
  }
  int tk_d16 = abs(two_j[1] - two_j[6]);
  int tk_s16 =     two_j[1] + two_j[6];
  int tk_d23 = abs(two_j[2] - two_j[3]);
  int tk_s23 =     two_j[2] + two_j[3];
  int tkmin5 = tk_d16 > tk_d23 ? tk_d16 : tk_d23;
  int tkmax5 = tk_s16 < tk_s23 ? tk_s16 : tk_s23;
  int tk_d48 = abs(two_j[4] - two_j[8]);
  int tk_s48  =    two_j[4] + two_j[8];
  tkmin5 = tkmin5 > tk_d48 ? tkmin5 : tk_d48;
  tkmax5 = tkmax5 < tk_s48 ? tkmax5 : tk_s48;
  int tkdiff5 = tkmax5 - tkmin5;
  if (tkdiff5 < tkdiff)   {
    tkmin  = tkmin5;
    tkmax  = tkmax5;
    tkdiff = tkdiff5;
    order  = 5; goto order5;
  }
} while (0);
switch (order) {
case 0:
order0:
  js[0].j.two_j1[0] = two_j[0];
  js[0].j.two_j2[0] = two_j[4];
  js[0].j.two_j4[0] = two_j[5];
  js[0].j.two_j5[0] = two_j[6];
  js[0].j.two_j6[0] = two_j[3];
  js[0].j.two_j1[1] = two_j[2];
  js[0].j.two_j2[1] = two_j[7];
  js[0].j.two_j4[1] = two_j[6];
  js[0].j.two_j5[1] = two_j[5];
  js[0].j.two_j6[1] = two_j[8];
  js[0].j.two_j1[2] = two_j[0];
  js[0].j.two_j2[2] = two_j[4];
  js[0].j.two_j4[2] = two_j[7];
  js[0].j.two_j5[2] = two_j[2];
  js[0].j.two_j6[2] = two_j[1];
  break;
case 1:
order1:
  js[0].j.two_j1[0] = two_j[0];
  js[0].j.two_j2[0] = two_j[5];
  js[0].j.two_j4[0] = two_j[4];
  js[0].j.two_j5[0] = two_j[6];
  js[0].j.two_j6[0] = two_j[3];
  js[0].j.two_j1[1] = two_j[1];
  js[0].j.two_j2[1] = two_j[8];
  js[0].j.two_j4[1] = two_j[6];
  js[0].j.two_j5[1] = two_j[4];
  js[0].j.two_j6[1] = two_j[7];
  js[0].j.two_j1[2] = two_j[0];
  js[0].j.two_j2[2] = two_j[5];
  js[0].j.two_j4[2] = two_j[8];
  js[0].j.two_j5[2] = two_j[1];
  js[0].j.two_j6[2] = two_j[2];
  sign ^= 1;  break;
case 2:
order2:
  js[0].j.two_j1[0] = two_j[0];
  js[0].j.two_j2[0] = two_j[7];
  js[0].j.two_j4[0] = two_j[4];
  js[0].j.two_j5[0] = two_j[2];
  js[0].j.two_j6[0] = two_j[1];
  js[0].j.two_j1[1] = two_j[3];
  js[0].j.two_j2[1] = two_j[8];
  js[0].j.two_j4[1] = two_j[2];
  js[0].j.two_j5[1] = two_j[4];
  js[0].j.two_j6[1] = two_j[5];
  js[0].j.two_j1[2] = two_j[0];
  js[0].j.two_j2[2] = two_j[7];
  js[0].j.two_j4[2] = two_j[8];
  js[0].j.two_j5[2] = two_j[3];
  js[0].j.two_j6[2] = two_j[6];
  sign ^= 1;  break;
case 3:
order3:
  js[0].j.two_j1[0] = two_j[0];
  js[0].j.two_j2[0] = two_j[8];
  js[0].j.two_j4[0] = two_j[5];
  js[0].j.two_j5[0] = two_j[1];
  js[0].j.two_j6[0] = two_j[2];
  js[0].j.two_j1[1] = two_j[3];
  js[0].j.two_j2[1] = two_j[7];
  js[0].j.two_j4[1] = two_j[1];
  js[0].j.two_j5[1] = two_j[5];
  js[0].j.two_j6[1] = two_j[4];
  js[0].j.two_j1[2] = two_j[0];
  js[0].j.two_j2[2] = two_j[8];
  js[0].j.two_j4[2] = two_j[7];
  js[0].j.two_j5[2] = two_j[3];
  js[0].j.two_j6[2] = two_j[6];
  break;
case 4:
order4:
  js[0].j.two_j1[0] = two_j[1];
  js[0].j.two_j2[0] = two_j[3];
  js[0].j.two_j4[0] = two_j[5];
  js[0].j.two_j5[0] = two_j[7];
  js[0].j.two_j6[0] = two_j[4];
  js[0].j.two_j1[1] = two_j[2];
  js[0].j.two_j2[1] = two_j[6];
  js[0].j.two_j4[1] = two_j[7];
  js[0].j.two_j5[1] = two_j[5];
  js[0].j.two_j6[1] = two_j[8];
  js[0].j.two_j1[2] = two_j[1];
  js[0].j.two_j2[2] = two_j[3];
  js[0].j.two_j4[2] = two_j[6];
  js[0].j.two_j5[2] = two_j[2];
  js[0].j.two_j6[2] = two_j[0];
  sign ^= 1;  break;
case 5:
order5:
  js[0].j.two_j1[0] = two_j[1];
  js[0].j.two_j2[0] = two_j[6];
  js[0].j.two_j4[0] = two_j[3];
  js[0].j.two_j5[0] = two_j[2];
  js[0].j.two_j6[0] = two_j[0];
  js[0].j.two_j1[1] = two_j[4];
  js[0].j.two_j2[1] = two_j[8];
  js[0].j.two_j4[1] = two_j[2];
  js[0].j.two_j5[1] = two_j[3];
  js[0].j.two_j6[1] = two_j[5];
  js[0].j.two_j1[2] = two_j[1];
  js[0].j.two_j2[2] = two_j[6];
  js[0].j.two_j4[2] = two_j[8];
  js[0].j.two_j5[2] = two_j[4];
  js[0].j.two_j6[2] = two_j[7];
  break;
}
