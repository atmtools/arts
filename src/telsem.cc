/* Copyright (C) 2017
   Oliver Lemke <olemke@core-dump.info>
   Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>

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
   USA.
*/

/*!
  \file   telsem.cc

  \brief  This file contains functions to handle TELSEM 2 atlas data.
*/

#include <cmath>
#include <utility>
#include "telsem.h"

extern Numeric EARTH_RADIUS;
extern Numeric DEG2RAD;
extern Numeric PI;

////////////////////////////////////////////////////////////////////////////////
// TelsemAtlas Class
////////////////////////////////////////////////////////////////////////////////

void TelsemAtlas::read(std::istream& is)
{
    name = "ssmi_mean_emis_climato";
    nchan = 7;
    dlat = 0.25;
    is >> ndat;
    emis.resize(ndat, nchan);
    emis = NAN;
    emis_err.resize(ndat, nchan);
    emis_err = NAN;
    classes1.resize(ndat);
    classes1 = -1;
    classes2.resize(ndat);
    classes2 = -1;
    cellnums.resize(ndat);
    cellnums = -1;

    equare();

    Index cellnum;
    Index class1;
    Index class2;
    Vector ssmi(2 * nchan);
    Index ipos = -1;
    for (Index j = 0; j < ndat; j++)
    {
        is >> cellnum;
        if (is.fail())
            throw std::runtime_error("Error reading cellnum.");
        for (Index nssmi = 0; nssmi < 2 * nchan; nssmi++)
        {
            is >> ssmi[nssmi];
            if (is.fail())
                throw std::runtime_error("Error reading emissivity.");
        }

        is >> class1 >> class2;
        if (is.fail())
            throw std::runtime_error("Error reading classes.");
        if (class1 > 0 && class2 > 0 && ipos < ndat)
        {
            ipos++;
            for (Index i = 0; i < nchan; i++)
            {
                emis(ipos, i) = ssmi[i];
                emis_err(ipos, i) = std::sqrt(ssmi[nchan + i]);
            }
            cellnums[ipos] = cellnum;
            classes1[ipos] = class1;
            classes2[ipos] = class2;
        }
    }
    telsem_calc_correspondence();
}

void TelsemAtlas::equare()
{
    Index maxlat = static_cast<Index>(floor(180.0 / dlat));

    ncells.resize(maxlat);
    firstcells.resize(maxlat);

    // Total number of cells.
    Index totcel = 0.0;

    Numeric rcelat = dlat * PI / 180.0;

    Numeric hezon = EARTH_RADIUS * sin(rcelat);
    Numeric aezon = 2.0 * PI * EARTH_RADIUS * hezon;
    Numeric aecell = aezon * dlat / 360.0;

    for (Index i = 0; i < maxlat / 2; ++i) {
        Numeric xlatb = static_cast<Numeric>(i) * dlat;
        Numeric xlate = xlatb + dlat;
        Numeric rlatb = DEG2RAD * xlatb;
        Numeric rlate = DEG2RAD * xlate;
        Numeric htb = EARTH_RADIUS * sin(rlatb);
        Numeric hte = EARTH_RADIUS * sin(rlate);
        Numeric htzone = hte - htb;
        Numeric azone  = 2.0 * PI * EARTH_RADIUS * htzone;
        Numeric rcells = azone / aecell;
        Index icellr = static_cast<Index>(floor(rcells + 0.5));

        totcel += 2 * icellr;

        Index lat1 = i + maxlat / 2;
        Index lat2 = maxlat / 2 - 1 - i;
        ncells[lat1] = icellr;
        ncells[lat2] = icellr;
    }

    firstcells[0] = 0;
    for (Index i = 1; i < maxlat; ++i) {
        firstcells[i] = firstcells[i - 1] + ncells[i];
    }
}

void TelsemAtlas::telsem_calc_correspondence()
{
    correspondence.resize(660066);
    correspondence = -1;
    for (Index j = 0; j < ndat; j++)
    {
        correspondence[cellnums[j]] = j;
    }
}

Index TelsemAtlas::calc_cellnum(Numeric lat,
                                Numeric lon) const
{
    Index cellnum = 0;
    Index ilat = static_cast<Index>(floor((lat + 90.0) / dlat));
    Index ilon = static_cast<Index>(
        floor(lon / 360.0 * static_cast<Numeric>(ncells[ilat]))
        ) + 1;
    for (Index i = 0; i < ilat; ++i) {
        cellnum += ncells[i];
    }
    cellnum += ilon;
    return cellnum;
}

std::pair<Numeric, Numeric> TelsemAtlas::get_coordinates(Index cellnum)
{

    Index index_lat_max = static_cast<Index>(std::floor(180.0 / dlat));
    Index index_lat = -1;
    Index index_lon = -1;
    if (cellnum >= firstcells[index_lat_max-1]) {
        index_lat = index_lat_max;
        index_lon = cellnum - firstcells[index_lat_max];
    } else {
        for (Index i = 0; i < index_lat_max; ++i) {
            if ((cellnum >= firstcells[i]) && (cellnum < firstcells[i + 1])) {
                index_lat = i;
                index_lon = cellnum - firstcells[i];
            }
        }
    }
    Numeric lat = (static_cast<Numeric>(index_lat) - 0.5) * dlat - 90.0;
    Numeric lon = (static_cast<Numeric>(index_lon) - 0.5)
        * (360.0 / static_cast<Numeric>(ncells[index_lat]));
    return std::make_pair(lat, lon);
}

Numeric TelsemAtlas::interp_freq2(Numeric emiss19,
                                  Numeric emiss37,
                                  Numeric emiss85,
                                  Numeric f,
                                  Index class2) const
{
    Numeric a = 0.0;
    Numeric b = 0.0;
    Numeric c = 0.0;
    Numeric emiss = 0.0;
    if (f <= 19.35) {
        a = 1;
        b = 0.0;
        c = 0.0;
        emiss = emiss19;
    } else if ((19.35 < f) && (f <= 37.0)) {
        a = (37.0 - f) / (37.0 - 19.35);
        b = (f - 19.35) / (37.0 - 19.35);
        c = 0.0;
        emiss = a * emiss19 + b * emiss37;
    } else if ((f > 37.0) && (f < 85.5)) {
        a = 0;
        b = (85.5 - f) / (85.5 - 37.0);
        c = (f - 37.0) / (85.5 - 37.0);
        emiss = b * emiss37 + c * emiss85;
    } else if (85.5 <= f) {
        a = 0;
        b = 0;
        c = 1.0;
        emiss = emiss85;
        if ((class2 > 9) && (class2 < 14) && (emiss85 > emiss37)) {
            if (f <= 150.0) {
                emiss = emiss85 + (f - 85.5) * (emiss85 - emiss37) / (85.5 - 37.0)
                    * rapport43_32[class2 - 10];
            } else if ((f > 150.0) && (f <= 190.0)) {
                emiss = emiss85 + (150.0 - 85.5) * (emiss85 - emiss37) / (85.5 - 37.0)
                    * rapport43_32[class2 - 10];
                emiss += (f - 150.0) * (emiss - emiss85) / (150.0 - 85.5)
                    * rapport54_43[class2 - 10];
            } else if (f > 190.0) {
                emiss = emiss85 + (150.0 - 85.5) * (emiss85 - emiss37) / (85.5 - 37.0)
                    * rapport43_32[class2 - 10];
                emiss += (190.0 - 150.0) * (emiss - emiss85) / (150.0 - 85.5)
                    * rapport54_43[class2 - 10];
            }
            if (emiss > 1.0) {
                emiss = 1.0;
            }
        }
    }
    return emiss;
}

std::pair<Numeric, Numeric> TelsemAtlas::emis_interp(Numeric theta,
                                                     Numeric freq,
                                                     Index class1,
                                                     Index class2,
                                                     const ConstVectorView &ev,
                                                     const ConstVectorView &eh) const
{
    Vector emiss_scal_h(3);
    Vector emiss_scal_v(3);

    for (Index i = 0; i < 3; ++i) {
        Numeric e0 = a0_k0[i + (class1 - 1)]
                   + a0_k1[i + (class1 - 1)] * ev[i]
                   + a0_k2[i + (class1 - 1)] * eh[i];

        Numeric a0 = a0_eveh[i+ (class1 - 1) * 3];
        Numeric a1 = a1_eveh[i+ (class1 - 1) * 3];
        Numeric a2 = a2_eveh[i+ (class1 - 1) * 3];
        Numeric a3 = a3_eveh[i+ (class1 - 1) * 3];
        Numeric b0 = b0_eveh[i+ (class1 - 1) * 3];
        Numeric b1 = b1_eveh[i+ (class1 - 1) * 3];
        Numeric b2 = b2_eveh[i+ (class1 - 1) * 3];
        Numeric b3 = b3_eveh[i+ (class1 - 1) * 3];

        Numeric s1_v = (theta - 53.0) / -53.0 * (e0 - a0) / a0;
        Numeric em53_v = a3 * pow(53.0, 3.0) + a2 * pow(53.0, 2.0) + a1 * 53.0 + a0;
        Numeric s2_v = theta / 53.0 * (ev[i] - em53_v) / em53_v;
        Numeric s_v = 1.0 + s1_v + s2_v;

        Numeric emtheta_v = a3 * pow(theta, 3) + a2 * pow(theta, 2) + a1 * theta + a0;
        emiss_scal_v[i] = s_v * emtheta_v;

        Numeric s1_h = (theta - 53.0) / -53.0 * (e0 - b0) / b0;
        Numeric em53_h = b3 * pow(53.0, 3.0) + b2 * pow(53.0, 2.0) + b1 * 53.0 + b0;
        Numeric s2_h = theta / 53.0 * (eh[i] - em53_h) / em53_h;
        Numeric s_h = 1.0 + s1_h + s2_h;

        Numeric emtheta_h = b3 * pow(theta, 3) + b2 * pow(theta, 2) + b1 * theta + b0;
        emiss_scal_h[i] = s_h * emtheta_h;

    }

    Numeric emiss_h = interp_freq2(emiss_scal_h[0],
                                   emiss_scal_h[1],
                                   emiss_scal_h[2],
                                   freq,
                                   class2);
    Numeric emiss_v = interp_freq2(emiss_scal_v[0],
                                   emiss_scal_v[1],
                                   emiss_scal_v[2],
                                   freq,
                                   class2);

    if (emiss_v < emiss_h) {
        emiss_v = 0.5 * (emiss_v + emiss_h);
        emiss_h = emiss_v;
    }
    return std::make_pair(emiss_h, emiss_v);
}

std::ostream& operator<<(std::ostream& os, const TelsemAtlas& ta)
{
    os << ta.name << std::endl;

    return os;
}

////////////////////////////////////////////////////////////////////////////////
// Regression Coefficients
////////////////////////////////////////////////////////////////////////////////

const std::array<Numeric, 30> TelsemAtlas::a0_k0 = {
    0.11509, 0.091535, 0.34796, 0.10525, 0.16627, 0.24434,
    0.29217, 0.23809, 0.28954, 0.17516, 0.19459, 0.28697,
    0.10521, 0.12126, 0.30278, 0.18212, 0.19625, 0.14551,
    -0.19202, 0.5411, 0.03739, 0.10292, 0.5486, -0.058937,
    -0.022672, 0.44492, -0.058448, -0.33894, -0.17621, 0.14742
};

const std::array<Numeric, 30> TelsemAtlas::a0_k1 = {
    0.61168, 0.59095, 0.7918, 0.60271, 0.69213, 0.62218,
    0.32728, 0.34334, 0.37062, 0.51217, 0.4491, 0.50101,
    0.48913, 0.41932, 0.29734, 0.64474, 0.30637, 0.031107,
    1.0405, 0.17538, 1.3215, 0.61819, 0.31298, 1.7218,
    0.87761, 0.47583, 1.2583, 1.0959, 0.92842, 0.51033
};

const std::array<Numeric, 30> TelsemAtlas::a0_k2 = {
    0.26726, 0.32033, -0.14778, 0.28547, 0.13592, 0.13193,
    0.37178, 0.41813, 0.33875, 0.30203, 0.35479, 0.20189,
    0.40663, 0.47493, 0.40668, 0.14811, 0.52382, 0.86634,
    0.14286, 0.27164, -0.37947, 0.2737, 0.12001, -0.67315,
    0.13492, 0.065463, -0.19316, 0.24905, 0.25475, 0.34637
};

const std::array<Numeric, 30> TelsemAtlas::a0_eveh = {
    0.9592599869E+00, 0.9565299749E+00, 0.9511899948E+00,
    0.9560700059E+00, 0.9541199803E+00, 0.9483199716E+00,
    0.9461100101E+00, 0.9439799786E+00, 0.9387800097E+00,
    0.9317600131E+00, 0.9289000034E+00, 0.9236800075E+00,
    0.9208700061E+00, 0.9190599918E+00, 0.9105200171E+00,
    0.9162799716E+00, 0.8937299848E+00, 0.8014699817E+00,
    0.9570500255E+00, 0.9213600159E+00, 0.7893999815E+00,
    0.9639400244E+00, 0.9530599713E+00, 0.8850200176E+00,
    0.9685299993E+00, 0.9622600079E+00, 0.9118800163E+00,
    0.8997200131E+00, 0.9012699723E+00, 0.9107499719E+00
};

const std::array<Numeric, 30> TelsemAtlas::a1_eveh = {
    0.3627802414E-07, -0.7778328204E-08, 0.4396108011E-07,
    0.2503205394E-06, 0.1996262995E-06, 0.2929977541E-06,
    0.4190530660E-06, 0.3655744649E-06, 0.3519195673E-06,
    0.5574374313E-06, 0.5273076340E-06, 0.5376484182E-06,
    0.1026844529E-05, 0.9679998811E-06, 0.8616486866E-06,
    0.3180800832E-06, 0.2886778532E-06, 0.2310362675E-06,
    -0.1118036366E-06, -0.1502856577E-06, 0.4842232926E-07,
    -0.8410978580E-08, -0.3478669441E-07, 0.2209441590E-06,
    0.2485776633E-06, 0.1800235907E-06, 0.2510202251E-06,
    0.2687000915E-06, 0.1740325644E-06, 0.3562134339E-06
};

const std::array<Numeric, 30> TelsemAtlas::a2_eveh = {
    0.3067140824E-05, 0.2520012231E-05, 0.4831396382E-05,
    0.8213598448E-05, 0.7378375358E-05, 0.1022081960E-04,
    0.1225889173E-04, 0.1165553113E-04, 0.1188659007E-04,
    0.1693615741E-04, 0.1648317448E-04, 0.1715818144E-04,
    0.2744720041E-04, 0.2642072104E-04, 0.2671847506E-04,
    0.1349592094E-04, 0.1261523357E-04, 0.5447756394E-05,
    0.2064244654E-05, 0.1919016057E-06, 0.5940860319E-06,
    0.5334760772E-05, 0.4130339221E-05, 0.4104662821E-05,
    0.6530796327E-05, 0.5727014013E-05, 0.7451782039E-05,
    0.1071246970E-04, 0.9539280654E-05, 0.1034286015E-04
};

const std::array<Numeric, 30> TelsemAtlas::a3_eveh = {
    -0.2004991551E-07,-0.6895366056E-07, -0.2047409282E-06,
    -0.7322448425E-07, -0.1273002681E-06, -0.2729916844E-06,
    -0.9421125213E-07, -0.1683332300E-06, -0.2726891637E-06,
    -0.1317753799E-06, -0.2107972250E-06, -0.3556060904E-06,
    -0.1889465580E-06, -0.2757958271E-06, -0.4909850304E-06,
    0.7339644004E-08, -0.4058669560E-06, -0.4146343997E-06,
    0.6170279931E-07, -0.1998567996E-06, -0.4713119139E-07,
    -0.1361754887E-07, -0.1765622955E-06, -0.2348146637E-06,
    -0.3901189061E-07, -0.1305666189E-06, -0.1533838798E-06,
    -.2679148992E-07, -0.4441960044E-07, -0.1815613899E-06
};

const std::array<Numeric, 30> TelsemAtlas::b0_eveh = {
    0.9592599869E+00, 0.9565299749E+00, 0.9511899948E+00,
    0.9560700059E+00, 0.9541199803E+00, 0.9483199716E+00,
    0.9461100101E+00, 0.9439799786E+00, 0.9387800097E+00,
    0.9317600131E+00, 0.9289000034E+00, 0.9236800075E+00,
    0.9208700061E+00, 0.9190599918E+00, 0.9105200171E+00,
    0.9162799716E+00, 0.8937299848E+00, 0.8014699817E+00,
    0.9570500255E+00, 0.9213600159E+00, 0.7893999815E+00,
    0.9639400244E+00, 0.9530599713E+00, 0.8850200176E+00,
    0.9685299993E+00, 0.9622600079E+00, 0.9118800163E+00,
    0.8997200131E+00,0.9012699723E+00,0.9107499719E+00
};

const std::array<Numeric, 30> TelsemAtlas::b1_eveh = {
    0.3626608347E-07, -0.7786279177E-08, 0.4393379172E-07,
    0.2502746099E-06, 0.1995944388E-06, 0.2929554341E-06,
    0.4189516289E-06, 0.3655020180E-06, 0.3518483140E-06,
    0.5572838404E-06, 0.5271903092E-06, 0.5375342766E-06,
    0.1026605219E-05, 0.9677979733E-06, 0.8614680951E-06,
    0.3179358714E-06, 0.2884899004E-06, 0.2308632219E-06,
    -0.1118781370E-06, -0.1503948681E-06, 0.4834672396E-07,
    -0.8455684153E-08, -0.3485171618E-07, 0.2208606134E-06,
    0.2485595019E-06, 0.1799959364E-06, 0.2509846695E-06,
    0.2686167306E-06, 0.1739760478E-06, 0.3561317214E-06
};

const std::array<Numeric, 30> TelsemAtlas::b2_eveh = {
    0.3065537157E-05, 0.2518960400E-05, 0.4829731552E-05,
    0.8209894986E-05, 0.7375769655E-05, 0.1021809931E-04,
    0.1225203869E-04, 0.1165053800E-04, 0.1188218721E-04,
    0.1692612022E-04, 0.1647546378E-04, 0.1715117833E-04,
    0.2743142431E-04, 0.2640772436E-04, 0.2670711910E-04,
    0.1348545720E-04, 0.1260529825E-04, 0.5439695997E-05,
    0.2058213340E-05, 0.1860650656E-06, 0.5898303925E-06,
    0.5330772183E-05, 0.4126528893E-05, 0.4100859314E-05,
    0.6528573977E-05, 0.5725009032E-05, 0.7449450095E-05,
    0.1070590315E-04, 0.9534271157E-05, 0.1033751869E-04
};

const std::array<Numeric, 30> TelsemAtlas::b3_eveh = {
    -0.1370247134E-06, -0.1436897747E-06, -0.2954870411E-06,
    -0.3118435643E-06, -0.2916583242E-06, -0.4311032171E-06,
    -0.5048401022E-06, -0.4662823869E-06, -0.5206445053E-06,
    -0.7210980471E-06, -0.6662896794E-06, -0.7548637200E-06,
    -0.1110204039E-05, -0.1030801400E-05, -0.1140921199E-05,
    -0.6330818110E-06, -0.9186441048E-06, -0.7947813856E-06,
    -0.3242539890E-06, -0.5027602583E-06, -0.2777987334E-06,
    -0.2747250676E-06, -0.3811997260E-06, -0.4102405455E-06,
    -0.1994112324E-06, -0.2555484855E-06, -0.2842682534E-06,
    -0.4413041665E-06, -0.3717419474E-06, -0.4975536854E-06
};

const std::array<Numeric, 4> TelsemAtlas::rapport43_32 = {0.62, 0.37, 0.46, 0.63};
const std::array<Numeric, 4> TelsemAtlas::rapport54_43 = {0.30, 0.60, 0.47, 0.35};
