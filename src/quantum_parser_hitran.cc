/* Copyright (C) 2015, The ARTS Developers.

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

/** \file

  Parser for quantum numbers from spectroscopic catalogs.

  \author Oliver Lemke
*/

#include "quantum_parser_hitran.h"
#include "global_data.h"
#include "absorption.h"
#include "messages.h"


// Parsing functions for fields
void parse_space(Rational& qn, String& s);
void parse_a1_hitran(Rational& qn, String& s);
void parse_a1_br_hitran(Rational& qn, String& s);
void parse_a2_hitran(Rational& qn, String& s);
void parse_a3_hitran(Rational& qn, String& s);
void parse_a4_hitran(Rational& qn, String& s);
void parse_a5_hitran(Rational& qn, String& s);
void parse_i1_hitran(Rational& qn, String& s);
void parse_i2_hitran(Rational& qn, String& s);
void parse_i3_hitran(Rational& qn, String& s);
void parse_f51_hitran(Rational& qn, String& s);

// Postprocessing functions for calculation of implicit quantum numbers
void postprocess_group5_hitran(QuantumNumberRecord& qnr, const Index species);



QuantumParserHITRAN2004::QuantumParserHITRAN2004()
{
    using global_data::species_data;

#define SKIP_X_SPACES(container, nspaces) \
    container.push_back_n(QuantumFieldDescription(QN_FINAL_ENTRY, parse_space), nspaces)

    // HITRAN Classes
    mclass.resize(10);

    // Class 1
    {
        Array<QuantumFieldDescription>& this_class = mclass[CI_CLASS1];
        SKIP_X_SPACES(this_class, 13);
        this_class.push_back(QuantumFieldDescription(QN_v1, parse_i2_hitran));
    }

    // Class 2
    {
        Array<QuantumFieldDescription>& this_class = mclass[CI_CLASS2];
        SKIP_X_SPACES(this_class, 12);
        this_class.push_back(QuantumFieldDescription(QN_X,  parse_a1_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v1, parse_i2_hitran));
    }

    // Class 3
    {
        Array<QuantumFieldDescription>& this_class = mclass[CI_CLASS3];
        SKIP_X_SPACES(this_class, 7);
        this_class.push_back(QuantumFieldDescription(QN_X,  parse_a1_hitran));
        this_class.push_back(QuantumFieldDescription(QN_i,  parse_a3_hitran));
        SKIP_X_SPACES(this_class, 2);
        this_class.push_back(QuantumFieldDescription(QN_v1, parse_i2_hitran));
    }

    // Class 4
    {
        Array<QuantumFieldDescription>& this_class = mclass[CI_CLASS4];
        SKIP_X_SPACES(this_class, 7);
        this_class.push_back(QuantumFieldDescription(QN_v1, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v2, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_l2, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v3, parse_i2_hitran));
    }

    // Class 5
    {
        Array<QuantumFieldDescription>& this_class = mclass[CI_CLASS5];
        SKIP_X_SPACES(this_class, 6);
        this_class.push_back(QuantumFieldDescription(QN_v1, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v2, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_l2, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v3, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_r,  parse_i1_hitran));
    }

    // Class 6
    {
        Array<QuantumFieldDescription>& this_class = mclass[CI_CLASS6];
        SKIP_X_SPACES(this_class, 9);
        this_class.push_back(QuantumFieldDescription(QN_v1, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v2, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v3, parse_i2_hitran));
    }

    // Class 7
    {
        Array<QuantumFieldDescription>& this_class = mclass[CI_CLASS7];
        this_class.push_back(QuantumFieldDescription(QN_v1, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v2, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v3, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v4, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v5, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_l,  parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_pm, parse_a1_hitran));
        this_class.push_back(QuantumFieldDescription(QN_r,  parse_i1_hitran));
        this_class.push_back(QuantumFieldDescription(QN_S_global, parse_a1_hitran));
    }

    // Class 8
    {
        Array<QuantumFieldDescription>& this_class = mclass[CI_CLASS8];
        SKIP_X_SPACES(this_class, 5);
        this_class.push_back(QuantumFieldDescription(QN_v1, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v2, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v3, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v4, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_S_global, parse_i2_hitran));
    }

    // Class 9
    {
        Array<QuantumFieldDescription>& this_class = mclass[CI_CLASS9];
        SKIP_X_SPACES(this_class, 3);
        this_class.push_back(QuantumFieldDescription(QN_v1, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v2, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v3, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v4, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v5, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v6, parse_i2_hitran));
    }

    // Class 10
    {
        Array<QuantumFieldDescription>& this_class = mclass[CI_CLASS10];
        SKIP_X_SPACES(this_class, 3);
        this_class.push_back(QuantumFieldDescription(QN_v1, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v2, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v3, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v4, parse_i2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_n_global, parse_a2_hitran));
        this_class.push_back(QuantumFieldDescription(QN_C,  parse_a2_hitran));
    }


    // HITRAN Groups
    mgroup.resize(6);

    // Group 1
    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP1].upper;
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_K,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Kc,  parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a5_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_hitran));
    }

    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP1].lower;
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_K,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Kc,  parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a5_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_hitran));
    }

    // Group 2
    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP2].upper;
        SKIP_X_SPACES(this_group, 10);
        this_group.push_back(QuantumFieldDescription(QN_F, parse_a5_hitran));
    }

    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP2].lower;
        SKIP_X_SPACES(this_group, 5);
        this_group.push_back(QuantumFieldDescription(QN_dJ,  parse_a1_hitran));
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a5_hitran));
    }

    // Group 3
    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP3].upper;
        SKIP_X_SPACES(this_group, 2);
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_C,   parse_a2_hitran));
        this_group.push_back(QuantumFieldDescription(QN_alpha,  parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a5_hitran));
    }

    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP3].lower;
        SKIP_X_SPACES(this_group, 2);
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_K,   parse_a2_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Kc,  parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a5_hitran));
    }

    // Group 4
    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP4].upper;
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_K,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_l,   parse_i2_hitran));
        this_group.push_back(QuantumFieldDescription(QN_C,   parse_a2_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a4_hitran));
    }

    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP4].lower;
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_K,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_l,   parse_i2_hitran));
        this_group.push_back(QuantumFieldDescription(QN_C,   parse_a2_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a4_hitran));
    }

    // Group 5
    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP5].upper;
        SKIP_X_SPACES(this_group, 10);
        this_group.push_back(QuantumFieldDescription(QN_F, parse_a5_hitran));
    }

    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP5].lower;
        SKIP_X_SPACES(this_group, 1);
        this_group.push_back(QuantumFieldDescription(QN_dN,  parse_a1_br_hitran));
        this_group.push_back(QuantumFieldDescription(QN_N,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_dJ,  parse_a1_br_hitran));
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a5_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_hitran));
    }

    // Group 6
    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP6].upper;
        SKIP_X_SPACES(this_group, 10);
        this_group.push_back(QuantumFieldDescription(QN_F, parse_a5_hitran));
    }

    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP6].lower;
        SKIP_X_SPACES(this_group, 3);
        this_group.push_back(QuantumFieldDescription(QN_dJ,  parse_a1_hitran));
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_f51_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a5_hitran));
    }

    mspecies.resize(species_data.nelem());

//    SetClassGroup("H2O", CI_CLASS6, GI_GROUP1);
    SetClassGroup("O2", CI_CLASS2, GI_GROUP5);

#undef SKIP_X_SPACES

}


void QuantumParserHITRAN2004::Parse(QuantumNumberRecord& qnr,
                                    const String& quantum_string,
                                    const Index species) const
{
    const QuantumClassGroup& qcg = mspecies[species];
    const Index qclass = mspecies[species].iclass;
    const Index qgroup = mspecies[species].igroup;

    if (qcg.iclass == -1 || qcg.igroup == -1)
        return;

    String qstr;

    qstr = quantum_string.substr(0,15);
    for (Index i = 0; i < mclass[qclass].nelem(); i++)
    {
        mclass[qclass][i].Parse(qnr.Upper(), qstr);
    }

    qstr = quantum_string.substr(15, 15);
    for (Index i = 0; i < mclass[qclass].nelem(); i++)
    {
        mclass[qclass][i].Parse(qnr.Lower(), qstr);
    }

    qstr = quantum_string.substr(30, 15);
    for (Index i = 0; i < mgroup[qgroup].upper.nelem(); i++)
    {
        mgroup[qgroup].upper[i].Parse(qnr.Upper(), qstr);
    }

    qstr = quantum_string.substr(45, 15);
    for (Index i = 0; i < mgroup[qgroup].lower.nelem(); i++)
    {
        mgroup[qgroup].lower[i].Parse(qnr.Lower(), qstr);
    }

    if (qgroup == GI_GROUP5)
        postprocess_group5_hitran(qnr, species);
}


void QuantumParserHITRAN2004::SetClassGroup(const String& species_name,
                                                       const ClassIds iclass,
                                                       const GroupIds igroup)
{
    Index species = species_index_from_species_name(species_name);
    mspecies[species].iclass = iclass;
    mspecies[species].igroup = igroup;
}


/////////////////////
/// Parsing functions
/////////////////////

void parse_space(Rational& qn, String& s)
{
    s.erase(0, 1);
    qn = RATIONAL_UNDEFINED;
}


void parse_a1_hitran(Rational& qn, String& s)
{
    s.erase(0, 1);
    qn = RATIONAL_UNDEFINED;
}


void parse_a1_br_hitran(Rational& qn, String& s)
{
    qn = 'Q' - s[0];
    if (qn > 2)
        throw std::runtime_error("Error parsing quantum number Br");
    s.erase(0, 1);
}


void parse_a2_hitran(Rational& qn, String& s)
{
    s.erase(0, 2);
    qn = RATIONAL_UNDEFINED;
}


void parse_a3_hitran(Rational& qn, String& s)
{
    s.erase(0, 3);
    qn = RATIONAL_UNDEFINED;
}


void parse_a4_hitran(Rational& qn, String& s)
{
    s.erase(0, 4);
    qn = RATIONAL_UNDEFINED;
}


void parse_a5_hitran(Rational& qn, String& s)
{
    String qnf = s.substr(0, 5);
    qn = RATIONAL_UNDEFINED;

    qnf.trim();
    if (qnf.nelem())
    {
        ArrayOfString as;
        qnf.split(as, ".");
        if (as.nelem() == 2)
        {
            Index nom;
            char* endptr;

            nom = strtol(as[0].c_str(), &endptr, 10);
            if (endptr != as[0].c_str()+as[0].nelem())
            {
                throw std::runtime_error("Error parsing quantum number of type A5");
            }

            if (as[1] == "5")
                qn = Rational(nom * 2 + 1, 2);
            else if (as[1] == "0")
                qn = nom;
            else
            {
                throw std::runtime_error("Error parsing quantum number of type A5");
            }
        }
    }
    s.erase(0, 5);
}


void parse_i1_hitran(Rational& qn, String& s)
{
    Index i;
    extract(i, s, 1);
    qn = i;
}


void parse_i2_hitran(Rational& qn, String& s)
{
    Index i;
    extract(i, s, 2);
    cout << s.substr(0, 2);
    qn = i;
}


void parse_i3_hitran(Rational& qn, String& s)
{
    Index i;
    extract(i, s, 3);
    qn = i;
}


void parse_f51_hitran(Rational& qn, String& s)
{
    s.erase(0, 7);
    qn = RATIONAL_UNDEFINED;
}



/////////////////////////////
/// Post-processing functions
/////////////////////////////

void postprocess_group5_hitran(QuantumNumberRecord& qnr, const Index /* species */)
{
    qnr.SetUpper(QN_N, qnr.Lower(QN_N) - qnr.Lower(QN_dN));
    qnr.SetUpper(QN_J, qnr.Lower(QN_J) - qnr.Lower(QN_dJ));

    // We don't need dN and dJ after this point
    qnr.SetLower(QN_dN, RATIONAL_UNDEFINED);
    qnr.SetLower(QN_dJ, RATIONAL_UNDEFINED);
}
