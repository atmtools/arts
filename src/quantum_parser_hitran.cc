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
void parse_space(Rational& qn, String& s, const Index species);
void parse_a1_br_hitran(Rational& qn, String& s, const Index species);
void parse_a1_pm_hitran(Rational& qn, String& s, const Index species);
void parse_a1_sym_hitran(Rational& qn, String& s, const Index species);
void parse_a1_s_hitran(Rational& qn, String& s, const Index species);
void parse_a1_x_hitran(Rational& qn, String& s, const Index species);
void parse_a2_hitran(Rational& qn, String& s, const Index species);
void parse_a3_hitran(Rational& qn, String& s, const Index species);
void parse_a4_hitran(Rational& qn, String& s, const Index species);
void parse_a5_hitran(Rational& qn, String& s, const Index species);
void parse_i1_hitran(Rational& qn, String& s, const Index species);
void parse_i2_hitran(Rational& qn, String& s, const Index species);
void parse_i3_hitran(Rational& qn, String& s, const Index species);
void parse_f51_hitran(Rational& qn, String& s, const Index species);

// Postprocessing functions for calculation of implicit quantum numbers
void postprocess_group1_hitran(QuantumNumberRecord& qnr, const Index species);
void postprocess_group2_hitran(QuantumNumberRecord& qnr, const Index species);
void postprocess_group5_hitran(QuantumNumberRecord& qnr, const Index species);
void postprocess_group6_hitran(QuantumNumberRecord& qnr, const Index species);
void postprocess_group6oh_hitran(QuantumNumberRecord& qnr, const Index species);



QuantumParserHITRAN2004::QuantumParserHITRAN2004()
{
    using global_data::species_data;

#define SKIP_X_SPACES(container, nspaces) \
    container.push_back_n(QuantumFieldDescription(QN_FINAL_ENTRY, parse_space), nspaces)

    // HITRAN Classes
    mclass.resize(CI_FINAL);

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
        this_class.push_back(QuantumFieldDescription(QN_X,  parse_a1_x_hitran));
        this_class.push_back(QuantumFieldDescription(QN_v1, parse_i2_hitran));
    }

    // Class 3
    {
        Array<QuantumFieldDescription>& this_class = mclass[CI_CLASS3];
        SKIP_X_SPACES(this_class, 7);
        this_class.push_back(QuantumFieldDescription(QN_X,     parse_a1_x_hitran));
        this_class.push_back(QuantumFieldDescription(QN_Omega, parse_a3_hitran));
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
        this_class.push_back(QuantumFieldDescription(QN_pm, parse_a1_pm_hitran));
        this_class.push_back(QuantumFieldDescription(QN_r,  parse_i1_hitran));
        this_class.push_back(QuantumFieldDescription(QN_S_global, parse_a1_s_hitran));
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
    mgroup.resize(GI_FINAL);

    // Group 1
    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP1].upper;
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Ka,  parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Kc,  parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a5_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_sym_hitran));
    }

    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP1].lower;
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Ka,  parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Kc,  parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a5_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_sym_hitran));
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
        this_group.push_back(QuantumFieldDescription(QN_dJ,  parse_a1_br_hitran));
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_sym_hitran));
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
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_sym_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a4_hitran));
    }

    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP4].lower;
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_K,   parse_i3_hitran));
        this_group.push_back(QuantumFieldDescription(QN_l,   parse_i2_hitran));
        this_group.push_back(QuantumFieldDescription(QN_C,   parse_a2_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_sym_hitran));
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
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_sym_hitran));
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
        this_group.push_back(QuantumFieldDescription(QN_dJ,  parse_a1_br_hitran));
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_f51_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_sym_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a5_hitran));
    }

    // Group 6 (OH)
    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP6OH].upper;
        SKIP_X_SPACES(this_group, 10);
        this_group.push_back(QuantumFieldDescription(QN_F, parse_a5_hitran));
    }

    {
        Array<QuantumFieldDescription>& this_group = mgroup[GI_GROUP6OH].lower;
        SKIP_X_SPACES(this_group, 1);
        this_group.push_back(QuantumFieldDescription(QN_dN,  parse_a1_br_hitran));
        this_group.push_back(QuantumFieldDescription(QN_dJ,  parse_a1_br_hitran));
        this_group.push_back(QuantumFieldDescription(QN_J,   parse_f51_hitran));
        this_group.push_back(QuantumFieldDescription(QN_Sym, parse_a1_sym_hitran));
        this_group.push_back(QuantumFieldDescription(QN_F,   parse_a5_hitran));
    }
    
    mspecies.resize(species_data.nelem());

    SetClassGroup("CO2", CI_CLASS5, GI_GROUP2);
    SetClassGroup("NO",  CI_CLASS3, GI_UNDEFINED); // GROUP6
    SetClassGroup("O2",  CI_CLASS2, GI_GROUP5);
    SetClassGroup("ClO", CI_CLASS3, GI_UNDEFINED); // GROUP6
    SetClassGroup("H2O", CI_CLASS6, GI_GROUP1);
    SetClassGroup("HO2", CI_CLASS6, GI_GROUP1);
    SetClassGroup("NO2", CI_CLASS6, GI_GROUP1);
    SetClassGroup("O3" , CI_CLASS6, GI_GROUP1);
    SetClassGroup("OH",  CI_CLASS3, GI_UNDEFINED); // GROUP6OH

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
        mclass[qclass][i].Parse(qnr.Upper(), qstr, species);
    }

    qstr = quantum_string.substr(15, 15);
    for (Index i = 0; i < mclass[qclass].nelem(); i++)
    {
        mclass[qclass][i].Parse(qnr.Lower(), qstr, species);
    }

    if (qgroup != GI_UNDEFINED)
    {
        qstr = quantum_string.substr(30, 15);
        for (Index i = 0; i < mgroup[qgroup].upper.nelem(); i++)
        {
            mgroup[qgroup].upper[i].Parse(qnr.Upper(), qstr, species);
        }
    }

    qstr = quantum_string.substr(45, 15);
    if (qgroup != GI_UNDEFINED)
    {
        for (Index i = 0; i < mgroup[qgroup].lower.nelem(); i++)
        {
            mgroup[qgroup].lower[i].Parse(qnr.Lower(), qstr, species);
        }
    }

    switch (qgroup)
    {
        case GI_GROUP1: postprocess_group1_hitran(qnr, species); break;
        case GI_GROUP2: postprocess_group2_hitran(qnr, species); break;
        case GI_GROUP5: postprocess_group5_hitran(qnr, species); break;
        case GI_GROUP6: postprocess_group6_hitran(qnr, species); break;
        case GI_GROUP6OH: postprocess_group6oh_hitran(qnr, species); break;
        default: break;
    }
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

void parse_space(Rational& qn, String& s, const Index /* species */)
{
    s.erase(0, 1);
    qn = RATIONAL_UNDEFINED;
}


void parse_a1_br_hitran(Rational& qn, String& s, const Index /* species */)
{
    qn = 'Q' - s[0];
    if (qn > 4)
        throw std::runtime_error("Error parsing quantum number Br");
    s.erase(0, 1);
}


void parse_a1_pm_hitran(Rational& qn, String& s, const Index /* species */)
{
    // FIXME: How to handle plus/minus?
    s.erase(0, 1);
    qn = RATIONAL_UNDEFINED;
}


void parse_a1_sym_hitran(Rational& qn, String& s, const Index species)
{
    qn = RATIONAL_UNDEFINED;
    const char ch = s[0];

     if ((ch == '+' || ch == '-')
         && (species == species_index_from_species_name("HO2")
             || species == species_index_from_species_name("NO2")))
    {
        if (ch == '+')
            qn = Rational(1, 2);
        else if (ch == '-')
            qn = Rational(-1 , 2);
        else
            throw std::runtime_error("Error parsing quantum number Sym");
    }
    else if ((ch == 'd' || ch == 'q')
             && (species == species_index_from_species_name("O2")
                 || species == species_index_from_species_name("N2")))
    {
        // FIXME: How to deal with symmetry parameter?
    }
    else if (ch =='g'
             && species == species_index_from_species_name("O2"))
    {
        // FIXME: What does 'g' mean? Docs only talk about d and q for O2
    }
    else if (ch == 'e' || ch == 'f')
    {
        // FIXME: How to handle linedoubling values?
    }
    else if (ch == '+' || ch == '-')
    {
        // FIXME: How to handle symmetry parameter?
    }
    else if (ch != ' ')
        throw std::runtime_error("Error parsing quantum number Sym");

    s.erase(0, 1);
}


void parse_a1_s_hitran(Rational& qn, String& s, const Index /* species */)
{
    // FIXME: How to handle S?
    s.erase(0, 1);
    qn = RATIONAL_UNDEFINED;
}


void parse_a1_x_hitran(Rational& qn, String& s, const Index /* species */)
{
    // FIXME: How to handle X?
    s.erase(0, 1);
    qn = RATIONAL_UNDEFINED;
}


void parse_a2_hitran(Rational& qn, String& s, const Index /* species */)
{
    // FIXME OLE
    s.erase(0, 2);
    qn = RATIONAL_UNDEFINED;
}


void parse_a3_hitran(Rational& qn, String& s, const Index /* species */)
{
    extract(qn, s, 3);
}


void parse_a4_hitran(Rational& qn, String& s, const Index /* species */)
{
    // FIXME OLE
    s.erase(0, 4);
    qn = RATIONAL_UNDEFINED;
}


void parse_a5_hitran(Rational& qn, String& s, const Index /* species */ )
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


void parse_i1_hitran(Rational& qn, String& s, const Index /* species */)
{
    Index i;
    extract(i, s, 1);
    qn = i;
}


void parse_i2_hitran(Rational& qn, String& s, const Index /* species */)
{
    Index i;
    extract(i, s, 2);
    qn = i;
}


void parse_i3_hitran(Rational& qn, String& s, const Index /* species */)
{
    Index i;
    extract(i, s, 3);
    qn = i;
}


void parse_f51_hitran(Rational& qn, String& s, const Index /* species */)
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
                throw std::runtime_error("Error parsing quantum number of type F5.1");
            }

            if (as[1] == "5")
                qn = Rational(nom * 2 + 1, 2);
            else if (as[1] == "0")
                qn = nom;
            else
            {
                throw std::runtime_error("Error parsing quantum number of type F5.1");
            }
            s.erase(0, 5);
        }
        else
        {
            bool valid_num = true;
            for (Index i = 0; valid_num && i < qnf.nelem(); i++)
                if (qnf[i] < '0' || qnf[i] > '9')
                    valid_num = false;
            if (!valid_num)
                throw std::runtime_error("Error parsing quantum number of type F5.1");

            Index q;
            extract(q, s, qnf.nelem());
        }
    }
}



/////////////////////////////
/// Post-processing functions
/////////////////////////////

void postprocess_group1_hitran(QuantumNumberRecord& qnr, const Index species)
{
    if (species == species_index_from_species_name("HO2")
        || species == species_index_from_species_name("NO2"))
    {
        // For these species, N is stored in field J, so we copy it here
        // to the correct quantum number
        qnr.SetUpper(QN_N, qnr.Upper(QN_J));
        qnr.SetLower(QN_N, qnr.Lower(QN_J));
        // Now we can calculate the correct J by using the Symmetry quantum
        // number. However, it does not represent Symmetry for these
        // species, but +1/2 or -1/2
        qnr.SetUpper(QN_J, qnr.Upper(QN_N) + qnr.Upper(QN_Sym));
        qnr.SetLower(QN_J, qnr.Lower(QN_N) + qnr.Lower(QN_Sym));

        // We don't need Sym after this point
        qnr.SetUpper(QN_Sym, RATIONAL_UNDEFINED);
        qnr.SetLower(QN_Sym, RATIONAL_UNDEFINED);
    }
}


void postprocess_group2_hitran(QuantumNumberRecord& qnr, const Index /* species */)
{
    qnr.SetUpper(QN_J, qnr.Lower(QN_J) - qnr.Lower(QN_dJ));

    // We don't need dN and dJ after this point
    qnr.SetLower(QN_dJ, RATIONAL_UNDEFINED);
}


void postprocess_group5_hitran(QuantumNumberRecord& qnr, const Index /* species */)
{
    qnr.SetUpper(QN_N, qnr.Lower(QN_N) - qnr.Lower(QN_dN));
    qnr.SetUpper(QN_J, qnr.Lower(QN_J) - qnr.Lower(QN_dJ));

    // We don't need dN and dJ after this point
    qnr.SetLower(QN_dN, RATIONAL_UNDEFINED);
    qnr.SetLower(QN_dJ, RATIONAL_UNDEFINED);
}


void postprocess_group6_hitran(QuantumNumberRecord& qnr, const Index species)
{
    qnr.SetUpper(QN_J, qnr.Lower(QN_J) - qnr.Lower(QN_dJ));

    if (species == species_index_from_species_name("NO"))
    {
        // This might also be true for ClO, but that is hidden in HITRAN,
        // it might be the opposite order of the plus-minus for ClO
        if (qnr.Upper(QN_Omega).Nom() == 3)
            qnr.SetUpper(QN_N, qnr.Upper(QN_J) + Rational(1,2));
        else if (qnr.Upper(QN_Omega).Nom() == 1)
            qnr.SetUpper(QN_N, qnr.Upper(QN_J) - Rational(1,2));

        if (qnr.Lower(QN_Omega).Nom() == 3)
            qnr.SetLower(QN_N, qnr.Upper(QN_J) + Rational(1,2));
        else if (qnr.Lower(QN_Omega).Nom() == 1)
            qnr.SetLower(QN_N, qnr.Upper(QN_J) - Rational(1,2));
    }

    // We don't need dJ after this point
    qnr.SetLower(QN_dJ, RATIONAL_UNDEFINED);
}


void postprocess_group6oh_hitran(QuantumNumberRecord& qnr, const Index /* species */)
{
    qnr.SetUpper(QN_J, qnr.Lower(QN_J) - qnr.Lower(QN_dJ));
    // FIXME: It is not 100% how to calculate N for OH, therefore
    // it is ignored for now

    // We don't need dN and dJ after this point
    qnr.SetLower(QN_dN, RATIONAL_UNDEFINED);
    qnr.SetLower(QN_dJ, RATIONAL_UNDEFINED);
}
