/* Copyright (C) 2019
 Richard Larsson <larsson@mps.mpg.de>

 
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

/** Contains the absorption namespace
 * @file   m_absorptionlines.cc
 * @author Richard Larsson
 * @date   2019-09-11
 * 
 * @brief  Contains the user interaction with absorption lines
 **/

#include "absorptionlines.h"
#include "auto_md.h"
#include "file.h"

void ReadHITRAN(ArrayOfAbsorptionLines& abs_lines2,
                const String& hitran_file,
                const Numeric& fmax,
                const Verbosity&)
{
  ifstream is;
  open_input_file(is, hitran_file);
  
  std::vector<Absorption::SingleLineExternal> v;
  v.resize(0);
  
  bool go_on = true;
  Index n = 0;
  while (go_on) {
    v.push_back(Absorption::ReadFromHitran2004Stream(is));
    
    if(v.back().bad) {
      v.pop_back();
      go_on = false;
    }
    else if(v.back().line.F0() > fmax) {
      v.pop_back();
      go_on = false;
    }
    else if(v.back().quantumidentity.Species() < 0 or v.back().quantumidentity.Isotopologue() < 0) {
      v.pop_back();
      go_on = true;
      n++;
    }
  }
  
  auto x = Absorption::split_list_of_external_lines(v, {QuantumNumberType::J}, {QuantumNumberType::v1});
  abs_lines2.resize(0);
  abs_lines2.reserve(x.size());
  for(auto& lines: x)
    abs_lines2.push_back(lines);
  
}

void abs_linesWriteSplitXML(const ArrayOfAbsorptionLines& abs_lines2,
                            const String& basename,
                            const Verbosity& verbosity)
{
  std::map<String, int> names;

  String true_basename = basename;
  if (not(true_basename.back() == '.' or true_basename.back() == '/'))
    true_basename += '.';

  for (auto& lines : abs_lines2) {
    auto name = lines.SpeciesName();
    const String fname = true_basename + name;

    WriteXML("ascii", lines,
             fname + '.' + std::to_string(names[name]++) + ".xml",
             0, "", "", "", verbosity);
    
    std::cout << lines.MetaData() << '\n';
  }
}
