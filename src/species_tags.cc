#include <cfloat>

#include "species_tags.h"

namespace Species {
Tag::Tag(String def) : spec_ind(-1), lower_freq(-1), upper_freq(-1), type(TagType::Plain), cia_2nd_species(Species::FINAL), cia_dataset_index(-1) {
  // Save input string for error messages:
  String def_original = def;

  // Name of species and isotopologue (aux variables):
  String name, isoname;
  // Aux index:
  Index n;

  // We cannot set a default value for the isotopologue, because the
  // default should be `ALL' and the value for `ALL' depends on the
  // species.

  // Extract the species name:
  n = def.find('-');  // find the '-'
  if (n != def.npos) {
    name = def.substr(0, n);  // Extract before '-'
    def.erase(0, n + 1);      // Remove from def
  } else {
    // n==def.npos means that def does not contain a '-'. In that
    // case we assume that it contains just a species name and
    // nothing else
    name = def;
    def = "";
  }
  
  // Remove whitespace
  name.trim();

  // Obtain species index from species name.
  // (This will also remove possible whitespace.)
  const Species spec = fromShortName(name);
  ARTS_USER_ERROR_IF (not good_enum(spec), "Bad species name: ", name, " extracted from: ", def_original)

  // Check if species name contains the special tag for
  // Faraday Rotation
  if (spec == Species::free_electrons) {
    constexpr Index ind = find_species_index(Species::free_electrons);
    spec_ind = ind;
    type = TagType::FreeElectrons;
    return;
  }

  // Check if species name contains the special tag for
  // Particles
  if (spec == Species::particles) {
    constexpr Index ind = find_species_index(Species::particles);
    spec_ind = ind;
    type = TagType::Particles;
    return;
  }

  // Set "all" species per default by leaving the joker in
  spec_ind = find_species_index(spec);
  ARTS_USER_ERROR_IF(spec_ind < 0, "Bad species extracted: ", spec, " from ", def_original)
  if (0 == def.nelem()) {
    return;
  }

  // Extract the isotopologue name/Zeeman flag:
  n = def.find('-');  // find the '-'
  if (n != def.npos) {
    isoname = def.substr(0, n);  // Extract before '-'
    def.erase(0, n + 1);         // Remove from def

    if ("Z" == isoname) {
      type = TagType::Zeeman;
      // Zeeman flag was present, now extract the isotopologue name:
      n = def.find('-');  // find the '-'
      if (n != def.npos) {
        isoname = def.substr(0, n);  // Extract before '-'
        def.erase(0, n + 1);         // Remove from def
      } else {
        // n==def.npos means that def does not contain a '-'. In that
        // case we assume that it contains just the isotopologue name and
        // nothing else.
        isoname = def;
        def = "";
      }
    }

    if ("HXSEC" == isoname) {
      type = TagType::HitranXsec;
      // Hitran Xsec flag was present, now extract the isotopologue name:
      n = def.find('-');  // find the '-'
      if (n != def.npos) {
        isoname = def.substr(0, n);  // Extract before '-'
        def.erase(0, n + 1);         // Remove from def
      } else {
        // n==def.npos means that def does not contain a '-'. In that
        // case we assume that it contains just the isotopologue name and
        // nothing else.
        isoname = def;
        def = "";
      }
    }
  } else {
    // n==def.npos means that def does not contain a '-'. In that
    // case we assume that it contains just the isotopologue name or
    // Zeeman flag and nothing else.
    isoname = def;
    def = "";
    if ("Z" == isoname) {
      type = TagType::Zeeman;
      // This means that there is nothing else to parse. Apparently
      // the user wants all isotopologues and no frequency limits.
      return;
    }
    if ("HXSEC" == isoname) {
      type = TagType::HitranXsec;
      // This means that there is nothing else to parse. Apparently
      // the user wants all isotopologues and no frequency limits.
      return;
    }
  }

  if ("*" == isoname) {
    // The user wants all isotopologues. This already set
  }
  else if ("CIA" == isoname)  // Check for "cia":
  {
    // The user wants this to use the CIA catalog:
    type = TagType::Cia;

    // We have to read in the second species, and the dataset index
    n = def.find('-');  // find the '-'

    ARTS_USER_ERROR_IF (n == def.npos,
      "Invalid species tag ", def_original, ".\n"
      "I am missing a minus sign (and a dataset index after that.)")

    String otherspec = def.substr(0, n);  // Extract before '-'
    def.erase(0, n + 1);                  // Remove from def

    cia_2nd_species = fromShortName(otherspec);
    ARTS_USER_ERROR_IF (not good_enum(cia_2nd_species), "Bad species name: ", otherspec, " extracted from: ", def_original)

    // Convert remaining def to dataset index.

    // Check that everything remaining is just numbers.
    for (Index i = 0; i < def.nelem(); ++i)
      ARTS_USER_ERROR_IF (!isdigit(def[i]),
          "Invalid species tag ", def_original, ".\n"
          "The tag should end with a dataset index")

    // Do the conversion from string to index:
    std::istringstream is(def);
    is >> cia_dataset_index;

    def = "";
  } else {
    spec_ind = find_species_index(spec, isoname);

    // Check if we found a matching isotopologue:
    ARTS_USER_ERROR_IF (spec_ind < 0,
        "Isotopologue ", isoname, " is not a valid isotopologue or "
        "absorption model for species ", name, ".\n"
        "Valid options are:\n", isotopologues_names(spec))

    // Check if the found isotopologue represents a predefined model
    // (continuum or full absorption model) and set the type accordingly:
    if (!isdigit(isoname[0])) type = TagType::Predefined;
  }

  if (0 == def.nelem()) {
    // This means that there is nothing else to parse. Apparently
    // the user wants no frequency limits.  Frequency defaults are
    // already set, so we can return directly.

    return;
  }

  ARTS_USER_ERROR_IF (def[0] != '*' && !isdigit(def[0]),
    "Expected frequency limits, but got \"", def, "\"")

  // Look for the two frequency limits:

  // Extract first frequency
  n = def.find('-');  // find the '-'
  if (n != def.npos) {
    // Frequency as a String:
    String fname;
    fname = def.substr(0, n);  // Extract before '-'
    def.erase(0, n + 1);       // Remove from def

    // Check for joker:
    if ("*" == fname) {
      // The default for lower_freq is already -1, meaning `ALL'.
      // So there is nothing to do here.
    } else {
      ARTS_USER_ERROR_IF (!isdigit(fname[0]),
        "Expected frequency limit, but got \"", fname, "\"")
        // Convert to Numeric:
      char* endptr;
      lower_freq = strtod(fname.c_str(), &endptr);
      ARTS_USER_ERROR_IF (endptr != fname.c_str() + fname.nelem(),
          "Error parsing frequency limit \"", fname, "\"")
    }
  } else {
    // n==def.npos means that def does not contain a '-'. In this
    // case that is not allowed!
    ARTS_USER_ERROR (
        "You must either specify both frequency limits\n"
        "(at least with jokers), or none.");
  }

  // Now there should only be the upper frequency left in def.
  // Check for joker:
  if ("*" == def) {
    // The default for upper_freq is already -1, meaning `ALL'.
    // So there is nothing to do here.
  } else {
    ARTS_USER_ERROR_IF (!isdigit(def[0]),
      "Expected frequency limit, but got \"", def, "\"")
    // Convert to Numeric:
    char* endptr;
    upper_freq = strtod(def.c_str(), &endptr);
    ARTS_USER_ERROR_IF (endptr != def.c_str() + def.nelem(),
      "Error parsing frequency limit \"", def, "\"")
  }
}

//! Return the full name of the tag.
/*! 
 * Examples:
 * 
 * \verbatim
 * O3-*-*-*         : All O3 lines
 * O3-nl            : O3, but without any lines
 * O3-666-*-*       : All O3-666 lines
 * O3-*-500e9-501e9 : All O3 lines between 500 and 501 GHz.
 * \endverbatim
 * 
 * \return The tag name as a string.
 */
String Tag::Name() const {
  std::ostringstream os;
  
  // First the species name:
  os << toShortName(Isotopologues[spec_ind].spec) << "-";
  
  // Is this a CIA tag?
  if (type == TagType::Cia) {
    os << "CIA-" << toShortName(cia_2nd_species) << "-" << cia_dataset_index;
    
  } else if (type == TagType::FreeElectrons || type == TagType::Particles) {
    os << toShortName(Isotopologues[spec_ind].spec);
  }
  // Hitran Xsec flag.
  else if (type == TagType::HitranXsec) {
    os << "HXSEC";
  } else {
    // Zeeman flag.
    if (type == TagType::Zeeman) os << "Z-";
    
    // Now the isotopologue. Can be a single isotopologue or ALL.
    os << Isotopologues[spec_ind].isotname << '-';
    
    // Now the frequency limits, if there are any. For this we first
    // need to determine the floating point precision.
    
    // Determine the precision, depending on whether Numeric is double
    // or float:
    int precision;
    #ifdef USE_FLOAT
    precision = FLT_DIG;
    #else
    #ifdef USE_DOUBLE
    precision = DBL_DIG;
    #else
    #error Numeric must be double or float
    #endif
    #endif
    
    if (0 > lower_freq) {
      // lower_freq < 0 means no lower limit.
      os << "*-";
    } else {
      os << std::setprecision(precision);
      os << lower_freq << "-";
    }
    
    if (0 > upper_freq) {
      // upper_freq < 0 means no upper limit.
      os << "*";
    } else {
      os << std::setprecision(precision);
      os << upper_freq;
    }
  }
  return os.str();
}
}

ArrayOfSpeciesTag2::ArrayOfSpeciesTag2(String names) {
  // There can be a comma separated list of tag definitions, so we
  // need to break the String apart at the commas.
  ArrayOfString tag_def;
  
  bool go_on = true;
  String these_names = names;
  while (go_on) {
    //          Index n = find_first( these_names, ',' );
    Index n = these_names.find(',');
    if (n == these_names.npos)  // Value npos indicates not found.
    {
      // There are no more commas.
      //              cout << "these_names: (" << these_names << ")\n";
      tag_def.push_back(these_names);
      go_on = false;
    } else {
      tag_def.push_back(these_names.substr(0, n));
      these_names.erase(0, n + 1);
    }
  }
  // tag_def now holds the different tag Strings for this group.
  
  // Set size to zero, in case the method has been called before.
  resize(0);

  for (Index s = 0; s < tag_def.nelem(); ++s) {
    Species::Tag this_tag(tag_def[s]);
    
    // Safety checks:
    if (s > 0) {
      // Tags inside a group must belong to the same species.
      ARTS_USER_ERROR_IF (front().Isotopologue().spec != this_tag.Isotopologue().spec,
                          "Tags in a tag group must belong to the same species.");
      
      // Zeeman tags and plain line by line tags must not be mixed. (Because
      // there can be only one line list per tag group.)
      ARTS_USER_ERROR_IF (((front().type == Species::TagType::Zeeman) &&
                          (this_tag.type == Species::TagType::Plain)) ||
                          ((front().type == Species::TagType::Plain) &&
                          (this_tag.type == Species::TagType::Zeeman)),
                                              "Zeeman tags and plain line-by-line tags must "
                                              "not be mixed in the same tag group.");
    }
    
    push_back(this_tag);
  }
}
