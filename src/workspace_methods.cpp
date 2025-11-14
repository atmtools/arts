#include "workspace_methods.h"

#include <mystring.h>

#include <exception>
#include <format>
#include <limits>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string>

#include "workspace_agendas.h"
#include "workspace_meta_methods.h"
#include "workspace_variables.h"

#if defined(__clang__)
#pragma clang optimize off
#endif

bool WorkspaceMethodInternalRecord::has_any() const {
  const auto cmp = Cmp::eq("Any");
  return std::ranges::any_of(gout_type, cmp) +
         std::ranges::any_of(gin_type, cmp);
}

bool WorkspaceMethodInternalRecord::has_overloads() const {
  for (auto& str : gout_type) {
    if (std::any_of(str.begin(), str.end(), Cmp::eq(','))) return true;
  }

  for (auto& str : gin_type) {
    if (std::any_of(str.begin(), str.end(), Cmp::eq(','))) return true;
  }

  return false;
}

std::string WorkspaceMethodInternalRecord::docstring() const try {
  std::string doc = std::format("/** {}\n", desc);
  if (pass_workspace) doc += "  @param[in] ws Workspace reference\n";

  for (auto& str : out) {
    if (std::any_of(
            in.begin(), in.end(), [&str](auto& var) { return str == var; }))
      doc += std::format("  @param[inout] {} As WSV\n", str);
    else
      doc += std::format("  @param[out] {} As WSV\n", str);
  }

  for (std::size_t i = 0; i < gout.size(); i++) {
    doc += std::format("  @param[out] {} {}\n", gout[i], gout_desc[i]);
  }

  for (auto& str : in) {
    if (std::any_of(
            out.begin(), out.end(), [&str](auto& var) { return str == var; }))
      continue;
    doc += std::format("  @param[in] {} As WSV\n", str);
  }

  for (std::size_t i = 0; i < gin.size(); i++) {
    doc += std::format("  @param[in] {} {}\n", gin[i], gin_desc[i]);
  }

  doc += " */";

  return doc;
} catch (std::exception& e) {
  throw std::runtime_error("Error in meta-function docstring():\n\n" +
                           std::string(e.what()));
}

std::vector<std::vector<std::string>>
WorkspaceMethodInternalRecord::generic_overloads() const {
  const auto cmp = Cmp::eq(',');

  std::vector<std::vector<std::string>> genvar;

  for (auto& str : gout_type) {
    if (std::any_of(str.begin(), str.end(), cmp))
      genvar.push_back(split(str, ","));
    else
      genvar.push_back(std::vector<std::string>{str});
  }

  for (auto& str : gin_type) {
    if (std::any_of(str.begin(), str.end(), cmp))
      genvar.push_back(split(str, ","));
    else
      genvar.push_back(std::vector<std::string>{str});
  }

  for (auto& v : genvar) {
    for (auto& s : v) {
      trim(s);
    }
  }

  return genvar;
}

std::string WorkspaceMethodInternalRecord::header(const std::string& name,
                                                  int overload) const try {
  const auto& wsv = internal_workspace_variables();
  const auto& wsa = internal_workspace_agendas();

  const std::string spaces(name.size() + return_type.size() + 2, ' ');

  const auto overloads = generic_overloads();
  int GVAR             = 0;

  std::string doc{};

  if (name.size() > 7 and wsa.contains(name.substr(0, name.size() - 7)) and
      name.substr(name.size() - 7) == "Execute") {
    return doc;
  }

  if (has_any()) {
    doc += "template <WorkspaceGroup T>\n";
  }

  doc += std::format("{} {}(", return_type, name);

  bool first = true;
  if (pass_workspace) {
    doc   += "const Workspace& ws";
    first  = false;
  }

  for (auto& str : out) {
    auto wsv_ptr = wsv.find(str);
    if (wsv_ptr != wsv.end()) {
      doc += std::format(
          "{}{}& {}", comma(first, spaces), wsv_ptr->second.type, str);
      continue;
    }

    auto wsa_ptr = wsa.find(str);
    if (wsa_ptr != wsa.end()) {
      doc += std::format("{}Agenda& {}", comma(first, spaces), str);
      continue;
    }

    throw std::runtime_error("WorkspaceMethodInternalRecord::header " + name);
  }

  for (const auto& i : gout) {
    doc += std::format(
        "{}{}& {}",
        comma(first, spaces),
        any_is_typename(overloads[GVAR][std::min<int>(
            overload, static_cast<int>(overloads[GVAR].size() - 1))]),
        i);
    GVAR++;
  }

  for (auto& str : in) {
    if (std::any_of(
            out.begin(), out.end(), [&str](auto& var) { return str == var; }))
      continue;

    auto wsv_ptr = wsv.find(str);
    if (wsv_ptr != wsv.end()) {
      doc += std::format(
          "{}const {}& {}", comma(first, spaces), wsv_ptr->second.type, str);
      continue;
    }

    throw std::runtime_error("WorkspaceMethodInternalRecord::header " + name);
  }

  for (const auto& i : gin) {
    doc += std::format(
        "{}const {}& {}",
        comma(first, spaces),
        any_is_typename(overloads[GVAR][std::min<int>(
            overload, static_cast<int>(overloads[GVAR].size() - 1))]),
        i);
    GVAR++;
  }

  doc += ")";

  return doc;
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format(R"(Error in meta-function header("{}", {}):

{}
)",
                  name,
                  overload,
                  std::string_view(e.what())));
}

int WorkspaceMethodInternalRecord::count_overloads() const try {
  const auto overloads = generic_overloads();
  int g                = 1;
  for (auto& x : overloads) {
    int ng = static_cast<int>(x.size());
    if (ng == 1) continue;
    if (g != ng and g != 1)
      throw std::runtime_error("Inconsistent number of overloads");
    g = ng;
  }
  return g;
} catch (std::exception& e) {
  throw std::runtime_error("Error in meta-function count_overloads():\n\n" +
                           std::string(e.what()));
}

std::string WorkspaceMethodInternalRecord::call(const std::string& name) const
    try {
  const std::string spaces(return_type.size() + 2, ' ');

  std::string doc = std::format("  {}(\n      ", name);

  bool first = true;
  if (pass_workspace) {
    doc   += "ws";
    first  = false;
  }

  for (auto& str : out) {
    doc += std::format("{}{}", comma(first, spaces), str);
  }

  for (const auto& i : gout) {
    doc += std::format("{}{}", comma(first, spaces), i);
  }

  for (auto& str : in) {
    if (std::any_of(
            out.begin(), out.end(), [&str](auto& var) { return str == var; }))
      continue;

    doc += std::format("{}{}", comma(first, spaces), str);
  }

  for (const auto& i : gin) {
    doc += std::format("{}{}", comma(first, spaces), i);
  }

  doc += ");";

  return doc;
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("Cannot create call of method \"{}\":\n\n{}",
                  name,
                  std::string_view(e.what())));
}

namespace {
void add_agenda_methods(
    std::unordered_map<std::string, WorkspaceMethodInternalRecord>& wsm_data) {
  for (auto& [agname, ag] : internal_workspace_agendas()) {
    if (wsm_data.contains(agname + "Execute") or
        wsm_data.contains(agname + "Set") or
        wsm_data.contains(agname + "ExecuteOperator") or
        wsm_data.contains(agname + "SetOperator"))
      throw std::runtime_error(std::format(R"(

The following methods are predefined.  Do not manually define these:

{0}Execute
{0}ExecuteOperator
{0}Set
{0}SetOperator

Remove the manual definition of these methods from workspace_methods.cpp.
)",
                                           agname));

    std::vector<std::string> input = ag.input;
    input.push_back(agname);
    std::vector<std::string> input_op = ag.input;
    input_op.push_back(agname + "_operator");

    wsm_data[agname + "Execute"] = {
        .desc        = "Executes *" + agname + "*, see it for more details\n",
        .author      = {"``Automatically Generated``"},
        .return_type = "Workspace",
        .return_desc = "The internal workspace, cleaned from all input/output.",
        .out         = ag.output,
        .in          = input,
        .pass_workspace = true,
    };

    wsm_data[agname + "ExecuteOperator"] = {
        .desc = std::format(
            "Executes an operator emulating *{0}*, see it, and also *{0}Operator*, for more details\n",
            agname),
        .author    = {"``Automatically Generated``"},
        .out       = ag.output,
        .in        = ag.input,
        .gin       = {agname + "_operator"},
        .gin_type  = {agname + "Operator"},
        .gin_value = {std::nullopt},
        .gin_desc  = {"Operator for *" + agname + "*"},
    };

    wsm_data[agname + "SetOperator"] = {
        .desc = std::format(
            "Set *{0}* to exclusively use provided external operator.  See *{0}Operator* for more details.\n",
            agname),
        .author    = {"``Automatically Generated``"},
        .out       = {agname},
        .gin       = {"f"},
        .gin_type  = {agname + "Operator"},
        .gin_value = {std::nullopt},
        .gin_desc  = {"Operator for *" + agname + "*"},
    };

    if (ag.enum_options.empty()) continue;

    wsm_data[agname + "Set"] = {
        .desc      = "Set *" + agname + "* to a specific predefined option\n",
        .author    = {"``Automatically Generated``"},
        .out       = {agname},
        .gin       = {"option"},
        .gin_type  = {"String"},
        .gin_value = {ag.enum_default.empty()
                          ? std::nullopt
                          : std::optional<String>(ag.enum_default)},
        .gin_desc  = {"Choice of generated agenda"},
    };
  }
}

void fix(
    std::unordered_map<std::string, WorkspaceMethodInternalRecord>& wsm_data) {
  for (auto& wsmr : wsm_data | std::views::values) {
    if (wsmr.desc.back() != '\n') wsmr.desc += '\n';
  }
}

std::unordered_map<std::string, WorkspaceMethodInternalRecord>
internal_workspace_methods_create() try {
  std::unordered_map<std::string, WorkspaceMethodInternalRecord> wsm_data;

  wsm_data["Ignore"] = {
      .desc      = R"--(Ignore a workspace variable.

This method is handy for use in agendas in order to suppress warnings
about unused input workspace variables. What it does is: Nothing!
In other words, it just ignores the variable it is called on.

This method can ignore any workspace variable you want.
)--",
      .author    = {"Stefan Buehler"},
      .gin       = {"input"},
      .gin_type  = {"Any"},
      .gin_value = {std::nullopt},
      .gin_desc  = {R"--(Variable to be ignored.)--"},
  };

  wsm_data["ReadXML"] = {
      .desc      = R"--(Reads a workspace variable from an XML file.

This method can read variables of any group.

If the given filename does not exist, this method will
also look for files with an added .xml, .xml.gz and .gz extension.

.. note::
    ARTS groups, including those that are not workspace groups,
    generally have two methods called ``fromxml`` and ``readxml``
    that you can access directly from python.  It is often more
    convenient and better to use these methods directly instead
    of using this workspace method.  It exists mainly for
    completeness and for use in agendas.
)--",
      .author    = {"Oliver Lemke"},
      .gout      = {"output"},
      .gout_type = {"Any"},
      .gout_desc = {R"--(Variable to be read.)--"},
      .gin       = {"filename"},
      .gin_type  = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc  = {R"--(Name of the XML file.)--"},
  };

  wsm_data["ReadXMLIndexed"] = {
      .desc      = R"--(As *ReadXML*, but reads indexed file names.

The variable is read from a file with name::

   <filename>.<file_index>.xml.

where <file_index> is the value of ``file_index``.

This means that ``filename`` shall here not include the .xml
extension.

.. note::
    ARTS groups, including those that are not workspace groups,
    generally have two methods called ``fromxml`` and ``readxml``
    that you can access directly from python.  It is often more
    convenient and better to use these methods directly instead
    of using this workspace method.  It exists mainly for
    completeness and for use in agendas.
)--",
      .author    = {"Oliver Lemke"},
      .gout      = {"output"},
      .gout_type = {"Any"},
      .gout_desc = {R"--(Workspace variable to be read.)--"},
      .gin       = {"file_index", "filename", "digits"},
      .gin_type  = {"Index", "String", "Index"},
      .gin_value = {std::nullopt, std::nullopt, Index{0}},
      .gin_desc =
          {R"--(Index of the file to read.)--",
           R"--(File name. See above.)--",
           R"--(Equalize the widths of all numbers by padding with zeros as necessary. 0 means no padding (default).)--"},
  };

  wsm_data["abs_xfit_dataReadSpeciesSplitCatalog"] = {
      .desc      = R"--(Reads HITRAN Crosssection coefficients

Reads coefficient files for HITRAN Xsec species defined
in *abs_species*.

.. tip::
    A common and perhaps more convenient alternative to this method
    is *ReadCatalogData*.

    If you have downloaded
    the ARTS catalog data - ``arts-cat-data`` - and set the environment
    variable ``ARTS_DATA_PATH`` to point to the location of this data,
    you can use that method to automagically read the data more easily
    than calling this method directly.
)--",
      .author    = {"Oliver Lemke"},
      .out       = {"abs_xfit_data"},
      .in        = {"abs_species"},
      .gin       = {"basename", "ignore_missing"},
      .gin_type  = {"String", "Index"},
      .gin_value = {std::nullopt, Index{0}},
      .gin_desc  = {R"--(Basepath to the files)--",
                    R"--(Ignore missing files (0: no, 1: yes))--"},
  };

  wsm_data["Touch"] = {
      .desc      = R"--(As *Ignore* but for agenda output.

This method is handy for use in agendas in order to suppress
warnings about not-produced output workspace variables.

What it does, in case the variable is initialized already, is:
Nothing!
In case the variable is not yet initialized, it is default initialized.
)--",
      .author    = {"Oliver Lemke"},
      .gout      = {"input"},
      .gout_type = {"Any"},
      .gout_desc = {R"--(Variable to do nothing with.)--"},
  };

  wsm_data["WignerInit"] = {
      .desc      = R"--(Initialize the Wigner tables

The default values take about 1 Gb memory.

The static data is kept in an external library and is therefore
only available inside ARTS.  Nevertheless, this must be set by
the application because any default value might be too small or
too large for the needs of any one application.

We rely on the Wigner library by :cite:t:`johansson2016`.  This method
allows mimicking the initialization of that library.

.. tip::

    If you are starved for memory, you can call *WignerUnload* after
    you are done with all Wigner calculations.  This will free the memory.
)--",
      .author    = {"Richard Larsson"},
      .gin       = {"fast_wigner_stored_symbols",
                    "largest_wigner_symbol_parameter",
                    "symbol_type"},
      .gin_type  = {"Index", "Index", "Index"},
      .gin_value = {Index{20000000}, Index{250}, Index{6}},
      .gin_desc =
          {R"--(Number of stored symbols possible before replacements)--",
           R"--(Largest symbol used for initializing factorials (e.g., largest J or L))--",
           "Type of symbol (3 or 6)"},
  };

  wsm_data["WignerUnload"] = {
      .desc =
          R"--(Unloads the Wigner tables from static data (see *WignerInit*)
)--",
      .author = {"Richard Larsson"},
  };

  wsm_data["WriteBuiltinPartitionFunctionsXML"] = {
      .desc      = R"--(Writes all the builtin partition functions to file.

All available partition functions are written to files in the select format
in the select directory

The temperature will be linearly spaced between [Tlow, Tupp] with N values

See *FileType* for valid ``output_file_format``.
)--",
      .author    = {"Richard Larsson"},
      .gin       = {"output_file_format", "dir", "Tlow", "Tupp", "N"},
      .gin_type  = {"String", "String", "Numeric", "Numeric", "Index"},
      .gin_value = {String{"ascii"},
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt},
      .gin_desc  = {"The format of the output",
                    R"--(The directory to write the data towards)--",
                    R"--(The lowest temperature)--",
                    R"--(The highest temperature)--",
                    R"--(The number of temperature points)--"},
  };

  wsm_data["WriteXML"] = {
      .desc      = R"--(Writes a workspace variable to an XML file.

This method can write variables of any group.

If no_clobber is set to 1, an increasing number will be
appended to the filename if the file already exists.

See *FileType* for valid ``output_file_format``.

.. note::
    ARTS groups, including those that are not workspace groups,
    generally have a method called ``savexml``
    that you can access directly from python.  It is often more
    convenient and better to use this method directly instead
    of using this workspace method.  It exists mainly for
    completeness and for use in agendas.
)--",
      .author    = {"Oliver Lemke"},
      .gin       = {"output_file_format", "input", "filename", "no_clobber"},
      .gin_type  = {"String", "Any", "String", "Index"},
      .gin_value = {String("ascii"), std::nullopt, std::nullopt, Index{0}},
      .gin_desc =
          {"The format of the output",
           R"--(Variable to be saved.)--",
           R"--(Name of the XML file.)--",
           R"--(0: Overwrite existing files, 1: Use unique filenames)--"},
  };

  wsm_data["WriteXMLIndexed"] = {
      .desc   = R"--(As *WriteXML*, but creates indexed file names.

The variable is written to a file with name::

  <filename>.<file_index>.xml.

where <file_index> is the value of ``file_index``.

This means that ``filename`` shall here not include the .xml
extension.

See *FileType* for valid ``output_file_format``.

.. note::
    ARTS groups, including those that are not workspace groups,
    generally have a method called ``savexml``
    that you can access directly from python.  It is often more
    convenient and better to use this method directly instead
    of using this workspace method.  It exists mainly for
    completeness and for use in agendas.
)--",
      .author = {"Patrick Eriksson, Oliver Lemke"},
      .gin =
          {"output_file_format", "file_index", "input", "filename", "digits"},
      .gin_type = {"String", "Index", "Any", "String", "Index"},
      .gin_value =
          {String("ascii"), std::nullopt, std::nullopt, std::nullopt, Index{0}},
      .gin_desc =
          {"The format of the output",
           R"--(Index number for files.)--",
           R"--(Workspace variable to be saved.)--",
           R"--(File name. See above.)--",
           R"--(Equalize the widths of all numbers by padding with zeros as necessary. 0 means no padding (default).)--"},
  };

  wsm_data["abs_cia_dataReadFromCIA"] = {
      .desc =
          R"--(Read data from a CIA data file for all CIA molecules defined
in *abs_species*.

The units in the HITRAN file are:

- Frequency: :math:`\textrm{cm}^{-1}`
- Binary absorption cross-section: :math:`\textrm{cm}^{5} \, \textrm{molec}^{-2}`

Upon reading we convert this to the ARTS internal SI units 
of Hz and :math:`\textrm{m}^{5} \, \textrm{molec}^{-2}`.
)--",
      .author    = {"Oliver Lemke"},
      .out       = {"abs_cia_data"},
      .in        = {"abs_species"},
      .gin       = {"catalogpath"},
      .gin_type  = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc  = {R"--(Path to the CIA catalog directory.)--"},
  };

  wsm_data["abs_cia_dataReadFromXML"] = {
      .desc =
          R"--(Read data from a CIA XML file and check that all CIA tags defined
in *abs_species* are present in the file.

The units of the data are described in *abs_cia_dataReadFromCIA*.
)--",
      .author    = {"Oliver Lemke"},
      .out       = {"abs_cia_data"},
      .in        = {"abs_species"},
      .gin       = {"filename"},
      .gin_type  = {"String"},
      .gin_value = {String("")},
      .gin_desc  = {R"--(Name of the XML file.)--"},
  };

  wsm_data["abs_cia_dataReadSpeciesSplitCatalog"] = {
      .desc      = R"--(Reads a species split CIA dataset.

The file names are expected to be of the form:

   <basename><Spec1>-CIA-<Spec2>.xml

where <Spec1> and <Spec2> are the *SpeciesEnum* names of the two species
involved in the CIA.

.. tip::
    A common and perhaps more convenient alternative to this method
    is *ReadCatalogData*.

    If you have downloaded
    the ARTS catalog data - ``arts-cat-data`` - and set the environment
    variable ``ARTS_DATA_PATH`` to point to the location of this data,
    you can use that method to automagically read the data more easily
    than calling this method directly.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_cia_data"},
      .in        = {"abs_species"},
      .gin       = {"basename", "ignore_missing"},
      .gin_type  = {"String", "Index"},
      .gin_value = {std::nullopt, Index{0}},
      .gin_desc =
          {R"--(The path to the split catalog files)--",
           R"--(Flag to continue in case nothing is found [0 throws, 1 continues])--"},
  };

  wsm_data["abs_predef_dataReadSpeciesSplitCatalog"] = {
      .desc =
          R"--(Reads *abs_predef_data* catalog but only for *abs_species*


The file names are expected to be of the form:

   <basename><Spec>-<Model>.xml

where <Spec> is the *SpeciesEnum* names of the species
and <Model> is the model name.  See *abs_speciesSet* for more
information on how to define a species with a predefined model.

If ``name_missing`` is true, missing models are set to named model, which is the most
common form of a predefined model.

.. tip::
    A common and perhaps more convenient alternative to this method
    is *ReadCatalogData*.

    If you have downloaded
    the ARTS catalog data - ``arts-cat-data`` - and set the environment
    variable ``ARTS_DATA_PATH`` to point to the location of this data,
    you can use that method to automagically read the data more easily
    than calling this method directly.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_predef_data"},
      .in        = {"abs_species"},
      .gin       = {"basename", "name_missing", "ignore_missing"},
      .gin_type  = {"String", "Index", "Index"},
      .gin_value = {std::nullopt, Index{1}, Index{0}},
      .gin_desc =
          {R"--(The path to the split catalog files)--",
           R"--(Flag to name models that are missing)--",
           R"--(Flag to otherwise (if not name_missing is true) ignore missing models)--"},
  };

  wsm_data["abs_speciesDefineAll"] = {
      .desc   = R"--(Sets *abs_species* to contain all species in ARTS
)--",
      .author = {"Richard Larsson"},
      .out    = {"abs_species"},
  };

  wsm_data["abs_speciesSet"] = {
      .desc      = R"--(Set *abs_species* to the named species.

The species that are defined by this method are used in various
file reading routines to populate both atmosphere- and absorption-related
data variables.  They select what data is required,

A tag begins with a valid *SpeciesEnum*.  The rest of the tag
is optional.

There 3 single tags and 2 combinatory tags.
These are named:

#. Joker.  Example ``"H2O"``.  Selects all Normal Isotopologue of the provided *SpeciesEnum*.
#. Normal Isotopologue.  Example ``"H2O-161"``.  Selects a specific Normal Isotopologue.
#. Predefined Model.  Example ``"H2O-PWR2022"``.  Selects a specific Predefined Model.
   For more information on Predefined Models, see *spectral_propmatAddPredefined*.
#. CIA.  Example ``"H2O-CIA-H2O"``.  Selects collusion-induced absorption between the two species.
   Any two *SpeciesEnum* in combination is valid.  The reverse combination is also valid and unique.
#. XFIT.  Example ``"H2O-XFIT"``.  Selects using cross-section fits for the species.

.. tip::
    We provide data to help make use of these tags in ``arts-cat-data``. Far
    from all use-cases are covered, but enough to help you get started.  See
    :func:`~pyarts3.data.download_arts_cat_data` for help to download and set it up.
)--",
      .author    = {"Stefan Buehler", "Richard Larsson"},
      .out       = {"abs_species"},
      .gin       = {"species"},
      .gin_type  = {"ArrayOfString"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {R"--(Specify one String for each tag group that you want to create. Inside the String, separate the tags by commas (plus optional blanks).)--"},
  };

  wsm_data["atm_fieldIGRF"] = {
      .desc      = R"--(Use IGRF to compute the magnetic field at each point.

The IGRF model is a model of the Earth's magnetic field. It is based on
spherical harmonics and is only valid for a limited time period.

The IGRF model is available via :cite:t:`Alken2021`.

.. note::
    The IGRF model is added as a functional object to the atmospheric field.
    If you need to retrieve the magnetic field, you must convert it to another type,
    use *atm_fieldSchmidthFieldFromIGRF*, which have been designed
    to support retrievals of the magnetic field via Legendre coefficients, or
    use *atm_fieldAbsoluteMagneticField*, which allow returning the magnitude
    of the magnetic field.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field"},
      .gin       = {"time"},
      .gin_type  = {"Time"},
      .gin_value = {Time{}},
      .gin_desc  = {"Time of data to use"},
  };

  wsm_data["atm_fieldSchmidthFieldFromIGRF"] = {
      .desc =
          R"--(For forward calculations, this should be similar to *atm_fieldIGRF*.

What it does different is that it is 1) not a direct computations matching the IGRF field,
instead averaging the Legendre coefficient matrices.

What it does very different is that it supports retrievals
of the magnetic field Legendre coefficients.

This is very much a WIP and not well tested.

The IGRF model is available via :cite:t:`Alken2021`.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field"},
      .gin       = {"time"},
      .gin_type  = {"Time"},
      .gin_value = {Time{}},
      .gin_desc  = {"Time of IGRF data to use"},
  };

  wsm_data["atm_fieldAbsoluteMagneticField"] = {
      .desc = R"--(Set the magnetic field to use the magnitude field functional.

The input field must be a *GeodeticField3* for all three parameters
to call this method.

The main purpose of this method is to retrieve the magnitude rather than the vector field.
)--",
      .author = {"Richard Larsson"},
      .out    = {"atm_field"},
      .in     = {"atm_field"},
  };

  wsm_data["atm_fieldAbsoluteWindField"] = {
      .desc   = R"--(Set the wind field to use the magnitude field functional.

The input field must be a *GeodeticField3* for all three parameters
to call this method.

The main purpose of this method is to retrieve the magnitude rather than the vector field.
)--",
      .author = {"Richard Larsson"},
      .out    = {"atm_field"},
      .in     = {"atm_field"},
  };

  wsm_data["freq_gridFitNonLTE"] = {
      .desc      = R"(Frequency grid useful for *atm_profileFitNonLTE*.

This method creates a frequency grid around the line-center of each
absorption line in the *abs_bands* variable.

.. note::
  For all sorted absorption line centers :math:`f_i`, the following should be 
  true:

  .. math::

    f_{i -1} + \delta f < f_i < f_{i + 1} - \delta f

  That is, the frequency ranges are not allowed to overlap.
)",
      .author    = {"Richard Larsson"},
      .out       = {"freq_grid"},
      .in        = {"abs_bands"},
      .gin       = {"df", "nf"},
      .gin_type  = {"Numeric", "Index"},
      .gin_value = {std::nullopt, Index{401}},
      .gin_desc =
          {R"--(Frequency grid around the line-center.  The range will cover :math:`f_i` \pm \left(1-\delta f\right)\forall i` of each absorption line :math:`i`, where this variable is :math:`\delta f`.)--",
           R"--(Number of frequency points per line.  The step between frequency grid points will be :math:`2\frac{\delta f}{N - 1}`, where this is :math:`N`.)--"},
  };

  wsm_data["freq_gridFromSingleFrequency"] = {
      .desc =
          R"(Composition method, creates a frequency grid from a single frequency.
)",
      .author = {"Richard Larsson"},
      .out    = {"freq_grid"},
      .in     = {"freq"},
  };

  wsm_data["atm_fieldFromProfile"] = {
      .desc =
          R"--(Sets the atmospheric field to be the 1D atmospheric profile.

The top of the atmosphere is the last value of the altitude grid.

All atmospheric points in the profile must contain the same set of parameters.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_profile", "alt_grid"},
      .gin       = {"altitude_extrapolation"},
      .gin_type  = {"InterpolationExtrapolation"},
      .gin_value = {InterpolationExtrapolation::Linear},
      .gin_desc  = {"Extrapolation method along the altitude grid"},
  };

  wsm_data["atm_profileExtract"] = {
      .desc =
          R"--(Extract an atmospheric profile from the atmospheric field.
)--",
      .author = {"Richard Larsson"},
      .out    = {"atm_profile"},
      .in     = {"atm_field", "alt_grid", "lat", "lon"},
  };

  wsm_data["atm_profileFromGrid"] = {
      .desc =
          R"--(Extract an atmospheric profile and its grids.

The key is used to find a *GeodeticField3* in the atmospheric field.  Its grids
must form a profile.  The profile is extracted and returned.  The grids are
returned as well.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_profile", "alt_grid", "lat", "lon"},
      .in        = {"atm_field"},
      .gin       = {"key"},
      .gin_type  = {"AtmKey"},
      .gin_value = {AtmKey::t},
      .gin_desc = {"Key to find the *GeodeticField3* in the atmospheric field"},
  };

  wsm_data["atm_profileFitNonLTE"] = {
      .desc =
          R"--(Fits non-LTE distributions to the level data.

The spectral flux is computed from the pseudo-2D assumption.

This method fits non-LTE distributions to the level data in the
atmospheric field.  It only works for absorption band data that
is separated by single-lines-per-band, and will produce nonsense
for overlapping line data.  If the lines overlap, the method will
keep introducing more-and-more energy into the system, meaning
that the method will not converge or turn to some extreme stable
state.

The statistical equilibrium equation is given by finding valid set of energy level distribution
:math:`n` such that for all valid energy level combination of upper levels :math:`i` and
lower levels :math:`j` the rate of change is zero for some :math:`n` that satisfies the equation

.. math::

    \frac{d n_i}{dt} =
        \sum_{j > i} \left[ n_j A_{ji} - \left( n_i B_{ij} - n_j B_{ji} \right) J_{ij} \right]
      - \sum_{j < i} \left[ n_i A_{ij} - \left( n_j B_{ji} - n_i B_{ij} \right) J_{ij} \right]
      + \sum_{j} \left[ n_j C_{ji} - n_i C_{ij} \right],

where :math:`A_{ij}` is the spontaneous emission rate, :math:`B_{ij}` is the
stimulated emission rate, :math:`B_{ij}` is the photon absorption rate,
:math:`J_{ij}` is the line-integrated flux, and :math:`C_{ij}`
is the collisional rate.

Generally, you need :math:`n` to compute :math:`J_{ij}`, making the problem non-linear.
Thus an iterative process is used to find the solution.  The iteration is considered
converged when the relative change in the energy level distribution is below the
convergence criterion.  Alternatively, the iteration is halted if the iteration count limit
is breached.

The method used here is based on :cite:t:`Yamada2018`
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_profile"},
      .in        = {"atm_profile",
                    "abs_bands",
                    "spectral_propmat_agenda",
                    "surf_field",
                    "freq_grid",
                    "alt_grid",
                    "lat",
                    "lon"},
      .gin       = {"collision_data",
                    "levels",
                    "pol",
                    "azimuth",
                    "dza",
                    "convergence_limit",
                    "iteration_limit",
                    "consider_limb"},
      .gin_type  = {"QuantumIdentifierGriddedField1Map",
                    "ArrayOfQuantumLevelIdentifier",
                    "Stokvec",
                    "Numeric",
                    "Numeric",
                    "Numeric",
                    "Index",
                    "Index"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    Stokvec{1.0, 0.0, 0.0, 0.0},
                    Numeric{0.0},
                    Numeric{5.0},
                    Numeric{1e-6},
                    Index{100},
                    Index{1}},
      .gin_desc =
          {"Collision data for the transitions - for :math:`C_{ij}` and :math:`C_{ji}`",
           "The order of the energy levels",
           "The polarization selection vector (use the default unless you know what you are doing)",
           "The azimuth of the radiation",
           "The zenith angle limit for the internal call to *za_gridProfilePseudo2D*",
           "Convergence criterion for the energy level distribution",
           "Maximum number of iterations",
           "Whether to add extra limb points in *za_gridProfilePseudo2D*"},
      .pass_workspace = true,
  };

  wsm_data["atm_fieldInit"] = {
      .desc =
          R"--(Initialize the atmospheric field with some altitude and isotopologue ratios

See *IsoRatioOption* for valid ``default_isotopologue``.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .gin       = {"toa", "default_isotopologue"},
      .gin_type  = {"Numeric", "String"},
      .gin_value = {std::nullopt, String{"Builtin"}},
      .gin_desc  = {R"--(Top of atmosphere altitude [m].)--",
                    "Default option for the isotopologue ratios"},
  };

  wsm_data["atm_pointInit"] = {
      .desc =
          R"--(Initialize an atmospheric point with some isotopologue ratios

See *IsoRatioOption* for valid ``default_isotopologue``.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_point"},
      .gin       = {"default_isotopologue"},
      .gin_type  = {"String"},
      .gin_value = {String{"Builtin"}},
      .gin_desc  = {"Default option for the isotopologue ratios"},
  };

  wsm_data["spectral_tramat_bkgFromPathPropagationBack"] = {
      .desc =
          R"--(Sets *spectral_tramat_bkg* to back of *spectral_tramat_cumulative_path*.

This is purely compositional and it is better to use pure python code if need this functionality
in your own control-flow.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_tramat_bkg"},
      .in     = {"spectral_tramat_cumulative_path"},
  };

  wsm_data["spectral_tramat_bkgFromPathPropagationFront"] = {
      .desc =
          R"--(Sets *spectral_tramat_bkg* to front of *spectral_tramat_cumulative_path*.

This is purely compositional and it is better to use pure python code if need this functionality
in your own control-flow.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_tramat_bkg"},
      .in     = {"spectral_tramat_cumulative_path"},
  };

  wsm_data["ray_path_zeeman_magnetic_fieldFromPath"] = {
      .desc      = R"--(Sets a path of Zeeman effect magnetic field properties.

This will return a list of magnetic field properties along the path.
The magnetic properties in Zeeman coordinates are the absolute strength [H],
the angle between the magnetic field and the line of sight [theta], and the
the rotation of the magnetic field in the plane perpendicular to the line of
sight [eta].

This is mostly a convenience method to allow diagnostic plotting of the
magnetic field along the path.
)--",
      .author    = {"Richard Larsson"},
      .gout      = {"ray_path_zeeman_magnetic_field"},
      .gout_type = {"ArrayOfVector3"},
      .gout_desc = {R"--(Along-the-path [H, theta, eta])--"},
      .in        = {"ray_path", "atm_path"},
  };

  wsm_data["abs_ecs_dataAddMakarov2020"] = {
      .desc   = R"--(Sets the O2-66 microwave band data for ECS.

This is based on the work of :cite:t:`Makarov2020`.
)--",
      .author = {"Richard Larsson"},
      .out    = {"abs_ecs_data"},
      .in     = {"abs_ecs_data"},
  };

  wsm_data["abs_ecs_dataAddTran2011"] = {
      .desc   = R"--(Sets the CO2-626, CO2-628, and CO2-636 band data for ECS.

Sets CO2 species.

This is based on the work of :cite:t:`Tran2011`.
)--",
      .author = {"Richard Larsson"},
      .out    = {"abs_ecs_data"},
      .in     = {"abs_ecs_data"},
  };

  wsm_data["abs_ecs_dataAddRodrigues1997"] = {
      .desc   = R"--(Sets the CO2-626, CO2-628, and CO2-636 band data for ECS.

Sets N2 and O2 species.

This is based on the work of :cite:t:`Rodrigues1997`.
)--",
      .author = {"Richard Larsson"},
      .out    = {"abs_ecs_data"},
      .in     = {"abs_ecs_data"},
  };

  wsm_data["abs_ecs_dataInit"] = {
      .desc   = R"--(Resets/initializes the ECS data.
)--",
      .author = {"Richard Larsson"},
      .out    = {"abs_ecs_data"},
  };

  wsm_data["abs_ecs_dataAddMeanAir"] = {
      .desc      = R"--(Sets ECS data for air from other data if available.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_ecs_data"},
      .in        = {"abs_ecs_data"},
      .gin       = {"vmrs", "species"},
      .gin_type  = {"Vector", "ArrayOfSpeciesEnum"},
      .gin_value = {std::nullopt, std::nullopt},
      .gin_desc  = {R"--(VMRs of air species)--", R"--(Air species)--"},
  };

  wsm_data["freq_gridWindShift"] = {
      .desc =
          R"--(Applies wind shift to the *freq_grid* for the local frequency grid.

Also sets *freq_wind_shift_jac*.

If the wind is 0 or nan, the *freq_grid* remains unchanged.
)--",
      .author = {"Richard Larsson"},
      .out    = {"freq_grid", "freq_wind_shift_jac"},
      .in     = {"freq_grid", "atm_point", "ray_point"},
  };

  wsm_data["freqWindShift"] = {
      .desc =
          R"--(Same as *freq_gridWindShift* but for single frequency values.
)--",
      .author = {"Richard Larsson"},
      .out    = {"freq", "freq_wind_shift_jac"},
      .in     = {"freq", "atm_point", "ray_point"},
  };

  wsm_data["atm_pathFromProfile"] = {
      .desc =
          R"--(Set ``atm_path = atm_profile``.

This is purely compositional and it is better to use pure python code if need this functionality
in your own control-flow.
)--",
      .author = {"Richard Larsson"},
      .out    = {"atm_path"},
      .in     = {"atm_profile"},
  };

  wsm_data["atm_pathFromPath"] = {
      .desc   = R"--(Gets the atmospheric points along the path.
)--",
      .author = {"Richard Larsson"},
      .out    = {"atm_path"},
      .in     = {"ray_path", "atm_field"},
  };

  wsm_data["spectral_tramat_cumulative_pathFromPath"] = {
      .desc =
          R"--(Sets *spectral_tramat_cumulative_path* by forward iteration of *spectral_tramat_path*
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_tramat_cumulative_path"},
      .in     = {"spectral_tramat_path"},
  };

  wsm_data["freq_grid_pathFromPath"] = {
      .desc   = R"--(Gets the frequency grids along the path.

The derivative transformation for wind parameters is also returned.

See *spectral_propmat_jacWindFix* for use of the wind shift data.
)--",
      .author = {"Richard Larsson"},
      .out    = {"freq_grid_path", "freq_wind_shift_jac_path"},
      .in     = {"freq_grid", "ray_path", "atm_path"},
  };

  wsm_data["spectral_propmat_pathFromPath"] = {
      .desc =
          R"--(Gets the propagation matrix and non-LTE source term along the path.

The calculations are in parallel if the program is not in parallel already.

Also outputs the *freq_grid_path* as a side effect (of wind).
)--",
      .author         = {"Richard Larsson"},
      .out            = {"spectral_propmat_path",
                         "spectral_srcvec_nlte_path",
                         "spectral_propmat_jac_path",
                         "spectral_srcvec_nlte_jac_path"},
      .in             = {"spectral_propmat_agenda",
                         "freq_grid_path",
                         "freq_wind_shift_jac_path",
                         "jac_targets",
                         "ray_path",
                         "atm_path"},
      .pass_workspace = true,
  };

  wsm_data["spectral_propmat_path_species_splitFromPath"] = {
      .desc =
          R"--(As *spectral_propmat_pathFromPath* but the output is split between the species in the
*select_species_list*.

The outer dimension of the output arrays are the size of the species selection list.  The inner dimensions
are as per *spectral_propmat_pathFromPath*.
)--",
      .author    = {"Richard Larsson"},
      .gout      = {"spectral_propmat_path_species_split",
                    "spectral_srcvec_nlte_path_species_split",
                    "spectral_propmat_jac_path_species_split",
                    "spectral_srcvec_nlte_jac_path_species_split"},
      .gout_type = {"ArrayOfArrayOfPropmatVector",
                    "ArrayOfArrayOfStokvecVector",
                    "ArrayOfArrayOfPropmatMatrix",
                    "ArrayOfArrayOfStokvecMatrix"},
      .gout_desc =
          {R"--(Propagation matrix for selected species)--",
           R"--(Non-LTE source vector for selected species)--",
           R"--(Jacobian of propagation matrix for selected species)--",
           R"--(Jacobian of non-LTE source vector for selected species)--"},
      .in             = {"spectral_propmat_agenda",
                         "freq_grid_path",
                         "freq_wind_shift_jac_path",
                         "jac_targets",
                         "ray_path",
                         "atm_path",
                         "select_species_list"},
      .pass_workspace = true,
  };

  wsm_data["spectral_propmat_scat_pathFromPath"] = {
      .desc =
          R"--(Gets the propagation matrix for scattering along the path.

The calculations are in parallel if the program is not in parallel already.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_propmat_scat_path"},
      .in     = {"spectral_propmat_scat_agenda", "freq_grid_path", "atm_path"},
      .pass_workspace = true,
  };

  wsm_data["spectral_propmat_pathAddScattering"] = {
      .desc =
          R"--(Adds the scattering part of the propagation matrix to the rest along the path.

The calculations are in parallel if the program is not in parallel already.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_propmat_path"},
      .in     = {"spectral_propmat_path", "spectral_propmat_scat_path"},
  };

  wsm_data["spectral_rad_srcvec_pathAddScattering"] = {
      .desc =
          R"--(Adds the scattering part of the propagation matrix to the rest along the path.

The calculations are in parallel if the program is not in parallel already.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad_srcvec_path"},
      .in     = {"spectral_rad_srcvec_path",
                 "spectral_rad_scat_path",
                 "spectral_propmat_path"},
  };

  wsm_data["spectral_rad_srcvec_pathFromPropmat"] = {
      .desc   = R"--(Gets the source term along the path.

Per *Stokvec* element (single frequency, single temperature, single derivative target),
the source term is computed as:

.. math::
  \vec{J} = B(T, f) + \mathbf{K}^{-1} \vec{S}

and the Jacobian is computed as:

.. math::
  \frac{\partial \vec{J}}{\partial x} =
    \frac{\partial B(T, f)}{\partial x} -
    \mathbf{K}^{-1} \left(
    \frac{\partial \mathbf{K}}{\partial x} \mathbf{K}^{-1} -
    \frac{\partial \vec{S}}{\partial x}  \right)

where:

.. list-table::
  :widths: auto
  :align: left
  :header-rows: 1

  * - Variable
    - Extracted from ARTS parameter
    - Meaning
  * - :math:`\vec{J}`
    - *spectral_rad_srcvec_path*
    - The spectral radiance source term along the path.
  * - :math:`B(T, f)`
    - ``None`` - this is computed locally
    - The Planck function at the temperature and frequency.
  * - :math:`\mathbf{K}`
    - *spectral_propmat_path*
    - The propagation matrix along the path.
  * - :math:`\vec{S}`
    - *spectral_srcvec_nlte_path*
    - The non-LTE source vector along the path.
  * - :math:`\frac{\partial \vec{J}}{\partial x}`
    - *spectral_rad_srcvec_jac_path*
    - The Jacobian of the spectral radiance source term with respect to the
      *jac_targets*.
  * - :math:`\frac{\partial B(T, f)}{\partial x}`
    - ``None`` - this is computed locally
    - The Jacobian of the Planck function with respect to the *jac_targets*.  Only tempertature is supported.
  * - :math:`\frac{\partial \mathbf{K}}{\partial x}`
    - *spectral_propmat_jac_path*
    - The Jacobian of the propagation matrix with respect to the *jac_targets*.
  * - :math:`\frac{\partial \vec{S}}{\partial x}`
    - *spectral_srcvec_nlte_jac_path*
    - The Jacobian of the non-LTE source vector with respect to the *jac_targets*.
  * - :math:`x`
    - *jac_targets*
    - The targets for the Jacobian computation.
  * - :math:`T`
    - *atm_path*
    - The temperature at the atmospheric point along the path.
  * - :math:`f`
    - *freq_grid_path*
    - The frequency grid at the atmospheric point along the path.

The output dimensions are:

- *spectral_rad_srcvec_path*: *ray_path* x *freq_grid*
- *spectral_rad_srcvec_jac_path*: *ray_path* x *freq_grid* x *jac_targets* (target count)
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad_srcvec_path", "spectral_rad_srcvec_jac_path"},
      .in     = {"spectral_propmat_path",
                 "spectral_srcvec_nlte_path",
                 "spectral_propmat_jac_path",
                 "spectral_srcvec_nlte_jac_path",
                 "freq_grid_path",
                 "atm_path",
                 "jac_targets"},
  };

  wsm_data["abs_predef_dataAddWaterMTCKD400"] = {
      .desc      = R"--(Sets the data for MT CKD 4.0 Water model

Note that the vectors must have the same length, and that wavenumbers must be growing
at a constant rate.  The minimum length is 4.

Note also that as this is predefined model data, the units of the values of the vectors
must be as described by each vector.

This is based on the works cited here: https://hitran.org/mtckd/

.. note::
    The method itself is implemented from scratch.  Using any version of
    data after version 4.0 is supported by this method - all that changes
    are the values of the vectors. 
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_predef_data"},
      .in        = {"abs_predef_data"},
      .gin       = {"ref_temp",
                    "ref_press",
                    "ref_h2o_vmr",
                    "self_absco_ref",
                    "for_absco_ref",
                    "wavenumbers",
                    "self_texp"},
      .gin_type  = {"Numeric",
                    "Numeric",
                    "Numeric",
                    "Vector",
                    "Vector",
                    "Vector",
                    "Vector"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt},
      .gin_desc  = {R"--(Reference temperature)--",
                    R"--(Reference pressure)--",
                    R"--(Reference volume mixing ratio of water)--",
                    R"--(Self absorption [1/(cm-1 molecules/cm^2])--",
                    R"--(Foreign absorption [1/(cm-1 molecules/cm^2)])--",
                    R"--(Wavenumbers [cm-1])--",
                    R"--(Self temperature exponent [-])--"},
  };

  wsm_data["abs_predef_dataInit"] = {
      .desc   = R"--(Initialize the predefined model data
)--",
      .author = {"Richard Larsson"},
      .out    = {"abs_predef_data"},
  };

  wsm_data["spectral_propmatAddCIA"] = {
      .desc =
          R"--(Add absorption coefficients for HITRAN collision induced absorption (CIA).

This interpolates the cross sections from *abs_cia_data*.  If too few
temperature grid-points are available, its polynomial order of interpolation
decreases to the maximum allowed.  Otherwise, both frequency and temperature
are interpolated using third order polynomials.

Given that the interpolation is :math:`\vec{x}_{ij} = f\left(T, \vec{f}\right)`,
where the *atm_point* temperature is :math:`T` and :math:`f` is the
*freq_grid*, :math:`i` is an index into the *abs_cia_data* and :math:`j`
is an index into the underlying :class:`~pyarts3.arts.CIARecord` data structure,
the absorption coefficient from CIA is given by

.. math::

    vec{\alpha}_\mathbf{CIA} = \sum_i n_{i,0} n_{i,1} \sum_j vec{x}_{ij}

where :math:`n_{i,0}` and :math:`n_{i,1}` are number densities
of the two species involved in the CIA.

The input ``T_extrapolfac`` sets a limit on the interpolation along
the temperature grid of the data.  If the temperature grid is
:math:`[T_0, T_1, \cdots, T_{n-1}, T_n]`, then this method throws
an error if

.. math::

    T < T_0 - \Delta T_e \left(T_1 - T_0\right)

or

.. math::

    T > T_n + \Delta T_e \left(T_n - T_{n-1}\right)

where :math:`\Delta T_e` is the extrapolation factor given by
``T_extrapolfac``.  If this happens to you and you believe you can
use the computations anyways, set ``T_extrapolfac`` to a very large
or infinite value.

The frequnecy grid interpolation is limited to within the range
of the available data. Any point in *freq_grid* outside this
range will simply be ignored.  The frequency interpolation can
thus not fail.

.. note::
    ``ignore_errors`` can be set to 1 to suppress runtime errors, but any
    error will result in NaN values in the output.  This is useful if you
    want to debug your results, but not if you want to use them.
)--",
      .author    = {"Stefan Buehler, Oliver Lemke"},
      .out       = {"spectral_propmat", "spectral_propmat_jac"},
      .in        = {"spectral_propmat",
                    "spectral_propmat_jac",
                    "select_species",
                    "jac_targets",
                    "freq_grid",
                    "atm_point",
                    "abs_cia_data"},
      .gin       = {"T_extrapolfac", "ignore_errors"},
      .gin_type  = {"Numeric", "Index"},
      .gin_value = {Numeric{0.5}, Index{0}},
      .gin_desc =
          {R"--(Temperature extrapolation factor (relative to grid spacing). :math:`\Delta T_e` in text above.)--",
           R"--(Set to 1 to suppress runtime errors (and return NAN values instead).)--"},
  };

  wsm_data["spectral_propmatAddFaraday"] = {
      .desc   = R"--(Calculates absorption matrix describing Faraday rotation.

Faraday rotation is a change of polarization state of an
electromagnetic wave propagating through charged matter by
interaction with a magnetic field. Hence, this method requires that
the magnetic field is non-zero and that the electron density is
held by *atm_point* (*SpeciesEnum* : ``free_electrons``).

Faraday rotation affects Stokes parameters 2 and 3 (but not
intensity!).

Like all *spectral_propmat*-modifying methods, the method is additive,
i.e., does not overwrite the propagation matrix *spectral_propmat*,
but adds further contributions.
)--",
      .author = {"Patrick Eriksson"},
      .out    = {"spectral_propmat", "spectral_propmat_jac"},
      .in     = {"spectral_propmat",
                 "spectral_propmat_jac",
                 "freq_grid",
                 "select_species",
                 "jac_targets",
                 "atm_point",
                 "ray_point"},
  };

  wsm_data["spectral_propmatAddPredefined"] = {
      .desc =
          R"--(Adds all of the predefined models in *abs_species* to the spectral_propmat

Only supports temperature and wind speed derivatives

.. |br| raw:: html

   <br />

Available models

.. list-table::
  :width: 0.5
  :widths: auto
  :align: left
  :header-rows: 1

  * - Model name
    - Description and limitations
    - Reference(s)

  * - ``H2O-ForeignContCKDMT320``
    - Foreign continua. |br|
      Use water cutoff of 25 cm-1 and ``H2O-SelfContCKDMT320``.
    - :cite:p:`MTCKD`

  * - ``H2O-ForeignContCKDMT350``
    - MT CKD 3.5 foreign continua. |br|
      Use water cutoff of 25 cm-1 and ``H2O-SelfContCKDMT350``.
    - :cite:p:`MTCKD`

  * - ``H2O-ForeignContCKDMT400``
    - MT CKD 4 foreign continua. |br|
      Use water cutoff of 25 cm-1 and ``H2O-SelfContCKDMT350``. |br|
      Requires an external data source.
    - :cite:t:`Mlawer2012` and :cite:p:`MTCKD`

  * - ``H2O-ForeignContStandardType``
    - Water microwave continua for foreign species.
    - :cite:t:`pwr:98`

  * - ``H2O-MPM89``
    - Microwave water absorption model.
    - :cite:t:`liebe:89`

  * - ``H2O-PWR2021``
    - Microwave water absorption model |br|
      developed by P.W. Rosenkranz.
    - :cite:t:`RosenkranzSoftware`.

  * - ``H2O-PWR2022``
    - Microwave water absorption model |br|
      developed by P.W. Rosenkranz.
    - :cite:t:`RosenkranzSoftware`.

  * - ``H2O-PWR98``
    - Microwave water absorption model.
    - :cite:t:`pwr:98`

  * - ``H2O-SelfContCKDMT320``
    - Self continua. |br|
      Use water cutoff of 25 cm-1 and ``H2O-ForeignContCKDMT320``.
    - :cite:p:`MTCKD`

  * - ``H2O-SelfContCKDMT350``
    - MT CKD 3.5 self continua. |br|
      Use water cutoff of 25 cm-1 and ``H2O-ForeignContCKDMT350``.
    - :cite:p:`MTCKD`

  * - ``H2O-SelfContCKDMT400``
    - MT CKD 4 self continua. |br|
      Use water cutoff of 25 cm-1 and ``H2O-SelfContCKDMT350``. |br|
      Requires an external data source.
    - :cite:t:`Mlawer2012` and :cite:p:`MTCKD`

  * - ``H2O-SelfContStandardType``
    - Water microwave continua for self.
    - :cite:t:`pwr:98`

  * - ``CO2-CKDMT252``
    - MT CKD absorption for CO2 version 2.52.
    - :cite:p:`MTCKD`

  * - ``O2-CIAfunCKDMT100``
    - CIA for oxygen from MT CKD.
    - :cite:p:`Thibault1997` and :cite:p:`MTCKD`

  * - ``O2-MPM2020``
    - 60 GHz and 118 GHz lines only. |br|
      Do not include the v 0 0 oxygen band manually.
    - :cite:t:`Makarov2020`

  * - ``O2-MPM89``
    - Oxygen microwave absorption model.
    - :cite:t:`liebeetal:93`

  * - ``O2-PWR2021``
    - Oxygen microwave absorption model |br|
      developed by P.W. Rosenkranz.
    - :cite:t:`pwr:98`

  * - ``O2-PWR2022``
    - Oxygen microwave absorption model |br|
      developed by P.W. Rosenkranz.
    - :cite:t:`pwr:98`

  * - ``O2-PWR98``
    - Oxygen microwave absorption model.
    - :cite:t:`pwr:93` and :cite:t:`liebeetal:92` and |br|
      M.J. Schwartz, Ph.D. thesis, M.I.T. (1997) and |br|
      :cite:t:`Rothman1998`.

  * - ``O2-SelfContStandardType``
    - Microwave continua term.
    - :cite:t:`pwr:93` and :cite:t:`liebeetal:93`

  * - ``O2-TRE05``
    - Oxygen microwave absorption model.
    - :cite:t:`liebeetal:93` and :cite:t:`tretyakov05:_60-ghz_jms`

  * - ``O2-v0v0CKDMT100``
    - CKD_MT 1.00 implementation of oxygen |br|
      collision-induced fundamental model.
    - :cite:p:`MTCKD` and :cite:p:`Mate1999` and :cite:t:`Mlawer1998`

  * - ``O2-v1v0CKDMT100``
    - MT CKD.
    - :cite:p:`MTCKD` and :cite:t:`Mlawer1998`

  * - ``O2-visCKDMT252``
    - MT CKD.
    - :cite:p:`MTCKD` and :cite:t:`Greenblatt1990`

  * - ``N2-CIAfunCKDMT252``
    - MT CKD.
    - :cite:p:`MTCKD` and :cite:t:`Lafferty1996`

  * - ``N2-CIArotCKDMT252``
    - MT CKD.
    - :cite:p:`MTCKD` and :cite:t:`borysow:86`

  * - ``N2-SelfContMPM93``
    - Microwave nitrogen absorption continua |br|
      from MPM93 model.
    - :cite:t:`liebeetal:93`

  * - ``N2-SelfContPWR2021``
    - Microwave nitrogen absorption continua |br|
      developed by P.W. Rosenkranz. |br|
      This includes O2-N2 and O2-O2 CIA (only applicable to Earth).
    - :cite:t:`pwr:98`

  * - ``N2-SelfContStandardType``
    - Microwave nitrogen absorption continua.
    - :cite:t:`pwr:93`

  * - ``liquidcloud-ELL07``
    - Water droplet absorption.
    - :cite:t:`Ellison2007`
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_propmat", "spectral_propmat_jac"},
      .in     = {"spectral_propmat",
                 "spectral_propmat_jac",
                 "abs_predef_data",
                 "select_species",
                 "jac_targets",
                 "freq_grid",
                 "atm_point"},
  };

  wsm_data["spectral_propmatAddXsecFit"] = {
      .desc =
          R"--(Calculate absorption cross sections per tag group for HITRAN xsec species.

This broadens the cross section data from *abs_xfit_data* and
interpolates it onto the current *freq_grid*.
)--",
      .author    = {"Oliver Lemke"},
      .out       = {"spectral_propmat", "spectral_propmat_jac"},
      .in        = {"spectral_propmat",
                    "spectral_propmat_jac",
                    "select_species",
                    "jac_targets",
                    "freq_grid",
                    "atm_point",
                    "abs_xfit_data"},
      .gin       = {"force_p", "force_t"},
      .gin_type  = {"Numeric", "Numeric"},
      .gin_value = {Numeric{-1}, Numeric{-1}},
      .gin_desc  = {R"--(Positive value forces constant pressure [Pa].)--",
                    R"--(Positive value forces constant temperature [K].)--"},
  };

  wsm_data["spectral_propmat_scatSpectralInit"] = {
      .desc =
          R"--(Initialize *spectral_propmat_scat* and co to zeroes.

This method must be used inside *spectral_propmat_scat_spectral_agenda* and then be called first.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_propmat_scat",
                 "spectral_absvec_scat",
                 "spectral_phamat_spectral"},
      .in     = {"freq_grid", "legendre_degree"},
  };

  wsm_data["legendre_degreeFromDisortSettings"] = {
      .desc =
          R"--(Sets *legendre_degree* to *disort_settings* ``legendre_polynomial_dimension``

Method is purely for convenience and composition.
)--",
      .author = {"Richard Larsson"},
      .out    = {"legendre_degree"},
      .in     = {"disort_settings"},
  };

  wsm_data["spectral_propmat_scatAddSpectralScatteringSpeciesTRO"] = {
      .desc =
          R"--(Adds *scattering_species* results for totally random oriented spectral calculations to
*spectral_propmat_scat* and co.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_propmat_scat",
                 "spectral_absvec_scat",
                 "spectral_phamat_spectral"},
      .in     = {"spectral_propmat_scat",
                 "spectral_absvec_scat",
                 "spectral_phamat_spectral",
                 "freq_grid",
                 "atm_point",
                 "scattering_species"},
  };

  wsm_data["spectral_propmat_scat_pathFromSpectralAgenda"] = {
      .desc =
          R"--(Compute *spectral_propmat_scat_path* and co for a path.
)--",
      .author         = {"Richard Larsson"},
      .out            = {"spectral_propmat_scat_path",
                         "spectral_absvec_scat_path",
                         "spectral_phamat_spectral_path"},
      .in             = {"freq_grid_path",
                         "atm_path",
                         "legendre_degree",
                         "spectral_propmat_scat_spectral_agenda"},
      .pass_workspace = true,
  };

  wsm_data["spectral_propmat_scatInit"] = {
      .desc =
          R"--(Initialize *spectral_propmat_scat* to zeroes.

This method must be used inside *spectral_propmat_scat_agenda* and then be called first.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_propmat_scat"},
      .in     = {"freq_grid"},
  };

  wsm_data["spectral_propmat_scatAirSimple"] = {
      .desc =
          R"--(Add simple air to *spectral_propmat_scat*.
)--",
      .author = {"Jon Petersen", "Richard Larsson"},
      .out    = {"spectral_propmat_scat"},
      .in     = {"spectral_propmat_scat", "freq_grid", "atm_point"},
  };

  wsm_data["spectral_propmatInit"] = {
      .desc =
          R"--(Initialize *spectral_propmat*, *spectral_srcvec_nlte*, and their derivatives to zeroes.

This method must be used inside *spectral_propmat_agenda* and then be called first.
)--",
      .author = {"Oliver Lemke", "Richard Larsson"},
      .out    = {"spectral_propmat",
                 "spectral_srcvec_nlte",
                 "spectral_propmat_jac",
                 "spectral_srcvec_nlte_jac"},
      .in     = {"jac_targets", "freq_grid"},
  };

  wsm_data["surf_fieldPlanet"] = {
      .desc =
          R"--(Initialize the surface field with the ellipsoid of a planet.

See *PlanetOrMoonType* for valid ``option``.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"surf_field"},
      .gin       = {"option", "surf_elevation"},
      .gin_type  = {"String", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0.0}},
      .gin_desc  = {R"--(Choice of planet or moon)--",
                    "Surface elevation over the full field"},
  };

  wsm_data["surf_fieldEarth"] = {
      .desc      = R"--(Earth reference ellipsoids.

The reference ellipsoid is set to model the Earth.

See *EarthEllipsoid* for valid ``model``
)--",
      .author    = {"Patrick Eriksson"},
      .out       = {"surf_field"},
      .gin       = {"model", "surf_elevation"},
      .gin_type  = {"String", "Numeric"},
      .gin_value = {String("Sphere"), Numeric{0.0}},
      .gin_desc  = {R"--(Model ellipsoid to use. Options listed above.)--",
                    "Surface elevation over the full field"},
  };

  wsm_data["surf_fieldEuropa"] = {
      .desc      = R"--(Europa reference ellipsoids.

The reference ellipsoid is set to model the Europa.

See *EuropaEllipsoid* for valid ``model``.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"surf_field"},
      .gin       = {"model", "surf_elevation"},
      .gin_type  = {"String", "Numeric"},
      .gin_value = {String("Sphere"), Numeric{0.0}},
      .gin_desc  = {R"--(Model ellipsoid to use. Options listed above.)--",
                    "Surface elevation over the full field"},
  };

  wsm_data["surf_fieldGanymede"] = {
      .desc      = R"--(Ganymede reference ellipsoids.

See *GanymedeEllipsoid* for valid ``model``.
)--",
      .author    = {"Takayoshi Yamada"},
      .out       = {"surf_field"},
      .gin       = {"model", "surf_elevation"},
      .gin_type  = {"String", "Numeric"},
      .gin_value = {String("Sphere"), Numeric{0.0}},
      .gin_desc  = {R"--(Model ellipsoid to use. Options listed above.)--",
                    "Surface elevation over the full field"},
  };

  wsm_data["surf_fieldInit"] = {
      .desc      = R"--(Manual setting of the reference ellipsoid.

The two values of the reference ellipsoid are set manually. The two
arguments correspond directly to first and second element of
reference ellipsoid.
)--",
      .author    = {"Patrick Eriksson"},
      .out       = {"surf_field"},
      .gin       = {"a", "b", "surf_elevation"},
      .gin_type  = {"Numeric", "Numeric", "Numeric"},
      .gin_value = {std::nullopt, std::nullopt, Numeric{0.0}},
      .gin_desc  = {R"--(Average or equatorial radius.)--",
                    R"--(Average or polar radius.)--",
                    "Surface elevation over the full field"},
  };

  wsm_data["surf_fieldIo"] = {
      .desc      = R"--(Io reference ellipsoids.

The reference ellipsoid is set to model the Io.

See *IoEllipsoid* for valid ``model``.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"surf_field"},
      .gin       = {"model", "surf_elevation"},
      .gin_type  = {"String", "Numeric"},
      .gin_value = {String("Sphere"), Numeric{0.0}},
      .gin_desc  = {R"--(Model ellipsoid to use. Options listed above.)--",
                    "Surface elevation over the full field"},
  };

  wsm_data["surf_fieldJupiter"] = {
      .desc      = R"--(Jupiter reference ellipsoids.

The reference ellipsoid is set to model the Jupiter.

See *JupiterEllipsoid* for valid ``model``.
)--",
      .author    = {"Patrick Eriksson"},
      .out       = {"surf_field"},
      .gin       = {"model", "surf_elevation"},
      .gin_type  = {"String", "Numeric"},
      .gin_value = {String("Sphere"), Numeric{0.0}},
      .gin_desc  = {R"--(Model ellipsoid to use. Options listed above.)--",
                    "Surface elevation over the full field"},
  };

  wsm_data["surf_fieldMars"] = {
      .desc      = R"--(Mars reference ellipsoids.

The reference ellipsoid is set to model the Mars.

See *MarsEllipsoid* for valid ``model``.
)--",
      .author    = {"Patrick Eriksson"},
      .out       = {"surf_field"},
      .gin       = {"model", "surf_elevation"},
      .gin_type  = {"String", "Numeric"},
      .gin_value = {String("Sphere"), Numeric{0.0}},
      .gin_desc  = {R"--(Model ellipsoid to use. Options listed above.)--",
                    "Surface elevation over the full field"},
  };

  wsm_data["surf_fieldMoon"] = {
      .desc      = R"--(Moon reference ellipsoids.

The reference ellipsoid is set to model the Moon.

See *MoonEllipsoid* for valid ``model``.
)--",
      .author    = {"Patrick Eriksson"},
      .out       = {"surf_field"},
      .gin       = {"model", "surf_elevation"},
      .gin_type  = {"String", "Numeric"},
      .gin_value = {String("Sphere"), Numeric{0.0}},
      .gin_desc  = {R"--(Model ellipsoid to use. Options listed above.)--",
                    "Surface elevation over the full field"},
  };

  wsm_data["surf_fieldVenus"] = {
      .desc      = R"--(Venus reference ellipsoids.

The reference ellipsoid is set to model the Venus.

See *VenusEllipsoid* for valid ``model``.
)--",
      .author    = {"Patrick Eriksson"},
      .out       = {"surf_field"},
      .gin       = {"model", "surf_elevation"},
      .gin_type  = {"String", "Numeric"},
      .gin_value = {String("Sphere"), Numeric{0.0}},
      .gin_desc  = {R"--(Model ellipsoid to use. Options listed above.)--",
                    "Surface elevation over the full field"},
  };

  wsm_data["water_equivalent_pressure_operatorMK05"] = {
      .desc =
          R"--(Calculate equivalent water pressure according to Murphy and Koop, 2005

Default is setting the saturation pressure to the one with respect
to water at temperatures >= 0C, and to the one with respect to ice
for <0C. The GIN ``only_liquid`` allows you to apply the liquid value
at all temperatures.

The saturation pressure with respect to liquid and ice water is
calculated according to Eq. 10 and 7, respectively, of:
Murphy, D. M., & Koop, T. (2005). Review of the vapour pressures of
ice and supercooled water for atmospheric applications. Quarterly
Journal of the Royal Meteorological Society, 131(608), 1539-1565.
)--",
      .author    = {"Patrick Eriksson", "Richard Larsson"},
      .out       = {"water_equivalent_pressure_operator"},
      .gin       = {"only_liquid"},
      .gin_type  = {"Index"},
      .gin_value = {Index{0}},
      .gin_desc =
          {"Set to 1 to use liquid saturation pressure at all temperatures"},
  };

  wsm_data["atm_fieldHydrostaticPressure"] = {
      .desc      = R"-x-(Add the hydrostatic pressure to the atmospheric field
    
The field must already be able to compute temperature as a function of
altitude, latitude, and longitude.

If a positive ``fixed_specific_gas_constant`` is not provided,
the field must also consist of correct volume mixing ratios so
that the mass of an average molecule can be computed.

The first altitude in *alt_grid*
is used as the altitude of the ``p0`` grid.  The extrapolation
outside of this range simply uses the formalism of  the select ``hydrostatic_option``.

.. note::
  The gradient changes only at the grid points of the *alt_grid*.
  Please make it dense enough to avoid missing features.  A recommendation
  is to extract the *alt_grid* directly from the temperature field.

  Also be aware that missing VMRs for important species, e.g.,
  :math:`\textrm{N}_2` or :math:`\textrm{O}_2`, will lead to incorrect results.
  The mean molecular mass need these VMRs, so if the VMRs are missing,
  the pressure will not be correct either.
)-x-",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field", "gravity_operator", "alt_grid"},
      .gin       = {"p0",
                    "fixed_specific_gas_constant",
                    "fixed_atmospheric_temperature",
                    "hydrostatic_option"},
      .gin_type  = {"GeodeticField2,Numeric", "Numeric", "Numeric", "String"},
      .gin_value = {std::nullopt,
                    Numeric{-1},
                    Numeric{-1},
                    String{"HydrostaticEquation"}},
      .gin_desc =
          {"Lowest altitude pressure field.  :math:`P_0` above.",
           "Specific gas constant if larger than 0",
           "Constant atmospheric temprature if larger than 0",
           "Computational option for levels.  See *HydrostaticPressureOption* for valid options."},
  };

  wsm_data["gravity_operatorCentralMass"] = {
      .desc =
          R"-x-(Sets a gravity operator from the gravitational constant and the mass of the planet

Gets the ellispoid from *surf_field*
)-x-",
      .author    = {"Richard Larsson"},
      .out       = {"gravity_operator"},
      .in        = {"surf_field"},
      .gin       = {"mass"},
      .gin_type  = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {"Gravitation constant so that the gravity at radius :math:`r` is :math:`GM / r^2`"},
  };

  wsm_data["spectral_flux_profileFromPathField"] = {
      .desc           = R"--(Computes the spectral flux from a field of paths.
)--",
      .author         = {"Richard Larsson"},
      .out            = {"spectral_flux_profile"},
      .in             = {"ray_path_field",
                         "atm_field",
                         "spectral_propmat_agenda",
                         "spectral_rad_space_agenda",
                         "spectral_rad_surface_agenda",
                         "surf_field",
                         "subsurf_field",
                         "freq_grid",
                         "alt_grid"},
      .pass_workspace = true,
  };

  wsm_data["spectral_flux_profileFromSpectralRadianceField"] = {
      .desc =
          R"--(Computes the spectral flux.  The input field must be a profile.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"spectral_flux_profile"},
      .in        = {"spectral_rad_field"},
      .gin       = {"pol"},
      .gin_type  = {"Stokvec"},
      .gin_value = {Stokvec{1.0, 0.0, 0.0, 0.0}},
      .gin_desc  = {"Polarization vector for the spectral flux profile"},
  };

  wsm_data["flux_profileIntegrate"] = {
      .desc      = R"--(Computes the spectral flux
)--",
      .author    = {"Richard Larsson"},
      .gout      = {"flux_profile"},
      .gout_type = {"Vector"},
      .gout_desc = {"The spectral flux profile"},
      .in        = {"spectral_flux_profile", "freq_grid"},
  };

  wsm_data["nlte_line_flux_profileIntegrate"] = {
      .desc =
          R"--(Integrate the spectral flux profile to get the line non-LTE flux
)--",
      .author = {"Richard Larsson"},
      .out    = {"nlte_line_flux_profile"},
      .in = {"spectral_flux_profile", "abs_bands", "atm_profile", "freq_grid"},
  };

  wsm_data["single_radClearskyEmissionPropagation"] = {
      .desc =
          R"--(Computes the spectral radiance for a single frequency using clear-sky emission propagation.

The path is built based on current optical properties and the radiative transfer
equation is solved along the path.  This means that the path is not precomputed
but built on-the-fly, allowing per-frequency refraction.
)--",
      .author = {"Richard Larsson"},
      .out    = {"single_rad", "single_rad_jac", "ray_path"},
      .in     = {"atm_field",
                 "freq",
                 "jac_targets",
                 "single_rad_space_agenda",
                 "single_rad_surface_agenda",
                 "single_propmat_agenda",
                 "ray_point_back_propagation_agenda",
                 "subsurf_field",
                 "surf_field",
                 "obs_pos",
                 "obs_los",
                 "max_stepsize"},
      .gin = {"polarization", "max_tau", "cutoff_tau", "hse_derivative", "N"},
      .gin_type  = {"Propmat", "Numeric", "Numeric", "Index", "Index"},
      .gin_value = {Propmat{0, 0, 0, 0, 0, 0, 0},
                    Numeric{0.01},
                    Numeric{14.0},
                    Index{0},
                    Index{1}},
      .gin_desc =
          {"Delta of the dispersion in polarizized form.  "
           "The dot-product of this and the propagation matrix is added to the internal *single_dispersion* variable.",
           "The maximum optical thickness per step, min of local *Propmat* A divided by ``max_tau`` "
           "and *max_stepsize* is passed to *ray_point_back_propagation_agendaExecute*.  "
           "Note that this is an approximation that will fail for highly non-linear absorption profiles.  "
           "As implemented, it takes too long steps if going from low to high absorption, and "
           "too short steps when going from high to low absorption.  See it as an approximation.",
           "Cutoff optical thickness for terminating the integration, "
           "computed as total *Propmat* A times distance.  "
           "If exceeded, the atmosphere is considered opaque and the temperature at that coordinate "
           "is used for the background radiation.  "
           "If not exceeded, the actual background is considered.  "
           "Note that errors will be large if exp(-``cutoff_tau``) is not small.",
           "Flag to compute the hypsometric distance derivatives",
           "Number of points to reserve in the ray path"},
      .pass_workspace = true,
  };

  wsm_data["spectral_radClearskyEmissionFrequencyDependentPropagation"] = {
      .desc =
          R"--(Wraps *single_radClearskyEmissionPropagation* for a vector of frequencies.
)--",
      .author         = {"Richard Larsson"},
      .out            = {"spectral_rad", "spectral_rad_jac"},
      .gout           = {"ray_paths"},
      .gout_type      = {"ArrayOfArrayOfPropagationPathPoint"},
      .gout_desc      = {"The ray paths for each frequency"},
      .in             = {"atm_field",
                         "freq_grid",
                         "jac_targets",
                         "single_rad_space_agenda",
                         "single_rad_surface_agenda",
                         "single_propmat_agenda",
                         "ray_point_back_propagation_agenda",
                         "subsurf_field",
                         "surf_field",
                         "obs_pos",
                         "obs_los",
                         "max_stepsize"},
      .pass_workspace = true,
  };
  wsm_data["spectral_radClearskyEmissionFrequencyDependentPropagation"].gin =
      wsm_data["single_radClearskyEmissionPropagation"].gin;
  wsm_data["spectral_radClearskyEmissionFrequencyDependentPropagation"]
      .gin_type = wsm_data["single_radClearskyEmissionPropagation"].gin_type;
  wsm_data["spectral_radClearskyEmissionFrequencyDependentPropagation"]
      .gin_value = wsm_data["single_radClearskyEmissionPropagation"].gin_value;
  wsm_data["spectral_radClearskyEmissionFrequencyDependentPropagation"]
      .gin_desc = wsm_data["single_radClearskyEmissionPropagation"].gin_desc;

  wsm_data["single_radFromVector"] = {
      .desc =
          R"--(Composition method to extract a single spectral radiance from a vector.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"single_rad", "single_rad_jac"},
      .in        = {"spectral_rad", "spectral_rad_jac"},
      .gin       = {"index"},
      .gin_type  = {"Index"},
      .gin_value = {Index{0}},
      .gin_desc  = {"Index into the first dimension of the spectral radiance"},
  };

  wsm_data["spectral_rad_bkgAgendasAtEndOfPath"] = {
      .desc           = R"--(Computes the background radiation.

This method checks the back of the ray path and calls
*spectral_rad_space_agenda* and *spectral_rad_surface_agenda* as needed.

If the back of the path is still in the atmosphere, an
error is raised.  As is it if the background position is unknown.

.. tip::
    To access subsurface emission, the *spectral_rad_surface_agenda* must be able to handle
    subsurface emission.  It will likely need the surface emission agenda as well,
    but that would be circular.  Therefore, please consider using *spectral_rad_closed_surface_agenda*
    to compute the surface emission without invoking recursion.
)--",
      .author         = {"Richard Larsson"},
      .out            = {"spectral_rad_bkg", "spectral_rad_bkg_jac"},
      .in             = {"freq_grid",
                         "jac_targets",
                         "ray_point",
                         "surf_field",
                         "subsurf_field",
                         "spectral_rad_space_agenda",
                         "spectral_rad_surface_agenda"},
      .pass_workspace = true,
  };

  wsm_data["spectral_radDefaultTransmission"] = {
      .desc =
          R"--(Sets default *spectral_rad* and *spectral_rad_jac* for transmission.

The Jacobian variable is all 0s, the background is [1 0 0 0] everywhere
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad", "spectral_rad_jac"},
      .in     = {"freq_grid", "jac_targets"},
  };

  wsm_data["spectral_radSetToBackground"] = {
      .desc =
          R"--(Set the spectral radiance to the background values.

.. note::
    This method is compositional and useful inside other agendas or meta-methods.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad", "spectral_rad_jac"},
      .in     = {"spectral_rad_bkg", "spectral_rad_bkg_jac"},
  };

  wsm_data["spectral_radSinglePathEmissionFrequencyLoop"] = {
      .desc =
          R"--(Computes the spectral radiance by looping over frequencies for a single path.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"spectral_rad", "spectral_rad_jac"},
      .in        = {"spectral_rad",
                    "spectral_rad_jac",
                    "jac_targets",
                    "ray_path",
                    "freq_grid_path",
                    "atm_path",
                    "spectral_propmat_path",
                    "spectral_srcvec_nlte_path",
                    "spectral_propmat_jac_path",
                    "spectral_srcvec_nlte_jac_path",
                    "surf_field",
                    "atm_field"},
      .gin       = {"hse_derivative"},
      .gin_type  = {"Index"},
      .gin_value = {Index{0}},
      .gin_desc  = {"Flag to compute the hypsometric distance derivatives"},
  };

  wsm_data["spectral_radUniformCosmicBackground"] = {
      .desc =
          R"--(Background spectral radiance is from a uniform cosmic background temperature.

.. math::
    \vec{I} = \left[ \begin{array}{c} \frac{2h\vec{\nu}^3}{c^2} \frac{1}{e^{\frac{h\vec{\nu}}{kT_c}} - 1} \\ 0 \\ 0 \\ 0 \end{array}\right],

where :math:`T_c` is the cosmic microwave background temperature.  :math:`\nu` is the frequency grid.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad"},
      .in     = {"freq_grid"},
  };

  wsm_data["spectral_radSurfaceBlackbody"] = {
      .desc =
          R"--(Set surface spectral radiance from Planck function of the surface temperature.

.. math::
    \vec{I} = \left[ \begin{array}{c} \frac{2h\vec{\nu}^3}{c^2} \frac{1}{e^{\frac{h\vec{\nu}}{kT_s}} - 1} \\ 0 \\ 0 \\ 0 \end{array}\right],

where :math:`T_s` is the surface temperature extracted at the input *ray_point*.  :math:`\nu` is the frequency grid.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad", "spectral_rad_jac"},
      .in     = {"freq_grid", "surf_field", "jac_targets", "ray_point"},
  };

  wsm_data["spectral_radSurfaceReflectance"] = {
      .desc =
          R"--(Set surface spectral radiance to use sub-surface emission and Fresnel reflectance.

The input path point must be close to the surface.

The *spectral_rad_closed_surface_agenda* should produce the surface emission,
though pure surface emission is fine.

The surface field must contain the surface refractive index.
The refractive index lives under the *SurfacePropertyTag* key "scalar refractive index".
)--",
      .author         = {"Richard Larsson"},
      .out            = {"spectral_rad", "spectral_rad_jac"},
      .in             = {"freq_grid",
                         "atm_field",
                         "surf_field",
                         "subsurf_field",
                         "jac_targets",
                         "ray_point",
                         "spectral_rad_observer_agenda",
                         "spectral_rad_closed_surface_agenda",
                         "spectral_surf_refl_agenda"},
      .pass_workspace = true,
  };

  wsm_data["spectral_rad_jacAddSensorJacobianPerturbations"] = {
      .desc   = R"--(Adds sensor properties to the *spectral_rad_jac*.

This is done via perturbation based on the input delta values to the sensor
Jacobian targets and a callback to *spectral_rad_observer_agenda* with
a modified *jac_targets*, making it safe to use this method inside
*spectral_rad_observer_agenda*.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad_jac"},
      .in =
          {
              "spectral_rad_jac",
              "spectral_rad",
              "measurement_sensor",
              "freq_grid",
              "jac_targets",
              "obs_pos",
              "obs_los",
              "atm_field",
              "surf_field",
              "subsurf_field",
              "spectral_rad_observer_agenda",
          },
      .pass_workspace = true,
  };

  wsm_data["spectral_rad_jacEmpty"] = {
      .desc   = R"--(Set the radiation derivative to empty.

Size : (*jac_targets*, *freq_grid*)
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad_jac"},
      .in     = {"freq_grid", "jac_targets"},
  };

  wsm_data["spectral_rad_jacFromBackground"] = {
      .desc =
          R"--(Sets *spectral_rad_jac* from the background values.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad_jac"},
      .in     = {"spectral_rad_bkg_jac", "spectral_tramat_bkg"},
  };

  wsm_data["spectral_rad_jacAddPathPropagation"] = {
      .desc =
          R"--(Adds the propagation variables to *spectral_rad_jac*.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad_jac"},
      .in     = {"spectral_rad_jac",
                 "spectral_rad_jac_path",
                 "jac_targets",
                 "atm_field",
                 "ray_path"},
  };

  wsm_data["spectral_rad_transform_operatorSet"] = {
      .desc =
          "Creates a *SpectralRadianceTransformOperator* from a *SpectralRadianceUnitType*.\n",
      .author    = {"Richard Larsson"},
      .out       = {"spectral_rad_transform_operator"},
      .gin       = {"option"},
      .gin_type  = {"SpectralRadianceUnitType"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {"The type of the spectral radiance transform operator to create. See *SpectralRadianceUnitType* for valid values and what they do."},
  };

  wsm_data["spectral_radApplyUnit"] = {
      .desc   = R"(Applies a unit to *spectral_rad*, returning a new field

See *SpectralRadianceUnitType* for valid use cases and limitations.

This effectively wraps the local creation of a *SpectralRadianceTransformOperator* call.

.. warning::
  This is a destructive method.  Any use of it means that it is undefined behavior
  to use *spectral_rad* or *spectral_rad_jac* in future methods.
)",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad", "spectral_rad_jac"},
      .in =
          {
              "spectral_rad",
              "spectral_rad_jac",
              "freq_grid",
              "ray_point",
              "spectral_rad_transform_operator",
          },
  };

  wsm_data["spectral_radApplyForwardUnit"] = {
      .desc =
          R"(Helper to call *spectral_radApplyUnit* when you do not have *spectral_rad_jac*.
)",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad"},
      .in =
          {
              "spectral_rad",
              "freq_grid",
              "ray_point",
              "spectral_rad_transform_operator",
          },
  };

  wsm_data["spectral_surf_reflFlatRealFresnel"] = {
      .desc =
          R"--(Set the surface reflectance to the flat real Fresnel reflectance

.. math::
    \begin{array}{lcr}
        \theta_2 &=& \arcsin\left(\frac{\Re{n_1}}{\Re{n_2}}\sin{\theta_1}\right)\\[5pt]
        R_v &=& \frac{n_2\cos\left(\theta_1\right) -
                      n_1\cos\left(\theta_2\right)}
                     {n_2\cos\left(\theta_1\right) +
                      n_1\cos\left(\theta_2\right)}\\[5pt]
        R_h &=& \frac{n_1\cos\left(\theta_1\right) -
                      n_2\cos\left(\theta_2\right)}
                     {n_1\cos\left(\theta_1\right) +
                      n_2\cos\left(\theta_2\right)},
    \end{array}

where :math:`\theta_1` is the angle of incidence, :math:`\theta_2` is the angle of refraction, and
:math:`n_1` and :math:`n_2` are the refractive indices of the two media.

We get :math:`n_1` and :math:`\theta_1` from the *ray_point* and extracts
:math:`n_2` from the *surf_field* parameter ``"scalar refractive index"``.

The reflectance matrix is

.. math::
    \mathbf{R} = \frac{1}{2}\left[
    \begin{array}{cccc}
            R_v\overline{R_v} + R_h\overline{R_h} & R_v\overline{R_v} - R_h\overline{R_h} & 0 & 0 \\
            R_v\overline{R_v} - R_h\overline{R_h} & R_v\overline{R_v} + R_h\overline{R_h} & 0 & 0 \\
            0 & 0 & \Re\left(R_h\overline{R_v} + R_v\overline{R_h}\right) & \Im\left(R_h\overline{R_v} - R_v\overline{R_h}\right) \\
            0 & 0 & \Im\left(R_v\overline{R_h} - R_h\overline{R_v}\right) & \Re\left(R_h\overline{R_v} + R_v\overline{R_h}\right) \\
    \end{array}\right]
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_surf_refl", "spectral_surf_refl_jac"},
      .in     = {"freq_grid", "surf_field", "ray_point", "jac_targets"},
  };

  wsm_data["spectral_surf_reflFlatScalar"] = {
      .desc =
          R"--(Set the surface reflectance to the flat real Fresnel reflectance

We get :math:`r` from the *surf_field* parameter ``"flat scalar reflectance"``.

The reflectance matrix is

.. math::
    \mathbf{R} = \left[
    \begin{array}{cccc}
            r&0&0&0\\
            0&r&0&0\\
            0&0&r&0\\
            0&0&0&r\\
    \end{array}\right]
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_surf_refl", "spectral_surf_refl_jac"},
      .in     = {"freq_grid", "surf_field", "ray_point", "jac_targets"},
  };

  wsm_data["spectral_propmat_jacWindFix"] = {
      .desc   = R"--(Fix for the wind field derivative.

The *spectral_propmat_agenda* will set the wind derivatives to 
those of the frequency derivative if this method is not used.  This
will cause the wind field to be treated as a frequency derivative,
meaning no *OEM* or other functionality that requires the Jacobian
matrix to be calculated will work.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_propmat_jac", "spectral_srcvec_nlte_jac"},
      .in     = {"spectral_propmat_jac",
                 "spectral_srcvec_nlte_jac",
                 "freq_grid",
                 "jac_targets",
                 "freq_wind_shift_jac"},
  };

  wsm_data["spectral_propmatAddLines"] = {
      .desc      = R"--(Add line-by-line absorption to the propagation matrix.

See :doc:`concept.absorption.lbl` for details.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"spectral_propmat",
                    "spectral_srcvec_nlte",
                    "spectral_propmat_jac",
                    "spectral_srcvec_nlte_jac"},
      .in        = {"spectral_propmat",
                    "spectral_srcvec_nlte",
                    "spectral_propmat_jac",
                    "spectral_srcvec_nlte_jac",
                    "freq_grid",
                    "jac_targets",
                    "select_species",
                    "abs_bands",
                    "abs_ecs_data",
                    "atm_point",
                    "ray_point"},
      .gin       = {"no_negative_absorption"},
      .gin_type  = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc =
          {"Turn off to allow individual absorbers to have negative absorption"},
  };

  wsm_data["spectral_propmatAddVoigtLTE"] = {
      .desc      = R"--(Add line-by-line absorption to the propagation matrix.

See :doc:`concept.absorption.lbl` for details.

This is only for LTE lines in Voigt.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"spectral_propmat", "spectral_propmat_jac"},
      .gout      = {"dispersion", "dispersion_jacobian"},
      .gout_type = {"Vector", "Matrix"},
      .gout_desc =
          {"Dispersion vector - only the main component (i.e., -imag(A) of the *Propmat*)",
           "Dispersion Jacobian matrix - only the main component (i.e., -imag(A) of the *Propmat*)"},
      .in        = {"spectral_propmat",
                    "spectral_propmat_jac",
                    "freq_grid",
                    "jac_targets",
                    "select_species",
                    "abs_bands",
                    "atm_point",
                    "ray_point"},
      .gin       = {"no_negative_absorption"},
      .gin_type  = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc =
          {"Turn off to allow individual absorbers to have negative absorption"},
  };

  wsm_data["single_propmatInit"] = {
      .desc   = R"--(Initialize single-point propagation matrix fields.
)--",
      .author = {"Richard Larsson"},
      .out    = {"single_propmat",
                 "single_propmat_jac",
                 "single_nlte_srcvec",
                 "single_nlte_srcvec_jac",
                 "single_dispersion",
                 "single_dispersion_jac"},
      .in     = {"jac_targets"}};

  wsm_data["single_propmatAddVoigtLTE"] = {
      .desc      = R"--(Add line-by-line absorption to the propagation matrix.

See :doc:`concept.absorption.lbl` for details.

This is only for LTE lines in Voigt.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"single_propmat",
                    "single_propmat_jac",
                    "single_dispersion",
                    "single_dispersion_jac"},
      .in        = {"single_propmat",
                    "single_propmat_jac",
                    "single_dispersion",
                    "single_dispersion_jac",
                    "freq",
                    "jac_targets",
                    "select_species",
                    "abs_bands",
                    "atm_point",
                    "ray_point"},
      .gin       = {"no_negative_absorption"},
      .gin_type  = {"Index"},
      .gin_value = {Index{1}},
      .gin_desc =
          {"Turn off to allow individual absorbers to have negative absorption"},
  };

  wsm_data["spectral_propmatAddLookup"] = {
      .desc     = R"--(Add line-by-line absorption to the propagation matrix.

See :doc:`concept.absorption.lookup` for details.
)--",
      .author   = {"Richard Larsson"},
      .out      = {"spectral_propmat", "spectral_propmat_jac"},
      .in       = {"spectral_propmat",
                   "spectral_propmat_jac",
                   "freq_grid",
                   "jac_targets",
                   "select_species",
                   "abs_lookup_data",
                   "atm_point"},
      .gin      = {"no_negative_absorption",
                   "p_interp_order",
                   "t_interp_order",
                   "water_interp_order",
                   "f_interp_order",
                   "extpolfac"},
      .gin_type = {"Index", "Index", "Index", "Index", "Index", "Numeric"},
      .gin_value =
          {Index{1}, Index{7}, Index{7}, Index{7}, Index{7}, Numeric{0.5}},
      .gin_desc =
          {"Turn off to allow individual absorbers to have negative absorption",
           "Interpolation order for pressure",
           "Interpolation order for temperature",
           "Interpolation order for water vapor",
           "Interpolation order for frequency",
           "Extrapolation factor"},
  };

  wsm_data["jac_targetsToggleRelativeHumidityAtmTarget"] = {
      .desc   = R"--(Toggles relative humidity or absolute retrievals.

If the target is in relative humidity mode, it becomes absolute.
If the target is not in relative humidity mode, it becomes relative humidity.

Overwrites all other functional toggles.
)--",
      .author = {"Richard Larsson"},
      .out    = {"jac_targets"},
      .in  = {"jac_targets", "atm_field", "water_equivalent_pressure_operator"},
      .gin = {"key", "nonnegative"},
      .gin_type =
          {"AtmKey,SpeciesEnum,SpeciesIsotope,QuantumLevelIdentifier,ScatteringSpeciesProperty",
           "Index"},
      .gin_value = {std::nullopt, Index{1}},
      .gin_desc  = {"Key to toggle",
                    "Whether or not to zero-out negative values"},
  };

  wsm_data["jac_targetsToggleRelativeAtmTarget"] = {
      .desc   = R"--(Toggles relative or absolute retrievals.

If the target is in relative mode, it becomes absolute.
If the target is not in relative mode, it becomes relative.

Overwrites all other functional toggles.
)--",
      .author = {"Richard Larsson"},
      .out    = {"jac_targets"},
      .in     = {"jac_targets", "atm_field"},
      .gin    = {"key"},
      .gin_type =
          {"AtmKey,SpeciesEnum,SpeciesIsotope,QuantumLevelIdentifier,ScatteringSpeciesProperty"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"Key to toggle"},
  };

  wsm_data["jac_targetsToggleLogarithmicAtmTarget"] = {
      .desc   = R"--(Toggles logarithmic or absolute retrievals.

If the target is in logarithmic mode, it becomes absolute.
If the target is not in logarithmic mode, it becomes logarithmic.

Overwrites all other functional toggles.
)--",
      .author = {"Richard Larsson"},
      .out    = {"jac_targets"},
      .in     = {"jac_targets", "atm_field"},
      .gin    = {"key"},
      .gin_type =
          {"AtmKey,SpeciesEnum,SpeciesIsotope,QuantumLevelIdentifier,ScatteringSpeciesProperty"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"Key to toggle"},
  };

  wsm_data["jac_targetsToggleLogRelAtmTarget"] = {
      .desc   = R"--(Toggles logarithmic/relative or absolute retrievals.

This means to take the logarithm of the relative value.

If the target is in logarithmic/relative mode, it becomes absolute.
If the target is not in logarithmic/relative mode, it becomes logarithmic/relative.

Overwrites all other functional toggles.
)--",
      .author = {"Richard Larsson"},
      .out    = {"jac_targets"},
      .in     = {"jac_targets", "atm_field"},
      .gin    = {"key"},
      .gin_type =
          {"AtmKey,SpeciesEnum,SpeciesIsotope,QuantumLevelIdentifier,ScatteringSpeciesProperty"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"Key to toggle"},
  };

  wsm_data["jac_targetsToggleRelativeSurfaceTarget"] = {
      .desc      = R"--(Toggles relative or absolute retrievals.

If the target is in relative mode, it becomes absolute.
If the target is not in relative mode, it becomes relative.

Overwrites all other functional toggles.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets", "surf_field"},
      .gin       = {"key"},
      .gin_type  = {"SurfaceKey,SurfacePropertyTag"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"Key to toggle"},
  };

  wsm_data["jac_targetsToggleLogarithmicSurfaceTarget"] = {
      .desc      = R"--(Toggles logarithmic or absolute retrievals.

If the target is in logarithmic mode, it becomes absolute.
If the target is not in logarithmic mode, it becomes logarithmic.

Overwrites all other functional toggles.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets", "surf_field"},
      .gin       = {"key"},
      .gin_type  = {"SurfaceKey,SurfacePropertyTag"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"Key to toggle"},
  };

  wsm_data["jac_targetsToggleLogRelSurfaceTarget"] = {
      .desc      = R"--(Toggles logarithmic/relative or absolute retrievals.

This means to take the logarithm of the relative value.

If the target is in logarithmic/relative mode, it becomes absolute.
If the target is not in logarithmic/relative mode, it becomes logarithmic/relative.

Overwrites all other functional toggles.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets", "surf_field"},
      .gin       = {"key"},
      .gin_type  = {"SurfaceKey,SurfacePropertyTag"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"Key to toggle"},
  };

  wsm_data["jac_targetsToggleRelativeSubsurfaceTarget"] = {
      .desc      = R"--(Toggles relative or absolute retrievals.

If the target is in relative mode, it becomes absolute.
If the target is not in relative mode, it becomes relative.

Overwrites all other functional toggles.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets", "subsurf_field"},
      .gin       = {"key"},
      .gin_type  = {"SubsurfaceKey"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"Key to toggle"},
  };

  wsm_data["jac_targetsToggleLogarithmicSubsurfaceTarget"] = {
      .desc      = R"--(Toggles logarithmic or absolute retrievals.

If the target is in logarithmic mode, it becomes absolute.
If the target is not in logarithmic mode, it becomes logarithmic.

Overwrites all other functional toggles.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets", "subsurf_field"},
      .gin       = {"key"},
      .gin_type  = {"SubsurfaceKey"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"Key to toggle"},
  };

  wsm_data["jac_targetsToggleLogRelSubsurfaceTarget"] = {
      .desc      = R"--(Toggles logarithmic/relative or absolute retrievals.

This means to take the logarithm of the relative value.

If the target is in logarithmic/relative mode, it becomes absolute.
If the target is not in logarithmic/relative mode, it becomes logarithmic/relative.

Overwrites all other functional toggles.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets", "subsurf_field"},
      .gin       = {"key"},
      .gin_type  = {"SubsurfaceKey"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"Key to toggle"},
  };

  wsm_data["jac_targetsConditionalClear"] = {
      .desc   = R"--(Clears *jac_targets* if *do_jac* evaluates false.
)--",
      .author = {"Richard Larsson"},
      .out    = {"jac_targets"},
      .in     = {"jac_targets", "do_jac"},
  };

  wsm_data["jac_targetsInit"] = {
      .desc   = R"--(Initialize or reset the *jac_targets*.
)--",
      .author = {"Richard Larsson"},
      .out    = {"jac_targets"},
  };

  wsm_data["jac_targetsOff"] = {
      .desc   = R"--(Turns off *jac_targets*
)--",
      .author = {"Richard Larsson"},
      .out    = {"jac_targets"},
  };

  wsm_data["jac_targetsFinalize"] = {
      .desc   = R"--(Finalize *jac_targets*.

The finalization computes the size of the required *model_state_vector*.
It is thus necessary if any *OEM* or other functionality that requires the
building of an actual Jacobian matrix.
)--",
      .author = {"Richard Larsson"},
      .out    = {"jac_targets"},
      .in     = {"jac_targets",
                 "atm_field",
                 "surf_field",
                 "subsurf_field",
                 "abs_bands",
                 "measurement_sensor"},
  };

  const auto jac2ret = [&wsm_data](const std::string& name) {
    auto v  = wsm_data[name];
    v.desc += std::format(R"(
This method wraps *{}* together with adding the covariance matrices,
to the *covariance_matrix_diagonal_blocks*, which are required to perform *OEM*.

The input covariance matrices must fit the size of the later computed model state
represented by the *jac_targets*.  The covariance matrix inverse 
)",
                          name);
    v.out.insert(v.out.begin() + 1, "covariance_matrix_diagonal_blocks");
    v.in.insert(v.in.begin() + 1, "covariance_matrix_diagonal_blocks");
    v.gin.insert(v.gin.end(), "matrix");
    v.gin.insert(v.gin.end(), "inverse");
    v.gin_type.insert(v.gin_type.end(), "BlockMatrix");
    v.gin_type.insert(v.gin_type.end(), "BlockMatrix");
    v.gin_value.insert(v.gin_value.end(), std::nullopt);
    v.gin_value.insert(v.gin_value.end(), BlockMatrix{});
    v.gin_desc.insert(v.gin_desc.end(), "The covariance diagonal block matrix");
    v.gin_desc.insert(v.gin_desc.end(),
                      "The inverse covariance diagonal block matrix");
    return v;
  };

  wsm_data["jac_targetsAddErrorPolyFit"] = {
      .desc =
          R"--(Set a measurement error to polynomial fit.

This is a generic error that is simply added to *measurement_vector* as if

.. math::

    y = y_0 + \epsilon(p_0,\; p_1,\; \cdots,\; p_n),

where :math:`y` represents *measurement_vector* and :math:`y_0` is the measurement vector without any errors)

Order 0 means constant: :math:`y = y_0 + a`

Order 1 means linear:   :math:`y = y_0 + a + b t`

and so on.  The derivatives that are added to the *model_state_vector* are
those with regards to a, b, etc..

.. note::

    The rule for the ``sensor_elem`` GIN is a bit complex.  Generally, methods such
    as *measurement_sensorAddSimple* will simply add a single unique frequency grid
    to all the different *SensorObsel* that they add to the *measurement_sensor*.
    The GIN ``sensor_elem`` is 0 for the first unique frequency grid, 1 for the second,
    and so on.  See *ArrayOfSensorObsel* member methods in python for help identifying
    and manipulating how many unique frequency grids are available in *measurement_sensor*.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets", "measurement_sensor"},
      .gin       = {"t", "sensor_elem", "polyorder"},
      .gin_type  = {"Vector", "Index", "Index"},
      .gin_value = {std::nullopt, std::nullopt, Index{0}},
      .gin_desc =
          {"The grid of :math:`y`.  As :math:`t` above.",
           "The sensor element whose frequency grid to use",
           "The order of the polynomial fit.  Maximum :math:`n` above."},
  };
  wsm_data["RetrievalAddErrorPolyFit"] = jac2ret("jac_targetsAddErrorPolyFit");

  wsm_data["jac_targetsAddSensorFrequencyPolyOffset"] = {
      .desc =
          R"--(Set sensor frequency derivative to use polynomial fitting offset

Order 0 means constant: :math:`f := f_0 + a`

Order 1 means linear:   :math:`f := f_0 + a + b f_0`

and so on.  The derivatives that are added to the *model_state_vector* are
those with regards to a, b, etc..

.. note::

    The rule for the ``sensor_elem`` GIN is a bit complex.  Generally, methods such
    as *measurement_sensorAddSimple* will simply add a single unique frequency grid
    to all the different *SensorObsel* that they add to the *measurement_sensor*.
    The GIN ``sensor_elem`` is 0 for the first unique frequency grid, 1 for the second,
    and so on.  See *ArrayOfSensorObsel* member methods in python for help identifying
    and manipulating how many unique frequency grids are available in *measurement_sensor*.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets", "measurement_sensor"},
      .gin       = {"d", "sensor_elem", "polyorder"},
      .gin_type  = {"Numeric", "Index", "Index"},
      .gin_value = {Numeric{0.1}, std::nullopt, Index{0}},
      .gin_desc =
          {"The perturbation used in methods that cannot compute derivatives analytically",
           "The sensor element whose frequency grid to use",
           "The order of the polynomial fit"},
  };
  wsm_data["RetrievalAddSensorFrequencyPolyOffset"] =
      jac2ret("jac_targetsAddSensorFrequencyPolyOffset");

  wsm_data["jac_targetsAddTemperature"] = {
      .desc      = R"--(Set temperature derivative.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets"},
      .gin       = {"d"},
      .gin_type  = {"Numeric"},
      .gin_value = {Numeric{0.1}},
      .gin_desc =
          {"The perturbation used in methods that cannot compute derivatives analytically"},
  };
  wsm_data["RetrievalAddTemperature"] = jac2ret("jac_targetsAddTemperature");

  wsm_data["jac_targetsAddPressure"] = {
      .desc      = R"--(Set pressure derivative.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets"},
      .gin       = {"d"},
      .gin_type  = {"Numeric"},
      .gin_value = {Numeric{0.1}},
      .gin_desc =
          {"The perturbation used in methods that cannot compute derivatives analytically"},
  };
  wsm_data["RetrievalAddPressure"] = jac2ret("jac_targetsAddPressure");

  wsm_data["jac_targetsAddMagneticField"] = {
      .desc      = R"--(Set magnetic field derivative.

See *FieldComponent* for valid ``component``.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets"},
      .gin       = {"component", "d"},
      .gin_type  = {"String", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0.1}},
      .gin_desc =
          {"The component to use [u, v, w]",
           "The perturbation used in methods that cannot compute derivatives analytically"},
  };
  wsm_data["RetrievalAddMagneticField"] =
      jac2ret("jac_targetsAddMagneticField");

  wsm_data["jac_targetsAddOverlappingMagneticField"] = {
      .desc   = R"--(Set magnetic field derivative for overlapping fields.

An overlapping field means that the derivative is computed but that the
x-component of the jacobian is at the same position as another Jacobian
target.

The reason for this method is that it allows representing the (signed) absolute
magnetic field derivative as a combination of the three magnetic field components.

To call this method, you first have added 1 component of the magnetic field
derivative, and then you call this method to add the second and third component.
)--",
      .author = {"Richard Larsson"},
      .out    = {"jac_targets"},
      .in     = {"jac_targets"},
  };
  wsm_data["RetrievalAddOverlappingMagneticField"] =
      jac2ret("jac_targetsAddOverlappingMagneticField");

  wsm_data["jac_targetsAddWindField"] = {
      .desc      = R"--(Set wind field derivative.

Note that the derivatives from methods that takes the freqeuncy will return
their derivatives as if these were frequency derivatives.

See *FieldComponent* for valid ``component``
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets"},
      .gin       = {"component", "d"},
      .gin_type  = {"String", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0.1}},
      .gin_desc =
          {"The component to use [u, v, w]",
           "The perturbation used in methods that cannot compute derivatives analytically"},
  };
  wsm_data["RetrievalAddWindField"] = jac2ret("jac_targetsAddWindField");

  wsm_data["jac_targetsAddOverlappingWindField"] = {
      .desc   = R"--(Set wind field derivative for overlapping fields.

An overlapping field means that the derivative is computed but that the
x-component of the jacobian is at the same position as another Jacobian
target.

The reason for this method is that it allows representing the (signed) absolute
wind speed derivative as a combination of the three wind field components.

To call this method, you first have added 1 component of the wind field
derivative, and then you call this method to add the second and third component.
)--",
      .author = {"Richard Larsson"},
      .out    = {"jac_targets"},
      .in     = {"jac_targets"},
  };
  wsm_data["RetrievalAddOverlappingWindField"] =
      jac2ret("jac_targetsAddOverlappingWindField");

  wsm_data["jac_targetsAddSpeciesVMR"] = {
      .desc      = R"--(Set volume mixing ratio derivative.

See *SpeciesEnum* for valid ``species``
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets"},
      .gin       = {"species", "d"},
      .gin_type  = {"SpeciesEnum", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0.1}},
      .gin_desc =
          {"The species of interest",
           "The perturbation used in methods that cannot compute derivatives analytically"},
  };
  wsm_data["RetrievalAddSpeciesVMR"] = jac2ret("jac_targetsAddSpeciesVMR");

  wsm_data["jac_targetsAddAtmosphere"] = {
      .desc      = R"--(Sets an atmospheric target.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets"},
      .gin       = {"target", "d"},
      .gin_type  = {"AtmKey,SpeciesEnum,SpeciesIsotope,QuantumLevelIdentifier",
                    "Numeric"},
      .gin_value = {std::nullopt, Numeric{0.1}},
      .gin_desc =
          {"The target of interest",
           "The perturbation used in methods that cannot compute derivatives analytically"},
  };
  wsm_data["RetrievalAddAtmosphere"] = jac2ret("jac_targetsAddAtmosphere");

  wsm_data["jac_targetsAddSurface"] = {
      .desc      = R"--(Sets a surface target
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets"},
      .gin       = {"target", "d"},
      .gin_type  = {"SurfaceKey,SurfacePropertyTag", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0.1}},
      .gin_desc =
          {"The target of interest",
           "The perturbation used in methods that cannot compute derivatives analytically"},
  };
  wsm_data["RetrievalAddSurface"] = jac2ret("jac_targetsAddSurface");

  wsm_data["jac_targetsAddSubsurface"] = {
      .desc      = R"--(Sets a subsurface target
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets"},
      .gin       = {"target", "d"},
      .gin_type  = {"SubsurfaceKey,SubsurfacePropertyTag", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0.1}},
      .gin_desc =
          {"The target of interest",
           "The perturbation used in methods that cannot compute derivatives analytically"},
  };
  wsm_data["RetrievalAddSubsurface"] = jac2ret("jac_targetsAddSubsurface");

  wsm_data["jac_targetsAddSpeciesIsotopologueRatio"] = {
      .desc      = R"--(Set isotopologue ratio derivative

See *SpeciesIsotope* for valid ``species``
)--",
      .author    = {"Richard Larsson"},
      .out       = {"jac_targets"},
      .in        = {"jac_targets"},
      .gin       = {"species", "d"},
      .gin_type  = {"SpeciesIsotope", "Numeric"},
      .gin_value = {std::nullopt, Numeric{0.1}},
      .gin_desc =
          {"The species isotopologue of interest",
           "The perturbation used in methods that cannot compute derivatives analytically"},
  };
  wsm_data["RetrievalAddSpeciesIsotopologueRatio"] =
      jac2ret("jac_targetsAddSpeciesIsotopologueRatio");

  wsm_data["abs_bandsReadHITRAN"] = {
      .desc =
          R"--(Reads HITRAN data from a file.

The HITRAN file is assumed sorted in frequency, with each line record
filling up one line of text.

If the full 160-char line record is consumed without reaching the end of the line,
qns' and qns'' are assumed appended with default HITRANonline format.

You may pass an inclusive frequency range to limit what is read.  This will limit
the data read to the range [fmin, fmax].  All data before fmin is limited to parsing
just up until the frequency, and the database is returned if the fmax frequency is
exceeded.

The optional parameter ``einstein_coefficient`` is used to indicate if it is to
be computed from the line strength, or simply read from the Hitran data.

.. warning::
   Several HITRAN lines has Einstein coefficients that will not reproduce the results
   of pure line strength simulations.  If the option is set to read the Einstein
   coefficicent ("A") instead of computing it ("S") the program will throw an error
   if missing data is encountered.
   
.. note::
   For the computed Einstein coefficients, if the upper
   degeneracy is missing, it will be set to either - (2J+1) or -1 if J is not a local
   quantum number.  Note that this will also make the Einstein coefficient negative.
   This should not affect the simulation, but it is a warning that the data is not
   complete.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_bands"},
      .in        = {},
      .gin       = {"file",
                    "frequency_range",
                    "line_strength_option",
                    "compute_zeeman_parameters"},
      .gin_type  = {"String", "Vector2", "String", "Index"},
      .gin_value = {std::nullopt,
                    Vector2{-std::numeric_limits<Numeric>::infinity(),
                            std::numeric_limits<Numeric>::infinity()},
                    String{"S"},
                    Index{1}},
      .gin_desc =
          {"Filename",
           "Frequency range selection",
           "Whether the Hitran line strenght or the Hitran Einstein coefficient is used, the latter has historically been less reliable",
           "Compute the Zeeman parameters from the HITRAN data (will not activate Zeeman calculations, this must be done manually afterwards)"},
  };

  wsm_data["abs_bandsLineMixingAdaptation"] = {
      .desc =
          R"--(Adapts select band to use ordered Line mixing coefficients.

This is an experimental feature and might not work.

The computations of line mixing are done on the grid of temperatures provided.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_bands"},
      .in        = {"abs_bands", "abs_ecs_data", "atm_point"},
      .gin       = {"temperatures",
                    "band_key",
                    "rosenkranz_fit_order",
                    "polynomial_fit_degree"},
      .gin_type  = {"AscendingGrid", "QuantumIdentifier", "Index", "Index"},
      .gin_value = {std::nullopt, std::nullopt, Index{1}, Index{3}},
      .gin_desc =
          {"The temperatures to use for the internal fitting",
           "The band to adapt",
           "The degree of Rosenkranz coefficients (1 for just fitting y, 2 for fitting also g and dv)",
           "The highest order of the polynomial fit (2 means square, 3 means cubic, etc)"},
  };

  wsm_data["abs_bandsSelectFrequencyByLine"] = {
      .desc =
          R"--(Remove all lines that strictly falls outside a frequency range

Also remove bands whose lines are all removed.
)--",
      .author    = {"Richard Larsson", "Oliver Lemke"},
      .out       = {"abs_bands"},
      .in        = {"abs_bands"},
      .gin       = {"fmin", "fmax"},
      .gin_type  = {"Numeric", "Numeric"},
      .gin_value = {-std::numeric_limits<Numeric>::infinity(),
                    std::numeric_limits<Numeric>::infinity()},
      .gin_desc  = {"Minimum frequency to keep", "Maximum frequency to keep"},
  };

  wsm_data["abs_bandsSelectFrequencyByBand"] = {
      .desc =
          R"--(Remove all bands whose lines all strictly falls outside a frequency range
)--",
      .author    = {"Richard Larsson", "Oliver Lemke"},
      .out       = {"abs_bands"},
      .in        = {"abs_bands"},
      .gin       = {"fmin", "fmax"},
      .gin_type  = {"Numeric", "Numeric"},
      .gin_value = {-std::numeric_limits<Numeric>::infinity(),
                    std::numeric_limits<Numeric>::infinity()},
      .gin_desc  = {"Minimum frequency to keep", "Maximum frequency to keep"},
  };

  wsm_data["abs_bandsKeepID"] = {
      .desc      = R"--(Keeps first band of ID

If ``line`` is positive, also keep only the line of this index
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_bands"},
      .in        = {"abs_bands"},
      .gin       = {"id", "line"},
      .gin_type  = {"QuantumIdentifier", "Index"},
      .gin_value = {std::nullopt, Index{-1}},
      .gin_desc  = {"Band to keep", "Line to keep (if positive)"},
  };

  wsm_data["abs_lookup_dataInit"] = {
      .desc =
          R"--(Initialize an empty lookup table.
)--",
      .author = {"Richard Larsson"},
      .out    = {"abs_lookup_data"},
  };

  wsm_data["abs_lookup_dataPrecompute"] = {
      .desc =
          R"--(Precompute the lookup table for a single species, adding it to the map.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_lookup_data"},
      .in        = {"abs_lookup_data",
                    "atm_profile",
                    "freq_grid",
                    "abs_bands",
                    "abs_ecs_data",
                    "select_species"},
      .gin       = {"temperature_perturbation", "water_perturbation"},
      .gin_type  = {"AscendingGrid", "AscendingGrid"},
      .gin_value = {AscendingGrid{}, AscendingGrid{}},
      .gin_desc =
          {"Temperature perturbation to use for the lookup table",
           "Water vapor perturbation to use for the lookup table (makes the species nonlinear)"},
  };

  wsm_data["abs_lookup_dataPrecomputeAll"] = {
      .desc =
          R"--(Compute the lookup table for all species in *abs_bands*.

Wraps *abs_lookup_dataPrecompute* for each species, passing ``water_perturbation`` along
for those species that are ``water_affected_species``.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_lookup_data"},
      .in        = {"abs_lookup_data",
                    "atm_profile",
                    "freq_grid",
                    "abs_bands",
                    "abs_ecs_data"},
      .gin       = {"temperature_perturbation",
                    "water_perturbation",
                    "water_affected_species"},
      .gin_type  = {"AscendingGrid", "AscendingGrid", "ArrayOfSpeciesEnum"},
      .gin_value = {AscendingGrid{}, AscendingGrid{}, ArrayOfSpeciesEnum{}},
      .gin_desc =
          {"Temperature perturbation to use for the lookup table",
           "Water vapor perturbation to use for the lookup table",
           "A list of absorption species that are affected by water vapor perturbations nonlinearly"},
  };

  wsm_data["abs_lookup_dataFromProfiles"] = {
      .desc =
          R"--(Compute the lookup table for all species in *abs_bands*.

Wraps *abs_lookup_dataPrecomputeAll* after creating a simple
*atm_path* from the input data.

Unlike *abs_lookup_dataPrecomputeAll*, this method will initialize
*abs_lookup_data*
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_lookup_data"},
      .in        = {"freq_grid", "abs_bands", "abs_ecs_data"},
      .gin       = {"pressure_profile",
                    "temperature_profile",
                    "vmr_profiles",
                    "temperature_perturbation",
                    "water_perturbation",
                    "water_affected_species",
                    "default_isotopologue_ratios"},
      .gin_type  = {"DescendingGrid",
                    "Vector",
                    "SpeciesEnumVectors",
                    "AscendingGrid",
                    "AscendingGrid",
                    "ArrayOfSpeciesEnum",
                    "String"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    AscendingGrid{},
                    AscendingGrid{},
                    ArrayOfSpeciesEnum{},
                    String{"Builtin"}},
      .gin_desc =
          {"Pressure profile [Pa]",
           "Temperature profile [K]",
           "Volume mixing ratio profiles {SpeciesEnum: [VMR]}",
           "Temperature perturbation to use for the lookup table",
           "Water vapor perturbation to use for the lookup table",
           "A list of absorption species that are affected by water vapor perturbations nonlinearly",
           "Default isotopologue ratio option to initialize the *AtmPoint* with"},
  };

  wsm_data["abs_lookup_dataSimpleWide"] = {
      .desc =
          R"--(Set up a simple wide lookup table for all species in *abs_bands*.

This method simply computes the profiles for Earth-like atmospheres (by defaults)
and pass them into *abs_lookup_dataFromProfiles*.

The pressure range is set up logarithmically and all other ranges are set linearly.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_lookup_data"},
      .in        = {"freq_grid", "abs_bands", "abs_ecs_data"},
      .gin       = {"water_affected_species",
                    "pressure_range",
                    "temperature_range",
                    "water_vmr_range",
                    "isoratio_option",
                    "vmr_value",
                    "atmospheric_steps",
                    "temperature_perturbation_steps",
                    "water_vmr_perturbation_steps"},
      .gin_type  = {"ArrayOfSpeciesEnum",
                    "Vector2",
                    "Vector2",
                    "Vector2",
                    "String",
                    "Numeric",
                    "Index",
                    "Index",
                    "Index"},
      .gin_value = {ArrayOfSpeciesEnum{},
                    Vector2{1e-2, 1100e2},
                    Vector2{150, 350},
                    Vector2{1e-4, 0.15},
                    String{"Builtin"},
                    Numeric{1e-9},
                    Index{80},
                    Index{15},
                    Index{15}},
      .gin_desc =
          {"A list of absorption species that are affected by water vapor perturbations nonlinearly",
           "Pressure range to consider - in increasing order [Pa]",
           "Temperature range to consider - in increasing order [K]",
           "Water VMR range to consider - in increasing order [vmr]",
           "Default isotopologue ratio option to initialize the *AtmPoint* with",
           "The VMR to use for the self-value broadening",
           "Number of steps in the atmospheric profile",
           "Number of steps in the temperature perturbation",
           "Number of steps in the water vapor perturbation"},
  };

  wsm_data["abs_bandsReadSpeciesSplitCatalog"] = {
      .desc      = R"--(Reads all species in *abs_species* from a basename

basename follows the standard ARTS rules.
For example if *abs_species* contains only ``H2O-161``, then a
basename of ``"lbl"`` will read the file ``"lbl.H2O-161.xml"``, and a
basename of ``"lbl/"`` will read the file ``"lbl/H2O-161.xml"``.

ignore_missing is a boolean that indicates if the method should
ignore missing files or not.  If set to true, the method will
ignore missing files and continue.  If set to false, the method
will throw an error if any file is missing.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_bands"},
      .in        = {"abs_species"},
      .gin       = {"basename", "ignore_missing"},
      .gin_type  = {"String", "Index"},
      .gin_value = {std::nullopt, Index{0}},
      .gin_desc  = {"Absolute or relative path to the directory",
                    "Ignore missing files instead of throwing an error"},
  };

  wsm_data["abs_bandsReadSpeciesSplitARTSCAT"] = {
      .desc =
          R"--(Same as *abs_bandsReadSpeciesSplitCatalog* but for reading the old ARTSCAT format.

One key difference is that ARTSCAT were often stored in a single file per *SpeciesEnum*
rather than per *SpeciesIsotope*, so the optional argument ``pure_species`` is available
and evaluates to true by default to use *SpeciesEnum*.  Switch this off to use the
*SpeciesIsotope* instead.

Note that ARTSCAT does not support many of the features of the modern line catalog format.
This reading routine is intended for use-as-is of the produced *abs_bands*.  Or after
minor changes, like removing absorption lines outside of some frequency span.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_bands"},
      .in        = {"abs_species"},
      .gin       = {"basename", "ignore_missing", "pure_species"},
      .gin_type  = {"String", "Index", "Index"},
      .gin_value = {std::nullopt, Index{0}, Index{1}},
      .gin_desc =
          {"Absolute or relative path to the directory",
           "Flag to ignore missing files instead of throwing an error.",
           "Flag that when true uses *SpeciesEnum* instead of *SpeciesIsotope* for file names."},
  };

  wsm_data["abs_bandsReadSplit"] = {
      .desc =
          R"--(Reads all xml-files in a given directory and puts them into *abs_bands*.

.. note::
    The ``dir`` path has to be absolute or relative to the working path, the environment
    variables are not considered
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_bands"},
      .gin       = {"dir"},
      .gin_type  = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"Absolute or relative path to the directory"},
  };

  wsm_data["abs_bandsSaveSplit"] = {
      .desc      = R"--(Saves all bands in *abs_bands* to a directory

This will create the directory if it does not exist.  It will also create
subdirectories that are the short-form of the isotopologue names.  The bands
will be stored as H2O-161.xml, H2O-162.xml, O2-66.xml, and so on

The ``dir`` path has to be absolute or relative to the working path, the environment
variables are not considered
)--",
      .author    = {"Richard Larsson"},
      .in        = {"abs_bands"},
      .gin       = {"dir"},
      .gin_type  = {"String"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"Absolute or relative path to the directory"},
  };

  wsm_data["ray_pathAddGeometricGridCrossings"] = {
      .desc =
          R"--(Fill the path with with points that crosses the grid of the atmspheric field.

The atmospheric field parameter must be gridded.  Only grids with size() > 1 are considered.

Points are added where the ray path crosses any of the three grids in pure geometrical manner.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"ray_path"},
      .in        = {"ray_path", "atm_field", "surf_field"},
      .gin       = {"atm_key"},
      .gin_type  = {"AtmKey"},
      .gin_value = {AtmKey::t},
      .gin_desc =
          {"The atmospheric field key for which the grid is expected if adding grid crossings is desired"},
  };

  wsm_data["ray_pathFillGeometricHalfStep"] = {
      .desc =
          R"--(Fill the path with geometric step points.

If two path points are more than *max_stepsize* apart, additional points are added
at half the distance between these two points.

This process is repeated until there are no more neighboring points for which the premise is true.
)--",
      .author = {"Richard Larsson"},
      .out    = {"ray_path"},
      .in     = {"ray_path", "surf_field", "max_stepsize"},
  };

  wsm_data["ray_pathFillGeometricStepwise"] = {
      .desc =
          R"--(Fill the path with geometric step points.

If two path points are more than *max_stepsize* apart, additional points are added
by propagating one of the points towards the other with a step length of *max_stepsize*.

This process is repeated until there are no more neighboring points for which the premise is true.
)--",
      .author = {"Richard Larsson"},
      .out    = {"ray_path"},
      .in     = {"ray_path", "surf_field", "max_stepsize"},
  };

  wsm_data["ray_pathFixUpdownAzimuth"] = {
      .desc =
          R"--(Fix azimuth angle errors that can occur for 180 and 0 degrees zenith.

These only matter for polarized radiative transfer.
)--",
      .author = {"Richard Larsson"},
      .out    = {"ray_path"},
      .in     = {"ray_path"},
  };

  wsm_data["ray_pathAddLimbPoint"] = {
      .desc =
          R"--(Add the limb point to the ray path
)--",
      .author = {"Richard Larsson"},
      .out    = {"ray_path"},
      .in     = {"ray_path", "surf_field"},
  };

  wsm_data["ray_pathRemoveNonAtm"] = {
      .desc =
          R"--(Remove non-atmospheric points to the ray path
)--",
      .author = {"Richard Larsson"},
      .out    = {"ray_path"},
      .in     = {"ray_path"},
  };

  wsm_data["ray_pathInit"] = {
      .desc =
          R"--(Initialize the ray path with a single point.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"ray_path"},
      .in        = {"atm_field", "surf_field"},
      .gin       = {"pos", "los", "as_sensor"},
      .gin_type  = {"Vector3", "Vector2", "Index"},
      .gin_value = {std::nullopt, std::nullopt, Index{1}},
      .gin_desc =
          {"The start position",
           "The start line-of-sight",
           "Whether or not the position is the sensor position or the observer position"},
  };

  wsm_data["ray_pathRemoveNearby"] = {
      .desc =
          R"--(Remove points that are too close to each other.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"ray_path"},
      .in        = {"ray_path", "surf_field"},
      .gin       = {"min_distance", "first"},
      .gin_type  = {"Numeric", "Index"},
      .gin_value = {std::nullopt, Index{0}},
      .gin_desc  = {"The minimum distance between points",
                    "Whether to remove the first or second point"},
  };

  wsm_data["ray_pathRemoveNonGeometricGridCrossings"] = {
      .desc =
          R"--(Remove all non-geometric grid crossings from the ray path.

The atmospheric field parameter must be gridded. All points overlapping with any of the three
grids are kept.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"ray_path"},
      .in        = {"ray_path", "atm_field"},
      .gin       = {"atm_key"},
      .gin_type  = {"AtmKey"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"The atmospheric key"},
  };

  wsm_data["ray_pathSetGeometricExtremes"] = {
      .desc =
          R"--(Add the geometric extremes to the ray path.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"ray_path"},
      .in        = {"ray_path", "atm_field", "surf_field"},
      .gin       = {"surf_search_accuracy", "surf_safe_search"},
      .gin_type  = {"Numeric", "Index"},
      .gin_value = {Numeric{0.1}, Index{1}},
      .gin_desc =
          {"The accuracy within which the surface intersection is counted as a hit",
           "Whether or not to search for the surface intersection in a safer but slower manner"},
  };

  wsm_data["ray_path_observer_agendaSetGeometric"] = {
      .desc =
          R"--(Set *ray_path_observer_agenda* from programmable geometric settings.

The default settings essentially call the default settings for *ray_pathGeometric*.

Options:

- ``max_step_option``: Choose the maximum distance between two points. Set *max_stepsize* for the distance.
- ``surf_search_accuracy`` and ``surf_safe_search``: The accuracy to search for
  surface intersections and whether or not to do it at all. 
- ``remove_nearby`` and ``remove_nearby_first``: The minimum distance between points, ignored if 0 or less.
  The second option tells which point to remove if they are too close.
- ``atm_key`` and ``add_crossings`` and ``remove_non_crossings``: The atmospheric field key for which the
  grid is expected if adding grid crossings is desired.  The other two options tell whether to add all grid
  points or remove non-crossings.  The removal happens after the filling of the path.
- ``fix_updown_azimuth``: Fix the azimuth angle when looking at 0 or 180 degrees.
- ``add_limb``: Add the limb point.
- ``remove_non_atm``: Remove points in space or in the subsurface.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"ray_path_observer_agenda"},
      .gin       = {"max_step_option",
                    "surf_search_accuracy",
                    "remove_nearby",
                    "atm_key",
                    "surf_safe_search",
                    "remove_nearby_first",
                    "add_crossings",
                    "remove_non_crossings",
                    "fix_updown_azimuth",
                    "add_limb",
                    "remove_non_atm"},
      .gin_type  = {"String",
                    "Numeric",
                    "Numeric",
                    "AtmKey",
                    "Index",
                    "Index",
                    "Index",
                    "Index",
                    "Index",
                    "Index",
                    "Index"},
      .gin_value = {String{"step"},
                    Numeric{0.1},
                    Numeric{0.0},
                    AtmKey::t,
                    Index{1},
                    Index{1},
                    Index{0},
                    Index{0},
                    Index{1},
                    Index{0},
                    Index{1}},
      .gin_desc =
          {"Option for max stepping.  See *ray_path_observer_agendaSetGeometricMaxStep*",
           "The accuracy to search for surface intersections",
           "The minimum distance between points, ignroed if 0 or less",
           "The atmospheric field key for which the grid is expected if adding grid crossings is desired",
           "Whether or not to search for the surface intersection in a safer but slower manner",
           "Which point (first or second) to remove if they are too close",
           "Add all grid crossings",
           "Remove non-crossings",
           "Fix the azimuth angle when looking at 0 or 180 degrees",
           "Add the limb point",
           "Remove non-atmospheric points"},
  };

  wsm_data["ray_pathGeometricUplooking"] = {
      .desc =
          R"--(Wraps *ray_pathGeometric* for straight uplooking paths from the surface altitude at the position
)--",
      .author = {"Richard Larsson"},
      .out    = {"ray_path"},
      .in     = {"atm_field", "surf_field", "lat", "lon", "max_stepsize"},
  };

  wsm_data["ray_pathGeometricDownlooking"] = {
      .desc =
          R"--(Wraps *ray_pathGeometric* for straight downlooking paths from the top-of-the-atmosphere altitude
)--",
      .author = {"Richard Larsson"},
      .out    = {"ray_path"},
      .in     = {"atm_field", "surf_field", "lat", "lon", "max_stepsize"},
  };

  wsm_data["ray_pathGeometric"] = {
      .desc      = R"--(Get a geometric radiation path

The path is defined by the origo and the line of sight.

The ``pos`` is either at the end or at the beginning of the path depending 
on the ``as_observer`` flag.  A value that evaluates to true means that it is
at the end of the path.  If ``as_observer`` is true, the ``los`` is therefore
looking backwards along the path.  Basically, ``as_observer`` true means that
``pos`` and ``los`` behaves as sensor pos and los.

The *max_stepsize* is the maximum step length in meters.  The path is first
created between the two extremes of either space and/or surface.  Afterwards,
there are additional points added every *max_stepsize* meters between these
points until no more fits (the last step is shorter or exactly *max_stepsize*).

Upon closing the method, the following options are available to modify
the output:

If ``add_limb`` is true, the limb point is added to the path at the end.  It
is computed using bisections to ensure that the zenith angle of the tangent
point is as close to 90 degrees as it can numerically be.

If ``remove_non_atm`` is true, all points that are not in the atmosphere are
removed.  It is recommended to remove these points as multiple methods will
either perform poorly or not at all with these points present.

If ``fix_updown_azimuth`` is true, the azimuthal angle of the path is
fixed to the initial azimuthal angle of the path.  Because calculations
of the azimuth angle makes use of IEEE atan2, some paths may produce
bad angles if this is turned off.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"ray_path"},
      .in        = {"atm_field", "surf_field", "max_stepsize"},
      .gin       = {"pos",
                    "los",
                    "surf_search_accuracy",
                    "as_observer",
                    "add_limb",
                    "remove_non_atm",
                    "fix_updown_azimuth",
                    "surf_safe_search"},
      .gin_type  = {"Vector3",
                    "Vector2",
                    "Numeric",
                    "Index",
                    "Index",
                    "Index",
                    "Index",
                    "Index"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    Numeric{0.1},
                    Index{1},
                    Index{0},
                    Index{1},
                    Index{1},
                    Index{1}},
      .gin_desc =
          {"The origo of the radiation path",
           "The line of sight of the radiation path",
           "The accuracy within which the surface intersection is counted as a hit",
           "Whether or not the path is as seen by the sensor or by the radiation (see text)",
           "Wheter or not to add the limb point",
           "Wheter or not to keep only atmospheric points",
           "Whether or not to attempt fix a potential issue with the path azimuthal angle",
           "Whether or not to search for the surface intersection in a safer but slower manner"},
  };

  wsm_data["ray_pathFromPointAndDepth"] = {
      .desc      = R"--(Create a depth profile ray path from a point.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"ray_path"},
      .in        = {"ray_point"},
      .gin       = {"depth_profile"},
      .gin_type  = {"DescendingGrid"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"List of depths"},
  };

  wsm_data["spectral_rad_operatorClearsky1D"] = {
      .desc           = R"--(Set up a 1D spectral radiance operator

The operator is set up to compute the spectral radiance at any point as seen from
a 1D atmospheric profile.

This method will share line-by-line,cross-section, collision-induced absorption, and
predefined model data with the workspace (if they exist already when this method is
called).
)--",
      .author         = {"Richard Larsson"},
      .out            = {"spectral_rad_operator"},
      .in             = {"atm_field", "surf_field", "alt_grid", "lat", "lon"},
      .gin            = {"cia_extrapolation", "cia_robust"},
      .gin_type       = {"Numeric", "Index"},
      .gin_value      = {Numeric{0.0}, Index{0}},
      .gin_desc       = {"The extrapolation distance for cia",
                         "The robustness of the cia extrapolation"},
      .pass_workspace = true,
  };

  wsm_data["spectral_rad_fieldProfilePseudo2D"] = {
      .desc =
          R"--(Computes the spectral radiance field assuming a profile and a pseudo-2D path.

A profile is defined as having space blackbody emission from the top and surface temperature
blackbody emissision from the surface.

Limb paths are only considered when the zenith angle misses the next lower level using the
same mechanism as in *za_gridProfilePseudo2D*.
)--",
      .author         = {"Richard Larsson"},
      .out            = {"spectral_rad_field"},
      .in             = {"spectral_propmat_agenda",
                         "atm_profile",
                         "surf_field",
                         "freq_grid",
                         "za_grid",
                         "alt_grid",
                         "lat",
                         "lon"},
      .gin            = {"azimuth"},
      .gin_type       = {"Numeric"},
      .gin_value      = {Numeric{0}},
      .gin_desc       = {"The azimuth"},
      .pass_workspace = true,
  };

  wsm_data["za_gridProfilePseudo2D"] = {
      .desc =
          R"--(A custom zenith grid for *spectral_rad_fieldProfilePseudo2D*
)--",
      .author    = {"Richard Larsson"},
      .out       = {"za_grid"},
      .in        = {"surf_field", "alt_grid", "lat", "lon"},
      .gin       = {"dza", "azimuth", "consider_limb"},
      .gin_type  = {"Numeric", "Numeric", "Index"},
      .gin_value = {Numeric{1}, Numeric{0}, Index{1}},
      .gin_desc  = {"The zenith grid max step size",
                    "The azimuth",
                    "Whether or not special care is given to the limb"},
  };

  wsm_data["spectral_rad_fieldFromOperatorPlanarGeometric"] = {
      .desc =
          R"--(Computes the spectral radiance field assuming planar geometric paths

A planar geometric path is just defined by a 1D atmospheric profile.  If the
*spectral_rad_operator* contains more than one latitude and/or longitude
point, their altitude profiles are treated independently.

Limitations:

- The zenith grid is not allowed to contain the value 90 degrees.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"spectral_rad_field"},
      .in        = {"spectral_rad_operator", "freq_grid", "za_grid"},
      .gin       = {"azimuth_grid"},
      .gin_type  = {"AzimuthGrid"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"The azimuth grid"},
  };

  wsm_data["spectral_rad_fieldFromOperatorPath"] = {
      .desc =
          R"--(Computes the spectral radiance field using *ray_path_observer_agenda*.

Each point is in computed individually, so there will be
zenith x azimuth x altitude x latitude x longitude x frequency number of calculations.
The positional arguments are taken from *spectral_rad_operator*.

If the code is not already in parallel operation mode when this method is called,
the first 5 dimensions are computed in parallel.
)--",
      .author         = {"Richard Larsson"},
      .out            = {"spectral_rad_field"},
      .in             = {"spectral_rad_operator",
                         "ray_path_observer_agenda",
                         "freq_grid",
                         "za_grid"},
      .gin            = {"azimuth_grid"},
      .gin_type       = {"AzimuthGrid"},
      .gin_value      = {std::nullopt},
      .gin_desc       = {"The azimuth grid"},
      .pass_workspace = true,
  };

  wsm_data["atm_fieldAppendBaseData"] = {
      .desc      = R"--(Append base data to the atmospheric field

This will look at the valid ``basename`` for files matching base
data.  The base data file names are of the form

- "<...>t.xml"
- "<...>p.xml"
- "<...>wind_u.xml"
- "<...>wind_v.xml"
- "<...>wind_w.xml"
- "<...>mag_u.xml"
- "<...>mag_v.xml"
- "<...>mag_w.xml"

If any of these files are found, they are appended to the atmospheric field.

See *InterpolationExtrapolation* for valid ``extrapolation``.

See *MissingFieldComponentError* for valid ``deal_with_field_component``.

The ``replace_existing`` is used to determine if the data should be replaced if it already
exists in the atmospheric field.

The ``allow_missing_pressure`` and ``allow_missing_temperature`` are used to determine
if the method should throw if the pressure or temperature is missing.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field"},
      .gin       = {"basename",
                    "extrapolation",
                    "deal_with_field_component",
                    "replace_existing",
                    "allow_missing_pressure",
                    "allow_missing_temperature"},
      .gin_type  = {"String", "String", "String", "Index", "Index", "Index"},
      .gin_value = {std::nullopt,
                    String{"Linear"},
                    String{"Throw"},
                    Index{1},
                    Index{0},
                    Index{0}},
      .gin_desc  = {"The base name of the files",
                    "The extrapolation to use.",
                    "How to deal with the field component.",
                    "Whether or not to replace existing data",
                    "Whether or not to allow missing pressure data",
                    "Whether or not to allow missing temperature data"},
  };

  wsm_data["atm_fieldAppendLineSpeciesData"] = {
      .desc =
          R"--(Append species data to the atmospheric field based on line data.

This will look at the valid ``basename`` for files matching base
data.  The base data file names are of the short-name form: "species.xml" (e.g., "H2O.xml").
See *SpeciesEnum* for valid short names.

See *InterpolationExtrapolation* for valid ``extrapolation``.

The ``missing_is_zero`` sets missing data to zero.

The ``replace_existing`` is used to determine if the data should be replaced if it already
exists in the atmospheric field.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field", "abs_bands"},
      .gin       = {"basename",
                    "extrapolation",
                    "missing_is_zero",
                    "replace_existing"},
      .gin_type  = {"String", "String", "Index", "Index"},
      .gin_value = {std::nullopt, String{"Linear"}, Index{0}, Index{0}},
      .gin_desc  = {"The base name of the files",
                    "The extrapolation to use.",
                    "Whether or not to zero-out missing data",
                    "Whether or not to replace existing data"},
  };

  wsm_data["atm_fieldAppendLineIsotopologueData"] = {
      .desc =
          R"--(Append isotopologue ratio data to the atmospheric field based on line data.

This will look at the valid ``basename`` for files matching base
data.  The base data file names are of the form: "species-n.xml" (e.g., "H2O-161.xml").
See *abs_speciesSet* for valid isotopologue names.

See *InterpolationExtrapolation* for valid ``extrapolation``.

The ``missing_is_zero`` sets missing data to zero.

The ``replace_existing`` is used to determine if the data should be replaced if it already
exists in the atmospheric field.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field", "abs_bands"},
      .gin       = {"basename",
                    "extrapolation",
                    "missing_is_zero",
                    "replace_existing"},
      .gin_type  = {"String", "String", "Index", "Index"},
      .gin_value = {std::nullopt, String{"Linear"}, Index{0}, Index{0}},
      .gin_desc  = {"The base name of the files",
                    "The extrapolation to use.",
                    "Whether or not to zero-out missing data",
                    "Whether or not to replace existing data"},
  };

  wsm_data["atm_fieldAppendLineLevelData"] = {
      .desc =
          R"--(Append NLTE data to the atmospheric field based on line data.

This will look at the valid ``basename`` for files matching base
data.  The base data file names are of the form: "species-n QN1 N1 N1 QN2 N2 N2.xml" (e.g., "O2-66 J 1 1 N 0 0.xml").
See *abs_speciesSet* for valid isotopologue names and
*QuantumLevelIdentifier* for constructing quantum numbers identifiers.

See *InterpolationExtrapolation* for valid ``extrapolation``.

The ``missing_is_zero`` sets missing data to zero.

The ``replace_existing`` is used to determine if the data should be replaced if it already
exists in the atmospheric field.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field", "abs_bands"},
      .gin       = {"basename",
                    "extrapolation",
                    "missing_is_zero",
                    "replace_existing"},
      .gin_type  = {"String", "String", "Index", "Index"},
      .gin_value = {std::nullopt, String{"Linear"}, Index{0}, Index{0}},
      .gin_desc  = {"The base name of the files",
                    "The extrapolation to use.",
                    "Whether or not to zero-out missing data",
                    "Whether or not to replace existing data"},
  };

  wsm_data["atm_fieldAppendTagsSpeciesData"] = {
      .desc =
          R"--(Append species data to the atmospheric field based on *abs_species*.

This will look at the valid ``basename`` for files matching base
data.  The base data file names are of the short-name form: "species.xml" (e.g., "H2O.xml").
See *SpeciesEnum* for valid short names.

See *InterpolationExtrapolation* for valid ``extrapolation``.

The ``missing_is_zero`` sets missing data to zero.

The ``replace_existing`` is used to determine if the data should be replaced if it already
exists in the atmospheric field.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field", "abs_species"},
      .gin       = {"basename",
                    "extrapolation",
                    "missing_is_zero",
                    "replace_existing"},
      .gin_type  = {"String", "String", "Index", "Index"},
      .gin_value = {std::nullopt, String{"Linear"}, Index{0}, Index{0}},
      .gin_desc  = {"The base name of the files",
                    "The extrapolation to use.",
                    "Whether or not to zero-out missing data",
                    "Whether or not to replace existing data"},
  };

  wsm_data["atm_fieldAppendCIASpeciesData"] = {
      .desc =
          R"--(Append species data to the atmospheric field based on collision-induced absorption data.

This will look at the valid ``basename`` for files matching base
data.  The base data file names are of the short-name form: "species1.xml" "species2.xml" (e.g., "H2O.xml" "CO2.xml").
See *SpeciesEnum* for valid short names.

See *InterpolationExtrapolation* for valid ``extrapolation``.

The ``missing_is_zero`` sets missing data to zero.

The ``replace_existing`` is used to determine if the data should be replaced if it already
exists in the atmospheric field.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field", "abs_cia_data"},
      .gin       = {"basename",
                    "extrapolation",
                    "missing_is_zero",
                    "replace_existing"},
      .gin_type  = {"String", "String", "Index", "Index"},
      .gin_value = {std::nullopt, String{"Linear"}, Index{0}, Index{0}},
      .gin_desc  = {"The base name of the files",
                    "The extrapolation to use.",
                    "Whether or not to zero-out missing data",
                    "Whether or not to replace existing data"},
  };

  wsm_data["atm_fieldAppendXsecSpeciesData"] = {
      .desc =
          R"--(Append species data to the atmospheric field based on absorption cross-section fit data.

This will look at the valid ``basename`` for files matching base
data.  The base data file names are of the short-name form: "species.xml" (e.g., "H2O.xml").
See *SpeciesEnum* for valid short names.

See *InterpolationExtrapolation* for valid ``extrapolation``.

The ``missing_is_zero`` sets missing data to zero.

The ``replace_existing`` is used to determine if the data should be replaced if it already
exists in the atmospheric field.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field", "abs_xfit_data"},
      .gin       = {"basename",
                    "extrapolation",
                    "missing_is_zero",
                    "replace_existing"},
      .gin_type  = {"String", "String", "Index", "Index"},
      .gin_value = {std::nullopt, String{"Linear"}, Index{0}, Index{0}},
      .gin_desc  = {"The base name of the files",
                    "The extrapolation to use.",
                    "Whether or not to zero-out missing data",
                    "Whether or not to replace existing data"},
  };

  wsm_data["atm_fieldAppendLookupTableSpeciesData"] = {
      .desc =
          R"--(Append species data to the atmospheric field based on absorption lookup table data.

This will look at the valid ``basename`` for files matching base
data.  The base data file names are of the short-name form: "species.xml" (e.g., "H2O.xml").
See *SpeciesEnum* for valid short names.

See *InterpolationExtrapolation* for valid ``extrapolation``.

The ``missing_is_zero`` sets missing data to zero.

The ``replace_existing`` is used to determine if the data should be replaced if it already
exists in the atmospheric field.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field", "abs_lookup_data"},
      .gin       = {"basename",
                    "extrapolation",
                    "missing_is_zero",
                    "replace_existing"},
      .gin_type  = {"String", "String", "Index", "Index"},
      .gin_value = {std::nullopt, String{"Linear"}, Index{0}, Index{0}},
      .gin_desc  = {"The base name of the files",
                    "The extrapolation to use.",
                    "Whether or not to zero-out missing data",
                    "Whether or not to replace existing data"},
  };

  wsm_data["atm_fieldAppendPredefSpeciesData"] = {
      .desc =
          R"--(Append species data to the atmospheric field based on absorption predefined model data.

This will look at the valid ``basename`` for files matching base
data.  The base data file names are of the short-name form: "species-MODEL.xml" (e.g., "H2O-ForeignContCKDMT400.xml").
See *SpeciesEnum* for valid short names.  Will also append H2O VMR if available as some predefined models requires it.

See *InterpolationExtrapolation* for valid ``extrapolation``.

The ``missing_is_zero`` sets missing data to zero.

The ``replace_existing`` is used to determine if the data should be replaced if it already
exists in the atmospheric field.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field", "abs_predef_data"},
      .gin       = {"basename",
                    "extrapolation",
                    "missing_is_zero",
                    "replace_existing"},
      .gin_type  = {"String", "String", "Index", "Index"},
      .gin_value = {std::nullopt, String{"Linear"}, Index{0}, Index{0}},
      .gin_desc  = {"The base name of the files",
                    "The extrapolation to use.",
                    "Whether or not to zero-out missing data",
                    "Whether or not to replace existing data"},
  };

  wsm_data["atm_fieldAppendAuto"] = {
      .desc =
          R"--(Append data to the atmospheric field based on available absorption data.

It is recommended to use *atm_fieldRead* rather than this method directly.

This method scans available data and calls (in order) the methods below if that
data is available on the workspace.  It is not possible to reproduce this
method call by manually calling each method below because that would require
defining the relevant data fields.

Wraps:

- *atm_fieldAppendLineSpeciesData* if the workspace contains *abs_bands*
- *atm_fieldAppendLineIsotopologueData* if ``load_isot`` is true and if the workspace contains *abs_bands*
- *atm_fieldAppendLineLevelData* if ``load_nlte`` is true and if the workspace contains *abs_bands*
- *atm_fieldAppendTagsSpeciesData* if the workspace contains *abs_species*
- *atm_fieldAppendLookupTableSpeciesData* if the workspace contains *abs_lookup_data*
- *atm_fieldAppendCIASpeciesData* if the workspace contains *abs_cia_data*
- *atm_fieldAppendXsecSpeciesData* if the workspace contains *abs_xfit_data*
- *atm_fieldAppendPredefSpeciesData* if the workspace contains *abs_predef_data*

See these individually for more details.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field"},
      .gin       = {"basename",
                    "extrapolation",
                    "missing_is_zero",
                    "replace_existing",
                    "load_isot",
                    "load_nlte"},
      .gin_type  = {"String", "String", "Index", "Index", "Index", "Index"},
      .gin_value = {std::nullopt,
                    String{"Linear"},
                    Index{0},
                    Index{0},
                    Index{0},
                    Index{0}},
      .gin_desc  = {"The base name of the files",
                    "The extrapolation to use.",
                    "Whether or not to zero-out missing data",
                    "Whether or not to replace existing data",
                    "Whether or not to load isotopologue data",
                    "Whether or not to load NLTE data"},
      .pass_workspace = true,
  };

  wsm_data["ReadCatalogData"] = {
      .desc =
          R"--(Reads split catalog data from a folder structure similar to ``arts-cat-data``

Wraps:

- *abs_bandsReadSpeciesSplitCatalog* with "lines/" added to ``basename``
- *abs_cia_dataReadSpeciesSplitCatalog* with "cia/" added to ``basename``
- *abs_xfit_dataReadSpeciesSplitCatalog* with "xsec/" added to ``basename``
- *abs_predef_dataReadSpeciesSplitCatalog* with "predef/" added to ``basename`` and ``name_missing`` = 1
)--",
      .author = {"Richard Larsson"},
      .out = {"abs_predef_data", "abs_xfit_data", "abs_cia_data", "abs_bands"},
      .in  = {"abs_species"},
      .gin = {"basename", "ignore_missing"},
      .gin_type  = {"String", "Index"},
      .gin_value = {String{}, Index{0}},
      .gin_desc  = {"Absolute or relative path to the data",
                    "Ignore missing files instead of throwing an error"},
  };

  wsm_data["abs_bandsSetZeeman"] = {
      .desc = R"--(Set the Zeeman splitting for lines within the frequency range

See *SpeciesIsotope* for valid ``species``
)--",
      .author    = {"Richard Larsson"},
      .out       = {"abs_bands"},
      .in        = {"abs_bands"},
      .gin       = {"species", "fmin", "fmax", "on"},
      .gin_type  = {"SpeciesIsotope", "Numeric", "Numeric", "Index"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, Index{1}},
      .gin_desc  = {"Isotopologue of the species",
                    "Minimum line frequency to set Zeeman splitting for",
                    "Maximum line frequency to set Zeeman splitting for",
                    "On or off"},
  };

  wsm_data["ray_pointPastGeometric"] = {
      .desc =
          R"--(Gets the previous geometric point along *ray_path*
)--",
      .author    = {"Richard Larsson"},
      .out       = {"ray_point"},
      .in        = {"ray_path", "atm_field", "surf_field", "max_stepsize"},
      .gin       = {"surf_search_accuracy", "surf_safe_search"},
      .gin_type  = {"Numeric", "Index"},
      .gin_value = {Numeric{0.1}, Index{1}},
      .gin_desc =
          {"The accuracy within which the surface intersection is counted as a hit",
           "Whether or not to search for the surface intersection in a safer but slower manner"},
  };

  wsm_data["ray_pointPastRefractive"] = wsm_data["ray_pointPastGeometric"];
  wsm_data["ray_pointPastRefractive"].desc =
      R"--(Gets the previous refractive point along *ray_path*

This basically wraps *ray_pointPastGeometric* but sets the zenith angle
to the refracted zenith angle, i.e.,

.. math::
    
    \theta_{refracted} = \arcsin\left(\frac{n_{current}}{n_{next}}\sin(\theta_{current})\right)
)--";
  wsm_data["ray_pointPastRefractive"].in.emplace_back("single_dispersion");

  wsm_data["ray_pointBackground"] = {
      .desc =
          R"--(Sets *ray_point* to the expected background point of *ray_path*
)--",
      .author = {"Richard Larsson"},
      .out    = {"ray_point"},
      .in     = {"ray_path"},
  };

  wsm_data["ray_pointForeground"] = {
      .desc =
          R"--(Sets *ray_point* to the expected foreground point of *ray_path*
)--",
      .author = {"Richard Larsson"},
      .out    = {"ray_point"},
      .in     = {"ray_path"},
  };

  wsm_data["ray_pointLowestFromPath"] = {
      .desc =
          R"(Sets *ray_point* to the lowest altitude point of *ray_path*.
)",
      .author = {"Richard Larsson"},
      .out    = {"ray_point"},
      .in     = {"ray_path"},
  };

  wsm_data["ray_pointHighestFromPath"] = {
      .desc =
          R"(Sets *ray_point* to the highest altitude point of *ray_path*.
)",
      .author = {"Richard Larsson"},
      .out    = {"ray_point"},
      .in     = {"ray_path"},
  };

  wsm_data["measurement_vectorFromOperatorPath"] = {
      .desc =
          R"--(Sets measurement vector by looping over all sensor elements

The core calculations happens inside the *spectral_rad_operator*.
)--",
      .author         = {"Richard Larsson"},
      .out            = {"measurement_vector"},
      .in             = {"measurement_sensor",
                         "spectral_rad_operator",
                         "ray_path_observer_agenda"},
      .pass_workspace = true,
  };

  wsm_data["measurement_vectorFromSensor"] = {
      .desc =
          R"--(Sets measurement vector by looping over all sensor elements

The core calculations happens inside the *spectral_rad_observer_agenda*.
)--",
      .author         = {"Richard Larsson"},
      .out            = {"measurement_vector", "measurement_jacobian"},
      .in             = {"measurement_sensor",
                         "jac_targets",
                         "atm_field",
                         "surf_field",
                         "subsurf_field",
                         "spectral_rad_transform_operator",
                         "spectral_rad_observer_agenda"},
      .pass_workspace = true,
  };

  wsm_data["measurement_jacobianAtmosphereTransformation"] = {
      .desc   = "Applies transformations to the atmospheric state Jacobian\n",
      .author = {"Richard Larsson"},
      .out    = {"measurement_jacobian"},
      .in     = {"measurement_jacobian",
                 "model_state_vector",
                 "atm_field",
                 "jac_targets"},
  };

  wsm_data["measurement_jacobianSurfaceTransformation"] = {
      .desc   = "Applies transformations to the surface state Jacobian\n",
      .author = {"Richard Larsson"},
      .out    = {"measurement_jacobian"},
      .in     = {"measurement_jacobian",
                 "model_state_vector",
                 "surf_field",
                 "jac_targets"},
  };

  wsm_data["measurement_jacobianSubsurfaceTransformation"] = {
      .desc   = "Applies transformations to the subsurface state Jacobian\n",
      .author = {"Richard Larsson"},
      .out    = {"measurement_jacobian"},
      .in     = {"measurement_jacobian",
                 "model_state_vector",
                 "subsurf_field",
                 "jac_targets"},
  };

  wsm_data["measurement_jacobianBandTransformation"] = {
      .desc   = "Applies transformations to the line-by-line state Jacobian\n",
      .author = {"Richard Larsson"},
      .out    = {"measurement_jacobian"},
      .in     = {"measurement_jacobian",
                 "model_state_vector",
                 "abs_bands",
                 "jac_targets"},
  };

  wsm_data["measurement_jacobianSensorTransformation"] = {
      .desc =
          "Applies transformations to the measurement sensor state Jacobian\n",
      .author = {"Richard Larsson"},
      .out    = {"measurement_jacobian"},
      .in     = {"measurement_jacobian",
                 "model_state_vector",
                 "measurement_sensor",
                 "jac_targets"},
  };

  wsm_data["measurement_sensorFromModelState"] = {
      .desc =
          R"--(Update *measurement_sensor* from *model_state_vector*.
)--",
      .author = {"Richard Larsson"},
      .out    = {"measurement_sensor"},
      .in     = {"measurement_sensor", "model_state_vector", "jac_targets"},
  };

  wsm_data["measurement_sensorInit"] = {
      .desc =
          R"--(Initialize *measurement_sensor* to empty.
)--",
      .author = {"Richard Larsson"},
      .out    = {"measurement_sensor"},
  };

  wsm_data["measurement_sensorAddSimple"] = {
      .desc =
          R"--(Adds a sensor with a dirac channel opening around the frequency grid.

All elements share position, line-of-sight, and frequency grid.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"measurement_sensor"},
      .in        = {"measurement_sensor", "freq_grid"},
      .gin       = {"pos", "los", "pol"},
      .gin_type  = {"Vector3", "Vector2", "Stokvec"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    rtepack::to_stokvec(PolarizationChoice::I)},
      .gin_desc =
          {"A position [alt, lat, lon]",
           "A line of sight [zenith, azimuth]",
           "The polarization whos dot-product with the spectral radiance becomes the measurement"},
  };

  wsm_data["measurement_sensorAddGaussianZenith"] = {
      .desc =
          R"--(Add a sensor to *measurement_sensor* that has a Gaussian zenith response.
)--",
      .author   = {"Richard Larsson"},
      .out      = {"measurement_sensor"},
      .in       = {"measurement_sensor", "freq_grid"},
      .gin      = {"pos", "los", "pol", "dza_grid", "std_za"},
      .gin_type = {"Vector3", "Vector2", "Stokvec", "AscendingGrid", "Numeric"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    rtepack::to_stokvec(PolarizationChoice::I),
                    std::nullopt,
                    std::nullopt},
      .gin_desc =
          {"A position [alt, lat, lon]",
           "A line of sight [zenith, azimuth]",
           "The polarization whos dot-product with the spectral radiance becomes the measurement",
           "The delta zenith grid for the Gaussian response",
           "The standard deviation for the Gaussian response"},
  };

  wsm_data["measurement_sensorAddSimpleGaussian"] = {
      .desc =
          R"--(Adds a sensor with a Gaussian channel opening around the frequency grid.

All elements share position, line-of-sight, and frequency grid.

Note that this means you only get "half" a Gaussian channel for the outermost channels.

The I component's distribution is normalized to 1 or 0 by itself, while
the Q, U, and V components' hypotenuse are normalized to 1 or 0 together.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"measurement_sensor"},
      .in        = {"measurement_sensor", "freq_grid"},
      .gin       = {"std", "pos", "los", "pol"},
      .gin_type  = {"Numeric", "Vector3", "Vector2", "Stokvec"},
      .gin_value = {std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    rtepack::to_stokvec(PolarizationChoice::I)},
      .gin_desc =
          {"The standard deviations of the channels",
           "A position [alt, lat, lon]",
           "A line of sight [zenith, azimuth]",
           "The polarization whos dot-product with the spectral radiance becomes the measurement"},
  };

  wsm_data["measurement_sensorAddVectorGaussian"] =
      wsm_data["measurement_sensorAddSimpleGaussian"];
  wsm_data["measurement_sensorAddVectorGaussian"].gin_type[0] = "Vector";

  wsm_data["measurement_sensorAddRawSensor"] = {
      .desc =
          R"--(Adds sensor elements from a raw perturbation of the sensor

The perturbation is a gridded field of up to 6-dimensions.  The input
frequency grid determines how many elements are added to the sensor.
The cartesian perturbation is added to the sensor's position, line of sight,
and frequency grid.

The order of the dimensions are:

1. Frequency (``"df"``)
2. Zenith angle (``"dza"``)
3. Azimuth angle (``"daa"``)
4. Altitude (``"dalt"``)
5. Latitude (``"dlat"``)
6. Longitude (``"dlon"``)

The quoted strings must be used as the grid names of the gridded field.

.. note::
    It is OK to have fewer than 6 dimensions, the missing dimensions
    will be assumed to have a size of 1.  Since the data is exhaustive,
    the missing dimensions will not affect the output.  What does matter
    is that the order of the dimensions do not change from the one above.
)--",
      .author = {"Richard Larsson"},
      .out    = {"measurement_sensor"},
      .in     = {"measurement_sensor", "freq_grid"},
      .gin    = {"pos", "los", "raw_sensor_perturbation", "normalize"},
      .gin_type =
          {"Vector3",
           "Vector2",
           "StokvecSortedGriddedField1,StokvecSortedGriddedField2,StokvecSortedGriddedField3,"
           "StokvecSortedGriddedField4,StokvecSortedGriddedField5,StokvecSortedGriddedField6,"
           "SortedGriddedField1,SortedGriddedField2,SortedGriddedField3,"
           "SortedGriddedField4,SortedGriddedField5,SortedGriddedField6",
           "Index"},
      .gin_value = {std::nullopt, std::nullopt, std::nullopt, Index{0}},
      .gin_desc =
          {"A position [alt, lat, lon]",
           "A line of sight [zenith, azimuth]",
           "The sensor perturbation grid",
           "Whether or not to normalize the perturbation to 1.0 for each element"},
  };

  wsm_data["sun_pathFromObserverAgenda"] = {
      .desc =
          R"--(Find a path that hits the sun if possible

The algorithm finds the pair of angles with the least error in regards to angular zenith and 
azimuth offset from the sun.  It uses this pair of angles to compute said path.  The algorithm
is iterative.  It first finds the geometric pair of angles pointing at the sun.  It then
computes the path, using the space-facing path point's pointing offset relative to the sun
to change the angles in the four directions (up, left, right, down) until it finds a better
solution.  If no better solution is found, the algorithm it refines the angular search to half
for every level of refinement above 1, it then stops.

Note that special care is taken to eliminate surface intersections so that part of the sun may
still be hit if it is above the horizon.  If the sun is entirerly below the horizon, the path
will point close to the horizon.

The two control parameters are the ``angle_cut`` and ``just_hit``.  The ``angle_cut`` is the limit
in degrees to which the algorithm should search for a better solution.  The ``just_hit`` is a flag
that just returns the first time a path hits the sun.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"sun_path"},
      .in        = {"surf_field", "ray_path_observer_agenda", "sun"},
      .gin       = {"pos", "angle_cut", "refinement", "just_hit"},
      .gin_type  = {"Vector3", "Numeric", "Index", "Index"},
      .gin_value = {std::nullopt, Numeric{0.0}, Index{1}, Index{0}},
      .gin_desc =
          {"An observer position [alt, lat, lon]",
           "The angle delta-cutoff in the iterative solver [0.0, ...]",
           "The refinement of the search algorithm (twice the power of this is the resultion)",
           "Whether or not it is enough to just hit the sun or if better accuracy is needed"},
      .pass_workspace = true,
  };

  wsm_data["ray_path_suns_pathFromPathObserver"] = {
      .desc =
          R"--(Wraps *sun_pathFromObserverAgenda* for all paths to all suns.
)--",
      .author = {"Richard Larsson"},
      .out    = {"ray_path_suns_path"},
      .in     = {"surf_field", "ray_path_observer_agenda", "ray_path", "suns"},
      .gin    = {"angle_cut", "refinement", "just_hit"},
      .gin_type  = {"Numeric", "Index", "Index"},
      .gin_value = {Numeric{0.0}, Index{1}, Index{0}},
      .gin_desc =
          {"The angle delta-cutoff in the iterative solver [0.0, ...]",
           "The refinement of the search algorithm (twice the power of this is the resultion)",
           "Whether or not it is enough to just hit the sun or if better accuracy is needed"},
      .pass_workspace = true,
  };

  wsm_data["spectral_radSunsOrCosmicBackground"] = {
      .desc =
          R"--(Get the spectral radiance of a sun or of the cosmic background if no sun is hit.

.. note::
    Only the first sun is used if multiple suns are defined, so it is advantageous to
    have sorted *suns* by distance before running this code.  If you only have one sun,
    this is of course not an issue but you could consider using *spectral_radSunOrCosmicBackground* instead.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad"},
      .in     = {"freq_grid", "ray_point", "suns", "surf_field"},
  };

  wsm_data["spectral_radSunOrCosmicBackground"] = {
      .desc =
          R"--(Get the spectral radiance of a sun or of the cosmic background if the sun is not hit.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad"},
      .in     = {"freq_grid", "sun_path", "sun", "surf_field"},
  };

  wsm_data["sunBlackbody"] = {
      .desc =
          R"--(Set *sun* to blackbody.

.. note::
    For a Sol-like sun there are huge differences in the UV-range
    between the actual sun spectrum and the blackbody spectrum
    with the effective temperature of the sun. The blackbody sun
    strongly overestimates the UV radiation.
)--",
      .author    = {"Jon Petersen", "Richard Larsson"},
      .out       = {"sun"},
      .in        = {"freq_grid", "lat", "lon"},
      .gin       = {"radius", "distance", "temperature"},
      .gin_type  = {"Numeric", "Numeric", "Numeric"},
      .gin_value = {6.963242e8, 1.495978707e11, 5772.0},
      .gin_desc =
          {"The radius of the sun in meter. "
           "Default is the radius of our sun. ",
           "The average distance between the sun and the planet in meter. "
           "Default value is set to 1 a.u. ",
           "The effective temperature of the suns photosphere in Kelvin. "
           "Default is the temperature of our sun - 5772 Kelvin "},
  };

  wsm_data["sunFromGrid"] = {
      .desc =
          R"(Extracts a sun spectrum from a field of such data.
          
The method allows to obtain the sun spectrum by
interpolation from a field of such data. 
The sun spectrum is expected to be stored as
the irradiance at the suns photosphere.

Unit:

- GriddedField2: [W m-2 Hz-1]

- Vector *freq_grid* [Hz]
- Vector ``stokes_dim`` [1]

Dimensions: [*freq_grid*, stokes_dim]

This method performs an interpolation onto the *freq_grid*.
The point of *freq_grid* that are outside the data frequency grid
are initialized according to planck's law of the temperature variable.
Hence, a temperature of 0 means 0s the edges of the *freq_grid*.
)",
      .author   = {"Jon Petersen", "Richard Larsson"},
      .out      = {"sun"},
      .in       = {"freq_grid", "lat", "lon"},
      .gin      = {"sun_spectrum_raw",
                   "radius",
                   "distance",
                   "temperature",
                   "description"},
      .gin_type = {"GriddedField2", "Numeric", "Numeric", "Numeric", "String"},
      .gin_value =
          {std::nullopt, 6.963242e8, 1.495978707e11, 5772.0, String{"A sun"}},
      .gin_desc =
          {"A raw spectrum",
           "The radius of the sun in meter. Default is the radius of our sun. ",
           "The average distance between the sun and the planet in meter. Default value is set to 1 a.u. ",
           "The effective temperature of the suns photosphere in Kelvin. Default is the temperature of our sun - 5772 Kelvin ",
           "A description of the sun."},
  };

  wsm_data["spectral_rad_scat_pathSunsFirstOrderRayleigh"] = {
      .desc           = R"--(Add *suns* to *spectral_rad_srcvec_path*.
)--",
      .author         = {"Richard Larsson"},
      .out            = {"spectral_rad_scat_path"},
      .in             = {"spectral_propmat_scat_path",
                         "ray_path",
                         "ray_path_suns_path",
                         "suns",
                         "jac_targets",
                         "freq_grid",
                         "atm_field",
                         "surf_field",
                         "spectral_propmat_agenda"},
      .gin            = {"depolarization_factor", "hse_derivative"},
      .gin_type       = {"Numeric", "Index"},
      .gin_value      = {Numeric{0.0}, Index{0}},
      .gin_desc       = {R"--(The depolarization factor to use.)--",
                         "Flag to compute the hypsometric distance derivatives"},
      .pass_workspace = true,
  };

  wsm_data["atm_fieldFromModelState"] = {
      .desc   = R"--(Sets *atm_field* to the state of the model.
)--",
      .author = {"Richard Larsson"},
      .out    = {"atm_field"},
      .in     = {"atm_field", "model_state_vector", "jac_targets"},
  };

  wsm_data["surf_fieldFromModelState"] = {
      .desc   = R"--(Sets *surf_field* to the state of the model.
)--",
      .author = {"Richard Larsson"},
      .out    = {"surf_field"},
      .in     = {"surf_field", "model_state_vector", "jac_targets"},
  };

  wsm_data["subsurf_fieldFromModelState"] = {
      .desc   = R"--(Sets *subsurf_field* to the state of the model.
)--",
      .author = {"Richard Larsson"},
      .out    = {"subsurf_field"},
      .in     = {"subsurf_field", "model_state_vector", "jac_targets"},
  };

  wsm_data["abs_bandsFromModelState"] = {
      .desc   = R"--(Sets *abs_bands* to the state of the model.
)--",
      .author = {"Richard Larsson"},
      .out    = {"abs_bands"},
      .in     = {"abs_bands", "model_state_vector", "jac_targets"},
  };

  wsm_data["model_state_vectorInit"] = {
      .desc =
          R"--(Sets *model_state_vector* to the size *jac_targets* demand.
)--",
      .author = {"Richard Larsson"},
      .out    = {"model_state_vector"},
      .in     = {"jac_targets"},
  };

  wsm_data["model_state_vectorPerturbations"] = {
      .desc =
          R"--(Sets *model_state_vector* to the size *jac_targets* demand.

Then fills it with the perturbations from the *jac_targets*.
)--",
      .author = {"Richard Larsson"},
      .out    = {"model_state_vector"},
      .in     = {"jac_targets"},
  };

  wsm_data["model_state_vectorFromSensor"] = {
      .desc   = R"--(Sets *model_state_vector*'s sensor part.
)--",
      .author = {"Richard Larsson"},
      .out    = {"model_state_vector"},
      .in     = {"model_state_vector", "measurement_sensor", "jac_targets"},
  };

  wsm_data["model_state_vectorFromAtmosphere"] = {
      .desc   = R"--(Sets *model_state_vector*'s atmospheric part.
)--",
      .author = {"Richard Larsson"},
      .out    = {"model_state_vector"},
      .in     = {"model_state_vector", "atm_field", "jac_targets"},
  };

  wsm_data["model_state_vectorFromSurface"] = {
      .desc   = R"--(Sets *model_state_vector*'s surface part.
)--",
      .author = {"Richard Larsson"},
      .out    = {"model_state_vector"},
      .in     = {"model_state_vector", "surf_field", "jac_targets"},
  };

  wsm_data["model_state_vectorFromSubsurface"] = {
      .desc   = R"--(Sets *model_state_vector*'s subsurface part.
)--",
      .author = {"Richard Larsson"},
      .out    = {"model_state_vector"},
      .in     = {"model_state_vector", "subsurf_field", "jac_targets"},
  };

  wsm_data["model_state_vectorFromBands"] = {
      .desc   = R"--(Sets *model_state_vector*'s absorption line part.
)--",
      .author = {"Richard Larsson"},
      .out    = {"model_state_vector"},
      .in     = {"model_state_vector", "abs_bands", "jac_targets"},
  };

  wsm_data["spectral_tramat_pathFromPath"] = {
      .desc      = R"--(Gets the transmission matrix in layers along the path.

The assumption is that each path variable forms a layer from the 
ray path.  So there is a reduction in size by one.  A demand therefore
is that there are at least 2 points in the path.

The derivatives first dimensions are also 2, the first for the derivative wrt
the level before and one for the level after.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"spectral_tramat_path", "spectral_tramat_jac_path"},
      .in        = {"spectral_propmat_path",
                    "spectral_propmat_jac_path",
                    "ray_path",
                    "atm_path",
                    "surf_field",
                    "jac_targets"},
      .gin       = {"hse_derivative"},
      .gin_type  = {"Index"},
      .gin_value = {Index{0}},
      .gin_desc  = {"Flag to compute the hypsometric distance derivatives"},
  };

  wsm_data["spectral_radStepByStepEmission"] = {
      .desc   = R"--(Gets the spectral radiance from the path.

This uses a step-by-step solver to propagate background radiation along the path.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad", "spectral_rad_jac_path"},
      .in     = {"spectral_tramat_path",
                 "spectral_tramat_cumulative_path",
                 "spectral_tramat_jac_path",
                 "spectral_rad_srcvec_path",
                 "spectral_rad_srcvec_jac_path",
                 "spectral_rad_bkg"},
  };

  wsm_data["spectral_radCumulativeTransmission"] = {
      .desc   = R"--(Gets the spectral radiance from the path transmission.

Also get the Jacobian of the spectral radiance with regards to the
path parameters.
)--",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad", "spectral_rad_jac_path"},
      .in     = {"spectral_tramat_path",
                 "spectral_tramat_cumulative_path",
                 "spectral_tramat_jac_path",
                 "spectral_rad_bkg"},
  };

  wsm_data["OEM"] = {
      .desc   = R"(Inversion by the so called optimal estimation method (OEM).

Work in progress ...

The cost function to minimise, including a normalisation with length
of *measurement_vector*, is:

.. math::
    \chi^2 = \chi^2_y + \chi^2_x

where:

.. math::
    \chi^2_y = \frac{1}{m} \left(\vec{y}-\vec{y}_f\right)^\top \mathbf{S}_\epsilon^{-1} \left(\vec{y}-\vec{y}_f\right)

.. math::
    \chi^2_x = \frac{1}{m} \left(\vec{x}-\vec{x}_a\right)^\top \mathbf{S}_a^{-1} \left(\vec{x}-\vec{x}_a\right)

where:

.. list-table::
  :widths: auto
  :align: left
  :header-rows: 1

  * - Variable
    - ARTS parameter
    - Meaning
  * - :math:`\vec{x}`
    - *model_state_vector*
    - The model state vector.  All model states that are :emphasis:`allowed` to vary.
  * - :math:`\vec{x}_a`
    - *model_state_vector_apriori*
    - The a priori model state vector.
  * - :math:`\vec{y}`
    - *measurement_vector*
    - The measurement vector. This is the measurement that the OEM is trying to fit.
  * - :math:`\vec{y}_f`
    - *measurement_vector_fitted*
    - The fitted measurement vector.  The simulated measurement vector for the model state vector.
  * - :math:`\mathbf{S}_\epsilon`
    - *measurement_vector_error_covariance_matrix*
    - The error covariance matrix of the measurement vector.
  * - :math:`\mathbf{S}_a`
    - *model_state_covariance_matrix*
    - The a priori covariance matrix of the model state vector.

The current implementation provides 3 methods for the minimization of
the cost functional: Linear, Gauss-Newton and Levenberg-Marquardt.
The Gauss-Newton minimizer attempts to find a minimum solution by 
fitting a quadratic function to the cost functional. The linear minimizer
is a special case of the Gauss-Newton method, since for a linear forward
model the exact solution of the minimization problem is obtained after
the first step. The Levenberg-Marquardt method adaptively constrains the
search region for the next iteration step by means of the so-called gamma-factor.
This makes the method more suitable for strongly non-linear problems.
If the gamma-factor is 0, Levenberg-Marquardt and Gauss-Newton method
are identical. Each minimization method (li,gn,lm) has an indirect
variant (li_cg,gn_cg,lm_cg), which uses the conjugate gradient solver
for the linear system that has to be solved in each minimzation step.
This of advantage for very large problems, that would otherwise require
the computation of expensive matrix products.

Description of the special input arguments:

    - ``method``:

        - ``"li"``: A linear problem is assumed and a single iteration is performed.
        - ``"li_cg"``: A linear problem is assumed and solved using the CG solver.
        - ``"gn"``: Non-linear, with Gauss-Newton iteration scheme.
        - ``"gn_cg"``: Non-linear, with Gauss-Newton and conjugate gradient solver.
        - ``"lm"``: Non-linear, with Levenberg-Marquardt (LM) iteration scheme.
        - ``"lm_cg"``: Non-linear, with Levenberg-Marquardt (LM) iteration scheme and conjugate gradient solver.

    - ``max_start_cost``:

      No inversion is done if the cost matching the a priori state is above
      this value. If set to a negative value, all values are accepted.
      This argument also controls if the start cost is calculated. If
      set to <= 0, the start cost in ``oem_diagnostics`` is set to NaN
      when using "li" and "gn".
    
    - ``x_norm``:

      A normalisation vector for *model_state_vector*. A normalisation of *model_state_vector* can be needed
      due to limited numerical precision. If this vector is set to be empty
      no normalisation is done (defualt case). Otherwise, this must be a
      vector with same length as *model_state_vector*, just having values above zero.
      Elementwise division between *model_state_vector* and ``x_norm`` (x./x_norm) shall give
      a vector where all values are in the order of unity. Maybe the best
      way to set ``x_norm`` is x_norm = sqrt( diag( Sx ) ).

    - ``max_iter``:

      Maximum number of iterations to perform. No effect for "li".

    - ``stop_dx``:\n"

      Iteration stop criterion. The criterion used is the same as given in Rodgers\' "Inverse Methods for Atmospheric Sounding"

    - ``lm_ga_settings``:

      Settings controlling the gamma factor, part of the "LM" method.
      This is a vector of length 6, having the elements (0-based index):

            0. Start value.
            1. Fractional decrease after succesfull iteration.
            2. Fractional increase after unsuccessful iteration.
            3. Maximum allowed value. If the value is passed, the inversion is halted.
            4. Lower treshold. If the threshold is passed, gamma is set to zero. If gamma must be increased from zero, gamma is set to this value.
            5. Gamma limit. This is an additional stop criterion. Convergence is not considered until there has been one succesful iteration having a gamma <= this value.
      
      The default setting triggers an error if "lm" is selected.

    - ``clear matrices``:

      With this flag set to 1, *measurement_jacobian* and *measurement_gain_matrix* are returned as empty matrices.

    - ``display_progress``:

      Controls if there is any screen output. The overall report level is ignored by this WSM.
)",
      .author = {"Patrick Eriksson"},
      .out    = {"model_state_vector",
                 "measurement_vector_fitted",
                 "measurement_jacobian",
                 "atm_field",
                 "abs_bands",
                 "measurement_sensor",
                 "surf_field",
                 "subsurf_field",
                 "measurement_gain_matrix"},
      .gout   = {"oem_diagnostics", "lm_ga_history", "errors"},
      .gout_type = {"Vector", "Vector", "ArrayOfString"},
      .gout_desc =
          {"Basic diagnostics of an OEM type inversion",
           "The series of gamma values for a Marquardt-levenberg inversion",
           "Errors encountered during OEM execution"},
      .in        = {"model_state_vector",
                    "measurement_vector_fitted",
                    "measurement_jacobian",
                    "atm_field",
                    "abs_bands",
                    "measurement_sensor",
                    "surf_field",
                    "subsurf_field",
                    "jac_targets",
                    "model_state_vector_apriori",
                    "model_state_covariance_matrix",
                    "measurement_vector",
                    "measurement_vector_error_covariance_matrix",
                    "inversion_iterate_agenda"},
      .gin       = {"method",
                    "max_start_cost",
                    "model_state_covariance_matrix_normalization",
                    "max_iter",
                    "stop_dx",
                    "lm_ga_settings",
                    "clear_matrices",
                    "display_progress"},
      .gin_type  = {"String",
                    "Numeric",
                    "Vector",
                    "Index",
                    "Numeric",
                    "Vector",
                    "Index",
                    "Index"},
      .gin_value = {std::nullopt,
                    Numeric{std::numeric_limits<Numeric>::infinity()},
                    Vector{},
                    Index{10},
                    Numeric{0.01},
                    Vector{},
                    Index{0},
                    Index{0}},
      .gin_desc =
          {"Iteration method. For this and all options below, see further above",
           "Maximum allowed value of cost function at start",
           "Normalisation of Sx",
           "Maximum number of iterations",
           "Stop criterion for iterative inversions",
           "Settings associated with the ga factor of the LM method",
           "An option to save memory",
           "Flag to control if inversion diagnostics shall be printed on the screen"},
      .pass_workspace = true,
  };

  wsm_data["measurement_vector_error_covariance_matrix_observation_systemCalc"] = {
      .desc =
          R"(Calculates the covariance matrix describing the error due to uncertainties in the observation system.

The uncertainties of the observation system are
described by *measurement_vector_error_covariance_matrix*,
which must be set by the user to include the
relevant contributions from the measurement and the forward model.

Prerequisite for the calculation of
``measurement_vector_error_covariance_matrix_observation_system`` is a successful *OEM*
computation where also the gain matrix has been computed.
)",
      .author = {"Simon Pfreundschuh"},
      .gout = {"measurement_vector_error_covariance_matrix_observation_system"},
      .gout_type = {"Matrix"},
      .gout_desc =
          {"Covariance matrix describing the retrieval error due to uncertainties of the observation system."},
      .in = {"measurement_gain_matrix",
             "measurement_vector_error_covariance_matrix"},
  };

  wsm_data["model_state_covariance_matrix_smoothing_errorCalc"] = {
      .desc =
          R"(Calculates the covariance matrix describing the error due to smoothing.
          
The calculation of ``model_state_covariance_matrix_smoothing_error``
also requires the averaging kernel matrix *measurement_averaging_kernel*
to be computed after a successful OEM calculation.
)",
      .author    = {"Simon Pfreundschuh"},
      .gout      = {"model_state_covariance_matrix_smoothing_error"},
      .gout_type = {"Matrix"},
      .gout_desc =
          {"Covariance matrix describing the retrieval error due to smoothing."},
      .in = {"measurement_averaging_kernel", "model_state_covariance_matrix"},
  };

  wsm_data["measurement_averaging_kernelCalc"] = {
      .desc =
          R"(Calculate the averaging kernel matrix.

This is done by describing the sensitivity of the
*OEM* retrieval with respect to the true state of the system. A prerequisite
for the calculation of the averaging kernel matrix is a successful *OEM*
calculation in which the *measurement_jacobian* and the gain matrix *measurement_gain_matrix* have been calculated.
)",
      .author = {"Simon Pfreundschuh"},
      .out    = {"measurement_averaging_kernel"},
      .in     = {"measurement_gain_matrix", "measurement_jacobian"},
  };

  wsm_data["model_state_vector_aprioriFromState"] = {
      .desc =
          R"(Sets the a priori state of the model state vector to the current state.
)",
      .author = {"Richard Larsson"},
      .out    = {"model_state_vector_apriori"},
      .in     = {"model_state_vector"},
  };

  wsm_data["measurement_vector_fittedFromMeasurement"] = {
      .desc =
          R"(Sets the fitted measurement vector to the current measurement vector.
)",
      .author = {"Richard Larsson"},
      .out    = {"measurement_vector_fitted"},
      .in     = {"measurement_vector"},
  };

  wsm_data["model_state_covariance_matrixInit"] = {
      .desc =
          R"(Initialises the model state covariance matrix to the identity matrix.
)",
      .author = {"Richard Larsson"},
      .out    = {"model_state_covariance_matrix"},
  };

  wsm_data["model_state_covariance_matrixAddSpeciesVMR"] = {
      .desc =
          R"(Set a species model state covariance matrix element.
)",
      .author    = {"Richard Larsson"},
      .out       = {"model_state_covariance_matrix"},
      .in        = {"model_state_covariance_matrix", "jac_targets"},
      .gin       = {"species", "matrix", "inverse"},
      .gin_type  = {"SpeciesEnum", "BlockMatrix", "BlockMatrix"},
      .gin_value = {std::nullopt, std::nullopt, BlockMatrix{}},
      .gin_desc  = {"The species to set the covariance matrix for",
                    "The covariance diagoinal block matrix",
                    "The inverse covariance diagoinal block matrix"},
  };

  wsm_data["measurement_vector_errorFromModelState"] = {
      .desc =
          R"(Set the error and its Jacobian from the state of the model.
)",
      .author = {"Richard Larsson"},
      .out    = {"measurement_vector_error", "measurement_jacobian_error"},
      .in     = {"measurement_sensor", "jac_targets", "model_state_vector"},
  };

  wsm_data["measurement_vectorConditionalAddError"] = {
      .desc =
          R"(Add the measurement error to the measurement.  Conditionally, also to the Jacobian.
)",
      .author = {"Richard Larsson"},
      .out    = {"measurement_vector", "measurement_jacobian"},
      .in     = {"measurement_vector",
                 "measurement_jacobian",
                 "measurement_vector_error",
                 "measurement_jacobian_error",
                 "do_jac"},
  };

  wsm_data["measurement_vector_error_covariance_matrixConstant"] = {
      .desc =
          R"(Sets a constant measurement vector error covariance matrix.
)",
      .author    = {"Richard Larsson"},
      .out       = {"measurement_vector_error_covariance_matrix"},
      .in        = {"measurement_sensor"},
      .gin       = {"value"},
      .gin_type  = {"Numeric"},
      .gin_value = {std::nullopt},
      .gin_desc  = {"The value of the covariance matrix diagonal"},
  };

  wsm_data["scattering_speciesInit"] = {
      .desc   = R"(Initialize scattering species.
)",
      .author = {"Richard Larsson"},
      .out    = {"scattering_species"},
  };

  wsm_data["disort_settingsDownwellingObserver"] = {
      .desc =
          R"(Set the downwelling boundary condition using the observer agenda.
)",
      .author    = {"Richard Larsson"},
      .out       = {"disort_settings"},
      .in        = {"disort_settings",
                    "freq_grid",
                    "ray_path",
                    "atm_field",
                    "surf_field",
                    "subsurf_field",
                    "spectral_rad_observer_agenda"},
      .gin       = {"pol"},
      .gin_type  = {"Stokvec"},
      .gin_value = {Stokvec{1.0, 0.0, 0.0, 0.0}},
      .gin_desc =
          {"The polarization state to select.  The dot-product of this and *spectral_rad* is used."},
      .pass_workspace = true,
  };

  wsm_data["spectral_radSubsurfaceDisortEmissionWithJacobian"] = {
      .desc           = R"--(Gets the spectral radiance from the path.

The Jacobian is computed by perturbations.  Sensor and absorption data are
not considered as part of the perturbations.

The method wraps calling *spectral_radSubsurfaceDisortEmission* by perturbing
*model_state_vector* for Jacobian calculations using *model_state_vectorPerturbations*.
)--",
      .author         = {"Richard Larsson"},
      .out            = {"spectral_rad", "spectral_rad_jac"},
      .in             = {"freq_grid",
                         "atm_field",
                         "surf_field",
                         "subsurf_field",
                         "jac_targets",
                         "ray_point",
                         "disort_quadrature_dimension",
                         "disort_fourier_mode_dimension",
                         "disort_legendre_polynomial_dimension",
                         "disort_settings_agenda",
                         "disort_settings_downwelling_wrapper_agenda"},
      .gin            = {wsm_data["ray_pathFromPointAndDepth"].gin[0]},
      .gin_type       = {wsm_data["ray_pathFromPointAndDepth"].gin_type[0]},
      .gin_value      = {wsm_data["ray_pathFromPointAndDepth"].gin_value[0]},
      .gin_desc       = {wsm_data["ray_pathFromPointAndDepth"].gin_desc[0]},
      .pass_workspace = true,
  };

  wsm_data["disort_settingsOpticalThicknessFromPath"] = {
      .desc      = R"(Get optical thickness from path.
)",
      .author    = {"Richard Larsson"},
      .out       = {"disort_settings"},
      .in        = {"disort_settings", "ray_path", "spectral_propmat_path"},
      .gin       = {"min_optical_depth"},
      .gin_type  = {"Numeric"},
      .gin_value = {Numeric{1e-11}},
      .gin_desc =
          {"The minimum increase in optical thickness per level.  The DISORT algorithm employed is numerically unstable if the change between levels is too small."},
  };

  wsm_data["disort_settingsSubsurfaceScalarAbsorption"] =
      wsm_data["disort_settingsOpticalThicknessFromPath"];
  wsm_data["disort_settingsSubsurfaceScalarAbsorption"].desc =
      "Get optical thickness from subsurface path.\n";
  wsm_data["disort_settingsSubsurfaceScalarAbsorption"].in[2] =
      "subsurf_profile";

  wsm_data["disort_settings_agendaSetup"] = {
      .desc =
          R"--(Setup for Disort standard calculations.

This method allows setting up *disort_settings_agenda* by named options.
A description of the options is given below.
)--",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings_agenda"},
      .gin    = {"layer_emission_setting",
                 "scattering_setting",
                 "space_setting",
                 "sun_setting",
                 "surf_setting",
                 "surf_lambertian_value",
                 wsm_data["disort_settingsOpticalThicknessFromPath"].gin[0]},
      .gin_type =
          {"String",
           "String",
           "String",
           "String",
           "String",
           "Vector",
           wsm_data["disort_settingsOpticalThicknessFromPath"].gin_type[0]},
      .gin_value =
          {String{"LinearInTau"},
           String{"None"},
           String{"CosmicMicrowaveBackgroundRadiation"},
           String{"None"},
           String{"Thermal"},
           Vector{},
           wsm_data["disort_settingsOpticalThicknessFromPath"].gin_value[0]},
      .gin_desc =
          {"Layer emission settings",
           "Scattering settings",
           "Space settings",
           "Sun settings",
           "Surface settings",
           "Surface lambertian value (must be the size of the frequency grid; used only when surface is set to a Lambertian variant)",
           wsm_data["disort_settingsOpticalThicknessFromPath"].gin_desc[0]},
  };

  wsm_data["disort_settings_agendaSubsurfaceSetup"] = {
      .desc =
          R"--(Setup for Disort subsurface calculations.

This method allows setting up *disort_settings_agenda* by named options.
A description of the options is given below.
)--",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings_agenda"},
      .gin    = {"sun_setting",
                 wsm_data["disort_settingsSubsurfaceScalarAbsorption"].gin[0],
                 "fading_bottom"},
      .gin_type =
          {"String",
           wsm_data["disort_settingsSubsurfaceScalarAbsorption"].gin_type[0],
           "Index"},
      .gin_value =
          {String{"None"},
           wsm_data["disort_settingsSubsurfaceScalarAbsorption"].gin_value[0],
           Index{1}},
      .gin_desc =
          {"Sun settings",
           wsm_data["disort_settingsSubsurfaceScalarAbsorption"].gin_desc[0],
           "If true, the bottom layer has no input from below (i.e., no emission or transmission from below)."},
  };

  wsm_data["disort_settingsLegendreCoefficientsFromPath"] = {
      .desc   = R"(Sets the legendre coefficients from the path-variable.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings", "spectral_phamat_spectral_path"},
  };

  wsm_data["disort_settingsSingleScatteringAlbedoFromPath"] = {
      .desc   = R"(Sets the single scattering albedo from the path-variable.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings",
                 "spectral_propmat_path",
                 "spectral_propmat_scat_path",
                 "spectral_absvec_scat_path"},
  };

  wsm_data["disort_settingsNoSun"] = {
      .desc   = R"(Turns off solar radiation in Disort calculations.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings"},
  };

  wsm_data["disort_settingsSetSun"] = {
      .desc =
          R"--(Uses Set the FOV to the sun input for Disort calculations.
)--",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in = {"disort_settings", "freq_grid", "surf_field", "sun", "ray_point"},
  };

  wsm_data["disort_settingsNoLayerThermalEmission"] = {
      .desc   = R"(Turns off source radiation in Disort calculations.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings"},
  };

  wsm_data["disort_settingsLayerThermalEmissionLinearInTau"] = {
      .desc =
          R"(Use a source function that changes linearly in optical thickness.

Note that you must have set the optical thickness before calling this.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings", "atm_path", "freq_grid"},
  };

  wsm_data["disort_settingsSubsurfaceLayerThermalEmissionLinearInTau"] = {
      .desc =
          R"(Use a source function that changes linearly in optical thickness.

Note that you must have set the optical thickness before calling this.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings", "subsurf_profile", "freq_grid"},
  };

  wsm_data["disort_settingsLayerNonThermalEmissionLinearInTau"] = {
      .desc =
          R"(Same as *disort_settingsLayerThermalEmissionLinearInTau* but considers non-LTE

This is WIP and should not be used.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings",
                 "atm_path",
                 "spectral_propmat_path",
                 "spectral_srcvec_nlte_path",
                 "freq_grid"},
  };

  wsm_data["disort_settingsNoSpaceEmission"] = {
      .desc =
          R"(Turns off boundary condition from space for Disort calculations.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings"},
  };

  wsm_data["disort_settingsCosmicMicrowaveBackgroundRadiation"] = {
      .desc =
          R"(Space radiation into Disort is isotropic cosmic background radiation.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings", "freq_grid"},
  };

  wsm_data["disort_settingsSubsurfaceEmissionByTemperature"] = {
      .desc =
          R"(Subsurface boundary emission into Disort is based on temperature.

Sets both upper and lower bounds.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings", "freq_grid", "subsurf_profile"},
  };

  wsm_data["disort_settingsNoSurfaceEmission"] = {
      .desc = R"(Turns boundary condition from surface for Disort calculations.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings"},
  };

  wsm_data["disort_settingsSurfaceEmissionByTemperature"] = {
      .desc =
          R"(Surface radiation into Disort is isotropic from surface temperature.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings", "freq_grid", "ray_point", "surf_field"},
  };

  wsm_data["disort_settingsNoLegendre"] = {
      .desc   = R"(Turns off Legendre coefficients in Disort calculations.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings"},
  };

  wsm_data["disort_settingsNoFractionalScattering"] = {
      .desc   = R"(Turns off fractional scattering in Disort calculations.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings"},
  };

  wsm_data["disort_settingsNoSingleScatteringAlbedo"] = {
      .desc   = R"(Turns off single albedo scattering in Disort calculations.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings"},
  };

  wsm_data["disort_settingsSubsurfaceScalarSingleScatteringAlbedo"] = {
      .desc   = R"(Turns off single albedo scattering in Disort calculations.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings", "subsurf_profile"},
  };

  wsm_data["disort_settingsNoSurfaceScattering"] = {
      .desc   = R"(Turns off BDRF in Disort calculations.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"disort_settings"},
  };

  wsm_data["disort_settingsSurfaceLambertian"] = {
      .desc      = R"(Sets the surface to Lambertian.
)",
      .author    = {"Richard Larsson"},
      .out       = {"disort_settings"},
      .in        = {"disort_settings"},
      .gin       = {"value"},
      .gin_type  = {"Numeric,Vector"},
      .gin_value = {std::nullopt},
      .gin_desc =
          {"The value of the BDRF in all directions (Numeric for constant, Vector for spectral)"},
  };

  wsm_data["disort_settingsInit"] = {
      .desc   = R"(Perform Disort calculations for spectral flux.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_settings"},
      .in     = {"freq_grid",
                 "ray_path",
                 "disort_quadrature_dimension",
                 "disort_legendre_polynomial_dimension",
                 "disort_fourier_mode_dimension"},
  };

  wsm_data["disort_spectral_rad_fieldApplyUnit"] = {
      .desc   = R"(Convert units of the Disort spectral radiance field.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_spectral_rad_field"},
      .in     = {"disort_spectral_rad_field",
                 "ray_point",
                 "spectral_rad_transform_operator"},
  };

  wsm_data["disort_spectral_rad_fieldCalc"] = {
      .desc      = R"(Perform Disort calculations for spectral radiance.
)",
      .author    = {"Richard Larsson"},
      .out       = {"disort_spectral_rad_field", "disort_quadrature"},
      .in        = {"disort_settings"},
      .gin       = {"azimuth_grid"},
      .gin_type  = {"AzimuthGrid"},
      .gin_value = {AzimuthGrid{{0.0}}},
      .gin_desc  = {"The azimuthal angles"},
  };

#ifdef ENABLE_CDISORT
  wsm_data["disort_spectral_rad_fieldCalcCdisort"] = {
      .desc      = R"(Perform CDisort calculations for spectral radiance.

CDisort is only included for testing and comparisons with our own disort implementation.
)",
      .author    = {"Oliver Lemke"},
      .out       = {"disort_spectral_rad_field", "disort_quadrature"},
      .in        = {"disort_settings",
                    "atm_path",
                    "freq_grid_path",
                    "ray_path",
                    "surf_field"},
      .gin       = {"azimuth_grid"},
      .gin_type  = {"AzimuthGrid"},
      .gin_value = {AzimuthGrid{{0.0}}},
      .gin_desc  = {"The azimuthal angles"},
  };
#endif

  wsm_data["disort_spectral_flux_fieldCalc"] = {
      .desc   = R"(Perform Disort calculations for spectral flux.
)",
      .author = {"Richard Larsson"},
      .out    = {"disort_spectral_flux_field"},
      .in     = {"disort_settings"},
  };

  wsm_data["SpectralFluxDisort"] = {
      .desc      = R"(Integrate Disort spectral radiance.
)",
      .author    = {"Richard Larsson"},
      .gout      = {"spectral_flux_field_up", "spectral_flux_field_down"},
      .gout_type = {"Matrix", "Matrix"},
      .gout_desc = {"Upward spectral flux field",
                    "Downward spectral flux field"},
      .in        = {"disort_spectral_flux_field"},
  };

  wsm_data["spectral_radIntegrateDisort"] = {
      .desc   = R"(Integrate Disort spectral radiance.
)",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad"},
      .in     = {"disort_spectral_rad_field", "disort_quadrature"},
  };

  wsm_data["spectral_radFromDisort"] = {
      .desc =
          R"(Extract spectral radiance from the Disort field at the ray path point.
)",
      .author = {"Richard Larsson"},
      .out    = {"spectral_rad"},
      .in     = {"disort_spectral_rad_field", "ray_point"},
  };

  wsm_data["RetrievalInit"] = {
      .desc   = R"(Initialize the retrieval setup.
)",
      .author = {"Richard Larsson"},
      .out    = {"jac_targets",
                 "model_state_covariance_matrix",
                 "covariance_matrix_diagonal_blocks"},
  };

  wsm_data["RetrievalFinalizeDiagonal"] = {
      .desc   = R"(Finalize the retrieval setup.

See *jac_targetsFinalize* for more information.
)",
      .author = {"Richard Larsson"},
      .out    = {"model_state_covariance_matrix", "jac_targets"},
      .in     = {"jac_targets",
                 "covariance_matrix_diagonal_blocks",
                 "atm_field",
                 "surf_field",
                 "subsurf_field",
                 "abs_bands",
                 "measurement_sensor"},
  };

  wsm_data["abs_bandsSetNonLTE"] = {
      .desc =
          R"--(Set all bands to use non-LTE calculations.
)--",
      .author = {"Richard Larsson"},
      .out    = {"abs_bands"},
      .in     = {"abs_bands"},
  };

  wsm_data["atm_fieldInitializeNonLTE"] = {
      .desc =
          R"--(Initialize the non-LTE atmospheric field from the LTE temperature field.

Note that the bands have to be 1-line long to work.

This is because of how non-LTE is implemented in ARTS.
)--",
      .author    = {"Richard Larsson"},
      .out       = {"atm_field"},
      .in        = {"atm_field", "abs_bands"},
      .gin       = {"normalization"},
      .gin_type  = {"Numeric"},
      .gin_value = {Numeric{0.0}},
      .gin_desc =
          {"Normalization factor for the non-LTE field - all species of same isotopologue will be summed to this value (non-positive means no normalization)"},
  };

  wsm_data["ray_path_observersFieldProfilePseudo2D"] = {
      .desc =
          R"(Get a list of observer positions and line of sights to represent observing all angles of a profile.

Three observer types are added:

- Downward looking.  At the top-of-atmosphere, cover [za+e, 180] degrees zenith.
- Limb looking.  At top of the atmosphere, cover [90, za-e] degrees zenith.
- Upward looking.  At the surface, cover [0, 90] degrees zenith.

Here za is the surface tangent zenith angle from the top of the atmosphere. e indicates
the smallest possible numerical offset from that angle in the signed direction.

.. note::

    Each position has their zenith angle coverage linearly separated in degrees.
    To avoid the top-of-atmosphere limb singularity and bottom of atmosphere limb
    overlap, the limb zentih angle grid is divided into `nlimb+1` segments.
    The 90 degree angle is then discarded.

.. seealso::

    :func:`pyarts3.plots.ArrayOfPropagationPathPoint.plot` for a visualization of the geometry.
)",
      .author = {"Richard Larsson"},
      .out    = {"ray_path_observers"},
      .in =
          {"atm_field", "surf_field", "ray_path_observer_agenda", "lat", "lon"},
      .gin       = {"azimuth", "nup", "nlimb", "ndown"},
      .gin_type  = {"Numeric", "Index", "Index", "Index"},
      .gin_value = {Numeric{0.0}, std::nullopt, std::nullopt, std::nullopt},
      .gin_desc  = {"Azimuth angle for the observer",
                    "Number of upward looking observers (min 2)",
                    "Number of limb looking observers (min 2)",
                    "Number of downward looking observers (min 2)"},
      .pass_workspace = true,
  };

  wsm_data["ray_path_observersFluxProfile"] = {
      .desc =
          R"(Add :math:`n` observers per altitude point.

The number :math:`n` must be uneven and larger than 2.
)",
      .author    = {"Richard Larsson"},
      .out       = {"ray_path_observers"},
      .in        = {"atm_field"},
      .gin       = {"azimuth", "n", "atm_key"},
      .gin_type  = {"Numeric", "Index", "AtmKey"},
      .gin_value = {Numeric{0.0}, std::nullopt, AtmKey::t},
      .gin_desc =
          {"Azimuth angle for the observer",
           "Number of limb looking observers (min 2).  :math:`n` above.",
           "The altitude profile key in the atmosphere"},
  };

  wsm_data["ray_path_fieldFluxProfile"] = {
      .desc =
          R"(Adds observers that covers all zenith angles for each altitude point.

By default, up-looking from surface, downlooking from top of atmosphere and limb looking
just hitting the surface and just missing the surface are added.

In addition to these, all up-looking ppoints will have additional observers for max ``dza``
resolution and all downlooking points will have additional observers for max ``dza``
resolution.

Additional work is requires if proper coverage of the limb is required
)",
      .author         = {"Richard Larsson"},
      .out            = {"ray_path_field"},
      .in             = {"atm_field", "ray_path_observer_agenda"},
      .gin            = {"azimuth", "dza", "atm_key"},
      .gin_type       = {"Numeric", "Numeric", "AtmKey"},
      .gin_value      = {Numeric{0.0}, Numeric{180.0}, AtmKey::t},
      .gin_desc       = {"Azimuth angle for the observer",
                         "The minimum step coverage in zenith angles",
                         "The altitude profile key in the atmosphere"},
      .pass_workspace = true,
  };

  wsm_data["ray_path_fieldFromObserverAgenda"] = {
      .desc =
          R"(Create a ray path field from a set of observers.
)",
      .author         = {"Richard Larsson"},
      .out            = {"ray_path_field"},
      .in             = {"ray_path_observers", "ray_path_observer_agenda"},
      .pass_workspace = true,
  };

  wsm_data["subsurf_profileFromPath"] = {
      .desc =
          R"(Extract a subsurface profile from a ray path.
)",
      .author = {"Richard Larsson"},
      .out    = {"subsurf_profile"},
      .in     = {"subsurf_field", "ray_path"},
  };

  /* 
  MANUAL ENTRIES BELOW THIS POINT MESSES WITH THE LOGIC OF WSM_DATA
  */

  add_agenda_methods(wsm_data);

  {  // spectral_propmat_agendaAuto
    const std::string meta = "spectral_propmat_agendaAuto";
    if (wsm_data.contains(meta)) {
      throw std::runtime_error(std::format(
          "Method name collision: \"{}\" already exists in the workspace methods",
          meta));
    }

    std::vector<std::string> gin{"use_abs_lookup_data"};
    std::vector<std::string> gin_type{"Index"};
    std::vector<std::string> gin_desc{
        "Whether or not to use the lookup table instead of pure line-by-line calculations"};
    std::vector<std::optional<Wsv>> gin_value{Index{0}};
    std::vector<std::string> gout{};
    std::vector<std::string> gout_type{};
    std::vector<std::string> gout_desc{};

    for (auto& method : {"spectral_propmatInit",
                         "spectral_propmatAddCIA",
                         "spectral_propmatAddLines",
                         "spectral_propmatAddFaraday",
                         "spectral_propmatAddXsecFit",
                         "spectral_propmatAddPredefined",
                         "spectral_propmatAddLookup"}) {
      try {
        const auto& x = wsm_data.at(method);

        for (Size i = 0; i < x.gin.size(); i++) {
          if (auto ptr = std::ranges::find(gin, x.gin[i]); ptr != gin.end()) {
            gin_desc[ptr - gin.begin()] +=
                std::format(", *{}*", std::string_view{method});
          } else {
            gin.emplace_back(x.gin[i]);
            gin_type.emplace_back(x.gin_type[i]);
            gin_desc.emplace_back(
                std::format("See *{}*", std::string_view{method}));
            gin_value.emplace_back(x.gin_value[i]);
          }
        }

        for (Size i = 0; i < x.gout.size(); i++) {
          if (auto ptr = std::ranges::find(gout, x.gout[i]);
              ptr != gout.end()) {
            gout_desc[ptr - gout.begin()] +=
                std::format(", *{}*", std::string_view{method});
          } else {
            gout.emplace_back(x.gout[i]);
            gout_type.emplace_back(x.gout_type[i]);
            gout_desc.emplace_back(
                std::format("See *{}*", std::string_view{method}));
          }
        }
      } catch (std::out_of_range&) {
        throw std::runtime_error(std::format(
            R"(Missing method: "{}",
cannot generate automatic method signature for "{}"
)",
            method,
            meta));
      }
    }

    wsm_data[meta] = {
        .desc =
            R"--(Sets the *spectral_propmat_agenda* automatically from absorption data and species tag meta information.

The following methods are considered for addition to the agenda:

- *spectral_propmatAddCIA*
- *spectral_propmatAddLines*
- *spectral_propmatAddFaraday*
- *spectral_propmatAddXsecFit*
- *spectral_propmatAddPredefined*

If ``use_abs_lookup_data`` evaluates to true, lookup table
calculations, via *spectral_propmatAddLookup*, are used instead of *spectral_propmatAddLines*.

Note that the signature of this method changes depending on the input methods.  This is important
because several generic input parameters are used in the methods.  Please see the individual methods
for more information.
)--",
        .author    = {"Richard Larsson"},
        .out       = {"spectral_propmat_agenda"},
        .gout      = gout,
        .gout_type = gout_type,
        .gout_desc = gout_desc,
        .in        = {"abs_species", "abs_bands"},
        .gin       = gin,
        .gin_type  = gin_type,
        .gin_value = gin_value,
        .gin_desc  = gin_desc,
    };
  }

  /*
  LEAVE THIS LAST AS IT REQUIRES THE DATA ABOVE TO FUNCTION
  */
  for (auto& m : internal_meta_methods()) {
    if (wsm_data.contains(m.name)) {
      throw std::runtime_error(std::format(
          "Method name collision: \"{}\" already exists in the workspace methods",
          m.name));
    }
    wsm_data[m.name] = m.create(wsm_data);
  }

  fix(wsm_data);
  return wsm_data;
} catch (std::exception& e) {
  throw std::runtime_error(std::format("Cannot create workspace methods:\n\n{}",
                                       std::string_view(e.what())));
}
}  // namespace

const std::unordered_map<std::string, WorkspaceMethodInternalRecord>&
internal_workspace_methods() {
  static const auto wsm_data = internal_workspace_methods_create();
  return wsm_data;
}
