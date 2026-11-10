Adding new species to ARTS
==========================

The species, isotopologues, and predefined models that ARTS support are
defined at compile time.  This document describes the layout and the
process of expanding the list of supported species.  The process is
different for isotopologues and predefined models.  Both requires that
the species already exist.

To define words here:

1. A species is the molecule or atom of interest, e.g., H2O, O3, CO2, etc.
2. An isotopologue is a specific isotopic variant of a species, e.g.,
   H2O-161, H2O-181, O3-666, etc.
3. A predefined model is a set of parameters and equations that describe
   the behavior of a species under certain conditions, e.g., H2O-PWR98,
   CO2-CKDMT252, etc.

All three of these are defined in individual XML files that can be read
by a compiled version of ARTS.
The process of adding new information generally involves adding new
XML files to the appropriate
directories and recompiling ARTS to make use of the new information.
The details of the XML file format and the compilation process are
described in the following sections.

.. tip::
   The easist way add new data of these is to
   copy an existing files and modify the
   relevant information.  There are a lot of
   tests run during the CI process to ensure that
   the basic structure of the files is correct.

To add a species
----------------

Species are defined in ``arts-cat-data/species/``.
They must be named in an index-ordered manner, e.g., ``0001.H2O.xml``,
``0002.CO2.xml``, etc.
Missing numbers are allowed but not recommended -
it is used during the CI process.

The file format is a :class:`~pyarts3.arts.SpeciesEnumInfo` XML file.
The information inside the file is the identifying index of the species,
the short name that humans are expected to use, and a long name that is used
primarily to define the value of the species in
the :class:`~pyarts3.arts.SpeciesEnum` class.
The short name is used in the ARTS XML files to refer to the species.

The reason these are different is that some species such
as ``NO+`` cannot be used
as C++/Python identifiers, making it cumbersome to write
automatic code generation from an enumeration.

To add an isotopologue
----------------------

Isotopologues are defined in ``arts-cat-data/isotopologues/``.
They must be named by the species they contain and the isotopologue
number that identifies the specific isotopologue,
e.g., ``H2O-161.xml``, ``H2O-181.xml``, etc.

The file format is a :class:`~pyarts3.arts.SpeciesIsotopologueInfo` XML file.
The content of these files is the short name of the species,
the isotopologue index,
the mass in atomic mass units of the isotoplogue,
the default isotopologue ratio in
Earth's atmosphere, and the degeneracy of the molecular lines.

You are not done yet, however, as you also need to add the isotopologue to
the partition function files, and possibly to
the HITRAN and/or JPL compatibility files.

To add a partition function
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The partition function files are defined
in ``arts-cat-data/partition-functions/``.
They must also be named by the species they contain and the isotopologue
number that identifies the specific isotopologue,
e.g., ``H2O-161.xml``, ``H2O-181.xml``, etc.

The partition function files are
a :class:`~pyarts3.arts.PartitionFunctionsData` XML file.
This is a pure data file that contains a tag for the format of the
partition function data, and a simple :class:`~pyarts3.arts.Matrix`
to hold the data.  The supported formats are defined by
the :class:`~pyarts3.arts.PartitionFunctionsType` enumeration.
See either type for more information on the supported
formats and the expected data layout.

To add support for HITRAN
^^^^^^^^^^^^^^^^^^^^^^^^^

This is required if you want to read a HITRAN line list containing the
isotopologue you want to add.

Supported HITRAN isotopologues are defined in ``arts-cat-data/hitran/``.
They must be named by the isotopologue the add
support for, e.g., ``H2O-161.xml``, ``H2O-181.xml``, etc.

The file format is a :class:`~pyarts3.arts.HitranSpeciesInfo` XML file.
The content of these files is the full name of the isotopologue,
the index of the species in the HITRAN database, and the tag character
of the isotopologue in the HITRAN database, and the HITRAN-defined
isotopologue ratio.
The latter is required to compute the Einstein A coefficients
from the HITRAN line strengths.

.. tip::

   If you need quantum numbers, use the HITRAN online database
   to append ``qns'`` and ``qns''`` at the end of the ``par-line``
   format.

To add support for JPL
^^^^^^^^^^^^^^^^^^^^^^

This is required if you want to read a JPL line list containing the
isotopologue you want to add.

Supported JPL isotopologues are defined in ``arts-cat-data/jpl/``.
They must be named by the six-digit JPL species tag they add support
for, e.g., ``017002.xml``, ``018003.xml``, etc.
Several tags may point to the same isotopologue, e.g., for
entries of different vibrational states.

The file format is a :class:`~pyarts3.arts.JplSpeciesInfo` XML file.
The content of these files is the JPL species tag,
the full name of the isotopologue, whether the quantum number
format is supported (currently unused), the JPL partition function
at the reference temperature, and the reference temperature of
the line strengths in the JPL database.
Both of the latter are given in the documentation of the JPL entry.

The JPL partition function is used to compute the Einstein coefficient
from the JPL intensity.  JPL degeneracies and partition functions often
leave out factors that are common to all levels,
such as nuclear spin statistics.
The upper state degeneracy is therefore scaled by the ratio of the partition
function of ARTS to that of JPL at the reference temperature, so that ARTS
reproduces the JPL intensity at the reference temperature.
At other temperatures, the line strengths follow
the partition function of ARTS.
An isotopologue that only exists in the JPL catalog should use the
JPL partition function table as its partition function, with the
``LogInterp`` type, so that its line strengths follow JPL at all temperatures.

.. admonition:: FIXME
   :class: danger

   The quantum numbers are currently completely ignored.  We would
   be happy to add support for them but will likely require external
   help to do so.

To add a predefined model
-------------------------

The predefined models are not defined in the
``arts-cat-data/`` directory, but is part of the
ARTS source code tree.  They are
defined in ``arts/src/core/spec/predef-models/``.
The must be named by the species they contain and the model name, e.g.,
``H2O-PWR98.xml``, ``CO2-CKDMT252.xml``, etc.

The file format is also
a :class:`~pyarts3.arts.SpeciesIsotopologueInfo` XML file,
but only the short name of the species and the isotopologue names are used.
