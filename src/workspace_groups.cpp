#include "workspace_groups.h"

std::unordered_map<std::string, WorkspaceGroupRecord>
internal_workspace_groups_creator() {
  std::unordered_map<std::string, WorkspaceGroupRecord> wsg_data;

  wsg_data["AbsorptionLines"] = {
      .file = "absorptionlines.h",
      .desc =
          "Contains line-by-line absorption information for a number of related absorption lines\n",
  };

  wsg_data["AbsorptionBand"] = {
      .file = "lbl.h",
      .desc =
          "Contains all information about bands of related absorption lines\n",
  };

  wsg_data["AbsorptionBands"] = {
      .file     = "lbl.h",
      .desc     = "A map of *QuantumIdentifier* to *AbsorptionBand*\n",
      .map_type = true,
  };

  wsg_data["SpeciesEnumVectors"] = {
      .file     = "atm.h",
      .desc     = "A map of *SpeciesEnum* to *Vector*\n",
      .map_type = true,
  };

  wsg_data["Agenda"] = {
      .file = "workspace_agenda_class.h",
      .desc = "Describes a set of function calls and variable definitions\n",
  };

  wsg_data["Any"] = {
      .file = "supergeneric.h",
      .desc =
          "Meta type for when methods can take any argument (avoid manual use)\n",
  };

  wsg_data["ArrayOfAbsorptionLines"] = {
      .file = "absorptionlines.h",
      .desc = "A list of *AbsorptionLines*\n",
  };

  wsg_data["ArrayOfArrayOfAbsorptionLines"] = {
      .file = "absorptionlines.h",
      .desc = "A list of *ArrayOfAbsorptionLines*\n",
  };

  wsg_data["ArrayOfAgenda"] = {
      .file = "workspace_agenda_class.h",
      .desc = "A list of *Agenda*\n",
  };

  wsg_data["ArrayOfArrayOfGriddedField1"] = {
      .file = "matpack.h",
      .desc = "A list of *ArrayOfGriddedField1*\n",
  };

  wsg_data["ArrayOfArrayOfGriddedField2"] = {
      .file = "matpack.h",
      .desc = "A list of *ArrayOfGriddedField2*\n",
  };

  wsg_data["ArrayOfArrayOfGriddedField3"] = {
      .file = "matpack.h",
      .desc = "A list of *ArrayOfGriddedField3*\n",
  };

  wsg_data["ArrayOfArrayOfIndex"] = {
      .file = "matpack.h",
      .desc = "A list of *ArrayOfIndex*\n",
  };

  wsg_data["ArrayOfArrayOfMatrix"] = {
      .file = "matpack.h",
      .desc = "A list of *ArrayOfMatrix*\n",
  };

  wsg_data["ArrayOfArrayOfScatteringMetaData"] = {
      .file = "optproperties.h",
      .desc = "A list of *ArrayOfScatteringMetaData*\n",
  };

  wsg_data["ArrayOfArrayOfSingleScatteringData"] = {
      .file = "optproperties.h",
      .desc = "A list of *ArrayOfSingleScatteringData*\n",
  };

  wsg_data["ArrayOfArrayOfSpeciesTag"] = {
      .file = "species_tags.h",
      .desc = "A list of *ArrayOfSpeciesTag*\n",
  };

  wsg_data["ArrayOfArrayOfString"] = {
      .file = "mystring.h",
      .desc = "A list of *ArrayOfString*\n",
  };

  wsg_data["ArrayOfArrayOfTensor3"] = {
      .file = "matpack.h",
      .desc = "A list of *ArrayOfTensor3*\n",
  };

  wsg_data["ArrayOfArrayOfTensor6"] = {
      .file = "matpack.h",
      .desc = "A list of *ArrayOfTensor6*\n",
  };

  wsg_data["ArrayOfArrayOfTime"] = {
      .file = "artstime.h",
      .desc = "A list of *ArrayOfTime*\n",
  };

  wsg_data["ArrayOfArrayOfVector"] = {
      .file = "matpack.h",
      .desc = "A list of *ArrayOfVector*\n",
  };

  wsg_data["ArrayOfAtmPoint"] = {
      .file = "atm.h",
      .desc = "A list of *AtmPoint*\n",
  };

  wsg_data["ArrayOfCIARecord"] = {
      .file = "cia.h",
      .desc = "A list of *CIARecord*\n",
  };

  wsg_data["ArrayOfGriddedField1"] = {
      .file = "matpack.h",
      .desc = "A list of *GriddedField1*\n",
  };

  wsg_data["ArrayOfGriddedField2"] = {
      .file = "matpack.h",
      .desc = "A list of *GriddedField2*\n",
  };

  wsg_data["ArrayOfGriddedField1Named"] = {
      .file = "matpack.h",
      .desc = "A list of *GriddedField1Named*\n",
  };

  wsg_data["ArrayOfNamedGriddedField2"] = {
      .file = "matpack.h",
      .desc = "A list of *NamedGriddedField2*\n",
  };

  wsg_data["ArrayOfGriddedField3"] = {
      .file = "matpack.h",
      .desc = "A list of *GriddedField3*\n",
  };

  wsg_data["ArrayOfGriddedField4"] = {
      .file = "matpack.h",
      .desc = "A list of *GriddedField4*\n",
  };

  wsg_data["ArrayOfIndex"] = {
      .file = "matpack.h",
      .desc = "A list of *Index*\n",
  };

  wsg_data["ArrayOfMatrix"] = {
      .file = "matpack.h",
      .desc = "A list of *Matrix*\n",
  };

  wsg_data["ArrayOfQuantumIdentifier"] = {
      .file = "quantum_numbers.h",
      .desc = "A list of *QuantumIdentifier*\n",
  };

  wsg_data["ArrayOfScatteringSpecies"] = {
      .file = "scattering/scattering_species.h",
      .desc = "Represents species of scattering paritlces in the atmosphere.",
  };

  wsg_data["ScatteringSpeciesProperty"] = {
      .file = "scattering/properties.h",
      .desc = "Meta data for scattering spefcies.",
};

  wsg_data["ArrayOfScatteringMetaData"] = {
      .file = "optproperties.h",
      .desc = "A list of *ScatteringMetaData*\n",
  };

  wsg_data["ArrayOfSingleScatteringData"] = {
      .file = "optproperties.h",
      .desc = "A list of *SingleScatteringData*\n",
  };

  wsg_data["SpeciesTag"] = {
      .file = "species_tags.h",
      .desc = R"(A tagged absorption species

These tags include the species and a lot of optional information
about the isotopologue, the absorption scheme, and the frequency limits
)",
  };

  wsg_data["ArrayOfSpeciesTag"] = {
      .file = "species_tags.h",
      .desc = R"--(A list of *SpeciesTag*
)--",
  };

  wsg_data["SurfaceKey"] = {
      .file = "enumsSurfaceKey.h",
      .desc = R"--(A surface key
)--",
  };

  wsg_data["SurfaceTypeTag"] = {
      .file = "surf.h",
      .desc = R"--(A surface type
)--",
  };

  wsg_data["SurfacePropertyTag"] = {
      .file = "surf.h",
      .desc = R"--(A surface property
)--",
  };

  wsg_data["ArrayOfSpeciesEnum"] = {
      .file = "species.h",
      .desc = R"--(A list of *SpeciesEnum*
)--",
  };

  wsg_data["SpeciesEnum"] = {
      .file = "enumsSpeciesEnum.h",
      .desc = R"--(An atmospheric species
)--",
  };

  wsg_data["AtmKey"] = {
      .file = "enumsAtmKey.h",
      .desc = R"--(An atmospheric key
)--",
  };

  wsg_data["LineShapeModelVariable"] = {
      .file = "enumsLineShapeModelVariable.h",
      .desc = R"--(A line shape model parameter
)--",
  };

  wsg_data["LineByLineVariable"] = {
      .file = "enumsLineByLineVariable.h",
      .desc = R"--(An line-by-line variable parameter
)--",
  };

  wsg_data["ArrayOfSparse"] = {
      .file = "matpack_sparse.h",
      .desc = "A list of *Sparse*\n",
  };

  wsg_data["Sun"] = {
      .file = "sun.h",
      .desc = R"-x-(A single sun.
          
Each sun is described by a struct with its spectrum, radius
distance from center of planet to center of sun,
temperature (if possible), latitude in the sky of the planet,
longitude in the sky of the planet and the type)-x-",
  };

  wsg_data["ArrayOfSun"] = {
      .file = "sun.h",
      .desc = "A list of *Sun*\n",
  };

  wsg_data["ArrayOfString"] = {
      .file = "mystring.h",
      .desc = "A list of *String*\n",
  };

  wsg_data["ArrayOfTelsemAtlas"] = {
      .file = "matpack.h",
      .desc = "A list of *TelsemAtlas*\n",
  };

  wsg_data["ArrayOfTensor3"] = {
      .file = "matpack.h",
      .desc = "A list of *Tensor3*\n",
  };

  wsg_data["ArrayOfTensor4"] = {
      .file = "matpack.h",
      .desc = "A list of *Tensor4*\n",
  };

  wsg_data["ArrayOfTensor5"] = {
      .file = "matpack.h",
      .desc = "A list of *Tensor5*\n",
  };

  wsg_data["ArrayOfTensor6"] = {
      .file = "matpack.h",
      .desc = "A list of *Tensor6*\n",
  };

  wsg_data["ArrayOfTensor7"] = {
      .file = "matpack.h",
      .desc = "A list of *Tensor7*\n",
  };

  wsg_data["ArrayOfTime"] = {
      .file = "matpack.h",
      .desc = "A list of *Time*\n",
  };

  wsg_data["ArrayOfVector"] = {
      .file = "matpack.h",
      .desc = "A list of *Vector*\n",
  };

  wsg_data["XsecRecord"] = {
      .file = "xsec_fit.h",
      .desc = "A single cross-section record",
  };

  wsg_data["ArrayOfXsecRecord"] = {
      .file = "xsec_fit.h",
      .desc =
          R"--(A list of cross-section records

These cross-section records contains information about the valid temperature and
pressure ranges as well as well as the fitting coefficients used to compute
and interpolate the cross-section to other temperatures and pressures
)--",
  };

  wsg_data["AtmKey"] = {
      .file = "enumsAtmKey.h",
      .desc = R"--(A key for atmospheric data
)--",
  };

  wsg_data["AtmField"] = {
      .file = "atm.h",
      .desc = R"--(An atmospheric field
)--",
  };

  wsg_data["AtmData"] = {
      .file = "atm.h",
      .desc = R"--(An atmospheric point
)--",
  };

  wsg_data["AtmPoint"] = {
      .file = "atm.h",
      .desc = R"--(An atmospheric point
)--",
  };

  wsg_data["CIARecord"] = {
      .file = "cia.h",
      .desc =
          R"--(Contains information to compute collision induced absorption for a pair of species

Holds an the record data in a gridded field with grids of temperature and frequency in
units of m^5 molec^(-2)
)--",
  };

  wsg_data["CallbackOperator"] = {
      .file = "callback.h",
      .desc = "Used to inject custom code into *Agenda*\n",
  };

  wsg_data["BlockMatrix"] = {
      .file = "covariance_matrix.h",
      .desc = "A block matrix for the covariance matrix\n",
  };

  wsg_data["CovarianceMatrix"] = {
      .file = "covariance_matrix.h",
      .desc = "A covariance matrix\n",
  };

  wsg_data["GasAbsLookup"] = {
      .file = "gas_abs_lookup.h",
      .desc = R"--(An absorption lookup table

This class holds an absorption lookup table, as well as all
information that is necessary to use the table to extract
absorption
)--",
  };

  wsg_data["GriddedField1"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 1 dimensional gridded set of *Numeric* data

The grid is 1 *Vector*

Both the data and the grid may be named)--",
  };

  wsg_data["GriddedField2"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 2 dimensional gridded set of *Numeric* data

The grid is a combination of 2 *Vector*

Both the data and the grid may be named
)--",
  };

  wsg_data["ComplexGriddedField2"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 2 dimensional gridded set of complex data

The grid is a combination of 2 *Vector*

Both the data and the grid may be named
)--",
  };

  wsg_data["GriddedField1Named"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 3 dimensional gridded set of *Numeric* data

The grid is a combination of 1 *Vector* and 1 *ArrayOfString*

Both the data and the grid may be named
)--",
  };

  wsg_data["NamedGriddedField2"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 3 dimensional gridded set of *Numeric* data

The grid is a combination of 1 *ArrayOfString* and 2 *Vector*

Both the data and the grid may be named
)--",
  };

  wsg_data["GriddedField3"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 3 dimensional gridded set of *Numeric* data

The grid is a combination of 3 *Vector*

Both the data and the grid may be named
)--",
  };

  wsg_data["NamedGriddedField3"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 4 dimensional gridded set of *Numeric* data

The grid is a combination of 1 *ArrayOfString* and 3 *Vector*

Both the data and the grid may be named
)--",
  };

  wsg_data["GriddedField4"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 4 dimensional gridded set of *Numeric* data

The grid is a combination of 4 *Vector*

Both the data and the grid may be named
)--",
  };

  wsg_data["GriddedField5"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 5 dimensional gridded set  of *Numeric* data

The grid is a combination of 5 *Vector*

Both the data and the grid may be named
)--",
  };

  wsg_data["GriddedField6"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 6 dimensional gridded set of *Numeric* data

The grid is a combination of 6 *Vector*

Both the data and the grid may be named
)--",
  };

  wsg_data["Index"] = {
      .file       = "matpack.h",
      .desc       = "A 64 bit signed integer type\n",
      .value_type = true,
  };

  wsg_data["LinemixingEcsData"] = {
      .file = "lbl.h",
      .desc =
          R"--(A map of line mixing data
)--",
  };

  wsg_data["MCAntenna"] = {
      .file = "mc_antenna.h",
      .desc = "An antenna object used by ``MCGeneral``\n",
  };

  wsg_data["Matrix"] = {
      .file = "matpack.h",
      .desc = "A 2 dimensional array of *Numeric*\n",
  };

  wsg_data["DisortBDRF"] = {
      .file = "disort.h",
      .desc = "A bidirectional reflectance function\n",
  };

  wsg_data["MatrixOfDisortBDRF"] = {
      .file = "disort.h",
      .desc = "A 2 dimensional array of *DisortBDRF*\n",
  };

  wsg_data["Numeric"] = {
      .file       = "matpack.h",
      .desc       = "IEEE 754 binary64 floating point number\n",
      .value_type = true,
  };

  wsg_data["PredefinedModelData"] = {
      .file = "predef.h",
      .desc =
          R"--(Contains any data required for a predefined model
)--",
  };

  wsg_data["QuantumIdentifier"] = {
      .file = "quantum_numbers.h",
      .desc =
          R"--(An ID for an absorption species state

It contains information about the species and a set of quantum numbers
and can thus be used to identify one of the following:

1. a species
2. an isotopologue of a species
3. an absorption band of an isotopologue
4. an absorption line of an isotopologue
5. the energy level of absorption band(s) of an isotopologue
6. the energy level of absorption line(s) of an isotopologue
)--",
  };

  wsg_data["Rational"] = {
      .file = "matpack.h",
      .desc = "Holds a rational number as two *Index* n / d\n",
  };

  wsg_data["ScatteringMetaData"] = {
      .file = "optproperties.h",
      .desc = "Holds meta data about the scattering\n",
  };

  wsg_data["SingleScatteringData"] = {
      .file = "optproperties.h",
      .desc = "Holds single scattering data\n",
  };

  wsg_data["Sparse"] = {
      .file = "matpack_sparse.h",
      .desc = "A sparse version of *Matrix*\n",
  };

  wsg_data["String"] = {
      .file       = "mystring.h",
      .desc       = "Basic string type\n",
      .value_type = true,
  };

  wsg_data["SurfaceField"] = {
      .file = "surf.h",
      .desc =
          R"--(A surface field that keeps relevant surface parameters
)--",
  };

  wsg_data["SurfacePoint"] = {
      .file = "surf.h",
      .desc =
          R"--(A surface point.

Keeps point values for the surface, including the local normal vector.
)--",
  };

  wsg_data["TelsemAtlas"] = {
      .file = "telsem.h",
      .desc = R"--(A telsem atlas

Represents a Telsem2 atlas containing land surface microwave emissivities.
Since the Atlas contains emissivities only for land surfaces, the data is
stored in a sparse format.
 
The emissivities are represented on an equal area grid and numbered
sequentially starting with the first latitude band at -90 degrees and
moving up to 90 degrees.

The correspondance array contains the data indices for each cellnumber
if it is contained in the Atlas and NAN otherwise.
)--",
  };

  wsg_data["Tensor3"] = {
      .file = "matpack.h",
      .desc = "A 3 dimensional array of *Numeric*\n",
  };

  wsg_data["Tensor4"] = {
      .file = "matpack.h",
      .desc = "A 4 dimensional array of *Numeric*\n",
  };

  wsg_data["Tensor5"] = {
      .file = "matpack.h",
      .desc = "A 5 dimensional array of *Numeric*\n",
  };

  wsg_data["Tensor6"] = {
      .file = "matpack.h",
      .desc = "A 6 dimensional array of *Numeric*\n",
  };

  wsg_data["Tensor7"] = {
      .file = "matpack.h",
      .desc = "A 7 dimensional array of *Numeric*\n",
  };

  wsg_data["Time"] = {
      .file = "artstime.h",
      .desc = R"(Represents a time stamp
)",
  };

  wsg_data["TessemNN"] = {
      .file = "tessem.h",
      .desc = "Data required by TESSEM to calculate surface emissivity\n",
  };

  wsg_data["Vector"] = {
      .file = "matpack.h",
      .desc = "A 1 dimensional array of *Numeric*\n",
  };

  wsg_data["VibrationalEnergyLevels"] = {
      .file = "nlte.h",
      .desc = "A map of vibrational energy levels for NLTE calculations\n",
  };

  // rtepack types
  wsg_data["Propmat"] = {
      .file = "rtepack.h",
      .desc = R"--(A single propagation matrix.

Due to the properties of a propagation matrix, only 7 independents need be stored.
The propagation matrix is thus represented as:

.. math::
    \begin {array} {rrrr}
    A & B & C & D \\
    B & A & U & V \\
    C &-U & A & W \\
    D &-V &-W & A
    \end {array}

This type is related to *Stokvec* in that its first 4 elements are the same as
the first 4 elements of *Stokvec* for pure clearsky radiative transfers.

This type is also related to *Muelmat* because it is computed often as the exponent
of this term multiplied by a negative distance.
)--",
  };

  wsg_data["Muelmat"] = {
      .file = "rtepack.h",
      .desc = "A single Mueller 4x4 matrix.\n",
  };

  wsg_data["Stokvec"] = {
      .file = "rtepack.h",
      .desc = "A single Stokes vector (of length 4).\n",
  };

  wsg_data["PropmatVector"] = {
      .file = "rtepack.h",
      .desc = "A vector of *Propmat*.\n",
  };

  wsg_data["MuelmatVector"] = {
      .file = "rtepack.h",
      .desc = "A vector of *Muelmat*.\n",
  };

  wsg_data["StokvecVector"] = {
      .file = "rtepack.h",
      .desc = "A vector of *Stokvec*.\n",
  };

  wsg_data["PropmatMatrix"] = {
      .file = "rtepack.h",
      .desc = "A matrix of *Propmat*.\n",
  };

  wsg_data["MuelmatMatrix"] = {
      .file = "rtepack.h",
      .desc = "A matrix of *Muelmat*.\n",
  };

  wsg_data["MuelmatTensor3"] = {
      .file = "rtepack.h",
      .desc = "A *Tensor3* of *Muelmat*.\n",
  };

  wsg_data["StokvecMatrix"] = {
      .file = "rtepack.h",
      .desc = "A matrix of *Stokvec*.\n",
  };

  wsg_data["StokvecTensor3"] = {
      .file = "rtepack.h",
      .desc = "A *Tensor3* but of *Stokvec*.\n",
  };

  wsg_data["StokvecTensor4"] = {
      .file = "rtepack.h",
      .desc = "A *Tensor4* but of *Stokvec*.\n",
  };

  wsg_data["StokvecTensor5"] = {
      .file = "rtepack.h",
      .desc = "A *Tensor5* but of *Stokvec*.\n",
  };

  wsg_data["StokvecTensor6"] = {
      .file = "rtepack.h",
      .desc = "A *Tensor6* but of *Stokvec*.\n",
  };

  wsg_data["StokvecGriddedField6"] = {
      .file = "rtepack.h",
      .desc = R"--(A 6-dimensional gridof *Stokvec*.
)--",
  };

  wsg_data["ArrayOfPropmatVector"] = {
      .file = "rtepack.h",
      .desc = "A list of *PropmatVector*.\n",
  };

  wsg_data["ArrayOfMuelmatVector"] = {
      .file = "rtepack.h",
      .desc = "A list of *MuelmatVector*.\n",
  };

  wsg_data["ArrayOfStokvecVector"] = {
      .file = "rtepack.h",
      .desc = "A list of *StokvecVector*.\n",
  };

  wsg_data["ArrayOfPropmatMatrix"] = {
      .file = "rtepack.h",
      .desc = "A list of *PropmatMatrix*.\n",
  };

  wsg_data["ArrayOfMuelmatMatrix"] = {
      .file = "rtepack.h",
      .desc = "A list of *MuelmatMatrix*.\n",
  };

  wsg_data["ArrayOfMuelmatTensor3"] = {
      .file = "rtepack.h",
      .desc = "A list of *MuelmatTensor3*.\n",
  };

  wsg_data["ArrayOfStokvecMatrix"] = {
      .file = "rtepack.h",
      .desc = "A list of *StokvecMatrix*.\n",
  };

  wsg_data["ArrayOfStokvecTensor3"] = {
      .file = "rtepack.h",
      .desc = "A list of *StokvecTensor3*.\n",
  };

  wsg_data["ArrayOfArrayOfPropmatVector"] = {
      .file = "rtepack.h",
      .desc = "A list of *ArrayOfPropmatVector*.\n",
  };

  wsg_data["ArrayOfArrayOfMuelmatVector"] = {
      .file = "rtepack.h",
      .desc = "A list of *ArrayOfMuelmatVector*.\n",
  };

  wsg_data["ArrayOfArrayOfStokvecVector"] = {
      .file = "rtepack.h",
      .desc = "A list of *ArrayOfStokvecVector*.\n",
  };

  wsg_data["ArrayOfArrayOfPropmatMatrix"] = {
      .file = "rtepack.h",
      .desc = "A list of *ArrayOfPropmatMatrix*.\n",
  };

  wsg_data["ArrayOfArrayOfMuelmatMatrix"] = {
      .file = "rtepack.h",
      .desc = "A list of *ArrayOfMuelmatMatrix*.\n",
  };

  wsg_data["ArrayOfArrayOfStokvecMatrix"] = {
      .file = "rtepack.h",
      .desc = "A list of *ArrayOfStokvecMatrix*.\n",
  };

  wsg_data["NumericUnaryOperator"] = {
      .file = "operators.h",
      .desc = R"--(A simple functional type.

This type will work as a function pointer that takes a single *Numeric*
to produce another *Numeric*.
)--",
  };

  wsg_data["NumericBinaryOperator"] = {
      .file = "operators.h",
      .desc = R"--(A simple functional type.

This type will work as a function pointer that takes two *Numeric*
to produce another *Numeric*.
)--",
  };

  wsg_data["NumericTernaryOperator"] = {
      .file = "operators.h",
      .desc = R"--(A simple functional type.

This type will work as a function pointer that takes three *Numeric*
to produce a single *Numeric*.
)--",
  };

  wsg_data["JacobianTargets"] = {
      .file = "jacobian.h",
      .desc = R"--(A list of targets for use in Jacobian Matrix calculations

This type flags the type of calculations that should be performed
when computing the Jacobian matrix or partial derivatives.
)--",
  };

  wsg_data["JacobianTargetType"] = {
      .file = "jacobian.h",
      .desc = R"--(A type of target for use in Jacobian Matrix calculations
)--",
  };

  wsg_data["JacobianTargetsDiagonalCovarianceMatrixMap"] = {
      .file = "retrieval_target.h",
      .desc =
          R"--(A map target types to matrix and inverse matrix pairs of *BlockMatrix*

The intended use of this type is to store required *BlockMatrix* objects so that
the user-interface for setting up retrieval targets can be simplified.
)--",
  };

  wsg_data["PropagationPathPoint"] = {
      .file = "path_point.h",
      .desc = R"--(A simple path-point of a propagation path

This point describes the origin and line-of-sight of the tracked
radiation.
)--",
  };

  wsg_data["ArrayOfPropagationPathPoint"] = {
      .file = "path_point.h",
      .desc = "A list of *PropagationPathPoint*.\n",
  };

  wsg_data["ArrayOfArrayOfPropagationPathPoint"] = {
      .file = "path_point.h",
      .desc = "A list of *ArrayOfPropagationPathPoint*.\n",
  };

  wsg_data["ArrayOfArrayOfArrayOfPropagationPathPoint"] = {
      .file = "path_point.h",
      .desc = "A list of *ArrayOfArrayOfPropagationPathPoint*.\n",
  };

  wsg_data["Vector3"] = {
      .file = "matpack.h",
      .desc = "A fixed-size 3D version of *Vector*.\n",
  };

  wsg_data["Vector2"] = {
      .file = "matpack.h",
      .desc = "A fixed-size 2D version of *Vector*.\n",
  };

  wsg_data["ArrayOfVector3"] = {
      .file = "matpack.h",
      .desc = "A list of *Vector3*\n",
  };

  wsg_data["ArrayOfVector2"] = {
      .file = "matpack.h",
      .desc = "A list of *Vector2*\n",
  };

  wsg_data["DescendingGrid"] = {
      .file = "matpack.h",
      .desc = "A sorted grid of always Descending values.\n",
  };

  wsg_data["AscendingGrid"] = {
      .file = "matpack.h",
      .desc = "A sorted grid of always ascending values.\n",
  };

  wsg_data["ArrayOfAscendingGrid"] = {
      .file = "matpack.h",
      .desc = "A list of *AscendingGrid*.\n",
  };

  wsg_data["SpectralRadianceOperator"] = {
      .file = "fwd.h",
      .desc = R"--(An operator for getting the *spectral_radiance*

An object of this type can be called with a frequency, position and
line-of-sight to get the corresponding spectral radiance.
)--",
  };

  wsg_data["SpeciesIsotope"] = {
      .file = "isotopologues.h",
      .desc = "Contains descriptions about an isotope.\n",
  };

  wsg_data["ArrayOfSpeciesIsotope"] = {
      .file = "isotopologues.h",
      .desc = "List of *SpeciesIsotope*.\n",
  };

  wsg_data["SensorPosLos"] = {
      .file = "obsel.h",
      .desc = "A position and line-of-sight of a sensor.\n",
  };

  wsg_data["SensorPosLosVector"] = {
      .file = "obsel.h",
      .desc = "Vector of *SensorPosLos*.\n",
  };

  wsg_data["SensorObsel"] = {
      .file = "obsel.h",
      .desc = R"(A single observation element.

Expected use of this type is to generate the measurement vector
of a sensor, where this observation element represent the readout
from that sensor in a convenient unit (commonly Kelvin or 
W sr :math:`^{-1}m :math:`^{-2}Hz :math:`^{-1}`, but not exclusively)

It deals with averaging the frequency grid sampled by a sensor element
and the transmission of the sensor system onto the sampling device, as
well as the sampling device's polarization response.
)",
  };

  wsg_data["ArrayOfSensorObsel"] = {
      .file = "obsel.h",
      .desc = "List of *SensorObsel*.\n",
  };

  wsg_data["DisortSettings"] = {
      .file = "disort.h",
      .desc = "The settings required to run Disort.\n",
  };

  wsg_data["AbsorptionLookupTable"] = {
      .file = "lookup_map.h",
      .desc = "A table of lookup calculations.\n",
  };

  wsg_data["AbsorptionLookupTables"] = {
      .file = "lookup_map.h",
      .desc = "A map of tables of of lookup calculations.\n",
      .map_type = true,
  };

  return wsg_data;
}

const std::unordered_map<std::string, WorkspaceGroupRecord>&
internal_workspace_groups() {
  static const auto wsg_data = internal_workspace_groups_creator();
  return wsg_data;
}
