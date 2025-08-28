#include "workspace_group_friends.h"

#include <arts_options.h>

#include <algorithm>
#include <stdexcept>

#include "workspace_groups.h"

namespace {
std::unordered_map<std::string, WorkspaceGroupRecord> group_friends_internal() {
  std::unordered_map<std::string, WorkspaceGroupRecord> wsg_data;

  wsg_data["AbsorptionBand"] = {
      .file = "lbl.h",
      .desc =
          R"(Contains all information about bands of related absorption lines.

This information includes

#. A list of :class:`~pyarts3.arts.AbsorptionLine`.

#. The line shape profile model.  See *LineByLineLineshape* for available line shape profiles.

#. The frequency cutoff value in [Hz] and type.  See *LineByLineCutoffType* for available cutoff types.
)",
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

  wsg_data["XsecRecord"] = {
      .file = "xsec_fit.h",
      .desc = R"(A single cross-section record

These cross-section records contains information about the valid temperature and
pressure ranges as well as well as the fitting coefficients used to compute
and interpolate the cross-section to other temperatures and pressures
)",
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

  wsg_data["SurfacePoint"] = {
      .file = "surf.h",
      .desc =
          R"--(A surface point.

This keeps two things:

#. The local normal vector.

#. A map of the same keys as *SurfaceField* but towards *Numeric* data.
   It is similar to how *AtmPoint* is towards *AtmField*.
)--",
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

  wsg_data["Specmat"] = {
      .file = "rtepack.h",
      .desc = "A single Complex Mueller 4x4 matrix.\n",
  };

  wsg_data["SpeciesTag"] = {
      .file = "species_tags.h",
      .desc = R"(A tagged absorption species

These tags are used to help ARTS identify the species
so that reading routines can find the correct data files.
)",
  };

  wsg_data["Sparse"] = {
      .file = "matpack_sparse.h",
      .desc = "A sparse version of *Matrix*\n",
  };

  wsg_data["ScatteringMetaData"] = {
      .file = "optproperties.h",
      .desc = "Holds meta data about the scattering\n",
  };

  wsg_data["SingleScatteringData"] = {
      .file = "optproperties.h",
      .desc = "Holds single scattering data\n",
  };

  wsg_data["SensorPosLos"] = {
      .file = "obsel.h",
      .desc = "A position and line-of-sight of a sensor.\n",
  };

  wsg_data["SensorPosLosVector"] = {
      .file = "obsel.h",
      .desc = "Vector of *SensorPosLos*.\n",
  };

  wsg_data["Rational"] = {
      .file = "matpack.h",
      .desc = "Holds a rational number as two *Index* n / d\n",
  };

  wsg_data["QuantumIdentifierNumericMap"] = {
      .file = "lbl.h",
      .desc =
          R"--(A map from *QuantumIdentifier* to *Numeric*.
)--",
      .map_type = true};

  wsg_data["PredefinedModelDataVariant"] = {
      .file = "predef.h",
      .desc =
          R"--(One of the following:

#. :class:`~pyarts3.arts.predef.PredefinedModelDataName`
#. :class:`~pyarts3.arts.predef.PredefinedModelDataWaterDataMTCKD4`
)--",
  };

  wsg_data["PairOfBlockMatrix"] = {
      .file = "retrieval_target.h",
      .desc = R"--(A pair of *BlockMatrix* objects)--",
  };

  wsg_data["NumericBinaryOperator"] = {
      .file = "operators.h",
      .desc = R"--(A simple functional type.

This type will work as a function pointer that takes two *Numeric*
to produce another *Numeric*.

.. math::

    m = f(x, y)

)--",
  };

  wsg_data["MuelmatMatrix"] = {
      .file = "rtepack.h",
      .desc = "A matrix of *Muelmat*.\n",
  };

  wsg_data["MatrixOfDisortBDRF"] = {
      .file = "disort.h",
      .desc = "A 2 dimensional array of *DisortBDRF*\n",
  };

  wsg_data["LinemixingSingleEcsData"] = {
      .file = "lbl.h",
      .desc =
          R"--(Data for single species ECS.
)--",
  };

  wsg_data["LinemixingSpeciesEcsData"] = {
      .file = "lbl.h",
      .desc =
          R"--(A map from *SpeciesEnum* to *LinemixingSingleEcsData*
)--",
      .map_type = true,
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

  wsg_data["DisortBDRF"] = {
      .file = "disort.h",
      .desc = "A bidirectional reflectance function\n",
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

  wsg_data["ArrayOfSubsurfacePoint"] = {
      .file = "subsurface.h",
      .desc =
          R"--(A list of *SubsurfacePoint*
)--",
  };

  wsg_data["ArrayOfStokvecTensor3"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *StokvecTensor3*.\n",
      .array_depth = 1,
  };

  wsg_data["ArrayOfSpeciesTag"] = {
      .file = "species_tags.h",
      .desc = R"--(A list of *SpeciesTag*
)--",
  };

  wsg_data["ArrayOfSpeciesIsotope"] = {
      .file = "isotopologues.h",
      .desc = "List of *SpeciesIsotope*.\n",
  };

  wsg_data["ArrayOfSparse"] = {
      .file = "matpack_sparse.h",
      .desc = "A list of *Sparse*\n",
  };

  wsg_data["ArrayOfNamedGriddedField2"] = {
      .file = "matpack.h",
      .desc = "A list of *NamedGriddedField2*\n",
  };

  wsg_data["ArrayOfMuelmatMatrix"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *MuelmatMatrix*.\n",
      .array_depth = 1,
  };

  wsg_data["ArrayOfGriddedField3"] = {
      .file = "matpack.h",
      .desc = "A list of *GriddedField3*\n",
  };

  wsg_data["ArrayOfGriddedField4"] = {
      .file = "matpack.h",
      .desc = "A list of *GriddedField4*\n",
  };

  wsg_data["ArrayOfGriddedField2"] = {
      .file = "matpack.h",
      .desc = "A list of *GriddedField2*\n",
  };

  wsg_data["ArrayOfGriddedField1"] = {
      .file = "matpack.h",
      .desc = "A list of *GriddedField1*\n",
  };

  wsg_data["AbsorptionLookupTable"] = {
      .file = "lookup_map.h",
      .desc = R"(A table of lookup calculations.

Effectively holds a *Tensor4* of pre-computed cross-section data.
This table is used to interpole to a pressure, temperature, water vmr, and frequency grid.
)",
  };

  wsg_data["ArrayOfArrayOfMuelmatMatrix"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *ArrayOfMuelmatMatrix*.\n",
      .array_depth = 2,
  };

  wsg_data["ArrayOfArrayOfMuelmatVector"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *ArrayOfMuelmatVector*.\n",
      .array_depth = 2,
  };

  wsg_data["ArrayOfArrayOfVector"] = {
      .file = "matpack.h",
      .desc = "A list of *ArrayOfVector*\n",
  };

  wsg_data["ArrayOfGriddedField1Named"] = {
      .file = "matpack.h",
      .desc = "A list of *GriddedField1Named*\n",
  };

  wsg_data["ArrayOfMatrix"] = {
      .file = "matpack.h",
      .desc = "A list of *Matrix*\n",
  };

  wsg_data["ArrayOfQuantumIdentifier"] = {
      .file = "quantum.h",
      .desc = "A list of *QuantumIdentifier*\n",
  };

  wsg_data["ArrayOfScatteringMetaData"] = {
      .file = "optproperties.h",
      .desc = "A list of *ScatteringMetaData*\n",
  };

  wsg_data["ArrayOfSingleScatteringData"] = {
      .file = "optproperties.h",
      .desc = "A list of *SingleScatteringData*\n",
  };

  wsg_data["ArrayOfTime"] = {
      .file = "matpack.h",
      .desc = "A list of *Time*\n",
  };

  wsg_data["ArrayOfVector"] = {
      .file = "matpack.h",
      .desc = "A list of *Vector*\n",
  };

  wsg_data["AtmData"] = {
      .file = "atm.h",
      .desc = R"--(An atmospheric data field.

Each atmospheric data field can be mapped to a single altitude-latitude-longitude coordinate,
producing a *Numeric* value at that point in the atmosphere.

It holds essentially two things:

#. Rules for how to extrapolate the data in the six directions of the atmosphere (up, down, north, south, east, west).

#. A data field.  This data field is one of the following types:

    #. *Numeric* - The data field is constant in the atmosphere.
       Cannot consider the extrapolation rules as there is no grid.

    #. *SortedGriddedField3* - The grids are altitude, latitude, longitude.
       Will consider the extrapolation rules but otherwise performs linear interpolation between all points.
       With the additional rule that longitude is considerd cyclic around [-180, 180).

    #. *NumericTernaryOperator* - The data field has functional form.
       Cannot consider the extrapolation rules as there is no grid.
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

  wsg_data["JacobianTargetType"] = {
      .file = "jacobian.h",
      .desc = R"--(A type of target for use in Jacobian Matrix calculations
)--",
  };

  wsg_data["SubsurfacePoint"] = {
      .file = "subsurface.h",
      .desc =

          R"--(A surface point.

This keeps a map of the same keys as *SubsurfaceField* but towards *Numeric* data.
It is similar to how *AtmPoint* is towards *AtmField*.
)--",
  };

  wsg_data["SensorObsel"] = {
      .file = "obsel.h",
      .desc = R"(A single observation element.

This should result in a single element in a *measurement_vector*.

Expected use of this type is to generate the measurement vector
of a sensor, where this observation element represent the readout
from that sensor in a convenient unit (commonly Kelvin or 
W sr :math:`^{-1}` m :math:`^{-2}` Hz :math:`^{-1}`, but not exclusively)

It deals with averaging the frequency grid sampled by a sensor element
and the transmission of the sensor system onto the sampling device, as
well as the sampling device's polarization response.

.. note::
   Multiple *SensorObsel* can be used to represent a single sensor
   with multiple channels, such as a sensor with multiple detectors
   or a sensor with multiple frequency channels.  They then can conveniently
   share much of their grids, but have different weighting functions
   and/or different sampling devices.
)",
  };

  wsg_data["Propmat"] = {
      .file = "rtepack.h",
      .desc = R"--(A single propagation matrix.

Due to the properties of a propagation matrix, only 7 independents need be stored.
The propagation matrix is thus represented as:

.. math::
    \mathbf{K} = \left[ \begin {array} {rrrr}
    A & B & C & D \\
    B & A & U & V \\
    C &-U & A & W \\
    D &-V &-W & A
    \end {array} \right]

This type is related to *Stokvec* in that its first 4 elements are the same as
the first 4 elements of *Stokvec* for pure clearsky radiative transfers.

This type is also related to *Muelmat* because it is computed often as the exponent
of this term multiplied by a negative distance.
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

  wsg_data["MuelmatTensor3"] = {
      .file = "rtepack.h",
      .desc = "A *Tensor3* of *Muelmat*.\n",
  };

  wsg_data["Muelmat"] = {
      .file = "rtepack.h",
      .desc = "A single Mueller 4x4 matrix.\n",
  };

  wsg_data["ArrayOfVector2"] = {
      .file = "matpack.h",
      .desc = "A list of *Vector2*\n",
  };

  wsg_data["GriddedField1"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 1 dimensional gridded set of *Numeric* data

The grid is 1 *Vector*

Both the data and the grid may be named)--",
  };

  wsg_data["GriddedField1Named"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 3 dimensional gridded set of *Numeric* data

The grid is a combination of 1 *Vector* and 1 *ArrayOfString*

Both the data and the grid may be named
)--",
  };

  wsg_data["GeodeticField3"] = {
      .file = "rtepack.h",
      .desc = R"--(A 3-dimensional gridof *Numeric*.

The grids are *altitude_grid* x *latitude_grid* x *longitude_grid*.
The types are *AscendingGrid* x *LatGrid* x *LonGrid*.
)--",
  };

  for (auto& g : internal_options()) {
    if (wsg_data.find(g.name) != wsg_data.end())
      throw std::runtime_error(
          "Duplicate workspace group name (name is reserved as options-group): " +
          g.name);

    wsg_data[g.name] = {
        .file = "enums.h",
        .desc = g.docs(),
    };
  }

  for (const auto& [name, g] : internal_workspace_groups()) {
    if (auto ptr = wsg_data.find(name); ptr != wsg_data.end())
      wsg_data.erase(ptr);
  }
  return wsg_data;
}
}  // namespace

const std::unordered_map<std::string, WorkspaceGroupRecord>&
workspace_group_friends() {
  static const auto friends = group_friends_internal();
  return friends;
}
