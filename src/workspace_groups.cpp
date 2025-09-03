#include "workspace_groups.h"

#include <arts_options.h>

#include <algorithm>
#include <stdexcept>

#include "workspace_agendas.h"

void add_select_options(
    std::unordered_map<std::string, WorkspaceGroupRecord>& wsg_data,
    const std::vector<std::string>& select_options) {
  const auto& options = internal_options();

  for (const auto& opt : select_options) {
    auto it = std::ranges::find_if(
        options, [&opt](const auto& o) { return o.name == opt; });

    if (it == options.end())
      throw std::runtime_error("Option " + opt +
                               " not found in internal options.");

    wsg_data[opt] = {
        .file = "enums.h",
        .desc = it->docs(),
    };
  }
}

void agenda_operators(
    std::unordered_map<std::string, WorkspaceGroupRecord>& wsg_data) {
  for (auto& [name, ag] : internal_workspace_agendas()) {
    wsg_data[name + "Operator"] = {
        .file = "auto_agenda_operators.h",
        .desc =
            ag.desc +
            "\n\nThis is the free-form customization point of the agenda *" +
            name +
            "*\n\n.. note::\n    Output constraints must be followed\n\n",
    };
  }
}

std::unordered_map<std::string, WorkspaceGroupRecord>
internal_workspace_groups_creator() {
  std::unordered_map<std::string, WorkspaceGroupRecord> wsg_data;

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
      .desc = R"(Describes a set of function calls and variable definitions

Agendas are effectively constrained callback methods.
They define a set of inputs and outputs.  These variables
are always local to a method that uses the agenda.
In addition to these local variables, the agenda has full access to copy or share any workspace variables.
However, by tracking how named workspace variables are used in an agenda, ARTS ensures that the agenda
does not change the global workspace while minimizing the number of variables that need to be copied or shared.
)",
  };

  wsg_data["Any"] = {
      .file = "supergeneric.h",
      .desc =
          "Meta type for any workspace group (see :doc:`workspace.groups`)\n",
  };

  wsg_data["ArrayOfArrayOfSpeciesTag"] = {
      .file = "species_tags.h",
      .desc = "A list of *ArrayOfSpeciesTag*\n",
  };

  wsg_data["ArrayOfAtmPoint"] = {
      .file = "atm.h",
      .desc = "A list of *AtmPoint*\n",
  };

  wsg_data["ArrayOfCIARecord"] = {
      .file = "cia.h",
      .desc = R"(A list of *CIARecord*
)",
  };

  wsg_data["ArrayOfIndex"] = {
      .file = "matpack.h",
      .desc = "A list of *Index*\n",
  };

  wsg_data["ArrayOfQuantumLevelIdentifier"] = {
      .file = "quantum.h",
      .desc = "A list of *QuantumLevelIdentifier*\n",
  };

  wsg_data["ArrayOfScatteringSpecies"] = {
      .file = "scattering/scattering_species.h",
      .desc = "Represents species of scattering paritlces in the atmosphere.",
  };

  wsg_data["ScatteringSpeciesProperty"] = {
      .file = "scattering/properties.h",
      .desc = "Meta data for scattering spefcies.",
  };

  wsg_data["SurfacePropertyTag"] = {
      .file = "surf.h",
      .desc = R"--(A surface property.

These tags are part of the keys that can be used to access a *SurfaceField*.
They are completely free-form and currently not used by ARTS internally.
Instead, they offer a customization point for users to define their own
surface properties and ensures we can access them in a consistent way.
)--",
  };

  wsg_data["ArrayOfSpeciesEnum"] = {
      .file = "species.h",
      .desc = R"--(A list of *SpeciesEnum*
)--",
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

  wsg_data["ArrayOfXsecRecord"] = {
      .file = "xsec_fit.h",
      .desc =
          R"--(A list of *XsecRecord*
)--",
  };

  wsg_data["AtmField"] = {
      .file = "atm.h",
      .desc = R"--(An atmospheric field.

An atmospheric field holds two things:

#. The top of the atmosphere altitude, which is the altitude at which the atmosphere ends.  Unit: m

#. A map of *AtmData*.  The available types of keys are:

   #. *AtmKey*

   #. *SpeciesEnum*

   #. *SpeciesIsotope*

   #. *QuantumIdentifier*

   #. *ScatteringSpeciesProperty*

   See each key for more information on what type of data it allows holding.
)--",
  };

  wsg_data["AtmPoint"] = {
      .file = "atm.h",
      .desc = R"--(An atmospheric point.

This can be thought of as the sampling all the *AtmData* of and *AtmField*
at a single altitude-latitude-longitude coordinate.
It, like *AtmField* also acts like a map.  They keys are the same as for *AtmField*. However,
the values are simply the *Numeric* data at that point in the atmosphere.
)--",
  };

  wsg_data["CIARecord"] = {
      .file = "cia.h",
      .desc =
          R"--(Contains information to compute collision induced absorption (CIA) for a pair of species.

A record holds a list of *GriddedField2* objects,
each which describe a separate band of absorption with dimensions of
frequency times temperature.  The *Matrix* objects in the *GriddedField2*
simply holds data to be interpolated in the frequency and temperature
with physical units of [m :math:`^5` per molecule :math:`^2`].
)--",
  };

  wsg_data["CallbackOperator"] = {
      .file = "callback.h",
      .desc = R"(Used to inject custom code into *Agenda*

This completely breaks the type system of ARTS and should only be used.
You are on your own when things go wrong with this.
)",
  };

  wsg_data["BlockMatrix"] = {
      .file = "covariance_matrix.h",
      .desc =
          R"(The data for a single :class:`~pyarts3.arts.Block`, likely part of a *CovarianceMatrix*.

This holds either a *Matrix* or a *Sparse* matrix.
)",
  };

  wsg_data["CovarianceMatrix"] = {
      .file = "covariance_matrix.h",
      .desc =
          R"(A covariance matrix is a square matrix that describes the covariance of some property.  

Please see the different workspace variables of this type for more information.

In ARTS, this square matrix is represented by two lists of :class:`~pyarts3.arts.Block`.
These are used to give both the covariance matrix and the inverse covariance matrix.
The block-structure allows for efficient storage and computation of the covariance matrix.
)",
  };

  wsg_data["GriddedField2"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 2 dimensional gridded set of *Numeric* data

The grid is a combination of 2 *Vector*

Both the data and the grid may be named
)--",
  };

  wsg_data["QuantumIdentifierVectorMap"] = {
      .file = "lbl.h",
      .desc =
          R"--(A map from *QuantumIdentifier* to *Vector*.
)--",
      .map_type = true};

  wsg_data["QuantumIdentifierGriddedField1Map"] = {
      .file = "lbl.h",
      .desc =
          R"--(A map from *QuantumIdentifier* to *GriddedField1*.
)--",
      .map_type = true};

  wsg_data["Index"] = {
      .file       = "matpack.h",
      .desc       = "A 64 bit signed integer type\n",
      .value_type = true,
  };

  wsg_data["LinemixingEcsData"] = {
      .file = "lbl.h",
      .desc =
          R"--(A map from *SpeciesIsotope* to *LinemixingSpeciesEcsData*
)--",
      .map_type = true,
  };

  wsg_data["Matrix"] = {
      .file = "matpack.h",
      .desc = "A 2 dimensional array of *Numeric*\n",
  };

  wsg_data["Numeric"] = {
      .file       = "matpack.h",
      .desc       = "IEEE 754 binary64 floating point number\n",
      .value_type = true,
  };

  wsg_data["PredefinedModelData"] = {
      .file = "predef.h",
      .desc =
          R"--(A map of *SpeciesIsotope* to *PredefinedModelDataVariant*
)--",
      .map_type = true,
  };

  wsg_data["QuantumIdentifier"] = {
      .file = "quantum.h",
      .desc =
          R"--(An ID for an absorption species state

It contains upper and lower level information of a quantum state.

It can identify:

1. a species
2. an isotopologue of a species
3. an absorption band of an isotopologue
4. an absorption line of an isotopologue
)--",
  };

  wsg_data["QuantumLevelIdentifier"] = {
      .file = "quantum.h",
      .desc =
          R"--(An ID for an absorption species state

It contains the level information of a quantum state.

It can identify:

1. a species
2. an isotopologue of a species
3. the energy level of absorption band(s) of an isotopologue
4. the energy level of absorption line(s) of an isotopologue
)--",
  };

  wsg_data["String"] = {
      .file       = "mystring.h",
      .desc       = "Basic string type\n",
      .value_type = true,
  };

  wsg_data["SurfaceField"] = {
      .file = "surf.h",
      .desc =
          R"--(A surface field.

A surface field effectively holds two things:

#. A *Vector2* of the ellipsoid.  a and b parameters.  Unit: m

#. A map of :class:`~pyarts3.arts.SurfaceData`.  The available types of keys are:

   #. *SurfaceKey*

   #. *SurfacePropertyTag*

   See each key for more information on what type of data it allows holding.
)--",
  };

  wsg_data["SubsurfaceField"] = {
      .file = "subsurface.h",
      .desc =
          R"--(A sub-surface field.

A sub-surface field effectively holds two things:

#. A *Numeric* of the deepest depth of the subsurface.  Unit: m

#. A map of :class:`~pyarts3.arts.SubsurfaceData`.  The available types of keys are:

   #. *SubsurfaceKey*

   See each key for more information on what type of data it allows holding.
)--",
  };

  wsg_data["Time"] = {
      .file = "artstime.h",
      .desc = R"(Represents a time stamp
)",
  };

  wsg_data["Vector"] = {
      .file = "matpack.h",
      .desc = "A 1 dimensional array of *Numeric*\n",
  };

  wsg_data["Stokvec"] = {
      .file = "rtepack.h",
      .desc = R"(A single Stokes vector (of length 4).

.. math::
    \vec{I} = \left[ \begin {array} {r}
      I \\ Q \\ U \\ V
    \end {array} \right]

The Stokes vector is used to represent the state of polarization of light.

The components of the Stokes vector are:

#. :math:`I` - the total intensity of the light

#. :math:`Q` - the difference in intensity between horizontally and vertically polarized light

#. :math:`U` - the difference in intensity between light polarized at +45 degrees and -45 degrees

#. :math:`V` - the difference in intensity between right and left circularly polarized light
)",
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

  wsg_data["SpecmatMatrix"] = {
      .file = "rtepack.h",
      .desc = "A matrix of *Muelmat*.\n",
  };

  wsg_data["StokvecMatrix"] = {
      .file = "rtepack.h",
      .desc = "A matrix of *Stokvec*.\n",
  };

  wsg_data["StokvecSortedGriddedField1"] = {
      .file = "rtepack.h",
      .desc = R"--(A 1-dimensional grid of *Stokvec*.

The grids are 1 *AscendingGrid*.
)--",
  };

  wsg_data["StokvecSortedGriddedField2"] = {
      .file = "rtepack.h",
      .desc = R"--(A 2-dimensional grid of *Stokvec*.

The grids are 2 *AscendingGrid*.
)--",
  };

  wsg_data["StokvecSortedGriddedField3"] = {
      .file = "rtepack.h",
      .desc = R"--(A 3-dimensional grid of *Stokvec*.

The grids are 3 *AscendingGrid*.
)--",
  };

  wsg_data["StokvecSortedGriddedField4"] = {
      .file = "rtepack.h",
      .desc = R"--(A 4-dimensional grid of *Stokvec*.

The grids are 4 *AscendingGrid*.
)--",
  };

  wsg_data["StokvecSortedGriddedField5"] = {
      .file = "rtepack.h",
      .desc = R"--(A 5-dimensional grid of *Stokvec*.

The grids are 5 *AscendingGrid*.
)--",
  };

  wsg_data["StokvecSortedGriddedField6"] = {
      .file = "rtepack.h",
      .desc = R"--(A 6-dimensional gridof *Stokvec*.

The grids are 6 *AscendingGrid*.
)--",
  };

  wsg_data["GriddedSpectralField6"] = {
      .file = "rtepack.h",
      .desc = R"--(A 6-dimensional gridof *Stokvec*.

The grids are altitude x latitude x longitude x zenith x azimuth x frequency of types
*AscendingGrid* x *LatGrid* x *LonGrid* x *ZenithGrid* x *AzimuthGrid* x *AscendingGrid* x.
)--",
  };

  wsg_data["ArrayOfPropmatVector"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *PropmatVector*.\n",
      .array_depth = 1,
  };

  wsg_data["ArrayOfMuelmatVector"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *MuelmatVector*.\n",
      .array_depth = 1,
  };

  wsg_data["ArrayOfStokvecVector"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *StokvecVector*.\n",
      .array_depth = 1,
  };

  wsg_data["ArrayOfPropmatMatrix"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *PropmatMatrix*.\n",
      .array_depth = 1,
  };

  wsg_data["ArrayOfSpecmatMatrix"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *SpecmatMatrix*.\n",
      .array_depth = 1,
  };

  wsg_data["ArrayOfMuelmatTensor3"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *MuelmatTensor3*.\n",
      .array_depth = 1,
  };

  wsg_data["ArrayOfStokvecMatrix"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *StokvecMatrix*.\n",
      .array_depth = 1,
  };

  wsg_data["ArrayOfArrayOfPropmatVector"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *ArrayOfPropmatVector*.\n",
      .array_depth = 2,
  };

  wsg_data["ArrayOfArrayOfStokvecVector"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *ArrayOfStokvecVector*.\n",
      .array_depth = 2,
  };

  wsg_data["ArrayOfArrayOfPropmatMatrix"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *ArrayOfPropmatMatrix*.\n",
      .array_depth = 2,
  };

  wsg_data["ArrayOfArrayOfStokvecMatrix"] = {
      .file        = "rtepack.h",
      .desc        = "A list of *ArrayOfStokvecMatrix*.\n",
      .array_depth = 2,
  };

  wsg_data["NumericUnaryOperator"] = {
      .file = "operators.h",
      .desc = R"--(A simple functional type.

This type will work as a function pointer that takes a single *Numeric*
to produce another *Numeric*.

.. math::

    m = f(x)

)--",
  };

  wsg_data["NumericTernaryOperator"] = {
      .file = "operators.h",
      .desc = R"--(A simple functional type.

This type will work as a function pointer that takes three *Numeric*
to produce a single *Numeric*.

.. math::

    m = f(x, y, z)

)--",
  };

  wsg_data["JacobianTargets"] = {
      .file = "jacobian.h",
      .desc = R"--(A list of targets for use in Jacobian Matrix calculations

This type flags the type of calculations that should be performed
when computing the Jacobian matrix or partial derivatives.
)--",
  };

  wsg_data["JacobianTargetsDiagonalCovarianceMatrixMap"] = {
      .file = "retrieval_target.h",
      .desc =
          R"--(A map target types to matrix and inverse matrix pairs of *BlockMatrix*

The intended use of this type is to store required *BlockMatrix* objects so that
the user-interface for setting up retrieval targets can be simplified.
)--",
      .map_type = true,
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

  wsg_data["DescendingGrid"] = {
      .file = "matpack.h",
      .desc = "A sorted *Vector* of always descending values.\n",
  };

  wsg_data["AscendingGrid"] = {
      .file = "matpack.h",
      .desc = "A sorted *Vector* of always ascending values.\n",
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
      .desc = "Contains name and data about an isotope.\n",
  };

  wsg_data["ArrayOfSensorObsel"] = {
      .file = "obsel.h",
      .desc = "List of *SensorObsel*.\n",
  };

  wsg_data["DisortFlux"] = {
      .file = "disort.h",
      .desc = R"(The flux result variable for Disort.

#. *AscendingGrid* frequency grid
#. *DescendingGrid* level altitude grid
#. *Matrix* upwelling flux
#. *Matrix* diffuse downwelling flux
#. *Matrix* direct downwelling flux
)",
  };

  wsg_data["DisortRadiance"] = {
      .file = "disort.h",
      .desc = R"(The radiance result variable for Disort.

#. *AscendingGrid* frequency grid
#. *DescendingGrid* level altitude grid
#. *AzimuthGrid* azimuth grid
#. *ZenithGrid* zenith grid
#. *Tensor4* radiance data
)",
  };

  wsg_data["DisortSettings"] = {
      .file = "disort.h",
      .desc = R"(The settings required to run Disort.

#. *Index* Quadrature dimension
#. *Index* Legendre order
#. *Index* Fourier order
#. *Index* The frequency grid
#. *Index* The level altitude grid
#. *Vector* Solar azimuth angles
#. *Vector* Solar zenith angles
#. *Vector* Solar source
#. *MatrixOfDisortBDRF* Bi-directional reflectance distribution functions
#. *Matrix* Optical thicknesses
#. *Matrix* Single scattering albedos
#. *Matrix* Fractional scattering
#. *Tensor3* Source function polynomial
#. *Tensor3* Legendre polynomial coefficients
#. *Tensor3* Positive boundary condition
#. *Tensor3* Negative boundary condition
)",
  };

  wsg_data["AbsorptionLookupTables"] = {
      .file     = "lookup_map.h",
      .desc     = "A map of *SpeciesEnum* to *AbsorptionLookupTable*.\n",
      .map_type = true,
  };

  wsg_data["ZenithGriddedField1"] = {
      .file = "matpack.h",
      .desc = R"--(A 1-dimensional gridof *Numeric*.

The grids are 1 *ZenithGrid*.
)--",
  };

  wsg_data["SortedGriddedField1"] = {
      .file = "matpack.h",
      .desc = R"--(A 1-dimensional gridof *Numeric*.

The grids are 1 *AscendingGrid*.
)--",
  };

  wsg_data["SortedGriddedField2"] = {
      .file = "matpack.h",
      .desc = R"--(A 2-dimensional gridof *Numeric*.

The grids are 2 *AscendingGrid*.
)--",
  };

  wsg_data["GeodeticField2"] = {
      .file = "matpack.h",
      .desc = R"--(A 2-dimensional gridof *Numeric*.

The grids are *latitude_grid* x *longitude_grid*.
The types are *LatGrid* x *LonGrid*.
)--",
  };

  wsg_data["SortedGriddedField3"] = {
      .file = "matpack.h",
      .desc = R"--(A 3-dimensional gridof *Numeric*.

The grids are 3 *AscendingGrid*.
)--",
  };

  wsg_data["SortedGriddedField4"] = {
      .file = "matpack.h",
      .desc = R"--(A 4-dimensional gridof *Numeric*.

The grids are 4 *AscendingGrid*.
)--",
  };

  wsg_data["SortedGriddedField5"] = {
      .file = "matpack.h",
      .desc = R"--(A 5-dimensional gridof *Numeric*.

The grids are 5 *AscendingGrid*.
)--",
  };

  wsg_data["SortedGriddedField6"] = {
      .file = "matpack.h",
      .desc = R"--(A 6-dimensional gridof *Numeric*.

The grids are 6 *AscendingGrid*.
)--",
  };

  wsg_data["ZenithGrid"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 1-dimensional vector of *Numeric* that are guaranteed to be within the range [0, 180].
)--",
  };

  wsg_data["AzimuthGrid"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 1-dimensional vector of *Numeric* that are guaranteed to be within the range [0, 360).
)--",
  };

  wsg_data["LonGrid"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 1-dimensional vector of *Numeric* that are guaranteed to be within the range [-180, 180).
)--",
  };

  wsg_data["LatGrid"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 1-dimensional vector of *Numeric* that are guaranteed to be within the range [-90, 90].
)--",
  };

  wsg_data["SpectralRadianceTransformOperator"] = {
      .file = "spectral_radiance_transform_operator.h",
      .desc =
          R"--(Transformation of *spectral_radiance* and *spectral_radiance_jacobian*

This type of transformation should be used limitedly.  It is useful
as the last step before creating a *measurement_vector*, as it is
used for in *measurement_vectorFromSensor* or just before displaying
data in a plotting routine.  It will destroy the *spectral_radiance*
and *spectral_radiance_jacobian* and replace them with the transformed
values.  They can likely not be reused for further calculations.

Parameters
----------
spectral_radiance : StokvecVector
    As WSV *spectral_radiance* **[INOUT]**
spectral_radiance_jacobian : StokvecMatrix
    As WSV *spectral_radiance_jacobian* **[INOUT]**
frequency_grid : AscendingGrid
    As WSV *frequency_grid* **[IN]**
ray_path_point : PropagationPathPoint
    As WSV *ray_path_point* **[IN]**
)--",
  };

  add_select_options(wsg_data,
                     {
                         "InterpolationExtrapolation",
                         "SpeciesEnum",
                         "AtmKey",
                         "SurfaceKey",
                         "SubsurfaceKey",
                         "SpectralRadianceUnitType",
                     });

  agenda_operators(wsg_data);

  return wsg_data;
}

const std::unordered_map<std::string, WorkspaceGroupRecord>&
internal_workspace_groups() {
  static const auto wsg_data = internal_workspace_groups_creator();
  return wsg_data;
}
