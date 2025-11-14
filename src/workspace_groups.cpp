#include "workspace_groups.h"

#include <arts_options.h>

#include <algorithm>
#include <stdexcept>
#include <string>

#include "workspace_agendas.h"

void add_arrays_of(
    std::unordered_map<std::string, WorkspaceGroupRecord>& wsg_data,
    const std::vector<std::string>& types,
    std::vector<std::string> extra_headers) {
  for (const auto& type : types) {
    auto& v = wsg_data["ArrayOf" + type] = {
        .file = "vector",
        .desc = "A list of *" + type + "*\n",
    };
    if (not extra_headers.empty()) {
      v.file = extra_headers.back();
      extra_headers.pop_back();
    }
  }
}

namespace {
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

  wsg_data["ArrayOfScatteringSpecies"] = {
      .file = "scattering/scattering_species.h",
      .desc = "Represents species of scattering particles in the atmosphere.",
  };

  wsg_data["ScatteringSpeciesProperty"] = {
      .file = "scattering/properties.h",
      .desc = R"(Meta data for scattering spefcies.

This is used to identify atmospheric data that are required for scattering calculations.

It is a combination of free-form strings and a *ParticulateProperty* - you need
to see the specific scattering models/methods for what type of data is required.
)",
  };

  wsg_data["SurfacePropertyTag"] = {
      .file = "surf.h",
      .desc = R"--(A surface property.

These tags are part of the keys that can be used to access a *SurfaceField* or *SurfacePoint*.
They are completely free-form and currently not used by ARTS internally.
Instead, they offer a customization point for users to define their own
surface properties and ensures we can access them in a consistent way.

Please see individual surface models/methods for keys that are relevant to run them.
They will generally throw an error if you lack the data.
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

  wsg_data["AtmField"] = {
      .file = "atm.h",
      .desc = R"--(An atmospheric field.

An atmospheric field holds two things:

#. The top of the atmosphere altitude, which is the altitude at which the atmosphere ends.  Unit: m

#. A virtual map of *AtmData*.  The available types of keys are:

   #. *AtmKey* - holds basic data temperature, pressure, wind, and magnetic field.

   #. *SpeciesEnum* - holds volume mixing ratios of absorption species.  Some species might have a different unit.

   #. *SpeciesIsotope* - holds isotopologue ratios of absorption species.  Defaults to built-in values

   #. *QuantumLevelIdentifier* - holds Non-LTE data - i.e., the direct ratios of upper and lower states (which are instead computed on-the-fly when LTE conditions are assumed).

   #. *ScatteringSpeciesProperty* - holds data required for scattering calculations.  The units and type of data are free-form and depend on the scattering model/method.
)--",
  };

  wsg_data["AtmPoint"] = {
      .file = "atm.h",
      .desc = R"--(An atmospheric point.

This can be thought of as the sampling all the *AtmData* of and *AtmField*
at a single altitude-latitude-longitude coordinate.
It, like *AtmField* also acts like a map.  They keys are the same as for *AtmField*. However,
the values are simply the *Numeric* data at that point in the atmosphere.

Always holds an acceptable representation of temperature, pressure, wind, and magnetic field,
although these might be 0 or NaN.

See *AtmField* for the type of data that the atmospheric point can
contain.
)--",
  };

  wsg_data["CallbackOperator"] = {
      .file = "callback.h",
      .desc = R"(Used to inject custom code into *Agenda*.

.. warning::
    This completely breaks the type system of ARTS and should only be used.
    You are on your own when things go wrong with this.
)",
  };

  wsg_data["BlockMatrix"] = {
      .file = "covariance_matrix.h",
      .desc =
          R"(The data for a single *Block*, likely part of a *CovarianceMatrix*.

This holds either a shared *Matrix* or a shared *Sparse* matrix.
)",
  };

  wsg_data["CovarianceMatrix"] = {
      .file = "covariance_matrix.h",
      .desc =
          R"(A covariance matrix is a square matrix that describes the covariance of some property.  

Please see the different workspace variables of this type for more information.

In ARTS, this square matrix is represented by two lists of *Block*.
These are used to give both the covariance matrix and the inverse covariance matrix.
The block-structure allows for efficient storage and computation of the covariance matrix.
)",
  };

  wsg_data["GriddedField2"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 2 dimensional gridded set of *Numeric* data

The grid is a combination of 2 *Vector*

Both the data and the grid may be named.  The grids are not sorted.
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
      .desc = R"(A 2 dimensional array of *Numeric*.

The python mapping allows treating this as a same rank :class:`~numpy.ndarray` in python.
)",
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

When used in the context of the atmosphere via *AtmField* or *AtmPoint*,
this is used to store Non-LTE data - i.e., the direct ratios of upper and
lower states (which are instead computed on-the-fly when LTE conditions
are assumed).
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

#. A map of *SurfaceData*.  The available types of keys are:

   #. *SurfaceKey* - holds basic surface data like elevation and temperature.

   #. *SurfacePropertyTag* - holds free-form surface properties.  The type of data is free-form and depends on the surface model/method.
)--",
  };

  wsg_data["SubsurfaceField"] = {
      .file = "subsurf.h",
      .desc =
          R"--(A sub-surface field.

A sub-surface field effectively holds two things:

#. A *Numeric* of the deepest depth of the subsurface.  Unit: m

#. A map of *SubsurfaceData*.  The available types of keys are:

   #. *SubsurfaceKey*

   #. *SubsurfacePropertyTag*

   See each key for more information on what type of data it allows holding.
)--",
  };

  wsg_data["SubsurfacePropertyTag"] = {
      .file = "subsurf.h",
      .desc = R"--(A custom property tag for subsurface fields data.

These tags are part of the keys that can be used to access a *SubsurfaceField* or *SubsurfacePoint*.
They are completely free-form and currently not used by ARTS internally.
Instead, they offer a customization point for users to define their own
subsurface properties and ensures we can access them in a consistent way.

Please see individual subsurface models/methods for keys that are relevant to run them.
They will generally throw an error if you lack the data.
)--",
  };

  wsg_data["Time"] = {
      .file = "artstime.h",
      .desc = R"(Represents a time stamp
)",
  };

  wsg_data["Vector"] = {
      .file = "matpack.h",
      .desc = R"(A 1 dimensional array of *Numeric*.

The python mapping allows treating this as a same rank :class:`~numpy.ndarray` in python.
)",
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

#. :math:`I` - the total intensity of the light.

#. :math:`Q` - the difference in intensity between horizontally and vertically polarized light.

#. :math:`U` - the difference in intensity between light polarized at +45 degrees and -45 degrees.

#. :math:`V` - the difference in intensity between right and left circularly polarized light.

The python mapping allows treating this 4-long :class:`~numpy.ndarray` in python.
)",
  };

  wsg_data["PropmatVector"] = {
      .file = "rtepack.h",
      .desc = R"(A vector of *Propmat*.

The python mapping allows treating this as a 2-dimensional :class:`~numpy.ndarray` with size 7 as columns.
)",
  };

  wsg_data["MuelmatVector"] = {
      .file = "rtepack.h",
      .desc = R"(A vector of *Muelmat*.

The python mapping allows treating this as a 3-dimensional :class:`~numpy.ndarray` with size 4x4 as rows and columns.
)",
  };

  wsg_data["MuelmatMatrix"] = {
      .file = "rtepack.h",
      .desc = R"(A matrix of *Muelmat*..

The python mapping allows treating this as a 4-dimensional :class:`~numpy.ndarray` with size 4x4 as rows and columns.
)",
  };

  wsg_data["StokvecVector"] = {
      .file = "rtepack.h",
      .desc = R"(A vector of *Stokvec*.

The python mapping allows treating this as a 2-dimensional :class:`~numpy.ndarray` with size 4 as columns.
)",
  };

  wsg_data["PropmatMatrix"] = {
      .file = "rtepack.h",
      .desc = R"(A matrix of *Propmat*.

The python mapping allows treating this as a 3-dimensional :class:`~numpy.ndarray` with size 7 as columns.
)",
  };

  wsg_data["SpecmatMatrix"] = {
      .file = "rtepack.h",
      .desc = R"(A matrix of *Muelmat*.

The python mapping allows treating this as a 4-dimensional :class:`~numpy.ndarray` with size 4x4 as rows and columns.
)",
  };

  wsg_data["StokvecMatrix"] = {
      .file = "rtepack.h",
      .desc = R"(A matrix of *Stokvec*.

The python mapping allows treating this as a 3-dimensional :class:`~numpy.ndarray` with size 4 for columns.
)",
  };

  wsg_data["StokvecSortedGriddedField1"] = {
      .file = "rtepack.h",
      .desc = R"--(A 1-dimensional grid of *Stokvec*.

The grids are 1 *AscendingGrid*.  The grids are fully sorted.
)--",
  };

  wsg_data["StokvecSortedGriddedField2"] = {
      .file = "rtepack.h",
      .desc = R"--(A 2-dimensional grid of *Stokvec*.

The grids are 2 *AscendingGrid*.  The grids are fully sorted.
)--",
  };

  wsg_data["StokvecSortedGriddedField3"] = {
      .file = "rtepack.h",
      .desc = R"--(A 3-dimensional grid of *Stokvec*.

The grids are 3 *AscendingGrid*.  The grids are fully sorted.
)--",
  };

  wsg_data["StokvecSortedGriddedField4"] = {
      .file = "rtepack.h",
      .desc = R"--(A 4-dimensional grid of *Stokvec*.

The grids are 4 *AscendingGrid*.  The grids are fully sorted.
)--",
  };

  wsg_data["StokvecSortedGriddedField5"] = {
      .file = "rtepack.h",
      .desc = R"--(A 5-dimensional grid of *Stokvec*.

The grids are 5 *AscendingGrid*.  The grids are fully sorted.
)--",
  };

  wsg_data["StokvecSortedGriddedField6"] = {
      .file = "rtepack.h",
      .desc = R"--(A 6-dimensional grid of *Stokvec*.

The grids are 6 *AscendingGrid*.  The grids are fully sorted.
)--",
  };

  wsg_data["GriddedSpectralField6"] = {
      .file = "rtepack.h",
      .desc = R"--(A 6-dimensional grid of *Stokvec*.

The grids are altitude x latitude x longitude x zenith x azimuth x frequency of types
*AscendingGrid* x *LatGrid* x *LonGrid* x *ZenithGrid* x *AzimuthGrid* x *AscendingGrid* x.
The grids are fully sorted.
)--",
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

  wsg_data["Vector3"] = {
      .file = "matpack.h",
      .desc = R"(A fixed-size 3D version of *Vector*.

The python mapping allows treating this as a 3-long :class:`~numpy.ndarray`.
)",
  };

  wsg_data["Vector2"] = {
      .file = "matpack.h",
      .desc = R"(A fixed-size 2D version of *Vector*.

The python mapping allows treating this as a 3-long :class:`~numpy.ndarray`.
)",
  };

  wsg_data["DescendingGrid"] = {
      .file = "matpack.h",
      .desc = R"(A sorted *Vector* of always descending values.

The python mapping allows treating this as a :class:`~numpy.ndarray`.
But because it has to be sorted in descending order,
modifying the values are not allowed.
)",
  };

  wsg_data["AscendingGrid"] = {
      .file = "matpack.h",
      .desc = R"(A sorted *Vector* of always ascending values.

The python mapping allows treating this as a :class:`~numpy.ndarray`.
But because it has to be sorted in ascending order,
modifying the values are not allowed.
)",
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
      .desc = R"(Contains name and data about an isotope.

This is used to identify a specific isotope in a species.
The allowed values for the isotope are predefined, see *abs_speciesSet* for available species.

For *PredefinedModelData*, this identifies the predefined model by name
and return any associated data.

For *LinemixingEcsData*, this identifies the linemixing ECS data by name
and returns any associated data.

For *AtmField* and *AtmPoint*, this identifies the isotopologue ratio by name
and returns any associated data.
)",
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
      .desc = R"--(A 1-dimensional grid of *Numeric*.

The grids are 1 *ZenithGrid*.  This grid is sorted.
)--",
  };

  wsg_data["SortedGriddedField1"] = {
      .file = "matpack.h",
      .desc = R"--(A 1-dimensional grid of *Numeric*.

The grids are 1 *AscendingGrid*.  This grid is sorted.
)--",
  };

  wsg_data["SortedGriddedField2"] = {
      .file = "matpack.h",
      .desc = R"--(A 2-dimensional grid of *Numeric*.

The grids are 2 *AscendingGrid*.  The grids are fully sorted.
)--",
  };

  wsg_data["GeodeticField2"] = {
      .file = "matpack.h",
      .desc = R"--(A 2-dimensional grid of *Numeric*.

The grids are *lat_grid* x *lon_grid*.
The types are *LatGrid* x *LonGrid*.  The grids are fully sorted.
)--",
  };

  wsg_data["SortedGriddedField3"] = {
      .file = "matpack.h",
      .desc = R"--(A 3-dimensional grid of *Numeric*.

The grids are 3 *AscendingGrid*.  The grids are fully sorted.
)--",
  };

  wsg_data["SortedGriddedField4"] = {
      .file = "matpack.h",
      .desc = R"--(A 4-dimensional grid of *Numeric*.

The grids are 4 *AscendingGrid*.  The grids are fully sorted.
)--",
  };

  wsg_data["SortedGriddedField5"] = {
      .file = "matpack.h",
      .desc = R"--(A 5-dimensional grid of *Numeric*.

The grids are 5 *AscendingGrid*.  The grids are fully sorted.
)--",
  };

  wsg_data["SortedGriddedField6"] = {
      .file = "matpack.h",
      .desc = R"--(A 6-dimensional grid of *Numeric*.

The grids are 6 *AscendingGrid*.  The grids are fully sorted.
)--",
  };

  wsg_data["ZenithGrid"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 1-dimensional vector of *Numeric* that are guaranteed to be within the range [0, 180].

In addition, the values are sorted in ascending order.

The python mapping allows treating this as a :class:`~numpy.ndarray`.  But because it has to be sorted in ascending order,
modifying the values are not allowed.
)--",
  };

  wsg_data["AzimuthGrid"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 1-dimensional vector of *Numeric* that are guaranteed to be within the range [0, 360).

In addition, the values are sorted in ascending order.

The python mapping allows treating this as a :class:`~numpy.ndarray`.  But because it has to be sorted in ascending order,
modifying the values are not allowed.
)--",
  };

  wsg_data["LonGrid"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 1-dimensional vector of *Numeric* that are guaranteed to be within the range [-180, 180).

In addition, the values are sorted in ascending order.

The python mapping allows treating this as a :class:`~numpy.ndarray`.  But because it has to be sorted in ascending order,
modifying the values are not allowed.
)--",
  };

  wsg_data["LatGrid"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 1-dimensional vector of *Numeric* that are guaranteed to be within the range [-90, 90].

In addition, the values are sorted in ascending order.

The python mapping allows treating this as a :class:`~numpy.ndarray`.  But because it has to be sorted in ascending order,
modifying the values are not allowed.
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
freq_grid : AscendingGrid
    As WSV *freq_grid* **[IN]**
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

  agenda_operators(wsg_data);

  add_arrays_of(wsg_data,
                {
                    "AtmPoint",
                    "SensorObsel",
                    "SpeciesTag",
                    "CIARecord",
                    "SpeciesEnum",
                    "QuantumLevelIdentifier",
                    "XsecRecord",
                    "Sun",
                    "String",
                    "SubsurfacePoint",
                    "PropmatVector",
                    "MuelmatVector",
                    "StokvecVector",
                    "PropmatMatrix",
                    "SpecmatMatrix",
                    "MuelmatTensor3",
                    "StokvecMatrix",
                    "ArrayOfPropmatVector",
                    "ArrayOfStokvecVector",
                    "ArrayOfPropmatMatrix",
                    "ArrayOfStokvecMatrix",
                    "PropagationPathPoint",
                    "ArrayOfPropagationPathPoint",
                    "ArrayOfArrayOfPropagationPathPoint",
                    "Vector3",
                    "AscendingGrid",
                },
                {
                    "optproperties.h",
                });

  return wsg_data;
}
}  // namespace

const std::unordered_map<std::string, WorkspaceGroupRecord>&
internal_workspace_groups() {
  static const auto wsg_data = internal_workspace_groups_creator();
  return wsg_data;
}
