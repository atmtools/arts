#include "workspace_group_friends.h"

#include <arts_options.h>

#include <stdexcept>

#include "workspace_groups.h"

namespace {
std::unordered_map<std::string, WorkspaceGroupRecord> group_friends_internal() {
  std::unordered_map<std::string, WorkspaceGroupRecord> wsg_data;

  wsg_data["Block"] = {
      .file = "covariance_matrix.h",
      .desc =
          R"(This holds a *BlockMatrix* and two *Range* objects and two indices to indicate the position in the *CovarianceMatrix*.
)",
  };

  wsg_data["AbsorptionBand"] = {
      .file = "lbl.h",
      .desc =
          R"(Contains information about a band of related absorption lines.

This information includes

#. A list of *AbsorptionLine*.

#. The line shape profile model.  See *LineByLineLineshape* for available line shape profiles.

#. The frequency cutoff value in [Hz] and type.  See *LineByLineCutoffType* for available cutoff types.

.. note::
    This type does not know about the species that the absorption band/lines belongs to.
    This is why it is often required to keep the *AbsorptionBands* object around.
)",
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

  wsg_data["ZeemanLineModel"] = {
      .file = "lbl.h",
      .desc = R"(Contains information about the Zeeman effect for a single line.

This includes

#. A boolean to indicate if the Zeeman effect is on or off.
#. The upper and lower level statistical weights as *Numeric*.
)",
  };

  wsg_data["QuantumState"] = {
      .file = "lbl.h",
      .desc =
          R"(A map of *QuantumNumberType* to quantum number values in an object with an upper and lower state.
)",
  };

  wsg_data["TemperatureModel"] = {
      .file = "lbl.h",
      .desc =
          R"(Contains information about how a line shape model variable depends on temperature.

See *LineShapeModelType* for available temperature model types.
)",
  };

  wsg_data["LineShapeSpeciesModel"] = {
      .file = "lbl.h",
      .desc =
          R"(Contains information about the line shape for a single line for a single species (or broadener).

This includes information is in the form of a map from *LineShapeModelVariable*
to *TemperatureModel*.
)",
  };

  wsg_data["LineShapeModel"] = {
      .file = "lbl.h",
      .desc = R"(Contains information about the line shape for a single line.

This includes

#. A boolean to indicate if the line-by-line calculations are per species.  This is an experimental feature and should normally be false.
#. The reference temperature as *Numeric*.
#. A map of *LineShapeSpeciesModel* objects for each species in the mixture.
)",
  };

  wsg_data["AbsorptionLine"] = {
      .file = "lbl.h",
      .desc =
          R"(Contains information about related absorption lines.

This type is often found as a member of *AbsorptionBand*.  It contains

#. Einstein A coefficient as *Numeric*.

#. Lower level energy as *Numeric*.

#. Upper level statistical weight as *Numeric*.

#. Lower level statistical weight as *Numeric*.

#. Zeeman model as *ZeemanLineModel*.

#. Line shape model as *LineShapeModel*.

#. Quantum numbers local to this line as *QuantumState*.


.. note::
    This type does not know about the species that the absorption band/lines belongs to.
    This is why it is often required to keep the *AbsorptionBands* object around.
)",
  };

  wsg_data["XsecRecord"] = {
      .file = "xsec_fit.h",
      .desc = R"(A single cross-section record.

These cross-section records contains information about the valid temperature and
pressure ranges as well as well as the fitting coefficients used to compute
and interpolate the cross-section to other temperatures and pressures
)",
  };

  wsg_data["Tensor3"] = {
      .file = "matpack.h",
      .desc = R"(A 3 dimensional array of *Numeric*.

The python mapping allows treating this as a same rank :class:`~numpy.ndarray` in python.
)",
  };

  wsg_data["Tensor4"] = {
      .file = "matpack.h",
      .desc = R"(A 4 dimensional array of *Numeric*.

The python mapping allows treating this as a same rank :class:`~numpy.ndarray` in python.
)",
  };

  wsg_data["Tensor5"] = {
      .file = "matpack.h",
      .desc = R"(A 5 dimensional array of *Numeric*.

The python mapping allows treating this as a same rank :class:`~numpy.ndarray` in python.
)",
  };

  wsg_data["Tensor6"] = {
      .file = "matpack.h",
      .desc = R"(A 6 dimensional array of *Numeric*.

The python mapping allows treating this as a same rank :class:`~numpy.ndarray` in python.
)",
  };

  wsg_data["Tensor7"] = {.file = "matpack.h",
                         .desc = R"(A 7 dimensional array of *Numeric*.

The python mapping allows treating this as a same rank :class:`~numpy.ndarray` in python.
)"};

  wsg_data["SurfacePoint"] = {
      .file = "surf.h",
      .desc =
          R"--(A surface point.

This keeps four things:

#. The elevation in meters as *Numeric*.

#. The temperature in Kelvin as *Numeric*.

#. The local normal vector.

#. A map of the same keys as *SurfaceField* (bar those in *SurfaceKey* that are extracted as above) but towards *Numeric* data.

.. note::
    It is required to keep the *SurfaceField* around if the reference ellipsoid is required.
)--",
  };

  wsg_data["StokvecTensor3"] = {
      .file = "rtepack.h",
      .desc = R"(A *Tensor3* but holds *Stokvec*.

When converted to a :class:`~numpy.ndarray` this will look
like a 4-dimensional array with the last dimension of size 4.
)",
  };

  wsg_data["StokvecTensor4"] = {
      .file = "rtepack.h",
      .desc = R"(A *Tensor4* but holds *Stokvec*.

When converted to a :class:`~numpy.ndarray` this will look
like a 5-dimensional array with the last dimension of size 4.
)",
  };

  wsg_data["StokvecTensor5"] = {
      .file = "rtepack.h",
      .desc = R"(A *Tensor5* but holds *Stokvec*.

When converted to a :class:`~numpy.ndarray` this will look
like a 6-dimensional array with the last dimension of size 4.
)",
  };

  wsg_data["StokvecTensor6"] = {
      .file = "rtepack.h",
      .desc = R"(A *Tensor6* but holds *Stokvec*.

When converted to a :class:`~numpy.ndarray` this will look
like a 7-dimensional array with the last dimension of size 4.
)",
  };

  wsg_data["Specmat"] = {
      .file = "rtepack.h",
      .desc = R"(A single Complex Mueller 4x4 matrix.
)",
  };

  wsg_data["SpeciesTag"] = {
      .file = "species_tags.h",
      .desc = R"(A tagged absorption species

These tags are used to help ARTS identify the species
so that reading routines can find the correct data files.

.. note::
    In previous versions of ARTS, this type had computational meaning.
    This is no longer the case and the type is now only used to help with file IO.
)",
  };

  wsg_data["Range"] = {
      .file = "matpack.h",
      .desc =
          R"(A data type to index a contiguous range inside the multidimensional structures; e.g., *Vector*, *Matrix*, etc.
)",
  };

  wsg_data["SurfaceData"] = {
      .file = "subsurface.h",
      .desc = R"(A data structure for surface field information.

This includes:

#. A data field that is one of *Numeric*, *GeodeticField2*, or *NumericBinaryOperator*.
#. *InterpolationExtrapolation* flags indicating how to extrapolate the data field in its interpolation routine.
)",
  };

  wsg_data["SubsurfaceData"] = {
      .file = "subsurface.h",
      .desc = R"(A data structure for subsurface field information.

This includes:

#. A data field that is one of *Numeric*, *GeodeticField3*, or *NumericTernaryOperator*.
#. *InterpolationExtrapolation* flags indicating how to extrapolate the data field in its interpolation routine.
)",
  };

  wsg_data["SubsurfacePoint"] = {
      .file = "subsurface.h",
      .desc = R"--(A subsurface point.

This keeps three things:

#. The subsurface temperature as *Numeric*.
#. The subsurface density as *Numeric*.
#. A map of the same keys as *SubsurfaceField* (bar those in *SubsurfaceKey* that are extracted as above) but towards *Numeric* data.
)--",
  };

  wsg_data["Sparse"] = {
      .file = "matpack_sparse.h",
      .desc = R"(A sparse version of *Matrix*

This maps to :class:`scipy.sparse.csr_matrix`.
)",
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
      .desc = R"(A position and line-of-sight of a sensor.

This maps to a 5-long :class:`numpy.ndarray` with the first three values
being the position in geodetics [altitude in meters, latitude in degrees, longitude in degrees]
nd the last two values being the line-of-sight
in terms of zenith and azimuth angles [both degrees].
)",
  };

  wsg_data["SensorPosLosVector"] = {
      .file = "obsel.h",
      .desc = R"(Vector of *SensorPosLos*.

This maps to a 2-dimensional :class:`numpy.ndarray` with shape (N, 5)
where N is the number of *SensorPosLos* in the vector.
)",
  };

  wsg_data["Rational"] = {
      .file = "matpack.h",
      .desc = R"(Holds a rational number as two *Index*

.. math::

    x = \frac{a}{b}

where :math:`a` is the numerator and :math:`b` is the denominator.
)",
  };

  wsg_data["QuantumIdentifierNumericMap"] = {
      .file = "lbl.h",
      .desc =
          R"--(A map from *QuantumIdentifier* to *Numeric*.
)--",
      .map_type = true,
  };

  wsg_data["PredefinedModelDataVariant"] = {
      .file = "predef.h",
      .desc =
          R"--(One of the following:

#. :class:`~pyarts3.arts.predef.PredefinedModelDataName`
#. :class:`~pyarts3.arts.predef.PredefinedModelDataWaterDataMTCKD4`

This is used when using predefined models to allow for different types of data
input.  Several types of predefined models have this data built into the code
and will use it directly but must live in the *absorption_predefined_model_data* as a :class:`~pyarts3.arts.predef.PredefinedModelDataName`.
)--",
  };

  wsg_data["PairOfBlockMatrix"] = {
      .file = "retrieval_target.h",
      .desc = R"--(A pair of *BlockMatrix* objects.
)--",
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

Both the data and the grid may be named.  The grids are not sorted.
)--",
  };

  wsg_data["GriddedField5"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 5 dimensional gridded set  of *Numeric* data

The grid is a combination of 5 *Vector*

Both the data and the grid may be named.  The grids are not sorted.
)--",
  };

  wsg_data["GriddedField6"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 6 dimensional gridded set of *Numeric* data

The grid is a combination of 6 *Vector*

Both the data and the grid may be named.  The grids are not sorted.
)--",
  };

  wsg_data["DisortBDRF"] = {
      .file = "disort.h",
      .desc = "A bidirectional reflectance function\n",
  };

  wsg_data["AbsorptionLookupTable"] = {
      .file = "lookup_map.h",
      .desc = R"(A table of lookup calculations.

Effectively holds a *Tensor4* of pre-computed cross-section data.
This table is used to interpole to a pressure, temperature, water vmr, and frequency grid.
)",
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

    #. *GeodeticField3* - The grids are altitude, latitude, longitude.
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

Both the data and the grid may be named.  The grids are not sorted.
)--",
  };

  wsg_data["JacobianTargetType"] = {
      .file = "jacobian.h",
      .desc = R"--(A type of target for use in Jacobian Matrix calculations.

Common for all targets is that they map data to and from the model state to
the model field.  That is they can transform a *Vector* to values in, e.g.,
an *AtmField*, *SurfaceField*, etc., and vice versa they can transform the fields
to values in a *Vector*.

See *model_state_vector* and method involving it for more information.
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

  wsg_data["NamedGriddedField2"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 3 dimensional gridded set of *Numeric* data

The grid is a combination of 1 *ArrayOfString* and 2 *Vector*

Both the data and the grid may be named.  The grids are not sorted.
)--",
  };

  wsg_data["GriddedField3"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 3 dimensional gridded set of *Numeric* data

The grid is a combination of 3 *Vector*

Both the data and the grid may be named.  The grids are not sorted.
)--",
  };

  wsg_data["NamedGriddedField3"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 4 dimensional gridded set of *Numeric* data

The grid is a combination of 1 *ArrayOfString* and 3 *Vector*

Both the data and the grid may be named.  The grids are not sorted.
)--",
  };

  wsg_data["MuelmatTensor3"] = {
      .file = "rtepack.h",
      .desc = R"(A *Tensor3* of *Muelmat*.

When converted to a :class:`~numpy.ndarray` this will look
like a 5-dimensional array with the last two dimensions of size 4.
)",
  };

  wsg_data["Muelmat"] = {
      .file = "rtepack.h",
      .desc = "A single Mueller 4x4 matrix.\n",
  };

  wsg_data["GriddedField1"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 1 dimensional gridded set of *Numeric* data

The grid is 1 *Vector*

Both the data and the grid may be named.  The grid is not sorted.
)--",
  };

  wsg_data["GriddedField1Named"] = {
      .file = "matpack.h",
      .desc =
          R"--(A 3 dimensional gridded set of *Numeric* data

The grid is a combination of 1 *Vector* and 1 *ArrayOfString*

Both the data and the grid may be named.  The grids are not sorted.
)--",
  };

  wsg_data["GeodeticField3"] = {
      .file = "rtepack.h",
      .desc = R"--(A 3-dimensional gridof *Numeric*.

The grids are *alt_grid* x *latitude_grid* x *longitude_grid*.
The types are *AscendingGrid* x *LatGrid* x *LonGrid*.  The grids are all sorted.
)--",
  };

  add_arrays_of(wsg_data,
                {
                    "Index",
                    "ArrayOfGriddedField1",
                    "ArrayOfGriddedField2",
                    "ArrayOfGriddedField3",
                    "ArrayOfIndex",
                    "ArrayOfMatrix",
                    "ArrayOfScatteringMetaData",
                    "ArrayOfSingleScatteringData",
                    "ArrayOfString",
                    "ArrayOfTensor3",
                    "ArrayOfTensor6",
                    "ArrayOfTime",
                    "Tensor3",
                    "Tensor4",
                    "Tensor5",
                    "Tensor6",
                    "Tensor7",
                    "SubsurfacePoint",
                    "StokvecTensor3",
                    "SpeciesIsotope",
                    "Sparse",
                    "NamedGriddedField2",
                    "MuelmatMatrix",
                    "GriddedField3",
                    "GriddedField4",
                    "GriddedField2",
                    "GriddedField1",
                    "ArrayOfMuelmatMatrix",
                    "ArrayOfMuelmatVector",
                    "ArrayOfVector",
                    "GriddedField1Named",
                    "Matrix",
                    "QuantumIdentifier",
                    "ScatteringMetaData",
                    "SingleScatteringData",
                    "Time",
                    "Vector",
                    "Vector2",
                },
                {});

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
