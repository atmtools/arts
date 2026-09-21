#include <debug.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/vector.h>
#include <python_interface.h>
#include <rtepack.h>
#include <rtepack_surface.h>

#include "hpy_arts.h"
#include "hpy_vector.h"

namespace Python {
// Defined in py_rtepack_stokvec.cpp, py_rtepack_propmat.cpp,
// py_rtepack_muelmat.cpp and py_rtepack_specmat.cpp respectively.
// py_rtepack.cpp used to bind all of these types in a single ~760-line
// function; that made it one of the most expensive translation units in the
// whole build (both in compile time and peak RAM), since nanobind's binding
// machinery is instantiated independently for every bound type. Splitting
// the (independent) groups across a handful of TUs lets the compiler
// process them in parallel and keeps any single TU's memory footprint down,
// without changing what gets bound.
void py_rtepack_stokvec(py::module_ &m);
void py_rtepack_propmat(py::module_ &m);
void py_rtepack_muelmat(py::module_ &m);
void py_rtepack_specmat(py::module_ &m);

void py_rtepack(py::module_ &m) try {
  // These must run first: the array-of-array bindings and the plain arrays
  // below (ArrayOfPropmatVector, etc.) require the element types' nanobind
  // bindings to already be registered.
  py_rtepack_stokvec(m);
  py_rtepack_propmat(m);
  py_rtepack_muelmat(m);
  py_rtepack_specmat(m);

  auto a1 = py::bind_vector<ArrayOfPropmatVector, py::rv_policy::reference_internal>(m, "ArrayOfPropmatVector");
  generic_interface(a1);
  vector_interface(a1);
  auto a2 =
      py::bind_vector<ArrayOfArrayOfPropmatVector, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfPropmatVector");
  generic_interface(a2);
  vector_interface(a2);
  auto a3 = py::bind_vector<ArrayOfPropmatMatrix, py::rv_policy::reference_internal>(m, "ArrayOfPropmatMatrix");
  generic_interface(a3);
  vector_interface(a3);
  auto a4 =
      py::bind_vector<ArrayOfArrayOfPropmatMatrix, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfPropmatMatrix");
  generic_interface(a4);
  vector_interface(a4);

  auto b1 = py::bind_vector<ArrayOfMuelmatVector, py::rv_policy::reference_internal>(m, "ArrayOfMuelmatVector");
  b1.def(
      "is_polarized",
      [](const ArrayOfMuelmatVector &m) {
        return stdr::any_of(m, [](const MuelmatVector &mm) {
          return stdr::any_of(mm, [](const Muelmat &m) { return m.is_polarized(); });
        });
      },
      "Check if the Mueller matrix represents a polarized state.");
  generic_interface(b1);
  vector_interface(b1);
  auto b2 =
      py::bind_vector<ArrayOfArrayOfMuelmatVector, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfMuelmatVector");
  generic_interface(b2);
  vector_interface(b2);
  auto b3 = py::bind_vector<ArrayOfMuelmatMatrix, py::rv_policy::reference_internal>(m, "ArrayOfMuelmatMatrix");
  generic_interface(b3);
  vector_interface(b3);
  auto b4 =
      py::bind_vector<ArrayOfArrayOfMuelmatMatrix, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfMuelmatMatrix");
  generic_interface(b4);
  vector_interface(b4);
  auto b5 = py::bind_vector<ArrayOfMuelmatTensor3, py::rv_policy::reference_internal>(m, "ArrayOfMuelmatTensor3");
  generic_interface(b5);
  vector_interface(b5);

  auto c1 = py::bind_vector<ArrayOfStokvecVector, py::rv_policy::reference_internal>(m, "ArrayOfStokvecVector");
  generic_interface(c1);
  vector_interface(c1);
  auto c2 =
      py::bind_vector<ArrayOfArrayOfStokvecVector, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfStokvecVector");
  generic_interface(c2);
  vector_interface(c2);
  auto c3 = py::bind_vector<ArrayOfStokvecMatrix, py::rv_policy::reference_internal>(m, "ArrayOfStokvecMatrix");
  generic_interface(c3);
  vector_interface(c3);
  auto c4 =
      py::bind_vector<ArrayOfArrayOfStokvecMatrix, py::rv_policy::reference_internal>(m, "ArrayOfArrayOfStokvecMatrix");
  generic_interface(c4);
  vector_interface(c4);
  auto c5 = py::bind_vector<ArrayOfStokvecTensor3, py::rv_policy::reference_internal>(m, "ArrayOfStokvecTensor3");
  generic_interface(c5);
  vector_interface(c5);

  auto d3 = py::bind_vector<ArrayOfSpecmatMatrix, py::rv_policy::reference_internal>(m, "ArrayOfSpecmatMatrix");
  generic_interface(d3);
  vector_interface(d3);

  auto rtepack  = m.def_submodule("rtepack");
  rtepack.doc() = "Interface to some of the core RTE functionality";

  rtepack.def(
      "two_level_radiative_transfer",
      [](const ArrayOfMuelmatVector  &Ts,
         const ArrayOfMuelmatTensor3 &dTs,
         const ArrayOfStokvecVector  &Js,
         const ArrayOfStokvecMatrix  &dJs,
         const StokvecVector         &I0) {
        StokvecVector        I;
        ArrayOfStokvecMatrix dI;

        const auto Pi = forward_cumulative_transmission(Ts);
        rtepack::two_level_linear_emission_step_by_step_full(I, dI, Ts, Pi, dTs, Js, dJs, I0);

        return std::pair{I, dI};
      },
      "Ts"_a,
      "dTs"_a,
      "Js"_a,
      "dJs"_a,
      "I0"_a,
      "Returns the two-level radiative transfer of the input matrices");

  auto tramat = py::class_<TransmittanceMatrix>(m, "TransmittanceMatrix");
  generic_interface(tramat);
  tramat.def_rw(
      "T",
      &TransmittanceMatrix::T,
      "The transmittance Mueller matrix; shape [nf, np] if defined\n\n.. :class:`~pyarts3.arts.MuelmatMatrix");
  tramat.def_rw("T_diag_m1",
                &TransmittanceMatrix::T_diag_m1,
                "The accurately cached diagonal of T-I; shape [nf, np] if defined\n\n.. "
                ":class:`~pyarts3.arts.StokvecMatrix");
  tramat.def_rw(
      "dT",
      &TransmittanceMatrix::dT,
      "The derivative of the transmittance Mueller matrix; shape [2, nf, np, nq] if defined\n\n.. :class:`~pyarts3.arts.MuelmatTensor4");
  tramat.def_rw(
      "P",
      &TransmittanceMatrix::P,
      "The cumulative from background transmittance Mueller matrix; shape [nf, np] if defined\n\n.. :class:`~pyarts3.arts.MuelmatMatrix");
  tramat.def_rw(
      "L",
      &TransmittanceMatrix::L,
      "The linear evolution Mueller matrix; shape [nf, np] if defined\n\n.. :class:`~pyarts3.arts.MuelmatMatrix");
  tramat.def_rw("L_diag_m1",
                &TransmittanceMatrix::L_diag_m1,
                "The accurately cached diagonal of L-I; shape [nf, np] if defined\n\n.. "
                ":class:`~pyarts3.arts.StokvecMatrix");
  tramat.def_rw(
      "dL",
      &TransmittanceMatrix::dL,
      "The derivative of the linear evolution Mueller matrix; shape [2, nf, np, nq] if defined\n\n.. :class:`~pyarts3.arts.MuelmatTensor4");
  tramat.def_rw("opt",
                &TransmittanceMatrix::option,
                "The option for the transmittance matrix calculations\n\n.. :class:`~pyarts3.arts.TransmittanceOption");

  auto srcvec = py::class_<SourceVector>(m, "SourceVector");
  generic_interface(srcvec);
  srcvec.def_rw("J", &SourceVector::J, "The source vectors; shape [nf, np]\n\n.. :class:`~pyarts3.arts.StokvecMatrix");
  srcvec.def_rw("dJ",
                &SourceVector::dJ,
                "The derivatives of the source vectors; shape [nf, np, nq]\n\n.. :class:`~pyarts3.arts.StokvecTensor3");

  auto rp  = m.def_submodule("rtepack");
  rp.doc() = "Module for RTEPACK functionality";

  auto tr = py::class_<rtepack::tran>(rp, "tran");
  generic_interface(tr);
  tr.doc() = "Class for computing the transmission Mueller matrix and its derivative";
  tr.def(py::init<Propmat, Propmat, Numeric>(), "k1"_a, "k2"_a, "r"_a)
      .def("__call__",
           static_cast<Muelmat (rtepack::tran::*)() const noexcept>(&rtepack::tran::operator()),
           "Returns the Mueller matrix")
      .def("deriv",
           static_cast<Muelmat (rtepack::tran::*)(
               const Muelmat &, const Propmat &, const Propmat &, const Propmat &, Numeric, Numeric) const>(
               &rtepack::tran::deriv),
           "t"_a,
           "k1"_a,
           "k2"_a,
           "dk"_a,
           "r"_a,
           "dr"_a,
           "Returns the derivative of the Mueller matrix")
      .def("linsrc",
           static_cast<Muelmat (rtepack::tran::*)() const noexcept>(&rtepack::tran::linsrc),
           "Returns the linear-in-tau evolve operator")
      .def("linsrc_deriv",
           static_cast<Muelmat (rtepack::tran::*)(const Propmat &, Numeric, Numeric) const>(
               &rtepack::tran::linsrc_deriv),
           "dk"_a,
           "r"_a,
           "dr"_a,
           "Returns the derivative of the linear-in-tau evolve operator")
      .def("expm1", &rtepack::tran::expm1, "Returns the Mueller matrix minus the identity matrix");

  rp.def("sqrt",
         &rtepack::sqrt,
         "K"_a,
         R"(Returns the square root of the Propagation matrix as a spectral matrix

.. math::
    \mathrm{K} = \mathrm{S} \mathrm{S},

The return value of this method is the spectral matrix :math:`\mathrm{S}`.

Parameters
----------
K : Specmat
    The propagation matrix

Returns
-------
Specmat
    The square root of the Propagation matrix as a spectral matrix.
)")
      .def("logK",
           &rtepack::logK,
           "m"_a,
           R"(Returns the logarithm of the Mueller matrix as a Propagation matrix

This only works for a Mueller matrix that has been generated from a Propagation matrix using

.. math::
    \mathrm{M} = e^{-\mathrm{K} r},

where :math:`\mathrm{M}` is the Mueller matrix, :math:`\mathrm{K}` is the Propagation matrix, and :math:`r` is the distance.

The return value of this method is then not the original Propagation matrix, but rather the logarithm of the Mueller matrix
as if it were generated from a Propagation matrix using the same equation.  That is

.. math::
    -\mathrm{K}r = \log(\mathrm{M})

Parameters
----------
m : Muelmat
    The Mueller matrix to compute the logarithm of.

Returns
-------
Propmat
    The logarithm of the Mueller matrix as a Propagation matrix.
)")
      .def("specular_reflected_direction",
           &rtepack::specular_reflected_direction,
           "k_inc"_a,
           "n_surface"_a,
           R"(Return the direction of the specularly reflected propagation vector.)")
      .def("fresnel_reflectance",
           &rtepack::fresnel_reflectance,
           "Rv"_a,
           "Rh"_a,
           R"(Return the Fresnel Mueller matrix for the complex amplitude
coefficients `Rv` and `Rh`.)")
      .def("fresnel_reflectance_specular",
           &rtepack::fresnel_reflectance_specular,
           "Rv"_a,
           "Rh"_a,
           "k_inc"_a,
           "n_surface"_a,
           R"(Return the Fresnel Mueller matrix for specular reflection. `k_inc`
is the incident propagation vector toward the surface and `n_surface` is the
outward surface normal.)")
      .def("fresnel_reflectance_nonspecular",
           &rtepack::fresnel_reflectance_nonspecular,
           "Rv"_a,
           "Rh"_a,
           "k_inc"_a,
           "k_out"_a,
           "n_surface"_a,
           R"(Return the Fresnel Mueller matrix for non-specular reflection with
independent incident and outgoing directions.)");

  rp.def("specular_reflected_direction",
         &rtepack::specular_reflected_direction,
         "k_inc"_a,
         "n_surface"_a,
         R"(Return the specular reflection direction.

Parameters
----------
k_inc : ~pyarts3.arts.Vector3
  Unit propagation vector of the incident beam (toward surface)
n_surface : ~pyarts3.arts.Vector3
  Outward unit surface normal

Returns
-------
k_ref : ~pyarts3.arts.Vector3
  Unit propagation vector of the specularly reflected beam (away from surface)
)");

  rp.def(
      "nonspecular_radiance_from_patches",
      [](const std::vector<Vector2> &coords,
         const StokvecVector        &sources,
         const Stokvec              &J,
         Complex                     Rv,
         Complex                     Rh,
         Vector2                     pos,
         Numeric                     h_pos,
         const Vector3              &n_surface,
         const Vector3              &k_out,
         Vector2                     ellipsoid,
         const GeodeticField2       &hfield) {
        return rtepack::nonspecular_radiance_from_patches(
            coords, sources, J, Rv, Rh, pos, h_pos, n_surface, k_out, ellipsoid, hfield);
      },
      "coords"_a,
      "sources"_a,
      "J"_a,
      "Rv"_a,
      "Rh"_a,
      "pos"_a,
      "h_pos"_a,
      "n_surface"_a,
      "k_out"_a,
      "ellipsoid"_a,
      "hfield"_a,
      R"(Accumulate non-specular scattered radiance from visible surface patches.

Integrates the reflectance contribution from each visible patch j:

    L_out = J + (1/pi) sum_j  R(k_j, k_out) . L_j . cos(theta_P) . dOmega_j

where dOmega_j is the solid angle of patch j as seen from the scatter point.

Parameters
----------
coords : list of ~pyarts3.arts.Vector2
  Visible (lat, lon) pairs [degrees] - typically from visible_coordinates()
sources : ~pyarts3.arts.StokvecVector
  Source radiance at each visible coordinate
J : ~pyarts3.arts.Stokvec
  Thermal emission Stokes vector at the scatter point
Rv : complex
  Fresnel amplitude for vertical polarisation
Rh : complex
  Fresnel amplitude for horizontal polarisation
pos : ~pyarts3.arts.Vector2
  Scatter-point (lat, lon) [degrees]
h_pos : float
  Scatter-point height [m]
n_surface : ~pyarts3.arts.Vector3
  Outward unit normal at the scatter point (ECEF)
k_out : ~pyarts3.arts.Vector3
  Unit propagation direction toward the sensor (ECEF)
ellipsoid : ~pyarts3.arts.Vector2
  Ellipsoid semi-axes (a, b) [m]
hfield : ~pyarts3.arts.GeodeticField2
  Height field - provides grid spacing and patch heights [m]

Returns
-------
L_out : ~pyarts3.arts.Stokvec
  Outgoing Stokes vector toward the sensor
)");

} catch (std::exception &e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize rtepack\n{}", e.what()));
}
}  // namespace Python
