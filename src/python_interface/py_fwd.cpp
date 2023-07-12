#include <fwd.h>
#include <py_auto_interface.h>
#include <pybind11/attr.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "debug.h"
#include "py_macros.h"

#define SingleModel(varname, FullName, SubName, docstr)                     \
  auto varname =                                                            \
      py::class_<FullName>(fwd, SubName)                                    \
          .def(py::init(                                                    \
                   [](Numeric T,                                            \
                      Numeric P,                                            \
                      const SpeciesIsotopologueRatios& isotopologue_ratios, \
                      const ArrayOfArrayOfSpeciesTag& allspecs,             \
                      const Vector& allvmrs,                                \
                      const AbsorptionLines& band,                          \
                      Index line) {                                         \
                     ARTS_USER_ERROR_IF(line >= band.NumLines(),            \
                                        "Out of bounds line index");        \
                     ARTS_USER_ERROR_IF(allspecs.nelem() != allvmrs.size(), \
                                        "Size mismatch");                   \
                     return new FullName{T,                                 \
                                         P,                                 \
                                         isotopologue_ratios,               \
                                         allspecs,                          \
                                         allvmrs,                           \
                                         band,                              \
                                         line};                             \
                   }),                                                      \
               "From values")                                               \
          .def(                                                             \
              "at",                                                         \
              py::vectorize(&FullName::at<>),                               \
              py::arg("f"),                                                 \
              py::doc(                                                      \
                  R"--(This times f * (1 - exp(hf / kT)) is the absorption coeff\
\
Parameters\
----------\
f : float or ~numpy.ndarray\
    Frequency in Hz\
\
Returns\
-------\
complex or ~numpy.ndarray\
    The absorption coefficient\
)--"));                                                                     \
  varname.doc() = docstr;

#define BandModel(varname, FullName, SubName, docstr)                       \
  auto varname =                                                            \
      py::class_<FullName>(fwd, SubName)                                    \
          .def(py::init(                                                    \
                   [](Numeric T,                                            \
                      Numeric P,                                            \
                      const SpeciesIsotopologueRatios& isotopologue_ratios, \
                      const ArrayOfArrayOfSpeciesTag& allspecs,             \
                      const Vector& allvmrs,                                \
                      const ArrayOfArrayOfAbsorptionLines& specbands) {     \
                     ARTS_USER_ERROR_IF(allspecs.nelem() != allvmrs.size(), \
                                        "Size mismatch");                   \
                     return new FullName{T,                                 \
                                         P,                                 \
                                         isotopologue_ratios,               \
                                         allspecs,                          \
                                         allvmrs,                           \
                                         specbands};                        \
                   }),                                                      \
               "From values")                                               \
          .def("size", &FullName::size, "Return number of lines")           \
          .def(                                                             \
              "at",                                                         \
              [](const FullName& self, const Numeric& f) {                  \
                return self.at(f);                                          \
              },                                                            \
              py::arg("f"),                                                 \
              py::doc(R"--(This returns the complex absorption coefficient\
\
Parameters\
----------\
f : float\
    Frequency in Hz\
\
Returns\
-------\
F : ~pyarts.arts.ComplexVector\
    The absorption coefficients\
)--"))                                                                      \
          .def(                                                             \
              "at",                                                         \
              [](const FullName& self, const Vector& f) {                   \
                return self.at(f);                                          \
              },                                                            \
              py::arg("f"),                                                 \
              py::doc(R"--(This returns the complex absorption coefficient\
\
Parameters\
----------\
f : ~pyarts.arts.Vector\
    Frequency in Hz\
\
Returns\
-------\
F : ~pyarts.arts.ComplexVector\
    The absorption coefficients\
)--"));                                                                     \
  varname.doc() = docstr;

namespace Python {
void py_fwd(py::module_& m) {
  auto fwd = m.def_submodule("fwd");

  SingleModel(
      lbl_mtckd_single,
      fwd::lbl::mtckd::single,
      "lbl_mtckd_single",
      "This is a single line model fitting the MT_CKD model w/o line mixing.");
  SingleModel(
      lbl_mtckd_single_lm,
      fwd::lbl::mtckd::single_lm,
      "lbl_mtckd_single_lm",
      "This is a single line model fitting the MT_CKD model w/ line mixing.");

  BandModel(lbl_mtckd_band,
            fwd::lbl::mtckd::band,
            "lbl_mtckd_band",
            "This is a full model fitting the MT_CKD model w/o line mixing.");
  BandModel(lbl_mtckd_band_lm,
            fwd::lbl::mtckd::band_lm,
            "lbl_mtckd_band_lm",
            "(This is a full model fitting the MT_CKD model w/ line mixing.");
  BandModel(lbl_full,
            fwd::lbl::full,
            "lbl",
            "This is a collection of line-by-line models.");
  lbl_full.def(
      "at_par",
      [](const fwd::lbl::full& self, const Vector& f) {
        return self.at_par(f);
      },
      "Parallel version of :func:`at`");

  py::class_<ForwardAbsorption>(m, "ForwardAbsorption");
  py::class_<ForwardRadiance>(m, "ForwardRadiance")
      .def(py::init<>(), "From nothing")
      .PythonInterfaceCopyValue(ForwardRadiance)
      .PythonInterfaceWorkspaceVariableConversion(ForwardRadiance)
      .PythonInterfaceFileIO(ForwardRadiance)
      .PythonInterfaceBasicRepresentation(ForwardRadiance)
      .def(
          "planar",
          [](const ForwardRadiance& fwd_prof, Numeric f, Numeric za) {
            return fwd_prof.planar(f, za);
          },
          py::arg("f"),
          py::arg("za"),
          py::doc(R"--(Plane parallel computer of radiance profile

Parameters
----------
f : float
    Frequency [Hz]
za : float
    Zenith angle [degrees]

Returns
-------
prof : ~pyarts.arts.Vector
    Profile of radiation
)--"))
      .def(
          "planar_par",
          [](const ForwardRadiance& fwd_prof, const Vector& f, Numeric za) {
            return fwd_prof.planar_par(f, za);
          },
          py::arg("f"),
          py::arg("za"),
          py::doc(R"--(Plane parallel computer of radiance profile

Parameters
----------
f :  ~pyarts.arts.Vector
    Frequency [Hz]
za : float
    Zenith angle [degrees]

Returns
-------
prof : ~pyarts.arts.Matrix
    Profile of radiation
)--"))
      .PythonInterfaceWorkspaceDocumentation(ForwardRadiance);

  py::class_<ForwardIrradiance>(m, "ForwardIrradiance")
      .def(py::init<>(), "From nothing")
      .def(py::init<ForwardRadiance>(), "From forward radiance")
      .PythonInterfaceCopyValue(ForwardIrradiance)
      .PythonInterfaceWorkspaceVariableConversion(ForwardIrradiance)
      .PythonInterfaceFileIO(ForwardIrradiance)
      .PythonInterfaceBasicRepresentation(ForwardIrradiance)
      .def_property_readonly(
          "altitude",
          [](const ForwardIrradiance& fwd_prof) {
            return Vector(fwd_prof.rad.altitude);
          },
          py::doc("The altitudes the irradiance is computed at"))
      .def(
          "planar",
          [](const ForwardIrradiance& fwd_prof, Numeric f, const Vector& za) {
            return fwd_prof.planar(f, za);
          },
          py::arg("f"),
          py::arg("za"),
          py::doc(R"--(Plane parallel computer of radiance profile

Parameters
----------
f : float
    Frequency [Hz]
za : float
    Zenith angle [degrees]

Returns
-------
prof : ~pyarts.arts.Vector
    Profile of radiation
)--"))
      .PythonInterfaceWorkspaceDocumentation(ForwardIrradiance);
}
}  // namespace Python