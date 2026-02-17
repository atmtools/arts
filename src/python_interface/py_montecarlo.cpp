#include <mc_antenna.h>
#include <python_interface.h>

namespace Python {
void py_montecarlo(py::module_& m) try {
  py::class_<MCAntenna> mca(m, "MCAntenna");

  mca.def(py::init<>())
      .def_rw("atype",
              &MCAntenna::atype,
              "The type of antenna pattern to use\n\n.. :class:`AntennaType`")
      .def_rw("sigma_aa",
              &MCAntenna::sigma_aa,
              "The spread of azimith to use\n\n.. :class:`Numeric`")
      .def_rw("sigma_za",
              &MCAntenna::sigma_za,
              "The spread of zenith to use\n\n.. :class:`Numeric`")
      .def_rw("aa_grid",
              &MCAntenna::aa_grid,
              "The azimuth grid\n\n.. :class:`Vector`")
      .def_rw("za_grid",
              &MCAntenna::za_grid,
              "The zenith grid\n\n.. :class:`Vector`")
      .def_rw("G_lookup",
              &MCAntenna::G_lookup,
              "The lookup table for the antenna gain\n\n.. :class:`Matrix`")
      .def("set_pencil_beam",
           &MCAntenna::set_pencil_beam,
           "Set the antenna pattern to a pencil beam")
      .def("set_gaussian",
           &MCAntenna::set_gaussian,
           R"(Set the antenna pattern to a Gaussian

Parameters
----------
za_sigma : Numeric
    The standard deviation of the Gaussian in the zenith direction
aa_sigma : Numeric
    The standard deviation of the Gaussian in the azimuth direction
)",
           "za_sigma"_a,
           "aa_sigma"_a)
      .def("set_gaussian_fwhm",
           &MCAntenna::set_gaussian_fwhm,
           R"(Set the antenna pattern to a Gaussian

Parameters
----------
za_fwhm : Numeric
    The full width half maximum of the Gaussian in the zenith direction
aa_fwhm : Numeric
    The full width half maximum of the Gaussian in the azimuth direction
)",
           "za_fwhm"_a,
           "aa_fwhm"_a)
      .def(
          "set_lookup",
          [](MCAntenna& self,
             const Vector& za,
             const Vector& aa,
             const Matrix& G) { self.set_lookup(za, aa, G); },
          R"(Set the antenna pattern from a lookup table

Parameters
----------
za : Vector
    The zenith angle grid for the lookup table
aa : Vector
    The azimuth angle grid for the lookup table
G : Matrix
    The gain values for the lookup table, where G[i, j] corresponds to za[i] and aa[j])",
          "za"_a,
          "aa"_a,
          "G"_a)
      .def("return_los",
           &MCAntenna::return_los,
           R"(Return the line of sight for a given azimuth and zenith angle
Parameters
----------
R_return : Matrix33
    The rotation matrix from the return path to the ENU frame
R_enu2ant : Matrix33
    The rotation matrix from the ENU frame to the antenna frame
Returns
-------
Numeric
    The line of sight for the given angles, normalized by the antenna gain at the boresight
)",
           "R_return"_a,
           "R_enu2ant"_a)
      .def("draw_los",
           &MCAntenna::draw_los,
           R"(Draw a random line of sight for the current antenna pattern
Parameters
----------
rng : RandomNumberGenerator
    A random number generator to use for sampling the antenna pattern
R_ant2enu : Matrix33
    The rotation matrix from the antenna frame to the ENU frame
bore_sight_los : Vector2
    The line of sight for the boresight of the antenna, in terms of zenith and azimuth angles
Returns
-------
(sampled_rte_los, R_los) : Tuple[Vector2, Matrix33]
    sampled_rte_los : Vector2
        The sampled line of sight in terms of zenith and azimuth angles
    R_los : Matrix33
        The rotation matrix corresponding to the sampled line of sight
)",
           "rng"_a,
           "R_ant2enu"_a,
           "bore_sight_los"_a);

  mca.doc() = "Monte Carlo Antenna pattern class";
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize montecarlo\n{}", e.what()));
}
}  // namespace Python