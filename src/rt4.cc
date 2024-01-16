/*!
  \file   rt4.cc
  \author Jana Mendrok <jana.mendrok@gmail.com>
  \date   2016-05-24
  
  \brief  Contains functions related to application of scattering solver RT4.
  
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts_constants.h"
#include "arts_conversions.h"
#include "config.h"

#ifdef ENABLE_RT4

#include <workspace.h>

#include <complex>
#include <stdexcept>

#include "check_input.h"
#include "interpolation.h"
#include "m_xml.h"
#include "optproperties.h"
#include "physics_funcs.h"
#include "rt4.h"
#include "special_interp.h"

inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);
inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;
inline constexpr Numeric COSMIC_BG_TEMP=Constant::cosmic_microwave_background_temperature;

void check_rt4_input(  // Output
    Index& nhstreams,
    Index& nhza,
    Index& nummu,
    // Input
    const Index& cloudbox_on,
    const Index& atmfields_checked,
    const Index& atmgeom_checked,
    const Index& cloudbox_checked,
    const Index& scat_data_checked,
    const ArrayOfIndex& cloudbox_limits,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    
    const Index& nstreams,
    const String& quad_type,
    const Index& add_straight_angles,
    const Index& pnd_ncols) {
  // Don't do anything if there's no cloudbox defined.
  ARTS_USER_ERROR_IF (!cloudbox_on,
        "Cloudbox is off, no scattering calculations to be"
        " performed.");

  ARTS_USER_ERROR_IF (atmfields_checked != 1,
        "The atmospheric fields must be flagged to have"
        " passed a consistency check (atmfields_checked=1).");
  ARTS_USER_ERROR_IF (atmgeom_checked != 1,
        "The atmospheric geometry must be flagged to have"
        " passed a consistency check (atmgeom_checked=1).");
  ARTS_USER_ERROR_IF (cloudbox_checked != 1,
        "The cloudbox must be flagged to have"
        " passed a consistency check (cloudbox_checked=1).");
  ARTS_USER_ERROR_IF (scat_data_checked != 1,
        "The scat_data must be flagged to have "
        "passed a consistency check (scat_data_checked=1).");

  ARTS_USER_ERROR_IF (4 < 0 || 4 > 2,
        "For running RT4, the dimension of stokes vector"
        " must be 1 or 2.\n");

  ARTS_USER_ERROR_IF (cloudbox_limits[0] != 0,
    "RT4 calculations currently only possible with"
    " lower cloudbox limit\n"
    "at 0th atmospheric level"
    " (assumes surface there, ignoring z_surface).\n")

  ARTS_USER_ERROR_IF (cloudbox_limits.size() != 6,
        "*cloudbox_limits* is a vector which contains the"
        " upper and lower limit of the cloud for all"
        " atmospheric dimensions. So its dimension must"
        " be 6");

  ARTS_USER_ERROR_IF (scat_data.empty(),
        "No single scattering data present.\n"
        "See documentation of WSV *scat_data* for options to"
        " make single scattering data available.\n");

  ARTS_USER_ERROR_IF (pnd_ncols != 1,
        "*pnd_field* is not 1D! \n"
        "RT4 can only be used for 1D!\n");

  ARTS_USER_ERROR_IF (quad_type.length() > 1,
    "Input parameter *quad_type* not allowed to contain more than a"
    " single character.\n"
    "Yours has ", quad_type.length(), ".\n")

  if (quad_type == "D" || quad_type == "G") {
    if (add_straight_angles)
      nhza = 1;
    else
      nhza = 0;
  } else if (quad_type == "L") {
    nhza = 0;
  } else {
    ARTS_USER_ERROR (
      "Unknown quadrature type: ", quad_type,
      ".\nOnly D, G, and L allowed.\n")
  }

  // RT4 actually uses number of angles in single hemisphere. However, we don't
  // want a bunch of different approaches used in the interface, so we apply the
  // DISORT (and ARTS) way of total number of angles here. Hence, we have to
  // ensure here that total number is even.
  ARTS_USER_ERROR_IF (nstreams / 2 * 2 != nstreams,
    "RT4 requires an even number of streams, but yours is ", nstreams,
    ".\n")
  
  nhstreams = nstreams / 2;
  // nummu is the total number of angles in one hemisphere, including both
  // the quadrature angles and the "extra" angles.
  nummu = nhstreams + nhza;

  // RT4 can only completely or azimuthally randomly oriented particles.
  bool no_arb_ori = true;
  for (Size i_ss = 0; i_ss < scat_data.size(); i_ss++)
    for (Size i_se = 0; i_se < scat_data[i_ss].size(); i_se++)
      if (scat_data[i_ss][i_se].ptype != PTYPE_TOTAL_RND &&
          scat_data[i_ss][i_se].ptype != PTYPE_AZIMUTH_RND)
        no_arb_ori = false;
  
  ARTS_USER_ERROR_IF (!no_arb_ori,
    "RT4 can only handle scattering elements of type ", PTYPE_TOTAL_RND,
    " (", PTypeToString(PTYPE_TOTAL_RND), ") and\n",
    PTYPE_AZIMUTH_RND, " (", PTypeToString(PTYPE_AZIMUTH_RND),
    "),\n"
    "but at least one element of other type (", PTYPE_GENERAL, "=",
    PTypeToString(PTYPE_GENERAL), ") is present.\n")

  //------------- end of checks ---------------------------------------
}

void get_quad_angles(  // Output
    VectorView mu_values,
    VectorView quad_weights,
    Vector& za_grid,
    Vector& aa_grid,
    // Input
    const String& quad_type,
    const Index& nhstreams,
    const Index& nhza,
    const Index& nummu) {
  if (quad_type == "D") {
    double_gauss_quadrature_(
        nhstreams, mu_values.unsafe_data_handle(), quad_weights.unsafe_data_handle());
  } else if (quad_type == "G") {
    gauss_legendre_quadrature_(
        nhstreams, mu_values.unsafe_data_handle(), quad_weights.unsafe_data_handle());
  } else  //if( quad_type=="L" )
  {
    lobatto_quadrature_(
        nhstreams, mu_values.unsafe_data_handle(), quad_weights.unsafe_data_handle());
  }

  // Set "extra" angle (at 0deg) if quad_type!="L" && !add_straight_angles
  if (nhza > 0) mu_values[nhstreams] = 1.;

  // FIXME: we should be able to avoid setting za_grid here in one way,
  // and resetting in another before leaving the WSM. This, however, requires
  // rearranging the angle order and angle assignment in the RT4-SSP prep
  // routines.
  za_grid.resize(2 * nummu);
  for (Index imu = 0; imu < nummu; imu++) {
    za_grid[imu] = acos(mu_values[imu]) * RAD2DEG;
    za_grid[nummu + imu] = 180. - za_grid[imu];
  }
  aa_grid.resize(1);
  aa_grid[0] = 0.;
}

void get_rt4surf_props(  // Output
    Vector& ground_albedo,
    Tensor3& ground_reflec,
    ComplexVector& ground_index,
    // Input
    ConstVectorView f_grid,
    const String& ground_type,
    const Numeric& surface_skin_t,
    ConstVectorView surface_scalar_reflectivity,
    ConstTensor3View surface_reflectivity,
    const GriddedField3& surface_complex_refr_index) {
  ARTS_USER_ERROR_IF (surface_skin_t < 0. || surface_skin_t > 1000.,
    "Surface temperature is set to ", surface_skin_t, " K,\n"
    "which is not considered a meaningful value.\n")

  const Index nf = f_grid.size();

  if (ground_type == "L")  // RT4's proprietary Lambertian
  {
    ARTS_USER_ERROR_IF (min(surface_scalar_reflectivity) < 0 ||
                        max(surface_scalar_reflectivity) > 1,
          "All values in *surface_scalar_reflectivity* must be inside [0,1].")

    // surface albedo
    if (surface_scalar_reflectivity.size() == f_grid.size())
      ground_albedo = surface_scalar_reflectivity;
    else if (surface_scalar_reflectivity.size() == 1)
      ground_albedo += surface_scalar_reflectivity[0];
    else {
      ARTS_USER_ERROR_IF(true,
      "For Lambertian surface reflection, the number of elements in\n"
      "*surface_scalar_reflectivity* needs to match the length of\n"
      "*f_grid* or be 1."
      "\n length of *f_grid* : ", f_grid.size(),
      "\n length of *surface_scalar_reflectivity* : ",
      surface_scalar_reflectivity.size(), "\n")
    }

  } else if (ground_type == "S")  // RT4's 'proprietary' Specular
  {
    const Index ref_sto = surface_reflectivity.nrows();

    chk_if_in_range("surface_reflectivity's 4", ref_sto, 1, 4);
    ARTS_USER_ERROR_IF (ref_sto != surface_reflectivity.ncols(),
      "The number of rows and columnss in *surface_reflectivity*\n"
      "must match each other.")

    ARTS_USER_ERROR_IF (min(surface_reflectivity(joker, 0, 0)) < 0 ||
                        max(surface_reflectivity(joker, 0, 0)) > 1,
          "All r11 values in *surface_reflectivity* must be inside [0,1].");

    // surface reflectivity
    if (surface_reflectivity.npages() == f_grid.size())
      if (ref_sto < 4)
        ground_reflec(joker, Range(0, ref_sto), Range(0, ref_sto)) =
            surface_reflectivity;
      else
        ground_reflec = surface_reflectivity(
            joker, Range(0, 4), Range(0, 4));
    else if (surface_reflectivity.npages() == 1)
      if (ref_sto < 4)
        for (Index f_index = 0; f_index < nf; f_index++)
          ground_reflec(f_index, Range(0, ref_sto), Range(0, ref_sto)) +=
              surface_reflectivity(0, joker, joker);
      else
        for (Index f_index = 0; f_index < nf; f_index++)
          ground_reflec(f_index, joker, joker) += surface_reflectivity(
              0, Range(0, 4), Range(0, 4));
    else {
      ARTS_USER_ERROR (
        "For specular surface reflection, the number of elements in\n"
        "*surface_reflectivity* needs to match the length of\n"
        "*f_grid* or be 1."
        "\n length of *f_grid* : ", f_grid.size(),
        "\n length of *surface_reflectivity* : ",
        surface_reflectivity.npages(), "\n")
    }

  } else if (ground_type == "F")  // RT4's proprietary Fresnel
  {
    //extract/interpolate from GriddedField
    Matrix n_real(nf, 1), n_imag(nf, 1);
    complex_n_interp(n_real,
                     n_imag,
                     surface_complex_refr_index,
                     "surface_complex_refr_index",
                     f_grid,
                     Vector(1, surface_skin_t));
    //ground_index = Complex(n_real(joker,0),n_imag(joker,0));
    for (Index f_index = 0; f_index < nf; f_index++) {
      ground_index[f_index] = Complex(n_real(f_index, 0), n_imag(f_index, 0));
    }
  } else {
    ARTS_USER_ERROR ( "Unknown surface type.\n")
  }
}

void za_grid_adjust(  // Output
    Vector& za_grid,
    // Input
    ConstVectorView mu_values,
    const Index& nummu) {
  for (Index j = 0; j < nummu; j++) {
    za_grid[nummu - 1 - j] = acos(mu_values[j]) * RAD2DEG;
    za_grid[nummu + j] = 180. - acos(mu_values[j]) * RAD2DEG;
    //cout << "setting scat_za[" << nummu-1-j << "]=" << za_grid[nummu-1-j]
    //     << " and [" << nummu+j << "]=" << za_grid[nummu+j]
    //     << " from mu[" << j << "]=" << mu_values[j] << "\n";
  }
}

void par_optpropCalc(Tensor5View emis_vector,
                     Tensor6View extinct_matrix,
                     //VectorView scatlayers,
                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                     const Vector& za_grid,
                     const Index& f_index,
                     ConstMatrixView pnd_profiles,
                     ConstVectorView t_profile,
                     const ArrayOfIndex& cloudbox_limits) {
  // Initialization
  extinct_matrix = 0.;
  emis_vector = 0.;

  const Index Np_cloud = pnd_profiles.ncols();
  const Index nummu = za_grid.size() / 2;

  ARTS_ASSERT(emis_vector.nbooks() == Np_cloud - 1);
  ARTS_ASSERT(extinct_matrix.nshelves() == Np_cloud - 1);

  // preparing input data
  Vector T_array{t_profile[Range(cloudbox_limits[0], Np_cloud)]};
  Matrix dir_array(za_grid.size(), 2, 0.);
  dir_array(joker, 0) = za_grid;

  // making output containers
  ArrayOfArrayOfTensor5 ext_mat_Nse;
  ArrayOfArrayOfTensor4 abs_vec_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;
  ArrayOfTensor5 ext_mat_ssbulk;
  ArrayOfTensor4 abs_vec_ssbulk;
  ArrayOfIndex ptype_ssbulk;
  Tensor5 ext_mat_bulk;
  Tensor4 abs_vec_bulk;
  Index ptype_bulk;

  opt_prop_NScatElems(ext_mat_Nse,
                      abs_vec_Nse,
                      ptypes_Nse,
                      t_ok,
                      scat_data,
                      T_array,
                      dir_array,
                      f_index);
  opt_prop_ScatSpecBulk(ext_mat_ssbulk,
                        abs_vec_ssbulk,
                        ptype_ssbulk,
                        ext_mat_Nse,
                        abs_vec_Nse,
                        ptypes_Nse,
                        pnd_profiles(joker, joker),
                        t_ok);
  opt_prop_Bulk(ext_mat_bulk,
                abs_vec_bulk,
                ptype_bulk,
                ext_mat_ssbulk,
                abs_vec_ssbulk,
                ptype_ssbulk);

  // Calculate layer averaged extinction and absorption and sort into RT4-format
  // data tensors
  for (Index ipc = 0; ipc < Np_cloud - 1; ipc++) {
    for (Index fi = 0; fi < abs_vec_bulk.nbooks(); fi++) {
      for (Index imu = 0; imu < nummu; imu++) {
        for (Index ist1 = 0; ist1 < 4; ist1++) {
          for (Index ist2 = 0; ist2 < 4; ist2++) {
            extinct_matrix(fi, ipc, 0, imu, ist2, ist1) =
                .5 * (ext_mat_bulk(fi, ipc, imu, ist1, ist2) +
                      ext_mat_bulk(fi, ipc + 1, imu, ist1, ist2));
            extinct_matrix(fi, ipc, 1, imu, ist2, ist1) =
                .5 * (ext_mat_bulk(fi, ipc, nummu + imu, ist1, ist2) +
                      ext_mat_bulk(fi, ipc + 1, nummu + imu, ist1, ist2));
          }
          emis_vector(fi, ipc, 0, imu, ist1) =
              .5 * (abs_vec_bulk(fi, ipc, imu, ist1) +
                    abs_vec_bulk(fi, ipc + 1, imu, ist1));
          emis_vector(fi, ipc, 1, imu, ist1) =
              .5 * (abs_vec_bulk(fi, ipc, nummu + imu, ist1) +
                    abs_vec_bulk(fi, ipc + 1, nummu + imu, ist1));
        }
      }
    }
  }
}

void sca_optpropCalc(  //Output
    Tensor6View scatter_matrix,
    Index& pfct_failed,
    //Input
    ConstTensor4View emis_vector,
    ConstTensor5View extinct_matrix,
    const Index& f_index,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    ConstMatrixView pnd_profiles,
    
    const Vector& za_grid,
    ConstVectorView quad_weights,
    const String& pfct_method,
    const Index& pfct_aa_grid_size,
    const Numeric& pfct_threshold,
    const Index& auto_inc_nstreams) {
  // FIXME: this whole funtions needs revision/optimization regarding
  // - temperature dependence (using new-type scat_data)
  // - using redundancies in sca_mat data at least for totally random
  //   orientation particles (are there some for az. random, too? like
  //   upper/lower hemisphere equiv?) - we might at least have a flag
  //   (transported down from calling function?) whether we deal exclusively
  //   with totally random orient particles (and then take some shortcuts...)

  // FIXME: do we have numerical issues, too, here in case of tiny pnd? check
  // with Patrick's Disort-issue case.

  // Initialization
  scatter_matrix = 0.;

  const Index N_se = pnd_profiles.nrows();
  const Index Np_cloud = pnd_profiles.ncols();
  const Index nza_rt = za_grid.size();

  ARTS_ASSERT(scatter_matrix.nvitrines() == Np_cloud - 1);

  // Check that total number of scattering elements in scat_data and pnd_profiles
  // agree.
  // FIXME: this should be done earlier. in calling method. outside freq- and
  // other possible loops. rather assert than runtime error.
  ARTS_USER_ERROR_IF (TotalNumberOfElements(scat_data) != N_se,
    "Total number of scattering elements in scat_data (",
    TotalNumberOfElements(scat_data), ") and pnd_profiles (", N_se,
    ") disagree.")

  ARTS_USER_ERROR_IF (pfct_aa_grid_size < 2,
    "Azimuth grid size for scatt matrix extraction"
    " (*pfct_aa_grid_size*) must be >1.\n"
    "Yours is ", pfct_aa_grid_size, ".\n")
  
  Vector aa_grid;
  nlinspace(aa_grid, 0, 180, pfct_aa_grid_size);

  Index i_se_flat = 0;
  Tensor5 sca_mat(N_se, nza_rt, nza_rt, 4, 4, 0.);
  Matrix ext_fixT_spt(N_se, nza_rt, 0.), abs_fixT_spt(N_se, nza_rt, 0.);

  // Precalculate azimuth integration weights for totally randomly oriented
  // (they are only determined by pfct_aa_grid_size)
  Numeric daa_totrand =
      1. / float(pfct_aa_grid_size - 1);  // 2*180./360./(pfct_aa_grid_size-1)

  // first we extract Z at one T, integrate the azimuth data at each
  // za_inc/za_sca combi (to get the Fourier series 0.th mode), then interpolate
  // to the mu/mu' combis we need in RT.
  //
  // FIXME: are the stokes component assignments correct? for totally random?
  // and azimuthally random? (as they are done differently...)
  for (Size i_ss = 0; i_ss < scat_data.size(); i_ss++) {
    for (Size i_se = 0; i_se < scat_data[i_ss].size(); i_se++) {
      SingleScatteringData ssd = scat_data[i_ss][i_se];
      Index this_f_index = (ssd.pha_mat_data.nlibraries() == 1 ? 0 : f_index);
      Index i_pfct;
      if (pfct_method == "low")
        i_pfct = 0;
      else if (pfct_method == "high")
        i_pfct = ssd.T_grid.size() - 1;
      else if( pfct_method=="median" )
        i_pfct = ssd.T_grid.size() / 2;
      else {
        ARTS_USER_ERROR_IF(true,
          "*pfct_method* must be \"low\", \"high\" or \"median\"." );
      }

      if (ssd.ptype == PTYPE_TOTAL_RND) {
        Matrix pha_mat(4, 4, 0.);
        for (Index iza = 0; iza < nza_rt; iza++) {
          for (Index sza = 0; sza < nza_rt; sza++) {
            Matrix pha_mat_int(4, 4, 0.);
            for (Index saa = 0; saa < pfct_aa_grid_size; saa++) {
              pha_matTransform(
                  pha_mat(joker, joker),
                  ssd.pha_mat_data(
                      this_f_index, i_pfct, joker, joker, joker, joker, joker),
                  ssd.za_grid,
                  ssd.aa_grid,
                  ssd.ptype,
                  sza,
                  saa,
                  iza,
                  0,
                  za_grid,
                  aa_grid);

              if (saa == 0 || saa == pfct_aa_grid_size - 1)
                pha_mat *= (daa_totrand / 2.);
              else
                pha_mat *= daa_totrand;
              pha_mat_int += pha_mat;
            }
            sca_mat(i_se_flat, iza, sza, joker, joker) = pha_mat_int;
          }
          ext_fixT_spt(i_se_flat, iza) =
              ssd.ext_mat_data(this_f_index, i_pfct, 0, 0, 0);
          abs_fixT_spt(i_se_flat, iza) =
              ssd.abs_vec_data(this_f_index, i_pfct, 0, 0, 0);
        }
      } else if (ssd.ptype == PTYPE_AZIMUTH_RND) {
        Index nza_se = ssd.za_grid.size();
        Index naa_se = ssd.aa_grid.size();
        Tensor4 pha_mat_int(nza_se, nza_se, 4, 4, 0.);
        ConstVectorView za_datagrid = ssd.za_grid;
        ConstVectorView aa_datagrid = ssd.aa_grid;
        ARTS_ASSERT(aa_datagrid[0] == 0.);
        ARTS_ASSERT(aa_datagrid[naa_se - 1] == 180.);
        Vector daa(naa_se);

        // Precalculate azimuth integration weights for this azimuthally
        // randomly oriented scat element
        // (need to do this per scat element as ssd.aa_grid is scat
        //  element specific (might change between elements) and need to do
        //  this on actual grid instead of grid number since the grid can,
        //  at least theoretically be non-equidistant)
        daa[0] = (aa_datagrid[1] - aa_datagrid[0]) / 360.;
        for (Index saa = 1; saa < naa_se - 1; saa++)
          daa[saa] = (aa_datagrid[saa + 1] - aa_datagrid[saa - 1]) / 360.;
        daa[naa_se - 1] =
            (aa_datagrid[naa_se - 1] - aa_datagrid[naa_se - 2]) / 360.;

        // first, extracting the phase matrix at the scatt elements own
        // polar angle grid, deriving their respective azimuthal (Fourier
        // series) 0-mode
        for (Index iza = 0; iza < nza_se; iza++)
          for (Index sza = 0; sza < nza_se; sza++) {
            for (Index saa = 0; saa < naa_se; saa++) {
              for (Index ist1 = 0; ist1 < 4; ist1++)
                for (Index ist2 = 0; ist2 < 4; ist2++)
                  pha_mat_int(sza, iza, ist1, ist2) +=
                      daa[saa] * ssd.pha_mat_data(this_f_index,
                                                  i_pfct,
                                                  sza,
                                                  saa,
                                                  iza,
                                                  0,
                                                  ist1 * 4 + ist2);
            }
          }

        // second, interpolating the extracted azimuthal mode to the RT4
        // solver polar angles
        for (Index iza = 0; iza < nza_rt; iza++) {
          for (Index sza = 0; sza < nza_rt; sza++) {
            GridPos za_sca_gp;
            GridPos za_inc_gp;
            Vector itw(4);
            Matrix pha_mat_lab(4, 4, 0.);
            Numeric za_sca = za_grid[sza];
            Numeric za_inc = za_grid[iza];

            gridpos(za_inc_gp, za_datagrid, za_inc);
            gridpos(za_sca_gp, za_datagrid, za_sca);
            interpweights(itw, za_sca_gp, za_inc_gp);

            for (Index ist1 = 0; ist1 < 4; ist1++)
              for (Index ist2 = 0; ist2 < 4; ist2++) {
                pha_mat_lab(ist1, ist2) =
                    interp(itw,
                           pha_mat_int(Range(joker), Range(joker), ist1, ist2),
                           za_sca_gp,
                           za_inc_gp);
                //if (ist1+ist2==1)
                //  pha_mat_lab(ist1,ist2) *= -1.;
              }

            sca_mat(i_se_flat, iza, sza, joker, joker) = pha_mat_lab;
          }
          ext_fixT_spt(i_se_flat, iza) =
              ssd.ext_mat_data(this_f_index, i_pfct, iza, 0, 0);
          abs_fixT_spt(i_se_flat, iza) =
              ssd.abs_vec_data(this_f_index, i_pfct, iza, 0, 0);
        }
      } else {
        ARTS_USER_ERROR ( "Unsuitable particle type encountered.")
      }
      i_se_flat++;
    }
  }

  ARTS_ASSERT(i_se_flat == N_se);
  // now we sum up the Z(mu,mu') over the scattering elements weighted by the
  // pnd_field data, deriving Z(z,mu,mu') and sorting this into
  // scatter_matrix
  // FIXME: it seems that at least for totally random, corresponding upward and
  // downward directions have exactly the same (0,0) elements. Can we use this
  // to reduce calc efforts (for sca and its normalization as well as for ext
  // and abs?)
  Index nummu = nza_rt / 2;
  for (Index scat_p_index_local = 0; scat_p_index_local < Np_cloud - 1;
       scat_p_index_local++) {
    Vector ext_fixT(nza_rt, 0.), abs_fixT(nza_rt, 0.);

    for (Index i_se = 0; i_se < N_se; i_se++) {
      Numeric pnd_mean = 0.5 * (pnd_profiles(i_se, scat_p_index_local + 1) +
                                pnd_profiles(i_se, scat_p_index_local));
      if (pnd_mean != 0.)
        for (Index iza = 0; iza < nummu; iza++) {
          // JM171003: not clear to me anymore why this check. if pnd_mean
          // is non-zero, then extinction should also be non-zero by
          // default?
          //if ( (extinct_matrix(scat_p_index_local,0,iza,0,0)+
          //      extinct_matrix(scat_p_index_local,1,iza,0,0)) > 0. )
          for (Index sza = 0; sza < nummu; sza++) {
            for (Index ist1 = 0; ist1 < 4; ist1++)
              for (Index ist2 = 0; ist2 < 4; ist2++) {
                // we can't use 4 jokers here since '*' doesn't
                // exist for Num*MatView. Also, order of stokes matrix
                // dimensions is inverted here (aka scat matrix is
                // transposed).
                scatter_matrix(scat_p_index_local, 0, iza, ist2, sza, ist1) +=
                    pnd_mean * sca_mat(i_se, iza, sza, ist1, ist2);
                scatter_matrix(scat_p_index_local, 1, iza, ist2, sza, ist1) +=
                    pnd_mean * sca_mat(i_se, nummu + iza, sza, ist1, ist2);
                scatter_matrix(scat_p_index_local, 2, iza, ist2, sza, ist1) +=
                    pnd_mean * sca_mat(i_se, iza, nummu + sza, ist1, ist2);
                scatter_matrix(scat_p_index_local, 3, iza, ist2, sza, ist1) +=
                    pnd_mean *
                    sca_mat(i_se, nummu + iza, nummu + sza, ist1, ist2);
              }
          }

          ext_fixT[iza] += pnd_mean * ext_fixT_spt(i_se, iza);
          abs_fixT[iza] += pnd_mean * abs_fixT_spt(i_se, iza);
        }
    }

    for (Index iza = 0; iza < nummu; iza++) {
      for (Index ih = 0; ih < 2; ih++) {
        if (extinct_matrix(scat_p_index_local, ih, iza, 0, 0) > 0.) {
          Numeric sca_mat_integ = 0.;

          // We need to calculate the nominal values for the fixed T, we
          // used above in the pha_mat extraction. Only this tells us whether
          // angular resolution is sufficient.
          //
          //Numeric ext_nom = extinct_matrix(scat_p_index_local,ih,iza,0,0);
          //Numeric sca_nom = ext_nom-emis_vector(scat_p_index_local,ih,iza,0);
          Numeric ext_nom = ext_fixT[iza];
          Numeric sca_nom = ext_nom - abs_fixT[iza];
          Numeric w0_nom = sca_nom / ext_nom;
          ARTS_ASSERT(w0_nom >= 0.);

          for (Index sza = 0; sza < nummu; sza++) {
            //                SUM2 = SUM2 + 2.0D0*PI*QUAD_WEIGHTS(J2)*
            //     .                 ( SCATTER_MATRIX(1,J2,1,J1, L,TSL)
            //     .                 + SCATTER_MATRIX(1,J2,1,J1, L+2,TSL) )
            sca_mat_integ +=
                quad_weights[sza] *
                (scatter_matrix(scat_p_index_local, ih, iza, 0, sza, 0) +
                 scatter_matrix(scat_p_index_local, ih + 2, iza, 0, sza, 0));
          }

          // compare integrated scatt matrix with ext-abs for respective
          // incident polar angle - consistently with scat_dataCheck, we do
          // this in terms of albedo deviation (since PFCT deviations at
          // small albedos matter less than those at high albedos)
          //            SUM1 = EMIS_VECTOR(1,J1,L,TSL)-EXTINCT_MATRIX(1,1,J1,L,TSL)
          Numeric w0_act = 2. * PI * sca_mat_integ / ext_nom;
          Numeric pfct_norm = 2. * PI * sca_mat_integ / sca_nom;
          /*
              cout << "sca_mat norm deviates " << 1e2*abs(1.-pfct_norm) << "%"
                   << " (" << abs(w0_act-w0_nom) << " in albedo).\n";
              */
          Numeric sca_nom_paropt =
              extinct_matrix(scat_p_index_local, ih, iza, 0, 0) -
              emis_vector(scat_p_index_local, ih, iza, 0);
          //Numeric w0_nom_paropt = sca_nom_paropt /
          //  extinct_matrix(scat_p_index_local,ih,iza,0,0);
          /*
              cout << "scat_p=" << scat_p_index_local
                   << ", iza=" << iza << ", hem=" << ih << "\n";
              cout << "  scaopt (paropt) w0_act= " << w0_act
                   << ", w0_nom = " << w0_nom
                   << " (" << w0_nom_paropt
                   << "), diff=" << w0_act-w0_nom
                   << " (" << w0_nom-w0_nom_paropt << ").\n";
              */

          if (abs(w0_act - w0_nom) > pfct_threshold) {
            if (pfct_failed >= 0) {
              if (auto_inc_nstreams) {
                pfct_failed = 1;
                return;
              } else {
                ARTS_USER_ERROR_IF(true,
                  "Bulk scattering matrix normalization deviates significantly\n"
                  "from expected value (", 1e2 * abs(1. - pfct_norm),
                  "%,"
                  " resulting in albedo deviation of ",
                  abs(w0_act - w0_nom), ").\n"
                  "Something seems wrong with your scattering data "
                  " (did you run *scat_dataCheck*?)\n"
                  "or your RT4 setup (try increasing *nstreams* and in case"
                  " of randomly oriented particles possibly also"
                  " pfct_aa_grid_size).");
              }
            }
          } else if (abs(w0_act - w0_nom) > pfct_threshold * 0.1 ||
                     abs(1. - pfct_norm) > 1e-2) {
          }
          // rescale scattering matrix to expected (0,0) value (and scale all
          // other elements accordingly)
          //
          // However, here we should not use the pfct_norm based on the
          // deviation from the fixed-temperature ext and abs. Instead, for
          // energy conservation reasons, this needs to be consistent with
          // extinct_matrix and emis_vector. Hence, recalculate pfct_norm
          // from them.
          //cout << "  scaopt (paropt) pfct_norm dev = " << 1e2*pfct_norm-1e2;
          pfct_norm = 2. * PI * sca_mat_integ / sca_nom_paropt;
          //cout << " (" << 1e2*pfct_norm-1e2 << ")%.\n";
          //
          // FIXME: not fully clear whether applying the same rescaling
          // factor to all stokes elements is the correct way to do. check
          // out, e.g., Vasilieva (JQSRT 2006) for better approaches.
          // FIXME: rather rescale Z(0,0) only as we should have less issues
          // for other Z elements as they are less peaked/featured.
          scatter_matrix(scat_p_index_local, ih, iza, joker, joker, joker) /=
              pfct_norm;
          scatter_matrix(
              scat_p_index_local, ih + 2, iza, joker, joker, joker) /=
              pfct_norm;
          //if (scat_p_index_local==49)
        }
      }
    }
  }
}

void rt4_test(Tensor4& out_rad,
              const String& datapath) {
  //emissivity.resize(4);
  //reflectivity.resize(4);

  Index nstokes = 2;
  Index nummu = 8;
  Index nuummu = 0;
  Numeric max_delta_tau = 1.0E-6;
  String quad_type = "L";
  Numeric ground_temp = 300.;
  String ground_type = "L";
  Numeric ground_albedo = 0.05;
  Complex ground_index;
  Numeric sky_temp = 0.;
  Numeric wavelength = 880.;
  //Index noutlevels=1;
  //ArrayOfIndex outlevels(1);
  //outlevels[0]=1;

  Vector height, temperatures, gas_extinct;
  Tensor5 sca_data;
  Tensor4 ext_data;
  Tensor3 abs_data;
  ReadXML(height, datapath + "z.xml");
  ReadXML(temperatures, datapath + "T.xml");
  ReadXML(gas_extinct, datapath + "abs_gas.xml");
  ReadXML(abs_data, datapath + "abs_par.xml");
  ReadXML(ext_data, datapath + "ext_par.xml");
  ReadXML(sca_data, datapath + "sca_par.xml");
  Index num_layers = height.size() - 1;
  Index num_scatlayers = 3;
  Vector scatlayers(num_layers, 0.);
  scatlayers[3] = 1.;
  scatlayers[4] = 2.;
  scatlayers[5] = 3.;

  // the read in sca/ext/abs_data is the complete set (and it's in the wrong
  // order for passing it directly to radtrano). before handing over to
  // fortran, we need to reduce it to the number of stokes elements to be
  // used. we can't use views here as all data needs to be continuous in
  // memory; that is, we have to explicitly copy the data we need.
  Tensor6 scatter_matrix(num_scatlayers, 4, nummu, nstokes, nummu, nstokes);
  for (Index ii = 0; ii < 4; ii++)
    for (Index ij = 0; ij < nummu; ij++)
      for (Index ik = 0; ik < nstokes; ik++)
        for (Index il = 0; il < nummu; il++)
          for (Index im = 0; im < nstokes; im++)
            for (Index in = 0; in < num_scatlayers; in++)
              scatter_matrix(in, ii, ij, ik, il, im) =
                  sca_data(im, il, ik, ij, ii);
  Tensor5 extinct_matrix(num_scatlayers, 2, nummu, nstokes, nstokes);
  for (Index ii = 0; ii < 2; ii++)
    for (Index ij = 0; ij < nummu; ij++)
      for (Index ik = 0; ik < nstokes; ik++)
        for (Index il = 0; il < nstokes; il++)
          for (Index im = 0; im < num_scatlayers; im++)
            extinct_matrix(im, ii, ij, ik, il) = ext_data(il, ik, ij, ii);
  Tensor4 emis_vector(num_scatlayers, 2, nummu, nstokes);
  for (Index ii = 0; ii < 2; ii++)
    for (Index ij = 0; ij < nummu; ij++)
      for (Index ik = 0; ik < nstokes; ik++)
        for (Index il = 0; il < num_scatlayers; il++)
          emis_vector(il, ii, ij, ik) = abs_data(ik, ij, ii);

  // dummy parameters necessary due to modified, flexible surface handling
  Tensor4 surf_refl_mat(nummu, nstokes, nummu, nstokes, 0.);
  Matrix surf_emis_vec(nummu, nstokes, 0.);
  Matrix ground_reflec(nstokes, nstokes, 0.);

  // Output variables
  Vector mu_values(nummu);
  Tensor3 up_rad(num_layers + 1, nummu, nstokes, 0.);
  Tensor3 down_rad(num_layers + 1, nummu, nstokes, 0.);

  radtrano_(nstokes,
            nummu,
            nuummu,
            max_delta_tau,
            quad_type.c_str(),
            ground_temp,
            ground_type.c_str(),
            ground_albedo,
            ground_index,
            ground_reflec.unsafe_data_handle(),
            surf_refl_mat.unsafe_data_handle(),
            surf_emis_vec.unsafe_data_handle(),
            sky_temp,
            wavelength,
            num_layers,
            height.unsafe_data_handle(),
            temperatures.unsafe_data_handle(),
            gas_extinct.unsafe_data_handle(),
            num_scatlayers,
            scatlayers.unsafe_data_handle(),
            extinct_matrix.unsafe_data_handle(),
            emis_vector.unsafe_data_handle(),
            scatter_matrix.unsafe_data_handle(),
            //noutlevels,
            //outlevels.unsafe_data_handle(),
            mu_values.unsafe_data_handle(),
            up_rad.unsafe_data_handle(),
            down_rad.unsafe_data_handle());

  //so far, output is in
  //    units W/m^2 um sr
  //    dimensions [numlayers+1,nummu,nstokes]
  //    sorted from high to low (altitudes) and 0 to |1| (mu)
  //WriteXML( "ascii", up_rad, "up_rad.xml", 0, "up_rad", "", "" );
  //WriteXML( "ascii", down_rad, "down_rad.xml", 0, "down_rad", "", "" );

  //to be able to compare with RT4 reference results, reshape output into
  //RT4-output type table (specifically, resort up_rad such that it runs from
  //zenith welling to horizontal, thus forms a continuous angle grid with
  //down_rad. if later changing up_rad/down_rad sorting such that it is in
  //line with cloudbox_field, then this has to be adapted as well...
  out_rad.resize(num_layers + 1, 2, nummu, nstokes);
  for (Index ii = 0; ii < nummu; ii++)
    out_rad(joker, 0, ii, joker) = up_rad(joker, nummu - 1 - ii, joker);
  //out_rad(joker,0,joker,joker) = up_rad;
  out_rad(joker, 1, joker, joker) = down_rad;
  //WriteXML( "ascii", out_rad, "out_rad.xml", 0, "out_rad", "", "" );
}

#endif /* ENABLE_RT4 */
