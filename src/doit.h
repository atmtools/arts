/*===========================================================================
  ===  File description 
  ===========================================================================*/

/*!
  \file   doit.h
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   2003-06-03
  
  \brief  Radiative transfer in cloudbox.
  
  This file contains functions related to the radiative transfer in the 
  cloudbox using the DOIT method.
*/

#ifndef doit_h
#define doit_h

#include <workspace.h>

//! Solves monochromatic VRTE for an atmospheric slab with constant conditions.
/*!
    The function can be used for cloudbox calculations.

    The function is best explained by considering a homogenous layer. That is,
    the physical conditions inside the layer are constant. In reality they
    are not constant, so in practical all coefficients have to be averaged
    before calling this function.
    Total extinction and absorption inside the layer are described by
    *ext_mat_av* and *abs_vec_av* respectively,
    the blackbdody radiation of the layer is given by *rte_planck_value*
    and the propagation path length through the layer is *lstep*.

    There is an additional scattering source term in the
    VRTE, the scattering integral term. For this function a constant
    scattering term is assumed. The radiative transfer step is only a part
    the iterative solution of the scattering problem, for more
    information consider AUG. In the clearsky case this variable has to be
    set to 0.

    When calling the function, the vector *stokes_vec* shall contain the
    Stokes vector for the incoming radiation. The function returns this
    vector, then containing the outgoing radiation on the other side of the
    layer.

    The function performs the calculations differently depending on the
    conditions to improve the speed. There are three cases: <br>
       1. The matrix ext_mat_gas is diagonal (unpolarised absorption). <br>
       2. The total general case.

    \param[in,out]  stokes_vec A Stokes vector.
    \param[out]     trans_mat Transmission matrix of slab.
    \param[in]      ext_mat_av Averaged extinction matrix.
    \param[in]      abs_vec_av Averaged absorption vector.
    \param[in]      sca_vec_av averaged scattering vector.
    \param[in]      lstep The length of the RTE step.
    \param[in]      rtp_planck_value Blackbody radiation.
    \param[in]      trans_is_precalc FIXMEDOC

    \author Richard Larsson,
    \date   2017-08-14
*/
void rte_step_doit_replacement(  //Output and Input:
    Stokvec& stokes_vec,
    Muelmat& trans_mat,
    //Input
    const Propmat& ext_mat_av,
    const Stokvec& abs_vec_av,
    const Stokvec& sca_vec_av,
    const Numeric& lstep,
    const Numeric& rtp_planck_value,
    const bool& trans_is_precalc=false);

//! Calculate ext_mat, abs_vec for all points inside the cloudbox for one
//  propagation direction.
/*!
  sca_vec can be obtained from the workspace variable doit_scat_field.
  As we need the average for each layer, it makes sense to calculte
  the coefficients once and store them in an array instead of
  calculating at each point the coefficient of the point above and
  the point below.

  \param[out]   ws Current Workspace
  \param[out]   ext_mat_field extinction matrix field
  \param[out]   abs_vec_field absorption vector field
  \param[in]    spt_calc_agenda Agenda for calculation of single scattering properties
  \param[in]    za_index Indices for
  \param[in]    aa_index    propagation direction
  \param[in]    cloudbox_limits Cloudbox limits.
  \param[in]    t_field Temperature field
  \param[in]    pnd_field Particle number density field.

  \author Claudia Emde
  \date 2002-06-03
*/
void cloud_fieldsCalc(const Workspace& ws,
                      // Output:
                      Tensor5View ext_mat_field,
                      Tensor4View abs_vec_field,
                      // Input:
                      const Agenda& spt_calc_agenda,
                      const Index& za_index,
                      const Index& aa_index,
                      const ArrayOfIndex& cloudbox_limits,
                      ConstTensor3View t_field,
                      ConstTensor4View pnd_field);

//! Convergence acceleration
/*!
 This function accelarate the convergence of the doit iteration by extrapolation
 of the doit_i_mono from the previous three iteration steps.
 The acceleration method is called Ng-Acceleration and was develop by Ng (1974)

 \param[out]    cloudbox_field_mono Radiation field in cloudbox
 \param[in]     acceleration_input Array of the previous three iteration steps
 \param[in]     accelerated Index wether to accelerate only the intensity or the
                whole Stokes Vector

 \author Jakob Doerr
 \date 2016-04-26

*/
void cloudbox_field_ngAcceleration(  //Output
    Tensor6& cloudbox_field_mono,
    //Input
    const ArrayOfTensor6& acceleration_input,
    const Index& accelerated);

//! Optimize the zenith angle grid
/*!
  This method optimizes the zenith angle grid. For optimization it uses the
  interpolation method given by *scat_za_interp* (0 - linear interpolation,
  1 - polynomial interpolation).
  As input it needs the intensity field calculated on a very fine zenith angle
  grid (*za_grid_fine*). The function picks out as many grid points as required
  to achieve the required accuracy (*acc* [%]). This methods optimizes only
  the intensity (first Stokes component) for 1D cases (first latitude and
  longitude of the intensity field.
  (Could be modified to optimize all Stokes components at the same time, if we
  don't want to use the clearsky field for grid optimuzation.)

  \param[out]   za_grid_opt Optimized zenith angle grid.
  \param[out]   cloudbox_field_opt Optimized intensity field.
  \param[in]    za_grid_fine Fine zenith angle grid.
  \param[in]    cloudbox_field_mono Radiation field calculated on a very fine za grid.
  \param[in]    acc Accuracy of optimization [%].
  \param[in]    scat_za_interp Interpolation method.

  \author Claudia Emde
  \date 2004-04-05
*/
void za_gridOpt(  //Output:
    Vector& za_grid_opt,
    Matrix& i_field_opt,
    // Input
    ConstVectorView za_grid_fine,
    ConstTensor6View i_field,
    const Numeric& acc,
    const Index& scat_za_interp);

//! Normalization of scattered field
/*!
  Calculate the scattered extinction field and apply the
  derived correction factor to doit_scat_field.

  Only 1D is supported.

  \param[in,out] ws Current workspace
  \param[in,out] doit_scat_field Scattered field
  \param[in]     cloudbox_field_mono Radiation field
  \param[in]     cloudbox_limits Cloudbox limits
  \param[in]     spt_calc_agenda Calculates single scattering properties
  \param[in]     atmopshere_dim Dimension of atmosphere
  \param[in]     za_grid Zenith angle grid
  \param[in]     aa_grid Azimuth angle grid
  \param[in]     pnd_field pnd field
  \param[in]     t_field Atmospheric temperature field
  \param[in]     norm_error_threshold  Normalization error threshold
  \param[in]     norm_debug Flag for normalization debug output

  \author Oliver Lemke
  \date 2013-01-17
*/
void doit_scat_fieldNormalize(const Workspace& ws,
                              Tensor6& doit_scat_field,
                              const Tensor6& cloudbox_field_mono,
                              const ArrayOfIndex& cloudbox_limits,
                              const Agenda& spt_calc_agenda,
                              const Vector& za_grid,
                              const Vector& aa_grid,
                              const Tensor4& pnd_field,
                              const Tensor3& t_field,
                              const Numeric& norm_error_threshold,
                              const Index& norm_debug);

#endif  //doit_h
