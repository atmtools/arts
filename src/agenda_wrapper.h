/*!
  \file   agenda_wrappers.h
  \author Simon Pfreundschuh <simonpf@chalmers.se>
  \date   Thu October 06 17:09:00 2016

  \brief Wrappers for the inversion_iterate_agenda in order to them as invlib
  forward models.
*/
#ifndef agenda_wrapper_h
#define agenda_wrapper_h

#include "oem.h"

//! Wrapper class for forward model.
/*!
  Wrapper class for the inversion_iterate_agendaExecute function to implement
  the forward model interface used by the non-linear oem function in oem.cc.
  The object is constructed with the pointers to the variables used as arguments
  for the function and then simply forwards the calls made to
  ForwardModel::evaluate() and ForwardModel::evaluate_jacobian() to
  inversion_iterate_agendaExecute.

 */
class AgendaWrapper {
 public:
  const unsigned int m, n;
  OEMMatrixReference jacobian;
  OEMVector yi;

  //! Create inversion_iterate_agendaExecute wrapper.
  /*!
  Initializes the wrapper object for the inversion_iterate_agendaExecute
  method. The object forwards the evaluate() and evaluate_jacobian() calls
  made by the iterative OEM methods to inversion_iterate_agendaExecute using
  the arguments provided to the constructor.

  \param ws_ Pointer to the workspace argument of the agenda execution function.
  function.
  \param inversion_iterate_agenda_ Pointer to the x argument of the agenda
  execution function.

*/
  AgendaWrapper(Workspace *ws_,
                unsigned int m_,
                unsigned int n_,
                Matrix &jacobian_,
                Vector &yi_,
                const Agenda *inversion_iterate_agenda_)
      : m(m_),
        n(n_),
        jacobian(jacobian_),
        yi(yi_),
        ws(ws_),
        inversion_iterate_agenda(inversion_iterate_agenda_),
        reuse_jacobian((jacobian_.nrows() != 0) && (jacobian_.ncols() != 0) &&
                       (yi_.nelem() != 0)),
        iteration_counter(0) {}

  AgendaWrapper(const AgendaWrapper &) = delete;
  AgendaWrapper(AgendaWrapper &&) = delete;
  AgendaWrapper &operator=(const AgendaWrapper &) = delete;
  AgendaWrapper &operator=(AgendaWrapper &&) = delete;

  //! Evaluate forward model and compute Jacobian.
  /*!

  Forwards the call to evaluate_jacobian() and evaluate() that is made by
  Gauss-Newton and Levenberg-Marquardt OEM methods using the variables pointed
  to by the pointers provided to the constructor as arguments.

  \param[out] y The measurement vector y = K(x) for the current state vector x
  as computed by the forward model.
  \param[out] J The Jacobian Ki=d/dx(K(x)) of the forward model.
  \param[in] x The current state vector x.
*/
  OEMMatrixReference Jacobian(const OEMVector &xi, OEMVector &yi_) {
    if (!reuse_jacobian) {
      inversion_iterate_agendaExecute(
          *ws, yi, jacobian, xi, 1, 0, *inversion_iterate_agenda);
      yi_ = yi;
      iteration_counter += 1;
    } else {
      reuse_jacobian = false;
      yi_ = yi;
    }
    return jacobian;
  }

  //! Evaluate forward model.
  /*!

  Forwards the call to evaluate that is made by Gauss-Newton and
  Levenberg-Marquardt OEM methods to the function pointers provided.

  \param[out] y The measurement vector y = K(x) for the current state vector x.
  \param[in] x The current state vector x.
*/
  OEMVector evaluate(const OEMVector &xi) {
    if (!reuse_jacobian) {
      Matrix dummy;
      inversion_iterate_agendaExecute(
          *ws, yi, dummy, xi, 0, iteration_counter, *inversion_iterate_agenda);
    } else {
      reuse_jacobian = false;
    }
    return yi;
  }

 private:
  Workspace *ws;
  const Agenda *inversion_iterate_agenda;
  bool reuse_jacobian;
  unsigned int iteration_counter;
};

#endif  // agenda_wrappers_h
