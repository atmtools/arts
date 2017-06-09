/*!
  \file   oem.h
  \author Simon Pfreundschuh <simonpf@chalmers.se>
  \date   Fri March 25 15:53:54 2016

  \brief Contains type aliases and interface classes for the
  invlib library.
*/

#ifndef oem_h
#define oem_h

#include <type_traits>

#include "invlib/algebra.h"
#include "invlib/map.h"
#include "invlib/optimization.h"
#include "invlib/algebra/solvers.h"
#include "invlib/algebra/precision_matrix.h"
#include "invlib/interfaces/arts_wrapper.h"
#include "invlib/profiling/timer.h"

////////////////////////////////////////////////////////////////////////////////
//  Type Aliases
////////////////////////////////////////////////////////////////////////////////

using OEMVector          = invlib::Vector<invlib::Timer<ArtsVector>>;
using OEMMatrix          = invlib::Matrix<ArtsMatrix>;
using OEMMatrixReference = invlib::Matrix<ArtsMatrixReference<Matrix>>;
using OEMSparse          = invlib::Matrix<invlib::Timer<ArtsMatrixReference<const Sparse>>>;

using Identity        = invlib::MatrixIdentity<OEMMatrix>;
using PrecisionMatrix = invlib::PrecisionMatrix<OEMMatrix>;
using PrecisionSparse = invlib::PrecisionMatrix<OEMSparse>;

// OEM types.

using invlib::Formulation;

// Standard Form.

template <typename ForwardModel>
using OEM_S_S = invlib::MAP<ForwardModel, OEMMatrix, Sparse,
    Sparse, OEMVector, Formulation::STANDARD>;

template <typename ForwardModel>
using OEM_PS_PS = invlib::MAP<ForwardModel, OEMMatrix, PrecisionSparse,
    PrecisionSparse, OEMVector, Formulation::STANDARD>;

template <typename ForwardModel>
using OEM_D_D = invlib::MAP<ForwardModel, OEMMatrix, OEMMatrix,
    OEMMatrix, OEMVector, Formulation::STANDARD>;

template <typename ForwardModel>
using OEM_PD_PD = invlib::MAP<ForwardModel, OEMMatrix, PrecisionMatrix,
    PrecisionMatrix, OEMVector, Formulation::STANDARD>;

// N-Form.

template <typename ForwardModel>
using OEM_NFORM_S_S = invlib::MAP<ForwardModel, OEMMatrix, Sparse,
    Sparse, OEMVector, Formulation::NFORM>;

template <typename ForwardModel>
using OEM_NFORM_PS_PS = invlib::MAP<ForwardModel, OEMMatrix, PrecisionSparse,
    PrecisionSparse, OEMVector, Formulation::NFORM>;

template <typename ForwardModel>
using OEM_NFORM_D_D = invlib::MAP<ForwardModel, OEMMatrix, OEMMatrix,
    OEMMatrix, OEMVector, Formulation::NFORM>;

template <typename ForwardModel>
using OEM_NFORM_PD_PD = invlib::MAP<ForwardModel, OEMMatrix, PrecisionMatrix,
    PrecisionMatrix, OEMVector, Formulation::NFORM>;

// M-Form.

template <typename ForwardModel>
using OEM_MFORM_S_S = invlib::MAP<ForwardModel, OEMMatrix, Sparse,
    Sparse, OEMVector, Formulation::MFORM>;

template <typename ForwardModel>
using OEM_MFORM_PS_PS = invlib::MAP<ForwardModel, OEMMatrix, PrecisionSparse,
    PrecisionSparse, OEMVector, Formulation::MFORM>;

template <typename ForwardModel>
using OEM_MFORM_D_D = invlib::MAP<ForwardModel, OEMMatrix, OEMMatrix,
    OEMMatrix, OEMVector, Formulation::MFORM>;

template <typename ForwardModel>
using OEM_MFORM_PD_PD = invlib::MAP<ForwardModel, OEMMatrix, PrecisionMatrix,
    PrecisionMatrix, OEMVector, Formulation::MFORM>;

// Solvers.

using Std     = invlib::Standard;
using CG      = invlib::ConjugateGradient<>;
template <typename TransformationMatrixType, typename SolverType>
class NormalizingSolver;

// Normed solvers. Currently uses dense matrices to store the normalization
// so this could be optimized.

template <typename SolverType = invlib::Standard>
using Normed  = NormalizingSolver<OEMMatrix, SolverType>;

// Optimization Methods.

using GN          = invlib::GaussNewton<Numeric, Normed<>>;
using GN_CG       = invlib::GaussNewton<Numeric, Normed<CG>>;
using LM_D        = invlib::LevenbergMarquardt<Numeric, OEMMatrix, Normed<>>;
using LM_CG_D     = invlib::LevenbergMarquardt<Numeric, OEMMatrix, Normed<CG>>;
using LM_Sparse_D = invlib::LevenbergMarquardt<Numeric, OEMMatrix, Normed<>>;
using LM_S        = invlib::LevenbergMarquardt<Numeric, OEMSparse, Normed<>>;
using LM_CG_S     = invlib::LevenbergMarquardt<Numeric, OEMSparse, Normed<CG>>;
using LM_I        = invlib::LevenbergMarquardt<Numeric, Identity,  Normed<>>;
using LM_CG_I     = invlib::LevenbergMarquardt<Numeric, Identity,  Normed<CG>>;
using LM_Sparse_S = invlib::LevenbergMarquardt<Numeric, OEMSparse, Normed<>>;

////////////////////////////////////////////////////////////////////////////////
//  Normalizing Solver
////////////////////////////////////////////////////////////////////////////////

/**
 * Invlib Solver class, that wraps around a given solver and transforms the linear
 * system from left and right with the given transformation matrix. This is used
 * to implement the normalization from qpack.
 */
template
<
    typename TransformationMatrixType,
    typename SolverType = invlib::Standard
>
class NormalizingSolver : SolverType
{
public:
    NormalizingSolver(const SolverType & s)
    : SolverType(s), apply(false), trans() {}

    template<typename ...Params>
    NormalizingSolver(const TransformationMatrixType &t, bool a, Params ...params)
    : SolverType(params...), apply(a), trans(t) {}

    template <typename MatrixType, typename VectorType>
    auto solve(const MatrixType & A, const VectorType & v)
    -> typename VectorType::ResultType
    {
        typename VectorType::ResultType w;
        if (apply) {
            typename VectorType::ResultType vv = trans * v;
            auto && ww = SolverType::solve(trans * A * trans, vv);
            w = trans * ww;
        } else {
            w = SolverType::solve(A, v);
        }
        return w;
    }

private:
    const bool apply = false;
    const TransformationMatrixType & trans;
};

////////////////////////////////////////////////////////////////////////////////
//  Custom Log Class
////////////////////////////////////////////////////////////////////////////////

template<typename T>
struct OptimizerLog;

/*
  Type trait specifying the log output for the Levenberg-Marquardt method.
   The member name contains the method name that is printed to the screen.
   The functions header() and log() specify what is printed in the last
   columns of the iteration table in the header and the line for each step,
   respectively.
*/
template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
    struct OptimizerLog<invlib::LevenbergMarquardt<RealType, DampingMatrix, Solver>>
{
    static constexpr auto name = "Levenberg-Marquardt";

    static std::string header()
    {
        std::string out = "Gamma Factor";
        return out;
    }

    static std::string log(
        const invlib::LevenbergMarquardt<RealType, DampingMatrix, Solver> &g,
        Vector & gamma_history, size_t i
        )
    {
        std::string lambda = std::to_string(g.get_lambda());
        std::string out(15 - lambda.size(), ' ');
        out += lambda;
        gamma_history[i] = g.get_lambda();
        return out;
    }

};

/*
  Specialization for Gauss-Newton method.
*/
template
<
typename RealType,
typename Solver
>
    struct OptimizerLog<invlib::GaussNewton<RealType, Solver>>
{
    static constexpr auto name = "Gauss-Newton";

    static std::string header()
    {
        return "";
    }

    static std::string log(
        const invlib::GaussNewton<RealType, Solver> &,
        Vector &, size_t
        )
    {
        return "";
    }

};

/*
  Specialized log type for specifying the output generated by invlib.
  The init function is called at the beginning of the compute method
  of an invlib MAP<...> object. Here this will print out the general
  info and the table header.
  In each iteration step, the step() method is called which prints
  a line for the iteration table and logs the gamma history in if
  the optimization method is the Levenberg-Marquardt method.
  Finally, finalize is called which finalizes the iteration table
  and prints out summarizing information to the screen.
*/
template
<
invlib::LogType type
>
class ArtsLog
{

public:

    ArtsLog(unsigned int v, Vector & g, bool l = false)
        : verbosity(v), gamma_history(g), linear(l),
        finalized(false) {}

    ~ArtsLog()
    {
        if (!finalized)
        {
            std::cout << invlib::separator() << std::endl << std::endl;
            std::cout << "Error during OEM computation." << std::endl;
            std::cout << std::endl;
            std::cout << invlib::center("----") << std::endl;
            std::cout << std::endl;
        }
    }

    template <typename... Params>
    void init(Params... params)
    {
        auto tuple = std::make_tuple(params...);

        auto & y       =  std::get<4>(tuple);
        scaling_factor =  1.0 / static_cast<Numeric>(y.nelem());
        std::cout << std::endl;
        std::cout << invlib::center("MAP Computation") << std::endl;

        // Print formulation.
        int formulation = static_cast<int>(std::get<6>(tuple));
        switch (formulation)
        {
        case 0:
                std::cout << "Formulation: Standard" << std::endl;
                break;
        case 1:
                std::cout << "Formulation: N-Form" << std::endl;
                break;

        case 2:
                std::cout << "Formulation: M-Form" << std::endl;
                break;
        }

        // Print optimization method.
        using OptimizationType =
            typename std::tuple_element<5, decltype(tuple)>::type;
        std::cout << "Method:      " << invlib::OptimizerLog<OptimizationType>::name;
        std::cout << std::endl;

        if (verbosity >= 1)
        {
            std::cout << std::endl;
            std::cout << std::setw(5) << "Step" << std::setw(15) << "Total Cost";
            std::cout << std::setw(15) << "x-Cost" << std::setw(15) << "y-Cost";
            std::cout << std::setw(15) << "Conv. Crit.";
            std::cout << std::setw(15) << OptimizerLog<OptimizationType>::header();
            std::cout << std::endl << invlib::separator() << std::endl;
        }
    }

    template <typename... Params>
    void step(Params... params)
    {
        if (verbosity >= 1)
        {
            auto tuple = std::make_tuple(params...);
            using OptimizationType =
                typename std::tuple_element<5, decltype(tuple)>::type;

            std::cout<< std::setw(5)  << std::get<0>(tuple);
            std::cout<< std::setw(15) << scaling_factor * std::get<1>(tuple);
            std::cout<< std::setw(15) << scaling_factor * std::get<2>(tuple);
            std::cout<< std::setw(15) << scaling_factor * std::get<3>(tuple);

            if (std::isnan(std::get<4>(tuple)))
            {
                std::cout<< std::setw(15) << " ";
            } else {
                std::cout<< std::setw(15) << std::get<4>(tuple);
            }
            std::cout<< OptimizerLog<OptimizationType>::log(std::get<5>(tuple), gamma_history, std::get<0>(tuple));
            std::cout << std::endl;
        }
    }

    template <typename... Params>
    void finalize(Params... params)
    {
        if (verbosity >= 1)
        {
            std::cout << invlib::separator() << std::endl;

            auto tuple = std::make_tuple(params...);
            std::cout << std::endl;

            std::cout << "Total number of steps:            ";
            std::cout << std::get<1>(tuple) << std::endl;
            std::cout << "Final scaled cost function value: ";
            std::cout << std::get<2>(tuple) * scaling_factor << std::endl;

            bool converged = std::get<0>(tuple);
            if (converged)
            {
                std::cout << "OEM computation converged." << std::endl;
            } else if (linear) {
                std::cout << "Linear OEM computation finished." << std::endl;
            } else {
                std::cout << "OEM computation DID NOT converge!" << std::endl;
            }
        }

        finalized = true;

    }

    template <typename... Params>
    void time(Params... params)
    {
        if (verbosity >= 1)
        {
            auto tuple = std::make_tuple(params...);
            std::cout << std::endl;
            std::cout << "Elapsed Time for Retrieval:                       ";
            std::cout << std::get<0>(tuple) << std::endl;
            std::cout << "Time in inversion_iterate Agenda (No Jacobian):   ";
            std::cout << std::get<1>(tuple) << std::endl;
            std::cout << "Time in inversion_iterate Agenda (With Jacobian): ";
            std::cout << std::get<2>(tuple) << std::endl;

        }

        std::cout << std::endl;
        std::cout << invlib::center("----") << std::endl;
        std::cout << std::endl;

    }


private:

    int      verbosity;
    Vector & gamma_history;
    Numeric  scaling_factor;
    bool     linear, finalized;

};

////////////////////////////////////////////////////////////////////////////////
//  Exception Handling
////////////////////////////////////////////////////////////////////////////////

/* If an error occurs in the MAP computation, invlib throws a nested exception,
 * which contains the original exception. This function extracts the error messages
 * and transforms them to a vector of strings.
 */
template <typename E>
std::vector<std::string> handle_nested_exception(
    const E & e,
    int level = 0
    )
{
    const std::exception * re;
    std::vector<std::string> errors{};

    re = dynamic_cast<const std::exception *>(&e);
    if (re)
    {
        std::string s{};

        // If invlib level, extend error description.
        if (level == 0)
        {
            s = "Run-time error in oem computation: ";
        }

        s += re->what();
        errors.push_back(s);
    }

    try {
        std::rethrow_if_nested(e);
    } catch(const std::exception & ne) {
        std::vector<std::string> sv(
            handle_nested_exception(ne, level + 1)
            );
        errors.insert(errors.end(), sv.begin(), sv.end());
    } catch(...) {}
    return errors;
}

#endif // oem_h
