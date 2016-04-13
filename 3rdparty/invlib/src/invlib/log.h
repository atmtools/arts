#ifndef LOG_H
#define LOG_H

#include <tuple>
#include <string>
#include <iostream>
#include <iomanip>

#include "invlib/profiling/timer.h"

namespace invlib
{

// ------------------------ //
//    Forward Declarations  //
// ------------------------ //

template
<
typename RealType,
typename Solver
>
class GaussNewton;

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
class LevenbergMarquardt;

// ------------------ //
//    Standard Log    //
// ------------------ //

enum class LogType {MAP, OPT_GN, OPT_LM, SUB, SOL_CG};

template
<
LogType type
>
class StandardLog
{

public:

    StandardLog(unsigned int v) : verbosity(v) {}

    template <typename... Params>
    void init(Params... params) {}

    template <typename... Params>
    void step(Params... params) {}

    template <typename... Params>
    void finalize(Params... params) {}

    template <typename... Params>
    void time(Params... params) {}

private:

    int verbosity;
};

// ---------------------- //
//  Formatting Functions  //
// ---------------------- //

std::string center(const std::string &s, int width = 80)
{
    auto padding_length = (width - s.size()) / 2;
    std::string padding(padding_length, ' ');
    std::string centered(padding);
    centered += s;
    centered += padding;
    return centered;
}

std::string separator(int width = 80)
{
    std::string separator(width, '-');
    return separator;
}
// ------------------- //
//      Type Names     //
// ------------------- //

template<typename T>
struct OptimizerLog;

template
<
typename RealType,
typename DampingMatrix,
typename Solver
>
struct OptimizerLog<LevenbergMarquardt<RealType, DampingMatrix, Solver>>
{
    static constexpr auto name = "Levenberg-Marquardt";

    static std::string header()
    {
        std::string lambda = "\u03BB";
        std::string out((15 - lambda.size()) / 2, ' ');
        out += lambda;
        return out;
    }

    static std::string log(const LevenbergMarquardt<RealType, DampingMatrix, Solver> &g)
    {
        std::string lambda = std::to_string(g.get_lambda());
        std::string out(15 - lambda.size(), ' ');
        out += lambda;
        return out;
    }

};

template
<
typename RealType,
typename Solver
>
struct OptimizerLog<GaussNewton<RealType, Solver>>
{
    static constexpr auto name = "Gauss-Newton";

    static std::string header()
    {
        return "";
    }

    static std::string log(const GaussNewton<RealType, Solver> &)
    {
        return "";
    }

};

// ---------------------- //
//     MAP Class          //
// ---------------------- //

template<>
template<typename... Params>
void StandardLog<LogType::MAP>::init(Params... params)
{
    auto tuple = std::make_tuple(params...);
    if (verbosity >= 1)
    {
        std::cout << std::endl;
        std::cout << center("MAP Computation") << std::endl;

        // Print formulation.
        int formulation = static_cast<int>(std::get<0>(tuple));
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
            typename std::tuple_element<1, decltype(tuple)>::type;
        std::cout << "Method:      " << OptimizerLog<OptimizationType>::name;
        std::cout << std::endl;

        if (verbosity >= 2)
        {
            std::cout << std::endl;
            std::cout << std::setw(5) << "Step" << std::setw(15) << "Total Cost";
            std::cout << std::setw(15) << "x-cost" << std::setw(15) << "y-cost";
            std::cout << std::setw(15) << "conv. crit.";
            std::cout << OptimizerLog<OptimizationType>::header();
            std::cout << std::endl << separator() << std::endl;
        }
    }
}

template<>
template<typename... Params>
void StandardLog<LogType::MAP>::step(Params... params)
{

    if (verbosity >= 2)
    {
        auto tuple = std::make_tuple(params...);
        using OptimizationType =
            typename std::tuple_element<5, decltype(tuple)>::type;

        std::cout<< std::setw(5) << std::get<0>(tuple);
        std::cout<< std::setw(15) << std::get<1>(tuple);
        std::cout<< std::setw(15) << std::get<2>(tuple);
        std::cout<< std::setw(15) << std::get<3>(tuple);
        std::cout<< std::setw(15) << std::get<4>(tuple);
        std::cout<< OptimizerLog<OptimizationType>::log(std::get<5>(tuple));
        std::cout << std::endl;
    }
}

template<>
template<typename... Params>
void StandardLog<LogType::MAP>::finalize(Params... params)
{
    if (verbosity >= 1)
    {
        if (verbosity >= 2)
            std::cout << separator() << std::endl;

        auto tuple = std::make_tuple(params...);
        std::cout << std::endl;

        std::cout << "Total number of steps: ";
        std::cout << std::get<1>(tuple) << std::endl;
        std::cout << "Final cost function value: ";
        std::cout << std::get<2>(tuple) << std::endl;

        bool converged = std::get<0>(tuple);
        if (converged)
        {
            std::cout << "MAP Computation converged." << std::endl;
        }
        else
        {
            std::cout << "MAP Computation NOT converged!" << std::endl;
        }


    }
}

template<>
template<typename... Params>
void StandardLog<LogType::MAP>::time(Params... params)
{
    if (verbosity >= 1)
    {
        auto tuple = std::make_tuple(params...);
        std::cout << std::endl;
        std::cout << "Total time           : ";
        std::cout << std::get<0>(tuple) << std::endl;
        std::cout << "Time in evaluate(...): ";
        std::cout << std::get<1>(tuple) << std::endl;
        std::cout << "Time in Jacobian(...): ";
        std::cout << std::get<2>(tuple) << std::endl;

        auto t1 = multiply_mm_time.count();
        auto t2 = multiply_mv_time.count();
        auto t3 = solve_time.count();
        auto t4 = invert_time.count();

        if (t1 > 0.0 || t2 > 0.0 || t3 > 0.0 || t4 > 0.0)
        {
            std::cout << std::endl;
            std::cout << "Time in MM multiply(...): ";
            std::cout << t1 << std::endl;
            std::cout << "Time in MV multiply(...): ";
            std::cout << t2 << std::endl;
            std::cout << "Time in solve(...): ";
            std::cout << t3 << std::endl;
            std::cout << "Time in invert(...): ";
            std::cout << t4 << std::endl;
        }
    }
}
// -------------------- //
//     MAP Class        //
// -------------------- //

template<>
template<typename... Params>
void StandardLog<LogType::SOL_CG>::init(Params... params)
{
    auto tuple = std::make_tuple(params...);
    if (verbosity >= 1)
    {
        std::cout << std::endl;
        std::cout << "CG Solver:" << std::endl;
        std::cout << "\tTolerance:             " << std::get<0>(tuple) << std::endl;
        std::cout << "\tInitial Residual Norm: " << std::get<1>(tuple) << std::endl;
        std::cout << "\tRight-hand side Norm:  " << std::get<2>(tuple) << std::endl;
    }
}

template<>
template<typename... Params>
void StandardLog<LogType::SOL_CG>::step(Params... params)
{

    if (verbosity >= 1)
    {
        auto tuple = std::make_tuple(params...);
        std::cout<< "Step " << std::setw(5) << std::get<0>(tuple) << ", ";
        std::cout<< "Normalized Residual: " << std::get<1>(tuple) << std::endl;
    }
}

template<>
template<typename... Params>
void StandardLog<LogType::SOL_CG>::finalize(Params... params)
{
    if (verbosity >= 1)
    {
        auto tuple = std::make_tuple(params...);
        std::cout << "Conjugate Gradient method converged after ";
        std::cout << std::get<0>(tuple) << " steps." << std::endl << std::endl;
    }
}

}      // namespace invlib

#endif // LOG_H
