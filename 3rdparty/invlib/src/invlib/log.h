#ifndef LOG_H
#define LOG_H

#include <tuple>
#include <string>
#include <iostream>
#include <iomanip>

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

enum class LogType {MAP, OPT_GN, OPT_LM, SUB};

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
        std::string out(15 - lambda.size(), ' ');
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
//     Log Functions      //
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
            typename std::tuple_element<4, decltype(tuple)>::type;

        std::cout<< std::setw(5) << std::get<0>(tuple);
        std::cout<< std::setw(15) << std::get<1>(tuple);
        std::cout<< std::setw(15) << std::get<2>(tuple);
        std::cout<< std::setw(15) << std::get<3>(tuple);
        std::cout<< OptimizerLog<OptimizationType>::log(std::get<4>(tuple));
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

}      // namespace invlib

#endif // LOG_H
