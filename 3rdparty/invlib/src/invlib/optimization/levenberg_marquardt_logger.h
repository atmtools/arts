#ifndef OPTIMIZATION_LEVENBERG_MARQUARDT_LOGGER_H
#define OPTIMIZATION_LEVENBERG_MARQUARDT_LOGGER_H

#include <iostream>
#include <iomanip>

namespace invlib
{

enum class Verbosity {SILENT, VERBOSE};

template <Verbosity V, std::ostream &stream> class LevenbergMarquardtLogger;

template<std::ostream &stream>
class LevenbergMarquardtLogger<Verbosity::SILENT, stream>
{
public:
    template <typename... Params>
    static void init(Params... params) {}

    template <typename... Params>
    static void step(Params... params) {}

    template <typename... Params>
    static void finalize(Params... params) {}
};

template<std::ostream &stream>
class LevenbergMarquardtLogger<Verbosity::VERBOSE, stream>
{
public:

    void separator(unsigned int length)
    {
        for (unsigned int i = 0; i < length; i++)
            stream << "-";
        stream << std::endl;
    }

    void init()
    {
        separator(30) << std::endl;
        stream << "LevenbergMarquardt Minimization" << std::endl << std::endl;
        stream << "Step      Cost      \u03bb" << std::endl;
        separator(30);
    }

    template <typename T>
        static void step(const T& t)
    {
        stream << std::setw(10) << t.step_count;
        stream << std::setw(10) << t.current_cost;
        stream << std::setw(10) << t.lambda << std::endl;
    }

    template <typename... Params>
    static void finalize(Params... params) {}
};

}      // namespace invlib

#endif // LEVENBERG_MARQUARDT_LOGGER_H
