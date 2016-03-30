#ifndef OPTIMIZATION_LOG
#define OPTIMIZATION_LOG

namespace invlib
{

enum class Verbosity {SILENT, VERBOSE};

template
<
Verbosity V,
template <typename... Args> class T,
typename... Args
>
class Logger;

template
<
template <typename... Args> class T,
typename... Args
>
class Logger<Verbosity::SILENT, T, Args...>
{
    template <typename... Params>
    static void init(Params... params) {}

    template <typename... Params>
    static void step(Params... params) {}

    template <typename... Params>
    static void finalize(Params... params) {}
};

}

#endif // LOG
