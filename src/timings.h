/*!
  \file   timer.h
  \author Simon Pfreundschuh <simonpf@chalmers.se>
  \date   Fri Apr 17 16:17:54 2015

  \brief Timings class.
*/

#ifndef timings_h
#define timings_h

#include "omp.h"
#include "arts.h"
#include <iostream>

//! Timings Class
/*!

  Simple class that handles a variable number of timers and computes
  the total time. The Timings class uses the OMP library to measure the
  runtime and therefore only works correctly if OMP is available.

*/
class Timings
{
public:

    friend std::ostream& operator<<( std::ostream &os, const Timings &timer );

    //! Simple Constructor
    /*!
      Constructs a Timings object without any timers.
     */
    Timings()
        : running(), times(), stamps(), names()
    {
        ntimers = 0;
        nrunning = 0;
        stamp = 0.0;
        total_time = 0.0;
    }

    //! Add timer.
    /*!
      Add timer with given name to Timings object. The name is used in the
      timer output.

      \param[in] name The name of the time used in the timer output.
     */
    Index add_timer( String name )
    {

        ntimers++;
        running.push_back( false );
        times.push_back( 0.0 );
        names.push_back( name );
        stamps.push_back( 0.0 );

        return (Index) ntimers - 1;
    }

    //! Get Wall Time
    /*!
      Returns the current wall time computed using omp_get_wtime(). If no OpenMP
      library is available, the functions return 0.0.

      \return The current wall time if an OpenMP implementation is available, 0.0
      otherwise.
    */
    Numeric get_time()
    {
        #ifdef OMP
        return omp_get_wtime();
        #else
        return 0.0;
        #endif
    }

    //! Start/Stop Timings
    /*!
      Starts or stops the timer with the given index. The indices of the timers
      reflect the order in which the timers have been added to the Timings object
      starting with 0 for the first timer. If the given timer is not running
      currently, it is started. If it is running, it is stopped and the passed
      time added to the total time of the timer.
    */
    void mark( Index i )
    {
        if (i < ntimers)
        {
            if (running[i])
            {
                Numeric t = get_time();
                times[i] += t - stamps[i];
                running[i] = false;

                nrunning--;
                if (nrunning == 0)
                    total_time += get_time() - stamp;

            }
            else
            {

                if (nrunning == 0)
                    stamp = get_time();

                nrunning++;

                stamps[i] = get_time();
                running[i] = true;

            }
        }
    }


private:

    Index ntimers, nrunning;
    Numeric total_time, stamp;
    std::vector<bool> running;
    std::vector<Numeric> times, stamps;
    std::vector<String> names;

};

#endif // timer_h
