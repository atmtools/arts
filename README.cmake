Here are the steps to use cmake to build arts.

Create a build directory in the arts directory (or anywhere else) and run
cmake to configure arts there and then build arts:

mkdir build
cd build
cmake ..
make

With cmake, you don't have to go into the src directory to only build the arts
executable. Just run 'make arts' inside your build directory.

If you have a multicore machine, don't forget to use the -j option to speed up
the compilation:

make -jX

Where X is the number of parallel build processes.
X=(Number of CPUs)+1 gives you usually the fastest compilation time.

To build a release version without assertions or debugging symbols use:

cmake -DCMAKE_BUILD_TYPE=Release ..
make

To switch back to the debug version use:

cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make

This is also the default configuration if you run cmake without
options in an empty build directory.

You can also disable certain features:

Disable assertions:
cmake -DNO_ASSERT=1 ..

Disable OpenMP:
cmake -DNO_OPENMP=1 ..

Disable NetCDF:
cmake -DNO_NETCDF=1 ..

Disable warnings being treated as errors:
cmake -DNO_WERROR=1 ..


If you want to compile with the Intel compiler, start with an empty build
directory and run:

CC=icc CXX=icpc cmake ..


