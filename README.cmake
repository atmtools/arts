Here are the steps to use cmake to build arts.

Create a build directory in the arts directory (or anywhere else) and run
cmake to configure arts there and then build arts:

mkdir build
cd build
cmake ..
make -j2

To build a release version without assertions and debugging symbols use:

cmake -DCMAKE_BUILD_TYPE=Release ..
make -j2

To switch back to the debug version use:

cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make -j2

You can also disable certain features:

Disable assertions:
cmake -DNO_ASSERT=1 ..

Disable assertions:
cmake -DNO_OPENMP=1 ..

Disable warnings being treated as errors:
cmake -DNO_WERROR=1 ..


If you want to compile with the Intel compiler, start with an empty build
directory and run:

CC=icc CXX=icpc cmake ..


