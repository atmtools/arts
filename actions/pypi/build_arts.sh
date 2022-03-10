#!/bin/sh

# Get arts source
mkdir -p arts && cd arts && git init .
git remote add origin https://github.com/${GITHUB_REPOSITORY}
git fetch --prune --depth=1 origin

# Build pyarts
git checkout --force ${GITHUB_REF#refs/heads/}; mkdir -p build; cd build;
cmake3 -DCMAKE_BUILD_TYPE=Release -DENABLE_FORTRAN=1 -DBLAS_blas_LIBRARY=/usr/lib64/atlas/libtatlas.so -DLAPACK_lapack_LIBRARY=/usr/lib64/atlas/libtatlas.so ..
make -j2 pyarts

# Packaging
cd python
python3 setup.py sdist bdist_wheel
auditwheel repair dist/pyarts*.whl
python3 -m twine upload wheelhouse/pyarts*.whl -u __token__ -p $INPUT_PYPI_ACCESS
