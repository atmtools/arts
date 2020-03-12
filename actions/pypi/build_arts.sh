#!/bin/sh

# Get arts source
git clone https://github.com/simonpf/arts/

# Build pyarts
cd arts; git checkout python_integration; mkdir build; cd build;
cmake3 -DENABLE_FORTRAN=1 -DBLAS_blas_LIBRARY=/usr/lib64/atlas/libtatlas.so -DLAPACK_lapack_LIBRARY=/usr/lib64/atlas/libtatlas.so ..
make pyarts

# Packaging
cd python
python3 setup.py sdist bdist_wheel
auditwheel repair dist/pyarts*.whl
python3 -m twine upload wheelhouse/pyarts*.whl -u __token__ -p ${pypi_access}
