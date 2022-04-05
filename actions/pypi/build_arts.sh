#!/bin/sh

PYTHONDIR=/opt/python/$INPUT_PYTHON_VERSION
echo "PYTHONDIR=$PYTHONDIR"
if [[ ! -f "${PYTHONDIR}/bin/python3" ]]; then
    echo "Python directory ${PYTHONDIR} not found"
    exit 1
fi
export PATH=$PYTHONDIR/bin:$PATH

python3 -m pip install --upgrade pip
python3 -m pip install --upgrade setuptools wheel auditwheel twine docutils lark-parser matplotlib netCDF4 numpy pytest scipy xarray cmake

# Get arts source
mkdir -p arts && cd arts && git init .
git remote add origin https://github.com/${GITHUB_REPOSITORY}
git fetch --prune --depth=1 origin

# Build pyarts
git checkout --force ${GITHUB_REF#refs/heads/}; mkdir -p build; cd build;
cmake -DCMAKE_PREFIX_PATH=$PYTHONDIR -DCMAKE_BUILD_TYPE=Release -DENABLE_FORTRAN=1 -DBLAS_blas_LIBRARY=/usr/lib64/atlas/libtatlas.so -DNUM_PYARTS_WSM=2 -DNUM_PYARTS_WSV=1 -DNUM_PYARTS_WSC=1 -DNUM_PYARTS_WSG=1 ..
make -j2 arts
make -j1 pyarts

# Packaging
cd python
python3 setup.py sdist bdist_wheel
auditwheel repair dist/pyarts*.whl
ls -lh dist/pyarts*.whl

# Upload to PyPi
if [[ ${GITHUB_REPOSITORY} == "atmtools/arts" ]]; then
    python3 -m twine upload wheelhouse/pyarts*.whl -u __token__ -p $INPUT_PYPI_ACCESS
fi

