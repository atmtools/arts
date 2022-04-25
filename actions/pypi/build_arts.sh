#!/bin/sh

PYTHONDIR=/opt/python/$INPUT_PYTHON_VERSION
echo "PYTHONDIR=$PYTHONDIR"
echo "GITHUB_REPOSITORY=@${GITHUB_REPOSITORY}@"
if [[ ! -f "${PYTHONDIR}/bin/python3" ]]; then
    echo "Python directory ${PYTHONDIR} not found"
    exit 1
fi
export PATH=$PYTHONDIR/bin:$PATH

python3 -m pip install --upgrade pip
python3 -m pip install --upgrade setuptools build wheel auditwheel twine docutils lark-parser matplotlib netCDF4 numpy pytest scipy xarray cmake

# Get arts source
mkdir -p arts && cd arts && git init .
git remote add origin https://github.com/${GITHUB_REPOSITORY}
git fetch --prune --depth=1 origin

# Build pyarts
git checkout --force ${GITHUB_REF#refs/heads/}; mkdir -p build; cd build;
cmake -DCMAKE_PREFIX_PATH=$PYTHONDIR -DCMAKE_BUILD_TYPE=Release -DENABLE_FORTRAN=1 -DBLAS_blas_LIBRARY=/usr/lib64/atlas/libtatlas.so -DNUM_PYARTS_WSM=2 -DNUM_PYARTS_WSV=1 -DNUM_PYARTS_WSC=1 -DNUM_PYARTS_WSG=1 ..
echo "########## CMakeCache.txt ##########"
cat CMakeCache.txt
echo "########## CMakeCache.txt ##########"
make -j2 arts
make -j1 pyarts
echo "########## Check Python version ##########"
make check-pyversion

# Packaging
cd python
python3 -m build
auditwheel repair dist/pyarts*.whl
ls -lh dist/pyarts*.whl

# Upload to PyPi
if [[ ${GITHUB_REPOSITORY} == "atmtools/arts" ]]; then
    python3 -m twine upload --verbose wheelhouse/pyarts*.whl -u __token__ -p $INPUT_PYPI_ACCESS
#elif [[ ${GITHUB_REPOSITORY} == "olemke/arts" ]]; then
#    python3 -m twine upload --verbose --repository testpypi wheelhouse/pyarts*.whl -u __token__ -p $INPUT_PYPI_ACCESS
fi

