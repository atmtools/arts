""" PyARTS, the Python interface for ARTS

PyARTS provides:

- full ARTS library for direct access from Python
- full scripting interface for ARTS to replace controlfiles
- reading and writing of ARTS XML data files
"""

import logging
import sys
import subprocess
import shutil

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import remove
from os.path import abspath, dirname, isfile, join

import builtins

DOCLINES = (__doc__ or "").split("\n")
__version__ = open(join("@ARTS_SRC_DIR@", "VERSION")).read().strip()

VERSION_TUPLE = __version__.split(".")
STABLE = int(VERSION_TUPLE[1]) % 2 == 0

here = abspath(dirname(__file__))

try:
    arts_libname = "libarts_api.so"
    lib_path = join("@ARTS_BINARY_DIR@", "src", arts_libname)
    if isfile(join("pyarts", "workspace", arts_libname)):
        remove(join("pyarts", "workspace", arts_libname))
    shutil.copy(lib_path, join("pyarts", "workspace"))
except:
    raise Exception(
        "Could not find ARTS API, which is required for the Python "
        "interface. Please make sure the installation was "
        "successful."
    )


setup(
    name="pyarts",
    version=__version__ + ("" if STABLE else ".dev0"),
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    url="https://github.com/atmtools/arts",
    author="The Typhon developers",
    author_email="arts.mi@lists.uni-hamburg.de",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.7",
    ],
    packages=find_packages(exclude=["contrib", "doc", "tests*"]),
    python_requires=">=3.6",
    install_requires=[
        "docutils",
        "matplotlib>=1.4",
        "netCDF4>=1.1.1",
        "numpy>=1.13",
        "scipy>=0.15.1",
        "setuptools>=0.7.2",
    ],
    extras_require={
        "docs": ["sphinx_rtd_theme"],
        "tests": ["pytest", "pint", "gdal",],
    },
    package_data={"": ["*.so"],},
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    include_package_data=True,
)
