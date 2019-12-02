"""
arts: The Python interface for ARTS
"""

import logging
import sys
import subprocess

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os.path import (abspath, dirname, join)

import builtins

DOCLINES = (__doc__ or '').split("\n")
__version__ = open(join("@ARTS_SRC_DIR@", 'VERSION')).read().strip()


here = abspath(dirname(__file__))

setup(
    name='arts',
    version=__version__,
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    url='https://github.com/atmtools/arts',
    author='The Typhon developers',
    author_email='arts.mi@lists.uni-hamburg.de',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
    ],

    packages=find_packages(exclude=['contrib', 'doc', 'tests*']),
    python_requires='>=3.6',

    install_requires=[
        'docutils',
        'imageio',
        'matplotlib>=1.4',
        'netCDF4>=1.1.1',
        'numba',
        'numexpr',
        'numpy>=1.13',
        'scipy>=0.15.1',
        'setuptools>=0.7.2',
    ],
    extras_require={
        'docs': ['sphinx_rtd_theme'],
        'tests': [
            'pytest',
            'pint',
            'gdal',
        ],
    },
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    package_data={
       'lib': ['/home/simonpf/build/arts/libarts_api.so'],
    },
    include_package_data=True,
)
