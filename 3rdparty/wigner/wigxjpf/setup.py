#!/usr/bin/env python

from setuptools import setup, find_packages

try:
    import pycparser as dummy
except ImportError:
    print('##############################################\n'
          '!        ERROR! Caught import error.         !\n'
          '!      Please, install pycparser first:      !\n'
          '!                                            !\n'
          '!            pip install pycparser           !\n'
          '!                                            !\n'
          '##############################################\n\n'
          'For details, see '
          'https://github.com/pypa/setuptools/issues/391#issuecomment-202919511'
          '\n\n')

setup(name='pywigxjpf',
      version='1.11',
      description='Wigner Symbols',
      long_description="""
WIGXJPF evaluates Wigner 3j, 6j and 9j symbols
accurately using prime factorisation and multi-word integer arithmetic.

Known issue
-----------

Please install pycparser first.

(See https://github.com/pypa/setuptools/issues/391#issuecomment-202919511 )

Inline usage information
------------------------

Available inline::

    import pywigxjpf as wig
    help(wig)                # For interfaces.
    help(wig.pywigxjpf)      # For usage information.

Library usage
-------------

The python interface to wigxjpf uses cffi.

Defines seven functions::

    wig_table_init(max_two_j,wigner_type)
    wig_table_free()
    wig_temp_init(max_two_j)
    wig_temp_free()

    wig3jj(jj1,jj2,jj3, mm1,mm2,mm3)
    wig6jj(jj1,jj2,jj3, jj4,jj5,jj6)
    wig9jj(jj1,jj2,jj3, jj4,jj5,jj6, jj7,jj8,jj9)

Note that the arguments are to be given as integers, with
twice the numeric value (this is what jj tries to indicate).
I.e. half-integer arguments will be passed as odd integers.

The two init functions must be called before evaluating any
symbol.

In addition, interfaces that take an array with the arguments
are also provided::

    wig3jj_array([jj1,jj2,jj3, mm1,mm2,mm3])
    wig6jj_array([jj1,jj2,jj3, jj4,jj5,jj6])
    wig9jj_array([jj1,jj2,jj3, jj4,jj5,jj6, jj7,jj8,jj9])
""",
      long_description_content_type='text/x-rst',
      license = "GPLv3",
      author='C. Forssen and H. T. Johansson',
      url='http://fy.chalmers.se/subatom/wigxjpf',
      packages=['pywigxjpf'],
      # package_data={'pywigxjpf': ['../lib/libwigxjpf_shared.*']},
      setup_requires=['numpy',"cffi>=1.0.0"],
      cffi_modules=["pywigxjpf/pywigxjpf_ffi_builder.py:ffibuilder"],
      install_requires=['numpy',"cffi>=1.0.0"],
     )
