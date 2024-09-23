# Python setup

from distutils.core import setup, Extension
from numpy import get_include

def getver():
    with open('CHANGELOG.md') as f:
        for l in f:
            s=l.find('* v')
            if s>=0:
                return l[s+3:].split()[0]
    return '3.5.2'

numpy_inc = get_include()		#  NumPy include path.
objs = "sht_init.o sht_kernels_a.o sht_kernels_s.o sht_odd_nlat.o sht_fly.o"
shtns_o = objs.split()			# transform to list of objects
libdir = "/usr/local"
if len(libdir) == 0:
	libdir = []
else:
	libdir = [libdir+"/lib"]
cargs = '-std=c99  -DSHTNS_VER="' + getver() +'"'
largs = ""
libs = "-lfftw3 -lm "
libslist = libs.replace('-l','').split()	# transform to list of libraries

shtns_module = Extension('_shtns', sources=['shtns_numpy_wrap.c'],
	extra_objects=shtns_o, depends=shtns_o,
	extra_compile_args=cargs.split(),
	extra_link_args=largs.split(),
	library_dirs=libdir,
	libraries=libslist,
	include_dirs=[numpy_inc])

setup(name='SHTns',
	version=getver(),
	description='High performance Spherical Harmonic Transform',
	license='CeCILL',
	author='Nathanael Schaeffer',
	author_email='nathanael.schaeffer@univ-grenoble-alpes.fr',
	url='https://bitbucket.org/nschaeff/shtns',
	ext_modules=[shtns_module],
	py_modules=["shtns"],
	requires=["numpy"],
	)
