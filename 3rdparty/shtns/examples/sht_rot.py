
from pylab import *
import shtns

Lmax = 13  # degree of fields to rotate. Note that arbitrary rotations require Mmax=Lmax

sh = shtns.sht(Lmax)         # Spherical Harmonic transform object
rot = shtns.rotation(Lmax)   # corresponding rotation object


# Define the rotation:
alpha = pi/2
beta = pi/2
gamma = 0
rot.set_angles_ZYZ(alpha, beta, gamma)   # defines rotation by giving Euler angles (alpha, beta, gamma), around axes Z,Y and Z respectively.
# or alternatively:
# rot.set_angles_ZXZ(alpha, beta, gamma)  # same, but around axes Z, X and Z
# rot.set_angle_axis(theta, Vx, Vy, Vz)   # same, but rotation of angle "theta" around vector (Vx,Vy,Vz)

### If you like, you can get the Wigner-d matrix https://en.wikipedia.org/wiki/Wigner_D-matrix#Wigner_(small)_d-matrix
# W3 = rot.wigner_d_matrix(3)  # Wigner-d rotation matrix for degree l=3


Qlm = sh.spec_array()   # empty spectral field
Qlm[1] = 1.  # a value somewhere

Rlm = rot.apply_real(Qlm)   # Rlm is the field Qlm rotated by the rotation defined above.
print('rotated by (%g,%g,%g) degrees (ZYZ)' % (rot.alpha*180/pi, rot.beta*180/pi, rot.gamma*180/pi))

sh.set_grid(nlat=64, nphi=128)
q = sh.synth(Qlm)
r = sh.synth(Rlm)

# plot the fields if possible
# inspired by: https://scipython.com/book/chapter-8-scipy/examples/visualizing-the-spherical-harmonics/

try:
	from matplotlib import cm, colors
	from mpl_toolkits.mplot3d import Axes3D
except:
	print("matplotlib required to plot.")
	exit()


phi = arange(sh.nphi)*2*pi/sh.nphi
theta = arccos(sh.cos_theta)
phi, theta = meshgrid(phi, theta)
# The Cartesian coordinates of the unit sphere
x = sin(theta) * cos(phi)
y = sin(theta) * sin(phi)
z = cos(theta)

# normalize to 0 - 1
fmin, fmax = q.min(), q.max()
q = (q-fmin)/(fmax-fmin)
r = (r-fmin)/(fmax-fmin)

fig = plt.figure(figsize=plt.figaspect(0.5))
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.seismic(q))
ax1.set_title('original')
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.seismic(r))
ax2.set_title('rotated')
for ax in (ax1, ax2):
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
show()
