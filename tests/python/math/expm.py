import pyarts3 as pa
import numpy as np
import scipy as sp

a = 1e-1
b = 1e-2
c = 1e-3
d = 1e-4
u = 1e-5
v = 1e-6
w = 1e-7

k = pa.arts.Propmat([a,b,c,d,u,v,w])

x = -1
a = np.array(k.exp(x))

b = sp.linalg.expm(x*np.array(k.as_matrix()))

x = np.allclose(a, b, rtol=1e-6, atol=1e-6)

assert x, "Matrix exponential of Propmat is not equal to scipy.linalg.expm"
