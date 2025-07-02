import numpy as np
import matplotlib.pyplot as plt

h2o_vmr_profile_parameters = [1e3, 1e-8, 1e5, 1e-2]

b = (
    np.log10(h2o_vmr_profile_parameters[3]) - np.log10(h2o_vmr_profile_parameters[1])
) / (np.log10(h2o_vmr_profile_parameters[2]) - np.log10(h2o_vmr_profile_parameters[0]))
a = np.log10(h2o_vmr_profile_parameters[1]) - b * np.log10(
    h2o_vmr_profile_parameters[0]
)
a10 = np.pow(10, a)

nump = 100
abs_p = np.linspace(1, 100000, num=nump)
vmr_h2o_default_profile = np.zeros(nump)

for i in range(len(vmr_h2o_default_profile)):
    if abs_p[i] < h2o_vmr_profile_parameters[0]:
        vmr_h2o_default_profile[i] = h2o_vmr_profile_parameters[1]
    else:
        vmr_h2o_default_profile[i] = a10 * np.pow(abs_p[i], b)

old = np.ones(nump) * 1e-3

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.semilogy(vmr_h2o_default_profile, abs_p / 100, label="New VMR")
ax1.semilogy(old, abs_p / 100, label="Old VMR")
ax1.invert_yaxis()
ax1.set_xlabel("VMR")
ax1.set_ylabel("Pressure [hPa]")
ax1.legend()
ax2.loglog(vmr_h2o_default_profile, abs_p / 100, label="New VMR")
ax2.loglog(old, abs_p / 100, label="Old VMR")
ax2.invert_yaxis()
ax2.set_xlabel("VMR")
ax2.legend()
plt.savefig("vmr.png")
plt.show()
