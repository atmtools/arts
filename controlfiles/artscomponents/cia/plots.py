import typhon
import numpy as np
import matplotlib.pyplot as plt

yCIA = typhon.arts.xml.load("tests/propmat.xml")[0].data[0][0][:, 0]
yHe = typhon.arts.xml.load("tests/propmat.xml")[1].data[0][0][:, 0]
dy = typhon.arts.xml.load("tests/dpropmat.xml")
y = yCIA + yHe
f = np.linspace(3e9, 30e12, len(yCIA)) / 1e9

plt.plot(f, yCIA)
plt.plot(f, yHe)
plt.xlabel("Frequency [GHz]")
plt.ylabel("XSEC [···]")
plt.title("Temperature Derivative")
plt.legend(["CIA-tags", "He-line tag"])
plt.show()

plt.subplot(2, 2, 1)
pyCIA = typhon.arts.xml.load("tests/propmat-dT.xml")[0].data[0][0][:, 0]
pyHe = typhon.arts.xml.load("tests/propmat-dT.xml")[1].data[0][0][:, 0]
py = pyCIA + pyHe
plt.plot(f, dy[0].data[0][0][:, 0], f, (py - y) / 0.0001)
plt.xlabel("Frequency [GHz]")
plt.ylabel("Derivative [···]")
plt.title("Temperature Derivative")

plt.subplot(2, 2, 2)
pyCIA = typhon.arts.xml.load("tests/propmat-dvmr-H2.xml")[0].data[0][0][:, 0]
pyHe = typhon.arts.xml.load("tests/propmat-dvmr-H2.xml")[1].data[0][0][:, 0]
py = pyCIA + pyHe
plt.plot(f, dy[1].data[0][0][:, 0], f, (py - y) / 0.00001)
plt.xlabel("Frequency [GHz]")
plt.ylabel("Derivative [···]")
plt.title("H$_2$ VMR Derivative")

plt.subplot(2, 2, 3)
pyCIA = typhon.arts.xml.load("tests/propmat-dvmr-He.xml")[0].data[0][0][:, 0]
pyHe = typhon.arts.xml.load("tests/propmat-dvmr-He.xml")[1].data[0][0][:, 0]
py = pyCIA + pyHe
plt.plot(f, dy[2].data[0][0][:, 0], f, (py - y) / 0.00001)
plt.xlabel("Frequency [GHz]")
plt.ylabel("Derivative [···]")
plt.title("He VMR Derivative")

plt.subplot(2, 2, 4)
pyCIA = typhon.arts.xml.load("tests/propmat-df.xml")[0].data[0][0][:, 0]
pyHe = typhon.arts.xml.load("tests/propmat-df.xml")[1].data[0][0][:, 0]
py = pyCIA + pyHe
plt.plot(f, dy[3].data[0][0][:, 0], f, (py - y) / 1e3)
plt.xlabel("Frequency [GHz]")
plt.ylabel("Derivative [···]")
plt.title("Frequency Derivative")

plt.tight_layout()
plt.show()
