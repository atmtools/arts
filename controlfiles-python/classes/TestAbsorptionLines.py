import numpy as np
import pyarts.pyarts_cpp as cxx
import test_functions as test

lineshapemodel = cxx.LineShapeModel([cxx.LineShapeSingleSpeciesModel(
    G0=cxx.LineShapeModelParameters("T1", 10000, 0.8))])

line = cxx.AbsorptionSingleLine(F0=1e9,
                                I0=1e-20,
                                glow=5,
                                gupp=7,
                                zeeman=cxx.ZeemanModel(1e-4, 1e-3),
                                localquanta="J 3 2",
                                lineshape=lineshapemodel)

x = cxx.AbsorptionLines(selfbroadening=False,
                        bathbroadening=True,
                        broadeningspecies=["AIR"],
                        quantumidentity="H2O-161",
                        lines=[line],
                        lineshapetype="VP")

assert x.ok

test.io(x, delete=True)

assert x.ok

# Test that line shapes can be computed (Check that the
# Lorentz half-width we compute is actually the half-width)
VMR = [1]
T = 250
P = 1e4
ls = cxx.LineShapeCalculator(x, 0, T, P, VMR)
assert np.isclose(ls.F(line.F0 + x.LineShapeOutput(0, T, P,
                  VMR).G0).real / ls.F(line.F0).real, 0.5)
