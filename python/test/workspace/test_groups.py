"""
Test handling of workspace groups of the Python interface.

Each group has to have a test{GROUP} function in the main test class
"""


import os
import uuid
import pickle

import numpy as np
import pyarts.arts as cxx
from pyarts.workspace import Workspace, arts_agenda
from copy import copy, deepcopy


class test:
    @staticmethod
    def io(x, fname=None, delete=False):
        if fname is None:
            fname = str(uuid.uuid4()) + ".xml"
        try:
            x.savexml(fname)
            x.readxml(fname)
        finally:
            if os.path.exists(fname) and delete:
                os.remove(fname)
        return True

    @staticmethod
    def array(x):
        assert len(x) == 1
        x.append(x[0])
        assert len(x) == 2
        x[0] = x[1]
        x.pop()
        assert len(x) == 1
        y = type(x)([x[0]])
        assert len(y) == 1
        return True

    @staticmethod
    def array_of_array(x):
        assert len(x) > 0
        y = type(x)([[x[0][0]]])
        assert len(y) == 1
        assert len(y[0]) == 1
        return True

    @staticmethod
    def shape_match(x, y):
        x = x.shape
        assert len(x) == len(y)
        for i in range(len(x)):
            assert x[i] == y[i]
        return True


list_of_groups = list(cxx.globals.workspace_groups().keys())
special_groups = ["CallbackFunction", "Any", "SpectralRadianceProfileOperator"]


class TestGroups:
    def testAbsorptionLines(self):
        lineshapemodel = cxx.LineShapeModel(
            [
                cxx.LineShapeSingleSpeciesModel(
                    G0=cxx.LineShapeModelParameters("T1", 10000, 0.8)
                )
            ]
        )

        line = cxx.AbsorptionSingleLine(
            F0=1e9,
            I0=1e-20,
            glow=5,
            gupp=7,
            zeeman=cxx.ZeemanModel(1e-4, 1e-3),
            localquanta="J 3 2",
            lineshape=lineshapemodel,
        )

        x = cxx.AbsorptionLines(
            selfbroadening=False,
            bathbroadening=True,
            broadeningspecies=["AIR"],
            quantumidentity="H2O-161",
            lines=[line],
            lineshapetype="VP",
        )

        assert x.ok

        test.io(x, delete=True)

        assert x.ok

    def testAny(self):
        cxx.Any()

    def testArrayOfAbsorptionLines(self):
        lineshapemodel = cxx.LineShapeModel(
            [
                cxx.LineShapeSingleSpeciesModel(
                    G0=cxx.LineShapeModelParameters("T1", 10000, 0.8)
                )
            ]
        )

        line = cxx.AbsorptionSingleLine(
            F0=1e9,
            I0=1e-20,
            glow=5,
            gupp=7,
            zeeman=cxx.ZeemanModel(1e-4, 1e-3),
            localquanta="J 3 2",
            lineshape=lineshapemodel,
        )

        x = cxx.AbsorptionLines(
            selfbroadening=False,
            bathbroadening=True,
            broadeningspecies=["AIR"],
            quantumidentity="H2O-161",
            lines=[line],
            lineshapetype="VP",
        )

        x = cxx.ArrayOfAbsorptionLines([x])

        test.io(x, delete=True)
        test.array(x)

    def testArrayOfAgenda(self):
        y = cxx.Agenda()
        x = cxx.ArrayOfAgenda([y])

        test.array(x)

    def testArrayOfArrayOfAbsorptionLines(self):
        lineshapemodel = cxx.LineShapeModel(
            [
                cxx.LineShapeSingleSpeciesModel(
                    G0=cxx.LineShapeModelParameters("T1", 10000, 0.8)
                )
            ]
        )

        line = cxx.AbsorptionSingleLine(
            F0=1e9,
            I0=1e-20,
            glow=5,
            gupp=7,
            zeeman=cxx.ZeemanModel(1e-4, 1e-3),
            localquanta="J 3 2",
            lineshape=lineshapemodel,
        )

        x = cxx.AbsorptionLines(
            selfbroadening=False,
            bathbroadening=True,
            broadeningspecies=["AIR"],
            quantumidentity="H2O-161",
            lines=[line],
            lineshapetype="VP",
        )

        x = cxx.ArrayOfAbsorptionLines([x])

        x = cxx.ArrayOfArrayOfAbsorptionLines([x])

        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

    def testArrayOfArrayOfGriddedField1(self):
        x = cxx.ArrayOfArrayOfGriddedField1(
            [cxx.ArrayOfGriddedField1([cxx.GriddedField1()])]
        )

        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

    def testArrayOfArrayOfGriddedField2(self):
        x = cxx.ArrayOfArrayOfGriddedField2(
            [cxx.ArrayOfGriddedField2([cxx.GriddedField2()])]
        )

        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

    def testArrayOfArrayOfGriddedField3(self):
        x = cxx.ArrayOfArrayOfGriddedField3(
            [cxx.ArrayOfGriddedField3([cxx.GriddedField3()])]
        )

        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

    def testArrayOfArrayOfIndex(self):
        x = cxx.ArrayOfArrayOfIndex([[1, 2, 3]])
        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

        x = cxx.ArrayOfArrayOfIndex([[1, 2, 3], [1, 2]])
        x = cxx.ArrayOfArrayOfIndex([np.array([1, 2, 3]), [1, 2]])
        x = cxx.ArrayOfArrayOfIndex(np.zeros(shape=(3, 6), dtype=int))
        assert np.all(np.array(x) == 0)

        np.array(x[0], copy=False)[:] = 1
        assert not np.all(np.array(x) == 0)

    def testArrayOfArrayOfMatrix(self):
        x = cxx.ArrayOfArrayOfMatrix([[[[1, 2, 3]]]])
        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

        x = cxx.ArrayOfArrayOfMatrix(np.zeros(shape=(3, 3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x[0][0])[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x[0][0], copy=False)[:] = 1
        assert not np.all(np.array(x) == 0)

    def testArrayOfArrayOfScatteringMetaData(self):
        x = cxx.ArrayOfArrayOfScatteringMetaData()
        test.io(x, delete=True)

    def testArrayOfArrayOfSingleScatteringData(self):
        x = cxx.ArrayOfArrayOfSingleScatteringData()
        test.io(x, delete=True)

    def testArrayOfArrayOfSpeciesTag(self):
        x = cxx.ArrayOfArrayOfSpeciesTag(["H2O,H2O-PWR98"])

        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

        x = cxx.ArrayOfArrayOfSpeciesTag(["H2O", "H2O-PWR98"])
        assert len(x) == 2

    def testArrayOfArrayOfString(self):
        x = cxx.ArrayOfArrayOfString([["OI"]])
        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

    def testArrayOfArrayOfTensor3(self):
        x = cxx.ArrayOfArrayOfTensor3([[[[[1, 2, 3]]]]])
        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

        x = cxx.ArrayOfArrayOfTensor3(np.zeros(shape=(3, 3, 3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x[0][0])[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x[0][0], copy=False)[:] = 1
        assert not np.all(np.array(x) == 0)

    def testArrayOfArrayOfTensor6(self):
        x = cxx.ArrayOfArrayOfTensor6([[[[[[[[1, 2, 3]]]]]]]])
        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

        x = cxx.ArrayOfArrayOfTensor6(np.zeros(shape=(3, 3, 3, 3, 3, 3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x[0][0])[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x[0][0], copy=False)[:] = 1
        assert not np.all(np.array(x) == 0)

    def testArrayOfArrayOfTime(self):
        x = cxx.ArrayOfArrayOfTime([cxx.ArrayOfTime([cxx.Time()])])
        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

    def testArrayOfArrayOfVector(self):
        x = cxx.ArrayOfArrayOfVector([[[1, 2, 3]]])
        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

        x = cxx.ArrayOfArrayOfVector(np.zeros(shape=(3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x[0][0])[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x[0][0], copy=False)[:] = 1
        assert not np.all(np.array(x) == 0)

    def testArrayOfCIARecord(self):
        x = cxx.ArrayOfCIARecord()
        test.io(x, delete=True)

    def testArrayOfGriddedField1(self):
        x = cxx.ArrayOfGriddedField1([cxx.GriddedField1()])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfGriddedField2(self):
        x = cxx.ArrayOfGriddedField2([cxx.GriddedField2()])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfGriddedField3(self):
        x = cxx.ArrayOfGriddedField3([cxx.GriddedField3()])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfGriddedField4(self):
        x = cxx.ArrayOfGriddedField4([cxx.GriddedField4()])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfIndex(self):
        x = cxx.ArrayOfIndex([1])
        test.io(x, delete=True)
        test.array(x)

        x = cxx.ArrayOfIndex([1, 2, 3, 4])
        x = cxx.ArrayOfIndex(np.zeros(shape=(5), dtype=int))
        assert np.all(np.array(x) == 0)

        np.array(x, copy=False)[:] = 1
        assert np.all(np.array(x) == 1)

    def testArrayOfMatrix(self):
        x = cxx.ArrayOfMatrix([[[1, 2, 3]]])
        test.io(x, delete=True)
        test.array(x)

        x = cxx.ArrayOfMatrix(np.zeros(shape=(3, 3, 3), dtype=int))
        assert np.all(np.array(x) == 0)

        np.array(x[0], copy=False)[:] = 1
        assert not np.all(np.array(x) == 1)

    def testArrayOfQuantumIdentifier(self):
        x = cxx.ArrayOfQuantumIdentifier(["H2O-161 J 1 1"])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfScatteringMetaData(self):
        x = cxx.ArrayOfScatteringMetaData([cxx.ScatteringMetaData()])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfSingleScatteringData(self):
        x = cxx.ArrayOfSingleScatteringData([cxx.SingleScatteringData()])
        # test.io(x, delete=True)
        test.array(x)

    def testArrayOfSparse(self):
        x = cxx.ArrayOfSparse([cxx.Sparse()])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfSpeciesTag(self):
        x = cxx.ArrayOfSpeciesTag("H2O")
        x = cxx.ArrayOfSpeciesTag(["H2O"])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfSun(self):
        ws = Workspace()
        ws.frequency_grid = [1e9, 2e9, 3e9]

        x = cxx.ArrayOfSun([cxx.Sun()])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfString(self):
        x = cxx.ArrayOfString(["OI"])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfTensor3(self):
        x = cxx.ArrayOfTensor3([[[[1, 2, 3]]]])
        test.io(x, delete=True)
        test.array(x)

        x = cxx.ArrayOfTensor3(np.zeros(shape=(3, 3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x[0])[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x[0], copy=False)[:] = 1
        assert not np.all(np.array(x) == 0)

    def testArrayOfTensor4(self):
        x = cxx.ArrayOfTensor4([[[[[1, 2, 3]]]]])
        test.io(x, delete=True)
        test.array(x)

        x = cxx.ArrayOfTensor4(np.zeros(shape=(3, 3, 3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x[0])[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x[0], copy=False)[:] = 1
        assert not np.all(np.array(x) == 0)

    def testArrayOfTensor5(self):
        x = cxx.ArrayOfTensor5([[[[[[1, 2, 3]]]]]])
        test.io(x, delete=True)
        test.array(x)

        x = cxx.ArrayOfTensor5(np.zeros(shape=(3, 3, 3, 3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x[0])[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x[0], copy=False)[:] = 1
        assert not np.all(np.array(x) == 0)

    def testArrayOfTensor6(self):
        x = cxx.ArrayOfTensor6([[[[[[[1, 2, 3]]]]]]])
        test.io(x, delete=True)
        test.array(x)

        x = cxx.ArrayOfTensor6(np.zeros(shape=(3, 3, 3, 3, 3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x[0])[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x[0], copy=False)[:] = 1
        assert not np.all(np.array(x) == 0)

    def testArrayOfTensor7(self):
        x = cxx.ArrayOfTensor7([[[[[[[[1, 2, 3]]]]]]]])
        test.io(x, delete=True)
        test.array(x)

        x = cxx.ArrayOfTensor7(np.zeros(shape=(3, 3, 3, 3, 3, 3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x[0])[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x[0], copy=False)[:] = 1
        assert not np.all(np.array(x) == 0)

    def testArrayOfTime(self):
        import datetime as datetime

        x = cxx.ArrayOfTime(["2017-01-01 15:30:20"])
        test.io(x, delete=True)
        test.array(x)

        x = cxx.ArrayOfTime([f"2020-01-{x} 00:00:00" for x in range(10, 30)])

        assert x.as_datetime[-1] == datetime.datetime(2020, 1, 29, 0, 0)

    def testArrayOfVector(self):
        x = cxx.ArrayOfVector([[1, 2, 3]])
        test.io(x, delete=True)
        test.array(x)

        x = cxx.ArrayOfVector(np.zeros(shape=(3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x[0])[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x[0], copy=False)[:] = 1
        assert not np.all(np.array(x) == 0)

    def testArrayOfXsecRecord(self):
        x = cxx.ArrayOfXsecRecord([cxx.XsecRecord()])
        test.io(x, delete=True)
        test.array(x)

    def testCIARecord(self):
        x = cxx.CIARecord()
        test.io(x, delete=True)

    def testCovarianceMatrix(self):
        x = cxx.CovarianceMatrix()
        test.io(x, delete=True)

    def testGasAbsLookup(self):
        x = cxx.GasAbsLookup()
        test.io(x, delete=True)

    def testGriddedField1(self):
        x = cxx.GriddedField1()
        test.io(x, delete=True)

    def testGriddedField2(self):
        x = cxx.GriddedField2()
        test.io(x, delete=True)

    def testGriddedField3(self):
        x = cxx.GriddedField3()
        test.io(x, delete=True)

    def testGriddedField4(self):
        x = cxx.GriddedField4()
        test.io(x, delete=True)

    def testGriddedField5(self):
        x = cxx.GriddedField5()
        test.io(x, delete=True)

    def testGriddedField6(self):
        x = cxx.GriddedField6()
        test.io(x, delete=True)

    def testIndex(self):
        x = cxx.Index(0)
        test.io(x, delete=True)

        assert x == 0
        x = x + 3
        assert x == 3

    def testMatrix(self):
        x = cxx.Matrix([[1, 2, 3]])
        test.io(x, delete=True)

        y = cxx.Matrix([[1], [2], [3]])
        test.shape_match(x @ y, [1, 1])
        test.shape_match(y @ x, [3, 3])

        z = cxx.Vector([1, 2, 3])
        test.shape_match((y @ x) @ z, [3])

        x = cxx.Matrix(np.zeros(shape=(3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x)[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x, copy=False)[:] = 1
        assert np.all(np.array(x) == 1)

        x += 1
        assert np.all(np.array(x) == 2)

        x *= 2
        assert np.all(np.array(x) == 4)

    def testNumeric(self):
        x = cxx.Numeric(0)
        test.io(x, delete=True)

        assert x == 0.0, f"{x} cannot evaluate as equal to 0."
        assert x == 0, f"{x} cannot evaluate as equal to 0"
        x = x + 2
        assert x == 2

    def testPredefinedModelData(self):
        x = cxx.PredefinedModelData()
        test.io(x, delete=True)

    def testQuantumIdentifier(self):
        x = cxx.QuantumIdentifier("O2-66 v 0 0")
        test.io(x, delete=True)

    def testRational(self):
        x = cxx.Rational()
        test.io(x, delete=True)

        y = x + 1
        x += 1
        assert y == x

        y = x * 2
        x *= 2
        assert y == x

        y = x / 3
        x /= 3
        assert y == x

        y = x - 4
        x -= 4
        assert y == x

        assert not (x != y)

        assert x <= y

        assert x >= y

        assert x < 0

        assert x > -4

    def testScatteringMetaData(self):
        x = cxx.ScatteringMetaData()
        test.io(x, delete=True)

    def testSparse(self):
        x = cxx.Sparse()
        test.io(x, delete=True)

    def testString(self):
        x = cxx.String("ho")
        test.io(x, delete=True)

        assert "ho" == x
        assert "h" == x[0]
        # assert hash("ho") == hash(x) - FIXME: breaks in nanobind

    def testTelsemAtlas(self):
        x = cxx.TelsemAtlas()
        test.io(x, delete=True)

    def testTensor3(self):
        x = cxx.Tensor3([[[1, 2, 3]]])
        test.io(x, delete=True)

        x = cxx.Tensor3(np.zeros(shape=(3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x)[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x, copy=False)[:] = 1
        assert np.all(np.array(x) == 1)

        x += 1
        assert np.all(np.array(x) == 2)

        x *= 2
        assert np.all(np.array(x) == 4)

    def testTensor4(self):
        x = cxx.Tensor4([[[[1, 2, 3]]]])
        test.io(x, delete=True)

        x = cxx.Tensor4(np.zeros(shape=(3, 3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x)[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x, copy=False)[:] = 1
        assert np.all(np.array(x) == 1)

        x += 1
        assert np.all(np.array(x) == 2)

        x *= 2
        assert np.all(np.array(x) == 4)

    def testTensor5(self):
        x = cxx.Tensor5([[[[[1, 2, 3]]]]])
        test.io(x, delete=True)

        x = cxx.Tensor5(np.zeros(shape=(3, 3, 3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x)[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x, copy=False)[:] = 1
        assert np.all(np.array(x) == 1)

        x += 1
        assert np.all(np.array(x) == 2)

        x *= 2
        assert np.all(np.array(x) == 4)

    def testTensor6(self):
        x = cxx.Tensor6([[[[[[1, 2, 3]]]]]])
        test.io(x, delete=True)

        x = cxx.Tensor6(np.zeros(shape=(3, 3, 3, 3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x)[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x, copy=False)[:] = 1
        assert np.all(np.array(x) == 1)

        x += 1
        assert np.all(np.array(x) == 2)

        x *= 2
        assert np.all(np.array(x) == 4)

    def testTensor7(self):
        x = cxx.Tensor7([[[[[[[1, 2, 3]]]]]]])
        test.io(x, delete=True)

        x = cxx.Tensor7(np.zeros(shape=(3, 3, 3, 3, 3, 3, 3)))
        assert np.all(np.array(x) == 0)

        np.array(x)[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x, copy=False)[:] = 1
        assert np.all(np.array(x) == 1)

        x += 1
        assert np.all(np.array(x) == 2)

        x *= 2
        assert np.all(np.array(x) == 4)

    def testTime(self):
        x = cxx.Time()
        test.io(x, delete=True)

    def testVector(self):
        x = cxx.Vector([1, 2, 3])
        test.io(x, delete=True)

        x = cxx.Vector(np.zeros(shape=(3)))
        assert np.all(np.array(x) == 0)

        np.array(x)[:] = 1
        assert np.all(np.array(x) == 0)

        np.array(x, copy=False)[:] = 1
        assert np.all(np.array(x) == 1)

        x += 1
        assert np.all(np.array(x) == 2)

        x *= 2
        assert np.all(np.array(x) == 4)

        assert x @ x == 48

    def testAtmField(self):
        x = cxx.AtmField()
        x["wind_u"] = 3.0
        x[cxx.SpeciesEnum("O2")] = cxx.GriddedField3(
            grids=[[1, 2, 3], [2, 3, 4, 5, 6, 7, 8], [5, 6, 7]],
            data=np.random.rand(3, 7, 3),
        )

        def test_fun(a, b, c):
            return float(a + b + c)

        x[cxx.SpeciesEnum("N2")] = test_fun
        test.io(x, delete=True)
        x[cxx.SpeciesEnum("N2")] = test_fun
        x.top_of_atmosphere = 2.5
        test.io(x, delete=True)

    def testAtmPoint(self):
        x = cxx.AtmPoint()
        x["wind_u"] = 3.0
        x[cxx.SpeciesEnum("O2")] = 0.21
        x[cxx.SpeciesEnum("N2")] = 0.79
        test.io(x, delete=True)

    def testVibrationalEnergyLevels(self):
        x = cxx.VibrationalEnergyLevels({"H2O-161": 3, "O2-66": 0.21})
        x["H2O-162"] = 5.0
        x["O2-68"] = 0.0021
        test.io(x, delete=True)

    def testNumericUnaryOperator(self):
        def sqrt(n):
            return np.sqrt(n)

        f = cxx.NumericUnaryOperator(sqrt)
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        assert np.allclose(np.sqrt(x) - f(x), 0)

    def testNumericTernaryOperator(self):
        def add(x, y, z):
            return x + y + z

        f = cxx.NumericTernaryOperator(add)
        x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
        assert np.allclose(x + x + x - f(x, x, x), 0)

    def test_xml(self):
        ignore_groups = [
            "CallbackOperator",
            "NumericUnaryOperator",
            "NumericTernaryOperator",
            "SingleScatteringData",
            "JacobianTargets",
            "SpectralRadianceOperator",
            "DisortBDRF",
            "MatrixOfDisortBDRF",
        ]

        groups = list(cxx.globals.workspace_groups().keys())
        groups.sort()

        for group in ignore_groups:
            if group not in groups:
                raise Exception(f"Ignored group {group} is not in workspace")

        fail = []
        for group in groups:
            if group in ignore_groups:
                print(f"Skipping group {group}")
                continue

            try:
                t = eval(f"cxx.{group}")
                print(f"Testing group {group}")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, does not exist:\n{e}\n\n")
                continue

            try:
                x = t()
                print(f"  Can init as {group}()")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, cannot init:\n{e}\n\n")
                continue

            try:
                x.savexml("test.xml")
                print(f"  Can save as xml ({group}.savexml())")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, cannot save:\n{e}\n\n")
                continue

            try:
                x.readxml("test.xml")
                print(f"  Can read as xml ({group}.readxml())")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, cannot read:\n{e}\n\n")
                continue

            try:
                t.fromxml("test.xml")
                print(f"  Can init ax xml ({group}.fromxml())")
            except Exception as e:
                fail.append(
                    f"\n\nFailed {group}, cannot init from file:\n{e}\n\n"
                )
                continue

        if len(fail) > 0:
            for g in fail:
                print(g)
            raise Exception(f"There are {len(fail)} tests failing")

    def test_print(self):
        ignore_groups = []

        groups = list(cxx.globals.workspace_groups().keys())
        groups.sort()

        for group in ignore_groups:
            if group not in groups:
                raise Exception(f"Ignored group {group} is not in workspace")

        fail = []
        for group in groups:
            if group in ignore_groups:
                print(f"Skipping group {group}")
                continue

            try:
                t = eval(f"cxx.{group}")
                print(f"Testing group {group}")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, does not exist:\n{e}\n\n")
                continue

            try:
                x = t()
                print(f"  Can init as {group}()")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, cannot init:\n{e}\n\n")
                continue

            try:
                print(x)
                print(f"  Can print (print({group}()))")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, cannot print:\n{e}\n\n")
                continue

            if "\0" in str(x):
                fail.append(f"\n\nFailed {group}, \\0 in str()")

            if "\0" in repr(x):
                fail.append(f"\n\nFailed {group}, \\0 in repr()")

            if "\0" in format(x):
                fail.append(f"\n\nFailed {group}, \\0 in format()")

            try:
                assert isinstance(
                    repr(x), str
                ), f"repr not a str but {type(repr(x))}"
                print(f"  Can repr (isinstance(repr({group}()), str))")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, cannot repr:\n{e}\n\n")
                continue

        if len(fail) > 0:
            for g in fail:
                print(g)
            raise Exception(f"There are {len(fail)} tests failing")

    def test_copy(self):
        ignore_groups = ["NumericUnaryOperator", "NumericTernaryOperator"]

        groups = list(cxx.globals.workspace_groups().keys())
        groups.sort()

        for group in ignore_groups:
            if group not in groups:
                raise Exception(f"Ignored group {group} is not in workspace")

        fail = []
        for group in groups:
            if group in ignore_groups:
                print(f"Skipping group {group}")
                continue

            try:
                t = eval(f"cxx.{group}")
                print(f"Testing group {group}")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, does not exist:\n{e}\n\n")
                continue

            try:
                x = t()
                print(f"  Can init as {group}()")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, cannot init:\n{e}\n\n")
                continue

            try:
                y = t(x)
                print(f"  Can copy init (x = {group}({group}()))")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, cannot copy init:\n{e}\n\n")
                continue

            try:
                y = copy(x)
                print(f"  Can copy init (x = copy({group}()))")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, cannot copy:\n{e}\n\n")
                continue

            try:
                y = deepcopy(x)
                print(f"  Can copy init (x = deepcopy({group}()))")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, cannot deepcopy:\n{e}\n\n")
                continue

        if len(fail) > 0:
            for g in fail:
                print(g)
            raise Exception(f"There are {len(fail)} tests failing")

    def test_construct_empty(self):
        ignore_groups = []

        groups = list(cxx.globals.workspace_groups().keys())
        groups.sort()

        fail = []
        for group in groups:
            if group in ignore_groups:
                print(f"Skipping group {group}")
                continue

            try:
                t = eval(f"cxx.{group}")
                print(f"Testing group {group}")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, does not exist:\n{e}\n\n")
                continue

            try:
                x = t()
                print(f"  Can init as {group}()")
            except Exception as e:
                fail.append(f"\n\nFailed {group}, cannot init:\n{e}\n\n")
                continue

        if len(fail) > 0:
            for g in fail:
                print(g)
            raise Exception(f"There are {len(fail)} tests failing")


if __name__ == "__main__":
    x = TestGroups()
    x.test_print()
