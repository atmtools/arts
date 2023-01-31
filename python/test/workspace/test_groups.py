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


list_of_groups = [x.name for x in cxx.get_wsv_groups()]
special_groups = ["CallbackFunction", "Any"]


class TestGroups:
    def test_if_test_exist(self):
        for g in list_of_groups:
            assert hasattr(self, f"test{g}"), f"Lacking test for {g}"

    def test_if_group_fulfill_basic_contracts(self):
        ws = Workspace()

        for g in list_of_groups:
            if g in special_groups:
                continue

            if g == "Agenda":
                x = eval("cxx.{}(ws)".format(g))
            else:
                x = eval("cxx.{}()".format(g))
            eval("ws.{}Create('{}__1', '', x)".format(g, g))
            exec("ws.{}__2 = x".format(g))
            assert hasattr(
                x, "readxml"
            ), f"Workspace group {g} lacks an xml reading routine!"
            assert hasattr(
                x, "savexml"
            ), f"Workspace group {g} lacks an xml saving routine!"
            assert x.__doc__ != "Missing description" and len(
                x.__doc__
            ), f"{g} lacks documentation"

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

        # Test that line shapes can be computed (Check that the
        # Lorentz half-width we compute is actually the half-width)
        VMR = [1]
        T = 250
        P = 1e4
        ls = cxx.LineShapeCalculator(x, 0, T, P, VMR)
        assert np.isclose(
            ls.F(line.F0 + x.LineShapeOutput(0, T, P, VMR).G0).real
            / ls.F(line.F0).real,
            0.5,
        )

    def testAgenda(self):
        ws = Workspace()
        ws.x = [4]

        x = cxx.Agenda(ws)

        def test(ws):
            ws.x = [2]
            print("in python:", ws.x.value)
            ws.x = [1]

        # Note that the actual output order here depends on your stream buffer
        # the values counting down is all that matters for
        x.add_workspace_method("Print", ws.x, 0)
        x.add_workspace_method("VectorSet", ws.x, [3])
        x.add_workspace_method("Print", ws.x, 0)
        x.add_callback_method(test)
        x.add_workspace_method("Print", ws.x, 0)
        x.add_workspace_method("Print", "Done!", 0)
        x.name = "test_agenda"
        x.check(ws)
        x.execute(ws)

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
        ws = Workspace()
        y = cxx.Agenda(ws)
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

        x = cxx.ArrayOfArrayOfAbsorptionLines(
            [x])

        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

    def testArrayOfArrayOfGriddedField1(self):
        x = cxx.ArrayOfArrayOfGriddedField1([cxx.ArrayOfGriddedField1([cxx.GriddedField1()])])

        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

    def testArrayOfArrayOfGriddedField2(self):
        x = cxx.ArrayOfArrayOfGriddedField2([cxx.ArrayOfGriddedField2([cxx.GriddedField2()])])

        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

    def testArrayOfArrayOfGriddedField3(self):
        x = cxx.ArrayOfArrayOfGriddedField3([cxx.ArrayOfGriddedField3([cxx.GriddedField3()])])

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

    def testArrayOfArrayOfPropagationMatrix(self):
        x = cxx.ArrayOfArrayOfPropagationMatrix(
            [cxx.ArrayOfPropagationMatrix([cxx.PropagationMatrix()])])

        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

    def testArrayOfArrayOfRadiationVector(self):
        x = cxx.ArrayOfArrayOfRadiationVector(
            [cxx.ArrayOfRadiationVector([cxx.RadiationVector()])])

        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

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

    def testArrayOfArrayOfStokesVector(self):
        x = cxx.ArrayOfArrayOfStokesVector([cxx.ArrayOfStokesVector([cxx.StokesVector()])])

        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

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
        x = cxx.ArrayOfArrayOfTime(1, cxx.ArrayOfTime(1, cxx.Time()))
        test.io(x, delete=True)
        test.array(x)
        test.array_of_array(x)

    def testArrayOfArrayOfTransmissionMatrix(self):
        x = cxx.ArrayOfArrayOfTransmissionMatrix(
            [cxx.ArrayOfTransmissionMatrix([cxx.TransmissionMatrix()])])

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

    def testArrayOfJacobianTarget(self):
        x = cxx.ArrayOfJacobianTarget(1, cxx.JacobianTarget())
        # test.io(x, delete=True)
        test.array(x)

    def testArrayOfMatrix(self):
        x = cxx.ArrayOfMatrix([[[1, 2, 3]]])
        test.io(x, delete=True)
        test.array(x)

        x = cxx.ArrayOfMatrix(np.zeros(shape=(3, 3, 3), dtype=int))
        assert np.all(np.array(x) == 0)

        np.array(x[0], copy=False)[:] = 1
        assert not np.all(np.array(x) == 1)

    def testArrayOfPpath(self):
        x = cxx.ArrayOfPpath(1, cxx.Ppath())
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfPropagationMatrix(self):
        x = cxx.ArrayOfPropagationMatrix([cxx.PropagationMatrix()])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfQuantumIdentifier(self):
        x = cxx.ArrayOfQuantumIdentifier(["H2O-161 J 1 1"])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfRadiationVector(self):
        x = cxx.ArrayOfRadiationVector([cxx.RadiationVector()])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfRetrievalQuantity(self):
        x = cxx.ArrayOfRetrievalQuantity(1, cxx.RetrievalQuantity())
        # test.io(x, delete=True)
        test.array(x)

    def testArrayOfScatteringMetaData(self):
        x = cxx.ArrayOfScatteringMetaData(1, cxx.ScatteringMetaData())
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfSingleScatteringData(self):
        x = cxx.ArrayOfSingleScatteringData(1, cxx.SingleScatteringData())
        # test.io(x, delete=True)
        test.array(x)

    def testArrayOfSparse(self):
        x = cxx.ArrayOfSparse(1, cxx.Sparse())
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfSpeciesTag(self):
        x = cxx.ArrayOfSpeciesTag("H2O")
        x = cxx.ArrayOfSpeciesTag(["H2O"])
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfSun(self):
        ws = Workspace()
        ws.stokes_dim = 1
        ws.f_grid = [1e9, 2e9, 3e9]

        ws.sunsAddSingleBlackbody(
            radius=20, distance=2000, temperature=5000, latitude=10, longitude=45
        )

        sun = ws.suns.value[0]

        assert sun.radius == 20
        assert sun.distance == 2000
        assert sun.latitude == 10
        assert sun.longitude == 45
        assert np.isclose(sun.spectrum[0, 0], 4.82602e-18, atol=1e-25)
        assert np.isclose(sun.spectrum[1, 0], 1.93040e-17, atol=1e-25)
        assert np.isclose(sun.spectrum[2, 0], 4.34338e-17, atol=1e-25)

        x = cxx.ArrayOfSun(1, cxx.Sun())
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfStokesVector(self):
        x = cxx.ArrayOfStokesVector(1, cxx.StokesVector())
        test.io(x, delete=True)
        test.array(x)

    def testArrayOfString(self):
        x = cxx.ArrayOfString(["OI"])
        test.io(x, delete=True)
        test.array(x)

        x = cxx.ArrayOfString(["OI", cxx.String("AI")])

    def testArrayOfTelsemAtlas(self):
        cxx.ArrayOfTelsemAtlas()

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

        x = cxx.ArrayOfTime(1, "2017-01-01 15:30:20")
        test.io(x, delete=True)
        test.array(x)

        x = cxx.ArrayOfTime([f"2020-01-{x} 00:00:00" for x in range(10, 30)])

        assert x.as_datetime[-1] == datetime.datetime(2020, 1, 29, 0, 0)

    def testArrayOfTransmissionMatrix(self):
        x = cxx.ArrayOfTransmissionMatrix([cxx.TransmissionMatrix()])
        test.io(x, delete=True)
        test.array(x)

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
        x = cxx.ArrayOfXsecRecord(1, cxx.XsecRecord())
        test.io(x, delete=True)
        test.array(x)

    def testCallbackFunction(self):
        def oi(ws):
            print("oi")
            ws.atmosphere_dim = 3

        x = cxx.CallbackFunction(oi)

        ws = Workspace()
        assert not ws.atmosphere_dim.init
        x(ws)
        assert ws.atmosphere_dim.init
        assert ws.atmosphere_dim.value.value == 3

    def testCIARecord(self):
        x = cxx.CIARecord()
        test.io(x, delete=True)

    def testCovarianceMatrix(self):
        x = cxx.CovarianceMatrix()
        test.io(x, delete=True)

    def testEnergyLevelMap(self):
        x = cxx.EnergyLevelMap()
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

    def testGridPos(self):
        x = cxx.GridPos()
        test.io(x, delete=True)

    def testHitranRelaxationMatrixData(self):
        x = cxx.HitranRelaxationMatrixData()

        test.io(x, delete=True)

    def testIndex(self):
        x = cxx.Index(0)
        test.io(x, delete=True)

        assert x.value == 0
        x.value = x + 3
        assert x.value == 3

    def testJacobianTarget(self):
        cxx.JacobianTarget()
        # test.io(x, delete=True)

    def testMapOfErrorCorrectedSuddenData(self):
        x = cxx.MapOfErrorCorrectedSuddenData()
        test.io(x, delete=True)

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

    def testMCAntenna(self):
        cxx.MCAntenna()
        # test.io(x, delete=True)

    def testNumeric(self):
        x = cxx.Numeric(0)
        test.io(x, delete=True)

        assert x.value == 0
        x.value = x + 2
        assert x.value == 2

    def testPpath(self):
        x = cxx.Ppath()
        test.io(x, delete=True)

    def testPredefinedModelData(self):
        x = cxx.PredefinedModelData()
        test.io(x, delete=True)

    def testPropagationMatrix(self):
        x = cxx.PropagationMatrix(100, 4, 4, 4)
        test.io(x, delete=True)

        assert np.all(np.array(x.data) == 0)

        np.array(x.data)[:] = 1
        assert np.all(np.array(x.data) == 0)

        np.array(x.data, copy=False)[:] = 1
        assert np.all(np.array(x.data) == 1)

    def testQuantumIdentifier(self):
        x = cxx.QuantumIdentifier("O2-66 v 0 0")
        test.io(x, delete=True)

    def testRadiationVector(self):
        x = cxx.RadiationVector(10, 4)
        test.io(x, delete=True)

        assert np.all(x[0] == 0)

        x[0][:] = 2
        assert np.all(x[0] == 2)

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

    def testSingleScatteringData(self):
        cxx.SingleScatteringData()
        # test.io(x, delete=True)

    def testSparse(self):
        x = cxx.Sparse()
        test.io(x, delete=True)

    def testSpeciesIsotopologueRatios(self):
        x = cxx.SpeciesIsotopologueRatios()
        test.io(x, delete=True)

    def testStokesVector(self):
        x = cxx.StokesVector(50, 4)
        test.io(x, delete=True)

        assert np.all(np.array(x.data) == 0)

        np.array(x.data)[:] = 1
        assert np.all(np.array(x.data) == 0)

        np.array(x.data, copy=False)[:] = 1
        assert np.all(np.array(x.data) == 1)

    def testString(self):
        x = cxx.String("ho")
        test.io(x, delete=True)

        assert "ho" == x
        assert "h" == x[0]

        x[1] = "i"
        assert "hi" == x
        assert hash("hi") == hash(x)

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

    def testTessemNN(self):
        cxx.TessemNN()
        # test.io(x, delete=True)

    def testTime(self):
        x = cxx.Time()
        test.io(x, delete=True)

    def testTimer(self):
        cxx.Timer()

    def testTransmissionMatrix(self):
        x = cxx.TransmissionMatrix(10, 4)
        test.io(x, delete=True)

        assert np.all(x[0] == np.diag([1, 1, 1, 1]))

        x[0][:] = 2
        assert np.all(x[0] == 2)

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

    def testVerbosity(self):
        cxx.Verbosity()
    
    def testAtmField(self):
        x = cxx.AtmField()
        x["wind_u"] = 3.
        x[cxx.ArrayOfSpeciesTag("O2")] = cxx.GriddedField3(
            [[1, 2, 3], [2, 3, 4, 5, 6, 7, 8], [5, 6, 7]], np.random.rand(3, 7, 3))

        def test_fun(a, b, c): return float(a+b+c)
        x[cxx.ArrayOfSpeciesTag("N2")] = test_fun
        test.io(x, delete=True)
        x[cxx.ArrayOfSpeciesTag("N2")] = test_fun
        x.top_of_atmosphere = 2.5
        x.regularize([1, 2, 3], [2, 3, 4, 5, 6, 7, 8], [5, 6, 7])
        test.io(x, delete=True)
        assert [3, 7, 3] == x.regularized_shape()
        
    def test_pickle(self):
        ws = Workspace()
        x = list(cxx.get_WsvGroupMap().keys())

        for i in range(len(x)):
            if x[i] == "CallbackFunction":
                continue
            elif x[i] == "Agenda":
                continue  # Cannot unpickle a standalone Agenda
            elif x[i] == "ArrayOfAgenda":
                continue  # Cannot unpickle a standalone ArrayOfAgenda
            else:
                exec(f"ws.v{i} = cxx.{x[i]}()")

            print(f"Pickling the workspace after adding a {x[i]}")

            with open("test.pcl", "wb") as f:
                pickle.dump(ws, f)

            with open("test.pcl", "rb") as f:
                ws2 = pickle.load(f)

        assert ws.number_of_initialized_variables() == ws2.number_of_initialized_variables(), \
            "Must be able to fully init a workspace containing no CallbackFunction"

        ws.testing = cxx.Index(3)

        @arts_agenda(ws=ws, set_agenda=True)
        def test_agenda(ws):
            ws.Print(ws.testing, 0)

        with open("ws.pcl", "wb") as f:
            pickle.dump(ws, f)

        with open("ws.pcl", "rb") as f:
            ws2 = pickle.load(f)

        print(ws)
        print()
        print(ws2)

        print("Predicted working output:\n\n3\n3\n3\n2\n\nActual output (if streams are synchronized):\n")
        ws.test_agenda.value.execute(ws)
        ws2.test_agenda.value.execute(ws2)
        ws2.testing = 2
        ws.test_agenda.value.execute(ws)
        ws2.test_agenda.value.execute(ws2)


if __name__ == "__main__":
    x = TestGroups()
    b = x.testAtmField()
