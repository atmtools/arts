# -*- coding: utf-8 -*-

"""This module contains classes within ARTS.
"""

# Full interface including internal manipulation
from pyarts.classes.AbsorptionLines import AbsorptionLines, ArrayOfAbsorptionLines, ArrayOfArrayOfAbsorptionLines
from pyarts.classes.Agenda import Agenda, ArrayOfAgenda
from pyarts.classes.Any import Any
from pyarts.classes.BasicTypes import Numeric, Index, ArrayOfIndex, ArrayOfArrayOfIndex, String, ArrayOfString, ArrayOfArrayOfString
from pyarts.classes.CIARecord import CIARecord, ArrayOfCIARecord
from pyarts.classes.EnergyLevelMap import EnergyLevelMap
from pyarts.classes.GasAbsLookup import GasAbsLookup
from pyarts.classes.GriddedField1 import GriddedField1, ArrayOfGriddedField1, ArrayOfArrayOfGriddedField1
from pyarts.classes.GriddedField2 import GriddedField2, ArrayOfGriddedField2, ArrayOfArrayOfGriddedField2
from pyarts.classes.GriddedField3 import GriddedField3, ArrayOfGriddedField3, ArrayOfArrayOfGriddedField3
from pyarts.classes.GriddedField4 import GriddedField4, ArrayOfGriddedField4
from pyarts.classes.GriddedField5 import GriddedField5
from pyarts.classes.GriddedField6 import GriddedField6
from pyarts.classes.GridPos import GridPos
from pyarts.classes.Matrix import Matrix, ArrayOfMatrix, ArrayOfArrayOfMatrix
from pyarts.classes.MCAntenna import MCAntenna
from pyarts.classes.Ppath import Ppath, ArrayOfPpath
from pyarts.classes.PropagationMatrix import PropagationMatrix, ArrayOfPropagationMatrix, ArrayOfArrayOfPropagationMatrix
from pyarts.classes.QuantumIdentifier import QuantumIdentifier, ArrayOfQuantumIdentifier
from pyarts.classes.RadiationVector import RadiationVector, ArrayOfRadiationVector, ArrayOfArrayOfRadiationVector
from pyarts.classes.Rational import Rational
from pyarts.classes.RetrievalQuantity import ArrayOfRetrievalQuantity
from pyarts.classes.ScatteringMetaData import ScatteringMetaData, ArrayOfScatteringMetaData, ArrayOfArrayOfScatteringMetaData
from pyarts.classes.SingleScatteringData import SingleScatteringData, ArrayOfSingleScatteringData, ArrayOfArrayOfSingleScatteringData
from pyarts.classes.SpeciesAuxData import SpeciesAuxData
from pyarts.classes.SpeciesTag import ArrayOfArrayOfSpeciesTag
from pyarts.classes.StokesVector import StokesVector, ArrayOfStokesVector, ArrayOfArrayOfStokesVector
from pyarts.classes.Tensor3 import Tensor3, ArrayOfTensor3, ArrayOfArrayOfTensor3
from pyarts.classes.Tensor4 import Tensor4, ArrayOfTensor4
from pyarts.classes.Tensor5 import Tensor5, ArrayOfTensor5
from pyarts.classes.Tensor6 import Tensor6, ArrayOfTensor6, ArrayOfArrayOfTensor6
from pyarts.classes.Tensor7 import Tensor7, ArrayOfTensor7
from pyarts.classes.TessemNN import TessemNN
from pyarts.classes.Timer import Timer
from pyarts.classes.TransmissionMatrix import TransmissionMatrix, ArrayOfTransmissionMatrix, ArrayOfArrayOfTransmissionMatrix
from pyarts.classes.Vector import Vector, ArrayOfVector, ArrayOfArrayOfVector
from pyarts.classes.Verbosity import Verbosity
from pyarts.classes.XsecRecord import ArrayOfXsecRecord

# No interface but creation and IO
from pyarts.classes.CovarianceMatrix import CovarianceMatrix
from pyarts.classes.Sparse import Sparse, ArrayOfSparse
from pyarts.classes.TelsemAtlas import TelsemAtlas, ArrayOfTelsemAtlas


# Attempt at conversions
import pyarts
import ctypes as c
from pyarts.workspace.api import arts_api as lib

def from_workspace(x):
    """ Converts a workspace variable to a python class

    Input:
        x:
            A workspace variable (pyarts.workspace.WorkspaceVariable)

    Output:
        A modifyable ARTS variable (pyarts.workspace.eval(x.group))
    """
    if isinstance(x, pyarts.workspace.WorkspaceVariable):
        v = lib.get_variable_data_pointer(x.ws.ptr, x.ws_id, x.group_id)
        typ = eval(x.group)

        if v is None:
            raise RuntimeError("Expects non-nullptr")

        return typ(c.c_void_p(v))
    else:
        raise TypeError("Expects WorkspaceVariable")


lib.get_variable_data_pointer.restype = c.c_void_p
lib.get_variable_data_pointer.argtypes = [c.c_void_p, c.c_long]
