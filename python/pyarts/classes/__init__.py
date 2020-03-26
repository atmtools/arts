# -*- coding: utf-8 -*-

"""This module contains classes within ARTS.
"""

# Full interface including internal manipulation
from pyarts.classes.AbsorptionLines import AbsorptionLines, ArrayOfAbsorptionLines, ArrayOfArrayOfAbsorptionLines
from pyarts.classes.Any import Any
from pyarts.classes.BasicTypes import Numeric, Index, ArrayOfIndex, ArrayOfArrayOfIndex, String, ArrayOfString, ArrayOfArrayOfString
from pyarts.classes.CIARecord import CIARecord, ArrayOfCIARecord
from pyarts.classes.EnergyLevelMap import EnergyLevelMap
from pyarts.classes.GriddedField1 import GriddedField1, ArrayOfGriddedField1, ArrayOfArrayOfGriddedField1
from pyarts.classes.GriddedField2 import GriddedField2, ArrayOfGriddedField2, ArrayOfArrayOfGriddedField2
from pyarts.classes.GriddedField3 import GriddedField3, ArrayOfGriddedField3, ArrayOfArrayOfGriddedField3
from pyarts.classes.GriddedField4 import GriddedField4, ArrayOfGriddedField4
from pyarts.classes.GriddedField5 import GriddedField5
from pyarts.classes.GriddedField6 import GriddedField6
from pyarts.classes.GridPos import GridPos
from pyarts.classes.Matrix import Matrix, ArrayOfMatrix, ArrayOfArrayOfMatrix
from pyarts.classes.Ppath import Ppath, ArrayOfPpath
from pyarts.classes.PropagationMatrix import PropagationMatrix, ArrayOfPropagationMatrix, ArrayOfArrayOfPropagationMatrix
from pyarts.classes.QuantumIdentifier import QuantumIdentifier, ArrayOfQuantumIdentifier
from pyarts.classes.RadiationVector import RadiationVector, ArrayOfRadiationVector, ArrayOfArrayOfRadiationVector
from pyarts.classes.Rational import Rational
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
from pyarts.classes.TransmissionMatrix import TransmissionMatrix, ArrayOfTransmissionMatrix, ArrayOfArrayOfTransmissionMatrix
from pyarts.classes.Vector import Vector, ArrayOfVector, ArrayOfArrayOfVector
from pyarts.classes.Verbosity import Verbosity
from pyarts.classes.XsecRecord import ArrayOfXsecRecord

# No interface but creation and IO
from pyarts.classes.Agenda import Agenda, ArrayOfAgenda
from pyarts.classes.CovarianceMatrix import CovarianceMatrix
from pyarts.classes.GasAbsLookup import GasAbsLookup
from pyarts.classes.MCAntenna import MCAntenna
from pyarts.classes.RetrievalQuantity import ArrayOfRetrievalQuantity
from pyarts.classes.Sparse import Sparse, ArrayOfSparse
from pyarts.classes.TelsemAtlas import TelsemAtlas, ArrayOfTelsemAtlas
from pyarts.classes.Timer import Timer
