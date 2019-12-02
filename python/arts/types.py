# -*- coding: utf-8 -*-

"""Collection of all ARTS types."""

from .griddedfield import (GriddedField1,
                           GriddedField2,
                           GriddedField3,
                           GriddedField4,
                           GriddedField5,
                           GriddedField6,
                           )
from .covariancematrix import (CovarianceMatrix)
from .scattering import (SingleScatteringData,
                         ScatteringMetaData,
                        )
from .retrieval import (RetrievalQuantity)
from .catalogues import (ArrayOfLineRecord,
                         CIARecord,
                         GasAbsLookup,
                         LineMixingRecord,
                         QuantumIdentifier,
                         QuantumNumberRecord,
                         QuantumNumbers,
                         Sparse,
                         SpeciesAuxData,
                         SpeciesTag,
                         PropagationMatrix,
                         StokesVector,
                         Ppath,
                         GridPos,
                         )
from .internals import (ARTSCAT5,
                        Rational,
                        PressureBroadening,
                        LineMixing,
                        PartitionFunctions,
                        )
from .xsec import (XsecRecord)

__all__ = []

classes = {
    'ArrayOfLineRecord': ArrayOfLineRecord,
    'CIARecord': CIARecord,
    'CovarianceMatrix': CovarianceMatrix,
    'GasAbsLookup': GasAbsLookup,
    'GriddedField1': GriddedField1,
    'GriddedField2': GriddedField2,
    'GriddedField3': GriddedField3,
    'GriddedField4': GriddedField4,
    'GriddedField5': GriddedField5,
    'GriddedField6': GriddedField6,
    'LineMixingRecord': LineMixingRecord,
    'QuantumIdentifier': QuantumIdentifier,
    'QuantumNumberRecord': QuantumNumberRecord,
    'QuantumNumbers': QuantumNumbers,
    'RetrievalQuantity': RetrievalQuantity,
    'ScatteringMetaData': ScatteringMetaData,
    'SingleScatteringData': SingleScatteringData,
    'Sparse': Sparse,
    'SpeciesAuxData': SpeciesAuxData,
    'SpeciesTag': SpeciesTag,
    'ARTSCAT5': ARTSCAT5,
    'Rational': Rational,
    'PressureBroadening': PressureBroadening,
    'LineMixing': LineMixing,
    'PartitionFunctions': PartitionFunctions,
    'XsecRecord': XsecRecord,
    'PropagationMatrix': PropagationMatrix,
    'StokesVector': StokesVector,
    'Ppath': Ppath,
    'GridPos': GridPos,
}
