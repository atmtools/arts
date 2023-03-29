# -*- coding: utf-8 -*-

"""Collection of all ARTS types."""

from .covariancematrix import (CovarianceMatrix)
from .scattering import (SingleScatteringData,
                         ScatteringMetaData,
                        )
from .catalogues import (Sparse
                         )

__all__ = []

classes = {
    'CovarianceMatrix': CovarianceMatrix,
    'ScatteringMetaData': ScatteringMetaData,
    'SingleScatteringData': SingleScatteringData,
    'Sparse': Sparse,
}
