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
from .catalogues import (Sparse
                         )

__all__ = []

classes = {
    'CovarianceMatrix': CovarianceMatrix,
    'GriddedField1': GriddedField1,
    'GriddedField2': GriddedField2,
    'GriddedField3': GriddedField3,
    'GriddedField4': GriddedField4,
    'GriddedField5': GriddedField5,
    'GriddedField6': GriddedField6,
    'ScatteringMetaData': ScatteringMetaData,
    'SingleScatteringData': SingleScatteringData,
    'Sparse': Sparse,
}
