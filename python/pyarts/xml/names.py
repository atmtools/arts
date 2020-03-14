# -*- coding: utf-8 -*-

__all__ = ['dimension_names', 'tensor_names', 'complex_tensor_names',
           'basic_types']

# Source: ARTS developer guide, section 3.4
dimension_names = [
    'ncols',
    'nrows',
    'npages',
    'nbooks',
    'nshelves',
    'nvitrines',
    'nlibraries']

tensor_names = [
    'Vector', 'Matrix', 'Tensor3', 'Tensor4', 'Tensor5', 'Tensor6', 'Tensor7']

complex_tensor_names = [
    'ComplexVector', 'ComplexMatrix', 'ComplexTensor3', 'ComplexTensor4',
    'ComplexTensor5', 'ComplexTensor6', 'ComplexTensor7']

basic_types = {
    'tuple': 'Array',
    'list': 'Array',
    'int': 'Index',
    'int8': 'Index',
    'int16': 'Index',
    'int32': 'Index',
    'int64': 'Index',
    'float': 'Numeric',
    'float16': 'Numeric',
    'float32': 'Numeric',
    'float64': 'Numeric',
    'float128': 'Numeric',
    'str': 'String',
    'str_': 'String',
    'NoneType': None,
}
