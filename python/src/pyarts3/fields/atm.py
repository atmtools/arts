import numpy as np
import pyarts3 as pa
from copy import deepcopy
from collections.abc import Callable

__all__ = ['from_table']


def from_table(table: np.ndarray, field_names: list[str], convert: dict[str, Callable] = dict(), row_major: bool = True, atm_in=None):
    """
    Convert a table of an atmospheric profile to an atmospheric field.

    The table must contain an 'altitude' column.  The altitude column is
    used to create the gridded data.

    Example:
    >>> import numpy as np
    >>> table = np.array([[1, 2], [3, 4], [5, 6]])  # 3 rows, 2 columns
    >>> field_names = ['temperature', 'pressure', 'altitude']
    >>> convert = {'temperature': lambda x: x * 100}
    >>> from_table(table, field_names, convert)
    {"top_of_atmosphere": -1.7976931348623157e+308 "Base": 2 "SpeciesEnum": 0 "SpeciesIsotope": 188 "QuantumLevelIdentifier": 0 "ScatteringSpeciesProperty": 0}

    Parameters
    ----------
    table : np.ndarray
        The table to convert.  Should be a 2D array x so that x[i] is a parameter
    field_names : list[str]
        The names of the fields to convert.  Should be a list of strings that converts into an atmospheric field key.
    convert : dict[str, Callable]
        A dictionary of conversion functions to apply to each field.  No key means no conversion.
    row_major : bool
        If True, the table is assumed to be in row-major order.  If False, the table is transposed.
    atm_in : pa.arts.AtmField, optional
        An existing atmospheric field to add the new fields to.  If None, a new atmospheric field is created by default construction.
    """

    if atm_in is None:
        atm = pa.arts.AtmField()
    else:
        atm = deepcopy(atm_in)

    if not row_major:
        table = np.array(table).T

    field = pa.arts.GriddedField3(data=np.zeros((table.shape[1], 1, 1)),
                                  grids=[np.zeros((table.shape[1])), [0.], [0.]],
                                  grid_names=['altitude', 'latitude', 'longitude'])

    found = False
    for data, field_name in zip(table, field_names):
        if field_name == 'altitude':
            field.grids = [data, [0.], [0.]]
            found = True
            break
    assert found, "No altitude field found in field_names"

    for data, field_name in zip(table, field_names):
        if field_name != 'altitude':
            key = pa.arts.AtmKey(field_name)
            field.data[:, 0, 0] = convert.get(field_name, lambda x: x)(data)
            field.dataname = field_name
            atm[key] = field

    return atm
