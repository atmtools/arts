# PyPI GitHub action

This action builds binary wheels of the pyarts Python package
on a manylinux2014 dist and pushes them to the Python Package
Index (PyPI).

## Input

The action expects the PyPI access token to be provided 
as input.

## Example usage

uses: actions/pypi
with:
pypi_access: ${{secrets.pypi_access}}
