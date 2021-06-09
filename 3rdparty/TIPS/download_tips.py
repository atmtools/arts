#!/bin/bash
#
# Download TIPS 2021 data from Zenodo.
#
# Sometimes the download fails for no apparent reason. If that happens, wait
# a few seconds and rerun the download script. It will continue from the point
# where the error occurred.
#
# Needs the zenodo_get module: pip install --user zenodo_get

python -m zenodo_get -R 5 -m -o TIPS2021 10.5281/zenodo.4708098

