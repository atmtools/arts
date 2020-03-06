"""
Simple script to convert .arts controlfiles to Python.
"""
import argparse
import os
import glob

import arts
from arts.parser import convert_to_python

parser = argparse.ArgumentParser(description='Convert ARTS controlfile to Python.')
parser.add_argument('files', metavar='files', type=str, nargs='+',
                   help='Single filename of base path for recursive conversion.')
args = parser.parse_args()
files = args.files

if len(files) == 1:
    if os.path.isdir(files[0]):
        files = glob.glob(os.path.join(files[0], "**/*.arts"), recursive=True)

for f in files:
    path, filename = os.path.split(f)
    filename_out, _ = os.path.splitext(filename)
    filename_out += ".py"
    f_out = os.path.join(path, filename_out)

    print("Converting {} to {}.".format(filename, filename_out))
    try:
        convert_to_python(f, f_out)
    except Exception as e:
        print("Error converting {}:".format(f))
        print(e)

