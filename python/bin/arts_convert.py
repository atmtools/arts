"""
Simple script to convert .arts controlfiles to Python.
"""
import argparse
import glob
import os
import sys

from pyarts.parser import convert_to_python

parser = argparse.ArgumentParser(description='Convert ARTS controlfile to Python.')
parser.add_argument('--nofail', action='store_true', help='Ignore failed conversions.')
parser.add_argument('--dry-run', action='store_true', help='Simulate conversion.')
parser.add_argument('files',
                    metavar='files',
                    type=str,
                    nargs='+',
                    help='Single filename of base path for recursive conversion.')
parser.add_argument('-o',
                    type=str,
                    default="",
                    help='Filename of base path for output files.')
args = parser.parse_args()
files = args.files
output = args.o

#
# Collect input files
#

if len(files) == 1 and os.path.isdir(files[0]):
    base_path = files[0]
    files = glob.glob(os.path.join(files[0], "**/*.arts"), recursive=True)
else:
    base_path = os.getcwd()

#
# Handle output
#

if len(files) > 1 and not os.path.isdir(output):
    sys.exit("If multiple input files are given output must be a directory.")
if len(files) > 1 and not all([os.path.relpath(f) for f in files]):
    sys.exit("If a directory is given as output, all input files must be given "
             "as relative paths.")

if os.path.isdir(output):

    def make_output(input):
        relpath = os.path.relpath(input, base_path)
        file_out = os.path.join(output, os.path.splitext(relpath)[0] + ".py")
        path_out = os.path.dirname(file_out)
        if not os.path.exists(path_out):
            os.makedirs(path_out)
        return file_out
else:

    def make_output(input):
        return os.path.splitext(input)[0] + ".py"


print("Running arts_convert ...")
converted = 0
for f in files:
    path, filename = os.path.split(f)
    filename_out = make_output(f)
    try:
        if args.dry_run:
            print(f"Would write {filename_out}")
            continue
        convert_to_python(f, filename_out)
        print(f"Wrote {filename_out}")
        converted += 1
    except Exception as e:
        print(f"Error converting {filename}:", e.args[0].splitlines()[-1])

if not args.dry_run and not args.nofail and converted < len(files):
    print(f"arts_convert done. "
          f"Failed to convert {len(files)-converted} out of {len(files)} files.")
    sys.exit(1)
else:
    print(f"arts_convert done. Converted {converted} out of {len(files)} files.")
