"""
Simple script to recursively run all *.py test in directory tree.
"""
import argparse
import os
import glob
from subprocess import Popen, PIPE, PIPE
from io import StringIO
from tqdm import tqdm

import arts
from arts.parser import convert_to_python

parser = argparse.ArgumentParser(description='Run ARTS test scripts.')
parser.add_argument('basedir', metavar='basedir', type=str, nargs=1,
                   help='Base directory containing the test scripts.')
args = parser.parse_args()
basedir = args.basedir[0]

files = glob.glob(os.path.join(basedir, "**/*.py"), recursive=True)

print("Found {} test files.\n\n".format(len(files)))

failed = []
for f in tqdm(files):

    dirname = os.path.dirname(f)
    p = Popen(["python", f], stdout=PIPE, stderr=PIPE, cwd=dirname)
    stdout, stderr = p.communicate()
    if not p.returncode == 0:
        failed += [f]

print("Passed {} out of {} tests.".format(len(files) - len(failed), len(failed)))
print("\nThe following files failed:\n")
for f in failed:
    print("\t" + f)

