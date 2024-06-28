import os
import sys
import urllib.request

nargs = len(sys.argv)

xml = "https://arts.mi.uni-hamburg.de/svn/rt/arts-xml-data/trunk/"
cat = "https://arts.mi.uni-hamburg.de/svn/rt/arts-cat-data/trunk/"

if sys.argv[1] == "xml":
    baseurl = xml
    basedir = "arts-xml-data"
elif sys.argv[1] == "cat":
    baseurl = cat
    basedir = "arts-cat-data"
else:
    baseurl = ""
    basedir = "custom-data"

for file in sys.argv[2:]:
    f = f"{basedir}/{file}"
    if not os.path.exists(f):
        os.makedirs(os.path.split(f)[0], exist_ok=True)
        print(f"Downloading {f}")
        urllib.request.urlretrieve(f"{baseurl}{file}", f)
        print(f"Downloaded {f}")
    else:
        print(f"Skipping {file}, exists at {os.path.abspath(f)}")
