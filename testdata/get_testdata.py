import logging
import os
import sys
import time
import urllib.request
import urllib.error

if os.environ.get("DEBUG_DOWNLOADING") is not None:
    import http.client

    http.client.HTTPConnection.debuglevel = 1

    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s %(levelname)s %(message)s"
    )
else:
    logging.basicConfig(level=logging.INFO, format="%(message)s")


def download_file(url, destination, nretry=5):
    """Download a file from a URL to a destination file."""
    backoff_factor = 1
    status_forcelist = [500, 502, 503, 504]

    for attempt in range(nretry + 1):
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req) as response:
                with open(destination, "wb") as file:
                    while True:
                        chunk = response.read(8192)
                        if not chunk:
                            break
                        file.write(chunk)
            return True
        except urllib.error.HTTPError as err:
            if err.code in status_forcelist and attempt < nretry:
                time.sleep(backoff_factor * (2 ** attempt))
                continue
            logging.error(f"Download failed: {err}")
            return False
        except urllib.error.URLError as err:
            if attempt < nretry:
                time.sleep(backoff_factor * (2 ** attempt))
                continue
            logging.error(f"Download failed: {err}")
            return False
        except Exception as err:
            logging.error(f"Download failed: {err}")
            return False
    return False


xml = "https://gitlab.rrz.uni-hamburg.de/atmtools/arts-xml-data/-/raw/main/"
cat = "https://gitlab.rrz.uni-hamburg.de/atmtools/arts-cat-data/-/raw/main/"
xml_fallback = "https://arts.mi.uni-hamburg.de/svn/rt/arts-xml-data/trunk/"
cat_fallback = "https://arts.mi.uni-hamburg.de/svn/rt/arts-cat-data/trunk/"

file_list_path = sys.argv[1]
filename = os.path.basename(file_list_path)

if filename.startswith("xml"):
    baseurl = xml
    baseurl_fallback = xml_fallback
    basedir = "arts-xml-data"
elif filename.startswith("cat"):
    baseurl = cat
    baseurl_fallback = cat_fallback
    basedir = "arts-cat-data"
else:
    raise RuntimeError(f"Unknown file list {file_list_path}")

if len(sys.argv) > 2:
    basedir = os.path.join(sys.argv[2], basedir)

with open(file_list_path, "r") as f:
    files = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]

for file in files:
    f = f"{basedir}/{file}"
    if not os.path.exists(f):
        os.makedirs(os.path.split(f)[0], exist_ok=True)
        if not download_file(f"{baseurl}{file}", f):
            logging.warning(
                f"Failed to download file {file} from {baseurl}, trying fallback URL {baseurl_fallback}"
            )
            baseurl = baseurl_fallback
            if not download_file(f"{baseurl_fallback}{file}", f):
                raise RuntimeError(f"Failed to download file {file}")
        logging.info(f"Downloaded {file}")
    else:
        logging.info(f"Skipping {file}, exists at {os.path.abspath(f)}")
