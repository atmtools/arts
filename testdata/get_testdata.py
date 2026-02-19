import logging
import os
import sys
import requests
import requests.adapters
from requests.adapters import HTTPAdapter, Retry

if os.environ.get("DEBUG_DOWNLOADING") is not None:
    import http.client

    http.client.HTTPConnection.debuglevel = 1

    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s %(levelname)s %(message)s"
    )
    logging.getLogger("urllib3").setLevel(logging.DEBUG)
    logging.getLogger("urllib3").propagate = True
else:
    logging.basicConfig(level=logging.INFO, format="%(message)s")


def download_file(url, destination, nretry=5):
    """Download a file from a URL to a destination file."""
    try:
        # Send a GET request to the URL
        s = requests.Session()
        retries = Retry(
            total=nretry, backoff_factor=1, status_forcelist=[500, 502, 503, 504]
        )
        s.mount("https://", HTTPAdapter(max_retries=retries))
        response = s.get(url)

        # Check if the request was successful
        response.raise_for_status()  # Raise an error for bad responses (4xx and 5xx)

        # Open the destination file in write binary mode
        with open(destination, "wb") as file:
            # Write the content to the file in chunks
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)

        return True
    except requests.exceptions.RequestException as err:
        logging.error(f"Download failed: {err}")
        return False


nargs = len(sys.argv)

xml = "https://gitlab.rrz.uni-hamburg.de/atmtools/arts-xml-data/-/raw/main/"
cat = "https://gitlab.rrz.uni-hamburg.de/atmtools/arts-cat-data/-/raw/main/"
xml_fallback = "https://arts.mi.uni-hamburg.de/svn/rt/arts-xml-data/trunk/"
cat_fallback = "https://arts.mi.uni-hamburg.de/svn/rt/arts-cat-data/trunk/"

if sys.argv[1] == "xml":
    baseurl = xml
    baseurl_fallback = xml_fallback
    basedir = "arts-xml-data"
elif sys.argv[1] == "cat":
    baseurl = cat
    baseurl_fallback = cat_fallback
    basedir = "arts-cat-data"
else:
    baseurl = ""
    basedir = "custom-data"

for file in sys.argv[2:]:
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
