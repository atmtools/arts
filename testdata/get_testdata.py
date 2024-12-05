import os
import sys
import requests
import requests.adapters
from requests.adapters import HTTPAdapter, Retry
from time import sleep


def download_file(url, destination):
    """Download a file from a URL to a destination file."""
    try:
        # Send a GET request to the URL
        s = requests.Session()
        retries = Retry(
            total=5, backoff_factor=0.1, status_forcelist=[500, 502, 503, 504]
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
    except requests.exceptions.RequestException:
        return False


def download_file_retry(url, destination, retries=30, delay=10):
    """Download a file from a URL to a destination file with retries."""
    x = 0
    for i in range(retries):
        sleep(x)
        x = delay
        if download_file(url, destination):
            print(f"Downloaded {destination}")
            return
        print(f"Failed to download from {url}, sleeping for {delay} seconds")
    raise RuntimeError(f"Failed to download file {destination} after {retries} retries")


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
        download_file_retry(f"{baseurl}{file}", f)
    else:
        print(f"Skipping {file}, exists at {os.path.abspath(f)}")
