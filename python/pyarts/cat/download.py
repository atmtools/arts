import os
import urllib.request
import zipfile


def _download_and_extract(url, extract_dir=".", verbose=False):
    """
    Downloads a zip file from a given URL and extracts it to a specified directory.

    Parameters:
        url (str): The URL of the zip file to download.
        extract_dir (str, optional): The directory to extract the zip file to. Defaults to ".".
        verbose (bool, optional): Whether to print info messages. Defaults to False.
    """
    if verbose:
        print(f"Downloading {url}")
    try:
        zip_path, _ = urllib.request.urlretrieve(url)
    except urllib.request.HTTPError as e:
        raise RuntimeError(f"Failed to download {url}: {e}")
    with zipfile.ZipFile(zip_path, "r") as f:
        f.extractall(extract_dir)


def _get_latest_stable_version(version):
    """
    Decrement version to previous stable version.

    Parameters:
        version (str): Version string in format "major.minor.micro".
    """
    major, minor, micro = map(int, version.split("."))
    return f"{major}.{minor}.{micro - 1 if micro % 2 else micro}"


def retrieve(download_dir=None, version=None, verbose=False):
    """
    Download and extract the ARTS XML and catalog data files from github.

    This function sets the ARTS data search path to the downloaded
    data directories, so that ARTS can find the required data files.

    If the environment variable ARTS_DATA_PATH or ARTS_INCLUDE_PATH is set,
    it is assumed that the user wants to use their own catalog locations
    and this function does nothing.

    Parameters:
        download_dir (str, optional):
            The directory where the data files will be stored.
            If not provided, the default is `~/.cache/arts`.
        version (str, optional):
            The version of ARTS to download the data files for. The default
            is the version of the currently installed pyarts package.
        verbose (bool, optional):
            Whether to print info messages. Defaults to False.
    """
    if version is None:
        from pyarts import __version__

        version = _get_latest_stable_version(__version__)
        if version != __version__:
            if verbose:
                print(
                    f"No downloadable catalogs are available for development versions of ARTS.\n"
                    f"Downloading the latest stable version {version} of the catalogs instead.\n"
                    f"Check out the catalogs with svn if you need the latest developmen version."
                )

    if download_dir is None:
        download_dir = os.path.join(os.getenv("HOME"), ".cache", "arts")

    GITHUB_URL = f"https://github.com/atmtools/arts/releases/download/v{version}/"
    artscatdata = "arts-cat-data-" + version
    artsxmldata = "arts-xml-data-" + version
    if os.getenv("ARTS_DATA_PATH"):
        if verbose:
            print("Skipping download, environment variable ARTS_DATA_PATH already set.")
        return
    if os.getenv("ARTS_INCLUDE_PATH"):
        if verbose:
            print("Skipping download, environment variable ARTS_INCLUDE_PATH already set.")
        return

    def retrieve_catalog(catname):
        if not os.path.exists(os.path.join(download_dir, catname)):
            os.makedirs(download_dir, exist_ok=True)
            _download_and_extract(GITHUB_URL + catname + ".zip",
                                  download_dir,
                                  verbose=verbose)
        elif verbose:
            print(f"Skipping download, data already exists at {download_dir}/{catname}")

    retrieve_catalog(artsxmldata)
    retrieve_catalog(artscatdata)

    from pyarts.arts import globals
    globals.parameters.datapath = [
        ".",
        os.path.join(download_dir, artsxmldata),
        os.path.join(download_dir, artscatdata),
    ]
