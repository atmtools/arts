# -*- coding: utf-8 -*-
"""Common functions for the pyarts package."""

import os
import urllib.request
import zipfile
from pyarts.arts.globals import parameters

from tqdm import tqdm


def download(data=("xml", "cat"), download_dir=None, verbose=False, **kwargs):
    """
    Download and extract data files.

    This function sets the ARTS data search path to the downloaded
    data directories, so that ARTS can find the required data files.

    If the environment variable ARTS_DATA_PATH or ARTS_INCLUDE_PATH is set,
    it is assumed that the user wants to use their own catalog locations
    and this function does nothing.

    Parameters:
        data (Tuple[str]):
            List of data types to download. Possible values are:
            - ``"xml"``: arts-xml-data package
            - ``"cat"``: arts-cat-data package
        download_dir (str, optional):
            The directory where the data files will be stored.
            If not provided, the default is ``~/.cache/arts/``.
        verbose (bool, optional):
            Whether to print info messages. Defaults to False.
        **kwargs:
            Additional options passed on to the data specific download functions, e.g.
            ``version`` for ARTS catalogs.
    """
    if os.getenv("ARTS_DATA_PATH"):
        if verbose:
            print("Skipping download, environment variable ARTS_DATA_PATH already set.")
        return
    if os.getenv("ARTS_INCLUDE_PATH"):
        if verbose:
            print(
                "Skipping download, environment variable ARTS_INCLUDE_PATH already set."
            )
        return

    if download_dir is None:
        download_dir = os.path.join(os.path.expanduser("~"), ".cache", "arts")

    datadirs = []
    for d in data:
        match d:
            case "xml":
                datadirs.append(
                    download_arts_xml_data(download_dir, verbose=verbose, **kwargs)
                )
            case "cat":
                datadirs.append(
                    download_arts_cat_data(download_dir, verbose=verbose, **kwargs)
                )
            case _:
                raise RuntimeError(f'Unknown download data type "{d}"')

    return datadirs


class _DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        return self.update(b * bsize - self.n)


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
        with _DownloadProgressBar(
            unit="B", unit_scale=True, miniters=1, desc=url.split("/")[-1]
        ) as t:
            zip_path, _ = urllib.request.urlretrieve(url, reporthook=t.update_to)
    except urllib.request.HTTPError as e:
        raise RuntimeError(f"Failed to download {url}: {e}")
    with zipfile.ZipFile(zip_path, "r") as f:
        f.extractall(extract_dir)


def _download_arts_data(catname, download_dir=None, version=None, verbose=False):
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
            ggkIf not provided, the default is `~/.cache/arts`.
        version (str, optional):
            The version of ARTS to download the data files for. The default
            is the version of the currently installed pyarts package.
        verbose (bool, optional):
            Whether to print info messages. Defaults to False.
    """
    if os.getenv("ARTS_DATA_PATH"):
        if verbose:
            print("Skipping download, environment variable ARTS_DATA_PATH already set.")
        return
    if os.getenv("ARTS_INCLUDE_PATH"):
        if verbose:
            print(
                "Skipping download, environment variable ARTS_INCLUDE_PATH already set."
            )
        return

    if download_dir is None:
        download_dir = os.path.join(os.path.expanduser("~"), ".cache", "arts")

    if version is None:
        from pyarts import __version__

        version = __version__

    GITHUB_URL = f"https://github.com/atmtools/arts/releases/download/v{version}/"
    if int(version[-1]) % 2:
        raise RuntimeError(
            f"Version {version} is not a release version.\n"
            f"Please check out the current catalogs with svn instead."
        )

    catname = catname + "-" + version
    catdir = os.path.join(download_dir, catname)
    if not os.path.exists(catdir):
        os.makedirs(download_dir, exist_ok=True)
        _download_and_extract(
            GITHUB_URL + catname + ".zip", download_dir, verbose=verbose
        )

    parameters.datapath.append(catdir)

    return catdir


def download_arts_xml_data(download_dir=None, version=None, verbose=False):
    """
    Download and extract the ARTS XML data files from github.

    This function sets the ARTS data search path to the downloaded
    data directories, so that ARTS can find the required data files.

    If the environment variable ARTS_DATA_PATH or ARTS_INCLUDE_PATH is set,
    it is assumed that the user wants to use their own catalog locations
    and this function does nothing.

    Parameters:
        download_dir (str, optional):
            The directory where the data files will be stored.
            ggkIf not provided, the default is `~/.cache/arts`.
        version (str, optional):
            The version of ARTS to download the data files for. The default
            is the version of the currently installed pyarts package.
        verbose (bool, optional):
            Whether to print info messages. Defaults to False.
    """
    return _download_arts_data("arts-xml-data", download_dir, version, verbose)


def download_arts_cat_data(download_dir=None, version=None, verbose=False):
    """
    Download and extract the ARTS catalog data files from github.

    This function sets the ARTS data search path to the downloaded
    data directories, so that ARTS can find the required data files.

    If the environment variable ARTS_DATA_PATH or ARTS_INCLUDE_PATH is set,
    it is assumed that the user wants to use their own catalog locations
    and this function does nothing.

    Parameters:
        download_dir (str, optional):
            The directory where the data files will be stored.
            ggkIf not provided, the default is `~/.cache/arts`.
        version (str, optional):
            The version of ARTS to download the data files for. The default
            is the version of the currently installed pyarts package.
        verbose (bool, optional):
            Whether to print info messages. Defaults to False.
    """
    return _download_arts_data("arts-cat-data", download_dir, version, verbose)
