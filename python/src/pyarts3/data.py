# -*- coding: utf-8 -*-
"""Common functions for the pyarts package."""

import os
import urllib.request
import zipfile
from pyarts3.arts.globals import parameters
import pyarts3 as pyarts
import numpy as np
from tqdm import tqdm
import xarray


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
        from pyarts3 import __version__

        version = __version__

    USE_SVN = False
    GITHUB_URL = f"https://github.com/atmtools/arts/releases/download/v{version}/"
    if int(version[-1]) % 2:
        if "CI" not in os.environ:
            USE_SVN = True
        else:
            raise RuntimeError(
                f"Version {version} is not a release version.\n"
                f"CI environment detected, SVN checkout prevented."
            )

    fullcatname = catname + "-" + version if not USE_SVN else catname + "-trunk"
    catdir = os.path.join(download_dir, fullcatname)
    if not os.path.exists(catdir):
        os.makedirs(download_dir, exist_ok=True)
        if USE_SVN:
            import subprocess

            print(f"Checking out {catname} from SVN")
            svn_path = f"https://radiativetransfer.org/svn/rt/{catname}/trunk/"
            command = ["svn", "checkout", svn_path, catdir]
            try:
                subprocess.run(command)
            except FileNotFoundError as e:
                raise FileNotFoundError(
                    "svn command not found, cannot check out from SVN"
                ) from e
        else:
            _download_and_extract(
                GITHUB_URL + fullcatname + ".zip", download_dir, verbose=verbose
            )

    if USE_SVN:
        print(
            f"Using {catname} from SVN, update manually if needed: `svn update {catdir}`"
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


def to_atmospheric_field(
    data: xarray.Dataset,
    remap: None | dict[str, str] = None,
    ignore: None | list[str] = None,
    *,
    atm: None | pyarts.arts.AtmField = None,
) -> pyarts.arts.AtmField:
    """
    Populates a ~pyarts3.arts.AtmField from an xarray Dataset-like structure

    Parameters
    ----------
    data : xarray.Dataset
        A dataset.  Coordinates should contain 'alt', 'lat', and 'lon'.
        All other keys() should be assignable to the atm-field by name.
    remap : None | dict[str, str], optional
        All names in this optional dict are renamed when accessing the dataset.
        For example, if the altitude-grid is called 'Alt' instead of 'alt', the
        dict should contain {..., 'alt': 'Alt', ...}.
    ignore : None | list[str], optional
        Ignore keys listed from assignment into the atmospheric field.
    atm : None | pyarts3.arts.AtmField
        The default atmospheric field to use.  Defaults to None to use default-
        constructed atmospheric field object.

    Returns
    -------
    atm : pyarts3.arts.AtmField
        An atmospheric field

    """
    if ignore is None:
        ignore = []

    if remap is None:
        remap = {}

    alt = (
        getattr(data, remap.get("alt", "alt")).data.flatten()
        if "alt" not in ignore
        else np.array([0])
    )
    lat = (
        getattr(data, remap.get("lat", "lat")).data.flatten()
        if "lat" not in ignore
        else np.array([0])
    )
    lon = (
        getattr(data, remap.get("lon", "lon")).data.flatten()
        if "lon" not in ignore
        else np.array([0])
    )

    GF3 = pyarts.arts.GeodeticField3(
        name="Generic",
        data=np.zeros((alt.size, lat.size, lon.size)),
        grid_names=["Altitude", "Latitude", "Longitude"],
        grids=[alt, lat, lon],
    )

    if atm is None:
        atm = pyarts.arts.AtmField()

    if hasattr(data, "top_of_atmosphere"):
        atm.top_of_atmosphere = data.top_of_atmosphere
    else:
        atm.top_of_atmosphere = max(alt)

    for k in data.keys():
        try:
            if k in ignore:
                continue

            k = remap.get(k, k)
            kstr = str(k)

            GF3.dataname = kstr
            np.asarray(GF3.data).flat[:] = getattr(data, kstr).data.flat[:]

            atm[k] = GF3

        except Exception as e:
            raise Exception(f"{e}\n\nFailed to assign key to atm: {k}")

    return atm


def to_absorption_species(
    atm_field: pyarts.arts.AtmField,
) -> pyarts.arts.ArrayOfArrayOfSpeciesTag:
    """Scans the ARTS data path for species relevant to the given atmospheric field.

    The scan is done over files in an arts-cat-data like directory structure.

    Args:
        atm_field (pyarts3.arts.AtmField): A relevant atmospheric field.

    Returns:
        pyarts3.arts.ArrayOfArrayOfSpeciesTag: All found species tags.
        The intent is that this is enough information to use pyarts3.Workspace.ReadCatalogData
    """
    species = atm_field.species_keys()

    out = []
    for spec in species:
        out.append(f"{spec}")

        if pyarts.arts.file.find_xml(f"xsec/{spec}-XFIT") is not None:
            out.append(f"{spec}-XFIT")

        for spec2 in species:
            if pyarts.arts.file.find_xml(f"cia/{spec}-CIA-{spec2}"):
                out.append(f"{spec}-CIA-{spec2}")

            if pyarts.arts.file.find_xml(f"cia/{spec2}-CIA-{spec}"):
                out.append(f"{spec2}-CIA-{spec}")

        if pyarts.arts.file.find_xml(f"cia/{spec}-CIA-{spec}") is not None:
            out.append(f"{spec}-CIA-{spec}")

        if spec == pyarts.arts.SpeciesEnum.Water:
            out.append("H2O-ForeignContCKDMT400")
            out.append("H2O-SelfContCKDMT400")
        elif spec == pyarts.arts.SpeciesEnum.CarbonDioxide:
            out.append("CO2-CKDMT252")

    return pyarts.arts.ArrayOfArrayOfSpeciesTag(np.unique(out))


def xarray_open_dataset(filename_or_obj, *args, **kwargs):
    """Wraps xarray.open_dataset to search for files in the current and in ARTS' data path.

    All args and kwargs are passed on to xarray.open_dataset directly.  Any FileNotFoundError
    will be caught and just raised if no path works.

    Args:
        filename_or_obj (_type_): _description_

    Raises:
        FileNotFoundError: _description_

    Returns:
        _type_: _description_
    """

    try:
        return xarray.open_dataset(filename_or_obj, *args, **kwargs)
    except FileNotFoundError:
        pass

    for p in parameters.datapath:
        try:
            return xarray.open_dataset(
                os.path.join(p, filename_or_obj), *args, **kwargs
            )
        except FileNotFoundError:
            pass

    raise FileNotFoundError(
        f"File not found in ARTS data path ({parameters.datapath}): {filename_or_obj}"
    )
