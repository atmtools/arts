from pathlib import Path
from ftplib import FTP
import os
import sys
import urllib.request

nargs = len(sys.argv)

XML = "https://arts.mi.uni-hamburg.de/svn/rt/arts-xml-data/trunk/"
CAT = "https://arts.mi.uni-hamburg.de/svn/rt/arts-cat-data/trunk/"
SSDB = (
    "ftp-projects.cen.uni-hamburg.de",
    "arts/ArtsScatDbase/v1.1.0/SSD/TotallyRandom"
)


def http_download(url, local_dir, path):
    """
    Download files from a http server.

    Args:
        url: String specifying the URL of the HTTP server.
        local_dir: Path object specifying the local directory to store the downloaded files in.
        path: A path specifying the location of the file to download on the
            server. This will be interpreted as a path relative to 'local_dir'
            to determine the output file to which the downloaded data will
            be written.
    """
    dest_path = local_dir / path
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {path}")
    urllib.request.urlretrieve(f"{url}{path}", dest_path)


def ftp_download(url, local_dir, path):
    """
    Download files from FTP server.

    Args:
        url: String specifying the URL of the FTP server or a tuple
             ``(url, base_path)`` specifying a base directory on the FTP
             server. In the latter case, the 'path' argument will be
             interpreted as a relative path with respect to ``base_path``.
        local_dir: Path object specifying the local directory to store the
            downloaded files in.
        path: A path specifying the location of the file to download on the
            server. This will be interpreted as a path relative to 'local_dir'
            to determine the output file to which the downloaded data will
            be written.
    """
    path = Path(path)
    dest_path = local_dir / path
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(url, tuple):
        url, base_dir = url
    else:
        base_dir = Path("")
    ftp = FTP(url)
    ftp.login()
    ftp.cwd(str((base_dir / path).parent))
    print(f"Downloading {path}")
    with open(dest_path, "wb") as output:
        ftp.retrbinary(f"RETR {path.name}", output.write)


if nargs > 3:
    workdir = Path.cwd()
    try:
        os.chdir(sys.argv[1])
        if sys.argv[2] == "xml":
            url = XML
            dest = Path("arts-xml-data")
            download = http_download
        elif sys.argv[2] == "cat":
            url = CAT
            dest = Path("arts-cat-data")
            download = http_download
        elif sys.argv[2] == "ssdb":
            url = SSDB
            dest = Path("arts-ssdb-data")
            download = ftp_download
        else:
            url = ""
            dest = Path("custom-data")
            download = http_download

        for path in sys.argv[3:]:
            output_file = dest / path
            if not output_file.exists():
                download(url, dest, path)
            else:
                print(f"Skipping {path}, exists at {output_file.absolute()}")
    finally:
        os.chdir(workdir)

