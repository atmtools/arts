import io
import logging
import os
import shutil
import sys
import time
import urllib.error
import urllib.request
import zipfile
from html.parser import HTMLParser

if os.environ.get("DEBUG_DOWNLOADING") is not None:
    import http.client

    http.client.HTTPConnection.debuglevel = 1

    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s %(levelname)s %(message)s"
    )
else:
    logging.basicConfig(level=logging.INFO, format="%(message)s")


def download_file(url, destination, nretry=1):
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
                time.sleep(backoff_factor * (2**attempt))
                continue
            logging.error(f"Download failed: {err}")
            return False
        except urllib.error.URLError as err:
            if attempt < nretry:
                time.sleep(backoff_factor * (2**attempt))
                continue
            logging.error(f"Download failed: {err}")
            return False
        except Exception as err:
            logging.error(f"Download failed: {err}")
            return False
    return False


def download_to_memory(url, nretry=1):
    """Download data from a URL into a BytesIO buffer."""
    backoff_factor = 1
    status_forcelist = [500, 502, 503, 504]

    for attempt in range(nretry + 1):
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req) as response:
                buf = io.BytesIO()
                while True:
                    chunk = response.read(8192)
                    if not chunk:
                        break
                    buf.write(chunk)
                buf.seek(0)
                return buf
        except urllib.error.HTTPError as err:
            if err.code in status_forcelist and attempt < nretry:
                time.sleep(backoff_factor * (2**attempt))
                continue
            logging.error(f"Download failed: {err}")
            return None
        except urllib.error.URLError as err:
            if attempt < nretry:
                time.sleep(backoff_factor * (2**attempt))
                continue
            logging.error(f"Download failed: {err}")
            return None
        except Exception as err:
            logging.error(f"Download failed: {err}")
            return None


class SVNIndexParser(HTMLParser):
    """Parse Apache SVN directory listing HTML for file links."""

    def __init__(self):
        super().__init__()
        self.entries = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            attrs_dict = dict(attrs)
            href = attrs_dict.get("href", "")
            if (
                href != "../"
                and "?" not in href
                and "#" not in href
                and "://" not in href
            ):
                self.entries.append(href)

    def error(self, message):
        pass


def list_svn_dir(url):
    """Return list of entry names from an SVN directory listing URL."""
    index_url = url if url.endswith("/") else (url + "/")
    try:
        req = urllib.request.Request(index_url)
        with urllib.request.urlopen(req) as response:
            html = response.read().decode("utf-8", errors="replace")
    except Exception as e:
        logging.warning(f"Failed to list SVN directory {index_url}: {e}")
        return None

    parser = SVNIndexParser()
    parser.feed(html)
    return parser.entries


def list_directory_svn(svn_base_url, path):
    """Recursively list files under a directory via SVN HTML fallback."""
    url = svn_base_url + path
    entries = list_svn_dir(url)
    if entries is None:
        return None

    files = []
    for entry in entries:
        if entry.endswith("/"):
            subfiles = list_directory_svn(svn_base_url, path + entry)
            if subfiles is None:
                return None
            files.extend(subfiles)
        else:
            files.append(path + entry)
    return files


def download_folder_zip(rel_path, basedir, baseurl):
    """Download a folder as a ZIP archive and extract in-place."""
    dest_dir = os.path.join(basedir, rel_path.rstrip("/"))
    if os.path.isdir(dest_dir):
        logging.info(
            f"Skipping {rel_path}, directory exists at {os.path.abspath(dest_dir)}"
        )
        return True

    raw_marker = "/-/raw/"
    idx = baseurl.find(raw_marker)
    if idx == -1:
        logging.error(f"Cannot derive archive URL from baseurl {baseurl}")
        return False

    project_base = baseurl[:idx]
    repo_name = project_base[project_base.rfind("/") + 1 :]
    archive_url = (
        f"{project_base}/-/archive/main/{repo_name}-main.zip"
        f"?ref_type=heads&path={rel_path.rstrip('/')}"
    )

    zip_buffer = download_to_memory(archive_url)
    if zip_buffer is None:
        return False

    try:
        with zipfile.ZipFile(zip_buffer) as z:
            # Detect the actual archive root prefix dynamically
            top_levels = {
                name.split("/")[0] + "/" for name in z.namelist() if "/" in name
            }
            if len(top_levels) != 1:
                raise RuntimeError(
                    f"Unexpected ZIP structure, top levels: {top_levels}"
                )
            archive_root = next(iter(top_levels))
            target_prefix = archive_root + rel_path.rstrip("/") + "/"

            for member in z.namelist():
                if not member.startswith(target_prefix):
                    continue
                if member.endswith("/"):
                    continue
                inner = member[len(target_prefix) :]
                if not inner:
                    continue
                final_path = os.path.join(dest_dir, inner)
                os.makedirs(os.path.dirname(final_path), exist_ok=True)
                with z.open(member) as src, open(final_path, "wb") as dst:
                    shutil.copyfileobj(src, dst)
        logging.info(f"Downloaded and extracted {rel_path}")
        return True
    except Exception as e:
        logging.error(f"Failed to extract ZIP for {rel_path}: {e}")
        return False


def download_item(rel_path, basedir, primary_url, fallback_url):
    """Download a single file if it doesn't already exist."""
    global use_fallback
    f = f"{basedir}/{rel_path}"
    if os.path.exists(f):
        logging.info(f"Skipping {rel_path}, exists at {os.path.abspath(f)}")
        return True

    os.makedirs(os.path.split(f)[0], exist_ok=True)

    if not use_fallback:
        if download_file(f"{primary_url}{rel_path}", f):
            logging.info(f"Downloaded {rel_path}")
            return True

        logging.warning(
            f"Failed to download file {rel_path} from {primary_url}, "
            f"trying fallback URL {fallback_url}"
        )
        use_fallback = True

    if download_file(f"{fallback_url}{rel_path}", f):
        logging.info(f"Downloaded {rel_path} (fallback)")
        return True

    return False


# If any download fails, use the fallback URL for all remaining downloads
use_fallback = False

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
    items = [
        line.strip() for line in f if line.strip() and not line.strip().startswith("#")
    ]

for item in items:
    logging.info(f"Processing {item}...")

    if item.endswith("/"):
        if not use_fallback:
            if download_folder_zip(item, basedir, baseurl):
                continue
            logging.warning(
                f"Failed to download folder {item} from GitLab as ZIP, "
                f"trying fallback URL {baseurl_fallback}"
            )
            use_fallback = True

        files_to_download = list_directory_svn(baseurl_fallback, item)
        if files_to_download is None:
            raise RuntimeError(f"Failed to list directory {item}")
        for rel_path in files_to_download:
            if not download_item(rel_path, basedir, baseurl, baseurl_fallback):
                raise RuntimeError(f"Failed to download file {rel_path}")
    else:
        if not download_item(item, basedir, baseurl, baseurl_fallback):
            raise RuntimeError(f"Failed to download file {item}")
