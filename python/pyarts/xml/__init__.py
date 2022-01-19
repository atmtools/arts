# -*- coding: utf-8 -*-

"""This module provides functionality for reading and writing ARTS XML files.
"""

import gzip
import glob
import itertools
import os
from os.path import isfile, join, basename, splitext, dirname
import tempfile

from . import read
from . import write

__all__ = [
    'load',
    'save',
    'load_directory',
    'load_indexed',
    'make_binary',
    'make_directory_binary',
    'update',
    'update_directory'
]


def save(var, filename, precision='.7e', format='ascii', comment=None,
         parents=False):
    """Save a variable to an ARTS XML file.

    Args:
        var: Variable to be stored.
        filename (str): Name of output XML file.
            If the name ends in .gz, the file is compressed on the fly.
        precision (str): Format for output precision.
        format (str): Output format: 'ascii' (default) or 'binary'.
        comment (str): Comment string included in a tag above data.
        parents (bool): Create missing parent directories.

    Note:
        Python's gzip module is extremely slow in writing. Consider
        compressing files manually after writing them normally.

    Example:
        >>> x = numpy.array([1.,2.,3.])
        >>> pyarts.xml.save(x, 'myvector.xml')

    """
    if parents:
        os.makedirs(dirname(filename), exist_ok=True)

    if filename.endswith('.gz'):
        if format != 'ascii':
            raise RuntimeError(
                'For zipped files, the output format must be "ascii"')
        xmlopen = gzip.open
    else:
        xmlopen = open

    if hasattr(var, "savexml"):
        if format == "ascii" and filename.endswith('.gz'):
            format = "zascii"
        var.savexml(filename, format)
        return

    with xmlopen(filename, mode='wt', encoding='UTF-8') as fp:
        if format == 'binary':
            with open(filename + '.bin', mode='wb') as binaryfp:
                axw = write.ARTSXMLWriter(fp, precision=precision,
                                          binaryfp=binaryfp)
                axw.write_header()
                if comment is not None:
                    axw.write_comment(comment)
                axw.write_xml(var)
                axw.write_footer()
        elif format == 'ascii':
            axw = write.ARTSXMLWriter(fp, precision=precision)
            axw.write_header()
            if comment is not None:
                axw.write_comment(comment)
            axw.write_xml(var)
            axw.write_footer()
        else:
            raise RuntimeError('Unknown output format "{}".'.format(format))


def load(filename, search_arts_path=True):
    """Load a variable from an ARTS XML file.

    The input file can be either a plain or gzipped XML file.

    By default, the current directory, ARTS_INCLUDE_PATH and ARTS_DATA_PATH
    environment variables are searched (in that order) if the passed filename is
    a relative path.

    Args:
        filename (str): Name of ARTS XML file.
        search_arts_path (bool): Set to False to ignore ARTS search paths.

    Returns:
        Data from the XML file. Type depends on data in file.

    Example:
        >>> pyarts.xml.load('tests/reference/matrix.xml')
        array([[ 0.,  1.],
               [ 2.,  3.]])

    """

    if search_arts_path and not os.path.isabs(filename):
        # Use dict to remove duplicate paths from list and preserve order
        for path in dict.fromkeys(
            [""] + os.environ.get("ARTS_INCLUDE_PATH", "").split(":") +
                os.environ.get("ARTS_DATA_PATH", "").split(":")).keys():
            checkfile = os.path.join(path, filename)
            if isfile(checkfile) or isfile(checkfile + ".gz"):
                filename = checkfile
                break

    # If file is not found, try the gzipped version.
    if not isfile(filename):
        if not isfile(filename + '.gz'):
            raise FileNotFoundError("No such file: '{}'".format(filename))
        else:
            filename += '.gz'

    if filename.endswith('.gz'):
        xmlopen = gzip.open
    else:
        xmlopen = open

    binaryfilename = filename + '.bin'
    artstag = None
    try:
        with xmlopen(filename, 'rb') as fp:
            if isfile(binaryfilename):
                with open(binaryfilename, 'rb',) as binaryfp:
                    root = read.parse(fp, binaryfp).getroot()
                    artstag = root[0].tag
                    return root.value()
            else:
                root = read.parse(fp).getroot()
                artstag = root[0].tag
                return root.value()
    except RuntimeError:
        if artstag == "Array" and len(root) and 'type' in root[0].attrib:
            artstag = f"ArrayOf{root[0].attrib['type']}"  # Fix for arrays
        
        import pyarts.classes as native_classes
        if hasattr(native_classes, artstag) and hasattr(
                getattr(native_classes, artstag), "readxml"):
            ret = getattr(native_classes, artstag)()
            ret.readxml(filename)
            return ret
        raise


def load_directory(directory, exclude=None):
    """Load all XML files in a given directory.

    Search given directory  for files with ``.xml`` or ``.xml.gz`` extension
    and try to load them using :func:`load`.

    Parameters:
        directory (str): Path to the directory.
        exclude (Container[str]): Filenames to exclude.

    Returns:
        dict: Filenames without extension are keys for the file content.

    Example:
        Load all files in ``foo`` except for the lookup table in
        ``abs_lookup.xml.``

        >>> load_directory('foo', exclude=['abs_lookup.xml'])
    """
    def includefile(f):
        """Check if to include file."""
        return basename(f) not in exclude if exclude is not None else True

    def stripext(f):
        """Strip the extension of a filename."""
        return splitext(f)[0]

    # Create a generator yielding all XML files to load (not excluded).
    xmlfiles = filter(includefile, glob.iglob(join(directory, '*.xml')))

    # Remove extension from zipped files to keep dictionary keys clean.
    # The `load` function looks for zipped files anyway.
    gzfiles = filter(includefile, glob.iglob(join(directory, '*.xml.gz')))
    gzfiles = map(stripext, gzfiles)

    # Store XML file contents in a dictionary, using the filename as key.
    return {stripext(basename(f)): load(f)
            for f in itertools.chain(xmlfiles, gzfiles)}


def load_indexed(filename):
    """Load all indexed XML files matching the given filename.

    The function searches all files matching the pattern
    ``<filename>.<file_index>.xml`` or ``<filename>.<file_index>.xml.gz``.

    A list with the loaded file contents is returned. The list indices are
    equivalent to the file indices.

    Parameters:
        filename (str): Filename.

    Returns:
        list: List of file contents.

    Example:
        Load all files matching the pattern ``foo.<file_index>.xml``.

        >>> load_indexed_xml('foo')

    """
    iidx = -2  # Relative position of fileindex in splitted filename.

    # Get all files matching the indexed filename format.
    files = glob.glob('{}.*.xml'.format(filename))

    # If no files are found, try the gzipped version.
    if len(files) == 0:
        files = glob.glob('{}.*.xml.gz'.format(filename))
        iidx = -3  # Correct fileindex position for gzipped files.

    # Extract indices from filenames.
    maxindex = max(int(x.split('.')[iidx]) for x in files)

    # Pre-allocate a list according to the maximum index found.
    ret = (maxindex + 1) * [None]

    # Fill list with file contents (file index matching list index).
    for f in files:
        findex = int(f.split('.')[iidx])
        ret[findex] = load(f)

    return ret


def make_binary(filename, out='', absolute_out=False, parents=True):
    """Loads xml-file at filename and saves it back in binary format

    Parameters:
        filename (str): Filename path.
        out (str): Path to save the binary.  Empty causes overwrite of file.
        absolute_out (bool): If true, then write file to out-path rather than
            to the relative path out.  Does nothing if file is in the working
            folder and out is relative.
        parents (bool): Create missing parent directories.

    Returns:
        str: Path to the created binary file.

    Example:
        Load t_field.xml and save it back as binary it as ./binary/t_field.xml
        and ./binary/t_field.bin

        >>> make_binary('t_field.xml', out='binary')
        'binary/t_field.xml'
    """

    xml_data = load(filename)
    if absolute_out:
        outfile = join(out, basename(filename))
    else:
        outfile = join(dirname(filename), out, basename(filename))

    save(xml_data, outfile, format='binary', parents=parents)

    return outfile


def make_directory_binary(directory, out='', absolute_out=False, parents=True):
    """Loads xml-files in directory and saves them back in binary format

    Parameters:
        directory (str): Directory path.
        out (str): Path to save the binary.
        absolute_out (bool): If true, then write file to out-path rather than
            to the relative path out.  Does nothing if file is in the working
            folder and out is relative.
        parents (bool): Create missing parent directories.

    Returns:
        list[str]: Paths to the created binary files.

    Example:
        Load arts-xml-data/spectroscopy/cia/hitran2011/ and save it back as
        binary it at arts-xml-data-binary/spectroscopy/cia/hitran2011/

        >>> make_directory_binary('arts-xml-data/spectroscopy/cia/hitran2011',
            out='arts-xml-data-binary/spectroscopy/cia/hitran2011',
            absolute_out=True)
        ['arts-xml-data-binary/spectroscopy/cia/hitran2011/hitran_cia2012_adapted.xml']
    """

    directory_of_xmls = load_directory(directory)
    outfiles = []  # Empty list to store output filepaths.

    if absolute_out:
        get_outfile = join(out, '{entry}.xml')
    else:
        get_outfile = join(directory, out, '{entry}.xml')

    for entry in directory_of_xmls:
        outfile = get_outfile.format(entry=entry)
        save(directory_of_xmls[entry],
             outfile,
             format='binary',
             parents=parents)
        outfiles.append(outfile)

    return outfiles


def update(filename, precision='%g'):
    """ Updates a file to the latest version of ARTS
    
    Wraps load()+save() from higher up in this file
    
    Only works for absolute paths.  To ensure the path is absolute,
    os.path.abspath is applied on filename as the first operation
    
    Attempts to store the file in the same format as it was read by, clobbering
    the original file.
    
    Parameters:
        filename (str): Filename path.
        precision (str): Format for output precision.
    """
    filename = os.path.abspath(filename)
    
    # Test file format so that the '.gz' is caught primarily
    if not os.path.isfile(filename) and os.path.isfile(filename + '.gz'):
        filename += '.gz'
    
    # Find the format for saving, it must be understood
    tascii = filename.endswith('.xml')
    zascii = filename.endswith('.gz')
    binary = os.path.isfile(filename + '.bin')
    format = 'binary' if binary else 'ascii'  # load will deal with zascii
    if not (tascii or zascii):
        raise RuntimeError(f'Must end with .xml or .gz, reads: {filename}')
    
    n = next(tempfile._get_candidate_names())
    fn = filename + f".pyarts.tmpfile.{n}.{'xml.gz' if zascii else 'xml'}"
    while os.path.isfile(fn) or os.path.isfile(fn + '.bin'):
        fn += f".pyarts.tmpfile.{n}.{'xml.gz' if zascii else 'xml'}"
    
    # load+save the file to temporary names
    save(load(filename, False), fn, precision=precision, format=format)
    
    # Move the files back
    os.rename(fn, filename)
    if binary: os.rename(fn + '.bin', filename + '.bin')

def update_directory(directory, precision='%g'):
    """ Update all files in a directory
    
    Wraps update() from higher up in this file for all files ending with
    .xml in the given directory
    
    Only works for absolute paths.  To ensure the path is absolute,
    os.path.abspath is applied on directory as the first operation
    
    There is a subset of .xml files that cannot be read by standard ARTS but
    requires specialized functions from within ARTS to be read.  For example,
    the old Artscat-N format line catalog files.
    
    Parameters:
        directory (str): Directory path.
        precision (str): Format for output precision.
        
    Returns:
        A dict of files with failed conversions and their error representations
    """
    directory = os.path.abspath(directory)
    
    out = {}
    for file in os.listdir(directory):
        if file.endswith('.xml'):
            fn = os.path.join(directory, file)
            try:
                update(fn, precision=precision)
            except Exception as e:
                out[fn] = repr(e)
    return out
