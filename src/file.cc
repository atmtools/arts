/* Copyright (C) 2000-2012 Stefan Buehler <sbuehler@ltu.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
   \file  file.cc

   This file contains basic functions to handle ASCII data files.

   \author Patrick Eriksson
   \date 2000-10-28 
*/

////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "arts.h"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>

#include "array.h"
#include "file.h"
#include "parameters.h"

////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

/**
   Gives the default file name for the ASCII formats.

   The default name is only used if the file name is empty.

   \param   filename Output:     file name
   \param    varname      variable name

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void filename_ascii(String& filename, const String& varname) {
  if ("" == filename) {
    filename = varname + ".aa";
  }
}

////////////////////////////////////////////////////////////////////////////
//   Functions to open and read ASCII files
////////////////////////////////////////////////////////////////////////////

/**
   Open a file for writing. If the file cannot be opened, the
   exception IOError is thrown. 
   @param     file File pointer 
   @param     name Name of the file to open
   @author    Stefan Buehler
   @version   1
   @exception ios_base::failure Could for example mean that the
                      directory is read only. */
void open_output_file(ofstream& file, const std::string_view name) {
  String ename = add_basedir(name);

  try {
    // Tell the stream that it should throw exceptions.
    // Badbit means that the entire stream is corrupted, failbit means
    // that the last operation has failed, but the stream is still
    // valid. We don't want either to happen!
    // FIXME: This does not yet work in  egcs-2.91.66, try again later.
    file.exceptions(ios::badbit | ios::failbit);

    file.open(ename.c_str());

    // See if the file is ok.
    // FIXME: This should not be necessary anymore in the future, when
    // g++ stream exceptions work properly. (In that case we would not
    // get here if there really was a problem, because of the exception
    // thrown by open().)
  } catch (const std::exception& e) {
    ARTS_USER_ERROR(
        "Cannot open output file: ",
        ename,
        "\nMaybe you don't have write access to the directory or the file?");
  }
}

/**
 Closes the file. If it is empty, the file is deleted.
 @param     file File pointer 
 @author    Oliver Lemke
 @exception ios_base::failure Could for example mean that the
 directory is read only. */
void cleanup_output_file(ofstream& file, const std::string_view name) {
  if (file.is_open()) {
    streampos fpos = file.tellp();
    file.close();
    if (!fpos) std::filesystem::remove(expand_path(name).c_str());
  }
}

/**
   Open a file for reading. If the file cannot be opened, the
   exception IOError is thrown. 
   @param     file File pointer 
   @param     name Name of the file to open
   @author    Stefan Buehler
   @version   1
   @exception ios_base::failure Somehow the file cannot be opened. */
void open_input_file(ifstream& file, const std::string_view name) {
  String ename{expand_path(name)};

  // Command line parameters which give us the include search path.
  extern const Parameters parameters;
  ArrayOfString allpaths = parameters.includepath;
  allpaths.insert(
      allpaths.end(), parameters.datapath.begin(), parameters.datapath.end());

  ArrayOfString matching_files;
  find_file(matching_files, ename, allpaths);

  if (matching_files.nelem()) ename = matching_files[0];

  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted.
  // On the other hand, end of file will not lead to an exception, you
  // have to check this manually!
  file.exceptions(ios::badbit);

  // c_str explicitly converts to c String.
  file.open(ename.c_str());

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly.
  ARTS_USER_ERROR_IF(!file,
                     "Cannot open input file: ",
                     ename,
                     "\nMaybe the file does not exist?");
}

/**
   Read an ASCII stream and append the contents to the String array
   text.  TEXT IS NOT OVERWRITTEN, BUT APPENDED!

   @param text Output. The contents fo the file
   @param is Stream from which to read
   @exception IOError Some error occured during the read
   @version   1
   @author Stefan Buehler */
ArrayOfString read_text_from_stream(istream& is) {
  ArrayOfString text;
  String linebuffer;

  // Read as long as `is' is good.
  // Contary to what I understood from the book, the explicit check
  // for eof is necessary here, otherwise the last line is read twice
  // if it is not terminated by a newline character!
  while (is && is.good() && !is.eof()) {
    // Read line from file into linebuffer:
    getline(is, linebuffer);

    // Append to end of text:
    text.push_back(linebuffer);
  }

  // Check for error:
  // FIXME: This should not be necessary anymore when stream
  // exceptions work properly.
  ARTS_USER_ERROR_IF(!is.eof(), "Read Error. Last line read:\n", linebuffer);

  return text;
}

/**
   Reads an ASCII file and appends the contents to the String vector
   text. This uses the function @see read_text_from_stream. TEXT IS
   NOT OVERWRITTEN, BUT APPENDED!  

   \param text  Output. The contents fo the file
   \param name  Name of file to read
   \exception IOError
   \version   1
   \author Stefan Buehler */
ArrayOfString read_text_from_file(const std::string_view name) {
  ArrayOfString text{};
  ifstream ifs;

  // Open input stream:
  open_input_file(ifs, name);
  // No need to check for error, because open_input_file throws a
  // runtime_error with an appropriate error message.

  // Read the text from the stream. Here we catch the exception,
  // because then we can issue a nicer error message that includes the
  // filename.
  try {
    text = read_text_from_stream(ifs);
  } catch (const std::runtime_error& x) {
    ARTS_USER_ERROR("Error reading file: ", name, "\n", x.what());
  }

  return text;
}

/**
    Replace all occurances of `what' in `s' with `with'.

    @param s Output. The String to act on.
    @param what The String to replace.
    @param with The replacement.

    @author Stefan Buehler */
void replace_all(String& s, const std::string_view what, const std::string_view with) {
  Index j = s.find(what);
  while (j != s.npos) {
    s.replace(j, 1, with);
    j = s.find(what, j + with.size());
  }
}

/**
  Checks if there is exactly one newline character
  at the end of the string.

  @param s The String to check.

  @return Error code (0=ok, 1=empty, 2=missing, 3=extra newline

  @author Oliver Lemke
*/
int check_newline(const std::string_view s) {
  String d = s;
  int result = 0;

  // Remove all whitespaces except \n
  replace_all(d, " ", "");
  replace_all(d, "\t", "");
  replace_all(d, "\r", "");

  const char* cp = d.c_str();
  while (*cp == '\n') cp++;

  if (!(*cp)) result = 1;

  if (!result && d[d.length() - 1] != '\n')
    result = 2;
  else if (!result && d.length() > 2 && d[d.length() - 1] == '\n' &&
           d[d.length() - 2] == '\n')
    result = 3;

  return result;
}

/**
  Checks if the given file exists.

  @param filename File to check.

  @return Error code (true = file exists, false = file doesn't exist)

  @author Oliver Lemke
*/
bool file_exists(const std::string_view filename) {
  return std::filesystem::exists(filename) &&
         !std::filesystem::is_directory(filename);
}

/**
  Searches through paths for a file with a matching name.

  If the filename starts with '/', the search path is ignored.

  @param[in,out] matches     Matching files are appended to this list.
  @param[in]     filename    File to find.
  @param[in]     paths       List of paths to look in for the file.
  @param[in]     extensions  List of extensions to add to base filename.

  @return True if matches were found, else false

  @author Oliver Lemke
*/
bool find_file(ArrayOfString& matches,
               const std::string_view filename,
               const ArrayOfString& paths,
               const ArrayOfString& extensions) {
  bool exists = false;
  std::string efilename{expand_path(filename)};

  // filename contains full path
  if (!paths.nelem() || std::filesystem::path(efilename).is_absolute()) {
    for (const auto& ext : extensions) {
      const String fullpath{efilename + ext};

      if (file_exists(fullpath)) {
        if (std::find(matches.begin(), matches.end(), fullpath) ==
            matches.end())
          matches.push_back(fullpath);
        exists = true;
      }
    }
  }
  // filename contains no or relative path
  else {
    for (const auto& path : paths) {
      for (const auto& ext : extensions) {
        const String fullpath{expand_path(path) + "/" + efilename + ext};

        if (file_exists(fullpath)) {
          if (std::find(matches.begin(), matches.end(), fullpath) ==
              matches.end())
            matches.push_back(fullpath);
          exists = true;
        }
      }
    }
  }

  return exists;
}

/**
  Find an xml file.
 
  If it doesn't exist in the current directory, also search the include
  and data path. Also tests if a compressed version exists.

  The filename will be modified to contain the full path to the found match.

  @param[in,out] filename   File to check.

  @throw runtime_error if file is not found.

  @author Oliver Lemke
*/
void find_xml_file(String& filename) {
  // Command line parameters which give us the include search path.
  extern const Parameters parameters;
  ArrayOfString allpaths = parameters.includepath;
  allpaths.insert(
      allpaths.end(), parameters.datapath.begin(), parameters.datapath.end());

  ArrayOfString matching_files;
  find_file(matching_files, filename, allpaths, {"", ".xml", ".gz", ".xml.gz"});

  ARTS_USER_ERROR_IF(!matching_files.nelem(),
                     "Cannot find input file: ",
                     filename,
                     "\nSearch path: ",
                     allpaths);

  filename = matching_files[0];
}

/** As find_xml_file but does not throw in the main body
 * 
 * The filename will be modified to contain the full path to the found match
 * if there is a match
 * 
 * @param[in,out] filename   File to check.
 * 
 * @return If file is found
 * 
 * @author Oliver Lemke
 */
bool find_xml_file_existence(String& filename) {
  // Command line parameters which give us the include search path.
  extern const Parameters parameters;
  ArrayOfString allpaths = parameters.includepath;
  allpaths.insert(
    allpaths.end(), parameters.datapath.begin(), parameters.datapath.end());
  
  ArrayOfString matching_files;
  find_file(matching_files, filename, allpaths, {"", ".xml", ".gz", ".xml.gz"});
  
  if (matching_files.nelem()) {
    filename = matching_files[0];
    return true;
  }

  return false;
}

/*!
 Expands the ~ to home directory location in given path.
 
 \param[in] path  String with path
 
 \return Expanded path.
 
 \author Oliver Lemke              
 \date   2010-04-30
 */
String expand_path(String path) {
  if ((path.length() == 1 && path[0] == '~') ||
      (path.length() > 1 && path[0] == '~' && path[1] == '/')) {
    path = String(std::getenv("HOME")) + String(path, 1);
  }

  return path;
}

/*!
 Adds base directory to the given path if it is a relative path.
 
 \param[in] path  String with path
 
 \return Expanded path.
 
 \author Oliver Lemke              
 \date   2012-05-01
 */
String add_basedir(const std::string_view path) {
  extern const Parameters parameters;
  String expanded_path{expand_path(path)};

  if (parameters.outdir.nelem() && path.length() && path[0] != '/') {
    expanded_path = parameters.outdir + '/' + expanded_path;
  }

  return expanded_path;
}

//! Return the parent directory of a path.
/** 
 Extracts the parent directory part of the given path.
 
 \param[out]  dirname  Parent directory of path
 \param[in]   path     Path
 
 \author Oliver Lemke
 */
String get_dirname(const std::string_view path) {
  String dirname{};
  if (!path.length()) return dirname;

  ArrayOfString fileparts;
  String{path}.split(fileparts, "/");
  if (path[0] == '/') dirname = "/";
  if (fileparts.nelem() > 1) {
    for (Index i = 0; i < fileparts.nelem() - 1; i++) {
      dirname += fileparts[i];
      if (i < fileparts.nelem() - 2) dirname += "/";
    }
  }

  return dirname;
}

//! Return list of files in directory.
/**
 Returns list of all files in the given path.

 \param[out] files    List of files
 \param[in]  dirname  Parent directory of path

 \author Oliver Lemke
 */
ArrayOfString list_directory(const std::string_view dirname) {
  ARTS_USER_ERROR_IF(!std::filesystem::is_directory(dirname),
                     "Error accessing directory: ",
                     dirname)

  ArrayOfString files{};
  for (const auto& filename :
       std::filesystem::directory_iterator{dirname}) {
    files.push_back(filename.path().string());
  }

  return files;
}

/** Make filename unique.

 Checks if a file (or a gzipped version of it) with the given name already
 exists und appends a unique number to the filename if necessary

 \param[in] filename   Filename
 \param[in]     extension  Optional, number is inserted before the extension
 \return Unique filename

 \author Oliver Lemke
 */

String make_filename_unique(const std::string_view filename, const String& extension) {
  String basename = filename;
  String extensionname;

  // Split filename into base and extension
  if (extension.length()) {
    size_t pos = filename.rfind(extension);
    if (pos == filename.length() - extension.length()) {
      basename = filename.substr(0, filename.length() - extension.length());
      extensionname = extension;
    }
  }

  Index filenumber = 0;
  ostringstream newfilename;
  newfilename << basename << extensionname;

  while (file_exists(newfilename.str()) ||
         file_exists(newfilename.str() + ".gz")) {
    filenumber++;
    newfilename.str("");
    newfilename << basename << "." << filenumber << extensionname;
  }

  return newfilename.str();
}
