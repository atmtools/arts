"""
This example shows how to load data from a file into a workspace variable.
"""

# Import the module
import pyarts


# Create a workspace
ws = pyarts.Workspace()


"""
Before we load data into a workspace variable, we need to know about
one of the tricks that ARTS uses to make the file IO more flexible.

The trick is that ARTS does not just look for files in the folder
that you specify, but also in a number of other folders.  This is
useful since a lot of radiative transfer data is the same regardless
of the scenario that you are simulating.  For example, the absorption
lines of H2O are the same regardless of whether you are simulating
a clear sky or a cloudy sky.

The folders that ARTS looks in are called search paths.  You can
print the search paths that ARTS uses by doing:
"""
print ("ARTS search paths:", pyarts.arts.globals.parameters.datapath)


"""
The results of this command will depend on your environment variables
at the time of first import of pyarts during the scripting.  The default
search paths can be changed by setting the ARTS_DATA_PATH environment
variable.  For example, if you want to add a folder called "my_data"
to the search paths, you can do:

export ARTS_DATA_PATH=/my_data

before starting python.  If you want to add multiple folders, you can
separate them with a colon:
export ARTS_DATA_PATH=/my_data:/my_other_data

Check the documentation for your shell if you want to make this change
permanent.

It is always possible to add additional search paths from within python:
"""
pyarts.arts.globals.parameters.datapath.append("/my_data")
print ("ARTS search paths:", pyarts.arts.globals.parameters.datapath)


"""
The following examples assumes that you have a copy of the arts-cat-data
and arts-xml-data repositories in the search paths.  If you do not have
these, you can download them from: https://www.radiativetransfer.org/tools/

It demonstrates different ways that you can load data from a file into
the workspace.  There is no "best" way to do this, it depends on the
context, and what you find most readable.
"""

# Call the WorkspaceVariable member method "readxml" to load data from file
ws.example_line_list = pyarts.arts.ArrayOfAbsorptionLines()  # Create an empty object
ws.example_line_list.readxml("lines/O2-66.xml")

# Call the Workspace Method "ReadXML" to load data from file
ws.example_griddedfield3 = pyarts.arts.GriddedField3()
ws.ReadXML(ws.example_griddedfield3, "planets/Earth/Fascod/tropical/tropical.t.xml")

# Call the pure python function "pyarts.xml.load" to load data from file
ws.example_griddedfield2 = pyarts.arts.GriddedField2()
ws.example_griddedfield2 = pyarts.xml.load("star/Sun/solar_spectrum_May_2004.xml")

"""
This example is a bit special.  The pyarts.xml.load function creates an
instance of the type that it finds in the file, and then loads the data into
that instance.  This means that you can also just duck type the variable:

>>> ws.example_griddedfield2_2 = pyarts.xml.load("star/Sun/solar_spectrum_May_2004.xml")
"""

# TESTING
# AS THIS FILE IS RUN TO TEST ARTS, WE NEED TO CHECK THAT
# THE CONTENT OF THE VARIABLES ARE GOOD.
assert pyarts.arts.globals.parameters.datapath.index("/my_data")
# END TESTING
