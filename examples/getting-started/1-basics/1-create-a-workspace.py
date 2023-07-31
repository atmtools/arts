# The first step is to import the pyarts package. This will make all
# pyarts modules available to you.
import pyarts

# The Workspace class is the common entry point for all pyarts simulations.
# It is used to interface with the ARTS C++ library and to manage the
# simulation state.
ws = pyarts.Workspace()

# That's it. You have created a Workspace. You can now use it to run
# simulations. For example, you can run the ARTS Hello World example:
ws.Print("Hello, World!", 0)
