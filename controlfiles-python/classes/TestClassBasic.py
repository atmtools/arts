"""
Note to developer encountering this error.  All workspace classes in ARTS
should be available in pyarts.  You have probably added a new group.  Unless
your group is very special it should support the following operations


import pyarts.pyarts_cpp as cxx
ws = cxx.Workspace()

x = cxx.Group()
ws.GroupCreate("Group__1", '', x)
ws.Group__2 = x
cxx.Group(ws.Group__1)
cxx.Group(ws.Group__2)
assert hasattr(x, "readxml")
assert hasattr(x, "loadxml")
assert hasattr(x, "__copy__")
assert hasattr(x, "__deepcopy__")
assert hasattr(x, "__str__")
assert hasattr(x, "__repr__")
"""

import pyarts.pyarts_cpp as cxx
ws = cxx.Workspace()

# Special groups
special_groups = ["Any"]

# All groups
list_of_groups = cxx.get_wsv_group_names()

# Test initialization
for g in list_of_groups:
    if g in special_groups: 
        print("Ignoring {}\n".format(g))
        continue

    try:
        print("Running tests for", g)
        
        print("Trying to defaul init")
        x = eval("cxx.{}()".format(g))
        
        print("Trying create workspace variable")
        eval("ws.{}Create('{}__1', '', x)".format(g, g))
        exec("ws.{}__2 = x".format(g))
        
        print("Trying create fron workspace variable")
        eval("cxx.{}(ws.{}__1)".format(g, g))
        eval("cxx.{}(ws.{}__2)".format(g, g))
        
        print("Checking File IO")
        assert hasattr(x, "readxml")
        assert hasattr(x, "savexml")
        
        print("Checking copy constructor")
        assert hasattr(x, "__copy__")
        assert hasattr(x, "__deepcopy__")
        
        print("Checking string and representation")
        assert hasattr(x, "__str__")
        assert hasattr(x, "__repr__")
        
        print("Success!\n")
    except:
        raise ImportError("Incomplete pyarts interface for {} in pyarts.classes".format(g))

""" Each class must be tested

If you encounter this error, add a Test{GROUP}.py for the group that is missing
in the same folder as this file is located
"""
if __name__ == "__main__":
    import os
    dir = os.path.dirname(__file__)
    cmake_path=os.path.join(os.path.dirname(dir), "CMakeLists.txt")
    cmake = open(cmake_path, 'r').read()
    
    for x in list_of_groups:
        if not os.path.exists(os.path.join(dir, f"Test{x}.py")):
            raise RuntimeError(f"Must have Test{x}.py in {dir}")
            
        run = f"arts_test_run_pyfile(fast classes/Test{x}.py)"
        if not run in cmake:
            raise RuntimeError(f"Must have \"{run}\" in {cmake_path}")
            
        