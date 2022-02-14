import pyarts


list_of_groups = pyarts.classes.list_of_workspace_variables()
for group in list_of_groups:
    if group == "CallbackFunction": continue
    try:
        if group in ["Index", "Numeric", "String"]:
            eval("pyarts.classes.{}(1)".format(group.val))
        else:
            eval("pyarts.classes.{}()".format(group.val))
    except:
        raise ImportError("Cannot use {} in pyarts.classes".format(group))


"""
Note to developer encountering this error.  All workspace classes in ARTS
should be available in pyarts.  You have probably added a new group.  Unless
your group is special and have different initialization in different build
modes (as indicated by the list of special groups above) you need to be able
to evaluate "pyarts.classes.YOURGROUP()" to show that the new group works in
pyarts.

The easiest way to add your new group is to look at one of the python/classes/*
examples and copy it (python class + C-API).  Furthermore, in 
python/classes/__init__.py, your new class should be imported.  Any subclass
you want can be supported in this manner as well.  All classes added this way
should preferably have a TestYOURGROUP.py file in controlfiles-python/classes/,
Subclasses can be tested here as well.

In C++, you will also need to create the C-interface.  There are plenty of
helper macros available in arts_api_classes.h/cc to facilitate this addition,
but depending on your class you might have to do more to make some of the
internals work.  Be careful, since setting a python class with "void *"
pointers might mean that it is a true reference, we can only actually delete
the C-data when we know it is not such a reference.  Functions that return
"void *" pointers must thus manually manage memory, which is a great annoyance.

Anyways... With hope, and take care!
"""
