from inspect import isclass, ismodule
import pyarts
import sys


class DocData:
    def __init__(self, c, f, m):
        self.classes = c
        self.funcs = f
        self.modules = m

    def print_classes(self, level=0):
        self.classes.sort()
        if len(self.classes) == 0: return ""
        
        out = f"""
Classes
{'-'*7 if level == 0 else '^'*7}

.. autosummary::
   :toctree: stubs

"""
           
        for name in self.classes:
            out  += f"   {name}\n"
        return out


    def print_funcs(self, level=0):
        self.funcs.sort()
        if len(self.funcs) == 0: return ""
        
        out = f"""
Functions
{'-'*9 if level == 0 else '^'*9}

.. autosummary::
   :toctree: stubs

"""
           
        for name in self.funcs:
            out  += f"   {name}\n"
        return out

    def print_submodules(self):
        keys = list(self.modules.keys())
        keys.sort()
        
        out = ""
        for submodule in keys:
            docdat = self.modules[submodule]
            if len(docdat.modules): print(f"""WARNING:
{submodule} has submodules

This is not considered by the auto-generated documentation""", file=sys.stderr)
            
            out += f"""
pyarts.arts.{submodule}
------------{'-'*len(submodule)}

.. automodule:: pyarts.arts.{submodule}

.. currentmodule:: pyarts.arts.{submodule}

""" + docdat.print_classes(1) + "\n\n" + docdat.print_funcs(1) + "\n\n" 
        return out
    
    def __str__(self):
        return """pyarts.arts
===========

.. automodule:: pyarts.arts

.. currentmodule:: pyarts.arts

""" + self.print_classes() + "\n\n" + self.print_funcs() + "\n\n" + self.print_submodules()


def doc_split(mod):
    classes = []
    functions = []
    submodules = {}
    for item in dir(mod):
        if item[0] == "_:" or item[-1] == "_" or  ":" in item: continue
        
        attr = getattr(mod, item)
        if isclass(attr): classes.append(item)
        elif ismodule(attr): submodules[item] = doc_split(attr)
        else: functions.append(item)
    
    return DocData(classes, functions, submodules)

print(doc_split(pyarts.arts))



