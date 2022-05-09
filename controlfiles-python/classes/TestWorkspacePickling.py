#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 16:31:41 2022

@author: u237023
"""

import pickle
import pyarts
import pyarts.arts as cxx

ws = pyarts.workspace.Workspace()
x = list(cxx.get_WsvGroupMap().keys())

for i in range(len(x)):
    print(f"Create {x[i]} on workspace", end='; ')
    exec(f"ws.v{i} = cxx.{x[i]}()")
    
    print(f"pickling the workspace", end='; ')
    pickle.dump(ws, open("test.pcl", 'wb'))
    
    print(f"unpickling the workspace")
    ws2 = pickle.load(open("test.pcl", 'rb'))

x = 4
assert ws.number_of_initialized_variables() - \
    ws2.number_of_initialized_variables() == x, \
        f"""
There should be {x} more initd vars because there are not picklable by design:
    - Agenda (needs workspace variables in fixed positions)
    - ArrayOfAgenda (as above)
    - CallbackFunction (needs code)
    - ArrayOfRetrievalQuantity (always set together with agenda)
"""
