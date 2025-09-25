Guide
=====

Here you find an overview of how different components of ARTS are intended to work.

There are several sections below that explain in some detail how to understand important
data structures and concepts in ARTS.  These sections focus on important topics of 
atmospheric radiative transfer that a user needs to understand to use ARTS successfully,
such as how the surface and the atmosphere are defined and represented in ARTS. 

The last subsection below is a collection of concepts that are required to understand
the background of some of the design decisions in ARTS.  It does not cover
implementation details, but it is useful to understand how different components
of ARTS are intended to work together.  These are generic topics that contains
theoretical constructs that are often implicitly used in ARTS to allow different
methods to work together on the assumption that they all represent a coherent
theoretical framework.

The subsections below are given in no particular order.

.. toctree::
   :maxdepth: 2
   
   user.atmospheric_field
   user.surface_field
   user.subsurface_field
   user.sensors
   
   concepts
