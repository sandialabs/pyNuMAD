.. _reference:


Reference
============

This page contains definitions for various terminology and abbreviations
used throughout pyNuMAD documentation and code. 

Terminology
-----------


Airfoil: A unitless (0-1) 2D outline of the desired aerodynamic shape (i.e., NACA-63-214, DU99-W-405, etc). 
These are defined at every span location (different from ispan).

Profiles: Usually there are not enough airfoils defined along the span to create a smooth 3D blade geometry. 
Profiles are unitless airfoil geometries at each blade interpolated span (ispan) location.

Station: A station is an airfoil object at a specified span location.

Stack: A stack of material plys specified by types and thicknesses.

Ply: A layer of material with defined properties, layup angle, and thickness.

Camber: A 2D collection of points midway between the HP curve and LP curve of an airfoil or profile.

Keypoint: Specified locations along the profile circumference used to define where there is a change in stack definition.

Component: TODO

Shearweb: TODO

Span: TODO

Chord: TODO

Twist: TODO

Prebend: TODO

Sweep: TODO

TE Type: TODO


Abbreviations
-------------

BOM: Bill of Materials

OML: Outer Mold Line

TE: Trailing Edge

LE: Leading Edge

PS: Pressure Side (high pressure)

SS: Suction Side (low pressure)