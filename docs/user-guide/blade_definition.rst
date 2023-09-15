.. _bladedefinition:

Blade Object
=================

The fundamental class of pyNuMAD is the Blade class.
Each blade houses a collection of subobjects as attributes 
which organize the various parameters and data associated blade
in a logical fashion. In what follows, each of these primary
attributes of the blade are explained at a high-level. For
more information, please refer to the API documentation.


Definition
------------

A Definition object provides attributes for the basic design
of a wind turbine blade. Typically this is populated by
a yaml file or an excel file, however it is possible to build a blade
from scratch by manually assigning all of the attributes.
Once a definition object has been loaded in either manually or from
a file, a user can make additional modifications to the blade, such as
adding additional stations or changing material assignments, before
generating downstream data.

In pyNuMAD, a blade is defined with the ``Definition`` object, or blade
object for short.S
many of the properties are parameterized by spanwise location. Refer to
X for a complete listing of ``BladeDef`` properties.

First and foremost there are *stations*. A station is an airfoil at a
specified span location. The airfoil is partitioned by *keypoints*,
shown in X. Various blade properties such as ``blade.leband``,
``blade.teband``, ``blade.sparcapwidth``, and ``blade.sparcapoffset`` help to
position the keypoints precisely. For example, ``blade.leband`` is the
arclength from the *le* keypoint to the keypoint *a*. *Regions* are
defined between the keypoints as listed in X. An adjacent
station helps define these regions as areas. Spanwise lines emanating
from each keypoint are connected to the corresponding keypoints on an
adjacent station; thus bounding the region with four curves. A suffix of
either HP or LP is added to each region name to distinguish regions on
the high pressure surface verses the low pressure surface. Other airfoil
properties and external blade shape data are defined with the ``AirfoilDef``
class and the ``StationDef`` object respectively.
Usually, the number of stations defined needs to be supplemented for
with interpolated stations.

Material properties, layup information, and thicknesses and widths are
additionally defined in the ``MaterialDef``, ``StackDef``, and ``ComponentDef`` respectively.
Refer to the :ref:`classDefs` for more information.

Geometry
--------

Typically a blade definition does not contain
high enough fidelity data for creating a mesh, so
pyNuMAD performs additional interpolation to
create a more detailed geometry. The data generated
from this process is stored in a Geometry object

Keypoints
----------

.. TODO

.. _bladeKeyPoints:
.. figure:: /_static/images/bladeKeyPoints.png
   :width: 6.5in
   :height: 2.23056in

   Relative locations of the blade keypoints.
   
   
.. _defineRegions:
.. table:: Region definition by keypoints (TE-Trailing edge, LE-leading edge)

    +----------------------------------+-----------------------------------+
    | Region Name                      | Bounding Keypoints                |
    +==================================+===================================+
    | LE                               | le & a                            |
    +----------------------------------+-----------------------------------+
    | LE Panel                         | a & b                             |
    +----------------------------------+-----------------------------------+
    | Spar                             | b & c                             |
    +----------------------------------+-----------------------------------+
    | TE Panel                         | c & d                             |
    +----------------------------------+-----------------------------------+
    | TE REINF                         | d & e                             |
    +----------------------------------+-----------------------------------+
    | TE Flatback                      | e & te                            |
    +----------------------------------+-----------------------------------+
