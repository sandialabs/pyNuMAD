.. _blade-overview:

Blade Overview
==============

The fundamental class of pyNuMAD is the Blade class.
Each blade object houses a collection of subobjects as attributes 
which organize the various parameters and data of the blade
in a logical fashion. In what follows, each of these primary
attributes of the blade are explained at a high-level. For
more information, please refer to the API documentation.
The below figure illustrates the dependencies between various
blade subobjects. For example, a definition object is required
to generate geometry, keypoints, bill of material, and material database.

.. _blade-tree:
.. figure:: /_static/images/blade_attribute_tree.png

   Dependency tree for blade subobjects

Definition (blade.definition)
-----------------------------

A Definition object provides attributes for the basic design
of a wind turbine blade - i.e. where the blade is *defined*. 
Typically this is populated by
a yaml file or an excel file, however it is possible to build a blade
from scratch by manually assigning all of the necessary attributes.
Once a definition object has been loaded in either manually or from
a file, a user can make additional modifications to the blade, such as
adding additional station locations or changing material assignments, before
generating downstream data.

Many of the attributes in Definition are parameterized by spanwise location.
For example, *stations* are airfoils at specified span locations. 
Other airfoil properties and external blade shape data are 
defined with the ``Airfoil``
class and the ``Station`` object respectively, and are stored in ``definition.stations``.
Material properties, layup information, and thicknesses and widths are
defined in the ``Material`` and ``Component`` classes 
and stored in ``blade.materials`` and ``blade.components``.

Geometry (blade.geometry)
---------------------------

Typically, the blade definition does not contain
high enough fidelity data for creating a mesh, so
pyNuMAD performs additional interpolation to
create a more detailed geometry. The Geometry class generates
and stores the interpoloated geometry. 

Keypoints (blade.keypoints)
---------------------------

The Keypoints class generates and organizes data related
to keypoints.
Airfoils are partitioned by keypoints,
as shown in :numref:`keypoints-fig`. Various definition properties such as 
``definition.leband``,
``definition.teband``, ``definition.sparcapwidth``, 
and ``definition.sparcapoffset`` help to
position the keypoints precisely. For example, ``definition.leband`` is the
arclength from the *le* keypoint to the keypoint *a*. *Regions* are
defined between the keypoints as listed in :numref:`define-regions`.
Adjacent stations help to define these regions as areas. Spanwise lines emanating
from each keypoint are connected to the corresponding keypoints on an
adjacent station; thus bounding the region with four curves. A suffix of
either HP or LP is added to each region name to distinguish regions on
the high pressure surface verses the low pressure surface. 


.. _keypoints-fig:
.. figure:: /_static/images/keypoints.png

   Keypoint locations
   
   
.. _define-regions:
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




StackDB (blade.stackdb)
-----------------------
.. TODO

BillOfMaterials (blade.bom)
---------------------------
The BillOfMaterials class organizes material types and quantities in
a blade. However, this class is not currently
used in pyNuMAD analyses and only exists for legacy reasons.


MaterialDB (blade.materialdb)
-----------------------------
The MaterialDB class is another class for organizing materials. 
However, this class also is not currently
used in pyNuMAD analyses and only exists for legacy reasons.
