.. _overview:

Overview
=======================


`pyNuMAD (Python Numerical Manufacturing And Design) <https://github.com/sandialabs/pyNuMAD>`_ simplifies the process of creating 
structural models of wind turbine blades at various fidelities. The tool manages information such as blade geometry,
material properties, layups, and layup placement. By interfacing with the following codes:
-   `ANSYS Mechanical® <http://www.ansys.com/>`__ 
-  `Abaqus <https://www.3ds.com/products-services/simulia/products/abaqus/>`__
-  `Cubit <https://cubit.sandia.gov/>`__
-  `VABS <https://analyswift.com/vabs-cross-sectional-analysis-tool-for-composite-beams/>`__
-  `BeamDyn <https://openfast.readthedocs.io/en/dev/source/user/beamdyn/index.html>`__

pyNuMAD can create beam, shell, or solid models. 

The code was created by converting select functionalities from the well-known NuMAD tool. The move to Python was motivated
by increasing code accessibility as well as to aid with integration with other WETO tools. Sandia developers have transitioned
the NuMAD repository to a static code. The graphical user interphase is currently unsupported. However new features have been implemented.
Namely, 

- In-house mesher for pure shell, pure solid, or shell + solid adhesive
- Cubit mesher for contiguous solid meshes and cross-sectional analysis
 
Supported input files:

- yaml files from the `IEA Ontology <https://windio.readthedocs.io/en/latest/>`__
- Excel files from the legacy version of NuMAD 


   

