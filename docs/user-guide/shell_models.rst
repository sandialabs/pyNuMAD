Making Shell Models
===================


Currently, shell models can be made using ANSYS or Abaqus.  Examples can be found in pyNuMAD/examples:

`ansys_analysis.py` runs an analysis using ANSYS with a model generated by the in-house mesher.

`write_abaqus_shell_model.py` generates an abaqus input file using the in-house mesher, which can be imported and analyzed using Abaqus CAE.

For shell models, the in-house mesher takes an option for whether to include the trailing edge adhesive, meshed with solid elements.  If included, the output gives constraint equations to tie the motion of the adhesive together with the blade.  write_abaqus_shell_model.py demonstrates the access and usage of this.