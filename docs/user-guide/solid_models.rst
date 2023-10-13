===================
Making Solid Models
===================


Pure solid model (explicitly descritezed sandwich panels) 
=========================================================

Pure solid models can be made in Sierra SM, Sierra SD or Abaqus.
Currenly only SD is partially supported since it does not allow for 
spatially varying material orientations. 

Continuous meshes in Sierra
----------------------------

#. The first step is to make a Genesis mesh file and the associated files with
   :py:func:`~pynumad.analysis.cubit.make_blade.cubit_make_solid_blade`. This 
   will create the following files

      * euler. Binary file with euler angles for material orientations
      * mat_ori.py
      * {wt_name}.g. Genesis mesh file. 

#. Make the Sierra input files with :py:func:`~pynumad.analysis.make_models.write_sierra_model`. 

   This creates the following files:

      * sm.i and/or sd.i: Sierra SM or SD input file

    Be sure to have ran the following import statement

    .. code-block:: python

       from pynumad.analysis.make_models import write_sierra_model



#. Apply the spatially varying material orientations to the Genisis mesh. This requires SEACAS 
   to be installed. One HPWS, you can issue `module load seacas`. Then run mat_ori.py from 
   a terminal. 

   .. Note:: 
      Reliance on SEACAS for placing the material orientations in the Genesis file is
      a temporary solution until this capability is added to Cubit.

#. Finally run Sierra SM with: launch -n 10 adagio -i sm.i, where n is the 
   number of CPUs.


An example called `cubit_solid.py` exists in the examples folder.

Discontinuous meshes in Abaqus
------------------------------

Layered solid model (homogenized sandwich panels)
==================================================

This capability does not yet exist.


