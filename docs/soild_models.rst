Solid Models
=============


Pure solid model (explicitly descritezed sandwich panels) 
---------------------------------------------------------

Pure solid models can in Sierra SM or SD. After meshing 
cross sections with :py:func:`~pynumad.analysis.cubit.make_blade.cubit_make_solid_blade`, call on 
:py:func:`~pynumad.analysis.make_models.write_sierra_model`. Be sure to issue the following 
import statement

.. code-block:: python

   from pynumad.analysis.make_models import write_sierra_model

An example called `cubit_solid.py` exists in the examples folder.




Layered solid model (homogenized sandwich panels)
-------------------------------------------------

This capability does not yet exist.


