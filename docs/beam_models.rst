.. _meshing:

Beam Models
==================================

Currently, only BeamDyn powered by VABS beam properties is supported. After meshing 
cross sections with :py:func:`~pynumad.analysis.cubit.make_blade.cubit_make_cross_sections`, call on 
:py:func:`~pynumad.analysis.make_models.write_beam_model`. Be sure to issue the following 
import statement

.. code-block:: python

   from pynumad.analysis.make_models import write_beam_model


An example called `cubit_beam.py` exists in the examples folder.