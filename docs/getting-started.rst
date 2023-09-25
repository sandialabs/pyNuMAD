.. _getting-started:

Getting Started
===============

Initializing a blade object
---------------------------

After a successful installation, pyNuMAD can be imported in a python
session in the usual way::

    import pynumad

Next, to initialize an empty blade object run::

    blade = pynumad.Blade()

To populate the blade with a yaml or Excel file you can run::

    blade.load_yaml("path/to/yaml")
    blade.load_excel("path/to/excel")

or initialize the blade with the path::

    blade = pynumad.Blade("path/to/yaml")
    blade = pynumad.Blade("path/to/excel")


Interfacing with other software
--------------------------------

Depending on your needs, you might need to add installation paths for third-party software such as ANSYS, Cubit, etc...
To do so, locate the file named `src/pynumad/software_paths.json` and modify it accordingly. Note that these paths 
are optional if you only need to read in a blade as described above and/or using the in-house mesher.