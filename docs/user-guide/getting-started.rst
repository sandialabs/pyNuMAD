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

**To Configure Cubit**
After downloading the main Cubit folder and placing it anywhere on your machine:

* Locate `cubit.py` in the main Cubit folder and add the full path to the "cubit" entry in `src/pynumad/software_paths.json`.
On Linux this may look something like: "/some/place/Cubit-16.10/bin/". If Mac it may be: "/Applications/Cubit.app/Contents/MacOS/".

* Locate `PyCubed_Main.py` in the main Cubit folder and add the full path to the "cubitEnhancements" entry in `src/pynumad/software_paths.json`.
On Linux this may look something like: "/some/place/Cubit-16.10/bin/python3/lib/python3.7/site-packages/". On Mac it may be: "/Applications/Cubit.app/Contents/MacOS/python3/lib/site-packages/".

**To Configure ANSYS**
Locate the `mapdl` executable in your ANSYS installation folder. Then place the name of the executable along with its path in the "ansys_path"
entry of `src/pynumad/software_paths.json`. For example, it may look like: "/some/where/ansys_inc/v222/ansys/bin/mapdl"