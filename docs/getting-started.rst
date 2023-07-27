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

To populate the blade with a yaml or excel file you can run::

    blade.load_yaml("path/to/yaml")
    blade.load_excel("path/to/excel")

or initialize the blade with the path::

    blade = pynumad.Blade("path/to/yaml")
    blade = pynumad.Blade("path/to/excel")