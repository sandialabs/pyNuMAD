.. _intallation:

Installation 
============

Download pyNuMAD
----------------

The pyNuMAD source code is hosted on the `pyNuMAD GitHub repository <https://github.com/sandialabs/pyNuMAD>`_. 
pyNuMAD users are recommended to clone the Github repository.
Cloning the repository allows users to easily pull the latest updates to the pyNuMAD source code.
These updates may improve the code's speed, accuracy and add additional functionality or advanced features.

.. TODO: this section doesn't exist
.. Developers who wish to contribute to pyNuMAD should see the corresponding Developer :ref:`dev-getting-started` section.

To download pyNuMAD using `git <https://git-scm.com/>`_, type the following in a git interface:: 

    git clone https://github.com/sandialabs/pyNuMAD

Installation
------------

After downloading the source, pyNuMAD can be installed by running
the following command in the root of the repository::

    pip install -e .

Currently, it is necessary to use `-e` for a local installation so that certain data files are not lost in the installation process.

Developers are recommended to install using the instructions on
:ref:`contributing<contributing>` page.
