.. _gettingStarted:

Getting Started 
================

Software Requirements
----------------------
NuMAD is developed in MATLAB, and requires the following software packages:

===========================================================================================	=============================
**Required Packages**        						    			**Oldest Compatible Version**
`PYTHON (Required) <https://www.python.org/>_`                        v3.?
`FAST (Optional) <https://www.nrel.gov/wind/nwtc/fastv7.html>`_ 	    			v7.02.00d
`ANSYS (Optional) <https://www.ansys.com/>`_	    			    			R12
===========================================================================================	=============================


.. Unsupported IECWind, ADAMS
 

Download NuMAD
----------------

The NuMAD source code is hosted on the `NuMAD GitHub repository <https://github.com/sandialabs/NuMAD>`_. 
NuMAD users are recommended to clone the Github repository.
Cloning the repository allows users to easily pull the latest updates to the NuMAD source code.
These updates may improve the code's speed, accuracy and add additional functionality or advanced features.

.. TODO: this section doesn't exist
.. Developers who wish to contribute to NuMAD should see the corresponding Developer :ref:`dev-getting-started` section.

 
To install NuMAD using `git <https://git-scm.com/>`_, type the following in a git interface:: 

    >> git clone https://github.com/sandialabs/NuMAD.git

The local copy of NuMAD can easily be updated to the latest version of the 
code hosted on the GitHub repository by using the pull command:: 

    >> git pull

For users who are new to git, it is recommended to go through examples on 
`GitHub <https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github>`_ 
or other sources while getting started. 
If you have problems downloading or installing please see the :ref:`troubleshooting` page.

.. Note::
    Users may also download a static version of NuMAD from the latest tagged 
    `NuMAD Release <https://github.com/sandialabs/NuMAD/releases>`_.  This is 
    the easiest way to obtain the NuMAD code, however it is more difficult to 
    manually download future updates.


.. _user-install:

Install NuMAD
---------------

Once you have downloaded the NuMAD source code, take the following steps to 
install the NuMAD code. The directory where the NuMAD code is contained is 
referred to as ``$NuMAD`` (e.g. ``C:/User/Documents/GitHub/NuMAD``). 



Installation Steps
~~~~~~~~~~~~~~~~~~

1.    Install `MATLAB <https://www.mathworks.com/products/matlab.html>`_
2.    Open ``addNumadPaths.m`` in an editor and set the ``numadPath`` to point to the source directory in the repo.
3.    If you plan to use ANSYS or other analysis programs, update the corresponding path variable to the actual path of each executable you plan to use. Save and close ``addNumadPaths.m.``
4.    Move the ``addNumadPaths.m`` file to the MATLAB directory of your user account (e.g. ``C:\Users\username\Documents\MATLAB``)

Every time a new MATLAB session is started type ``addNumadPaths`` to be able to run NuMAD.




