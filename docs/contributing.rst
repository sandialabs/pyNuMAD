.. _contributing:


Contributing
============

Thank you for considering contributing to pyNuMAD. 
We welcome contributions from the community in the form of bug fixes, feature
enhancements, documentation updates, etc. All contributions are processed
through pull-requests or issues on GitHub. Please follow these guidelines for contributing.

Reporting issues and bugs
-------------------------

This section guides you through the process of submitting an issue for pyNuMAD.
To report issues or bugs please `create a new
issue <https://github.com/sandialabs/pyNuMAD/issues/new>`_ on GitHub.
	
Following these guidelines will help maintainers understand your issue,
reproduce the behavior, and develop a fix in an expedient fashion. Before
submitting your bug report, please perform a cursory
search to see if the problem has been already reported. If it has been reported, and the
issue is still open, add a comment to the existing issue instead of opening a
new issue.

Tips for effective bug reporting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Use a clear descriptive title for the issue

- Describe the steps to reproduce the problem, the behavior you observed after
  following the steps, and the expected behavior

- Provide the SHA ID of the git commit that you are using

- For runtime errors, provide a function call stack


Submitting pull-requests
^^^^^^^^^^^^^^^^^^^^^^^^

Contributions can take the form of bug fixes, feature enhancements,
documentation updates. All updates to the repository are managed via `pull
requests <https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/proposing-changes-to-your-work-with-pull-requests>`_.
One of the easiest ways to get started is by looking at `open
issues <https://github.com/sandialabs/pyNuMAD/issues>`_ and contributing fixes,
enhancements that address those issues. If your code contribution involves large
changes or additions to the codebase, we recommend opening an issue first and
discussing your proposed changes with the core development team to ensure that
your efforts are well directed, and so that your submission can be reviewed and
merged seamlessly by the maintenance team.

Guidelines for preparing and submitting pull-requests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Use a clear descriptive title for your pull-requests

- Describe if your submission is a bugfix, documentation update, or a feature
  enhancement. Provide a concise description of your proposed changes. 
  
- Provide references to open issues, if applicable, to provide the necessary
  context to understand your pull request
  
- Make sure that your pull-request merges cleanly with the `main` branch of
  pyNuMAD. When working on a feature, always create your feature branch off of
  the latest `main` commit
  
- Ensure that the code compiles without warnings, (leave for later? the unit tests and regression
  tests all pass without errors, and the documentation builds properly with your
  modifications)
  
- New physics models and code enhancements should be accompanied with relevant
  updates to the documentation, supported by necessary verification and
  validation, as well as unit tests and regression tests
  
  
Once a pull-request is submitted you will iterate with pyNuMAD maintainers
until your changes are in an acceptable state and can be merged in. You can push
addditional commits to the branch used to create the pull-request to reflect the
feedback from maintainers and users of the code.


Coding conventions
^^^^^^^^^^^^^^^^^^

These are the conventions we require. Note that these have not necessarily been followed in the past but moving forward the will:

- We indent using four spaces (soft tabs)
- Moving forward we ALWAYS put spaces after list items and method parameters (`[1, 2, 3]`, not `[1,2,3]`), around operators (`x = 1`, not `x=1`), and around hash arrows.
- We use camelCase when naming functions and variables
- This is open source software. Consider the people who will read your code, and make it look nice for them.



Developer Installation
----------------------

To maintain a local installation, developers should use the following commands::
    
    git clone https://github.com/sandialabs/pyNuMAD
    cd pynumad
    pip install -e .


Running Common Developer Tasks
------------------------------

Many common developer tasks have been implemented through nox for convenient and consistent
results. The following subsections describe the various tasks available through nox.

Testing
^^^^^^^
To run tests locally, run::

    nox -s tests

at the root of the repository. Note that this requires the installation
of nox.

Formatting
^^^^^^^^^^
To check if your code complies to the repositories standards run::

  nox -s style


Black

Building Docs
^^^^^^^^^^^^^

To build docs locally, navigate to ``pynumad/docs`` and run::

  nox docs

After building, the static html files can be found in ``_build/html``.

Serving Docs
^^^^^^^^^^^^

Documentation
-------------

Conventions
^^^^^^^^^^^

Docstrings
^^^^^^^^^^

The documentation for pyNuMAD adheres to NumPy style docstrings. Not only does this
help to keep a consistent style, but it is also necessary for the API documentation
to be parsed and displayed correctly. For an example of what this should look like::

    def func(arg1, arg2):
    """Summary line.

    Extended description of function.

    Parameters
    ----------
    arg1 : int
        Description of arg1
    arg2 : str
        Description of arg2

    Returns
    -------
    bool
        Description of return value

    """
    return True

Additional examples can be found in the 
`napoleon documentation <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`_.
The following boilerplate can be copy-pasted into the top of a function definition
to help get things started::

    """Summary line.

    Extended description of function.

    Parameters
    ----------

    Returns
    -------


    """
