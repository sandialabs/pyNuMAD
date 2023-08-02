# -*- coding: utf-8 -*-
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# -- Project information -----------------------------------------------------
project = 'PyNuMAD'
copyright = '2023, Kirk Bonney'
author = 'Kirk Bonney'
version = '1.0'
release = '1.0'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    # 'sphinx.ext.viewcode', # commenting out for now b/c bad render width
    'sphinx.ext.napoleon',
    'sphinxcontrib.bibtex',
]
napoleon_use_rtype = False
viewcode_import = True
numpydoc_show_class_members = True
numpydoc_show_inherited_class_members = False
numpydoc_class_members_toctree = False
autodoc_member_order = 'bysource'
autoclass_content = 'both'
bibtex_bibfiles = ['refs/publications.bib','refs/conclusion.bib']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
language = "en"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '_user']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_style = 'css/my_style.css'

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'PyNuMADdoc'