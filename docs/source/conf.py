# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
# sys.path.insert(0, os.path.abspath('../../pdb2sql'))

# -- Project information -----------------------------------------------------

project = 'pdb2sql'
copyright = '2019, Netherlands eScience Center'
author = 'Nicolas Renaud'

# The full version, including alpha/beta/rc tags
release = '0.2.2'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosectionlabel',
    'IPython.sphinxext.ipython_directive',
    'IPython.sphinxext.ipython_console_highlighting',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'default'
# html_theme = 'sphinx_rtd_theme'
html_theme = 'nature'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- Options for extensions -------------------------------------------------

# Only the __init__ method’s docstring is inserted.
autoclass_content = 'init'
# order members by source code order
autodoc_member_order = 'bysource'
# Disable docstring inheritance
autodoc_inherit_docstrings = False
# mock the packges that is not avaiable on your machine
# autodoc_mock_imports = ['cython', 'sqlalchemy', 'matplotlib',
#                         'numpy', 'schema', 'tqdm', 'pandas']

# napoleon
napoleon_numpy_docstring = False
napoleon_use_rtype = False

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# intersphinx
intersphinx_mapping = {
    'python': ('https://docs.python.org/', None),
    'numpy': ('http://docs.scipy.org/doc/numpy/', None)
}

## autosummary
# Make _autosummary files and include them
autosummary_generate = True
# autosummary_imported_members = True


## ipython
ipython_warning_is_error = False
ipython_execlines = [
    "import numpy as np",
    'pd.options.display.encoding="utf8"',
]


header = """\

.. ipython:: python
   :suppress:

   import os
   os.chdir(r'{}')
""".format(
    os.path.dirname(os.path.dirname(__file__))
)