# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('../lib'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'GELATO'
copyright = '2025, Interstellar Technologies Inc.'
author = 'Interstellar Technologies Inc.'
release = '0.8.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
    'sphinx.ext.intersphinx',
    'breathe',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Extension configuration -------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
autodoc_member_order = 'bysource'

# Mock imports for modules that can't be imported during doc build
autodoc_mock_imports = [
    'user_constraints',
    'dynamics_c',
    'coordinate_c',
    'utils_c',
    'USStandardAtmosphere_c',
    'IIP_c',
]

# Suppress warnings for duplicate object descriptions
suppress_warnings = ['autodoc.import_object', 'ref.duplicate']

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
}

# Breathe Configuration for C++ Documentation
breathe_projects = {
    "GELATO_CPP": "_build/doxygen/xml"
}
breathe_default_project = "GELATO_CPP"

# -- GitHub Pages configuration ----------------------------------------------
# Base URL for GitHub Pages (will be set by sphinx.ext.githubpages)
html_baseurl = 'https://istellartech.github.io/GELATO/'
html_context = {
    "display_github": True,
    "github_user": "istellartech",
    "github_repo": "GELATO",
    "github_version": "master",
    "conf_py_path": "/docs/",
}

# Additional HTML theme options
html_theme_options = {
    'navigation_depth': 4,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'includehidden': True,
    'titles_only': False
}
