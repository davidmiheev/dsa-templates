# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
import inspect
import sphinx_copybutton

sys.path.append('../..')

import stack
import all_subs, bin_search, trie, bfs, dfs, dijkstra
import dynam_prog, backtrack, union_find, fenwick, min_span_tree
import top_sort

project = 'dsa-templates'
copyright = '2024, David Mikheev'
author = 'David Mikheev'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx_copybutton',
    'sphinx_rtd_theme'
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

add_module_names = False
add_function_parentheses = True

# Enable Google style docstrings
napoleon_google_docstring = True

# Enable NumPy style docstrings
napoleon_numpy_docstring = True

# Set to True to use Sphinx’s native docstring section style
napoleon_use_param = True
napoleon_use_rtype = True

# Set to True to create a list of class methods in the order they appear in the source code
napoleon_use_ivar = False


latex_engine = 'xelatex'
latex_elements = {
    'preamble': r"""
\usepackage{graphicx}
\usepackage{amssymb, amsmath, amsthm, mathrsfs, mathbbol}
\usepackage{fontspec, mathspec}
"""
}
