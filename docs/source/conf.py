# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SHARP Pack'
version = '2.0'
copyright = '2025, SHARP Pack Developers'
author = 'Dil Limbu'
release = '2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        'sphinx_rtd_theme',
        'sphinxcontrib.bibtex'
        ]

templates_path = ['_templates']
bibtex_bibfiles = ['ref.bib']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'alabaster'
html_show_sourcelink = False
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_favicon = '_static/rpsh.png'
