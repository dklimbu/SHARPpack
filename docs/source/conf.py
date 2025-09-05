# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SHARP Pack'
version = '2.0'
copyright = '2025, SHARP Pack Developers'
author = ''
release = version

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
numfig = True
html_show_sourcelink = False
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_favicon = '_static/figures/rpsh.png'


# -- Options for MathJax -------------------------------------------------
mathjax3_config = {
    "tex": {
        "macros": {
            "angstrom": r"\unicode{8491}",   # Å (Angström symbol, U+212B)
        }
    }
}

# For LaTeX/PDF build
latex_elements = {
    'preamble': r'''
        \newcommand{\angstrom}{\text{\AA}}  % Use LaTeX's Å symbol
    ''',
    'extraclassoptions': 'openany,oneside',
    'releasename': "Version",
}

#for latexpdf to remove blank page
#latex_elements = {
#  'extraclassoptions': 'openany,oneside'
#}


