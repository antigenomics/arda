"""Sphinx configuration for arda."""

import os
import sys

sys.path.insert(0, os.path.abspath("../src"))

project = "arda"
author = "Mikhail Shugay"
copyright = "2026, Mikhail Shugay"
release = "2.0.1"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
]

# Google/NumPy docstrings via napoleon.
napoleon_google_docstring = True
napoleon_numpy_docstring = False

autosummary_generate = False
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {"members": True, "undoc-members": True, "show-inheritance": True}
# Heavy/optional deps that need not import at doc-build time.
autodoc_mock_imports = ["polars", "typer", "requests", "arda._markup"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Served as a project page under the org custom domain (antigenomics.github.io
# -> docs.isalgo.dev), so arda's docs live at docs.isalgo.dev/arda/.
html_baseurl = "https://docs.isalgo.dev/arda/"
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "github_url": "https://github.com/antigenomics/arda",
    "show_nav_level": 2,
}
