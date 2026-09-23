"""Sphinx configuration for arda."""

import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath("../src"))

project = "arda"
author = "Mikhail Shugay"
copyright = "2026, Mikhail Shugay"
# Read from the package source rather than a literal: this said 2.10.0 while the package was
# 2.22.0, i.e. twelve releases of docs that named the wrong version on every page. Parsed, not
# imported -- the docs job installs sphinx and polars, not arda, so `import arda` would fail.
release = re.search(
    r'__version__ = "([^"]+)"',
    (Path(__file__).parent.parent / "src/arda/__init__.py").read_text(),
).group(1)
version = release

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
autodoc_mock_imports = ["typer", "requests", "arda._markup"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Served as a project page under the org custom domain (antigenomics.github.io
# -> docs.isalgo.dev), so arda's docs live at docs.isalgo.dev/arda/.
html_baseurl = "https://docs.isalgo.dev/arda/"
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
# The version is in the brand, on every page: which release a page documents is the first thing a
# reader needs and the last thing anyone remembers to write down.
html_title = f"arda {release}"
html_theme_options = {
    "github_url": "https://github.com/antigenomics/arda",
    "show_prev_next": False,
    # The sidebar was one flat toctree in source order -- fifteen full page titles, no grouping,
    # nothing marking where the reader is. `index.rst` now carries five captioned toctrees; these
    # options make the theme render that structure instead of flattening it.
    "show_nav_level": 2,        # open each caption's pages, do not collapse to the caption alone
    "navigation_depth": 3,      # let a page's own sections show under it
    "collapse_navigation": False,
    "header_links_before_dropdown": 4,  # Installation, Usage, Samples, Recipes; rest under More
    "show_toc_level": 2,        # the right-hand on-page TOC: subsections too, not just top level
}

# Left sidebar: the whole captioned tree on every page (see `_templates/site-nav.html`). The
# stock `sidebar-nav-bs` renders only the CHILDREN of the current top-level page, and every page
# here is a top-level sibling -- so it emitted a bare "Section Navigation" heading with nothing
# under it.
html_sidebars = {
    "**": ["site-nav"],
    "index": [],  # the landing page has its own card grid; a duplicate tree adds nothing
}
