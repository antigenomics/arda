"""Force a platform-tagged, non-purelib wheel.

The package is pure Python *source*, but it carries a compiled executable, so a
``py3-none-any`` wheel would be installed on every platform and be wrong on all but one.
setuptools decides purity from the presence of extension modules, and there are none here —
hence the explicit override. ``build_wheel.py`` supplies the tag via ``$ARDA_MMSEQS_PLAT``.
"""

from setuptools import setup
from setuptools.dist import Distribution


class BinaryDistribution(Distribution):
    def has_ext_modules(self) -> bool:  # noqa: D102 — setuptools hook
        return True

    def is_pure(self) -> bool:  # noqa: D102 — setuptools hook
        return False


setup(distclass=BinaryDistribution)
