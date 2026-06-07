Installation
============

``arda`` uses a dedicated conda environment for the MMseqs2 binary and the C++
toolchain; the package itself installs with pip and builds a small C++ extension.

Bootstrap
---------

.. code-block:: bash

   bash setup.sh
   conda activate arda

``setup.sh`` flags:

* ``--no-conda`` — use the already-active environment instead of creating ``arda``.
* ``--build-db`` — rebuild the reference database after install (needs IgBLAST).
* ``--tests`` — run the fast unit + synthetic suites.

What gets installed
-------------------

* The ``arda`` conda env (Python, ``mmseqs2``, a C++ compiler, perl).
* The latest IgBLAST release into ``bin/`` (gitignored) — only needed to rebuild
  references, not at annotation time.
* The ``arda`` package + the ``arda._markup`` C++ extension (editable install).

The committed ``database/vdj/<organism>/`` references mean most users do not need
to build anything.
