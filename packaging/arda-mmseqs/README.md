# arda-mmseqs

A static [MMseqs2](https://github.com/soedinglab/MMseqs2) binary, packaged as a wheel.

Install it and `arda` finds mmseqs with nothing else to set up:

```bash
pip install 'arda-mapper[mmseqs]'
```

Plain `pip install arda-mapper` still works — it downloads the same binary on first use.
This package only removes that first-run download (and the need for a network at all).

## Why a separate distribution

pip cannot be asked to *prefer* one wheel of a project over another: build tags are ranked,
not chosen. So "arda with a bundled binary" and "arda without one" cannot be two variants of
`arda-mapper`; the binary has to live in its own distribution that an extra depends on. This
is the shape `cmake`, `ninja`, `ruff` and `patchelf` all use on PyPI.

## Why arda still checks the binary

An mmseqs index is only reusable by the release that compiled it. `arda.mmseqs.mmseqs_binary()`
therefore version-matches every candidate against the precompiled indexes in `database/`, and
this wheel is just the first candidate it tries — a stale bundle gets rejected like any other.
That is why the version pin on this package does not have to be exact to be safe.

## Building

```bash
python build_wheel.py                                                   # this platform
python build_wheel.py --plat manylinux_2_17_x86_64 \
                      --asset mmseqs-linux-avx2.tar.gz                  # cross-build
```

One wheel per platform. The asset table is arda's own (`arda._mmseqs_fetch.default_asset`), so
the bundled build and the auto-fetched build can never disagree.

## Licence

MMseqs2 is MIT-licensed (Steinegger & Söding); the upstream notice ships with the wheel as
`LICENSE.mmseqs2`, as MIT requires. This packaging is MIT too.
