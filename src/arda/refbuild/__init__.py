"""Offline reference-database construction (Phase 1).

Downloads IMGT germline references, enumerates in-frame V(D)J rearrangements,
runs IgBLAST to extract region markup, translates to protein, and writes the
curated ``database/vdj/<species>/`` artifacts.
"""
