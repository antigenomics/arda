"""Runtime mapping & markup transfer (Phase 2).

Maps input sequences against the curated reference DB with MMseqs2 and projects
reference region markup onto each query via the compiled ``arda._markup`` hot path.

TODO (D-segment mapping): scaffolds are V*J only, so FR/CDR coordinates and
v_call/j_call are assigned but ``d_call`` is not. To add it, after V/J transfer,
align the CDR3 interior (between the projected CDR3 start and end) against a D
germline DB and emit ``d_call`` + ``d_sequence_{start,end}``. Must handle **double
D-D junctions** (D-D fusions in IGH/TRD) — emit a second ``d2_call`` when two D
segments are present. See memory/markup-transfer.md.
"""
