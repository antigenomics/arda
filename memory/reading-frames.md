# Reading frames: the subtle part of scaffold assembly

## V frame is NOT always 0
We assumed IMGT V-REGION starts at codon 1 (frame 0). **False for partial-5′
alleles** — e.g. human `IGHV5-51*05` (245 nt, "partial in 5′ and in 3′") reads
cleanly only in frame 2. The IMGT header field 8 ("codon start") is **always 1**
and useless here; the real signal is the `partial in 5'` flag (field ~14).

## Fix: auto-detect coding frame
Rather than parse fragile IMGT fields, `translate.detect_coding_frame` picks the
frame (0/1/2) with the fewest stop codons. For a real germline V exactly one
frame is stop-free and it is unique (verified across alleles). We trim the
leading `frame` nt so every V — and thus the whole scaffold — reads in frame 0.

## J frame comes from the aux file
`bin/optional_file/<organism>_gl.aux` column 2 = 0-based "first coding frame start
position". We pad the V–J junction with N so `(len(V_trimmed) + n_pad + jframe) %
3 == 0`. This put the canonical FR4 `WGQGTLVTVSS` exactly at the J translation —
the sanity check that the frame logic is correct.

## Validation signal
After the fix, ~18/20 human IGH smoke scaffolds were productive with correct
canonical FR4. A few J alleles still yield non-productive/incomplete markup; those
scaffolds are dropped from the reference and counted in `build.log`.
