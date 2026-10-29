"""Major-only derivatives of immutable references; query evidence is never cached."""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import tempfile

from ..alleles import major_allele_map
from .._locking import build_lock
from ..paths import data_dir
from ..refbuild.imgt import read_fasta


def major_reference(ref):
    """Restrict searchable templates to major alleles, retaining gene ambiguity."""
    from .reference import Reference, _load_markup
    from ..refbuild.segments import build_segment_reference

    base = ref.target_fasta.parent
    selected = major_allele_map(ref.organism)
    keep = set(selected.values())
    names = ["markup.tsv", "markup.aa.tsv", "alleles.fasta", "alleles.aa.fasta", "combinations.tsv"]
    digest = hashlib.sha256(b"major-reference-v1" + json.dumps(selected, sort_keys=True).encode())
    for name in names:
        digest.update((base / name).read_bytes())
    out = data_dir() / "major_references" / digest.hexdigest()
    out.parent.mkdir(parents=True, exist_ok=True)
    with build_lock(out.parent / (out.name + ".lock"), done=lambda: (out / "READY").exists()) as ours:
        if ours:
            with tempfile.TemporaryDirectory(dir=out.parent) as tmp:
                tmp = Path(tmp)
                ids = set()
                for name in ("markup.tsv", "markup.aa.tsv", "combinations.tsv"):
                    with (base / name).open() as handle, (tmp / name).open("w") as target:
                        reader = csv.DictReader(handle, delimiter="\t")
                        writer = csv.DictWriter(target, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
                        writer.writeheader()
                        for row in reader:
                            call_keys = ("v_calls", "j_calls") if name == "combinations.tsv" else ("v_call", "j_call")
                            eligible = True
                            for key in call_keys:
                                old = row.get(key, "")
                                row[key] = ",".join(a for a in old.split(",") if a in keep)
                                eligible &= not old or bool(row[key])
                            if eligible:
                                writer.writerow(row)
                                ids.add(row["scaffold_id"])
                for name in ("alleles.fasta", "alleles.aa.fasta"):
                    with (tmp / name).open("w") as handle:
                        for header, sequence in read_fasta(base / name):
                            if header in ids:
                                handle.write(f">{header}\n{sequence}\n")
                build_segment_reference(ref.organism, out_dir=tmp)
                (tmp / "READY").write_text("major-only reference\n")
                os.rename(tmp, out)
    entries = {}
    markup = "markup.aa.tsv" if ref.seqtype == "aa" else "markup.tsv"
    _load_markup(out / markup, entries)
    if ref.seqtype == "nt":
        _load_markup(out / "segments.markup.tsv", entries)
    result = Reference(ref.organism, ref.seqtype, out / ref.target_fasta.name, entries,
                       {locus: [(a, s) for a, s in alleles if a in keep] for locus, alleles in ref.d_germlines.items()},
                       {key: anchor for key, anchor in ref.anchors.items() if key[1] in keep},
                       major_alleles=selected)
    return result
