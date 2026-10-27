"""Measure mapping batches, wall time and sampled process-tree resources.

Run from the repository root with its Python environment. Uses the committed
paired fixture unless --r1/--r2 are supplied. Fresh mapping outputs are temporary;
only the JSON measurements survive. This measures mapping, not the whole pipeline.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time


WORKER = r'''
import json, os, sys
os.environ["POLARS_MAX_THREADS"] = "1"
from arda import mmseqs
original = mmseqs.run
trace = sys.argv.pop(1)
def run(args, **kw):
    row = {"command": args[0]}
    if "--threads" in args:
        row["threads"] = int(args[args.index("--threads") + 1])
    if args[0] == "createdb":
        with open(args[1]) as fh:
            row["queries"] = sum(line.startswith(">") for line in fh)
    with open(trace, "a") as fh:
        fh.write(json.dumps(row) + "\n")
    return original(args, **kw)
mmseqs.run = run
from arda.cli import main
main()
'''


def sample_tree(pid):
    snapshot = subprocess.run(["ps", "-axo", "pid=,ppid=,rss=,comm="],
                              capture_output=True, text=True, check=True, timeout=10)
    rows = [line.split(None, 3) for line in snapshot.stdout.splitlines()]
    rows = {int(r[0]): (int(r[1]), int(r[2]), r[3]) for r in rows if len(r) == 4}
    family = {pid}
    while True:
        found = {p for p, (parent, _, _) in rows.items() if parent in family}
        if found <= family:
            break
        family |= found
    family &= rows.keys()
    if sys.platform == "darwin":
        result = subprocess.run(["ps", "-M", "-p", ",".join(map(str, family))],
                                capture_output=True, text=True, timeout=10)
        threads = max(0, len(result.stdout.splitlines()) - 1)
    else:
        threads = 0
        for p in family:
            try:
                threads += len(os.listdir(f"/proc/{p}/task"))
            except FileNotFoundError:
                pass
    return {"rss_kib": sum(rows[p][1] for p in family), "live_threads": threads,
            "mmseqs_processes": sum(Path(rows[p][2]).name == "mmseqs" for p in family)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1", type=Path, default=Path("tests/data/rnaseq_real/reads_1.fq.gz"))
    ap.add_argument("--r2", type=Path, default=Path("tests/data/rnaseq_real/reads_2.fq.gz"))
    ap.add_argument("--mode", choices=("bulk", "amplicon"), default="bulk")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--timeout", type=float, default=180)
    ap.add_argument("--limit", type=int, help="Maximum input fragments per run (ARDA's existing limit).")
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if args.repeats < 1:
        ap.error("--repeats must be positive")
    from arda import __version__
    from arda._threads import resolve_threads

    measurements, digest = [], None
    for budget in dict.fromkeys((1, min(2, resolve_threads()), 0)):
        for rep in range(args.repeats):
            with tempfile.TemporaryDirectory(prefix="arda_parallel_") as directory:
                root = Path(directory)
                out, report, trace = (root / name for name in ("out.tsv", "report.json", "trace.jsonl"))
                cmd = [sys.executable, "-c", WORKER, str(trace), "map", "--r1", str(args.r1),
                       "--r2", str(args.r2), "-o", str(out), "--report", str(report),
                       "--threads", str(budget)]
                cmd += (["--prefilter"] if args.mode == "bulk" else
                        ["--two-pass", "--fast-segments", "--v-only-on-segment"])
                if args.limit is not None:
                    cmd += ["--limit", str(args.limit)]
                peak = dict(rss_kib=0, live_threads=0, mmseqs_processes=0)
                start = time.perf_counter()
                with (root / "run.log").open("w+") as log:
                    proc = subprocess.Popen(cmd, stdout=log, stderr=log, start_new_session=True)
                    try:
                        while proc.poll() is None:
                            if time.perf_counter() - start > args.timeout:
                                raise TimeoutError(f"Mapping exceeded {args.timeout}s")
                            for key, value in sample_tree(proc.pid).items():
                                peak[key] = max(peak[key], value)
                            time.sleep(0.1)
                        if proc.returncode:
                            log.seek(0)
                            raise RuntimeError(log.read())
                    finally:
                        if proc.poll() is None:
                            import signal

                            os.killpg(proc.pid, signal.SIGTERM)
                            proc.wait(timeout=15)
                current = hashlib.sha256(out.read_bytes()).hexdigest()
                if digest is not None and current != digest:
                    raise AssertionError("AIRR output changed across thread counts/repetitions")
                digest = current
                result = json.loads(report.read_text())
                measurements.append(dict(requested_threads=budget, threads=result["threads"],
                                         repetition=rep + 1, wall_seconds=time.perf_counter() - start,
                                         total_reads=result["total_reads"], mapped_reads=result["mapped_reads"],
                                         digest=current, peak_sampled=peak,
                                         commands=[json.loads(s) for s in trace.read_text().splitlines()]))
                args.output.write_text(json.dumps(dict(arda=__version__, mode=args.mode,
                    sampling_interval_seconds=0.1, measurements=measurements), indent=2) + "\n")
                print(json.dumps({k: v for k, v in measurements[-1].items() if k != "commands"}), flush=True)


if __name__ == "__main__":
    main()
